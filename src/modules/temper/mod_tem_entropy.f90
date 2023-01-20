!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



module mod_tem_entropy

  use def_kintyp,         only : ip,rp
  use def_parame
  use def_master
  use def_domain
  use def_temper
  use mod_gradie
  use mod_ker_proper 
  use mod_solver,         only : solver_lumped_mass_system
  use mod_solver,         only : solver_explicit
  use mod_memory,         only : memory_alloca
  use mod_memory,         only : memory_deallo
  use mod_communications, only : PAR_MAX
  use mod_elmgeo

  implicit none

  real(rp),   pointer :: bt(:)         ! Energy RHS
  real(rp),   pointer :: tt(:)         ! Temperature
  real(rp),   pointer :: tt0(:)        ! Old temperature 
  real(rp),   pointer :: rt(:,:)         ! Energy residual
  !real(rp),   pointer :: rt(:)         ! Energy residual
  real(rp),   pointer :: Mass(:)       ! lumped mass
  real(rp),   pointer :: Tim(:)        ! projection dt/(rho*cp)  or dt/rho
  real(rp),   pointer :: entro(:,:)    ! entropy
  real(rp),   pointer :: grads(:,:)    ! gradient of entropy
  real(rp),   pointer :: wmeanl(:,:)   ! local wmean 
  type(r1p),  pointer :: envit_gp(:)   ! save entropy viscosity for postprocessing
  real(rp)            :: dt
!  real(rp)            :: eminmax(2)
  real (rp)           :: a(5)
  real (rp)           :: b(5,5)
  real (rp)           :: sigma(5)
  type(soltyp), pointer :: solve1(:)
  type(soltyp), pointer :: solve2(:)

  private

  public :: tem_entropy_solution
  public :: tem_entropy_viscosity
  public :: tem_entropy_postprocess
  public :: tem_entropy_initialization
  public :: tem_entropy_memory

contains 

  subroutine tem_entropy_initialization()

    nullify(rt)
    nullify(Tim)
    nullify(tt)
    nullify(tt0)
    nullify(entro)
    nullify(grads)
    nullify(envit_gp)
    
  end subroutine tem_entropy_initialization
  
  subroutine tem_entropy_memory()

    integer(ip) :: ipoin,ielem,pelty,pgaus

    if(kfl_coupl(ID_NASTAL,ID_CHEMIC)==1) then
       wmeanl => wmean
    else
       call memory_alloca(mem_modul(1:2,modul),'WMEANL' ,'tem_entropy_allocate',wmeanl, max(1_ip,npoin), 2_ip)
    endif

    call memory_alloca(mem_modul(1:2,modul),'TT'      ,'tem_entropy_allocate', tt    ,  max(1_ip,npoin))
    call memory_alloca(mem_modul(1:2,modul),'TT0'     ,'tem_entropy_allocate', tt0   ,  max(1_ip,npoin))
    call memory_alloca(mem_modul(1:2,modul),'RT'      ,'tem_entropy_allocate', rt    ,  max(1_ip,npoin),5_ip)
    !call memory_alloca(mem_modul(1:2,modul),'RT'      ,'tem_entropy_allocate', rt    ,  max(1_ip,npoin))
    call memory_alloca(mem_modul(1:2,modul),'TIM'     ,'tem_entropy_allocate', Tim   ,  max(1_ip,npoin))
    call memory_alloca(mem_modul(1:2,modul),'ENTRO'   ,'tem_entropy_allocate', entro ,  max(1_ip,npoin),  2_ip)
    call memory_alloca(mem_modul(1:2,modul),'GRADS'   ,'tem_entropy_allocate', grads ,  ndime          ,  max(1_ip,npoin))

    call memory_alloca(mem_modul(1:2,modul),'ENVIT_GP','tem_entropy_allocate', envit_gp, max(1_ip,nelem))
    do ielem = 1,nelem
       pelty = ltype(ielem)
       pgaus = ngaus(pelty)
       call memory_alloca(mem_modul(1:2,modul),'ENVIT_GP % A','tem_entropy_allocate', envit_gp(ielem)%a,pgaus)
    end do

    if(kfl_coupl(ID_NASTAL,ID_CHEMIC)==0) then
       do ipoin=1,npoin
          wmeanl(ipoin,:) = 1.0_rp
       enddo
    endif

    dt = 1.0_rp
    
    solve1 => solve(1:)
    solve2 => solve(2:)

    !
    ! RK2 table
    !
    !!a(4) = 1.0_rp

    !!b(4,3) = 1.0_rp

    !!sigma(3) = 0.5_rp
    !!sigma(4) = 0.5_rp

    !
    ! RK4 table
    !
    a(2) = 0.5_rp
    a(3) = 0.5_rp
    a(4) = 1.0_rp

    b(2,1) = 0.5_rp 
    b(3,2) = 0.5_rp
    b(4,3) = 1.0_rp

    sigma(1) = 1.0_rp/6.0_rp 
    sigma(2) = 1.0_rp/3.0_rp 
    sigma(3) = 1.0_rp/3.0_rp 
    sigma(4) = 1.0_rp/6.0_rp 

  end subroutine tem_entropy_memory

  subroutine tem_entropy_eval(istep)

     integer(ip), intent(in) :: istep
     integer(ip)             :: ipoin
     integer(ip)             :: KFL_ENTROPY_TEM_LOCAL

     !
     ! Get time steps
     !
     dt = 1.0_rp/dtinv
     do ipoin = 1,npoin
        dt_rho_cp_tem(ipoin) = 0.0_rp
     end do
     !
     ! Get RHS 
     !
     KFL_ENTROPY_TEM_LOCAL = kfl_entropy_tem
     kfl_entropy_tem = 0_ip
     call tem_matrix()
     kfl_entropy_tem = KFL_ENTROPY_TEM_LOCAL

     if ( INOTEMPTY ) then
        if ( kfl_rhs_scal_tem == 0 ) then
           call solver_lumped_mass_system(1_ip,dt_rho_cp_tem)
        end if

        if ( kfl_rhs_scal_tem == 0 ) then
           do ipoin = 1,npoin
              Tim(ipoin) = dt_rho_cp_tem(ipoin) / dt
              rt(ipoin,istep) = bt(ipoin)
           end do
        else if ( kfl_rhs_scal_tem > 0 ) then
           do ipoin = 1,npoin
              Tim(ipoin) = 1.0_rp
              rt(ipoin,istep) = bt(ipoin)
           end do
        end if
     end if

     !
     ! Exchange RHS:
     !
     call rhsmod(1_ip,rt)

     if(istep == 4_ip) then
        call tem_entropy_solution_sij(4_ip,sigma,1.0_rp)
     else
        call tem_entropy_solution_sij(istep,b(istep+1,:),a(istep+1))
     endif

  end subroutine tem_entropy_eval

  subroutine tem_entropy_solution_sij(istep,weight1,weight2)

     integer(ip) , intent(in) :: istep
     real(rp) ,    intent(in) :: weight1(5)
     real(rp) ,    intent(in) :: weight2
     integer(ip)              :: ipoin,jstep
     real(rp)                 :: aux2

     do ipoin = 1,npoin
        aux2 = 0.0_rp
        do jstep=1,istep
           aux2   = aux2 + rt(ipoin,jstep)*weight1(jstep)
        end do
        tt(ipoin) = Tim(ipoin)*aux2*dt
     end do
     !
     ! Solve system
     !
     call solver_explicit(solve1,tt,EXCHANGE=.false.,solve_consistent=solve2)
     do ipoin = 1,npoin
        tt(ipoin) = tt0(ipoin) + tt(ipoin)
     end do
     !
     ! Dirichet bbcc 
     !
     do ipoin = 1,npoin
        if( kfl_fixno_tem(1,ipoin) > 0 ) &
             tt(ipoin) = bvess_tem(1,ipoin,1)
     end do

  end subroutine tem_entropy_solution_sij

  subroutine tem_entropy_solution()

    integer(ip) :: ipoin
    real(rp)    :: tt_max
!    real(rp)    :: entr_max, entr_min
    logical(rp), parameter :: vel_square=.false.

    !-----------------------------------------------------------------
    !
    ! Allocate memory if necessary
    !
    !-----------------------------------------------------------------

    bt   => rhsid
    Mass => vmass

    tt_max = 0.0_rp

    do ipoin=1,npoin
       tt0(ipoin)   = unkno(ipoin) ! current unknown 
       tt(ipoin)    = unkno(ipoin) ! prediction of unknown
    enddo 

    !
    ! Runge Kutta stage
    !
    if (kfl_disable_entpre == 0) then
       call tem_entropy_eval(1_ip)
       call tem_endite(ITASK_INNITE)
       call tem_entropy_eval(2_ip)
       call tem_endite(ITASK_INNITE)
       call tem_entropy_eval(3_ip)
       call tem_endite(ITASK_INNITE)
       call tem_entropy_eval(4_ip)
    end if

    do ipoin=1,npoin
       tt(ipoin)       = unkno(ipoin)
       unkno(ipoin)    = tt0(ipoin)
       therm(ipoin,1)  = tt0(ipoin)
    enddo

    do ipoin = 1,npoin
       tt_max = max(abs(tt0(ipoin)),tt_max) 
    end do
    call PAR_MAX(tt_max)

    !entr_max = 0.0_rp
    !entr_min = 1.0e10_rp

    if (.not. vel_square) then

       !
       ! Entropy solution
       !
       do ipoin = 1,npoin
          entro(ipoin, 2) = cosh(tt0(ipoin)/tt_max)  ! Entropy values at current step
          entro(ipoin, 1) = cosh(tt(ipoin)/tt_max)   ! Entropy values at prediction step
          !entro(ipoin, 2) = tt0(ipoin)**2 
          !entro(ipoin, 1) = tt(ipoin)**2  
       end do
       
    else
       !
       ! Use kinetic energy
       !
       do ipoin = 1,npoin
          tt(ipoin) = 0.5_rp*dot_product(advec(1:ndime,ipoin,1),advec(1:ndime,ipoin,1)) 
       end do
       
       do ipoin = 1,npoin
          entro(ipoin, 2) = tt0(ipoin)
          entro(ipoin, 1) = tt(ipoin)
       end do
       
    end if

    !do ipoin = 1,npoin
    !   entr_max = max(abs(entro(ipoin, 1)),entr_max) 
    !   entr_min = min(abs(entro(ipoin, 1)),entr_min) 
    !end do

    !eminmax(1) =  entr_max
    !eminmax(2) = -entr_min
    !call PAR_MAX(2_ip,eminmax)
    !entr_max = eminmax(1)
    !entr_min = -eminmax(2)

    !do ipoin = 1,npoin
    !   entro(ipoin, 2) = entro(ipoin,2)/abs(entr_max - entr_min+1.0e-14_rp) 
    !   entro(ipoin, 1) = entro(ipoin,1)/abs(entr_max - entr_min+1.0e-14_rp)  
    !end do

    !
    ! Compute gradients
    !
    if( INOTEMPTY ) call gradie(entro(:,2),grads)

    if (vel_square) then
       !
       ! Save old step:
       !
       do ipoin = 1,npoin
          tt0(ipoin) = tt(ipoin)
       end do
       
    end if
    
  end subroutine tem_entropy_solution

  !!subroutine tem_entropy_solution()
  !!  integer(ip) :: ipoin
  !!  integer(ip) :: dummi    
  !!  real(rp)    :: rho(2), cp(2), garho,entr_max, entr_min
  !!  logical(rp), parameter :: vel_square=.false.
  !!  integer(ip)            :: KFL_ENTROPY_TEM_LOCAL

  !!  !-----------------------------------------------------------------
  !!  !
  !!  ! Allocate memory if necessary
  !!  !
  !!  !-----------------------------------------------------------------

  !!  bt   => rhsid
  !!  Mass => vmass

  !!  entr_max = 0.0_rp
  !!  entr_min = 1.0e10_rp

  !!  if (.not. vel_square) then
  !!     !
  !!     ! Copy unknown
  !!     !
  !!     do ipoin=1,npoin
  !!        tt0(ipoin)   = unkno(ipoin) ! current unknown 
  !!        tt(ipoin)    = unkno(ipoin) ! prediction of unknown
  !!     enddo

  !!     !
  !!     ! Get time steps
  !!     !
  !!     dt = 1.0_rp/dtinv
  !!     do ipoin = 1,npoin
  !!        dt_rho_cp_tem(ipoin) = 0.0_rp
  !!     end do
  !!     !
  !!     ! Get RHS 
  !!     !
  !!     KFL_ENTROPY_TEM_LOCAL = kfl_entropy_tem
  !!     kfl_entropy_tem = 0_ip
  !!     call tem_matrix()
  !!     kfl_entropy_tem = KFL_ENTROPY_TEM_LOCAL

  !!     call solver_lumped_mass_system(1_ip,dt_rho_cp_tem)

  !!     do ipoin = 1,npoin
  !!        Tim(ipoin) = dt_rho_cp_tem(ipoin) / dt
  !!        rt(ipoin) = bt(ipoin)
  !!     end do
  !!     !
  !!     ! Exchange RHS:
  !!     !
  !!     call rhsmod(1_ip,rt)
  !!     !
  !!     ! Do an Euler forward step:
  !!     !
  !!     do ipoin = 1,npoin
  !!        tt(ipoin) = tt(ipoin) + ( rt(ipoin)*Tim(ipoin) )* dt / Mass(ipoin)
  !!     end do

  !!     !
  !!     ! Dirichet bbcc 
  !!     !
  !!     do ipoin = 1,npoin
  !!        if( kfl_fixno_tem(1,ipoin) > 0 ) &
  !!             tt(ipoin) = bvess_tem(1,ipoin,1)
  !!     end do

  !!     !
  !!     ! Entropy solution
  !!     !
  !!     do ipoin = 1,npoin
  !!        !! call ker_proper('DENSI','IPOIN',ipoin,dummi,rho)
  !!        !! call ker_proper('SPHEA','IPOIN',ipoin,dummi,cp)
  !!        !! !garho = (prthe(1)/gasco) * (1.0_rp/tt(ipoin))
  !!        !! garho = cp(1)*(prthe(1)/gasco) * (wmeanl(ipoin,1)/tt(ipoin))
  !!        !! rho(1) = rho(1)*cp(1)
  !!        !! entro(ipoin, 2) = (rho(1)/0.4_rp)*log(prthe(1)/rho(1)**1.4_rp)  
  !!        !! entro(ipoin, 1) = (garho/0.4_rp)*log(prthe(1)/garho**1.4_rp)  

  !!        entro(ipoin, 2) = tt0(ipoin)**2 
  !!        entro(ipoin, 1) = tt(ipoin)**2  
  !!     end do
  !!     
  !!  else
  !!     !
  !!     ! Use kinetic energy
  !!     !
  !!     do ipoin = 1,npoin
  !!        tt(ipoin) = 0.5_rp*dot_product(advec(1:ndime,ipoin,1),advec(1:ndime,ipoin,1)) 
  !!     end do
  !!     
  !!     do ipoin = 1,npoin
  !!        entro(ipoin, 2) = tt0(ipoin)
  !!        entro(ipoin, 1) = tt(ipoin)
  !!     end do
  !!     
  !!  end if

  !!  do ipoin = 1,npoin
  !!     entr_max = max(abs(entro(ipoin, 1)),entr_max) 
  !!     entr_min = min(abs(entro(ipoin, 1)),entr_min) 
  !!  end do

  !!  eminmax(1) =  entr_max
  !!  eminmax(2) = -entr_min
  !!  call PAR_MAX(2_ip,eminmax)
  !!  entr_max = eminmax(1)
  !!  entr_min = -eminmax(2)

  !!  do ipoin = 1,npoin
  !!     entro(ipoin, 2) = entro(ipoin,2)/abs(entr_max - entr_min+1.0e-14_rp) 
  !!     entro(ipoin, 1) = entro(ipoin,1)/abs(entr_max - entr_min+1.0e-14_rp)  
  !!  end do
  !!  !
  !!  ! Compute gradients
  !!  !
  !!  if( INOTEMPTY ) call gradie(entro(:,2),grads)

  !!  if (vel_square) then
  !!     !
  !!     ! Save old step:
  !!     !
  !!     do ipoin = 1,npoin
  !!        tt0(ipoin) = tt(ipoin)
  !!     end do
  !!     
  !!  end if
  !!  
  !!end subroutine tem_entropy_solution


 subroutine tem_entropy_viscosity(ielem,pnode,pgaus,igaui,igauf,gpsha,gpcar,elvel,gpden,hleng,gpvol,gpdif) 
    integer(ip), intent(in)    :: ielem,pnode,pgaus,igaui,igauf
    real(rp),    intent(in)    :: gpsha(mnode,mgaus)
    real(rp),    intent(in)    :: gpcar(ndime,mnode,mgaus)
    real(rp),    intent(in)    :: elvel(ndime,mnode)
    real(rp),    intent(in)    :: gpden(mgaus)
    real(rp),    intent(in)    :: hleng(ndime)
    real(rp),    intent(in)    :: gpvol(mgaus)
    real(rp),    intent(inout) :: gpdif(mgaus)

    integer(ip)               :: igaus, inode, ipoin, ielty, p_order
    real(rp)                  :: elSres(mnode)
    real(rp)                  :: betae, ve, h, gpenv

    ielty = ltype(ielem)
    p_order = element_type(ielty)%order

    elSres = 0.0_rp
    do inode=1,pnode
       ipoin = lnods(inode,ielem)
       elSres(inode) = elSres(inode) + (entro(ipoin,1) - entro(ipoin,2))*dtinv
       elSres(inode) = elSres(inode) + dot_product(elvel(1:ndime,inode),grads(1:ndime,ipoin))
    end do


    ve = 0.0_rp
    betae = 0.0_rp
    do inode = 1,pnode
       ipoin = lnods(inode,ielem)
       ve = max(abs(elSres(inode)),ve) 
       betae = max(sqrt(dot_product(elvel(1:ndime,inode),elvel(1:ndime,inode))), betae)
    end do

    h = 1e15_rp
    h = min(hleng(1),h)
    h = min(hleng(2),h)
    if (ndime == 3) then
       h = min(hleng(3),h)
    end if
    h = h/dble(p_order)

    !
    ! Add viscosity
    !
    do igaus=igaui,igauf
       ! gpenv = 1.0_rp*min(0.5_rp*betae*h*gpden(igaus),dtinv*gpden(igaus)*ve*h**2_ip)
       ! gpenv = 0.5_rp*betae*h*gpden(igaus)
       !gpenv = min( 0.5_rp*betae*h*gpden(igaus), 1.0_rp*gpden(igaus)*ve*h**2_ip )
       gpenv = min( 0.5_rp*betae*h*gpden(igaus), 10.0_rp*gpden(igaus)*ve*h**2_ip )
       
       envit_gp(ielem) % a(igaus)   = gpenv
       gpdif(igaus)                 = gpdif(igaus) + gpenv 
    end do

 end subroutine tem_entropy_viscosity


 subroutine tem_entropy_postprocess(envit)

   real(rp),    intent(inout)  :: envit(:)
    integer(ip)                 :: ipoin

    !
    ! Project Gauss point values to nodes
    !
    if (associated(envit_gp)) then
       call smooth (envit_gp, envit)      
    else
       do ipoin=1,npoin
          envit(ipoin) = 0.0_rp
       enddo
    endif
 end subroutine tem_entropy_postprocess


 end module mod_tem_entropy
