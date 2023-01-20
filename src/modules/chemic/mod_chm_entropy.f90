!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



module mod_chm_entropy

  use def_kintyp,         only : ip, rp, r1p
  use def_master,         only : INOTMASTER, ITASK_INNITE, soltyp, dtinv, rhsid, unkno, mem_modul, modul, solve
  use def_domain,         only : ndime, mgaus, mnode, npoin, nelem, lnods, ngaus, vmass, ltype
  use def_chemic,         only : nclas_chm, kfl_disable_entpre_chm, kfl_spray_chm, dt_rho_chm, dt_chm, kfl_fixno_chm, bvess_chm
  use mod_gradie,         only : gradie
  use mod_solver,         only : solver_lumped_mass_system
  use mod_solver,         only : solver_explicit
  use mod_memory,         only : memory_alloca
  use mod_memory,         only : memory_deallo
  use mod_communications, only : PAR_MAX, PAR_MIN

  implicit none

  real(rp),   pointer :: bt(:)        ! Scalars RHS
  real(rp),   pointer :: tt(:)        ! Scalars
  real(rp),   pointer :: tt0(:,:)       ! Scalars at previous step
  !real(rp),   pointee :: rt(:,:)      ! Scalars residual
  real(rp),   pointer :: rt(:,:,:)      ! Scalars residual
  real(rp),   pointer :: aux_tt(:,:)  ! Unknown residuals
  real(rp),   pointer :: Mass(:)      ! lumped mass
  real(rp),   pointer :: Tim(:,:)     ! projection dt / (rho)
  real(rp),   pointer :: entro(:,:,:) ! entropy
  real(rp),   pointer :: grads(:,:,:) ! gradient of entropy
  type(r1p),  pointer :: envic_gp(:,:)! save entropy viscosity for postprocessing
  real(rp)            :: dt
  real(rp)            :: a(5)
  real(rp)            :: b(5,5)
  real(rp)            :: sigma(5)
  type(soltyp), pointer :: solve1(:)
  type(soltyp), pointer :: solve2(:)

  private

  public :: chm_entropy_solution
  public :: chm_entropy_viscosity
  public :: chm_entropy_postprocess
  public :: chm_entropy_initialization
  public :: chm_entropy_memory

contains

   !! NEW STUFF FOR RK !!

   subroutine chm_entropy_initialization()

      nullify(rt)
      nullify(aux_tt)
      nullify(Tim)
      nullify(tt)
      nullify(tt0)
      nullify(entro)
      nullify(grads)
      nullify(envic_gp)

   end subroutine chm_entropy_initialization

   subroutine chm_entropy_memory()

      integer(ip)       :: ielem,pelty,pgaus,iclas

      !
      ! Allocate variables for this module
      !
      call memory_alloca(mem_modul(1:2,modul),'TT'     ,'chm_entropy_allocate',tt   , max(nclas_chm,npoin*nclas_chm))
      call memory_alloca(mem_modul(1:2,modul),'TT0'    ,'chm_entropy_allocate',tt0  , nclas_chm,npoin)
      call memory_alloca(mem_modul(1:2,modul),'RT'     ,'chm_entropy_allocate',rt   , nclas_chm, max(1_ip,npoin),5_ip)
      call memory_alloca(mem_modul(1:2,modul),'AUX_TT' ,'chm_entropy_allocate'      , aux_tt,nclas_chm,npoin)
      call memory_alloca(mem_modul(1:2,modul),'TIM'    ,'chm_entropy_allocate',Tim  , nclas_chm, max(1_ip,npoin))
      call memory_alloca(mem_modul(1:2,modul),'ENTRO'  ,'chm_entropy_allocate',entro, nclas_chm, max(1_ip,npoin), 2_ip)
      call memory_alloca(mem_modul(1:2,modul),'GRADS'  ,'chm_entropy_allocate',grads, nclas_chm, ndime, max(1_ip,npoin))

      call memory_alloca(mem_modul(1:2,modul),'ENVIC_GP','chm_entropy_allocate', envic_gp, nclas_chm, max(1_ip,nelem))
      do ielem = 1,nelem
         pelty = ltype(ielem)
         pgaus = ngaus(pelty)
         do iclas = 1,nclas_chm
            call memory_alloca(mem_modul(1:2,modul),'ENVIC_GP % A','chm_entropy_allocate', envic_gp(iclas,ielem)%a, pgaus)
         enddo
      end do

      solve1 => solve(1:)
      solve2 => solve(2:)

      !tt(1:npoin*nclas_chm) = unkno(1:npoin*nclas_chm)
      !do ii=1,npoin*nclas_chm
      !   tt(ii) = unkno(ii)
      !enddo

      dt = 1.0_rp

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

   end subroutine chm_entropy_memory

   subroutine chm_entropy_eval(iclai,iclaf,istep)

      integer(ip) , intent(in) :: iclai
      integer(ip) , intent(in) :: iclaf
      integer(ip) , intent(in) :: istep
      integer(ip)              :: ipoin,iclass,kpoin,iclaf_gas,iclai_liq

      external                 :: chm_matrix
      external                 :: chm_partis
      external                 :: rhsmod

      do ipoin = 1,npoin
         dt_rho_chm(ipoin) = 0.0_rp
      end do

      if ( kfl_spray_chm /= 0_ip ) then
         do ipoin = 1,npoin
            dt_chm(ipoin) = 0.0_rp
         end do
      end if

      if( INOTMASTER ) call chm_matrix()

      call solver_lumped_mass_system(1_ip,dt_rho_chm)
      if ( kfl_spray_chm /= 0_ip ) &
           call solver_lumped_mass_system(1_ip,dt_chm)

      !
      ! Coupling with partis
      !
      call chm_partis()

      !
      ! Get density
      !
      if ( kfl_spray_chm /= 0_ip ) then

         iclaf_gas = iclaf     - 2
         iclai_liq = iclaf_gas + 1

         do ipoin = 1,npoin
            do iclass = iclai,iclaf_gas
               Tim(iclass,ipoin) = dt_rho_chm(ipoin)/dt
            end do
            do iclass = iclai_liq,iclaf
               Tim(iclass,ipoin) = dt_chm(ipoin)/dt
            end do
         end do
      else
         do ipoin = 1,npoin
            do iclass = iclai,iclaf
               Tim(iclass,ipoin) = dt_rho_chm(ipoin)/dt
            end do
         end do
      end if

      !
      ! Get RHS
      !
      kpoin = 0
      do ipoin = 1,npoin
         do iclass = iclai,iclaf
            kpoin             = kpoin + 1
            rt(iclass,ipoin,istep) = bt(kpoin)
         end do
      end do

      !
      ! Exchange RHS:
      !
      call rhsmod(nclas_chm,rt(:,:,istep))

      !
      ! Solve substeps
      !
      if(istep == 4_ip) then
         call chm_entropy_solution_sij(iclai,iclaf,4_ip,sigma)
      else
         call chm_entropy_solution_sij(iclai,iclaf,istep,b(istep+1,:))
      endif

   end subroutine chm_entropy_eval

   subroutine chm_entropy_solution_sij(iclai,iclaf,istep,weight1)

      integer(ip) , intent(in) :: iclai
      integer(ip) , intent(in) :: iclaf
      integer(ip) , intent(in) :: istep
      real(rp) ,    intent(in) :: weight1(5)
      integer(ip)              :: ipoin,jstep,iclass,kpoin
      real(rp)                 :: aux2

      kpoin = 0
      do ipoin = 1,npoin
         do iclass = iclai,iclaf
            aux2 = 0.0_rp
            kpoin     = kpoin + 1
            do jstep=1,istep
               aux2   = aux2 + rt(iclass,ipoin,jstep)*weight1(jstep)
            end do
            !tt(kpoin) = tt0(iclass,ipoin) + ( Tim(iclass,ipoin,1)*aux2*dt(1) ) / Mass(ipoin)
            aux_tt(iclass,ipoin) = Tim(iclass,ipoin)*aux2*dt
         end do
      end do
      !
      ! Solve system
      !
      call solver_explicit(solve1,aux_tt,EXCHANGE=.false.,solve_consistent=solve2)

      kpoin = 0
      do ipoin = 1,npoin
         do iclass = iclai,iclaf
         kpoin = kpoin+1
         tt(kpoin) = tt0(iclass,ipoin) + aux_tt(iclass,ipoin)
         end do
      end do

      !
      ! Dirichet bbcc
      !
      kpoin = 0
      do ipoin = 1,npoin
         do iclass = iclai,iclaf
            kpoin = kpoin + 1
            if( kfl_fixno_chm(iclass,ipoin) > 0 ) &
                 tt(kpoin) = bvess_chm(iclass,ipoin)
         end do
      end do

   end subroutine chm_entropy_solution_sij

   subroutine chm_entropy_solution()

      integer(ip) :: ipoin,iclas
      integer(ip) :: iclass, iclai, iclaf, kpoin
      real(rp)    :: entr_max, entr_min

      external    :: chm_endite

      bt       => rhsid
      Mass     => vmass

      iclai = 1_ip
      iclaf = nclas_chm

      !
      ! Copy unknown
      !
      kpoin = 0
      do ipoin = 1,npoin
         do iclas = iclai,iclaf
            kpoin = kpoin + 1
            tt(kpoin) = unkno(kpoin)
            tt0(iclas,ipoin) = tt(kpoin)
         end do
      end do

      !
      ! RK substeps
      !
      write(*,'("--| ENTROPY PREDICTION FLAG IS: ",i0)') kfl_disable_entpre_chm
      if (kfl_disable_entpre_chm == 0) then
         call chm_entropy_eval(iclai,iclaf,1_ip)
         call chm_endite(ITASK_INNITE)
         call chm_entropy_eval(iclai,iclaf,2_ip)
         call chm_endite(ITASK_INNITE)
         call chm_entropy_eval(iclai,iclaf,3_ip)
         call chm_endite(ITASK_INNITE)
         call chm_entropy_eval(iclai,iclaf,4_ip)
      end if

      !
      ! Reset unknown
      !
      kpoin = 0
      do ipoin = 1,npoin
         do iclas = iclai,iclaf
            kpoin = kpoin + 1
            tt(kpoin) = unkno(kpoin)
            unkno(kpoin) = tt0(iclas,ipoin)
         end do
      end do

      !
      ! Entropy solution
      !
      kpoin = 0
      do ipoin = 1,npoin
         do iclass = iclai,iclaf
            kpoin           = kpoin + 1
            entro(iclass, ipoin, 2) = tt0(iclass,ipoin)**2
            entro(iclass, ipoin, 1) = tt(kpoin)**2
         enddo
      end do

      !
      ! Compute entropy function and gradient
      !
      do iclass = iclai,iclaf
         entr_max = 0.0_rp
         entr_min = 1.0e10_rp
         do ipoin = 1,npoin
            entr_max = max(abs(entro(iclass,ipoin, 1)),entr_max)
            entr_min = min(abs(entro(iclass,ipoin, 1)),entr_min)
         end do

         call PAR_MAX(entr_max)
         call PAR_MIN(entr_min)

         do ipoin = 1,npoin
            entro(iclass,ipoin, 2) = entro(iclass,ipoin,2)/abs(entr_max - entr_min+1.0e-14_rp)
            entro(iclass,ipoin, 1) = entro(iclass,ipoin,1)/abs(entr_max - entr_min+1.0e-14_rp)
         end do


          !
          ! Compute gradients
          !
          call gradie(entro(iclass,:,2), grads(iclass,:,:))
      enddo

   end subroutine chm_entropy_solution

   !! END OF NEW STUFF !!

   !! OLD VERSION !!

  !subroutine chm_entropy_allocate()

  !  integer(ip), save :: ipass = 0
  !  integer(ip)       :: ii,ielem,pelty,pgaus,iclas


  !  if( ipass == 0 ) then

  !     ipass = 1
  !     nullify(rt)
  !     nullify(Tim)
  !     nullify(tt)
  !     nullify(tt0)
  !     nullify(entro)
  !     nullify(grads)
  !     nullify(envic_gp)

  !     bt       => rhsid
  !     Mass     => vmass

  !     if( kfl_coupl_chm == 2 ) then
  !        !
  !        ! Allocate variables for this module
  !        !
  !        call memory_alloca(mem_modul(1:2,modul),'TT'     ,'chm_entropy_allocate',tt   , max(nclas_chm,npoin*nclas_chm))
  !        call memory_alloca(mem_modul(1:2,modul),'TT0'    ,'chm_entropy_allocate',tt0  , max(nclas_chm,npoin*nclas_chm))
  !        call memory_alloca(mem_modul(1:2,modul),'RT'     ,'chm_entropy_allocate',rt   , nclas_chm, max(1_ip,npoin))
  !        call memory_alloca(mem_modul(1:2,modul),'TIM'    ,'chm_entropy_allocate',Tim  , nclas_chm, max(1_ip,npoin))
  !        call memory_alloca(mem_modul(1:2,modul),'ENTRO'  ,'chm_entropy_allocate',entro, nclas_chm, max(1_ip,npoin), 2_ip)
  !        call memory_alloca(mem_modul(1:2,modul),'GRADS'  ,'chm_entropy_allocate',grads, nclas_chm, ndime, max(1_ip,npoin))

  !        call memory_alloca(mem_modul(1:2,modul),'ENVIC_GP','chm_entropy_allocate', envic_gp, nclas_chm, max(1_ip,nelem))
  !        do ielem = 1,nelem
  !           pelty = ltype(ielem)
  !           pgaus = ngaus(pelty)
  !           do iclas = 1,nclas_chm
  !              call memory_alloca(mem_modul(1:2,modul),'ENVIC_GP % A','chm_entropy_allocate', envic_gp(iclas,ielem)%a, pgaus)
  !           enddo
  !        end do
  !     else
  !        call runend('CHM_ENTROPY_ALLOCATE: Entropy viscosity is only implemented for monolithic coupling.')
  !     endif


  !     tt(1:npoin*nclas_chm) = unkno(1:npoin*nclas_chm)
  !     do ii=1,npoin*nclas_chm
  !        tt(ii) = unkno(ii)
  !     enddo

  !     dt = 1.0_rp
  !  end if

  !end subroutine chm_entropy_allocate


  !subroutine chm_entropy_solution()
  !  integer(ip) :: ipoin,ii
  !  integer(ip) :: iclaf_gas, iclai_liq, iclass, iclai, iclaf, kpoin
  !  integer(ip) :: dummi
  !  real(rp)    :: rho(2), garho,entr_max, entr_min
  !  logical(rp), parameter :: vel_square=.false.
  !
  !  !-----------------------------------------------------------------
  !  !
  !  ! Allocate memory if necessary
  !  !
  !  !-----------------------------------------------------------------

  !  iclai = 1_ip
  !  iclaf = nclas_chm

  !  call chm_entropy_allocate()

  !  if (.not. vel_square) then
  !      !
  !      ! Copy unknown
  !      !
  !      do ii=1,npoin*nclas_chm
  !         tt0(ii) = unkno(ii) ! current unknown
  !         tt(ii)  = unkno(ii) ! prediction of unknown
  !      enddo
  !
  !      !
  !      ! Get time steps
  !      !
  !      dt = 1.0_rp/dtinv
  !      dt_rho_chm = 0.0_rp
  !      if ( kfl_spray_chm /= 0_ip ) &
  !          dt_chm = 0.0_rp

  !      !
  !      ! Get RHS
  !      !
  !      kfl_entropy_chm = 0_ip
  !      call chm_matrix()
  !      kfl_entropy_chm = 1_ip
  !  endif


  !  if( .not.vel_square ) then
  !     call solver_lumped_mass_system(1_ip,dt_rho_chm)
  !     if ( kfl_spray_chm /= 0_ip ) &
  !        call solver_lumped_mass_system(1_ip,dt_chm)

  !     !
  !     ! Coupling with partis
  !     !
  !     call chm_partis()

  !     !
  !     ! Get density:
  !     !
  !     if ( kfl_spray_chm /= 0_ip ) then

  !        iclaf_gas = iclaf     - 2
  !        iclai_liq = iclaf_gas + 1

  !        do ipoin = 1,npoin
  !           do iclass = iclai,iclaf_gas
  !              Tim(iclass,ipoin) = dt_rho_chm(ipoin)/dt
  !           end do
  !           do iclass = iclai_liq,iclaf
  !              Tim(iclass,ipoin) = dt_chm(ipoin)/dt
  !           end do
  !        end do

  !     else
  !        do ipoin = 1,npoin
  !           do iclass = iclai,iclaf
  !              Tim(iclass,ipoin) = dt_rho_chm(ipoin)/dt
  !           end do
  !        end do
  !     end if
  !
  !
  !     !
  !     ! Get right hand side:
  !     !
  !     kpoin = 0
  !     do ipoin = 1,npoin
  !        do iclass = iclai,iclaf
  !           kpoin            = kpoin + 1
  !           rt(iclass,ipoin) = bt(kpoin)
  !        end do
  !     end do

  !
  !     !
  !     ! Exchange RHS:
  !     !
  !     call rhsmod(nclas_chm,rt)
  !  end if

  !
  !  if( INOTMASTER) then
  !     if (vel_square) then
  !         !
  !         ! Use kinetic energy
  !         !
  !         kpoin = 0
  !         do ipoin = 1,npoin
  !            do iclass = iclai,iclaf
  !               kpoin           = kpoin + 1
  !               tt(kpoin) = 0.5_rp*dot_product(advec(1:ndime,ipoin,1),advec(1:ndime,ipoin,1))
  !            enddo
  !         end do
  !         !
  !         ! Entropy solution
  !         !
  !         kpoin = 0
  !         do ipoin = 1,npoin
  !            do iclass = iclai,iclaf
  !               kpoin           = kpoin + 1
  !               entro(iclass, ipoin, 2) = tt0(kpoin)
  !               entro(iclass, ipoin, 1) = tt(kpoin)
  !            enddo
  !         end do

  !     else
  !         !
  !         ! Do an Euler forward step:
  !         !
  !         kpoin = 0
  !         do ipoin = 1,npoin
  !            do iclass = iclai,iclaf
  !               kpoin     = kpoin + 1
  !               tt(kpoin) = tt(kpoin) + ( rt(iclass,ipoin)*Tim(iclass,ipoin) )* dt / Mass(ipoin)
  !            enddo
  !         end do

  !         !
  !         ! Dirichet bbcc
  !         !
  !         kpoin = 0
  !         do ipoin = 1,npoin
  !            do iclass = iclai,iclaf
  !               kpoin            = kpoin + 1
  !               if( kfl_fixno_chm(iclass,ipoin) > 0 ) &
  !                  tt(kpoin) = bvess_chm(iclass,ipoin)
  !           enddo
  !         end do

  !         !
  !         ! Entropy solution
  !         !
  !         kpoin = 0
  !         do ipoin = 1,npoin
  !            do iclass = iclai,iclaf
  !               kpoin           = kpoin + 1
  !               entro(iclass, ipoin, 2) = tt0(kpoin)**2
  !               entro(iclass, ipoin, 1) = tt(kpoin)**2
  !            enddo
  !         end do
  !     endif

  !  end if

  !  do iclass = iclai,iclaf
  !     entr_max = 0.0_rp
  !     entr_min = 1.0e10_rp
  !     do ipoin = 1,npoin
  !        entr_max = max(abs(entro(iclass,ipoin, 1)),entr_max)
  !        entr_min = min(abs(entro(iclass,ipoin, 1)),entr_min)
  !     end do
  !
  !     call PAR_MAX(entr_max)
  !     call PAR_MIN(entr_min)

  !     do ipoin = 1,npoin
  !        entro(iclass,ipoin, 2) = entro(iclass,ipoin,2)/abs(entr_max - entr_min+1.0e-14_rp)
  !        entro(iclass,ipoin, 1) = entro(iclass,ipoin,1)/abs(entr_max - entr_min+1.0e-14_rp)
  !     end do


  !      !
  !      ! Compute gradients
  !      !
  !      call gradie(entro(iclass,:,2), grads(iclass,:,:))
  !  enddo
  !
  !
  !  if (vel_square) then
  !      !
  !      ! Save old step:
  !      !
  !      tt0(1:npoin*nclas_chm) = tt(1:npoin*nclas_chm)
  !  endif
  !end subroutine chm_entropy_solution

 subroutine chm_entropy_viscosity(ielem,pnode,pgaus,igaui,igauf,iclas,elvel,gpden,hleng,del_gpdif)
    integer(ip), intent(in)    :: ielem,pnode,pgaus,igaui,igauf,iclas
    real(rp),    intent(in)    :: elvel(ndime,mnode)
    real(rp),    intent(in)    :: gpden(mgaus)
    real(rp),    intent(in)    :: hleng(ndime)
    real(rp),    intent(inout) :: del_gpdif(mgaus,nclas_chm)

    integer(ip)               :: igaus,inode, ipoin
    real(rp)                  :: Smax, Smin, Savg, elSres(mnode)
    real(rp)                  :: betae, ve, h, gpenv


    del_gpdif(1:pgaus,iclas) = 0.0_rp

    elSres = 0.0_rp
    do inode=1,pnode
       ipoin = lnods(inode,ielem)
       elSres(inode) = elSres(inode) + (entro(iclas,ipoin,1) - entro(iclas,ipoin,2))*dtinv
       elSres(inode) = elSres(inode) + dot_product(elvel(1:ndime,inode),grads(iclas,1:ndime,ipoin))
    end do


   Smin = 1.0e10_rp
   Smax = -1.0e10_rp
   Savg = 0.0_rp
   do inode=1,pnode
      Smax = max(Smax, elSres(inode))
      Smin = max(Smin, elSres(inode))
   !   Savg = Savg + entro(ipoin,1)
   end do


    ve = 0.0_rp
    betae = 0.0_rp
    do inode = 1,pnode
       ipoin = lnods(inode,ielem)
       ve = max(abs(elSres(inode)),ve)
       !ve = max(abs(elSres(inode))/abs(Smax - Smin+1.0e-14_rp),ve)
       betae = max(sqrt(dot_product(elvel(1:ndime,inode),elvel(1:ndime,inode))), betae)
    end do

    h = 1e15_rp
    h = min(hleng(1),h)
    h = min(hleng(2),h)
    if (ndime == 3) then
       h = min(hleng(3),h)
    end if

    do igaus=igaui,igauf
       ! gpenv = 1.0_rp*min(0.5_rp*betae*h*gpden(igaus),dtinv*gpden(igaus)*ve*h**2_ip)
       ! gpenv = 0.5_rp*betae*h*gpden(igaus)
       gpenv = min( 0.5_rp*betae*h*gpden(igaus), 1.0_rp*gpden(igaus)*ve*h**2_ip )

       envic_gp(iclas,ielem) % a(igaus)   = gpenv
       del_gpdif(igaus,iclas)             = del_gpdif(igaus,iclas) + gpenv
    end do

 end subroutine chm_entropy_viscosity


 subroutine chm_entropy_postprocess(itask, envic)
    use def_domain,        only  : npoin
    use mod_memory,        only  : memory_alloca
    use mod_memory,        only  : memory_deallo

    integer(ip), intent(in)     :: itask
    real(rp),    intent(out)    :: envic(:)
    real(rp),   pointer         :: temp_vi(:)
    integer(ip)                 :: ipoin, iclas

    external                    :: smooth

    if (associated(envic_gp)) then
       select case(itask)
       case(0_ip)
          !
          ! Take maximum of values
          !
          nullify(temp_vi)
          call memory_alloca(mem_modul(1:2,modul),'TEMP_VI' ,'chm_entropy_postprocess',temp_vi, max(1_ip,npoin))

          do ipoin=1,npoin
             envic(ipoin) = 0.0_rp
          enddo

          do iclas=1,nclas_chm
             call smooth (envic_gp(iclas,:), temp_vi)
             do ipoin=1,npoin
                envic(ipoin) = max(temp_vi(ipoin),envic(ipoin))
             enddo
          enddo

          call memory_deallo(mem_modul(1:2,modul),'TEMP_VI' ,'chm_entropy_postprocess',temp_vi)

       case default
          !
          ! iclas = itask
          ! Project Gauss point values to nodes
          !
          call smooth (envic_gp(itask,:), envic)
       end select
    else
       do ipoin=1,npoin
          envic(ipoin) = 0.0_rp
       enddo
    endif


 end subroutine chm_entropy_postprocess


 end module mod_chm_entropy
