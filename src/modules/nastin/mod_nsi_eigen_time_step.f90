!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!




!-----------------------------------------------------------------------
!> @addtogroup Nastin
!> @{
!> @file    mod_nsi_eigen_time_step.f90
!> @author  Guillaume Houzeaux and Oriol Lehmkuhl
!> @brief   Critical time step
!> @details Critical time step based on eigenvalues
!-----------------------------------------------------------------------

module mod_nsi_eigen_time_step

#include "def_vector_size.inc"
  use def_kintyp,              only : ip,rp
  use def_parame,              only : pi
  use def_master,              only : INOTEMPTY,veloc,mem_modul
  use def_master,              only : velom,modul
  use def_domain,              only : ndime,ntens,mnode,elmar,npoin,lnods
  use def_domain,              only : vmass,mnode,coord,lgaus,ltype,nnode
  use def_domain,              only : nelem,lnnod
  use def_domain,              only : ompss_domain   ! Required by OMPSS
  use def_domain,              only : ompss_domains
  use mod_ker_proper,          only : ker_proper
  use mod_communications,      only : PAR_MAX
  use def_nastin,              only : fcons_nsi, kfl_algor_nsi
  use mod_memory,              only : memory_alloca
  use mod_memory,              only : memory_deallo
  use mod_element_integration, only : element_shape_function_derivatives_jacobian
  use mod_element_integration, only : element_shape_function_derivatives_jacobian_vector_2
  use mod_elmgeo_vector,       only : elmgeo_cartesian_derivatives_jacobian_vectorized_cpu
  use mod_parall,              only : par_omp_nelem_chunk
  use mod_parall,              only : num_subd_par
  use mod_parall,              only : num_pack_par
  use mod_parall,              only : list_elements_par
  use mod_parall,              only : typ_list_elements_par
  use mod_parall,              only : num_subd_cpu
  use mod_parall,              only : num_pack_cpu
  use mod_parall,              only : list_elements_cpu       
 
  implicit none
  private

  real(rp), pointer :: convective_pos_term(:,:)
  real(rp), pointer :: convective_neg_term(:,:)
  real(rp), pointer :: viscous_term(:,:)
  
  public :: nsi_eigen_time_step_all

contains
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-07-27
  !> @brief   Compute the time step based on matrix eigenvalues
  !> @details Compartible with vectorization, OmpSs and OpenMP
  !> 
  !-----------------------------------------------------------------------

  subroutine nsi_eigen_time_step_all(dt,dtp)

    real(rp),    optional,          intent(out)   :: dt
    real(rp),    optional, pointer, intent(inout) :: dtp(:)
    integer(ip)                                   :: isubd,ipack,ielem
    integer(ip)                                   :: pnode,pgaus
#ifdef ALYA_OMPSS
    integer(ip)                                   :: jsubd
    integer(ip)                                   :: num_neigh
#endif
    integer(ip)                                   :: ipoin,idime,VECTOR_DIM
    integer(ip)                                   :: num_subd
    integer(ip),                 pointer          :: num_pack(:)
    type(typ_list_elements_par), pointer          :: list_elements(:)
    real(rp)                                      :: convective_max(ndime)
    real(rp)                                      :: viscous_max(ndime)
    real(rp)                                      :: fact0,fact1,phi

    nullify(convective_pos_term)
    nullify(convective_neg_term)
    nullify(viscous_term)

    call memory_alloca(mem_modul(1:2,modul),'CONVECTIVE_POS_TERM','nsi_eigen_time_step_element_operations',convective_pos_term,ndime,npoin)
    call memory_alloca(mem_modul(1:2,modul),'CONVECTIVE_NEG_TERM','nsi_eigen_time_step_element_operations',convective_neg_term,ndime,npoin)
    call memory_alloca(mem_modul(1:2,modul),'VISCOUS_TERM'       ,'nsi_eigen_time_step_element_operations',viscous_term       ,ndime,npoin)

    !if( VECTOR_SIZE == VECTOR_SIZE_CPU ) then
       num_subd      =  num_subd_par
       num_pack      => num_pack_par
       list_elements => list_elements_par
    !else
    !   num_subd      =  num_subd_cpu 
    !   num_pack      => num_pack_cpu
    !   list_elements => list_elements_cpu
    !end if

    do isubd = 1,num_subd

#ifdef ALYA_OMPSS
       
       num_neigh = size(ompss_domains(isubd) % neighbours,KIND=ip)

       !-----------------------------------------------------------------------------
       !$OMP TASK         COMMUTATIVE(                                              &
       !$OMP              [ompss_domains(ompss_domains(isubd) % neighbours(jsubd)), &
       !$OMP              jsubd = 1,num_neigh] ) PRIORITY(num_neigh)                &
       !$OMP FIRSTPRIVATE ( num_neigh,jsubd,isubd )                                 &
       !$OMP SHARED       ( ompss_domains )                                         &
       !-----------------------------------------------------------------------------
#else
       !-----------------------------------------------------------------------------
       !$OMP PARALLEL DO                                                            &
       !$OMP SCHEDULE     ( DYNAMIC , par_omp_nelem_chunk )                         &
       !$OMP SHARED       ( isubd,par_omp_nelem_chunk )                             &
       !-----------------------------------------------------------------------------
#endif
       !-----------------------------------------------------------------------------
       !$OMP DEFAULT      ( SHARED )                                                &
       !$OMP PRIVATE      ( ipack,pnode,pgaus,ielem,VECTOR_DIM )                    &
       !$OMP SHARED       ( list_elements,num_pack,lnnod,lgaus )                     
       !-----------------------------------------------------------------------------
       
       do ipack = 1,num_pack(isubd)

          ielem      = list_elements(isubd) % packs(ipack) % l(1)                 ! Select first element
          pnode      = lnnod(ielem)                                               ! Number of nodes
          pgaus      = lgaus(ielem)                                               ! Number of Gauss points
          VECTOR_DIM = size(list_elements(isubd) % packs(ipack) % l,kind=ip)

          call nsi_eigen_time_step_element_operations(&
               VECTOR_DIM,pnode,pgaus,&
               list_elements(isubd) % packs(ipack) % l)
       end do

#ifdef ALYA_OMPSS
       !$OMP END TASK
#else
       !$OMP END PARALLEL DO
#endif
    end do
#ifdef ALYA_OMPSS
    !$OMP  TASKWAIT
#endif
    !
    ! Project
    ! 
    if( INOTEMPTY ) then
       call rhsmod(ndime,convective_pos_term)
       call rhsmod(ndime,convective_neg_term)
       call rhsmod(ndime,viscous_term)
       do ipoin = 1,npoin
          do idime = 1,ndime
             convective_pos_term(idime,ipoin) = ( convective_pos_term(idime,ipoin) + convective_neg_term(idime,ipoin) ) / vmass(ipoin)
             viscous_term(idime,ipoin)        = viscous_term(idime,ipoin) / vmass(ipoin)
          end do
       end do
       !
       ! Max over nodes
       !
       do idime = 1,ndime
          convective_max(idime) = maxval(convective_pos_term(idime,:))
          viscous_max(idime)    = maxval(viscous_term(idime,:))
       end do
    else
       convective_max = 0.0_rp
       viscous_max    = 0.0_rp
    end if

    if( present(dt) ) then
       !
       ! Global time step
       !
       call PAR_MAX(ndime,convective_max)
       call PAR_MAX(ndime,viscous_max)       
       fact0 = 2.0_rp * maxval(viscous_max) 
       fact1 = maxval(convective_max)    
       if(abs(fact1) < 1e-15_rp ) then
          phi = 0.0_rp
       else
          phi = pi*0.5_rp - atan(fact0/fact1)
       endif
       if(kfl_algor_nsi == 3 .and. abs(fact1)>1e-10_rp) then
          phi = 0.5_rp*pi
          dt = nsi_eigen_time_step_Topt(phi) / sqrt(fact1*fact1) 
       else
          dt = nsi_eigen_time_step_Topt_static(phi) / sqrt(fact0*fact0 + fact1*fact1) 
       end if
       
    else if ( present(dtp) .and. INOTEMPTY ) then
       !
       ! Nodal time step
       !
       do ipoin = 1,npoin
          fact0 = 2.0_rp * maxval(viscous_term(:,ipoin)) 
          fact1 = maxval(convective_pos_term(:,ipoin))
          if(abs(fact1) < 1e-15_rp ) then
             phi = 0.0_rp
          else
             phi = pi*0.5_rp - atan(fact0/fact1)
          endif
          if(kfl_algor_nsi == 3 .and. abs(fact1)>1e-10_rp) then
             phi        = 0.5_rp*pi
             dtp(ipoin) = nsi_eigen_time_step_Topt(phi) / sqrt(fact1*fact1)
          else
             print*,phi,nsi_eigen_time_step_Topt_static(phi)
             dtp(ipoin) = nsi_eigen_time_step_Topt_static(phi) / sqrt(fact0*fact0 + fact1*fact1)           
          end if
       end do
    end if
    
    call memory_deallo(mem_modul(1:2,modul),'CONVECTIVE_POS_TERM','nsi_eigen_time_step_element_operations',convective_pos_term)
    call memory_deallo(mem_modul(1:2,modul),'CONVECTIVE_NEG_TERM','nsi_eigen_time_step_element_operations',convective_neg_term)
    call memory_deallo(mem_modul(1:2,modul),'VISCOUS_TERM'       ,'nsi_eigen_time_step_element_operations',viscous_term)

  end subroutine nsi_eigen_time_step_all  

  function nsi_eigen_time_step_Kopt(phi) 

    real(rp), intent(in) :: phi
    real(rp)             :: nsi_eigen_time_step_Kopt
    real(rp)             :: c6,c7,c8,c9,c10,c11,c12,c13,c14
    real(rp)             :: phi1,phi2,phi3,k1,k2

    c6   =  2403400.0_rp
    c7   = -5018490.0_rp
    c8   =  2620140.0_rp
    c9   =  2945.0_rp
    c10  = -6665.76_rp
    c11  =  3790.54_rp
    c12  =  4.80513_rp
    c13  = -16.9473_rp
    c14  =  15.0155_rp
    phi1 =  atan(164.0_rp/99.0_rp)
    phi2 =  pi/3.0_rp
    phi3 = (3.0_rp/5.0_rp)*(3.0_rp/5.0_rp)*pi
    k1   =  0.73782212_rp 
    k2   =  0.44660387_rp

    nsi_eigen_time_step_Kopt = 0.0_rp

    if      (phi>=0.0_rp .and. phi<=phi1  ) then
       nsi_eigen_time_step_Kopt = 1.0_rp

    else if (phi> phi1   .and. phi<=phi2  ) then
      nsi_eigen_time_step_Kopt = nsi_eigen_time_step_G(phi,c6 ,c7 ,c8 ,phi1,phi2  ,1.0_rp,k1 )

    else if (phi> phi2   .and. phi<=phi3  ) then
      nsi_eigen_time_step_Kopt = nsi_eigen_time_step_G(phi,c9,c10,c11,phi2,phi3  ,k1 ,k2 )

    else if (phi> phi3   .and. phi<=pi/2.0_rp) then
      nsi_eigen_time_step_Kopt = nsi_eigen_time_step_G(phi,c12,c13,c14,phi3,pi/2.0_rp,k2 ,0.0_rp)

    else 
       call runend('TE HA PETAO!')
    end if

  end function nsi_eigen_time_step_Kopt

  function nsi_eigen_time_step_G( x, a, b, c, x0, x1, f0, f1) 
     real(rp), intent(in) :: x, a, b, c, x0, x1, f0, f1
     real(rp)             :: nsi_eigen_time_step_G

     nsi_eigen_time_step_G = (a*x*x+b*x+c) * nsi_eigen_time_step_Q(x,x0,x1) + nsi_eigen_time_step_L(x,x0,x1,f0,f1)

  end function nsi_eigen_time_step_G

  function nsi_eigen_time_step_Q( x, x0, x1) 
     real(rp), intent(in) :: x, x0, x1
     real(rp)             :: nsi_eigen_time_step_Q

     nsi_eigen_time_step_Q = ((x-x0)*(x-x1))

  end function nsi_eigen_time_step_Q

  function nsi_eigen_time_step_L( x, x0, x1, f0, f1) 
     real(rp), intent(in) :: x, x0, x1, f0, f1
     real(rp)             :: nsi_eigen_time_step_L

     !if( abs(x1-x0)<1.0e-12 ) print*,'c=',x0,x1
     nsi_eigen_time_step_L = (f0+(x-x0)*(f1-f0)/(x1-x0))

  end function nsi_eigen_time_step_L

  function nsi_eigen_time_step_Topt(phi)

    real(rp), intent(in) :: phi
    real(rp)             :: nsi_eigen_time_step_Topt
    real(rp)             :: c1,c2,c3,c4,c5
    real(rp)             :: phi1,t1

    c1   =  0.0647998_rp
    c2   = -0.386022_rp
    c3   =  3.72945_rp
    c4   = -9.38143_rp
    c5   =  7.06574_rp
    phi1 = atan(164.0_rp/99.0_rp)
    t1   =  0.9302468_rp

    nsi_eigen_time_step_Topt = 0.0_rp

    if      (phi>=0.0_rp .and. phi<=phi1  ) then
      nsi_eigen_time_step_Topt = nsi_eigen_time_step_G(phi,0.0_rp,c1,c2,0.0_rp ,phi1,4.0_rp/3.0_rp,t1 )

    else if (phi> phi1   .and. phi<=pi/2.0_rp) then
      nsi_eigen_time_step_Topt = nsi_eigen_time_step_G(phi,c3,c4,c5,phi1,pi/2.0_rp,t1,1.0_rp)

    else 
       call runend('TE HA PETAO!')
    end if

  end function nsi_eigen_time_step_Topt

  function nsi_eigen_time_step_Topt_static(phi)

    real(rp), intent(in) :: phi
    real(rp)             :: nsi_eigen_time_step_Topt_static
    real(rp)             :: a1,a2,a3,a4
    real(rp)             :: b1,b2,b3,b4
    real(rp)             :: c1,c2,c3,c4
    real(rp)             :: f0,f1,f2,f3,f4
    real(rp)             :: phi0,phi1,phi2,phi3,phi4


    a1   = -0.111497_rp
    b1   =  0.0246904_rp
    c1   = -0.114765_rp

    a2   = -36.1939_rp
    b2   =  90.0298_rp
    c2   = -58.1013_rp

    a3   = -22671.5_rp
    b3   =  68271.6_rp
    c3   = -51439.8_rp

    a4   =  0.0_rp
    b4   =  0.0_rp
    c4   =  0.0_rp 

    phi0 =  0.0_rp
    f0   =  1.0_rp

    phi1 =  2.0_rp*pi/5.0_rp
    f1   =  0.837063_rp

    phi2 =  19.0_rp*pi/40.0_rp
    f2   =  0.60948_rp

    phi3 =  997.0_rp*pi/2000.0_rp
    f3   =  0.261504_rp

    phi4 =  pi/2.0_rp
    f4   =  1.32e-04_rp

    nsi_eigen_time_step_Topt_static = 0.0_rp

    !if(kfl_paral<=0) print*,'b=',phi,phi0,phi1,phi2,phi3,phi4

    if      (phi>=phi0 .and. phi<=phi1  ) then
      nsi_eigen_time_step_Topt_static = nsi_eigen_time_step_G(phi,a1,b1,c1,phi0,phi1,f0,f1)

    else if (phi> phi1   .and. phi<=phi2) then
      nsi_eigen_time_step_Topt_static = nsi_eigen_time_step_G(phi,a2,b2,c2,phi1,phi2,f1,f2)

    else if (phi> phi2   .and. phi<=phi3) then
      nsi_eigen_time_step_Topt_static = nsi_eigen_time_step_G(phi,a3,b3,c3,phi2,phi3,f2,f3)

    else if (phi> phi3   .and. phi<=phi4) then
      nsi_eigen_time_step_Topt_static = nsi_eigen_time_step_G(phi,a4,b4,c4,phi3,phi4,f3,f4)

    else
       call runend('TE HA PETAO!')
    end if

  end function nsi_eigen_time_step_Topt_static

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-07-27
  !> @brief   Element operations to compute time step
  !> @details Assemble the nodal convective and diffusive terms.
  !>          The convective term is split into a negative and positive
  !>          part to be able to take the row absolute value
  !>          ELCOND = \int_V u.grad(Nj) Ni dV
  !>          ELVISC = \int_V mu*grad(Nj)*grad(Ni) dV
  !>          Their units are therefore [ m^3 * rho / s ]
  !>
  !>          There are two methods:
  !>          - At center of gravity
  !>          - Using the user-defined intergration rule
  !> 
  !-----------------------------------------------------------------------
  
  subroutine nsi_eigen_time_step_element_operations(VECTOR_DIM,pnode,pgaus,list_elements)

    integer(ip),          intent(in) :: VECTOR_DIM
    integer(ip),          intent(in) :: pnode
    integer(ip),          intent(in) :: pgaus
    integer(ip), pointer, intent(in) :: list_elements(:)
    integer(ip)                      :: ielem,pelty,inode,ipoin,jdime
    integer(ip)                      :: idime,ivect,plapl,kgaus,dummi
    integer(ip)                      :: jnode,ielem0(1),igaus
    real(rp)                         :: elcod(VECTOR_SIZE,ndime,mnode)
    real(rp)                         :: elvel(VECTOR_SIZE,ndime,mnode)
    real(rp)                         :: elvem(VECTOR_SIZE,ndime,mnode)
    real(rp)                         :: elconv_pos(VECTOR_SIZE,ndime,mnode)
    real(rp)                         :: elconv_neg(VECTOR_SIZE,ndime,mnode)
    real(rp)                         :: elvisc(VECTOR_SIZE,ndime,mnode)
    real(rp)                         :: fact0(VECTOR_SIZE)
    real(rp)                         :: fact1(VECTOR_SIZE) 
    real(rp)                         :: fact2(VECTOR_SIZE) 
    real(rp)                         :: gpsha(VECTOR_SIZE,pnode,pgaus)        ! N
    real(rp)                         :: gpcar(VECTOR_SIZE,ndime,mnode,pgaus)  ! dN/dxi
    real(rp)                         :: gpvol(VECTOR_SIZE,pgaus)              ! w*|J|, |J|
    real(rp)                         :: gpden(VECTOR_SIZE,pgaus)
    real(rp)                         :: gpvis(VECTOR_SIZE,pgaus)
    real(rp)                         :: gpnut(VECTOR_SIZE,pgaus)
    real(rp)                         :: agrau(VECTOR_SIZE,pnode)
    real(rp)                         :: gpvel(VECTOR_SIZE,ndime)
    real(rp)                         :: gpvem(VECTOR_SIZE,ndime)
    real(rp)                         :: xjaci(VECTOR_SIZE,ndime,ndime,pgaus)
    integer(ip)                      :: list_elements_p(VECTOR_SIZE)          ! List of elements (always positive)

#define DEF_VECT 1:VECTOR_SIZE_CPU
    !
    ! Initialization
    !    
    elconv_pos = 0.0_rp
    elconv_neg = 0.0_rp
    elvisc     = 0.0_rp
    plapl      = 0
    kgaus      = 1 !pgaus
    pelty      = ltype(list_elements(1))

    xjaci(:,:,:,:) = 0.0_rp

    !
    ! List of elements. Put last non-zero element in the list to
    ! avoid zero arrays
    !    
    list_elements_p = list_elements
    ielem0 = minloc(list_elements,list_elements==0)
    if(ielem0(1) > 1 ) list_elements_p(ielem0(1):VECTOR_SIZE) = list_elements(ielem0(1)-1)

    !--------------------------------------------------------------------
    !
    ! Gather
    !
    !--------------------------------------------------------------------

    do ivect = 1,VECTOR_SIZE
       ielem = list_elements_p(ivect)
       do inode = 1,pnode
          ipoin = lnods(inode,ielem)
          elcod(ivect,1:ndime,inode) = coord(1:ndime,ipoin)
          elvel(ivect,1:ndime,inode) = veloc(1:ndime,ipoin,1)
       end do
    end do
    if( associated(velom) ) then
       do ivect = 1,VECTOR_SIZE
          ielem = list_elements_p(ivect)
          do inode = 1,pnode
             ipoin = lnods(inode,ielem)
             elvem(ivect,1:ndime,inode) = velom(1:ndime,ipoin)
          end do
       end do
    end if

    !--------------------------------------------------------------------
    !
    ! Assembly
    !
    !--------------------------------------------------------------------
    !
    ! Cartesian derivatives, Hessian matrix and volume: GPCAR, GPVOL
    !
    if( kgaus == 1 ) then
       call elmgeo_cartesian_derivatives_jacobian_vectorized_cpu(&
            ndime,mnode,pnode,kgaus,elcod,elmar(pelty) % dercg,xjaci,gpcar,gpvol)
       do igaus = 1,kgaus
          gpvol(DEF_VECT,igaus) = gpvol(DEF_VECT,igaus) * elmar(pelty) % weicg
          do inode = 1,pnode
             gpsha(DEF_VECT,inode,igaus) = elmar(pelty) % shacg(inode)
          end do
       end do
    else
       call elmgeo_cartesian_derivatives_jacobian_vectorized_cpu(&
            ndime,mnode,pnode,kgaus,elcod,elmar(pelty) % deriv,xjaci,gpcar,gpvol)
       do igaus = 1,kgaus
          gpvol(DEF_VECT,igaus) = gpvol(DEF_VECT,igaus) * elmar(pelty) % weigp(igaus)
          do inode = 1,pnode
             gpsha(DEF_VECT,inode,igaus) = elmar(pelty) % shape(inode,igaus)
          end do
       end do
    end if
    !
    ! Properties
    !
    if( kgaus == 1 ) then
       call ker_proper('DENSI','COG  ',dummi,list_elements_p,gpden,pnode,kgaus,gpsha,gpcar,VECTOR_DIM=VECTOR_DIM)
       call ker_proper('VISCO','COG  ',dummi,list_elements_p,gpvis,pnode,kgaus,gpsha,gpcar,VECTOR_DIM=VECTOR_DIM)
       !call ker_proper('TURBU','COG  ',dummi,list_elements_p,gpnut,pnode,kgaus,gpsha,gpcar,VECTOR_DIM=VECTOR_DIM)
       gpnut(DEF_VECT,:) = 0.0_rp
    else
       call ker_proper('DENSI','PGAUS',dummi,list_elements_p,gpden,pnode,kgaus,gpsha,gpcar,VECTOR_DIM=VECTOR_DIM)  
       call ker_proper('VISCO','PGAUS',dummi,list_elements_p,gpvis,pnode,kgaus,gpsha,gpcar,VECTOR_DIM=VECTOR_DIM)         
       !call ker_proper('TURBU','PGAUS',dummi,list_elements_p,gpnut,pnode,kgaus,gpsha,gpcar,VECTOR_DIM=VECTOR_DIM)         
       gpnut(DEF_VECT,:) = 0.0_rp
    end if
    !
    ! a.grad(Ni)
    !
    do igaus = 1,kgaus
       gpvis(DEF_VECT,igaus) = gpvis(DEF_VECT,igaus) + gpden(DEF_VECT,igaus) * gpnut(DEF_VECT,igaus) 
       gpvel(DEF_VECT,:)     = 0.0_rp
       gpvem(DEF_VECT,:)     = 0.0_rp
       agrau(DEF_VECT,:)     = 0.0_rp 
       do inode = 1,pnode
          do idime = 1,ndime
             gpvel(DEF_VECT,idime) = gpvel(DEF_VECT,idime) + elvel(DEF_VECT,idime,inode) * gpsha(DEF_VECT,inode,igaus)
          end do
       end do
       if( associated(velom) ) then
          do inode = 1,pnode
             do idime = 1,ndime
                gpvem(DEF_VECT,idime) = gpvem(DEF_VECT,idime) + elvem(DEF_VECT,idime,inode) * gpsha(DEF_VECT,inode,igaus)
             end do
          end do          
       end if
       do inode = 1,pnode
          do idime = 1,ndime
             agrau(DEF_VECT,inode) = agrau(DEF_VECT,inode) + (gpvel(DEF_VECT,idime)-gpvem(DEF_VECT,idime)) * gpcar(DEF_VECT,idime,inode,igaus)
          end do
       end do
       !
       ! Viscous and convective term
       !
       do inode = 1,pnode
          do jnode = 1,pnode
             if( jnode /= inode ) then
                fact0(DEF_VECT) = gpvol(DEF_VECT,igaus) * agrau(DEF_VECT,jnode) * gpsha(DEF_VECT,inode,igaus)
                do idime = 1,ndime
                   elconv_neg(DEF_VECT,idime,inode) = elconv_neg(DEF_VECT,idime,inode) - min(fact0(DEF_VECT),0.0_rp) ! Avoid if: only negative fact0
                   elconv_pos(DEF_VECT,idime,inode) = elconv_pos(DEF_VECT,idime,inode) + max(fact0(DEF_VECT),0.0_rp) ! Avoid if: only positive fact0
                end do
             end if
          end do
          do idime = 1,ndime
             do jdime = 1,ndime
                elvisc(DEF_VECT,idime,inode) = elvisc(DEF_VECT,idime,inode) &
                     + abs( gpvol(DEF_VECT,igaus) * gpvis(DEF_VECT,igaus) / gpden(DEF_VECT,igaus) &
                     * gpcar(DEF_VECT,jdime,inode,igaus) * gpcar(DEF_VECT,jdime,inode,igaus))
             end do
          end do
       end do
       !
       ! Skew-symmetric
       !
       if( fcons_nsi > 0.1_rp ) then   
          fact0(DEF_VECT) = fcons_nsi * gpvol(DEF_VECT,igaus)
          do inode = 1,pnode
             do jnode = 1,pnode
                if( inode /= jnode ) then
                   fact1(DEF_VECT) = agrau(DEF_VECT,inode) * gpsha(DEF_VECT,jnode,igaus) * fact0(DEF_VECT)
                   fact2(DEF_VECT) = agrau(DEF_VECT,jnode) * gpsha(DEF_VECT,inode,igaus) * fact0(DEF_VECT)
                   do idime = 1,ndime
                      elconv_neg(DEF_VECT,idime,inode) = elconv_neg(DEF_VECT,idime,inode) - min(fact1(DEF_VECT),0.0_rp) - min(fact2(DEF_VECT),0.0_rp)
                      elconv_pos(DEF_VECT,idime,inode) = elconv_pos(DEF_VECT,idime,inode) + max(fact1(DEF_VECT),0.0_rp) + max(fact2(DEF_VECT),0.0_rp) 
                   end do
                end if
             end do
          end do
       end if
    end do

    !--------------------------------------------------------------------
    !
    ! Scatter
    !
    !--------------------------------------------------------------------

    do ivect = 1,VECTOR_SIZE
       ielem = list_elements(ivect)
       if( ielem > 0 ) then
          if( ltype(ielem) > 0 ) then
             do inode = 1,pnode
                ipoin = lnods(inode,ielem)
                do idime = 1,ndime
#ifdef NO_COLORING                
                   !$OMP ATOMIC
#endif
                   convective_pos_term(idime,ipoin) = convective_pos_term(idime,ipoin) + elconv_pos(ivect,idime,inode)
#ifdef NO_COLORING
                   !$OMP ATOMIC
#endif
                   convective_neg_term(idime,ipoin) = convective_neg_term(idime,ipoin) + elconv_neg(ivect,idime,inode)
#ifdef NO_COLORING
                   !$OMP ATOMIC
#endif
                   viscous_term(idime,ipoin)        = viscous_term(idime,ipoin)        + elvisc(ivect,idime,inode)
                end do
             end do
          end if
       end if
    end do

  end subroutine nsi_eigen_time_step_element_operations
  
end module mod_nsi_eigen_time_step
!> @}
