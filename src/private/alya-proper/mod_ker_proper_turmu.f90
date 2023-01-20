!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> 
!> @author  houzeaux
!> @date    2020-11-19
!> @brief   Turbulent properties
!> @details Module for element loop to compute properties
!> 
!-----------------------------------------------------------------------

module mod_ker_proper_turmu

#include "def_vector_size.inc"
  use def_master, only : kfl_paral
  use def_vectorization
  use def_kintyp
  use def_domain
  use def_kermod
  use mod_ker_proper_generic, only : elgat_typ
  use mod_ker_proper_generic, only : ker_proper_value
  use mod_ker_proper_generic, only : ker_proper_gradient
  use mod_ker_proper_generic, only : ker_proper_product
  use mod_ker_proper_generic, only : ker_proper_heaviside
  use mod_maths_solver,       only : maths_eigen_3x3_symmetric_matrix
  use mod_ker_proper,         only : ker_proper
  use mod_frivel,             only : frivel
       
#ifdef OPENACCHHH
!$acc routine(ker_proper_gradient_vector) seq
#endif

  implicit none
  private

  integer(ip), parameter :: VECTOR_DIM2 = 1
  real(rp),    parameter :: eps         = epsilon(1.0_rp)
  real(rp),    parameter :: pow(3)      = (/ 1.0_rp , 0.5_rp , 0.3333333_rp /)

#ifdef OPENACCHHH
#define DEF_VECT ivect

#define DEF_VECT_IN ivect:ivect

#else

#define DEF_VECT_IN 1:VECTOR_SIZE

#define DEF_VECT 1:VECTOR_SIZE
#endif

#define DEF_VECT2 1:VECTOR_DIM2
  !
  ! Wall viscosity
  !
  public :: ker_proper_turmu_walvi                  
  !
  ! VREMAN turbulence model
  !
  public :: ker_proper_turmu_vreman_init             
  public :: ker_proper_turmu_vreman_operations       
  !
  ! SMAGORINSKY turbulence model
  !
  public :: ker_proper_turmu_smagorinsky_init       
  public :: ker_proper_turmu_smagorinsky_operations 
  !
  ! SIGMA turbulence model
  !
  public :: ker_proper_turmu_sigma_init       
  public :: ker_proper_turmu_sigma_operations 
  !
  ! WALE turbulence model
  !
  public :: ker_proper_turmu_wale_init       
  public :: ker_proper_turmu_wale_operations 
  !
  ! ILSA turbulence model
  !
  public :: ker_proper_turmu_ilsa_init       
  public :: ker_proper_turmu_ilsa_operations 
  !
  ! TKE-SGS turbulence model
  !
  public :: ker_proper_turmu_tkesgs_init       
  public :: ker_proper_turmu_tkesgs_operations 
 
contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-11-13
  !> @brief   Wall viscosity
  !> @details Variables involving the instantaneous velocity:
  !>          gpgnavv => elgnavv => gpvis_nsw
  !>
  !>          gpvis_nsw(ivect,:) = ( ... ) / abs(elgnavvt(ivect))
  !>
  !-----------------------------------------------------------------------

  subroutine ker_proper_turmu_walvi(ielem,pnode,pgaus,gpcar,elvel,gpvis,gpden,gpvis_nsw)

    integer(ip), intent(in)    :: ielem
    integer(ip), intent(in)    :: pnode
    integer(ip), intent(in)    :: pgaus
    real(rp),    intent(in)    :: gpcar(VECTOR_DIM2,ndime,mnode,pgaus)              ! dN/dxi
    real(rp),    intent(in)    :: elvel(VECTOR_DIM2,ndime,pnode)
    real(rp),    intent(in)    :: gpvis(VECTOR_DIM2,pgaus)
    real(rp),    intent(in)    :: gpden(VECTOR_DIM2,pgaus)
    real(rp),    intent(inout) :: gpvis_nsw(VECTOR_DIM2,pgaus)
    integer(ip)                :: ivect
    real(rp)                   :: elnnsw(VECTOR_DIM2,ndime)
    real(rp)                   :: gpgnavv(VECTOR_DIM2,ndime,pgaus)  ! Gauss Point Gradient in Normal dir of AVerage Velocity
    real(rp)                   :: elgnavv(VECTOR_DIM2,ndime)        ! ELement Gradient in Normal dir of AVerage Velocity
    real(rp)                   :: elgnavvt(VECTOR_DIM2)             ! ELement Gradient in Normal dir of AVerage Tangent Velocity
    real(rp)                   :: elywal(VECTOR_DIM2)
    real(rp)                   :: auxvi(VECTOR_DIM2)
    real(rp)                   :: auxde(VECTOR_DIM2)
    real(rp)                   :: avelavv(VECTOR_DIM2,ndime)        ! AVerage ELement AVerage Velocity
    real(rp)                   :: auxi(VECTOR_DIM2)
    real(rp)                   :: velfr(VECTOR_DIM2)
    real(rp)                   :: rough_aux,vikin,tveno
    integer(ip)                :: idime,jdime,ibopo
    integer(ip)                :: inode,igaus,ipoin
    integer(ip)                :: kount1
    real(rp)                   :: elnor(VECTOR_DIM2)

    elnor   = 0.0_rp
    avelavv = 0.0_rp

    do ivect = 1,VECTOR_DIM2
       !ielem = list_elements(ivect)
       if( ielem > 0 ) then
          if( kfl_nswel_ker(ielem) > 0_ip ) then
             elnnsw(ivect,1:ndime) = normal_nsw_ker(1:ndime,kfl_nswel_ker(ielem)) 
             elnor(ivect)          = dot_product(elnnsw(ivect,1:ndime),elnnsw(ivect,1:ndime))
          else
             elnnsw(ivect,:)       = 0.0_rp
             elnor(ivect)          = 0.0_rp
          end if
          if( kfl_delta == 0 ) then
             elywal(ivect) = delta_dom
          else
             elywal(ivect) = ywale(ielem)
          end if
          avelavv(ivect,1:ndime) = lnsw_exch(ielem) % velav(1:ndime)
       end if
    end do

    !----------------------------------------------------------------------
    !
    ! Gauss point values
    !
    !----------------------------------------------------------------------
    !
    ! GPGNAVV = dj uav_i nj =  dj N_I_i nj  Uav_I  ! gauss point gradient in normal direction of the average velocity
    !
    gpgnavv(DEF_VECT2,:,:) = 0.0_rp

    do igaus = 1,pgaus
       do inode = 1,pnode
          do idime = 1,ndime
             do jdime = 1,ndime
                gpgnavv(DEF_VECT2,idime,igaus) = gpgnavv(DEF_VECT2,idime,igaus) + elvel(DEF_VECT2,idime,inode) * &
                     gpcar(DEF_VECT2,jdime,inode,igaus) * elnnsw(DEF_VECT2,jdime) 
             end do
          end do
       end do
    end do
    !
    ! Obtain the average over gauss points of the normal gradient of the average velocity for all gauss points
    ! I could merge this loop with the upper one and obtain elgnavv directly
    ! Also average density and viscosity
    !
    elgnavv = 0.0_rp
    auxvi   = 0.0_rp
    auxde   = 0.0_rp
    do igaus = 1,pgaus
       do idime = 1,ndime
          elgnavv(DEF_VECT2,idime)  = elgnavv(DEF_VECT2,idime)  + gpgnavv(DEF_VECT2,idime,igaus)
       end do
       auxvi(DEF_VECT2)  = auxvi(DEF_VECT2) + gpvis(DEF_VECT2,igaus)
       auxde(DEF_VECT2)  = auxde(DEF_VECT2) + gpden(DEF_VECT2,igaus)
    end do
    do idime = 1,ndime
       elgnavv(DEF_VECT2,idime) = elgnavv(DEF_VECT2,idime) / real(pgaus,rp)
    end do
    auxvi(DEF_VECT2)  = auxvi(DEF_VECT2)  / real(pgaus,rp)
    auxde(DEF_VECT2)  = auxde(DEF_VECT2)  / real(pgaus,rp)
    !
    ! Substract normal component to keep only tangential one
    !
    !
    auxi(DEF_VECT2) = 0.0_rp
    do idime=1,ndime
       auxi(DEF_VECT2) = auxi(DEF_VECT2) + elgnavv(DEF_VECT2,idime) * elnnsw(DEF_VECT2,idime)
    end do
    do idime=1,ndime
       elgnavv(DEF_VECT2,idime) = elgnavv(DEF_VECT2,idime) - auxi(DEF_VECT2) * elnnsw(DEF_VECT2,idime)
    end do
    !
    ! obtain modulus of tangential component
    !
    auxi(DEF_VECT2) = 0.0_rp
    do idime=1,ndime
       auxi(DEF_VECT2) = auxi(DEF_VECT2) + elgnavv(DEF_VECT2,idime) * elgnavv(DEF_VECT2,idime)
    end do
    elgnavvt(DEF_VECT2) = sqrt(auxi(DEF_VECT2))
    gpvis_nsw = 0.0_rp
    !
    ! Obtain tange from wall_law.
    ! Compute U*: VELFR   this part is not vectorized for the moment  I would need a vector frivel , not difficult
    ! also Time average of (mu+mut) d u_t / dn  - to be usad later
    !
    do ivect = 1,VECTOR_DIM2
       !ielem = list_elements(ivect)
       if( ielem > 0_ip ) then
          if( kfl_nswel_ker(ielem) > 0 .and. elnor(ivect) > eps ) then
             vikin = auxvi(ivect) / auxde(ivect)

             tveno = sqrt(dot_product(avelavv(ivect,1:ndime),avelavv(ivect,1:ndime)))   ! |u_tan-u_fix_tan|
             !call vecnor(avelavv(ivect,1:ndime),ndime,tveno,2_ip)        

             if( kfl_rough == 0 ) then
                !
                ! Constant roughness
                !
                rough_aux = rough_dom

             else if( kfl_rough > 0 ) then
                !
                ! Roughness
                !
                rough_aux = 0.0_rp
                kount1 = 0_ip
                do inode = 1,pnode
                   ipoin = lnods(inode,ielem)
                   ibopo = lpoty(ipoin)
                   if( ibopo /= 0 ) then
                      kount1 = kount1 + 1_ip
                      rough_aux = rough_aux + rough(ipoin)
                   end if
                end do
                rough_aux = rough_aux / real(kount1,rp)

             else
                !
                ! No roughness
                !
                rough_aux = 0.0_rp

             end if

             call frivel(kfl_ustar,elywal(ivect),rough_aux,tveno,vikin,velfr(ivect))      ! u*

             gpvis_nsw(ivect,1) = ( velfr(ivect) * velfr(ivect) * auxde(ivect) ) / abs(elgnavvt(ivect)+1.0e-30_rp)
             do igaus = 2,pgaus
                gpvis_nsw(ivect,igaus) = gpvis_nsw(ivect,1) 
             end do

          end if
       end if
    end do

  end subroutine ker_proper_turmu_walvi

  !-----------------------------------------------------------------------
  !> 
  !> @author  lehmkuhl
  !> @date    2020-11-19
  !> @brief   VREMAN
  !> @details VREMAN turbulence model
  !> 
  !-----------------------------------------------------------------------

  subroutine ker_proper_turmu_vreman_init(elgat)

    type(elgat_typ), intent(inout) :: elgat

    elgat % kfl_elcod = .true.
    elgat % kfl_elvel = .true.

  end subroutine ker_proper_turmu_vreman_init

  subroutine ker_proper_turmu_vreman_operations(&
       VECTOR_DIM,list_elements_loc,pelty,pnode,pgaus,gpsha,gpcar,&
       hleng,elgat,param,gpnut,grnut)

    integer(ip),     intent(in)    :: VECTOR_DIM
    integer(ip),     intent(in)    :: list_elements_loc(VECTOR_DIM) 
    integer(ip),     intent(in)    :: pelty
    integer(ip),     intent(in)    :: pnode
    integer(ip),     intent(in)    :: pgaus
    real(rp),        intent(in)    :: gpsha(pnode,pgaus)
    real(rp),        intent(in)    :: gpcar(VECTOR_DIM,ndime,mnode,pgaus)
    real(rp),        intent(in)    :: hleng(VECTOR_DIM,ndime)
    type(elgat_typ), intent(in)    :: elgat
    real(rp),        intent(in)    :: param(mlapa_ker)
    real(rp),        intent(inout) :: gpnut(VECTOR_DIM,pgaus)
    real(rp),        intent(inout) :: grnut(VECTOR_DIM,ndime,pgaus)
    integer(ip)                    :: inode,igaus,idime,jdime,kdime,ivect
    integer(ip)                    :: VECTOR_DIM_IN
    real(rp)                       :: gpgve(VECTOR_DIM,ndime,ndime)
    real(rp)                       :: alpha(VECTOR_DIM)
    real(rp)                       :: Bbeta(VECTOR_DIM)
    real(rp)                       :: xmile(VECTOR_DIM)
    real(rp)                       :: G__ij(VECTOR_DIM,3,3)

#ifdef OPENACCHHH
    VECTOR_DIM_IN = 1
    !$acc enter data create (G__ij, gpgve(1:VECTOR_DIM,1:ndime,1:ndime), xmile, alpha, Bbeta, gpnut       ) &
    !$acc            copyin ( elgat, elgat%elvel(1:VECTOR_DIM,1:ndime,1:pnode), param) 
#else
    VECTOR_DIM_IN = VECTOR_DIM
#endif

#ifdef OPENACCHHH
    !$acc parallel loop gang vector default(present)
    do ivect = 1,VECTOR_DIM
#endif
       xmile = ker_proper_product(VECTOR_DIM,ndime,ivect,hleng)

       xmile(DEF_VECT) = xmile(DEF_VECT) ** pow(ndime)

       do igaus = 1,pgaus

          gpgve(DEF_VECT,:,:) = 0.0_rp
          do inode = 1,pnode
             do idime = 1,ndime
                do jdime = 1,ndime
                   gpgve(DEF_VECT,jdime,idime) =  gpgve(DEF_VECT,jdime,idime) &
                        + elgat%elvel(DEF_VECT,idime,inode) * gpcar(DEF_VECT,jdime,inode,igaus)
                end do
             end do
          end do
 
          G__ij(DEF_VECT,:,:) = 0.0_rp
          !
          ! G = g^T*g
          !
          do idime = 1_ip,ndime
             do jdime = 1_ip,ndime
                do kdime = 1_ip,ndime
                   G__ij(DEF_VECT,idime,jdime) = G__ij(DEF_VECT,idime,jdime) &
                        + (gpgve(DEF_VECT,idime,kdime)*gpgve(DEF_VECT,jdime,kdime)*xmile(DEF_VECT)*xmile(DEF_VECT))
                end do
             end do
          end do

          alpha(DEF_VECT) = (G__ij(DEF_VECT,1_ip,1_ip) + G__ij(DEF_VECT,2_ip,2_ip) &
               &           + G__ij(DEF_VECT,3_ip,3_ip))/(xmile(DEF_VECT)*xmile(DEF_VECT))
          Bbeta(DEF_VECT) =  G__ij(DEF_VECT,1_ip,1_ip) * G__ij(DEF_VECT,2_ip,2_ip) &
               &           + G__ij(DEF_VECT,2_ip,2_ip) * G__ij(DEF_VECT,3_ip,3_ip) &
               &           + G__ij(DEF_VECT,3_ip,3_ip) * G__ij(DEF_VECT,1_ip,1_ip) &
               &           - G__ij(DEF_VECT,1_ip,2_ip) * G__ij(DEF_VECT,1_ip,2_ip) &
               &           - G__ij(DEF_VECT,2_ip,3_ip) * G__ij(DEF_VECT,2_ip,3_ip) &
               &           - G__ij(DEF_VECT,1_ip,3_ip) * G__ij(DEF_VECT,1_ip,3_ip)
          !
          ! If ALPHA < 10 * EPS, the GPNUT = 0
          !
          call ker_proper_heaviside(VECTOR_DIM_IN,alpha(DEF_VECT_IN),10.0_rp*eps, gpnut(DEF_VECT_IN,igaus)) 
          gpnut(DEF_VECT,igaus) = gpnut(DEF_VECT,igaus) * param(1) * sqrt (max( Bbeta(DEF_VECT)/max(alpha(DEF_VECT),eps),eps) )

     end do

#ifdef OPENACCHHH
    end do
    !$acc end parallel loop

    !$acc update self ( gpnut(1:VECTOR_DIM,1:pgaus) )

    !$acc exit data delete(G__ij, gpgve, xmile, alpha, Bbeta, &
    !$acc                     elgat, elgat%elvel, param, gpnut) 

#endif

  end subroutine ker_proper_turmu_vreman_operations

  !-----------------------------------------------------------------------
  !> 
  !> @author  lehmkuhl
  !> @date    2020-11-19
  !> @brief   SMAGORINSKY
  !> @details SMAGORINSKY turbulence model
  !> 
  !-----------------------------------------------------------------------

  subroutine ker_proper_turmu_smagorinsky_init(elgat)

    type(elgat_typ), intent(inout) :: elgat

    elgat % kfl_elcod = .true.
    elgat % kfl_elvel = .true.

  end subroutine ker_proper_turmu_smagorinsky_init

  subroutine ker_proper_turmu_smagorinsky_operations(&
       VECTOR_DIM,list_elements_loc,pelty,pnode,pgaus,gpsha,gpcar,&
       hleng,elgat,param,gpnut,grnut)

    integer(ip),     intent(in)    :: VECTOR_DIM
    integer(ip),     intent(in)    :: list_elements_loc(VECTOR_DIM) 
    integer(ip),     intent(in)    :: pelty
    integer(ip),     intent(in)    :: pnode
    integer(ip),     intent(in)    :: pgaus
    real(rp),        intent(in)    :: gpsha(pnode,pgaus)
    real(rp),        intent(in)    :: gpcar(VECTOR_DIM,ndime,mnode,pgaus)
    real(rp),        intent(in)    :: hleng(VECTOR_DIM,ndime)
    type(elgat_typ), intent(in)    :: elgat
    real(rp),        intent(in)    :: param(mlapa_ker)
    real(rp),        intent(inout) :: gpnut(VECTOR_DIM,pgaus)
    real(rp),        intent(inout) :: grnut(VECTOR_DIM,ndime,pgaus)
    integer(ip)                    :: VECTOR_DIM_IN
    integer(ip)                    :: inode,igaus,idime,jdime,ivect
    real(rp)                       :: gpgve(VECTOR_DIM,ndime,ndime)
    real(rp)                       :: seci4(VECTOR_DIM)
    real(rp)                       :: xmile(VECTOR_DIM)
 
#ifdef OPENACCHHH
    VECTOR_DIM_IN = 1
    !$acc enter data create (gpgve(1:VECTOR_DIM,1:ndime,1:ndime), xmile, seci4, gpnut ) &
    !$acc            copyin ( param,  &
    !$acc                    elgat, elgat%elvel(1:VECTOR_DIM,1:ndime,1:pnode), param ) 
#else
    VECTOR_DIM_IN = VECTOR_DIM
#endif

   
#ifdef OPENACCHHH
    !$acc parallel loop gang vector default(present)
    do ivect = 1,VECTOR_DIM
#endif
      xmile = ker_proper_product(VECTOR_DIM,ndime,ivect,hleng)
      xmile(DEF_VECT) = xmile(DEF_VECT) ** pow(ndime)


       do igaus = 1,pgaus

          gpgve(DEF_VECT,:,:) = 0.0_rp
          do inode = 1,pnode
             do idime = 1,ndime
                do jdime = 1,ndime
                   gpgve(DEF_VECT,jdime,idime) =  gpgve(DEF_VECT,jdime,idime) &
                        + elgat%elvel(DEF_VECT,idime,inode) * gpcar(DEF_VECT,jdime,inode,igaus)
                end do
             end do
          end do
        
          !
          ! 2 S_ij : S_ij
          !
          seci4(DEF_VECT) = 0.0_rp
          do idime = 1,ndime          
             do jdime = 1,ndime         
                seci4(DEF_VECT) = seci4(DEF_VECT) + gpgve(DEF_VECT,idime,jdime) &
                     &            * (gpgve(DEF_VECT,idime,jdime) + gpgve(DEF_VECT,jdime,idime))
             end do
          end do
          gpnut(DEF_VECT,igaus) = param(1)*xmile(DEF_VECT)*xmile(DEF_VECT)*sqrt(abs(seci4(DEF_VECT)))
       end do

#ifdef OPENACCHHH
    end do
    !$acc end parallel loop

    !$acc update self ( gpnut )

    !$acc exit data delete (gpgve, xmile, seci4, &
    !$acc                   param, gpnut, elgat, elgat%elvel, param  ) 

#endif

  end subroutine ker_proper_turmu_smagorinsky_operations

  !-----------------------------------------------------------------------
  !> 
  !> @author  lehmkuhl
  !> @date    2020-11-19
  !> @brief   SIGMA
  !> @details SIGMA turbulence model. Smart values to test:
  !>
  !>           gpgve               =  0.0_rp
  !>           gpgve(DEF_VECT,1,1) =  1.0_rp
  !>           gpgve(DEF_VECT,1,2) =  2.0_rp
  !>           gpgve(DEF_VECT,1,3) =  3.0_rp
  !>           gpgve(DEF_VECT,2,1) = -4.0_rp
  !>           gpgve(DEF_VECT,2,2) = -5.0_rp
  !>           gpgve(DEF_VECT,2,3) =  6.0_rp
  !>           gpgve(DEF_VECT,3,1) = -7.0_rp
  !>           gpgve(DEF_VECT,3,2) =  8.0_rp
  !>           gpgve(DEF_VECT,3,3) =  9.0_rp
  !> 
  !-----------------------------------------------------------------------

  subroutine ker_proper_turmu_sigma_init(elgat)

    type(elgat_typ), intent(inout) :: elgat

    elgat % kfl_elcod = .true.
    elgat % kfl_elvel = .true.

  end subroutine ker_proper_turmu_sigma_init

  subroutine ker_proper_turmu_sigma_operations(&
       VECTOR_DIM,list_elements_loc,pelty,pnode,pgaus,gpsha,gpcar,&
       hleng,elgat,param,gpnut,grnut)

    integer(ip),     intent(in)    :: VECTOR_DIM
    integer(ip),     intent(in)    :: list_elements_loc(VECTOR_DIM) 
    integer(ip),     intent(in)    :: pelty
    integer(ip),     intent(in)    :: pnode
    integer(ip),     intent(in)    :: pgaus
    real(rp),        intent(in)    :: gpsha(pnode,pgaus)
    real(rp),        intent(in)    :: gpcar(VECTOR_DIM,ndime,mnode,pgaus)
    real(rp),        intent(in)    :: hleng(VECTOR_DIM,ndime)
    type(elgat_typ), intent(in)    :: elgat
    real(rp),        intent(in)    :: param(mlapa_ker)
    real(rp),        intent(inout) :: gpnut(VECTOR_DIM,pgaus)
    real(rp),        intent(inout) :: grnut(VECTOR_DIM,ndime,pgaus)
    integer(ip)                    :: VECTOR_DIM_IN
    integer(ip)                    :: inode,igaus,idime,jdime,kdime,ivect
    real(rp)                       :: gpgve(VECTOR_DIM,ndime,ndime)
    real(rp)                       :: G__ij(VECTOR_DIM,3,3)
    real(rp)                       :: G_val(VECTOR_DIM,3)
    real(rp)                       :: xmile(VECTOR_DIM)
    real(rp)                       :: sigma(VECTOR_DIM,3)
    real(rp)                       :: denom(VECTOR_DIM)
    real(rp)                       :: vdumm(3,3)

#ifdef OPENACCHHH
    VECTOR_DIM_IN = 1
    !$acc enter data create (G__ij, G_val, gpgve(1:VECTOR_DIM,1:ndime,1:ndime), xmile, sigma, denom, gpnut ) &
    !$acc            copyin (elgat, elgat%elvel(1:VECTOR_DIM,1:ndime,1:pnode), param  ) 
#else
    VECTOR_DIM_IN = VECTOR_DIM
#endif


#ifdef OPENACCHHH
    !$acc parallel loop gang vector default(present)
    do ivect = 1,VECTOR_DIM
#endif

      xmile = ker_proper_product(VECTOR_DIM,ndime,ivect,hleng)
      xmile(DEF_VECT) = xmile(DEF_VECT) ** pow(ndime)


       do igaus = 1,pgaus

          gpgve(DEF_VECT,:,:) = 0.0_rp
          do inode = 1,pnode
             do idime = 1,ndime
                do jdime = 1,ndime
                   gpgve(DEF_VECT,jdime,idime) =  gpgve(DEF_VECT,jdime,idime) &
                        + elgat%elvel(DEF_VECT,idime,inode) * gpcar(DEF_VECT,jdime,inode,igaus)
                end do
             end do
          end do
 
          G__ij(DEF_VECT,:,:) = 0.0_rp                      ! G = g^T*g = g_ki * g_kj 
          G_val(DEF_VECT,:)   = 0.0_rp                      ! Note gpgve(kdime,idime) is g_ik = du_i/dx_k from  Nicoud
          
          do idime = 1,ndime
             do jdime = 1,ndime
                do kdime = 1,ndime
                   G__ij(DEF_VECT,idime,jdime) = G__ij(DEF_VECT,idime,jdime) &
                        + gpgve(DEF_VECT,idime,kdime)*gpgve(DEF_VECT,jdime,kdime)
                end do
             end do
          end do
          
#ifndef OPENACCHHH
          do ivect = 1,VECTOR_DIM
             call maths_eigen_3x3_symmetric_matrix(G__ij(ivect,:,:),G_val(ivect,:),vdumm)
          end do
#endif
          sigma(DEF_VECT,1) = sqrt(max(G_val(DEF_VECT,1),0.0_rp))
          sigma(DEF_VECT,2) = sqrt(max(G_val(DEF_VECT,2),0.0_rp))
          sigma(DEF_VECT,3) = sqrt(max(G_val(DEF_VECT,3),0.0_rp))
          denom(DEF_VECT)   = sigma(DEF_VECT,1) * sigma(DEF_VECT,1)
          
          call ker_proper_heaviside(VECTOR_DIM_IN,denom(DEF_VECT_IN),1.0e-10_rp,gpnut(DEF_VECT_IN,igaus))

          gpnut(DEF_VECT,igaus) = gpnut(DEF_VECT,igaus) * param(1) * xmile(DEF_VECT) * xmile(DEF_VECT) &
               & * ( sigma(DEF_VECT,3) * ( sigma(DEF_VECT,1) - sigma(DEF_VECT,2) )  &
               & * ( sigma(DEF_VECT,2) - sigma(DEF_VECT,3) ) ) / ( max(denom(DEF_VECT),eps) )

       end do
     
#ifdef OPENACCHHH
    end do
    !$acc end parallel loop

    !$acc update self ( gpnut )
    !$acc exit data delete (G__ij, G_val, gpgve, xmile, sigma, denom, &
    !$acc                   gpnut, elgat, elgat%elvel, param ) 


#endif

  end subroutine ker_proper_turmu_sigma_operations

  !-----------------------------------------------------------------------
  !> 
  !> @author  lehmkuhl
  !> @date    2020-11-19
  !> @brief   WALE
  !> @details WALE turbulence model
  !> 
  !-----------------------------------------------------------------------

  subroutine ker_proper_turmu_wale_init(elgat)

    type(elgat_typ), intent(inout) :: elgat

    elgat % kfl_elcod = .true.
    elgat % kfl_elvel = .true.

  end subroutine ker_proper_turmu_wale_init

  subroutine ker_proper_turmu_wale_operations(&
       VECTOR_DIM,list_elements_loc,pelty,pnode,pgaus,gpsha,gpcar,&
       hleng,elgat,param,gpnut,grnut)

    integer(ip),     intent(in)    :: VECTOR_DIM
    integer(ip),     intent(in)    :: list_elements_loc(VECTOR_DIM) 
    integer(ip),     intent(in)    :: pelty
    integer(ip),     intent(in)    :: pnode
    integer(ip),     intent(in)    :: pgaus
    real(rp),        intent(in)    :: gpsha(pnode,pgaus)
    real(rp),        intent(in)    :: gpcar(VECTOR_DIM,ndime,mnode,pgaus)
    real(rp),        intent(in)    :: hleng(VECTOR_DIM,ndime)
    type(elgat_typ), intent(in)    :: elgat
    real(rp),        intent(in)    :: param(mlapa_ker)
    real(rp),        intent(inout) :: gpnut(VECTOR_DIM,pgaus)
    real(rp),        intent(inout) :: grnut(VECTOR_DIM,ndime,pgaus)
    integer(ip)                    :: inode,igaus,idime,jdime,kdime,ivect
    integer(ip)                    :: VECTOR_DIM_IN
    real(rp)                       :: gpgve(VECTOR_DIM,ndime,ndime)
    real(rp)                       :: g2_kk(VECTOR_DIM)
    real(rp)                       :: g2_ij(VECTOR_DIM,3,3)
    real(rp)                       :: seci4(VECTOR_DIM)
    real(rp)                       :: seci5(VECTOR_DIM)
    real(rp)                       :: sd_ij(VECTOR_DIM)
    real(rp)                       :: denom(VECTOR_DIM)
    real(rp)                       :: xmile(VECTOR_DIM)

#ifdef OPENACCHHH
    VECTOR_DIM_IN = 1
    !$acc enter data create (gpgve(1:VECTOR_DIM,1:ndime,1:ndime), xmile, g2_kk, g2_ij, seci4,         &
    !$acc                    seci5, sd_ij, denom, gpnut                     )  &
    !$acc            copyin (elgat, elgat%elvel(1:VECTOR_DIM,1:ndime,1:pnode), param ) 
#else
    VECTOR_DIM_IN = VECTOR_DIM
#endif


#ifdef OPENACCHHH
    !$acc parallel loop gang vector default(present)
    do ivect = 1,VECTOR_DIM
#endif
       
      xmile = ker_proper_product(VECTOR_DIM,ndime,ivect,hleng)
      xmile(DEF_VECT) = xmile(DEF_VECT) ** pow(ndime)


       do igaus = 1,pgaus
          !
          ! WALE  - Wall-adapting local eddy-viscosity - Nicoud-Ducros 99
          ! 
          ! Computation of the square of the velocity gradient g2_ij
          !
          gpgve(DEF_VECT,:,:) = 0.0_rp
          do inode = 1,pnode
             do idime = 1,ndime
                do jdime = 1,ndime
                   gpgve(DEF_VECT,jdime,idime) =  gpgve(DEF_VECT,jdime,idime) &
                        + elgat%elvel(DEF_VECT,idime,inode) * gpcar(DEF_VECT,jdime,inode,igaus)
                end do
             end do
          end do
  
          g2_ij(DEF_VECT,:,:) = 0.0_rp          
          do idime = 1,ndime
             do jdime = 1,ndime
                do kdime = 1,ndime
                   g2_ij(DEF_VECT,idime,jdime) = g2_ij(DEF_VECT,idime,jdime) &
                        + gpgve(DEF_VECT,idime,kdime)*gpgve(DEF_VECT,kdime,jdime)
                end do
             end do
          end do
          
          g2_kk(DEF_VECT) = g2_ij(DEF_VECT,1,1) + g2_ij(DEF_VECT,2,2) + g2_ij(DEF_VECT,3,3)
          
          seci4(DEF_VECT) = 0.0_rp
          seci5(DEF_VECT) = 0.0_rp
          sd_ij(DEF_VECT) = 0.0_rp
          denom(DEF_VECT) = 0.0_rp
          
          do idime = 1,ndime
             do jdime = 1,ndime
                seci4(DEF_VECT) = seci4(DEF_VECT) + 0.5_rp*gpgve(DEF_VECT,idime,jdime) * &          ! S_ij : S_ij
                     (gpgve(DEF_VECT,idime,jdime) + gpgve(DEF_VECT,jdime,idime)) 
                sd_ij(DEF_VECT) = 0.5_rp * ( g2_ij(DEF_VECT,idime,jdime) + g2_ij(DEF_VECT,jdime,idime) ) 
                seci5(DEF_VECT) = seci5(DEF_VECT) + sd_ij(DEF_VECT)*sd_ij(DEF_VECT)                 ! Sd_ij:Sd_ij
             end do
             
          end do
          
          seci5(DEF_VECT)       = seci5(DEF_VECT) - g2_kk(DEF_VECT)*g2_kk(DEF_VECT)/3.0_rp         
          denom(DEF_VECT)       = max ( eps , seci4(DEF_VECT)**2.5_rp + seci5(DEF_VECT)**1.25_rp )
          gpnut(DEF_VECT,igaus) = param(1) * xmile(DEF_VECT) * xmile(DEF_VECT) * (seci5(DEF_VECT)**1.5_rp)/denom(DEF_VECT)

       end do
     
#ifdef OPENACCHHH
    end do
    !$acc end parallel loop

    !$acc update self ( gpnut )
 
    !$acc exit data delete (gpgve, xmile, g2_kk, g2_ij, seci4, seci5,  &
    !$acc                    sd_ij, denom, hleng, gpnut, elgat, elgat%elvel, param ) 

#endif

  end subroutine ker_proper_turmu_wale_operations

 !-----------------------------------------------------------------------
  !> 
  !> @author  Matias
  !> @date    2020-11-19
  !> @brief   TKE-SGS
  !> @details TKE-SGS turbulence model
  !> 
  !-----------------------------------------------------------------------
  
  subroutine ker_proper_turmu_tkesgs_init(elgat)

    type(elgat_typ), intent(inout) :: elgat

    elgat % kfl_elvel = .true.
    elgat % kfl_eltem = .true.
    elgat % kfl_eltur = .true.     
    elgat % kfl_elwal = .true.

  end subroutine ker_proper_turmu_tkesgs_init

  subroutine ker_proper_turmu_tkesgs_operations(&
       VECTOR_DIM,list_elements_loc,pelty,pnode,pgaus,gpsha,gpcar,&
       hleng,elgat,param,gpnut,grnut)


    use def_master, only  : tempe
    integer(ip),     intent(in)    :: VECTOR_DIM
    integer(ip),     intent(in)    :: list_elements_loc(VECTOR_DIM) 
    integer(ip),     intent(in)    :: pelty
    integer(ip),     intent(in)    :: pnode
    integer(ip),     intent(in)    :: pgaus
    real(rp),        intent(in)    :: gpsha(pnode,pgaus)
    real(rp),        intent(in)    :: gpcar(VECTOR_DIM,ndime,pnode,pgaus)
    real(rp),        intent(in)    :: hleng(VECTOR_DIM,ndime)
    type(elgat_typ), intent(in)    :: elgat
    real(rp),        intent(in)    :: param(mlapa_ker)
    real(rp),        intent(inout) :: gpnut(VECTOR_DIM,pgaus)
    real(rp),        intent(inout) :: grnut(VECTOR_DIM,ndime,pgaus)
    integer(ip)                    :: VECTOR_DIM_IN
    integer(ip)                    :: inode,igaus,idime,ivect
#ifdef OPENACCHHH
    real(rp)                       :: gpgve(VECTOR_DIM,ndime,ndime)
#endif
    real(rp)                       :: gptke(VECTOR_DIM), gpwal(VECTOR_DIM)
    real(rp)                       :: length(VECTOR_DIM)
    real(rp)                       :: xmile(VECTOR_DIM)
    real(rp)                       :: teref, grtez(VECTOR_DIM)
    real(rp)                       :: grtem(VECTOR_DIM, ndime), toler(VECTOR_DIM)

    toler(DEF_VECT) =1.0e-8
#ifdef OPENACCHHH
    VECTOR_DIM_IN = 1
    !$acc enter data create (gpgve(1:VECTOR_DIM,1:ndime,1:ndime), xmile, gpnut ) &
    !$acc            copyin ( param,  &
    !$acc                    elgat, elgat%eltur(1:VECTOR_DIM,1:pnode), &
    !$acc                    elgat%eltem(1:VECTOR_DIM,1:pnode),elgat%elwal(1:VECTOR_DIM,1:pnode), param ) 
#else
    VECTOR_DIM_IN = VECTOR_DIM
#endif

   
#ifdef OPENACCHHH
    !$acc parallel loop gang vector default(present)
    loop_ivect:do ivect = 1,VECTOR_DIM
#endif
       ! xmile= (h1*h2*h3)^(1/3) ! can be greater than volume when not using heahedras
      xmile = ker_proper_product(VECTOR_DIM,ndime,ivect,hleng)
      xmile(DEF_VECT) = xmile(DEF_VECT) ** pow(ndime)

      
      ! vertical gr of tempe
      teref = 290.0_rp
      do igaus =1, pgaus
         ! TKE at gauss point
         gptke = 0.0_rp
         do inode = 1,pnode
            gptke(DEF_VECT) = gptke(DEF_VECT) &
                 + elgat%eltur(DEF_VECT,1,inode) * gpsha(inode, igaus)
         end do
         
         if (associated(tempe)) then ! thermally coupled
            !
            ! Compute termal things
            !
            ! Temperature gradient opossed to gravity
            grtem = 0.0_rp ! tempe gradient
            grtez = 0.0_rp ! tempe gradient projected in direction opossed to gravity
            gpwal=  0.0_rp
            do inode = 1,pnode
               do idime = 1,ndime
                  grtem(DEF_VECT,idime) = grtem(DEF_VECT, idime) + gpcar(DEF_VECT,idime,inode, igaus) * elgat%eltem(DEF_VECT,inode)
               end do
               gpwal(DEF_VECT) =  gpwal(DEF_VECT) + elgat%elwal(DEF_VECT,inode)*gpsha(inode, igaus)
            end do
            do idime = 1,ndime
               grtez(DEF_VECT) =  grtez(DEF_VECT) - grtem(DEF_VECT,idime)  * gravi(idime) ! gravi is a unit vector 
            end do
        
            ! Calculation of mixing length for stable flows
!!$            
            where (gpwal.gt.70.0_rp.and.grtez.gt.toler(1)) 
               length = max (min(xmile, 0.76_rp*sqrt(gptke*teref/(grtez*grnor))),xmile*0.01_rp)
            elsewhere
               length = xmile
            end where
           
         else ! neutral problem, not thermaly coupled
            length = xmile
         end if
         
         gpnut(DEF_VECT, igaus) = 0.1_rp*length(DEF_VECT)*sqrt(gptke(DEF_VECT))
         ! ALERT, Prandtl number is missing
      end do
      
#ifdef OPENACCHHH
    end do loop_ivect
    !$acc end parallel loop

    !$acc update self ( gpnut )

    !$acc exit data delete (gpgve, xmile, &
    !$acc                   param, gpnut, elgat, elgat%eltur, elgat%eltem,elgat%elwal, param  ) 

#endif

  end subroutine ker_proper_turmu_tkesgs_operations

  !-----------------------------------------------------------------------
  !> 
  !> @author  lehmkuhl
  !> @date    2020-11-19
  !> @brief   ILSA
  !> @details ILSA turbulence model
  !> 
  !-----------------------------------------------------------------------

  subroutine ker_proper_turmu_ilsa_init(elgat)

    type(elgat_typ), intent(inout) :: elgat

    elgat % kfl_elcod = .true.
    elgat % kfl_elvel = .true.
    elgat % kfl_elwal = .true.

  end subroutine ker_proper_turmu_ilsa_init

  subroutine ker_proper_turmu_ilsa_operations(&
       VECTOR_DIM,list_elements_loc,pelty,pnode,pgaus,gpsha,gpcar,&
       hleng,elgat,param,gpnut,grnut)

    use mod_ker_ILSA_sgsmodel
    integer(ip),     intent(in)    :: VECTOR_DIM
    integer(ip),     intent(in)    :: list_elements_loc(VECTOR_DIM) 
    integer(ip),     intent(in)    :: pelty
    integer(ip),     intent(in)    :: pnode
    integer(ip),     intent(in)    :: pgaus
    real(rp),        intent(in)    :: gpsha(pnode,pgaus)
    real(rp),        intent(in)    :: gpcar(VECTOR_DIM,ndime,mnode,pgaus)
    real(rp),        intent(in)    :: hleng(VECTOR_DIM,ndime)
    type(elgat_typ), intent(in)    :: elgat
    real(rp),        intent(in)    :: param(mlapa_ker)
    real(rp),        intent(inout) :: gpnut(VECTOR_DIM,pgaus)
    real(rp),        intent(inout) :: grnut(VECTOR_DIM,ndime,pgaus)
    integer(ip)                    :: VECTOR_DIM_IN,ivect
    integer(ip)                    :: igaus,dummi
    real(rp)                       :: gpgve(VECTOR_DIM,ndime,ndime,pgaus)
    real(rp)                       :: gpwal(VECTOR_DIM,pgaus)
    real(rp)                       :: gpden(VECTOR_DIM,pgaus)
    real(rp)                       :: gpvis(VECTOR_DIM,pgaus)
    real(rp)                       :: gpvel(VECTOR_DIM,ndime,pgaus)
    real(rp)                       :: gpnuILSA(VECTOR_DIM,pgaus)
 
#ifdef OPENACCHHH
    VECTOR_DIM_IN = 1
#else
    VECTOR_DIM_IN = VECTOR_DIM
#endif

   
#ifdef OPENACCHHH
    do ivect = 1,VECTOR_DIM
#endif
       
       do igaus = 1,pgaus
          !
          ! ILSA  - Wall-adapting local eddy-viscosity - Nicoud-Ducros 99
          ! 
          ! Computation of the square of the velocity gradient g2_ij
          !
          call ker_proper('DENSI','PGAUS',dummi,list_elements_loc,gpden,pnode,pgaus,VECTOR_DIM=VECTOR_DIM)
          call ker_proper('VISCO','PGAUS',dummi,list_elements_loc,gpvis,pnode,pgaus,VECTOR_DIM=VECTOR_DIM)
          call ker_proper_gradient(VECTOR_DIM_IN,mnode,pnode,pgaus,gpcar(DEF_VECT_IN,:,:,:),elgat % elvel(DEF_VECT_IN,:,:),gpgve(DEF_VECT_IN,:,:,:))
 
          gpwal                    = ker_proper_value   (VECTOR_DIM,pnode,pgaus,gpsha,elgat % elwal)
          gpvel                    = ker_proper_value   (VECTOR_DIM,pnode,pgaus,gpsha,elgat % elvel)
          gpnuILSA(DEF_VECT,igaus) = gpvis(DEF_VECT,igaus)/gpden(DEF_VECT,igaus)

       end do
     
#ifdef OPENACCHHH
    end do
#endif

#ifndef OPENACCHHH
    do ivect = 1,VECTOR_DIM
#endif       
       call ker_ILSA_sgs_viscosity(&
            pgaus,1_ip,pgaus,gpnuILSA(ivect,:),gpnut(ivect,:),gpvel(ivect,:,:),&
            gpgve(ivect,:,:,:),list_elements_loc(ivect),param(1),& 
            param(2),param(3),hleng(ivect,:),gpwal(ivect,:))
          
#ifndef OPENACCHHH
    end do
#endif

  end subroutine ker_proper_turmu_ilsa_operations

end module mod_ker_proper_turmu
!> @}
