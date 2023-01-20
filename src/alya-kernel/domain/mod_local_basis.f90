!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------------
!> @addtogroup Domain
!> @{
!> @file    mod_local_basis.f90
!> @author  Gerard Guillamet
!> @date    June, 2018
!>          - Subroutine written
!> @brief   ToolBox for coordinate system transformations and material axes
!>
!> @details ToolBox for coordinate system transformations and material axes
!>
!>          \verbatim
!>          This toolBox includes useful transformations using coordinate
!>          systems
!>          \endverbatim
!>
!> @}
!------------------------------------------------------------------------------

module mod_local_basis

  use def_kintyp,         only : ip, rp, lg
  use def_master,         only : LOCAL_BASIS_CARTESIAN
  use def_master,         only : LOCAL_BASIS_CYLINDRICAL
  use def_master,         only : LOCAL_BASIS_SPHERICAL
  use def_master,         only : LOCAL_BASIS_CONSTANT
  use def_master,         only : ISLAVE
  use def_master,         only : IMASTER
  use mod_maths,          only : maths_vectorial_product
  use mod_maths,          only : maths_normalize_vector
  use mod_maths,          only : maths_matrix_vector_multiplication
  use mod_maths,          only : maths_matrix_transpose_vector_multiplication
  use def_domain,         only : lobas,num_lobas
  use mod_domain,         only : domain_memory_allocate
  use mod_domain,         only : domain_memory_deallocate
  use def_domain,         only : ndime,npoin, lpoty,exnor,coord
  use def_domain,         only : skcos,xfiel
  use mod_maths,          only : maths_vectorial_product
  use mod_maths,          only : maths_local_orthonormal_basis
  use mod_maths,          only : maths_normalize_vector
  use mod_maths_geometry, only : maths_distance_point_segment
  use mod_exchange,       only : exchange_init
  use mod_exchange,       only : exchange_add
  use mod_exchange,       only : exchange_end
  implicit none
  private

  interface local_basis_global_to_local
     module procedure &
          local_basis_global_to_local_s,&
          local_basis_global_to_local_0,&
          local_basis_global_to_local_1,&
          local_basis_global_to_local_2,&
          local_basis_global_to_local_3
  end interface local_basis_global_to_local

  interface local_basis_local_to_global
     module procedure &
          local_basis_local_to_global_s,&
          local_basis_local_to_global_0,&
          local_basis_local_to_global_1,&
          local_basis_local_to_global_2,&
          local_basis_local_to_global_3
  end interface local_basis_local_to_global

  interface local_basis_matrix
     module procedure &
       local_basis_matrix_0,&
       local_basis_matrix_1
  end interface local_basis_matrix

  public  :: local_basis
  public  :: local_basis_local_to_global
  public  :: local_basis_global_to_local
  public  :: local_basis_parall
  public  :: local_basis_matrix

contains

  !-----------------------------------------------------------------------
  !>
  !> @author  houzeaux
  !> @date    2020-03-30
  !> @brief   Coordinate system
  !> @details Broacast data of coordinate system
  !>
  !-----------------------------------------------------------------------

  subroutine local_basis_parall()

    integer(ip) :: ii

    if( ISLAVE ) call domain_memory_allocate('LOBAS',NUMBER1=num_lobas)

    call exchange_init()
    do ii = 1,num_lobas
       call exchange_add(lobas(ii) % type)
       call exchange_add(lobas(ii) % param)
    end do
    call exchange_end()
    
  end subroutine local_basis_parall

  !----------------------------------------------------------------------------
  !> @author  Gerard Guillamet
  !> @date    August, 2018
  !> @brief   Build rotation matrix by using different system of coordinates
  !> @details
  !>
  !>          \verbatim
  !>
  !>           Nomenclature:
  !>           X',Y',Z' : Local cartesian axes
  !>           R        : Radial axis
  !>           T        : Tangencial axis
  !>           Z_       : Axial axis
  !>           N        : Exterior normal
  !>           T1,T2    : First and sencond tangencial axes
  !>
  !>                                               (Rotation matrix
  !>                                                  colum order)
  !>           - Cartesian system   (User defined)    X'- Y'- Z'
  !>
  !>           - Cylindrical system (User defined)    R - T - Zc
  !>
  !>           - Spherical system   (User defined)
  !>
  !>           - Exterior normal    (Alya)            N - T1 - T2
  !>
  !>          \endverbatim
  !>
  !> @warning To implement a function of a projection of a point to continuum 
  !>          line given 2 points. In this case the segment is finite.
  !----------------------------------------------------------------------------

  subroutine local_basis()

    integer(ip) :: ipoin,ibopo,icsys
    real(rp)    :: rho
    real(rp)    :: N(ndime)
    real(rp)    :: R(ndime), Z_(ndime), T(ndime)
    real(rp)    :: a(ndime), b(ndime), proje(ndime), dista
    real(rp)    :: localbasis(ndime,ndime)

    do icsys = 1,num_lobas

       call domain_memory_allocate('LOBAS % LOCAL_BASIS',NUMBER1=icsys)

       do ipoin = 1, npoin

          ibopo      = lpoty(ipoin)
          localbasis = 0.0_rp

          if( ibopo > 0 ) then

             if( lobas(icsys) % type == LOCAL_BASIS_CARTESIAN ) then
                !
                ! Cartesian LOCAL_BASIS
                !
                localbasis(:,1) = lobas(icsys) % param(  ndime+1:ndime*2) - lobas(icsys) % param(1)
                localbasis(:,2) = lobas(icsys) % param(ndime*2+1:ndime*3) - lobas(icsys) % param(2)
                call maths_normalize_vector(ndime,localbasis(:,1),rho)
                call maths_normalize_vector(ndime,localbasis(:,2),rho)
                !
                ! N = P x S
                if( ndime == 3_ip ) then
                   call maths_vectorial_product(localbasis(:,1), localbasis(:,2), N(:), 3_ip)
                   localbasis(:,3) = N(:)
                end if

             else if( lobas(icsys) % type == LOCAL_BASIS_CONSTANT ) then
                !
                ! Constant
                !
                localbasis(:,1) = lobas(icsys) % param(1:ndime)
                localbasis(:,2) = lobas(icsys) % param(ndime+1:2*ndime)
                if( ndime == 3 ) localbasis(:,3) = lobas(icsys) % param(2*ndime+1:3*ndime)

             else if( lobas(icsys) % type == LOCAL_BASIS_CYLINDRICAL ) then
                !
                ! Cylindrical LOCAL_BASIS (Z_ = Local Axial axis)
                !
                if( ndime == 2_ip ) then
                   ! Points defining the axial axis
                   a(1:ndime) = lobas(icsys) % param(1:ndime) ! Center
                   !
                   ! Radial axis
                   R(1) = coord(1,ipoin) - a(1)
                   R(2) = coord(2,ipoin) - a(2)
                   call maths_normalize_vector(ndime,R(:),rho)
                   ! Tangencial axis
                   T(1) = -R(2)
                   T(2) =  R(1)
                   !
                   ! Local basis
                   localbasis(:,1) = R(:)
                   localbasis(:,2) = T(:)

                else if( ndime == 3_ip ) then
                   ! Points defining the axial axis
                   a(1:ndime) = lobas(icsys) % param(1:3) 
                   b(1:ndime) = lobas(icsys) % param(4:6)
                   !
                   ! Axial axis
                   Z_(1) = lobas(icsys) % param(4) - lobas(icsys) % param(1)
                   Z_(2) = lobas(icsys) % param(5) - lobas(icsys) % param(2)
                   Z_(3) = lobas(icsys) % param(6) - lobas(icsys) % param(3)
                   call maths_normalize_vector(ndime,Z_(:),rho)
                   ! Radial axis
                   call maths_distance_point_segment(a,b,coord(:,ipoin),proje,dista)
                   R(1) = coord(1,ipoin) - proje(1)
                   R(2) = coord(2,ipoin) - proje(2)
                   R(3) = coord(3,ipoin) - proje(3)
                   call maths_normalize_vector(ndime,R(:),rho)
                   !
                   ! Local basis
                   localbasis(:,1) = R(:)
                   localbasis(:,3) = Z_(:)
                   ! T = Z x R
                   call maths_vectorial_product(Z_(:), R(:), T(:), 3_ip)
                   call maths_normalize_vector(ndime,T(:),rho)
                   localbasis(:,2) = T(:)

                end if

             else if ( lobas(icsys) % type == LOCAL_BASIS_SPHERICAL ) then
                !
                ! Spherical LOCAL_BASIS
                !
                call runend('MOD_LOCAL_BASIS: NOT IMPLEMENTED')

             else

                call runend('MOD_LOCAL_BASIS: DOES NOT EXIST')

             end if
             !
             ! Built rotation matrix
             !
             lobas(icsys) % local_basis(:,:,ibopo) = localbasis(:,:)
          end if

       end do

    end do

  end subroutine local_basis

  !-----------------------------------------------------------------------
  !>
  !> @author  houzeaux
  !> @date    2019-09-02
  !> @brief   Rotation of velocity
  !> @details This routine rotates the nodal velocities using the appropiate
  !>          rotation matrix:
  !>          ITASK=1 ... From global to local
  !>          ITASK=2 ... From local to global
  !>          Modifications need to be done only if there exist image nodes or
  !>          boundary conditions in skew systems.
  !>
  !-----------------------------------------------------------------------

  subroutine local_basis_local_to_global_s(kfl_fixrs,uu,NODE)

    integer(ip),                    intent(in)    :: kfl_fixrs
    real(rp),                       intent(inout) :: uu(*)
    integer(ip),          optional, intent(in)    :: NODE
    integer(ip), pointer                          :: kfl_fixrs_loc(:)

    if( kfl_fixrs /= 0 ) then
       allocate(kfl_fixrs_loc(1))
       kfl_fixrs_loc(1) = kfl_fixrs
       call local_basis_rotation_RP(2_ip,kfl_fixrs_loc,uu,NODE=NODE)
       deallocate(kfl_fixrs_loc)
    end if

  end subroutine local_basis_local_to_global_s

  subroutine local_basis_local_to_global_0(kdime,kfl_fixrs,uu,skcos_in,NODE)

    integer(ip),                    intent(in)    :: kdime
    integer(ip), pointer,           intent(in)    :: kfl_fixrs(:)
    real(rp),                       intent(inout) :: uu(*)
    real(rp),    pointer, optional, intent(in)    :: skcos_in(:,:,:)
    integer(ip),          optional, intent(in)    :: NODE

    call local_basis_rotation_RP(2_ip,kfl_fixrs,uu,skcos_in,NODE)

  end subroutine local_basis_local_to_global_0

  subroutine local_basis_local_to_global_1(kfl_fixrs,uu,skcos_in,NODE)

    integer(ip), pointer,           intent(in)    :: kfl_fixrs(:)
    real(rp),    pointer,           intent(inout) :: uu(:)
    real(rp),    pointer, optional, intent(in)    :: skcos_in(:,:,:)
    integer(ip),          optional, intent(in)    :: NODE

    if( associated(uu) ) then
       call local_basis_rotation_RP(2_ip,kfl_fixrs,uu,skcos_in,NODE)
    end if

  end subroutine local_basis_local_to_global_1

  subroutine local_basis_local_to_global_2(kfl_fixrs,uu,skcos_in,NODE)

    integer(ip), pointer,           intent(in)    :: kfl_fixrs(:)
    real(rp),    pointer,           intent(inout) :: uu(:,:)
    real(rp)   , pointer, optional, intent(in)    :: skcos_in(:,:,:)
    integer(ip),          optional, intent(in)    :: NODE

    if( associated(uu) ) then
       call local_basis_rotation_RP(2_ip,kfl_fixrs,uu,skcos_in,NODE)
    end if

  end subroutine local_basis_local_to_global_2

  subroutine local_basis_local_to_global_3(kfl_fixrs,uu,skcos_in,LAST_COMPONENT,NODE)

    integer(ip), pointer,           intent(in)    :: kfl_fixrs(:)
    real(rp),    pointer,           intent(inout) :: uu(:,:,:)
    real(rp)   , pointer, optional, intent(in)    :: skcos_in(:,:,:)
    integer(ip),          optional, intent(in)    :: LAST_COMPONENT
    integer(ip),          optional, intent(in)    :: NODE

    if( associated(uu) ) then
       if( present(LAST_COMPONENT) ) then
          call local_basis_rotation_RP(2_ip,kfl_fixrs,uu(:,:,LAST_COMPONENT),skcos_in,NODE)
       else
          call local_basis_rotation_RP(2_ip,kfl_fixrs,uu,skcos_in,NODE)
       end if
    end if

  end subroutine local_basis_local_to_global_3

  subroutine local_basis_global_to_local_s(kfl_fixrs,uu,NODE)

    integer(ip),                    intent(in)    :: kfl_fixrs
    real(rp),                       intent(inout) :: uu(*)
    integer(ip),          optional, intent(in)    :: NODE
    integer(ip), pointer                          :: kfl_fixrs_loc(:)

    if( kfl_fixrs /= 0 ) then
       allocate(kfl_fixrs_loc(1))
       kfl_fixrs_loc(1) = kfl_fixrs
       call local_basis_rotation_RP(1_ip,kfl_fixrs_loc,uu,NODE=NODE)
       deallocate(kfl_fixrs_loc)
    end if

  end subroutine local_basis_global_to_local_s

  subroutine local_basis_global_to_local_0(kdime,kfl_fixrs,uu,skcos_in,NODE)

    integer(ip), pointer,           intent(in)    :: kfl_fixrs(:)
    integer(ip),                    intent(in)    :: kdime
    real(rp),                       intent(inout) :: uu(*)
    real(rp),    pointer, optional, intent(in)    :: skcos_in(:,:,:)
    integer(ip),          optional, intent(in)    :: NODE

    call local_basis_rotation_RP(1_ip,kfl_fixrs,uu,skcos_in,NODE)

  end subroutine local_basis_global_to_local_0

  subroutine local_basis_global_to_local_1(kfl_fixrs,uu,skcos_in,NODE)

    integer(ip), pointer,           intent(in)    :: kfl_fixrs(:)
    real(rp),    pointer,           intent(inout) :: uu(:)
    real(rp),    pointer, optional, intent(in)    :: skcos_in(:,:,:)
    integer(ip),          optional, intent(in)    :: NODE

    if( associated(uu) ) then
       call local_basis_rotation_RP(1_ip,kfl_fixrs,uu,skcos_in,NODE)
    end if

  end subroutine local_basis_global_to_local_1

  subroutine local_basis_global_to_local_2(kfl_fixrs,uu,skcos_in,NODE)

    integer(ip), pointer,           intent(in)    :: kfl_fixrs(:)
    real(rp),    pointer,           intent(inout) :: uu(:,:)
    real(rp)   , pointer, optional, intent(in)    :: skcos_in(:,:,:)
    integer(ip),          optional, intent(in)    :: NODE

    if( associated(uu) ) then
       call local_basis_rotation_RP(1_ip,kfl_fixrs,uu,skcos_in,NODE)
    end if

  end subroutine local_basis_global_to_local_2

  subroutine local_basis_global_to_local_3(kfl_fixrs,uu,skcos_in,LAST_COMPONENT,NODE)

    integer(ip), pointer,           intent(in)    :: kfl_fixrs(:)
    real(rp),    pointer,           intent(inout) :: uu(:,:,:)
    real(rp)   , pointer, optional, intent(in)    :: skcos_in(:,:,:)
    integer(ip),          optional, intent(in)    :: LAST_COMPONENT
    integer(ip),          optional, intent(in)    :: NODE

    if( associated(uu) ) then
       if( present(LAST_COMPONENT) ) then
          call local_basis_rotation_RP(1_ip,kfl_fixrs,uu(:,:,LAST_COMPONENT),skcos_in,NODE)
       else
          call local_basis_rotation_RP(1_ip,kfl_fixrs,uu,skcos_in,NODE)
       end if
    end if

  end subroutine local_basis_global_to_local_3

  subroutine local_basis_rotation_RP(itask,kfl_fixrs,uu,skcos_in,NODE)

    integer(ip),                    intent(in)    :: itask
    integer(ip), pointer,           intent(in)    :: kfl_fixrs(:)
    real(rp),                       intent(inout) :: uu(*)
    real(rp)   , pointer, optional, intent(in)    :: skcos_in(:,:,:)
    integer(ip),          optional, intent(in)    :: NODE
    integer(ip)                                   :: ipoin,ibopo,iroty,itotv
    integer(ip)                                   :: idime,jdime,kdime,kpoin
    integer(ip)                                   :: ipoii,ipoif
    real(rp)                                      :: worma(3),worve(3)
    real(rp),    pointer                          :: rotma(:,:)
    real(rp),    target                           :: rotma_w(ndime,ndime)

    if( associated(kfl_fixrs) ) then

       if( present(NODE) ) then
          ipoii = NODE
          ipoif = NODE
       else
          ipoii = 1
          ipoif = npoin
       end if

       do ipoin = ipoii,ipoif

          ibopo = lpoty(ipoin)
          if( present(NODE) ) then
             kpoin = 1
          else
             kpoin = ipoin
          end if
          iroty = kfl_fixrs(kpoin)

          if( ibopo > 0 .and. iroty /= 0 ) then

             if(      iroty == -1 ) then
                !
                ! EXNOR
                !
                rotma => exnor(:,:,ibopo)

             else if( iroty == -2 ) then
                !
                ! SKCOS_IN
                !
                if( present(skcos_in) ) then
                   rotma => skcos_in(:,:,ibopo)
                else
                   call runend('MOD_LOCAL_BASIS: SKCOS_IN IS MISSING')
                end if

             else if( iroty == -3 ) then
                !
                ! SKCOS
                !
                rotma => skcos(:,:,ibopo)

             else if( iroty > 1000 ) then
                !
                ! LOBAS % LOCAL_BASIS
                !
                rotma => lobas(iroty-1000) % local_basis(:,:,ibopo)

             else if( iroty > 0 ) then
                !
                ! FIELD
                !
                kdime = 0
                do jdime = 1,ndime
                   do idime = 1,ndime
                      kdime = kdime + 1
                      rotma_w(idime,jdime) = xfiel(iroty) % a(kdime,ipoin,1)
                   end do
                end do
                rotma => rotma_w

             end if

             if( itask == 1 ) then
                !
                ! Global to local
                !
                itotv          = (kpoin-1)*ndime
                worve(1:ndime) = uu(itotv+1:itotv+ndime)
                call maths_matrix_transpose_vector_multiplication(ndime,ndime,rotma,worve,worma)
                uu(itotv+1:itotv+ndime) = worma(1:ndime)

             else if( itask == 2 ) then
                !
                ! Local to global
                !
                itotv          = (kpoin-1)*ndime
                worve(1:ndime) = uu(itotv+1:itotv+ndime)
                call maths_matrix_vector_multiplication(ndime,ndime,rotma,worve,worma)
                uu(itotv+1:itotv+ndime) = worma(1:ndime)

             end if

          end if
       end do

    end if

  end subroutine local_basis_rotation_RP

  !-----------------------------------------------------------------------
  !>
  !> @author  houzeaux
  !> @date    2020-03-31
  !> @brief   Matrix
  !> @details Get rotation matrix
  !>
  !-----------------------------------------------------------------------

  subroutine local_basis_matrix_0(ipoin,ibopo,iroty,rotma)

    integer(ip), intent(in)  :: ipoin
    integer(ip), intent(in)  :: ibopo
    integer(ip), intent(in)  :: iroty
    real(rp),    intent(out) :: rotma(ndime,ndime)

    call local_basis_matrix_RP(ipoin,ibopo,iroty,rotma)

  end subroutine local_basis_matrix_0

  subroutine local_basis_matrix_1(ipoin,lpoty_loc,kfl_fixrs,rotma)

    integer(ip),          intent(in)    :: ipoin
    integer(ip), pointer, intent(in)    :: lpoty_loc(:)
    integer(ip), pointer, intent(in)    :: kfl_fixrs(:)
    real(rp),             intent(out)   :: rotma(ndime,ndime)
    integer(ip)                         :: iroty,ibopo

    iroty = kfl_fixrs(ipoin)
    ibopo = lpoty_loc(ipoin)
    call local_basis_matrix_RP(ipoin,ibopo,iroty,rotma)

  end subroutine local_basis_matrix_1

  subroutine local_basis_matrix_RP(ipoin,ibopo,iroty,rotma)

    integer(ip), intent(in)  :: ipoin
    integer(ip), intent(in)  :: iroty
    integer(ip), intent(in)  :: ibopo
    real(rp),    intent(out) :: rotma(ndime,ndime)
    integer(ip)              :: idime,jdime,kdime

    if( ibopo > 0 .and. iroty /= 0 ) then

       if(      iroty == -1 ) then
          !
          ! EXNOR
          !
          rotma(1:ndime,1:ndime) = exnor(1:ndime,1:ndime,ibopo)

       else if( iroty == -2 ) then
          !
          ! SKCOS_IN
          !
          call runend('MOD_LOCAL_BASIS: SKCOS_IN IS MISSING')

       else if( iroty == -3 ) then
          !
          ! SKCOS
          !
          rotma(1:ndime,1:ndime) = skcos(1:ndime,1:ndime,ibopo)

       else if( iroty > 1000 ) then
          !
          ! LOBAS % LOCAL_BASIS
          !
          rotma(1:ndime,1:ndime) = lobas(iroty-1000) % local_basis(1:ndime,1:ndime,ibopo)

       else if( iroty > 0 ) then
          !
          ! FIELD
          !
          kdime = 0
          do jdime = 1,ndime
             do idime = 1,ndime
                kdime = kdime + 1
                rotma(idime,jdime) = xfiel(iroty) % a(kdime,ipoin,1)
             end do
          end do

       end if

    end if

  end subroutine local_basis_matrix_RP

end module mod_local_basis

