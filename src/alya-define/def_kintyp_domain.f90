!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Kinds_and_types
!> @{
!> @file    def_kintyp_domain.g90
!> @author  houzeaux
!> @date    2020-04-04
!> @brief   Functions
!> @details Communications
!-----------------------------------------------------------------------

module def_kintyp_domain

  use def_kintyp_basic, only : ip,rp,i1p
  use def_kintyp_dims,  only : nelty 
  use def_kintyp_dims,  only : mfiel 
 
  !----------------------------------------------------------------------
  !
  ! Element data base type
  !
  !----------------------------------------------------------------------

  type ::  elm
     ! Isoparametric
     integer(ip)          :: pgaus            ! Number of Gauss points
     real(rp),    pointer :: shape(:,:)       ! pnode,pgaus
     real(rp),    pointer :: deriv(:,:,:)     ! ndime,pnode,pgaus
     real(rp),    pointer :: heslo(:,:,:)     ! ntens,pnode,pgaus
     ! User integration rule
     real(rp),    pointer :: weigp(:)
     real(rp),    pointer :: shaga(:,:)
     real(rp),    pointer :: posgp(:,:)
     ! Bubble
     real(rp),    pointer :: shape_bub(:,:)   ! 1,pgaus
     real(rp),    pointer :: deriv_bub(:,:,:) ! ndime,1,pgaus
     real(rp),    pointer :: heslo_bub(:,:,:) ! ntens,1,pgaus
     ! Center of gravity
     real(rp),    pointer :: shacg(:)
     real(rp),    pointer :: dercg(:,:)
     real(rp),    pointer :: hescg(:,:)
     real(rp)             :: weicg
     ! Close rule
     real(rp),    pointer :: shapc(:,:)
     real(rp),    pointer :: deric(:,:,:)
     real(rp),    pointer :: heslc(:,:,:)
     real(rp),    pointer :: weigc(:)
     ! IB integration rule
     real(rp),    pointer :: shaib(:,:)
     real(rp),    pointer :: derib(:,:,:)
     real(rp),    pointer :: weiib(:)
  end type elm

  type elmgp
     real(rp),    pointer :: gpvol(:)
     real(rp),    pointer :: gpcar(:,:,:)
     real(rp),    pointer :: gphes(:,:,:)
     real(rp),    pointer :: hleng(:)
     real(rp),    pointer :: tragl(:,:)
  end type elmgp

  type elm_cloud
     integer(ip)          :: pgaus            ! Number of Gauss points
     real(rp),    pointer :: shape(:,:)       ! Shape function (pnode,pgaus)
     real(rp),    pointer :: deriv(:,:,:)     ! Derivative (ndime,pnode,pgaus)
     real(rp),    pointer :: heslo(:,:,:)     ! Hessian (ntens,pnode,pgaus)
     real(rp),    pointer :: weigp(:)         ! Gauss points weights
     real(rp),    pointer :: posgp(:,:)       ! Gauss points positions in canonical element
  end type elm_cloud
  
  !----------------------------------------------------------------------
  !
  ! Mesh
  !
  !----------------------------------------------------------------------

  type cell

     real(rp)      :: coor(3,2)
     integer(ip)   :: neigh(6)
     integer(ip)   :: level
     integer(ip)   :: marked
     real(rp)      :: rsize

  end type cell
  
   type typ_lobas
     real(rp),   pointer :: local_basis(:,:,:)
     integer(ip)         :: type
     real(rp)            :: param(9)
  end type typ_lobas

  !----------------------------------------------------------------------
  !
  ! Element bin
  !
  !----------------------------------------------------------------------

  type typ_element_bin
     real(rp)             :: comin(3)
     real(rp)             :: comax(3)
     integer(ip)          :: boxes(3)
     integer(ip), pointer :: bin_size(:,:,:)
     type(i1p),   pointer :: list_elements(:,:,:)
  end type typ_element_bin

  !----------------------------------------------------------------------
  !
  ! OMPSS_DOMAIN
  !
  !----------------------------------------------------------------------

  type ompss_domain
     integer(ip), pointer :: neighbours(:)
     integer(ip), pointer :: elements(:)
     integer(ip)          :: neighIdx
     integer(ip)          :: elemIdx
  end type ompss_domain

  !------------------------------------------------------------------------
  !
  ! KD-TREE
  !
  !------------------------------------------------------------------------

  type netyp
     integer(ip)          :: nelem             ! Number of elements for a node
     integer(ip), pointer :: eleme(:)          ! Ids of elements
     integer(ip), pointer :: ltype(:)          ! Types of elements
  end type netyp
  
end module def_kintyp_domain
!> @}
