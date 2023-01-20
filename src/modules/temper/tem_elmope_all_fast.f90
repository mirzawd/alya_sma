!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tem_elmope_all_fast(order)
  !------------------------------------------------------------------------
  !****f* Temper/tem_elmop2
  ! NAME 
  !    tem_elmop2
  ! DESCRIPTION
  !    ORDER=1:
  !      Temperature equation, elemental operations:
  !      1. Compute elemental matrix and RHS 
  !      2. Impose Dirichlet boundary conditions
  !      3. Assemble them
  !    ORDER=4:
  !      Update the subgrid scale
  ! USES
  ! USED BY
  !    tem_matrix
  !------------------------------------------------------------------------

#include "def_vector_size.inc"
  use def_master
  use def_domain
  use def_temper
  use def_kermod
  !##
  use mod_parall,                 only : num_subd_par
  use mod_parall,                 only : num_pack_par
  use mod_parall,                 only : list_elements_par
  use mod_parall,                 only : typ_list_elements_par
  use mod_tem_fast

  implicit none

  integer(ip), intent(in) :: order                     ! =2: compute SGS only

  integer(ip) :: ielem                     ! Indices and dimensions
  integer(ip) :: pnode
  integer(ip) :: pgaus

  real(rp)    :: dtmin
  integer(ip) :: isubd, ipack, num_subd
  integer(ip),                 pointer :: num_pack(:)
  type(typ_list_elements_par), pointer :: list_elements(:)
  integer(ip) :: VECTOR_DIM

#ifdef EVENT
  call mpitrace_user_function(1)
#endif

  VECTOR_DIM    =  VECTOR_SIZE
  num_subd      =  num_subd_par
  num_pack      => num_pack_par
  list_elements => list_elements_par

  subds: do isubd = 1,num_subd
    packs: do ipack = 1,num_pack(isubd)
        !
        ! Element dimensions
        !
        ielem = list_elements(isubd) % packs(ipack) % l(1)  ! Select first element
        pnode = lnnod(ielem)                                ! Number of nodes
        pgaus = lgaus(ielem)                                ! Number of Gauss points

        call mod_tem_element_operations_fast(VECTOR_DIM,pnode,pgaus,list_elements(isubd) % packs(ipack) % l,order)

    end do packs
  end do subds

  if( kfl_reset == 0 ) then
     if( dtmin /= 0.0_rp ) then
        dtcri_tem = min(dtcri_tem,dtmin)
     end if
  end if

#ifdef EVENT
  call mpitrace_user_function(0)
#endif

end subroutine tem_elmope_all_fast

