!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    sld_elmope_vect_expl.f90
!> @author  Adria Quintanas-Corominas
!> @date    November 2021
!> @brief   Vectorised elemental loop
!> @details Vectorised elemental loop 
!> @}
!------------------------------------------------------------------

subroutine sld_vect_elmope_expl(itask)

  use def_kintyp,                  only : ip,rp
  use def_domain,                  only : lnnod, lelch
  use def_elmtyp,                  only : ELINT, ELFEM
  use def_domain,                  only : lgaus
  use mod_parall,                  only : num_subd_par
  use mod_parall,                  only : num_pack_par
  use mod_parall,                  only : list_elements_par
  use mod_parall,                  only : typ_list_elements_par
  use mod_sld_vect_assembly_elfem, only : sld_vect_element_operations_ELFEM_expl
  use mod_sld_vect_assembly_elint, only : sld_vect_element_operations_ELINT

  implicit none

  integer(ip), intent(in)              :: itask

  integer(ip)                          :: isubd,ipack,ielem
  integer(ip)                          :: pnode,pgaus,pelch
  integer(ip)                          :: num_subd
  integer(ip),                 pointer :: num_pack(:)
  type(typ_list_elements_par), pointer :: list_elements(:)
  real(rp)                             :: time_detail(10)

  num_subd      =  num_subd_par
  num_pack      => num_pack_par
  list_elements => list_elements_par

  time_detail = 0.0_rp
  !
  ! Elemental assembly 
  ! 
  subds: do isubd = 1, num_subd
     packs: do ipack = 1, num_pack(isubd)

        ielem = list_elements(isubd) % packs(ipack) % l(1)  ! Select first element
        pnode = lnnod(ielem)                                ! Number of nodes
        pgaus = lgaus(ielem)                                ! Number of Gauss points
        pelch = lelch(ielem)

        ! It is assumed that the elements are also agruped according to the lelch
        if(     lelch(ielem) == ELFEM ) then
           ! Standard elements
           call sld_vect_element_operations_ELFEM_expl(itask,pnode,pgaus,list_elements(isubd) % packs(ipack) % l)

        elseif( lelch(ielem) == ELINT ) then
           ! Interface elements
           call sld_vect_element_operations_ELINT(itask,pnode,list_elements(isubd) % packs(ipack) % l)

        endif

     enddo packs
  enddo subds

end subroutine sld_vect_elmope_expl
