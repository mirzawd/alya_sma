!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Levels
!> @{
!> @file    lev_output_interface_mesh.f90
!> @author  houzeaux
!> @date    2020-09-07
!> @brief   Interface mesh
!> @details Output interface mesh
!> @} 
!-----------------------------------------------------------------------

subroutine lev_output_interface_mesh()

  use def_kintyp_basic,      only : ip,rp
  use def_kintyp_mesh_basic, only : mesh_type_basic
  use mod_mesh_type_basic,   only : mesh_type_basic_output
  use mod_mesh_type_basic,   only : mesh_type_basic_parallel
  use def_master,            only : ittim
  use def_master,            only : fleve
  use def_master,            only : mem_modul
  use def_master,            only : modul
  use def_master,            only : kfl_paral
  use def_kermod,            only : ndivi
  use def_domain,            only : meshe
  use mod_parall,            only : PAR_COMM_MY_CODE
  use mod_strings,           only : integer_to_string
  use def_levels,            only : npp_inter_lev
  use mod_communications,    only : PAR_SUM
  use mod_messages,          only : messages_live
  
  implicit none

  type(mesh_type_basic)         :: mesh
  real(rp),             pointer :: dista(:)
  integer(ip)                   :: nelem_sum
  
  if( npp_inter_lev /= 0 ) then
     if( mod(ittim,npp_inter_lev) == 0 ) then
        !
        ! Create mesh
        !
        dista => fleve(:,1)
        call mesh % init()
        call mesh % cut_level(meshe(ndivi),dista,mem_modul(1:2,modul))

        nelem_sum = mesh % nelem
        call PAR_SUM(nelem_sum)
        if( nelem_sum > 0 ) then
           call messages_live('LEVELS: OUTPUT INTERFACE MESH, # ELEMENTS= '//integer_to_string(nelem_sum))
           call mesh_type_basic_parallel(mesh,PAR_COMM_MY_CODE,kfl_paral)
           mesh % name = 'LEVEL-'//integer_to_string(ittim)
           call mesh_type_basic_output(mesh)
        end if
        !
        ! Deallocate 
        !
        call mesh % deallo()
     end if
  end if

end subroutine lev_output_interface_mesh
