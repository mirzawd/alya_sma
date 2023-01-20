!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Domain
!> @{
!> @file    domain_destructor.f90
!> @author  houzeaux
!> @date    2019-06-11
!> @brief   Destroy domain
!> @details Reinitialize all domain arrays
!> @} 
!-----------------------------------------------------------------------

subroutine domain_destructor()

  use def_kintyp
  use def_master
  use def_domain
  use def_kermod,            only : ielse
  use mod_elsest,            only : elsest_deallocate
  use mod_htable,            only : htades
  use mod_memory,            only : memory_deallo
  use mod_memory,            only : memory_resize, memory_alloca, memory_deallo
  use mod_par_bin_structure, only : par_bin_structure_deallocate
  use mod_mesh_type,         only : mesh_type_deallocate
  use mod_redistribution,    only : redistribution_domain,redistribution_initialization
  use mod_mesh_type,         only : mesh_type_update_last_mesh
  use mod_messages,          only : messages_live
  use mod_output,            only : output_file_names
  use mod_iofile,            only : iofile_open_unit
  use mod_iofile,            only : iofile_close_unit
  use mod_graphs,            only : graphs_elepoi_deallocate
  use mod_domain,            only : domain_memory_deallocate
  use mod_renumbering_nodes, only : permR_nodes
  use mod_std
  implicit none
  !
  ! Others
  !
  if( INOTMASTER ) then
     call elsest_deallocate(ielse)
  end if
  call par_bin_structure_deallocate()
  !
  ! Domain
  !
  call graphs_elepoi_deallocate(pelpo,  lelpo,  PELPO_NAME='PELPO',  LELPO_NAME='LELPO',memor=memor_dom)
  call graphs_elepoi_deallocate(pelpo_2,lelpo_2,PELPO_NAME='PELPO_2',LELPO_NAME='LELPO_2',memor=memor_dom)

  call memory_deallo(memor_dom,'C_DOM'      ,'domain_destructor' , c_dom       )
  call memory_deallo(memor_dom,'R_DOM'      ,'domain_destructor' , r_dom       )
  call memory_deallo(memor_dom,'LPOTY'      ,'domain_destructor' , lpoty       )
  call memory_deallo(memor_dom,'LNLEV'      ,'domain_destructor' , lnlev       )
  call memory_deallo(memor_dom,'LELEV'      ,'domain_destructor' , lelev       )
  call memory_deallo(memor_dom,'LBLEV'      ,'domain_destructor' , lblev       )
  call memory_deallo(memor_dom,'PERMR_NODES','domain_destructor' , permR_nodes )
  
  call domain_memory_deallocate('LNNOD')
  call domain_memory_deallocate('LNNOB')
  call domain_memory_deallocate('LGAUS')
  
  call domain_memory_deallocate('LPERI')
  
  call domain_memory_deallocate('LNSEC')
  call domain_memory_deallocate('LESEC')
  call domain_memory_deallocate('LBSEC')

  call domain_memory_deallocate('LBONO')
  !
  ! Hash table for global numbering
  !
  call htades( htable_lninv_loc )
  !
  ! Solver
  !
  nzdom     = 0
  nzdom_own = 0
  nzmat     = 0                 
  nzrhs     = 0
  nzpre     = 0
  !
  ! No halo computed yet!
  !
  npoin_2   = npoin
  nelem_2   = nelem 
  nboun_2   = nboun
  !
  ! Current module is kernel
  !
  modul = 0
  call moddef(9_ip) 

end subroutine domain_destructor
