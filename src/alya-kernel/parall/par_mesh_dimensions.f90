!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Parall
!> @{
!> @file    par_mesh_dimensions.f90
!> @author  houzeaux
!> @date    2019-06-13
!> @brief   Mesh dimensions
!> @details Mesh dimensions
!> @} 
!-----------------------------------------------------------------------

subroutine par_mesh_dimensions()

  use def_kintyp,         only : ip    
  use def_master,         only : IMASTER
  use def_master,         only : ISEQUEN
  use def_master,         only : npoin_par
  use def_master,         only : nelem_par
  use def_master,         only : nboun_par
  use def_master,         only : npart
  use def_master,         only : npoi1
  use def_master,         only : npoi2
  use def_master,         only : npoi3
  use def_domain,         only : npoin
  use def_domain,         only : nelem
  use def_domain,         only : nboun
  use def_domain,         only : npoin_origi
  use def_domain,         only : npoin_total
  use def_domain,         only : nelem_total
  use def_domain,         only : nboun_total
  use mod_communications, only : PAR_GATHER
  use mod_communications, only : PAR_SUM
  use mod_parall,         only : PAR_CODE_SIZE
  use mod_parall,         only : par_memor
  use mod_memory,         only : memory_alloca

  implicit none

  integer(ip), pointer :: npoin_par0(:)
  integer(ip), pointer :: nelem_par0(:)
  integer(ip), pointer :: nboun_par0(:)

  nullify(npoin_par0,nelem_par0,nboun_par0)
  
  if( ISEQUEN ) then
     
     npoin_origi = npoin
     npoin_total = npoin 
     nelem_total = nelem 
     nboun_total = nboun

  else

     npoin_origi = npoi1+(npoi3-npoi2)+1 ! Cannot use npoin_own because mesh has not been renumbered
     npoin_total = npoin 
     nelem_total = nelem 
     nboun_total = nboun
     call PAR_SUM(npoin_origi)
     call PAR_SUM(npoin_total)
     call PAR_SUM(nelem_total)
     call PAR_SUM(nboun_total)

     if( IMASTER ) then
        allocate(npoin_par0(0:PAR_CODE_SIZE-1))
        allocate(nelem_par0(0:PAR_CODE_SIZE-1))
        allocate(nboun_par0(0:PAR_CODE_SIZE-1))
        if( .not. associated(npoin_par) ) call memory_alloca(par_memor,'NPOIN_PAR' ,'par_interface_exchange' , npoin_par  , npart )
        if( .not. associated(nelem_par) ) call memory_alloca(par_memor,'NELEM_PAR' ,'par_interface_exchange' , nelem_par  , npart )
        if( .not. associated(nboun_par) ) call memory_alloca(par_memor,'NBOUN_PAR' ,'par_interface_exchange' , nboun_par  , npart )
     end if
     call PAR_GATHER(npoin,npoin_par0)
     call PAR_GATHER(nelem,nelem_par0)
     call PAR_GATHER(nboun,nboun_par0) 
     if( IMASTER ) then
        npoin_par(1:PAR_CODE_SIZE-1) = npoin_par0(1:PAR_CODE_SIZE-1) 
        nelem_par(1:PAR_CODE_SIZE-1) = nelem_par0(1:PAR_CODE_SIZE-1) 
        nboun_par(1:PAR_CODE_SIZE-1) = nboun_par0(1:PAR_CODE_SIZE-1) 
        deallocate(npoin_par0) 
        deallocate(nelem_par0)
        deallocate(nboun_par0)
     end if
  end if

end subroutine par_mesh_dimensions
