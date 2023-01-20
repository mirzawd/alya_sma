!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine par_senbcs()
  !------------------------------------------------------------------------
  !****f* Parall/par_senbcs
  ! NAME
  !    par_senbcs
  ! DESCRIPTION
  !    Send boundary conditions to slaves  
  ! OUTPUT
  ! USED BY
  !    Domain
  !***
  !------------------------------------------------------------------------
  use def_kintyp
  use def_parame
  use def_parall
  use def_domain
  use def_master
  use mod_memchk
  use mod_domain, only : domain_memory_allocate
  use mod_domain, only : domain_memory_deallocate
  implicit none  

  strin = 'par_senbcs'
  strre = 'par_senbcs'
  strch = 'par_senbcs'

  !----------------------------------------------------------------------
  !
  ! Codes
  !
  !----------------------------------------------------------------------
  !
  ! Slaves allocate memory
  !
  if( ISLAVE ) then
     nboun_2 = nboun
     call domain_memory_allocate('KFL_CODNO')
     call domain_memory_allocate('KFL_CODBO')
  end if
  !
  ! Gather code arrays
  !
  if( kfl_icodn > 0 ) call par_parari('GAT',NPOIN_TYPE,mcono*npoin,kfl_codno)
  if( kfl_icodb > 0 ) call par_parari('GAT',NBOUN_TYPE,nboun,kfl_codbo)
  !
  ! Master deallocates memory
  !
  if( IMASTER ) then
     call domain_memory_deallocate('KFL_CODNO')
     call domain_memory_deallocate('KFL_CODBO')
  end if

end subroutine par_senbcs
