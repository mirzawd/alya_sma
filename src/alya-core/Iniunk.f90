!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine Iniunk()
  !-----------------------------------------------------------------------
  !****f* master/Iniunk
  ! NAME
  !    Iniunk
  ! DESCRIPTION
  !    This routine ask the modules to compute their initial condition
  ! USES
  !    moduls
  ! USED BY
  !    Alya
  !***
  !-----------------------------------------------------------------------
  
  use def_kintyp,          only : ip,rp
  use def_master,          only : modul
  use def_master,          only : iblok
  use def_master,          only : nblok
  use def_master,          only : mitim
  use def_master,          only : kfl_gotim
  use def_master,          only : cutim
  use def_master,          only : lun_outpu
  use def_master,          only : ID_KERMOD
  use def_master,          only : ITASK_BEFORE
  use def_master,          only : ITASK_AFTER
  use def_master,          only : ITASK_INIUNK
  use mod_coupling_driver, only : COU_DRIVER
  use mod_communications,  only : PAR_BARRIER
  use mod_outfor,          only : outfor
  use mod_messages,        only : messages_live
  use mod_mass_matrix,     only : mass_matrix_consistent
  use mod_moduls,          only : moduls 
  implicit none
  !
  ! Call modules
  !
  modul = ID_KERMOD
  call COU_DRIVER(ITASK_BEFORE,ITASK_INIUNK)
  call Kermod(-ITASK_INIUNK)
  modul = ID_KERMOD
  call COU_DRIVER(ITASK_AFTER,ITASK_INIUNK)
  
  modul = 0 

  call messages_live('INITIAL SOLUTION','START SECTION')

  do iblok = 1,nblok
     call moduls(ITASK_INIUNK)
  end do
  !
  ! No time step is performed
  !
  if( mitim == 0 )  then
     kfl_gotim = 0
     cutim     = 0.0_rp
  end if
  call messages_live('INITIAL SOLUTION','END SECTION') 
  !
  ! Start iterating!
  !
  call outfor(91_ip,lun_outpu,' ')

end subroutine Iniunk
