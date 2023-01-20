!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine par_inidat()
  !------------------------------------------------------------------------
  !****f* Parall/par_inidat
  ! NAME
  !    par_inidat
  ! DESCRIPTION
  !    Send initial data to slaves using MPI
  ! OUTPUT
  ! USED BY
  !    Parall
  !***
  !------------------------------------------------------------------------
  use def_parame
  use def_master
  use def_parall
  use def_inpout
  use mod_communications, only : PAR_BROADCAST
  implicit none
  integer(ip), target :: kfl_varia_tmp(6)

  if( ISEQUEN .and. kfl_ptask == 0 ) then
     !
     ! Enable Master to perform preprocess
     !
     kfl_paral=0
     call vocabu(-1_ip,0_ip,0_ip) 

  else if( ( IMASTER .and. ( kfl_ptask /= 0 .or. nproc_par > 1 ) ) .or. ISLAVE ) then
     !
     ! Exchange type of restart file (ASCII/binary), file hierarchy
     ! and file prefix
     !
     call vocabu(-1_ip,0_ip,0_ip) 
     kfl_varia_tmp(1) =  kfl_ascii_par
     kfl_varia_tmp(2) =  kfl_fileh_par
     kfl_varia_tmp(3) =  nsire_par
     kfl_varia_tmp(5) =  kfl_async
     kfl_varia_tmp(6) =  npart
     call PAR_BROADCAST(6_ip,kfl_varia_tmp)
     kfl_ascii_par    =  kfl_varia_tmp(1)
     kfl_fileh_par    =  kfl_varia_tmp(2)
     nsire_par        =  kfl_varia_tmp(3)
     kfl_async        =  kfl_varia_tmp(5) 
     npart            =  kfl_varia_tmp(6) 
     call vocabu(-1_ip,0_ip,0_ip)
  end if

end subroutine par_inidat
 
