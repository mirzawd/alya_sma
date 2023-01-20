!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine errors(itask,ierro,iwarn,messa)
  !------------------------------------------------------------------------
  !****f* master/errors
  ! NAME 
  !    errors
  ! DESCRIPTION
  !    ITASK == 0 ... Error from master
  !    ITASK == 1 ... From master
  !    ITASK == 2 ... From master
  !    ITASK == 2 ... From module
  ! USES
  ! USED BY
  !***
  !------------------------------------------------------------------------
  use def_master
  use def_domain
  use def_solver
  use mod_communications, only : PAR_BROADCAST
  use mod_outfor,         only : outfor
  implicit none
  integer(ip),  intent(in)    :: itask
  integer(ip),  intent(inout) :: ierro
  integer(ip),  intent(inout) :: iwarn
  character(*), intent(in)    :: messa
  integer(ip)                 :: kfl_ptask_old
  character(20)               :: werro
  
  kfl_ptask_old = kfl_ptask
  kfl_ptask     = 1

  if( itask == 3 ) then

     !-------------------------------------------------------------------
     !
     ! Called by modules to treat errors and/or warnings
     !
     !-------------------------------------------------------------------

     call PAR_BROADCAST(ierro)
     call PAR_BROADCAST(iwarn)
     !
     ! Warning
     !
     if( iwarn /= 0 ) call outfor( 3_ip,momod(modul)%lun_outpu,' ')
     !
     ! Stop
     !
     werro = intost(ierro)
     if( ierro /= 0 ) call outfor(-4_ip,momod(modul)%lun_outpu,trim(werro))

  else

     !-------------------------------------------------------------------
     !
     ! Errors detected by Kernel
     !
     !-------------------------------------------------------------------

     call PAR_BROADCAST(ierro)

     if( ierro /= 0 ) then
        if( itask == 1 ) then        
           call outfor(  4_ip,lun_outpu,trim(messa))
        else if( itask == 2 ) then        
           call outfor(-47_ip,lun_outpu,trim(messa))
        end if
     end if
  end if

  kfl_ptask = kfl_ptask_old 

end subroutine errors
