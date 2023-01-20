!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



module mod_log

  use def_kintyp,            only : ip,rp
  use def_master,            only : lun_outpu  
  use mod_iofile,            only : iofile_flush_unit
  
  implicit none

  public :: log_runend_error, log_end_of_analysis

contains

  subroutine log_runend_error(message)
    character(len=*), intent(in) :: message
    integer(4) :: lun_outpu4
    lun_outpu4 = int(lun_outpu,4)
    write(lun_outpu4,100) trim(message)
    call iofile_flush_unit(lun_outpu4)

100   format(//,&
         & 5x,//,&
         & 5x,'--------------------------------------------------------',/,&
         & 5x,/,&
         & 5x,'|- AN ERROR HAS BEEN DETECTED:',/,/,&
         & 9x,a,/)
  end subroutine log_runend_error
  
  subroutine log_end_of_analysis()
    integer(4) :: lun_outpu4
    lun_outpu4 = int(lun_outpu,4)
    write(lun_outpu4,101)
    call iofile_flush_unit(lun_outpu4)

101 format(//,&
         & 5x,'|- END OF ANALYSIS',/)
  end subroutine log_end_of_analysis

end module mod_log
