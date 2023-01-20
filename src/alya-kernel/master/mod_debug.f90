!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



module mod_debug

  use def_kintyp
  use mod_par_tools, only : PAR_GLOBAL_TO_LOCAL_NODE
  use def_master,    only : INOTMASTER,IMASTER
  use def_master,    only : kfl_paral,lninv_loc,igene
  use def_domain,    only : coord
  use mod_communications
  implicit none

  character(5) :: wopos_tmp(3)

contains

  subroutine debug_interface()

    use def_domain, only : npoin
    use def_master, only : npoi1,npoi2,npoi3

    implicit none
    integer(ip)          :: ipoin,ierro1,ierro2
    integer(ip), pointer :: itest(:)

    ierro1 = 0
    ierro2 = 0
    nullify(itest)

    if( INOTMASTER ) then

       allocate(itest(npoin))
       itest = 0
       do ipoin = 1,npoi1
          itest(ipoin) = 1
       end do
       do ipoin = 1,npoi2,npoi3
          itest(ipoin) = 1
       end do

       call PAR_INTERFACE_NODE_EXCHANGE(coord,'DIFF')
       call PAR_INTERFACE_NODE_EXCHANGE(itest,'SUM')

       ierro1 = minval(itest)
       ierro2 = maxval(itest)

       deallocate(itest)

    end if

    call PAR_MIN(ierro1)
    call PAR_MAX(ierro2)

    if( ierro1 /= 0 ) then
       if( IMASTER ) write(*,*) 'MIN VALUE= ',ierro1
       call runend('ERROR WITH INTERFACE')
    end if

    if( ierro2 /= 1 ) then
       if( IMASTER ) write(*,*) 'MAX VALUE= ',ierro2
       call runend('ERROR WITH INTERFACE')
    end if

  end subroutine debug_interface

end module mod_debug
