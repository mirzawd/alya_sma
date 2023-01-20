!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



#ifdef CMAKE

subroutine infrev()
  use def_master
  use gitinfo
  implicit none
  character(150) :: log_file
  integer, parameter :: u=1457845

  log_file = adjustl(trim(namda))//"-git.log"

  open(u, file=log_file, action="write")
  write(u,*) "Remote URL: "//trim(GIT_REMOTE)
  write(u,*) "Revision: "//trim(GIT_REVISION)
  write(u,*) "Branch: "//trim(GIT_BRANCH)
  close(u)

end subroutine infrev

#endif
