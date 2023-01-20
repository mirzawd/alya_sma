!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine par_senprt_int(ndim1,ivari)

  use def_master
  implicit none
  integer(ip), intent(in) :: ndim1
  integer(ip), target     :: ivari(ndim1)

  npari =  ndim1
  parin => ivari

  if( PART_AND_WRITE() ) then
     call par_sendin()
  else
     call par_receiv()
  end if

  npari = 0
  nullify(parin)
  
end subroutine par_senprt_int
