!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine par_senprt_rea(ndim1,rvari)

  use def_master
  implicit none
  integer(ip), intent(in) :: ndim1
  real(rp),    target     :: rvari(ndim1)

  nparr =  ndim1
  parre => rvari

  if( PART_AND_WRITE() ) then
     call par_sendin()
  else
     call par_receiv()
  end if

  nparr = 0
  nullify(parre)
  
end subroutine par_senprt_rea
