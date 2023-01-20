!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine iexcha(intva)
  !-----------------------------------------------------------------------
  !****f* iexcha
  ! NAME
  !    iexcha
  ! DESCRIPTION
  !    This routine exchange integer data individually
  ! USES
  ! USED BY
  !    nsi_sendat
  !    nsa_sendat
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  implicit none
  integer(ip) :: intva

  npari=npari+1

  if(parii==2) then
     if(kfl_ptask==1) then
        if(kfl_paral<=0) parin(npari) = intva
        if(kfl_paral>=1) intva        = parin(npari)
     else if(kfl_ptask==0) then
        if(kfl_paral<=0) parin(npari) = intva
     else if(kfl_ptask==2) then
        intva = parin(npari) 
     end if
  end if

end subroutine iexcha
