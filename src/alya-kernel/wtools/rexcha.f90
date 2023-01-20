!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine rexcha(reava)
  !-----------------------------------------------------------------------
  !****f* rexcha
  ! NAME
  !    rexcha
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
  real(rp) :: reava

  nparr=nparr+1

  if(parii==2) then
     if(kfl_ptask==1) then
        if(kfl_paral<=0) parre(nparr) = reava
        if(kfl_paral>=1) reava        = parre(nparr)
     else if(kfl_ptask==0) then
        if(kfl_paral<=0) parre(nparr) = reava
     else if(kfl_ptask==2) then
        reava = parre(nparr) 
     end if
  end if

end subroutine rexcha

