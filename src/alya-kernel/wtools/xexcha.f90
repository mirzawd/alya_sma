!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine xexcha(cpxva)
  !-----------------------------------------------------------------------
  !****f* xexcha
  ! NAME
  !    iexcha
  ! DESCRIPTION
  !    This routine exchange complex data individually.
  ! USES
  ! USED BY
  !    nsi_sendat
  !    nsa_sendat
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  implicit none
  complex(rp) :: cpxva

  nparx=nparx+1

  if(parii==2) then
     if(kfl_ptask==1) then
        if(kfl_paral<=0) parcx(nparx) = cpxva
        if(kfl_paral>=1) cpxva        = parcx(nparx)
     else if(kfl_ptask==0) then
        if(kfl_paral<=0) parcx(nparx) = cpxva
     else if(kfl_ptask==2) then
        cpxva = parcx(nparx) 
     end if
  end if

end subroutine xexcha
