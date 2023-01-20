!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine par_bounda()
!-------------------------------------------------------------------------------
!****f* parall/par_bounda
! NAME
!    par_bounda
! DESCRIPTION
!    
! INPUT
!    lnpar_par 
!    xadj 
!    adj
! OUTPUT
!    part
! USED BY
!    par_arrays
!***
!-------------------------------------------------------------------------------
  use def_kintyp
  use def_domain
  use def_master
  use def_parall
  implicit none
  integer(ip) :: iboun, ielem, domai

  nboun_par(1:npart_par) = 0
  nboun_total = 0
  if(kfl_bouel==1) then
     do iboun= 1, nboun
        ielem            = lelbo(iboun)
        if( ielem == 0 ) then
           call runend('PAR_BOUNDA: YOU HAVE NOT PROVIDED THE BOUNDARY/ELEMENT CONNECTIVITY FILE (LELBO)')
        end if
        domai            = lepar_par(ielem)
        lbpar_par(iboun) = domai
        nboun_par(domai) = nboun_par(domai) + 1
     enddo
     do domai = 1,npart_par
        nboun_total = nboun_total + nboun_par(domai)
     end do
  else if(nboun>0) then
     call runend('PARALL: MUST GIVE BOUNDARY/ELEMENT CONNECTIVITY')
  end if

end subroutine par_bounda
