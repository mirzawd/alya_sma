!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine chm_bouope()
  !-------------------------------------------------------------------------
  !****f* emic/chm_bouope
  ! NAME
  !    chm_bouope
  ! DESCRIPTION
  !    Boundary operations
  ! USES
  ! USED BY
  !    chm_matrix
  !***
  !------------------------------------------------------------------------
  use def_master,     only : INOTMASTER
  use def_domain,     only : nboun
  use def_kintyp,     only : ip
  use def_chemic,     only : iclaf_chm, iclai_chm, kfl_robin_chm, kfl_fixbo_chm
  implicit none
  integer(ip) :: iboun,kfl_robin
  integer(ip) :: iclas

  external    :: runend

  kfl_robin = 0

  if( INOTMASTER .and. kfl_robin_chm /= 0 ) then
     iboun     = 0
     do while( iboun < nboun )
        iboun = iboun + 1
        iclas = iclai_chm - 1
        do while( iclas < iclaf_chm )
           iclas = iclas + 1
           if( kfl_fixbo_chm(iclas,iboun) == 2 ) then
              kfl_robin = 1
              iclas     = iclaf_chm
              iboun     = nboun
           end if
        end do
     end do

     if( kfl_robin == 1 ) then
       call runend('ROBIN-TYPE BCs ARE NOT CODED IN CHEMIC')
     end if

  else  !! Master

     if( kfl_robin == 1 ) then
       call runend('ROBIN-TYPE BCs ARE NOT CODED IN CHEMIC')
     end if

  end if

end subroutine chm_bouope
