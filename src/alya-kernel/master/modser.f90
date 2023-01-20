!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine modser
  !------------------------------------------------------------------------
  !****f* master/outerr
  ! NAME 
  !    outerr
  ! DESCRIPTION
  !    This routine checks the order of modules iterations, and compute
  !    number of modules and services
  ! OUTPUT
  !    LMORD(IORDE,IBLOK) ... Module number to run at IORDE position
  !                           within block IBLOK
  !    NMODU ................ Number of modules
  !    NSERV ................ Number of services
  ! USES
  ! USED BY
  !    Reapro
  !***
  !------------------------------------------------------------------------
  use def_kintyp
  use def_master
  implicit none
  integer(ip) :: imodu,iorde
  !
  ! Automatic ordering of the modules: LMORD
  !
  if( lmord(1,1) == 0 ) then
     nblok = 1     
     iorde = 0
     do imodu = 1,mmodu-1
        if( kfl_modul(imodu) /= 0 ) then
           iorde = iorde + 1
           lmord(iorde,1) = imodu
        end if
     end do
  else
     
  end if
  !
  ! Count the number of modules used: NMODU
  !
  nmodu = 0
  do imodu = 1,mmodu
     if( kfl_modul(imodu) /= 0 ) nmodu = nmodu + 1
  end do

end subroutine modser
