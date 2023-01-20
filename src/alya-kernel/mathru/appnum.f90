!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine appnum(inume,filna)
  !-----------------------------------------------------------------------
  !****f* outrut/appnum
  ! NAME
  !    appnum
  ! DESCRIPTION
  !    Append a number to a file name:
  !    name.rst <= name-00000**.rst
  ! USES
  ! USED BY
  !***
  !-----------------------------------------------------------------------
  use def_kintyp
  use def_parame
  use def_master
  use def_domain
  use mod_iofile
  use def_postpr
  implicit none
  integer(ip),    intent(in)    :: inume
  character(150), intent(inout) :: filna
  character(7)                  :: wnume
  character(6)                  :: wexte
  integer(ip)                   :: ileng

  if( inume < 10 ) then
     write(wnume,'(a,i1)') '000000',inume
  else if( inume < 100 ) then
     write(wnume,'(a,i2)') '00000',inume
  else if( inume < 1000 ) then
     write(wnume,'(a,i3)') '0000',inume
  else if( inume < 10000 ) then
     write(wnume,'(a,i4)') '000',inume
  else if( inume < 100000 ) then
     write(wnume,'(a,i5)') '00',inume
  else if( inume < 1000000 ) then
     write(wnume,'(a,i6)') '0',inume
  else if( inume < 10000000 ) then
     write(wnume,'(i7)') inume
  else
     call runend('APPNUM: COULD NOT GENERATE FILE NAME')
  end if
  !
  ! Append number before extension
  !
  ileng = len(adjustl(trim(filna)))

  if( trim(filna(ileng-2:ileng-2)) == '.' ) then

     wexte = adjustl(trim(filna(ileng-2:ileng)))
     filna = adjustl(trim(filna(1:ileng-3)))//'-'//wnume//trim(wexte)

  else if( trim(filna(ileng-3:ileng-3)) == '.' ) then

     wexte = adjustl(trim(filna(ileng-3:ileng)))
     filna = adjustl(trim(filna(1:ileng-4)))//'-'//wnume//trim(wexte)

  else if( trim(filna(ileng-4:ileng-4)) == '.' ) then

     wexte = adjustl(trim(filna(ileng-4:ileng)))
     filna = adjustl(trim(filna(1:ileng-5)))//'-'//wnume//trim(wexte)

  else if( trim(filna(ileng-6:ileng-6)) == '.' ) then

     wexte = adjustl(trim(filna(ileng-6:ileng)))
     filna = adjustl(trim(filna(1:ileng-7)))//'-'//wnume//trim(wexte)

  end if

end subroutine appnum
