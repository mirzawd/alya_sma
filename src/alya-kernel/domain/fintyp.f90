!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine fintyp(ndime,nnode,ielty)
  !-----------------------------------------------------------------------
  !****f* Domain/fintyp
  ! NAME
  !    fintyp
  ! DESCRIPTION
  !    This routine defines the following derivated parameters:
  ! USED BY
  !    Domain
  !***
  !-----------------------------------------------------------------------

  ! This is VERY dangerous!!! For 3D elements, many have identical number
  ! of nodes, like TET20 and HEX20. We need to include topology as input as well.

  use      def_elmtyp
  use mod_outfor, only : outfor

  implicit none
  integer(ip), intent(in)  :: ndime,nnode
  integer(ip), intent(out) :: ielty

  ielty=0
  if( ndime == 0 ) then
     ielty=POINT
  else if( ndime == 1 ) then
     if(nnode == 2) then
        ielty = BAR02
     else if(nnode == 3) then
        ielty = BAR03
     else if(nnode == 4) then
        ielty = BAR04
     end if
  else if(ndime == 2) then
     if(nnode == 3) then
        ielty = TRI03
     else if(nnode == 4) then
        ielty = QUA04
     else if(nnode == 6) then
        ielty = TRI06
     else if(nnode == 8) then
        ielty = QUA08
     else if(nnode == 9) then
        ielty = QUA09 
     else if(nnode == 10) then
        ielty = TRI10
     else if(nnode == 16) then
        ielty = QUA16 
     end if
  else
     if(nnode == 3) then
        ielty = SHELL
     else if(nnode == 4) then
        ielty = TET04
     else if(nnode == 5) then
        ielty = PYR05
     else if(nnode == 6) then
        ielty = PEN06
     else if(nnode == 8) then
        ielty = HEX08
     else if(nnode == 10) then
        ielty = TET10
     else if(nnode == 14) then
        ielty = PYR14
     else if(nnode == 15) then
        ielty = PEN15
     else if(nnode == 18) then
        ielty = PEN18
     else if(nnode == 20) then
        !ielty = HEX20
        ielty = TET20
     else if(nnode == 27) then
        ielty = HEX27
     end if
  end if
  if(ielty == 0) call outfor(13_ip,nnode,' ')

end subroutine fintyp
