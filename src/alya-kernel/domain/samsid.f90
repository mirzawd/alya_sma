!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine samsid(orige,desti,poin1,poin2,ndime,istru)
  !-----------------------------------------------------------------------
  ! NAME
  !    samsid
  ! DESCRIPTION
  !    Determine if a point is inside a triangle using the same side technique
  ! USED BY
  !    insiib
  !***
  !----------------------------------------------------------------------- 

  use def_kintyp
  
  implicit   none
  integer(ip), intent(in)   :: ndime
  integer(ip), intent(out)  :: istru
  integer(ip)               :: idime
  real(rp),    intent(in)   :: poin1(ndime),poin2(ndime),orige(ndime),desti(ndime)
  real(rp)                  :: side1(3),side2(3),side3(3),cp1(3),cp2(3),resul
  
  side1(3) = 0.0_rp
  side2(3) = 0.0_rp 
  side3(3) = 0.0_rp 
  do idime=1,ndime
     side1(idime) = desti(idime) - orige(idime)
     side2(idime) = poin1(idime) - orige(idime)
     side3(idime) = poin2(idime) - orige(idime)
  end do
  
  call vecpro(side1,side2,cp1,3)
  call vecpro(side1,side3,cp2,3)
  
  resul=0.0_rp
  do idime = 1,ndime
     resul = resul + (cp1(idime) * cp2(idime))
  end do
  if (resul >= 0) then
     istru = 1
  else
     istru = 0
  end if
  
end subroutine samsid
