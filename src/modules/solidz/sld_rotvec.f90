!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine sld_rotvec(irotype,vecold,baslo,vecnew,ndime)
!**************************************************************************
!
!**** This subroutine performs rotations and inverse rotations on vecold,
!**** according to this:  vecnew(i)= baslo(i,j)^(irotype) * vecold(j).
!
!**************************************************************************
  use      def_kintyp  
  implicit none
  integer(ip)             :: irotype,ndime,irots,idime
  real(rp),    intent(in) :: baslo(ndime,ndime)
  real(rp)                :: vecnew(ndime),vecold(ndime),vecaux(3)
  
  if (irotype.eq.-1) then
     !...  vecaux= baslo * vecold
     do irots=1,ndime
        vecaux(irots)=0.0_rp
        do idime=1,ndime
           vecaux(irots)= vecaux(irots)+vecold(idime)*baslo(irots,idime)
        end do
     end do
  else if (irotype.eq.1) then
     !...  vecaux= baslo^(-1) * vecold
     do irots=1,ndime
        vecaux(irots)=0.0_rp
        do idime=1,ndime
           vecaux(irots)= vecaux(irots)+vecold(idime)*baslo(idime,irots)
        end do
     end do
  end if
  
  !...  vecnew <-- vecaux
  do idime=1,ndime
     vecnew(idime)=vecaux(idime)
  end do
  
end subroutine sld_rotvec
