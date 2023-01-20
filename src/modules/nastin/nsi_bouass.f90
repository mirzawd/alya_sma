!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_bouass(&
     pevat,ndofn,pnodb,lboel,gbsha,gbsur,tract,wmatr,&
     wrhsi,elmat,elrhs)
  !-----------------------------------------------------------------------
  !****f* Nastin/nsi_bouass
  ! NAME 
  !    nsi_bouass
  ! DESCRIPTION
  ! USES
  ! USED BY
  !    nsi_matrix
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only     :  ip,rp
  use def_domain, only     :  ndime
  implicit none
  integer(ip), intent(in)  :: pevat,ndofn,pnodb
  integer(ip), intent(in)  :: lboel(pnodb)
  real(rp),    intent(in)  :: wmatr(pevat,pevat),wrhsi(pevat)
  real(rp),    intent(in)  :: gbsha(pnodb),gbsur,tract(ndime)
  real(rp),    intent(inout) :: elmat(pevat,pevat),elrhs(pevat)
  integer(ip)              :: ievat,jevat,inodb,idime,ievab

  do ievat = 1,pevat
     do jevat = 1,pevat
        elmat(jevat,ievat) = elmat(jevat,ievat) + gbsur * wmatr(jevat,ievat)
     end do
     elrhs(ievat) = elrhs(ievat) + gbsur * wrhsi(ievat)
  end do

  do inodb = 1,pnodb
     ievab = (lboel(inodb)-1) * ndofn
     do idime = 1,ndime
        ievab = ievab + 1
        elrhs(ievab) = elrhs(ievab) + gbsur * tract(idime) * gbsha(inodb)
     end do
  end do

end subroutine nsi_bouass
