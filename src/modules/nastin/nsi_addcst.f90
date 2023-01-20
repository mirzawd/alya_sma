!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_addcst(ndofn,amatr_nsi,amatr)
  !-----------------------------------------------------------------------
  !****f* Nastin/nsi_addcst
  ! NAME 
  !    nsi_addcst
  ! DESCRIPTION
  !    This routine adds the constant matrix to the momentum
  !    global matrix
  ! USES
  ! USED BY
  !    nsi_matrix
  !***
  !-----------------------------------------------------------------------
  use def_kintyp,   only     :  ip,rp
  use def_domain,   only     :  npoin,ndime,r_dom
  use def_nastin,   only     :  kfl_savco_nsi,kfl_fixno_nsi
  implicit none
  integer(ip), intent(in)    :: ndofn
  real(rp),    intent(inout) :: amatr_nsi(*),amatr(ndofn,ndofn,*)
  integer(ip)                :: ipoin,izdom,idime

  if(kfl_savco_nsi==1) then
     !
     ! Impose boundary conditions on constant matrix
     !     
     do ipoin=1,npoin
        if(kfl_fixno_nsi(1,ipoin)/=0) then
           do izdom=r_dom(ipoin),r_dom(ipoin+1)-1
              amatr_nsi(izdom)=0.0_rp
           end do
        end if
     end do
     
  end if

  if(kfl_savco_nsi>=1) then
     !
     ! Add constant matrix
     !
     kfl_savco_nsi=2     
     do ipoin=1,npoin
        do izdom=r_dom(ipoin),r_dom(ipoin+1)-1
           do idime=1,ndime
              amatr(idime,idime,izdom)=amatr(idime,idime,izdom)&
                   +amatr_nsi(izdom)
           end do
        end do
     end do
  end if

end subroutine nsi_addcst
