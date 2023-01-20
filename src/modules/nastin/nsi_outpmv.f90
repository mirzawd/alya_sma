!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_outpmv(itask,istar,istop,value)
  !-----------------------------------------------------------------------
  !****f* Nastin/nsi_outpmv
  ! NAME 
  !    nsi_outpmv
  ! DESCRIPTION
  !    This routine outputs PMV and PPD
  ! USES
  ! USED BY
  !    nsi_turnon
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only     :  ip,rp
  use def_master, only     :  veloc,tempe
  use def_domain, only     :  ndime
  use def_nastin, only     :  cloth_nsi,metab_nsi,wetme_nsi,ambie_nsi,&
       &                      radia_nsi,relat_nsi
  implicit none
  integer(ip), intent(in)  :: itask,istar,istop
  real(rp),    intent(out) :: value(*)
  real(rp)                 :: velno,PMV,PPD,tempe_nsi
  integer(ip)              :: idime,istat,ipoin


  do ipoin=istar,istop          
     if(ambie_nsi==-999.0_rp) then
        tempe_nsi = tempe(ipoin,1)
     else
        tempe_nsi = ambie_nsi
     end if
     velno=0.0_rp
     do idime=1,ndime
        velno=velno+veloc(idime,ipoin,1)*veloc(idime,ipoin,1)
     end do
     velno=sqrt(velno)
     call nsi_pmvppd(&
          cloth_nsi,metab_nsi,wetme_nsi,tempe_nsi,radia_nsi,&
          velno,relat_nsi,PMV,PPD,istat)
     if(itask==-1) then
        value(1)     = PMV
     else if(itask==1) then
        value(ipoin) = PMV
     else if(itask==-2) then
        value(1)     = PPD
     else if(itask==2) then
        value(ipoin) = PPD
     end if
  end do

end subroutine nsi_outpmv
