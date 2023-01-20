!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_boumas(&
     pevat,pnodb,ndofn,lboel,gbsha,&
     gbden,baloc,wmatr)

  !-----------------------------------------------------------------------
  !****f* Nastin/nsi_boumas
  ! NAME 
  !    nsi_boumas
  ! DESCRIPTION
  !    This routine computes the surface integral that comes when integrating 
  !    by parts the mass equation for low mach
  !    int(rho vel·n qh)
  ! USES
  !   
  ! USED BY
  !    nsi_bouope
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only     :  ip,rp             
  use def_domain, only     :  ndime
  implicit none
  integer(ip), intent(in)  :: pnodb,pevat,ndofn
  integer(ip), intent(in)  :: lboel(pnodb)
  real(rp),    intent(in)  :: gbden
  real(rp),    intent(in)  :: gbsha(pnodb)
  real(rp),    intent(in)  :: baloc(ndime,ndime)
  real(rp),    intent(out) :: wmatr(pevat,pevat)
  integer(ip)              :: inodb,jnodb,jdime,ievab,jevab
  real(rp)                 :: fact0

  do inodb=1, pnodb
     ievab = (lboel(inodb) -1)*ndofn + ndime +1
     do jnodb = 1, pnodb
        fact0 = gbsha(inodb)*gbsha(jnodb)*gbden 
        jevab = (lboel(jnodb) -1)*ndofn
        do jdime =1, ndime    
           jevab = jevab +1
           wmatr(ievab, jevab) = wmatr (ievab, jevab) + fact0*baloc(jdime, ndime)
        end do
     end do
  end do

end subroutine nsi_boumas
