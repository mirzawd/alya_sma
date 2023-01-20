!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_bouopp(&
     lboel,lnodb,gbsha,baloc,wmatr,pnode,pnodb,ndofn,&
     pevat,shapp)
  !-----------------------------------------------------------------------
  !****f* Nastin/nsi_bouopp
  ! NAME 
  !    nsi_bouopp
  ! DESCRIPTION
  !    This routine computes the contribution to WMATR for the Navier-
  !    Stokes equations due to open boundaries. In this case, the term
  !    S.n (S = -p I + 2 mu Sym(grad(u)) being the Cauchy stress tensor) 
  !    is considered unknown.
  ! USES
  ! USED BY
  !    nsi_bouope 
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp
  use def_master, only       :  tempe
  use def_kermod, only       :  gasco
  use def_domain, only       :  ndime
  use def_nastin, only       :  kfl_regim_nsi
  implicit none
  integer(ip), intent(in)    :: pnode,pnodb,pevat,ndofn
  integer(ip), intent(in)    :: lnodb(pnodb),lboel(pnodb)
  real(rp),    intent(in)    :: gbsha(pnodb),shapp(pnode)
  real(rp),    intent(in)    :: baloc(ndime)
  real(rp),    intent(inout) :: wmatr(pevat,pevat)
  integer(ip)                :: inodb,ldime,ievab,jnods,jevat,ipoin
  real(rp)                   :: xmuit,prod1,gbtem

  if(kfl_regim_nsi==1) then
     !
     ! Compressible regime: unknown is density
     !
     gbtem=0.0_rp
     do inodb=1,pnodb
        ipoin=lnodb(inodb)
        gbtem=gbtem+gbsha(inodb)*tempe(ipoin,1)
     end do
     do ldime=1,ndime
        prod1=baloc(ldime)*gasco*gbtem
        do inodb=1,pnodb
           ievab=(lboel(inodb)-1)*ndofn+ldime
           xmuit=gbsha(inodb)*prod1
           do jnods=1,pnode
              jevat=jnods*ndofn
              wmatr(ievab,jevat)=wmatr(ievab,jevat)&
                   +xmuit*shapp(jnods)
           end do
        end do
     end do
     
  else
     !
     ! Contribution from the pressure term     
     !
     do ldime=1,ndime
        prod1=baloc(ldime)
        do inodb=1,pnodb
           ievab=(lboel(inodb)-1)*ndofn+ldime
           xmuit=gbsha(inodb)*prod1
           do jnods=1,pnode
              jevat=jnods*ndofn
              wmatr(ievab,jevat)=wmatr(ievab,jevat)&
                   +xmuit*shapp(jnods)
           end do
        end do
     end do
  end if

end subroutine nsi_bouopp
