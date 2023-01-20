!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tur_sasimu(ndime,pnode,gptur,gpden,gradk,gradw,gppr2,sasso)
  !------------------------------------------------------------------------
  !****f* Turbul/tur_sstkom
  ! NAME
  !   tur_sasimu
  ! DESCRIPTION
  !    Computes the source term for the scale adaptive simulation method
  !    (SST model)
  ! OUTPUT 
  !    SASSO
  ! USES
  ! USED BY
  !    tur_sstkom
  !***
  !------------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp
  use def_turbul, only       :  nturb_tur
  implicit none
  integer(ip), intent(in)    :: ndime,pnode
  real(rp),    intent(in)    :: gptur(nturb_tur)
  real(rp),    intent(out)   :: sasso
  real(rp),    intent(in)    :: gpden,gppr2
  real(rp),    intent(in)    :: gradk(3),gradw(3)
  integer(ip)                :: idime,inode
  real(rp)                   :: gpkin,gpome,kappa,zeta2,csas,sigphi
  real(rp)                   :: u2dash,grk,grw,gplap
  real(rp)                   :: lensca,vklen,lratio,argk,argw,kograd
  real(rp)                   :: comp1,comp2,compsum

  real(rp)                   :: elvel(3,pnode),gphes(3,pnode)
  integer(ip)                :: plapl

  ! turbulence variables
  gpkin = gptur(1)
  gpome = gptur(2)

  kappa  = 0.41_rp
  zeta2  = 3.51_rp
  csas   = 2.0_rp
  sigphi = 2.0_rp/3.0_rp

  u2dash = 0.0_rp
  grk    = 0.0_rp
  grw    = 0.0_rp

  if( plapl == 1 ) then
     if( ndime == 1 ) then
        do inode=1,pnode
           gplap=&
                gphes(1,inode)
           u2dash = u2dash + sqrt(              &
                (gplap*elvel(1,inode))**2.0_rp)
        end do
     else if( ndime == 2 ) then
        do inode=1,pnode
           gplap=&
                gphes(1,inode)&
                +gphes(2,inode)
           u2dash = u2dash + sqrt(              &
                (gplap*elvel(1,inode))**2.0_rp &
                +(gplap*elvel(2,inode))**2.0_rp)
        end do
     else
        do inode=1,pnode
           gplap=&
                gphes(1,inode)&
                +gphes(2,inode)&
                +gphes(3,inode)
           u2dash = u2dash + sqrt(              &
                (gplap*elvel(1,inode))**2.0_rp &
                +(gplap*elvel(2,inode))**2.0_rp &
                +(gplap*elvel(3,inode))**2.0_rp)
        end do
     end if
  else
     gplap=0.0_rp
     u2dash=0.0_rp
  end if

  do idime = 1,ndime
    grk = grk + gradk(idime)
    grw = grw + gradw(idime)
  end do

  lensca = sqrt(gpkin)/( (0.09_rp**0.25_rp) *gpome )
  vklen = kappa*2.0_rp*gppr2/u2dash
  lratio = lensca/vklen
  argk = 1.0_rp/( gpkin*gpkin )*grk*grk
  argw = 1.0_rp/( gpome*gpome )*grw*grw
  kograd = max( argk , argw )
  comp1 = gpden*zeta2*kappa*2.0_rp*gppr2*(lratio)**2.0_rp
  comp2 = csas*2.0_rp*gpden*gpkin/sigphi*kograd
  compsum = comp1-comp2
  sasso = max( compsum , 0.0_rp )

end subroutine tur_sasimu
