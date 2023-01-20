!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine ker_biflu2(&
     itask,pgaus,pnode,gpsha,gpcar,ellev,prope_wat,prope_air,&
     thicl,hleng,gppro,grpro)
  !-----------------------------------------------------------------------
  !****f* kermod/ker_biflu2
  ! NAME 
  !    ker_biflu2
  ! DESCRIPTION
  !    Compute a property using the bifluid model - the 2 indicates it includes de gradient
  ! USES
  !
  ! USED BY
  !    modules
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp  
  use def_parame, only       :  pi
  use def_domain, only       :  ndime,mnode
  implicit none
  integer(ip),  intent(in)   :: itask
  integer(ip),  intent(in)   :: pgaus
  integer(ip),  intent(in)   :: pnode
  real(rp),     intent(in)   :: gpsha(pnode,pgaus)
  real(rp),     intent(in)   :: gpcar(ndime,mnode,pgaus)
  real(rp),     intent(in)   :: ellev(pnode)
  real(rp),     intent(in)   :: prope_wat
  real(rp),     intent(in)   :: prope_air
  real(rp),     intent(in)   :: thicl
  real(rp),     intent(in)   :: hleng
  real(rp),     intent(out)  :: gppro(pgaus)
  real(rp),     intent(out)  :: grpro(ndime,pgaus)
  integer(ip)                :: igaus,inode,idime
  real(rp)                   :: phi,eps,f,grphi(3),ooeps,oovpi
 

  if( thicl < 1.0e-10_rp ) then

     do igaus = 1,pgaus
        phi = 0.0_rp
        do inode = 1,pnode
           phi = phi + ellev(inode) * gpsha(inode,igaus)
        end do
        if( phi >= 0.0_rp ) then
           gppro(igaus) = prope_wat
           grpro(1:ndime,igaus) = 0.0_rp
        else
           gppro(igaus) = prope_air
           grpro(1:ndime,igaus) = 0.0_rp
        end if
     end do


else

  if( thicl > 0.0_rp ) then
     eps = thicl
  else
     eps = -thicl * hleng
  end if

  ooeps = 1.0_rp/eps
  oovpi = 1.0_rp/pi

  do igaus = 1,pgaus
     phi = 0.0_rp
     do inode = 1,pnode
        phi = phi + ellev(inode) * gpsha(inode,igaus)
     end do

     if( phi < -eps ) then
        gppro(igaus)         = prope_air
        grpro(1:ndime,igaus) = 0.0_rp
     else if( phi > eps ) then
        gppro(igaus)         = prope_wat
        grpro(1:ndime,igaus) = 0.0_rp
     else 
        f = 0.5_rp*(1.0_rp+phi*ooeps+sin(pi*phi*ooeps)*oovpi)
        gppro(igaus) = prope_air + (prope_wat-prope_air)*f
        if( itask == 1 ) then
           grphi = 0.0_rp
           do inode = 1,pnode
              do idime = 1,ndime
                 grphi(idime) = grphi(idime) + ellev(inode) * gpcar(idime,inode,igaus)
              end do
           end do
           do idime = 1,ndime
              grpro(idime,igaus) = (prope_wat-prope_air)*0.5_rp*(grphi(idime)*ooeps)*(1.0_rp+cos(pi*phi*ooeps))
           end do
        end if
     end if

  end do
end if
end subroutine ker_biflu2
