!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tem_elmtes(&
     pnode,pgaus,plapl,gpsha,gpcar,gplap,gprea,gpadv,&
     gpgrd,gpcon,gpcod,gpstp,gpstt,ttemp)
  !-----------------------------------------------------------------------
  !****f* Temper/tem_elmtes
  ! NAME
  !   tem_elmtes
  ! DESCRIPTION
  !    Compute the adjoint operator
  !    -[rho*a+grad(k)].grad(v) -k*Lapl(v) -k/r*dv/dr+ r*v
  ! OUTPUT 
  !    GPADG ... Adjoint operator at Gauss point
  ! USES
  ! USED BY
  !    tem_elmope
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only     :  ip,rp 
  use def_domain, only     :  ndime,mnode,kfl_naxis
  use def_temper, only     :  kfl_grdif_tem,kfl_advec_tem,kfl_taust_tem
  implicit none
  integer(ip), intent(in)  :: pnode,pgaus,plapl
  real(rp),    intent(in)  :: gpsha(pnode,pgaus)
  real(rp),    intent(in)  :: gpcar(ndime,mnode,pgaus)
  real(rp),    intent(in)  :: gplap(pnode,pgaus),gpgrd(ndime,pgaus)
  real(rp),    intent(in)  :: gpcon(pgaus),gpcod(ndime,pgaus)
  real(rp),    intent(in)  :: gprea(pgaus)
  real(rp),    intent(in)  :: gpadv(pnode,pgaus)
  real(rp),    intent(in)  :: gpstp(pgaus),gpstt(pgaus)
  real(rp),    intent(out) :: ttemp(pnode,pgaus)
  integer(ip)              :: inode,idime,igaus
  real(rp)                 :: fact1

  if(kfl_taust_tem/=0) then
     !
     ! Source term: (tau^{-1}*tau'-tau'*s)*Ni
     !
     do igaus=1,pgaus
        do inode=1,pnode
           ttemp(inode,igaus)=gpsha(inode,igaus)&
                *(gpstt(igaus)-gpstp(igaus)*gprea(igaus))
        end do
     end do
     !
     ! Advection: tau'*rho*a.grad(Ni), rho=rho*cp
     !
     if(kfl_advec_tem/=0) then
        do igaus=1,pgaus
           do inode=1,pnode
              ttemp(inode,igaus)=ttemp(inode,igaus)&
                   +gpstp(igaus)*gpadv(inode,igaus)
           end do
        end do
     end if
     !
     ! tau'*grad(k).grad(Ni)
     !
     if(kfl_grdif_tem==1) then
        do igaus=1,pgaus
           do inode=1,pnode
              do idime=1,ndime
                 ttemp(inode,igaus)=ttemp(inode,igaus)&
                      +gpstp(igaus)*gpgrd(idime,igaus)&
                      *gpcar(idime,inode,igaus)
              end do
           end do
        end do
     end if
     !
     ! tau'*k*Lapl(Ni)
     !
     if(plapl==1) then     
        do igaus=1,pgaus
           do inode=1,pnode
              ttemp(inode,igaus)=ttemp(inode,igaus)&
                   +gpstp(igaus)*gplap(inode,igaus)
           end do
        end do
     end if
     !
     ! Axisymmetric: tau1'*k/r*dNi/dr
     !
     if(kfl_naxis==1) then     
        do igaus=1,pgaus
           fact1=gpstp(igaus)*gpcon(igaus)/gpcod(1,igaus)
           do inode=1,pnode
              ttemp(inode,igaus)=ttemp(inode,igaus)&
                   +fact1*gpcar(1,inode,igaus)
           end do
        end do
     end if

  end if

end subroutine tem_elmtes
