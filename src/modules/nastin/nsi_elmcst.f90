!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_elmcst(&
     pnode,pgaus,plapl,gplap,gpsha,gpcar,gpvol,gpvis,elmac)
  !-----------------------------------------------------------------------
  !****f* Nastin/nsi_elmcst
  ! NAME 
  !    nsi_elmcst
  ! DESCRIPTION
  !    Compute constant element matrix
  ! USES
  ! USED BY
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only     :  ip,rp
  use def_domain, only     :  ndime,mnode
  implicit none
  integer(ip), intent(in)  :: pnode,pgaus,plapl
  real(rp),    intent(in)  :: gpsha(pnode,pgaus)
  real(rp),    intent(in)  :: gpcar(ndime,mnode,pgaus),gpvol(pgaus)
  real(rp),    intent(in)  :: gpvis(pgaus),gplap(pnode,pgaus)
  real(rp),    intent(out) :: elmac(pnode,pnode)
  integer(ip)              :: igaus,inode,jnode
  real(rp)                 :: fact1,fact2
  !
  ! Initialization
  ! 
  do inode=1,pnode
     do jnode=1,pnode
        elmac(jnode,inode)=0.0_rp
     end do
  end do
  !
  ! ( mu*dui/dxk , dv/dxk ) 
  !
  if(ndime==2) then
     do igaus=1,pgaus
        fact2=gpvis(igaus)*gpvol(igaus)
        do inode=1,pnode
           do jnode=inode+1,pnode
              fact1=gpcar(1,inode,igaus)*gpcar(1,jnode,igaus)&
                   +gpcar(2,inode,igaus)*gpcar(2,jnode,igaus)
              fact1=fact1*fact2
              elmac(inode,jnode)=elmac(inode,jnode)+fact1
              elmac(jnode,inode)=elmac(jnode,inode)+fact1
           end do
           elmac(inode,inode)=elmac(inode,inode)&
                & +fact2*( gpcar(1,inode,igaus)*gpcar(1,inode,igaus)&  
                &         +gpcar(2,inode,igaus)*gpcar(2,inode,igaus))  
        end do
     end do
  else
     do igaus=1,pgaus
        fact2=gpvis(igaus)*gpvol(igaus)
        do inode=1,pnode
           do jnode=inode+1,pnode
              fact1=gpcar(1,inode,igaus)*gpcar(1,jnode,igaus)&
                   +gpcar(2,inode,igaus)*gpcar(2,jnode,igaus)&
                   +gpcar(3,inode,igaus)*gpcar(3,jnode,igaus)
              fact1=fact1*fact2
              elmac(inode,jnode)=elmac(inode,jnode)+fact1
              elmac(jnode,inode)=elmac(jnode,inode)+fact1
           end do
           elmac(inode,inode)=elmac(inode,inode)&
                & +fact2*(gpcar(1,inode,igaus)*gpcar(1,inode,igaus)&  
                &        +gpcar(2,inode,igaus)*gpcar(2,inode,igaus)&  
                &        +gpcar(3,inode,igaus)*gpcar(3,inode,igaus))  
        end do
     end do
  end if
  ! 
  ! Viscosity: ( mu*d^2ui/dxk^2, vj )
  !
  if(plapl==1) then
     do igaus=1,pgaus
        do inode=1,pnode 
           fact1=gpsha(inode,igaus)*gpvol(igaus)
           do jnode=1,pnode
              elmac(inode,jnode)=elmac(inode,jnode)&
                   +gplap(jnode,igaus)*fact1
           end do
        end do
     end do
  end if

end subroutine nsi_elmcst
