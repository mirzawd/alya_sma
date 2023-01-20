!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine lev_velfun(&
     kfl_advec_lev,ndime,pnode,lnods,coord,elvel)
  !-----------------------------------------------------------------------
  !****f* Levels/lev_velfun
  ! NAME
  !   lev_velfun
  ! DESCRIPTION
  !   Compute elvelity according to the function number 
  ! INPUT
  !   KFL_ADVEC_LEV ... Function
  !   COORD ........... Coordinates
  !   NDIME ........... Dimension
  !   PNODE ........... Number of nodes or Gauss points
  ! OUTPUT 
  !   ELVEL ........... Elvelity
  ! USES
  ! USED BY
  !    lev_elmope
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only     :  ip,rp
  use def_parame, only     :  pi 
  use def_master, only     :  veloc, ittim, dtinv
  implicit none
  integer(ip), intent(in)  :: kfl_advec_lev,ndime,pnode
  integer(ip), intent(in)  :: lnods(*)
  real(rp),    intent(in)  :: coord(ndime,pnode)
  real(rp),    intent(out) :: elvel(ndime,pnode)
  integer(ip)              :: ipoin,inode,idime
  real(rp)                 :: x,y,z,theta,p,t

  if(kfl_advec_lev==1) then
     do inode=1,pnode
        ipoin=lnods(inode)
        do idime=1,ndime
           elvel(idime,inode)=veloc(idime,ipoin,1)
        end do
     end do

  else if(kfl_advec_lev==2) then
     do inode=1,pnode
        x=coord(1,inode)
        y=coord(2,inode)
        elvel(1,inode)= 0.5_rp*(1.0_rp-x*x)*(1.0_rp+y)
        elvel(2,inode)=-0.5_rp*x*(4.0_rp-(1.0_rp+y)**2)
     end do

  else if(kfl_advec_lev==3) then
     do inode=1,pnode
        elvel(1,inode) =  1.0_rp
        elvel(2,inode) = -1.0_rp   
     end do

  else if(kfl_advec_lev==4) then
     do inode=1,pnode
        elvel(1,inode) = 1.0_rp
        elvel(2,inode) = 0.0_rp
     end do

  else if(kfl_advec_lev==5) then
     do inode=1,pnode
        theta            =  1.1071487  ! tan-1(2)
        elvel(1,inode) =  cos(theta)
        elvel(2,inode) = -sin(theta)
     end do

     !
     !  Zalesak test velocity field
     !
  else if(kfl_advec_lev==6) then
     do inode=1,pnode
        x=coord(1,inode)
        y=coord(2,inode)
        elvel(1,inode) = (0.5_rp-y)
        elvel(2,inode) = (x-0.5_rp)
     end do

  else if(kfl_advec_lev==7) then
     p=3.141592653589793_rp
     do inode=1,pnode
        x=coord(1,inode)
        y=coord(2,inode)
        elvel(1,inode)=-sin(p*x)*sin(p*x)*2*sin(p*y)*cos(p*y)
        elvel(2,inode)=sin(p*y)*sin(p*y)*2*sin(p*x)*cos(p*x)
     end do

  else if(kfl_advec_lev==8) then
     do inode=1,pnode
        x=coord(1,inode)
        y=coord(2,inode)
        elvel(1,inode)=(50_rp-y)/100_rp
        elvel(2,inode)=(x-50_rp)/100_rp
        elvel(3,inode)=0_rp
     end do

  else if(kfl_advec_lev==9) then
     do inode=1,pnode
        x=coord(1,inode)
        y=coord(2,inode)
        z=coord(3,inode)
        elvel(1,inode)=(50_rp-y)/100_rp
        elvel(2,inode)=(x-50_rp)/100_rp
        elvel(3,inode)=0_rp
     end do

  else if(kfl_advec_lev==10) then

     p=3.141592653589793_rp

     if(ittim==0) then
        t=ittim/dtinv
     else
        t=(ittim-1)/dtinv
     endif

     do inode=1,pnode

        x = coord(1,inode)
        y = coord(2,inode)
        z = coord(3,inode)
        elvel(1,inode)= 2*sin(pi*x)*sin(pi*x)*sin(2*pi*y)*sin(2*pi*z)*cos(pi*t/3.0)
        elvel(2,inode)= -sin(2*pi*x)*sin(pi*y)*sin(pi*y)*sin(2*pi*z)*cos(pi*t/3.0)
        elvel(3,inode)= -sin(2*pi*x)*sin(2*pi*y)*sin(pi*z)*sin(pi*z)*cos(pi*t/3.0)

     end do


  end if

end subroutine lev_velfun
