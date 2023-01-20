!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine velchl(pnode,elcod,elvel,chale,hleng)
  !-----------------------------------------------------------------------
  !****f* Domain/elmchl
  ! NAME
  !   elmchl
  ! DESCRIPTION
  !   Compute the element length in the flow direction
  ! OUTPUT
  !   CHALE
  ! USES
  ! USED BY
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only     :  ip,rp
  use def_domain, only     :  ndime
  implicit none
  integer(ip), intent(in)  :: pnode
  real(rp),    intent(in)  :: elvel(ndime,pnode),elcod(ndime,pnode)
  real(rp),    intent(in)  :: hleng(ndime)
  real(rp),    intent(out) :: chale(2)
  integer(ip)              :: idime,inode,jnode,knode
  real(rp)                 :: velcg(3),coocg(3),rnode
  real(rp)                 :: vx1,vy1,vx2,vy2,m,x1,y1,x2,y2
  real(rp)                 :: xmin1,xmin2,ymin1,ymin2,det,venor

  rnode=real(pnode,rp)
  do idime=1,ndime
     velcg(idime)=0.0_rp
     coocg(idime)=0.0_rp
     do inode=1,pnode
        velcg(idime)=velcg(idime)+elvel(idime,inode)
        coocg(idime)=coocg(idime)+elcod(idime,inode)
     end do
     velcg(idime)=velcg(idime)/rnode
     coocg(idime)=coocg(idime)/rnode
  end do
  venor=0.0_rp
  do idime=1,ndime
     venor=venor+velcg(idime)*velcg(idime)
  end do
  if(venor<1.0e-12_rp) then
     chale(1)=hleng(ndime)
     chale(2)=hleng(ndime)
  else
     if(ndime==2) then
        inode=0
        knode=0
        do while(inode<pnode)
           !
           ! Length in flow direction
           !
           inode=inode+1
           if(inode==pnode) then
              jnode=1
           else
              jnode=inode+1
           end if
           xmin1 = coocg(1)
           ymin1 = coocg(2)
           vx1   = velcg(1)
           vy1   = velcg(2)
           xmin2 = elcod(1,inode)
           ymin2 = elcod(2,inode)
           vx2   = elcod(1,jnode)-elcod(1,inode)
           vy2   = elcod(2,jnode)-elcod(2,inode)
           det   = vy2*vx1-vy1*vx2
           if(abs(det)>0.0_rp) then
              m = (vx1*(ymin1-ymin2)+vy1*(xmin2-xmin1))/det
              if(m>=0.0_rp.and.m<=1.0_rp) then
                 if(knode==0) then
                    x1    = xmin2+m*vx2
                    y1    = ymin2+m*vy2
                    knode = 1
                 else
                    x2    = xmin2+m*vx2
                    y2    = ymin2+m*vy2
                    knode = 2
                    inode = pnode
                 end if
              end if
           end if
        end do
        if(knode==2) then
           chale(1)=sqrt( (x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) )
        else
           call runend('VELCHL: COULD NOT COMPUTE DISTANCE')
        end if

        inode=0
        knode=0
        do while(inode<pnode)
           !
           ! Length in orthogonal direction to flow
           !
           inode=inode+1
           if(inode==pnode) then
              jnode=1
           else
              jnode=inode+1
           end if
           vx1   = -velcg(2)
           vy1   =  velcg(1)
           xmin2 =  elcod(1,inode)
           ymin2 =  elcod(2,inode)
           vx2   =  elcod(1,jnode)-elcod(1,inode)
           vy2   =  elcod(2,jnode)-elcod(2,inode)
           det   =  vy2*vx1-vy1*vx2
           if(abs(det)>0.0_rp) then
              m = (vx1*(ymin1-ymin2)+vy1*(xmin2-xmin1))/det
              if(m>=0.0_rp.and.m<=1.0_rp) then
                 if(knode==0) then
                    x1    = xmin2+m*vx2
                    y1    = ymin2+m*vy2
                    knode = 1
                 else
                    x2    = xmin2+m*vx2
                    y2    = ymin2+m*vy2
                    knode = 2
                    inode = pnode
                 end if
              end if
           end if
        end do
        if(knode==2) then
           chale(2)=sqrt( (x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) )
        else
           call runend('VELCHL: COULD NOT COMPUTE DISTANCE')
        end if
     else
        call runend('VELCHL: NOT AVAILABLE IN 3D')
     end if
  end if

end subroutine velchl
