!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine lev_iniun0()
  !-----------------------------------------------------------------------
  !****f* wavequ/lev_iniun0
  ! NAME 
  !    lev_iniun0
  ! DESCRIPTION
  !    This routine does part of what was done in lev_iniunk  but it is done before 
  !    so that ker_updpro(itask_iniunk) has fleve to be used in bifluid law
  ! USES
  !    
  ! USED BY
  !    lev_turnon
  !***
  !-----------------------------------------------------------------------
  use def_parame,           only :  ip,rp, pi
  use def_inpout
  use def_master,           only :  veloc,fleve,kfl_paral,INOTMASTER
  use def_domain
  use def_levels
  use mod_chktyp,           only : check_type
  use mod_chktyp,           only : check_type

  implicit none
  real(rp)                 :: fx
  real(rp)                 :: x,y,z,xc,yc,zc,rc,distx,disty,distz
  real(rp)                 :: xi,xj,a,b,c,d
  real(rp)                 :: vx,vy,vz
  real(rp)                 :: zmed,ztop,zbot
  real(rp)                 :: xmin_jump, xmax_jump, h1, h2, mramp, ylevel, ylevel_max
  integer(ip)              :: kfl_value,ipoin,i,j,inter

  if( INOTMASTER ) then

     !
     ! This part was previously in lev_iniunk - we put it here so that ker_updpro(itask_iniunk) has fleve to calculate 
     ! bifluid props. something similar may need to be done for other modules
     !
     !
     ! Default is air
     !
     do ipoin = 1,npoin
        fleve(ipoin,3) = -1000.0_rp
     end do

     if( kfl_inlev_lev < -1000 ) then
        !
        ! Take initial condition from a value function field: Guillaume para Simone
        ! 
        kfl_value = -kfl_inlev_lev-1000
        call check_type(xfiel,kfl_value,1_ip,npoin) ! Check if value function exist                      
        bvess_lev(1,1:npoin,1) = xfiel(kfl_value) % a(1,1:npoin,1)

        fleve(1:npoin,3)         = bvess_lev(1,1:npoin,1)

     else if ( kfl_inlev_lev == -1 ) then
        !
        ! Prescribed height at z
        !
        do ipoin=1,npoin
           fleve(ipoin,3)     = height_lev - coord(ndime,ipoin) + 1.0e-8_rp
           !fleve(ipoin,3)     = height_lev - coord(ndime,ipoin) 
           bvess_lev(1,ipoin,1) = fleve(ipoin,3) 
        end do


     else if(kfl_inlev_lev==-2) then
        !
        ! Prescribed height at y
        !
        do ipoin=1,npoin
           fleve(ipoin,3)     = height_lev - coord(2,ipoin) + 1.0e-8_rp
           !fleve(ipoin,3)     = height_lev - coord(2,ipoin) 
           bvess_lev(1,ipoin,1) = fleve(ipoin,3) 
        end do

     else if(kfl_inlev_lev==-3) then
        !
        ! Prescribed height at x
        !
        do ipoin=1,npoin
           fleve(ipoin,3)     = height_lev - coord(1,ipoin) + 1.0e-8_rp
           !fleve(ipoin,3)     = height_lev - coord(1,ipoin) 
           bvess_lev(1,ipoin,1) = fleve(ipoin,3) 
        end do

     else if(kfl_inlev_lev==-4) then
        !
        ! Prescribed height at z - inverse
        !
        do ipoin=1,npoin
           fleve(ipoin,3)     = -( height_lev - coord(ndime,ipoin) + 1.0e-8_rp )
           !fleve(ipoin,3)     = -( height_lev - coord(ndime,ipoin) )
           bvess_lev(1,ipoin,1) = fleve(ipoin,3) 
        end do


     else if(kfl_inlev_lev==-5) then
        !
        ! Prescribed height at y - inverse
        !
        do ipoin=1,npoin
           fleve(ipoin,3)     = -( height_lev - coord(2,ipoin) + 1.0e-8_rp )
           !fleve(ipoin,3)     = -( height_lev - coord(2,ipoin) )
           bvess_lev(1,ipoin,1) = fleve(ipoin,3) 
        end do

     else if(kfl_inlev_lev==-6) then
        !
        ! Prescribed height at x - inverse
        !
        do ipoin=1,npoin
           fleve(ipoin,3)     = -( height_lev - coord(1,ipoin) + 1.0e-8_rp )
           !fleve(ipoin,3)     = -( height_lev - coord(1,ipoin) )
           bvess_lev(1,ipoin,1) = fleve(ipoin,3) 
        end do

     else if(kfl_inlev_lev==1) then

        do ipoin=1,npoin
           fx=coord(1,ipoin)
           fleve(ipoin,3)= fx 
        end do

     else if(kfl_inlev_lev==2) then
        !
        !  Zalesak interface initialisation
        !
        xc=0.5_rp
        yc=0.75_rp
        rc=0.15_rp

        do ipoin=1,npoin

           x=coord(1,ipoin)
           y=coord(2,ipoin)

           if(sqrt((x-xc)*(x-xc)+(y-yc)*(y-yc))<=rc) then

              if(y<0.85_rp) then
                 if(x>=0.47_rp.and.x<=0.5_rp) then
                    fleve(ipoin,3)=-x+0.475_rp
                 else if(x<=0.53_rp.and.x>0.5_rp) then
                    fleve(ipoin,3)=x-0.525_rp
                 else 
                    fleve(ipoin,3)= -sqrt((x-xc)*(x-xc)+(y-yc)*(y-yc))+rc
                 endif
              endif

              if(y==0.85_rp) then
                 if(x>=0.48_rp.and.x<=0.52_rp) then
                    fleve(ipoin,3)=0._rp
                 else 
                    fleve(ipoin,3)=-sqrt((x-xc)*(x-xc)+(y-yc)*(y-yc))+rc
                 endif

              else if (y>0.85_rp) then
                 fleve(ipoin,3)=-sqrt((x-xc)*(x-xc)+(y-yc)*(y-yc))+rc
              endif

           else 
              fleve(ipoin,3)= -sqrt((x-xc)*(x-xc)+(y-yc)*(y-yc))+rc
           endif

           fleve(ipoin,3)=fleve(ipoin,3)+1e-08

        end do

     else if(kfl_inlev_lev==3) then

        !
        !  Vortex interface initialisation
        !
        xc=0.5_rp
        yc=0.75_rp
!!$           xc=0.25_rp
!!$           yc=0.5_rp
        rc=0.15_rp

        do ipoin=1,npoin

           x=coord(1,ipoin)
           y=coord(2,ipoin)
!!$              fleve(ipoin,3)= 2.0_rp*(-sqrt((x-xc)*(x-xc)+(y-yc)*(y-yc))+rc)+1e-08
           fleve(ipoin,3)= -sqrt((x-xc)*(x-xc)+(y-yc)*(y-yc))+rc+1e-08
!!$              fleve(ipoin,3)= x*x

        end do

     else if(kfl_inlev_lev==4) then
        !
        !  Rigid Body Rotation of Zalesak's Sphere
        !

        xc=50_rp
        yc=75_rp
        zc=50_rp
        rc=15_rp

        do ipoin=1,npoin

           x=coord(1,ipoin)
           y=coord(2,ipoin)
           z=coord(3,ipoin)

           if(sqrt((x-xc)*(x-xc)+(y-yc)*(y-yc)+(z-zc)*(z-zc))<=rc+1e-02) then

              if(y<85_rp) then
                 if(x>=47_rp.and.x<=50_rp) then
                    fleve(ipoin,3)=-x+47.5_rp
                 else if(x<=53_rp.and.x>50_rp) then
                    fleve(ipoin,3)=x-52.5_rp
                 else 
                    fleve(ipoin,3)= -sqrt((x-xc)*(x-xc)+(y-yc)*(y-yc)+(z-zc)*(z-zc))+rc
                 endif
              endif

              if(y==85_rp) then
                 if(x>=48_rp.and.x<=52_rp) then
                    fleve(ipoin,3)=0._rp
                 else 
                    fleve(ipoin,3)=-sqrt((x-xc)*(x-xc)+(y-yc)*(y-yc)+(z-zc)*(z-zc))+rc
                 endif

              else if (y>85_rp) then
                 fleve(ipoin,3)=-sqrt((x-xc)*(x-xc)+(y-yc)*(y-yc)+(z-zc)*(z-zc))+rc
              endif

           else 
              fleve(ipoin,3)= -sqrt((x-xc)*(x-xc)+(y-yc)*(y-yc)+(z-zc)*(z-zc))+rc
           endif
           fleve(ipoin,3)=fleve(ipoin,3)+1e-08 

        end do

     else if(kfl_inlev_lev==5) then
        !
        !  Three-Dimensional Deformation Field
        !
        xc=0.35_rp
        yc=0.35_rp
        zc=0.35_rp
        rc=0.15_rp


        do ipoin=1,npoin

           x=coord(1,ipoin)
           y=coord(2,ipoin)
           z=coord(3,ipoin)

!!$              fleve(ipoin,3)= 2.0_rp*(-sqrt((x-xc)*(x-xc)+(y-yc)*(y-yc)+(z-zc)*(z-zc))+rc)
           fleve(ipoin,3) = -sqrt((x-xc)*(x-xc)+(y-yc)*(y-yc)+(z-zc)*(z-zc))+rc
           fleve(ipoin,3) = fleve(ipoin,3)+1e-08 

        end do

     else if(kfl_inlev_lev==6) then

        print*,' kikou inlev ',6

        xc=50_rp
        yc=50_rp
        zc=50_rp
        rc=25_rp

        do ipoin=1,npoin

           x=coord(1,ipoin)
           y=coord(2,ipoin)
           z=coord(3,ipoin)

!!$             fleve(ipoin,3)= 2.0_rp*(-sqrt((x-xc)*(x-xc)+(y-yc)*(y-yc)+(z-zc)*(z-zc))+rc)
           fleve(ipoin,3)= -sqrt((x-xc)*(x-xc)+(y-yc)*(y-yc)+(z-zc)*(z-zc))+rc+1e-08

        end do

     else if(kfl_inlev_lev==7) then
        !
        ! Equilibrium interface
        !
        yc=1.0_rp
        do ipoin=1,npoin
           y=coord(2,ipoin)
           x=coord(1,ipoin)
           fleve(ipoin,3) = -(y-yc)+1e-08
           bvess_lev(1,ipoin,1) = -(y-yc)+1e-08
        end do

        if (.not.associated(veloc)) call runend('lev_iniunk: veloc not associated')

        vy=0.0_rp
        do ipoin=1,npoin
           if(kfl_fixno_lev(1,ipoin)==5_ip) then
              vx=0.0_rp
              veloc(1,ipoin,1)=vx
              veloc(1,ipoin,2)=vx
              veloc(1,ipoin,3)=vx
           else
              vx=1.776_rp
              veloc(1,ipoin,1)=vx
              veloc(1,ipoin,2)=vx
              veloc(1,ipoin,3)=vx
           endif
           veloc(2,ipoin,1)=vy
           veloc(2,ipoin,2)=vy
           veloc(2,ipoin,3)=vy
        enddo

     else if(kfl_inlev_lev==8) then

        !
        ! Column interface
        !
        yc=0.5_rp
        xc=0.25_rp
        !yc=0.10_rp
        !xc=0.25_rp
        a=1.0_rp
        b=yc-a*xc

        do ipoin=1,npoin
           y=coord(2,ipoin)
           x=coord(1,ipoin)

           distx=-(x-xc)
           disty=-(y-yc)
           c=a*x+b

           if(y<c) then
              fleve(ipoin,3)= distx+1e-08
           else
              fleve(ipoin,3)= disty+1e-08
           endif

        end do


        if(kfl_paral==-1) then
           ! We need a gauge for this case
           ncapt_lev=0
           do ipoin=1,npoin
              y=coord(2,ipoin)
              if(abs(y)<1e-06) then
                 ncapt_lev=ncapt_lev+1
              endif
           enddo
           allocate(capin_lev(ncapt_lev))
           ncapt_lev=0
           do ipoin=1,npoin
              x=coord(1,ipoin)
              y=coord(2,ipoin)
              if(abs(y)<1e-06) then
                 ncapt_lev=ncapt_lev+1
                 capin_lev(ncapt_lev)=ipoin
              endif
           enddo
           do i=1,ncapt_lev
              xi=coord(1,capin_lev(i))
              do j=i+1,ncapt_lev
                 xj=coord(1,capin_lev(j))
                 if(xj<xi) then
                    inter=capin_lev(i)
                    capin_lev(i)=capin_lev(j)
                    capin_lev(j)=inter
                 endif
              enddo
           enddo
        endif

     else if(kfl_inlev_lev==9) then
        !
        ! Column interface with obstacle
        !
        yc=0.292_rp
        xc=0.146_rp

        do ipoin=1,npoin
           y=coord(2,ipoin)
           x=coord(1,ipoin)

           distx=-(x-xc)
           disty=-(y-yc)
           if(y<x+xc) then
              fleve(ipoin,3)= distx+1e-08
           else
              fleve(ipoin,3)= disty+1e-08
           endif

        end do

     else if(kfl_inlev_lev==10) then
        !
        ! 2D Steady state flat surface with no initial velocity:
        ! to test HS equilibrium:
        !
        yc=0.25_rp

        do ipoin=1,npoin

           y=coord(2,ipoin)

           disty=-(y-yc)
           fleve(ipoin,3)       = disty + 1e-08
           bvess_lev(1,ipoin,1) = disty + 1e-08

        enddo

        vx=0.0_rp
        vy=0.0_rp

        if (.not.associated(veloc)) call runend('lev_iniunk: veloc not associated')

        do ipoin=1,npoin

           veloc(1,ipoin,1)=vx
           veloc(1,ipoin,2)=vx
           veloc(1,ipoin,3)=vx

           veloc(2,ipoin,1)=vy
           veloc(2,ipoin,2)=vy
           veloc(2,ipoin,3)=vy
        enddo


        do ipoin = 1,npoin
           fleve(ipoin,1) = fleve(ipoin,3)
           fleve(ipoin,2) = fleve(ipoin,3)
        end do

     else if(kfl_inlev_lev==11) then

        zc=0.55_rp
        xc=1.992_rp
        a=1.0_rp
        b=zc+xc

        do ipoin=1,npoin
           x=coord(1,ipoin)
           y=coord(2,ipoin)
           z=coord(3,ipoin)

           distx=(x-xc)
           distz=-(z-zc)

           c=-a*x+b

           if(z<c) then
              fleve(ipoin,3)= distx+1e-08
           else
              fleve(ipoin,3)= distz+1e-08
           endif

        end do


     else if(kfl_inlev_lev==12) then

        zc=0.5_rp
        xc=2.97_rp
        a=(zc)/(xc-3.22_rp)
        b=-a*3.22_rp

        do ipoin=1,npoin
           x=coord(1,ipoin)
           y=coord(2,ipoin)
           z=coord(3,ipoin)

           distx=(x-xc)
           distz=-(z-zc)

           if(z<a*x+b) then
              fleve(ipoin,3)= distx+1e-08
           else
              fleve(ipoin,3)= distz+1e-08
           endif


        end do

     else if(kfl_inlev_lev==13) then
        !
        ! Dam Break 3D   
        !
        if( ndime == 2 ) then

           yc=0.5_rp
           xc=0.5_rp
           a=1.0_rp
           b=yc-xc

           do ipoin=1,npoin
              x=coord(1,ipoin)
              y=coord(2,ipoin)

              distx=-(x-xc)
              disty=-(y-yc)

              c=a*x+b

              if(y<c) then
                 fleve(ipoin,3)= distx+1e-08
              else
                 fleve(ipoin,3)= disty+1e-08
              endif

           end do
        else
           zc=0.5_rp
           xc=0.25_rp
           a=1.0_rp
           b=zc-xc

           do ipoin=1,npoin
              x=coord(1,ipoin)
              y=coord(2,ipoin)
              z=coord(3,ipoin)

              distx=-(x-xc)
              distz=-(z-zc)

              c=a*x+b

              if(z<c) then
                 fleve(ipoin,3)= distx+1e-08
              else
                 fleve(ipoin,3)= distz+1e-08
              endif

           end do
        end if

     else if(kfl_inlev_lev==14) then
        !
        ! Dam Break 3D   
        !
        zc=0.55_rp
        xc=1.992_rp
        a=1.0_rp
        b=zc+xc

        do ipoin=1,npoin
           x=coord(1,ipoin)
           y=coord(2,ipoin)
           z=coord(3,ipoin)

           distx=(x-xc)
           distz=-(z-zc)

           c=-a*x+b

           if(z<c) then
              fleve(ipoin,3)= distx+1e-08
           else
              fleve(ipoin,3)= distz+1e-08
           endif

        end do

     else if(kfl_inlev_lev==15) then

        yc=0.55_rp
        xc=1.992_rp
        a=1.0_rp
        b=yc+xc

        do ipoin=1,npoin

           x=coord(1,ipoin)
           y=coord(2,ipoin)

           distx=(x-xc)
           disty=-(y-yc)

           c=-a*x+b

           if(y<c) then
              fleve(ipoin,3)= distx+1e-08
           else
              fleve(ipoin,3)= disty+1e-08
           endif


        end do

     else if(kfl_inlev_lev==17) then
        !
        ! Flow Stick   
        !
        zc=0.1_rp

        do ipoin=1,npoin

           x=coord(1,ipoin)
           y=coord(2,ipoin)
           z=coord(3,ipoin)

           distz=-(z-zc)
           fleve(ipoin,3) = distz+1e-08
           bvess_lev(1,ipoin,1) = distz+1e-08

        enddo

        vy=0.0_rp
        vz=0.0_rp

        if (.not.associated(veloc)) call runend('lev_iniunk: veloc not associated')

        do ipoin=1,npoin
           if(kfl_fixno_lev(1,ipoin)<7_ip) then
              vx=1.0_rp
              veloc(1,ipoin,1)=vx
              veloc(1,ipoin,2)=vx
              veloc(1,ipoin,3)=vx
           else
              vx=0.0_rp
              veloc(1,ipoin,1)=vx
              veloc(1,ipoin,2)=vx
              veloc(1,ipoin,3)=vx
           endif
           veloc(2,ipoin,1)=vy
           veloc(2,ipoin,2)=vy
           veloc(2,ipoin,3)=vy
           veloc(3,ipoin,1)=vz
           veloc(3,ipoin,2)=vz
           veloc(3,ipoin,3)=vz
        enddo


        ! Dam Obs
     else if(kfl_inlev_lev==18) then

        zc=0.55_rp
        xc=1.992_rp
        a=-1.0_rp
        b=zc-a*xc

        do ipoin=1,npoin
           x=coord(1,ipoin)
           y=coord(2,ipoin)
           z=coord(3,ipoin)

           distx=(x-xc)
           distz=-(z-zc)

           c=a*x+b

           if(z<c) then
              fleve(ipoin,3)= distx+1e-08
           else
              fleve(ipoin,3)= distz+1e-08
           endif


        enddo

     else if(kfl_inlev_lev==19) then

        zc=0.0_rp

        do ipoin=1,npoin

           x=coord(1,ipoin)
           y=coord(2,ipoin)
           z=coord(3,ipoin)

           zc = sin(2*x*3.14)*0.2+0.5

           distz=-(z-zc)
           fleve(ipoin,3) = distz+1e-08
           bvess_lev(1,ipoin,1) = distz+1e-08

        enddo

     else if(kfl_inlev_lev==20) then

        zc=0.5_rp
        xc=0.25_rp
        a=1.0_rp
        b=zc-a*xc

        do ipoin=1,npoin
           x=coord(1,ipoin)
           y=coord(2,ipoin)
           z=coord(3,ipoin)

           distx=-(x-xc)
           distz=-(z-zc)

           c=a*x+b

           if(z<c) then
              fleve(ipoin,3)= distx+1e-08
           else
              fleve(ipoin,3)= distz+1e-08
           endif


        enddo

     else if(kfl_inlev_lev==21) then

        !
        !  Vortex interface initialisation
        !
        xc=0.5_rp
        yc=0.75_rp
        rc=0.15_rp

        do ipoin=1,npoin

           x=coord(1,ipoin)
           y=coord(2,ipoin)
           fleve(ipoin,3)= 2.0_rp*(-sqrt((x-xc)*(x-xc)+(y-yc)*(y-yc))+rc)+1e-08

        end do

     else if(kfl_inlev_lev==22) then

        if(kfl_paral==1) print*,' kikou inlev ',22
        !
        !  Vortex interface initialisation
        !
        xc=50_rp
        yc=50_rp
        zc=50_rp
        rc=25_rp

        do ipoin=1,npoin

           x=coord(1,ipoin)
           y=coord(2,ipoin)
           z=coord(3,ipoin)

           fleve(ipoin,3)= 2.0_rp*(-sqrt((x-xc)*(x-xc)+(y-yc)*(y-yc)+(z-zc)*(z-zc))+rc)
!!$             fleve(ipoin,3)= -sqrt((x-xc)*(x-xc)+(y-yc)*(y-yc)+(z-zc)*(z-zc))+rc+1e-08

        end do

     else if(kfl_inlev_lev==23) then
        !
        !  Three-Dimensional Deformation Field 
        !
        xc=0.5_rp
        yc=0.5_rp
        zc=0.5_rp
        rc=0.25_rp

        do ipoin=1,npoin
           x=coord(1,ipoin)
           y=coord(2,ipoin)
           z=coord(3,ipoin)

           fleve(ipoin,3)= 2.0_rp*(-sqrt((x-xc)*(x-xc)+(y-yc)*(y-yc)+(z-zc)*(z-zc))+rc)
           fleve(ipoin,3) = fleve(ipoin,3)+1e-08 
        end do

     else if(kfl_inlev_lev==24) then

        !
        ! Column interface
        !
        yc=3.5_rp
        xc=3.5_rp
        a=1.0_rp
        b=yc-a*xc
        do ipoin=1,npoin
           y=coord(2,ipoin)
           x=coord(1,ipoin)

           distx=-(x-xc)
           disty=-(y-yc)
           c=a*x+b

           if(y<c) then
              fleve(ipoin,3)= distx+1e-08
           else
              fleve(ipoin,3)= disty+1e-08
           endif

        end do


        if(kfl_paral==-1) then
           ! We need a gauge for this case
           ncapt_lev=0
           do ipoin=1,npoin
              y=coord(2,ipoin)
              if(abs(y)<1e-06) then
                 ncapt_lev=ncapt_lev+1
              endif
           enddo
           allocate(capin_lev(ncapt_lev))
           ncapt_lev=0
           do ipoin=1,npoin
              x=coord(1,ipoin)
              y=coord(2,ipoin)
              if(abs(y)<1e-06) then
                 ncapt_lev=ncapt_lev+1
                 capin_lev(ncapt_lev)=ipoin
              endif
           enddo
           do i=1,ncapt_lev
              xi=coord(1,capin_lev(i))
              do j=i+1,ncapt_lev
                 xj=coord(1,capin_lev(j))
                 if(xj<xi) then
                    inter=capin_lev(i)
                    capin_lev(i)=capin_lev(j)
                    capin_lev(j)=inter
                 endif
              enddo
           enddo
        endif

     else if(kfl_inlev_lev==25) then

        !
        !  circle of radius rc arround xc,yc
        !
        xc=0.5_rp
        yc=0.5_rp
        rc=0.2_rp

        do ipoin=1,npoin

           x=coord(1,ipoin)
           y=coord(2,ipoin)
           fleve(ipoin,3)= -sqrt((x-xc)*(x-xc)+(y-yc)*(y-yc))+rc    !+1e-08
           bvess_lev(1,ipoin,1)=fleve(ipoin,3)
        end do

     else if(kfl_inlev_lev==26) then

        !
        !  circle of radius rc arround xc,yc
        !
        xc=0.0_rp
        yc=0.0_rp
        rc=0.0005_rp

        do ipoin=1,npoin

           x=coord(1,ipoin)
           y=coord(2,ipoin)
           fleve(ipoin,3)= -sqrt((x-xc)*(x-xc)+(y-yc)*(y-yc))+rc    !+1e-08
           bvess_lev(1,ipoin,1)=fleve(ipoin,3)
        end do

     else if(kfl_inlev_lev==27) then

        !
        !  for toilet filled with water between ztop and zbot
        !
        ztop=0.26
        zbot=-0.01
        zmed = (ztop+zbot)/2.0_rp
        do ipoin=1,npoin

           z=coord(3,ipoin)
           if(z>zmed) then
              fleve(ipoin,3)= ztop - z
           else
              fleve(ipoin,3)= z - zbot
           end if
           if(coord(2,ipoin) < -0.001_rp) fleve(ipoin,3) = min(fleve(ipoin,3),-0.04_rp) ! modification for the connecting tube hadrien has added

           bvess_lev(1,ipoin,1)=fleve(ipoin,3)
        end do

     else if(kfl_inlev_lev==28) then
        !
        ! Dam Break 3D MARIN
        !
        zc=0.55_rp
        xc=1.228_rp
        a=1.0_rp
        b=zc-xc

        do ipoin=1,npoin
           x=coord(1,ipoin)
           y=coord(2,ipoin)
           z=coord(3,ipoin)

           distx=-(x-xc)
           distz=-(z-zc)

           c=a*x+b

           if(z<c) then
              fleve(ipoin,3)= distx+1e-08
           else
              fleve(ipoin,3)= distz+1e-08
           endif

        end do

     else if(kfl_inlev_lev==29) then

        !
        !  sphere of radius rc arround xc,yc,zc
        !
        xc=0.0_rp
        yc=0.0_rp
        zc=0.0_rp
        rc=0.0005_rp

        do ipoin=1,npoin

           x=coord(1,ipoin)
           y=coord(2,ipoin)
           z=coord(3,ipoin)
           fleve(ipoin,3)= rc - sqrt((x-xc)**2+(y-yc)**2+(z-zc)**2)    !+1e-08
           bvess_lev(1,ipoin,1)=fleve(ipoin,3)
        end do

     else if(kfl_inlev_lev==30) then

        !
        !  Two bubbles of radius rc one on top of each other
        !
        xc=0.5_rp
        yc=0.5_rp
        rc=0.2_rp

        do ipoin=1,npoin

           x=coord(1,ipoin)
           y=coord(2,ipoin)

           if (y .lt. (yc+1.5*rc)) then !lower bubble
              fleve(ipoin,3)= -sqrt((x-xc)*(x-xc)+(y-yc)*(y-yc))+rc    !+1e-08
              bvess_lev(1,ipoin,1)=fleve(ipoin,3)
           else
              fleve(ipoin,3)= -sqrt((x-xc)*(x-xc)+(y-yc-3.0*rc)*(y-yc-3.0*rc))+rc    !+1e-08
              bvess_lev(1,ipoin,1)=fleve(ipoin,3)
           endif
        end do

     else if(kfl_inlev_lev==31) then
        !
        ! Equilibrium interface SM testing
        !
        xmin_jump = 0.35
        xmax_jump = 0.5
        h1        = 0.032
        h2        = 0.149953
        mramp = (h2 - h1)/(xmax_jump - xmin_jump);

        do ipoin=1,npoin
           y=coord(2,ipoin)
           x=coord(1,ipoin)

           if( x <= xmin_jump ) then !0.2950 is the x-coordinate where the jump starts

              fleve(ipoin,3)     = -(y-h1)+1e-08
              bvess_lev(1,ipoin,1) = -(y-h1)+1e-08

           else if( x > xmin_jump .and. x < xmax_jump ) then

              ylevel = min(mramp*(x - xmin_jump) + h1, h2)
              fleve(ipoin,3)     = -(y-ylevel)+1e-08
              bvess_lev(1,ipoin,1) = -(y-ylevel)+1e-08

           else ! if( x >= xmax_jump) then

              ylevel_max = h2;

              fleve(ipoin,3)     = -(y-ylevel_max)+1e-08
              bvess_lev(1,ipoin,1) = -(y-ylevel_max)+1e-08

           end if

        end do

        if (.not.associated(veloc)) call runend('lev_iniunk: veloc not associated')

        vy=0.0_rp
        do ipoin=1,npoin
           y=coord(2,ipoin)
           x=coord(1,ipoin)

           if( x <= xmin_jump ) then !0.2950 is the x-coordinate where the jump starts

              yc = h1
              if( y <= yc) then !0.2950 is the x-coordinate where the jump starts
                 vx               = 1.5_rp
                 veloc(1,ipoin,1) = vx
                 veloc(1,ipoin,2) = vx
                 veloc(1,ipoin,3) = vx
              else !air
                 vx               = 0.0_rp
                 veloc(1,ipoin,1) = vx
                 veloc(1,ipoin,2) = vx
                 veloc(1,ipoin,3) = vx
              end if

           else if( x > xmin_jump .and. x < xmax_jump ) then

              ylevel = min(mramp*(x - xmin_jump) + h1, h2)
              if( y <= ylevel) then !water under the ramp
                 vx               = 0.0 !1.5_rp*h1/ylevel  !continuity
                 veloc(1,ipoin,1) = vx
                 veloc(1,ipoin,2) = vx
                 veloc(1,ipoin,3) = vx

              else !air
                 vx               = 0.0
                 veloc(1,ipoin,1) = vx
                 veloc(1,ipoin,2) = vx
                 veloc(1,ipoin,3) = vx
              end if

           else
              ylevel_max = mramp*xmax_jump;

              yc = h2
              if( y <= yc) then
                 vx               = 0.643000 !U_out
                 veloc(1,ipoin,1) = vx
                 veloc(1,ipoin,2) = vx
                 veloc(1,ipoin,3) = vx
              else !air
                 vx               = 0.0
                 veloc(1,ipoin,1) = vx
                 veloc(1,ipoin,2) = vx
                 veloc(1,ipoin,3) = vx

              end if

           end if

           veloc(2,ipoin,1)=vy
           veloc(2,ipoin,2)=vy
           veloc(2,ipoin,3)=vy
        enddo

     else if(kfl_inlev_lev==32) then
        !
        ! Equilibrium interface SM testing
        !
        xmin_jump = 0.35
        xmax_jump = 1.0
        h1        = 0.07
        h2        = 0.35
        mramp = (h2 - h1)/(xmax_jump - xmin_jump);

        do ipoin=1,npoin
           y=coord(2,ipoin)
           x=coord(1,ipoin)

           if( x <= xmin_jump ) then !0.2950 is the x-coordinate where the jump starts
              fleve(ipoin,3)     = -(y-h1)+1e-08
              bvess_lev(1,ipoin,1) = -(y-h1)+1e-08

           else if( x > xmin_jump .and. x < xmax_jump ) then

              ylevel = min(mramp*(x - xmin_jump) + h1, h2)
              fleve(ipoin,3)     = -(y-ylevel)+1e-08
              bvess_lev(1,ipoin,1) = -(y-ylevel)+1e-08

           else ! if( x >= xmax_jump) then

              ylevel_max = h2;

              fleve(ipoin,3)     = -(y-ylevel_max)+1e-08
              bvess_lev(1,ipoin,1) = -(y-ylevel_max)+1e-08

           end if

        end do

        if (.not.associated(veloc)) call runend('lev_iniunk: veloc not associated')

        vy=0.0_rp
        do ipoin=1,npoin
           y=coord(2,ipoin)
           x=coord(1,ipoin)

           if( x <= xmin_jump ) then !0.2950 is the x-coordinate where the jump starts

              yc = h1
              if( y <= yc) then !0.2950 is the x-coordinate where the jump starts
                 vx               = 1.5_rp
                 veloc(1,ipoin,1) = vx
                 veloc(1,ipoin,2) = vx
                 veloc(1,ipoin,3) = vx
              else !air
                 vx               = 0.0_rp
                 veloc(1,ipoin,1) = vx
                 veloc(1,ipoin,2) = vx
                 veloc(1,ipoin,3) = vx
              end if

           else if( x > xmin_jump .and. x < xmax_jump ) then

              ylevel = min(mramp*(x - xmin_jump) + h1, h2)
              if( y <= ylevel) then !water under the ramp
                 vx               = 0.0 !1.5_rp*h1/ylevel  !continuity
                 veloc(1,ipoin,1) = vx
                 veloc(1,ipoin,2) = vx
                 veloc(1,ipoin,3) = vx

              else !air
                 vx               = 0.0
                 veloc(1,ipoin,1) = vx
                 veloc(1,ipoin,2) = vx
                 veloc(1,ipoin,3) = vx
              end if

           else
              ylevel_max = mramp*xmax_jump;

              yc = h2
              if( y <= yc) then
                 vx               = 0.643000 !U_out
                 veloc(1,ipoin,1) = vx
                 veloc(1,ipoin,2) = vx
                 veloc(1,ipoin,3) = vx
              else !air
                 vx               = 0.0
                 veloc(1,ipoin,1) = vx
                 veloc(1,ipoin,2) = vx
                 veloc(1,ipoin,3) = vx

              end if

           end if

           veloc(2,ipoin,1)=vy
           veloc(2,ipoin,2)=vy
           veloc(2,ipoin,3)=vy
        enddo

     else if(kfl_inlev_lev==33) then
        !
        ! Equilibrium interface SM testing
        !
        xmin_jump = 1.0_rp
        xmax_jump = 1.8_rp
        h1        = 0.032_rp
        h2        = 0.15_rp
        mramp = (h2 - h1)/(xmax_jump - xmin_jump);

        do ipoin=1,npoin
           y=coord(2,ipoin)
           x=coord(1,ipoin)

           if( x <= xmin_jump ) then !0.2950 is the x-coordinate where the jump starts
              fleve(ipoin,3)     = -(y-h1)+1e-08
              bvess_lev(1,ipoin,1) = -(y-h1)+1e-08

           else if( x > xmin_jump .and. x < xmax_jump ) then

              ylevel = min(mramp*(x - xmin_jump) + h1, h2)
              fleve(ipoin,3)     = -(y-ylevel)+1e-08
              bvess_lev(1,ipoin,1) = -(y-ylevel)+1e-08

           else ! if( x >= xmax_jump) then

              ylevel_max = h2;

              fleve(ipoin,3)     = -(y-ylevel_max)+1e-08
              bvess_lev(1,ipoin,1) = -(y-ylevel_max)+1e-08

           end if

        end do

     else if(kfl_inlev_lev==34) then
        !
        ! Equilibrium interface SM testing
        !
        xmin_jump = 0.3_rp
        xmax_jump = 0.5_rp
        h1        = 0.05_rp
        h2        = 0.138_rp

        do ipoin=1,npoin
           y=coord(2,ipoin)
           x=coord(1,ipoin)

           if( x <= xmin_jump ) then !0.2950 is the x-coordinate where the jump starts

              fleve(ipoin,3)     = -(y-h1)+1e-08
              bvess_lev(1,ipoin,1) = -(y-h1)+1e-08

           else if( x > xmin_jump ) then

              fleve(ipoin,3)     = -(y-h2)+1e-08
              bvess_lev(1,ipoin,1) = -(y-h2)+1e-08

           end if

        end do

     else if( kfl_inlev_lev == 35 ) then
        !
        ! Cylinder
        !
        xc=0.0_rp
        yc=0.0_rp
        rc=0.5_rp

        do ipoin = 1,npoin

           y = coord(2,ipoin)
           x = coord(1,ipoin)
           d = sqrt((x-xc)*(x-xc)+(y-yc)*(y-yc))

           if( d <=rc ) then
              !
              ! Air
              !
              fleve(ipoin,3)       = -d
              bvess_lev(1,ipoin,1) = -d

           else 
              !
              ! Water
              !
              fleve(ipoin,3)       =  d
              bvess_lev(1,ipoin,1) =  d

           end if

        end do

     end if

  end if

end subroutine lev_iniun0

