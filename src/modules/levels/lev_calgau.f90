!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine lev_calgau(trian)
  !-----------------------------------------------------------------------
  !****f* Levels/lev_calgau
  ! NAME 
  !    lev_calgau
  ! DESCRIPTION
  !    Compute Interface Gauges
  ! USES
  !    
  ! USED BY
  !    lev_outgau
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_levels
  use def_domain
  implicit none

  real(rp),    intent(in)  :: trian(ndime,ndime)
  integer (ip)             :: igaug, sames
  real(rp)                 :: x1,y1,z1,x2,y2,z2,x3,y3,z3,xp,yp,zp, s, t
  real(rp)                 :: pscal, cp1x, cp1y, cp1z, cp2x, cp2y, cp2z
  real(rp)                 :: deter,x2x1,y2y1,x3x1,y3y1,xx1,yy1

  do igaug = 1, npp_nbgau_lev

     if((typga_lev(igaug)== 3_ip).and.(findg_lev(igaug)/=1_ip)) then

        x1 = trian(1,1)
        y1 = trian(2,1)
        z1 = 0.0_rp

        x2 = trian(1,2)
        y2 = trian(2,2)
        z2 = 0.0_rp

        x3 = trian(1,3)
        y3 = trian(2,3)
        z3 = 0.0_rp

        xp = cogau_lev(1,igaug)
        yp = cogau_lev(2,igaug)
        zp = 0.0_rp

        sames=0_ip

        if((x1/=x2).or.(x1/=x3)) then
           if((y1/=y2).or.(y1/=y3)) then
              cp1x=(y2-y1)*(zp-z1)-(z2-z1)*(yp-y1)
              cp1y=(z2-z1)*(xp-x1)-(x2-x1)*(zp-z1)
              cp1z=(x2-x1)*(yp-y1)-(y2-y1)*(xp-x1)

              cp2x=(y2-y1)*(z3-z1)-(z2-z1)*(y3-y1)
              cp2y=(z2-z1)*(x3-x1)-(x2-x1)*(z3-z1)
              cp2z=(x2-x1)*(y3-y1)-(y2-y1)*(x3-x1)

              pscal=cp1x*cp2x+cp1y*cp2y+cp1z*cp2z      

              if(pscal>=0) then
                 sames=sames+1
              endif

              cp1x=(y3-y2)*(zp-z2)-(z3-z2)*(yp-y2)
              cp1y=(z3-z2)*(xp-x2)-(x3-x2)*(zp-z2)
              cp1z=(x3-x2)*(yp-y2)-(y3-y2)*(xp-x2)

              cp2x=(y3-y2)*(z1-z2)-(z3-z2)*(y1-y2)
              cp2y=(z3-z2)*(x1-x2)-(x3-x2)*(z1-z2)
              cp2z=(x3-x2)*(y1-y2)-(y3-y2)*(x1-x2)

              pscal=cp1x*cp2x+cp1y*cp2y+cp1z*cp2z      

              if(pscal>=0) then
                 sames=sames+1
              endif

              cp1x=(y1-y3)*(zp-z3)-(z1-z3)*(yp-y3)
              cp1y=(z1-z3)*(xp-x3)-(x1-x3)*(zp-z3)
              cp1z=(x1-x3)*(yp-y3)-(y1-y3)*(xp-x3)

              cp2x=(y1-y3)*(z2-z3)-(z1-z3)*(y2-y3)
              cp2y=(z1-z3)*(x2-x3)-(x1-x3)*(z2-z3)
              cp2z=(x1-x3)*(y2-y3)-(y1-y3)*(x2-x3)

              pscal=cp1x*cp2x+cp1y*cp2y+cp1z*cp2z      

              if(pscal>=0) then
                 sames=sames+1
              endif
           endif
        endif

        !    For P1 triangles we have:
        !    x = (1-s-t)*x1 + s*x2 + t*x3
        !    y = (1-s-t)*y1 + s*y2 + t*y3
        !
        !    This linear problem is solved for (s,t):
        !         (x3-x1)(y-y1) -(y3-y1)(x-x1)
        !    s =  ----------------------------
        !         (x3-x1)(y2-y1)-(y3-y1)(x2-x1)
        !
        !         (x-x1)(y2-y1) -(y-y1)(x2-x1)
        !    t =  ----------------------------
        !         (x3-x1)(y2-y1)-(y3-y1)(x2-x1)

        if(sames==3_ip) then

           findg_lev(igaug) = 1_ip
           ! find  gauss coordinates in parametric triangle
           x2x1     = x2-x1
           y2y1     = y2-y1
           x3x1     = x3-x1
           y3y1     = y3-y1

           deter = x3x1*y2y1-y3y1*x2x1
           if( deter /= 0.0_rp ) deter = 1.0_rp/deter
           xx1   = cogau_lev(1,igaug) - x1
           yy1   = cogau_lev(2,igaug) - y1
           s     = deter*(x3x1*yy1-y3y1*xx1)
           t     = deter*(y2y1*xx1-x2x1*yy1)
           valga_lev(igaug) = (1.0_rp -s -t) * trian(3,1)+ s* trian(3,2) + t* trian(3,3)

        endif

     else if (typga_lev(igaug)/= 3_ip) then

        print*,' only vertical type of gauge is implemented '
        stop  

     end if

  enddo

end subroutine lev_calgau
