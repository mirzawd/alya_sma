!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine pts_analytics(ielem, dt, coord_p, coord_kp1)
  
    use def_domain
    use def_kermod
    use def_master
    use mod_elmgeo
    use def_kintyp, only : ip,rp
    use mod_maths,  only : maths_invert_matrix
    implicit none
    !real(rp),    intent(in)  :: coord(*)
    !real(rp),    intent(in)  :: lnods(*)
    real(rp),    intent(in)  :: dt
    !integer(ip), intent(in)  :: ndime
    integer(ip), intent(in)  :: ielem 
    real(rp)                 :: detm
    real(rp)                 :: node1(ndime), node2(ndime), node3(ndime), node4(ndime)
    real(rp)                 :: veloc1(ndime), veloc2(ndime), veloc3(ndime), veloc4(ndime)
    real(rp),   intent(in)   :: coord_p(ndime)
    real(rp),   intent(out)  :: coord_kp1(ndime)
    real(rp)                 :: a_s, b_s, c_s, a_t, b_t, c_t, A_x, A_y, B_x, B_y, C_x, C_y
    real(rp)                 :: landa_p, landa_m, sq

    real(rp) :: denom
    real(rp) :: aa(ndime,ndime),bb(ndime)
    real(rp) :: tt(ndime,ndime),tinv(ndime,ndime)
    real(rp) :: a1,b1,c1
    real(rp) :: a2,b2,c2
    real(rp) :: a3,b3,c3
    real(rp) :: a,b,c,d
!    real(rp) :: t
    real(rp) :: l1,l2,xx(3),x0(3)

    !print *, 'x1,y1', coord(1:ndime, lnods(1, ielem))
    !print *, 'x2,y2', coord(1:ndime, lnods(2, ielem))
    !print *, 'x3,y3', coord(1:ndime, lnods(3, ielem))
    
    node1(1:ndime) = coord(1:ndime, lnods(1, ielem))
    node2(1:ndime) = coord(1:ndime, lnods(2, ielem)) 
    node3(1:ndime) = coord(1:ndime, lnods(3, ielem))
    if (ndime == 3) then
      node4(1:ndime) = coord(1:ndime, lnods(4, ielem))
    end if

    veloc1(1:ndime) = advec(1:ndime, lnods(1, ielem),1)
    veloc2(1:ndime) = advec(1:ndime, lnods(2, ielem),1)
    veloc3(1:ndime) = advec(1:ndime, lnods(3, ielem),1)
    if (ndime == 3) then
       veloc4(1:ndime) = advec(1:ndime, lnods(4, ielem),1)
    end if
    
    !veloc1(1)=1.0_rp;veloc1(2)=2.0_rp
    !veloc2(1)=1.0_rp;veloc2(2)=2.0_rp
    !veloc3(1)=1.0_rp;veloc3(2)=2.0_rp
    !x0(1:ndime) = (node1(1:ndime)+node2(1:ndime)+node3(1:ndime))/3.0_rp

    !veloc1=1.0_rp;veloc2=1.0_rp;veloc3=1.0_rp
    
    denom   = (node2(2)-node3(2))*(node1(1)-node3(1)) + (node3(1)-node2(1))*(node1(2)-node3(2)) 
    a1      = (node2(2)-node3(2)) / denom
    b1      = (node3(1)-node2(1)) / denom
    c1      = ( (node2(2)-node3(2))*(-node3(1)) + (node3(1)-node2(1))*(-node3(2)) ) / denom
    
    a2      = (node3(2)-node1(2)) / denom
    b2      = (node1(1)-node3(1)) / denom
    c2      = ( (node3(2)-node1(2))*(-node3(1)) + (node1(1)-node3(1))*(-node3(2)) ) / denom
    
    a3      = -a1-a2
    b3      = -b1-b2
    c3      = 1.0_rp-c1-c2


    aa(1,1) = a1*veloc1(1)+a2*veloc2(1)+a3*veloc3(1)
    aa(1,2) = b1*veloc1(1)+b2*veloc2(1)+b3*veloc3(1)
    
    aa(2,1) = a1*veloc1(2)+a2*veloc2(2)+a3*veloc3(2)
    aa(2,2) = b1*veloc1(2)+b2*veloc2(2)+b3*veloc3(2)

    bb(1)   = c1*veloc1(1)+c2*veloc2(1)+c3*veloc3(1)
    bb(2)   = c1*veloc1(2)+c2*veloc2(2)+c3*veloc3(2)

    a       = 1.0_rp
    b       = -(aa(1,1)+aa(2,2))
    c       = aa(1,1)*aa(2,2)-aa(1,2)*aa(2,1)
    d       = b**2 - 4.0_rp*a*c

    !t       = aa(1,1)+aa(2,2)
    !d       = aa(1,1)*aa(2,2)-aa(1,2)*aa(2,1)
    !l1      = 0.5_rp * t + sqrt(0.25_rp*t**2-d)
    !l1      = 0.5_rp * t - sqrt(0.25_rp*t**2-d)
    !xx(1)   = aa(1,1)*x0(1)+aa(1,2)*x0(2)+bb(1)
    !xx(2)   = aa(2,1)*x0(1)+aa(2,2)*x0(2)+bb(2)
    !print*,'x=',xx
    !stop
    if(d<0.0_rp) then
       d = -1.0_rp*d
       print*,'a=',a,b,c
       stop
    end if
       print*,'d=',a,b,c,d

    l1      = (-b-sqrt(d))/(2.0_rp*a)
    l2      = (-b+sqrt(d))/(2.0_rp*a)

    if( c /= 0.0_rp ) then
       tt(1,1) = l1-d
       tt(2,1) = c
       tt(1,2) = l2-d
       tt(2,2) = c
    else if( b /= 0.0_rp ) then
       tt(1,1) = b
       tt(2,1) = l1-a
       tt(1,2) = b
       tt(2,2) = l2-a
    else
       tt(1,1) = 1.0_rp
       tt(2,1) = 0.0_rp
       tt(1,2) = 0.0_rp
       tt(2,2) = 1.0_rp     
    end if
    tinv    = tt
    call maths_invert_matrix(ndime,tinv)
    !
    ! alpha = T^{-1} b
    ! x0'   = T^{-1} x0
    !
    a1           = tinv(1,1)*bb(1)+tinv(1,2)*bb(2)
    a2           = tinv(1,2)*bb(1)+tinv(2,2)*bb(2)
    x0(1)        = tinv(1,1)*coord_p(1)+tinv(1,2)*coord_p(2)
    x0(2)        = tinv(1,2)*coord_p(1)+tinv(2,2)*coord_p(2)
    print*,'A=',l1,l2

    xx(1)        = (x0(1)+a1/l1)*exp(-l1*dt)-a1/l1
    xx(2)        = (x0(2)+a2/l2)*exp(-l2*dt)-a2/l2
    print*,'B'
    !
    ! x = T x'
    !
    coord_kp1(1) = tt(1,1)*xx(1)+tt(1,2)*xx(2)
    coord_kp1(2) = tt(2,1)*xx(1)+tt(2,2)*xx(2)
    print*,'C'

    return
    !print *, 'veloc1', veloc1(1:ndime)
    !print *, 'veloc2', veloc2(1:ndime)
    !print *, 'veloc3', veloc3(1:ndime)
    ! Calculate the main Matrix determinant
    if (ndime == 2) then
      detm = (node2(1) - node1(1)) * (node3(2) - node1(2)) &
           - (node3(1) - node1(1)) * (node2(2) - node1(2))
      !print *, 'detm',  detm
  
    else if (ndime == 3) then
      detm = (node2(1) - node1(1)) * (node3(2) - node1(2)) * (node4(3) - node1(3)) &
           + (node3(1) - node1(1)) * (node4(2) - node1(2)) * (node2(3) - node1(3)) &
           + (node4(1) - node1(1)) * (node2(2) - node1(2)) * (node3(3) - node1(3)) &
           - (node4(1) - node1(1)) * (node3(2) - node1(2)) * (node2(3) - node1(3)) &
           - (node3(1) - node1(1)) * (node2(2) - node1(2)) * (node4(3) - node1(3)) &
           - (node2(1) - node1(1)) * (node4(2) - node1(2)) * (node3(3) - node1(3))
    end if
    ! s (= N_2) = a_s*x + b_s*y + c_s
    ! t (= N_3) = a_t*x + b_t*y + c_t
    a_s = (node3(2) - node1(2))/detm
    b_s = - (node3(1) - node1(1))/detm
    c_s = (node1(2) * (node3(1) - node1(1)) - node1(1) * (node3(2) - node1(2)))/detm
    a_t = - (node2(2) - node1(2))/detm
    b_t = (node2(1) - node1(1))/detm
    c_t = (node1(1) * (node1(2) - node1(2)) - node1(2) * (node2(1) - node1(1)))/detm
    !print *, 'a_s, b_s, c_s, a_t, b_t, c_t', a_s, b_s, c_s, a_t, b_t, c_t
    A_x = - (a_s + a_t) * veloc1(1) + a_s * veloc2(1) + a_t * veloc3(1)
    A_y = - (a_s + a_t) * veloc1(2) + a_s * veloc2(2) + a_t * veloc3(2)
    B_x = - (b_s + b_t) * veloc1(1) + b_s * veloc2(1) + b_t * veloc3(1)
    B_y = - (b_s + b_t) * veloc1(2) + b_s * veloc2(2) + b_t * veloc3(2)
    C_x = (1.0_rp - c_s - c_t) * veloc1(1) + c_s * veloc2(1) + c_t * veloc3(1)
    C_y = (1.0_rp - c_s - c_t) * veloc1(2) + c_s * veloc2(2) + c_t * veloc3(2)
    !print *, 'A_x, A_y, B_x, B_y, C_x, C_y', A_x, A_y, B_x, B_y, C_x, C_y
    ! Eigenvalues
    sq = 4.0_rp * A_y * B_x + (A_x - B_y)**2
    if (sq < 0.0_rp ) then
        print *, 'raiz negativa', sq
        sq = -1.0_rp*sq
        print *, 'nueva raiz', sq
    end if
    landa_p = 0.5_rp * (A_x + B_y + sqrt(sq))
    landa_m = 0.5_rp * (A_x + B_y - sqrt(sq))
    coord_kp1(1) = exp(dt * landa_m) * coord_p(1) + (C_x / landa_m) * (exp(dt * landa_m) - 1)
    coord_kp1(2) = exp(dt * landa_p) * coord_p(2) + (C_y / landa_p) * (exp(dt * landa_p) - 1)

end subroutine pts_analytics
