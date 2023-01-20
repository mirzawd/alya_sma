!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



module mod_mag_quadra

  use def_kintyp, only: ip, rp

  implicit none

  type stcqua
    integer(ip) :: nq
    real(rp), allocatable :: wq(:)
    real(rp), allocatable :: xq(:,:)
    real(rp), allocatable :: lq(:,:)
  end type stcqua

contains

  !##############################################################
  function quadraTetrahedron(N) result(quad)

    implicit none

    integer(ip), intent(in) :: N

    type (stcqua) :: quad

    select case (N)

    case (1_ip)
      quad%nq = 4_ip
      allocate(quad%wq(quad%nq), quad%xq(3, quad%nq))
      quad%wq = [ 1.0_rp, 1.0_rp, 1.0_rp, 1.0_rp ] / 24.0_rp
      quad%xq = reshape( [ 0.5854101966249685_rp, 0.1381966011250105_rp, 0.1381966011250105_rp, &
                           0.1381966011250105_rp, 0.1381966011250105_rp, 0.1381966011250105_rp, &
                           0.1381966011250105_rp, 0.1381966011250105_rp, 0.5854101966249685_rp, &
                           0.1381966011250105_rp, 0.5854101966249685_rp, 0.1381966011250105_rp &
                          ], [ 3, 4 ] )

    case (2_ip)
      quad%nq = 5_ip
      allocate(quad%wq(quad%nq), quad%xq(3, quad%nq))
      quad%wq = [-0.8_rp, 0.45_rp, 0.45_rp, 0.45_rp, 0.45_rp] / 6.0_rp
      quad%xq = reshape( [ 6.0_rp, 6.0_rp, 6.0_rp, &
                           12.0_rp, 4.0_rp, 4.0_rp, &
                           4.0_rp, 4.0_rp, 4.0_rp, &
                           4.0_rp, 4.0_rp, 12.0_rp, &
                           4.0_rp, 12.0_rp, 4.0_rp &
                         ], [3, 5] ) / 24.0_rp

    case (3_ip)
      quad%nq = 11_ip
      allocate(quad%wq(quad%nq), quad%xq(3, quad%nq))
      quad%wq = [ -0.0789333333333333_rp, 0.0457333333333333_rp, 0.0457333333333333_rp, 0.0457333333333333_rp, &
                   0.0457333333333333_rp, 0.1493333333333333_rp, 0.1493333333333333_rp, 0.1493333333333333_rp, &
                   0.1493333333333333_rp, 0.1493333333333333_rp, 0.1493333333333333_rp] / 6.0_rp
      quad%xq = reshape( [ 0.25_rp, 0.25_rp, 0.25_rp, &
                           0.7857142857142857_rp, 0.0714285714285714_rp, 0.0714285714285714_rp, &
                           0.0714285714285714_rp, 0.0714285714285714_rp, 0.0714285714285714_rp, &
                           0.0714285714285714_rp, 0.0714285714285714_rp, 0.7857142857142857_rp, &
                           0.0714285714285714_rp, 0.7857142857142857_rp, 0.0714285714285714_rp, &
                           0.1005964238332008_rp, 0.3994035761667992_rp, 0.3994035761667992_rp, &
                           0.3994035761667992_rp, 0.1005964238332008_rp, 0.3994035761667992_rp, &
                           0.3994035761667992_rp, 0.3994035761667992_rp, 0.1005964238332008_rp, &
                           0.3994035761667992_rp, 0.1005964238332008_rp, 0.1005964238332008_rp, &
                           0.1005964238332008_rp, 0.3994035761667992_rp, 0.1005964238332008_rp, &
                           0.1005964238332008_rp, 0.1005964238332008_rp, 0.3994035761667992_rp &
                         ], [3, 11] )

    case (4_ip)
      quad%nq = 15_ip
      allocate(quad%wq(quad%nq), quad%xq(3, quad%nq))
      quad%wq = [ 0.1817020685825351_rp, 0.0361607142857143_rp, 0.0361607142857143_rp, &
                  0.0361607142857143_rp, 0.0361607142857143_rp, 0.0698714945161738_rp, &
                  0.0698714945161738_rp, 0.0698714945161738_rp, 0.0698714945161738_rp, &
                  0.0656948493683187_rp, 0.0656948493683187_rp, 0.0656948493683187_rp, &
                  0.0656948493683187_rp, 0.0656948493683187_rp, 0.0656948493683187_rp ] / 6.0_rp
      quad%xq = reshape( [ 0.25_rp, 0.25_rp, 0.25_rp, &
                           0.0_rp, 1.0_rp / 3.0_rp, 1.0_rp / 3.0_rp, &
                           1.0_rp / 3.0_rp, 1.0_rp / 3.0_rp, 1.0_rp / 3.0_rp, &
                           1.0_rp / 3.0_rp, 1.0_rp / 3.0_rp, 0.0_rp, &
                           1.0_rp / 3.0_rp, 0.0_rp, 1.0_rp / 3.0_rp, &
                           0.7272727272727273_rp, 0.0909090909090909_rp, 0.0909090909090909_rp, &
                           0.0909090909090909_rp, 0.0909090909090909_rp, 0.0909090909090909_rp, &
                           0.0909090909090909_rp, 0.0909090909090909_rp, 0.7272727272727273_rp, &
                           0.0909090909090909_rp, 0.7272727272727273_rp, 0.0909090909090909_rp, &
                           0.4334498464263357_rp, 0.0665501535736643_rp, 0.0665501535736643_rp, &
                           0.0665501535736643_rp, 0.4334498464263357_rp, 0.0665501535736643_rp, &
                           0.0665501535736643_rp, 0.0665501535736643_rp, 0.4334498464263357_rp, &
                           0.0665501535736643_rp, 0.4334498464263357_rp, 0.4334498464263357_rp, &
                           0.4334498464263357_rp, 0.0665501535736643_rp, 0.4334498464263357_rp, &
                           0.4334498464263357_rp, 0.4334498464263357_rp, 0.0665501535736643_rp &
                         ], [3, 15] )

    case default
      call runend("quadraTetrahedron: error")

    end select

  end function quadraTetrahedron
  !##############################################################


  !##############################################################
  function quadraHexahedron(N) result(quad)

    implicit none

    integer(ip), intent(in) :: N
    
    type (stcqua) :: quad, quadLin

    integer(ip) :: i, j, k, s

    quadLin = quadraLine(N)

    quad % nq = N**3
    allocate(quad % wq(quad % nq), quad % xq(3, quad % nq))

    s = 0_ip
    do i = 1, N
      do j = 1, N
        do k = 1, N
          s = s + 1_ip
          quad % wq(s) = quadLin % wq(i) * quadLin % wq(j) * quadLin % wq(k)
          quad % xq(1:3, s) = [ quadLin % xq(1, i), quadLin % xq(1, j), quadLin % xq(1, k) ]
        end do
      end do
    end do

  end function quadraHexahedron
  !##############################################################


  !##############################################################
  function quadraLine(N) result(quad)

    implicit none

    integer(ip), intent(in) ::             N

    type (stcqua) ::                    quad

    ! Gauss quadrature points in the 1-D domain [-1, 1]

    quad%nq = N
    allocate(quad%wq(quad%nq), quad%xq(1, quad%nq))

    select case (N)

    case (1)
      quad%wq(1) = 2.0_rp

      quad%xq(1, 1) = 0.0_rp

    case (2)
      quad%wq(1) = 1.0_rp
      quad%wq(2) = 1.0_rp

      quad%xq(1, 1) = -1.0_rp / sqrt(3.0_rp)
      quad%xq(1, 2) = 1.0_rp / sqrt(3.0_rp)

    case (3)
      quad%wq(1) = 5.0_rp / 9.0_rp
      quad%wq(2) = 8.0_rp / 9.0_rp
      quad%wq(3) = 5.0_rp / 9.0_rp

      quad%xq(1, 1) = -sqrt(0.6_rp)
      quad%xq(1, 2) = 0.0_rp
      quad%xq(1, 3) = sqrt(0.6_rp)

    case (4)
      quad%wq(1) = 0.5_rp - sqrt(30.0_rp) / 36.0_rp
      quad%wq(2) = 0.5_rp + sqrt(30.0_rp) / 36.0_rp
      quad%wq(3) = 0.5_rp + sqrt(30.0_rp) / 36.0_rp
      quad%wq(4) = 0.5_rp - sqrt(30.0_rp) / 36.0_rp

      quad%xq(1, 1) = -sqrt((3.0_rp + 2.0_rp * sqrt(1.2_rp)) / 7.0_rp)
      quad%xq(1, 2) = -sqrt((3.0_rp - 2.0_rp * sqrt(1.2_rp)) / 7.0_rp)
      quad%xq(1, 3) = sqrt((3.0_rp - 2.0_rp * sqrt(1.2_rp)) / 7.0_rp)
      quad%xq(1, 4) = sqrt((3.0_rp + 2.0_rp * sqrt(1.2_rp)) / 7.0_rp)

    case (5)
      quad%wq(1) = (322.0_rp - 13.0_rp * sqrt(70.0_rp)) / 900.0_rp
      quad%wq(2) = (322.0_rp + 13.0_rp * sqrt(70.0_rp)) / 900.0_rp
      quad%wq(3) = 128.0_rp / 225.0_rp
      quad%wq(4) = (322.0_rp + 13.0_rp * sqrt(70.0_rp)) / 900.0_rp
      quad%wq(5) = (322.0_rp - 13.0_rp * sqrt(70.0_rp)) / 900.0_rp

      quad%xq(1, 1) = -sqrt(5.0_rp + 2.0_rp * sqrt(10.0_rp / 7.0_rp)) / 3.0_rp
      quad%xq(1, 2) = -sqrt(5.0_rp - 2.0_rp * sqrt(10.0_rp / 7.0_rp)) / 3.0_rp
      quad%xq(1, 3) = 0.0_rp
      quad%xq(1, 4) = sqrt(5.0_rp - 2.0_rp * sqrt(10.0_rp / 7.0_rp)) / 3.0_rp
      quad%xq(1, 5) = sqrt(5.0_rp + 2.0_rp * sqrt(10.0_rp / 7.0_rp)) / 3.0_rp

    case (6)
      quad%wq(1) = 0.171324492379170_rp
      quad%wq(2) = 0.360761573048139_rp
      quad%wq(3) = 0.467913934572691_rp
      quad%wq(4) = 0.467913934572691_rp
      quad%wq(5) = 0.360761573048139_rp
      quad%wq(6) = 0.171324492379170_rp

      quad%xq(1, 1) = -0.932469514203152_rp
      quad%xq(1, 2) = -0.661209386466265_rp
      quad%xq(1, 3) = -0.238619186083197_rp
      quad%xq(1, 4) = 0.238619186083197_rp
      quad%xq(1, 5) = 0.661209386466265_rp
      quad%xq(1, 6) = 0.932469514203152_rp

    case (7)
      quad%wq(1) = 0.129484966168870_rp
      quad%wq(2) = 0.279705391489277_rp
      quad%wq(3) = 0.381830050505119_rp
      quad%wq(4) = 0.417959183673469_rp
      quad%wq(5) = 0.381830050505119_rp
      quad%wq(6) = 0.279705391489277_rp
      quad%wq(7) = 0.129484966168870_rp

      quad%xq(1, 1) = -0.949107912342759_rp
      quad%xq(1, 2) = -0.741531185599394_rp
      quad%xq(1, 3) = -0.405845151377397_rp
      quad%xq(1, 4) = 0.0_rp
      quad%xq(1, 5) = 0.405845151377397_rp
      quad%xq(1, 6) = 0.741531185599394_rp
      quad%xq(1, 7) = 0.949107912342759_rp
  
    case default
      deallocate(quad%wq, quad%xq)
      call runend("quadraLine: error")

    end select

  end function quadraLine
  !##############################################################


  !##############################################################
  function quadraTriangle(N) result(quad)

    implicit none

    integer(ip), intent(in) ::              N

    type (stcqua) ::                      quad

    select case (N)

    case (1)
      ! 1 Gauss Node
      quad%nq = 1_ip
      allocate(quad%wq(quad%nq), quad%xq(2, quad%nq))
      quad%wq = [ 0.5_rp ]
      quad%xq = reshape( [ 1.0_rp, 1.0_rp ], [ 2, 1 ] ) / 3.0_rp

    case (2)
      ! 3 Gauss Nodes
      quad%nq = 3_ip
      allocate(quad%wq(quad%nq), quad%xq(2, quad%nq))
      quad%wq = [ 1.0_rp, 1.0_rp, 1.0_rp ] / 6.0_rp
      quad%xq = reshape( [ 0.0_rp, 1.0_rp, &
                                1.0_rp, 0.0_rp, &
                                1.0_rp, 1.0_rp ], &
                                [ 2, 3 ] ) / 2.0_rp

    case (3)
      ! 4 Gauss Nodes
      quad%nq = 4_ip
      allocate(quad%wq(quad%nq), quad%xq(2, quad%nq))
      quad%wq = [ -27.0_rp, 25.0_rp, 25.0_rp, 25.0_rp ] / (96.0_rp)
      quad%xq = reshape( [ 1.0_rp/3.0_rp, 1.0_rp/3.0_rp, &
                                1.0_rp/5.0_rp, 3.0_rp/5.0_rp, &
                                1.0_rp/5.0_rp, 1.0_rp/5.0_rp, &
                                3.0_rp/5.0_rp, 1.0_rp/5.0_rp ], &
                                [ 2, 4 ] )

    case (4)
      ! 7 Gauss Nodes
      quad%nq = 7_ip
      allocate(quad%wq(quad%nq), quad%xq(2, quad%nq))
      quad%wq = [ 1.0_rp / 40.0_rp, 1.0_rp / 15.0_rp, 1.0_rp / 40.0_rp, 1.0_rp / 15.0_rp, &
        1.0_rp / 40.0_rp, 1.0_rp / 15.0_rp, 9.0_rp / 40.0_rp ]
      quad%xq = reshape( [ 0.0_rp, 1.0_rp, &
                                0.0_rp, 0.5_rp, &
                                0.0_rp, 0.0_rp, &
                                0.5_rp, 0.0_rp, &
                                1.0_rp, 0.0_rp, &
                                0.5_rp, 0.5_rp, &
                                1.0_rp / 3.0_rp, 1.0_rp / 3.0_rp ], [ 2, 7 ] )

    case default
      call runend("quadraTriangle: error")

    end select

  end function quadraTriangle
  !##############################################################


  !##############################################################
  function quadraQuadrangle(N) result(quad)

    implicit none

    integer(ip), intent(in) ::              N

    type (stcqua) ::                      quad, quadLin

    integer(ip) ::     i, j, s

    quadLin = quadraLine(N)

    quad%nq = N**2
    allocate(quad%wq(quad%nq), quad%xq(2, quad%nq))

    s = 0_ip
    do i = 1, N
      do j = 1, N
        s = s + 1_ip
        quad % wq(s) = quadLin % wq(i) * quadLin % wq(j)
        quad % xq(1:2, s) = [ quadLin % xq(1, i), quadLin % xq(1, j) ]
      end do
    end do

  end function quadraQuadrangle
  !##############################################################

end module mod_mag_quadra
