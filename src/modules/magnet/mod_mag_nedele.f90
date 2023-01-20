!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



module mod_mag_nedele

  use def_kintyp, only: ip, rp
  use mod_mag_linalg, only: matdet, mag_vecmod
  use mod_mag_lagran, only: mag_lagdif

  implicit none

contains

  !##############################################################
  function mag_edgdof(pelty) result(numdof)

    implicit none

    integer(ip), intent(in) :: pelty

    integer(ip) :: numdof

    select case (pelty)

    case (1_ip)
      numdof = 0_ip

    case (10_ip)
      ! TRI03
      numdof = 3_ip

    case (12_ip)
      ! QUA04
      numdof = 4_ip

    case (30_ip)
      ! TET04
      numdof = 6_ip

    case (37_ip)
      ! HEX08
      numdof = 12_ip

    case default
      numdof = 0_ip
      call runend("mag_edgdof: unknown element type")

    end select

  end function mag_edgdof
  !##############################################################


  !##############################################################
  function mag_nedtri(n0, nodes) result(basis)

    implicit none

    real(rp), intent(in) :: n0(:)
    real(rp), intent(in) :: nodes(:,:)

    real(rp) :: basis(2, 3)
    real(rp), dimension(2) :: nA, nB, nC
    real(rp) :: x, y, m
    real(rp), dimension(3) :: l, a, b, c, d

    nA = nodes(1:2, 1)
    nB = nodes(1:2, 2)
    nC = nodes(1:2, 3)

    l(1) = sqrt( (nB(1) - nC(1))**2 + (nB(2) - nC(2))**2 )
    l(2) = sqrt( (nA(1) - nC(1))**2 + (nA(2) - nC(2))**2 )
    l(3) = sqrt( (nB(1) - nA(1))**2 + (nB(2) - nA(2))**2 )

    d(1) = (nA(1) - nC(1)) * (nB(2) - nC(2)) - (nB(1) - nC(1)) * (nA(2) - nC(2))
    a(1) = (nB(2) - nC(2)) / d(1)
    b(1) = (nC(1) - nB(1)) / d(1)
    c(1) = (nB(1) * nC(2) - nC(1) * nB(2)) / d(1)

    d(2) = (nB(1) - nA(1)) * (nC(2) - nA(2)) - (nC(1) - nA(1)) * (nB(2) - nA(2))
    a(2) = (nC(2) - nA(2)) / d(2)
    b(2) = (nA(1) - nC(1)) / d(2)
    c(2) = (nC(1) * nA(2) - nA(1) * nC(2)) / d(2)

    d(3) = -((nA(1) - nB(1)) * (nC(2) - nB(2)) - (nC(1) - nB(1)) * (nA(2) - nB(2)))
    a(3) = (nA(2) - nB(2)) / d(3)
    b(3) = (nB(1) - nA(1)) / d(3)
    c(3) = (nA(1) * nB(2) - nB(1) * nA(2)) / d(3)

    x = n0(1); y = n0(2)

    m = a(3) * b(2) - a(2) * b(3)
    basis(:, 1) = [ m * y + a(3) * c(2) - a(2) * c(3), -m * x + b(3) * c(2) - c(3) * b(2) ] * l(1)
    m = a(1) * b(3) - a(3) * b(1)
    basis(:, 2) = [ m * y + a(1) * c(3) - a(3) * c(1), -m * x + b(1) * c(3) - b(3) * c(1) ] * l(2)
    m = a(2) * b(1) - a(1) * b(2)
    basis(:, 3) = [ m * y + a(2) * c(1) - a(1) * c(2), -m * x + b(2) * c(1) - b(1) * c(2) ] * l(3)

  end function mag_nedtri
  !##############################################################


  !##############################################################
  function mag_nedrce(nodes) result(basis)

    implicit none

    real(rp), intent(in) :: &
      nodes(:, :)

    real(rp) :: &
      basis(2, 4),    &
      edgvec(2)

    integer(ip) :: &
      iedge,    &
      ipoin,    &
      jpoin

    do iedge = 1, 4
      ipoin = iedge
      jpoin = mod(ipoin, 4_ip) + 1_ip

      edgvec = nodes(1:2, jpoin) - nodes(1:2, ipoin)

      basis(:, iedge) = 0.5_rp * edgvec / sqrt( dot_product(edgvec, edgvec) )
    end do

  end function mag_nedrce
  !##############################################################


  !##############################################################
  function mag_nedrec(n0, nodes) result(basis)

    implicit none

    real(rp), intent(in) :: &
      n0(:),    &
      nodes(:, :)

    real(rp) :: &
      basis(2, 4),    &
      edgvec1(2),    &
      edgvec2(2),     &
      auxvec(2),      &
      len1,     &
      len2,   &
      weight1,        &
      weight2

    integer(ip) :: &
      iedge,    &
      ipoin1,    &
      jpoin1,   &
      ipoin2, &
      jpoin2

    do iedge = 1,4
!
!    jpoin2     v2
!      o<----------------o jpoin1 = ipoin2
!                       /^
!              auxvec  / |
!                     /  |
!                    o   | v1
!                  n0    |
!                        |
!                        |
!                        o ipoin1
!
      !
      ! Edge 1: ipoin1 ---> jpoin1 = ipoin2 ---> jpoin2
      !
      ipoin1 = iedge
      jpoin1 = mod(ipoin1, 4_ip) + 1_ip
      !
      edgvec1 = nodes(1:2, jpoin1) - nodes(1:2, ipoin1)
      !
      len1 = mag_vecmod(edgvec1)
      !
      ! Edge 2: ipoin2 ---> jpoin2
      !
      ipoin2 = jpoin1
      jpoin2 = mod(ipoin2, 4_ip) + 1_ip
      !
      edgvec2 = nodes(1:2, jpoin2) - nodes(1:2, ipoin2)
      !
      len2 = mag_vecmod(edgvec2)
      !
      auxvec = n0(1:2) - nodes(1:2, jpoin1)
      !
      weight1 = - dot_product(edgvec1, auxvec) / len1**2
      weight2 = dot_product(edgvec2, auxvec) / len2**2
      !
      if ( weight2 < 0.0_rp .or. weight2 > 1.0_rp &
        .or. weight1 < 0.0_rp .or. weight1 > 1.0_rp) call runend('mag_nedrec: point outside element')
      !
      basis(:, iedge) = edgvec1 * (1.0_rp - weight2)
      !
    end do
    !
  end function mag_nedrec
  !##############################################################


  !##############################################################
  subroutine mag_nedtet(nref, nodes, basis, curl, detJ)

    implicit none

    real(rp), intent(in) :: nref(:), nodes(:,:)

    real(rp), intent(out) :: basis(:,:), curl(:,:), detJ

    real(rp) :: &
      xi, eta, dseta,    &
      l(12),    &
      JT(3,3),    &
      grad_xi(3), grad_eta(3), grad_dseta(3),    &
      u(3), v(3), w(3)

    !
    ! Node coordinates in reference element
    !
    xi = nref(1); eta = nref(2); dseta = nref(3)
    !
    ! Element edge lengths
    !
    l(1) = mag_vecmod(nodes(:,2) - nodes(:,1))
    l(2) = mag_vecmod(nodes(:,3) - nodes(:,2))
    l(3) = mag_vecmod(nodes(:,1) - nodes(:,3))
    l(4) = mag_vecmod(nodes(:,4) - nodes(:,1))
    l(5) = mag_vecmod(nodes(:,4) - nodes(:,2))
    l(6) = mag_vecmod(nodes(:,4) - nodes(:,3))
    !
    ! Transposed Jacobian
    !
    JT = matmul(nodes(1:3, 1:4), mag_lagdif(nref(1:3), 30_ip))
    !
    ! Jacobian determinant
    !
    detJ = matdet(JT)
    !
    grad_xi = [JT(2,2) * JT(3,3) - JT(3,2) * JT(2,3), &
               JT(3,2) * JT(1,3) - JT(1,2) * JT(3,3), &
               JT(1,2) * JT(2,3) - JT(2,2) * JT(1,3)] / detJ
    !
    grad_eta = [JT(2,3) * JT(3,1) - JT(2,1) * JT(3,3), &
                JT(1,1) * JT(3,3) - JT(3,1) * JT(1,3), &
                JT(2,1) * JT(1,3) - JT(1,1) * JT(2,3)] / detJ
    !
    grad_dseta = [JT(2,1) * JT(3,2) - JT(2,2) * JT(3,1), &
                  JT(1,2) * JT(3,1) - JT(1,1) * JT(3,2), &
                  JT(1,1) * JT(2,2) - JT(1,2) * JT(2,1)] / detJ
    !
    basis(1:3,1) = l(1) * ( (1.0_rp - eta - dseta) * grad_xi + xi * grad_eta + xi * grad_dseta )
    basis(1:3,2) = l(2) * ( xi * grad_eta - eta * grad_xi )
    basis(1:3,3) = l(3) * ( -eta * grad_xi - (1.0_rp - xi - dseta) * grad_eta - eta * grad_dseta )
    basis(1:3,4) = l(4) * ( -dseta * grad_xi - dseta * grad_eta - (1.0_rp - xi - eta) * grad_dseta )
    basis(1:3,5) = l(5) * ( dseta * grad_xi - xi * grad_dseta )
    basis(1:3,6) = l(6) * ( dseta * grad_eta - eta * grad_dseta)
    !
    u = mag_vecpro(grad_xi, grad_eta)
    v = mag_vecpro(grad_eta, grad_dseta)
    w = mag_vecpro(grad_dseta, grad_xi)
    !
    curl(1:3,1) = 2.0_rp * l(1) * (u - w)
    curl(1:3,2) = 2.0_rp * l(2) * u
    curl(1:3,3) = 2.0_rp * l(3) * (u - v)
    curl(1:3,4) = 2.0_rp * l(4) * (-w + v)
    curl(1:3,5) = 2.0_rp * l(5) * w
    curl(1:3,6) = -2.0_rp * l(6) * v
    !
  end subroutine mag_nedtet
  !##############################################################


  !##############################################################
  subroutine mag_nedhex(nref, nodes, basis, curl, detJ)

    implicit none

    real(rp), intent(in) :: nref(:), nodes(:,:)

    real(rp), intent(out) :: basis(:,:), curl(:,:), detJ

    real(rp) :: &
      xi, eta, dseta,    &
      l(12),    &
      JT(3,3),    &
      grad_xi(3), grad_eta(3), grad_dseta(3),    &
      u(3), v(3), w(3)

    !
    ! Node coordinates in reference element
    !
    xi = nref(1); eta = nref(2); dseta = nref(3)
    !
    ! Element edge lengths
    !
    l(1) = mag_vecmod(nodes(:,2) - nodes(:,1))
    l(2) = mag_vecmod(nodes(:,4) - nodes(:,3))
    l(3) = mag_vecmod(nodes(:,6) - nodes(:,5))
    l(4) = mag_vecmod(nodes(:,8) - nodes(:,7))
    l(5) = mag_vecmod(nodes(:,4) - nodes(:,1))
    l(6) = mag_vecmod(nodes(:,3) - nodes(:,2))
    l(7) = mag_vecmod(nodes(:,8) - nodes(:,5))
    l(8) = mag_vecmod(nodes(:,7) - nodes(:,6))
    l(9) = mag_vecmod(nodes(:,5) - nodes(:,1))
    l(10) = mag_vecmod(nodes(:,6) - nodes(:,2))
    l(11) = mag_vecmod(nodes(:,8) - nodes(:,4))
    l(12) = mag_vecmod(nodes(:,7) - nodes(:,3))
    !
    ! Transposed Jacobian 
    !
    JT = matmul(nodes(1:3, 1:8), mag_lagdif(nref(1:3), 37_ip))
    !
    ! Jacobian determinant
    !
    detJ = matdet(JT)
    !
    grad_xi = [JT(2,2) * JT(3,3) - JT(3,2) * JT(2,3), &
               JT(3,2) * JT(1,3) - JT(1,2) * JT(3,3), &
               JT(1,2) * JT(2,3) - JT(2,2) * JT(1,3)] / detJ
    !
    grad_eta = [JT(2,3) * JT(3,1) - JT(2,1) * JT(3,3), &
                JT(1,1) * JT(3,3) - JT(3,1) * JT(1,3), &
                JT(2,1) * JT(1,3) - JT(1,1) * JT(2,3)] / detJ
    !
    grad_dseta = [JT(2,1) * JT(3,2) - JT(2,2) * JT(3,1), &
                  JT(1,2) * JT(3,1) - JT(1,1) * JT(3,2), &
                  JT(1,1) * JT(2,2) - JT(1,2) * JT(2,1)] / detJ
    !
    basis(1:3,1) = 0.125_rp * l(1) * (1.0_rp - eta) * (1.0_rp - dseta) * grad_xi
    basis(1:3,2) = 0.125_rp * l(2) * (1.0_rp + eta) * (1.0_rp - dseta) * grad_xi
    basis(1:3,3) = 0.125_rp * l(3) * (1.0_rp - eta) * (1.0_rp + dseta) * grad_xi
    basis(1:3,4) = 0.125_rp * l(4) * (1.0_rp + eta) * (1.0_rp + dseta) * grad_xi
    basis(1:3,5) = 0.125_rp * l(5) * (1.0_rp - xi) * (1.0_rp - dseta) * grad_eta
    basis(1:3,6) = 0.125_rp * l(6) * (1.0_rp + xi) * (1.0_rp - dseta) * grad_eta
    basis(1:3,7) = 0.125_rp * l(7) * (1.0_rp - xi) * (1.0_rp + dseta) * grad_eta
    basis(1:3,8) = 0.125_rp * l(8) * (1.0_rp + xi) * (1.0_rp + dseta) * grad_eta
    basis(1:3,9) = 0.125_rp * l(9) * (1.0_rp - xi) * (1.0_rp - eta) * grad_dseta
    basis(1:3,10) = 0.125_rp * l(10) * (1.0_rp + xi) * (1.0_rp - eta) * grad_dseta
    basis(1:3,11) = 0.125_rp * l(11) * (1.0_rp - xi) * (1.0_rp + eta) * grad_dseta
    basis(1:3,12) = 0.125_rp * l(12) * (1.0_rp + xi) * (1.0_rp + eta) * grad_dseta
    !
    u = mag_vecpro(grad_xi, grad_eta)
    v = mag_vecpro(grad_eta, grad_dseta)
    w = mag_vecpro(grad_dseta, grad_xi)
    !
    curl(1:3,1) = 0.125_rp * l(1) * ( (1.0_rp - dseta) * u - (1.0_rp - eta) * w )
    curl(1:3,2) = 0.125_rp * l(2) * ( - (1.0_rp - dseta) * u - (1.0_rp + eta) * w )
    curl(1:3,3) = 0.125_rp * l(3) * ( (1.0_rp + dseta) * u + (1.0_rp - eta) * w )
    curl(1:3,4) = 0.125_rp * l(4) * ( - (1.0_rp + dseta) * u + (1.0_rp + eta) * w )
    curl(1:3,5) = 0.125_rp * l(5) * ( - (1.0_rp - dseta) * u + (1.0_rp - xi) * v )
    curl(1:3,6) = 0.125_rp * l(6) * ( (1.0_rp - dseta) * u + (1.0_rp + xi) * v )
    curl(1:3,7) = 0.125_rp * l(7) * ( - (1.0_rp + dseta) * u - (1.0_rp - xi) * v )
    curl(1:3,8) = 0.125_rp * l(8) * ( (1.0_rp + dseta) * u - (1.0_rp + xi) * v )
    curl(1:3,9) = 0.125_rp * l(9) * ( (1.0_rp - eta) * w - (1.0_rp - xi) * v )
    curl(1:3,10) = 0.125_rp * l(10) * ( - (1.0_rp - eta) * w - (1.0_rp + xi) * v )
    curl(1:3,11) = 0.125_rp * l(11) * ( (1.0_rp + eta) * w + (1.0_rp - xi) * v )
    curl(1:3,12) = 0.125_rp * l(12) * ( - (1.0_rp + eta) * w + (1.0_rp + xi) * v )
    !
  end subroutine mag_nedhex
  !##############################################################


  !##############################################################
  subroutine mag_nedtri2(nref, nodes, basis, curl, detJ)

    implicit none

    real(rp), intent(in) :: nref(:), nodes(:,:)

    real(rp), intent(out) :: basis(:,:), curl(:), detJ

    real(rp) :: &
      xi, eta,    &
      l(3),    &
      grad_xi(2), grad_eta(2),    &
      JT(2,2)
    !
    ! Node coordinates in reference element
    !
    xi = nref(1); eta = nref(2)
    !
    ! Edge lengths
    !
    l(1) = mag_vecmod( nodes(1:2, 2) - nodes(1:2, 1) )
    l(2) = mag_vecmod( nodes(1:2, 3) - nodes(1:2, 2) )
    l(3) = mag_vecmod( nodes(1:2, 1) - nodes(1:2, 3) )
    !
    ! Transposed Jacobian
    !
    JT = matmul(nodes(1:2, 1:3), mag_lagdif(nref(1:2), 10_ip))
    !
    detJ = matdet(JT)
    !
    ! Gradients
    !
    grad_xi = [JT(2,2), -JT(1,2)] / detJ
    grad_eta = [-JT(2,1), JT(1,1)] / detJ
    !
    ! Basis
    ! 
    basis(1:2,1) = l(1) * ( (1.0_rp - eta) * grad_xi + xi * grad_eta )
    basis(1:2,2) = l(2) * ( xi * grad_eta - eta * grad_xi )
    basis(1:2,3) = l(3) * ( - eta * grad_xi - (1.0_rp - xi) * grad_eta )
    !
    ! Curls
    !
    curl(1:3) = l * 2.0_rp * (grad_xi(1) * grad_eta(2) - grad_xi(2) * grad_eta(1))
    !
  end subroutine mag_nedtri2
  !##############################################################


  !##############################################################
  subroutine mag_nedqua(nref, nodes, basis, curl, detJ)
    ! 
    ! This function receives the position in the reference element
    ! and the vertices of the quadrangle and returns the vector basis
    ! functions at nref
    !
    implicit none

    real(rp), intent(in) :: nref(:), nodes(:,:)

    real(rp), intent(out) :: basis(:,:), curl(:), detJ

    real(rp) :: &
      xi, eta,    &
      l(4),    &
      grad_xi(2), grad_eta(2),    &
      JT(2,2)
    !
    ! Node coordinates in reference element
    !
    xi = nref(1); eta = nref(2)
    !
    ! Element edge lengths
    !
    l(1) = mag_vecmod( nodes(:, 2) - nodes(:, 1) )
    l(2) = mag_vecmod( nodes(:, 4) - nodes(:, 3) )
    l(3) = mag_vecmod( nodes(:, 4) - nodes(:, 1) )
    l(4) = mag_vecmod( nodes(:, 3) - nodes(:, 2) )
    !
    ! Transposed Jacobian
    !
    JT = matmul(nodes(1:2, 1:4), mag_lagdif(nref(1:2), 12_ip))
    !
    ! Determinant of Jacobian
    !
    detJ = matdet(JT)
    !
    ! Gradients
    !
    grad_xi = [JT(2,2), -JT(1,2)] / detJ
    grad_eta = [-JT(2,1), JT(1,1)] / detJ
    !
    ! Basis
    !
    basis(1:2, 1) = 0.25_rp * (1 - eta) * grad_xi * l(1)
    basis(1:2, 2) = 0.25_rp * (1 + eta) * grad_xi * l(2)
    basis(1:2, 3) = 0.25_rp * (1 - xi) * grad_eta * l(3)
    basis(1:2, 4) = 0.25_rp * (1 + xi) * grad_eta * l(4)
    !
    ! Curls
    !
    curl(1:4) = 0.25_rp * l * (grad_xi(1) * grad_eta(2) - grad_xi(2) * grad_eta(1))
    curl(2) = -curl(2)
    curl(3) = -curl(3)
    !
  end subroutine mag_nedqua
  !##############################################################


  !##############################################################
  function mag_vecpro(u, v) result(w)

    implicit none

    real(rp), intent(in) :: &
      u(3),    &
      v(3)

    real(rp) :: &
      w(3)

    w(1) = u(2) * v(3) - u(3) * v(2)
    w(2) = u(3) * v(1) - u(1) * v(3)
    w(3) = u(1) * v(2) - u(2) * v(1)

  end function mag_vecpro
  !##############################################################


  !##############################################################
  function mag_vecpro2(u, v) result(w)

    implicit none

    real(rp), intent(in) :: &
      u(2),    &
      v(2)

    real(rp) :: &
      w

    w = u(1) * v(2) - u(2) * v(1)

  end function mag_vecpro2
  !##############################################################

end module mod_mag_nedele
