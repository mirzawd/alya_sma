!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



module mod_mag_linalg

  use def_kintyp, only: ip, rp

  implicit none

contains

  !##############################################################
  function mag_vecmod(a) result(na)

    implicit none

    real(rp), intent(in) ::    a(:)

    real(rp) :: na

    na = SQRT ( dot_product(a, a) )

  end function mag_vecmod
  !##############################################################

  
  !##############################################################
  function mag_vecvec(a, b) result(c)

    implicit none

    real(rp), intent(in) ::    a(:), b(:)

    real(rp) ::  c(size(a))

    integer(ip) ::  na, nb

    na = size(a); nb = size(b)

!    if (na /= nb) stop "mag_vecvec: vector dimensions do not match"

    if (na == 3_ip .and. nb == 3_ip) then
      c(1) = a(2) * b(3) - a(3) * b(2)
      c(2) = a(3) * b(1) - a(1) * b(3)
      c(3) = a(1) * b(2) - a(2) * b(1)
    elseif (na == 2_ip .and. nb == 3_ip) then
      c(1) = a(2) * b(3) 
      c(2) = - a(1) * b(3)
      c(3) = a(1) * b(2) - a(2) * b(1)
    elseif (na == 3_ip .and. nb == 2_ip) then
      c(1) = - a(3) * b(2)
      c(2) = a(3) * b(1)
      c(3) = a(1) * b(2) - a(2) * b(1)
    elseif (na == 2_ip .and. nb == 2_ip) then
      c(1) = 0.0_rp
      c(2) = 0.0_rp
      c(3) = a(1) * b(2) - a(2) * b(1)
    else
      stop "mag_vecvec: check vector dimensions"
    end if

  end function mag_vecvec
  !##############################################################
  

  !##############################################################
  function mag_matvec(A, x) result(y)

    implicit none

    real(rp), intent(in) ::             A(:,:)

    real(rp), intent(in) ::             x(:)

    real(rp) ::            y(size(A, 1))

    integer(ip) ::                      i, j, nr, nc

    nr = size(A, 1); nc = size(A, 2)

    do i = 1, nr
      y(i) = 0.0_rp
      do j = 1, nc
        y(i) = y(i) + A(i, j) * x(j)
      end do
    end do

  end function mag_matvec
  !##############################################################


  !##############################################################
  function mag_vecmat(x, A) result(y)

    implicit none

    real(rp), intent(in) ::             A(:,:)

    real(rp), intent(in) ::             x(:)

    real(rp) ::            y(size(A, 2))

    integer(ip) ::                      i, j, nr, nc

    nr = size(A, 1); nc = size(A, 2)

    do j = 1, nc
      y(j) = 0.0_rp
      do i = 1, nr
        y(j) = y(j) + x(i) * A(i, j)
      end do
    end do

  end function mag_vecmat
  !##############################################################


  !##############################################################
  subroutine lupp(A, P, STAT)

    implicit none

    real(rp), intent(inout) :: A(:,:)

    real(rp), intent(out) :: P(:,:)

    integer(ip) :: m, n, i, j, k

    integer, intent(out), optional :: STAT

    real(rp) :: piv

    real(rp), dimension(:) :: aswap(size(A, 2)), pswap(size(A, 2))

    ! Everything goes just fine...
    if (present(STAT)) STAT = 0

    m = size(A, 1); n = size(A, 2)

    if (m /= n) then
      ! Matrix is not square!
      if (present(STAT)) STAT = 1
    end if

    ! Permutation matrix
    P = eye(n)

    do k = 1, n

      j = k - 1 + maxindex(abs(A(k:n, k)))

      piv = A(j, k)

      if (piv == 0) then
        ! Matrix is singular
        if (present(STAT)) STAT = 2
        exit
      end if

      ! Swap rows
      aswap = A(k, :); A(k, :) = A(j, :); A(j, :) = aswap
      pswap = P(k, :); P(k, :) = P(j, :); P(j, :) = pswap

      ! Update rows in matrix
      do i = k+1, n
        A(i, k+1:n) = A(i, k+1:n) - A(i, k) / piv * A(k, k+1:n)
        A(i, k) = A(i, k) / piv
      end do

    end do

  end subroutine lupp
  !##############################################################


  !##############################################################
  subroutine lusolver(x, A, b, STAT)

    implicit none

    real(rp), intent(out) :: x(:)

    real(rp), intent(in) :: A(:,:)

    real(rp), intent(in) :: b(:)

    integer, intent(out), optional :: STAT

    real(rp), dimension(:) ::      y(size(x))

    real(rp), dimension(:,:) ::    LU(size(A,1), size(A,2)),     &
                                   P(size(A,1), size(A,2))

    integer(ip) ::                              m,      &
                                                n,      &
                                                i

    integer ::                                  STATUS

    m = size(A, 1); n = size(A, 2)

    ! LU decomposition
    LU = A
    call lupp(LU, P, STATUS)

    ! Re-order RHS
    y = matmul(P, b)

    ! Re-order RHS
    y = matmul(P, b)

    ! Solve lower triangular system
    do i = 1, n
      y(i) = y(i) - dot_product(LU(i, 1:i-1), y(1:i-1))
    end do

    x = y
    ! Solve upper triangular system
    do i = n, 1, -1
      x(i) = (x(i) - dot_product(LU(i, i+1:n), x(i+1:n))) / LU(i, i)
    end do

!    STAT = STATUS

  end subroutine lusolver
  !##############################################################


  !##############################################################
  function eye(n)

    implicit none

    integer(ip), intent(in) :: n

    real(rp), dimension(n, n) :: eye

    integer(ip) :: i

    eye = 0_rp
    do i = 1, n
      eye(i, i) = 1_rp
    end do

  end function eye
  !##############################################################


  !##############################################################
  function maxindex(v)

    implicit none

    real(rp), intent(in) :: v(:)

    integer(ip) :: maxindex, i, n
    real(rp) :: maxvalue

    n = size(v)

    if (n == 0) then
      maxindex = 0
      return
    end if

    maxvalue = v(1)
    maxindex = 1
    do i = 2, n
      if (v(i) > maxvalue) then
        maxvalue = v(2)
        maxindex = i
      end if
    end do

  end function maxindex
  !##############################################################


  !##############################################################
  function matdet(A) result(det)

    implicit none

    real(rp), intent(in) :: A(:,:)

    real(rp) :: det

    integer(ip) :: nr, nc

    nr = size(A, 1); nc = size(A, 2)

    if (nr /= nc) stop "matdet: number of rows and columns must be the same"

    if (nr == 1) then
      det = A(1,1)
    elseif (nr == 2) then
      det = A(1,1) * A(2,2) - A(1,2) * A(2,1)
    elseif (nr == 3) then
      det = A(1,1) * A(2,2) * A(3,3) + A(2,1) * A(3,2) * A(1,3) + A(1,2) * A(2,3) * A(3,1) &
          - A(1,3) * A(2,2) * A(3,1) - A(2,3) * A(3,2) * A(1,1) - A(1,2) * A(2,1) * A(3,3)
    else
      stop "matdet: to be programmed"
    end if

  end function matdet
  !##############################################################


  !##############################################################
!  function matinv3(A) result(invA)
!
!    implicit none
!
!    real(rp), intent(in) :: A(3,3)
!
!    real(rp) :: invA(3,3)
!
!    real(rp) :: detA
!
!    detA = matdet(A)
!
!    invA = 
!
!  end function matinv3
  !##############################################################
    
end module mod_mag_linalg
