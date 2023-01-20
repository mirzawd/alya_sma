!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Domain
!> @{
!> @file    def_chebyshev.f90
!> @author  houzeaux
!> @date    2020-04-24
!> @brief   Definitions used for Chebyshev
!> @details CDefinitions used for Chebyshev
!>   
!-----------------------------------------------------------------------

module def_chebyshev

  use def_kintyp_basic, only : ip,rp
  implicit none
  private  

  integer(ip), target    :: atoIJK3(64) = &
       [1,4,16,15,2,3,11,12,9,14,33,36,10,13,34,35,     &
       5,8,32,31,6,7,27,28,25,30,53,56,26,29,54,55,     &
       17,20,44,43,18,19,39,40,37,42,57,60,38,41,58,59, &
       21,24,52,51,22,23,47,48,45,50,61,64,46,49,62,63  ]

  integer(ip), target    :: atoIJ3(16) = &
       [1,4,12,11,2,3,7,8,5,10,13,16,6,9,14,15]

  public :: chebyshev_roots
  public :: lagrange_roots
  public :: DoubleTensorProduct
  public :: TripleTensorProduct
  public :: var_interpolate
  public :: atoIJK3
  public :: atoIJ3

contains

  !
  ! Computes the roots of Tn(cos(y)) = cos(ny)       !
  ! in the interval [-1,1], and includes the         !
  ! endpoints. Roots are given in the following      !
  ! order:                                           !
  !  - xi(1) = -1                                    !
  !  - xi(2) =  1                                    !
  !    xi(3:porder+1) = -cos(pi*[1:porder-1]/porder) !
  !

  pure subroutine chebyshev_roots(porder,xi_chb)

    integer(ip), intent(in)  :: porder
    real(rp),    intent(out) :: xi_chb(porder+1)

    !do i = 1,porder+1
    !   xi_chb(i) = -cos(v_pi*dble(i-1)/dble(porder))
    !end do

    select case ( porder )
    case ( 3_ip ) 
       xi_chb(1) = -1.0_rp
       xi_chb(2) = -0.447213595499958_rp
       xi_chb(3) =  0.447213595499958_rp
       xi_chb(4) =  1.0_rp
    end select

  end subroutine chebyshev_roots

  !
  ! Computes the equispaced loci for the Lagrange    !
  ! poly in the interval [-1,1], and includes the    !
  ! endpoints. Roots are given in the following      !
  ! order:                                           !
  !  - xi(1) = -1                                    !
  !  - xi(2) =  1                                    !
  !    xi(3:porder+1) = xi(1)+[1:N-1]*(2/N)          !
  ! 

  pure subroutine lagrange_roots(porder,xi_lag)


    integer(ip), intent(in)  :: porder
    real(rp),    intent(out) :: xi_lag(porder+1)
    integer(ip)              :: i

    do i = 1,porder+1
       xi_lag(i) = -1.0_rp+(2.0_rp*real(i-1,rp)/real(porder,rp))
    end do

  end subroutine lagrange_roots

  pure subroutine TripleTensorProduct(nnode,porder,xi_grid,s,t,z,N,dN,dlxigp_ip)

    integer(ip),           intent(in)     :: nnode
    integer(ip),           intent(in)     :: porder
    real(rp),              intent(in)     :: s, t, z, xi_grid(porder+1)
    real(rp),              intent(out)    :: N(nnode)
    real(rp),    optional, intent(out)    :: dN(3,nnode)
    real(rp),    optional, intent(out)    :: dlxigp_ip(3,porder+1)
    integer(ip)                           :: i, j, k, c
    real(rp),    dimension(porder+1)      :: lxi_ip, leta_ip, lzeta_ip
    real(rp),    dimension(porder+1)      :: dlxi_ip, dleta_ip, dlzeta_ip

    call eval_lagrangePoly     (porder,xi_grid,s,lxi_ip)
    call eval_lagrangePoly     (porder,xi_grid,t,leta_ip)
    call eval_lagrangePoly     (porder,xi_grid,z,lzeta_ip)
    call eval_lagrangePolyDeriv(porder,xi_grid,s,dlxi_ip)
    call eval_lagrangePolyDeriv(porder,xi_grid,t,dleta_ip)
    call eval_lagrangePolyDeriv(porder,xi_grid,z,dlzeta_ip)

    if(present(dlxigp_ip)) then
       dlxigp_ip(1,:) = dlxi_ip(:)
       dlxigp_ip(2,:) = dleta_ip(:)
       dlxigp_ip(3,:) = dlzeta_ip(:)
    end if

    c = 0
    do k = 1,porder+1
       do i = 1,porder+1
          do j = 1,porder+1
             c = c+1
             N(atoIJK3(c))    = lxi_ip(i)*leta_ip(j)*lzeta_ip(k)
             if( present(dN) ) then
                dN(1,atoIJK3(c)) = dlxi_ip(i)*leta_ip(j)*lzeta_ip(k)
                dN(2,atoIJK3(c)) = lxi_ip(i)*dleta_ip(j)*lzeta_ip(k)
                dN(3,atoIJK3(c)) = lxi_ip(i)*leta_ip(j)*dlzeta_ip(k)
             end if
          end do
       end do
    end do

  end subroutine TripleTensorProduct

  pure subroutine DoubleTensorProduct(nnode,porder,xi_grid,s,t,N,dN,dlxigp_ip)

    integer(ip),           intent(in)     :: nnode
    integer(ip),           intent(in)     :: porder
    real(rp),              intent(in)     :: s, t, xi_grid(porder+1)
    real(rp),              intent(out)    :: N(nnode)
    real(rp),    optional, intent(out)    :: dN(2,nnode)
    real(rp),    optional, intent(out)    :: dlxigp_ip(2,porder+1)
    integer(ip)                           :: i, j, c
    real(rp),    dimension(porder+1)      :: lxi_ip, leta_ip
    real(rp),    dimension(porder+1)      :: dlxi_ip, dleta_ip

    call eval_lagrangePoly     (porder,xi_grid,s,lxi_ip)
    call eval_lagrangePoly     (porder,xi_grid,t,leta_ip)
    call eval_lagrangePolyDeriv(porder,xi_grid,s,dlxi_ip)
    call eval_lagrangePolyDeriv(porder,xi_grid,t,dleta_ip)

    if(present(dlxigp_ip)) then
       dlxigp_ip(1,:) = dlxi_ip(:)
       dlxigp_ip(2,:) = dleta_ip(:)
    end if

    c = 0
    do i = 1,porder+1
       do j = 1,porder+1
          c = c+1
          N(atoIJ3(c))    = lxi_ip(i)*leta_ip(j)
          if( present(dN) ) then
             dN(1,atoIJ3(c)) = dlxi_ip(i)*leta_ip(j)
             dN(2,atoIJ3(c)) = lxi_ip(i)*dleta_ip(j)
          end if
       end do
    end do

  end subroutine DoubleTensorProduct

  pure subroutine eval_lagrangePoly(porder,xi,xi_p,l_ip)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Evaluates the Lagrange poly of order N at a      !
    ! location xi(p), where xi E [-1,1]. Returns an    !
    ! array with all possible values l_ip can assume.  !
    ! Remark that l_ip = 1 if p==i and 0 if            !
    ! p==j~=i.                                         !
    ! Expects a grid of the form:                      !
    !                                                  !
    ! -1                        1 xi                   !
    !  o----v----v----v----v----o                      !
    !  0    2    3    4    5    1 i                    !
    !                                                  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    integer(ip), intent(in)  :: porder
    real(rp),    intent(in)  :: xi(porder+1), xi_p
    real(rp),    intent(out) :: l_ip(porder+1)
    integer(ip)              :: i, j, lorder(porder+1)

    lorder(1) = 1
    lorder(2) = porder+1
    do i = 3,porder+1
       lorder(i) = i-1
    end do
    do i = 0,porder ! Nodal loop
       l_ip(i+1) = 1.0_rp
       do j = 0,porder ! Product series
          if (j .ne. (lorder(i+1)-1)) then
             l_ip(i+1) = l_ip(i+1)*((xi_p-xi(j+1))/(xi(lorder(i+1))-xi(j+1)))
          end if
       end do
    end do

  end subroutine eval_lagrangePoly

  !
  ! Evaluates the derivative of the Lagrange poly    !
  ! of order N at a location xi(p), where            !
  ! xi E [-1,1]. Returns an array with all possible  !
  ! values dl_ip can assume.                         !
  ! Expects a grid of the form:                      !
  !                                                  !
  ! -1                        1 xi                   !
  !  o----v----v----v----v----o                      !
  !  0    2    3    4    5    1 i                    !
  !                                                  !
  !

  pure subroutine eval_lagrangePolyDeriv(porder,xi,xi_p,dl_ip)

    integer(ip), intent(in)  :: porder
    real(rp),    intent(in)  :: xi(porder+1), xi_p
    real(rp),    intent(out) :: dl_ip(porder+1)
    integer(ip)              :: i, j, m, lorder(porder+1)
    real(rp)                 :: aux

    lorder(1) = 1
    lorder(2) = porder+1
    do i = 3,porder+1
       lorder(i) = i-1
    end do
    do i = 0,porder ! Nodal loop
       dl_ip(i+1) = 0.0_rp
       do j = 0,porder ! Product series
          aux = 1.0_rp
          if (j .ne. (lorder(i+1)-1)) then
             do m = 0,porder
                if (m .ne. j .and. m .ne. lorder(i+1)-1) then
                   aux = aux * (xi_p-xi(m+1))/(xi(lorder(i+1))-xi(m+1))
                end if
             end do
             dl_ip(i+1) = dl_ip(i+1) + (1.0_rp/(xi(lorder(i+1))-xi(j+1)))*aux
          end if
       end do
    end do

  end subroutine eval_lagrangePolyDeriv

  !
  ! Interpolates a variable using element shape      !
  ! functions N(xi,eta,zeta). Given a coordinate     !
  ! set X = (x,y,z), v(X) = N_a*v_a, where v_a are   !
  ! element nodal values.
  ! Expects a grid of the form:                      !
  !                                                  !
  ! -1                        1 xi                   !
  !  o----v----v----v----v----o                      !
  !  0    2    3    4    5    1 i                    !
  !                                                  !
  !

  pure subroutine var_interpolate(nnode,var,Neval,var_a)

    integer(ip), intent(in)  :: nnode
    real(rp),    intent(in)  :: var(nnode), Neval(nnode)
    real(rp),    intent(out) :: var_a
    integer(ip)              :: inode

    var_a = 0.0_rp
    do inode = 1,nnode
       var_a = var_a+Neval(inode)*var(inode)
    end do

  end subroutine var_interpolate

end module def_chebyshev
!> @}
