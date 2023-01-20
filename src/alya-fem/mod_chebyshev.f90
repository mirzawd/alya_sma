!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Kinds_and_types
!> @{
!> @file    mod_chebyshev.f90
!> @author  houzeaux
!> @date    2020-04-04
!> @brief   Iso-parametric arrays
!> @details Chebyshev shape function, derivative and Hessian
!-----------------------------------------------------------------------

module mod_chebyshev

  use def_kintyp_basic, only : ip,rp
  use def_chebyshev,    only : chebyshev_roots
  use def_chebyshev,    only : TripleTensorProduct
  use def_chebyshev,    only : DoubleTensorProduct

  private

  public :: cheby2
  public :: cheby3

contains

  !-----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @brief   3D shape functions
  !> @details Evaluate shape functions, derivative and Hessian
  !>          for 1-d continuos with 2 3 & 4 nodes
  !>
  !>          HEXAHEDRA:   64  nodes
  !
  !
  !-----------------------------------------------------------------------

  pure subroutine cheby2(x,nnode,N,dN,dN2,ierro)

    integer(ip), intent(in)            :: nnode
    real(rp),    intent(in)            :: x(2)
    real(rp),    intent(out)           :: N(nnode)
    real(rp),    intent(out)           :: dN(2,nnode)
    real(rp),    intent(out), optional :: dN2(3,nnode)
    integer(ip), intent(out), optional :: ierro
    integer(ip)                        :: ii,jj,ierr
    integer(ip), parameter             :: porder=3
    real(rp)                           :: xi, eta
    real(rp)                           :: dlxigp_ip(2,porder+1)
    real(rp)                           :: xi_grid(porder+1)

    ierr = 0
    xi   = x(1)
    eta  = x(2)

    if( present(dN2) ) then
       do ii = 1,nnode
          do jj = 1,3
             dN2(jj,ii) = 0.0_rp
          end do
       end do
    end if 
    do ii = 1,nnode
       do jj = 1,2
          dN(jj,ii) = 0.0_rp
       end do
    end do

    if( nnode == 16 ) then

       call chebyshev_roots    (porder,xi_grid)
       call DoubleTensorProduct(nnode,porder,xi_grid,xi,eta,N,dN,dlxigp_ip)

    else   

       ierr = 1

    end if

    if( present(ierro) ) ierro = ierr

  end subroutine cheby2
  
  !-----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @brief   3D shape functions
  !> @details Evaluate shape functions, derivative and Hessian
  !>          for 1-d continuos with 2 3 & 4 nodes
  !>
  !>          HEXAHEDRA:   64  nodes
  !
  !
  !-----------------------------------------------------------------------

  pure subroutine cheby3(x,nnode,N,dN,dN2,ierro)

    integer(ip), intent(in)            :: nnode
    real(rp),    intent(in)            :: x(3)
    real(rp),    intent(out)           :: N(nnode)
    real(rp),    intent(out)           :: dN(3,nnode)
    real(rp),    intent(out), optional :: dN2(6,nnode)
    integer(ip), intent(out), optional :: ierro
    integer(ip)                        :: ii,jj,ierr
    integer(ip), parameter             :: porder=3
    real(rp)                           :: xi, eta, zeta
    real(rp)                           :: dlxigp_ip(3,porder+1)
    real(rp)                           :: xi_grid(porder+1)

    ierr = 0
    xi   = x(1)
    eta  = x(2)
    zeta = x(3)

    if( present(dN2) ) then
       do ii = 1,nnode
          do jj = 1,6
             dN2(jj,ii) = 0.0_rp
          end do
       end do
    end if 
    do ii = 1,nnode
       do jj = 1,3
          dN(jj,ii) = 0.0_rp
       end do
    end do

    if( nnode == 64 ) then

       call chebyshev_roots(porder,xi_grid)
       call TripleTensorProduct(nnode,porder,xi_grid,xi,eta,zeta,N,dN,dlxigp_ip)

    else   

       ierr = 1

    end if

    if( present(ierro) ) ierro = ierr

  end subroutine cheby3

end module mod_chebyshev
!> @}
