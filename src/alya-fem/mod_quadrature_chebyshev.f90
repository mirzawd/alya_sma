!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Domain
!> @{
!> @file    mod_quadrature_chebyshev.f90
!> @author  houzeaux
!> @date    2020-04-24
!> @brief   Close quadrature
!> @details Close quadrature rule
!>   
!-----------------------------------------------------------------------

module mod_quadrature_chebyshev

  use def_kintyp_basic, only : ip,rp
  use def_chebyshev,    only : chebyshev_roots
  use def_chebyshev,    only : atoIJK3
  use def_chebyshev,    only : atoIJ3

  implicit none
  private
  integer(ip), parameter :: morder=10

  public :: cheby_qua

contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2022-02-09
  !> @brief   Closed rule for hex
  !> @details This routine sets up the integration constants of closed 
  !>          integration rules for brick elements:
  !>
  !>               NDIME = 1             NDIME = 2             NDIME = 3
  !>
  !>           NGAUS  EXACT POL.     NGAUS  EXACT POL.     NGAUS  EXACT POL. 
  !>           -----  ----------     -----  ----------     -----  ----------  
  !>             2      q1           2 x 2     q1          2x2x2     q1   
  !>             3      q3           3 x 3     q3          3x3x3     q3
  !>             4      q3           4 x 4     q3          4x4x4     q3
  !>                                   8       p2            20      p2 
  !>
  !-----------------------------------------------------------------------

  subroutine cheby_qua(ndime,ngaus,posgp,weigp,ierro)

    integer(ip), intent(in)  :: ndime,ngaus
    integer(ip), intent(out) :: ierro
    real(rp),    intent(out) :: posgp(ndime,ngaus),weigp(ngaus)
    real(rp)                 :: q, sl, s0, sr, ws, w0
    integer(ip)              :: inode, i, j, k, lorder(morder+1)
    integer(ip)              :: porder
    real(rp)                 :: xi(morder+1),w1d(morder+1)

    ierro=0
    !
    ! 2D 
    !
    if( ndime == 2 ) then

       select case ( ngaus )

       case ( 1_ip )

          posgp(1,1) = 0.0_rp
          posgp(1,2) = 0.0_rp
          posgp(1,3) = 0.0_rp
          weigp(1)   = 4.0_rp

       case ( 8_ip )

          q          = 1.0_rp/sqrt(3.0_rp)
          posgp(1,1) = -q ! xi
          posgp(1,2) = -q ! eta
          posgp(1,3) = -q ! zeta
          posgp(2,1) = -q
          posgp(2,2) = -q
          posgp(2,3) =  q
          posgp(3,1) =  q
          posgp(3,2) = -q
          posgp(3,3) =  q
          posgp(4,1) =  q
          posgp(4,2) = -q
          posgp(4,3) = -q
          posgp(5,1) = -q ! xi
          posgp(5,2) =  q ! eta
          posgp(5,3) = -q ! zeta
          posgp(6,1) = -q
          posgp(6,2) =  q
          posgp(6,3) =  q
          posgp(7,1) =  q
          posgp(7,2) =  q
          posgp(7,3) =  q
          posgp(8,1) =  q
          posgp(8,2) =  q
          posgp(8,3) = -q
          weigp(1) = 1.0_rp
          weigp(2) = 1.0_rp
          weigp(3) = 1.0_rp
          weigp(4) = 1.0_rp
          weigp(5) = 1.0_rp
          weigp(6) = 1.0_rp
          weigp(7) = 1.0_rp
          weigp(8) = 1.0_rp

       case ( 27_ip) 

          sl            = -sqrt(0.6_rp)
          s0            = 0.0_rp
          sr            = sqrt(0.6_rp)
          posgp( 1,1:3) = [sl,sl,sl]
          posgp( 2,1:3) = [sl,sl,s0]
          posgp( 3,1:3) = [sl,sl,sr]
          posgp( 4,1:3) = [s0,sl,sl]
          posgp( 5,1:3) = [s0,sl,s0]
          posgp( 6,1:3) = [s0,sl,sr]
          posgp( 7,1:3) = [sr,sl,sl]
          posgp( 8,1:3) = [sr,sl,s0]
          posgp( 9,1:3) = [sr,sl,sr]
          posgp(10,1:3) = [sl,s0,sl]
          posgp(11,1:3) = [sl,s0,s0]
          posgp(12,1:3) = [sl,s0,sr]
          posgp(13,1:3) = [s0,s0,sl]
          posgp(14,1:3) = [s0,s0,s0]
          posgp(15,1:3) = [s0,s0,sr]
          posgp(16,1:3) = [sr,s0,sl]
          posgp(17,1:3) = [sr,s0,s0]
          posgp(18,1:3) = [sr,s0,sr]
          posgp(19,1:3) = [sl,sr,sl]
          posgp(20,1:3) = [sl,sr,s0]
          posgp(21,1:3) = [sl,sr,sr]
          posgp(22,1:3) = [s0,sr,sl]
          posgp(23,1:3) = [s0,sr,s0]
          posgp(24,1:3) = [s0,sr,sr]
          posgp(25,1:3) = [sr,sr,sl]
          posgp(26,1:3) = [sr,sr,s0]
          posgp(27,1:3) = [sr,sr,sr]

          ws        = 5.0_rp/9.0_rp
          w0        = 8.0_rp/9.0_rp
          weigp( 1) = ws*ws*ws
          weigp( 2) = ws*ws*w0
          weigp( 3) = ws*ws*ws
          weigp( 4) = w0*ws*ws
          weigp( 5) = w0*ws*w0
          weigp( 6) = w0*ws*ws
          weigp( 7) = ws*ws*ws
          weigp( 8) = ws*ws*w0
          weigp( 9) = ws*ws*ws
          weigp(10) = ws*w0*ws
          weigp(11) = ws*w0*w0
          weigp(12) = ws*w0*ws
          weigp(13) = w0*w0*ws
          weigp(14) = w0*w0*w0
          weigp(15) = w0*w0*ws
          weigp(16) = ws*w0*ws
          weigp(17) = ws*w0*w0
          weigp(18) = ws*w0*ws
          weigp(19) = ws*ws*ws
          weigp(20) = ws*ws*w0
          weigp(21) = ws*ws*ws
          weigp(22) = w0*ws*ws
          weigp(23) = w0*ws*w0 
          weigp(24) = w0*ws*ws
          weigp(25) = ws*ws*ws
          weigp(26) = ws*ws*w0
          weigp(27) = ws*ws*ws

       case ( 16_ip )

          porder = 3
          call chebyshev_roots(porder,xi)
          lorder(1) = 1
          lorder(2) = porder+1
          do i = 3,porder+1
             lorder(i) = i-1
          end do

          if (porder == 3) then
             w1d(1:4) = [1.0_rp/6.0_rp, 5.0_rp/6.0_rp, 5.0_rp/6.0_rp, 1.0_rp/6.0_rp]
          else
             STOP 1 
          end if

          select case ( porder )
             
          case ( 3_ip )

             inode = 0
             do i = 1,porder+1
                do j = 1,porder+1
                   inode = inode + 1
                   posgp(1:2,atoIJ3(inode)) = [xi(lorder(i)), xi(lorder(j))]
                   weigp(atoIJ3(inode))     = w1d(lorder(i))*w1d(lorder(j))
                end do
             end do
             
          end select

       case default

          ierro = 1
          posgp = 0.0_rp
          weigp = 0.0_rp

       end select

    else if( ndime == 3 ) then

       select case ( ngaus )

       case ( 64_ip )

          porder = 3
          call chebyshev_roots(porder,xi)
          lorder(1) = 1
          lorder(2) = porder+1
          do i = 3,porder+1
             lorder(i) = i-1
          end do

          if (porder == 3) then
             w1d(1:4) = [1.0_rp/6.0_rp, 5.0_rp/6.0_rp, 5.0_rp/6.0_rp, 1.0_rp/6.0_rp]
          else
             STOP 1 
          end if

          select case ( porder )
             
          case ( 3_ip )
             
             inode = 0
             do k = 1,porder+1
                do i = 1,porder+1
                   do j = 1,porder+1
                      inode = inode + 1
                      posgp(1:3,atoIJK3(inode)) = [xi(lorder(i)), xi(lorder(j)), xi(lorder(k))]
                      weigp(atoIJK3(inode))     = w1d(lorder(i))*w1d(lorder(j))*w1d(lorder(k))
                   end do
                end do
             end do
             
          end select
          

       case default

          ierro = 1

       end select

    end if

  end subroutine cheby_qua

end module mod_quadrature_chebyshev
!> @}

