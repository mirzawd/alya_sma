!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Domain
!> @{
!> @file    mod_quadrature_gauss.f90
!> @author  houzeaux
!> @date    2020-04-24
!> @brief   Trapezoidal-rule quadrature
!> @details Trapezoidal-rule quadrature. This quadrature converges
!>          with 1/n^2 where n is the number of points in each direction.
!>          This is an upper bound.
!>   
!-----------------------------------------------------------------------

module mod_quadrature_trapezoidal

  use def_kintyp_basic, only : ip,rp
  implicit none
  private

  public :: trape_bar 
  public :: trape_qua    

contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-01-26
  !> @brief   Rule for bar elements
  !> @details This routine sets up the integration constants for
  !>          trapezoidal rule.
  !-----------------------------------------------------------------------

  subroutine trape_bar(ngaus,posgp,weigp,ierro)

    integer(ip), intent(in)  :: ngaus
    real(rp),    intent(out) :: posgp(ngaus)
    real(rp),    intent(out) :: weigp(ngaus)
    integer(ip), intent(out) :: ierro
    real(rp)                 :: delta
    integer(ip)              :: igaus

    ierro =  0

    select case ( ngaus )

    case ( 1_ip )

       posgp(1) =  0.0_rp
       weigp(1) =  2.0_rp

    case ( 2_ip )

       posgp(1) = -1.0_rp
       weigp(1) =  1.0_rp
       posgp(2) =  1.0_rp
       weigp(2) =  1.0_rp

    case default

       delta        =  2.0_rp/real(ngaus-1,rp)
       posgp(1)     = -1.0_rp
       weigp(1)     =  0.5_rp * delta
       posgp(ngaus) =  1.0_rp
       weigp(ngaus) =  0.5_rp * delta

       do igaus = 2,ngaus-1
          posgp(igaus) = -1.0_rp + real(igaus,rp) * delta
          weigp(igaus) = delta
       end do

    end select

  end subroutine trape_bar
     
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-01-26
  !> @brief   Rule for bar elements
  !> @details This routine sets up the integration constants for
  !>          trapezoidal rule.
  !-----------------------------------------------------------------------

  subroutine trape_qua(ndime,ngaus,posgp,weigp,ierro)

    integer(ip), intent(in)  :: ndime,ngaus
    real(rp),    intent(out) :: posgp(ndime,ngaus)
    real(rp),    intent(out) :: weigp(ngaus)
    integer(ip), intent(out) :: ierro
    real(rp)                 :: dx,dy,wx,wy
    real(rp)                 :: dz,wz
    integer(ip)              :: ii,jj,kk,nn,gg

    ierro =  0
    
    select case ( ngaus )

    case ( 1_ip )

       posgp(1:2,1) =  (/0.0_rp,0.0_rp/)
       weigp(1)     =  4.0_rp

    case default

       select case ( ndime )
          
       case ( 2_ip )
          
          nn = int(sqrt(real(ngaus,rp)),ip)
          if( nn*nn /= ngaus ) then
             ierro = 1
             return
          end if
          dx =  2.0_rp/real(nn-1,rp)
          dy =  2.0_rp/real(nn-1,rp)

          gg = 0

          do jj = 1,nn
             if( jj == 1 .or. jj == nn ) then
                wy = 0.5_rp
             else
                wy = 1.0_rp
             end if
             do ii = 1,nn
                if( ii == 1 .or. ii == nn ) then
                   wx = 0.5_rp
                else
                   wx = 1.0_rp
                end if
                gg = gg + 1
                posgp(1,gg) = -1.0_rp + real(ii-1,rp) * dx
                posgp(2,gg) = -1.0_rp + real(jj-1,rp) * dy
                weigp(gg)   =  wx * wy * dx * dy
             end do
          end do
          
       case ( 3_ip )
          
          nn = nint(real(ngaus,rp)**(1.0_rp/3.0_rp),ip)
          if( nn*nn*nn /= ngaus ) then
             ierro = 1
             return
          end if
          dx = 2.0_rp/real(nn-1,rp)
          dy = 2.0_rp/real(nn-1,rp)
          dz = 2.0_rp/real(nn-1,rp)
          gg = 0
          
          do kk = 1,nn
             if( kk == 1 .or. kk == nn ) then
                wz = 0.5_rp
             else
                wz = 1.0_rp
             end if
             do jj = 1,nn
                if( jj == 1 .or. jj == nn ) then
                   wy = 0.5_rp
                else
                   wy = 1.0_rp
                end if
                do ii = 1,nn
                   if( ii == 1 .or. ii == nn ) then
                      wx = 0.5_rp
                   else
                      wx = 1.0_rp
                   end if
                   gg          =  gg + 1
                   posgp(1,gg) = -1.0_rp + real(ii-1,rp) * dx
                   posgp(2,gg) = -1.0_rp + real(jj-1,rp) * dy
                   posgp(3,gg) = -1.0_rp + real(kk-1,rp) * dz
                   weigp(gg)   =  wx * wy * wz * dx * dy * dz
                end do
             end do
          end do

       end select
       
    end select

  end subroutine trape_qua
     
end module mod_quadrature_trapezoidal
!> @}
