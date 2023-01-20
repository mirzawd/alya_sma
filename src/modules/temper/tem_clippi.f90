!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Temper
!> @{
!> @file    tem_clippi.f90
!> @author  houzeaux
!> @date    2018-12-28
!> @brief   This routine performs several types of updates for the thermal variable (T or h)
!> @details Clipping
!> @}
!-----------------------------------------------------------------------

subroutine tem_clippi()

  use def_parame
  use def_master
  use def_domain
  use def_temper
  implicit none
  integer(ip) :: ipoin

  if (kfl_negat_tem==1) then
     do ipoin = 1,nunkn_tem
        if( unkno(ipoin) < negat_tem ) then
           unkno(ipoin) = negat_tem
        end if
     end do
  endif
  if (kfl_posit_tem==1) then
     do ipoin = 1,nunkn_tem
        if( unkno(ipoin) > posit_tem ) then
           unkno(ipoin) = posit_tem
        end if
     end do
  endif

end subroutine tem_clippi
