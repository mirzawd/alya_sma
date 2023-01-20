!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup CPU_Time
!> @{
!> @file    variability.f90
!> @author  houzeaux
!> @date    2019-04-02
!> @brief   Variability
!> @details Subroutine used for variability study
!> @} 
!-----------------------------------------------------------------------

subroutine variability(time1)

  use def_kintyp,         only : ip,rp
  use mod_communications, only : PAR_MAX
  use mod_communications, only : PAR_AVERAGE
 implicit none
  real(rp), intent(in) :: time1
  integer(ip)          :: rank_max
  real(rp)             :: timem,timea,timei,timef,timed
  integer(ip), save    :: ipass = 0

  call cputim(timei)
  timem = time1
  timea = time1
  call PAR_MAX(timem,rank_max_owner=rank_max)
  call PAR_AVERAGE(timea)

  call cputim(timef)
  timed = timef-timei
  timei = timef
  if( ipass > 0 ) write(90,10) timed,timea,timem,rank_max

  ipass = ipass + 1

10 format(3(1x,e13.6),i4)
  
end subroutine variability
