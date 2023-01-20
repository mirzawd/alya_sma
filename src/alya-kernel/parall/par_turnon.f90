!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Parall
!> @{
!> @file    par_turnon.f90
!> @author  houzeaux
!> @date    2020-05-08
!> @brief   Start Parall 
!> @details Turn on parall
!> @} 
!-----------------------------------------------------------------------

subroutine par_turnon()
  use def_kintyp_basic, only : rp
  use mod_outfor,       only : outfor
  use def_parall
  use mod_parall
  real(rp) :: time1,time2

  call cputim(time1)
  !
  ! Broadcast data
  !
  call par_inidat()  
  call par_sendat(1_ip)                                
  !
  ! Open files
  !
  call par_openfi(0_ip)
  call par_openfi(1_ip)
  call outfor(27_ip,lun_outpu_par,' ')
  
  call cputim(time2)
  cpu_paral(2)=time2-time1

end subroutine par_turnon
