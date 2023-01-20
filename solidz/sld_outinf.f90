!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    sld_outinf.f90
!> @author  Solidz team
!> @date    May 2020
!> @brief   
!> @details 
!> @}
!-----------------------------------------------------------------------

subroutine sld_outinf(itask)
  
  use def_kintyp_basic,      only : ip, rp, lg
  use def_master,            only : INOTSLAVE, kfl_rstar
  use mod_sld_cardiac_cycle, only : kfl_cardiac_cycle, sld_cardiac_cycle_write_res
  use def_solidz,            only : lun_sysnet_heart_res_sld, lun_sysnet_system_res_sld
  use mod_sld_post_reaction, only : kfl_preac_sld, sld_post_reaction_write_res
#ifndef PROPER_ELEM_PRIVATE_OFF
  use mod_ker_sysnet, only: kfl_sysnet, sysnet_write_alya_res
#endif
  
  implicit none
  
  integer(ip), intent(in) :: itask

  if( INOTSLAVE ) then

     select case(itask)

     case( 1_ip )
        !
        ! Write heading files
        !
        if( kfl_rstar /= 2 ) then
          if( kfl_cardiac_cycle ) call sld_cardiac_cycle_write_res(1_ip)
          if( kfl_preac_sld     ) call sld_post_reaction_write_res(1_ip)

#ifndef PROPER_ELEM_PRIVATE_OFF
          if(kfl_sysnet) call sysnet_write_alya_res(1_ip, lun_sysnet_heart_res_sld)     ! *-sysnetHeart.sld.res 
          if(kfl_sysnet) call sysnet_write_alya_res(2_ip, lun_sysnet_system_res_sld)    ! *-sysnetSystem.sld.res
#endif
       end if

     case( 2_ip )
        !
        ! Write results at the end of the step
        !
        if( kfl_cardiac_cycle ) call sld_cardiac_cycle_write_res(2_ip)
        if( kfl_preac_sld     ) call sld_post_reaction_write_res(2_ip)

#ifndef PROPER_ELEM_PRIVATE_OFF
        if(kfl_sysnet) call sysnet_write_alya_res(3_ip, lun_sysnet_heart_res_sld)
        if(kfl_sysnet) call sysnet_write_alya_res(4_ip, lun_sysnet_system_res_sld)
#endif

     end select

  end if

end subroutine sld_outinf

