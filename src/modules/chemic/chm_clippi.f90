!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Chemic
!> @{
!> @file    chm_clippi.f90
!> @author  houzeaux
!> @date    2018-12-28
!> @brief   Clipping
!> @details Clipping
!> @}
!-----------------------------------------------------------------------

subroutine chm_clippi()

!  use def_parame
  use def_master,      only     : unkno, therm
  use def_domain,      only     : npoin
  use def_kintyp,      only     : ip, rp
  use def_kermod, only          : lookup_fw
  use def_chemic,      only     : table_fw, kfl_model_chm, kfl_negat_chm, kfl_posit_chm, kfl_premix_chm, kfl_solve_cond_CMC_chm,&
                                  kfl_spray_chm, nclas_chm, Zs_CMC_chm, xZs_chm, xZr_chm, kfl_multimod_chm, kfl_tab_fw_chm_diff
  use mod_interp_tab,  only     : fw_scale_cont_var
  use mod_interp_tab,  only     : max_lookup_dim
  implicit none
  integer(ip)               :: ipoin,kpoin,idimt,pos
  real(rp)                  :: control(max_lookup_dim)   ! input of table lookup function
  real(rp)                  :: scale_control(max_lookup_dim)
  real(rp)                  :: lim_control(max_lookup_dim,2_ip)
  integer(ip)               :: ind(max_lookup_dim)

  select case(kfl_model_chm)

  case (1_ip,2_ip) ! At Begste

     !
     ! Prevent unde/over shoots
     !
     if ( kfl_negat_chm == 1 .or. kfl_posit_chm == 1) then

        !
        ! Unreacting spray
        !
        if (kfl_premix_chm == 1 .and. kfl_spray_chm /= 0_ip ) then

           kpoin = 0
           do ipoin = 1,npoin
              !
              ! Liquid volume fraction: phi_L
              !
              kpoin = (ipoin - 1) * nclas_chm + 3
              unkno(kpoin) = max(0.0_rp,min(1.0_rp,unkno(kpoin)))
           end do

        else

           kpoin = 0
           ind   = 1_ip
           do ipoin = 1,npoin
              kpoin = (ipoin - 1) * nclas_chm

              control = 0.0_rp
              if (kfl_multimod_chm == 1_ip) then
                 do idimt = 1, lookup_fw(kfl_tab_fw_chm_diff) % main_table % ndim
                    if (lookup_fw(kfl_tab_fw_chm_diff) % kfl_chm_control(idimt) > 0) then
                       !
                       ! >0: one of the conces
                       !
                       control(idimt) = unkno(kpoin+lookup_fw(kfl_tab_fw_chm_diff) % kfl_chm_control(idimt))
                    else
                       if (lookup_fw(kfl_tab_fw_chm_diff) % kfl_chm_control(idimt) == -1) then
                          !
                          ! -1: enthalpy
                          !
                          control(idimt) = therm(ipoin,1)
                       elseif (lookup_fw(kfl_tab_fw_chm_diff) % kfl_chm_control(idimt) == -2) then
                          !
                          ! -2: scalar dissipation rate
                          !
                          control(idimt) = xZr_chm(ipoin) + xZs_chm(ipoin)
                       endif
                    endif
                 enddo

                 call fw_scale_cont_var( control, scale_control, lim_control, lookup_fw(kfl_tab_fw_chm_diff), ind)

                 do idimt = 1, lookup_fw(kfl_tab_fw_chm_diff) % main_table % ndim
                    if (lookup_fw(kfl_tab_fw_chm_diff) % kfl_chm_control(idimt) > 0) then
                        unkno(kpoin+lookup_fw(kfl_tab_fw_chm_diff) % kfl_chm_control(idimt)) = max(lim_control(idimt,1),&
                        min(lim_control(idimt,2),unkno(kpoin+lookup_fw(kfl_tab_fw_chm_diff) % kfl_chm_control(idimt))))
                    endif
                 enddo
              else
                 do idimt = 1, table_fw % main_table % ndim
                    if (table_fw % kfl_chm_control(idimt) > 0) then
                       !
                       ! >0: one of the conces
                       !
                       control(idimt) = unkno(kpoin+table_fw % kfl_chm_control(idimt))
                    else
                       if (table_fw % kfl_chm_control(idimt) == -1) then
                          !
                          ! -1: enthalpy
                          !
                          control(idimt) = therm(ipoin,1)
                       elseif (table_fw % kfl_chm_control(idimt) == -2) then
                          !
                          ! -2: scalar dissipation rate
                          !
                          control(idimt) = xZr_chm(ipoin) + xZs_chm(ipoin)
                       endif
                    endif
                 enddo

                 call fw_scale_cont_var( control, scale_control, lim_control, table_fw, ind)

                 do idimt = 1, table_fw % main_table % ndim
                    if (table_fw % kfl_chm_control(idimt) > 0) then
                        unkno(kpoin+table_fw % kfl_chm_control(idimt)) = max(lim_control(idimt,1),min(lim_control(idimt,2),&
                            unkno(kpoin+table_fw % kfl_chm_control(idimt))))
                    endif
                 enddo
              end if 
           end do


           if (kfl_spray_chm /= 0_ip ) then
              kpoin = 0
              do ipoin = 1,npoin
                 !
                 ! Liquid volume fraction: phi_L
                 !
                 kpoin = (ipoin - 1) * nclas_chm + 5
                 unkno(kpoin) = max(0.0_rp,min(1.0_rp,unkno(kpoin)))
              end do
           end if

        end if

     end if


  case (4) ! Clipping when solving mixing variables with CMC

     if (kfl_solve_cond_CMC_chm == 0) then
        !
        ! Clip mixture fraction and its variance
        !
        do ipoin = 1, npoin
           pos = (ipoin-1_ip)*4_ip
           unkno(3_ip+pos) = max(min(unkno(3_ip+pos), Zs_CMC_chm), 0.0_rp)
           unkno(4_ip+pos) = max(min(unkno(4_ip+pos), &
                                 unkno(3_ip+pos)*(Zs_CMC_chm-unkno(3_ip+pos))), 0.0_rp)
           xZr_chm(ipoin)   = max(xZr_chm(ipoin), 0.0_rp)
           xZs_chm(ipoin)   = max(xZs_chm(ipoin), 0.0_rp)
        end do
     end if

  end select

end subroutine chm_clippi
