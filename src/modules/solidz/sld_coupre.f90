!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    sld_coupre.f90
!> @author  J.C. Cajas
!> @date    09/06/2016
!> @brief   This routine uses the temporal predictor for zonal coupling
!> @details This routine uses the temporal predictor for zonal coupling
!> @}
!-----------------------------------------------------------------------

subroutine sld_coupre()

  use def_kintyp, only :  ip, rp
  
  use def_coupli, only :  coupling_type
  use def_coupli, only :  mcoup
  use def_coupli, only :  UNKNOWN
  use def_coupli, only :  RESIDUAL

  use def_master, only :  itti2
  use def_master, only :  solve_sol
  use def_master, only :  modul
  use def_master, only :  current_code
  use def_master, only :  ID_NASTIN
  use def_master, only :  ID_SOLIDZ
  use def_master, only :  INOTMASTER
  use def_master, only :  gdepo

  use def_domain, only :  ndime

  use def_solidz, only :  kfl_fixno_sld
  use def_solidz, only :  kfl_gdepo
  use def_solidz, only :  bvess_sld

  implicit none 

  real(rp)             :: foref
  integer(ip)          :: icoup, ipoin, kpoin, idime, jdime
  !
  ! Use of the prediction in zonal coupling
  !
  if( mcoup > 0 .and. itti2 > 2_ip )then

     do icoup = 1_ip, mcoup

        if( current_code == coupling_type(icoup) % code_target .and. modul == coupling_type(icoup) % module_target .and. &
             coupling_type(icoup) % temporal_predictor == 1_ip )then

           !
           ! Predict Force 
           !
           if( coupling_type(icoup) % what == RESIDUAL )then
              !
              ! Predict Force from FSI
              !
              if( coupling_type(icoup) % module_source == ID_NASTIN )then
 
                 if( kfl_gdepo == 1_ip ) then
                    !
                    ! Push forward for the force coming from nastin
                    !
                    if( INOTMASTER ) then

                       do kpoin = 1, coupling_type(icoup) % wet % npoin_wet

                          ipoin = coupling_type(icoup) % wet % lpoin_wet(kpoin)

                          do idime = 1,ndime

                             foref = 0.0_rp

                             if( kfl_fixno_sld(idime,ipoin) /= 1 ) then

                                do jdime = 1,ndime

                                   foref = foref + gdepo(idime,jdime,ipoin) * coupling_type(icoup) % values_converged(jdime,kpoin,1_ip) 

                                end do

                                solve_sol(1) % bvnat(idime,ipoin) = foref

                             end if

                          end do

                       end do

                    end if  ! I NOT MASTER
                 else
                    if( INOTMASTER ) then
                       do kpoin = 1, coupling_type(icoup) % wet % npoin_wet
                          ipoin = coupling_type(icoup) % wet % lpoin_wet(kpoin)
                          do idime = 1,ndime
                             if( kfl_fixno_sld(idime,ipoin) /= 1 ) then
                                solve_sol(1) % bvnat(idime,ipoin) = coupling_type(icoup) % values_converged(idime,kpoin,1_ip)
                             end if
                          end do
                       end do
                    end if
                 end if

              else if(  coupling_type(icoup) % module_source == ID_SOLIDZ )then
                 ! 
                 ! Force coming from lagrangian framework (other solidz)
                 !
                 do kpoin = 1_ip, coupling_type(icoup) % wet % npoin_wet

                    ipoin = coupling_type(icoup) % wet % lpoin_wet(kpoin)
                    do idime = 1_ip, ndime
                       solve_sol(1) % bvnat(idime,ipoin) =  coupling_type(icoup) % values_converged(idime,kpoin,1_ip)
                    end do

                 end do
                 
              else 
                 !
                 ! Not recognized force source
                 !
                 call runend('SLD_COUPRE: Not recognized force source for prediction, check sld_begzon ')
                 
              end if     ! Forces

           else if( coupling_type(icoup) % what == UNKNOWN  )then
              ! 
              ! Predict displacements
              !
              if( INOTMASTER )then

                 do kpoin = 1_ip, coupling_type(icoup) % wet % npoin_wet
                    ipoin = coupling_type(icoup) % wet % lpoin_wet(kpoin)
                    do idime = 1_ip, ndime
                       bvess_sld(idime,ipoin,1_ip) = coupling_type(icoup) % values_converged(idime,kpoin,1_ip)
                    end do
                 end do

              end if

           end if        ! What (unknown, residual)

        end if           ! Code target

     end do              ! Coupling loop

  end if                 ! mcoup > 0

end subroutine sld_coupre
