!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Alefor
!> @{
!> @file    ale_coupre.f90
!> @author  J.C. Cajas
!> @date    09/06/2016
!> @brief   This routine uses the temporal predictor for zonal coupling
!> @details This routine uses the temporal predictor for zonal coupling
!> @} 
subroutine ale_coupre()
  !-----------------------------------------------------------------------
  !****f* Solidz/sld_coupre
  ! NAME 
  !    sld_coupre
  ! DESCRIPTION
  !    This routine uses the temporal predictor for zonal coupling
  ! USES
  ! USED BY
  !    sld_begste
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only :  ip, rp

  use def_coupli, only :  coupling_type
  use def_coupli, only :  mcoup
  use def_coupli, only :  UNKNOWN
  use def_coupli, only :  RESIDUAL

  use def_master, only :  itti2
  use def_master, only :  modul
  use def_master, only :  current_code
  use def_master, only :  ID_NASTIN
  use def_master, only :  ID_SOLIDZ
  use def_master, only :  INOTMASTER
  use def_master, only :  bvess_ale

  use def_domain, only :  ndime


  implicit none 

  integer(ip)          :: icoup, ipoin, kpoin, idime
  !
  ! Use of the prediction in zonal coupling
  !
  if( mcoup > 0 .and. itti2 > 2_ip )then

     do icoup = 1_ip, mcoup

        if( current_code == coupling_type(icoup) % code_target .and. modul == coupling_type(icoup) % module_target &
             .and. coupling_type(icoup) % temporal_predictor == 1_ip )then
           
           !
           ! Predict Force 
           !
           if( coupling_type(icoup) % what == RESIDUAL )then
              !
              ! Predict Force
              !
              call runend(' Seriously residual for ALE? Are you kidding me? ale_coupre')

           else if( coupling_type(icoup) % what == UNKNOWN  )then
              ! 
              ! Predict displacements
              !
              if( INOTMASTER )then

                 do kpoin = 1_ip, coupling_type(icoup) % wet % npoin_wet
                    ipoin = coupling_type(icoup) % wet % lpoin_wet(kpoin)
                    do idime=1,ndime
                       bvess_ale(idime,ipoin,1) = coupling_type(icoup) % values_converged(idime,kpoin,1_ip)
                    enddo
                 end do

              end if

           end if        ! What (unknown, residual)

        end if           ! Code target

     end do              ! Coupling loop

  end if                 ! mcoup > 0

end subroutine ale_coupre
