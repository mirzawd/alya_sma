!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine lev_sendat(order)
  !-----------------------------------------------------------------------
  !****f* Levels/lev_sendat
  ! NAME
  !    lev_sendat
  ! DESCRIPTION
  !    This routine exchange data 
  ! USES
  ! USED BY
  !    lev_turnon
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master 
  use def_kermod
  use def_solver
  use def_domain
  use def_levels
  use def_inpout
  use mod_memchk
  use mod_opebcs
  use mod_exchange,           only : exchange_init
  use mod_exchange,           only : exchange_add
  use mod_exchange,           only : exchange_end
  use mod_solver,             only : solver_parall
  use mod_output_postprocess, only : output_postprocess_parall
  implicit none
  integer(ip), intent(in) :: order
  integer(ip)             :: ji,ki
  integer(4)              :: istat

  select case (order)

  case(1_ip)     
     !
     ! Exchange data read in lev_reaphy, lev_reanut and lev_reaous
     !
     call exchange_init()
     !
     ! Exchange of lev_reaphy variables 
     !
     call exchange_add(kfl_inlev_lev)
     call exchange_add(kfl_advec_lev)
     call exchange_add(kfl_reave_lev)
     call exchange_add(nmate_lev)
     call exchange_add(thicl)                    ! belongs to defmaster. Shared with other modules, but only after itask_turnon.
     !
     ! Exchange of lev_reanut variables 
     !        
     call exchange_add(kfl_timet_lev)
     call exchange_add(kfl_tisch_lev)
     call exchange_add(kfl_timco_lev)
     call exchange_add(kfl_tiacc_lev)
     call exchange_add(kfl_normc_lev)
     call exchange_add(kfl_ellen_lev)
     call exchange_add(kfl_zonal_lev)
     call exchange_add(neule_lev)
     call exchange_add(miinn_lev)
     call exchange_add(inred_lev)
     call exchange_add(nfred_lev)
     call exchange_add(tyred_lev)
     call exchange_add(nstre_lev)
     call exchange_add(kfl_locre_lev)
     call exchange_add(nbitr_lev)
     call exchange_add(kfl_corvo_lev)
     call exchange_add(safet_lev) 
     call exchange_add(sstol_lev)
     call exchange_add(cotol_lev)
     call exchange_add(cpuit_lev)
     call exchange_add(supgp_lev)
     call exchange_add(npp_gauge_lev)           
     call exchange_add(npp_nbgau_lev)
     call exchange_add(npp_inter_lev)
     do ji=1,ngaug_lev
        call exchange_add(typga_lev(ji)) 
     end do
     do ki=1,ngaug_lev
        do ji=1,3
           call exchange_add(cogau_lev(ji,ki))
        end do
     end do
     !
     ! Exchange of lev_reabcs variables 
     !
     call exchange_add(kfl_inlev_lev) 
     call exchange_add(kfl_conbc_lev) 
     call exchange_add(height_lev)

     !
     ! Solver and postprocess
     !
     call solver_parall()
     call output_postprocess_parall()

     call exchange_end()

     !------------------------------------------------------------------- 
     !
     ! Variables read in reabcs
     !
     !------------------------------------------------------------------- 

     call boundary_conditions_exchange(tncod_lev)
     call boundary_conditions_exchange(tgcod_lev)

  case (6)


  end select

end subroutine lev_sendat
