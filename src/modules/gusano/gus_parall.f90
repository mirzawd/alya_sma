!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Gusano
!> @{
!> @file    gus_parall.f90
!> @author  houzeaux
!> @date    2020-10-20
!> @brief   Parall
!> @details Broadcast data to slaves
!> @} 
!-----------------------------------------------------------------------

subroutine gus_parall()

  use def_parame
  use def_master
  use def_domain
  use def_solver
  use def_gusano
  use def_inpout
  use mod_memory
  use mod_opebcs
  use mod_output_postprocess, only : output_postprocess_parall
  use mod_solver,             only : solver_parall
  use mod_exchange,           only : exchange_add
  use mod_exchange,           only : exchange_end
  implicit none

  !----------------------------------------------------------------
  !
  ! Exchange data read in gus_reaphy
  !
  !----------------------------------------------------------------
  
  call exchange_add(kfl_timei_gus)
  call exchange_add(kfl_regim_gus)
  call exchange_add(kfl_bendi_gus)
  call exchange_add(kfl_conve_gus)
  
  !----------------------------------------------------------------
  !
  ! Exchange data read in gus_reanut
  !
  !----------------------------------------------------------------
  
  call exchange_add(momod(modul) % miinn)
  call exchange_add(kfl_stabi_gus       )
  call exchange_add(kfl_algor_gus       )
  call exchange_add(cotol_gus           )
  call exchange_add(safet_gus           )
  call exchange_add(sstol_gus           )
  call solver_parall()
  
  !----------------------------------------------------------------
  !
  ! Exchange data read in gus_reaous
  !
  !----------------------------------------------------------------
  
  call output_postprocess_parall()
  call exchange_end()

  !----------------------------------------------------------------
  !
  ! Boundary conditions
  !
  !----------------------------------------------------------------
  
  call boundary_conditions_exchange(tncod_gus)  
  call boundary_conditions_exchange(tbcod_gus)

end subroutine gus_parall
