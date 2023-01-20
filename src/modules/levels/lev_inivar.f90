!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Domain
!> @{
!> @file    domvar.f90
!> @author  Guillaume Houzeaux
!> @date    01/10/2015
!> @brief   Define some variables
!> @details Initialize and define some domain variables
!>          ITASK = 1 ... Before reading the module's data
!>                = 2 ... After reading the data as some variables
!>                        may depend on them
!> @} 
!-----------------------------------------------------------------------
subroutine lev_inivar(itask)
  use def_master
  use def_solver
  use def_levels
  use def_domain
  implicit none
  integer(ip), intent(in) :: itask 

  select case ( itask )
  
  case ( 1_ip )
     !
     ! Nullfiy modules' pointers
     !
     nullify(flsex_lev)
     nullify(tncod_lev)
     nullify(tgcod_lev)
     nullify(kfl_fixno_lev)
     nullify(kfl_funno_lev)
     nullify(bvess_lev)
     nullify(lnodb_lev)
     nullify(lnodp_lev)
     nullify(lebsu_lev)
     nullify(lebsp_lev)
     nullify(psbel_lev)
     nullify(lsbel_lev)
     nullify(nredm_lev)
     nullify(capin_lev)
     nullify(coord_lev)
     nullify(coorp_lev)
     nullify(dista_lev)
     nullify(flev0_lev)
     nullify(norml_lev)
     nullify(normv_lev)
     nullify(icupt_lev)
     nullify(elcro_lev)
     nullify(walld_lev)
     !
     ! Postprocess Variable names and types
     !
     postp(1) % wopos (1,1) = 'LEVEL'
     postp(1) % wopos (2,1) = 'SCALA'

     postp(1) % wopos (1,2) = 'VELOC'
     postp(1) % wopos (2,2) = 'VECTO'

     postp(1) % wopos (1,3) = 'GRADL'
     postp(1) % wopos (2,3) = 'VECTO'

     postp(1) % wopos (1,4) = 'DISTA'
     postp(1) % wopos (2,4) = 'SCALA'
     
     postp(1) % wopos (1,5) = 'DISPM'
     postp(1) % wopos (2,5) = 'VECTO'

     postp(1) % wopos (1,6) = 'VELOM'
     postp(1) % wopos (2,6) = 'VECTO'
     !
     ! Solver
     !
     call soldef(-5_ip)
     solve(1) % wprob     = 'LEVEL_SET'
     solve(1) % kfl_solve = 1         

     solve(2) % wprob     = 'REDISTANTIATION'
     solve(2) % kfl_solve = 1      

     solve(3) % wprob     = 'POISSON'
     solve(3) % kfl_solve = 1     

     solve(4) % wprob     = 'MESHD'
     solve(4) % kfl_solve = 1   
     solve(4) % ndofn     = ndime
     
     solve(5) % wprob     = 'GENERALIZED_DISTANCE'
     solve(5) % kfl_solve = 1   
     solve(5) % ndofn     = 1
     !
     ! Others
     !
     cpuit_lev      = 0.0_rp                    ! CPU time per iteration
     kfl_alloc_lev  = 0                         ! Interface not allocated

  case ( 2_ip )
     !
     ! Problem is linear
     !
     kfl_normc_lev = 2                          ! L2 norm for convergence
     kfl_tiaor_lev = kfl_tiacc_lev              ! Time accuracy: save original value
!!     solve(2) = solve(1)                      ! Now sussman redist uses its own solver
     solve(2) % wprob     = 'REDISTANTIATION'
     solve(2) % kfl_solve = 1                   ! Output solver info
     solve(2) % lun_solve = 1451

  end select
     
end subroutine lev_inivar
