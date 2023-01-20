!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine ale_parall(order)
  !-----------------------------------------------------------------------
  !****f* Alefor/ale_parall
  ! NAME
  !    ale_parall
  ! DESCRIPTION
  !    This routine exchange data 
  ! USES
  ! USED BY
  !    ale_turnon
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_solver
  use def_alefor
  use def_inpout
  use mod_memchk
  use mod_opebcs
  use def_kermod,             only : number_space_time_function
  use mod_alefor,             only : alefor_memory_allocate
  use mod_exchange,           only : exchange_init
  use mod_exchange,           only : exchange_add
  use mod_exchange,           only : exchange_end
  use mod_solver,             only : solver_parall
  use mod_output_postprocess, only : output_postprocess_parall
  implicit none
  integer(ip), intent(in) :: order
  integer(ip)             :: ji
  integer(4)              :: istat

  if( ISEQUEN ) return

  select case (order)

  case(1_ip)    

     !------------------------------------------------------------------- 
     !
     ! Exchange data read in ale_reaphy, ale_reanut and ale_reaous
     !
     !------------------------------------------------------------------- 

     call exchange_init()
     call exchange_add(nrbod)
     call exchange_add(kfl_rigid_ale)
     call exchange_end()
     if( ISLAVE .and. kfl_rigid_ale == 1 ) call alefor_memory_allocate('RBBOU')

     call exchange_init()
     !
     ! Exchange of ale_reaphy, ale_reanut variables 
     !
     call exchange_add(kfl_smoot_ale)
     call exchange_add(kfl_timef_ale)
     call exchange_add(nsmoo_ale)
     call exchange_add(kfl_defor_ale)
     call exchange_add(ndefo_ale)
     call exchange_add(kfl_smobo_ale)
     call exchange_add(kfl_fixsm_ale)
     call exchange_add(nsmob_ale)
     call exchange_add(kfl_crist_ale)
     call exchange_add(kfl_foexo_ale)
     call exchange_add(kfl_disor_ale)
     call exchange_add(kfl_nforc_ale)
     call exchange_add(kfl_sensi_ale)
     call exchange_add(moddi_ale)
     call exchange_add(modvi_ale)
     call exchange_add(ansmo_ale)
     call exchange_add(resmo_ale)
     call exchange_add(defor_param_ale)
     !
     ! Variables for rigid body
     !
     if( kfl_rigid_ale == 1 ) call ale_sendat(1_ip)
     !
     ! Variables computed in ibm_readim and ibm_cderda
     !
     do ji=1,nelty
        call exchange_add(lexib(ji))
     end do
     do ji=1,nelty
        call exchange_add(ngaib(Ji))
     end do
     do ji=1,nelty
        call exchange_add(lruib(ji))
     end do

     call exchange_add(mnoib)
     call exchange_add(mnodi)
     call exchange_add(mgaib)
     !
     ! Solver and postprocess
     !
     call solver_parall()
     call output_postprocess_parall()

     call exchange_end()
     !
     ! Boundary conditions
     !
     call boundary_conditions_exchange(tncod_ale)
     call boundary_conditions_exchange(tgcod_ale)
     call boundary_conditions_exchange(tbcod_ale)

  end select

end subroutine ale_parall
