!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!




!-----------------------------------------------------------------------
!> @defgroup Exmedi
!> @{
!> @file    Exmedi.f90
!> @author  Mariani Vazques
!> @date    2018-12-30
!> @brief   Main Exmedi surboutine
!> @details This is the main subroutine for the Excitable Media module.
!>          A Bidomain Model is used for cardiac electrical propagation.
!>          Two different kinds of models are implemented for ionic currents:
!>          - No subcell models (approximate models): FitzHugh-Nagumo (FHN)
!>          - Subcell models: TenTusscher (TT), Luo-Rudy (LR), ...
!> @} 
!-----------------------------------------------------------------------

subroutine Exmedi(order)

  use      def_parame
  use      def_domain
  use      def_master
  use      def_solver
  use      def_exmedi

  implicit none
  integer(ip) :: order
  
  select case (order)
     
  case(ITASK_TURNON)
     call exm_turnon
  case(ITASK_BEGRUN)
     call exm_begrun
  case( ITASK_SOLMEM)
     call exm_solmem
  case(ITASK_INIUNK)
     call exm_iniunk
  case(ITASK_OUTPUT)
     call exm_output
  case(ITASK_TIMSTE) 
     call exm_timste
  case(ITASK_BEGSTE) 
     call exm_begste
  case(ITASK_DOITER)
     call exm_doiter
  case(ITASK_CONCOU)
     call exm_concou
  case(ITASK_CONBLK)
     return
!!     call exm_conblk   --> to be done when needed
  case(ITASK_ENDSTE)
     call exm_endste
  case(ITASK_REDIST)
     call exm_redist
  case(ITASK_TURNOF)
     call exm_turnof
  case( ITASK_READ_RESTART )
     call exm_restar(ITASK_READ_RESTART)
  case( ITASK_WRITE_RESTART )
     call exm_restar(ITASK_WRITE_RESTART)
     
  end select
  !
  ! Coupling
  ! 
  if( order > 1000 ) call exm_plugin(order-1000_ip) 
  
end subroutine Exmedi
      
