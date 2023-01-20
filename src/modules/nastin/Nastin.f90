!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @defgroup Nastin
!> Incompressible - Boussinesq - Low Mach Navier-Stokes equations
!> @{
!> @file    Nastin.f90
!> @date    10/10/1972
!> @author  Guillaume Houzeaux
!> @brief   Incompressible NSI main
!> @details Nastin: incompressible - Low Mach Navier-Stokes equations. Main subroutine
!> @}
!------------------------------------------------------------------------
subroutine Nastin(order)

  use      def_master
  use      def_nastin

  implicit none
  
  integer(ip), intent(in) :: order

  select case ( order )

  case( ITASK_TURNON        ) ; call nsi_turnon()
  case( ITASK_BEGRUN        ) ; call nsi_begrun()
  case( ITASK_SOLMEM        ) ; call nsi_solmem()
  case( ITASK_TIMSTE        ) ; call nsi_timste()
  case( ITASK_INIUNK        ) ; call nsi_iniunk()
  case( ITASK_BEGSTE        ) ; call nsi_begste()
  case( ITASK_DOITER        ) ; call nsi_doiter()
  case( ITASK_CONCOU        ) ; call nsi_concou()
  case( ITASK_CONBLK        ) ; call nsi_conblk()
  case( ITASK_ENDSTE        ) ; call nsi_endste()
  case( ITASK_OUTPUT        ) ; call nsi_output()
  case( ITASK_TURNOF        ) ; call nsi_turnof()
  case( ITASK_REDIST        ) ; call nsi_redist()
  case( ITASK_INTERP        ) ; call nsi_interp()
  case( ITASK_READ_RESTART  ) ; call nsi_restar(ITASK_READ_RESTART)
  case( ITASK_WRITE_RESTART ) ; call nsi_restar(ITASK_WRITE_RESTART)
  case( 1001:               ) ; call nsi_plugin(order-1000_ip) 
  end select
  
end subroutine Nastin

