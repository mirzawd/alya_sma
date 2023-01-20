!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Gusano
!> @{
!> @file    Gusano.f90
!> @author  houzeaux
!> @date    2020-10-19
!> @brief   Incompressible 1D model
!> @details Incompressible flow using a network model
!> @}
!------------------------------------------------------------------------
subroutine Gusano(order)

  use def_master
  use def_gusano

  implicit none
  integer(ip), intent(in) :: order

  select case ( order )

  case( ITASK_TURNON        ) ; call gus_turnon()
  case( ITASK_SOLMEM        ) ; call gus_solmem()
  case( ITASK_INIUNK        ) ; call gus_iniunk()
  case( ITASK_BEGSTE        ) ; call gus_begste()
  case( ITASK_DOITER        ) ; call gus_doiter()
  case( ITASK_ENDSTE        ) ; call gus_endste()
  case( ITASK_OUTPUT        ) ; call gus_output()
  case( ITASK_TIMSTE        ) ; call gus_timste()
  case( ITASK_BEGRUN        ) ; call gus_begrun()
!!$  case( ITASK_CONCOU        ) ; call gus_concou()
!!$  case( ITASK_CONBLK        ) ; call gus_conblk()
!!$  case( ITASK_TURNOF        ) ; call gus_turnof()
!!$  case( ITASK_REDIST        ) ; call gus_redist()
!!$  case( ITASK_INTERP        ) ; call gus_interp()
!!$  case( ITASK_READ_RESTART  ) ; call gus_restar(ITASK_READ_RESTART)
!!$  case( ITASK_WRITE_RESTART ) ; call gus_restar(ITASK_WRITE_RESTART)
  end select 
  !
  ! Coupling
  ! 
  if( order > 1000 ) call gus_plugin(order-1000_ip) 

end subroutine Gusano

