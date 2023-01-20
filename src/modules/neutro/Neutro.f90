!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup Neutro 
!> @{
!> @file    Neutro.f90
!> @date    29/03/2016
!> @author  Guillaume Houzeaux
!> @brief   Generic module
!> @details Generic module
!> @}
!------------------------------------------------------------------------

subroutine Neutro(order)

  use def_master
  implicit none
  integer(ip), intent(in) :: order
!   integer(ip)             :: icoup
  external :: neu_turnon, neu_iniunk, neu_begste, neu_doiter, neu_concou, neu_output

  select case ( order )

  case( ITASK_TURNON )
     call neu_turnon()
  case( ITASK_TIMSTE ) 
    ! call neu_timste()
  case( ITASK_INIUNK ) 
     call neu_iniunk() ! initial solution, initial valuo of neutr, if there are restartt too
  case( ITASK_BEGSTE )
     call neu_begste() ! initial soluton in temporal step (this is redundant with iniunk. To me is worg is there are restart)
  case( ITASK_DOITER )
     call neu_doiter() ! make the work In stationay case this is all!!
  case( ITASK_CONCOU )
     call neu_concou() ! Check convergence and otuput residual. Its nothing in case of not time funtionallity
  case( ITASK_CONBLK )
     !call neu_conblk() 
  case( ITASK_ENDSTE )
    ! call neu_endste()
  case( ITASK_OUTPUT )
     call neu_output() ! Initial solution, end of a time step and end of run. 
  case( ITASK_TURNOF )
     !call neu_turnof()
  end select
  !
  ! Coupling
  ! 
  !if( order > 1000 ) call neu_plugin(order-1000_ip) 
  
end subroutine Neutro

