!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine Temper(order)
  !-----------------------------------------------------------------------
  !****f* Temper/temper
  ! NAME
  !   Temper
  ! DESCRIPTION
  !   This routine deals with the temperature equation. The task done
  !   corresponds to the order given by the master.
  ! USES
  !    tem_turnon
  !    tem_timste
  !    tem_begste
  !    tem_doiter
  !    tem_concon
  !    tem_conblk
  !    tem_newmsh
  !    tem_endste
  !    tem_turnof
  ! USED BY
  !    Reapro
  !    Turnon
  !    Timste
  !    Begste
  !    Doiter
  !    Concon
  !    Conblk
  !    Newmsh
  !    Endste
  !    Turnof
  !***
  !-----------------------------------------------------------------------
  use def_master

  implicit none
  
  integer(ip), intent(in) :: order

  !
  ! Normal order
  !
  select case (order)

  case( ITASK_TURNON )
     call tem_turnon()
  case( ITASK_TIMSTE ) 
     call tem_timste()
  case( ITASK_INIUNK )
     call tem_iniunk()
  case( ITASK_BEGSTE )
     call tem_begste()
  case( ITASK_DOITER )
     call tem_doiter()
  case( ITASK_CONCOU )
     call tem_concou()
  case( ITASK_CONBLK )
     call tem_conblk()
  case( ITASK_ENDSTE )
     call tem_endste()
  case( ITASK_OUTPUT )
     call tem_output()
  case( ITASK_TURNOF )
     call tem_turnof()
  case( ITASK_SOLMEM )
     call tem_solmem()
  case( ITASK_BEGRUN )
     call tem_begrun()
  case( ITASK_REDIST )
     call tem_redist()
  case( ITASK_INTERP )
     call tem_interp()
  case( ITASK_READ_RESTART )
     call tem_restar(ITASK_READ_RESTART)
  case( ITASK_WRITE_RESTART )
     call tem_restar(ITASK_WRITE_RESTART)
  end select
  !
  ! Coupling
  !
  if( order > 1000 ) call tem_plugin(order-1000_ip) 

end subroutine Temper
