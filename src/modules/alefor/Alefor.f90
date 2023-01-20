!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!





subroutine Alefor(order)
!-----------------------------------------------------------------------
!****f* Alefor/alefor
! NAME
!   Alefor
! DESCRIPTION
!   This routine deals with the ALE formulation equation. The task done
!   corresponds to the order given by the master.
! USES
!    ale_turnon
!    ale_begste
!    ale_doiter
!    ale_gencon
!    ale_endste
!    ale_turnof
! USED BY
!    Turnon
!    Begste
!    Doiter
!    Gencon
!    Endste
!    Turnof
!***
!-----------------------------------------------------------------------
  use      def_parame
  use      def_master

  implicit none
  
  integer(ip) :: order

  select case (order)

  case(ITASK_TURNON)
     call ale_turnon()
  case(ITASK_SOLMEM)
     call ale_solmem()
  case(ITASK_BEGRUN)
     call ale_begrun()
  case(ITASK_TIMSTE) 
     !call ale_timste()
  case(ITASK_INIUNK) 
     call ale_iniunk()
  case(ITASK_BEGSTE) 
     call ale_begste()
  case(ITASK_BEGZON)
     call ale_begzon()
  case(ITASK_DOITER)
     call ale_doiter()
  case(ITASK_CONCOU)
     call ale_concou()
  case(ITASK_CONBLK)
     !call ale_conblk()
  case(ITASK_ENDSTE)
     call ale_endste()
  case(ITASK_OUTPUT)
     call ale_output()
  case(ITASK_TURNOF)
     call ale_turnof()
  case( ITASK_READ_RESTART )
     call ale_restar(ITASK_READ_RESTART)
  case( ITASK_WRITE_RESTART )
     call ale_restar(ITASK_WRITE_RESTART)
  case( ITASK_REDIST )
     call ale_redist()
  case( ITASK_INTERP )
     call ale_interp()
  end select
  !
  ! Coupling
  !
  if( order > 1000_ip ) call ale_plugin(order-1000) ! Compute and send 
 
end subroutine Alefor
      
