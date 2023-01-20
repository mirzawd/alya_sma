!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine Levels(order)
  !-----------------------------------------------------------------------
  !****f* Levels/Levels
  ! NAME
  !   Levels
  ! DESCRIPTION
  !   This routine deals with the level set advection equation. The task done
  !   corresponds to the order given by the master.
  !
  !   Created:     01 February 2007 
  !   by:          Anne-Cecile Lesage
  !   Responsible: Anne-Cecile Lesage  
  ! USES
  !    lev_turnon
  !    lev_timste
  !    lev_begste
  !    lev_doiter
  !    lev_concon
  !    lev_conblk
  !    lev_newmsh
  !    lev_endste
  !    lev_turnof
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
  use      def_master
  implicit none
  integer(ip), intent(in) :: order

  select case (order)

  case(ITASK_TURNON)
     call lev_turnon()
  case(ITASK_TIMSTE) 
     call lev_timste()
  case(ITASK_INIUNK) 
     call lev_iniunk()
  case(ITASK_BEGSTE) 
     call lev_begste()
  case(ITASK_DOITER)
     call lev_doiter()
  case(ITASK_CONCOU)
     call lev_concou()
  case(ITASK_CONBLK)
     call lev_conblk()
  case(ITASK_ENDSTE)
     call lev_endste()
  case(ITASK_OUTPUT)
     call lev_output()
  case(ITASK_TURNOF)
     call lev_turnof()

  end select

end subroutine Levels
