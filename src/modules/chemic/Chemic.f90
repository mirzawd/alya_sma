!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine Chemic(order)
  !-----------------------------------------------------------------------
  !****f* chemic/Chemic
  ! NAME
  !   Chemic
  ! DESCRIPTION
  !   This routine deals with the ADS equation. The task done
  !   corresponds to the order given by the master.
  ! USES
  !    chm_turnon
  !    chm_timste
  !    chm_begste
  !    chm_doiter
  !    chm_concon
  !    chm_conblk
  !    chm_newmsh
  !    chm_endste
  !    chm_turnof
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
  use def_master,     only : ITASK_TURNON, ITASK_INIUNK, ITASK_TIMSTE, ITASK_BEGSTE, ITASK_DOITER, ITASK_CONCOU,&
                             ITASK_CONBLK, ITASK_ENDSTE, ITASK_OUTPUT, ITASK_TURNOF, ITASK_SOLMEM, ITASK_REDIST,&
                             ITASK_READ_RESTART, ITASK_WRITE_RESTART
  use def_kintyp,      only: ip

  implicit none
  integer(ip), intent(in) :: order

  external                :: chm_turnon
  external                :: chm_iniunk
  external                :: chm_timste
  external                :: chm_begste
  external                :: chm_doiter
  external                :: chm_concou
  external                :: chm_conblk
  external                :: chm_endste
  external                :: chm_output
  external                :: chm_turnof
  external                :: chm_solmem
  external                :: chm_redist
  external                :: chm_restar
  external                :: chm_plugin

  select case (order)

  case(ITASK_TURNON)
     call chm_turnon()
  case(ITASK_INIUNK)
     call chm_iniunk()
  case(ITASK_TIMSTE)
     call chm_timste()
  case(ITASK_BEGSTE)
     call chm_begste()
  case(ITASK_DOITER)
     call chm_doiter()
  case(ITASK_CONCOU)
     call chm_concou()
  case(ITASK_CONBLK)
     call chm_conblk()
  case(ITASK_ENDSTE)
     call chm_endste()
  case(ITASK_OUTPUT)
     call chm_output()
  case(ITASK_TURNOF)
     call chm_turnof()
  case(ITASK_SOLMEM)
     call chm_solmem()
  case(ITASK_REDIST)
     call chm_redist()
  case(ITASK_READ_RESTART)
     call chm_restar(ITASK_READ_RESTART)
  case(ITASK_WRITE_RESTART)
     call chm_restar(ITASK_WRITE_RESTART)

  end select
  !
  ! Coupling
  !
  if( order > 1000 ) call chm_plugin(order-1000_ip)

end subroutine Chemic
