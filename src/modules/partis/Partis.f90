!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine Partis(order)
  !-----------------------------------------------------------------------
  !****f* partis/Partis
  ! NAME
  !   Partis
  ! DESCRIPTION
  !   This routine deals with the ADS equation. The task done
  !   corresponds to the order given by the master.
  ! USES
  !    pts_turnon
  !    pts_timste
  !    pts_begste
  !    pts_doiter
  !    pts_concon
  !    pts_conblk
  !    pts_newmsh
  !    pts_endste
  !    pts_turnof
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
!  integer,dimension(8) :: values
!  real(rp)       :: dt1,dt2
  
  select case (order)
 
  case(ITASK_TURNON)
     call pts_turnon()
  case(ITASK_INIUNK) 
     call pts_iniunk()
  case(ITASK_BEGRUN) 
     call pts_begrun()
  case(ITASK_TIMSTE) 
     !call pts_timste()
  case(ITASK_BEGSTE) 
     call pts_begste()
  case(ITASK_DOITER)
     call pts_doiter()
  case(ITASK_SOLMEM)
     call pts_solmem()
  case(ITASK_CONCOU) 
     !call pts_concou()
  case(ITASK_CONBLK)
     !call pts_conblk()
  case(ITASK_ENDSTE)
     call pts_endste()
  case(ITASK_REDIST)
     call pts_redist()
  case( ITASK_READ_RESTART )
     call pts_restar(ITASK_READ_RESTART)
  case( ITASK_WRITE_RESTART )
     call pts_restar(ITASK_WRITE_RESTART)
  case(ITASK_OUTPUT)
!    call date_and_time(VALUES=values)
!    dt1=values(5)*3600+values(6)*60+values(7)+values(8)*0.001    ! seconds

#ifdef DBPARTICLES
     call pts_output_parall_db()
     call pts_output()
#else
     call pts_output()
#endif          
!     call date_and_time(VALUES=values)
!     dt2=values(5)*3600+values(6)*60+values(7)+values(8)*0.001    ! seconds

!     write(*,*) 'time_sec_output',dt2 - dt1
  case(ITASK_TURNOF)
     !call pts_turnof()
  end select
  !
  ! Coupling
  ! 
  if( order > 1000 ) call pts_plugin(order-1000_ip) 

end subroutine Partis
