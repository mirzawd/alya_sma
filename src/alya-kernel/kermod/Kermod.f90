!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine Kermod(order)
  !-----------------------------------------------------------------------
  !****f* kermod/Kermod
  ! NAME 
  !    Kermod
  ! DESCRIPTION
  !    This routine deals with the incompressible NS equations.
  !    Kermod is monitorized for Paraver. 
  ! USES
  !    ker_turnon
  !    ker_timste
  !    ker_begste
  !    ker_doiter
  !    ker_concon
  !    ker_conblk
  !    ker_newmsh
  !    ker_endste
  !    ker_turnof
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
  use mod_messages, only : messages_live
  implicit none
  integer(ip), intent(in) :: order

  modul = mmodu
  call moddef(9_ip)

  select case ( order )

  case( -ITASK_TURNON )
     call messages_live(trim(namod(modul))//': READ DATA')
     call ker_turnon(1_ip)

  case( ITASK_BEGRUN )     
     call ker_begrun()
     
  case( ITASK_TURNON )
     call ker_turnon(2_ip)

  case( ITASK_TIMSTE ) 
     call ker_timste()

  case( -ITASK_INIUNK )
     call messages_live(trim(namod(modul))//': INITIAL SOLUTION')
     call ker_iniunk()

  case( ITASK_BEGSTE ) 
     call ker_begste()

  case( ITASK_DOITER )
     call ker_doiter()

  case(  ITASK_CONCOU )
     !call ker_concou()

  case(  ITASK_CONBLK )
     !call ker_conblk()

  case(  ITASK_ENDSTE )
     call ker_endste()

  case( -ITASK_OUTPUT )
     call ker_output()

  case( ITASK_TURNOF )
     !call ker_turnof()

  case( ITASK_SOLMEM )
     call ker_solmem()
     
  case( ITASK_REDIST )
     call ker_redist()
     
  case( ITASK_INTERP )
     call ker_interp()
     
  case( ITASK_READ_RESTART )
     call ker_restar(ITASK_READ_RESTART)
     
  case( ITASK_WRITE_RESTART )
     call ker_restar(ITASK_WRITE_RESTART)

  end select
  !
  ! Coupling
  ! 
  if( order > 1000 ) call ker_plugin(order-1000_ip) 

  modul = 0
  call moddef(9_ip)

end subroutine Kermod
