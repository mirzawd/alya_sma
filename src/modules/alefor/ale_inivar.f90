!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine ale_inivar(itask)
  !-----------------------------------------------------------------------
  !****f* Alefor/ale_inivar
  ! NAME 
  !    ale_inicar
  ! DESCRIPTION
  !    This routine initializes some variables
  !    ITASK=1 ... When starting the run (from Turnon)
  !    ITASK=2 ... First time step. This is needed as some variables 
  !                are not initialized before
  ! USES
  ! USED BY
  !    ale_turnon
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_alefor
  use def_solver
  use mod_iofile
  use mod_arrays, only : arrays_register
  use mod_alefor, only : alefor_initialization
  implicit none
  integer(ip), intent(in) :: itask

  select case ( itask )

  case ( 0_ip )
     !
     ! Nullify pointers
     !
     call alefor_initialization()
     !
     ! variable registration
     !
     call arrays_register((/'DISPM','VECTO','NPOIN','PRIMA'/),dispm    ,ENTITY_POSITION=2_ip,TIME_POSITION=3_ip)
     call arrays_register((/'VELOM','VECTO','NPOIN','PRIMA'/),velom    ,ENTITY_POSITION=2_ip)
     call arrays_register((/'COALE','VECTO','NPOIN','PRIMA'/),velom    ,ENTITY_POSITION=2_ip,TIME_POSITION=3_ip)
     call arrays_register((/'COORI','VECTO','NPOIN','PRIMA'/),coord_ori,ENTITY_POSITION=2_ip)
     call arrays_register((/'GROUP','SCALA','NPOIN','SECON'/))
     call arrays_register((/'BVESS','VECTO','NPOIN','SECON'/))
     !
     ! Solver
     ! Displacement is solved separately
     !     
     call soldef(-1_ip)
     solve(1) % wprob     = 'DISPLACEMENT'    ! Equation name
     solve(1) % kfl_solve = 1                 ! Output flag
     solve(1) % ndofn     = ndime             ! Dimension

  end select

end subroutine ale_inivar
