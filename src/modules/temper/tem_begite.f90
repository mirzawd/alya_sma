!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!




subroutine tem_begite
  !-----------------------------------------------------------------------
  !****f* Temper/tem_begite
  ! NAME 
  !    tem_begite
  ! DESCRIPTION
  !    This routine starts an internal iteration for the temperature
  !    equation
  ! USES
  !    tem_tittim
  !    tem_updbcs
  !    tem_inisol
  !    tem_updunk
  ! USED BY
  !    tem_doiter
  !***
  !-----------------------------------------------------------------------
  use      def_parame
  use      def_master
  use      def_domain
  use      def_temper
  use      mod_messages,        only : livinf
  use      mod_tem_therm_press, only : therm_press_update
  implicit none
  !
  ! Initializations
  !  
  kfl_goite_tem = 1 
  itinn(modul)  = 0
  if(itcou==1) call tem_tistep()
  call livinf(15_ip,' ',modul)
  !
  ! Update boundary conditions
  !
  call tem_updbcs(ITASK_BEGITE)
  !
  ! Set up the solver parameters for the temperature equation
  !
  call tem_inisol()
  !
  ! Set up the parameters for the optimization
  !
  call tem_iniopt()
  !
  ! Obtain the initial guess for inner iterations
  !
  call tem_updunk(ITASK_BEGITE)
  !
  ! Coupling
  !
  call tem_coupli(ITASK_BEGITE)
  !
  ! Low-Mach initializes thermodynamic pressure
  !
  if( itcou == 1 .and. ittim == 1 ) then 
     call therm_press_update(-1_ip)
  end if
  ! ********************
  !   Test for coupling 
  ! call tem_temexch(1_ip)
  ! ********************

end subroutine tem_begite
    
