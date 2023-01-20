!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tur_begite()
  !-----------------------------------------------------------------------
  !****f* Turbul/tur_begite
  ! NAME 
  !    tur_begite
  ! DESCRIPTION
  !    This routine starts an internal iteration for the turbulence
  !    equations
  ! USES
  !    tur_tittim
  !    tur_frivel
  !    tur_updbcs
  !    tur_inisol
  !    tur_updunk
  ! USED BY
  !    tur_doiter
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_turbul
  use mod_messages, only : livinf
  implicit none
  !
  ! Initializations
  !
  kfl_goite_tur = 1
  itinn(modul)  = 0
  if(itcou==1) call tur_tistep()
  call livinf(15_ip,' ',modul)
  !
  ! Adaptive b.c.
  !
  call tur_adapti()
  !
  ! Compute friction velocity on walls
  !
  call tur_frivel()
  !
  ! Compute second order velocity gradients
  !
  call tur_grave2()
  !
  ! Production term
  !
  call tur_produc()
  !
  ! Average velocity
  !
  call tur_bouave(1_ip)
  !
  ! Compute grad(sqrt(k))
  !
  call tur_grsqki()
  !
  ! Compute magnitude of vorticity or strain rate
  !
  call tur_vortic()
  !
  ! Update boundary conditions
  !
  call tur_updbcs(TUR_BEFORE_GLOBAL_ITERATION)
  !
  ! Set up the parameters for the optimization
  !
  call tur_iniopt()
  !
  ! Obtain the initial guess for inner iterations
  !
  call tur_updunk(ITASK_BEGITE)

end subroutine tur_begite
