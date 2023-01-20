!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    sld_solite.f90
!> @author  Mariano Vazquez
!> @date    16/11/1966
!> @brief   This routine solves an iteration of the module equations.
!> @details This routine solves an iteration of the module equations.
!> @}
!-----------------------------------------------------------------------

subroutine sld_solite()

  use def_kintyp,               only : ip, rp
  use def_master,               only : itinn, modul
  use mod_messages,             only : livinf
  use def_solidz,               only : cpu_ass_sol_sld
  use def_solidz,               only : kfl_timet_sld, kfl_tisch_sld
  use def_solidz,               only : kfl_rigid_sld
  use def_solidz,               only : SLD_EXPLICIT_SCHEME, SLD_IMPLICIT_SCHEME
  use def_solidz,               only : kfl_timet_sld
  use mod_sld_solution_methods, only : sld_implicit, sld_explicit, sld_explicit_tw
  use mod_sld_rbo,              only : sld_rbo_solrbo

  implicit none

  real(rp)    :: time1,time2

  !----------------------------------------------------------------------
  !
  ! Initialization
  !
  !----------------------------------------------------------------------
  !
  ! Update inner iteration counter
  !
  itinn(modul) = itinn(modul) + 1_ip
  !
  ! Init cpu time
  !
  call cputim(time1)

  !----------------------------------------------------------------------
  !
  ! Solve equations
  !
  !----------------------------------------------------------------------

  if( kfl_rigid_sld == 0_ip ) then
     !
     ! Deformable body
     !
     if( kfl_timet_sld == SLD_EXPLICIT_SCHEME ) then

        call livinf(165_ip,'  EXPLICITLY... ',0_ip)

        if(      kfl_tisch_sld == 1_ip ) then
           call sld_explicit()    ! Central differences scheme
        else if ( kfl_tisch_sld == 2_ip ) then
           call sld_explicit_tw() ! TW scheme
        end if

     else if( kfl_timet_sld == SLD_IMPLICIT_SCHEME ) then

        call livinf(165_ip,'  IMPLICITLY... ',0_ip)

        if( kfl_tisch_sld == 1_ip ) then
           call sld_implicit()    ! Beta-Newmark implicit scheme
        end if

     end if

  else
     !
     ! Rigid body
     !
     if( kfl_tisch_sld == 3_ip ) then

        ! Classical Runge-Kutta
        call livinf(165_ip,'  RUNGE-KUTTA... ',0_ip)
        call sld_rbo_solrbo()

     end if

  end if
  !
  ! End cpu time
  !
  call cputim(time2)
  cpu_ass_sol_sld(2) = time2 - time1

end subroutine sld_solite
