!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_solite()
  !-----------------------------------------------------------------------
  !****f* Nastin/nsi_solite
  ! NAME 
  !    nsi_solite
  ! DESCRIPTION
  !    This routine solves an iteration for the incompressible
  !    Navier-Stokes equations, using:
  !    - A Monolithic scheme
  !    - A block Gauss-Seidel scheme
  ! USES
  !    nsi_solmon
  !    nsi_solbgs
  !    nsi_rotunk
  ! USED BY
  !    nsi_doiter
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_nastin
  use mod_nsi_schur_complement,  only : nsi_schur_complement_solution
  use mod_nsi_bubble,            only : nsi_bubble_update
  use mod_solver,                only : solver_postprocess
  use mod_nsi_fractional_step,   only : nsi_fractional_step_solution
  use mod_nsi_multi_step_fs,     only : nsi_multi_step_fs_solution
  use mod_nsi_semi_implicit,     only : nsi_semi_implicit_solution
  use mod_nsi_immersed_boundary, only : nsi_ib_lagrange_multiplier_update
  use def_kermod,                only : kfl_adj_prob
  use mod_messages,              only : livinf
  use mod_messages,              only : messages_live

  implicit none
  integer(ip) :: kfl_linea_old
#ifdef outmateocoe
  integer(ip)        :: ipass_aux
#endif

  !----------------------------------------------------------------------
  !
  ! Initialization
  !
  !----------------------------------------------------------------------
  !
  ! Update inner iteration counter and write headings in the solver file.
  !
  itinn(modul) = itinn(modul) + 1
  ittot_nsi    = ittot_nsi + 1
  !
  ! Linearization
  !
  kfl_linea_old = kfl_linea_nsi
  if( itinn(modul) <= npica_nsi .and. kfl_adj_prob == 0 ) kfl_linea_nsi = 1
  !
  ! Update boundary conditions
  !
  call nsi_updbcs(ITASK_INNITE)
  !
  ! Initialize SGS residual
  !
  call nsi_solsgs(1_ip)

  !----------------------------------------------------------------------
  !
  ! Solve equations 
  !
  !----------------------------------------------------------------------

  if( NSI_MONOLITHIC ) then
     !
     ! Monolithic
     !
     call nsi_inisol(1_ip)
     call nsi_matrix()
#ifdef outmateocoe
     call nsi_matndof(rhsid,unkno,amatr,r_sol,c_sol,ndime+1,npoin,ipass_aux)
#endif
#ifdef solve_w_agmg
     call nsi_agmgsol(rhsid,unkno,amatr,r_sol,c_sol,ndime+1,npoin)
#else
     call solver(rhsid,unkno,amatr,pmatr) 
     
#endif       

  else if( NSI_FRACTIONAL_STEP ) then
     !
     ! Fractional step
     !
     if(kfl_tisch_nsi == 4) then ! Multi step Fractional Step
        call nsi_multi_step_fs_solution()
     else                        ! One step Fractional Step
        call nsi_fractional_step_solution()
     end if
     
  else if( NSI_SEMI_IMPLICIT ) then
     !
     ! Semi-implicit
     !
     call nsi_semi_implicit_solution()
  
  else if( NSI_SCHUR_COMPLEMENT ) then
     !
     ! Schur complement
     !
     call nsi_schur_complement_solution()

  end if
  !
  ! Solver postprocess
  !  
  call solver_postprocess(momod(modul) % solve,amatr,rhsid,unkno)
  !
  ! Linearization
  !
  kfl_linea_nsi = kfl_linea_old
  
  !----------------------------------------------------------------------
  !
  ! Immersed boundary method
  !
  !----------------------------------------------------------------------

  call nsi_ib_lagrange_multiplier_update()

  !if( INOTMASTER ) then
  !   call PAR_INTERFACE_NODE_EXCHANGE(ndime,unkno,'TAKE MIN')
  !   call PAR_INTERFACE_NODE_EXCHANGE(1_ip,unkno(ndbgs_nsi+1:),'TAKE MIN')
  !end if

end subroutine nsi_solite
