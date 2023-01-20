!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    mod_sld_solution_methods.f90
!> @author  Matias Rivero (matias.rivero@bsc.es)
!>          Adria Quintanas (adria.quintanas@udg.edu)
!>          Gerard Guillamet (gerard.guillamet@bsc.es)
!> @date    January, 2017
!> @brief   Solidz time integration schemes in order to compute a non-linear
!>          iteration.
!> @details Time integration schemes\n
!>
!>          - Explicit Centred Diferences scheme
!>          - Beta-Newmark Implicit scheme
!>          - Dissipative Tchamwa-Wielgosz Explicit scheme
!>
!>          References:\n
!>          T. Belytschko, W. K. Liu, B. Moran, K. I. Elkhodary
!>          Nonlinear Finite elements for Continua and Structures\n
!>          E. J. Barbero. Introduction to Composite Materials Design\n
!>
!>          L. Mahéo, V. Grolleau, G. Rio. Damping efficiency of the\n
!>          Tchamwa-Wielgosz explicit dissipative scheme under instantaneous\n
!>          loading conditions\n
!>
!>
!> @}
!------------------------------------------------------------------------

module mod_sld_solution_methods

  use def_kintyp,           only : ip, rp
  use def_master,           only : INOTMASTER, ITER_K, TIME_N
  use def_master,           only : INOTEMPTY
  use def_master,           only : displ, unkno
  use def_master,           only : momod, solve_sol, solve
  use def_master,           only : itinn, modul, dtime, mem_modul
  use def_master,           only : rhsid, amatr, pmatr
  use def_domain,           only : ndime, npoin
  use mod_solver,           only : solver_solve, solver_lumped_mass_system
  use mod_array_operations, only : array_operations_initialization
  use mod_memory,           only : memory_alloca
  use mod_memory,           only : memory_deallo
  use mod_local_basis,      only : local_basis_global_to_local
  use mod_local_basis,      only : local_basis_local_to_global
  use def_solidz,           only : ndofn_sld
  use def_solidz,           only : kfl_timei_sld
  use def_solidz,           only : kfl_fixno_sld, kfl_fixrs_sld
  use def_solidz,           only : kfl_local_sld, kfl_conta_sld
  use def_solidz,           only : kfl_bvesv_sld
  use def_solidz,           only : SLD_STATIC_PROBLEM, SLD_DYNAMIC_PROBLEM
  use def_solidz,           only : SLD_IMPLICIT_SCHEME, SLD_EXPLICIT_SCHEME
  use def_solidz,           only : dunkn_sld, veloc_sld, accel_sld, ddisp_sld
  use def_solidz,           only : unknotmp_sld, veloctmp_sld, bvess_sld
  use def_solidz,           only : tifac_sld, vmass_sld, last_iters_sld
  use def_solidz,           only : fintt_sld, finte_sld, fexte_sld, fextt_sld
  use mod_sld_energy,       only : sld_energy
  use mod_sld_bclocal,      only : sld_bclocal_displ, sld_bclocal_rhsid_fixity, sld_bclocal_roback
  use mod_sld_fe2,          only : fe2_set_strains, fe2_homogenize

  implicit none

  integer(ip), parameter  :: FIRSTITERATION = 1_ip
  integer(ip)             :: ipoin,idime,itott

  public                  :: &
       sld_explicit,         & ! Centred Differences Explicit scheme
       sld_explicit_tw,      & ! Dissipative Explicit Tchamwa-Wielgosz scheme
       sld_implicit,         & ! Newmark-Beta Implicit scheme
       sld_rk4,              & ! Runge-Kutta scheme
       sld_rk4_derivs

contains

  !------------------------------------------------------------------------------
  !> @author  Matias Rivero (matias.rivero@bsc.es)
  !> @author  Adria Quintanas (adria.quintanas@udg.edu)
  !> @author  Gerard Guillamet (gerard.guillamet@bsc.es)
  !> @date    January, 2016
  !>           - Subroutine written
  !> @brief   Centred differences Explicit time integration scheme
  !> @details Central Difference Method Box 6.1 (pag. 332)
  !>          from Belytschko book Second Edition
  !>
  !>          References:\n
  !>          T. Belytschko, W. K. Liu, B. Moran, K. I. Elkhodary
  !>          Nonlinear Finite Elements for Continua and Structures\n
  !>
  !>
  !------------------------------------------------------------------------------

  subroutine sld_explicit()

    use def_master,          only : cutim

    implicit none

    real(rp)          :: t_end, t_mid  !< t^{n+1}, t^{n+1/2}

    if ( INOTEMPTY ) then
       !
       ! Time updates: t^{n+1} and t^{n+1/2}
       !
       t_end = cutim + dtime
       t_mid = 0.5_rp*(cutim + t_end)
       !
       ! Local Axes
       !
       call local_basis_global_to_local(kfl_fixrs_sld,    displ,LAST_COMPONENT=TIME_N)
       call local_basis_global_to_local(kfl_fixrs_sld,veloc_sld,LAST_COMPONENT=TIME_N)
       call local_basis_global_to_local(kfl_fixrs_sld,accel_sld,LAST_COMPONENT=TIME_N)

       do ipoin = 1,npoin
          itott = (ipoin-1)*ndofn_sld
               
          do idime = 1,ndime
             !
             ! First partial update nodal velocities: {v}^{n+1/2}
             !
             ! - Free nodes
             veloctmp_sld(idime,ipoin) = veloc_sld(idime,ipoin,TIME_N) + (t_mid - cutim)*accel_sld(idime,ipoin,TIME_N)
             !
             ! Update nodal displacements: {d}^{n+1}
             !
             unkno(itott+idime) = displ(idime,ipoin,TIME_N) + dtime*veloctmp_sld(idime,ipoin)
             ! - Fixed displacements
             if ( kfl_fixno_sld(idime,ipoin) == 1_ip ) unkno(itott+idime) = bvess_sld(idime,ipoin,ITER_K)
             ! - Fixed displacements (PDN Contact nodes)
             if ( kfl_fixno_sld(1,ipoin)     == 3_ip ) then
                unkno(itott+idime) = bvess_sld(idime,ipoin,ITER_K) + dtime*veloctmp_sld(idime,ipoin)
             end if
             !
             displ(idime,ipoin,ITER_K) = unkno(itott+idime)
          end do

       end do
       !
       ! Local Axes
       !
       call local_basis_local_to_global(kfl_fixrs_sld,unkno                          )
       call local_basis_local_to_global(kfl_fixrs_sld,displ,    LAST_COMPONENT=ITER_K)
       call local_basis_local_to_global(kfl_fixrs_sld,displ,    LAST_COMPONENT=TIME_N)
       call local_basis_local_to_global(kfl_fixrs_sld,veloc_sld,LAST_COMPONENT=TIME_N)
       call local_basis_local_to_global(kfl_fixrs_sld,accel_sld,LAST_COMPONENT=TIME_N) 
       
    end if
    !
    ! FE2
    !
    call fe2_set_strains ()
    call fe2_homogenize ()
    !
    ! Right hand side - get forces
    !
    call sld_matrix(SLD_EXPLICIT_SCHEME)

    if ( INOTMASTER ) then
       if ( kfl_local_sld /= 0_ip .or. kfl_conta_sld /= 0_ip ) then
          !
          ! Fixity for PDN-contact
          !
          call sld_bclocal_rhsid_fixity()
       end if
    end if
    !
    ! Solution
    !   {a}^{n+1} = [M]^{1}*{f}, where f is the RHS
    !
    !   Using the solver:
    solve(1) % xdiag = 1.0_rp
    call array_operations_initialization(dunkn_sld)
    call solver_solve(momod(modul) % solve,amatr,rhsid,dunkn_sld,vmass_sld)
    !   Without solver: comment the above lines and uncomment the <**NO_SOLVER**>
    ! Save number of iterations
    last_iters_sld = solve_sol(1)%iters

    if ( INOTMASTER ) then
       !
       ! Updates
       !
       t_end = cutim + dtime
       t_mid = 0.5_rp*(cutim + t_end)
       do ipoin = 1,npoin
          itott = (ipoin - 1_ip)*ndofn_sld
          !
          ! Acceleration: {a}^{n+1}
          !
          !accel_sld(1:ndime,ipoin,ITER_K) = rhsid(itott+1:itott+ndime) / vmass_sld(ipoin) ! <**NO_SOLVER**>
          accel_sld(1:ndime,ipoin,ITER_K) = dunkn_sld(itott+1:itott+ndime)
          !
          ! Second partial updated nodal velocities: {v}^{n+1}
          !
          veloc_sld(1:ndime,ipoin,ITER_K) = veloctmp_sld(1:ndime,ipoin) + (t_end - t_mid)*accel_sld(1:ndime,ipoin,ITER_K)
          !
          ! Other updates
          ! - Update displacement increment
          ddisp_sld(1:ndime,ipoin,ITER_K) = displ(1:ndime,ipoin,ITER_K) - displ(1:ndime,ipoin,TIME_N)
          ! - Update internal and exteranl forces
          fintt_sld(1:ndime,ipoin,ITER_K) = finte_sld(itott+1:itott+ndime)
          fextt_sld(1:ndime,ipoin,ITER_K) = fexte_sld(itott+1:itott+ndime)
       end do
       !
       ! Local Axes
       !
       call local_basis_local_to_global(kfl_fixrs_sld,veloc_sld,LAST_COMPONENT=ITER_K)
       call local_basis_local_to_global(kfl_fixrs_sld,accel_sld,LAST_COMPONENT=ITER_K) 

    end if
    !
    ! Check energy balance at time step {n+1}
    !
    call sld_energy()

  end subroutine sld_explicit

  !------------------------------------------------------------------------------
  !> @author  Matias Rivero (matias.rivero@bsc.es)
  !>          Adria Quintanas (adria.quintanas@udg.edu)
  !>          Gerard Guillamet (gerard.guillamet@bsc.es)
  !> @date    January, 2016
  !> @brief   Beta-Newmark implicit time integration scheme
  !> @details Box 6.1 (pag. 313) from Belytschko book.
  !>
  !>          References:\n
  !>          T. Belytschko, W. K. Liu, B. Moran, K. I. Elkhodary
  !>          Nonlinear Finite elements for Continua and Structures\n
  !>
  !------------------------------------------------------------------------------

  subroutine sld_implicit()

    implicit none

    real(rp)     :: dt, dt2                           ! Time parameters
    real(rp)     :: nmkbeta, nmkgamma                 ! Newmark parameters
    !
    ! Newmark initialization
    !
    if ( INOTMASTER ) then

       ! Times
       dt  = dtime
       dt2 = dtime*dtime

       ! Newmark parameters
       nmkbeta = tifac_sld(1)
       nmkgamma= tifac_sld(2)

       ! Initial guess for the first iteration
       if ( itinn(modul) == FIRSTITERATION ) then

          ! Dynamic problem
          if ( kfl_timei_sld == SLD_DYNAMIC_PROBLEM ) then    ! in case of dynamic problem
             ! displacement, acceleration and velocities and Dirichlet boundary conditions
             do ipoin = 1, npoin
                itott = (ipoin-1)*ndofn_sld

                do idime=1,ndime
                   ! Free nodes
                   unknotmp_sld(itott+idime) = displ(idime,ipoin,TIME_N) + veloc_sld(idime,ipoin,TIME_N)*dt &
                        + (0.5_rp - nmkbeta)*accel_sld(idime,ipoin,TIME_N)*dt2
                   unkno(itott+idime) = unknotmp_sld(itott+idime)
                   displ(idime,ipoin,ITER_K) = unknotmp_sld(itott+idime)
                   veloctmp_sld(idime,ipoin) = veloc_sld(idime,ipoin,TIME_N) + (1.0_rp - nmkgamma)*accel_sld(idime,ipoin,TIME_N)*dt
                   ! Fixed nodes
                   if( kfl_fixno_sld(idime,ipoin) > 0 .and. kfl_fixno_sld(idime,ipoin) /= 3 ) then
                      unkno(itott+idime)            = bvess_sld(idime,ipoin,ITER_K)
                      displ(idime,ipoin,ITER_K)     = bvess_sld(idime,ipoin,ITER_K)
                   end if
                end do
                !
                ! Local Axes (Local --> Global)
                !
                call sld_bclocal_displ(ndofn_sld,ipoin)

             end do

          ! Static problem (or equilibrium problem)
          else if ( kfl_timei_sld == SLD_STATIC_PROBLEM ) then
             ! displacement and Dirichlet boundary conditions
             do ipoin = 1, npoin
                itott = (ipoin-1)*ndofn_sld

                do idime = 1, ndime
                   ! Free nodes
                   unkno(itott+idime) = displ(idime,ipoin,TIME_N)
                   ! Fixed nodes
                   if( kfl_fixno_sld(idime,ipoin) > 0 .and. kfl_fixno_sld(idime,ipoin) /= 3 ) then
                      unkno(itott+idime)            = bvess_sld(idime,ipoin,ITER_K)
                      displ(idime,ipoin,ITER_K)     = bvess_sld(idime,ipoin,ITER_K)
                   end if
                end do
                !
                ! Local Axes (Local --> Global)
                !
                call sld_bclocal_displ(ndofn_sld,ipoin)

             end do
          end if

       end if
       !
       ! Newmark corrections
       !
       if ( kfl_timei_sld == SLD_DYNAMIC_PROBLEM ) then
          do ipoin = 1,npoin
             itott = (ipoin-1)*ndofn_sld
             accel_sld(1:ndime,ipoin,ITER_K) = &
                  (1.0_rp / (nmkbeta*dt2)) * ( unkno(itott+1:itott+ndime) - unknotmp_sld(itott+1:itott+ndime) )
             veloc_sld(1:ndime,ipoin,ITER_K) = veloctmp_sld(1:ndime,ipoin) + nmkgamma*dt*accel_sld(1:ndime,ipoin,ITER_K)
          end do
       end if

    end if
    !
    ! FE2
    !
    call fe2_set_strains ()
    call fe2_homogenize ()
    !
    ! Assemble Jacobian matrix (J) and RHS (r)
    !
    call sld_matrix(SLD_IMPLICIT_SCHEME)
    !
    ! Solve
    !
    call array_operations_initialization(dunkn_sld)

    call solver_solve(momod(modul) % solve,amatr,rhsid,dunkn_sld)

    ! Save solver iterations
    last_iters_sld = solve_sol(1)%iters

    if ( INOTMASTER ) then
       !
       ! Correcting displacement and imposing boundary conditions, for both STATIC and DYNAMIC
       !
       do ipoin = 1,npoin
          itott = (ipoin-1)*ndofn_sld
          !
          ! Local Axes: correct rotation back (Local --> Global)
          !
          call sld_bclocal_roback(ndofn_sld,ipoin)
          !
          ! Updates
          !
          ! Update displacement
          unkno(itott+1:itott+ndime) = unkno(itott+1:itott+ndime) + dunkn_sld(itott+1:itott+ndime)
          ! Update displacement increment
          ddisp_sld(1:ndime,ipoin,ITER_K) = unkno(itott+1:itott+ndime) - displ(1:ndime,ipoin,TIME_N)
          ! Update global force vectors
          fintt_sld(1:ndime,ipoin,ITER_K) = finte_sld(itott+1:itott+ndime)
          fextt_sld(1:ndime,ipoin,ITER_K) = fexte_sld(itott+1:itott+ndime)
       end do
    end if
    !
    ! Check energy balance
    !  
    call sld_energy()

  end subroutine sld_implicit

  !------------------------------------------------------------------------------
  !> @author  Gerard Guillamet (gerard.guillamet@bsc.es)
  !> @date    April, 2019
  !> @brief   Dissipative Tchamwa-Wielgosz scheme
  !> @details Dissipative Tchamwa-Wielgosz scheme
  !>
  !>          References:\n
  !>          L. Mah\'{\e}o, V. Grolleau, G. Rio. Damping efficiency of the\n
  !>          Tchamwa-Wielgosz explicit dissipative scheme under instantaneous\n
  !>          loading conditions\n
  !>
  !------------------------------------------------------------------------------
  
  subroutine sld_explicit_tw()

    implicit none

    real(rp)          :: phi           !< Parameter controlling the TW scheme

    if( INOTMASTER ) then

       ! TW parameter
       phi = tifac_sld(4)
       !
       ! Local Axes
       !
       call local_basis_global_to_local(kfl_fixrs_sld,    displ,LAST_COMPONENT=TIME_N)
       call local_basis_global_to_local(kfl_fixrs_sld,veloc_sld,LAST_COMPONENT=TIME_N)
       call local_basis_global_to_local(kfl_fixrs_sld,accel_sld,LAST_COMPONENT=TIME_N)

       if( kfl_bvesv_sld == 0 ) then
          !
          ! Prescribed displacement
          !
          do ipoin = 1,npoin
             itott = (ipoin-1)*ndofn_sld                    
             do idime = 1,ndime
                !
                ! Update nodal velocities: {v}^{n+1}
                !
                ! - Free nodes
                veloc_sld(idime,ipoin,ITER_K) = veloc_sld(idime,ipoin,TIME_N) + dtime*accel_sld(idime,ipoin,TIME_N)
                !
                ! Update nodal displacements: {d}^{n+1}
                !
                unkno(itott+idime) = displ(idime,ipoin,TIME_N) &
                     + dtime*veloc_sld(idime,ipoin,TIME_N) + phi*(dtime**2)*accel_sld(idime,ipoin,TIME_N)
                !
                ! Prescribed displacement
                !
                if( kfl_fixno_sld(idime,ipoin) == 1_ip ) unkno(itott+idime) = bvess_sld(idime,ipoin,ITER_K)
                if( kfl_fixno_sld(1,ipoin)        == 3_ip ) then
                   unkno(itott+idime) = bvess_sld(idime,ipoin,ITER_K) &
                        + dtime*veloc_sld(idime,ipoin,TIME_N) + phi*(dtime**2)*accel_sld(idime,ipoin,TIME_N)
                end if
                displ(idime,ipoin,ITER_K) = unkno(itott+idime)
             end do
          end do
       else
          !
          ! Prescribed velocity
          !
          do ipoin = 1,npoin
             itott = (ipoin-1)*ndofn_sld                    
             do idime = 1,ndime
                !
                ! Update nodal velocities: {v}^{n+1}
                !
                ! - Free nodes
                veloc_sld(idime,ipoin,ITER_K) = veloc_sld(idime,ipoin,TIME_N) + dtime*accel_sld(idime,ipoin,TIME_N)
                !
                ! Update nodal displacements: {d}^{n+1}
                !
                unkno(itott+idime) = displ(idime,ipoin,TIME_N) &
                     + dtime*veloc_sld(idime,ipoin,TIME_N) + phi*(dtime**2)*accel_sld(idime,ipoin,TIME_N)
                !
                ! Prescribed velocities
                !
                if( kfl_fixno_sld(idime,ipoin) == 1_ip ) then
                   if( abs(bvess_sld(idime,ipoin,ITER_K)) > 0.0_rp ) then
                      veloc_sld(idime,ipoin,ITER_K) = bvess_sld(idime,ipoin,ITER_K)
                   else
                      unkno(itott+idime)            = 0.0_rp
                      veloc_sld(idime,ipoin,ITER_K) = 0.0_rp
                   end if
                end if
                displ(idime,ipoin,ITER_K) = unkno(itott+idime)
             end do
          end do
       end if
       !
       ! Local Axes
       !
       call local_basis_local_to_global(kfl_fixrs_sld,unkno                          )
       call local_basis_local_to_global(kfl_fixrs_sld,displ,    LAST_COMPONENT=ITER_K)
       call local_basis_local_to_global(kfl_fixrs_sld,displ,    LAST_COMPONENT=TIME_N)
       call local_basis_local_to_global(kfl_fixrs_sld,veloc_sld,LAST_COMPONENT=TIME_N)
       call local_basis_local_to_global(kfl_fixrs_sld,accel_sld,LAST_COMPONENT=TIME_N)

    end if
    !
    ! Right hand side - get forces
    !
    call sld_matrix(SLD_EXPLICIT_SCHEME)

    if( INOTMASTER ) then

       if ( kfl_conta_sld /= 0_ip ) then
          !
          ! Fixity for PDN-contact
          !
          call sld_bclocal_rhsid_fixity()

       end if

    end if
    !
    ! Solution
    !   {a}^{n+1} = [M]^{1}*{f}, where f is the RHS
    !
    solve(1) % xdiag = 1.0_rp
    call array_operations_initialization(dunkn_sld)
    call solver_solve(momod(modul) % solve,amatr,rhsid,dunkn_sld,vmass_sld)
    !   Without solver: comment the above lines and uncomment the <**NO_SOLVER**>
    ! Save number of iterations
    last_iters_sld = solve_sol(1)%iters

    if( INOTMASTER ) then
       if( kfl_bvesv_sld == 0 ) then
          !
          ! Prescribed displacement
          !
          do ipoin = 1,npoin
             itott = (ipoin - 1_ip)*ndofn_sld
             do idime = 1,ndime
                !
                ! Acceleration: {a}^{n+1}
                !
                !accel_sld(idime,ipoin,ITER_K) = rhsid(itott+idime) / vmass_sld(ipoin) ! <**NO_SOLVER**>
                accel_sld(idime,ipoin,ITER_K) = dunkn_sld(itott+idime)
                !
                ! Impose Dirichlet condition
                !
                if( kfl_fixno_sld(idime,ipoin) == 1_ip ) then
                   veloc_sld(idime,ipoin,ITER_K) = 0.0_rp
                   accel_sld(idime,ipoin,ITER_K) = 0.0_rp
                end if
                !
                ! Other updates
                ! - Update displacement increment
                ddisp_sld(idime,ipoin,ITER_K) = displ(idime,ipoin,ITER_K) - displ(idime,ipoin,TIME_N)
                ! - Update internal and exteranl forces
                fintt_sld(idime,ipoin,ITER_K) = finte_sld(itott+idime)
                fextt_sld(idime,ipoin,ITER_K) = fexte_sld(itott+idime)
             end do
          end do
       else
          !
          ! Prescribed velocity
          !
          do ipoin = 1,npoin
             itott = (ipoin - 1_ip)*ndofn_sld
             do idime = 1,ndime
                !
                ! Acceleration: {a}^{n+1}
                !
                !accel_sld(idime,ipoin,ITER_K) = rhsid(itott+idime) / vmass_sld(ipoin) ! <**NO_SOLVER**>
                accel_sld(idime,ipoin,ITER_K) = dunkn_sld(itott+idime)
                !
                ! Impose Dirichlet condition
                !
                if( kfl_fixno_sld(idime,ipoin) == 1_ip ) then
                   if( abs(bvess_sld(idime,ipoin,ITER_K)) > 0.0_rp ) then
                      accel_sld(idime,ipoin,ITER_K) = 0.0_rp
                   else
                      displ(idime,ipoin,ITER_K)     = 0.0_rp
                      veloc_sld(idime,ipoin,ITER_K) = 0.0_rp
                      accel_sld(idime,ipoin,ITER_K) = 0.0_rp
                   end if
                end if
                !
                ! Other updates
                ! - Update displacement increment
                ddisp_sld(idime,ipoin,ITER_K) = displ(idime,ipoin,ITER_K) - displ(idime,ipoin,TIME_N)
                ! - Update internal and exteranl forces
                fintt_sld(idime,ipoin,ITER_K) = finte_sld(itott+idime)
                fextt_sld(idime,ipoin,ITER_K) = fexte_sld(itott+idime)
             end do
          end do
       end if
       !
       ! Local Axes
       !
       call local_basis_local_to_global(kfl_fixrs_sld,veloc_sld,LAST_COMPONENT=ITER_K)
       call local_basis_local_to_global(kfl_fixrs_sld,accel_sld,LAST_COMPONENT=ITER_K)
       
    end if
    !
    ! Check energy balance at time step {n+1}
    !
    call sld_energy()

  end subroutine sld_explicit_tw

  !-----------------------------------------------------------------------
  !> @author  Gerard Guillamet
  !> @date    January, 2019
  !>          - Subroutine written
  !> @brief   RK4 scheme
  !> @details RK4 scheme
  !>
  !>          References:\n
  !>          W. H. Press, Saul A. Teukolsky, W. T. Vetterling and\n
  !>          Brian P. Flannery. Numerical Recipes in Fortran 77: The
  !>          Art of Scientific Computing\n
  !>
  !-----------------------------------------------------------------------

  subroutine sld_rk4(itask,dt,nderi,u)

    implicit none

    integer(ip), intent(in)    :: itask    !< What to solve
    integer(ip), intent(in)    :: nderi    !< Number of derivatives
    real(rp),    intent(in)    :: dt       !< Time step size
    real(rp),    intent(inout) :: u(nderi) !< Solution

    real(rp) :: k1(nderi), k2(nderi), k3(nderi), k4(nderi)
    real(rp) :: uk(nderi)
    !
    ! Initialize
    !
    uk = 0.0_rp
    k1 = 0.0_rp
    k2 = 0.0_rp
    k3 = 0.0_rp
    k4 = 0.0_rp
    !
    ! Calculate partial steps
    !
    call sld_rk4_derivs(itask,nderi,u(:),k1(:))

    uk(1:nderi) = u(1:nderi) + k1(1:nderi)*dt/2.0_rp
    call sld_rk4_derivs(itask,nderi,uk(:),k2(:))

    uk(1:nderi) = u(1:nderi) + k2(1:nderi)*dt/2.0_rp
    call sld_rk4_derivs(itask,nderi,uk(:),k3(:))

    uk(1:nderi) = u(1:nderi) + k3(1:nderi)*dt
    call sld_rk4_derivs(itask,nderi,uk(:),k4(:))
    !
    ! Combine partial steps
    !
    u(:) = u(:) + (k1(:) + 2.0_rp*(k2(:) + k3(:)) + k4(:))*dt/6.0_rp

  end subroutine sld_rk4

  !-----------------------------------------------------------------------
  !> @author  Gerard Guillamet
  !> @date    January, 2019
  !> @brief   Functions for Runge-Kutta scheme
  !> @details This routine includes the right hand side of n first-order
  !>          differential equations.
  !>
  !>          1.
  !>          2.
  !>          3.
  !>
  !-----------------------------------------------------------------------

  subroutine sld_rk4_derivs(itask,nderi,u,deriv)

    use def_solidz, only : grnor_sld, gravi_sld
    use def_solidz, only : rbmas_sld
    use def_solidz, only : rbfor_sld

    implicit none

    integer(ip), intent(in)  :: itask        !< Select case
    integer(ip), intent(in)  :: nderi        !< Number of derivatives
    real(rp),    intent(in)  :: u(nderi)     !< Value of dependent variable
    real(rp),    intent(out) :: deriv(nderi) !< Value of the derivative

    select case (itask)

    case( 1_ip )
       !
       ! Falling sphere:  m * dx^2/dt^2 = m * g - Fext
       ! dx/dt      = v
       ! dx^2/dt*2  = dv/dt = g - Fext/m
       !
       ! Linear velocities
       deriv(1:ndime) = u(ndime+1:ndime*2)
       ! Linear accelerations
       deriv(ndime+1:ndime*2) = grnor_sld*gravi_sld(1:ndime) - rbfor_sld(1:ndime)/rbmas_sld

    end select

  end subroutine sld_rk4_derivs

end module mod_sld_solution_methods
