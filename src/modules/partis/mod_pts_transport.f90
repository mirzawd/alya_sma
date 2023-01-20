!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Partis
!> @{
!> @file    mod_pts_transport.f90
!> @author  houzeaux
!> @date    2018-09-21
!> @brief   Module for particle transport
!> @details Contains all subroutine to compute the transport of particles
!>          
!-----------------------------------------------------------------------
module mod_pts_transport

  use def_master,              only : vorti
  use def_kintyp_comm,         only : comm_data_par
  use def_master,              only : momod,mem_modul,modul,igene
  use def_master,              only : cutim,dtime,dtinv,ittim
  use def_master,              only : nparr,zeror,pard1
  use def_master,              only : advec,gisca,CPU_ASSEMBLY
  use def_master,              only : npasr,velom
  use def_master,              only : ISLAVE,IMASTER,IPARALL,INOTEMPTY
  use def_master,              only : INOTMASTER,INOTSLAVE
  use def_master,              only : momentum_sink
  use def_master,              only : mass_sink
  use def_master,              only : heat_sink
  use def_master,              only : lelbf,leinv_loc,cpu_modul
  use def_master,              only : ID_NASTIN,ID_PARTIS
  use def_master,              only : vesgs,kfl_paral
  use def_parame,              only : pi
  use def_kermod,              only : gravi,grnor,kfl_detection,relse
  use def_domain,              only : lnods,coord,ltype,nnode,npoin
  use def_domain,              only : ngaus
  use def_domain,              only : lelel_2,pelel_2,nelem,leldo
  use def_domain,              only : ndime,mnode,mgaus,walln
  use def_domain,              only : element_bin
  use def_domain,              only : element_bin_boxes
  use def_domain,              only : walld
  use def_domain,              only : elmar
  use def_elmtyp,              only : TET04,TRI03
  use mod_ker_proper,          only : ker_proper
  use mod_elmgeo,              only : elmgeo_natural_coordinates
  use mod_elmgeo,              only : elmgeo_gauss_to_element
  use mod_memory,              only : memory_alloca
  use mod_memory,              only : memory_deallo
  use mod_maths,               only : maths_mapping_coord_to_3d
  use mod_parall,              only : PAR_CODE_SIZE
  use mod_communications,      only : PAR_ALLGATHER
  use mod_communications,      only : PAR_SEND_RECEIVE
  use mod_communications,      only : PAR_DEFINE_COMMUNICATOR
  use mod_communications,      only : PAR_COMM_RANK_AND_SIZE
  use mod_communications,      only : PAR_SUM
  use mod_communications,      only : PAR_ALLGATHER
  use mod_communications,      only : PAR_ALLGATHERV
  use mod_communications,      only : PAR_MAX
  use mod_communications,      only : PAR_MIN
  use mod_ker_detection,       only : ker_events_particle_not_converged
  use mod_maths,               only : maths_local_orthonormal_basis
  use mod_maths,               only : maths_vector_to_new_basis
  use mod_maths,               only : maths_vector_from_new_basis
  use mod_ker_timeline,        only : ker_timeline
  use mod_physics,             only : physics_drag_force
  use mod_physics,             only : physics_lift_force 
  use mod_physics,             only : physics_set_liquid_temperature
  use mod_physics,             only : physics_Cunningham
  use mod_pts_parallelization, only : pts_parallelization_migration
  use mod_messages,            only : livinf
  use mod_messages,            only : messages_live
  use mod_pts_particle,        only : pts_particle_mass
  use mod_pts_particle,        only : pts_particle_diameter
  use mod_pts_thermodynamic,   only : pts_thermodynamic_transport
  use mod_pts_thermodynamic,   only : pts_thermodynamic_properties
  use mod_iofile,              only : iofile_flush_unit
  use mod_timings,             only : timings_assembly
  use mod_random,              only : random_initialization
  use mod_random,              only : random_end
  use mod_physics,             only : physics_T_2_HCp
  use mod_output_postprocess,  only : output_postprocess_check_variable_postprocess
  use mod_projec,              only : projec_elements_to_nodes
  use mod_solver,              only : solver_lumped_mass_system
  use mod_strings,             only : integer_to_string
  use def_partis  
  use def_mpi
#include "def_mpi.inc"
  
  implicit none
  private

  real(rp) :: toler          ! Tolerance for Elsest
  real(rp) :: ti,tf          ! Time limits
  real(rp) :: strex,ovstr    ! Adaptive time step
  real(rp) :: g(3)           ! Gravity

  public :: pts_transport_particles
  public :: pts_transport_finalize

contains

  !------------------------------------------------------------------------
  !> @addtogroup Partis
  !> @{
  !> @name   Partis inner iteration
  !> @file    pts_soltie.f90
  !> @author  Guillaume Houzeaux
  !> @date    28/06/2012
  !> @brief   This routine integrates particle paths
  !> @details
  !>
  !>    ALGORITHM
  !>    ---------
  !>
  !>       We are at time n+1:
  !>
  !>       do while( no han llegado todas a su tiempo )
  !>
  !>          Loop over particles
  !>            do while( t < tf )
  !>              Start from last element: ielem
  !>              x^n+1 = x^n + dt * u
  !>              Look for neighbors of ielem
  !>              - if found
  !>                  Validate new position
  !>                  x^n = x^n+1
  !>                  Save element
  !>                  if( fringe element ) save particle to be sent
  !>              - else
  !>                  dt <= dt / 2
  !>              - endif
  !>              t = t + dt
  !>            end do
  !>          end do
  !>          Send what I have to send
  !>          Receive what I have to receive
  !>       end do
  !>
  !>    NEWMARK
  !>    -------
  !>
  !>    The Newmark scheme has a small inconsistency. The mag force is computed at k
  !>    and not k+1
  !>
  !>    v^{k+1}_NS => F_drag                             x^{k+1} = f(x^k, a^{k+1}, a^k)
  !>                           => F^{k+1} => a^{k+1} =>
  !>    x^k        => F_mag                              u^{k+1} = f(a^{k+1}, a^k)
  !>
  !>    x^{n+1} = x^n + u^n dt + 1/2 a^n dt^2 + 1/6 a^n dt^3
  !>    a^{n+1} = a^n + a'^n dt + 1/2 a''^n dt^2 =>
  !>    a'^n dt = ( a^{n+1} - a^n ) / dt =>
  !>
  !>    u^{n+1} = u^n + a^n dt + gamma dt ( a^{n+1} - a^n )
  !>    x^{n+1} = x^n + u^n dt + 1/2 a^n dt^2 + beta dt^2 ( a^{n+1} - a^n )
  !>
  !>    x^{n+1} = x^n + dt [ u^{n+1} - a^n dt - gamma * dt ( a^{n+1} - a^n ) ] + 1/2 a^n dt^2 + beta dt^2 ( a^{n+1} - a^n )
  !>    x^{n+1} = x^n + dt u^{n+1} - a^n dt^2 - gamma * dt^2 ( a^{n+1} - a^n ) ] + 1/2 a^n dt^2 + beta dt^2 ( a^{n+1} - a^n )
  !>    x^{n+1} = x^n + dt u^{n+1} + dt^2 ( - a^n  - gamma * ( a^{n+1} - a^n ) + 1/2 a^n + beta ( a^{n+1} - a^n ) ]
  !>    x^{n+1} = x^n + dt u^{n+1} + dt^2 [ (-1/2+gamma-beta) a^n  + (-gamma+beta)  a^{n+1} ]
  !>
  !>                            gamma   beta
  !>    Fox Goodwin             1 / 2   1 / 12    conditionnaly stable
  !>    Linear acceleration     1 / 2   1 / 6     conditionnaly stable
  !>    Average acceleration    1 / 2   1 / 4     unconditionnaly stable
  !>    Nosotros                0.75    0.390625  diffusive
  !>    Nosotros                1.0     0.5625    super diffusive
  !>    External                1 / 2   0         explicit
  !>
  !>    Relation between beta and gamma: beta = 0.25*(gamma+1/2)^2
  !>
  !>    RANDOM WALK
  !>    -----------
  !>
  !>    It is also possible to simulate the diffusion process by using a random walk
  !>    algorithm. Once the new position of the particle has been found, the
  !>    position is corrected as:
  !>
  !>    x^{n+1} <= x^{n+1} + eps * (2D*dt)^1/2,
  !>
  !>    where D is the diffusion coefficient [m^2/s] and eps follows a normal distribution.
  !>
  !>    The diffusion is given by
  !>           k*T
  !>    D = ---------
  !>        6*pi*mu*r
  !>
  !>    where T is temperature, k is Boltzmann's constant, mu viscosity and r ths solute radius.
  !>
  !>    The process is the following:
  !>
  !>    1. Generate two random numbers U1 and U2 in [0,1]
  !>    2. Perform a box-muller transformation to generate from U1 and U2 a normal distribution
  !>       eps1 = sqrt( -2 *log(U1) ) * cos(2*pi*U2)
  !>       eps2 = sqrt( -2 *log(U1) ) * sin(2*pi*U2)
  !>    3. Update the position
  !>       x <= x + eps1 * (2D*dt)^1/2
  !>       y <= y + eps1 * (2D*dt)^1/2
  !>
  !>    A simple test can be peformed to check the normal distribution:
  !>
  !>       The probablity P(r,t) is the normal defined as
  !>
  !>       P(r,t) = 1/(4*D*pi*t)^(d/2) * exp( -r^2/(4*D*t) )
  !>
  !>       for a d-dimensional diffusion problem.
  !>
  !>       k is related to the diffusion as k^2 = 2D and P satisfies the
  !>       2-dimensional Poisson equation:
  !>       dP/dt = D ( d^2P/dx^2 + d^2P/dy^2 )
  !>
  !>    1. Inject particles at 0.05,0.05 in a domain [0,0] x [0.01,0.01].
  !>    2. Choose zero convection velocity; D = 5.35*10^-5 [m^2/s];
  !>       dt = 10^-5 [s]. With dt = 10^-4 [s], statistics is bad.
  !>    3. At a given time, say t=0.01 [s], count the number of particles
  !>       with radius in ranges from [0:R] where R goes from 0 to 0.01.
  !>    4. Divide the results by the total number of particles.
  !>    5. Compare the results with the following 2D normal distribution:
  !>
  !>       f(R) = \int_0^{2pi} dtheta \int_0^R r P(r,t) dr
  !>            = 1 - exp( -R^2/(4*D*t) )
  !>
  !>       using \int x*exp(-c*x^2) dx = -1/2c * exp(-c*x^2)
  !>
  !>       Note that in 2D we have \int_0^2pi dtheta \int_0^infty r P(r,t) dr = 1.
  !>
  !>       In 3D, ther result is:
  !>
  !>       f(R) = \int_0^{2pi} dtheta \int_0^pi sin(phi) dphi \int_0^R r^2 P(r,t) dr
  !>            = 4*pi/(4*D*pi*t)^3/2 \int_0^R r^2 * exp( -r^2/(4*D*t) ) dr
  !>
  !>    INJECTION
  !>    ---------
  !>
  !>       KFL_MODLA ......... Model for Lagrangian transport (0=no particle,1=velocity,1=drag)
  !>       KFL_INJLA ......... Injection model
  !>       TINLA ............. Initial time of injection
  !>       TPELA ............. Time period of injection
  !>       MPALA ............. Maximum number of parameters
  !>       PARLA(MPALA) ...... Parameters for the injection
  !>
  !>    OTHERS
  !>    ------
  !>
  !>       MLAGR ......................... Max. number of particle in each subdomain
  !>       NLAGR ......................... Number of particles in each subdomain (just needed for info)
  !>       NLACC_PTS ..................... Total number of existing and disappeared particles in all subdomains
  !>       NLAGR_EXISTING_PTS ............ Total number of existing particles in all subdomains
  !>       NLAGR_NON_MIGRATING_PTS ....... Particles going from one subdomain to another
  !>       NLAGR_GOING_OUT_PTS ........... Particles deposited: going out of the computational domain
  !>       NLAGR_ZERO_TIME_PTS ........... Particles that vanish because of zero time step
  !>       NLAGR_DEPOSITED_PTS ........... Particles deposited, bue boundary not found
  !>       PELEL_2,LELEL_2 ............... Element connectivity linked list including
  !>                                       the neighbors element in parallel
  !>
  !>    LAGRTYP definition
  !>    ------------------
  !>
  !>       Transport of Lagrangian particles. For particle ilagr
  !>       lagrtyp(ilagr) % ilagr     =  i ......... Particle number
  !>       lagrtyp(ilagr) % kfl_exist = -1 ......... I am the owner
  !>                                  = -2 ......... Particle is deposited: go out of domain
  !>                                  = -3 ......... Particle vanishes: time step too small
  !>                                  = -4 ......... Particle is out of the flow
  !>                                  = -5 ......... Particle has just been injected
  !>                                  = -6 ......... Particle fully evaporated
  !>                                  =  i ......... Send particle to neighbor i
  !>                                  =  0 ......... No particle at that position ILAGR
  !>       lagrtyp(ilagr) % coord(3)  = x,y,z ...... Coordinates
  !>       lagrtyp(ilagr) % veloc(3)  = vx,vy,vz ... Old velocity
  !>       lagrtyp(ilagr) % t         = t .......... Current time
  !>
  !> @}
  !------------------------------------------------------------------------

  subroutine pts_transport_particles()
#ifdef _OPENMP
    use omp_lib
#endif
    implicit none
#ifdef ALYA_DLB
  include 'dlbf.h'
#endif
    integer(ip)                      :: ielem,ielel,jelem,dumm0,inode
    integer(ip)                      :: nneig,ilagr,ilagr_local
    real(rp)                         :: time1,time2
    !
    ! Parallelization
    !
    integer(ip)                      :: comcont
    MY_MPI_COMM                      :: PAR_COMM_TO_USE4
    integer(ip),         allocatable :: nlagr_migrating(:)
    type(comm_data_par), pointer     :: commu
    !
    ! Special modes
    !
    logical(lg)                      :: debugmode
    !
    ! Options
    !
    debugmode = .false.     
    !
    ! Random numbers
    !
    call random_initialization(broadcast_seed=.false.)
    !
    ! Define communication array
    !
    nullify(commu)
    if( IPARALL ) call PAR_DEFINE_COMMUNICATOR('IN MY ZONE',PAR_COMM_TO_USE4,commu)

    !----------------------------------------------------------------------
    !
    ! USEFUL FOR DEBUGGING IN PARALLEL:
    !
    ! Reorder elements linked list to go through the neighboring elements
    ! in the same order in sequential and parallel
    !
    !----------------------------------------------------------------------

    if( debugmode .and. INOTMASTER ) then
       do ielem = 1,nelem
          dumm0 = pelel_2(ielem+1) - pelel_2(ielem)
          call memgen(1_ip,dumm0,0_ip)
          inode = 0
          do ielel = pelel_2(ielem),pelel_2(ielem+1)-1
             jelem = lelel_2(ielel)
             inode = inode + 1
             gisca(inode) = leinv_loc(jelem)
          end do
          call heapsorti2(1_ip,inode,gisca,lelel_2(pelel_2(ielem)))
          call memgen(3_ip,dumm0,0_ip)
       end do
    end if

    !----------------------------------------------------------------------
    !
    ! Global number for particles: NLAGR, starting from NLACC_PTS
    !
    !----------------------------------------------------------------------
    !
    ! Allocate memory to communicate with my neighbors
    !
    if( IPARALL ) then
       nneig = commu % nneig
    else
       nneig = 0
    end if
    !
    ! Temporary array to send/receive to my neighbors
    !
    allocate(nlagr_migrating(max(1_ip,nneig)))

    !----------------------------------------------------------------------
    !
    ! Numerical constants and initialization
    !
    !----------------------------------------------------------------------

    toler              = relse(1)           ! Tolerance for element search
    time_transport_pts = 0.0_rp             ! Total CPU time
    particles_sent     = 0                  ! # of particles received
    particles_recv     = 0                  ! # of particles sent
    comm_loops_pts     = 0                  ! # communication loops
    comcont            = 0                  ! Communication loop active
    nlagr_existing_pts = 0                  ! Number of existing particles
    nlagr_local_pts    = 0                  ! Taken positions in particle type
    ti                 = cutim - dtime      ! Initial time
    tf                 = cutim              ! Final time: t^{n+1}
    ovstr              = 1.2_rp             ! Inverse stretching factor
    strex              = 1.0_rp / ovstr     ! Stretching factor
    g                  = grnor * gravi      ! Gravity g

    !----------------------------------------------------------------------
    !
    ! NLAGR_EXISTING_PTS= Number of existing particles (should be equal to nlacc_pts)
    !
    !----------------------------------------------------------------------

    do ilagr = 1,mlagr
       if( lagrtyp(ilagr) % kfl_exist == -1 ) then
          nlagr_existing_pts = nlagr_existing_pts + 1
          lagrtyp(ilagr) % t = ti
       end if
    end do

    call PAR_SUM(nlagr_existing_pts,'IN MY ZONE')

    call messages_live(&
         'TRANSPORT LAGRANGIAN PARTICLES= '//integer_to_string(nlagr_existing_pts),&
         'MODULE')

    call ker_timeline('INI_SOLVER',nlagr_existing_pts)

    !----------------------------------------------------------------------
    !
    ! Loop over communication iterations
    !
    !----------------------------------------------------------------------

    do while( comcont /= -1 )
       !
       ! Communcation variables
       !
       comcont                 = comcont + 1
       comm_loops_pts          = comm_loops_pts + 1
       nlagr_migrating         = 0

       call pts_compute_permutation()

       call cputim(time1)
#ifdef ALYA_DLB
  if( dlb_enable() < DLB_SUCCESS ) call par_livinf(17_ip,'DLB COULD NOT BE ENABLED',0_ip) 
#endif

       !----------------------------------------------------------------------
       !
       ! Loop over particles
       !
       !----------------------------------------------------------------------
       !
       !$OMP PARALLEL DO SCHEDULE (DYNAMIC,1000)                    &
       !$OMP DEFAULT      (NONE)                                    &
       !$OMP PRIVATE      (ilagr_local,ilagr)                       &
       !$OMP SHARED       (nlagr_local_pts,permu_nlagr_pts)         &
       !$OMP REDUCTION    (+:nlagr_non_migrating_pts,nlagr_migrating)
       !
       do ilagr_local = 1,nlagr_local_pts
          ilagr = permu_nlagr_pts(ilagr_local)
          call pts_transport_single_particle(ilagr,nlagr_non_migrating_pts,nlagr_migrating)
       end do
       !
       !$OMP END PARALLEL DO
       !
       call cputim(time2)
       time_transport_pts = time_transport_pts + time2 - time1
       !
       ! All reduce total number of non-crossing particles
       ! if NLAGR_EXISTING_PTS /= NLAGR_NON_MIGRATING_PTS, some particles migrate to other subdomains
       !
       call PAR_SUM(nlagr_non_migrating_pts,'IN MY ZONE')
#ifdef ALYA_DLB
  if( INOTMASTER ) then
     if( dlb_disable() < DLB_SUCCESS ) call par_livinf(17_ip,'DLB COULD NOT BE DISABLED',0_ip) 
  end if
#endif
       if( nlagr_existing_pts == nlagr_non_migrating_pts ) then
          !
          ! Communication is over
          !
          comcont = -1
       else if( ISLAVE ) then
          !
          ! Send/receive number of particles and recompute permutation array
          !
          call pts_parallelization_migration(commu,nlagr_migrating,particles_sent,particles_recv)
       end if

    end do

    deallocate(nlagr_migrating)
    call random_end()

    call timings_assembly(time_transport_pts,TYPE_OF_ASSEMBLY='PARTICLE')
    call ker_timeline('END_SOLVER',nlagr_existing_pts)    


  end subroutine pts_transport_particles

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-09-21
  !> @brief   Transport a particle
  !> @details Transport a partile with local ID ilagr
  !> 
  !-----------------------------------------------------------------------


  subroutine pts_transport_single_particle(ilagr,nlagr_non_migrating_pts,nlagr_migrating)

    integer(ip), intent(in)   :: ilagr                                    !< Particle number
    integer(ip), intent(out)  :: nlagr_non_migrating_pts                  !< Number of existing particles
    integer(ip), intent(out)  :: nlagr_migrating(:)                       !< Number of particle migrate to neighbors

    integer(ip)               :: ielem,dumm0,jelem,pelty,pnode,inode
    integer(ip)               :: peltj,pnodj
    integer(ip)               :: ipoin,idime,ifoun,ineig,itype,iboun
    integer(ip)               :: pgaus
    real(rp)                  :: elcod(ndime,mnode)
    real(rp)                  :: coloc(3),deriv(ndime,mnode),shapf(mnode)
    real(rp)                  :: dummr(2),coord_kp1(3),xxd(3)
    real(rp)                  :: t,dtc,hleng,venor,xfact
    real(rp)                  :: grafo,buofo,Re,dtg,mass
    real(rp)                  :: accel_kp1(3),veloc_kp1(3)                ! a^{k+1},u^{k+1},x^{k+1}
    real(rp)                  :: vefl1(3),vefl2(3)                        ! Fluid velocity at n,n^+1
    real(rp)                  :: v_fluid_k(3),dt_k
    real(rp)                  :: v_fluid_km1(3)
    real(rp)                  :: v_fluid_km2(3)
    real(rp)                  :: T_fluid_k

    integer(ip)               :: iwall,istat
    real(rp)                  :: xinte(3),dista

    real(rp)                  :: denfl,visfl,diame,Du,Cd                  ! Drag
    real(rp)                  :: turvi,gpnut(mgaus)                       ! Turbulent diffusion
    real(rp)                  :: spher
    real(rp)                  :: denpa,Cc,tau_p

    real(rp)                  :: eps(3),D,r                               ! Random walk

    integer(ip)               :: itint,itmax                              ! Time integration
    real(rp)                  :: alpha_str,h
    logical(lg)               :: accept_time_step
    logical(lg)               :: inscont

    real(rp)                  :: tau,Stk(2)                               ! Stk(1): instantaneous Stk, Stk(2): effective Stk
    real(rp)                  :: K                                        ! Saffman constant coefficient
    real(rp)                  :: urel(3),vorti_fl(3)                      ! relative velocity, voriticity
    real(rp)                  :: dMp                                      ! Mass change of particle 
    real(rp)                  :: conv_heat_flux                           ! Heat flux between dropet and fluid
    real(rp)                  :: H_vapor_chm                              ! Enthalpy of vaporized liquid at enthalpy of fuel in lookup table 
    real(rp)                  :: mass_kp1, tempe_kp1 
    real(rp)                  :: confl,sphfl,Dvg_m,Yv_surf,Yv_fluid_k
    real(rp)                  :: Therm_fluid_k,conce_seen(nclas_pts)
    real(rp)                  :: Pr_m, Sc_m, LK_m, xvap, w_nonf_seen
    logical(lg)               :: newmark
    !
    ! Update solution at k+1
    !
    ! k-2          k-1            k            k+1
    !  o-------------o-------------o-------------o------>
    !      dt^k-2         dt^k-1         dt^k
    !            
    itmax     = 200                                          ! Maximum number of time steps
    iwall     = 0                                            ! Do not touch wall
    inscont   = .true.                                       ! Continue...
    newmark   = .true.                                       ! Flag to set if newmark did not converge
    itint     = 0                                            ! Total number time step number (including non-accepted ones)

    itype     = lagrtyp(ilagr) % itype                       ! Particle type

    t         = lagrtyp(ilagr) % t                           ! Particle time
    dt_k      = lagrtyp(ilagr) % dt_k                        ! Particle time step guess: t^k+1 - t^k
    alpha_str = lagrtyp(ilagr) % stret                       ! Stretching

    if( parttyp(itype) % kfl_therm == 0 ) then
       denpa     = parttyp(itype) % denpa                    ! Particle density
    else
       call physics_set_liquid_temperature( parttyp(itype) % liq , lagrtyp(ilagr) % tempe_k)
       denpa     = parttyp(itype) % liq % rho                ! Particle density 
    endif
    spher     = parttyp(itype) % spher                       ! Particle sphericity
    mass      = pts_particle_mass(itype,ilagr)               ! Call FUNCTION to calculate particle mass 
    diame     = pts_particle_diameter(itype,ilagr,denpa)     ! Call FUNCTION to calculate particle diameter
    r         = 0.5_rp * diame                               ! Particle radius
    Cd        = 0.0_rp                                       ! Drag Coefficient inicalization
    Re        = 0.0_rp                                       ! Reynold's Particle initialization
    K         = 2.594_rp                                     ! Constant Saffman coefficient
    eps       = 0.0_rp                                       ! Brownian motion
    urel      = 0.0_rp                                       ! Relative velocity
    Stk       = 0.0_rp                                       ! Stokes number
    dista     = 0.0_rp                                       ! Distance to wall
    Cc        = physics_Cunningham(mean_free_path_pts,diame) ! Cunningham slip correction factor (1 if mean free path < 0 )
    veloc_kp1 = 0.0_rp                                       ! Particle velocity
    coord_kp1 = 0.0_rp                                       ! Particle coordinates
    v_fluid_k = 0.0_rp                                       ! Fluid velocity
    T_fluid_k = 0.0_rp                                       ! Fluid temperature

    !-------------------------------------------------------------------
    !
    ! Loop over time
    !
    ! ti               t <--dtk-->    tf
    ! o----------------|----------|---o
    ! n                              n+1
    ! <------------------------------->
    !              dtime
    !
    ! INSCONT is set to .false. in the following conditions:
    ! 1. Particle has zero time step
    ! 2. Particle go out of computational domain (deposited or outflow)
    ! 3. Particle should migrate to a neighboring subdomain
    !
    !-------------------------------------------------------------------

    do while( t < tf-dtmin_pts .and. inscont )

       ielem = lagrtyp(ilagr) % ielem
       hleng = hleng_pts(ielem)
       !
       ! Modify time step
       !
       if( parttyp(itype) % kfl_tstep < 0 ) then
          dt_k = max(dtmin_pts,min(dt_k*alpha_str,dtime_pts)) ! Prescribed time step (dt_k could have been lower to find first neighbor element)
       else
          dt_k = max(dtmin_pts,dt_k*alpha_str)                ! Adaptive time step (increase or decrease according to ALPHA_STR)
       end if
       dt_k   = min(dtime,dt_k)                               ! Time step cannot be higher than global time step
       dtg    = dt_k
       !
       ! Time and predicted time step
       !
       t      = lagrtyp(ilagr) % t
       !
       ! Synchronize with tf when dt is too large
       !
       if( t + dt_k >= tf-zeror ) dt_k = tf - t
       !
       ! Time step is too small, consider we have arrived
       !
       if( dt_k < zeror ) then
          inscont = .false.
          t       = tf
          goto 20
       end if
       !
       ! Advance in time
       !
       t     = t + dt_k
       itint = itint + 1
       !
       ! Particle is in element IELEM. Get value of shape function SHAPF in it
       ! Compute element length HLENG
       !
       pelty = ltype(ielem)
       pnode = nnode(pelty)
       pgaus = ngaus(pelty)
       do inode = 1,pnode
          ipoin = lnods(inode,ielem)
          elcod(1:ndime,inode) = coord(1:ndime,ipoin)
       end do
       call elmgeo_natural_coordinates(             &
            ndime,pelty,pnode,elcod,shapf,deriv,    &
            lagrtyp(ilagr) % coord,coloc,ifoun,toler)
       !
       ! A particle did not find its element! This is very strange
       !
       if( ifoun == 0 ) then
          lagrtyp(ilagr) % kfl_exist = -3
          inscont = .false.
          t       = tf
          print *, 'particle did not find its element a=',lagrtyp(ilagr) % ilagr
          print *, 'particle did not find its element b=',lagrtyp(ilagr) % coord(1:ndime)
          print *, 'particle did not find its element c=',leinv_loc(ielem),pelty,pnode
          goto 20
       end if
       !
       ! Interpolate fluid velocity uf(t^k,x^k) at old particle position
       !
       ! VEFL1 = u^n+1 at tf
       ! VEFL2 = u^n   at ti
       !
       ! ti              t^k            tf
       ! o----------------o-------------o
       ! VEFL2           x^k            VEFL1
       !
       !
       vefl1 = 0.0_rp
       vefl2 = 0.0_rp
       if( associated(advec) ) then
          do inode = 1,pnode
             ipoin = lnods(inode,ielem)
             vefl1(1:ndime) = vefl1(1:ndime) + shapf(inode) * advec(1:ndime,ipoin,1)
             vefl2(1:ndime) = vefl2(1:ndime) + shapf(inode) * advec(1:ndime,ipoin,3)
          end do
       end if
       if( associated(vesgs) .and. kfl_vesgs_pts /= 0 ) then
          call elmgeo_gauss_to_element(pnode,pgaus,shapf,elmar(pelty) % shaga,ndime,vesgs(ielem) % a(:,:,1),vefl1,INITIALIZE=.false.)
          call elmgeo_gauss_to_element(pnode,pgaus,shapf,elmar(pelty) % shaga,ndime,vesgs(ielem) % a(:,:,2),vefl2,INITIALIZE=.false.)
       end if
       xfact = dtinv * ( lagrtyp(ilagr) % t - ti )
       lagrtyp(ilagr) % v_fluid_k(1:ndime) = (1.0_rp-xfact) * vefl2(1:ndime) + xfact * vefl1(1:ndime)
       !
       !  Fluid velocity at previous times, a priori unknown at k+1
       !
       !   n     k-2     k-1     k     k+1     n+1
       !   o------|-------|------|------|-------o
       !   ti      dt_km2  dt_km1  dt_k         tf
       !
       !   uf^k   = uf(x^k,t^k)     = v_fluid_k(1:3)   = lagrtyp(ilagr) % v_fluid_k(1:3)
       !   uf^k-1 = uf(x^k-1,t^k-1) = v_fluid_km1(1:3) = lagrtyp(ilagr) % v_fluid_km1(1:3)
       !   uf^k-2 = uf(x^k-2,t^k-2) = v_fluid_km2(1:3) = lagrtyp(ilagr) % v_fluid_km2(1:3)
       !
       !   dt^k   = t^k+1 - t^k
       !   dt^k-1 = t^k   - t^k-1
       !   dt^k-2 = t^k-1 - t^k-2
       !
       v_fluid_k   = lagrtyp(ilagr) % v_fluid_k
       v_fluid_km1 = lagrtyp(ilagr) % v_fluid_km1
       v_fluid_km2 = lagrtyp(ilagr) % v_fluid_km2
       !
       ! Fluid properties
       !
       if (parttyp(itype) % kfl_therm /= 0) then
          call pts_thermodynamic_properties(itype, lagrtyp(ilagr) % tempe_k, ielem, pnode, lnods(1:pnode,ielem), shapf, &
               confl, sphfl, denfl, visfl, Dvg_m, Pr_m, Sc_m, LK_m, Yv_surf, Yv_fluid_k, &
               Therm_fluid_k, T_fluid_k, xvap, w_nonf_seen, conce_seen)
       else
          call ker_proper('VISCO','IGAUS',dumm0,ielem,visfl,pnode,1_ip,shapf)   ! Fluid viscosity
          call ker_proper('DENSI','IGAUS',dumm0,ielem,denfl,pnode,1_ip,shapf)   ! Fluid density
       endif
       !
       ! Turbulent diffusion
       !
       if (parttyp(itype) % kfl_turbu /= 0) then
          call ker_proper('TURBU','PGAUS',dumm0,ielem,gpnut,pnode,pgaus,elmar(pelty) % shape) ! nut
          call elmgeo_gauss_to_element(pnode,pgaus,shapf,elmar(pelty) % shaga,gpnut,turvi)
          visfl = visfl + denfl * turvi / parttyp(itype) % tursc 
       end if
       !
       ! Diffusion coefficent values
       !
       if( parttyp(itype) % kfl_brown /= 0 ) then
          D = parttyp(itype) % diffu
          call pts_brown(&
               & itint,ielem,itype,pnode,shapf,&
               & dt_k,r,visfl,denfl,Cc,D,eps(1),eps(2),eps(3))
       else
          D = 0.0_rp
       end if

       !----------------------------------------------------------
       !
       !     ----------------------------------------------
       !     UPDATE ACCELERATION, VELOCITY, POSITION AT K+1
       !     ----------------------------------------------
       !
       !----------------------------------------------------------

       if( kfl_exacs_pts /= 0 ) then
          !
          ! Exact acceleration
          !
          call pts_exacso(3_ip,t,accel_kp1,dummr,dummr)
          veloc_kp1(1:ndime) = lagrtyp(ilagr) % veloc(1:ndime) &
               + dt_k * ( (1.0_rp-gamma_pts) * accel_kp1(1:ndime) + gamma_pts * lagrtyp(ilagr) % accel(1:ndime))
          coord_kp1(1:ndime) = lagrtyp(ilagr) % coord(1:ndime) + dt_k * lagrtyp(ilagr) % veloc(1:ndime) &
               + dt_k * dt_k * ( 0.5_rp*lagrtyp(ilagr) % accel(1:ndime) &
               + beta_pts *(accel_kp1(1:ndime)-lagrtyp(ilagr) % accel(1:ndime)) )

       else if( parttyp(itype) % kfl_modla == 1 ) then
          !
          ! Pure transport
          !
          if( kfl_order_pts < 0 .and. ( pelty == TET04 .or. pelty == TRI03 ) ) then
             !
             ! Analytical integration
             !
             call pts_analytics(ielem, dt_k, lagrtyp(ilagr) % coord, coord_kp1)

          else
             !
             ! Adams-Bashforth
             !
             call pts_transport_Adams_Bashforth(&
                  lagrtyp(ilagr),dt_k,coord_kp1,veloc_kp1,&
                  accel_kp1)
          end if

       else if( parttyp(itype) % kfl_modla == 2) then
          !
          ! Force model
          !
          grafo     = real( parttyp(itype) % kfl_grafo, rp ) ! Gravity  force = 1.0
          buofo     = real( parttyp(itype) % kfl_buofo, rp ) ! Buoyancy force = 1.0
          veloc_kp1 = lagrtyp(ilagr) % veloc                 ! Initial guess
          !
          ! Vorticity useful for Saffman Mei force
          !
          if( parttyp(itype) % kfl_saffm /= 0 ) then
             vorti_fl(1:ndime) = 0.0_rp
             do inode = 1,pnode
                ipoin = lnods(inode,ielem)
                vorti_fl(1:ndime) = vorti_fl(1:ndime) + shapf(inode) * vorti(1:ndime,ipoin)
             end do
          end if
          !
          ! Distance to wall
          !
          if(parttyp(itype) % kfl_drafo == 8 ) then
             dista = dot_product(shapf(1:pnode),walld(lnods(1:pnode,ielem)))
          else
             dista = huge(1.0_rp)
          end if

          if(      parttyp(itype) % kfl_schem == 2 ) then
             !
             ! Runge-Kutta 4th-order
             !
             call pts_rungk4(                                                           &
                  ndime,parttyp(itype) % kfl_drafo,parttyp(itype) % kfl_extfo,          &
                  grafo,buofo,g,v_fluid_k,visfl,denfl,veloc_kp1,lagrtyp(ilagr) % coord, &
                  denpa,diame,spher,dista,t,dt_k,veloc_kp1,coord_kp1)

          else if(  parttyp(itype) % kfl_schem == 3 ) then
             !
             ! Runge-Kutta 4th/5th-order with time step control
             !
             call pts_rungk45(                                                          &
                  ndime,parttyp(itype) % kfl_drafo,parttyp(itype) % kfl_extfo,          &
                  grafo,buofo,g,v_fluid_k,visfl,denfl,veloc_kp1,lagrtyp(ilagr) % coord, &
                  denpa,diame,spher,dista,t,dt_k,veloc_kp1,xxd)

          else if( parttyp(itype) % kfl_schem == 1 ) then
             !
             ! Analytical time integration
             !
             call runend('PTS_TRANSPORT_SINGLE_PARTICLE: NOT CODED')

          else if( parttyp(itype) % kfl_schem == 0 ) then
             !
             ! Numerical time integration: Newmark + Newton-Raphson
             !
             call pts_transport_newmark_newton_raphson(&
                  lagrtyp(ilagr),t,dt_k,Cc,visfl,denfl,vorti_fl,&
                  eps,dista,diame,denpa,spher,grafo,buofo,&
                  coord_kp1,veloc_kp1,accel_kp1,Cd,Re,tau_p,&
                  Stk,istat)
             !
             ! Update new position or go back in time if Newton-Raphson not converged
             !
             if( istat /= 0 ) then
                alpha_str         = strex
                t                 = lagrtyp(ilagr) % t
                accept_time_step  = .false.
                newmark           = .false.  !temporarily set to throw an error in 20              
                goto 20
             end if
          end if
          call pts_stk_local(ielem,tau_p,lagrtyp(ilagr) % t_inject,Stk(1),Stk(2))

       end if

       !----------------------------------------------------------
       !
       !                  -----------
       !                  RANDOM WALK
       !                  -----------
       !
       !----------------------------------------------------------

       if( parttyp(itype) % kfl_brown == 1 ) then
          coord_kp1 = coord_kp1 + eps
       end if

       !----------------------------------------------------------
       !
       !                  ------------------
       !                  THERMONDYMIC MODEL
       !                  ------------------
       !
       ! Compute particle temperature and mass
       !
       !----------------------------------------------------------

       if( kfl_thermo_pts /= 0 ) then
          call pts_thermodynamic_transport(              &
               ilagr,lagrtyp(ilagr),Re,ielem,pnode,      &
               lnods(1:pnode,ielem),shapf,dt_k,mass_kp1, &
               tempe_kp1,conv_heat_flux)
       end if

       !----------------------------------------------------------
       !
       !                    -----------------------
       !                    SLIP AND BOUNCING WALLS
       !                    -----------------------
       !
       !----------------------------------------------------------
       call pts_transport_special_walls(&
            lagrtyp(ilagr),ielem,pnode,shapf,deriv,elcod, &
            diame,coord_kp1,veloc_kp1,accel_kp1)

       !----------------------------------------------------------
       !
       !                    ---------------
       !                    WALL DEPOSITION
       !                    ---------------
       !
       ! Check if trajectory crosses the wall
       ! IWALL = 0 ... stay in domain, slip, bouncing
       !       = 1 ... outflow, wall
       !
       !----------------------------------------------------------

       iwall = 0
       if( lboue_pts(ielem) > 0 .or. kfl_walld_pts > 0 ) then
          call pts_cross_wall(&
               ielem,pnode,ilagr,diame,hleng,toler,shapf,deriv,iwall,&
               iboun,coord_kp1,veloc_kp1,accel_kp1,xinte,t,dt_k)
          if( iwall == 1 ) then
             !
             ! Hit a wall: go out of time loop
             !
             inscont                         = .false.
             lagrtyp(ilagr) % iboun          = iboun
             lagrtyp(ilagr) % coord(1:ndime) = xinte(1:ndime)
             lagrtyp(ilagr) % t              = t
             goto 20

          end if
       end if

       !----------------------------------------------------------
       !
       !                  ------------------
       !                  ADAPTIVE TIME STEP
       !                  ------------------
       !
       !----------------------------------------------------------

       accept_time_step = .true.
       if( parttyp(itype) % kfl_tstep > 0 .and. parttyp(itype) % kfl_modla == 2 ) then
          call pts_transport_length_scale(&
               itype,diame,denpa,visfl,denfl,veloc_kp1,Du,hleng,h)
          call pts_adapti(&
               parttyp(itype) % kfl_tstep,h,coord_kp1,lagrtyp(ilagr)%coord,veloc_kp1,lagrtyp(ilagr)%veloc,&
               v_fluid_k,accel_kp1,lagrtyp(ilagr)%accel,tau,dt_k,lagrtyp(ilagr)%dt_km1,parttyp(itype) % safet,&
               alpha_str,accept_time_step,parttyp(itype) % kfl_modla)
          if( .not. accept_time_step ) then
             t = lagrtyp(ilagr) % t
             goto 20
          end if
       else
          alpha_str = ovstr
       end if
       !
       ! Critical time step based on element minimum length HLENG
       !
       venor = sqrt(dot_product(veloc_kp1(1:ndime),veloc_kp1(1:ndime)))
       if( venor /= 0.0_rp .or. D /= 0.0_rp ) then
          dtc = 1.0_rp / ( venor / hleng + 2.0_rp * D / hleng**2 )
       else
          dtc = 1.0e12_rp
       end if
       
       !----------------------------------------------------------
       !
       !                      -------------------
       !                      SEARCH HOST ELEMENT
       !                      -------------------
       !
       !----------------------------------------------------------
       !
       ! Start by checking in IELEM... Is it really likely?
       !
       call pts_transport_element_search(ielem,coord_kp1,ifoun,jelem,shapf)
       if( ifoun /= 0 ) then          
          !
          ! Residence time
          !
          if( kfl_resid_pts > 0 ) then
             if( jelem == ielem ) then
                !$OMP ATOMIC
                resid_pts(itype,ielem) = resid_pts(itype,ielem) + dt_k
             else
                !
                ! Intersection between segment and face should be computed in mod_elmgeo
                ! If JELEM is a ghost element, info will be further sent in pts_endite()
                !
                !$OMP ATOMIC
                resid_pts(itype,ielem) = resid_pts(itype,ielem) + 0.5_rp * dt_k
                !$OMP ATOMIC
                resid_pts(itype,jelem) = resid_pts(itype,jelem) + 0.5_rp * dt_k
             end if
          end if

          !----------------------------------------------------------
          !
          !                      -------------------
          !                      UPDATE PROPERTIES
          !                      -------------------
          !
          !----------------------------------------------------------


          if( jelem > nelem ) then
             !
             ! JELEM > NELEM: particle goes to neighboring subdomain INEIG
             !
             ineig                      = leldo(1,jelem-nelem)            ! Which neighbor holds this element
             lagrtyp(ilagr) % kfl_exist = ineig                           ! Neighbor
             inscont                    = .false.                         ! Go out of internal time loop
             nlagr_migrating(ineig)     = nlagr_migrating(ineig) + 1      ! Number of particles to send to neighbor INEIG
             lagrtyp(ilagr) % ielem     = leldo(2,jelem-nelem)            ! Local numbering of JELEM in my neighbor
          else
             !
             ! Element found: accept time step update new position
             !
             lagrtyp(ilagr) % ielem     = jelem                           ! JELEM remains in my subdomain
          end if

          lagrtyp(ilagr) % ittim        = lagrtyp(ilagr) % ittim + 1      ! One new time step
          lagrtyp(ilagr) % Cd           = Cd                              ! Accept drag
          lagrtyp(ilagr) % Re           = Re                              ! Accept Reynolds
          lagrtyp(ilagr) % Stk          = Stk                             ! Stokes number ! (Stk(1): instantaneous, Stk(2): effective)
          lagrtyp(ilagr) % stret        = alpha_str                       ! Stretching
          lagrtyp(ilagr) % t            = t                               ! Accept new time
          lagrtyp(ilagr) % dt_k         = dtg                             ! Guess for next time step k
          lagrtyp(ilagr) % dt_km2       = lagrtyp(ilagr) % dt_km1         ! Save time step k-2
          lagrtyp(ilagr) % dt_km1       = dt_k                            ! Save time step k-1
          lagrtyp(ilagr) % coord_km1    = lagrtyp(ilagr) % coord_k        ! Previous coordinates k-1
          lagrtyp(ilagr) % coord_k      = lagrtyp(ilagr) % coord          ! Previous coordinates k
          lagrtyp(ilagr) % coord        = coord_kp1                       ! Accept position k+1
          lagrtyp(ilagr) % veloc        = veloc_kp1                       ! Accept velocity k+1
          lagrtyp(ilagr) % accel        = accel_kp1                       ! Accept acceleration k+1

          if ( parttyp(itype) % kfl_therm /= 0) then
             lagrtyp(ilagr) % tempe_km1    = lagrtyp(ilagr) % tempe_k    ! Temperature at k+1
             lagrtyp(ilagr) % mass_km1     = lagrtyp(ilagr) % mass_k     ! Mass at k+1
             lagrtyp(ilagr) % tempe_k      = tempe_kp1                   ! Temperature at k+1
             lagrtyp(ilagr) % mass_k       = mass_kp1                    ! Mass at k+1
             !
             ! Recalculate postprocessing quantities
             !
             call physics_set_liquid_temperature( parttyp(itype) % liq , lagrtyp(ilagr) % tempe_k)
             lagrtyp(ilagr) % diam_k = pts_particle_diameter(itype,ilagr,parttyp(itype) % liq % rho)   ! Diameter at k+1
          endif
          lagrtyp(ilagr) % v_fluid_km2  = lagrtyp(ilagr) % v_fluid_km1    ! Previous fluid velocity k-2
          lagrtyp(ilagr) % v_fluid_km1  = lagrtyp(ilagr) % v_fluid_k      ! Previous fluid velocity k-1
          !
          ! Compute derivated parameters              
          !
          call pts_transport_derivated_parameters(lagrtyp(ilagr))




          !----------------------------------------------------------
          !
          !                      -------------------
          !                      COUPLING WITH CFD
          !                      -------------------
          !
          !----------------------------------------------------------

          !  
          ! Compute evaporation rate for coupling
          ! dm_p = m^k+1 - m^k
          !  
          if( kfl_thermo_pts == 0 ) then
             dMp = 0.0_rp
          else
             dMp = lagrtyp(ilagr) % mass_k - lagrtyp(ilagr) % mass_km1
          endif

          !
          ! Momentum exchange for coupling with CFD
          !
          if( kfl_momentum_sink_pts /= 0 ) then 

             peltj = ltype(jelem)
             pnodj = nnode(peltj)
             do inode = 1,pnodj
                ipoin = lnods(inode,jelem)
                if ((ipoin >= 1)) then
                   do idime = 1,ndime
                      !$OMP ATOMIC
                      momentum_sink(idime,ipoin) = momentum_sink(idime,ipoin)     &
                           &  - shapf(inode) / dtime * parttyp(itype) % n_drop *  & 
                           & ( mass *  lagrtyp(ilagr) % acced(idime) * dt_k +     &
                           &   lagrtyp(ilagr) % veloc(idime) * dMp )
                   end do
                else
                   print'(A,I6,A,I6,A,I6,A,I6,A,I6)','Element is problematic. jelem/nelem: ',jelem,'/',nelem,', lnods(inode,jelem)/npoin: ', lnods(inode,jelem),'/',npoin, ', ilagr: ', ilagr
                endif
             end do
             !print'(A,10(1X,E16.8))', 'shapf', shapf
             !print'(A,10(1X,E16.8))', 'dtime', dtime
             !print'(A,10(1X,E16.8))', 'mass', mass
             !print'(A,10(1X,E16.8))', 'lagrtyp(ilagr) % acced', lagrtyp(ilagr) % acced
             !print'(A,10(1X,E16.8))', 'dt_k', dt_k
             !print'(A,10(1X,E16.8))', 'lagrtyp(ilagr) % veloc', lagrtyp(ilagr) % veloc
             !print'(A,10(1X,E16.8))', 'dMp', dMp
          end if

          !
          ! Mass exchange for coupling with CFD
          !
          if( kfl_mass_sink_pts /= 0 ) then 

             peltj = ltype(jelem)
             pnodj = nnode(peltj)
             do inode = 1,pnodj
                ipoin = lnods(inode,jelem)
                if ((ipoin >= 1)) then
                   !$OMP ATOMIC
                   mass_sink(ipoin) = mass_sink(ipoin) - dMp * parttyp(itype) % n_drop / dtime * shapf(inode) 
                endif
             end do
             !print'(A,10(1X,E16.8, A))', 'dMp', dMp, ' = ', lagrtyp(ilagr) % mass_k, ' - ', lagrtyp(ilagr) % mass_km1
             !print'(A,10(1X,E16.8))', 'dtime', dtime
             !print'(A,10(1X,E16.8))', 'shapf', shapf
          end if

          !
          ! Heat exchange for coupling with CFD
          !
          if( kfl_heat_sink_pts /= 0 ) then

             if ( parttyp(itype) % kfl_therm == 2) then
                !  
                ! Compute enthalpy of vaporized mixture:  H_vapor_chm
                !  
                call physics_T_2_HCp(lagrtyp(ilagr) % tempe_k, parttyp(itype) % cpcoef_v_chm, H_vapor_chm, dummr(1))
             else   

                H_vapor_chm = parttyp(itype) % L_vapor
                call runend('In mod_pts_transport: Enthalpy of vaporized mixture to be computed for non-tabulated models')
             end if

             peltj = ltype(jelem)
             pnodj = nnode(peltj)
             do inode = 1,pnodj
                !             dE_drop     d(m cp T)     
                ! S_entha = - ------- = - ---------         
                !               dt          dt                   
                ipoin = lnods(inode,jelem)
                if ((ipoin >= 1)) then
                   !$OMP ATOMIC
                   heat_sink(ipoin) = heat_sink(ipoin) - ( conv_heat_flux * dt_k  + dMp * H_vapor_chm ) / dtime * shapf(inode) * parttyp(itype) % n_drop
                endif
             end do

             !print'(A,10(1X,E16.8))', 'conv_heat_flux', conv_heat_flux
             !print'(A,10(1X,E16.8))', 'dt_k', dt_k
             !print'(A,10(1X,E16.8))', 'dMp', dMp
             !print'(A,10(1X,E16.8))', 'H_vapor_chm', H_vapor_chm
             !print'(A,10(1X,E16.8))', 'dtime', dtime
             !print'(A,10(1X,E16.8))', 'shapf', shapf
          end if

       else
          !
          ! Element not found: reduce time step
          !
          t                             = lagrtyp(ilagr) % t              ! Go back to previous time
          alpha_str                     = strex                           ! Decrease time step
          accept_time_step              = .false.                         ! Do not accept time step
       end if
       !
       ! Very small dt: remove particle
       !
20     continue
       if( ( dt_k < 0.1_rp * dtmin_pts .or. itint > itmax ) .and. iwall /= 1 ) then
          inscont = .false.                     ! Go out of time loop
          if( lelbf(ielem) % n /= 0 ) then
             !lagrtyp(ilagr) % kfl_exist = PTS_PARTICLE_HITS_WALL  ! Assume it hits the wall
             lagrtyp(ilagr) % kfl_exist = PTS_PARTICLE_ZERO_TIME   ! It is lost
             print*,'A----la pierdo',lagrtyp(ilagr) % ilagr, lagrtyp(ilagr) % itype
             print*,'A----con coord-kk',coord_kp1
          else
             lagrtyp(ilagr) % kfl_exist = PTS_PARTICLE_ZERO_TIME   ! Assume it is lost for numerical reasons... :o(
            print*,'B----la pierdo',lagrtyp(ilagr) % ilagr, lagrtyp(ilagr) % itype
            print*,'B----con coord-kk',coord_kp1
        end if
          if( .not. newmark ) then                 
             write(*,'(A,I8,A,3(1X,E18.8),A,1X,E18.8)') 'Newmark did not converge after all iteraitons. Particle id=', lagrtyp(ilagr) % ilagr,'; coords=',lagrtyp(ilagr) % coord(1:3),'; diameter=',diame 
          end if
       end if

       newmark = .true.

    end do
    !
    ! Particle stay in current subdomain...
    !
    if( lagrtyp(ilagr) % kfl_exist < 0 ) nlagr_non_migrating_pts = nlagr_non_migrating_pts + 1

  end subroutine pts_transport_single_particle

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-09-25
  !> @brief   Finalize the transport
  !> @details Finalize the transport of particles and
  !>          1. Compute some useful statistics:
  !>             NLAGR_GOING_OUT_PTS = # of wall or outflow particles
  !>             NLAGR_ZERO_TIME_PTS = # particles with zero time step
  !>             NLAGR_DEPOSITED_PTS = # particles deposited without boundary
  !>             NLAGR_HITS_WALL_PTS = # particles deposited on a wall
  !>          2. Detect events
  !>          3. Accumulate deposited particle in DEPOB_PTS
  !> 
  !-----------------------------------------------------------------------

  subroutine pts_transport_finalize()

    integer(ip)                  :: ilagr,ielem,iboun,klagr,dummi(5),itype
    integer(ip)                  :: idime,ipoin
    integer(4)                   :: my_nlagr_4,ipart4,ndime4,my_rank4
    real(rp)                     :: time_max,time_ave,load_balance
    real(rp)                     :: mass_in_dom, enth_in_dom, massp, enthp, dummr
    real(rp)                     :: mom_in_dom(ndime), kine_in_dom, velop(ndime)
    real(rp)                     :: bound_box_min(ndime),bound_box_max(ndime) 
    integer(ip)                  :: nlagr_evaporated_pts
    integer(4),          pointer :: par_nlagr_4(:)
    real(rp),            pointer :: coord_nlagr_4(:,:)
    real(rp),            pointer :: my_coord_nlagr_4(:,:)
    real(rp),            pointer :: my_veloc_nlagr_4(:,:)
    real(rp),            pointer :: auxvar(:)
    real(rp),            pointer :: auxelm(:)
    real(rp)                     :: diame 

    nullify(par_nlagr_4)
    nullify(coord_nlagr_4)
    nullify(my_coord_nlagr_4)
    nullify(my_veloc_nlagr_4)

    !----------------------------------------------------------------------
    !
    ! Count deposited particles
    !
    !----------------------------------------------------------------------

    nlagr_going_out_pts = 0
    nlagr_zero_time_pts = 0
    nlagr_deposited_pts = 0
    nlagr_hits_wall_pts = 0
    nlagr_evaporated_pts= 0
    
    if( INOTMASTER ) then
       do ilagr = 1,mlagr
          !do ilagr_local = 1,nlagr_local_pts
          !ilagr = permu_nlagr_pts(ilagr_local)
          if( lagrtyp(ilagr) % kfl_exist == PTS_PARTICLE_HITS_WALL ) nlagr_hits_wall_pts = nlagr_hits_wall_pts + 1

          if( lagrtyp(ilagr) % kfl_exist == PTS_PARTICLE_HITS_WALL   .or. &
              lagrtyp(ilagr) % kfl_exist == PTS_PARTICLE_OUTFLOW     .or. &
              lagrtyp(ilagr) % kfl_exist == PTS_PARTICLE_EVAPORATED  .or. &
              lagrtyp(ilagr) % kfl_exist == PTS_PARTICLE_MOVING_MESH ) then
             !
             ! Deposited, outflow, or fully evaporated particle
             !
             nlagr_going_out_pts    = nlagr_going_out_pts + 1
             if( lagrtyp(ilagr) % boundary_set == -1 ) nlagr_deposited_pts = nlagr_deposited_pts + 1
             if (lagrtyp(ilagr) % kfl_exist == PTS_PARTICLE_EVAPORATED) nlagr_evaporated_pts = nlagr_evaporated_pts + 1

          else if( lagrtyp(ilagr) % kfl_exist == PTS_PARTICLE_ZERO_TIME ) then
             !
             ! Zero time step
             !
             nlagr_zero_time_pts = nlagr_zero_time_pts + 1

          end if
       end do
    end if
    
    !
    ! Gather thermodynamic properties
    !
    bound_box_min =  1.0e15
    bound_box_max = -1.0e15
    mass_in_dom   =  0.0_rp
    enth_in_dom   =  0.0_rp
    mom_in_dom    =  0.0_rp
    kine_in_dom   =  0.0_rp
    if( INOTEMPTY .and. kfl_thermo_pts /= 0 ) then
       
       if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVMIE') ) then
          nullify(auxelm)
          call memory_alloca(mem_modul(1:2,modul),'AUXELM','pts_transport_finalize',auxelm,nelem)
       endif
       do ilagr = 1,mlagr
          if( lagrtyp(ilagr) % kfl_exist == PTS_PARTICLE_EXISTS ) then
             itype = lagrtyp(ilagr) % itype 
             !
             ! Update liquid properties
             !
             call physics_set_liquid_temperature( parttyp(itype) % liq, lagrtyp(ilagr) % tempe_k ) 
             
             !
             ! Mass 
             !
             massp       = pts_particle_mass( lagrtyp(ilagr) % itype, ilagr) 

             !
             ! Enthalpy
             !
             if ( parttyp(itype) % kfl_therm == 2) then
                !  
                ! Compute enthalpy of vaporized mixture
                !  
                call physics_T_2_HCp(lagrtyp(ilagr) % tempe_k, parttyp(itype) % cpcoef_v_chm, enthp, dummr)
             else   

                enthp = parttyp(itype) % liq % cp * lagrtyp(ilagr) % tempe_k
             end if 
             enthp = enthp - parttyp(itype) % liq % Lv

             !
             ! Momentum and kinetic energy
             !
             velop(1:ndime) = lagrtyp(ilagr) % veloc(1:ndime)
             
             
             !
             ! Calculate total quantities subdomain
             !
             mass_in_dom = mass_in_dom + massp * parttyp(itype) % n_drop
             enth_in_dom = enth_in_dom + massp * enthp * parttyp(itype) % n_drop
             mom_in_dom  = mom_in_dom  + massp * velop * parttyp(itype) % n_drop
             kine_in_dom = kine_in_dom + 0.5_rp * massp * dot_product(velop,velop) * parttyp(itype) % n_drop
             
             !
             ! Calculate bounding box
             !
             do idime = 1,ndime
                bound_box_min(idime) = min(bound_box_min(idime),lagrtyp(ilagr) % coord(idime)) 
                bound_box_max(idime) = max(bound_box_max(idime),lagrtyp(ilagr) % coord(idime)) 
             enddo


             !
             ! Acummulate average Mie scattering if needed
             !
             if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVMIE') ) then
                ielem = lagrtyp(ilagr) % ielem
                diame = pts_particle_diameter(lagrtyp(ilagr) % itype,ilagr,parttyp(lagrtyp(ilagr) % itype) % liq % rho)
                if( ielem > 0 .and. ielem <= nelem  ) then
                   auxelm(ielem) = auxelm(ielem) +  parttyp(lagrtyp(ilagr) % itype) % n_drop * diame**2 
                end if
             end if

          end if
       end do

       !
       ! Summ avarage Mie 
       !
       if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVMIE') ) then
          nullify(auxvar)
          call memory_alloca(mem_modul(1:2,modul),'AUXVAR','pts_transport_finalize',auxvar,npoin)
          call projec_elements_to_nodes(auxelm,auxvar)

          call solver_lumped_mass_system(1_ip,auxvar,EXCHANGE=.false.)

          do ipoin = 1,npoin
             avg_mie_pts(ipoin) = avg_mie_pts(ipoin) + auxvar(ipoin)
          enddo
          call memory_deallo(mem_modul(1:2,modul),'AUXVAR','pts_transport_finalize',auxvar)
          call memory_deallo(mem_modul(1:2,modul),'AUXELM','pts_transport_finalize',auxelm)
       endif

    endif

    dummi(1) = nlagr_going_out_pts ! Out of computational domain
    dummi(2) = nlagr_zero_time_pts ! With zero time step
    dummi(3) = nlagr_deposited_pts ! Deposited but no boundary was found
    dummi(4) = nlagr_hits_wall_pts ! Hitting a wall
    dummi(5) = nlagr_evaporated_pts! Evaporated

    call PAR_SUM(5_ip,dummi,'IN MY ZONE')

    nlagr_going_out_pts = dummi(1)
    nlagr_zero_time_pts = dummi(2)
    nlagr_deposited_pts = dummi(3)
    nlagr_hits_wall_pts = dummi(4)
    nlagr_evaporated_pts= dummi(5)

    if( nlagr_going_out_pts > 0 ) call livinf(-9_ip,'PARTICLES GOING OUT OF COMPUTATIONAL DOMAIN= ',nlagr_going_out_pts)
    if( nlagr_zero_time_pts > 0 ) call livinf(-9_ip,'PARTICLES WITH ZERO TIME STEP= ',              nlagr_zero_time_pts)
    if( nlagr_deposited_pts > 0 ) call livinf(-9_ip,'PARTICLES DEPOSITED OUT OF BOUNDARY= ',        nlagr_deposited_pts)
    if( nlagr_hits_wall_pts > 0 ) call livinf(-9_ip,'PARTICLES HITTING A WALL= ',                   nlagr_hits_wall_pts)
    if( nlagr_evaporated_pts> 0 ) call livinf(-9_ip,'PARTICLES EVAPORATED=',                        nlagr_evaporated_pts)

    !----------------------------------------------------------------------
    !
    ! Timing
    !
    !----------------------------------------------------------------------

    cpu_modul(CPU_ASSEMBLY,modul) = cpu_modul(CPU_ASSEMBLY,modul) + time_transport_pts
    time_max = time_transport_pts
    time_ave = time_transport_pts

    call PAR_MAX(time_max)
    call PAR_SUM(time_ave)
    call PAR_MAX(particles_sent)
    call PAR_MAX(particles_recv)
    
    !
    ! Thermodynamic statistics:
    !
    if( kfl_thermo_pts /= 0 ) then
       call PAR_SUM(mass_in_dom)
       call PAR_SUM(enth_in_dom)
       call PAR_SUM(kine_in_dom)
       do idime = 1,ndime
          call PAR_SUM(mom_in_dom(idime))
          call PAR_MAX(bound_box_max(idime))
          call PAR_MIN(bound_box_min(idime))
       enddo
    endif


    time_ave     = time_ave / max(1.0_rp,real(PAR_CODE_SIZE-1,rp))
    load_balance = time_ave / max(zeror,time_max)

    !
    ! Output Convergence data
    !
    if( INOTSLAVE ) then
       if( kfl_thermo_pts /= 0 ) then
          !
          ! Output properties relevant for evaporating particles
          !
          write(momod(modul) % lun_conve,112) &
               cutim,               &
               nlacc_pts,           &
               nlagr_existing_pts,  &
               nlagr_going_out_pts, &
               nlagr_zero_time_pts, &
               time_max,            &
               time_ave,            &
               load_balance,        &
               particles_sent,      &
               particles_recv,      &
               comm_loops_pts,      &
               mass_in_dom,         &
               enth_in_dom,         &
               mom_in_dom,          &
               kine_in_dom,         &
               bound_box_min,       &
               bound_box_max
       else
          !
          ! Output usual properties
          !
          write(momod(modul) % lun_conve,111) &
               cutim,               &
               nlacc_pts,           &
               nlagr_existing_pts,  &
               nlagr_going_out_pts, &
               nlagr_zero_time_pts, &
               time_max,            &
               time_ave,            &
               load_balance,        &
               particles_sent,      &
               particles_recv,      &
               comm_loops_pts
       endif
       call iofile_flush_unit(momod(modul) % lun_conve)
    end if


    !----------------------------------------------------------------------
    !
    ! Save particles on deposition map
    !
    !----------------------------------------------------------------------

    if( INOTMASTER ) then
       if( kfl_depos_pts == 1 .or. kfl_depos_surface_pts /= 0) then
          do ilagr = 1,mlagr
             !do ilagr_local = 1,nlagr_local_pts
             !ilagr = permu_nlagr_pts(ilagr_local)
             if( lagrtyp(ilagr) % kfl_exist == PTS_PARTICLE_HITS_WALL ) then
                ielem = lagrtyp(ilagr) % ielem
                iboun = lagrtyp(ilagr) % iboun
                if( iboun > 0 ) then
                   !depoe_pts(lagrtyp(ilagr) % itype,ielem) = depoe_pts(lagrtyp(ilagr) % itype,ielem) + 1.0_rp
                   depob_pts(lagrtyp(ilagr) % itype,iboun) = depob_pts(lagrtyp(ilagr) % itype,iboun) + 1.0_rp
                end if
             end if
          end do
       end if
    end if

    !----------------------------------------------------------------------
    !
    ! Events: particles are lost
    !
    !----------------------------------------------------------------------

    my_nlagr_4 = int(nlagr_zero_time_pts,4)
    if( nlagr_zero_time_pts > 0 .and. kfl_detection /= 0 ) then

       if( IPARALL ) then
          !
          ! All gather particle coordinates
          !
          call memory_alloca(mem_modul(1:2,modul),'PAR_NLAGR_4','pts_solite',par_nlagr_4,int(PAR_CODE_SIZE,4),'INITIALIZE',0_4)
          call PAR_ALLGATHER(my_nlagr_4,par_nlagr_4,1_4,'IN MY CODE')
          if( INOTMASTER ) then
             call memory_alloca(mem_modul(1:2,modul),'MY_COORD_NLAGR_4','pts_solite',my_coord_nlagr_4,ndime,int(my_nlagr_4,ip),'INITIALIZE')
             call memory_alloca(mem_modul(1:2,modul),'MY_VELOC_NLAGR_4','pts_solite',my_veloc_nlagr_4,ndime,int(my_nlagr_4,ip),'INITIALIZE')
             my_nlagr_4 = 0
             do ilagr = 1,mlagr
                !do ilagr_local = 1,nlagr_local_pts
                !ilagr = permu_nlagr_pts(ilagr_local)
                if( lagrtyp(ilagr) % kfl_exist == PTS_PARTICLE_ZERO_TIME ) then
                   !if( lagrtyp(ilagr) % kfl_exist == -1 ) then
                   my_nlagr_4 = my_nlagr_4 + 1
                   my_coord_nlagr_4(1:ndime,my_nlagr_4) = lagrtyp(ilagr) % coord(1:ndime)
                   my_veloc_nlagr_4(1:ndime,my_nlagr_4) = lagrtyp(ilagr) % veloc(1:ndime)
                end if
             end do
          end if
          ndime4 = int(ndime,4)
          do ipart4 = 0,int(PAR_CODE_SIZE,4)-1_4
             par_nlagr_4(ipart4) = par_nlagr_4(ipart4) * ndime4
          end do
          call memory_alloca(mem_modul(1:2,modul),'COORD_NLAGR_4','pts_solite',coord_nlagr_4,ndime,nlagr_zero_time_pts,'INITIALIZE')
          call PAR_ALLGATHERV(my_coord_nlagr_4,coord_nlagr_4,par_nlagr_4,'IN MY CODE')
          !
          ! Output event
          !
          call PAR_COMM_RANK_AND_SIZE(my_rank4,wherein='IN MY ZONE')
          klagr = 0
          do ipart4 = 0,int(PAR_CODE_SIZE,4)-1_4
             par_nlagr_4(ipart4) = par_nlagr_4(ipart4) / ndime4
             do ilagr = 1,par_nlagr_4(ipart4)
                klagr = klagr + 1
                if( ipart4 == my_rank4 ) then
                   call ker_events_particle_not_converged(ipart4,coord_nlagr_4(1:ndime,klagr),ndime,advec,my_veloc_nlagr_4(1:ndime,klagr))
                else
                   call ker_events_particle_not_converged(ipart4,coord_nlagr_4(1:ndime,klagr),ndime,advec)
                end if
             end do
          end do

          call memory_deallo(mem_modul(1:2,modul),'COORD_NLAGR_4'   ,'pts_solite',coord_nlagr_4)
          call memory_deallo(mem_modul(1:2,modul),'MY_VELOC_NLAGR_4','pts_solite',my_veloc_nlagr_4)
          call memory_deallo(mem_modul(1:2,modul),'MY_COORD_NLAGR_4','pts_solite',my_coord_nlagr_4)
          call memory_deallo(mem_modul(1:2,modul),'PAR_NLAGR_4'     ,'pts_solite',par_nlagr_4)

       end if
    end if
    !
    ! Formats
    !
111 format((e12.6,2x,4(2x,i7),3(2x,e12.6),3(2x,i7)))
112 format((e12.6,2x,4(2x,i7),3(2x,e12.6),3(2x,i7),66(2x,e12.6)))

  end subroutine pts_transport_finalize

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-09-25
  !> @brief   Element search
  !> @details Look for host element for particle at position COORD_PK1
  !> 
  !-----------------------------------------------------------------------

  subroutine pts_transport_element_search(ielem,coord_kp1,ifoun,jelem,shapf)

    integer(ip), intent(in)  :: ielem          !< Initial guess
    real(rp),    intent(in)  :: coord_kp1(3)   !< Particle coordinates
    integer(ip), intent(out) :: ifoun          !< If element is found
    integer(ip), intent(out) :: jelem          !< Last element checked
    real(rp),    intent(out) :: shapf(mnode)   !< Element shape function
    integer(ip)              :: pelty,pnode
    integer(ip)              :: ipoin,inode
    integer(ip)              :: ii,jj,kk,ielel,melel
    real(rp)                 :: elcod(ndime,mnode)
    real(rp)                 :: coloc(3)
    real(rp)                 :: deriv(ndime,mnode)

    ifoun = 0
    jelem = ielem
    pelty = ltype(jelem)
    pnode = nnode(pelty)
    do inode = 1,pnode
       ipoin = lnods(inode,jelem)
       elcod(1:ndime,inode) = coord(1:ndime,ipoin)
    end do
    call elmgeo_natural_coordinates(          &
         ndime,pelty,pnode,elcod,shapf,deriv, &
         coord_kp1,coloc,ifoun,toler)

    if( ifoun == 0 ) then
       if( kfl_usbin_pts == 1 ) then
          !
          ! Use element neighboring bin to reduce the search
          !
          call maths_mapping_coord_to_3d(&
               ndime,element_bin_boxes,element_bin(ielem) % comin,&
               element_bin(ielem) % comax,coord_kp1,ii,jj,kk)
          if( ii*jj*kk /= 0 ) then
             ielel = 1
             melel = element_bin(ielem) % bin_size(ii,jj,kk)
             do while( ifoun == 0 .and. ielel <= melel )
                jelem = element_bin(ielem) % list_elements(ii,jj,kk) % l(ielel)
                if( jelem /= ielem ) then
                   pelty = ltype(jelem)
                   pnode = nnode(pelty)
                   do inode = 1,pnode
                      ipoin = lnods(inode,jelem)
                      elcod(1:ndime,inode) = coord(1:ndime,ipoin)
                   end do
                   call elmgeo_natural_coordinates(         &
                        ndime,pelty,pnode,elcod,shapf,deriv,&
                        coord_kp1,coloc,ifoun,toler)
                end if
                ielel = ielel + 1
             end do
          end if
       else
          !
          ! Look for host element JELEM in neigboring list of IELEM
          !
          ielel = pelel_2(ielem)
          melel = pelel_2(ielem+1)-1
          do while( ifoun == 0 .and. ielel <= melel )
             jelem = lelel_2(ielel)
             pelty = ltype(jelem)
             pnode = nnode(pelty)
             do inode = 1,pnode
                ipoin = lnods(inode,jelem)
                elcod(1:ndime,inode) = coord(1:ndime,ipoin)
             end do
             call elmgeo_natural_coordinates(          &
                  ndime,pelty,pnode,elcod,shapf,deriv, &
                  coord_kp1,coloc,ifoun,toler)
             ielel = ielel + 1
          end do
       end if
    end if

  end subroutine pts_transport_element_search

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-09-25
  !> @brief   Length scale
  !> @details Compute a length scale to normalize the error used
  !>          in the adaptive time step strategy
  !> 
  !-----------------------------------------------------------------------

  subroutine pts_transport_length_scale(itype,diame,denpa,visfl,denfl,veloc_kp1,Du,hleng,h)

    integer(ip), intent(in)  :: itype
    real(rp),    intent(in)  :: diame
    real(rp),    intent(in)  :: denpa
    real(rp),    intent(in)  :: visfl
    real(rp),    intent(in)  :: denfl
    real(rp),    intent(in)  :: veloc_kp1(3)
    real(rp),    intent(in)  :: Du
    real(rp),    intent(in)  :: hleng
    real(rp),    intent(out) :: h
    real(rp)                 :: uu

    if( parttyp(itype) % chale == -3.0_rp ) then
       !
       ! Tau, h = tau * |uf|
       !
       uu = dot_product(veloc_kp1(1:ndime),veloc_kp1(1:ndime))
       h  = denpa * diame * diame / ( 18.0_rp * visfl ) * (sqrt(uu)+zeror)

    else if( parttyp(itype) % chale == -2.0_rp ) then
       !
       ! Re = 1
       !
       h = visfl / ( denfl * ( Du + zeror ) )

    else if( parttyp(itype) % chale == -1.0_rp ) then
       !
       ! Mesh size
       !
       h = hleng

    else if( parttyp(itype) % chale == 0.0_rp ) then
       !
       ! Particle diameter
       !
       h = sqrt(diame)

    else
       !
       ! User-defined
       !
       h = parttyp(itype) % chale

    end if
  end subroutine pts_transport_length_scale

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-09-27
  !> @brief   Compute parameters
  !> @details Compute some derivated parameters of the particle
  !>          PARTICLE % DISTA   ... distance covered by particle
  !>          PARTICLE % COORD1D ... 1d coordinate of particle following
  !>                                 its path
  !> 
  !-----------------------------------------------------------------------

  subroutine pts_transport_derivated_parameters(particle)

    type(latyp), intent(inout) :: particle    !< Particle type
    real(rp)                   :: dista
    real(rp)                   :: xvect_old(3)
    real(rp)                   :: xvect_new(3)

    xvect_new = particle % coord   - particle % coord_k
    xvect_old = particle % coord_k - particle % coord_km1
    dista     = sqrt(dot_product(xvect_new,xvect_new))
    if( dot_product(xvect_new,xvect_old) < 0.0_rp .and. particle % ittim > 1 ) then
       particle % sign = -particle % sign 
    end if

    particle % dista   = particle % dista   + dista
    particle % coord1d = particle % coord1d + particle % sign * dista

  end subroutine pts_transport_derivated_parameters

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-09-27
  !> @brief   Treat special walls
  !> @details Treatment of slip and bouncing walls
  !> 
  !-----------------------------------------------------------------------

  subroutine pts_transport_special_walls(&
       particle,ielem,pnode,shapf,deriv,elcod,diame,coord_kp1,veloc_kp1,accel_kp1)

    use mod_elmgeo
    type(latyp), intent(inout) :: particle
    integer(ip), intent(in)    :: ielem              !< Element number
    integer(ip), intent(in)    :: pnode              !< Number of nodes
    real(rp),    intent(in)    :: shapf(pnode)       !< Shape function in element
    real(rp),    intent(in)    :: deriv(ndime,pnode) !< Shape function in element
    real(rp),    intent(in)    :: elcod(ndime,pnode) !< Shape function in element
    real(rp),    intent(in)    :: diame              !< Particle diameter
    real(rp),    intent(inout) :: coord_kp1(3)       !< Particle new coordinate
    real(rp),    intent(inout) :: veloc_kp1(3)       !< Particle new velocity
    real(rp),    intent(inout) :: accel_kp1(3)       !< Particle new acceleration
    integer(ip)                :: idime
    real(rp)                   :: ubas(3)
    real(rp)                   :: abas(3)
    real(rp)                   :: xbas(3)
    real(rp)                   :: xold(3)
    real(rp)                   :: umsh(3)
    real(rp)                   :: dista,bnorm,radius!,fmix
    real(rp)                   :: basis(ndime,ndime),rmax
    real(rp)                   :: xvect(3),xdotn
    integer(ip)                :: special_condition
    real(rp)                   :: gpcar(ndime,pnode)
    !
    ! Check if we touch a slip or bouncing wall
    !
    special_condition = -1
    radius            = 0.5_rp*diame
    rmax              = 1.2_rp*radius

    if( kfl_slip_wall_pts > 0 ) then
       dista = dot_product(shapf(1:pnode),walld_slip_pts(lnods(1:pnode,ielem)))
       if( dista <= rmax ) special_condition = PTS_SLIP_CONDITION
    end if
    if( kfl_bouncing_wall_pts > 0 ) then
       dista = dot_product(shapf(1:pnode),walld_bouncing_pts(lnods(1:pnode,ielem)))
       if( dista <= rmax ) special_condition = PTS_BOUNCING_CONDITION
    end if
    if( special_condition > 0 ) then
       !
       ! BASIS(1:NDIME,1)   = normal to the wall
       ! BASIS(1:NDIME,2-3) = other vectors of basis
       !
       !OJOif( associated(walld_slip_pts) ) then
       if( special_condition == PTS_SLIP_CONDITION) then
          call elmgeo_cartesian_derivatives(ndime,pnode,elcod,deriv,gpcar)
          do idime = 1,ndime
             basis(idime,1) = -dot_product(gpcar(idime,1:pnode),walld_slip_pts(lnods(1:pnode,ielem)))
          end do
       else
          call elmgeo_cartesian_derivatives(ndime,pnode,elcod,deriv,gpcar)
          do idime = 1,ndime
             basis(idime,1) = -dot_product(gpcar(idime,1:pnode),walld_bouncing_pts(lnods(1:pnode,ielem)))
          end do
          !do idime = 1,ndime
          ! basis(idime,1) =  dot_product(shapf(1:pnode),walln(idime,lnods(1:pnode,ielem)))
          !end do
       end if

       bnorm            = dot_product(basis(1:ndime,1),basis(1:ndime,1))
       basis(1:ndime,1) = basis(1:ndime,1) / sqrt(bnorm)
       call maths_local_orthonormal_basis(ndime,basis)
       !
       ! Check if we are going out the wall
       !
       xvect = coord_kp1 - particle % coord
       xdotn = dot_product(xvect(1:ndime),basis(1:ndime,1))
       if( xdotn > 0.0_rp ) then

          
          xbas = coord_kp1
          xold = particle % coord
          ubas = veloc_kp1
          abas = accel_kp1
          call maths_vector_to_new_basis(ndime,basis,xbas)
          call maths_vector_to_new_basis(ndime,basis,xold)
          call maths_vector_to_new_basis(ndime,basis,ubas)
          call maths_vector_to_new_basis(ndime,basis,abas)
          
          if( if_moving_mesh_pts ) then
             do idime = 1,ndime
                umsh(idime) = dot_product(shapf(1:pnode),velom(idime,lnods(1:pnode,ielem)))
             end do
             call maths_vector_to_new_basis(ndime,basis,umsh)
          else
             umsh = 0.0_rp
          end if
          if(      special_condition == PTS_BOUNCING_CONDITION ) then
             !
             ! Invert normal velocity
             !
             xbas(1) =  xold(1) 
             ubas(1) = -ubas(1)-umsh(1)
             abas(1) =  abas(1)
             
          else if( special_condition == PTS_SLIP_CONDITION     ) then
             !
             ! Cancel out normal velocity
             !
             !fmix    =  max(0.0_rp,min((dista-radius)/(rmax-radius),1.0_rp))
             xbas(1) =  xold(1)
             ubas(1) =  umsh(1)
             abas(1) =  0.0_rp !-abas(1)
          end if

          call maths_vector_from_new_basis(ndime,basis,xbas)
          call maths_vector_from_new_basis(ndime,basis,ubas)
          call maths_vector_from_new_basis(ndime,basis,abas)
          coord_kp1 = xbas
          veloc_kp1 = ubas
          accel_kp1 = abas

       end if
    end if

  end subroutine pts_transport_special_walls

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-10-06
  !> @brief   Adams Bashforth
  !> @details Adams Bashforth method to advance in time.
  !>          https://en.wikipedia.org/wiki/Linear_multistep_method#Two-step_Adams%E2%80%93Bashforth
  !>
  !>          This method is known as a two-step method. More precisely, is is known as the second-order
  !>          Adams-Bashforth method (or AB method) dating back to 1883
  !>          If dt is constant, it gives:
  !>          x^n+1 = x^n + 3/2*dt*u^n -1/2*dt*u^n-1
  !>         
  !>              k-1     k      k+1
  !>           o---|------|-------|--------o
  !>                <----> <------>
  !>                 dto      dt
  !>         
  !>         >          From Taylor series:
  !>          (1)  x^{k+1} = x^k + dt * uk + 1/2 a^k * dt^2 + O(dt^3)
  !>          (2)  u^{k-1} = u^k - a^k dto + O(dto^2)
  !>          where u^k = u(t^k,x^k)
  !>         
  !>          From (2) we have:
  !>          a^k = ( u^k - u ^{k-1} ) / dto + O(dto)
  !>          Substitute this in (1):
  !>          x^{k+1} = x^k + dt * uk + 1/2 ( u^k - u^{k-1} ) / dto * dt^2 + O(dto*dt^2) + O(dt^3)
  !>          x^{k+1} = x^k + ( dt + 1/2 dt^2 / dto ) u^k - 1/2 dt^2 / dto u^{k-1}
  !> 
  !-----------------------------------------------------------------------

  subroutine pts_transport_Adams_Bashforth(&
       particle,dt_k,coord_kp1,veloc_kp1,&
       accel_kp1)

    type(latyp), intent(in)  :: particle          !< Particle
    real(rp),    intent(in)  :: dt_k              !< Current time step
    real(rp),    intent(out) :: coord_kp1(3)      !< New coordinates
    real(rp),    intent(out) :: veloc_kp1(3)      !< New velocity
    real(rp),    intent(out) :: accel_kp1(3)      !< New acceleration
    real(rp)                 :: alpha_k,alpha_km1
    real(rp)                 :: alpha_km2,beta2
    real(rp)                 :: dt012,dt01,dt12
    real(rp)                 :: dt_km1,dt_km2

    dt_km1 = particle % dt_km1
    dt_km2 = particle % dt_km2

    if( abs(kfl_order_pts) == 2 ) then
       !
       ! Second order AB
       !
       beta2              =  dt_k / dt_km1
       alpha_k            =  dt_k * ( 1.0_rp + 0.5_rp * beta2 )
       alpha_km1          = -dt_k * 0.5_rp * beta2
       veloc_kp1(1:ndime) =  particle % v_fluid_k(1:ndime)
       accel_kp1(1:ndime) =  ( particle % v_fluid_k(1:ndime) - particle % v_fluid_km1(1:ndime) ) / dt_km1
       coord_kp1(1:ndime) =  particle % coord(1:ndime)   &
            &                + alpha_k   * particle % v_fluid_k(1:ndime)  &
            &                + alpha_km1 * particle % v_fluid_km1(1:ndime)

    else if( abs(kfl_order_pts) == 3 ) then
       !
       ! Third order AB
       !
       dt012              = dt_km2 + dt_km1 + dt_k
       dt01               = dt_km1 + dt_k
       dt12               = dt_km2 + dt_km1
       alpha_k            = ( 0.25_rp*(dt012*dt01*(dt012+dt01)-dt12*dt_km1*(dt12+dt_km1))&
            &               + 1.0_rp/12.0_rp*(-dt012**3-dt01**3+dt12**3+dt_km1**3))/( dt_km1*dt12)
       alpha_km1          = ( 0.25_rp*(dt012*dt_k *(dt012+dt_k))&
            &               + 1.0_rp/12.0_rp*(-dt012**3-dt_k**3 +dt12**3))/(-dt_km2*dt_km1)
       alpha_km2          = ( 0.25_rp*(dt01 *dt_k *(dt01 +dt_k))&
            &               + 1.0_rp/12.0_rp*(-dt01**3 -dt_k**3 +dt_km1**3))/( dt_km2*dt12)
       veloc_kp1(1:ndime) = particle % v_fluid_k(1:ndime)
       accel_kp1(1:ndime) = ( particle % v_fluid_k(1:ndime) - particle % v_fluid_km1(1:ndime) ) / dt_km1
       coord_kp1(1:ndime) = particle % coord(1:ndime)                     &
            &               + alpha_k   * particle % v_fluid_k(1:ndime)   &
            &               + alpha_km1 * particle % v_fluid_km1(1:ndime) &
            &               + alpha_km2 * particle % v_fluid_km2(1:ndime)

    else
       
       call runend('PTS_TRANSPORT_ADAMS_BASHOFRTH: WRONG ORDER')
       
    end if
    
  end subroutine pts_transport_Adams_Bashforth

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-10-06
  !> @brief   Newton F=ma
  !> @details Newmark/Newton-Raphson method to solve Newton's law
  !> 
  !-----------------------------------------------------------------------

  subroutine pts_transport_newmark_newton_raphson(&
       particle,t,dt_k,Cc,visfl,denfl,vorti_fl,&
       eps,dista,diame,denpa,spher,grafo,buofo,&
       coord_kp1,veloc_kp1,accel_kp1,Cd,Re,tau_p,&
       Stk,istat)

    type(latyp), intent(inout) :: particle          !< Particle
    real(rp),    intent(in)    :: t                 !< Current time 
    real(rp),    intent(in)    :: dt_k              !< Current time step
    real(rp),    intent(in)    :: Cc                !< Cunningham slip correction factor
    real(rp),    intent(in)    :: visfl             !< Fluid viscosity
    real(rp),    intent(in)    :: denfl             !< Fluid density
    real(rp),    intent(in)    :: vorti_fl(3)       !< Fluid vorticity
    real(rp),    intent(in)    :: eps(3)            !< Brownian motion
    real(rp),    intent(in)    :: dista             !< Distance to wall
    real(rp),    intent(in)    :: diame             !< Diameter
    real(rp),    intent(in)    :: denpa             !< Density
    real(rp),    intent(in)    :: spher             !< Sphericity
    real(rp),    intent(in)    :: grafo             !< Gravity force
    real(rp),    intent(in)    :: buofo             !< Buoyancy
    real(rp),    intent(out)   :: coord_kp1(3)      !< New coordinates
    real(rp),    intent(inout) :: veloc_kp1(3)      !< New velocity
    real(rp),    intent(inout) :: accel_kp1(3)      !< New acceleration
    real(rp),    intent(out)   :: Cd                !< Drag coefficient
    real(rp),    intent(out)   :: Re                !< Reynolds number
    real(rp),    intent(out)   :: tau_p             !< Tau particle
    real(rp),    intent(out)   :: Stk(2)            !< Stokes number
    integer(ip), intent(out)   :: istat             !< Status (converged/not converged)
    integer(ip)                :: niter,iiter,itype
    real(rp)                   :: xerro,urel(3)
    real(rp)                   :: Du,xnume
    real(rp)                   :: ff(3),df(3),dRedu(3)
    real(rp)                   :: dCddRe,deltu(3)
    real(rp)                   :: alpha,accel(3)
    real(rp)                   :: tauinv,tau,uf,wf
    real(rp)                   :: Cls,Res
    real(rp)                   :: nu,CdRe,xdeno
    !
    ! Initialization
    !
    niter = 100                  ! Max # of iterations
    xerro = 1.0_rp               ! Error
    iiter = 0                    ! Iteration number
    itype = particle % itype     ! Particle type
    nu    = visfl / denfl        ! Fluid kinematic viscosity
    accel = 0.0_rp               ! Constant acceleration
    !
    ! External force
    !
    if( parttyp(itype) % kfl_extfo /= 0 ) then
       call extefo(&
            parttyp(itype) % kfl_extfo,particle % coord,&
            denpa,spher,denfl,visfl,t,particle % accee)
       accel = accel + particle % accee
    end if
    !
    ! Gravity and buoyancy
    !
    particle % acceg = - g * buofo * denfl / denpa
    accel            = accel + particle % acceg + g * grafo                   
    !
    ! Brownian force
    !
    if( parttyp(itype) % kfl_brown == 2 ) then
       accel = accel + eps
    end if
    !
    ! Avoid doint an extra iteration
    !
    if( parttyp(itype) % kfl_drafo + parttyp(itype) % kfl_extfo == 0 ) niter = 1

    do while( iiter < niter .and. xerro > 1.0e-12_rp )

       iiter     =  iiter + 1
       urel      =  particle % v_fluid_k-veloc_kp1
       Du        =  sqrt(dot_product(urel,urel))+zeror
       df        = -1.0_rp
       accel_kp1 =  accel
       !
       ! Drag force:
       ! a^{n+1} = 1/tau * ( u_fluid - u )
       ! tau is referred to as relaxation time and tau = tau(u)
       ! tau =  ( rho_p * d^2 ) / ( 3/4 * mu * Cd * Re )
       !
       ! aD = rho * Ap / (2 cc) * Cd/mp * |u_fluid-u| (u_fluid-u)
       !
       if( parttyp(itype) % kfl_drafo /= 0 ) then
          if ( diame > 0.0_rp ) then
             call physics_drag_force(parttyp(itype) % kfl_drafo,Du,visfl,denfl,diame,Cd,Re,CdRe,dCddRe,spher,dista)                      
             alpha            =  0.75_rp * visfl / ( denpa * diame * diame * Cc )
             tauinv           =  alpha * CdRe
             tau              =  1.0_rp / ( tauinv + zeror )
             uf               =  sqrt(dot_product(particle % v_fluid_k(1:ndime),particle % v_fluid_k(1:ndime)))
             tau_p            =  (denpa * diame * diame)/(18.0_rp*visfl) 
             Stk              =  tau * uf / diame
             dRedu            = -diame / nu * urel / Du
             df               =  df - dt_k * gamma_pts * alpha * ( CdRe - urel * dRedu * dCddRe )
             particle % acced =  alpha * CdRe * urel
             accel_kp1        =  accel_kp1 + particle % acced
          else
             alpha            = 1.0_rp 
             tauinv           = 1.0_rp 
             tau              = 0.0_rp 
             uf               = 1.0_rp 
             tau_p            = 0.0_rp 
             Stk              = 1.0_rp 
             dRedu            = 0.0_rp 
             df               = 1.0_rp 
             particle % acced = 0.0_rp 
             accel_kp1        = 0.0_rp 
          endif
       end if
       !
       ! Saffman Mei
       !
       if( parttyp(itype) % kfl_saffm /= 0 ) then
          wf = sqrt(dot_product(vorti_fl(1:ndime),vorti_fl(1:ndime)))
          call physics_lift_force(parttyp(itype) % kfl_saffm,Du,wf,visfl,denfl,diame,dista,Cls,Re,Res) 
          accel_kp1(1) = accel_kp1(1) + denfl * pi/8.0_rp * diame * diame * diame * Cls * ( urel(2) * vorti_fl(3) - urel(3) * vorti_fl(2) )
          accel_kp1(2) = accel_kp1(2) + denfl * pi/8.0_rp * diame * diame * diame * Cls * ( urel(3) * vorti_fl(1) - urel(1) * vorti_fl(3) )
          accel_kp1(3) = accel_kp1(3) + denfl * pi/8.0_rp * diame * diame * diame * Cls * ( urel(1) * vorti_fl(2) - urel(2) * vorti_fl(1) )
       end if
       ff        =  ( -veloc_kp1 + particle % veloc + dt_k * ( (1.0_rp-gamma_pts)*particle % accel + gamma_pts*(accel_kp1) ) )
       deltu     = -ff / ( df + zeror )
       veloc_kp1 =  veloc_kp1 + deltu
       !
       ! Residual
       !
       xnume = dot_product(deltu(1:ndime),deltu(1:ndime))
       xdeno = dot_product(veloc_kp1(1:ndime),veloc_kp1(1:ndime))
       xnume = sqrt(xnume)
       xdeno = sqrt(xdeno) + zeror
       xerro = xnume / xdeno
    end do
    !
    ! Update new position or go back in time if not converged
    !
    if( iiter == niter .and. niter /= 1 ) then
       istat = 1
    else
       istat = 0
       coord_kp1 = particle % coord                          &
            &      + dt_k * particle % veloc                 &
            &      + dt_k * dt_k * ( 0.5_rp*particle % accel &
            &      + beta_pts *(accel_kp1-particle % accel) )
    end if

  end subroutine pts_transport_newmark_newton_raphson
  
end module mod_pts_transport
!> @}
