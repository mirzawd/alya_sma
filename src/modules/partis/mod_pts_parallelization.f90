!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Partis
!> @{
!> @file    mod_pts_parallelization.f90
!> @author  houzeaux
!> @date    2018-09-21
!> @brief   Parallelization subroutines
!> @details Some parallel functions for distributed memory
!>
!>          To add a postprocess variable:
!>          1. Check mvarp_pts is sufficently high in def_partis
!>          2. Add the variable in pack and unpack (option LIST_OF_PROPERTIES)
!>
!>          To add a variable to be migrated:
!>          1. Modify number_migrated_variables
!>          2. Add the variable in pack and unpack
!>
!>          Do not forget to add the same variables in pts_restar!
!>
!-----------------------------------------------------------------------

module mod_pts_parallelization

  use def_kintyp_basic,   only : ip,rp,lg
  use def_kintyp_comm,    only : comm_data_par
  use def_master,         only : kfl_paral
  use def_master,         only : modul
  use def_master,         only : mem_modul
  use mod_memory,         only : memory_alloca
  use mod_memory,         only : memory_deallo
  use mod_memory,         only : memory_resize
  use mod_memory,         only : memory_size
  use mod_communications, only : PAR_SEND_RECEIVE
  use def_partis,         only : nlagr_local_pts
  use def_partis,         only : nlagr_free_pts
  use def_partis,         only : permu_nlagr_pts
  use def_partis,         only : pts_name_to_variable_number
  use def_partis,         only : migrated_variables_pts
  use def_partis,         only : number_migrated_variables_pts
  use def_partis,         only : mvarp_pts
  use def_partis,         only : lagrtyp
  use def_partis,         only : latyp
  use def_partis,         only : mlagr
  use def_partis,         only : PTS_PARTICLE_EXISTS
  use def_partis,         only : kfl_thermo_pts
  use def_partis,         only : postprocess_var_pts
  use def_partis,         only : ntyla_pts 
  use def_partis,         only : parttyp

  implicit none

  private 

  interface pts_parallelization_unpack
     module procedure pts_parallelization_unpack_all,&
          &           pts_parallelization_unpack_array,&
          &           pts_parallelization_unpack_single
  end interface  pts_parallelization_unpack

  public :: pts_parallelization_migration
  public :: pts_parallelization_pack
  public :: pts_parallelization_unpack
  public ::pts_parallelization_initialization
  
contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-09-21
  !> @brief   Pack a type
  !> @details Pack some given variables particle type into an array
  !> 
  !-----------------------------------------------------------------------

  subroutine pts_parallelization_unpack_single(particle,ii,xx,LIST_OF_PROPERTIES,wname,status)

    type(latyp),           intent(inout) :: particle
    integer(ip),           intent(inout) :: ii
    real(rp),    pointer,  intent(in)    :: xx(:)
    logical(lg), pointer,  intent(in)    :: LIST_OF_PROPERTIES(:)
    character(*),          intent(in)    :: wname
    integer(ip), optional, intent(out)   :: status
    character(5)                         :: wname_array(1)

    wname_array = wname
    call pts_parallelization_unpack_array(particle,ii,xx,LIST_OF_PROPERTIES,wname_array,status)

  end subroutine pts_parallelization_unpack_single

  subroutine pts_parallelization_unpack_array(particle,ii,xx,LIST_OF_PROPERTIES,wname,status) 

    type(latyp),           intent(inout) :: particle
    integer(ip),           intent(inout) :: ii
    real(rp),    pointer,  intent(in)    :: xx(:)
    logical(lg), pointer,  intent(in)    :: LIST_OF_PROPERTIES(:)
    character(*),          intent(in)    :: wname(:)
    integer(ip), optional, intent(out)   :: status
    integer(ip)                          :: jj,kk,ll
    do kk = 1,size(wname)
       jj = pts_name_to_variable_number(trim(wname(kk)))
       if( jj == 0 ) then
          particle % ilagr = 0
          if( present(status) ) status = -1
       else
          ll = count(LIST_OF_PROPERTIES(1:jj)) + ii
          if(      trim(wname(kk)) == 'ILAGR' ) then
             particle % ilagr = int(xx(ll),ip)
          else if( trim(wname(kk)) == 'EXIST' ) then
             particle % kfl_exist = int(xx(ll),ip)
          else
             call runend('UNPACK NOT CODED YET')
          end if

       end if
    end do

  end subroutine pts_parallelization_unpack_array

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-09-21
  !> @brief   Pack a type
  !> @details Pack a particle type into an array
  !> 
  !-----------------------------------------------------------------------

  subroutine pts_parallelization_pack(particle,ii,xx,LIST_OF_PROPERTIES)

    type(latyp),          intent(in)    :: particle
    integer(ip),          intent(inout) :: ii
    real(rp),    pointer, intent(inout) :: xx(:)
    logical(lg),          intent(in)    :: LIST_OF_PROPERTIES(:)
    integer(ip)                         :: kk
    
    do kk = 1,size(LIST_OF_PROPERTIES)
       if( LIST_OF_PROPERTIES(kk) ) then
          select case ( kk )  
          case (  1 ) ; ii = ii + 1 ; xx(ii) = particle % t               
          case (  2 ) ; ii = ii + 1 ; xx(ii) = real(particle % ilagr,rp)
          case (  3 ) ; ii = ii + 1 ; xx(ii) = real(particle % itype,rp) 
          case (  4 ) ; ii = ii + 1 ; xx(ii) = real(particle % ielem,rp)
          case (  5 ) ; ii = ii + 1 ; xx(ii) = real(particle % kfl_exist,rp)
          case (  6 ) ; ii = ii + 1 ; xx(ii) = real(particle % ittim,rp)
          case (  7 ) ; ii = ii + 1 ; xx(ii) = real(particle % boundary_set,rp)                             
          case (  8 ) ; ii = ii + 1 ; xx(ii) = particle % coord(1)
          case (  9 ) ; ii = ii + 1 ; xx(ii) = particle % coord(2)
          case ( 10 ) ; ii = ii + 1 ; xx(ii) = particle % coord(3)
          case ( 11 ) ; ii = ii + 1 ; xx(ii) = particle % veloc(1)
          case ( 12 ) ; ii = ii + 1 ; xx(ii) = particle % veloc(2)
          case ( 13 ) ; ii = ii + 1 ; xx(ii) = particle % veloc(3)
          case ( 14 ) ; ii = ii + 1 ; xx(ii) = particle % accel(1)
          case ( 15 ) ; ii = ii + 1 ; xx(ii) = particle % accel(2)
          case ( 16 ) ; ii = ii + 1 ; xx(ii) = particle % accel(3)
          case ( 17 ) ; ii = ii + 1 ; xx(ii) = particle % dt_k
          case ( 18 ) ; ii = ii + 1 ; xx(ii) = particle % Cd
          case ( 19 ) ; ii = ii + 1 ; xx(ii) = particle % Stk(1)
          case ( 20 ) ; ii = ii + 1 ; xx(ii) = particle % Stk(2)
          case ( 21 ) ; ii = ii + 1 ; xx(ii) = particle % v_fluid_k(1)
          case ( 22 ) ; ii = ii + 1 ; xx(ii) = particle % v_fluid_k(2)
          case ( 23 ) ; ii = ii + 1 ; xx(ii) = particle % v_fluid_k(3)                
          case ( 24 ) ; ii = ii + 1 ; xx(ii) = particle % acced(1)
          case ( 25 ) ; ii = ii + 1 ; xx(ii) = particle % acced(2)
          case ( 26 ) ; ii = ii + 1 ; xx(ii) = particle % acced(3)
          case ( 27 ) ; ii = ii + 1 ; xx(ii) = particle % accee(1)
          case ( 28 ) ; ii = ii + 1 ; xx(ii) = particle % accee(2)
          case ( 29 ) ; ii = ii + 1 ; xx(ii) = particle % accee(3)
          case ( 30 ) ; ii = ii + 1 ; xx(ii) = particle % acceg(1)
          case ( 31 ) ; ii = ii + 1 ; xx(ii) = particle % acceg(2)
          case ( 32 ) ; ii = ii + 1 ; xx(ii) = particle % acceg(3)
          case ( 33 ) ; ii = ii + 1 ; xx(ii) = particle % stret
          case ( 34 ) ; ii = ii + 1 ; xx(ii) = particle % dt_km1
          case ( 35 ) ; ii = ii + 1 ; xx(ii) = particle % dt_km2
          case ( 36 ) ; ii = ii + 1 ; xx(ii) = particle % dtg
          case ( 37 ) ; ii = ii + 1 ; xx(ii) = particle % v_fluid_km1(1)
          case ( 38 ) ; ii = ii + 1 ; xx(ii) = particle % v_fluid_km1(2)
          case ( 39 ) ; ii = ii + 1 ; xx(ii) = particle % v_fluid_km1(3)
          case ( 40 ) ; ii = ii + 1 ; xx(ii) = particle % v_fluid_km2(1)
          case ( 41 ) ; ii = ii + 1 ; xx(ii) = particle % v_fluid_km2(2)
          case ( 42 ) ; ii = ii + 1 ; xx(ii) = particle % v_fluid_km2(3)
          case ( 43 ) ; ii = ii + 1 ; xx(ii) = particle % t_inject
          case ( 44 ) ; ii = ii + 1 ; xx(ii) = particle % coord_k(1)
          case ( 45 ) ; ii = ii + 1 ; xx(ii) = particle % coord_k(2)
          case ( 46 ) ; ii = ii + 1 ; xx(ii) = particle % coord_k(3)
          case ( 47 ) ; ii = ii + 1 ; xx(ii) = particle % coord_km1(1)
          case ( 48 ) ; ii = ii + 1 ; xx(ii) = particle % coord_km1(2)
          case ( 49 ) ; ii = ii + 1 ; xx(ii) = particle % coord_km1(3)
          case ( 50 ) ; ii = ii + 1 ; xx(ii) = particle % dista
          case ( 51 ) ; ii = ii + 1 ; xx(ii) = particle % coord1d
          case ( 52 ) ; ii = ii + 1 ; xx(ii) = particle % sign
          case ( 53 ) ; ii = ii + 1 ; xx(ii) = particle % tempe_k
          case ( 54 ) ; ii = ii + 1 ; xx(ii) = particle % tempe_km1
          case ( 55 ) ; ii = ii + 1 ; xx(ii) = particle % mass_k
          case ( 56 ) ; ii = ii + 1 ; xx(ii) = particle % mass_km1
          case ( 57 ) ; ii = ii + 1 ; xx(ii) = real(particle % mpi_rank,rp)
          case ( 58 ) ; ii = ii + 1 ; xx(ii) = particle % diam_k
          case ( 59 ) ; ii = ii + 1 ; xx(ii) = particle % Temp_fluid_k
          case ( 60 ) ; ii = ii + 1 ; xx(ii) = particle % Yvap_fluid_k
          case ( 61 ) ; ii = ii + 1 ; xx(ii) = particle % diam_0
          case ( 62 ) ; ii = ii + 1 ; xx(ii) = particle % mass_0
          case ( 63 ) ; ii = ii + 1 ; xx(ii) = particle % BM
          end select
       end if
    end do

  end subroutine pts_parallelization_pack

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-09-21
  !> @brief   Unpack an array
  !> @details Unpack an array into a a particle type
  !> 
  !-----------------------------------------------------------------------

  subroutine pts_parallelization_unpack_all(particle,ii,xx,LIST_OF_PROPERTIES) 

    type(latyp),          intent(inout) :: particle
    integer(ip),          intent(inout) :: ii
    real(rp),    pointer, intent(in)    :: xx(:)
    logical(lg),          intent(in)    :: LIST_OF_PROPERTIES(:)
    integer(ip)                         :: kk

       do kk = 1,size(LIST_OF_PROPERTIES)
          if( LIST_OF_PROPERTIES(kk) ) then
             select case ( kk )  
             case (  1 ) ; ii = ii + 1 ; particle % t               = xx(ii) 
             case (  2 ) ; ii = ii + 1 ; particle % ilagr           = int(xx(ii),ip) 
             case (  3 ) ; ii = ii + 1 ; particle % itype           = int(xx(ii),ip)
             case (  4 ) ; ii = ii + 1 ; particle % ielem           = int(xx(ii),ip)
             case (  5 ) ; ii = ii + 1 ; particle % kfl_exist       = int(xx(ii),ip)
             case (  6 ) ; ii = ii + 1 ; particle % ittim           = int(xx(ii),ip)
             case (  7 ) ; ii = ii + 1 ; particle % boundary_set    = int(xx(ii),ip)         
             case (  8 ) ; ii = ii + 1 ; particle % coord(1)        = xx(ii) 
             case (  9 ) ; ii = ii + 1 ; particle % coord(2)        = xx(ii) 
             case ( 10 ) ; ii = ii + 1 ; particle % coord(3)        = xx(ii) 
             case ( 11 ) ; ii = ii + 1 ; particle % veloc(1)        = xx(ii) 
             case ( 12 ) ; ii = ii + 1 ; particle % veloc(2)        = xx(ii) 
             case ( 13 ) ; ii = ii + 1 ; particle % veloc(3)        = xx(ii) 
             case ( 14 ) ; ii = ii + 1 ; particle % accel(1)        = xx(ii) 
             case ( 15 ) ; ii = ii + 1 ; particle % accel(2)        = xx(ii) 
             case ( 16 ) ; ii = ii + 1 ; particle % accel(3)        = xx(ii) 
             case ( 17 ) ; ii = ii + 1 ; particle % dt_k            = xx(ii) 
             case ( 18 ) ; ii = ii + 1 ; particle % Cd              = xx(ii) 
             case ( 19 ) ; ii = ii + 1 ; particle % Stk(1)          = xx(ii) 
             case ( 20 ) ; ii = ii + 1 ; particle % Stk(2)          = xx(ii) 
             case ( 21 ) ; ii = ii + 1 ; particle % v_fluid_k(1)    = xx(ii) 
             case ( 22 ) ; ii = ii + 1 ; particle % v_fluid_k(2)    = xx(ii) 
             case ( 23 ) ; ii = ii + 1 ; particle % v_fluid_k(3)    = xx(ii)       
             case ( 24 ) ; ii = ii + 1 ; particle % acced(1)        = xx(ii) 
             case ( 25 ) ; ii = ii + 1 ; particle % acced(2)        = xx(ii) 
             case ( 26 ) ; ii = ii + 1 ; particle % acced(3)        = xx(ii) 
             case ( 27 ) ; ii = ii + 1 ; particle % accee(1)        = xx(ii) 
             case ( 28 ) ; ii = ii + 1 ; particle % accee(2)        = xx(ii) 
             case ( 29 ) ; ii = ii + 1 ; particle % accee(3)        = xx(ii) 
             case ( 30 ) ; ii = ii + 1 ; particle % acceg(1)        = xx(ii) 
             case ( 31 ) ; ii = ii + 1 ; particle % acceg(2)        = xx(ii) 
             case ( 32 ) ; ii = ii + 1 ; particle % acceg(3)        = xx(ii) 
             case ( 33 ) ; ii = ii + 1 ; particle % stret           = xx(ii) 
             case ( 34 ) ; ii = ii + 1 ; particle % dt_km1          = xx(ii) 
             case ( 35 ) ; ii = ii + 1 ; particle % dt_km2          = xx(ii) 
             case ( 36 ) ; ii = ii + 1 ; particle % dtg             = xx(ii) 
             case ( 37 ) ; ii = ii + 1 ; particle % v_fluid_km1(1)  = xx(ii) 
             case ( 38 ) ; ii = ii + 1 ; particle % v_fluid_km1(2)  = xx(ii) 
             case ( 39 ) ; ii = ii + 1 ; particle % v_fluid_km1(3)  = xx(ii) 
             case ( 40 ) ; ii = ii + 1 ; particle % v_fluid_km2(1)  = xx(ii) 
             case ( 41 ) ; ii = ii + 1 ; particle % v_fluid_km2(2)  = xx(ii) 
             case ( 42 ) ; ii = ii + 1 ; particle % v_fluid_km2(3)  = xx(ii) 
             case ( 43 ) ; ii = ii + 1 ; particle % t_inject        = xx(ii) 
             case ( 44 ) ; ii = ii + 1 ; particle % coord_k(1)      = xx(ii)
             case ( 45 ) ; ii = ii + 1 ; particle % coord_k(2)      = xx(ii)
             case ( 46 ) ; ii = ii + 1 ; particle % coord_k(3)      = xx(ii)
             case ( 47 ) ; ii = ii + 1 ; particle % coord_km1(1)    = xx(ii)
             case ( 48 ) ; ii = ii + 1 ; particle % coord_km1(2)    = xx(ii)
             case ( 49 ) ; ii = ii + 1 ; particle % coord_km1(3)    = xx(ii)
             case ( 50 ) ; ii = ii + 1 ; particle % dista           = xx(ii)
             case ( 51 ) ; ii = ii + 1 ; particle % coord1d         = xx(ii)
             case ( 52 ) ; ii = ii + 1 ; particle % sign            = xx(ii)
             case ( 53 ) ; ii = ii + 1 ; particle % tempe_k         = xx(ii)  
             case ( 54 ) ; ii = ii + 1 ; particle % tempe_km1       = xx(ii)
             case ( 55 ) ; ii = ii + 1 ; particle % mass_k          = xx(ii)
             case ( 56 ) ; ii = ii + 1 ; particle % mass_km1        = xx(ii)
             case ( 57 ) ; ii = ii + 1 ; particle % mpi_rank        = int(xx(ii),ip)
             case ( 58 ) ; ii = ii + 1 ; particle % diam_k          = xx(ii)
             case ( 59 ) ; ii = ii + 1 ; particle % Temp_fluid_k    = xx(ii)
             case ( 60 ) ; ii = ii + 1 ; particle % Yvap_fluid_k    = xx(ii)
             case ( 61 ) ; ii = ii + 1 ; particle % diam_0          = xx(ii)
             case ( 62 ) ; ii = ii + 1 ; particle % mass_0          = xx(ii)
             case ( 63 ) ; ii = ii + 1 ; particle % BM              = xx(ii)
             end select
          end if
       end do

  end subroutine pts_parallelization_unpack_all

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-09-21
  !> @brief   Migrate particles
  !> @details Migrate particles to neighboring subdomains
  !> 
  !-----------------------------------------------------------------------

  subroutine pts_parallelization_migration(commu,nlagr_migrating,particles_sent,particles_recv)

    type(comm_data_par),          intent(in)    :: commu               !< Communicator
    integer(ip),                  intent(inout) :: nlagr_migrating(:)  !< Number of particles to send
    integer(ip),       optional,  intent(inout) :: particles_sent      !< Number of particles sent
    integer(ip),       optional,  intent(inout) :: particles_recv      !< Number of particle received
    integer(ip)                                 :: ineig,nlagr_migrate
    integer(ip)                                 :: new_size,ipars
    integer(ip)                                 :: nparr,ilagr_local
    integer(ip)                                 :: nlagr_last,icror
    integer(ip)                                 :: ifoun,nneig,iparr
    integer(ip)                                 :: npasr,ilagr
    integer(ip),         allocatable            :: nlagr_receiving(:)
    real(rp),            pointer                :: parrs_pts(:)
    real(rp),            pointer                :: parre_pts(:)
    integer(ip)                                 :: my_particles_sent   
    integer(ip)                                 :: my_particles_recv

    if( commu % nneig == 0 ) return
    ! 
    ! Nullify
    !
    nullify(parrs_pts)
    nullify(parre_pts)
    ! 
    ! Variables to migrate
    !
    nlagr_migrate = number_migrated_variables_pts
    !
    ! Exchange number of particles to send and receive
    !
    nneig = commu % nneig
    allocate(nlagr_receiving(nneig))
    do ineig = 1,nneig
       call PAR_SEND_RECEIVE(1_ip,1_ip,nlagr_migrating(ineig:ineig),nlagr_receiving(ineig:ineig),'IN MY ZONE',commu % neights(ineig) )
    end do
    !
    ! Count total [number of particles to receive] - [number of particles to send]
    ! UPS!!! I cannot take into account the send particles because the place is freed
    ! as we go through the neigbors
    !
    my_particles_sent = 0
    my_particles_recv = 0
    do ineig = 1,nneig
       my_particles_sent = my_particles_sent + nlagr_migrating(ineig) 
       my_particles_recv = my_particles_recv + nlagr_receiving(ineig)
    end do
    if( present(particles_sent) ) particles_sent = my_particles_sent
    if( present(particles_recv) ) particles_recv = my_particles_recv

    if( my_particles_recv > nlagr_free_pts ) then
       new_size = int(1.2_rp*real(mlagr+my_particles_recv-nlagr_free_pts,rp),ip)
       call pts_reallocate(new_size)
       if( size(permu_nlagr_pts) < mlagr ) then
          call memory_resize(mem_modul(1:2,modul),'PERMU_NLAGR_PTS','pts_solite',permu_nlagr_pts,mlagr)
       end if
    end if
    !
    ! Send/Recv crossing particles to neighbors: coord, tinit, ilagr, etc.
    ! This is the minimum information for the particle to be correctly tracked by
    ! the neighboring subdomain
    !
    do ineig = 1,nneig

       if( nlagr_migrating(ineig) > 0 ) then
          npasr = nlagr_migrating(ineig) * nlagr_migrate
          call memory_alloca(mem_modul(1:2,modul),'PARRS_PTS','pts_solite',parrs_pts,npasr,'DO_NOT_INITIALIZE')
          ipars = 0
          do ilagr_local = 1,nlagr_local_pts
             ilagr = permu_nlagr_pts(ilagr_local)
             if(ilagr/=0)then
               if( lagrtyp(ilagr) % kfl_exist == ineig ) then
                 call pts_parallelization_pack(lagrtyp(ilagr),ipars,parrs_pts,migrated_variables_pts)                      
                 lagrtyp(ilagr) % kfl_exist   = 0         ! Take off particle from my list
                 permu_nlagr_pts(ilagr_local) = 0         ! Take off particle from permutation (useless)
               end if
             end if
          end do
       end if
       
       if( nlagr_receiving(ineig) > 0 ) then
          nparr = nlagr_receiving(ineig) * nlagr_migrate
          call memory_alloca(mem_modul(1:2,modul),'PARRE_PTS','pts_solite',parre_pts,nparr)
       end if

       call PAR_SEND_RECEIVE(parrs_pts,parre_pts,'IN MY ZONE',commu % neights(ineig) )
       !
       ! Allocate new particles
       !
       ilagr      = 1
       iparr      = 0
       nlagr_last = 0

       do icror = 1,nlagr_receiving(ineig)
          !
          ! Look for a free space in LAGRTYP to save particle
          ! If there is no space, allocate more memory
          !
          ifoun = 0
          ilagr = nlagr_last
          loop_find_position: do while( ilagr < mlagr )
             ilagr = ilagr + 1
             if( lagrtyp(ilagr) % kfl_exist == 0 ) then
                ifoun = 1
                exit loop_find_position
             end if
          end do loop_find_position
          nlagr_last = ilagr
          if( ifoun == 0 ) then
             call runend('PTS_SOLITE: WE ARE IN TROUBLE!')
          end if
          call lagrtyp(ilagr) % init()
          lagrtyp(ilagr) % kfl_exist = -1
          lagrtyp(ilagr) % mpi_rank  = kfl_paral
          call pts_parallelization_unpack(lagrtyp(ilagr),iparr,parre_pts,migrated_variables_pts)
        
       end do
       !
       ! Deallocate
       !
       call memory_deallo(mem_modul(1:2,modul),'PARRE_PTS','pts_solite',parre_pts)
       call memory_deallo(mem_modul(1:2,modul),'PARRS_PTS','pts_solite',parrs_pts)
       nlagr_migrating(ineig) = 0

    end do
    !
    ! Recompute permutation
    !
    !call pts_compute_permutation()
    
    deallocate(nlagr_receiving)

  end subroutine pts_parallelization_migration

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-10-03
  !> @brief   Initialization
  !> @details INitialization useful for parallelization purpose 
  !> 
  !-----------------------------------------------------------------------

  subroutine pts_parallelization_initialization()
    integer(ip) :: ivari,itype 
    !
    ! Variables to migrate
    !
    migrated_variables_pts        = .false.
    migrated_variables_pts(pts_name_to_variable_number(    'T')) = .true.
    migrated_variables_pts(pts_name_to_variable_number('ILAGR')) = .true.
    migrated_variables_pts(pts_name_to_variable_number('ITYPE')) = .true.
    migrated_variables_pts(pts_name_to_variable_number('IELEM')) = .true.
    migrated_variables_pts(pts_name_to_variable_number('ITTIM')) = .true.
    migrated_variables_pts(pts_name_to_variable_number('COORX')) = .true.
    migrated_variables_pts(pts_name_to_variable_number('COORY')) = .true.
    migrated_variables_pts(pts_name_to_variable_number('COORZ')) = .true.
    migrated_variables_pts(pts_name_to_variable_number('VELOX')) = .true.
    migrated_variables_pts(pts_name_to_variable_number('VELOY')) = .true.
    migrated_variables_pts(pts_name_to_variable_number('VELOZ')) = .true.
    migrated_variables_pts(pts_name_to_variable_number('ACCEX')) = .true.
    migrated_variables_pts(pts_name_to_variable_number('ACCEY')) = .true.
    migrated_variables_pts(pts_name_to_variable_number('ACCEZ')) = .true.
    migrated_variables_pts(pts_name_to_variable_number(  'DTK')) = .true.
    migrated_variables_pts(pts_name_to_variable_number(   'CD')) = .true.
    migrated_variables_pts(pts_name_to_variable_number( 'STK1')) = .true. 
    migrated_variables_pts(pts_name_to_variable_number( 'STK2')) = .true. 
    migrated_variables_pts(pts_name_to_variable_number('VFM1X')) = .true.
    migrated_variables_pts(pts_name_to_variable_number('VFM1Y')) = .true.
    migrated_variables_pts(pts_name_to_variable_number('VFM1Z')) = .true.
    migrated_variables_pts(pts_name_to_variable_number('ACCDX')) = .true.
    migrated_variables_pts(pts_name_to_variable_number('ACCDY')) = .true.
    migrated_variables_pts(pts_name_to_variable_number('ACCDZ')) = .true.
    migrated_variables_pts(pts_name_to_variable_number('STRET')) = .true.
    migrated_variables_pts(pts_name_to_variable_number('DTKM1')) = .true.
    migrated_variables_pts(pts_name_to_variable_number('DTKM2')) = .true.
    migrated_variables_pts(pts_name_to_variable_number(  'DTG')) = .true.
    migrated_variables_pts(pts_name_to_variable_number('VFM1X')) = .true.
    migrated_variables_pts(pts_name_to_variable_number('VFM1Y')) = .true.
    migrated_variables_pts(pts_name_to_variable_number('VFM1Z')) = .true.
    migrated_variables_pts(pts_name_to_variable_number('VFM2X')) = .true.
    migrated_variables_pts(pts_name_to_variable_number('VFM2Y')) = .true.
    migrated_variables_pts(pts_name_to_variable_number('VFM2Z')) = .true.
    migrated_variables_pts(pts_name_to_variable_number('TINJE')) = .true.
    migrated_variables_pts(pts_name_to_variable_number( 'COKX')) = .true.
    migrated_variables_pts(pts_name_to_variable_number( 'COKY')) = .true.
    migrated_variables_pts(pts_name_to_variable_number( 'COKZ')) = .true.
    migrated_variables_pts(pts_name_to_variable_number('COK1X')) = .true.
    migrated_variables_pts(pts_name_to_variable_number('COK1Y')) = .true.
    migrated_variables_pts(pts_name_to_variable_number('COK1Z')) = .true.
    migrated_variables_pts(pts_name_to_variable_number('DISTA')) = .true.
    migrated_variables_pts(pts_name_to_variable_number('COO1D')) = .true.
    migrated_variables_pts(pts_name_to_variable_number( 'SIGN')) = .true.
    !
    ! Migrate initial diameter if it's an output
    ! or it's needed for th evaporation
    !
    ivari = pts_name_to_variable_number('DIAM0')
    if (postprocess_var_pts(ivari)) migrated_variables_pts(ivari) = .true.
    do itype = 1,ntyla_pts
       if (abs(parttyp(itype) % kfl_dmini) == 3) migrated_variables_pts(ivari) = .true.
    enddo
    !
    ! Migrate initial mass if it's an output
    ! or it's needed for th evaporation
    !
    ivari = pts_name_to_variable_number('MASS0')
    if (postprocess_var_pts(ivari)) migrated_variables_pts(ivari) = .true.
    do itype = 1,ntyla_pts
       if (abs(parttyp(itype) % kfl_dmini) == 2) migrated_variables_pts(ivari) = .true.
    enddo
   
    !
    ! Migrate mass or diameter if they are outputs
    !
    ivari = pts_name_to_variable_number('DIAMK')
    if (postprocess_var_pts(ivari)) migrated_variables_pts(ivari) = .true.
    ivari = pts_name_to_variable_number('MASSK')
    if (postprocess_var_pts(ivari)) migrated_variables_pts(ivari) = .true.


    if( kfl_thermo_pts /= 0 ) then
       migrated_variables_pts(pts_name_to_variable_number('TEMPK')) = .true. 
       migrated_variables_pts(pts_name_to_variable_number('TEKM1')) = .true.     
       migrated_variables_pts(pts_name_to_variable_number('MASSK')) = .true.
       migrated_variables_pts(pts_name_to_variable_number('MAKM1')) = .true.

       !
       ! Migrate some stuff that would not be recalculated before the output
       !
       ivari = pts_name_to_variable_number('TEMPF')
       if (postprocess_var_pts(ivari)) migrated_variables_pts(ivari) = .true.
       ivari = pts_name_to_variable_number('YVAPF')
       if (postprocess_var_pts(ivari)) migrated_variables_pts(ivari) = .true.
    end if
    
    number_migrated_variables_pts = count(migrated_variables_pts)
    
  end subroutine pts_parallelization_initialization
  
end module mod_pts_parallelization
!> @}

