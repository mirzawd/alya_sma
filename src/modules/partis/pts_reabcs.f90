!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Partis
!> @{
!> @file    pts_reabcs.f90
!> @author  Guillaume Houzeaux
!> @date    16/11/1966   
!> @brief   Read boundary conditions 
!> @details The boundary codes are:
!>          KFL_FIXBOU(IBOUN) = 0 ... Outflow
!>                              1 ... Wall
!>                              3 ... Slip
!>                              4 ... bounce 
!>
!>  To read prescribed velocity fields, these need to be specified in dom.dat and here in pts.dat. Periodic fields are supported. 
!>  Linear interpolation in time is also supported
!>  boundary_conditions
!>     VELOC = field_id
!> 
!>
!>  For bouncing boundaries, WALLD must not be used. Kernel needs only solver for normal propagation
!>
!>  For slip boundairies, WALL_NORMAL nees to be added to ker.dat, potentially for partis
!>  also WALL_DISTANCE needs to be added. Sample definitions can be like this:
!>     NUMERICAL_TREATMENT
!>       WALL_DISTANCE                                
!>         ALGEBRAIC_SOLVER    DEFLATED_CG, ITERA= 5000,TOLER= 1e-9, ADAPTIVE, RATIO=1e-9
!>         PRECONDITIONING     DIAGONAL               
!>         CODES, BOUNDARIES
!>           1 3
!>         END_CODES
!>       END_WALL_DISTANCE 
!>     
!>       WALL_NORMAL
!>         ALGEBRAIC_SOLVER    DEFLATED_CG, ITERA= 5000,TOLER= 1e-9, ADAPTIVE, RATIO=1e-9
!>         PRECONDITIONING     DIAGONAL               
!>         CODES, BOUNDARIES
!>           1 3
!>         END_CODES
!>       END_WALL_NORMAL
!>     END_NUMERICAL_TREATMENT
!>
!>   For everything else:
!>
!>       PHYSICAL_PROBLEM
!>         MAXIMUM int    $maximum number of particles to preallocate memory, if it is smaller than the number of injected particles, the memory will be reallocated
!>         TYPE= int  $1,2,...
!>           MODEL= FORCE
!>           FORCES= DRAG, GRAVITY, BUOYANCY
!>           DENSITY= real
!>           DIAMETER= real
!>         END_TYPE
!>       END_PHYSICAL_PROBLEM
!>
!>
!>       BOUNDARY_CONDITIONS
!>         INJECTOR= int
!>           GEOMETRY: injector_name, PARAM= injector dependent numbers
!>           TYPE:      int | ALL   $ type of the particles to inject (provided in the PHYSICAL_PROBLEM section
!>           DISTRIBUTION:         UNICA | UNIPO   $ spatial particle disrtibution, optional
!>           NPARTICLES: ASIS | [ MASS_FLOW:  CONSTANT= mdot | TIME_FUNCTION=NAMEFUN ] $ mass flow overwrites the number of particles given in the geometry  
!>           STOCA: ON  $random particles, otherwise nonrandom (no clue what is the nonrandom)
!>           INITIAL_TIME: real  $initial time of injection
!>           FINAL_TIME: real  $final time of injection
!>           PERIOD_TIME: real  $time step of injection
!>           VELOCITY: SPRAY, PARAM=velocity, angle_degrees
!>           SIZE_DISTRIBUTION: CONSTANT | UNIFORM PARAM= dmin, dmax | ROSIN_RAMMLER PARAM= dmin, dmax, dmean, nparam $ size distribution overwrites the diameter from the type 
!>         END_INJECTOR
!>         CODES, BOUNDARIES
!>           ...
!>           int1 1                                    $   int1=id of the boundary, 1 - wall
!>           ...
!>         END_CODES
!>       END_BOUNDARY_CONDITIONS
!>
!>       Injector types:
!>           GEOMETRY: CIRCLE, PARAMS=x_center, y_cenetr, z_center, radius, x_normal, y_normal, z_normal, number_of_particles
!>
!>           The normal has to point in the direction of injection.
!>
!>           The real number of particles injected by default is number_of_particles^2, 
!>           unless NPARTICLES: ASIS is specified, then number_of_particles is injected.
!>
!>           DISTRIBUTION parameter (works only for CIRCLE injector for 3D meshes) specifies 
!>           the spatial distrbution of the particles, by default it's uniform in polar
!>           coordinates UNIPO (as it was implmented). Specifying UNICA will produce uniform
!>           distribution in cartesian coordinates.
!>        
!>        Other injector types to be added.
!> @} 
!-----------------------------------------------------------------------

subroutine pts_reabcs()

  use def_kintyp
  use def_master
  use def_kermod
  use def_inpout
  use def_domain
  use def_partis
  use mod_opebcs
  use mod_memory
  use mod_ecoute,   only : ecoute
  use mod_opebcs,   only : boundary_conditions_read_boundary_codes
  use mod_opebcs,   only : opebcs_initialization_structure
  use mod_opebcs,   only : opebcs_initialization_variable
  use mod_messages, only : livinf
  use mod_ker_space_time_function
  use mod_pts_injection, only : PTS_INJ_SIZEDIST_CONST
  use mod_pts_injection, only : PTS_INJ_SIZEDIST_UNIFORM
  use mod_pts_injection, only : PTS_INJ_SIZEDIST_ROSINRAMMLER
  use mod_pts_injection, only : PTS_INJ_SIZEDIST_TABULATED
  use mod_pts_injection, only : PTS_INJ_SIZEDIST_FILED
  use mod_pts_injection, only : PTS_INJ_SIZEDIST_FILEM
  use mod_pts_injection, only : PTS_INJ_TEMPE_CONST
  use mod_pts_injection, only : PTS_INJ_TEMPE_FILE
  use mod_pts_injection, only : PTS_INJ_VELOC_ZERO    
  use mod_pts_injection, only : PTS_INJ_VELOC_FLUID   
  use mod_pts_injection, only : PTS_INJ_VELOC_NORMAL  
  use mod_pts_injection, only : PTS_INJ_VELOC_GAUSSIAN
  use mod_pts_injection, only : PTS_INJ_VELOC_CONIC   
  use mod_pts_injection, only : PTS_INJ_VELOC_SPRAY   
  use mod_pts_injection, only : PTS_INJ_VELOC_CONST
  use mod_pts_injection, only : PTS_INJ_VELOC_AXSPRAY   
  use mod_pts_injection, only : PTS_INJ_VELOC_FILE
  use mod_pts_injection, only : PTS_INJ_VELOC_SIZE_DEP_SPRAY 
  use mod_pts_injection, only : PTS_INJ_FLUCT_VELOC_CONST
  use mod_pts_injection, only : PTS_INJ_FLUCT_VELOC_UNIFORM
  use mod_pts_injection, only : PTS_INJ_FLUCT_VELOC_NORMAL
  use mod_pts_injection, only : PTS_INJ_FLOW_NMODIFIED 
  use mod_pts_injection, only : PTS_INJ_FLOW_NASIS     
  use mod_pts_injection, only : PTS_INJ_FLOW_MASSFLOW  
  use mod_pts_injection, only : PTS_INJ_FLOW_VOLUMEFLOW
  use mod_pts_injection, only : PTS_INJ_SPATDIST_UNICARTESIAN
  use mod_pts_injection, only : PTS_INJ_SPATDIST_UNIPOLAR    
  use mod_pts_injection, only : PTS_INJ_GEO_SQUARE    
  use mod_pts_injection, only : PTS_INJ_GEO_SPHERE    
  use mod_pts_injection, only : PTS_INJ_GEO_SEMISPHERE
  use mod_pts_injection, only : PTS_INJ_GEO_CIRCLE    
  use mod_pts_injection, only : PTS_INJ_GEO_RECTANGLE 
  use mod_pts_injection, only : PTS_INJ_GEO_POINT     
  use mod_pts_injection, only : PTS_INJ_GEO_CONE      
  use mod_pts_injection, only : PTS_INJ_GEO_RANDOM    
  use mod_pts_injection, only : PTS_INJ_GEO_SEGMENT   
  use mod_pts_injection, only : PTS_INJ_GEO_FILE      
  use mod_pts_injection, only : PTS_INJ_GEO_ANNULUS   
  use mod_pts_injection, only : pts_injection_memory
  implicit none
  character(len=8) :: fmt,num1,str2
  integer(ip)      :: iinj,ii,ifunc
  real(rp)         :: normlen
  character(5) :: wfname
  !
  ! Allocate memory for boundary codes
  !
  call opebcs_initialization_structure(1_ip,tbcod_pts)     ! Velocity 
  call opebcs_initialization_variable (1_ip,tbcod_pts)

  if( INOTSLAVE ) then

     kfl_imax_injector_pts          = 0
     ninj_pts                       = 0
     tinla_pts                      = 0.0_rp                   
     tfila_pts                      = 1.0e16_rp
     tpela_pts                      = 1.0e6_rp
     kfl_boundary_injection         = 0
     codbo_pts                      = 0_ip 
     !
     ! Reach section
     !
     call ecoute('pts_reabcs')
     do while( words(1) /= 'BOUND' )
        call ecoute('pts_reabcs')
     end do
     call ecoute('pts_reabcs')
     !
     ! Read boundary conditions field
     !
     do while( words(1) /= 'ENDBO' )

        if( words(1) == 'CODES'.and.exists('BOUND') ) then
           ! 
           ! User-defined codes on boundaries
           !          
           kfl_fixbo => kfl_fixbo_pts
           bvnat     => bvnat_pts
           tbcod     => tbcod_pts(1:)
           call boundary_conditions_read_boundary_codes('PARTICLES')
           !call reacod(2_ip)

        else if( words(1) == 'INJEC') then
           !
           ! Injection field
           !
           iinj = getint('INJEC',1_ip,'#INJECTOR')
           if( iinj > pts_minj ) call runend('PTS_REABCS: WRONG INJECTOR NUMBER, INCREASE PTS_MINJ IN DEF_PARTIS... SORRY FOR THAT!')
           !
           ! Keep track of highest injector index 
           !
           kfl_imax_injector_pts = max(kfl_imax_injector_pts,iinj)
           call ecoute('pts_reabcs')
           do while( words(1) /= 'ENDIN' )
              if( words(1) == 'GEOME' .or. words(1) == 'BOUND' ) then
                 !
                 ! Get boundary number
                 !
                 if( words(1) == 'BOUND' ) then
                    kfl_boundary_injection = 1_ip
                    codbo_pts(iinj) = getint('BOUND ',0_ip,'BOUNDARY CODE')
                 end if
                 !
                 ! Geometry of the injector
                 !
                 if( words(2) == 'SQUAR' ) then
                    injection_pts(iinj) % kfl_geometry = PTS_INJ_GEO_SQUARE
                    !
                    ! Initialize for old style input:
                    !
                    injection_pts(iinj) % num_part                 = int(param(2+7),ip)
                    injection_pts(iinj) % geo_coord_min(1:ndime)   = param(2+1:2+ndime)
                    injection_pts(iinj) % geo_coord_max(1:ndime)   = param(5+1:5+ndime)
                 else if( words(2) == 'SPHER' ) then
                    injection_pts(iinj) % kfl_geometry = PTS_INJ_GEO_SPHERE
                    !
                    ! Initialize for old style input:
                    !
                    injection_pts(iinj) % num_part                 = int(param(2+5),ip)
                    injection_pts(iinj) % geo_coord(1:ndime)       = param(2+1:2+ndime)
                    injection_pts(iinj) % geo_rad                  = param(2+4)
                 else if( words(2) == 'SEMIS' ) then
                    injection_pts(iinj) % kfl_geometry = PTS_INJ_GEO_SEMISPHERE
                    !
                    ! Initialize for old style input:
                    !
                    injection_pts(iinj) % num_part                 = int(param(2+8),ip)
                    injection_pts(iinj) % geo_coord(1:ndime)       = param(2+1:2+ndime)
                    injection_pts(iinj) % geo_rad                  = param(2+4)
                    injection_pts(iinj) % geo_normal(1:ndime)      = param(6+1:6+ndime)
                 else if( words(2) == 'CIRCL' ) then
                    injection_pts(iinj) % kfl_geometry = PTS_INJ_GEO_CIRCLE
                    !
                    ! Initialize for old style input:
                    !
                    injection_pts(iinj) % num_part                 = int(param(2+8),ip)
                    injection_pts(iinj) % geo_coord(1:ndime)       = param(2+1:2+ndime)
                    injection_pts(iinj) % geo_rad                  = param(2+4)
                    injection_pts(iinj) % geo_normal(1:ndime)      = param(6+1:6+ndime)
                 else if( words(2) == 'RECTA' ) then
                    injection_pts(iinj) % kfl_geometry = PTS_INJ_GEO_RECTANGLE
                    !
                    ! Initialize for old style input:
                    !
                    injection_pts(iinj) % num_part                 = int(param(2+10),ip)
                    injection_pts(iinj) % geo_coord1(1:ndime)      = param(2+1:2+ndime)
                    injection_pts(iinj) % geo_coord2(1:ndime)      = param(5+1:5+ndime)
                    injection_pts(iinj) % geo_coord3(1:ndime)      = param(8+1:8+ndime)
                 else if( words(2) == 'POINT' ) then
                    injection_pts(iinj) % kfl_geometry = PTS_INJ_GEO_POINT
                    !
                    ! Initialize for old style input:
                    !
                    injection_pts(iinj) % num_part                 = max(1_ip,int(param(2+4),ip))
                    injection_pts(iinj) % geo_coord(1:ndime)       = param(2+1:2+ndime)
                 else if( words(2) == 'CONE'  ) then
                    injection_pts(iinj) % kfl_geometry = PTS_INJ_GEO_CONE
                    !
                    ! Initialize for old style input:
                    !
                    injection_pts(iinj) % num_part                 = int(param(2+9),ip)
                    injection_pts(iinj) % geo_coord(1:ndime)       = param(2+1:2+ndime)
                    injection_pts(iinj) % geo_rad                  = param(2+4)
                    injection_pts(iinj) % geo_height               = param(2+5)
                    injection_pts(iinj) % geo_normal(1:ndime)      = param(7+1:7+ndime)
                 else if( words(2) == 'RANDO' ) then
                    injection_pts(iinj) % kfl_geometry = PTS_INJ_GEO_RANDOM
                    !
                    ! Initialize for old style input:
                    !
                    injection_pts(iinj) % num_part                 = int(param(2+7),ip)
                    injection_pts(iinj) % geo_coord_min(1:ndime)   = param(2+1:2+ndime)
                    injection_pts(iinj) % geo_coord_max(1:ndime)   = param(5+1:5+ndime)
                 else if( words(2) == 'SEGME' ) then
                    injection_pts(iinj) % kfl_geometry = PTS_INJ_GEO_SEGMENT
                    !
                    ! Initialize for old style input:
                    !
                    injection_pts(iinj) % num_part                 = int(param(2+7),ip)
                    injection_pts(iinj) % geo_coord_min(1:ndime)   = param(2+1:2+ndime)
                    injection_pts(iinj) % geo_coord_max(1:ndime)   = param(5+1:5+ndime)
                 else if( words(2) == 'ANNUL' ) then
                    injection_pts(iinj) % kfl_geometry = PTS_INJ_GEO_ANNULUS
                    !
                    ! Initialize for old style input:
                    !
                    injection_pts(iinj) % num_part                 = int(param(2+9),ip)
                    injection_pts(iinj) % geo_coord(1:ndime)       = param(2+1:2+ndime)
                    injection_pts(iinj) % geo_rad                  = param(2+4)
                    injection_pts(iinj) % geo_radmin               = param(2+5)
                    injection_pts(iinj) % geo_normal(1:ndime)      = param(7+1:7+ndime)
                 end if
                 
                 !
                 ! Input geometry parameters
                 ! 
                 if( exists('PARTI') ) injection_pts(iinj) % num_part         = getint('PARTI',0_ip,  '#PARTICLES')
                 if( exists('X    ') ) injection_pts(iinj) % geo_coord(1)     = getrea('X    ',0.0_rp,'Central coordinate x')
                 if( exists('Y    ') ) injection_pts(iinj) % geo_coord(2)     = getrea('Y    ',0.0_rp,'Central coordinate y')
                 if( exists('Z    ') ) injection_pts(iinj) % geo_coord(3)     = getrea('Z    ',0.0_rp,'Central coordinate z')
                 if( exists('XMIN ') ) injection_pts(iinj) % geo_coord_min(1) = getrea('XMIN ',0.0_rp,'Minimum coordinate x')
                 if( exists('YMIN ') ) injection_pts(iinj) % geo_coord_min(2) = getrea('YMIN ',0.0_rp,'Minimum coordinate y')
                 if( exists('ZMIN ') ) injection_pts(iinj) % geo_coord_min(3) = getrea('ZMIN ',0.0_rp,'Minimum coordinate z')
                 if( exists('XMAX ') ) injection_pts(iinj) % geo_coord_max(1) = getrea('XMAX ',0.0_rp,'Maximum coordinate x')
                 if( exists('YMAX ') ) injection_pts(iinj) % geo_coord_max(2) = getrea('YMAX ',0.0_rp,'Maximum coordinate y')
                 if( exists('ZMAX ') ) injection_pts(iinj) % geo_coord_max(3) = getrea('ZMAX ',0.0_rp,'Maximum coordinate z')
                 if( exists('NX   ') ) injection_pts(iinj) % geo_normal(1)    = getrea('NX   ',0.0_rp,'Normal vector x')
                 if( exists('NY   ') ) injection_pts(iinj) % geo_normal(2)    = getrea('NY   ',0.0_rp,'Normal vector y')
                 if( exists('NZ   ') ) injection_pts(iinj) % geo_normal(3)    = getrea('NZ   ',0.0_rp,'Normal vector z')
                 if( exists('RADIU') ) injection_pts(iinj) % geo_rad          = getrea('RADIU',0.0_rp,'Radius')
                 if( exists('MINRA') ) injection_pts(iinj) % geo_radmin       = getrea('MINRA',0.0_rp,'Minimum radius')
                 if( exists('HEIGH') ) injection_pts(iinj) % geo_height       = getrea('HEIGH',0.0_rp,'Height')

                 !
                 ! Normalize normal vector
                 !
                 normlen = sqrt(dot_product(injection_pts(iinj) % geo_normal,injection_pts(iinj) % geo_normal))
                 injection_pts(iinj) % geo_normal = injection_pts(iinj) % geo_normal/max(normlen,1.0e-16_rp) 


              else if ( words(1) == 'COORD' ) then
                 injection_pts(iinj) % num_part=getint('NUMBE',1_ip,'*PARTICLES')
                 call livinf(-9_ip,'INJECTED PARTICLES FROM FILE   ',injection_pts(iinj) % num_part )
                 injection_pts(iinj) % kfl_geometry = PTS_INJ_GEO_FILE
                 call pts_injection_memory(iinj,4_ip) 

                 call ecoute('pts_reabcs')
                 ii = 1
                 do while( words(1) /= 'ENDCO' )
                    if( ii > injection_pts(iinj) % num_part ) call runend('PTS_REABCS: WRONG NUMBER OF PARTICLES IN COORDINATE LIST')
                    injection_pts(iinj) % coord_list(1:ndime, ii ) = param(1:ndime)
                    call ecoute('pts_reabcs')
                    ii = ii + 1_ip
                 end do

              else if( words(1) == 'DISTR' ) then
                 if( words(2) == 'UNICA' ) then ! uniform distribution in cartesian coordinates
                    injection_pts(iinj) % kfl_geo_spatial_dist = PTS_INJ_SPATDIST_UNICARTESIAN
                 else if( words(2) == 'UNIPO') then ! uniform distribution in polar coordinates
                    injection_pts(iinj) % kfl_geo_spatial_dist = PTS_INJ_SPATDIST_UNIPOLAR
                 else
                    call runend("pts_reabcs: Particle distribution " // words(2) // " is not implemented")
                 end if
                  
              else if( words(1) == 'NPART' ) then
                 if ( words(2) == 'MASSF' .or. words(2) == 'VOLUM' ) then
                    !
                    ! Read type of mass flow 
                    !
                    if ( exists('MASSF') ) then
                       injection_pts(iinj) % kfl_flow = PTS_INJ_FLOW_MASSFLOW 
                    elseif ( exists('VOLUM') ) then
                       injection_pts(iinj) % kfl_flow = PTS_INJ_FLOW_VOLUMEFLOW
                    endif

                    if ( exists('CONST') ) injection_pts(iinj) % flow_rate = getrea('CONST',0.0_rp,'#Flow rate of injector')
                    if ( exists('FLOWR') ) injection_pts(iinj) % flow_rate = getrea('FLOWR',0.0_rp,'#Flow rate of injector')

                    if ( exists('SPACE') ) then
                       wfname                            = getcha('SPACE','NULL ','#Space/time Function name')
                       injection_pts(iinj) % kfl_flowfun = space_time_function_number(wfname)        
                    elseif ( exists('TIMEF') ) then
                       wfname                  = getcha('TIMEF','NULL ','#Time Function name')
                       do ifunc = 1,number_time_function
                          if( trim(wfname) == trim(time_function(ifunc) % name) ) then
                             injection_pts(iinj) % kfl_flowfun = -1_ip * ifunc
                          end if
                       end do
                       if (injection_pts(iinj) % kfl_flowfun == 0) &
                           call runend('pts_reabcs: time function: '//trim(wfname)//' is not defined in .ker.dat')
                    endif

                 elseif ( words(2) == 'ASIS ' ) then
                     injection_pts(iinj) % kfl_flow = PTS_INJ_FLOW_NASIS
                 else
                     injection_pts(iinj) % kfl_flow = PTS_INJ_FLOW_NMODIFIED
                 end if

              else if( (words(1) == 'STOCA') .or. (words(1) == 'STOCH') ) then
                 injection_pts(iinj) % kfl_random = 1
              else if( words(1) == 'PIPE' ) then
                 injection_pts(iinj) % kfl_random = 2

              else if( words(1) == 'TYPE ' ) then
                 !
                 ! On which type to apply
                 !
                 if( words(2) == 'ALL  ' ) then
                    injection_pts(iinj) % kfl_particle_type = 0
                 else
                    injection_pts(iinj) % kfl_particle_type =  getint('TYPE ',0_ip,'#TYPE ON WHICH TO APPLY INJECTOR')
                 end if
              
              else if( words(1) == 'INITI' ) then
                 !
                 ! Initial injection time
                 !
                 injection_pts(iinj) % time_initial = getrea('INITI',0.0_rp,'#INITIAL INJECTION TIME LAGRANGIAN PARTICLE')
              else if( words(1) == 'PERIO' ) then
                 !
                 ! Injection period
                 !
                 injection_pts(iinj) % time_period  = getrea('PERIO',0.0_rp,'#PERIOD INJECTION TIME LAGRANGIAN PARTICLE')
              else if( words(1) == 'FINAL' ) then
                 !
                 ! Final injection time
                 !
                 injection_pts(iinj) % time_final   = getrea('FINAL',1.0e16_rp,'#FINAL INJECTION TIME LAGRANGIAN PARTICLE')
                 
              else if( words(1) == 'VELOC' ) then
                 !
                 ! Injection velocity
                 !
                 if(     words(2) == 'ZERO ' ) then
                    injection_pts(iinj) % kfl_veloc = PTS_INJ_VELOC_ZERO
                 else if( words(2) == 'FLUID' ) then
                    injection_pts(iinj) % kfl_veloc = PTS_INJ_VELOC_FLUID
                 else if( words(2) == 'UNIFO' ) then
                    injection_pts(iinj) % kfl_veloc = PTS_INJ_VELOC_NORMAL
                    !
                    ! Initialize for old style input:
                    !
                    injection_pts(iinj) % vel_mag = param(3)
                 else if( words(2) == 'GAUSS' ) then
                    injection_pts(iinj) % kfl_veloc = PTS_INJ_VELOC_GAUSSIAN
                    !
                    ! Initialize for old style input:
                    !
                    injection_pts(iinj) % vel_mag   = param(3)
                    injection_pts(iinj) % vel_sigma = param(4)
                 else if( words(2) == 'CONIC' ) then
                    injection_pts(iinj) % kfl_veloc = PTS_INJ_VELOC_CONIC
                    !
                    ! Initialize for old style input:
                    !
                    injection_pts(iinj) % vel_mag   = param(3)
                    injection_pts(iinj) % vel_angle = param(4)
                 else if( words(2) == 'SPRAY' ) then
                    injection_pts(iinj) % kfl_veloc = PTS_INJ_VELOC_SPRAY
                    !
                    ! Initialize for old style input:
                    !
                    injection_pts(iinj) % vel_mag   = param(3)
                    injection_pts(iinj) % vel_angle = param(4)
                 else if( words(2) == 'CONST' ) then
                    injection_pts(iinj) % kfl_veloc = PTS_INJ_VELOC_CONST
                    !
                    ! Initialize for old style input:
                    !
                    injection_pts(iinj) % vel_vec(1) = param(3)
                    injection_pts(iinj) % vel_vec(2) = param(4)
                    injection_pts(iinj) % vel_vec(3) = param(5)
                 else if( words(2) == 'AXSPR' ) then
                    injection_pts(iinj) % kfl_veloc = PTS_INJ_VELOC_AXSPRAY
                    !
                    ! Initialize for old style input:
                    !
                    injection_pts(iinj) % vel_ax    = param(3)
                    injection_pts(iinj) % vel_angle = param(4)
                 else if( words(2) == 'SDSPR' ) then
                    injection_pts(iinj) % kfl_veloc = PTS_INJ_VELOC_SIZE_DEP_SPRAY
                    !
                    ! Initialize for old style input:
                    !
                    injection_pts(iinj) % vel_mag       = param(3)  ! largest droplets
                    injection_pts(iinj) % vel_diam_L    = param(4)
                    injection_pts(iinj) % vel_ang_max_L = param(5)  
                    injection_pts(iinj) % vel_ang_min_L = param(6)
                    injection_pts(iinj) % vel_diam_s    = param(7)
                    injection_pts(iinj) % vel_mag_s     = param(8)
                    injection_pts(iinj) % vel_ang_max_s = param(9)  
                    injection_pts(iinj) % vel_ang_min_s = param(10)
                 else if( words(2) == 'LIST ' ) then
                    injection_pts(iinj) % kfl_veloc = PTS_INJ_VELOC_FILE
                 end if

                 if ( exists('SPACE') ) then
                     wfname                 = getcha('SPACE','NULL ','#Space/time Function name')
                     injection_pts(iinj) % kfl_velfun = space_time_function_number(wfname)        
                 endif

                 if ( exists('VELOX') ) injection_pts(iinj) % vel_vec(1)   = getrea('VELOX',0.0_rp,'Droplet x velocity') 
                 if ( exists('VELOY') ) injection_pts(iinj) % vel_vec(2)   = getrea('VELOY',0.0_rp,'Droplet y velocity') 
                 if ( exists('VELOZ') ) injection_pts(iinj) % vel_vec(3)   = getrea('VELOZ',0.0_rp,'Droplet z velocity')  
                 if ( exists('VELMA') ) injection_pts(iinj) % vel_mag      = getrea('VELMA',0.0_rp,'Droplet velocity magnitude') 
                 if ( exists('VELAX') ) injection_pts(iinj) % vel_ax       = getrea('VELAX',0.0_rp,'Droplet axial velocity') 
                 if ( exists('SIGMA') ) injection_pts(iinj) % vel_sigma    = getrea('SIGMA',0.0_rp,'Droplet velocity sigma for Gaussian profile') 
                 if ( exists('ANGLE') ) injection_pts(iinj) % vel_angle    = getrea('ANGLE',0.0_rp,'Droplet velocity angle for cones') 

                 !
                 ! Size dependent velocity 
                 !
                 if ( exists('DIAML') ) injection_pts(iinj) % vel_diam_L    = getrea('DIAML',0.0_rp,'Diameter of large droplets') 
                 if ( exists('VMAGL') ) injection_pts(iinj) % vel_mag       = getrea('VMAGL',0.0_rp,'Velocity magnitude of large droplets') 
                 if ( exists('AMAXL') ) injection_pts(iinj) % vel_ang_max_L = getrea('AMAXL',0.0_rp,'Outer angle of large droplets') 
                 if ( exists('AMINL') ) injection_pts(iinj) % vel_ang_min_L = getrea('AMINL',0.0_rp,'Inner angle of large droplets') 
                 if ( exists('DIAMS') ) injection_pts(iinj) % vel_diam_s    = getrea('DIAMS',0.0_rp,'Diameter of small droplets') 
                 if ( exists('VMAGS') ) injection_pts(iinj) % vel_mag_s     = getrea('VMAGS',0.0_rp,'Velocity magnitude of small droplets') 
                 if ( exists('AMAXS') ) injection_pts(iinj) % vel_ang_max_s = getrea('AMAXS',0.0_rp,'Outer angle of small droplets') 
                 if ( exists('AMINS') ) injection_pts(iinj) % vel_ang_min_s = getrea('AMINS',0.0_rp,'Inner angle of small droplets') 


                 !
                 ! Velocity fluctuations
                 !
                 if ( exists('FLUCT') ) then
                    wfname = getcha('FLUCT','CONST','#Velocity fluctuation type')
                    if (wfname == "CONST") injection_pts(iinj) % kfl_fluct_veloc = PTS_INJ_FLUCT_VELOC_CONST 
                    if (wfname == "UNIFO") injection_pts(iinj) % kfl_fluct_veloc = PTS_INJ_FLUCT_VELOC_UNIFORM
                    if (wfname == "NORMA") injection_pts(iinj) % kfl_fluct_veloc = PTS_INJ_FLUCT_VELOC_NORMAL 
                    
                    if (     injection_pts(iinj) % kfl_fluct_veloc == PTS_INJ_FLUCT_VELOC_UNIFORM &
                        .or. injection_pts(iinj) % kfl_fluct_veloc == PTS_INJ_FLUCT_VELOC_NORMAL  ) then 
                       injection_pts(iinj) % fluct_vel_std = getrea('FLSTD',0.0_rp,'Normalized standard deviation around velocity magnitude.')
                    endif

                 endif

                 !
                 ! Read listed velocity values
                 !
                 if (injection_pts(iinj) % kfl_veloc == PTS_INJ_VELOC_FILE) then
                    call pts_injection_memory(iinj,3_ip) 
                    call ecoute('pts_reabcs')
                    ii = 1
                    do while( words(1) /= 'ENDVE' )
                       if( ii > injection_pts(iinj) % num_part ) call runend('PTS_REABCS: WRONG NUMBER OF PARTICLES IN VELOCITY LIST')
                       injection_pts(iinj) % vel_list(1:ndime, ii ) = param(1:ndime)
                       call ecoute('pts_reabcs')
                       ii = ii + 1_ip
                    end do
                 endif

              else if( words(1) == 'TEMPE' ) then
                 !
                 ! Injection tempreature
                 !
                 if( words(2) == 'CONST' ) then
                    injection_pts(iinj) % kfl_tempe     = PTS_INJ_TEMPE_CONST 
                    !
                    ! Initialize for old style input:
                    !
                    injection_pts(iinj) % tempe         = param(3)
                 elseif ( words(2) == 'LIST ' ) then
                    injection_pts(iinj) % kfl_tempe     = PTS_INJ_TEMPE_FILE 
                 end if

                 if ( exists('VALUE') ) then 
                     injection_pts(iinj) % tempe = getrea('VALUE',injection_pts(iinj) % tempe,'Droplet temperature') 
                 else if ( exists('CONST') ) then
                     injection_pts(iinj) % tempe = getrea('CONST',injection_pts(iinj) % tempe,'Droplet temperature') 
                 else if ( exists('TEMPE') ) then
                     injection_pts(iinj) % tempe = getrea('TEMPE',injection_pts(iinj) % tempe,'Droplet temperature') 
                 endif

                 !
                 ! Read listed temperature values
                 !
                 if (injection_pts(iinj) % kfl_tempe == PTS_INJ_TEMPE_FILE) then
                    call pts_injection_memory(iinj,2_ip) 
                    call ecoute('pts_reabcs')
                    ii = 1
                    do while( words(1) /= 'ENDTE' )
                       if( ii > injection_pts(iinj) % num_part ) call runend('PTS_REABCS: WRONG NUMBER OF PARTICLES IN TEMPERATURE LIST')
                       injection_pts(iinj) % tempe_list(ii)  = param(1)
                       call ecoute('pts_reabcs')
                       ii = ii + 1_ip
                    end do
                 endif
                 

              else if( words(1) == 'SIZED' ) then 
                 !
                 ! Size distribution
                 !
                 if(      words(2) == 'CONST' ) then
                    injection_pts(iinj) % kfl_size_dist     =  PTS_INJ_SIZEDIST_CONST
                 else if( words(2) == 'UNIFO' ) then
                    injection_pts(iinj) % kfl_size_dist     =  PTS_INJ_SIZEDIST_UNIFORM
                    !
                    ! Initialize for old style input:
                    !
                    injection_pts(iinj) % size_dmin = param(3)
                    injection_pts(iinj) % size_dmax = param(4)
                 else if( words(2) == 'ROSIN' ) then
                    injection_pts(iinj) % kfl_size_dist     =  PTS_INJ_SIZEDIST_ROSINRAMMLER
                    !
                    ! Initialize for old style input:
                    !
                    injection_pts(iinj) % size_dmin     = param(3)
                    injection_pts(iinj) % size_dmax     = param(4)
                    injection_pts(iinj) % size_rr_dbar  = param(5)
                    injection_pts(iinj) % size_rr_n     = param(6)
                 else if( words(2) == 'CDFFR' ) then
                    injection_pts(iinj) % kfl_size_dist     =  PTS_INJ_SIZEDIST_TABULATED
                    injection_pts(iinj) % kfl_size_lookupfw =  getint('CDFFR',1_ip,'#Index of cummulative density function framework') 
                 elseif ( words(2) == 'LISTD' ) then
                    injection_pts(iinj) % kfl_size_dist     = PTS_INJ_SIZEDIST_FILED 
                 elseif ( words(2) == 'LISTM' ) then
                    injection_pts(iinj) % kfl_size_dist     = PTS_INJ_SIZEDIST_FILEM
                 end if

                 !
                 ! Process named parameters if they are present
                 !
                 if ( exists('DIAME') ) injection_pts(iinj) % size_diame   = getrea('DIAME',0.0_rp,'Droplet diameter') 
                 if ( exists('DMINI') ) injection_pts(iinj) % size_dmin    = getrea('DMINI',0.0_rp,'Minimum droplet diameter') 
                 if ( exists('DMAXI') ) injection_pts(iinj) % size_dmax    = getrea('DMAXI',0.0_rp,'Maximum droplet diameter') 

                 if ( exists('DBARR') ) injection_pts(iinj) % size_rr_dbar = getrea('DBARR',0.0_rp,'Rosin-Rammler distribuition mean diameter parameter') 
                 if ( exists('NROSI') ) injection_pts(iinj) % size_rr_n    = getrea('NROSI',0.0_rp,'Rosin-Rammler distribuition width parameter') 

                 if ( injection_pts(iinj) % kfl_size_dist == PTS_INJ_SIZEDIST_ROSINRAMMLER ) then
                    !
                    ! Precompute Rossin-Rammler correction factor
                    !
                    injection_pts(iinj) % size_rr_dbar = injection_pts(iinj) % size_rr_dbar - injection_pts(iinj) % size_dmin  
                    injection_pts(iinj) % size_rr_K    = 1.0_rp - exp( -1.0_rp*((injection_pts(iinj) % size_dmax - injection_pts(iinj) % size_dmin)/injection_pts(iinj) % size_rr_dbar)**injection_pts(iinj) % size_rr_n )
                 endif

                 !
                 ! Read listed size values
                 !
                 if (injection_pts(iinj) % kfl_size_dist == PTS_INJ_SIZEDIST_FILED .or. &
                    &injection_pts(iinj) % kfl_size_dist == PTS_INJ_SIZEDIST_FILEM) then
                    call pts_injection_memory(iinj,1_ip) 
                    call ecoute('pts_reabcs')
                    ii = 1
                    do while( words(1) /= 'ENDSI' )
                       if( ii > injection_pts(iinj) % num_part ) call runend('PTS_REABCS: WRONG NUMBER OF PARTICLES IN SIZE LIST')
                       injection_pts(iinj) % size_list(ii)  = param(1)
                       call ecoute('pts_reabcs')
                       ii = ii + 1_ip
                    end do
                 endif


              end if
              call ecoute('pts_reabcs')
           end do

        else if( words(1)(1:4) == 'INJE') then
           !
           ! Inject particles with a saqure, sphere, semi-sphere or a circle
           !
           fmt = '(I1.1)'
           num1 = words(1)(5:5)
           read (num1, fmt) iinj
           kfl_imax_injector_pts = max(kfl_imax_injector_pts,iinj)
           if( iinj > 0 .and. iinj <= pts_minj) then
              if( words(2) == 'SQUAR' ) then
                 injection_pts(iinj) % kfl_geometry = PTS_INJ_GEO_SQUARE
                 !
                 ! Initialize for old style input:
                 !
                 injection_pts(iinj) % num_part                 = int(param(2+7),ip)
                 injection_pts(iinj) % geo_coord_min(1:ndime)   = param(2+1:2+ndime)
                 injection_pts(iinj) % geo_coord_max(1:ndime)   = param(5+1:5+ndime)
              else if( words(2) == 'SPHER' ) then
                 injection_pts(iinj) % kfl_geometry = PTS_INJ_GEO_SPHERE
                 !
                 ! Initialize for old style input:
                 !
                 injection_pts(iinj) % num_part                 = int(param(2+5),ip)
                 injection_pts(iinj) % geo_coord(1:ndime)       = param(2+1:2+ndime)
                 injection_pts(iinj) % geo_rad                  = param(2+4)
              else if( words(2) == 'SEMIS' ) then
                 injection_pts(iinj) % kfl_geometry = PTS_INJ_GEO_SEMISPHERE
                 !
                 ! Initialize for old style input:
                 !
                 injection_pts(iinj) % num_part                 = int(param(2+8),ip)
                 injection_pts(iinj) % geo_coord(1:ndime)       = param(2+1:2+ndime)
                 injection_pts(iinj) % geo_rad                  = param(2+4)
                 injection_pts(iinj) % geo_normal(1:ndime)      = param(6+1:6+ndime)
              else if( words(2) == 'CIRCL' ) then
                 injection_pts(iinj) % kfl_geometry = PTS_INJ_GEO_CIRCLE
                 !
                 ! Initialize for old style input:
                 !
                 injection_pts(iinj) % num_part                 = int(param(2+8),ip)
                 injection_pts(iinj) % geo_coord(1:ndime)       = param(2+1:2+ndime)
                 injection_pts(iinj) % geo_rad                  = param(2+4)
                 injection_pts(iinj) % geo_normal(1:ndime)      = param(6+1:6+ndime)
                 !
                 ! Normalize normal vector
                 !
                 normlen = sqrt(dot_product(injection_pts(iinj) % geo_normal,injection_pts(iinj) % geo_normal))
                 injection_pts(iinj) % geo_normal = injection_pts(iinj) % geo_normal/normlen 
              else if( words(2) == 'RECTA' ) then
                 injection_pts(iinj) % kfl_geometry = PTS_INJ_GEO_RECTANGLE
                 !
                 ! Initialize for old style input:
                 !
                 injection_pts(iinj) % num_part                 = int(param(2+10),ip)
                 injection_pts(iinj) % geo_coord1(1:ndime)      = param(2+1:2+ndime)
                 injection_pts(iinj) % geo_coord2(1:ndime)      = param(5+1:5+ndime)
                 injection_pts(iinj) % geo_coord3(1:ndime)      = param(8+1:8+ndime)
              else if( words(2) == 'POINT' ) then
                 injection_pts(iinj) % kfl_geometry = PTS_INJ_GEO_POINT
                 !
                 ! Initialize for old style input:
                 !
                 injection_pts(iinj) % num_part                 = max(1_ip,int(param(2+4),ip))
                 injection_pts(iinj) % geo_coord(1:ndime)       = param(2+1:2+ndime)
              else if( words(2) == 'CONE'  ) then
                 injection_pts(iinj) % kfl_geometry = PTS_INJ_GEO_CONE
                 !
                 ! Initialize for old style input:
                 !
                 injection_pts(iinj) % num_part                 = int(param(2+9),ip)
                 injection_pts(iinj) % geo_coord(1:ndime)       = param(2+1:2+ndime)
                 injection_pts(iinj) % geo_rad                  = param(2+4)
                 injection_pts(iinj) % geo_height               = param(2+5)
                 injection_pts(iinj) % geo_normal(1:ndime)      = param(7+1:7+ndime)
              else if( words(2) == 'RANDO' ) then
                 injection_pts(iinj) % kfl_geometry = PTS_INJ_GEO_RANDOM
                 !
                 ! Initialize for old style input:
                 !
                 injection_pts(iinj) % num_part                 = int(param(2+7),ip)
                 injection_pts(iinj) % geo_coord_min(1:ndime)   = param(2+1:2+ndime)
                 injection_pts(iinj) % geo_coord_max(1:ndime)   = param(5+1:5+ndime)
              else if( words(2) == 'SEGME' ) then
                 injection_pts(iinj) % kfl_geometry = PTS_INJ_GEO_SEGMENT
                 !
                 ! Initialize for old style input:
                 !
                 injection_pts(iinj) % num_part                 = int(param(2+7),ip)
                 injection_pts(iinj) % geo_coord_min(1:ndime)   = param(2+1:2+ndime)
                 injection_pts(iinj) % geo_coord_max(1:ndime)   = param(5+1:5+ndime)
              else if( words(2) == 'ANNUL' ) then
                 injection_pts(iinj) % kfl_geometry = PTS_INJ_GEO_ANNULUS
                 !
                 ! Initialize for old style input:
                 !
                 injection_pts(iinj) % num_part                 = int(param(2+9),ip)
                 injection_pts(iinj) % geo_coord(1:ndime)       = param(2+1:2+ndime)
                 injection_pts(iinj) % geo_rad                  = param(2+4)
                 injection_pts(iinj) % geo_radmin               = param(2+5)
                 injection_pts(iinj) % geo_normal(1:ndime)      = param(7+1:7+ndime)
                 !
                 ! Normalize normal vector
                 !
                 normlen = sqrt(dot_product(injection_pts(iinj) % geo_normal,injection_pts(iinj) % geo_normal))
                 injection_pts(iinj) % geo_normal = injection_pts(iinj) % geo_normal/normlen 
              end if
           else
              write(str2,fmt) pts_minj
              call runend("NUMBER OF INJECTIONS MUST BE BETWEEN 1 AND "//str2)
           end if


        else if( words(1) == 'INITI' ) then
           !
           ! Initial injection time
           !
           tinla_pts = getrea('INITI',1.0_rp,'#INITIAL INJECTION TIME LAGRANGIAN PARTICLE')
        else if( words(1) == 'FINAL' ) then
           !
           ! Default final injection time
           !
           tfila_pts = getrea('FINAL',1.0_rp,'#FINAL INJECTION TIME LAGRANGIAN PARTICLE')
        else if( words(1) == 'PERIO' ) then
           !
           ! Default injection period
           !
           tpela_pts = getrea('PERIO',1.0_rp,'#PERIOD INJECTION TIME LAGRANGIAN PARTICLE')
        else if( words(1) == 'VELOC' ) then
           !
           ! Default injection velocity
           !
           call runend('PTS_REABCS: OBSOLETE WAY OF DEFINING INJECTION VELOCITY')

        end if

        call ecoute('pts_reabcs')
     end do

  end if

end subroutine pts_reabcs
