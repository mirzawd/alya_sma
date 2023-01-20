!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine pts_reanut()
  !------------------------------------------------------------------------
  !****f* Partis/pts_reanut
  ! NAME 
  !    pts_reanut
  ! DESCRIPTION
  !    This routine reads data
  ! OUTPUT
  ! USES
  ! USED BY
  !    Reapro
  !***
  !------------------------------------------------------------------------
  use def_kintyp
  use def_master
  use def_kermod
  use def_inpout
  use def_domain
  use def_partis
  use mod_messages, only : messages_live
  use mod_ecoute, only :  ecoute
  implicit none
  integer(ip) :: ivari
  logical(lg) :: reset_time_step

  if( INOTSLAVE ) then

     kfl_adapt_pts =  0                  ! Constant time step
     kfl_usbin_pts =  0                  ! Use element bin
     kfl_order_pts =  2                  ! Second order AB
     kfl_walld_pts =  0                  ! Do not use walld for checking deposition
     kfl_vesgs_pts =  0                  ! Velocity sgs 
     kfl_thermo_timsch_pts =  1          ! Implicit time integration for mass and energy
     gamma_pts     =  0.5_rp             ! Newmark's gamma_pts constant
     beta_pts      =  1.0_rp / 4.0_rp    ! Newmark's beta constant
     chale_pts     = -1.0_rp             ! Characteristic length
     safet_pts     =  1.0_rp             ! Safety factor
     dtime_pts     =  1.0_rp             ! Minimum time step
     gamma_pts     =  0.75_rp            ! Newmark gamma
     beta_pts      =  0.390625_rp        ! Newmark beta
     
     reset_time_step = .false.

     ivari         = 1
     
     !-------------------------------------------------------------
     !
     ! Numerical problem
     !
     !-------------------------------------------------------------

     call ecoute('ker_reanut')
     do while( words(1) /= 'NUMER' )
        call ecoute('ker_reanut')
     end do
     call ecoute('ker_reanut')

     do while(words(1) /= 'ENDNU' )

        if(      words(1) == 'SUBGR' ) then
           !
           ! Subgris scale tracking
           !
           if( exists('VELOC') ) kfl_vesgs_pts = 1
           
        else if( words(1) == 'NEWMA' ) then
           !
           ! Time integration Neumark coefficients
           !
           gamma_pts = getrea('GAMMA',0.5_rp, '#GAMMA FOR NEWMARK')
           beta_pts  = getrea('BETA ',0.25_rp,'#BETA FOR NEWMARK')

        else if( words(1) == 'ORDER' ) then
           !
           ! Time integration order
           !
           kfl_order_pts = getint('ORDER',2_ip, '#TIME INTEGRATION ORDER')
           
        else if( words(1) == 'WALLD' ) then
           !
           ! PDE based wall distance
           !
           if( option('WALLD') )    kfl_walld_pts = 1
           if( words(2) == 'PDE  ') kfl_walld_pts = 1
           
        else if( words(1) == 'ELEME' ) then
           !
           ! ELement search strategy
           !
           if( words(2) == 'BIN  ' ) then
              kfl_usbin_pts = 1 
              !if( if_moving_mesh_pts == 1)  call runend('PARTIS: BECAUSE OF MOVING_MESH BIN ELEMENT SEARCH STRATEGY IS NOT AVAILABLE')

           else if( words(2) == 'NEIGH' ) then
              kfl_usbin_pts = 0
           end if
           
        else if( words(1) == 'SLIPW' ) then
           !
           ! Slip wall distance solver
           !
           ivari = 1
           call ecoute('pts_reanut')
           do while( words(1) /= 'ENDSL' )
              if( words(1) == 'ALGEB' ) then
                 solve_sol => solve(1:)
                 call reasol(1_ip)
              end if
              call ecoute('pts_reanut')
           end do
           
        else if( words(1) == 'BOUNC' ) then
           !
           ! Bouncing wall distance solver
           !
           ivari = 2
           call ecoute('pts_reanut')
           do while( words(1) /= 'ENDBO' )
              if( words(1) == 'ALGEB' ) then
                 solve_sol => solve(2:)
                 call reasol(1_ip)
              end if
              call ecoute('pts_reanut')
           end do
           
        else if( words(1) == 'ALGEB' ) then
           !
           ! Solver for slip wall
           !
           solve_sol => solve(ivari:)
           call reasol(1_ip)
           ivari = 1

        else if( words(1) == 'ADAPT' ) then
           !
           ! Time adaptation scheme
           !
           if( exists('CHARA') ) then
              if( exists('DIAME') ) then
                 chale_pts =  0.0_rp
              else if( exists('ELEME') ) then
                 chale_pts = -1.0_rp
              else if( exists('VISCO') ) then
                 chale_pts = -2.0_rp
              else if( exists('TAU  ') ) then
                 chale_pts = -3.0_rp
              else
                 chale_pts =  getrea('CHARA',1.0_rp,'#CHARACTERISTIC LENGTH')
              end if
           end if

           if( exists('SAFET') ) then
              safet_pts =  getrea('SAFET',1.0_rp,'#SAFETY FACTOR')
           end if

           if( words(2) == 'ERROR' ) then
              kfl_adapt_pts =  1
           else if( words(2) == 'POSIT' ) then
              kfl_adapt_pts =  2
           else if( words(2) == 'FROMC' ) then
              kfl_adapt_pts =  3
           else if( words(2) == 'VELOC' ) then
              kfl_adapt_pts =  4
           end if
                      
           reset_time_step = .true.
           
        else if( words(1) == 'TIMES' ) then
           !
           ! Time step. The time step of the particles will be the minimum between
           ! the extrnal one and DTMIN_PTS
           !
           kfl_adapt_pts = -1
           dtime_pts     = getrea('TIMES',1.0_rp,'#TIME STEP')
           reset_time_step = .true.

        else if( words(1) == 'THERM' ) then
           !
           ! Thermodynamic time scheme
           !
           if( exists('EXPL1') ) kfl_thermo_timsch_pts = 0
           if( exists('IMPL1') ) kfl_thermo_timsch_pts = 1
           
        end if

        call ecoute('ker_reanut')
     end do
     !
     ! Set all particle to global options
     !
     if( reset_time_step ) then
        call messages_live('TIME STEP STRATEGY IS GLOBAL AND OVERWRITES THE ONES OF THE TYPES','WARNING')
        parttyp(:) % kfl_tstep = kfl_adapt_pts
        parttyp(:) % dtime     = dtime_pts        
        parttyp(:) % safet     = safet_pts        
        parttyp(:) % chale     = chale_pts        
     end if
     
  end if
end subroutine pts_reanut
