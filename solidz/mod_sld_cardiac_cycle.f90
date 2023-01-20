!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-------------------------------------------------------------------------------
!> @addtogroup variablepressure boundary conditions
!> @{
!> @authors Adria Quintanas-Corominas : adria.quintanas@bsc.es
!> @authors Alfonso Santiago          : alfonso.santiago@bsc.es
!> @date    May, 2020
!> @}
!------------------------------------------------------------------------------
module mod_sld_cardiac_cycle
  ! ----------------------------------------
  use def_kintyp_basic,               only : ip, rp, lg
  use mod_bouder
  use mod_eccoupling,                 only : eccou_do_volume_integral
  use def_master,                     only : mem_modul, modul

  ! ----------------------------------------
  implicit none
  ! ----------------------------------------
  integer(ip),   private, parameter       :: max_cycles = 12_ip !TODO: make it caclulatable or user specified
  integer(ip),   private, parameter       :: max_valves = 10_ip !TODO: make calculatable
  integer(ip),   private, parameter       :: max_cavities = 4_ip !TODO: make calculatable
  integer(ip),   parameter                :: &
         PHASE_TRANSITION_DVOL = 1_ip,       &
         PHASE_TRANSITION_DCAL = 2_ip
  integer(ip),   private, parameter       :: &
         VOL_METHOD_PLANE_POINT   = 1_ip,    &
         VOL_METHOD_PLANE_NODESET = 2_ip,    &
         VOL_METHOD_DIVERGENCE    = 3_ip,    & !Volume calculation based on divergence theorem (slower)
         VOL_METHOD_MIXEDPROD     = 4_ip       !Volume calculation based on mixed product of vertices (faster) 

  ! ----------------------------------------
  integer(ip),   protected                :: ncavi = 0_ip       !number of cavities in the problem
  integer(ip),   protected                :: ncycles = 0_ip     !current number of cycles, counter
  logical(lg),   protected                :: kfl_cardiac_cycle  = .False.
  logical(lg),   protected                :: kfl_heart_cavity   = .False.
  logical(lg),   protected                :: kfl_pressure_cycle = .False.
  real(rp),      protected                :: volcai(3) = 0.0_rp !(t,c_prev,c_current) : t- timestep number at which calcium was calculated, c - previous and current calcium
  integer(ip),   protected                :: phase41_trans_method = PHASE_TRANSITION_DVOL
  ! ----------------------------------------

  type cavity_hole
     real(rp)                               :: centroid(3) !this gets modified on every timestep
     integer(ip)                            :: number_of_nodes !local to the partition
     integer(ip), dimension(:,:), pointer   :: edges => null() !local to the partition
  end type

  type heart_cavity
     logical(lg)                            :: initialized = .False.
     integer(ip)                            :: bouset ! codbo for the cavity 
     integer(ip), dimension(max_valves)     :: valset ! codbos for each of the valves
     integer(ip)                            :: nvalves ! number of valves calculated from .dat
     integer(ip), dimension(3)           :: plane_nodset
     real(rp),    dimension(3)           :: plane_origin
     real(rp),    dimension(3)           :: plane_normal
     real(rp),    dimension(3)           :: plane_vector
     real(rp),    dimension(3,3)         :: plane_points
     real(rp),    dimension(5)           :: volume
     real(rp),    dimension(3)           :: dvol
     real(rp),    dimension(3)           :: ddvol
     real(rp)                            :: ini_vol
     real(rp)                            :: end_dia_vol
     real(rp)                            :: end_sys_vol
     real(rp),    dimension(3)           :: prest = 0.0_rp
     real(rp),    dimension(3)           :: pres = 0.0_rp
     integer(ip)                         :: phase = 0_ip
     integer(ip)                         :: phase_counter = 0_ip
     integer(ip)                         :: volume_method = 0_ip
     !------------------------------
     ! Holes related variables
     type(cavity_hole), dimension(:), pointer :: valves
     integer(ip)                         :: hole_field ! node field id that has hole node labels
   contains
     procedure,   pass                   :: init => init_heart_cavity
  end type heart_cavity
  ! ----------------------------------------
  type windkessel
     real(rp)                            :: r
     real(rp)                            :: c
     real(rp)                            :: p_ini
     real(rp), dimension(3)              :: pres=0.0_rp
   contains
     procedure,   pass                   :: init => init_windkessel
  end type windkessel
  ! ----------------------------------------
  type cycle_control
     logical(lg)                         :: defined = .False.
     logical(lg)                         :: active = .False.
     integer(ip)                         :: cycle_id
     integer(ip)                         :: cavity_controled
     integer(ip)                         :: n_beats = 0_ip
     integer(ip)                         :: max_beats = 0_ip
     integer(ip)                         :: prevcycle = 0_ip
     integer(ip)                         :: nextcycle = 0_ip
     real(rp)                            :: last_phase_change=0.0_rp
     real(rp)                            :: tzero
     real(rp)                            :: tpstr
     real(rp)                            :: pzero
     real(rp)                            :: pstr0
     real(rp), dimension(2)              :: gains_contraction=0.0_rp
     type(windkessel)                    :: assoc_wdk
     real(rp), dimension(2)              :: gains_relaxation=0.0_rp
     real(rp)                            :: ppost
   contains
     procedure,   pass                   :: init => init_cycle_control
  end type cycle_control
  ! ----------------------------------------
  type(heart_cavity),  protected          :: cavities(max_cavities)
  type(cycle_control), protected          :: cycles(max_cycles)
  ! ----------------------------------------

  private :: compute_volume_plane
  private :: compute_volume_holes
  public  :: cardiac_cycle_redist

contains 

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-12-22
  !> @brief   Initialization
  !> @details Initialization of classes
  !> 
  !-----------------------------------------------------------------------

  subroutine init_cycle_control(self)
    class(cycle_control), intent(inout) :: self

    self % defined           = .False.
    self % active            = .False.
    self % cycle_id          = 0_ip
    self % cavity_controled  = 0_ip
    self % n_beats           = 0_ip
    self % max_beats         = 0_ip
    self % prevcycle         = 0_ip
    self % nextcycle         = 0_ip
    self % last_phase_change = 0.0_rp
    self % tzero             = 0.0_rp
    self % tpstr             = 0.0_rp
    self % pzero             = 0.0_rp
    self % pstr0             = 0.0_rp
    self % gains_contraction = 0.0_rp
    self % gains_relaxation  = 0.0_rp
    self % ppost             = 0.0_rp

    call self % assoc_wdk % init()

  end subroutine init_cycle_control

  subroutine init_heart_cavity(self)
    class(heart_cavity), intent(inout) :: self

    self % initialized   = .False.
    self % bouset        = 0_ip
    self % valset        = 0_ip
    self % plane_nodset  = 0_ip
    self % plane_origin  = 0.0_rp
    self % plane_normal  = 0.0_rp
    self % plane_vector  = 0.0_rp
    self % plane_points  = 0.0_rp
    self % volume        = 0.0_rp
    self % dvol          = 0.0_rp
    self % ddvol         = 0.0_rp
    self % ini_vol       = 0.0_rp
    self % end_dia_vol   = 0.0_rp
    self % end_sys_vol   = 0.0_rp
    self % prest         = 0.0_rp
    self % pres          = 0.0_rp
    self % phase         = 0_ip
    self % phase_counter = 0_ip
    self % nvalves       = 0_ip
    self % hole_field    = 0_ip  
    nullify( self % valves )
  end subroutine init_heart_cavity

  subroutine init_windkessel(self)
    class(windkessel), intent(inout) :: self

    self % r     = 0.0_rp
    self % c     = 0.0_rp
    self % p_ini = 0.0_rp
    self % pres  = 0.0_rp

  end subroutine init_windkessel

  ! ---------------------------------------------------------------------
  !> @author   Adria  Quintanas-Corominas (adria.quintanas@bsc.es)
  !> @author   Alfonso Santiago (alfonso.santiago@bsc.es)
  ! ---------------------------------------------------------------------
  subroutine sld_cardiac_cycle_initialization ( )

    integer(ip) :: icavi, icycle

    do icavi = 1,size(cavities)
       call cavities(icavi) % init()
    end do
    do icycle = 1,size(cycles)
       call cycles(icycle) % init()
    end do

  end subroutine sld_cardiac_cycle_initialization

  ! ---------------------------------------------------------------------
  !> @author   Adria  Quintanas-Corominas (adria.quintanas@bsc.es)
  !> @author   Alfonso Santiago (alfonso.santiago@bsc.es)
  ! ---------------------------------------------------------------------
  subroutine sld_cardiac_cycle_initialise_variables ( &
       )
    ! -------------------------------
    implicit none
    ! -------------------------------
    kfl_cardiac_cycle  = .False.
    kfl_heart_cavity   = .False.
    kfl_pressure_cycle = .False. 
    ncavi              = 0_ip            ! number of cavities (ventricles)
    ncycles            = 0_ip            ! number of cavities (ventricles)
    volcai                = 0.0_rp
    phase41_trans_method  = PHASE_TRANSITION_DVOL
    ! -------------------------------

  end subroutine sld_cardiac_cycle_initialise_variables

  ! ---------------------------------------------------------------------
  !> @author   Adria  Quintanas-Corominas (adria.quintanas@bsc.es)
  !> @author   Alfonso Santiago (alfonso.santiago@bsc.es)
  ! ---------------------------------------------------------------------
  subroutine sld_cardiac_cycle_write_res( &
       itask )
    use def_solidz, only            :  lun_carcy_res_sld
    use def_master, only            :  ittim, cutim, intost
    use def_master, only            :  TIME_N
    ! -------------------------------
    implicit none
    ! -------------------------------
    integer(ip), intent(in)        :: itask
    real(rp)                       :: auxr(max_cavities)
    integer(ip)                    :: icycle, cycid(max_cavities), i, icavi
    ! -------------------------------
    select case(itask)

    case( 1_ip )

       write(lun_carcy_res_sld, "(100a)") '# --| ALYA Cardiac Cycle Computations'
       write(lun_carcy_res_sld, "(100a)") '# --| '
       write(lun_carcy_res_sld, "(100a)") '# --|  1. Step increment'
       write(lun_carcy_res_sld, "(100a)") '# --|  2. Step time'
       write(lun_carcy_res_sld, "(100a)") '# --|  3. Overall calcium integral'


       i = 4_ip
       do icavi = 1_ip, max_cavities
         if ( cavities(icavi) % initialized ) then
            write(lun_carcy_res_sld,"(100a)") '# --|  '//trim(intost(i))//'. Cavity #'//trim(intost(icavi))
            i=i+1_ip
            write(lun_carcy_res_sld,"(100a)") '# --|  '//trim(intost(i))//'. Cavity #'//trim(intost(icavi))//': cycle'
            i=i+1_ip
            write(lun_carcy_res_sld,"(100a)") '# --|  '//trim(intost(i))//'. Cavity #'//trim(intost(icavi))//': phase'
            i=i+1_ip
            write(lun_carcy_res_sld,"(100a)") '# --|  '//trim(intost(i))//'. Cavity #'//trim(intost(icavi))//': volume'
            i=i+1_ip
            write(lun_carcy_res_sld,"(100a)") '# --|  '//trim(intost(i))//'. Cavity #'//trim(intost(icavi))//': cavity pressure'
            i=i+1_ip
            write(lun_carcy_res_sld,"(100a)") '# --|  '//trim(intost(i))//'. Cavity #'//trim(intost(icavi))//': outflow windkessel pressure'
            i=i+1_ip
         end if
       end do
            
       write(lun_carcy_res_sld,"(a)",advance='no') '#Iter Time Ca '
       do icavi = 1_ip, max_cavities
         if ( cavities(icavi) % initialized ) then
            write(lun_carcy_res_sld,"(a)",advance='no') 'C'//trim(intost(icavi))//' C'//trim(intost(icavi))//'_cycle C'&
               //trim(intost(icavi))//'_phase C'//trim(intost(icavi))//'_V C'//trim(intost(icavi))//'_P C'&
               //trim(intost(icavi))//'_outflow_P ' 
         end if
      end do
      write(lun_carcy_res_sld,"(a)") ''

    case( 2_ip )
       auxr(:)=0.0_rp
       cycid(:)=0_ip

       do icavi = 1_ip, max_cavities
         if( cavities(icavi) % initialized ) then 
            do icycle=1, max_cycles
               if(  cycles(icycle) % active .and. cycles(icycle) % cavity_controled .eq. icavi) then
                  auxr(icavi)  = cycles(icycle) % assoc_wdk % pres(TIME_N)
                  cycid(icavi) = cycles(icycle) % cycle_id
               endif
            enddo
         end if
       end do

       write(lun_carcy_res_sld, 101, advance='no')     &
            ittim,                                     &  ! Step increment
            cutim,                                     &   ! Step time
            volcai(3)                                 ! Calcium integral overall

       do icavi = 1_ip, max_cavities
         if ( cavities(icavi) % initialized ) then
            write(lun_carcy_res_sld, 102, advance='no')     &
               &   icavi,                                   &  ! Cavity #1
               &   cycid(icavi),                            &  ! Cavity #1: cycle
               &   cavities(icavi) % phase,                 &  ! Cavity #1: Phase.
               &   cavities(icavi) % volume(TIME_N),        &  ! Cavity #1: volume
               &   cavities(icavi) % pres(TIME_N),          &  ! Cavity #1: cavity pressure
               &   auxr(icavi)                                 ! Cavity #1: outflow pressure
         end if
       end do
       write(lun_carcy_res_sld, '(a)') '' 

       flush( lun_carcy_res_sld )

    end select
    ! -------------------------------
101 format(2x, i9, 2x, e16.8e3, 2x,  e16.8e3 )
102 format( (3(2x,i2) ,3(2x,e16.8e3)) )
    ! -------------------------------
  end subroutine sld_cardiac_cycle_write_res

  ! ---------------------------------------------------------------------
  !> @author   Adria  Quintanas-Corominas (adria.quintanas@bsc.es)
  !> @author   Alfonso Santiago (alfonso.santiago@bsc.es)
  ! ---------------------------------------------------------------------
  subroutine sld_cardiac_cycle_write_cvg( &
       itask )
    use def_solidz, only            :  lun_carcy_cvg_sld
    use def_master, only            :  itinn, ittim, modul, intost
    use def_master, only            :  ITER_K
    ! -------------------------------
    implicit none
    ! -------------------------------
    integer(ip), intent(in)        :: itask
    real(rp)                       :: auxr(max_cavities)
    integer(ip)                    :: icycle, cycid(max_cavities), i, icavi
    ! -------------------------------
    select case(itask)

    case( 1_ip )

       write(lun_carcy_cvg_sld,"(100a)") '# --| ALYA Cardiac Cycle Computations'
       write(lun_carcy_cvg_sld,"(100a)") '# --| '
       write(lun_carcy_cvg_sld,"(100a)") '# --|  1. Step increment'
       write(lun_carcy_cvg_sld,"(100a)") '# --|  2. Iteration number'

       i = 3_ip
       do icavi = 1_ip, max_cavities
         if ( cavities(icavi) % initialized ) then
            write(lun_carcy_cvg_sld,"(100a)") '# --|  '//trim(intost(i))//'. Cavity #'//trim(intost(icavi))
            i=i+1_ip
            write(lun_carcy_cvg_sld,"(100a)") '# --|  '//trim(intost(i))//'. Cavity #'//trim(intost(icavi))//': cycle'
            i=i+1_ip
            write(lun_carcy_cvg_sld,"(100a)") '# --|  '//trim(intost(i))//'. Cavity #'//trim(intost(icavi))//': phase'
            i=i+1_ip
            write(lun_carcy_cvg_sld,"(100a)") '# --|  '//trim(intost(i))//'. Cavity #'//trim(intost(icavi))//': volume'
            i=i+1_ip
            write(lun_carcy_cvg_sld,"(100a)") '# --|  '//trim(intost(i))//'. Cavity #'//trim(intost(icavi))//': cavity pressure'
            i=i+1_ip
            write(lun_carcy_cvg_sld,"(100a)") '# --|  '//trim(intost(i))//'. Cavity #'//trim(intost(icavi))//': outflow windkessel pressure'
            i=i+1_ip
         end if
       end do
       write(lun_carcy_cvg_sld,"(100a)") '# --|'
            
       write(lun_carcy_cvg_sld,"(a)",advance='no') '#Iter Time '
       do icavi = 1_ip, max_cavities
         if ( cavities(icavi) % initialized ) then
            write(lun_carcy_cvg_sld,"(a)",advance='no') 'C'//trim(intost(icavi))//' C'//trim(intost(icavi))//'_cycle C'&
               //trim(intost(icavi))//'_phase C'//trim(intost(icavi))//'_V C'//trim(intost(icavi))//'_P C'&
               //trim(intost(icavi))//'_outflow_P ' 
         end if
      end do
      write(lun_carcy_cvg_sld,"(100a)") ''


    case( 2_ip )
       auxr(:)=0.0_rp
       cycid(:)=0_ip

       do icavi = 1_ip,max_cavities
         if (cavities(icavi) % initialized ) then
            do icycle=1, max_cycles
               if(  cycles(icycle) % active .and. cycles(icycle) % cavity_controled .eq. icavi) then
                  auxr(icavi)  = cycles(icycle) % assoc_wdk % pres(ITER_K)
                  cycid(icavi) = cycles(icycle) % cycle_id
               endif
            enddo
         end if
       end do

       write(lun_carcy_cvg_sld, 101, advance = 'no') &
            &   ittim,                               &  ! Step increment
            &   itinn(modul)                            ! Iteration number

       do icavi=1_ip, max_cavities
         if( cavities(icavi) % initialized ) then
            write(lun_carcy_cvg_sld, 102, advance='no')                      &
                  &   icavi,                                   &  ! Cavity #1
                  &   cycid(icavi),                            &  ! Cavity #1: cycle
                  &   cavities(icavi) % phase,                 &  ! Cavity #1: Phase.
                  &   cavities(icavi) % volume(ITER_K),        &  ! Cavity #1: volume
                  &   cavities(icavi) % pres(ITER_K),          &  ! Cavity #1: cavity pressure
                  &   auxr(icavi)                                 ! Cavity #1: outflow pressure
         end if
       end do
       write(lun_carcy_cvg_sld, '(a)') ''

       flush( lun_carcy_cvg_sld )

    end select
    ! -------------------------------
101 format(2x, i9, 2x, i9 )
102 format( (3(2x,i2) ,3(2x,e16.8e3)) )
    ! -------------------------------
  end subroutine sld_cardiac_cycle_write_cvg

  ! ---------------------------------------------------------------------
  !> @author   Adria  Quintanas-Corominas (adria.quintanas@bsc.es)
  !> @author   Alfonso Santiago (alfonso.santiago@bsc.es)
  ! ---------------------------------------------------------------------
  subroutine sld_cardiac_cycle_read_data( )
    use def_master,           only :  intost
    use def_domain,           only :  ndime
    use def_inpout,           only :  words, param, exists, getint, getrea, getcha, nnpar
    use mod_ecoute,           only :  ecoute
    use mod_messages,         only :  livinf, messages_live
    use def_inpout,           only :  getint
    ! -------------------------------
    implicit none
    ! -------------------------------
    integer(ip)                    :: icavi, ii, ipar
    integer(ip)                    :: id_aux, cycle_transition_id
    real(rp)                       :: p1(3), p2(3), v1(3), nrm
    character(5)                   :: cycle_transition_method
    ! -------------------------------
    kfl_cardiac_cycle = .true.

    call ecoute('sld_cardiac_cycle_read_data')

    !CARDIAC_CYCLE section
    do while(words(1)/='ENDCA')
    

       if( words(1) == 'CAVIT' )then
          ncavi=ncavi+1
          kfl_heart_cavity = .true.

          !
          ! Get cavity identification and check it was initialized
          !
          icavi = int(param(1), kind=ip)
          if( cavities(icavi) % initialized )then
             call runend('MOD_SLD_CARDIAC_CYCLE: Two cavities with the same ID')
          else
             cavities(icavi) % initialized   = .True.
             cavities(icavi) % bouset        = 0_ip
             cavities(icavi) % valset        = 0_ip
             cavities(icavi) % plane_nodset  = 0_ip
             cavities(icavi) % plane_origin  = 0.0_rp 
             cavities(icavi) % plane_normal  = 0.0_rp
             cavities(icavi) % plane_vector  = 0.0_rp
             cavities(icavi) % volume        = 0.0_rp
             cavities(icavi) % volume_method = VOL_METHOD_MIXEDPROD !default
          end if

          call ecoute('sld_cardiac_cycle_read_data')

          do while(words(1)/='ENDCA')

              
             !
             ! Boundary defining the cavity
             !
             if (words(1)=='BOUND') then 
                cavities(icavi) % bouset = int(param(1), kind=ip)
             else if(  words(1) == 'PLANE' )then
                !
                ! Plane definition for open cavity
                !
                if( exists('POINT') )then

                   ! Plane defined by 2 points
                   p1(1:ndime) = param(2:4)
                   p2(1:ndime) = param(5:7)
                   v1(:)= p1(:) - p2(:) 
                   nrm = sqrt(v1(1)**2 + v1(2)**2 + v1(3)**2)
                   if( nrm < 0.0_rp ) call runend('MOD_SLD_CARDIAC_CYCLE: ||n|| of the plane is <= 0')
                   cavities(icavi) % plane_vector(:) = v1(:)/nrm
                   cavities(icavi) % volume_method = VOL_METHOD_PLANE_POINT

                elseif( exists('NODSE') .or. exists('SETNO') )then
                   ! Example usage:
                   ! - Calcuate volume simply using a vector in a hole closing plane
                   ! CAVITY: 1
                   !     BOUNDARY: 3 
                   !     PLANE NODSET 438 448
                   ! END_CAVITY
                   ! - Calcuate volume below the plane defined by 3 nodes. Nodes must be in clockwise order, when looking at the cavity
                   !   from the top. The triangles above the plane will not contribute to the volume 
                   ! CAVITY: 1
                   !     BOUNDARY: 3 
                   !     PLANE NODSET 438 448 465
                   ! END_CAVITY


                   ! Plane defined by a node set
                   if ( nnpar == 2 ) then
                        cavities(icavi) % plane_nodset = (/int(param(2),kind=ip), int(param(3),kind=ip), 0_ip/)
                   elseif( nnpar==3 ) then
                        cavities(icavi) % plane_nodset = (/int(param(2),kind=ip), int(param(3),kind=ip), int(param(4),kind=ip)/)
                   else
                        call runend("PLANE NODESET FOR CAVITY "//trim(intost(icavi))//" HAS TO HAVE 2 OR 3 VALUES. "//trim(intost(nnpar))//" VALUES ARE PROVIDED")
                   end if

                   cavities(icavi) % volume_method = VOL_METHOD_PLANE_NODESET

                elseif( exists('ORIGI') )then

                   ! Plane origin
                   cavities(icavi) % plane_origin = param(2:4)
                   

                else
                   ! Plane defined by the normal vector
                   call runend('MOD_SLD_CARDIAC_CYCLE: PLANE POINT / OIRIGIN / NODSET')

                endif

             else if (words(1)=='VALFI') then 
               cavities(icavi) % hole_field = getint('VALFI',0_ip,'*Node field with hole labels')
             else if (words(1)=='VALSE') then 
               ! calculate volume by closing automatically all holes (adding valve centroid) on the 
               ! intersection of the CODBOs in BOUNDARY and VALSET, each VALSET specified a element 
               ! ring marking the hole, e.g.
               ! CAVITY: 1
               !     BOUNDARY: 3 
               !     VALSET 4 5 6
               !     VOLMETHOD MIXED|DIVERGENCE (optional)
               ! END_CAVITY
               ! Here the holes will consist of the nodes on the intersection of codbo 3&4. 3&5, 3&6
               ! By default the volume is calculated as 1/6th of the sum of mixed products of vertices 
               ! of triangulated faces (VOLMETHOD MIXED). Alternatively (VOLMETHOD DIVERGENCE) calculated volume
               ! as a flux through the closed surface (slower method, advantages are not clear)
               !
               ! Alternatively
               ! CAVITY: 1
               !     BOUNDARY: 3 
               !     VALSET 4 5 6
               !     VOLMETHOD MIXED|DIVERGENCE (optional)
               !     VALFIELD 1
               ! END_CAVITY
               ! Will use field 1 defined on NODES, where cavity with CODBO 3 will have holes defined by nodes with 
               ! labes 4,5,6 in the field. Make sure you label all the holes. You can use create_hole_field.py in Utils/user

               if ( nnpar>max_valves ) then
                  call runend("SOLIDZ_CARDIAC_CYCLE: Number of holes in specified in parameteres "//trim(intost(nnpar))//&
                        " is larger than the maximum "// trim(intost(max_valves))//&
                        "for cavity "//trim(intost(icavi)))
               end if
               do ipar = 1_ip, nnpar
                  cavities(icavi) % valset(ipar) = int(param(ipar),ip)
               enddo               
               cavities(icavi) % nvalves = nnpar        
               
                  
             else if (words(1) == 'VOLME') then
               if(words(2)=='DIVER') then
                  cavities(icavi) % volume_method = VOL_METHOD_DIVERGENCE
               else if (words(2)=='MIXED') then
                  cavities(icavi) % volume_method = VOL_METHOD_MIXEDPROD
               else
                  call runend("Unknown volume calculation method "//trim(words(2))//" for cavity "//trim(intost(icavi)))
               end if

             end if
             
             call ecoute('sld_cardiac_cycle_read_data')
          end do

       else if( words(1) == 'VARIA') then
          kfl_pressure_cycle = .True.
          ncycles=ncycles+1
          do while( words(1) /= 'ENDVA' )
             call ecoute('sld_cardiac_cycle_read_data')
             if( words(1) == 'CYCLE') then
                id_aux=getint('CYCLE',0_ip,'#CYCLE_ID')
                cycles(id_aux) % cycle_id = id_aux
                cycles(id_aux) % defined  = .True.
                do while( words(1) /= 'ENDCY' )
                   if ( words(1) == 'FORCA') cycles(id_aux) %  cavity_controled = int(param(1),ip)
                   if ( words(1) == 'MAXBE') cycles(id_aux) %  max_beats = int(param(1),ip)
                   if ( words(1) == 'PREVC') cycles(id_aux) %  prevcycle = int(param(1),ip)
                   if ( words(1) == 'NEXTC') cycles(id_aux) %  nextcycle = int(param(1),ip)
                   if ( words(1) == 'TZERO') cycles(id_aux) %  tzero = param(1)
                   if ( words(1) == 'TPSTR') cycles(id_aux) %  tpstr = param(1)
                   if ( words(1) == 'GAINC') then
                      cycles(id_aux) %  gains_contraction(1) = param(1)
                      cycles(id_aux) %  gains_contraction(2) = param(2)
                   endif
                   if ( words(1) == 'PZERO') cycles(id_aux) %  pzero = param(1)
                   if ( words(1) == 'PSTR0') cycles(id_aux) %  pstr0 = param(1)
                   if ( words(1) == 'PART0') cycles(id_aux) %  assoc_wdk % p_ini = param(1)
                   if ( words(1) == 'CPRES') cycles(id_aux) %  assoc_wdk % c     = param(1)
                   if ( words(1) == 'RPRES') cycles(id_aux) %  assoc_wdk % r     = param(1)
                   if ( words(1) == 'GAINR') then
                      cycles(id_aux) %  gains_relaxation(1) = param(1)
                      cycles(id_aux) %  gains_relaxation(2) = param(2)
                   endif
                   if ( words(1) == 'PPOST') cycles(id_aux) %  ppost = param(1)

                   call ecoute('sld_cardiac_cycle_read_data')
                end do
             endif
          enddo
       else if( words(1) == 'TRANS' ) then
          ! Usage
          ! CARDIAC_CYCLE
          !    ...
          !    TRANSITION PHASE=41 METHOD=(DVOLUME | DCALCIUM)
          !    ...
          ! END_CARDIAC_CYCLE
          ! PHASE = phase id, 41 means 4 to 1. The only one programmed for now
          ! METHOD - method of phase trnasition: 
          !     DVOLUME  -- 4->1 transition happens when volume starts decreasing
          !     DCALCIUM -- 4->1 transition happens when calcium starts to increase 
          cycle_transition_id     = getint('PHASE',0_ip,'!Phase IDs')
          cycle_transition_method = getcha('METHO','12345','!Phase transition method')
          if ( cycle_transition_id .ne. 41_ip ) call runend("sld_cardiac_phases: Unknown cardiac cycle phase id: "//trim(intost(cycle_transition_id))) 
          if ( cycle_transition_method == 'DVOLU' ) then !delta volume
             phase41_trans_method = PHASE_TRANSITION_DVOL
             call messages_live("CARDIAC CYCLE: Phase 4->1 transition is based on volume decrease")
          elseif ( cycle_transition_method == 'DCALC' ) then !delta calcium
             phase41_trans_method = PHASE_TRANSITION_DCAL
             call messages_live("CARDIAC CYCLE: Phase 4->1 transition is based on calcium increase")
          else
             call runend("sld_cardiac_phases: Unknown cardiac cycle "//trim(intost(cycle_transition_id))//" transition method: "//trim(cycle_transition_method)) 
          end if
       end if

       call ecoute('sld_cardiac_cycle_read_data')
    end do
    ! ------------------------------

    if ( cavities(icavi) % valset(1)==0_ip ) then !if automatic hole closing is off
       !
       ! Check if the nodeset is valid
       !
       do icavi= 1, max_cavities
            if( cavities(icavi) % initialized )then
               do ii = 1, 2
                  if ( cavities(icavi) % plane_nodset(ii) == 0_ip ) then
                     call runend("VOLUME CALCULATION FOR CAVITY "//trim(intost(icavi))//" NODE NUMBER "//trim(intost(ii))//" IS UNDEFINED. 2 OR 3 NODES ARE NEEDED")
                  end if
               enddo
               if (cavities(icavi) % plane_nodset(3) == 0_ip) then
                  call messages_live("VOLUME CALCULATION FOR CAVITY "//trim(intost(icavi))//" USES 2 NODES")
               else
                  call messages_live("VOLUME CALCULATION FOR CAVITY "//trim(intost(icavi))//" USES 3 NODES")
               end if
            end if
       enddo
    else
       if( any(cavities(icavi) % plane_nodset(:) .ne. 0_ip ) ) then
          call runend("VOLUME CALCULATION FOR CAVITY "//trim(intost(icavi))//" USES 2 METHODS FOR VOLUME CALCULATION (valset and plane)")
       end if
    end if
       
       
  end subroutine sld_cardiac_cycle_read_data

  ! ---------------------------------------------------------------------
  !> @author   Adria Quintanas-Corominas (adria.quintanas@bsc.es)
  !> @author   Alfonso Santiago (alfonso.santiago@bsc.es)
  ! ---------------------------------------------------------------------
  subroutine sld_cardiac_cycle_send_data( &
       itask )
    ! -------------------------------
    use mod_exchange,           only : exchange_add
    implicit none
    integer(ip),  intent(in)       :: itask
    integer(ip)                    :: icavi, icycle
    ! ------------------------------- 
    select case(itask)

    case( 1_ip )

       call exchange_add(kfl_cardiac_cycle)
       call exchange_add(kfl_heart_cavity)
       call exchange_add(kfl_pressure_cycle)
       call exchange_add(ncavi)
       call exchange_add(ncycles)

    case( 2_ip )

       call exchange_add(volcai) 
       call exchange_add(phase41_trans_method)

       if( kfl_heart_cavity )then

          do icavi = 1,size(cavities)
             call exchange_add(cavities(icavi) % initialized )
             call exchange_add(cavities(icavi) % phase )
             call exchange_add(cavities(icavi) % bouset )
             call exchange_add(cavities(icavi) % plane_nodset)
             call exchange_add(cavities(icavi) % plane_vector)
             call exchange_add(cavities(icavi) % plane_origin)
             call exchange_add(cavities(icavi) % plane_normal)
             call exchange_add(cavities(icavi) % plane_points)
             call exchange_add(cavities(icavi) % valset)
             call exchange_add(cavities(icavi) % phase_counter)
             call exchange_add(cavities(icavi) % nvalves )
             call exchange_add(cavities(icavi) % volume_method)
             call exchange_add(cavities(icavi) % hole_field)
          end do

       end if

       if( kfl_pressure_cycle )then

          do icycle = 1,size(cycles)
             call exchange_add(cycles(icycle) % defined )
             call exchange_add(cycles(icycle) % active )
             call exchange_add(cycles(icycle) % cycle_id )
             call exchange_add(cycles(icycle) % cavity_controled )
             call exchange_add(cycles(icycle) % cavity_controled )
             call exchange_add(cycles(icycle) % n_beats )
             call exchange_add(cycles(icycle) % max_beats )
             call exchange_add(cycles(icycle) % prevcycle )
             call exchange_add(cycles(icycle) % nextcycle )
             call exchange_add(cycles(icycle) % tzero)
             call exchange_add(cycles(icycle) % tpstr)
             call exchange_add(cycles(icycle) % pzero)
             call exchange_add(cycles(icycle) % pstr0)
             call exchange_add(cycles(icycle) % gains_contraction)
             call exchange_add(cycles(icycle) % assoc_wdk % p_ini)
             call exchange_add(cycles(icycle) % assoc_wdk % c)
             call exchange_add(cycles(icycle) % assoc_wdk % r)
             call exchange_add(cycles(icycle) % gains_relaxation)
             call exchange_add(cycles(icycle) % ppost)
          end do

       endif

    end select
  end subroutine sld_cardiac_cycle_send_data


  ! ---------------------------------------------------------------------
  !> @author   Adria Quintanas-Corominas (adria.quintanas@bsc.es)
  !> @author   Alfonso Santiago (alfonso.santiago@bsc.es)
  ! ---------------------------------------------------------------------
  subroutine sld_cardiac_cycle_manage_restart()
    use mod_restart,  only : restart_add
    use mod_strings,  only : integer_to_string
    ! -------------------------------
    implicit none
    ! -------------------------------
    integer(ip)                    :: icavi, icycle
    ! -------------------------------

    do icavi = 1, max_cavities
       !call restart_add(cavities(icavi) % initialized ,'initialized' //integer_to_string(icavi))
       !call restart_add(cavities(icavi) % bouset      ,'bouset'      //integer_to_string(icavi))     
       !call restart_add(cavities(icavi) % plane_nodset,'plane_nodset'//integer_to_string(icavi))
       !call restart_add(cavities(icavi) % plane_origin,'plane_origin'//integer_to_string(icavi))
       !call restart_add(cavities(icavi) % plane_vector,'plane_vector'//integer_to_string(icavi))
       !call restart_add(cavities(icavi) % plane_points,'plane_points'//integer_to_string(icavi))
       call restart_add(cavities(icavi) % volume        ,'volume'      //integer_to_string(icavi))
       call restart_add(cavities(icavi) % dvol          ,'dvol'        //integer_to_string(icavi))
       call restart_add(cavities(icavi) % ddvol         ,'ddvol'       //integer_to_string(icavi))
       call restart_add(cavities(icavi) % ini_vol       ,'ini_vol'     //integer_to_string(icavi))
       call restart_add(cavities(icavi) % end_dia_vol   ,'end_dia_vol' //integer_to_string(icavi))
       call restart_add(cavities(icavi) % end_sys_vol   ,'end_sys_vol' //integer_to_string(icavi))
       call restart_add(cavities(icavi) % prest         ,'prest'       //integer_to_string(icavi))
       call restart_add(cavities(icavi) % pres          ,'pres'        //integer_to_string(icavi))
       call restart_add(cavities(icavi) % phase         ,'phase'       //integer_to_string(icavi))
       call restart_add(cavities(icavi) % phase_counter ,'phasecounter'//integer_to_string(icavi))
    end do

    do icycle = 1, max_cycles
       !call restart_add(cycles(icycle) % defined           ,'defined'          //integer_to_string(icycle))          
       !call restart_add(cycles(icycle) % active            ,'active'           //integer_to_string(icycle))         
       !call restart_add(cycles(icycle) % cycle_id          ,'cycle_id'         //integer_to_string(icycle))        
       !call restart_add(cycles(icycle) % cavity_controled  ,'cavity_controled' //integer_to_string(icycle))
       !call restart_add(cycles(icycle) % n_beats           ,'n_beats'          //integer_to_string(icycle))      
       !call restart_add(cycles(icycle) % max_beats         ,'max_beats'        //integer_to_string(icycle))       
       !call restart_add(cycles(icycle) % prevcycle         ,'prevcycle'        //integer_to_string(icycle))       
       !call restart_add(cycles(icycle) % nextcycle         ,'nextcycle'        //integer_to_string(icycle))        
       call restart_add(cycles(icycle) % last_phase_change ,'last_phase_change'//integer_to_string(icycle))
       !call restart_add(cycles(icycle) % tzero             ,'tzero'            //integer_to_string(icycle))          
       !call restart_add(cycles(icycle) % tpstr             ,'tpstr'            //integer_to_string(icycle))          
       !call restart_add(cycles(icycle) % pzero             ,'pzero'            //integer_to_string(icycle))          
       !call restart_add(cycles(icycle) % pstr0             ,'pstr0'            //integer_to_string(icycle))    
       !call restart_add(cycles(icycle) % gains_contraction ,'gains_contraction'//integer_to_string(icycle))
       !call restart_add(cycles(icycle) % gains_relaxation  ,'gains_relaxation' //integer_to_string(icycle)) 
       !call restart_add(cycles(icycle) % ppost             ,'ppost'            //integer_to_string(icycle))        
       call restart_add(cycles(icycle) % assoc_wdk % pres  ,'assoc_wdk % pres' //integer_to_string(icycle))
       !call restart_add(cycles(icycle) % assoc_wdk % p_ini ,'assoc_wdk % p_ini'//integer_to_string(icycle))
       !call restart_add(cycles(icycle) % assoc_wdk % r     ,'assoc_wdk % r'    //integer_to_string(icycle))
       !call restart_add(cycles(icycle) % assoc_wdk % c     ,'assoc_wdk % c'    //integer_to_string(icycle))
    end do
    call restart_add(volcai, 'volcai')

    ! -------------------------------
  end subroutine sld_cardiac_cycle_manage_restart

  ! ---------------------------------------------------------------------
  !> @author   Adria Quintanas-Corominas (adria.quintanas@bsc.es)
  !> @author   Alfonso Santiago (alfonso.santiago@bsc.es)
  !> @author   Mariano Vázquez (mariano.vazquez@bsc.es)
  ! ---------------------------------------------------------------------
  subroutine sld_cardiac_cycle_manage_cavity_pressure( &
       itask )
    !--------------------------------
    implicit none
    ! -------------------------------
    integer(ip),  intent(in)       :: itask
    integer(ip)                    :: icycle
    ! -------------------------------
    if( .not. kfl_pressure_cycle ) return
    do icycle=1, max_cycles
       if(cycles(icycle) % defined) then
          cycles(icycle) % active = .True.
          call sld_cardiac_phases(itask, cycles(icycle) % cycle_id)
       endif
    enddo
    ! -------------------------------
  endsubroutine sld_cardiac_cycle_manage_cavity_pressure

  ! ---------------------------------------------------------------------
  !> @author   Adria Quintanas-Corominas (adria.quintanas@bsc.es)
  !> @author   Alfonso Santiago (alfonso.santiago@bsc.es)
  !> @author   Mariano Vázquez (mariano.vazquez@bsc.es)
  ! ---------------------------------------------------------------------
  subroutine sld_cardiac_phases( &
       itask, id )
    use def_master,           only :  ittim, itinn, cutim, modul, dtime, intost
    use def_master,           only :  ITER_K, TIME_N, TIME_N_MINUS_1, ITASK_ENDINN, ITASK_ENDITE
    !--------------------------------
    implicit none
    ! -------------------------------
    integer(ip),  intent(in)       :: itask
    integer(ip),  intent(in)       :: id
    ! -------------------------------
    integer(ip)                    :: icavi
    integer(ip)                    :: dummi
    real(rp)                       :: t_init, t_pres, preload_press, prestr_press, c_wdk, r_wdk, ini_wdk, p_fill, fill_rate
    real(rp)                       :: dvol_aux
    real(rp)                       :: delta_phase_change
    real(rp)                       :: Cvol, Cvel          
    real(rp)                       :: gain_err_c, gain_derr_c, gain_err_r, gain_derr_r
    ! -------------------------------

    ! ------------------------< VARIABLE DEFINITION >-----------------------------------
    icavi         = cycles(id) %  cavity_controled
    t_init        = cycles(id) %  tzero
    t_pres        = cycles(id) %  tpstr
    preload_press = cycles(id) %  pzero
    prestr_press  = cycles(id) %  pstr0
    gain_err_c    = cycles(id) %  gains_contraction(1)
    gain_derr_c   = cycles(id) %  gains_contraction(2)
    ini_wdk       = cycles(id) %  assoc_wdk % p_ini
    c_wdk         = cycles(id) %  assoc_wdk % c
    r_wdk         = cycles(id) %  assoc_wdk % r
    gain_err_r    = cycles(id) %  gains_relaxation(1)
    gain_derr_r   = cycles(id) %  gains_relaxation(2)
    p_fill        = cycles(id) %  ppost

    ! -------< INITIALISATION. FIRST TIME ITERATION AND FIRST INNER ITERATION>----------
    ! ----------------------------------------------------------------------------------
    if ((ittim == 1) .and. (itinn(modul) == 1)) then
       cycles(id) % assoc_wdk % pres(ITER_K) = ini_wdk
       cycles(id) % assoc_wdk % pres(TIME_N) = ini_wdk
       cavities(icavi) % ini_vol  = cavities(icavi) % volume(ITER_K)
    endif



    select case(itask)

       ! -------------< WHAT HAPPENS AT THE END OF EACH INNER ITERATION>-----------------
       ! -------------------------------------------------------------------------------
    case (ITASK_ENDINN)

       ! Calculate volume increment and rate change
       cavities(icavi) % dvol(ITER_K) = cavities(icavi) % volume(ITER_K)   - cavities(icavi) % volume(TIME_N)
       cavities(icavi) % ddvol(ITER_K) = cavities(icavi) % dvol(ITER_K) / dtime

       ! -------------------< EVOLVE WINDKESSEL (aortic pressure)  MODEL >---------------------
       if(cavities(icavi) % phase .ne. 0) then
          dvol_aux=0.0_rp
          if(cavities(icavi) % phase .eq. 2) dvol_aux = cavities(icavi) % ddvol(ITER_K)
          cycles(id) % assoc_wdk % pres(ITER_K) = cycles(id) % assoc_wdk % pres(TIME_N) - (dtime / c_wdk) *  &
               (dvol_aux + (cycles(id) % assoc_wdk % pres(ITER_K) / r_wdk))
       endif

       select case(cavities(icavi) % phase)
          ! -----------------------< PHASES OF THE CYCLE  >------------------------------
          !
          ! Phases:
          !   
          !    0) Preload and prestress (phase = 0)
          !        * PARAMETERS:
          !            + pzero_sld
          !            + tzero_sld
          !            + pstr0_sld
          !            + tpstr_sld

          !    1) Isovolumetric contraction:  Keeps volume constant
          !        * PARAMETERS:
          !            + gains_contraction
          !    2) Ejection: Ventricular pressure equal to aorta
          !
          !    3) Isovolumetric relaxation: Keeps volume constant
          !        * PARAMETERS:
          !            + gains_relaxation
          !
          !    4) Filling: fills ventricle
          !        * PARAMETERS:
          !            + none
          !
          ! -------------------------------------------------------------------------------
          ! --------------------< FLAG=0. PRELOAD AND PRE-STRESS  >------------------------
       case(0_ip)
          ! ------------------< ADDING PRELOAD CONTRIBUTION >-------------------------
          cavities(icavi) % pres(ITER_K) = preload_press * (cutim / t_init)


          ! ------------------< ADDING PRESTRES CONTRIBUTION >-------------------------
          if(prestr_press .gt. 0.1_rp) then
             if (cutim < t_pres) then
                cavities(icavi) % prest(ITER_K) = prestr_press * (cutim / t_pres)
             else
                Cvol = cavities(icavi) % volume(ITER_K) / cavities(icavi) % prest(ITER_K)
                Cvel = 1.0_rp / 1.0_rp ! Acceptable ranges are from 1.0 to 1.5

                dvol_aux = cavities(icavi) % volume (ITER_K) - cavities(icavi) % ini_vol

                cavities(icavi) % prest(ITER_K) = cavities(icavi) % prest(TIME_N) + dvol_aux / Cvol + cavities(icavi) % ddvol(ITER_K) / Cvel

                if (cavities(icavi) % prest(ITER_K) < 1.0_rp) then
                   call runend('sld_windk: negative pre-stress!')
                   cavities(icavi) % prest(ITER_K) = 1.0_rp
                endif
             end if
          endif

          ! ------------------< SAVE LAST COMPUTED VALUE >-------------------------
          cavities(icavi) % end_dia_vol = cavities(icavi) % volume(ITER_K)
          ! ----------------------------------------------------------------------------

          ! -------------< FLAG=1. ISOVOLUMETRIC CONTRACTION  >------------------------
       case(1_ip)
          !gain_err = cavities(icavi) % volume(ITER_K) / cavities(icavi) % pres(ITER_K)

          dvol_aux = cavities(icavi) % volume (ITER_K) - cavities(icavi) % ini_vol
          gain_err_c = cavities(icavi) % volume (ITER_K) / cavities(icavi) % pres(ITER_K)

          cavities(icavi) % pres(ITER_K) = cavities(icavi) % pres(TIME_N) - gain_err_c * dvol_aux - gain_derr_c * cavities(icavi) % ddvol(ITER_K)
          if(cavities(icavi) % pres(ITER_K) .lt. cavities(icavi) % pres(TIME_N))then
             cavities(icavi) % pres(ITER_K) = cavities(icavi) % pres(TIME_N)
          endif

          ! Save new value in end diastolic in case oscilations or inertia make it larger
          if(cavities(icavi) % volume(TIME_N) .gt. cavities(icavi) % end_dia_vol)then
             cavities(icavi) % end_dia_vol = cavities(icavi) % volume(TIME_N)
          endif

          ! ----------------------------------------------------------------------------

          ! -------------< FLAG=2. EJECTION  >------------------------------------------
       case(2_ip)
          ! -------------------< VENTRICLE PRESSURE EQUAL TO WDK >---------------------
          cavities(icavi) % pres(ITER_K) = cycles(id) % assoc_wdk % pres(ITER_K)
          cavities(icavi) % end_sys_vol = cavities(icavi) % volume(ITER_K)     
          ! ----------------------------------------------------------------------------

          ! -------------< FLAG=3. ISOVOLUMETRIC RELAXATION  >---------------------------
       case(3_ip)
          ! Fixed pressure is overrided and replaced with this automatic value
          gain_err_r =  cavities(icavi) % pres(ITER_K)/cavities(icavi) % volume(ITER_K)

          dvol_aux = cavities(icavi) % volume(ITER_K) - cavities(icavi) % end_sys_vol

          cavities(icavi) % pres(ITER_K) = cavities(icavi) % pres(TIME_N) - gain_err_r*dvol_aux  - gain_derr_r * cavities(icavi) % ddvol(ITER_K)

          ! clipping works well with automatic gain
          if(cavities(icavi) % pres(ITER_K) .gt. cavities(icavi) % pres(TIME_N))then
             cavities(icavi) % pres(ITER_K) = cavities(icavi) % pres(TIME_N)
          endif

          ! ----------------------------------------------------------------------------

          ! -------------< FLAG=4. DIASTOLIC FILLING > ---------------------------------
       case(4_ip)

          fill_rate=(preload_press-p_fill)/(cavities(icavi) % ini_vol - cavities(icavi) % end_sys_vol)
          if(fill_rate.lt.-500.0_rp) fill_rate=-500.0_rp

          cavities(icavi) % pres(ITER_K) = cavities(icavi) % pres(TIME_N) + fill_rate* cavities(icavi) % dvol(ITER_K)

          if(cavities(icavi) % pres(ITER_K).lt.preload_press)then
             cavities(icavi) % pres(ITER_K)=preload_press
          endif

       end select


       ! -------------< WHAT HAPPENS AT THE END OF EACH TIME ITERATION>-----------------
       ! -------------------------------------------------------------------------------
    case (ITASK_ENDITE)

       ! -------------< VARIABLE UPDATE >-------------------------------------------   
       cavities(icavi) % dvol(TIME_N)  = cavities(icavi) % volume(TIME_N) - cavities(icavi) % volume(TIME_N_MINUS_1)
       cavities(icavi) % ddvol(TIME_N) = cavities(icavi) % dvol(TIME_N) / dtime
       cavities(icavi) % pres(TIME_N) = cavities(icavi) % pres(ITER_K)
       cavities(icavi) % prest(TIME_N) = cavities(icavi) % prest(ITER_K)
       cycles(id) % assoc_wdk % pres(TIME_N) = cycles(id) % assoc_wdk % pres(ITER_K)
       delta_phase_change= cutim-cycles(id) % last_phase_change

       if ( phase41_trans_method == PHASE_TRANSITION_DCAL ) then
         if ( ittim > int(volcai(1), kind=ip) ) then !calcualte only one time, not every time this is called for each cavity
            volcai(2) = volcai(3) ! prev = current old
            call eccou_do_volume_integral( dummi, "doal", volcai(3) )
            volcai(1) = real( ittim, kind=ip )
         end if
         !IF(INOTSLAVE) then
         !  WRITE(1997, "(5e12.5)") cutim, volcai(1), volcai(2), volcai(3), volcai(3) - volcai(2)
         !ENDIF
       end if

       ! ----------------------------------------------------------------------------


       ! -------------< PHASE TRANSITION >-------------------------------------------
       ! Phase description:
       !   -> 0:  initialisation
       !   -> 1:  isovol contraction
       !   -> 2:  ejection
       !   -> 3:  isovol relax
       !   -> 4:  filling
       !   
       ! -------------< preload(0)->isolov contraction(1) >--------------------------
       if (cutim .ge. t_init) then
         call sld_cardiac_phase_transition(icavi, id, delta_phase_change)
       end if

    end select
    ! -------------------------------
  end subroutine sld_cardiac_phases

  
  subroutine sld_cardiac_phase_transition(icavi, id, delta_phase_change)
      use def_master,           only :  cutim, intost
      use def_master,           only :  TIME_N


      implicit none
      integer(ip), intent(in)        :: icavi
      integer(ip), intent(in)        :: id
      real(rp),    intent(in)        :: delta_phase_change
      logical(lg)                    :: make_41_transition
      real(rp)                       :: p_fill

      p_fill = cycles(id) %  ppost


      ! -------------< PHASE TRANSITION >-------------------------------------------
      ! Phase description:
      !   -> 0:  initialisation
      !   -> 1:  isovol contraction
      !   -> 2:  ejection
      !   -> 3:  isovol relax
      !   -> 4:  filling
      !   
      ! -------------< preload(0)->isolov contraction(1) >--------------------------
      if ( cavities(icavi) % phase == 0_ip ) then
         cavities(icavi) % phase = 1_ip
         cycles(id) % last_phase_change=cutim
         cavities(icavi) % phase_counter = 0_ip
         ! --------------------------------------------------------------------------
         ! -------------< isovol contract(1)->ejection(2) >--------------------------
         ! CONDITIONS:
         !   * in phase 1
         !   * wdk pressure greater than windkessel pressure 
         !
      elseif ((cavities(icavi) % phase == 1_ip) .and. &
            (cavities(icavi) % pres(TIME_N) > cycles(id) % assoc_wdk % pres(TIME_N))) then

         cavities(icavi) % phase = 2_ip
         cycles(id) % last_phase_change=cutim
         cavities(icavi) % phase_counter = 0_ip

         ! --------------------------------------------------------------------------
         ! -------------< ejection(2)->isovol relax(2) >-----------------------------
         ! CONDITIONS:
         !   * in phase 2
         !   * Volume is increasing (geometry is relaxing)
         !   * It's still contracted (required for oscilations)
         !
      elseif (cavities(icavi) % phase == 2_ip .and. delta_phase_change.gt.0.05_rp) then 

        if ( (cavities(icavi) % dvol(TIME_N) > 0.00001_rp)  .and. &
             (cavities(icavi) % volume(TIME_N) < 0.99_rp * cavities(icavi) % end_dia_vol)) then
            cavities(icavi) % phase_counter = cavities(icavi) % phase_counter + 1_ip
        else 
            cavities(icavi) % phase_counter = 0_ip
        end if
     

        if (cavities(icavi) % phase_counter > 5_ip) then
            cavities(icavi) % phase = 3_ip
            cycles(id) % last_phase_change=cutim
            cavities(icavi) % phase_counter = 0_ip
        end if
        
         ! --------------------------------------------------------------------------
         ! -------------< isovol relax(3)->filling(4) >------------------------------
         ! CONDITIONS:
         !   * in phase 3
         !   * Pressure is under threashold
         !
      elseif ((cavities(icavi) % phase == 3_ip) .and. &
            (cavities(icavi) % pres(TIME_N) < p_fill))then

         cavities(icavi) % phase = 4
         cavities(icavi) % phase_counter = 0
         cycles(id) % last_phase_change=cutim
         ! --------------------------------------------------------------------------
         ! -------------< filling(4)->isovol contract(1) >---------------------------
         ! CONDITIONS:
         !   * in phase 4
         !   * Geometry is contracting
         !
      elseif ((cavities(icavi) % phase == 4_ip) .and. &
            delta_phase_change.gt.0.05_rp ) then 
   
         make_41_transition = .FALSE.
         if ( phase41_trans_method == PHASE_TRANSITION_DCAL ) then
            make_41_transition = (volcai(3) - volcai(2)) > 10*epsilon(1.0_rp)
         elseif( phase41_trans_method == PHASE_TRANSITION_DVOL ) then
            if (cavities(icavi) % dvol(TIME_N) < -0.000001_rp) then
               cavities(icavi) % phase_counter = cavities(icavi) % phase_counter + 1_ip
            else 
               cavities(icavi) % phase_counter = 0_ip
            end if

            make_41_transition = ( cavities(icavi) % phase_counter >5_ip )
         else
            call runend("sld_cardiac_phases: Unknown 4->1 cardiac cycle transition method: "//trim(intost(phase41_trans_method))) 
         end if

         if ( make_41_transition ) then
            cavities(icavi) % phase = 1_ip
            cycles(id) % last_phase_change=cutim
            cycles(id) % n_beats = cycles(id) % n_beats + 1_ip
            cavities(icavi) % phase_counter = 0_ip
         end if
      end if

  end subroutine sld_cardiac_phase_transition




  subroutine sld_cardiac_cycle_compute_cavity_volume( itask )
      use def_master, only : ITASK_ENDITE, TIME_N_MINUS_2, TIME_N_MINUS_1,TIME_N,ITER_K, ittim, intost, INOTMASTER
      implicit none
      integer(ip),  intent(in)       :: itask
      integer(ip)                    :: icavi, ivalve

      call calculate_valve_centroids()

      do icavi= 1, max_cavities
         if( cavities(icavi) % initialized )then
            if( cavities(icavi) % valset(1) == 0_ip ) then
               call compute_volume_plane(icavi)
            else               
               call compute_volume_holes(icavi)
            end if
         end if

         if ( cavities(icavi) % volume(ITER_K) < -epsilon(0.0_rp) ) then
            call runend( "VOLUME FOR CAVITY "//trim(intost(icavi))//" IS NEGATIVE" )
         end if
   

         !
         ! From solidz to master arrays
         ! 
         !update volume, only at the end of the iterations loop
         if (itask == ITASK_ENDITE) then
            cavities(icavi) % volume(TIME_N_MINUS_2) = cavities(icavi) % volume(TIME_N_MINUS_1)  !V^n-2
            cavities(icavi) % volume(TIME_N_MINUS_1) = cavities(icavi) % volume(TIME_N)          !V^n-1
            cavities(icavi) % volume(TIME_N)         = cavities(icavi) % volume(ITER_K)          !V^n
  
            if (ittim == 1) then
               cavities(icavi) % volume(TIME_N_MINUS_2) = cavities(icavi) % volume (TIME_N)
               cavities(icavi) % volume(TIME_N_MINUS_1) = cavities(icavi) % volume (TIME_N)      !initially V^n-1 = V^n = V^i 
            end if
         end if
  
      end do

      if(INOTMASTER) then
         do icavi= 1, max_cavities
            if( cavities(icavi) % initialized )then
               if( cavities(icavi) % valset(1) .ne. 0_ip ) then
                  do ivalve = 1_ip, cavities(icavi) % nvalves
                     cavities(icavi) % valves(ivalve) % centroid(:) = 0.0_rp
                  end do
               end if
            end if
         end do
      end if
      ! -------------------------------    
      !stop 1
  end subroutine


  ! ---------------------------------------------------------------------
  !> @author   Adria Quintanas-Corominas (adria.quintanas@bsc.es)
  !> @author   Alfonso Santiago (alfonso.santiago@bsc.es)
  !> @details  Calculation of the volume of a cavity cutted by the basal plane:
  !>           V = Int_{dOmega_0} (x \cdot e_1) (e_1 \cdot n) d dOmega_0    
  !>           where, x : current coordinates
  !>                  e : vector in the basal plane
  !>                  n : vector normal to the cavity surface
  !>           Refs: Levrero-Florencio et al (2020) DOI: 10.1016/j.cma.2019.112762
  !>                 Quarteroni et al. (2017)       DOI: 10.1016/j.cma.2016.05.031 
  ! ---------------------------------------------------------------------
   subroutine compute_volume_plane( icavi )
      use def_master,           only :  displ, intost
      use def_master,           only :  INOTMASTER, ITER_K, TIME_N, TIME_N_MINUS_1, TIME_N_MINUS_2, ITASK_ENDITE
      use def_domain,           only :  mnode, mnodb, ltypb, nnode, ngaus, lelbo, ltype, nboun, npoin_own
      use def_domain,           only :  elmar, ndimb, kfl_codbo, lnodb, coord, lnods
      use def_domain,           only :  htable_lninv_loc
      use mod_htable,           only :  htalid
      use mod_communications,   only :  PAR_SUM  
      use mod_maths_basic,      only :  maths_cross_product, maths_normalize_vector
      use mod_messages,         only :  messages_live

      !-----------------------------
      implicit none
      !-----------------------------
      ! -------------------------------
      integer(ip),  intent(in)       :: icavi
      ! -------------------------------
      integer(ip)                    :: ielem, ipoin, idime, iboun, igaub, inodb, ii
      integer(ip)                    :: pelty, pblty, pnodb, pgaub, pnode, inode
      real(rp)                       :: gpel(3), v1(3), v2(3)
      real(rp)                       :: baloc(3,3), bocod(3,mnodb), elcod(3,mnode), eucta
      logical(lg)                    :: add_vol_contribution
      ! -------------------------------
      !
      ! Initialise  
      !
      cavities(icavi) % volume(ITER_K) = 0.0_rp
      cavities(icavi) % plane_points(:,:) = 0.0_rp
    


      !
      ! Locate the nodes of the plane in the domain
      !
      if( INOTMASTER ) then
         do ii = 1, size( cavities(icavi) % plane_nodset, 1_ip, kind=ip )
            if ( cavities(icavi) % plane_nodset(ii) > 0_ip ) then
               ipoin = htalid(htable_lninv_loc,cavities(icavi) % plane_nodset(ii))
               if( ipoin > 0_ip .and. ipoin <= npoin_own ) then
                  cavities(icavi) % plane_points(:,ii) = coord(:,ipoin) + displ(:,ipoin,1)
               end if
            end if
         enddo
      end if

      ! mod_parall global_to_local
      do ii = 1, size( cavities(icavi) % plane_points, 2_ip, kind=ip )
         call PAR_SUM( 3_ip, cavities(icavi) % plane_points(:,ii) )
      enddo

      !
      ! Compute the new origin and normal vector of the plane
      !
      if( all(cavities(icavi) % plane_nodset(1:2) > 0_ip)  )then
         v1(:) = cavities(icavi) % plane_points(:,1) - cavities(icavi) % plane_points(:,2) 
         call maths_normalize_vector(3_ip, v1)
         cavities(icavi) % plane_vector(:) = v1(:)

         if ( cavities(icavi) % plane_nodset(3) > 0_ip ) then ! 3 nodes are provided, caculate plane normal
            ! Cavity view from top at the hole. Points clockwise on the boundary plane
            ! 
            !         v1
            !  p1 <--------- p2
            !   \           /
            !    \         / v2
            !     \       /
            !      \    |/_
            !         p3
            !
            ! v1 = p1-p2-
            ! v2 = p3-p2
            ! normal = v2 x v1, pointing inward.
            v2(:) =  cavities(icavi) % plane_points(:,3) - cavities(icavi) % plane_points(:,2)
            call maths_normalize_vector(3_ip, v2)
            cavities(icavi) % plane_normal(:) = maths_cross_product( v2(:), v1(:), 3_ip )
            cavities(icavi) % plane_origin(:) = cavities(icavi) % plane_points(:,1)
         end if
      endif

      !
      ! Compute the volume of the cavity
      !
      if( INOTMASTER ) then
         do iboun = 1, nboun
            if( kfl_codbo(iboun) == cavities(icavi) % bouset )then

               pblty=ltypb(iboun) 
               pnodb=nnode(pblty)
               pgaub=ngaus(pblty)
               ielem=lelbo(iboun)
               pelty=ltype(ielem)
               pnode=nnode(pelty)

               ! Gather bondary and element coordiantes
               do inodb=1,pnodb
                  ipoin=lnodb(inodb,iboun)
                  do idime=1, 3
                     bocod(idime,inodb) = coord(idime,ipoin) + displ(idime,ipoin,1)
                  end do
               end do

               do inode=1,pnode
                  ipoin=lnods(inode,ielem)
                  do idime=1, 3
                     elcod(idime,inode) = coord(idime,ipoin) + displ(idime,ipoin,1)
                  end do
               end do

               ! Divergence theorem
               do igaub = 1, pgaub

                  call bouder(pnodb,3_ip,ndimb,elmar(pblty)%deriv(:,:,igaub),bocod,baloc,eucta)    
                  call chenor(pnode,baloc,bocod,elcod)

                  gpel = 0.0_rp
                  do inodb = 1_ip,pnodb
                     do idime= 1_ip, 3_ip
                        gpel(idime) = gpel(idime) + elmar(pblty)%shape(inodb,igaub) * &
                           bocod(idime,inodb)
                     end do
                  end do

                  add_vol_contribution = .TRUE.

                  ! if the usr specified a 3 point plane, exclude the faces in the exterior of the cavity
                  if ( cavities(icavi) % plane_nodset(3) > 0_ip ) then ! 3 nodes are provided, check on which side of the plane is the point
                  if ( dot_product( gpel(:) - cavities(icavi) % plane_origin(:), cavities(icavi) % plane_normal(:) ) < 0 ) then
                     add_vol_contribution = .FALSE.
                  end if
                  end if

                  ! includes shifting the gpel to the coordinate system centered at one of the vector nodes
                  if (add_vol_contribution) then
                     cavities(icavi) % volume(ITER_K) =  cavities(icavi) % volume(ITER_K) &
                           &    - dot_product(gpel(:) - cavities(icavi) % plane_points(:,1), cavities(icavi) % plane_vector(:)) &
                           &    * dot_product(cavities(icavi) % plane_vector(:), baloc(:,3)) &
                           &    * elmar(pblty) % weigp(igaub) * eucta
                  end if

               enddo
            end if
         end do !iboun
      endif

      call PAR_SUM( cavities(icavi) % volume(ITER_K) )
   end subroutine 


   subroutine cardiac_cycle_redist()
      implicit none
      integer(ip)  :: icavi, ivalve

      do icavi= 1_ip, max_cavities
         IF( cavities(icavi) % initialized ) THEN
            if( cavities(icavi) % valset(1) > 0_ip ) then ! the volume is calculated by automatic hole closing
               !clean up memory if allocated (in case of repartition)
               if( associated( cavities(icavi) % valves ) ) then
                  do ivalve = 1_ip, size(cavities(icavi) % valves, kind=ip)
                     IF ( associated(cavities(icavi) % valves(ivalve) % edges) ) THEN 
                        deallocate( cavities(icavi) % valves(ivalve) % edges )
                        nullify( cavities(icavi) % valves(ivalve) % edges )
                     END IF
                  end do
                  
                  deallocate( cavities(icavi) % valves )
                  nullify( cavities(icavi) % valves )
               end if
         

               !Alocate and recalculate nodes and connectivity
               call init_cavity_hole_info(icavi)
            end if
         end if
      end do
   end subroutine

   ! Allocates arrays, identifies the hole nodes and connectivity
   ! To be run at the start and after repartitioning
   subroutine init_cavity_hole_info(icavi)
      use def_domain, only : kfl_codno, kfl_codbo, lnodb, nnode, ltypb, nboun
      use def_master, only : INOTMASTER, ITER_K
      implicit none
      integer(ip), intent(in) :: icavi
      integer(ip)             :: valve_num_nodes( cavities(icavi) % nvalves )
      integer(ip)             :: iboun, inodb, ivalve, pblty, pnodb


      cavities(icavi) % volume(ITER_K) = 0.0_rp

      if (INOTMASTER) then
         if( .not. associated( cavities(icavi) % valves ) ) THEN
            if( cavities(icavi) % nvalves>0_ip ) then
               allocate( cavities(icavi) % valves(cavities(icavi) % nvalves) )
               do ivalve = 1_ip, cavities(icavi) % nvalves
                  cavities(icavi) % valves(ivalve) % centroid(:) = 0.0_rp
                  cavities(icavi) % valves(ivalve) % number_of_nodes = 0_ip
               end do
            end if
         end if


         !Count number of nodes on each hole
         call identify_connectivity(.true.)
         do ivalve = 1_ip, cavities(icavi) % nvalves
            cavities(icavi) % valves(ivalve) % number_of_nodes = valve_num_nodes(ivalve)
         end do

         cavities(icavi) % volume(ITER_K) = 0.0_rp
         do ivalve = 1_ip, cavities(icavi) % nvalves
            IF (.NOT. associated(cavities(icavi) % valves(ivalve) % edges)) THEN 
               if ( cavities(icavi) % valves(ivalve) % number_of_nodes > 0_ip ) then
                  allocate( cavities(icavi) % valves(ivalve) % edges( 2_ip, &
                           cavities(icavi) % valves(ivalve) % number_of_nodes ) )
                  cavities(icavi) % valves(ivalve) % edges(:,:) = 0_ip
               end if
            END IF
         end do
         


         !---------------------------------------------------------
         ! Identify connectivity along the hole
         call identify_connectivity(.false.)

         ! end Identify connectivity along the hole
         !---------------------------------------------------------   
      end if

   contains
      subroutine identify_connectivity(count_only)
         implicit none
         logical(lg), intent(in) :: count_only
         integer(ip)             :: ipoin_1, ipoin_2

         !On a closed 1-d boundary, number of nodes == number of edges(line segments)
         valve_num_nodes(:) = 0_ip
         DO iboun = 1_ip, nboun !loop over boundary elements
            !I have to iterate over bouset rather than particular valset because of the field
            !When the node field is given, valset is not on the codbo anymore
            IF( kfl_codbo(iboun) == cavities(icavi) % bouset ) THEN
               DO ivalve = 1_ip, cavities(icavi) % nvalves 
                  if ( count_only .or. associated(cavities(icavi) % valves(ivalve) % edges) ) then

                     pblty=ltypb(iboun) 
                     pnodb=nnode(pblty)

                     !Iterate over all oriented edges and identify the ones belonging to the hole
                     do inodb=1_ip, pnodb
                        ipoin_1 = lnodb( inodb, iboun )
                        if( inodb<pnodb ) then !i think this if is less costly than modulo()
                           ipoin_2 = lnodb( inodb + 1_ip, iboun )
                        else
                           ipoin_2 = lnodb( 1_ip, iboun )
                        end if

                        if( node_on_hole(icavi, ipoin_1, cavities(icavi) % valset(ivalve), cavities(icavi) % bouset) .and. &
                            node_on_hole(icavi, ipoin_2, cavities(icavi) % valset(ivalve), cavities(icavi) % bouset) ) then

                           valve_num_nodes(ivalve) = valve_num_nodes(ivalve) + 1_ip 
                           if (.not. count_only) then
                              !if we ietarte over bouset, nodes need flipping
                              cavities(icavi) % valves(ivalve) % edges(:, valve_num_nodes(ivalve) ) = (/ipoin_2, ipoin_1/)
                              !if we ietarate over valset, no need to flip nodes
                              !cavities(icavi) % valves(ivalve) % edges(:, valve_num_nodes(ivalve) ) = (/ipoin_1, ipoin_2/)
                           end if
                        end if
                     end do
            
                  end if          
               END DO
            end if 
         END DO

      end subroutine

      logical(lg) pure function node_on_hole( icavi, ipoin, icode_valve, icode_cavity )
         !icode_valve -- has to be the code of the valve, the code you pass in valset
         !icode_cavity -- has to be the code you pass in bouset, the code assigned to the cavity
         use def_domain, only : xfiel
         implicit none
         integer(ip), intent(in) :: icavi, ipoin, icode_valve, icode_cavity
         integer(ip)             :: ifiel

         ifiel = cavities(icavi) % hole_field
         if( ifiel == 0_ip ) then
            !if no field is provided, look at the intersection of codnos
            node_on_hole = any( icode_valve == kfl_codno(:,ipoin) ) .and. any( icode_cavity == kfl_codno(:,ipoin) )
         else
            !if the field is provided, check that the valset code corresponds to the node code
            node_on_hole = ( icode_valve == int(xfiel(ifiel) % a(1,ipoin,1), kind=ip) )
         end if
      end function

   end subroutine
 

   subroutine calculate_valve_centroids()
      use def_master,           only : INOTMASTER, displ
      use def_domain,           only : coord, ndime
      use mod_communications,   only : PAR_SUM

      implicit none
      integer(ip) :: icavi, ivalve, IDX, ipoin, icoord
      real(rp)    :: centroid(ndime)
      integer(ip) :: nnode_valve


      if (INOTMASTER) then
         !-----------------------------------------
         ! calculate centroid contribution of each valve
         !
         do icavi =1_ip, max_cavities
            if( cavities(icavi) % initialized .and. cavities(icavi) % valset(1) .ne. 0_ip ) then
               DO ivalve = 1_ip, cavities(icavi) % nvalves
                  if(associated(cavities(icavi) % valves(ivalve) % edges)) then
                     DO IDX  = 1_ip, size( cavities(icavi) % valves(ivalve) % edges, 2_ip, kind=ip ) 
                        IF (cavities(icavi) % valves(ivalve) % edges (1,IDX)>0_ip)  THEN
                           !points are ordered around the hole in the connectivity array, 
                           !so in theory we can just grab the first from connectivity
                           !or alternatively all points, but then divide the result by 2
                           !do ii=1_ip,2_ip
                           !   ipoin = cavities(icavi) % valves(ivalve) % edges(ii,IDX) 
                           !   do icoord = 1,3_ip
                           !      cavities(icavi) % valves(ivalve) % centroid( icoord ) = &
                           !         cavities(icavi) % valves(ivalve) % centroid( icoord ) + &
                           !         coord(icoord,ipoin) + displ(icoord,ipoin,1)                               
                           !   end do
                           !end do
                           ipoin = cavities(icavi) % valves(ivalve) % edges(1,IDX) 
                           do icoord = 1,3_ip
                              cavities(icavi) % valves(ivalve) % centroid( icoord ) = &
                                 cavities(icavi) % valves(ivalve) % centroid( icoord ) + &
                                 coord(icoord,ipoin) + displ(icoord,ipoin,1)                               
                           end do
                        end if
                     end do
                  end if
               end do
            end if
         end do
      END IF !if INOTMASTER

      !-----------------------------------------
      ! calculate centroid of each valve from the contributions of the partitions
      !
      do icavi =1_ip, max_cavities
         if( cavities(icavi) % initialized ) then
            DO ivalve = 1_ip, cavities(icavi) % nvalves
               if (INOTMASTER) then 
                  nnode_valve = cavities(icavi) % valves(ivalve) % number_of_nodes
                  centroid(:) = cavities(icavi) % valves(ivalve) % centroid
               else 
                  nnode_valve = 0_ip
                  centroid(:) = 0.0_rp
               end if

               call PAR_SUM( nnode_valve )         
               call PAR_SUM( 3_ip, centroid )

               if( INOTMASTER ) then
                  cavities(icavi) % valves(ivalve) % centroid(:) = centroid / real( nnode_valve, kind = rp ) 
                  !cavities(icavi) % valves(ivalve) % centroid(:) = centroid / real( 2_ip * nnode_valve, kind = rp ) 
               end if
            END DO
         end if
      end do
      !
      ! end calculate centroid of each valve from the contributions of the partitions
      !-----------------------------------------      
   end subroutine


   subroutine compute_volume_holes( icavi )
      use def_master,           only :  displ, intost
      use def_master,           only :  INOTMASTER, ITER_K, TIME_N, TIME_N_MINUS_1, TIME_N_MINUS_2, ITASK_ENDITE
      use def_domain,           only :  mnode, mnodb, ltypb, nnode, ngaus, lelbo, ltype, nboun
      use def_domain,           only :  elmar, ndimb, kfl_codbo, lnodb, coord, lnods
      use def_domain,           only :  ndime
      use mod_htable,           only :  htalid
      use mod_messages,         only :  messages_live
      use mod_maths_basic,      only :  maths_cross_product, maths_normalize_vector
      use mod_communications,   only :  PAR_SUM,PAR_MAX
      
      !-----------------------------
      implicit none
      !-----------------------------
      ! -------------------------------
      integer(ip),  intent(in)                   :: icavi
      ! -------------------------------
      integer(ip)                                :: ielem, ipoin, idime, iboun, igaub, inodb, icoord
      integer(ip)                                :: pelty, pblty, pnodb, pgaub, pnode, inode
      real(rp)                                   :: gpel(3),r_1(3), r_2(3)
      real(rp)                                   :: baloc(3,3), bocod(3,mnodb), elcod(3,mnode), eucta
      integer(ip)                                :: idx,ipoint_1, ipoint_2, ivalve
      real(rp)                                   :: area
      real(rp)                                   :: norm
      real(rp)                                   :: centroid(ndime)

      ! -------------------------------
      !
      ! Initialise  
      !       
      !------------------------------------------------

      cavities(icavi) % volume(ITER_K) = 0.0_rp
    
      !TODO: Add check fior watertightness. The idea is that each node of all elements with a certain CODBO
      ! has to be visited even number of times. So just make an array or bools and toggle them every time a 
      ! triangle is created, at the end synchronize the array across partitions and make sure everything is 0
      ! I might need to talk to Guillaume to see how to do it best. Don't feel like keeping a big array of all 
      ! boundary nodes

      IF( INOTMASTER ) THEN
         !Calculate volume contribution by the triangles on the boundary

         DO iboun = 1_ip, nboun
            IF ( kfl_codbo(iboun) == cavities(icavi) % bouset  ) THEN ! 

               pblty=ltypb(iboun) 
               pnodb=nnode(pblty)
               pgaub=ngaus(pblty)
               ielem=lelbo(iboun)
               pelty=ltype(ielem)
               pnode=nnode(pelty)

               ! Gather bondary and element coordiantes
               do inodb=1_ip,pnodb
                  ipoin=lnodb(inodb,iboun)
                  do idime=1_ip, ndime
                     bocod(idime,inodb) = coord(idime,ipoin) + displ(idime,ipoin,1)
                  end do   
               end do

               select case ( cavities(icavi) % volume_method )
               case (VOL_METHOD_MIXEDPROD) 
                  !------------------------------------------
                  ! Volume calculation method based on mixed product
                  ! When all triangles have the same orientation, i-th triangle having vertices ai,bi,ci
                  ! The volume is 1/6 sum_i ( ai . (bi x ci) )
                  !
                  ! triangulate the element
                  do inodb=2_ip,pnodb-1                  
                     cavities(icavi) % volume(ITER_K) = cavities(icavi) % volume(ITER_K) - &
                        mixed_product( bocod(:,1), bocod(:,inodb), bocod(:,inodb+1_ip) )
                     !write(1977,*) bocod(:,1), bocod(:,inodb), bocod(:,inodb+1_ip)
                  end do
                  ! Volume calculation method based on mixed product
                  !------------------------------------------
               case (VOL_METHOD_DIVERGENCE) 
                  !------------------------------------------
                  ! Volume calculation method based on divergence theorem
                  do inode=1_ip,pnode
                     ipoin=lnods(inode,ielem)
                     do idime=1_ip, ndime
                        elcod(idime,inode) = coord(idime,ipoin) + displ(idime,ipoin,1)
                     end do
                  end do


                  ! Divergence theorem
                  DO igaub = 1_ip, pgaub

                     call bouder(pnodb,3_ip,ndimb,elmar(pblty)%deriv(:,:,igaub),bocod,baloc,eucta)    
                     call chenor(pnode,baloc,bocod,elcod)
                     gpel = 0.0_rp
                     do inodb = 1_ip,pnodb
                        do idime= 1_ip, ndime
                           gpel(idime) = gpel(idime) + elmar(pblty)%shape(inodb,igaub) * &
                              bocod(idime,inodb)
                        end do
                     end do

                     cavities(icavi) % volume(ITER_K) = cavities(icavi) % volume(ITER_K) &
                        - dot_product(gpel(:), baloc(:,3)) * elmar(pblty) % weigp(igaub) * eucta

                  END DO !DO igaub = 1, pgaub
                  ! End Volume calculation method based on divergence theorem
                  !------------------------------------------
               case default
                  call runend("Unknown volume calculation method "//trim(intost(cavities(icavi) % volume_method))//" for cavity "//trim(intost(icavi)))
               end select
            END IF
         END DO

         !Calculate volume contribution from the hole covers
         DO ivalve = 1_ip, cavities(icavi) % nvalves 
            if (associated(cavities(icavi) % valves(ivalve) % edges)) then
               DO IDX  = 1_ip, size( cavities(icavi) % valves(ivalve) % edges, 2_ip, kind=ip ) 
                  IF (cavities(icavi) % valves(ivalve) % edges(1,IDX)>0_ip)  THEN

                     ipoint_1 = cavities(icavi) % valves(ivalve) % edges(1,IDX)
                     ipoint_2 = cavities(icavi) % valves(ivalve) % edges(2,IDX)       

                     select case (cavities(icavi) % volume_method)
                     case ( VOL_METHOD_MIXEDPROD ) 
                        !------------------------------------------
                        ! Volume calculation based on mixed product
                        cavities(icavi) % volume(ITER_K) = cavities(icavi) % volume(ITER_K) - &
                           mixed_product( &
                              coord(:, ipoint_1) + displ(:, ipoint_1, 1), &
                              coord(:, ipoint_2) + displ(:, ipoint_2, 1), &
                              cavities(icavi) % valves(ivalve) % centroid(:) )
                        !write(1977,*) coord(:, ipoint_1) + displ(:, ipoint_1, 1), &
                        !   coord(:, ipoint_2) + displ(:, ipoint_2, 1), &
                        !   cavities(icavi) % valves(ivalve) % centroid(:)

                        ! end Volume calculation based on mixed product
                        !------------------------------------------
                     case (VOL_METHOD_DIVERGENCE) 
                        !------------------------------------------
                        ! Volume calculation based divergence theorem
                        do icoord=1_ip,ndime
                           r_1(icoord) = coord(icoord, ipoint_1) + displ(icoord, ipoint_1, 1) - &
                              cavities(icavi) % valves(ivalve) % centroid(icoord)
                           r_2(icoord) = coord(icoord, ipoint_2) + displ(icoord, ipoint_2, 1) - &
                              cavities(icavi) % valves(ivalve) % centroid(icoord)
                        end do
                        
                        !triangle centroid
                        centroid = (1.0_rp/3.0_rp)* &
                           (cavities(icavi) % valves(ivalve) % centroid(:) + &
                           coord(:, ipoint_1) + displ(:, ipoint_1, 1) + &
                           coord(:, ipoint_2) + displ(:, ipoint_2, 1))

                        !this normal will be pointing outwards from the cavity
                        cavities(icavi) % plane_normal(:) = maths_cross_product( r_1(:), r_2(:), 3_ip )
                        norm = sqrt(sum(cavities(icavi) % plane_normal(:)**2_ip))
                        area = 0.5_rp * norm

                        if( norm>epsilon(1.0_rp) ) then
                           cavities(icavi) % plane_normal(:) = cavities(icavi) % plane_normal(:) / norm
                        end if
                  
                        !here we flip the normal to point inwards
                        cavities(icavi) % volume(ITER_K) =  cavities(icavi) % volume(ITER_K) - &
                           dot_product(centroid(:), cavities(icavi) % plane_normal(:) ) * area

                        ! end Volume calculation based divergence theorem
                        !------------------------------------------
                     case default 
                        call runend("Unknown volume calculation method "//trim(intost(cavities(icavi) % volume_method))//" for cavity "//trim(intost(icavi)))
                     end select

                  END IF
               END DO
            end if
         END DO

      END IF

      call PAR_SUM(cavities(icavi) % volume(ITER_K))

      select case (cavities(icavi) % volume_method )
      case(VOL_METHOD_DIVERGENCE)
         cavities(icavi) % volume(ITER_K) = (1.0_rp/3.0_rp)*cavities(icavi) % volume(ITER_K)            
      case(VOL_METHOD_MIXEDPROD)
         cavities(icavi) % volume(ITER_K) = (1.0_rp/6.0_rp)*cavities(icavi) % volume(ITER_K)
      case default 
         call runend("Unknown volume calculation method "//trim(intost(cavities(icavi) % volume_method))//" for cavity "//trim(intost(icavi)))
      end select

   contains
      real(rp) pure function mixed_product(a, b, c)
         implicit none
         real(rp), intent(in) :: a(3),b(3),c(3)

         !a1 a2 a3
         !b1 b2 b3
         !c1 c2 c3
         mixed_product = &
            a(1)*( b(2)*c(3) - b(3)*c(2) ) + &
            a(2)*( b(3)*c(1) - b(1)*c(3) ) + &
            a(3)*( b(1)*c(2) - b(2)*c(1) )

      end function
   end subroutine 


  !--------------------------------------------------
  !
  ! Member acess for tests
  ! DO NOT USE OUTSIDE OF TESTS
  !
  subroutine test_set_cavity_phase( icavi, phase )
      implicit none
      integer(ip), intent(in) :: icavi, phase

      cavities(icavi) % phase = phase
  end subroutine test_set_cavity_phase

  subroutine test_set_cavity_phase_counter( icavi, value )
      implicit none
      integer(ip), intent(in) :: icavi, value

      cavities(icavi) % phase_counter = value
  end subroutine test_set_cavity_phase_counter

  subroutine test_set_phase41_transition_method( value )
      implicit none
      integer(ip), intent(in) :: value

      phase41_trans_method = value
  end subroutine test_set_phase41_transition_method

  subroutine test_set_volcai( value )
      implicit none
      real(rp), intent(in), dimension(3) :: value

      volcai(:) = value(:)
  end subroutine test_set_volcai

  subroutine test_set_cavity_dvol( icavi, value )
      use def_master, only: TIME_N
      implicit none
      integer(ip), intent(in) :: icavi
      real(rp), intent(in)    :: value

      cavities(icavi) % dvol(TIME_N) = value
  end subroutine test_set_cavity_dvol


end module mod_sld_cardiac_cycle
