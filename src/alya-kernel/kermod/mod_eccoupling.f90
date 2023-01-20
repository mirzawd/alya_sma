!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-------------------------------------------------------------------------------
!> @addtogroup Electromechanical coupling
!> @{
!> @authors Adria Quintanas-Corominas : adria.quintanas@bsc.es
!> @authors Alfonso Santiago          : alfonso.santiago@bsc.es
!> @author  Constantine Butakoff
!> @date    March, 2020
!> @}
!------------------------------------------------------------------------------
module mod_eccoupling
    ! =========================================================================
    use def_kintyp,   only                    :  ip, rp, lg
    use def_master,   only                    :  cutim, dtime
    use def_master,   only                    :  modul, kfl_paral
    use mod_parall,   only                    :  PAR_MY_UNIVERSE_RANK
    use mod_messages, only                    :  messages_live
    use def_domain,   only                    :  nmate
    use mod_biofibers, only                   :  kfl_biofibers, biofibers
    ! ----------------------------------------
    implicit none
    ! -< DIMENSIONS >-------------------------
    integer(ip), parameter                    :: nprop_ecc = 22_ip
    integer(ip), parameter                    :: nsdvs_ecc = 8_ip  
    integer(ip), parameter                    :: ncomp_ecc = 2_ip
    ! -< CELL MODELS >------------------------
    integer(ip), parameter                    :: EXMSLD_CELL_NOMODEL     =  0_ip
    integer(ip), parameter                    :: EXMSLD_CELL_FITZHUGH    =  1_ip
    integer(ip), parameter                    :: EXMSLD_CELL_FENTON      =  2_ip
    integer(ip), parameter                    :: EXMSLD_CELL_TENTUSCHER  =  3_ip
    integer(ip), parameter                    :: EXMSLD_CELL_TT2006      =  4_ip
    integer(ip), parameter                    :: EXMSLD_CELL_OHARA       =  5_ip
    integer(ip), parameter                    :: EXMSLD_CELL_OHARA_INAPA =  6_ip
    integer(ip), parameter                    :: EXMSLD_CELL_SCATRIA     =  7_ip
    integer(ip), parameter                    :: EXMSLD_CELL_SCVENTRI    =  8_ip        
    integer(ip), parameter                    :: EXMSLD_CELL_TORORD      =  9_ip        
    integer(ip), parameter                    :: EXMSLD_CELL_COURTE      = 10_ip     
    integer(ip), parameter                    :: EXMSLD_CELL_MAXMODELID  = EXMSLD_CELL_COURTE     

    ! -< ELECTRO-MECHANICAL MODELS >----------
    integer(ip), parameter                    :: EXMSLD_EMEC_NO_EMEC         = 0_ip
    integer(ip), parameter                    :: EXMSLD_EMEC_HUNTER          = 1_ip ! simplified hunter model
    integer(ip), parameter                    :: EXMSLD_EMEC_LAND            = 3_ip
    integer(ip), parameter                    :: EXMSLD_EMEC_LAND_BIDIR      = 4_ip
    integer(ip), parameter                    :: EXMSLD_EMEC_HUNTER_ORIGINAL = 5_ip ! original hunter model with several state variables
    ! -< ELECTRO-MECHANICAL PARAMETERS >------
    real(rp),    parameter, private           :: lamda_thr = 1.2_rp
    ! ----< GENERAL VARIABLES >---------------
    logical(lg),           protected          :: kfl_exmsld_ecc        = .false.
    logical(lg),           protected          :: kfl_exmsld_moncod_ecc = .false.
    logical(lg),           protected          :: kfl_exmsld_mulcod_ecc = .false.
    logical(lg),           protected          :: kfl_exmsld_lands_ecc  = .false.
    logical(lg),           protected          :: kfl_print             = .false.
    logical(lg),           protected          :: kfl_flags_init        = .false.
    logical(lg),           private            :: kfl_state_initialised = .false.
    ! ----< VARIABLES >-----------------------
    real(rp),    pointer                      :: troponin_ecc(:,:)   => null() 
    real(rp),    pointer, protected           :: props_ecc(:,:)      => null() 
    real(rp),    pointer, protected           :: ortho_ecc(:,:)      => null()
    real(rp),    pointer                      :: state_ecc(:,:,:)    => null()
    real(rp),    pointer                      :: strch_ecc(:,:,:)    => null()    
    real(rp),    pointer, protected           :: calcium_ecc(:)      => null() ! = vconc(1,:,1)  * 1000
    real(rp),    pointer, protected           :: displ_lagr_ecc(:,:,:) => null()
    real(rp),    pointer, protected           :: fiber_ecc(:,:)      => null()
    real(rp),    pointer, protected           :: fibe0_ecc(:,:)      => null()
    real(rp),    pointer, protected           :: celty_ecc(:)        => null()
    real(rp),    pointer, protected           :: cell_ca0_ecc(:,:)   => null()
    real(rp),    pointer, protected           :: cell_ca50_ecc(:,:)  => null() 
    real(rp),    pointer, private             :: cell_ca50_permate(:)=> null() 
    real(rp),    pointer, protected           :: cell_caX_ecc(:)     => null()
    real(rp),    pointer, protected           :: gp_Ca_ecc(:)        => null()
    real(rp),    pointer, protected           :: gp_Ca0_ecc(:)       => null()
    real(rp),    pointer, protected           :: gp_Ca50_ecc(:)      => null()
    real(rp),    pointer, protected           :: gp_props_ecc(:,:)   => null()
    integer(ip), pointer, protected           :: kfl_cellmod(:)      => null() ! Electrophysiology cell model  EXMEDI - SOLIDZ
    integer(ip), pointer, protected           :: kfl_eccty(:)        => null() ! Electromechanical coupling type (per material) EXMEDI - SOLIDZ
    logical(lg), pointer, protected           :: kfl_force_ca50_calculation(:) => null()
    integer(ip), pointer, private             :: kfl_ca50_provided(:) => null()
    integer(ip),          private             :: MPI_RANK_MASTER_EXM
    integer(ip),          private             :: MPI_RANK_MASTER_SLD
    real(rp),             private             :: coupling_delay_ecc = 0.0_rp  
    ! This should be made private
    real(rp)                                  :: volcai_ecc     = 0.0_rp
    integer(ip),          protected           :: ncelltypes_ecc = 1_ip

    ! ----------------------------------------
    abstract interface
        subroutine  get_active_traction( props, lmb, Ca, T, dTdlmb, sdv )
            import                            :: ip, rp, lg, dtime
            real(rp),          intent(in)     :: props(:)  ! Properties
            real(rp),          intent(in)     :: lmb       ! Lambda (stretches)
            real(rp),          intent(in)     :: Ca        ! Calcium
            real(rp),          intent(inout)  :: T         ! Active traction
            real(rp),          intent(inout)  :: dTdlmb    ! Derivative active traction w.r.t lamda
            real(rp),          intent(inout)  :: sdv(:,:)  ! State variables
        end subroutine get_active_traction
    end interface
    procedure(get_active_traction),  pointer , private :: get_active_traction_generic
    ! ----------------------------------------
    public :: exm_sld_ecc_initialize_troponin
    public :: exm_sld_ecc_interchange_troponin
    public :: eccou_show_warning
    public :: eccou_manage_variables
    public :: eccou_initialise_flags
    public :: eccou_assign_fields
    public :: eccou_interpolate_variables_at_gp
    public :: eccou_allocate_memory
    public :: eccou_read_data
    public :: eccou_send_data
    public :: eccou_manage_arrays
    public :: eccou_manage_restart
    public :: eccou_set_state_variables
    public :: eccou_save_log
    public :: eccou_set_cellmod
    public :: eccou_get_cellmod
    public :: eccou_coupling_active
    public :: eccou_do_volume_integral
    public :: eccou_set_ca
    public :: eccou_set_ca0
    public :: eccou_set_ca50
    public :: eccou_set_ncelltypes
    private :: calc_pC50_ref
    private :: eccou_ca50_setdefaults
    ! ----------------------------------------

    contains


    subroutine eccou_ca50_setdefaults()
      use def_master, only : intost
      implicit None
      integer(ip) :: icellt, imate

      !switch  ca50 to pc50 if required, after reading eccoupling section
      do imate = 1,nmate
        cell_ca50_ecc(:,imate) = cell_ca50_permate(imate)

        !if not provided, load ca50 from the defaults
        if (kfl_ca50_provided(imate) == 0_ip) then 
            if(      kfl_eccty(imate) == EXMSLD_EMEC_HUNTER ) then
                cell_ca50_ecc(:,imate) = props_ecc( 2,imate)
            else if( kfl_eccty(imate) == EXMSLD_EMEC_HUNTER_ORIGINAL )then
                cell_ca50_ecc(:,imate) = props_ecc( 14,imate) 
            else if( kfl_eccty(imate) == EXMSLD_EMEC_LAND .or. &
                    kfl_eccty(imate) == EXMSLD_EMEC_LAND_BIDIR )then
                cell_ca50_ecc(:,imate) = props_ecc( 3,imate)
            endif
        end if

        !if ca50 is provided and needs conversion to pc50_ref for the hunter model
        if ( (kfl_ca50_provided(imate) > 0_ip) ) then
            call messages_live("ECCOUPLING: USING CA50 FROM ker.dat TO CALCULATE pC50_ref FOR MATERIAL "//trim(intost(imate)))
            do icellt = 1,ncelltypes_ecc
                cell_ca50_ecc(icellt,imate) = calc_pC50_ref( cell_ca50_ecc(icellt,imate) )
            end do
        end if
      end do
    end subroutine


    subroutine eccou_set_ncelltypes(n)
      implicit none
      integer(ip), intent(in) :: n
      ncelltypes_ecc = n
    end subroutine

    ! Return .true. if the coupling delay is not active any more
    logical(lg) function eccou_coupling_active()
        use def_master, only : cutim
        implicit none

        eccou_coupling_active = ( cutim >= coupling_delay_ecc )
    end function eccou_coupling_active

    logical(lg) function eccou_imate_is_active(imate)
        use def_master,    only :  cutim
        integer(ip), intent(in) :: imate
        if( kfl_cellmod(imate) /= EXMSLD_CELL_NOMODEL .and. &
            kfl_eccty(imate) /= EXMSLD_EMEC_NO_EMEC   .and. &
            cutim >= coupling_delay_ecc )then
            eccou_imate_is_active = .true.
        else
            eccou_imate_is_active = .false.
        endif
    end function eccou_imate_is_active


    ! -------------------------------------------------------------------------
    !> @author  Adria Quintanas-Corominas (adria.quintanas@bsc.es)
    !> @author  Alfonso Santiago (alfonso.santiago@bsc.es)
    ! -------------------------------------------------------------------------
    subroutine eccou_set_cellmod(imate, cellmod_id)
        implicit none
        integer(ip), intent(in) :: imate, cellmod_id
        kfl_cellmod(imate) = cellmod_id
    end subroutine eccou_set_cellmod


    integer(ip) function eccou_get_cellmod(imate)
        implicit none
        integer(ip), intent(in) :: imate
        eccou_get_cellmod = kfl_cellmod(imate) 
    end function eccou_get_cellmod


    ! -------------------------------------------------------------------------
    !> @author  Adria Quintanas-Corominas (adria.quintanas@bsc.es)
    !> @author  Alfonso Santiago (alfonso.santiago@bsc.es)
    ! -------------------------------------------------------------------------
    subroutine eccou_show_warning()
        use mod_messages,         only :  livinf, messages_live
        ! -------------------------------
        call messages_live("------------------------------------------------------------","WARNING")
        call messages_live(" A deprecated syntax is used for COUPLING_SCHEME in SLD.DAT ",'WARNING')
        call messages_live(" Instead, use the following syntax in KER.DAT :             ","WARNING")
        call messages_live("   PHYSICAL_PROBLEM                                         ","WARNING")
        call messages_live("     FIBERS FIELD 1                  $ (mandatory ID_FIELD) ","WARNING")
        call messages_live("     CELL_TYPE FIELD 2               $ (mandatory ID_FIELD) ","WARNING")
        call messages_live("     ECCOUPLING                                             ","WARNING")
        call messages_live("         MATERIAL:    1              $ (mandatory)          ","WARNING")  
        call messages_live("         MODEL:       HUNTER | LAND  $ (mandatory)          ","WARNING")  
        call messages_live("     END_ECCOUPLING                                         ","WARNING")
        call messages_live("   END_PHYSICAL_PROBLEM                                     ","WARNING")
        call messages_live("------------------------------------------------------------","WARNING")
        call runend('SLD_REAPHY: READ THE PREVIOUS WARNINGS')
        ! -------------------------------
    end subroutine eccou_show_warning

    ! -------------------------------------------------------------------------
    !> @author  Adria Quintanas-Corominas (adria.quintanas@bsc.es)
    !> @author  Alfonso Santiago (alfonso.santiago@bsc.es)
    !> @TODO    In the mono-code coupling this routine is called two times.
    ! -------------------------------------------------------------------------
    subroutine eccou_initialise_flags()
        use def_master,            only : coupling
        use def_coupli,            only : mcoup, coupling_type
        ! -------------------------------
        implicit none
        ! -------------------------------
        integer(ip)                    :: ii
        ! -------------------------------
        ! Mono-code coupling?
         kfl_exmsld_moncod_ecc = .true.
         kfl_exmsld_mulcod_ecc = .false.
        ! Multi-code coupling?
        if( mcoup > 0_ip )then
            do ii = 1, mcoup
                if( coupling_type(ii) % variable == 'ECCOU' .or. &
                    coupling_type(ii) % variable == 'DISPL' .or. &
                    coupling_type(ii) % variable == 'ALEFO' )then
                    kfl_exmsld_moncod_ecc = .false.
                    kfl_exmsld_mulcod_ecc = .true.
                endif
            enddo
        endif
        ! -------------------------------
    end subroutine eccou_initialise_flags 

    ! -------------------------------------------------------------------------
    !> @author  Adria Quintanas-Corominas (adria.quintanas@bsc.es)
    !> @author  Alfonso Santiago (alfonso.santiago@bsc.es)
    ! -------------------------------------------------------------------------
    subroutine eccou_assign_fields()
        use def_domain, only : xfiel
        use def_kermod, only : kfl_celltype_fun
        use def_master, only : intost

        ! -------------------------------
        implicit none
        ! -------------------------------
        if( kfl_celltype_fun < 0_ip ) then
            call messages_live("ECCOUPLING IS USING FIELD "//trim(intost(-kfl_celltype_fun))//" TO READ HETEROGENEOUS CELL TYPE")
            celty_ecc => xfiel(-kfl_celltype_fun) % a(1,:,1)
        end if

        !if( kfl_fiber_long_fun    < 0_ip ) then 
        !    call messages_live("ECCOUPLING IS USING FIELD "//trim(intost(-kfl_fiber_long_fun))//" TO READ FIBERS")
        !    fibe0_ecc => xfiel(-kfl_fiber_long_fun) % a(:,:,1)
        !end if
        ! -------------------------------
    end subroutine eccou_assign_fields 

    ! -------------------------------------------------------------------------
    !> @author  Adria Quintanas-Corominas (adria.quintanas@bsc.es)
    !> @author  Alfonso Santiago (alfonso.santiago@bsc.es)
    ! -------------------------------------------------------------------------
    subroutine eccou_manage_arrays( &
        itask,  wtask )
        ! -------------------------------
        use def_domain,           only :  npoin, nelem, mgaus, ndime
        use mod_arrays,           only :  arrays, arrays_number, arrays_register
        ! -------------------------------
        implicit none
        ! -------------------------------
        integer(ip),      intent(in)           :: itask
        character(len=*), intent(in), optional :: wtask
        ! -------------------------------
        integer(ip),      parameter            :: REGISTER = 0_ip
        integer(ip),      parameter            :: DOJOB    = 1_ip
        ! -------------------------------
        select case(itask)

            case( REGISTER )
                call arrays_register( (/ 'CALCI', 'SCALA', 'NPOIN', 'SECON' /))
                call arrays_register( (/ 'EMSVD', 'SCALA', 'NELEM', 'PRIMA' /), state_ecc,    ENTITY_POSITION= 2_ip, TIME_POSITION= 3_ip )
                call arrays_register( (/ 'EMSTR', 'SCALA', 'NELEM', 'PRIMA' /), strch_ecc,    ENTITY_POSITION= 2_ip, TIME_POSITION= 3_ip )
                call arrays_register( (/ 'EMTRO', 'SCALA', 'NPOIN', 'PRIMA' /), troponin_ecc, ENTITY_POSITION= 1_ip, TIME_POSITION= 2_ip )
                call arrays_register( (/ 'EMFIB', 'SCALA', 'NELEM', 'PRIMA' /), fiber_ecc,    ENTITY_POSITION= 2_ip)

             case( DOJOB )
               call arrays( arrays_number('EMSVD'), wtask, state_ecc,    nsdvs_ecc*mgaus, nelem, ncomp_ecc)
               call arrays( arrays_number('EMSTR'), wtask, strch_ecc,    mgaus,           nelem, ncomp_ecc)
               call arrays( arrays_number('EMTRO'), wtask, troponin_ecc, npoin,           ncomp_ecc)
               ! TODO: Unify gpfib and fiber_ecc by creating a new way to register (nsdv, ngaus, nelem) arrays.
               if( wtask == 'WRITE RESTART') call set_fiber_ecc() 
               call arrays( arrays_number('EMFIB'), wtask, fiber_ecc,    ndime*mgaus,     nelem)
               if( wtask == 'READ RESTART' ) call get_fiber_ecc()

               ! To avoid re-initialise the state varibles
               if( wtask == 'READ RESTART' ) kfl_state_initialised = .true.

        end select
        ! -------------------------------
        contains 

            subroutine set_fiber_ecc()
                use def_master, only : gpfib
                use def_domain, only : ndime, nelem, mgaus
                implicit none
                integer(ip)         :: ielem, ii, jj, kk 
                do ielem = 1, nelem
                    do ii = 1, mgaus
                        do kk = 1, ndime 
                           jj = (ii - 1_ip)*nsdvs_ecc + kk
                           fiber_ecc(jj,ielem) = gpfib(kk,ii,ielem)
                        enddo
                    enddo
                enddo    
            end subroutine set_fiber_ecc

            subroutine get_fiber_ecc()
                use def_master, only : gpfib 
                use def_domain, only : ndime, nelem, mgaus
                implicit none
                integer(ip)         :: ielem, ii, jj, kk
                do ielem = 1, nelem
                    do ii = 1, mgaus
                        do kk = 1, ndime
                           jj = (ii - 1_ip)*nsdvs_ecc + kk
                           gpfib(kk,ii,ielem) = fiber_ecc(jj,ielem)
                        enddo
                    enddo
                enddo
            end subroutine get_fiber_ecc

    end subroutine eccou_manage_arrays

    ! ---------------------------------------------------------------------
    !> @author   Adria Quintanas-Corominas (adria.quintanas@bsc.es)
    !> @author   Alfonso Santiago (alfonso.santiago@bsc.es)
    !> @author   Guillaume Houzeaux (guillaume.houzeaux@bsc.es)
    ! ---------------------------------------------------------------------
    subroutine eccou_manage_restart( itask )
        use mod_restart, only : restart_add
        ! -------------------------------
        implicit none
        ! -------------------------------
        integer(ip), intent(in) :: itask

        if( ncelltypes_ecc > 0 .and. kfl_exmsld_ecc ) then
           call restart_add(cell_ca0_ecc, 'cell_ca0_ecc' )
           call restart_add(cell_ca50_ecc,'cell_ca50_ecc')
           call restart_add(cell_caX_ecc, 'cell_caX_ecc ')
        end if
        
    end subroutine eccou_manage_restart

    ! -------------------------------------------------------------------------
    !> @author  Adria Quintanas-Corominas (adria.quintanas@bsc.es)
    !> @author  Alfonso Santiago (alfonso.santiago@bsc.es)
    ! -------------------------------------------------------------------------
    subroutine eccou_allocate_memory( &
        itask )
        ! -------------------------------
        use def_master,           only :  modul, mem_modul, INOTMASTER, gpfib, displ
        use def_domain,           only :  npoin, nelem, mgaus, nmate, ndime
        use mod_memory,           only :  memory_alloca, memory_deallo
        ! -------------------------------
        implicit none
        ! -------------------------------
        integer(ip), intent(in)        :: itask
        integer(ip), parameter         :: KERNEL_READAT          = 0_ip
        integer(ip), parameter         :: KERNEL_MEMALL          = 1_ip
        integer(ip), parameter         :: KERNEL_MEMALL_NOCOUP   = 2_ip
        integer(ip), parameter         :: EXMEDI_MEMALL          = 3_ip
        integer(ip), parameter         :: KERNEL_READAT_FINALIZE = 4_ip
        ! -------------------------------
        select case(itask)
            
            case( KERNEL_READAT )
                ! Allocate memmory for the variables read in the ker.dat
                call memory_alloca(mem_modul(1:2,modul), 'PROPS_ECC'                 , 'mod_eccoupling', props_ecc    , nprop_ecc, nmate); props_ecc(:,:) = 0.0_rp
                call memory_alloca(mem_modul(1:2,modul), 'ORTHO_ECC'                 , 'mod_eccoupling', ortho_ecc    ,      2_ip, nmate); ortho_ecc(:,:) = 0.0_rp
                call memory_alloca(mem_modul(1:2,modul), 'kfl_ca50_provided'         , 'mod_eccoupling', kfl_ca50_provided , nmate); kfl_ca50_provided(:) = 0_ip
                call memory_alloca(mem_modul(1:2,modul), 'kfl_force_ca50_calculation', 'mod_eccoupling', kfl_force_ca50_calculation , nmate); kfl_force_ca50_calculation(:) = .FALSE.
                call memory_alloca(mem_modul(1:2,modul), 'kfl_cellmod'               , 'mod_eccoupling', kfl_cellmod  , nmate); kfl_cellmod(:) = EXMSLD_CELL_NOMODEL
                call memory_alloca(mem_modul(1:2,modul), 'kfl_eccty'                 , 'mod_eccoupling', kfl_eccty    , nmate); kfl_eccty(:) = real(EXMSLD_EMEC_NO_EMEC,rp)
                call memory_alloca(mem_modul(1:2,modul), 'cell_ca50_permate'         , 'mod_eccoupling', cell_ca50_permate, nmate); cell_ca50_permate(:) = 0.0_rp
                call memory_alloca(mem_modul(1:2,modul), 'CELL_CAX_ECC'              , 'mod_eccoupling', cell_caX_ecc , nmate); cell_caX_ecc(:) = 1.0_rp

                kfl_ca50_provided(:) = 0_ip

            case( KERNEL_READAT_FINALIZE )
                call memory_alloca(mem_modul(1:2,modul), 'CELL_CA0_ECC' , 'mod_eccoupling', cell_ca0_ecc , ncelltypes_ecc, nmate); cell_ca0_ecc(:,:) = 0.0_rp
                call memory_alloca(mem_modul(1:2,modul), 'CELL_CA50_ECC', 'mod_eccoupling', cell_ca50_ecc, ncelltypes_ecc, nmate); cell_ca50_ecc(:,:) = 0.0_rp

                call eccou_ca50_setdefaults()

            case( KERNEL_MEMALL )
                ! If vconc has been associated means that it is a monocoupling. Therefore, calcium_ecc only needs to point
                ! the correct position of vconc. Otherwise, the allocation of memmory for calcium_ecc is required for the 
                ! solidz size.  

              !TODO: this will fail during repartitioning, but I have no idea right now hot o deal with arrays per gauss point
                call memory_alloca(mem_modul(1:2,modul), 'CALCIUM_ECC'  , 'mod_eccoupling', calcium_ecc  , npoin)
                call memory_alloca(mem_modul(1:2,modul), 'GPFIB'        , 'mod_eccoupling', gpfib        , ndime , mgaus, nelem)
                call memory_alloca(mem_modul(1:2,modul), 'GP_CA_ECC'    , 'mod_eccoupling', gp_Ca_ecc    , mgaus)
                call memory_alloca(mem_modul(1:2,modul), 'GP_CA0_ECC'   , 'mod_eccoupling', gp_Ca0_ecc   , mgaus)
                call memory_alloca(mem_modul(1:2,modul), 'GP_CA50_ECC'  , 'mod_eccoupling', gp_Ca50_ecc  , mgaus)
                call memory_alloca(mem_modul(1:2,modul), 'GP_PROPS_ECC' , 'mod_eccoupling', gp_props_ecc , nprop_ecc, mgaus)
                ! Initialise 
                if( INOTMASTER )then
                    state_ecc(:,:,:) = 0.0_rp
                    strch_ecc(:,:,:) = 1.0_rp
                endif

            case( KERNEL_MEMALL_NOCOUP )
                ! 
                call memory_alloca(mem_modul(1:2,modul), 'kfl_cellmod' , 'mod_eccoupling', kfl_cellmod, nmate); kfl_cellmod(:) = EXMSLD_CELL_NOMODEL
                
            case( EXMEDI_MEMALL )
                ! EXMEDI when multi code coupling
                call memory_alloca(mem_modul(1:2,modul), 'DISPL_LAGR_ECC', 'mod_eccoupling', displ_lagr_ecc, ndime, npoin,1_ip)
                if( INOTMASTER )displ_lagr_ecc=0.0_rp
                displ=>displ_lagr_ecc

        end select
        ! -------------------------------
    end subroutine eccou_allocate_memory


    subroutine eccou_save_log()
        use def_master, only : momod, intost, INOTSLAVE
        ! -------------------------------
        implicit none
        ! -------------------------------
        integer(ip)                    :: imate, icell
        ! -------------------------------
        if (INOTSLAVE) then
            write( momod(modul)%lun_outpu, * ) NEW_LINE('a')
            write( momod(modul)%lun_outpu, * ) "INITIAL COUPLING PARAMETERS FOR EP PROBLEM" 
            do imate = 1,nmate
                write( momod(modul)%lun_outpu, * ) "   Material "//trim(intost(imate))
                if(      kfl_eccty(imate) == EXMSLD_EMEC_HUNTER  ) then
                    write( momod(modul)%lun_outpu, * ) "      Model    = HUNTER"
                    write( momod(modul)%lun_outpu, "(a,e16.8E3)" ) "      Cam      = ", props_ecc( 1,imate)
                    write( momod(modul)%lun_outpu, "(a,e16.8E3)" ) "      Time     = ", props_ecc( 3,imate) 
                    write( momod(modul)%lun_outpu, "(a,e16.8E3)" ) "      Eta      = ", props_ecc( 4,imate) 
                    write( momod(modul)%lun_outpu, "(a,e16.8E3)" ) "      Hillc    = ", props_ecc( 5,imate) 
                    write( momod(modul)%lun_outpu, "(a,e16.8E3)" ) "      k        = ", props_ecc( 6,imate) 
                    write( momod(modul)%lun_outpu, "(a,e16.8E3)" ) "      Traction = ", props_ecc( 7,imate) 

                else if (kfl_eccty(imate) == EXMSLD_EMEC_HUNTER_ORIGINAL )then
                    write( momod(modul)%lun_outpu, * ) "      Model    = HUNTER ORIGINAL"
                    write( momod(modul)%lun_outpu, "(a,e16.8E3)" ) "      rho_0    = ", props_ecc( 1,imate) 
                    write( momod(modul)%lun_outpu, "(a,e16.8E3)" ) "      rho_1    = ", props_ecc( 2,imate) 
                    write( momod(modul)%lun_outpu, "(a,e16.8E3)" ) "      Ca_max   = ", props_ecc( 4,imate) 
                    write( momod(modul)%lun_outpu, "(a,e16.8E3)" ) "      Ca_bmax  = ", props_ecc( 5,imate) 
                    write( momod(modul)%lun_outpu, "(a,e16.8E3)" ) "      tau_Ca   = ", props_ecc( 6,imate) 
                    write( momod(modul)%lun_outpu, "(a,e16.8E3)" ) "      gamma    = ", props_ecc( 8,imate) 
                    write( momod(modul)%lun_outpu, "(a,e16.8E3)" ) "      n_ref    = ", props_ecc( 9,imate) 
                    write( momod(modul)%lun_outpu, "(a,e16.8E3)" ) "      alpha_0  = ", props_ecc(10,imate) 
                    write( momod(modul)%lun_outpu, "(a,e16.8E3)" ) "      beta_0   = ", props_ecc(11,imate) 
                    write( momod(modul)%lun_outpu, "(a,e16.8E3)" ) "      beta_1   = ", props_ecc(12,imate) 
                    write( momod(modul)%lun_outpu, "(a,e16.8E3)" ) "      beta_2   = ", props_ecc(13,imate) 
                    write( momod(modul)%lun_outpu, "(a,e16.8E3)" ) "      T_ref    = ", props_ecc(15,imate) 
                    write( momod(modul)%lun_outpu, "(a,e16.8E3)" ) "      a        = ", props_ecc(16,imate) 
                    write( momod(modul)%lun_outpu, "(a,e16.8E3)" ) "      A_1      = ", props_ecc(17,imate) 
                    write( momod(modul)%lun_outpu, "(a,e16.8E3)" ) "      A_2      = ", props_ecc(18,imate) 
                    write( momod(modul)%lun_outpu, "(a,e16.8E3)" ) "      A_3      = ", props_ecc(19,imate) 
                    write( momod(modul)%lun_outpu, "(a,e16.8E3)" ) "      alpha_1  = ", props_ecc(20,imate) 
                    write( momod(modul)%lun_outpu, "(a,e16.8E3)" ) "      alpha_2  = ", props_ecc(21,imate) 
                    write( momod(modul)%lun_outpu, "(a,e16.8E3)" ) "      alpha_3  = ", props_ecc(22,imate) 

                else if( kfl_eccty(imate) == EXMSLD_EMEC_LAND .or. &
                           kfl_eccty(imate) == EXMSLD_EMEC_LAND_BIDIR )then

                    if (kfl_eccty(imate) == EXMSLD_EMEC_LAND_BIDIR) then
                        write( momod(modul)%lun_outpu, * ) "      Model    = LAND BIDIRECTIONAL"
                    else
                        write( momod(modul)%lun_outpu, * ) "      Model    = LAND"
                    end if

                    write( momod(modul)%lun_outpu, "(a,e16.8E3)" ) "      k_TRPN     = ", props_ecc( 1,imate) 
                    write( momod(modul)%lun_outpu, "(a,e16.8E3)" ) "      n_TRPN     = ", props_ecc( 2,imate) 
                    write( momod(modul)%lun_outpu, "(a,e16.8E3)" ) "      k_u        = ", props_ecc( 4,imate) 
                    write( momod(modul)%lun_outpu, "(a,e16.8E3)" ) "      n_tm       = ", props_ecc( 5,imate) 
                    write( momod(modul)%lun_outpu, "(a,e16.8E3)" ) "      TRPN_50    = ", props_ecc( 6,imate) 
                    write( momod(modul)%lun_outpu, "(a,e16.8E3)" ) "      k_uw       = ", props_ecc( 7,imate) 
                    write( momod(modul)%lun_outpu, "(a,e16.8E3)" ) "      k_ws       = ", props_ecc( 8,imate) 
                    write( momod(modul)%lun_outpu, "(a,e16.8E3)" ) "      r_w        = ", props_ecc( 9,imate) 
                    write( momod(modul)%lun_outpu, "(a,e16.8E3)" ) "      r_s        = ", props_ecc(10,imate) 
                    write( momod(modul)%lun_outpu, "(a,e16.8E3)" ) "      y_s        = ", props_ecc(11,imate) 
                    write( momod(modul)%lun_outpu, "(a,e16.8E3)" ) "      y_w        = ", props_ecc(12,imate) 
                    write( momod(modul)%lun_outpu, "(a,e16.8E3)" ) "      phi        = ", props_ecc(13,imate) 
                    write( momod(modul)%lun_outpu, "(a,e16.8E3)" ) "      A_eff      = ", props_ecc(14,imate) 
                    write( momod(modul)%lun_outpu, "(a,e16.8E3)" ) "      beta_0     = ", props_ecc(15,imate) 
                    write( momod(modul)%lun_outpu, "(a,e16.8E3)" ) "      beta_1     = ", props_ecc(16,imate) 
                    write( momod(modul)%lun_outpu, "(a,e16.8E3)" ) "      T_ref      = ", props_ecc(17,imate) 
                    write( momod(modul)%lun_outpu, "(a,e16.8E3)" ) "      Ca0        = ", props_ecc(18,imate) 
                    write( momod(modul)%lun_outpu, "(a,e16.8E3)" ) "      beta_2     = ", props_ecc(19,imate)
                    write( momod(modul)%lun_outpu, "(a,e16.8E3)" ) "      Cat50      = ", props_ecc(20,imate)
               
                else if( kfl_eccty(imate) == EXMSLD_EMEC_NO_EMEC )then
                    write( momod(modul)%lun_outpu, * ) "      Model    = NONE"
                else
                    call runend('ERROR : ECCOU_READ_DATA : Unknown coupling model')
                endif

                do icell = 1,ncelltypes_ecc
                    write( momod(modul)%lun_outpu, * ) "      Cell type "//trim(intost(icell))//" Ca0  = ", cell_ca0_ecc(icell,imate)
                    write( momod(modul)%lun_outpu, "(a,e16.8E3)" ) "      Cell type "//trim(intost(icell))//" Ca50 = ", 10.0_rp**(6.0_rp - cell_ca50_ecc(icell,imate))
                    write( momod(modul)%lun_outpu, "(a,e16.8E3)" ) "      Cell type "//trim(intost(icell))//" pC50_ref = ", cell_ca50_ecc(icell,imate)
                end do
            end do


            flush(momod(modul)%lun_outpu)
        end if
    end subroutine eccou_save_log


    ! -------------------------------------------------------------------------
    !> @author  Adria Quintanas-Corominas (adria.quintanas@bsc.es)
    !> @author  Alfonso Santiago (alfonso.santiago@bsc.es)
    ! -------------------------------------------------------------------------
    subroutine eccou_read_data()
        use def_master,           only :  intost, retost
        use def_inpout,           only :  words, exists, getint, getrea, getcha
        use mod_ecoute,           only :  ecoute
        use mod_messages,         only :  messages_live
        use mod_messages,         only :  livinf

        ! -------------------------------
        implicit none
        ! -------------------------------
        integer(ip)                    :: imate
        real(rp)                       :: timec_ecc
        ! -------------------------------
        call messages_live('   READING ELECTRO-MECANICHAL COUPLING (ECCOU)')
  
        ! Initialise some variables
        imate = 1_ip
        props_ecc(:,:)  = 0.0_rp
        ortho_ecc(:,:)  = 0.0_rp
        cell_caX_ecc(:) = 1.0_rp
        coupling_delay_ecc = 0.0_rp

        ! Initialise coupling flag
        kfl_exmsld_ecc = .true. 

        call ecoute('eccou_read_data')
        coupling_sld_exm: &
        do while( words(1) /= 'ENDEC' )

            if(      words(1) == 'MATER' )then
  
                imate = getint('MATER=',1_ip,'#Current material')
  
            else if( words(1) == 'DELAY' ) then
                coupling_delay_ecc = getrea("DELAY",0.0_rp,"!Coupling delay")
                call messages_live("ECCOUPLING: COUPLING WILL BE DELAYED UNTIL "//trim(retost(coupling_delay_ecc,'(F10.4)'))//"s","WARNING")

            else if( words(1) == 'CONST' .or. words(1) == 'MODEL' )then
                ! Select model according to the second word
                if(      words(2) == 'HUNTE' )then

                    if(  words(3) == 'ORIGI' )then
                        kfl_eccty(imate) = EXMSLD_EMEC_HUNTER_ORIGINAL
                        call set_default_HUNTER_ORIGINAL(imate)
                    else
                        kfl_eccty(imate) = EXMSLD_EMEC_HUNTER
                        call set_default_HUNTER(imate)
                    endif
  
                else if( words(2) == 'LAND ')then

                    kfl_exmsld_lands_ecc = .true.
                    if(  words(3) == 'BIDIR' )then
                        kfl_eccty(imate) = EXMSLD_EMEC_LAND_BIDIR
                        call set_default_LAND(imate)
                    else
                        kfl_eccty(imate) = EXMSLD_EMEC_LAND
                        call set_default_LAND(imate)
                    endif
                else if ( words(2) == 'NOMOD') then 
                    kfl_eccty(imate) = EXMSLD_EMEC_NO_EMEC
                else
                    call runend('eccou_read_data : NO ACTIVE MODEL DEFINED IN KER.DAT')
                endif
  
            else if( words(1) == 'CONTR' )then
                ! Control coupling factor
                if(      kfl_eccty(imate) == EXMSLD_EMEC_HUNTER )then
                    props_ecc(7,imate) = props_ecc(7,imate)   * getrea("CONTR",0.0_rp,"!Control coupling factor") 
                else if( kfl_eccty(imate) == EXMSLD_EMEC_HUNTER_ORIGINAL )then
                    props_ecc(15,imate) = props_ecc(15,imate) * getrea("CONTR",0.0_rp,"!Control coupling factor")
                else if( kfl_eccty(imate) == EXMSLD_EMEC_LAND .or. &
                         kfl_eccty(imate) == EXMSLD_EMEC_LAND_BIDIR )then
                    props_ecc(17,imate) = props_ecc(17,imate) * getrea("CONTR",0.0_rp,"!Control coupling factor")
                else
                    call runend('ERROR : ECCOU_READ_DATA : CONTR defined before MODEL')
                endif
          
            else if( words(1) == 'CAT50' )then
                ! Sensitivity reduction factor of Cat50 for Land model
                if(      kfl_eccty(imate) == EXMSLD_EMEC_LAND .or. & 
                         kfl_eccty(imate) == EXMSLD_EMEC_LAND_BIDIR )then
                     props_ecc(20,imate) = props_ecc(20,imate) * getrea("CAT50",1.0_rp,"!Sensitivity reduction factor of Cat50 for Land model") 
                else
                    call runend('ERROR : ECCOU_READ_DATA : CAT50 defined only for LAND MODEL')
                endif

  
            else if( words(1) == 'TIMEC' )then
                ! Time constant for the [Ca] fct
                timec_ecc = getrea("TIMEC",0.0_rp,"!Time constant for the [Ca] fct")
                
            else if( words(1) == 'HILLC' )then
                ! Hill coefficient
                if(      kfl_eccty(imate) == EXMSLD_EMEC_HUNTER  ) then
                    props_ecc( 5,imate) = getrea("HILLC",0.0_rp,"!Hill coefficient") 
                else if (kfl_eccty(imate) == EXMSLD_EMEC_HUNTER_ORIGINAL )then
                    !TODO: no idea which property of hunter original is hill coefficient
                else if( kfl_eccty(imate) == EXMSLD_EMEC_LAND .or. &
                         kfl_eccty(imate) == EXMSLD_EMEC_LAND_BIDIR )then
                    call runend('ERROR : ECCOU_READ_DATA : HILLC defined for Land model')
                else
                    call runend('ERROR : ECCOU_READ_DATA : CAL50 defined before MODEL')
                endif
  
            else if( words(1) == 'CAL50' )then
                ! Ca50 coefficient
                if( words(2) == 'AUTO ' ) then
                    kfl_force_ca50_calculation(imate) = .TRUE.
                else if( words(2) == 'MATER' ) then
                    kfl_ca50_provided(imate) = -getint('MATER=',0_ip,'#Current material')    

                    if ( kfl_ca50_provided(imate)==0_ip ) then 
                        call runend("ECCOUPLING MATERIAL "//trim(intost(imate))//": ERROR READING CA50 MATERIAL NUMBER")
                    end if
                else
                    kfl_ca50_provided(imate) = 1_ip
                    cell_ca50_permate(imate) = getrea("CAL50",0.0_rp,"!Ca50 coefficient, in solidz units (i.e. exmedi * 1000)")                    
                end if
            else if( words(1) == 'SCALE' )then
                ! Scale Ca50 ! TODO: Modify the name
                cell_caX_ecc(imate) = getrea("SCALE",0.0_rp,"!Scale Ca50")
  
            else if( words(1) == 'ORTK1' )then
                ! First coefficient for the orthotropic activation
                ortho_ecc( 2,imate) = getrea("ORTK1",0.0_rp,"!First coefficient for the orthotropic activation")
  
            else if( words(1) == 'ORTK2' )then
                ! Second coefficient for the orthotropic activation
                ortho_ecc( 3,imate) = getrea("ORTK2",0.0_rp,"!Second coefficient for the orthotropic activation")
            
            endif

            call ecoute('eccou_read_data')
        end do &
        coupling_sld_exm


        contains

            subroutine set_default_HUNTER(imate)
                implicit none
                integer(ip), intent(in) :: imate
                props_ecc( 1,imate) = 1.00_rp    ! Cam 
                props_ecc( 2,imate) = 0.50_rp    ! Ca50
                props_ecc( 3,imate) = 0.06_rp    ! Time
                props_ecc( 4,imate) = 1.45_rp    ! Eta
                props_ecc( 5,imate) = 3.00_rp    ! Hillc
                props_ecc( 6,imate) = 0.30_rp    ! k
                props_ecc( 7,imate) = 1.0e6_rp   ! Traction
                props_ecc( 8 ,imate) = 0.0_rp    ! Ca0 -- set later to the correct value
                props_ecc( 9 ,imate) = 1.95_rp   ! Beta_1
                props_ecc( 10,imate) = 0.31_rp   ! Beta_2
            end subroutine set_default_HUNTER

            subroutine set_default_HUNTER_ORIGINAL(imate)
                implicit none
                integer(ip), intent(in) :: imate
                props_ecc( 1,imate) = 100.0_rp   ! rho_0   
                props_ecc( 2,imate) = 163.0_rp   ! rho_1
                props_ecc( 3,imate) = 0.0_rp     ! Ca_0
                props_ecc( 4,imate) = 1.08_rp    ! Ca_max  
                props_ecc( 5,imate) = 1.26_rp    ! Ca_bmax in paper was 2.26 but goes to shit
                props_ecc( 6,imate) = 0.06_rp    ! tau_Ca  
                props_ecc( 8,imate) = 2.7_rp     ! gamma   
                props_ecc( 9,imate) = 6.9_rp     ! n_ref   
                props_ecc(10,imate) = 26.0_rp    ! alpha_0 
                props_ecc(11,imate) = 1.45_rp    ! beta_0  
                props_ecc(12,imate) = 1.95_rp    ! beta_1  
                props_ecc(13,imate) = 0.31_rp    ! beta_2  
                props_ecc(14,imate) = 0.5        ! 6.301_rp   ! pC50_ref, this corresponds to Ca50 of 0.5
                props_ecc(15,imate) = 100.0_rp   ! T_ref   
                props_ecc(16,imate) = 0.5_rp     ! a       
                props_ecc(17,imate) = 50.0_rp    ! A_1     
                props_ecc(18,imate) = 175.0_rp   ! A_2     
                props_ecc(19,imate) = 175.0_rp   ! A_3     
                props_ecc(20,imate) = 33.0_rp    ! alpha_1 
                props_ecc(21,imate) = 2850.0_rp  ! alpha_2 
                props_ecc(22,imate) = 2850.0_rp  ! alpha_3 
            end subroutine set_default_HUNTER_ORIGINAL

            subroutine set_default_LAND(imate)
                implicit none
                integer(ip), intent(in) :: imate
                props_ecc( 1,imate) = 100.0_rp  ! k_TRPN    
                props_ecc( 2,imate) = 2.0_rp    ! n_TRPN   
                props_ecc( 3,imate) = 0.805_rp  ! Ca50_ref  
                props_ecc( 4,imate) = 1000.0_rp ! k_u      
                props_ecc( 5,imate) = 5.0_rp    ! n_tm     
                props_ecc( 6,imate) = 0.35_rp   ! TRPN_50     
                props_ecc( 7,imate) = 182.0_rp  ! k_uw         
                props_ecc( 8,imate) = 12.0_rp   ! k_ws       
                props_ecc( 9,imate) = 0.5_rp    ! r_w        
                props_ecc(10,imate) = 0.25_rp   ! r_s         
                props_ecc(11,imate) = 0.0085_rp ! y_s          
                props_ecc(12,imate) = 0.615_rp  ! y_w          
                props_ecc(13,imate) = 2.23_rp   ! phi         
                props_ecc(14,imate) = 25.0_rp   ! A_eff       
                props_ecc(15,imate) = 2.3_rp    ! beta_0     
                props_ecc(16,imate) = -2.4_rp   ! beta_1      
                props_ecc(17,imate) = 120.0_rp  ! T_ref             
                props_ecc(18,imate) = 0.0_rp    ! Ca0, set later on
                props_ecc(19,imate) = 0.31_rp   ! beta_2
                props_ecc(20,imate) = 1.0_rp    ! Ca50 reducction factor
            end subroutine set_default_LAND

    end subroutine eccou_read_data

    ! -------------------------------------------------------------------------
    !> @author  Adria Quintanas-Corominas (adria.quintanas@bsc.es)
    !> @author  Alfonso Santiago (alfonso.santiago@bsc.es)
    !> @Details Interchange the flags and properties readed in ker.dat.
    !           called in ker_parall.f90
    ! -------------------------------------------------------------------------
    subroutine eccou_send_data(itask)
        use def_master
        use mod_exchange,          only : exchange_add
        use mod_exchange,          only : exchange_end
        use mod_exchange,          only : exchange_init

        ! -------------------------------
        implicit none
        ! -------------------------------
        integer(ip)                    :: itask
        ! -------------------------------
        if( ISEQUEN ) return

        select case(itask)
        case(0_ip)
          call exchange_init()
          call exchange_add(kfl_exmsld_ecc)
          call exchange_add(kfl_exmsld_moncod_ecc)
          call exchange_add(kfl_exmsld_mulcod_ecc)
          call exchange_add(kfl_exmsld_lands_ecc)
          call exchange_add(coupling_delay_ecc)
          call exchange_add(ncelltypes_ecc)
          call exchange_end()

          ! Only continue in case of coupling
          if( .not. kfl_exmsld_ecc ) return

          if( ISLAVE )then
              call eccou_allocate_memory(0_ip)
          endif
  
          call exchange_init()
          call exchange_add(kfl_cellmod)
          call exchange_add(kfl_eccty)
          call exchange_add(props_ecc)
          call exchange_add(cell_caX_ecc)
          call exchange_add(ortho_ecc)
          call exchange_add(kfl_force_ca50_calculation)
          call exchange_add(kfl_ca50_provided)
          call exchange_add(cell_ca50_permate)
          call exchange_end()
        case(1_ip) !after reading ker.dat
          if (kfl_exmsld_ecc) then
              call exchange_init()
              call exchange_add(cell_ca0_ecc)
              call exchange_add(cell_ca50_ecc)
              call exchange_end()
          end if
        case default
          call runend("eccou_send_data: Undefined task requested")
        end select

    end subroutine eccou_send_data

    subroutine eccou_send_cellmod()
        use def_master
        use mod_exchange,          only : exchange_add
        use mod_exchange,          only : exchange_end
        use mod_exchange,          only : exchange_init

        call exchange_init()
        call exchange_add(kfl_cellmod)
        call exchange_end()
        
    end subroutine

    ! -------------------------------------------------------------------------
    !> @author  Adria Quintanas-Corominas (adria.quintanas@bsc.es)
    !> @author  Alfonso Santiago (alfonso.santiago@bsc.es)
    ! -------------------------------------------------------------------------
    subroutine eccou_do_plugin( &
        icoup )
        use def_domain,            only : npoin
        use def_master,            only : modul, mem_modul, ID_EXMEDI, ID_SOLIDZ
        use mod_couplings,         only : COU_INTERPOLATE_NODAL_VALUES
        use mod_couplings,         only : COU_SET_FIXITY_ON_TARGET
        use mod_memory,            only : memory_deallo, memory_alloca
        use mod_communications,    only : par_barrier, par_broadcast, par_sum

        ! -------------------------------
        implicit none
        ! -------------------------------
        integer(ip), intent(in)         :: icoup
        logical(lg), save               :: first_time = .true.
        real(rp), pointer, dimension(:) :: aux
        ! -------------------------------
        if( first_time )then
            call get_master_ranks()
            call exchange_initilisation()
            call exchange_physics()
            nullify(aux)
        endif
        call memory_alloca(mem_modul(1:2,modul),'AUX','eccou_do_plugin',aux,npoin)
        if( modul == ID_SOLIDZ )then
            ! Recieve 
            call COU_INTERPOLATE_NODAL_VALUES(icoup,1_ip,calcium_ecc,aux)
        elseif( modul == ID_EXMEDI )then
            ! Send
            call COU_INTERPOLATE_NODAL_VALUES(icoup,1_ip,aux,calcium_ecc)
        endif
        call memory_deallo(mem_modul(1:2,modul),'AUX','eccou_do_plugin',aux)
        ! -------------------------------
    end subroutine eccou_do_plugin

    ! -------------------------------------------------------------------------
    !> @author  Adria Quintanas-Corominas (adria.quintanas@bsc.es)
    !> @author  Alfonso Santiago (alfonso.santiago@bsc.es)
    ! -------------------------------------------------------------------------
    subroutine eccou_do_volume_integral( &
        ipoin, itask, integral)
        use def_domain,            only : vmass, npoin
        use mod_communications,    only :  PAR_SUM,PAR_BARRIER
        use def_master,            only : ID_EXMEDI, ID_SOLIDZ, INOTMASTER
        use def_master,            only : npoi1, npoi2, npoi3

        implicit none
        integer(ip), intent(in)                :: ipoin
        character(len=4), intent(in)           :: itask
        integer(ip)                            :: i
        real(rp)                               :: integral


        select case(itask)

        case( "doal" )
            ! Take profit of the loop to integrate the calcium
            integral=-1.0_rp
            volcai_ecc=0.0_rp
            IF(INOTMASTER) then
                do i=1,npoin
                    if ((i<=npoi1) .OR. ((i>=npoi2) .AND. (i<=npoi3))) then
                        volcai_ecc = volcai_ecc + vmass(i) * calcium_ecc(i)
                    end if
                enddo
            ENDIF
            call PAR_SUM(volcai_ecc, "IN MY CODE")
            integral=volcai_ecc
        case( "node" )
            ! Take profit of the loop to integrate the calcium
            volcai_ecc = volcai_ecc + vmass(ipoin) * calcium_ecc(ipoin) 
        case( "doma" )
            ! Use the synch point to add volcai
            call PAR_SUM(volcai_ecc, "IN MY CODE")
        end select
    end subroutine eccou_do_volume_integral

    ! -------------------------------------------------------------------------
    !> @author  Adria Quintanas-Corominas (adria.quintanas@bsc.es)
    !> @author  Alfonso Santiago (alfonso.santiago@bsc.es)
    ! -------------------------------------------------------------------------
    subroutine eccou_do_lagrangian( &
        icoup )
        use mod_communications,    only :  PAR_SUM,PAR_BARRIER
        use mod_couplings,         only : COU_INTERPOLATE_NODAL_VALUES
        use def_domain,            only : npoin, ndime
        use def_master,            only : modul, mem_modul, ID_EXMEDI, ID_SOLIDZ, INOTMASTER
        use mod_memory,            only : memory_deallo, memory_alloca
        !use mod_couplings,         only : COU_SET_FIXITY_ON_TARGET
        !use mod_communications,    only : par_barrier, par_broadcast
        ! -------------------------------
        implicit none
        ! -------------------------------
        integer(ip), intent(in)         :: icoup
        logical(lg), save               :: first_time = .true.
        real(rp), pointer               :: aux(:,:)
        real(rp), pointer               :: xvalu(:,:)
        ! -------------------------------
        if( first_time )then
            nullify(aux)
            nullify(xvalu)
        endif
        !
        if(INOTMASTER)then
            call memory_alloca(mem_modul(1:2,modul),'AUX','eccou_do_plugin',aux,ndime,npoin)
            call memory_alloca(mem_modul(1:2,modul),'XVALU','eccou_do_plugin',xvalu,ndime,npoin)
        else
            call memory_alloca(mem_modul(1:2,modul),'AUX','eccou_do_plugin',aux,1_ip,1_ip)
            call memory_alloca(mem_modul(1:2,modul),'XVALU','eccou_do_plugin',xvalu,1_ip,1_ip)
        endif
        ! Recieve 
        call COU_INTERPOLATE_NODAL_VALUES(icoup,ndime,xvalu,aux)


       IF(INOTMASTER) displ_lagr_ecc(:,:,1)=xvalu(:,:)


        call memory_deallo(mem_modul(1:2,modul),'AUX','eccou_do_plugin',aux)
        call memory_deallo(mem_modul(1:2,modul),'XVALU','eccou_do_plugin',xvalu)

    end subroutine eccou_do_lagrangian


    subroutine get_master_ranks()
        use def_master,            only :  modul, IMASTER, ID_EXMEDI, ID_SOLIDZ
        use mod_parall,            only :  PAR_MY_UNIVERSE_RANK
        use mod_communications,    only :  PAR_SUM
        ! -------------------------------
        implicit none
        ! -------------------------------
        real(rp), pointer              :: ranks(:)
        ! -------------------------------
        allocate(ranks(2))
        ranks(:) = 0_ip
        if( modul == ID_SOLIDZ )then
            if( IMASTER )then
               ranks(1) = PAR_MY_UNIVERSE_RANK
            endif
        elseif( modul == ID_EXMEDI )then
            if( IMASTER )then
               ranks(2) = PAR_MY_UNIVERSE_RANK
            endif
        endif
        call PAR_SUM(ranks,'IN THE UNIVERSE')
        MPI_RANK_MASTER_SLD = int(ranks(1))
        MPI_RANK_MASTER_EXM = int(ranks(2))
        deallocate(ranks)
        ! -------------------------------
    end subroutine get_master_ranks


    subroutine exchange_initilisation()
        use mod_communications,    only :  PAR_BROADCAST
        ! -------------------------------
        implicit none
        real(rp), pointer              ::  p1(:), p2(:) 
        ! -------------------------------
        allocate(p1(nmate)); p1 = kfl_eccty(1:nmate)
        allocate(p2(nmate)); p2 = kfl_cellmod(1:nmate)
        call PAR_BROADCAST(nmate,p1,wherein='IN THE UNIVERSE',root_rank=MPI_RANK_MASTER_SLD)
        call PAR_BROADCAST(nmate,p2,wherein='IN THE UNIVERSE',root_rank=MPI_RANK_MASTER_EXM)
        kfl_eccty(1:nmate) = p1; deallocate(p1)
        kfl_cellmod(1:nmate) = p2; deallocate(p2)
        ! -------------------------------
    end subroutine exchange_initilisation


    subroutine exchange_physics()
        use mod_communications,    only :  PAR_BROADCAST
        ! -------------------------------
        implicit none
        ! -------------------------------
        call PAR_BROADCAST(ncelltypes_ecc,nmate,cell_ca0_ecc,wherein='IN THE UNIVERSE',root_rank=MPI_RANK_MASTER_EXM)
        call PAR_BROADCAST(ncelltypes_ecc,nmate,cell_ca50_ecc,wherein='IN THE UNIVERSE',root_rank=MPI_RANK_MASTER_EXM)
        ! -------------------------------
    end subroutine exchange_physics


    ! -------------------------------------------------------------------------
    !> @author  Adria Quintanas-Corominas (adria.quintanas@bsc.es)
    !> @author  Alfonso Santiago (alfonso.santiago@bsc.es)
    ! -------------------------------------------------------------------------
    subroutine eccou_manage_variables( &
        itask )
        use def_master,            only : INOTEMPTY, ITASK_ENDSTE
        ! -------------------------------
        implicit none
        ! -------------------------------
        integer(ip), intent(in)        :: itask
        ! -------------------------------
        if( INOTEMPTY )then
            select case(itask)

                case( ITASK_ENDSTE )
                    ! ITER_AUX = 2; ITER_K = 1
                    state_ecc(:,:,1) = state_ecc(:,:,2) 
                    strch_ecc(:,:,1) = strch_ecc(:,:,2)

            end select
        endif
        ! -------------------------------
    end subroutine eccou_manage_variables


    subroutine exm_sld_ecc_initialize_troponin( &
        )
        ! -------------------------------
        implicit none
        ! -------------------------------
        integer(ip)                    :: ii,jj
        ! -------------------------------
        if(associated(troponin_ecc))then
            do jj = 1, size(troponin_ecc,DIM=2_ip,KIND=ip)
                do ii = 1, size(troponin_ecc,DIM=1_ip,KIND=ip)
                   troponin_ecc(ii,jj) = 0.0_rp
                enddo
            enddo
        endif
        ! -------------------------------
    end subroutine exm_sld_ecc_initialize_troponin


    subroutine exm_sld_ecc_interchange_troponin( &
        )
        use def_domain,            only : vmass
        ! -------------------------------
        implicit none
        ! -------------------------------
        call rhsmod(1_ip,troponin_ecc(:,1))
        troponin_ecc(:,1) = troponin_ecc(:,1) / vmass(:)
        ! -------------------------------
    end subroutine exm_sld_ecc_interchange_troponin


    subroutine eccou_assemble_troponin( &
        ielem, pnode, pgaus,  gp_N, gp_vol )
        use def_domain,            only : lnods
        use mod_matrix,            only : matrix_assemble_element_RHS
        ! -------------------------------
        implicit none
        ! -------------------------------
        integer(ip), intent(in)        :: ielem
        integer(ip), intent(in)        :: pnode
        integer(ip), intent(in)        :: pgaus
        real(rp),    intent(in)        :: gp_vol(pgaus)
        real(rp),    intent(in)        :: gp_N(pnode,pgaus)
        ! -------------------------------
        integer(ip)                    :: inode, igaus, iposi
        real(rp)                       :: gp_trop(pgaus)
        real(rp)                       :: el_trop(pnode)
        ! -------------------------------
        ! Recover current value of the troponning the gauss point.
        ! It is computed into get_active_traction_LANDS, which is 
        ! called by solidz during the elemental loop assembly. 
        do igaus = 1, pgaus
            iposi = (igaus - 1_ip)*nsdvs_ecc + 3_ip 
            gp_trop(igaus) = state_ecc(iposi,ielem,2)        
        enddo

        ! Integrate troponin
        el_trop(:) = 0.0_rp
        do inode = 1, pnode
            do igaus = 1, pgaus
                el_trop(inode) = el_trop(inode) + gp_N(inode,igaus) * gp_vol(igaus) * gp_trop(igaus)
            end do
        end do

        ! Assemble troponin
        call matrix_assemble_element_RHS(1_ip,1_ip,pnode,lnods(:,ielem),el_trop(:),troponin_ecc(:,1))
        ! -------------------------------
    end subroutine eccou_assemble_troponin


    subroutine eccou_interpolate_variables_at_gp( &
        ielem )
        ! -------------------------------
        use def_domain,            only : lnnod, ltype, lnods, elmar, nodemat, ngaus
        use def_kermod,            only :  kfl_celltype_fun

        ! -------------------------------
        implicit none
        ! -------------------------------
        integer(ip), intent(in)        :: ielem
        ! -------------------------------
        integer(ip)                    :: pnode, pelty, pgaus, ipoin, igaus, inode, imate, icell, iprop
        integer(ip), pointer           :: lpoin(:)
        real(rp),    pointer           :: gp_N(:,:)
        ! -------------------------------
        ! Recover elemental information
        pnode = lnnod(ielem)
        pelty = ltype(ielem)
        pgaus = ngaus(pelty)
        lpoin => lnods(:,ielem)
        gp_N  => elmar(pelty)%shape

        ! Interpolate calcium and properties from nodes to gp
        gp_Ca_ecc = 0.0_rp
        gp_Ca0_ecc = 0.0_rp
        gp_Ca50_ecc = 0.0_rp
        gp_props_ecc = 0.0_rp
        do inode = 1,pnode
            ipoin = lpoin(inode)
            imate = nodemat(ipoin)

            !interpolate calcium only if it's calculated
            if ( kfl_cellmod(imate) /= EXMSLD_CELL_NOMODEL ) then
                if ( kfl_celltype_fun < 0_ip ) then
                    icell = int(celty_ecc(ipoin))
                else
                    icell = 1_ip
                end if

                ! calcium_ecc(:) is pointing the vconc(1,:,1)
                do igaus = 1, pgaus
                    gp_Ca_ecc(igaus) = gp_Ca_ecc(igaus) + gp_N(inode,igaus) * calcium_ecc(ipoin) 
                    gp_Ca0_ecc(igaus) = gp_Ca0_ecc(igaus) + gp_N(inode,igaus) * cell_ca0_ecc(icell,imate) 
                    ! in the following *1000 is not needed, for some reason
                    gp_Ca50_ecc(igaus) = gp_Ca50_ecc(igaus) + gp_N(inode,igaus) * cell_ca50_ecc(icell,imate) * cell_caX_ecc(imate)
                    do iprop = 1, nprop_ecc
                        gp_props_ecc(iprop,igaus) = gp_props_ecc(iprop,igaus) + props_ecc(iprop,imate) * gp_N(inode,igaus) 
                    enddo
                enddo          
            end if

        end do
        ! -------------------------------
    end subroutine eccou_interpolate_variables_at_gp


    ! -------------------------------------------------------------------------
    !> @author  Adria Quintanas-Corominas (adria.quintanas@bsc.es)
    !> @author  Alfonso Santiago (alfonso.santiago@bsc.es)
    !> @details Ref_A : Levrero-Florencio et al. 2019 (DOI: 10.1016/j.cma.2019.112762)
    ! -------------------------------------------------------------------------
    subroutine eccou_add_active_stress_and_moduli( &
        itask, ielem, imate, pnode, pgaus, gp_N, gp_F, gp_str, gp_pio, gp_dds, gp_tmo )
        ! -------------------------------
        implicit none
        ! -------------------------------
        integer(ip), intent(in)        :: itask
        integer(ip), intent(in)        :: ielem                   ! Number of element
        integer(ip), intent(in)        :: imate                   ! Number of material
        integer(ip), intent(in)        :: pnode                   ! Number of nodes
        integer(ip), intent(in)        :: pgaus                   ! Number of gauss points
        real(rp),    intent(in)        :: gp_N(pnode,pgaus)       ! Shape functions
        real(rp),    intent(in)        :: gp_F(3,3,pgaus)         ! Deformation Gradient at gauss points
        real(rp),    intent(inout)     :: gp_str(3,3,pgaus)       ! Second Piola-Kirchhoff at gauss points
        real(rp),    intent(inout)     :: gp_pio(3,3,pgaus)       ! First Piola-Kirchhoff at gauss points
        real(rp),    intent(inout)     :: gp_dds(3,3,3,3,pgaus)   ! Second elasticity tangent moduli at gauss points
        real(rp),    intent(inout)     :: gp_tmo(3,3,3,3,pgaus)   ! First elasticity tangent moduli at gauss points
        ! -------------------------------
        integer(ip)                    :: igaus                   ! Loop counters
        integer(ip)                    :: i, j, k, l, m, n, ii, jj 
        real(rp)                       :: i4f, i4s, i4n           ! Kinematics invariants
        real(rp)                       :: hf, hs, hn              ! Stretches at gauss points
        real(rp)                       :: statvar(10,2)           ! Local state variables
        real(rp)                       :: gp_T, gp_T_per, gp_dTdh ! Active tractions at the guass points
        real(rp)                       :: gp_C(3,3)               ! Right-Cauchy strain tensor
        real(rp)                       :: gp_S(3,3,pgaus)         ! Second Piola-Kirchhoff at gauss points
        real(rp)                       :: gp_PK1(3,3,pgaus)       ! First Piola-Kirchhoff at gauss points
        real(rp)                       :: gp_dSdE(3,3,3,3,pgaus)  ! Second elasticity tangent moduli at gauss points
        real(rp)                       :: gp_dPdF(3,3,3,3,pgaus)  ! First elasticity tangent moduli at gauss points
        real(rp)                       :: gp_fbr_0(3,pgaus)       ! Fiber direction in the reference config. at gauss points
        real(rp)                       :: gp_sht_0(3,pgaus)       ! Sheat direction in the reference config. at gauss points
        real(rp)                       :: gp_nrm_0(3,pgaus)       ! Norma direction in the reference config. at gauss points
        real(rp)                       :: fxfxfxf(3,3,3,3)        ! Auxiliar tensors 
        real(rp)                       :: sxsxsxs(3,3,3,3)
        real(rp)                       :: nxnxnxn(3,3,3,3)
        real(rp)                       :: sxsxfxf(3,3,3,3)
        real(rp)                       :: nxnxfxf(3,3,3,3)


        gp_dTdh = 0.0_rp

        ! -------------------------------
        ! Recover reference fiber at GP
        call biofibers % get_reference_fibers_at_gp( &
            ielem, pnode, pgaus, gp_N, gp_fbr_0, gp_sht_0, gp_nrm_0 ) 

        ! Select law according to the type of coupling
        select case( kfl_eccty(imate) )

            case( EXMSLD_EMEC_HUNTER )
                gp_props_ecc(2,1:pgaus) = gp_Ca50_ecc(1:pgaus)
                gp_props_ecc(8,1:pgaus) = gp_Ca0_ecc(1:pgaus)
                get_active_traction_generic => get_active_traction_HUNTER
            
            case( EXMSLD_EMEC_HUNTER_ORIGINAL )
                gp_props_ecc(3,1:pgaus) = gp_Ca0_ecc(1:pgaus)
                gp_props_ecc(14,1:pgaus) = gp_Ca50_ecc(1:pgaus) 
                get_active_traction_generic => get_active_traction_HUNTER_ORIGINAL 

            case( EXMSLD_EMEC_LAND, EXMSLD_EMEC_LAND_BIDIR )
                gp_props_ecc(3,1:pgaus) = gp_Ca50_ecc(1:pgaus)
                gp_props_ecc(18,1:pgaus) = gp_Ca0_ecc(1:pgaus)
                get_active_traction_generic => get_active_traction_LANDS
            
        end select

        ! Second elasticity tensors 
        loop_gp_S_dSdE: &
        do igaus = 1, pgaus

            ! Right Cauchy tensor
            gp_C(:,:) = TxM(gp_F(:,:,igaus),gp_F(:,:,igaus))

            ! Compute invatiants
            ! - i4_f = sqrt( f_0 · C f_0 )
            i4f = ROUNDOFF(VxMxV(gp_fbr_0(:,igaus),gp_C(:,:),gp_fbr_0(:,igaus)))
            i4s = ROUNDOFF(VxMxV(gp_sht_0(:,igaus),gp_C(:,:),gp_sht_0(:,igaus)))
            i4n = ROUNDOFF(VxMxV(gp_nrm_0(:,igaus),gp_C(:,:),gp_nrm_0(:,igaus)))
            ! - consistency condition ( if s0 and n0 are not defined can produce numerical
            !   errors. )
            if( i4s < 1.0e-16_rp ) i4s = 1.0_rp
            if( i4n < 1.0e-16_rp ) i4n = 1.0_rp

            ! Compute stretches
            ! - h_f = sqrt( i4_f )
            hf = sqrt(i4f)
            hs = sqrt(i4s)
            hn = sqrt(i4n)

            ! Auxiliary indices
            ii = (igaus - 1_ip)*nsdvs_ecc + 1_ip
            jj = ii + (nsdvs_ecc - 1_ip)
            ! Get state variables
            statvar(1:nsdvs_ecc,1)   = state_ecc(ii:jj,ielem,1) 
            statvar(  nsdvs_ecc+1,1) = strch_ecc(igaus,ielem,1)

            ! Compute active traction according to the electro-mechanical model
            call get_active_traction_generic( gp_props_ecc(:,igaus), hf + 1.0e-6_rp, gp_Ca_ecc(igaus), gp_T_per, gp_dTdh, statvar(:,:) )
            call get_active_traction_generic( gp_props_ecc(:,igaus), hf, gp_Ca_ecc(igaus), gp_T, gp_dTdh, statvar(:,:) )
            gp_dTdh = (gp_T_per - gp_T) / 1.0e-6_rp

            ! Set state variables
            state_ecc(ii:jj,ielem,2) = statvar(1:nsdvs_ecc,2)

            ! Set the stretch along the fiber direction
            strch_ecc(igaus,ielem,2) = hf

            ! Compute second Piola-Kirchoff stress tensor
            ! - Ref_A (Eq.2.11)
            gp_S(:,:,igaus) = 0.0_rp
            gp_S(:,:,igaus) = gp_T * ( VxV(gp_fbr_0(:,igaus),gp_fbr_0(:,igaus)) / i4f +  &
                                       VxV(gp_sht_0(:,igaus),gp_sht_0(:,igaus)) / i4s * ortho_ecc(1,imate) +  & ! Kort1
                                       VxV(gp_nrm_0(:,igaus),gp_nrm_0(:,igaus)) / i4n * ortho_ecc(2,imate) )    ! Kort2

            ! Compute second elasticity tensor
            ! - Auxiliar terms
            fxfxfxf = VxVxVxV(gp_fbr_0(:,igaus),gp_fbr_0(:,igaus),gp_fbr_0(:,igaus),gp_fbr_0(:,igaus))
            sxsxsxs = VxVxVxV(gp_sht_0(:,igaus),gp_sht_0(:,igaus),gp_sht_0(:,igaus),gp_sht_0(:,igaus))
            nxnxnxn = VxVxVxV(gp_nrm_0(:,igaus),gp_nrm_0(:,igaus),gp_nrm_0(:,igaus),gp_nrm_0(:,igaus))
            sxsxfxf = VxVxVxV(gp_sht_0(:,igaus),gp_sht_0(:,igaus),gp_fbr_0(:,igaus),gp_fbr_0(:,igaus))
            nxnxfxf = VxVxVxV(gp_nrm_0(:,igaus),gp_nrm_0(:,igaus),gp_fbr_0(:,igaus),gp_fbr_0(:,igaus))
            ! - Ref_A (Appendix F (Warning, possible mistake in the first summand in the reference))
            gp_dSdE(:,:,:,:,igaus) = 0.0_rp
            gp_dSdE(:,:,:,:,igaus) = ( gp_dTdh / hf**3 - 2.0_rp * gp_T / hf**4 ) * fxfxfxf + &
                                     ( gp_dTdh / (hf * hs**2) * sxsxfxf - 2.0_rp * gp_T * sxsxsxs / hs**4 ) * ortho_ecc(1,imate) + &
                                     ( gp_dTdh / (hf * hn**2) * nxnxfxf - 2.0_rp * gp_T * nxnxnxn / hn**4 ) * ortho_ecc(2,imate)

        end do &
        loop_gp_S_dSdE


        ! Frist elasticity tensors
        loop_gp_P_dPdF: &
        do igaus = 1, pgaus

            ! Compute first Piola-Kirchoff (P) stress tensor from the second one (S)
            ! - Ref_B (Eq.)
            gp_PK1(:,:,igaus) = 0.0_rp
            do k = 1, 3
                do j = 1, 3
                    do i = 1, 3
                        gp_PK1(j,k,igaus) = gp_PK1(j,k,igaus) + gp_F(j,i,igaus) * gp_S(i,k,igaus)
                    enddo
                enddo
            enddo

            ! Compute first elasticity tensor (dPdF) from the second one (dSdE = 2 dSdC)
            ! - Ref_B (Eq.)
            gp_dPdF(:,:,:,:,igaus) = 0.0_rp
            do i = 1, 3
                do j = 1, 3
                    do k = 1, 3 
                        do l = 1, 3
                            do m = 1, 3
                                do n = 1, 3
                                    gp_dPdF(i,j,k,l,igaus)= gp_dPdF(i,j,k,l,igaus) + gp_dSdE(m,j,n,l,igaus) * gp_F(i,m,igaus) * gp_F(k,n,igaus)
                                enddo
                            enddo
                            ! Equivalent to  gp_dPdF(i,j,k,l,igaus) = gp_dPdF(i,j,k,l,igaus) + I_3x3(i,k) * gp_S(j,l,igaus)
                            if( i == k ) gp_dPdF(i,j,k,l,igaus) = gp_dPdF(i,j,k,l,igaus) + gp_S(j,l,igaus)
                        enddo
                    enddo
                enddo
            enddo

        enddo &
        loop_gp_P_dPdF

        ! Add contributions
        gp_str = gp_str + gp_S
        gp_pio = gp_pio + gp_PK1
        if( itask  == 2 )then
            gp_dds = gp_dds + gp_dSdE
            gp_tmo = gp_tmo + gp_dPdF
        endif

        contains

            pure function TxM(A,B) result(R)
                implicit none
                real(rp), dimension(3,3), intent(in) :: A
                real(rp), dimension(3,3), intent(in) :: B
                real(rp), dimension(3,3)             :: R
                integer(ip)                          :: i, j, k
                R = 0.0_rp
                do k = 1, 3
                    do j = 1, 3
                        do i = 1, 3
                            R(i,j) = R(i,j) + A(k,i)*B(k,j)
                        enddo
                    enddo
                enddo
            end function TxM

            pure function VxMxV(a,M,b) result(r)
                implicit none
                real(rp), dimension(3),   intent(in) :: a
                real(rp), dimension(3,3), intent(in) :: M
                real(rp), dimension(3),   intent(in) :: b
                real(rp)                             :: r
                integer(ip)                          :: i, j
                r = 0.0_rp
                do j = 1, 3
                    do i = 1, 3
                        r = r + a(i)*M(i,j)*b(j)
                    enddo
                enddo
            end function VxMxV

            pure function ROUNDOFF(val) result(res)
                implicit none
                real(rp), intent(in)           :: val
                real(rp)                       :: res
                real(rp), parameter            :: thr = 1.0e-12_rp
                if( abs((val - 1.0_rp)) < thr )then
                    res = 1.0_rp
                else
                    res = val
                endif
            end function ROUNDOFF

            pure function VxV(a,b) result(R)
                implicit none
                real(rp), dimension(3),  intent(in) :: a
                real(rp), dimension(3),  intent(in) :: b
                real(rp), dimension(3,3)            :: R
                integer(ip)                         :: i, j
                do i = 1, 3
                    do j = 1, 3
                        R(i,j) = a(i)*b(j)
                    enddo
                enddo
            end function VxV

            pure function VxVxVxV(a,b,c,d) result(R)
                implicit none
                real(rp), dimension(3), intent(in) :: a
                real(rp), dimension(3), intent(in) :: b
                real(rp), dimension(3), intent(in) :: c
                real(rp), dimension(3), intent(in) :: d
                real(rp), dimension(3,3,3,3)       :: R
                integer(ip)                        :: i, j, k, l
                do i = 1, 3
                    do j = 1, 3
                        do k = 1, 3
                            do l = 1, 3
                                R(i,j,k,l) = a(i)*b(j)*c(k)*d(l)
                            enddo
                        enddo
                    enddo
                enddo
            end function VxVxVxV

    end subroutine eccou_add_active_stress_and_moduli


    ! -------------------------------------------------------------------------
    !> @author  Adria Quintanas-Corominas (adria.quintanas@bsc.es)
    !> @author  Alfonso Santiago (alfonso.santiago@bsc.es)
    !> @author  Pierre Lafortune 
    !> @details Simplified Hunter model as per 
    !           Ref_B : P. Lafortune et al. 2012, DOI: 10.1002/cnm.1494
    ! -------------------------------------------------------------------------
    subroutine get_active_traction_HUNTER( &
        props, lamda, Ca, tract, dtract, sdv )
        ! -------------------------------
        implicit none
        ! -------------------------------
        real(rp),          intent(in)    :: props(:)  ! Properties
        real(rp),          intent(in)    :: lamda     ! Lambda (stretches)
        real(rp),          intent(in)    :: Ca        ! Calcium
        real(rp),          intent(inout) :: tract     ! Active traction
        real(rp),          intent(inout) :: dtract    ! Derivative active traction w.r.t lamda
        real(rp),          intent(inout) :: sdv(:,:)  ! State variables
        ! -------------------------------
        real(rp)                       :: Ca0, Cam, Ca50, Ca2p
        real(rp)                       :: thau, eta, n, k, n_ref, pC50_ref
        real(rp)                       :: tmax, beta_1, beta_2, pC50
        real(rp)                       :: lamda_auxi       
        ! -------------------------------
        ! Assign properties
        Cam  = props(1)
        pC50_ref = props(2)
        thau = props(3)
        eta  = props(4)
        n_ref    = props(5)
        k    = props(6)
        tmax = props(7)
        Ca0  = props(8)
        beta_1   = props(9)
        beta_2   = props(10)

        ! Compute Ca2p according to the cellular model
        Ca2p = max(Ca - Ca0,0.0_rp)
        lamda_auxi = min(lamda_thr,lamda)        ! Tropomyosin kinetics (z)

        n    = n_ref*(1.0_rp + beta_1*(lamda_auxi - 1.0_rp))
        pC50 = pC50_ref*(1.0_rp + beta_2*(lamda_auxi - 1.0_rp))
        Ca50 = 10.0_rp**(6.0_rp - pC50) !pC50_ref needs to be calculated from Ca50 of exmedi with this formula
        
        !else !old behaviour
        !    n = n_ref
        !    Ca50 = pC50_ref
        !end if

        ! Compute active traction
        ! { Ref_B : Eq. 21 }
        tract = (Ca2p**n)*tmax*(1.0_rp + eta*(lamda - 1.0_rp))/((Ca2p**n) + (Ca50**n))

        ! Compute derivative active traction w.r.t. lambda along the fiber direction
        dtract = (Ca2p**n)*tmax*eta/((Ca2p**n) + (Ca50)**n)
        ! -------------------------------
    end subroutine get_active_traction_HUNTER


    ! -------------------------------------------------------------------------
    !> @author  Jazmin Aguado
    !> @author  Constantine Butakoff
    !> @author  Adria Quintanas-Corominas 
    !> @details Ref_A : This is the original model as per P.J. Hunter  et al. 
    !               1998, DOI: 10.1016/S0079-6107(98)00013-3, not simplified
    ! -------------------------------------------------------------------------
    subroutine get_active_traction_HUNTER_ORIGINAL( &
        props, lmb, Ca, tract, dtract, sdv )
        use def_master,           only :  dtime
!        use def_master,           only :  cutim
        ! -------------------------------
        implicit none
        ! -------------------------------
        real(rp),          intent(in)    :: props(:)  ! Properties
        real(rp),          intent(in)    :: lmb       ! Lambda (stretches)
        real(rp),          intent(in)    :: Ca        ! Calcium
        real(rp),          intent(inout) :: tract     ! Active traction
        real(rp),          intent(inout) :: dtract    ! Derivative active traction w.r.t lamda
        real(rp),          intent(inout) :: sdv(:,:)  ! State variables
        ! -------------------------------
        integer(ip), parameter         :: old = 1_ip, new = 2_ip
        ! -------------------------------
        real(rp)                       :: n, Ca_0, dlmbdt
        real(rp)                       :: rho_0, rho_1, Ca_max, Ca_bmax, tau_Ca, n_ref, alpha_0, &
        &                                 beta_1, beta_2, pC50_ref, T_ref, a, A_1, A_2, A_3, alpha_1,  &
        &                                 alpha_2, alpha_3, lmb_0, Ca_b, T_0, T_a, Q_1, Q_2, Q_3, Q,   &
        &                                 Ca_i, dCaidt, pC50, C50, dzdt, z, dT0dt, gamma, beta_0
        real(rp)                       :: lamda_auxi       
        ! -------------------------------
        ! Assign properties
        rho_0    =  props( 1)
        rho_1    =  props( 2)
        Ca_0     =  props( 3)
        Ca_max   =  props( 4)
        Ca_bmax  =  props( 5)
        tau_Ca   =  props( 6)
        gamma    =  props( 8)
        n_ref    =  props( 9)
        alpha_0  =  props(10)
        beta_0   =  props(11)
        beta_1   =  props(12)
        beta_2   =  props(13)
        pC50_ref =  props(14)
        T_ref    =  props(15)
        a        =  props(16)
        A_1      =  props(17)
        A_2      =  props(18)
        A_3      =  props(19)
        alpha_1  =  props(20)
        alpha_2  =  props(21)
        alpha_3  =  props(22)
        
        ! Recover values of the previous time step
        lmb_0 = sdv(1,old)
        Ca_b  = sdv(2,old)
        T_0   = sdv(3,old)
        T_a   = sdv(4,old)
        z     = sdv(5,old)
        Q_1   = sdv(6,old)
        Q_2   = sdv(7,old)
        Q_3   = sdv(8,old)

        ! Compute rate of lambda
        lamda_auxi = min(lamda_thr,lmb)        ! Tropomyosin kinetics (z)
        dlmbdt = (lmb - lmb_0)/dtime

        ! TnC - Ca++ binding kinetics
        ! Time-dependent intracellular concentration of Ca^{2+} {Eq.1}
        Ca_i   = max(Ca - Ca_0, 0.0_rp)!Ca_0 + (Ca_max - Ca_0) * (cutim / tau_Ca) * exp(1.0_rp - cutim / tau_Ca)
        ! Binding kinetics TnC-Ca^{2+} {Eq.2}
        T_0 = max(T_0, 1.0e-16_rp)
        dCaidt = rho_0 * Ca_i * (Ca_bmax - Ca_b) - rho_1 * (1.0_rp - T_a / (gamma * T_0)) * Ca_b
        Ca_b   = Ca_b + dtime * dCaidt

        if ( Ca_b>Ca_bmax ) then
            call runend("HUNTER_ORIGINAL: Something went terribly worng: Ca_b > Ca_bmax")
        end if

        if ((z < 0.0_rp) .or. (z > 1.0_rp)) then
            call runend("HUNTER_ORIGINAL: z is not in [0,1]")
        end if


        ! Tropomyosin kinetics (z)
        ! - Hill parameters (Eq.11) & (Eq.12)
        n    = n_ref*(1.0_rp + beta_1*(lamda_auxi - 1.0_rp))
        pC50 = pC50_ref*(1.0_rp + beta_2*(lamda_auxi - 1.0_rp))
        C50  = 10.0_rp**(6.0_rp - pC50) !pC50_ref needs to be calculated from Ca50 of exmedi with this formula
        ! - First order tropomsyosin kinetics (Eq.4)
        dzdt = alpha_0 * ((Ca_b / C50)**n * (1.0_rp - z) - z)
        z = z + dtime * dzdt

        ! Isometric tension (T0)
        !- Derivative of the isometric tension-length with respect of time:
        !        dT0/dt = dT0/dz*dz/dt + dT0/dlamda*dlamda/dt
        !    where T0 from (Eq. 9) and z  from (Eq. 5) for steady state conditions 
        !    or otherwise from (Eq. 4)
        dT0dt = T_ref * (1.0_rp + (lmb - 1.0_rp) * beta_0) * dzdt + (T_ref * beta_0 * z * dlmbdt) 
        T_0 = T_0 + dtime * dT0dt

        ! Internal state variables
        Q_1 = Q_1 * fexp(-alpha_1 * dtime) + A_1 * dlmbdt * dtime ! (1.0_rp / alpha_1) * (1.0_rp - fexp(-alpha_1 * dtime))
        Q_2 = Q_2 * fexp(-alpha_2 * dtime) + A_2 * dlmbdt * dtime ! (1.0_rp / alpha_2) * (1.0_rp - fexp(-alpha_2 * dtime))
        Q_3 = Q_3 * fexp(-alpha_3 * dtime) + A_3 * dlmbdt * dtime ! (1.0_rp / alpha_3) * (1.0_rp - fexp(-alpha_3 * dtime))
        Q  = Q_1 + Q_2 + Q_3

        ! Active traction
        T_a = max(T_0*(1.0_rp + a * Q)/(1.0_rp - Q),1.0e-16_rp)
  
        ! Convert to GS units
        tract  = T_a * 10000_rp

        ! Set current values
        sdv(1,new) = lmb
        sdv(2,new) = Ca_b
        sdv(3,new) = T_0
        sdv(4,new) = T_a
        sdv(5,new) = z
        sdv(6,new) = Q_1
        sdv(7,new) = Q_2
        sdv(8,new) = Q_3

        contains

           real(rp) pure function fexp(xval) 
               real(rp), intent(in) :: xval
               if( xval < -100_rp )then
                  fexp = 0.0_rp
               else
                  fexp = exp(xval)
               endif
               
           end function fexp

    end subroutine get_active_traction_HUNTER_ORIGINAL


    ! -------------------------------------------------------------------------
    !> @author  Francesc Levrero-Florencio
    !> @author  Adria Quintanas-Corominas (adria.quintanas@bsc.es)
    !> @author  Alfonso Santiago (alfonso.santiago@bsc.es)
    ! -------------------------------------------------------------------------
    subroutine get_active_traction_LANDS( &
        props, lamda, Ca, tract, dtract, sdv )
        use def_master,           only :  dtime
        ! -------------------------------
        implicit none
        ! -------------------------------
        real(rp),        intent(in)    :: props(:)  ! Properties
        real(rp),        intent(in)    :: lamda     ! Lambda (stretches)
        real(rp),        intent(in)    :: Ca        ! Calcium
        real(rp),        intent(inout) :: tract     ! Active traction
        real(rp),        intent(inout) :: dtract    ! Derivative active traction w.r.t lamda
        real(rp),        intent(inout) :: sdv(:,:)  ! State variables
        ! -------------------------------
        real(rp)                       :: k_TRPN, n_TRPN, TRPN_50, Ca50_ref
        real(rp)                       :: n_tm
        real(rp)                       :: k_u, k_u_2, k_uw, k_ws, k_wu, k_su
        real(rp)                       :: r_s, r_w
        real(rp)                       :: y_s, y_w, y_su, y_wu
        real(rp)                       :: zeta_s, zeta_w
        real(rp)                       :: c_s, c_w
        real(rp)                       :: phi
        real(rp)                       :: A_eff, A_s
        real(rp)                       :: beta_0, beta_1, beta_2, pC50, pC50_ref
        real(rp)                       :: T_ref
        real(rp)                       :: lamda_curr, lamda_auxi, lamda_prev, lamda_rate
        real(rp)                       :: B, S, W, h, Cai, Ca0, Ca50, CaTRPN, term
        integer(ip)                    :: ii
        integer(ip), parameter         :: old = 1_ip, new = 2_ip
        integer(ip), parameter         :: NR_ite_max = 100_ip
        real(rp),    parameter         :: NR_res_tol = 1.0e-10_rp
        integer(ip)                    :: NR_ite
        real(rp)                       :: NR_res_nrm
        real(rp),    dimension(6,2)    :: NR_sdv
        real(rp),    dimension(6)      :: NR_res, ODE_ydot
        real(rp),    dimension(6,6)    :: NR_jac
        real(rp),    dimension(6)      :: dRdlmd, dSdlmd
        real(rp)                       :: dhdlmd
        ! -------------------------------
        ! Assign properties
        k_TRPN   = props( 1)
        n_TRPN   = props( 2)
        pC50_ref = props( 3)
        k_u      = props( 4)
        n_tm     = props( 5)
        TRPN_50  = props( 6)
        k_uw     = props( 7)
        k_ws     = props( 8)
        r_w      = props( 9)
        r_s      = props(10)
        y_s      = props(11)
        y_w      = props(12)
        phi      = props(13)
        A_eff    = props(14)
        beta_0   = props(15)
        beta_1   = props(16)
        T_ref    = props(17)
        Ca0      = props(18)
        beta_2   = props(19)

        ! Compute lambda values
        lamda_curr = lamda
        lamda_auxi = min(lamda_thr,lamda_curr)
        lamda_prev = sdv(nsdvs_ecc+1,old)             ! Previous lambda from state variable
        lamda_rate = (lamda_curr - lamda_prev)/dtime  ! Lambda rate throught finite difference

        ! Compute Ca_50
        Cai      = Ca !max(Ca - Ca0,0.0_rp) 
        ! this ca50 handling is borrowed from Hunter or something like that, I can't remember now
        pC50     = pC50_ref*(1.0_rp + beta_2*(lamda_auxi - 1.0_rp))
        Ca50_ref = 10.0_rp**(6.0_rp - pC50)                
        Ca50     = Ca50_ref + beta_1*(lamda_auxi - 1.0_rp)

        ! Compute h
        h = max(0.0_rp, 1.0_rp + beta_0*(lamda_auxi + min(0.87_rp, lamda_auxi) - 1.87_rp))

        ! Compute dependent parameters 
        k_wu = k_uw*((1.0_rp/r_w)-1.0_rp)-k_ws
        k_su = k_ws*((1.0_rp/r_s)-1.0_rp)*r_w
        A_s = (A_eff*r_s)/((1.0_rp-r_s)*r_w+r_s)
        c_s = phi*k_ws*(1.0_rp-r_s)*r_w/r_s
        c_w = phi*k_uw*((1.0_rp-r_s)*(1.0_rp-r_w))/((1.0_rp-r_s)*r_w)
        k_u_2 = (k_u*(TRPN_50**n_tm))/(1.0_rp-r_s-(1.0_rp-r_s)*r_w)

        ! Loop for the Newton-Raphson of the Backward Euler integration scheme
        ! - initialization
        NR_sdv(1:6,old) = sdv(1:6,old)
        NR_sdv(1:6,new) = sdv(1:6,old)
        S      = NR_sdv(1,new)
        W      = NR_sdv(2,new)
        CaTRPN = NR_sdv(3,new)
        B      = NR_sdv(4,new)
        zeta_s = NR_sdv(5,new)
        zeta_w = NR_sdv(6,new)
        ! - loop
        NR_loop: do NR_ite = 1, NR_ite_max

            ! Regularisation of some state variable 
            ! - compute y_su
            y_su = 0.0_rp
            if(     zeta_s + 1.0_rp < 0.0_rp )then
                y_su = -y_s*(zeta_s + 1.0_rp )
            elseif( zeta_s + 1.0_rp > 1.0_rp )then
                y_su = y_s*zeta_s
            endif
            ! - compute y_wu
            y_wu = y_w*abs(zeta_w)

            ! Compute the RHS of the ODE system
            ODE_ydot(1) = k_ws*W - k_su*S - y_su*S
            ODE_ydot(2) = k_uw*(1.0_rp - B - S - W) - k_wu*W - k_ws*W - y_wu*W
            ODE_ydot(3) = k_TRPN*(((Cai/Ca50)**n_TRPN)*(1.0_rp - CaTRPN) - CaTRPN)
            ODE_ydot(4) = k_u_2*min((CaTRPN**(-n_tm/2.0_rp)),100.0_rp)*(1.0_rp - B - S - W) - &
            &   k_u*(CaTRPN**(n_tm/2.0_rp))*B
            ODE_ydot(5) = A_s*lamda_rate - c_s*zeta_s
            ODE_ydot(6) = A_s*lamda_rate - c_w*zeta_w

            ! Compute Jacobian of the linearised system of equations
            NR_jac(:,:) = 0.0_rp
            ! jac(1,:)
            NR_jac(1,1) = -k_su - y_su
            NR_jac(1,2) =  k_ws
            if(     zeta_s + 1.0_rp < 0.0_rp )then
                NR_jac(1,5) = -y_s
            elseif( zeta_s + 1.0_rp > 1.0_rp )then
                NR_jac(1,5) = y_s
            endif
            ! jac(2,:)
            NR_jac(2,1) = -k_uw
            NR_jac(2,2) = -k_uw - k_wu -k_ws - y_wu
            NR_jac(2,4) = -k_uw
            if( zeta_w < 0.0_rp )then
                NR_jac(2,6) = -y_w
            else
                NR_jac(2,6) = y_w
            endif
            ! jac(3,:)
            NR_jac(3,3) = -k_TRPN*(((Cai/Ca50)**n_TRPN) + 1.0_rp)
            ! jac(4,:)
            if (CaTRPN**(-n_tm/2.0_rp) <= 100.0_rp) then
                term = CaTRPN**(-n_tm/2.0_rp)
            else
                term = 100.0_rp
            end if
            NR_jac(4,1) = -k_u_2*term
            NR_jac(4,2) = -k_u_2*term
            if( CaTRPN**(-n_tm/2.0_rp) <= 100.0_rp )then
                NR_jac(4,3) = (-n_tm/2.0_rp)*(k_u_2*(1.0_rp - B - S - W)*(CaTRPN**(-(n_tm/2.0_rp) - 1.0_rp)) + &
                &   k_u*B*(CaTRPN**((n_tm/2.0_rp)-1.0_rp)))
            else
                NR_jac(4,3) = (-n_tm/2.0_rp)*k_u*B*(CaTRPN**((n_tm/2.0_rp)-1.0_rp))
            end if
            NR_jac(4,4) = -k_u_2*term - k_u*(CaTRPN**(n_tm/2.0_rp))
            ! jac(5,:)
            NR_jac(5,5) = -c_s
            ! jac(6,:)
            NR_jac(6,6) = -c_w
            ! jac = I_6x6 - jac()
            NR_jac(:,:) = -NR_jac(:,:)*dtime
            do ii = 1, 6
                NR_jac(ii,ii) = 1.0_rp + NR_jac(ii,ii)
            enddo

            ! Compute residual of the linearised system of equations 
            ! (previous solution are the state variables)
            NR_res(:) = NR_sdv(:,new) - NR_sdv(:,old) - dtime*ODE_ydot(:)

            ! Solve the system of the equations
            call invert(NR_jac,6_ip,6_ip)
            NR_sdv(:,new) = NR_sdv(:,new) - matmul(NR_jac(:,:),NR_res(:))
            S      = NR_sdv(1,new)
            W      = NR_sdv(2,new)
            CaTRPN = NR_sdv(3,new)
            B      = NR_sdv(4,new)
            zeta_s = NR_sdv(5,new)
            zeta_w = NR_sdv(6,new)


            ! Check the convergence
            ! - compute norm of the residual of the linearised system pf equations
            NR_res_nrm = 0.0_rp
            do ii = 1, 6
                NR_res_nrm = NR_res_nrm + NR_res(ii)**2
            enddo
            NR_res_nrm = sqrt(NR_res_nrm)
            ! - apply convergence criterion
            if( NR_res_nrm < NR_res_tol )then
                exit
            endif

        enddo NR_loop

        ! Compute the active tractions
        ! - S      = NR_sdv(1,2)
        ! - W      = NR_sdv(2,2)
        ! - CaTRN  = NR_sdv(3,2)
        ! - B      = NR_sdv(4,2)
        ! - zeta_s = NR_sdv(5,2)
        ! - zeta_w = NR_sdv(6,2)
        tract = h*(T_ref/r_s)*(NR_sdv(1,2)*(NR_sdv(5,2)+1.0_rp)+NR_sdv(2,2)*NR_sdv(6,2))
        ! - convert from kPa to GS units 
        tract = tract*10000.0_rp

        ! Compute the derivative of the active tractions w.r.t.  lambda
        ! - dR/dlamda
        dRdlmd(:) = 0.0_rp
        if( lamda_curr < lamda_thr )then
            dRdlmd(3) = beta_1*n_TRPN*k_TRPN*dtime*(1.0_rp - CaTRPN)* &
            &   ((Ca**n_TRPN)*(Ca50**(-n_TRPN - 1.0_rp)))
        endif
        dRdlmd(5) = - A_s
        dRdlmd(6) = - A_s
        ! - dSdlamda
        ! NR_jac = NR_jac_inv (see above)
        dSdlmd(:) = -matmul(NR_jac,dRdlmd)
        ! - dhdlambda
        dhdlmd = 0.0_rp
        if(     lamda_curr < ((1.87_rp*beta_0 - 1.0_rp)/(2.0_rp*beta_0)) )then
            dhdlmd = 0.0_rp
        elseif( lamda_curr < 0.87_rp) then
            dhdlmd = 2.0_rp*beta_0
        elseif( lamda_curr < 1.2_rp )then
            dhdlmd = beta_0
        end if
        ! - DtractDlmd
        dtract = (T_ref*10000_rp/r_s)*(dhdlmd*((zeta_s + 1.0_rp)*S + zeta_w*W) + &
        &   h*((dSdlmd(5)*S + zeta_s*dSdlmd(1)) + (dSdlmd(6)*W + zeta_w*dSdlmd(2))))

        ! Assign new state variables
        sdv(1:6,new) = NR_sdv(1:6,new)                      ! Computed through the Backeuler
        sdv(  nsdvs_ecc+1,new) = sdv( nsdvs_ecc+1,old)      !  We didnt modify lamda previus
        ! -------------------------------
    end subroutine get_active_traction_LANDS


    subroutine eccou_set_cell_type(celty)
        use def_kermod,           only :  kfl_celltype_fun
        implicit none
        real(rp), intent(in) :: celty(:)
        integer(ip)          :: ii

        if( kfl_celltype_fun < 0_ip ) return

        do ii = 1, size(celty)
            celty_ecc(ii) = real(celty(ii),rp)
        enddo
    end subroutine eccou_set_cell_type


    subroutine eccou_set_sdv_at_gp( ielem, gp_sdv, xx )
        implicit none
        integer(ip), intent(in)           :: ielem
        integer(ip), intent(in), optional :: xx 
        real(rp),    intent(in)           :: gp_sdv(:,:)
        integer(ip)                       :: ii, jj, igaus, icomp=1_ip
        if(  size(gp_sdv,dim=1) > nsdvs_ecc ) call runend('Trying to set more state than allocated')
        if(  present(xx) ) icomp = xx
        do igaus = 1, size(gp_sdv,dim=2)
            ii = (igaus - 1_ip) * nsdvs_ecc + 1_ip
            jj = ii + (nsdvs_ecc - 1_ip)
            state_ecc(ii:jj,ielem,icomp) = gp_sdv(:,igaus)
        enddo
    end subroutine eccou_set_sdv_at_gp


    subroutine eccou_get_sdv_at_gp( ielem, gp_sdv, xx )
        implicit none
        integer(ip), intent(in)           :: ielem
        integer(ip), intent(in), optional :: xx
        real(rp),    intent(inout)        :: gp_sdv(:,:)
        integer(ip)                       :: ii, jj, igaus, icomp=1
        if(  size(gp_sdv,dim=1) > nsdvs_ecc ) call runend('Trying to get more state than allocated')
        if( present(xx) ) icomp = xx
        gp_sdv = 0
        do igaus = 1, size(gp_sdv,dim=2)
            ii = (igaus - 1_ip) * nsdvs_ecc + 1_ip
            jj = ii + (nsdvs_ecc - 1_ip)
            gp_sdv(:,igaus) = state_ecc(ii:jj,ielem,icomp)
        enddo
    end subroutine eccou_get_sdv_at_gp


    subroutine eccou_set_ca( ca_per_node, scaling_coeff )
        ! to pass ca from exmedi to here
        implicit none
        integer(ip)                        :: i
        real(rp), dimension(:), intent(in) :: ca_per_node
        real(rp)                           :: scaling_coeff
        
        do i=1_ip, size(ca_per_node, kind=ip)
            calcium_ecc(i) = ca_per_node(i) * scaling_coeff
        end do
    end subroutine eccou_set_ca

    !subroutine eccou_get_ca( ca_per_node, scaling_coeff )
    !    ! to pass ca from here to exmedi
    !    implicit none
    !    integer(ip)                           :: i
    !    real(rp), dimension(:), intent(inout) :: ca_per_node
    !    real(rp)                              :: scaling_coeff
    !    
    !    do i=1_ip, size(ca_per_node, kind=ip)
    !        ca_per_node(i) = calcium_ecc(i) * scaling_coeff
    !    end do
    !end subroutine eccou_get_ca



    subroutine eccou_set_ca0(icell, imate, ca0)
        implicit none
        integer(ip), intent(in) :: icell, imate
        real(rp),    intent(in) :: ca0
        cell_ca0_ecc(icell,imate)  = ca0
    end subroutine eccou_set_ca0


    subroutine eccou_set_ca50(icell, imate, ca50)
        use def_master, only : intost
        implicit none
        integer(ip), intent(in) :: icell, imate
        real(rp),    intent(in) :: ca50

        if( kfl_ca50_provided(imate) .ne. 0_ip ) return

        if ( abs(ca50)<epsilon(1.0_rp) ) then
            call runend("ECCOUPLING::eccou_set_ca50: CA50 CANNOT BE 0")
        end if

        call messages_live("ECCOUPLING USING EXTERNAL CA50 TO CALCULATE pC50_Ref FOR MATERIAL "//trim(intost(imate))//" CELLTYPE "//trim(intost(icell)))
        cell_ca50_ecc(icell,imate)  = calc_pC50_ref(ca50)

    end subroutine eccou_set_ca50


    subroutine eccou_set_state_variables()
        use def_domain,   only :  nelem, ngaus, lmate, ltype
        implicit none
        integer(ip)            :: igaus, ii, jj, ielem, pelty, pgaus, imate
        if( kfl_state_initialised ) return

        do ielem = 1,nelem
            imate = lmate(ielem)
            pelty = ltype(ielem)
            pgaus = ngaus(pelty)
            if( kfl_eccty(imate) == EXMSLD_EMEC_HUNTER_ORIGINAL )then
                do igaus = 1,pgaus
                    ii = (igaus - 1_ip)*nsdvs_ecc + 1_ip
                    jj = ii + (nsdvs_ecc - 1_ip)
                    state_ecc(ii:jj,ielem,1) = get_state_HUNTER()
                end do
                state_ecc(:,ielem,2) = state_ecc(:,ielem,1)
            endif
        end do

        kfl_state_initialised = .true.

        contains 

           function get_state_HUNTER() result(statvar)
               real(rp) :: statvar(nsdvs_ecc)
               statvar = 0.0_rp
               statvar(1) = 1.0_rp    ! lmb_0
               statvar(3) = 1e-12_rp  ! T_0
               statvar(4) = 1e-12_rp  ! T_a
           end function

    end subroutine eccou_set_state_variables


    subroutine eccou_finalize_initial_conditions() 
        ! this should be called after all ca50 are calcuated and set
        use def_master, only : intost
        implicit none
        integer(ip) :: imate, ref_material, icellt

        do imate = 1,nmate
            if ( kfl_ca50_provided(imate) < 0_ip ) then
                if (kfl_ca50_provided(-1_ip*(kfl_ca50_provided(imate))) < 0_ip) then
                    call runend("ECCOUPLING: CAL50 FOR MATERIAL "//trim(intost(imate))//" HAS TOO MUCH RECURSION")
                end if

                ref_material = -1_ip * kfl_ca50_provided(imate)
                call messages_live("ECCOUPLING: USING CA50 FROM MATERIAL "//trim(intost(ref_material))//" FOR MATERIAL "//trim(intost(imate)) )

                do icellt = 1,ncelltypes_ecc
                    cell_ca50_ecc(icellt,imate) = cell_ca50_ecc(icellt,ref_material)
                end do
            end if
        end do
    end subroutine eccou_finalize_initial_conditions


    real(rp) function calc_pC50_ref( ca50 )
        implicit none
        real(rp) :: ca50

        if ( abs(ca50)<epsilon(1.0_rp) ) then
            call runend("ECCOUPLING::calc_pC50_ref: CA50 CANNOT BE 0")
        end if

        calc_pC50_ref = 6.0_rp-log10(ca50)  !ca50 comes in units of solidz
    end function calc_pC50_ref

end module mod_eccoupling
