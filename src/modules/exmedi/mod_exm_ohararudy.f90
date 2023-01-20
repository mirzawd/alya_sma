!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Exmedi
!> @{
!> @file    mod_exm_ohararudy.f90
!> @author  mixed
!> @date    2019-11-16
!> @brief   mod_exm_ohararudy
!> @details mod_exm_ohararudy
!> @}
!-----------------------------------------------------------------------
module mod_exm_ohararudy

    use def_master
    use def_exmedi
    use mod_eccoupling, only : kfl_exmsld_ecc, kfl_force_ca50_calculation 
    use mod_messages,   only : messages_live
    use mod_parall,     only : PAR_MY_CODE_RANK
    use mod_exm_drugs
    use mod_exm_memory
    use mod_memory,     only : memory_size
    use mod_exm_cellmodel
    use mod_exm_oharaprecalc
    use mod_eccoupling, only : EXMSLD_CELL_OHARA, EXMSLD_CELL_OHARA_INAPA, EXMSLD_EMEC_LAND, EXMSLD_EMEC_LAND_BIDIR, kfl_cellmod, kfl_eccty


    implicit none
    private

    !integer(ip), parameter, public :: & ! Minimum array lengths required by this model. exm_inivar can make these larger
    !   EXM_OHR_NVARS_vaulo_exm = 29_ip, &
    !   EXM_OHR_NVARS_vcolo_exm = 11_ip, &
    !   EXM_OHR_NVARS_vicel_exm = 26_ip
    integer(4), private, parameter :: state_ip = 8 !precision to store ints


    integer(ip), parameter  :: &
      !precalculated hardcoded initial conditions
      EXM_CELL_NORMAL_1000BCL             = 0_ip, & !ohara-rudy,    MYOCYTE=normal,        CYCLELENGTH=1000 (genderless)
      EXM_CELL_NORMAL_70BPM               = 1_ip, & !ohara-passini, MYOCYTE=normal,        CYCLELENGTH=857
      EXM_CELL_NORMAL_HUMAN_MALE_600BCL   = 2_ip, & !ohara-passini, MYOCYTE=normal male,   CYCLELENGTH=600
      EXM_CELL_NORMAL_HUMAN_FEMALE_600BCL = 3_ip, & !ohara-passini, MYOCYTE=normal female, CYCLELENGTH=600
      !forces recalculation of the initial conditions
      EXM_CELL_NORMAL_MODIFIED            = 4_ip  !modified

    integer(ip), parameter :: nviclo = 26_ip  ! number of elements in viclo_exm

    real(rp)              :: voltage_local               ! voltage local to this module


    type(ohara_precalc),     dimension(:), pointer :: precalc_per_imate => null()
    type(DRUG_CONDUCTANCES), dimension(:), pointer :: drug_effects_per_imate => null()


    real(rp), dimension(:,:), pointer :: &
      vaulo_exm => null(),&            ! Auxiliary variables for single cell solution
      vcolo_exm => null()              ! Concentrations for single cell solution
    real(rp), dimension(:),   pointer :: &
      viclo_exm => null()              ! Currents for single cell solution


    type, public :: ohara_constants
      real(rp) :: &
        gnal, gks, gk1, gna, gkr, pca, &
        nao, cao, ko, rgas, temp, farad, &
        leng, rad, vcell, ageo, acap, vmyo, vnsr, vjsr, vss, &
        kmcamk, acamk, bcamk, camko, kmkam, kmcam,&
        pkna, ahf, ahs, tauHL, gto, &
        aff, afs, kmn, k2n, zca, &
        gncx, pnak, gkatp, fkatp, gkb, cmdnmax, &
        pcap, pcana, pcak, pcanap, pcakp, &
        kna1, kna2, kna3, kasymm, wna, wca, wnaca, kcaon, kcaoff, qna, qca, kmcaact, &
        k1p, k1m, k2p, k2m, k3m, k4m, knai0, knao0, delta, &
        kki, kko, mgadp, mgatp, kmgatp, hbig, ep, khp, knap, kxkur, &
        pnab, pcab, gpca, kmcmdn, trpnmax, kmtrpn, bsrmax, &
        kmbsr, bslmax, kmbsl, csqnmax, kmcsqn, &
        dtimon, j, ccm
    end type

    type, public :: ohara_outputs
        real(rp) :: xioni, Inet, qnet
    end type

    public :: exm_oneohr

    !private :: exm_ohr_allocate_locals
    !private :: exm_ohr_deallocate_locals
    !private :: localvoltages2globalvolatges
    !private :: exm_oceohr
    !private :: exm_init_voltages
    !private :: localvoltages2globalvolatges

    !private :: save_ohara_state
    !private :: load_ohara_state

    public :: exm_ohara_constants_initialize
    public :: exm_ohara_execute_timestep
    public :: exm_ohara_ionicurrents
    public :: exm_ohara_init

contains


subroutine exm_ohara_init()
    ! initializes precalculate structure to be used in the fem part, not 0d initialization
    use def_domain,     only : nmate
    use mod_eccoupling, only : kfl_cellmod
    implicit none
    integer(ip) :: imate

    call memory_alloca(mem_modul(1:2,modul),'precalc_per_imate','exm_ohara_init', precalc_per_imate, nmate)     
    call memory_alloca(mem_modul(1:2,modul),'drug_effects_per_imate','exm_ohara_init', drug_effects_per_imate, nmate)     

    do imate = 1,nmate
        if (any(kfl_cellmod(imate) == (/EXMSLD_CELL_OHARA, EXMSLD_CELL_OHARA_INAPA/))) then
            if ( exm_drugs_ondrugs(imate) ) then
                drug_effects_per_imate(imate) = exm_drugs_calculate_effect(imate)
            else 
                drug_effects_per_imate(imate) % pca  = 0.0_rp
                drug_effects_per_imate(imate) % gkr  = 0.0_rp
                drug_effects_per_imate(imate) % gna  = 0.0_rp
                drug_effects_per_imate(imate) % gk1  = 0.0_rp
                drug_effects_per_imate(imate) % gnal = 0.0_rp
                drug_effects_per_imate(imate) % gks  = 0.0_rp        
                drug_effects_per_imate(imate) % Ito  = 0.0_rp        
            end if
    
            call precalc_per_imate(imate) % init( kfl_cellmod(imate), .TRUE. )   ! with safe exponent
        end if
    end do
end subroutine exm_ohara_init




character(100) function generate_state_filename(mat, icelltype)
    use def_master, only : namda
    implicit none    
    integer(ip), intent(in) :: mat, icelltype

    !generate_state_filename = adjustl(trim(namda))//'.cm_state.m'//trim(intost(mat))//&
    !                          "c"//trim(intost(icelltype))//&
    !                          "rp"//trim(intost(int(rp, kind=ip)))//&
    !                          "ip"//trim(intost(int(ip, kind=ip)))//".bin"

    generate_state_filename = adjustl(trim(namda))//'.cm_state.m'//trim(intost(mat))//&
                              "c"//trim(intost(icelltype))//&
                              "rp"//trim(intost(int(rp, kind=ip)))//&
                              "ip8.bin"
end function


integer(ip) function get_cellmodel_file_handle( mat, icelltype )
   implicit none
   integer(ip), intent(in) :: mat, icelltype

   get_cellmodel_file_handle = 1000000_ip + 100_ip*mat + icelltype   
   if (IPARALL) get_cellmodel_file_handle = get_cellmodel_file_handle + PAR_MY_CODE_RANK*10000

end function get_cellmodel_file_handle



subroutine save_ohara_state( mat, icelltype, ohara_type, outputs )
   !!! Make sure integers are saved with precicion 8
   use mod_restart,   only : restart_add, restart_ini, restart_end
   use def_master,    only : ITASK_WRITE_RESTART
   implicit none

   integer(ip), intent(in)                     :: mat, icelltype, ohara_type
   TYPE(CELL_MODEL_OUTPUTS), intent(inout)     :: outputs

   integer(ip)                                 :: post_file_handle
   integer(ip)                                 :: iostat

   post_file_handle = get_cellmodel_file_handle( mat, icelltype )
   
   !save ohara_type, BCL, conductances, drugs, vaulo_exm(:, 1), vcolo_exm(:, 1), voltage, outputs
   !save initial volatge_local, viclo_exm(:, 1), vcolo_exm(:, 1:3), vaulo_exm(:, 1:2)
   OPEN (unit=post_file_handle, file=trim(generate_state_filename(mat, icelltype)), form="unformatted", &
         access='stream', status="REPLACE", iostat = iostat)

   if (iostat == 0_ip) then
      !save drugs
      call exm_drugs_save_to_file( post_file_handle )

      write( post_file_handle ) int(ohara_type, kind=state_ip)

      !save kfl_exmsld_ecc, and if true, save kfl_force_ca50_calculation(mat) 
      write( post_file_handle ) kfl_exmsld_ecc
      if (kfl_exmsld_ecc) write( post_file_handle ) kfl_force_ca50_calculation(mat) 

      !save exm_drugs_ondrugs(mat)
      write( post_file_handle ) exm_drugs_ondrugs(mat)

      !save conductances
      write( post_file_handle ) ttparmate_exm(icelltype, :,  mat)

      !save moneclmate_exm(:,mat)  and kfl_hfmod_exm(mat), kfl_inaga_exm(mat)
      !write( post_file_handle ) moneclmate_exm(1,mat), moneclmate_exm(2,mat), kfl_hfmod_exm(mat), kfl_inaga_exm(mat)
      write( post_file_handle ) int(moneclmate_exm(1,mat), kind=state_ip), int(moneclmate_exm(2,mat), kind=state_ip)

      !save steady state parameters
      write( post_file_handle ) int(kfl_steadystate_variable(mat), kind=state_ip), steady_tol_cellm(mat), timestep_cellm(mat)

      !save vaulo_exm, vcolo_exm, elmlo_exm 
      write( post_file_handle ) vaulo_exm, vcolo_exm, voltage_local
      
      !save outputs
      write( post_file_handle ) outputs % S
      write( post_file_handle ) outputs % W
      write( post_file_handle ) outputs % CaTRPN
      write( post_file_handle ) outputs % B
      write( post_file_handle ) outputs % zeta_s
      write( post_file_handle ) outputs % zeta_w
      write( post_file_handle ) outputs % Ca50   
      write( post_file_handle ) outputs % Lambda 
      write( post_file_handle ) outputs % Vinit  
      write( post_file_handle ) int( outputs % nbeats, kind=state_ip )
      write( post_file_handle ) outputs % toler  
      write( post_file_handle ) outputs % rmse   
      write( post_file_handle ) int( outputs % success, kind=state_ip )
      write( post_file_handle ) outputs % dt


      CLOSE(post_file_handle)
   else
      call messages_live( "Failed to open "//trim(generate_state_filename(mat, icelltype))//" for writing", "WARNING" )
   end if

end subroutine save_ohara_state

logical(lg) function load_ohara_state( mat, icelltype, ohara_type, outputs )
   !!! Make sure integers are read with precicion 8
   use mod_restart, only   : restart_add, restart_ini, restart_end
   use def_master,  only   : ITASK_READ_RESTART
   implicit none
   integer(ip), intent(in)                     :: mat, icelltype, ohara_type
   TYPE(CELL_MODEL_OUTPUTS), intent(inout)     :: outputs
   integer(ip)                                 :: iostat
   integer(ip)                                 :: post_file_handle
   integer(state_ip)                           :: ivar, ivar1
!   integer(state_ip)                           :: ivar2, ivar3
   logical(lg)                                 :: lvar
   real(rp)                                    :: fvar, fvar1
   real(rp)                                    :: fbuffer( max(nvars_ttparmate_exm, nauxi_exm, nconc_exm) )



   load_ohara_state = .TRUE.

   post_file_handle = get_cellmodel_file_handle( mat, icelltype )

   OPEN (unit=post_file_handle, file=trim(generate_state_filename(mat, icelltype)), form="unformatted", &
         access='stream', status="OLD", iostat = iostat)


   if (iostat == 0_ip) then
      !load drugs
      load_ohara_state = load_ohara_state .and. exm_drugs_compare_to_file( post_file_handle )

      read( post_file_handle, iostat = iostat, err = 10 ) ivar
      load_ohara_state = load_ohara_state .and. ( ivar == ohara_type )

      !read kfl_exmsld_ecc, and if true, save kfl_force_ca50_calculation(mat) 
      read( post_file_handle, iostat = iostat, err = 10  ) lvar
      load_ohara_state = load_ohara_state .and. ( lvar .eqv. kfl_exmsld_ecc )
            
      if (kfl_exmsld_ecc) then
         read( post_file_handle, iostat = iostat, err = 10  ) lvar
         load_ohara_state = load_ohara_state .and. ( lvar .eqv. kfl_force_ca50_calculation(mat) )         
      end if

      !read exm_drugs_ondrugs(mat)
      read( post_file_handle, iostat = iostat, err = 10  ) lvar
      load_ohara_state = load_ohara_state .and. ( lvar .eqv. exm_drugs_ondrugs(mat) )         
      
      !moneclmate_exm : 2_ip x nmate  , integer
      !vaulo_exm : nauxi_exm x 3_ip  (:,1) real
      !vcolo_exm : nconc_exm x 3_ip  vcolo_exm(:, 1) real
      !elmlo_exm (2), real

      !read conductances
      read( post_file_handle, iostat = iostat, err = 10  ) fbuffer(1:nvars_ttparmate_exm)
      load_ohara_state = load_ohara_state .and. &
            all( abs( fbuffer(1:nvars_ttparmate_exm) - ttparmate_exm(icelltype, 1:nvars_ttparmate_exm,  mat) )<epsilon(1.0_rp) )

      !read moneclmate_exm(:,mat)  and kfl_hfmod_exm(mat), kfl_inaga_exm(mat)
      !read( post_file_handle, iostat = iostat, err = 10  ) ivar, ivar1, ivar2, ivar3
      read( post_file_handle, iostat = iostat, err = 10  ) ivar, ivar1
      load_ohara_state = load_ohara_state .and. ( ivar  == moneclmate_exm(1,mat) )
      load_ohara_state = load_ohara_state .and. ( ivar1 == moneclmate_exm(2,mat) )
      !load_ohara_state = load_ohara_state .and. ( ivar2 == kfl_hfmod_exm(mat)    )
      !load_ohara_state = load_ohara_state .and. ( ivar3 == kfl_inaga_exm(mat)    )

      !read steady state parameters
      read( post_file_handle, iostat = iostat, err = 10  ) ivar, fvar, fvar1
      load_ohara_state = load_ohara_state .and. ( ivar == kfl_steadystate_variable(mat) )
      load_ohara_state = load_ohara_state .and. ( abs(fvar - steady_tol_cellm(mat))<epsilon(1.0_rp) )
      load_ohara_state = load_ohara_state .and. ( abs(fvar - timestep_cellm  (mat))<epsilon(1.0_rp) )

      !if inputs match, proceed loading the state variables
      if ( load_ohara_state .and. (iostat == 0) ) then
         !read vaulo_exm, vcolo_exm, elmlo_exm 
         read( post_file_handle, iostat = iostat, err = 10  ) vaulo_exm, vcolo_exm, voltage_local
         
         !read outputs
         read( post_file_handle, iostat = iostat, err = 10  ) outputs % S
         read( post_file_handle, iostat = iostat, err = 10  ) outputs % W
         read( post_file_handle, iostat = iostat, err = 10  ) outputs % CaTRPN
         read( post_file_handle, iostat = iostat, err = 10  ) outputs % B
         read( post_file_handle, iostat = iostat, err = 10  ) outputs % zeta_s
         read( post_file_handle, iostat = iostat, err = 10  ) outputs % zeta_w
         read( post_file_handle, iostat = iostat, err = 10  ) outputs % Ca50   
         read( post_file_handle, iostat = iostat, err = 10  ) outputs % Lambda 
         read( post_file_handle, iostat = iostat, err = 10  ) outputs % Vinit  
         read( post_file_handle, iostat = iostat, err = 10  ) outputs % nbeats 
         read( post_file_handle, iostat = iostat, err = 10  ) outputs % toler  
         read( post_file_handle, iostat = iostat, err = 10  ) outputs % rmse   
         read( post_file_handle, iostat = iostat, err = 10  ) outputs % success   
         read( post_file_handle, iostat = iostat, err = 10  ) outputs % dt
      end if

10    if ( iostat .ne. 0 ) then
         load_ohara_state = .FALSE.
      end if

      CLOSE(post_file_handle)

   else
      load_ohara_state = .FALSE.
   end if
      
end function load_ohara_state


subroutine exm_ohr_allocate_locals()
   use mod_memory, only : memory_alloca
   implicit none

   call memory_alloca(mem_modul(1:2,modul),'vaulo_exm', 'mod_exm_ohararudy', vaulo_exm,  nauxi_exm, 2_ip) 
   call memory_alloca(mem_modul(1:2,modul),'vcolo_exm', 'mod_exm_ohararudy', vcolo_exm,  nconc_exm, 2_ip) 
   call memory_alloca(mem_modul(1:2,modul),'viclo_exm', 'mod_exm_ohararudy', viclo_exm,  nviclo)
   vaulo_exm = 0.0_rp
   vcolo_exm = 0.0_rp
   viclo_exm = 0.0_rp

end subroutine exm_ohr_allocate_locals


subroutine exm_ohr_deallocate_locals()
   use mod_memory, only : memory_deallo
   implicit none

   call memory_deallo(mem_modul(1:2,modul),'vaulo_exm', 'mod_exm_ohararudy', vaulo_exm) 
   call memory_deallo(mem_modul(1:2,modul),'vcolo_exm', 'mod_exm_ohararudy', vcolo_exm) 
   call memory_deallo(mem_modul(1:2,modul),'viclo_exm', 'mod_exm_ohararudy', viclo_exm) 
end subroutine exm_ohr_deallocate_locals






!------------------------------------------------------------------------
!> @addtogroup Exmedi
!> @{
!> @file    exm_oceohr.f90
!> @author  Jazmin Aguado-Sierra
!> @brief   Single cell run for Initial condition setup for Ohara-Rudy 2011 heterogeneous model
!> @date   16/NOV/1966
!> @details Runs a single cell simulation at th given frequency and pathologic conditions \n
!!    It performs single cell runs under normal, heart failure or drugs \n
!> @} 
!!-----------------------------------------------------------------------
subroutine exm_ohara_ionicurrents(ipoin,xioni,dioni,cai)

    use      def_parame
    use      def_master
    use      def_elmtyp
    use      def_domain
    use      def_exmedi
    use      def_kermod,     only : kfl_celltype_fun
    use      mod_eccoupling
    use      mod_exm_cellmodel
    use      mod_exm_drugs


    ! definition of variables
    implicit none
    integer(ip), intent(in) :: ipoin !< node
    real(rp), intent(out)   :: xioni   !< current
    real(rp), intent(out)   :: dioni !< current derivative
    real(rp), intent(out)   :: cai !< current calcium
    integer(ip)             :: imate,ituss_exm


    !old! real(rp)    :: vinf, xitaux, vaux0, vaux1, vaux2, vaux3 
    !old! real(rp)    :: vffrt, vfrt, ena, ek, eks, bt, a_rel, btp, a_relp
    !old! real(rp)    :: v1, v2 
    !old! real(rp)    :: a1, a2, a3, a4, b1, b2, b3, b4, k3p, k4p
    !old! real(rp)    :: h1, h2, h3, h4, h5, h6, h7, h8, h9, h10, h11, h12, x1, x2, x3, x4
    !old! real(rp)    :: k1, k2, k3, k4, k5, k6, k7, k8, e1, e2, e3, e4, k3pp, k4pp, knai
    !old! real(rp)    :: phical, phicana, phicak
    !old! real(rp)    :: hp, hh, finalp, ii, afcaf, afcas, fcap 
    !old! real(rp)    :: ksca  
    !old! real(rp)    :: hca, hna, allo, zna, jncxna, jncxca
    !old! real(rp)    :: knao, pp, zk
    !old! real(rp)    :: jnakna, jnakk, xkb, ipp, aif, ais, ikatp
    !old! real(rp)    :: fp, axrf, axrs, xr, rkr, rk1, fjupp 
    !old! real(rp)    :: fca, ff, ficalp
    !old! real(rp)    :: finap, fjrelp
    !old! real(rp)    :: jupnp, jupp, fitop
    !old! real(rp)    :: cai, cass, cansr, cajsr, nass, nai, ki, kss, camkt
    !old! real(rp)    :: rhsx1, rhsx2, rhsx, val0
    !old! real(rp)    :: HFjleak, INaCa, cai

    !old! real(rp)    :: vinf, xitaux, vaux0, vaux1, vaux2, vaux3, a2bas
    !old! real(rp)    ::  vffrt, vfrt, ena, ek, eks, bt, a_rel, btp, a_relp, ccm, Inet
    !old! real(rp)    :: pkna, farad,  v1, v2
    !old! real(rp)    ::  a1, a2, a3, a4, b1, b2, b3, b4, k1p, k2p, k3p, k4p, k1m, k2m, k3m, k4m
    !old! real(rp)    :: h1, h2, h3, h4, h5, h6, h7, h8, h9, h10, h11, h12, x1, x2, x3, x4
    !old! real(rp)    :: k1, k2, k3, k4, k5, k6, k7, k8, e1, e2, e3, e4, k3pp, k4pp
    !old! real(rp)    :: knai,knao
    !old! real(rp)    :: gna, gnal, gto, gks, gk1, gncx, gkb, gkr, gpca, phical, phicana, phicak, nao
    !old! real(rp)    :: ahf, ahs, hp, hh, finalp, ii, aff, afs, afcaf, afcas, fcap
    !old! real(rp)    :: zca, pca, pcap, pcak, pcanap, pcakp, ksca, kasymm
    !old! real(rp)    :: wna, wca, wnaca, kcaon, kcaoff, qna, qca, hca, hna, kmcaact, allo, zna, jncxna, jncxca
    !old! real(rp)    :: delta, kki, kko, mgadp, mgatp, kmgatp, ep, khp, knap, kxkur, pp, zk
    !old! real(rp)    :: jnakna, jnakk, pnak, xkb, pnab, ipp, aif, ais
    !old! real(rp)    :: fp, ko, axrf, axrs, xr, rkr, rk1, fjupp, pcab, pcana
    !old! real(rp)    :: cao, ff, fca, ficalp, finap,  fjrelp
    !old! real(rp)    :: rgas, temp, jupnp, jupp, fitop, hbig, camko, kmcam
    !old! real(rp)    :: cass, cansr, cajsr, nass, nai, ki, kss, camkt
    !old! real(rp)    :: rhsx1, rhsx2, rhsx, val0, bslmax, bsrmax,cmdnmax, csqnmax
    !old! real(rp)    :: leng, rad, vcell, ageo, acap, vmyo, vnsr, vjsr, vss, kmcamk, kmcmdn
    !old! real(rp)    :: acamk, bcamk, trpnmax, kmcsqn, kmtrpn, kmbsl, kmbsr, INaCa
    !old! real(rp)    :: dtimeEP
    !old! real(rp)    :: sac, ikatp, gkatp, fkatp, tauHL, HFjleak


    !old! real(rp) :: inf_results(EXM_OHR_INF_EQUATIONS), tau_results(EXM_OHR_TAU_EQUATIONS)
    !old! real(rp), parameter :: KNA1 = 15.0_rp
    !old! real(rp), parameter :: KNA2 = 5.0_rp
    !old! real(rp), parameter :: KNA3 = 88.12_rp
    !old! real(rp), parameter :: INV_KNA1 = 1.0_rp / KNA1
    !old! real(rp), parameter :: INV_KNA2 = 1.0_rp / KNA2
    !old! real(rp), parameter :: INV_KNA3 = 1.0_rp / KNA3
    !old! real(rp), parameter :: KNAI0 = 9.073_rp
    !old! real(rp), parameter :: KNAO0 = 27.78_rp
    logical :: has_land

    type(ohara_constants)   :: params
    type(ohara_outputs)     :: oo
    integer(ip)             :: model_type
    real(rp)                :: a2bas, sac

    if(kfl_celltype_fun == 0_ip) then !if there is no celltype field
        call runend("exm_ohara_ionicurrents: Cell type filed is required")
    end if


    if(INOTMASTER) then
        imate      = nodemat(ipoin)
        model_type = kfl_cellmod(imate)
            
        !old! dtimeEP=dtime * 1000.0_rp
        if (modab_exm == 0_ip) then ! no gradient
            a2bas = 1.0_rp
        else                        ! gradient
            a2bas = atbhe_exm(1,ipoin)
        end if

        ituss_exm = int(celty_exm(ipoin), kind=ip)  
        if ( all( model_type .ne. (/EXMSLD_CELL_OHARA, EXMSLD_CELL_OHARA_INAPA/)) ) call runend("exm_ohara_ionicurrents called for a node whose cell model is not OHara")
        if (imate==0) call runend('Point '//trim(intost(ipoin))//' (local numbering) has material 0')
        

        !!!!  DRUG DEFINITION TO ICAL, IKR, INA, IK1, INAL
        if ( exm_drugs_ondrugs(imate) ) then
            call exm_ohara_constants_initialize( params, ttparmate_exm(:, :, imate), timestep_cellm(imate), ituss_exm, model_type, drug_effects_per_imate(imate) )
        else 
            call exm_ohara_constants_initialize( params, ttparmate_exm(:, :, imate), timestep_cellm(imate), ituss_exm, model_type )
        endif
        params % dtimon = dtime * 1000.0_rp        

        !-----------------------------------------------------------------
        !
        ! Apex-base heterogeneity
        !
        params % gks    = params % gks * a2bas !add apex-base heterogeneity
        !
        !
        !-----------------------------------------------------------------
      
        sac=0.0_rp !initialize sac current
      
        has_land = .false.
        if (kfl_exmsld_ecc) then
            if (kfl_eccty(imate) == EXMSLD_EMEC_LAND .or. kfl_eccty(imate) == EXMSLD_EMEC_LAND_BIDIR) has_land = .TRUE.
        end if

        !-----------------------------------------
        !
        !call precalc % init (model_type, .TRUE.)
        call precalc_per_imate(imate) % set_voltage( elmag(ipoin,ITER_K) )
        !call times(1) % ini()
        call precalc_per_imate(imate) % calculate()
        !call times(1) % add()
        !
        !------------------------------------------
        oo % qnet   = qneto_exm(ipoin) !TODO: qnet and inet should be reset to 0 after each beat as in 0d model???
        oo % Inet   = 0.0_rp
        oo % xioni  = 0.0_rp

        call exm_ohara_execute_timestep( params,  ttparmate_exm(:, :, imate), model_type, has_land, ituss_exm, &
            vconc(:,ipoin,2), vconc(:,ipoin,1), vauxi_exm(:,ipoin,2), vauxi_exm(:,ipoin,1), vicel_exm(:,ipoin), &
            elmag(ipoin,ITER_K), precalc_per_imate(imate), oo ) 

        xioni = oo%xioni
        qneto_exm(ipoin) = oo % qnet

        if (has_land) then
            call exm_ohaland_calcium(params%kmcmdn, params%cmdnmax, params%vnsr, &
                params%vmyo, params%vss, params%acap, params%farad, params % dtimon, &
                ipoin, cai, sac, elmag(ipoin,ITER_K), imate)

            vconc(1,ipoin,1) = cai

            if (kfl_eccty(imate) == EXMSLD_EMEC_LAND_BIDIR) then
                xioni = xioni + sac
            end if
        end if

        

        vauxi_exm(:,ipoin,2)=vauxi_exm(:,ipoin,1)         
        dioni = 0.0_rp
        cai   = vconc(1,ipoin,1)
  end if

end subroutine exm_ohara_ionicurrents



!------------------------------------------------------------------------
!> @addtogroup Exmedi
!> @{
!> @file    exm_ohaland_calcium.f90
!> @author  Jazmin Aguado-Sierra and Francesc Levrero-Florencio
!> @brief   Single cell run for Initial condition setup for Ohara-Rudy 2011 heterogeneous model
!> @date   16/NOV/1966
!> @details Runs a single cell simulation at th given frequency and pathologic conditions \n
!!    It performs single cell runs under normal, heart failure or drugs \n
!> @} 
!!-----------------------------------------------------------------------
subroutine exm_ohaland_calcium(kmcmdn,cmdnmax,vnsr,vmyo,vss,acap,farad,dtimeEP,ipoin,cai,sac,volt,imate)
    !subroutine exm_ohaland_calcium(kmcmdn,kmtrpn,cmdnmax,vnsr,vmyo,vss,acap,farad,dtimeEP,ipoin,cai,sac,volt)
    
      use      def_parame
      use      def_master
      use      def_elmtyp
      use      def_domain
      use      def_exmedi
      use      mod_eccoupling, only: troponin_ecc, props_ecc
    
      implicit none
    
      real(rp), intent(in) :: kmcmdn,cmdnmax,vnsr,vmyo,vss,acap,farad,dtimeEP,volt
      !real(rp), intent(in) :: kmtrpn
      real(rp), intent(out) :: cai, sac
      real(rp) :: vaux1,vaux3,vaux4,rhsx1,rhsx2,rhsx,val0,k_TRPN,n_TRPN,Ca50,Ca50_ref,beta_1,lambda,  &
       C_tens(3,3),dCaTRPN,lambda0,gsac,esac,Ca50_fact
    
      integer(ip) , intent(in) :: ipoin, imate
      integer(ip) :: idime,jdime,kdime
    
      ! Declare some parameters of the Land model
      k_TRPN =   props_ecc(1 , imate) * 0.001_rp !0.1_rp
      n_TRPN =   props_ecc(2 , imate) !2.0_rp
      Ca50_ref = props_ecc(3 , imate)/1000.0_rp !0.805_rp, convert from solidz units
      beta_1 =   props_ecc(16, imate) !-2.4_rp
      Ca50_fact = props_ecc(20, imate)!1.0_rp by defaults
    
      ! Recover the right Cauchy-Green strain tensor
      C_tens = 0.0_rp
      do idime = 1,ndime
        do jdime = 1,ndime
          do kdime = 1,ndime
            C_tens(idime,jdime) = C_tens(idime,jdime)+gdepo(kdime,idime,ipoin)*                   &
             gdepo(kdime,jdime,ipoin)
          end do
        end do
      end do
    
      ! Recover the stretch in the fibre direction
      lambda = 0.0_rp
      do idime = 1,ndime
        do jdime = 1,ndime
          lambda = lambda+fiber(idime,ipoin)*C_tens(idime,jdime)*fiber(jdime,ipoin)
        end do
      end do
      lambda = SQRT(lambda)
    
      ! Calculate calcium sensitivity
      lambda0 = lambda
      lambda = MIN(1.2_rp,lambda)
      Ca50 = (Ca50_ref+beta_1*(lambda-1.0_rp))*Ca50_fact
    
    
      ! Calcium calcium in the coupled model
      vaux1 = (kmcmdn+vconc(1,ipoin,2)) * (kmcmdn+vconc(1,ipoin,2))
      vaux3 = 1.0_rp/(1.0_rp + (cmdnmax*kmcmdn/vaux1))
      rhsx1 = vicel_exm(14,ipoin) + vicel_exm(13,ipoin) - (2.0_rp*vicel_exm(8,ipoin))
      vaux4 = (vconc(1,ipoin,2)/(Ca50_ref*Ca50_fact))**n_TRPN * (1.0_rp-(troponin_ecc(ipoin,1)*0.001_rp))
      dCaTRPN = k_TRPN *( vaux4 - (troponin_ecc(ipoin,1)*0.001_rp))
      rhsx2 = -(vicel_exm(18,ipoin)*vnsr/vmyo) + (vicel_exm(15,ipoin)*vss/vmyo) - (0.07_rp * dCaTRPN)
      rhsx = vaux3 *(-(rhsx1 *acap/(2.0_rp*farad*vmyo)) + rhsx2)
      val0 = vconc(1,ipoin,2)
      cai = val0 + dtimeEP*rhsx
    
      ! Calculate SAC current
      if (lambda0 >= 1.0_rp) then
        gsac = 0.1_rp
        esac = 0.0_rp
        sac = gsac*((lambda0-1.0_rp)/(1.1_rp-1.0_rp))*(volt-esac)
      else
        sac = 0.0_rp  
      end if
    
    end subroutine exm_ohaland_calcium



!------------------------------------------------------------------------
!> @addtogroup Exmedi
!> @{
!> @file    exm_oneohr.f90
!> @date    12/04/2013
!> @author  Jazmin Aguado-Sierra
!> @brief   Sets initial conditions for Ohara-Rudy 2011 model
!> @details Runs for changes of cycle length and drug administration \n
!> @}
!------------------------------------------------------------------------
! icelltype -- epi, endo, mid (1-3)
! Returns:
! success_status == 0_ip -- sucess, >0 -- error, see EXM_ONEOHR_ERROR* for the meaning
   subroutine exm_oneohr(mat, icelltype, outputs)
      implicit none
      integer(ip), intent(in)                   :: mat, icelltype
      TYPE(CELL_MODEL_OUTPUTS), intent(inout)   :: outputs

      logical(lg)                               :: has_land
      integer(ip)                               :: ohara_type

      call exm_ohr_allocate_locals()

      outputs % success = EXM_CELLMODEL_NOTINITIALIZED

      ohara_type = kfl_cellmod(mat)

      has_land = .FALSE.

      if (kfl_exmsld_ecc) then
         if (kfl_eccty(mat) == EXMSLD_EMEC_LAND .or. kfl_eccty(mat) == EXMSLD_EMEC_LAND_BIDIR) has_land = .TRUE.
      end if


!---------------------------------------------------------------------
!
!      This piece of code allows using hardcoded initial conditions. disabled for now, since the model saves the state after running into .bin files      
!
!---------------------------------------------------------------------
!     
!      if (kfl_hfmodmate_exm(mat) == EXM_MYOCYTE_NORMAL) then
!
!         if ((moneclmate_exm(2, mat) == 1000) .and. ( .not. exm_drugs_ondrugs(mat) ) .and. (ohara_type == EXMSLD_CELL_OHARA) ) then  
!            !if it's Normal and at 1000 bpm, not ina passini, taken from the paper
!            call messages_live("EXMEDI LOADING HARDCODED CELL MODEL FOR MATERIAL "//trim(intost(mat))//": 1000ms, NOT INA_PASSINI")
!            call exm_init_voltages(EXM_CELL_NORMAL_1000BCL, mat, icelltype, outputs)
!
!            if ( kfl_exmsld_ecc ) then
!               call messages_live("MATERIAL "//trim(intost(mat))//" CELL MODEL: 1000ms, NOT INA_PASSINI CANNOT PASS CA50 TO THE MODELS. MAKE SURE YOU SET CAL50","WARNING")
!            end if
!            outputs % success = EXM_CELLMODEL_LOADED
!            
!         else if ((moneclmate_exm(2, mat) == 857) .and. ( .not.  exm_drugs_ondrugs(mat) ) .and. (ohara_type == EXMSLD_CELL_OHARA_INAPA)) then   
!            !Normal at 70 bpm
!            call messages_live("EXMEDI LOADING HARDCODED CELL MODEL FOR MATERIAL "//trim(intost(mat))//": 70bpm, INA_PASSINI")
!            call exm_init_voltages(EXM_CELL_NORMAL_70BPM, mat, icelltype, outputs)
!            call localvoltages2globalvolatges(mat, icelltype)
!            outputs % success = EXM_CELLMODEL_LOADED
!
!         else if((moneclmate_exm(2,mat) == 600) .and. (kfl_hfmod_exm(mat) == EXM_MYOCYTE_MALE) .and.( .not. exm_drugs_ondrugs(mat) ) .and. (ohara_type == EXMSLD_CELL_OHARA_INAPA)) then  !if it's a Normal MALE HUMAN  AT 600BCL
!            call messages_live("EXMEDI LOADING HARDCODED CELL MODEL FOR MATERIAL "//trim(intost(mat))//": 600ms, MALE, INA_PASSINI")
!            call exm_init_voltages(EXM_CELL_NORMAL_HUMAN_MALE_600BCL, mat, icelltype, outputs)
!            call localvoltages2globalvolatges(mat, icelltype)
!
!            if ( kfl_exmsld_ecc ) then
!               call messages_live("MATERIAL "//trim(intost(mat))//" CELL MODEL: 600ms, MALE, INA_PASSINI CANNOT PASS CA50 TO THE MODELS. MAKE SURE YOU SET CAL50","WARNING")
!            end if
!            outputs % success = EXM_CELLMODEL_LOADED
!
!         else if((moneclmate_exm(2,mat) == 600) .and. (kfl_hfmod_exm(mat) == EXM_MYOCYTE_FEMALE) .and.( .not. exm_drugs_ondrugs(mat) ) .and. (ohara_type == EXMSLD_CELL_OHARA_INAPA)) then  !if it's a Normal FEMALE HUMAN at 600BCL
!            call messages_live("EXMEDI LOADING HARDCODED CELL MODEL FOR MATERIAL "//trim(intost(mat))//": 600ms, FEMALE, INA_PASSINI")
!            call exm_init_voltages(EXM_CELL_NORMAL_HUMAN_FEMALE_600BCL, mat, icelltype, outputs)
!            call localvoltages2globalvolatges(mat, icelltype)
!
!            if ( kfl_exmsld_ecc ) then
!               call messages_live("MATERIAL "//trim(intost(mat))//" CELL MODEL: 600ms, FEMALE, INA_PASSINI CANNOT PASS CA50 TO THE MODELS. MAKE SURE YOU SET CAL50","WARNING")
!            end if
!
!            outputs % success = EXM_CELLMODEL_LOADED
!         end if
!      end if

      if ( outputs % success == EXM_CELLMODEL_NOTINITIALIZED ) then ! includes (kfl_hfmodmate_exm(mat) == EXM_MYOCYTE_MODIFIED)

            call exm_init_voltages(EXM_CELL_NORMAL_MODIFIED, mat, icelltype, outputs)

            if (kfl_user_specified_celltypes_exm(mat, icelltype) .eq. 1_ip) then
               if( .not. load_ohara_state( mat, icelltype, ohara_type, outputs ) ) then
                  call exm_oceohr(mat, icelltype, ohara_type, has_land, outputs)
                  call save_ohara_state( mat, icelltype, ohara_type, outputs )
               else
                  !call messages_live("Successfully loaded O'Hara-Passini state for material "//trim(intost(mat))//" celltype "//trim(intost(icelltype)),"REPORT")
                  outputs % success = EXM_CELLMODEL_LOADED
               end if
            end if

            call localvoltages2globalvolatges(mat, icelltype)

      end if

      call exm_ohr_deallocate_locals()

   end subroutine exm_oneohr

   subroutine localvoltages2globalvolatges(mat, ituss_exm)
      implicit none

      integer(ip), intent(in) :: mat, ituss_exm

      if (ituss_exm == EXM_CELLTYPE_EPI) then  !epi
         vauxi_exm_initial(:, ituss_exm, mat) = vaulo_exm(:, 1)
         vconc_initial(:, ituss_exm, mat) = vcolo_exm(:, 1)
         vminimate_exm(ituss_exm, mat) = voltage_local
      else if (ituss_exm == EXM_CELLTYPE_ENDO) then  !endo
         vauxi_exm_initial(:, ituss_exm, mat) = vaulo_exm(:, 1)
         vconc_initial(:, ituss_exm, mat) = vcolo_exm(:, 1)
         vminimate_exm(ituss_exm, mat) = voltage_local
      else if (ituss_exm == EXM_CELLTYPE_MID) then  !mid
         vauxi_exm_initial(:, ituss_exm, mat) = vaulo_exm(:, 1)
         vconc_initial(:, ituss_exm, mat) = vcolo_exm(:, 1)
         vminimate_exm(ituss_exm, mat) = voltage_local
      end if
   end subroutine



   subroutine exm_ohara_constants_initialize( c, material_params, timestep, celltype, model_type, drug_effects )
      !material_params = ttparmate_exm(:, :, imate), 
      implicit none
      type(ohara_constants),    intent(inout) :: c
      integer(ip),              intent(in)    :: celltype, model_type
      real(rp), dimension(:,:), intent(in)    :: material_params
      real(rp),                 intent(in)    :: timestep
      type(DRUG_CONDUCTANCES), optional       :: drug_effects

      c % gna  = 75.0_rp
      c % pca  = 0.0001_rp
      c % gkr  = 0.046_rp
      c % gk1  = 0.1908_rp
      c % gnal = 0.0075_rp
      c % gks  = 0.0034_rp
      c % gto  = 0.02_rp 

      if ( present(drug_effects) ) then
         c % gna  = c % gna  * drug_effects % gna 
         c % pca  = c % pca  * drug_effects % pca 
         c % gkr  = c % gkr  * drug_effects % gkr 
         c % gk1  = c % gk1  * drug_effects % gk1 
         c % gnal = c % gnal * drug_effects % gnal
         c % gks  = c % gks  * drug_effects % gks   
         c % gto  = c % gto  * drug_effects % Ito
      end if

      !extracellular ionic concentrations
      c % nao = 140.0_rp * material_params(1, 13)
      c % cao = 1.8_rp   * material_params(2, 13)
      c % ko  = 5.4_rp   * material_params(3, 13)
   
      !physical constants
      c % rgas  = 8314.0_rp
      c % temp  = 310.0_rp
      c % farad = 96485.0_rp
   
      !cell geometry
      c % leng  = 0.01_rp
      c % rad   = 0.0011_rp
      c % vcell = 1000.0_rp * 3.14_rp * c % rad * c % rad * c % leng
      c % ageo  = 2.0_rp    * 3.14_rp * c % rad * c % rad + 2.0_rp * 3.14_rp * c % rad * c % leng
      c % acap  = 2.0_rp    * c % ageo
      c % vmyo  = 0.68_rp   * c % vcell
      c % vnsr  = 0.0552_rp * c % vcell
      c % vjsr  = 0.0048_rp * c % vcell
      c % vss   = 0.02_rp   * c % vcell
   
      !%camk constants
      c % kmcamk = 0.15_rp
      c % acamk  = 0.05_rp
      c % bcamk  = 0.00068_rp
      c % camko  = 0.05_rp
      c % kmcam  = 0.0015_rp
   
      !Current constants
      c % pkna = 0.01833_rp
      !!!!!!   calculate ina
      c % ahf  = 0.99_rp
      c % ahs  = 1.0_rp - c % ahf
      !gnal = 0.0075_rp

      c % tauHL= 200.0_rp 
   
      !!! calculate ical
      !!!!!! calculate ff
      c % aff = 0.6_rp
      c % afs = 1.0_rp - c % aff
      !!!!! calculate icana and icak
      c % kmn = 0.002_rp
      c % k2n = 1000.0_rp
      c % zca = 2.0_rp
   
      !!  inal
      !!!  calculate IKr
   
      !!!  calculate IKs
      !gks = 0.0034_rp
      !!!! calculate IK1
   
      c % gncx  = 0.0008_rp
      c % pnak  = 30.0_rp
      !!IKATP Channel for ischemia
      c % gkatp = 0.05_rp
      c % fkatp = 1.0_rp
      !!calculate IKb
      c % gkb   = 0.003_rp
   
      !%calcium buffer constants
      c % cmdnmax = 0.05_rp

      if (celltype == EXM_CELLTYPE_EPI) then !!epi
         c % gnal    = c % gnal    *0.6_rp * material_params(3, 6 )
         c % gto     = c % gto     *4.0_rp * material_params(3, 1 )
         c % pca     = c % pca     *1.2_rp * material_params(3, 9 )  !epi
         c % gkr     = c % gkr     *1.3_rp * material_params(3, 4 )
         c % gk1     = c % gk1     *1.2_rp * material_params(3, 3 )
         c % gncx    = c % gncx    *1.1_rp * material_params(3, 7 )
         c % pnak    = c % pnak    *0.9_rp * material_params(3, 10)
         c % gkb     = c % gkb     *0.6_rp * material_params(3, 8 )
         c % cmdnmax = c % cmdnmax *1.3_rp * material_params(3, 11)
         c % gks     = c % gks     *1.4_rp * material_params(3, 2 )
         c % gna     = c % gna             * material_params(3, 5 )
         c % fkatp   = c % fkatp           * material_params(3, 14)
         c % tauHL   =            200.0_rp * material_params(3, 17)

      else if (celltype == EXM_CELLTYPE_MID) then !!MID
         c % gto     = c % gto    *4.0_rp * material_params(2, 1 )
         c % pca     = c % pca    *2.5_rp * material_params(2, 9 )  !mid
         c % gkr     = c % gkr    *0.8_rp * material_params(2, 4 )
         c % gks     = c % gks            * material_params(2, 2 )
         c % gk1     = c % gk1    *1.3_rp * material_params(2, 3 )
         c % gncx    = c % gncx   *1.4_rp * material_params(2, 7 )
         c % pnak    = c % pnak   *0.7_rp * material_params(2, 10)
         c % gna     = c % gna            * material_params(2, 5 )
         c % gnal    = c % gnal           * material_params(2, 6 )
         c % gkb     = c % gkb            * material_params(2, 8 )
         c % cmdnmax = c % cmdnmax        * material_params(2, 11)
         c % fkatp   = c % fkatp          * material_params(2, 14)
         c % tauHL   =           200.0_rp * material_params(2, 17)
         if ( model_type == EXMSLD_CELL_OHARA_INAPA ) then
            c % pca = (c % pca*1.8_rp)/2.5_rp   !mid !Minchole modification
         end if
      else if (celltype == EXM_CELLTYPE_ENDO) then  !ENDO
         c % gto     = c % gto              * material_params(1, 1 )
         c % pca     = c % pca     * 1.0_rp * material_params(1, 9 )  !endo
         c % gkr     = c % gkr              * material_params(1, 4 )
         c % gks     = c % gks              * material_params(1, 2 )
         c % gk1     = c % gk1              * material_params(1, 3 )
         c % gncx    = c % gncx             * material_params(1, 7 )
         c % pnak    = c % pnak             * material_params(1, 10)
         c % gna     = c % gna              * material_params(1, 5 )
         c % gnal    = c % gnal             * material_params(1, 6 )
         c % gkb     = c % gkb              * material_params(1, 8 )
         c % cmdnmax = c % cmdnmax          * material_params(1, 11)
         c % fkatp   = c % fkatp            * material_params(1, 14)
         c % tauHL   =             200.0_rp * material_params(1, 17)
      end if

      c % pcap   = 1.1_rp       * c % pca
      c % pcana  = 0.00125_rp   * c % pca
      c % pcak   = 0.0003574_rp * c % pca
      c % pcanap = 0.00125_rp   * c % pcap
      c % pcakp  = 0.0003574_rp * c % pcap
   
      !%calculate inaca_i
      c % kna1    = 15.0_rp
      c % kna2    = 5.0_rp
      c % kna3    = 88.12_rp
      c % kasymm  = 12.5_rp
      c % wna     = 60000.0_rp
      c % wca     = 60000.0_rp
      c % wnaca   = 5000.0_rp
      c % kcaon   = 1500000.0_rp
      c % kcaoff  = 5000.0_rp
      c % qna     = 0.5224_rp
      c % qca     = 0.1670_rp
      c % kmcaact = 0.000150_rp !% and to calculate inaca_ss
   
      !%calculate inak
      c % k1p    = 949.5_rp
      c % k1m    = 182.4_rp
      c % k2p    = 687.2_rp
      c % k2m    = 39.4_rp
      c % k3m    = 79300.0_rp
      c % k4m    = 40.0_rp
      c % knai0  = 9.073_rp
      c % knao0  = 27.78_rp
      c % delta  = -0.1550_rp
      c % kki    = 0.5_rp
      c % kko    = 0.3582_rp
      c % mgadp  = 0.05_rp
      c % mgatp  = 9.8_rp
      c % kmgatp = 0.0000001698_rp
      c % hbig   = 0.00000010_rp
      c % ep     = 4.2_rp
      c % khp    = 0.0000001698_rp
      c % knap   = 224.0_rp
      c % kxkur  = 292.0_rp
   
      !%calculate inab
      c % pnab = 0.000000000375_rp
   
      !%calculate icab
      c % pcab = 0.000000025_rp
   
      !%calculate ipca
      c % gpca    = 0.0005_rp
      c % kmcmdn  = 0.00238_rp
      c % trpnmax = 0.07_rp
      c % kmtrpn  = 0.0005_rp
      c % bsrmax  = 0.047_rp
      c % kmbsr   = 0.00087_rp
      c % bslmax  = 1.124_rp
      c % kmbsl   = 0.0087_rp
      c % csqnmax = 10.0_rp
      c % kmcsqn  = 0.8_rp
   
      if( timestep >= 0.0_rp ) then
        c % dtimon = timestep * 1000.0_rp   !pass to miliseconds
      else
        c % dtimon = 0.005_rp   !0.005_rp  0.0015_rp
      end if

      c % ccm    = 1.0_rp
   
      !c % Inet   = 0.0_rp
      !c % qnet   = 0.0_rp
      !c % xioni  = 0.0_rp

   end subroutine exm_ohara_constants_initialize


! Single cell run for Initial condition setup for Ohara-Rudy 2011 heterogeneous model
! Runs a single cell simulation at th given frequency and pathologic conditions \n
! It performs single cell runs under normal, heart failure or drugs \n
!-----------------------------------------------------------------------
!mat - material id
!ituss_exm -  3-Epicardium, 1-endocardium, 2-midmyocardium
!model_type - EXMSLD_CELL_OHARA, EXM_OCEOHR_LAND
!extra_outputs - real 1D array preallocated
!Returns: 0 - sucess, >0 - error, <0 not executed
!If has_land=TRUE returns land parameters in outputs
!else land_parameters are not set
   subroutine exm_oceohr(mat, ituss_exm, model_type, has_land, outputs)
      use mod_exm_cellmodel_convergence
      use def_master, only : namda

      ! definition of variables
      implicit none
      integer(ip), intent(in)                     :: mat, ituss_exm, model_type
      logical(lg), intent(in)                     :: has_land
      TYPE(CELL_MODEL_OUTPUTS), intent(inout)     :: outputs
   
   
      integer(ip) :: i, nsamples_per_beat
      real(rp)    :: Ca_max, Ca_min, aux1, dur, stim
      integer(ip) :: nbeat, itim, j
      real(rp)    :: batim, atim

      !LAND variables
      real(rp)    :: k_TRPN, n_TRPN, Ca50_ref, k_u, n_tm, TRPN_50, k_uw, &    
                    k_ws, r_w, r_s, y_s, y_w, phi, A_eff, beta_0, beta_1, &  
                    T_ref, k_wu, k_su, A_s, c_s, c_w, k_u_2, &
                    lambda, lambda_rate, Ca50, &
                    S, W, CaTRPN, B, zeta_s, zeta_w,&
                    vaux1, vaux3, rhsx1, rhsx2, rhsx,&
                    y_su, y_wu, ydot(6)




      !number of voltages actually saved in the array, needed to verify that there was no error
      integer(ip) :: ss_array_current_id
      real(rp)    :: ss_epsilon ! threshold squared to declare steady state when comparing two beats
      real(rp)    :: ss_min_voltage_threshold
   
      integer(ip)                       :: post_frequecy, post_file_handle, dump_file_handle
      logical                           :: save_variables
      
      type(ConvergenceChecker)          :: convergence
      type(ohara_constants)             :: params

      type(ohara_precalc)               :: precalc
      type(ohara_outputs)               :: oo

      !to save calcium and voltage and other variables
      save_variables = ( kfl_save_convergence_cellmodel == 1_ip )
      post_frequecy = 200_ip
   
      post_file_handle = -1_ip
      
      !sanity check
      SELECT CASE (model_type)
      CASE (EXMSLD_CELL_OHARA)
         ! nothing for now
      CASE (EXMSLD_CELL_OHARA_INAPA)
         ! nothing for now
      CASE DEFAULT
         call runend("mod_exm_oceohr: Unknown cell model: "//trim(intost(model_type)))
      END SELECT
   
      call precalc % init(model_type, .FALSE.)

   
      !write file header
      if (save_variables) then
         post_file_handle = get_cellmodel_file_handle( mat, ituss_exm )

         OPEN (unit=post_file_handle, file=adjustl(trim(namda))//'.cm_0d.m'//trim(intost(mat))//"c"//trim(intost(ituss_exm))//".csv", form="FORMATTED", status="REPLACE")
   
   
         if (.not. has_land) then
            WRITE (post_file_handle,*) "Global_Time Beat_Time Calcium Voltage QNet"        ! running ohara without land
         else
            WRITE (post_file_handle,*) "Global_Time Beat_Time Calcium Voltage CaTRPN Qnet" ! running ohara with land
         end if
            
      end if
   
      !!!!  DRUG DEFINITION TO ICAL, IKR, INA, IK1, INAL
      if ( exm_drugs_ondrugs(mat) ) then
         call exm_ohara_constants_initialize( params, ttparmate_exm(:, :, mat), timestep_cellm(mat), ituss_exm, model_type, drug_effects_per_imate(mat) )
      else 
         call exm_ohara_constants_initialize( params, ttparmate_exm(:, :, mat), timestep_cellm(mat), ituss_exm, model_type )
      endif
   
      oo % Inet   = 0.0_rp
      oo % qnet   = 0.0_rp
      oo % xioni  = 0.0_rp


      !---------------------------------------
      !
      !   Land constants
      !
      !------------------------------------------
      k_TRPN   = 0.1_rp
      n_TRPN   = 2.0_rp
      Ca50_ref = 0.805_rp !solidz uses different units, Ca in solidz = Ca in exmedi * 1000
      k_u      = 1.0_rp
      n_tm     = 5.0_rp
      TRPN_50  = 0.35_rp
      k_uw     = 0.182_rp
      k_ws     = 0.012_rp
      r_w      = 0.5_rp
      r_s      = 0.25_rp
      y_s      = 0.0085_rp
      y_w      = 0.615_rp
      phi      = 2.23_rp
      A_eff    = 25.0_rp
      beta_0   = 2.3_rp
      beta_1   = -2.4_rp
      T_ref    = 120.0_rp
      k_wu     = k_uw * ((1.0_rp/r_w) - 1.0_rp) - k_ws
      k_su     = k_ws * ((1.0_rp/r_s) - 1.0_rp) * r_w
      A_s      = (A_eff * r_s)/((1.0_rp - r_s)*r_w + r_s)
      c_s      = phi * k_ws * (1.0_rp - r_s) * r_w/r_s
      c_w      = phi * k_uw * ((1.0_rp - r_s) * (1.0_rp - r_w))/((1.0_rp - r_s) * r_w)
      k_u_2    = (k_u * (TRPN_50**n_tm))/(1.0_rp - r_s - (1.0_rp - r_s) * r_w)

      !Declare initial variables
      lambda      = 1.0_rp
      lambda_rate = 0.0_rp
      Ca50        = Ca50_ref + beta_1*(lambda - 1.0_rp)

      ! Declare initial values of state variables
      S      = 0.0_rp
      W      = 0.0_rp
      CaTRPN = 0.000001_rp
      B      = 1.0_rp
      zeta_s = 0.0_rp
      zeta_w = 0.0_rp

      !---------------------------------------
      !
      !   End of Land constants
      !
      !------------------------------------------



      nbeat  = 0_ip
      itim   = 0_ip
      j      = 0_ip
      atim   = 0.0_rp
      batim  = 0.0_rp
      stim   = 0.0_rp

   
      !----------------------------------------------
      !
      ! begin: Initialize variables for steady state detection
      !
      !----------------------------------------------
      !allocate the arrays to store two full beats
      nsamples_per_beat = nint(moneclmate_exm(2, mat)/ params % dtimon, kind=ip)
   
      if (steady_tol_cellm(mat) >= 0) then
         ss_epsilon = steady_tol_cellm(mat)
      else
         select case (kfl_steadystate_variable(mat))
         case (EXM_CELL_STEADY_VOLTAGE)
            ss_epsilon = 1.0_rp    !this is maximum allowed  mean voltage difference
         case (EXM_CELL_STEADY_CALCIUM)
            ss_epsilon = 1.0e-7_rp 
         case default
            call runend("EXM_OCEOHR: Unknown name of the variable to determine the cell convergence.")
         end select
      end if
      call convergence % Init( ss_epsilon, nsamples_per_beat )
      call convergence % SetMaxSteps( 3_ip )
   
      ss_min_voltage_threshold = -75.0_rp !the min volytage within the beat has to be lower than this, to consider that we are in steady state
      !----------------------------------------------
      !
      ! end: Initialize variables for steady state detection
      !
      !----------------------------------------------
   
      outputs % success = EXM_CELLMODEL_NOTCONVERGED !return error by default
   
      ! moneclmate_exm(1, mat) - number of beats
      ! moneclmate_exm(2, mat) - cycle length
      do while (batim < (moneclmate_exm(1, mat)*moneclmate_exm(2, mat)))
         Ca_max = 0.0_rp !reset Ca max for the next beat
         Ca_min = 1.0e+30_rp !reset
   
         aux1 = 0.0_rp
         atim = params % dtimon * aux1
   
         !for steady state verification
         ss_array_current_id = 1_ip
         oo % qnet = 0.0_rp
         oo % Inet = 0.0_rp
         do while (atim < moneclmate_exm(2, mat))
            !if (atim < moneclmate_exm(2) )then
            aux1 = aux1 + 1.0_rp
            atim = params % dtimon * aux1
            dur  = 0.5_rp ! duration time (ms)
            if (atim < dur) then
               ! External stimulus current I_stim [pA/pF]
               stim = -80.0_rp     ! Value of I_stim current [pA/pF]
            else
               stim = 0.0_rp
            end if

            !add vcolo(,2) as input, and vcolo(,1) as output, pass inout viclo(:,1)
            !add vaulo_in, vaulo_out(:,:) as inout
            call exm_ohara_execute_timestep( params,  ttparmate_exm(:, :, mat), model_type, has_land, ituss_exm, &
                           vcolo_exm(:, 2), vcolo_exm(:, 1), vaulo_exm(:, 2), vaulo_exm(:, 1), viclo_exm, voltage_local, precalc, oo ) 
               
            voltage_local = voltage_local + params % dtimon * ( -(oo % xioni + stim) ) ! *params % farad/params % acap


            if( has_land ) then
                vaux1 = (params % kmcmdn + vcolo_exm(1, 2) )*(params % kmcmdn + vcolo_exm(1, 2) )
                vaux3 = 1.0_rp/(1.0_rp + (params % cmdnmax*params % kmcmdn/vaux1))
                rhsx1 = viclo_exm(14) + viclo_exm(13) - (2.0_rp*viclo_exm(8))
                rhsx2 = -(viclo_exm(18) * params % vnsr / params % vmyo) + (viclo_exm(15)*params % vss/params % vmyo) -&
                         (k_TRPN * (((vcolo_exm(1, 2)*1000.0_rp / Ca50)**n_TRPN)*(1.0_rp-CaTRPN)-&
                         CaTRPN)) * params % trpnmax

                rhsx = vaux3*(-(rhsx1*params % acap/(2.0_rp*params % farad*params % vmyo)) + rhsx2)
                vcolo_exm(1, 1) = vcolo_exm(1, 2) + params % dtimon*rhsx             !!! cai
            else 
                if ( vcolo_exm(1, 1) < 0.0_rp ) call runend("Land not coupled but exm_ohara_execute_timestep did not calculate Cai")
            end if


            !------------------------------------------------------
            ! 
            !   Ca50 stuff
            !
            if ( Ca_max < vcolo_exm(1, 1) * 1000.0_rp) then 
               Ca_max = vcolo_exm(1, 1) * 1000.0_rp !change of units
            elseif ( Ca_min > vcolo_exm(1, 1) * 1000.0_rp) then
               Ca_min = vcolo_exm(1, 1) * 1000.0_rp !change of units
            endif
            ! 
            !   Ca50 stuff
            !
            !------------------------------------------------------
   

            !------------------------------------------------------
            ! 
            !   Land stuff
            !
            !----------------------------------------------------
            if (has_land) then
       
                ! State variable dependent parameters (regularised parameters)
                if ((zeta_s + 1.0_rp) < 0.0_rp) then
                   y_su = y_s*(-zeta_s - 1.0_rp)
                elseif ((zeta_s + 1.0_rp) > 1.0_rp) then
                   y_su = y_s*zeta_s
                else
                   y_su = 0.0_rp
                end if
                y_wu = y_w*ABS(zeta_w)
       
                ! Define the RHS of the ODE system of Land
                ydot(1) = k_ws*W - k_su*S - y_su*S
                ydot(2) = k_uw*(1.0_rp - B - S - W) - k_wu*W - k_ws*W - y_wu*W
                ydot(3) = k_TRPN*(((vcolo_exm(1, 1)*1000.0_rp/Ca50)**n_TRPN)*(1.0_rp - CaTRPN) - CaTRPN)
                ydot(4) = k_u_2*MIN((CaTRPN**(-n_tm/2.0_rp)), 100.0_rp)*(1.0_rp - B - S - W) - k_u*(CaTRPN**(n_tm/2.0_rp))*B
                ydot(5) = A_s*lambda_rate - c_s*zeta_s
                ydot(6) = A_s*lambda_rate - c_w*zeta_w
       
                ! Update variables of Land
                S      = S      + params % dtimon*ydot(1)
                W      = W      + params % dtimon*ydot(2)
                CaTRPN = CaTRPN + params % dtimon*ydot(3)
                B      = B      + params % dtimon*ydot(4)
                zeta_s = zeta_s + params % dtimon*ydot(5)
                zeta_w = zeta_w + params % dtimon*ydot(6)                
            end if 
            !------------------------------------------------------
            ! 
            !   End of Land
            !
            !----------------------------------------------------

            vcolo_exm(:, 2) = vcolo_exm(:, 1)
            vaulo_exm(:, 2) = vaulo_exm(:, 1)



   
            itim = itim + 1_ip
   
            if (save_variables) then
               if (j >= post_frequecy) then  !it was 10 for saving
   
                  if (.not. has_land ) then
                     write (post_file_handle, "(5(e16.8E3, ' '))") itim* params % dtimon, atim, vcolo_exm(1, 1), voltage_local, oo % qnet
                  else 
                     write (post_file_handle, "(6(e16.8E3, ' '))") itim* params % dtimon, atim, vcolo_exm(1, 1), voltage_local, CaTRPN, oo % qnet
                  end if
                     
                  j = 0
               end if
               j = j + 1
            end if !save_variables
   
            
            select case (kfl_steadystate_variable(mat))
            case (EXM_CELL_STEADY_VOLTAGE)
               convergence % array (ss_array_current_id) = voltage_local
            case (EXM_CELL_STEADY_CALCIUM)
               convergence % array (ss_array_current_id) = vcolo_exm(1, 1)
            case default
               call runend("EXM_OCEOHR: Unknown name of the variable to determine the cell convergence.")
            end select
   
            ss_array_current_id = ss_array_current_id + 1_ip
   
   
         end do
   
         batim = itim * params % dtimon
   
         !Update Ca50
         if( kfl_exmsld_ecc ) then ! TODO: the condition (.not. has_land) should disappear eventually
            if ( (.not. has_land) .or. ( kfl_force_ca50_calculation(mat) ) ) then
               Ca50 = (Ca_min+Ca_max) * 0.5_rp 
            end if
         else
            Ca50 = (Ca_min+Ca_max) * 0.5_rp 
         end if

         !--------------------------------------------
         !
         !  begin: check if we reached steady state. If so, break
         !
         !----------------------------------------
         if (ss_array_current_id - 1 .NE. nsamples_per_beat) then
            call runend("EXM_OCEOHR: Something went wrong in the loop of one beat. The number of timesteps executed is incorrect.")
         end if
   
         if ( convergence % converged() ) then
            outputs % success = EXM_CELLMODEL_CONVERGED !success
            exit
         end if
   
   
         call convergence % EndIteration()         
         !--------------------------------------------
         !
         !  end: check if we reached steady state
         !
         !----------------------------------------
   
      end do
   
      !------------------------------------------------
      !
      ! Save stats
      !
      !----------------------------------------------
      outputs % nbeats = itim/nsamples_per_beat
      outputs % toler = convergence % GetTolerance()
      outputs % rmse  = convergence % GetMetricValue()
   
      if (save_variables) close (post_file_handle)
      
   
      !-----------------------------------------------
      !
      !   Extra outputs
      !
      outputs % S      = S      ! land only
      outputs % W      = W      ! land only
      outputs % CaTRPN = CaTRPN * 1000.0_rp ! land only change from ms to s
      outputs % B      = B      ! land only
      outputs % zeta_s = zeta_s ! land only
      outputs % zeta_w = zeta_w ! land only
      outputs % Lambda = 1.0_rp ! not set here, introduced in TOR  
      outputs % Ca50   = Ca50 / 1000.0_rp   ! general. Denormalize, so that it has the same units as Ca0, and can be renormalized again in mod_exm_sld_eccoupling.f90
      outputs % dt     = params % dtimon / 1000.0_rp! ms -> s 
      !
      !   Extra outputs end
      !
      !-----------------------------------------------------
   
   
      !--------------------------------------------------------
      ! 
      ! This fragment dumps variables for predefining the conditions
      ! 
      !--------------------------------------------------------
      if ( kfl_save_init_cellmodel==1_ip ) then
         dump_file_handle = 2000000 + 100*mat + ituss_exm
   
         if (IPARALL) dump_file_handle = dump_file_handle + PAR_MY_CODE_RANK*10000 
   
         OPEN (unit=dump_file_handle, file='ohara_dump.m'//trim(intost(mat))//"c"//trim(intost(ituss_exm))//".csv", form="FORMATTED", status="REPLACE")
   
         WRITE(dump_file_handle,*) "voltage_local = "//trim(retost(voltage_local))//"_rp"
   
         WRITE(dump_file_handle,*) "viclo_exm(1:26) = 0.0_rp"
         do i=1,memory_size(vcolo_exm, 1_ip)
            WRITE(dump_file_handle,*) "vcolo_exm("//trim(intost(i))//", 1:2) = "//trim(retost(vcolo_exm(i, 1)))//"_rp"
         end do
         do i=1,memory_size(vaulo_exm, 1_ip)
            WRITE(dump_file_handle,*) "vaulo_exm("//trim(intost(i))//", 1:2) = "//trim(retost(vaulo_exm(i, 1)))//"_rp"
         end do

         WRITE(dump_file_handle,*) "outputs % S      = "//trim(retost(outputs % S     ))//"_rp"
         WRITE(dump_file_handle,*) "outputs % W      = "//trim(retost(outputs % W     ))//"_rp" 
         WRITE(dump_file_handle,*) "outputs % CaTRPN = "//trim(retost(outputs % CaTRPN))//"_rp" 
         WRITE(dump_file_handle,*) "outputs % B      = "//trim(retost(outputs % B     ))//"_rp" 
         WRITE(dump_file_handle,*) "outputs % zeta_s = "//trim(retost(outputs % zeta_s))//"_rp" 
         WRITE(dump_file_handle,*) "outputs % zeta_w = "//trim(retost(outputs % zeta_w))//"_rp" 
         WRITE(dump_file_handle,*) "outputs % Lambda = "//trim(retost(outputs % Lambda))//"_rp" 
         WRITE(dump_file_handle,*) "outputs % Ca50   = "//trim(retost(outputs % Ca50  ))//"_rp" 
         WRITE(dump_file_handle,*) "outputs % nbeats = "//trim(intost(outputs % nbeats))//"_ip"
         WRITE(dump_file_handle,*) "outputs % toler  = "//trim(retost(outputs % toler ))//"_rp"
         WRITE(dump_file_handle,*) "outputs % toler  = "//trim(retost(outputs % toler ))//"_rp"
         WRITE(dump_file_handle,*) "outputs % rmse   = "//trim(retost(outputs % rmse  ))//"_rp"
         WRITE(dump_file_handle,*) "outputs % dt     = "//trim(retost(outputs % dt    ))//"_rp"

         CLOSE(dump_file_handle)
      end if
   
   
   end subroutine exm_oceohr
   
   
   !modifies params (Inet, qnet, S, W, CaTRPN, B, zeta_s, zeta_w) , viclo
   !output_params.xioni is calculated, Inet and qnet are updated, so make sure they are initialized before
   subroutine exm_ohara_execute_timestep( params, material_params, model_type, has_land, cell_type, &
                                          vcolo_in, vcolo_out, vaulo_in, vaulo_out, viclo, &
                                          voltage, precalc, output_params )
      use mod_memory, only : memory_size
      implicit none
      type(ohara_constants),     intent(in)      :: params
      type(ohara_outputs),       intent(inout)   :: output_params
      real(rp), dimension(:),    intent(in)      :: vcolo_in 
      real(rp), dimension(:),    intent(out)     :: vcolo_out
      real(rp), dimension(:),    intent(inout)   :: viclo
      real(rp), dimension(:),    intent(in)      :: vaulo_in
      real(rp), dimension(:),    intent(out)     :: vaulo_out
      integer(ip),               intent(in)      :: model_type
      integer(ip),               intent(in)      :: cell_type
      logical(lg),               intent(in)      :: has_land
      real(rp), dimension(:,:),  intent(in)      :: material_params ! = ttparmate(:,:,imate)
      real(rp),                  intent(in)      :: voltage
      type(ohara_precalc),       intent(inout)   :: precalc

      integer(ip) :: i
      real(rp)    :: vinf, xitaux, vaux0, vaux1, vaux2, vaux3 
      real(rp)    :: vffrt, vfrt, ena, ek, eks, bt, a_rel, btp, a_relp
      real(rp)    :: v1, v2 
      real(rp)    :: a1, a2, a3, a4, b1, b2, b3, b4, k3p, k4p
      real(rp)    :: h1, h2, h3, h4, h5, h6, h7, h8, h9, h10, h11, h12, x1, x2, x3, x4
      real(rp)    :: k1, k2, k3, k4, k5, k6, k7, k8, e1, e2, e3, e4, k3pp, k4pp, knai
      real(rp)    :: phical, phicana, phicak
      real(rp)    :: hp, hh, finalp, ii, afcaf, afcas, fcap 
      real(rp)    :: ksca  
      real(rp)    :: hca, hna, allo, zna, jncxna, jncxca
      real(rp)    :: knao, pp, zk
      real(rp)    :: jnakna, jnakk, xkb, ipp, aif, ais, ikatp
      real(rp)    :: fp, axrf, axrs, xr, rkr, rk1, fjupp 
      real(rp)    :: fca, ff, ficalp
      real(rp)    :: finap, fjrelp
      real(rp)    :: jupnp, jupp, fitop
      real(rp)    :: cai, cass, cansr, cajsr, nass, nai, ki, kss, camkt
      real(rp)    :: rhsx1, rhsx2, rhsx, val0
      real(rp)    :: HFjleak, INaCa

      !real(rp)    :: delepi, ta, devel, recov, fss, vaux4, tif,  tis, tj, iss, anca, dnca
      !real(rp)    :: tff, tfcaf, ths, tm, jss, km2n
   

      if ( size(vcolo_in , kind=ip) .ne. nconc_exm ) call runend("exm_ohara_execute_timestep: vcolo_in has incorrect number of elements") 
      if ( size(vcolo_out, kind=ip) .ne. nconc_exm ) call runend("exm_ohara_execute_timestep: vcolo_out has incorrect number of elements") 
      if ( size(vaulo_in , kind=ip) .ne. nauxi_exm ) call runend("exm_ohara_execute_timestep: vaulo_in has incorrect number of elements") 
      if ( size(vaulo_out, kind=ip) .ne. nauxi_exm ) call runend("exm_ohara_execute_timestep: vaulo_out has incorrect number of elements") 
      if ( size(viclo    , kind=ip) .ne. nviclo    ) call runend("exm_ohara_execute_timestep: viclo has incorrect number of elements") 

      !set voltage for the precalculation and precalculate Vinf and Taus
      call precalc % set_voltage( voltage )


      !% START OF CEIOHR CURRENT CALCULATIONS
      vaux1 = params % rgas  * params % temp / params % farad
      ena   = vaux1 * log(  params % nao / vcolo_in(5))
      ek    = vaux1 * log(  params % ko  / vcolo_in(3))
      eks   = vaux1 * log(( params % ko + params % pkna * params % nao ) / (vcolo_in(3) + params % pkna*vcolo_in(5)))

      !%convenient shorthand calculations
      vffrt = voltage * params % farad * params % farad/(params % rgas*params % temp)
      vfrt  = voltage * params % farad/(params % rgas*params % temp)

        !!! First calculate the auxiliaries, Currents, then concentrations and then calculate new voltage,

        !!! Update CaMKa
      vaux2 = 1.0_rp/(1.0_rp + (params % kmcam/vcolo_in(6)))
      viclo(26) = vaux2*params % camko*(1.0_rp - vcolo_in(11))
      viclo(22) = viclo(26) + vcolo_in(11)  ! camka

      !old-duplicate! viclo(26) = vaux2*params % camko*(1.0_rp - vcolo_in(11))
      !old-duplicate! viclo(22) = viclo(26) + vcolo_in(11)  ! camka


      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      !!!! CALCULATE GATES
      !call times(2) % ini()

      call precalc % calculate_gates(vaulo_in, vaulo_out, cell_type, params%ko, vcolo_in(6), params%dtimon, params%tauHL, params%kmn, params % k2n)
      !call times(2) % add()


        !!!!!!   calculate INA
      !calculate m gate
      !old! vinf            = exp(-(voltage + 39.57_rp)/9.871_rp)
      !old! vinf            = 1.0_rp/(1.0_rp + vinf)
      !old! vaux1           = 8.552_rp*exp(-(voltage + 77.42_rp)/5.955_rp)
      !old! xitaux          = vaux1 + (6.765_rp*exp((voltage + 11.64_rp)/34.77_rp))
      !old! tm              = 1.0_rp/xitaux
      !old! vaux0           = vaulo_in(1)
      !old! vaux3           = vinf - (vinf - vaux0)*exp(-params % dtimon/tm)
      !old! vaulo_out(1) = vaux3 ! value of variable m

        !!!!!!! hf
      !old! if ( model_type == EXMSLD_CELL_OHARA_INAPA ) then
      !old!    vinf   = exp((voltage + 78.50_rp)/6.22_rp) !modified from Passini et al. with LMS no reg(original: exp((voltage + 82.90_rp) / 6.086_rp))
      !old!    vinf   = 1.0_rp/(1.0_rp + vinf)
      !old!    vaux1  = exp(-(voltage + 3.8875_rp)/7.8579_rp) !modified from Passini et al. with LMS no reg(original: exp(-(voltage + 1.196_rp) / 6.285_rp))
      !old!    vaux2  = exp((voltage - 0.4963_rp)/9.1843_rp) !modified from Passini et al. with LMS no reg(original: exp((voltage + 0.5096_rp) / 20.27_rp))
      !old!    xitaux = (0.000003686_rp*vaux1) + (16.0_rp*vaux2) !modified from Passini et al. with LMS no reg(original: (0.00001432_rp * vaux1) + (6.149_rp * vaux2) )
      !old! else
      !old!    vinf   = exp((voltage + 82.90_rp)/6.086_rp) !original O'Hara-Rudy formulation
      !old!    vinf   = 1.0_rp/(1.0_rp + vinf)
      !old!    vaux1  = exp(-(voltage + 1.196_rp)/6.285_rp) !original O'Hara-Rudy formulation
      !old!    vaux2  = exp((voltage + 0.5096_rp)/20.27_rp) !origina O'Hara-Rudy formulation
      !old!    xitaux = (0.00001432_rp*vaux1) + (6.149_rp*vaux2) !original O'Hara-Rudy formulation
      !old! end if
      !old! xitaux          = 1.0_rp/xitaux
      !old! vaux0           = vaulo_in(2)
      !old! vaux3           = vinf - (vinf - vaux0)*exp(-params % dtimon/xitaux)
      !old! vaulo_out(2) = vaux3      ! value of variable hf
      
      !old! jss = vinf
      !old!   !!!! hs
      !old! vaux1           = exp(-(voltage + 17.95_rp)/28.05_rp)
      !old! vaux2           = exp((voltage + 5.730_rp)/56.66_rp)
      !old! xitaux          = (0.009794_rp*vaux1) + (0.3343_rp*vaux2)
      !old! ths             = 1.0_rp/xitaux
      !old! vaux0           = vaulo_in(3)
      !old! vaux3           = vinf - (vinf - vaux0)*exp(-params % dtimon/ths)
      !old! vaulo_out(3) = vaux3      ! value of variable hs

        !!!!!  params % j
      !old! if ( model_type == EXMSLD_CELL_OHARA_INAPA ) then
      !old!    vaux1  = exp(-(voltage + 116.7258_rp)/7.6005_rp) !modified from Passini et al. with LMS no reg(original: exp(-(v+100.6)/8.281))
      !old!    vaux2  = exp((voltage + 6.2719_rp)/9.0358_rp) !modified from Passini et al. with LMS no reg(original:  exp((voltage + 0.9941_rp) / 38.45_rp))
      !old!    xitaux = (0.8628_rp*vaux1) + (1.1096_rp*vaux2) !modified from Passini et al. with LMS no reg(original: (0.02136_rp * vaux1) + (0.3052_rp * vaux2))
      !old!    tj     = 4.8590_rp + (1.0_rp/xitaux) !modified from Passini et al. with LMS no reg(original: 2.038_rp + (1.0_rp / xitaux))
      !old! else
      !old!    vaux1  = exp(-(voltage + 100.6_rp)/8.281_rp) !original O'Hara-Rudy formulation
      !old!    vaux2  = exp((voltage + 0.9941_rp)/38.45_rp) !original O'Hara-Rudy formulation
      !old!    xitaux = (0.02136_rp*vaux1) + (0.3052_rp*vaux2) !original O'Hara-Rudy formulation
      !old!    tj     = 2.038_rp + (1.0_rp/xitaux) !original O'Hara-Rudy formulation
      !old! end if
      !old! vaux0           = vaulo_in(4)
      !old! vaux3           = jss - (jss - vaux0)*exp(-params % dtimon/tj)
      !old! vaulo_out(4) = vaux3      ! value of variable params % j

        !!!!  hsp
      !old! if ( model_type == EXMSLD_CELL_OHARA_INAPA ) then
      !old!    vinf = 1.0_rp/(1.0_rp + exp((voltage + 84.7_rp)/6.22_rp)) !modified from Passini et al. (original:   1.0_rp / (1.0_rp + exp((voltage + 89.1_rp) / 6.086_rp)))
      !old! else
      !old!    vinf = 1.0_rp/(1.0_rp + exp((voltage + 89.1_rp)/6.086_rp)) !original:   1.0_rp / (1.0_rp + exp((voltage + 89.1_rp) / 6.086_rp)))
      !old! end if
      !old! xitaux          = 3.0_rp*ths
      !old! vaux0           = vaulo_in(5)
      !old! vaux3           = vinf - (vinf - vaux0)*exp(-params % dtimon/xitaux)
      !old! vaulo_out(5) = vaux3      ! value of variable hsp

        !!! jp
      !old! xitaux          = 1.46_rp*tj
      !old! vaux0           = vaulo_in(6)
      !old! vaux3           = jss - (jss - vaux0)*exp(-params % dtimon/xitaux)
      !old! vaulo_out(6) = vaux3      ! value of variable jp

      !Calculate INA
             !!!!  hp
      hp       = (params % ahf*vaulo_out(2)) + (params % ahs*vaulo_out(5))
      vaux1    = params % ahf*vaulo_out(2)
      hh       = vaux1 + (params % ahs*vaulo_out(3))  !h
      finap    = 1.0_rp/(1.0_rp + (params % kmcamk/viclo(22)))
      vaux1    = params % gna*(voltage - ena)
      vaux1    = vaux1*vaulo_out(1)*vaulo_out(1)*vaulo_out(1)
      vaux2    = (1.0_rp - finap)*hh*vaulo_out(4)
      vaux2    = vaux2 + (finap*hp*vaulo_out(6))
      viclo(1) = vaux1*vaux2  !!!!  ina

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !%calculate inal
        !!!!!!  ml
      !old! vinf            = 1.0_rp + exp(-(voltage + 42.85_rp)/5.264_rp)
      !old! vinf            = 1.0_rp/vinf
      !old! vaux0           = vaulo_in(7)
      !old! vaux3           = vinf - (vinf - vaux0)*exp(-params % dtimon/tm)
      !old! vaulo_out(7) = vaux3      ! value of variable ml

        !!!!! hl
      !old! vinf = 1.0_rp + exp((voltage + 87.61_rp)/7.488_rp)
      !old! vinf = 1.0_rp/vinf
      !old! vaux0 = vaulo_in(8)
      !old! vaux3 = vinf - (vinf - vaux0)*exp(-params % dtimon/params % tauHL)
      !old! vaulo_out(8) = vaux3      ! value of variable hl

        !!!!  hlp
      !old! vinf            = 1.0_rp + exp((voltage + 93.81_rp)/7.488_rp)
      !old! vinf            = 1.0_rp/vinf
      !old! xitaux          = 3.0_rp*200.0_rp
      !old! vaux0           = vaulo_in(9)
      !old! vaux3           = vinf - (vinf - vaux0)*exp(-params % dtimon/xitaux)
      !old! vaulo_out(9) = vaux3      ! value of variable hlp

          !!!  Calculate inal current
      finalp   = 1.0_rp + (params % kmcamk/viclo(22))
      finalp   = 1.0_rp/finalp
      vaux1    = params % gnal*(voltage - ena)*vaulo_out(7)
      vaux2    = (1.0_rp - finalp)*vaulo_out(8) + (finalp*vaulo_out(9))
      viclo(2) = vaux1*vaux2

        !!!!! calculate ito
        !!!!!!  calculate variable a
      !old! vinf             = 1.0_rp + exp(-(voltage - 14.34_rp)/14.82_rp)
      !old! vinf             = 1.0_rp/vinf
      !old! vaux1            = 1.0_rp/(1.2089_rp*(1.0_rp + exp(-(voltage - 18.4099_rp)/29.3814_rp)))
      !old! vaux2            = 3.5_rp/(1.0_rp + exp((voltage + 100.0_rp)/29.3814_rp))
      !old! ta               = 1.0515_rp/(vaux1 + vaux2)
      !old! vaux0            = vaulo_in(10)
      !old! vaux3            = vinf - (vinf - vaux0)*exp(-params % dtimon/ta)
      !old! vaulo_out(10) = vaux3      ! value of variable a

        !!!!! calculate ifast  iss
      !old! vinf = 1.0_rp + exp((voltage + 43.94_rp)/5.711_rp)
      !old! iss  = 1.0_rp/vinf
      !old! if (cell_type == EXM_CELLTYPE_EPI) then
      !old!    delepi = (1.0_rp + exp((voltage + 70.0_rp)/5.0_rp))
      !old!    delepi = (0.95_rp/delepi)
      !old!    delepi = 1.0_rp - delepi
      !old! else
      !old!    delepi = 1.0_rp
      !old! end if
      !old! vaux1            = 0.08004_rp*exp((voltage + 50.0_rp)/16.59_rp) !tif
      !old! vaux2            = 0.3933_rp*exp(-(voltage + 100.0_rp)/100.0_rp)
      !old! vaux3            = vaux1 + vaux2
      !old! xitaux           = 4.562_rp + (1.0_rp/vaux3)
      !old! xitaux           = delepi*xitaux
      !old! tif              = xitaux
      !old! vaux0            = vaulo_in(11)
      !old! vaux3            = iss - (iss - vaux0)*exp(-params % dtimon/tif)
      !old! vaulo_out(11) = vaux3      ! value of variable ifast

        !!!!!!  calculate islow  tis
      !old! vaux1            = 0.001416_rp*exp(-(voltage + 96.52_rp)/59.05_rp)
      !old! vaux2            = 0.00000001780_rp*exp((voltage + 114.1_rp)/8.079_rp)
      !old! vaux3            = vaux1 + vaux2
      !old! xitaux           = 23.62_rp + (1.0_rp/vaux3)
      !old! xitaux           = delepi*xitaux
      !old! tis              = xitaux
      !old! vaux0            = vaulo_in(12)
      !old! vaux3            = iss - (iss - vaux0)*exp(-params % dtimon/tis)
      !old! vaulo_out(12) = vaux3      ! value of variable islow

        !!!!  calculate ap  (assp)
      !old! vinf             = 1.0_rp + exp(-(voltage - 24.34_rp)/14.82_rp)
      !old! vinf             = 1.0_rp/vinf
      !old! vaux0            = vaulo_in(13)
      !old! vaux3            = vinf - (vinf - vaux0)*exp(-params % dtimon/ta)
      !old! vaulo_out(13) = vaux3      ! value of variable ap

        !!!!!!!!!! calculate ifp (dti_develop)
      !old! vaux1 = exp((voltage - 167.4_rp)/15.89_rp)
      !old! vaux2 = exp(-(voltage - 12.23_rp)/0.2154_rp)
      !old! devel = vaux1 + vaux2
      !old! devel = 1.354_rp + (0.0001_rp/devel)

      !old! vaux1            = exp((voltage + 70.0_rp)/20.0_rp)
      !old! vaux1            = 1.0_rp + vaux1
      !old! recov            = 1.0_rp - (0.5_rp/vaux1)
      !old! xitaux           = devel*recov*tif
      !old! vaux0            = vaulo_in(14)
      !old! vaux3            = iss - (iss - vaux0)*exp(-params % dtimon/xitaux)
      !old! vaulo_out(14) = vaux3      ! value of variable ifp

        !!!!! calculate isp
      !old! xitaux           = devel*recov*tis
      !old! vaux0            = vaulo_in(15)
      !old! vaux3            = iss - (iss - vaux0)*exp(-params % dtimon/xitaux)
      !old! vaulo_out(15) = vaux3      ! value of variable isp

             !!!!!!  Ito current
             !!!! calculate aif, ais and ii
      !old! vaux1    =  exp((voltage - 213.6_rp)/151.2_rp)      
      !old! aif      = 1.0_rp/(1.0_rp + vaux1)
      aif      = precalc % get_inf(OHR_A_I_FAST)
      ais      = 1.0_rp - aif
      ii       = aif*vaulo_out(11) + (ais*vaulo_out(12))
      ipp      = aif*vaulo_out(14) + (ais*vaulo_out(15))
      vaux1    = params % kmcamk/viclo(22)
      fitop    = 1.0_rp/(1.0_rp + vaux1)
      vaux2    = (fitop*vaulo_out(13)*ipp)
      vaux1    = ((1.0_rp - fitop)*vaulo_out(10)*ii) + vaux2
      vaux1    = params % gto*(voltage - ek)*vaux1
      viclo(3) = vaux1        !!! ito

        !!! calculate ical, icana, icak
        !%calculate ical
        !!!!  calculate gate d
      !old! vaux1            = exp(-(voltage + 3.940_rp)/4.230_rp)
      !old! vinf             = 1.0_rp/(1.0_rp + vaux1)
      !old! vaux1            = exp(0.09_rp*(voltage + 14.0_rp))
      !old! vaux2            = exp(-0.05_rp*(voltage + 6.0_rp))
      !old! xitaux           = vaux1 + vaux2
      !old! xitaux           = 0.6_rp + (1.0_rp/xitaux)
      !old! vaux0            = vaulo_in(16)
      !old! vaux3            = vinf - (vinf - vaux0)*exp(-params % dtimon/xitaux)
      !old! vaulo_out(16) = vaux3      ! value of gate d

        !!!! calculate gate ff
      !old! vaux1            = exp((voltage + 19.58_rp)/3.696_rp)
      !old! vinf             = 1.0_rp/(1.0_rp + vaux1)
      !old! fss              = vinf
      !old! vaux1            = 0.0045_rp*exp((voltage + 20.0_rp)/10.0_rp)
      !old! vaux2            = 0.0045_rp*exp(-(voltage + 20.0_rp)/10.0_rp)
      !old! tff              = 7.0_rp + (1.0_rp/(vaux1 + vaux2))
      !old! vaux0            = vaulo_in(17)
      !old! vaux3            = fss - (fss - vaux0)*exp(-params % dtimon/tff)
      !old! vaulo_out(17) = vaux3      ! value of gate ff

        !!!!!!  calculate gate fs
      !old! vaux1            = 0.000035_rp*exp((voltage + 5.0_rp)/6.0_rp)
      !old! vaux2            = 0.000035_rp*exp(-(voltage + 5.0_rp)/4.0_rp)
      !old! xitaux           = 1000.0_rp + (1.0_rp/(vaux1 + vaux2))
      !old! vaux0            = vaulo_in(18)
      !old! vaux3            = fss - (fss - vaux0)*exp(-params % dtimon/xitaux)
      !old! vaulo_out(18) = vaux3      ! value of variable fs

        !!!!! calculate fcass
             !!fcass=fss
      !old! vaux1            = (0.04_rp*exp((voltage - 4.0_rp)/7.0_rp))
      !old! vaux2            = (0.04_rp*exp(-(voltage - 4.0_rp)/7.0_rp))
      !old! tfcaf            = 7.0_rp + (1.0_rp/(vaux1 + vaux2))
      !old! vaux0            = vaulo_in(19)
      !old! vaux3            = fss - (fss - vaux0)*exp(-params % dtimon/tfcaf)
      !old! vaulo_out(19) = vaux3      ! value of variable fcaf

        !!!! calculate gate fcas
      !old! vaux1            = (0.00012_rp*exp(voltage/7.0_rp))
      !old! vaux2            = (0.00012_rp*exp(-voltage/3.0_rp))
      !old! xitaux           = 100.0_rp + (1.0_rp/(vaux1 + vaux2))
      !old! vaux0            = vaulo_in(20)
      !old! vaux3            = fss - (fss - vaux0)*exp(-params % dtimon/xitaux)
      !old! vaulo_out(20) = vaux3      ! value of variable fcas

        !!! calculate gate jca
      !old! xitaux           = 75.0_rp
      !old! vaux0            = vaulo_in(21)
      !old! vaux3            = fss - (fss - vaux0)*exp(-params % dtimon/xitaux)
      !old! vaulo_out(21) = vaux3      ! value of variable jca

        !!! calculate gate ffp
      !old! xitaux           = 2.5_rp*tff
      !old! vaux0            = vaulo_in(23)
      !old! vaux3            = fss - (fss - vaux0)*exp(-params % dtimon/xitaux)
      !old! vaulo_out(23) = vaux3      ! value of variable ffp

        !!!! calculate gate fcafp
      !old! xitaux           = 2.5_rp*tfcaf
      !old! vaux0            = vaulo_in(24)
      !old! vaux3            = fss - (fss - vaux0)*exp(-params % dtimon/xitaux)
      !old! vaulo_out(24) = vaux3      ! value of variable ffp   

        !!!!! calculate icana and icak
      !old! km2n             = vaulo_out(21)*1.0_rp
      !old! vaux1            = params % kmn/vcolo_in(6)
      !old! vaux4            = (1.0_rp + vaux1)*(1.0_rp + vaux1)*(1.0_rp + vaux1)*(1.0_rp + vaux1)
      !old! vaux2            = (params % k2n/km2n) + vaux4
      !old! anca             = 1.0_rp/vaux2
      !old! dnca             = (anca*params % k2n/km2n)
      !old! vaux0            = vaulo_in(22)
      !old! vaux3            = dnca - (dnca - vaux0)*exp(-params % dtimon*km2n)
      !old! vaulo_out(22)    = vaux3

        !!!!!! calculate ff  = F
      ff    = (params % aff*vaulo_out(17)) + (params % afs*vaulo_out(18))   !value of variable ff
            !!! calculate fca
      !old! vaux1 = 1.0_rp + exp((voltage - 10.0_rp)/10.0_rp)
      !old! afcaf = 0.3_rp + (0.6_rp/vaux1)
      afcaf = precalc % get_inf(OHR_A_F_CA_FAST) 
      afcas = 1.0_rp - afcaf
      fca   = (afcaf*vaulo_out(19)) + (afcas*vaulo_out(20))
        !!!!  calculate fcap
      fcap  = afcaf*vaulo_out(24) + (afcas*vaulo_out(20))
        !!!!! calculate fp
      fp      = params % aff*vaulo_out(23) + (params % afs*vaulo_out(18))
      vaux2   = exp(2.0_rp*vfrt)
      vaux1   = exp(1.0_rp*vfrt)
      vaux3   = 1.0_rp/(vaux2 - 1.0_rp)
      phical  = 4.0_rp*vffrt*((vcolo_in(6)*vaux2) - (0.341_rp*params % cao))*vaux3
      vaux3   = 1.0_rp/(vaux1 - 1.0_rp)
      phicana = 1.0_rp*vffrt*((0.75_rp*vcolo_in(2)*vaux1) - (0.75_rp*params % nao))*vaux3
      phicak  = 1.0_rp*vffrt*((0.75_rp*vcolo_in(4)*vaux1) - (0.75_rp*params % ko))*vaux3
      ficalp  = 1.0_rp/(1.0_rp + (params % kmcamk/viclo(22)))

        !!!! calculate ical current
      vaux1 = (fp*(1.0_rp - vaulo_out(22)) + (vaulo_out(21)*fcap*vaulo_out(22)))
      vaux1 = ficalp*params % pcap*phical*vaulo_out(16)*vaux1
      vaux2 = ff*(1.0_rp - vaulo_out(22)) + (vaulo_out(21)*fca*vaulo_out(22))
      vaux3 = (1.0_rp - ficalp)*params % pca*phical*vaulo_out(16)
      viclo(4) = (vaux3*vaux2) + vaux1  !!! ical


        !!!! calculate icana current
      vaux1 = (fp*(1.0_rp - vaulo_out(22)) + (vaulo_out(21)*fcap*vaulo_out(22)))
      vaux1 = vaux1*ficalp*params % pcanap*phicana*vaulo_out(16)
      vaux2 = (ff*(1.0_rp - vaulo_out(22))) + (vaulo_out(21)*fca*vaulo_out(22))
      vaux2 = vaux2*(1.0_rp - ficalp)*params % pcana*phicana*vaulo_out(16)
      viclo(24) = vaux2 + vaux1  !! icana

        !!!!! calculate icak
      vaux1 = fp*(1.0_rp - vaulo_out(22)) + (vaulo_out(21)*fcap*vaulo_out(22))
      vaux1 = vaux1*ficalp*params % pcakp*phicak*vaulo_out(16)
      vaux2 = ff*(1.0_rp - vaulo_out(22)) + (vaulo_out(21)*fca*vaulo_out(22))
      vaux2 = vaux2*(1.0_rp - ficalp)*params % pcak*phicak*vaulo_out(16)
      viclo(25) = vaux2 + vaux1   !!! icak

        !!!!  calculate ikr GATES  XRF
      !old! vinf   = exp(-(voltage + 8.337_rp)/6.789_rp)
      !old! vinf   = 1.0_rp/(1.0_rp + vinf)
      !old! vaux1  = 0.00004123_rp*exp(-(voltage - 47.78_rp)/20.38_rp)
      !old! vaux2  = 0.3652_rp*exp((voltage - 31.66_rp)/3.869_rp)
      !old! xitaux = 12.98_rp + (1.0_rp/(vaux1 + vaux2))
      !old! vaux0  = vaulo_in(25)
      !old! vaux3  = vinf - (vinf - vaux0)*exp(-params % dtimon/xitaux)
      !old! vaulo_out(25) = vaux3      ! value of variable xrf
        !! GATE XRS
      !old! vaux1  = 0.06629_rp*exp((voltage - 34.70_rp)/7.355_rp)
      !old! vaux2  = 0.00001128_rp*exp(-(voltage - 29.74_rp)/25.94_rp)
      !old! xitaux = 1.865_rp + (1.0_rp/(vaux1 + vaux2))
      !old! vaux0  = vaulo_in(26)
      !old! vaux3  = vinf - (vinf - vaux0)*exp(-params % dtimon/xitaux)
      !old! vaulo_out(26) = vaux3      ! value of variable xrs

      !%calculate ikr CURRENT
      !old! vaux1 = exp((voltage + 54.81_rp)/38.21_rp)
      !old! axrf  = 1.0_rp/(1.0_rp + vaux1)
      axrf  = precalc % get_inf(OHR_A_XR_FAST) !1.0_rp/(1.0_rp + vaux1)
      axrs  = 1.0_rp - axrf
      xr    = axrf*vaulo_out(25) + (axrs*vaulo_out(26))
      !old! vaux1 = 1.0_rp + exp((voltage + 55.0_rp)/75.0_rp)
      !old! vaux2 = 1.0_rp + exp((voltage - 10.0_rp)/30.0_rp)
      !old! vaux3 = vaux1*vaux2
      !old! rkr   = 1.0_rp/vaux3
      rkr   = precalc % get_inf(OHR_R_KR_1) * precalc % get_inf(OHR_R_KR_2)
      vaux1 = sqrt(params % ko/5.4_rp)*(voltage - ek)
      viclo(5) = params % gkr*vaux1*xr*rkr   !!! ikr

        !!! %calculate iks GATES
      !old! vinf = exp(-(voltage + 11.60_rp)/8.932_rp)
      !old! vinf = 1.0_rp/(1.0_rp + vinf)
      !old! vaux1 = 0.0002326_rp*exp((voltage + 48.28_rp)/17.80_rp)
      !old! vaux2 = 0.001292_rp*exp(-(voltage + 210.0_rp)/230.0_rp)
      !old! xitaux = 817.3_rp + (1.0_rp/(vaux1 + vaux2))
      !old! vaux0 = vaulo_in(27)
      !old! vaux3 = vinf - (vinf - vaux0)*exp(-params % dtimon/xitaux)
      !old! vaulo_out(27) = vaux3      ! value of variable xs1

      !vinf = xs2ss = xs1ss
      !old! vaux1 = 0.01_rp*exp((voltage - 50.0_rp)/20.0_rp)
      !old! vaux2 = 0.0193_rp*exp(-(voltage + 66.54_rp)/31.0_rp)
      !old! xitaux = 1.0_rp/(vaux1 + vaux2)
      !old! vaux0 = vaulo_in(28)
      !old! vaux3 = vinf - (vinf - vaux0)*exp(-params % dtimon/xitaux)
      !old! vaulo_out(28) = vaux3      ! value of variable xs2

        !!!!  calculate gate xk1 GATES
      !old! vaux1 = 1.0_rp/(1.5692_rp*params % ko + 3.8115_rp)
      !old! vaux2 = voltage + 2.5538_rp*params % ko + 144.59_rp
      !old! vinf = exp(-vaux2*vaux1)
      !old! vinf = 1.0_rp/(1.0_rp + vinf)
      !old! vaux1 = exp(-(voltage + 127.2_rp)/20.36_rp)
      !old! vaux2 = exp((voltage + 236.8_rp)/69.33_rp)
      !old! xitaux = 122.2_rp/(vaux1 + vaux2)
      !old! vaux0 = vaulo_in(29)
      !old! vaux3 = vinf - (vinf - vaux0)*exp(-params % dtimon/xitaux)
      !old! vaulo_out(29) = vaux3      ! value of variable xk1

      !print *,abs(vaulo_out - vaulo_tmp )
      !stop 1

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!
      !!  Concentrations and everything else
      !!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         !!!! calculate iks CURRENTS
      ksca = (0.000038_rp/vcolo_in(1))**(1.4_rp) !!(7.0_rp/5.0_rp)
      ksca = 1.0_rp + (0.6_rp/(1.0_rp + ksca))
      vaux1 = params % gks*ksca*vaulo_out(27)*vaulo_out(28)
      viclo(6) = vaux1*(voltage - eks)


             !!!! calculate ik1 CURRENT
      rk1 = 1.0_rp + safe_exp_cond((voltage + 105.8_rp - 2.6_rp*params % ko)/9.493_rp, precalc%get_use_safe_exp())
      rk1 = 1.0_rp/rk1
      vaux1 = params % gk1*sqrt(params % ko)*rk1*vaulo_out(29)
      viclo(7) = vaux1*(voltage - ek) !!! ik1
      
      !!!! calculate IkATP CURRENT for ISCHEMIA modelling
      ikatp=params % gkatp*params % fkatp *(params % ko/5.4_rp)**0.24_rp * (voltage-ek)
      viclo(23)=ikatp
      
      !%calculate INaCa CURRENT
      hca    = exp(params % qca*vfrt)
      hna    = exp(params % qna*vfrt)
      h1     = 1.0_rp + (vcolo_in(5)/params % kna3)*(1.0_rp + hna)
      h2     = (vcolo_in(5)*hna)/(params % kna3*h1)
      h3     = 1.0_rp/h1
      h4     = 1.0_rp + vcolo_in(5)/params % kna1*(1.0_rp + (vcolo_in(5)/params % kna2))
      h5     = vcolo_in(5)*vcolo_in(5)/(h4*params % kna1*params % kna2)
      h6     = 1.0_rp/h4
      h7     = 1.0_rp + params % nao/params % kna3*(1.0_rp + (1.0_rp/hna))
      h8     = params % nao/(params % kna3*hna*h7)
      h9     = 1.0_rp/h7
      h10    = params % kasymm + 1.0_rp + params % nao/params % kna1*(1.0_rp + (params % nao/params % kna2))
      h11    = params % nao*params % nao/(h10*params % kna1*params % kna2)
      h12    = 1.0_rp/h10
      k1     = h12*params % cao*params % kcaon
      k2     = params % kcaoff
      k3p    = h9*params % wca
      k3pp   = h8*params % wnaca
      k3     = k3p + k3pp
      k4p    = h3*params % wca/hca
      k4pp   = h2*params % wnaca
      k4     = k4p + k4pp
      k5     = params % kcaoff
      k6     = h6*vcolo_in(1)*params % kcaon
      k7     = h5*h2*params % wna
      k8     = h8*h11*params % wna
      x1     = (k2*k4*(k7 + k6)) + (k5*k7*(k2 + k3))
      x2     = (k1*k7*(k4 + k5)) + (k4*k6*(k1 + k8))
      x3     = (k1*k3*(k7 + k6)) + (k8*k6*(k2 + k3))
      x4     = (k2*k8*(k4 + k5)) + (k3*k5*(k1 + k8))
      e1     = x1/(x1 + x2 + x3 + x4)
      e2     = x2/(x1 + x2 + x3 + x4)
      e3     = x3/(x1 + x2 + x3 + x4)
      e4     = x4/(x1 + x2 + x3 + x4)
      allo   = (params % kmcaact/vcolo_in(1))*(params % kmcaact/vcolo_in(1))
      allo   = 1.0_rp/(1.0_rp + allo)
      zna    = 1.0_rp
      jncxna = (3.0_rp*(e4*k7 - e1*k8)) + (e3*k4pp) - (e2*k3pp)
      jncxca = (e2*k2) - (e1*k1)
      vaux1  = (zna*jncxna) + (params % zca*jncxca)
      viclo(8) = 0.8_rp*params % gncx*allo*vaux1  !!! inaca_i Current

      !%calculate inaca_ss
      h1     = 1.0_rp + vcolo_in(2)/params % kna3*(1.0_rp + hna)
      h2     = (vcolo_in(2)*hna)/(params % kna3*h1)
      h3     = 1.0_rp/h1
      h4     = 1.0_rp + vcolo_in(2)/params % kna1*(1.0_rp + (vcolo_in(2)/params % kna2))
      h5     = vcolo_in(2)*vcolo_in(2)/(h4*params % kna1*params % kna2)
      h6     = 1.0_rp/h4
      h7     = 1.0_rp + (params % nao/params % kna3)*(1.0_rp + (1.0_rp/hna))
      h8     = params % nao/(params % kna3*hna*h7)
      h9     = 1.0_rp/h7
      h10    = params % kasymm + 1.0_rp + params % nao/params % kna1*(1.0_rp + (params % nao/params % kna2))
      h11    = params % nao*params % nao/(h10*params % kna1*params % kna2)
      h12    = 1.0_rp/h10
      k1     = h12*params % cao*params % kcaon
      k2     = params % kcaoff
      k3p    = h9*params % wca
      k3pp   = h8*params % wnaca
      k3     = k3p + k3pp
      k4p    = h3*params % wca/hca
      k4pp   = h2*params % wnaca
      k4     = k4p + k4pp
      k5     = params % kcaoff
      k6     = h6*vcolo_in(6)*params % kcaon
      k7     = h5*h2*params % wna
      k8     = h8*h11*params % wna
      x1     = (k2*k4*(k7 + k6)) + (k5*k7*(k2 + k3))
      x2     = (k1*k7*(k4 + k5)) + (k4*k6*(k1 + k8))
      x3     = (k1*k3*(k7 + k6)) + (k8*k6*(k2 + k3))
      x4     = (k2*k8*(k4 + k5)) + (k3*k5*(k1 + k8))
      e1     = x1/(x1 + x2 + x3 + x4)
      e2     = x2/(x1 + x2 + x3 + x4)
      e3     = x3/(x1 + x2 + x3 + x4)
      e4     = x4/(x1 + x2 + x3 + x4)
      allo   = (params % kmcaact/vcolo_in(6))*(params % kmcaact/vcolo_in(6))
      allo   = 1.0_rp/(1.0_rp + allo)
      jncxna = (3.0_rp*(e4*k7 - e1*k8)) + (e3*k4pp) - (e2*k3pp)
      jncxca = (e2*k2) - (e1*k1)
      viclo(9) = 0.2_rp*params % gncx*allo*((zna*jncxna) + (params % zca*jncxca))   !!! inaca_ss

      INaCa = viclo(8) + viclo(9)

      !%calculate inak
      k3p    = 1899.0_rp
      k4p    = 639.0_rp
      knai   = params % knai0*exp(params % delta*vfrt/3.0_rp)
      knao   = params % knao0*exp((1.0_rp - params % delta)*vfrt/3.0_rp)
      pp     = params % ep/(1.0_rp + params % hbig/params % khp + (vcolo_in(5)/params % knap) + (vcolo_in(3)/params % kxkur))
      vaux1  = (vcolo_in(5)/knai)*(vcolo_in(5)/knai)*(vcolo_in(5)/knai)
      vaux2  = (1.0_rp + (vcolo_in(5)/knai))*(1.0_rp + (vcolo_in(5)/knai))*(1.0_rp + (vcolo_in(5)/knai))
      vaux3  = (1.0_rp + (vcolo_in(3)/params % kki))*(1.0_rp + (vcolo_in(3)/params % kki))
      a1     = (params % k1p*vaux1)/(vaux2 + vaux3 - 1.0_rp)
      b1     = params % k1m*params % mgadp
      a2     = params % k2p
      vaux1  = (params % nao/knao)*(params % nao/knao)*(params % nao/knao)
      vaux2  = (1.0_rp + (params % nao/knao))*(1.0_rp + (params % nao/knao))*(1.0_rp + (params % nao/knao))
      vaux3  = (1.0_rp + (params % ko/params % kko))*(1.0_rp + (params % ko/params % kko))
      b2     = (params % k2m*vaux1)/(vaux2 + vaux3 - 1.0_rp)
      vaux1  = (params % ko/params % kko)*(params % ko/params % kko)
      a3     = (k3p*vaux1)/(vaux2 + vaux3 - 1.0_rp)
      b3     = (params % k3m*pp*params % hbig)/(1.0_rp + (params % mgatp/params % kmgatp))
      a4     = (k4p*params % mgatp/params % kmgatp)/(1.0_rp + (params % mgatp/params % kmgatp))
      vaux1  = (vcolo_in(3)/params % kki)*(vcolo_in(3)/params % kki)
      vaux2  = (1.0_rp + (vcolo_in(5)/knai))*(1.0_rp + (vcolo_in(5)/knai))*(1.0_rp + (vcolo_in(5)/knai))
      vaux3  = (1.0_rp + (vcolo_in(3)/params % kki))*(1.0_rp + (vcolo_in(3)/params % kki))
      b4     = (params % k4m*vaux1)/(vaux2 + vaux3 - 1.0_rp)
      x1     = (a4*a1*a2) + (b2*b4*b3) + (a2*b4*b3) + (b3*a1*a2)
      x2     = (b2*b1*b4) + (a1*a2*a3) + (a3*b1*b4) + (a2*a3*b4)
      x3     = (a2*a3*a4) + (b3*b2*b1) + (b2*b1*a4) + (a3*a4*b1)
      x4     = (b4*b3*b2) + (a3*a4*a1) + (b2*a4*a1) + (b3*b2*a1)
      e1     = x1/(x1 + x2 + x3 + x4)
      e2     = x2/(x1 + x2 + x3 + x4)
      e3     = x3/(x1 + x2 + x3 + x4)
      e4     = x4/(x1 + x2 + x3 + x4)
      zk     = 1.0_rp
      jnakna = 3.0_rp*((e1*a3) - (e2*b3))
      jnakk  = 2.0_rp*((e4*b1) - (e3*a1))
      viclo(10) = params % pnak*(zna*jnakna + zk*jnakk)  !!!! inak

      !%calculate ikb CURRENT
      !old! xkb = 1.0_rp/(1.0_rp + exp(-(voltage - 14.48_rp)/18.34_rp))
      xkb = precalc % get_inf(OHR_X_KB)
      viclo(11) = params % gkb*xkb*(voltage - ek)  !!!! ikb

      !%calculate inab CURRENT
      vaux1 = params % pnab*vffrt*(vcolo_in(5)*exp(vfrt) - params % nao)
      viclo(12) = vaux1/(exp(vfrt) - 1.0_rp)

      !%calculate icab CURRENT
      vaux1 = (vcolo_in(1)*exp(2.0_rp*vfrt) - (0.341_rp*params % cao))/(exp(2.0_rp*vfrt) - 1.0_rp)
      viclo(13) = params % pcab*4.0_rp*vffrt*vaux1   !!!!icab

      !%calculate ipca CURRENT
      viclo(14) = params % gpca*vcolo_in(1)/(0.0005_rp + vcolo_in(1))  !!! ipca

      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           !!  VOLTAGE INTEGRATION
      v1 = 0.0_rp
      v2 = 0.0_rp
      do i = 1, 14
         v1 = v1 + viclo(i)
      end do
      do i = 23, 25
         v2 = v2 + viclo(i)
      end do
      do i=2,7 
        output_params % Inet = output_params % Inet + viclo(i)
      end do
      output_params % xioni = v1+v2    !+sac

      output_params % Inet = output_params % Inet + viclo(23)
      output_params % qnet = output_params % qnet + (output_params % Inet * params % dtimon)

      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !!! Update CaMKa
      !old-duplicate! vaux2 = 1.0_rp/(1.0_rp + ( params % kmcam/vcolo_in(6)))
      !old-duplicate! viclo(26) = vaux2*params % camko*(1.0_rp - vcolo_in(11))
      !old-duplicate! viclo(22) = viclo(26) + vcolo_in(11)  ! camka

      !%update camk
      rhsx1 = params % acamk*viclo(26)*(viclo(26) + vcolo_in(11))
      rhsx = rhsx1 - (params % bcamk*vcolo_in(11))
      val0 = vcolo_in(11)
      camkt = val0 + params % dtimon*rhsx
      !k1 = rhsx
      !k2 = rhsx + 0.5_rp * params % dtimon* k1
      !k3 = rhsx + 0.5_rp * params % dtimon * k2
      !k4 = rhsx + params % dtimon * k3
      !camkt = val0 + ((params % dtimon / (6.0_rp)) *(k1 + 2.0_rp * k2 + 2.0_rp * k3 + k4))
      vcolo_out(11) = camkt      ! value of camkt concentration

      !%calculate diffusion fluxes
      viclo(16) = (vcolo_in(2) - vcolo_in(5))/2.0_rp    !!! jdiffna
      viclo(17) = (vcolo_in(4) - vcolo_in(3))/2.0_rp    !!! jdiffk
      viclo(15) = (vcolo_in(6) - vcolo_in(1))/0.2_rp    !!! jdiff
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      !%calculate ryanodione receptor calcium induced calcium release from the jsr
      bt = 4.75_rp
      a_rel = 0.5_rp*bt
      vaux1 = (1.5_rp/vcolo_in(8))**8.0_rp
      vinf = a_rel*(-viclo(4))/(1.0_rp + vaux1)

      if     (cell_type == EXM_CELLTYPE_MID) then
         vinf = vinf * 1.7_rp * material_params(2, 16)
      elseif (cell_type == EXM_CELLTYPE_ENDO) then
         vinf = vinf          * material_params(1, 16)
      elseif (cell_type == EXM_CELLTYPE_EPI) then
         vinf = vinf          * material_params(3, 16)
      end if
      xitaux = bt/(1.0_rp + 0.0123_rp/vcolo_in(8))

      if (xitaux < 0.001_rp) then  !fixed bug was wrongly 0.005
         xitaux = 0.001_rp    !fixed bug was wrongly 0.005
      end if

      vaux0 = vcolo_in(9)
      !vaux1 = (vinf - vaux0) * xitaux;
      !vaux2 = (vinf - (vaux0 + 0.5_rp * params % dtimon * vaux1)) * xitaux;
      !vaux3 = (vinf - (vaux0 + 0.5_rp * params % dtimon * vaux2)) * xitaux;
      !vaux4 = (vinf - (vaux0 + params % dtimon * vaux3)) * xitaux;
      !vaux5 = vaux0 + (params % dtimon / 6.0_rp) * (vaux1 + 2.0_rp * vaux2 + 2.0_rp * vaux3 + vaux4);
      !vaux3 = vaux5
      vaux3 = vinf - (vinf - vaux0)*safe_exp_cond(-params % dtimon/xitaux, precalc%get_use_safe_exp())
      vcolo_out(9) = vaux3      ! value of variable jrelnp


             !!  JRELP
      btp = 1.25_rp*bt
      a_relp = 0.5_rp*btp
      vaux1 = (1.5_rp/vcolo_in(8))**8.0_rp
      vinf = a_relp*(-viclo(4))/(1.0_rp + vaux1)
      if (cell_type == EXM_CELLTYPE_MID) then
         vinf = vinf*1.7_rp
      end if
      xitaux = btp/(1.0_rp + (0.0123_rp/vcolo_in(8)))

      if (xitaux < 0.001_rp) then    !fixed bug was wrongly 0.005
         xitaux = 0.001_rp   !fixed bug was wrongly 0.005
      end if
      vaux0 = vcolo_in(10)
      !vaux1 = (vinf - vaux0) * xitaux;
      !vaux2 = (vinf - (vaux0 + 0.5_rp * params % dtimon * vaux1)) * xitaux;
      !vaux3 = (vinf - (vaux0 + 0.5_rp * params % dtimon * vaux2)) * xitaux;
      !vaux4 = (vinf - (vaux0 + params % dtimon * vaux3)) * xitaux;
      !vaux5 = vaux0 + (params % dtimon / 6.0_rp) * (vaux1 + 2.0_rp * vaux2 + 2.0_rp * vaux3 + vaux4);
      !vaux3 = vaux5
      vaux3 = vinf - (vinf - vaux0)*safe_exp_cond(-params % dtimon/xitaux,  precalc%get_use_safe_exp())
      vcolo_out(10) = vaux3      ! value of jrelp

      fjrelp = 1.0_rp/(1.0_rp + (params % kmcamk/viclo(22)))
      viclo(21) = (1.0_rp - fjrelp)*vcolo_out(9) + (fjrelp*vcolo_out(10)) !!! jrel

      !%calculate serca pump, ca uptake flux
      jupnp = (0.004375_rp*vcolo_in(1)/(vcolo_in(1) + 0.00092_rp))
      vaux1 = 1.0_rp/(vcolo_in(1) + 0.00092_rp - 0.00017_rp)
      jupp = 2.75_rp*0.004375_rp*vcolo_in(1)*vaux1
      if      (cell_type == EXM_CELLTYPE_EPI) then
         jupnp   = jupnp * 1.3_rp * material_params(3, 12)
         jupp    = jupp  * 1.3_rp * material_params(3, 12)
         HFjleak =                  material_params(3, 15)
      else if (cell_type == EXM_CELLTYPE_ENDO) then  !fixed this input 
         jupnp   = jupnp * material_params(1, 12)
         jupp    = jupp  * material_params(1, 12)
         HFjleak =         material_params(1, 15)
      else if (cell_type == EXM_CELLTYPE_MID) then
         jupnp   = jupnp * material_params(2, 12)
         jupp    = jupp  * material_params(2, 12)
         HFjleak =         material_params(2, 15)
      end if
      fjupp = (1.0_rp/(1.0_rp + params % kmcamk/viclo(22)))
      viclo(19) = (0.0039375_rp*vcolo_in(7)/15.0_rp )* HFjleak  !!!! jleak
      viclo(18) = (1.0_rp - fjupp)*jupnp + fjupp*jupp - viclo(19)  !!! jup

      !%calculate tranlocation flux
      viclo(20) = (vcolo_in(7) - vcolo_in(8))/100.0_rp  !!!! jtr

                     !!!!!  CONCENTRATIONS!!!!

      !%update intracellular concentrations, using buffers for cai, cass, cajsr
             !! calculate na current
      vaux1 = viclo(1) + viclo(2) + viclo(12)
      vaux2 = 3.0_rp*viclo(8) + 3.0_rp*viclo(10)
      vaux3 = -(vaux1 + vaux2)*params % acap/(params % farad*params % vmyo)
      rhsx = vaux3 + (viclo(16)*params % vss/params % vmyo)
      val0 = vcolo_in(5)
      nai = val0 + params % dtimon*rhsx
      !k1 = rhsx
      !k2 = rhsx + 0.5_rp * params % dtimon* k1
      !k3 = rhsx + 0.5_rp * params % dtimon * k2
      !k4 = rhsx + params % dtimon * k3
      !nai = val0 + (params % dtimon / 6.0_rp) *(k1 + 2.0_rp * k2 + 2.0_rp * k3 + k4)
      vcolo_out(5) = nai


             !!! calculate na current in subspace ss
      vaux1 = (viclo(24) + 3.0_rp*viclo(9))*params % acap/(params % farad*params % vss)
      rhsx = -vaux1 - viclo(16)
      val0 = vcolo_in(2)
      nass = val0 + params % dtimon*rhsx
      !k1 = rhsx
      !k2 = rhsx + 0.5_rp * params % dtimon* k1
      !k3 = rhsx + 0.5_rp * params % dtimon * k2
      !k4 = rhsx + params % dtimon * k3
      !nass = val0 + (params % dtimon / 6.0_rp) *(k1 + 2.0_rp * k2 + 2.0_rp * k3 + k4)
      vcolo_out(2) = nass             !!! nass

              !!! calculate k current
      vaux1 = viclo(3) + viclo(5) + viclo(6) + viclo(23)
      vaux2 = viclo(7) + viclo(11) - (2.0_rp*viclo(10))
      vaux3 = (viclo(17)*params % vss/params % vmyo)
      rhsx = -((vaux1 + vaux2)*params % acap/(params % farad*params % vmyo)) + vaux3
      val0 = vcolo_in(3)
      ki = val0 + params % dtimon*rhsx
      !k1 = rhsx
      !k2 = rhsx + 0.5_rp * params % dtimon* k1
      !k3 = rhsx + 0.5_rp * params % dtimon * k2
      !k4 = rhsx + params % dtimon * k3
      !ki = val0 + (params % dtimon / 6.0_rp) *(k1 + 2.0_rp * k2 + 2.0_rp * k3 + k4)
      vcolo_out(3) = ki             !!! ki

             !!!!  calculate k current in the subspace ss
      rhsx = -(viclo(25)*params % acap/(params % farad*params % vss)) - viclo(17)
      val0 = vcolo_in(4)
      kss = val0 + params % dtimon*rhsx
      !k1 = rhsx
      !k2 = rhsx + 0.5_rp * params % dtimon* k1
      !k3 = rhsx + 0.5_rp * params % dtimon * k2
      !k4 = rhsx + params % dtimon * k3
      !kss = val0 + (params % dtimon / 6.0_rp) *(k1 + 2.0_rp * k2 + 2.0_rp * k3 + k4)
      vcolo_out(4) = kss             !!! kss

      !!!  calculate ca current (Cai)

      if( .not. has_land ) then
          vaux1 = (params % kmcmdn + vcolo_in(1))*(params % kmcmdn + vcolo_in(1))
          vaux2 = (params % kmtrpn + vcolo_in(1))*(params % kmtrpn + vcolo_in(1))
          vaux3 = 1.0_rp/(1.0_rp + (params % cmdnmax*params % kmcmdn/vaux1) + (params % trpnmax*params % kmtrpn/vaux2))
          rhsx1 = viclo(14) + viclo(13) - (2.0_rp*viclo(8))
          rhsx2 = -(viclo(18)*params % vnsr/params % vmyo) + (viclo(15)*params % vss/params % vmyo)

          rhsx = vaux3*(-(rhsx1*params % acap/(2.0_rp*params % farad*params % vmyo)) + rhsx2)
          val0 = vcolo_in(1)
          cai = val0 + params % dtimon*rhsx
          !k1 = rhsx
          !k2 = rhsx + 0.5_rp * params % dtimon* k1
          !k3 = rhsx + 0.5_rp * params % dtimon * k2
          !k4 = rhsx + params % dtimon * k3
          !cai = val0 + (params % dtimon / 6.0_rp) *(k1 + 2.0_rp * k2 + 2.0_rp * k3 + k4)
          vcolo_out(1) = cai             !!! cai
      else
          vcolo_out(1) = -1.0_rp !LAND part is outside, make sure Cai gets calculated
      end if
      !moved_out! else (if has_land)
      !moved_out!    vaux1 = (params % kmcmdn + vcolo_in(1))*(params % kmcmdn + vcolo_in(1))
      !moved_out!    vaux3 = 1.0_rp/(1.0_rp + (params % cmdnmax*params % kmcmdn/vaux1))
      !moved_out!    rhsx1 = viclo(14) + viclo(13) - (2.0_rp*viclo(8))
      !moved_out!    rhsx2 = -(viclo(18)*params % vnsr/params % vmyo) + (viclo(15)*params % vss/params % vmyo) - (params % k_TRPN*(((vcolo_in(1)*1000.0_rp / params % Ca50)**params % n_TRPN)*(1.0_rp-params % CaTRPN)-params % CaTRPN))*params % trpnmax
      !moved_out! end if
      !moved_out!    
      !moved_out! rhsx = vaux3*(-(rhsx1*params % acap/(2.0_rp*params % farad*params % vmyo)) + rhsx2)
      !moved_out! val0 = vcolo_in(1)
      !moved_out! cai = val0 + params % dtimon*rhsx
      !moved_out! !k1 = rhsx
      !moved_out! !k2 = rhsx + 0.5_rp * params % dtimon* k1
      !moved_out! !k3 = rhsx + 0.5_rp * params % dtimon * k2
      !moved_out! !k4 = rhsx + params % dtimon * k3
      !moved_out! !cai = val0 + (params % dtimon / 6.0_rp) *(k1 + 2.0_rp * k2 + 2.0_rp * k3 + k4)
      !moved_out! vcolo_out(1) = cai             !!! cai

             !!!! calculate ca current in the subspace ss
      vaux1 = (params % kmbsr + vcolo_in(6))*(params % kmbsr + vcolo_in(6))
      vaux2 = (params % kmbsl + vcolo_in(6))*(params % kmbsl + vcolo_in(6))
      vaux3 = 1.0_rp/(1.0_rp + (params % bsrmax*params % kmbsr/vaux1) + (params % bslmax*params % kmbsl/vaux2)) !Bcass
      rhsx1 = viclo(4) - (2.0_rp*viclo(9))
      rhsx2 = viclo(21)*params % vjsr/params % vss
      rhsx = vaux3*(-rhsx1*params % acap/(2.0_rp*params % farad*params % vss) + rhsx2 - viclo(15))
      val0 = vcolo_in(6)
      cass = val0 + params % dtimon*rhsx
      !k1 = rhsx
      !k2 = rhsx + 0.5_rp * params % dtimon* k1
      !k3 = rhsx + 0.5_rp * params % dtimon * k2
      !k4 = rhsx + params % dtimon * k3
      !cass = val0 + (params % dtimon / 6.0_rp) *(k1 + 2.0_rp * k2 + 2.0_rp * k3 + k4)
      vcolo_out(6) = cass             !!! cass

             !!!! calculate ca current in the sarcoplasmic reticulum nsr
      rhsx1 = viclo(20)*params % vjsr/params % vnsr
      rhsx = viclo(18) - rhsx1
      val0 = vcolo_in(7)
      cansr = val0 + params % dtimon*rhsx
      !k1 = rhsx
      !k2 = rhsx + 0.5_rp * params % dtimon* k1
      !k3 = rhsx + 0.5_rp * params % dtimon * k2
      !k4 = rhsx + params % dtimon * k3
      !cansr = val0 + (params % dtimon / 6.0_rp) *(k1 + 2.0_rp * k2 + 2.0_rp * k3 + k4)
      vcolo_out(7) = cansr             !!! cansr

             !!!! calculate ca current in the junctional sarcoplasmic reticulum jsr
      vaux1 = (params % kmcsqn + vcolo_in(8))*(params % kmcsqn + vcolo_in(8))
      vaux3 = params % csqnmax*params % kmcsqn/vaux1
      vaux2 = 1.0_rp/(1.0_rp + vaux3)
      rhsx = vaux2*(viclo(20) - viclo(21))
      val0 = vcolo_in(8)
      cajsr = val0 + params % dtimon*rhsx
      !k1 = rhsx
      !k2 = rhsx + 0.5_rp * params % dtimon* k1
      !k3 = rhsx + 0.5_rp * params % dtimon * k2
      !k4 = rhsx + params % dtimon * k3
      !cajsr = val0 + ((params % dtimon / (6.0_rp)) *(k1 + 2.0_rp * k2 + 2.0_rp * k3 + k4))
      vcolo_out(8) = cajsr             !!! cajsr
   
   contains 
#include "exm_safe_exp.f90.inc"


   end subroutine exm_ohara_execute_timestep



   subroutine exm_init_voltages(condition_type, mat, celltype, outputs)
      implicit none

      integer(ip), intent(in) :: condition_type, mat, celltype
      TYPE(CELL_MODEL_OUTPUTS), intent(inout)        :: outputs


      select case (condition_type)

!      case (EXM_CELL_NORMAL_HUMAN_FEMALE_600BCL)
!         select case (celltype)
!         case (EXM_CELLTYPE_ENDO) !1     
!            voltage_local     =   -87.664902277317367_rp
!            viclo_exm(1:26)   = 0.0_rp
!            vcolo_exm(1, 1:2) =    1.1455758338192559E-004_rp
!            vcolo_exm(2, 1:2) =    8.8679915666602636_rp
!            vcolo_exm(3, 1:2) =    142.60572798713852_rp
!            vcolo_exm(4, 1:2) =    142.60568727724311_rp
!            vcolo_exm(5, 1:2) =    8.8678427399199471_rp
!            vcolo_exm(6, 1:2) =    1.1587317953261430E-004_rp
!            vcolo_exm(7, 1:2) =    3.1369878324189133_rp
!            vcolo_exm(8, 1:2) =    2.7944805944253157_rp
!            vcolo_exm(9, 1:2) =    7.0702906508163200E-007_rp
!            vcolo_exm(10, 1:2) =    8.8375189003228634E-007_rp
!            vcolo_exm(11, 1:2) =    7.4513326488172979E-002_rp
!            !vcolo_exm(12, 1:2) =    0.0000000000000000_rp
!            vaulo_exm(1, 1:2) =    7.5972512633682120E-003_rp
!            vaulo_exm(2, 1:2) =   0.81351764109286584_rp
!            vaulo_exm(3, 1:2) =   0.81353058798150524_rp
!            vaulo_exm(4, 1:2) =   0.80819140775794773_rp
!            vaulo_exm(5, 1:2) =   0.61670107627034498_rp
!            vaulo_exm(6, 1:2) =   0.78906362841872835_rp
!            vaulo_exm(7, 1:2) =    2.0070623504288491E-004_rp
!            vaulo_exm(8, 1:2) =   0.41228785270410750_rp
!            vaulo_exm(9, 1:2) =   0.19128492327140989_rp
!            vaulo_exm(10, 1:2) =    1.0241152945904416E-003_rp
!            vaulo_exm(11, 1:2) =   0.99952689235110237_rp
!            vaulo_exm(12, 1:2) =   0.32466518417485002_rp
!            vaulo_exm(13, 1:2) =    5.2182028268431757E-004_rp
!            vaulo_exm(14, 1:2) =   0.99952692611944138_rp
!            vaulo_exm(15, 1:2) =   0.36076385857210003_rp
!            vaulo_exm(16, 1:2) =    2.5350254970419230E-009_rp
!            vaulo_exm(17, 1:2) =   0.99999998999729311_rp
!            vaulo_exm(18, 1:2) =   0.77728679107429577_rp
!            vaulo_exm(19, 1:2) =   0.99999998999759832_rp
!            vaulo_exm(20, 1:2) =   0.98697365138801163_rp
!            vaulo_exm(21, 1:2) =   0.99421320450114170_rp
!            vaulo_exm(22, 1:2) =    9.0113765701275959E-003_rp
!            vaulo_exm(23, 1:2) =   0.99999997975008115_rp
!            vaulo_exm(24, 1:2) =   0.99999998956811254_rp
!            vaulo_exm(25, 1:2) =    7.3225421282804371E-004_rp
!            vaulo_exm(26, 1:2) =   0.70374586064651934_rp
!            vaulo_exm(27, 1:2) =   0.46041257143396341_rp
!            vaulo_exm(28, 1:2) =    2.0270889409482307E-004_rp
!            vaulo_exm(29, 1:2) =   0.99684836820524714_rp
!            !vaulo_exm(30, 1:2) =    0.0000000000000000_rp         
!         case (EXM_CELLTYPE_MID)  !2 
!            voltage_local =   -87.435195284762415_rp
!            viclo_exm(1:26) = 0.0_rp
!            vcolo_exm(1, 1:2) =    1.1890657170525055E-004_rp
!            vcolo_exm(2, 1:2) =    9.7706696760397378_rp
!            vcolo_exm(3, 1:2) =    141.66483952218840_rp
!            vcolo_exm(4, 1:2) =    141.66479112719344_rp
!            vcolo_exm(5, 1:2) =    9.7704648272293344_rp
!            vcolo_exm(6, 1:2) =    1.1898783498170927E-004_rp
!            vcolo_exm(7, 1:2) =    2.8021201980086876_rp
!            vcolo_exm(8, 1:2) =    2.4236072150570278_rp
!            vcolo_exm(9, 1:2) =    1.4071451498516887E-006_rp
!            vcolo_exm(10, 1:2) =    1.7586168867274698E-006_rp
!            vcolo_exm(11, 1:2) =    7.8962722885925360E-002_rp
!            !vcolo_exm(12, 1:2) =    0.0000000000000000_rp
!            vaulo_exm(1, 1:2) =    7.7747291376501189E-003_rp
!            vaulo_exm(2, 1:2) =   0.80784540409976191_rp
!            vaulo_exm(3, 1:2) =   0.80786105547820475_rp
!            vaulo_exm(4, 1:2) =   0.80265157891836159_rp
!            vaulo_exm(5, 1:2) =   0.60792453114967615_rp
!            vaulo_exm(6, 1:2) =   0.78403552237781493_rp
!            vaulo_exm(7, 1:2) =    2.0965657866126548E-004_rp
!            vaulo_exm(8, 1:2) =   0.41110392493267284_rp
!            vaulo_exm(9, 1:2) =   0.19178946886523357_rp
!            vaulo_exm(10, 1:2) =    1.0400968637168772E-003_rp
!            vaulo_exm(11, 1:2) =   0.99950747422536579_rp
!            vaulo_exm(12, 1:2) =   0.32788959685532149_rp
!            vaulo_exm(13, 1:2) =    5.2996757393124953E-004_rp
!            vaulo_exm(14, 1:2) =   0.99950751019835571_rp
!            vaulo_exm(15, 1:2) =   0.36509196261881782_rp
!            vaulo_exm(16, 1:2) =    2.6765011517336192E-009_rp
!            vaulo_exm(17, 1:2) =   0.99999998935554446_rp
!            vaulo_exm(18, 1:2) =   0.78131589613574182_rp
!            vaulo_exm(19, 1:2) =   0.99999998935589263_rp
!            vaulo_exm(20, 1:2) =   0.98788614956431220_rp
!            vaulo_exm(21, 1:2) =   0.99431744834580349_rp
!            vaulo_exm(22, 1:2) =    9.9542465724058414E-003_rp
!            vaulo_exm(23, 1:2) =   0.99999998133288193_rp
!            vaulo_exm(24, 1:2) =   0.99999998895728170_rp
!            vaulo_exm(25, 1:2) =    6.6213956799707088E-004_rp
!            vaulo_exm(26, 1:2) =   0.69998596202985708_rp
!            vaulo_exm(27, 1:2) =   0.45745549722925427_rp
!            vaulo_exm(28, 1:2) =    2.0754361583136973E-004_rp
!            vaulo_exm(29, 1:2) =   0.99690661867951635_rp
!            !vaulo_exm(30, 1:2) =    0.0000000000000000_rp               
!         case (EXM_CELLTYPE_EPI)  !3 
!            voltage_local =   -87.397959426295145_rp
!            viclo_exm(1:26) = 0.0_rp
!            vcolo_exm(1, 1:2) =    8.5433179446188953E-005_rp
!            vcolo_exm(2, 1:2) =    9.7033871421276938_rp
!            vcolo_exm(3, 1:2) =    140.63418994101886_rp
!            vcolo_exm(4, 1:2) =    140.63414316826149_rp
!            vcolo_exm(5, 1:2) =    9.7032770674922091_rp
!            vcolo_exm(6, 1:2) =    8.3969893961209863E-005_rp
!            vcolo_exm(7, 1:2) =    8.7338011361931187_rp
!            vcolo_exm(8, 1:2) =    8.5090253044155926_rp
!            vcolo_exm(9, 1:2) =    9.4348216776708040E-007_rp
!            vcolo_exm(10, 1:2) =    1.1792687199169907E-006_rp
!            vcolo_exm(11, 1:2) =   0.23304109107375579_rp
!            !vcolo_exm(12, 1:2) =    0.0000000000000000_rp
!            vaulo_exm(1, 1:2) =    7.8038762211742538E-003_rp
!            vaulo_exm(2, 1:2) =   0.80697892118129511_rp
!            vaulo_exm(3, 1:2) =   0.80698048877554496_rp
!            vaulo_exm(4, 1:2) =   0.80427282308871484_rp
!            vaulo_exm(5, 1:2) =   0.60674144785235762_rp
!            vaulo_exm(6, 1:2) =   0.79118681882009789_rp
!            vaulo_exm(7, 1:2) =    2.1114419773127173E-004_rp
!            vaulo_exm(8, 1:2) =   0.42822745461488432_rp
!            vaulo_exm(9, 1:2) =   0.20738103535881042_rp
!            vaulo_exm(10, 1:2) =    1.0426904028929566E-003_rp
!            vaulo_exm(11, 1:2) =   0.99950453328143118_rp
!            vaulo_exm(12, 1:2) =   0.99424423356861691_rp
!            vaulo_exm(13, 1:2) =    5.3128975397284521E-004_rp
!            vaulo_exm(14, 1:2) =   0.99950453350151214_rp
!            vaulo_exm(15, 1:2) =   0.99701906127760564_rp
!            vaulo_exm(16, 1:2) =    2.7000139951840619E-009_rp
!            vaulo_exm(17, 1:2) =   0.99999998925616540_rp
!            vaulo_exm(18, 1:2) =   0.81022448363846356_rp
!            vaulo_exm(19, 1:2) =   0.99999998925619915_rp
!            vaulo_exm(20, 1:2) =   0.99151690213156662_rp
!            vaulo_exm(21, 1:2) =   0.99630964501327157_rp
!            vaulo_exm(22, 1:2) =    2.6414450407093136E-003_rp
!            vaulo_exm(23, 1:2) =   0.99999998823676040_rp
!            vaulo_exm(24, 1:2) =   0.99999998919167776_rp
!            vaulo_exm(25, 1:2) =    3.0381511578187663E-004_rp
!            vaulo_exm(26, 1:2) =   0.67685433163727815_rp
!            vaulo_exm(27, 1:2) =   0.41836968833467597_rp
!            vaulo_exm(28, 1:2) =    2.0676096444089193E-004_rp
!            vaulo_exm(29, 1:2) =   0.99691451138133580_rp
!            !vaulo_exm(30, 1:2) =    0.0000000000000000_rp                      
!         case default
!            call runend("EXM_INIT_VOLATGES: Undefined cell type.")
!         end select
!
!      case (EXM_CELL_NORMAL_HUMAN_MALE_600BCL)
!         select case (celltype)
!         case (EXM_CELLTYPE_ENDO) !1 
!            voltage_local =   -87.840677417569623_rp
!            viclo_exm(1:26) = 0.0_rp
!            vcolo_exm(1, 1:2) =    1.0642364619202648E-004_rp
!            vcolo_exm(2, 1:2) =    7.9772077683146829_rp
!            vcolo_exm(3, 1:2) =    143.86444519547052_rp
!            vcolo_exm(4, 1:2) =    143.86440869454694_rp
!            vcolo_exm(5, 1:2) =    7.9770757936138228_rp
!            vcolo_exm(6, 1:2) =    1.0772077419226267E-004_rp
!            vcolo_exm(7, 1:2) =    1.9994980946461740_rp
!            vcolo_exm(8, 1:2) =    1.7061582377871256_rp
!            vcolo_exm(9, 1:2) =    3.1469432182109148E-007_rp
!            vcolo_exm(10, 1:2) =    3.9256396771398118E-007_rp
!            vcolo_exm(11, 1:2) =    4.0081050009406739E-002_rp
!            !vcolo_exm(12, 1:2) =    0.0000000000000000_rp
!            vaulo_exm(1, 1:2) =    7.4641620415121365E-003_rp
!            vaulo_exm(2, 1:2) =   0.81778093747987457_rp
!            vaulo_exm(3, 1:2) =   0.81778996374016499_rp
!            vaulo_exm(4, 1:2) =   0.81525534740532923_rp
!            vaulo_exm(5, 1:2) =   0.62341439722652625_rp
!            vaulo_exm(6, 1:2) =   0.80380141909773406_rp
!            vaulo_exm(7, 1:2) =    1.9411613137063226E-004_rp
!            vaulo_exm(8, 1:2) =   0.43868600477080355_rp
!            vaulo_exm(9, 1:2) =   0.21407703205078979_rp
!            vaulo_exm(10, 1:2) =    1.0120484465101275E-003_rp
!            vaulo_exm(11, 1:2) =   0.99954127668039483_rp
!            vaulo_exm(12, 1:2) =   0.35456205922719886_rp
!            vaulo_exm(13, 1:2) =    5.1566877333477810E-004_rp
!            vaulo_exm(14, 1:2) =   0.99954130276542819_rp
!            vaulo_exm(15, 1:2) =   0.39798819179101697_rp
!            vaulo_exm(16, 1:2) =    2.4318143636470063E-009_rp
!            vaulo_exm(17, 1:2) =   0.99999999046347776_rp
!            vaulo_exm(18, 1:2) =   0.81829227729265408_rp
!            vaulo_exm(19, 1:2) =   0.99999999046370502_rp
!            vaulo_exm(20, 1:2) =   0.99159904776531749_rp
!            vaulo_exm(21, 1:2) =   0.99607918538745932_rp
!            vaulo_exm(22, 1:2) =    6.8365035674538763E-003_rp
!            vaulo_exm(23, 1:2) =   0.99999998908994181_rp
!            vaulo_exm(24, 1:2) =   0.99999999037131027_rp
!            vaulo_exm(25, 1:2) =    2.9908608383753262E-004_rp
!            vaulo_exm(26, 1:2) =   0.67433044628262184_rp
!            vaulo_exm(27, 1:2) =   0.41460644563979088_rp
!            vaulo_exm(28, 1:2) =    1.9695576450942084E-004_rp
!            vaulo_exm(29, 1:2) =   0.99680275895078962_rp
!            !vaulo_exm(30, 1:2) =    0.0000000000000000_rp            
!         case (EXM_CELLTYPE_MID)  !2    
!            voltage_local =   -87.435195284762415_rp
!            viclo_exm(1:26) = 0.0_rp
!            vcolo_exm(1, 1:2) =    1.1890657170525055E-004_rp
!            vcolo_exm(2, 1:2) =    9.7706696760397378_rp
!            vcolo_exm(3, 1:2) =    141.66483952218840_rp
!            vcolo_exm(4, 1:2) =    141.66479112719344_rp
!            vcolo_exm(5, 1:2) =    9.7704648272293344_rp
!            vcolo_exm(6, 1:2) =    1.1898783498170927E-004_rp
!            vcolo_exm(7, 1:2) =    2.8021201980086876_rp
!            vcolo_exm(8, 1:2) =    2.4236072150570278_rp
!            vcolo_exm(9, 1:2) =    1.4071451498516887E-006_rp
!            vcolo_exm(10, 1:2) =    1.7586168867274698E-006_rp
!            vcolo_exm(11, 1:2) =    7.8962722885925360E-002_rp
!            !vcolo_exm(12, 1:2) =    0.0000000000000000_rp
!            vaulo_exm(1, 1:2) =    7.7747291376501189E-003_rp
!            vaulo_exm(2, 1:2) =   0.80784540409976191_rp
!            vaulo_exm(3, 1:2) =   0.80786105547820475_rp
!            vaulo_exm(4, 1:2) =   0.80265157891836159_rp
!            vaulo_exm(5, 1:2) =   0.60792453114967615_rp
!            vaulo_exm(6, 1:2) =   0.78403552237781493_rp
!            vaulo_exm(7, 1:2) =    2.0965657866126548E-004_rp
!            vaulo_exm(8, 1:2) =   0.41110392493267284_rp
!            vaulo_exm(9, 1:2) =   0.19178946886523357_rp
!            vaulo_exm(10, 1:2) =    1.0400968637168772E-003_rp
!            vaulo_exm(11, 1:2) =   0.99950747422536579_rp
!            vaulo_exm(12, 1:2) =   0.32788959685532149_rp
!            vaulo_exm(13, 1:2) =    5.2996757393124953E-004_rp
!            vaulo_exm(14, 1:2) =   0.99950751019835571_rp
!            vaulo_exm(15, 1:2) =   0.36509196261881782_rp
!            vaulo_exm(16, 1:2) =    2.6765011517336192E-009_rp
!            vaulo_exm(17, 1:2) =   0.99999998935554446_rp
!            vaulo_exm(18, 1:2) =   0.78131589613574182_rp
!            vaulo_exm(19, 1:2) =   0.99999998935589263_rp
!            vaulo_exm(20, 1:2) =   0.98788614956431220_rp
!            vaulo_exm(21, 1:2) =   0.99431744834580349_rp
!            vaulo_exm(22, 1:2) =    9.9542465724058414E-003_rp
!            vaulo_exm(23, 1:2) =   0.99999998133288193_rp
!            vaulo_exm(24, 1:2) =   0.99999998895728170_rp
!            vaulo_exm(25, 1:2) =    6.6213956799707088E-004_rp
!            vaulo_exm(26, 1:2) =   0.69998596202985708_rp
!            vaulo_exm(27, 1:2) =   0.45745549722925427_rp
!            vaulo_exm(28, 1:2) =    2.0754361583136973E-004_rp
!            vaulo_exm(29, 1:2) =   0.99690661867951635_rp
!            !vaulo_exm(30, 1:2) =    0.0000000000000000_rp               
!         case (EXM_CELLTYPE_EPI)  !3 
!            voltage_local =   -87.708936866906583_rp
!            viclo_exm(1:26) = 0.0_rp
!            vcolo_exm(1, 1:2) =    8.0028671995560340E-005_rp
!            vcolo_exm(2, 1:2) =    8.6136328499931931_rp
!            vcolo_exm(3, 1:2) =    142.87867795811988_rp
!            vcolo_exm(4, 1:2) =    142.87863578222584_rp
!            vcolo_exm(5, 1:2) =    8.6135333361932869_rp
!            vcolo_exm(6, 1:2) =    8.0347466740107238E-005_rp
!            vcolo_exm(7, 1:2) =    3.3957368575898395_rp
!            vcolo_exm(8, 1:2) =    3.1113557277334025_rp
!            vcolo_exm(9, 1:2) =    4.7491800922715200E-007_rp
!            vcolo_exm(10, 1:2) =    5.9362015321547805E-007_rp
!            vcolo_exm(11, 1:2) =    3.9107271067814031E-002_rp
!            !vcolo_exm(12, 1:2) =    0.0000000000000000_rp
!            vaulo_exm(1, 1:2) =    7.5636859241911452E-003_rp
!            vaulo_exm(2, 1:2) =   0.81463650323285075_rp
!            vaulo_exm(3, 1:2) =   0.81463987634369694_rp
!            vaulo_exm(4, 1:2) =   0.81361739286268686_rp
!            vaulo_exm(5, 1:2) =   0.61856297179956743_rp
!            vaulo_exm(6, 1:2) =   0.80730247582796777_rp
!            vaulo_exm(7, 1:2) =    1.9903435708260969E-004_rp
!            vaulo_exm(8, 1:2) =   0.45892752188085889_rp
!            vaulo_exm(9, 1:2) =   0.23701050019239650_rp
!            vaulo_exm(10, 1:2) =    1.0210647092482490E-003_rp
!            vaulo_exm(11, 1:2) =   0.99953077560081349_rp
!            vaulo_exm(12, 1:2) =   0.99757663568230093_rp
!            vaulo_exm(13, 1:2) =    5.2026513031970978E-004_rp
!            vaulo_exm(14, 1:2) =   0.99953077622510722_rp
!            vaulo_exm(15, 1:2) =   0.99877342602002361_rp
!            vaulo_exm(16, 1:2) =    2.5086639233833685E-009_rp
!            vaulo_exm(17, 1:2) =   0.99999999012175589_rp
!            vaulo_exm(18, 1:2) =   0.86178725892506836_rp
!            vaulo_exm(19, 1:2) =   0.99999999012183993_rp
!            vaulo_exm(20, 1:2) =   0.99587532884666397_rp
!            vaulo_exm(21, 1:2) =   0.99788907700093699_rp
!            vaulo_exm(22, 1:2) =    2.2320899641454662E-003_rp
!            vaulo_exm(23, 1:2) =   0.99999999004661755_rp
!            vaulo_exm(24, 1:2) =   0.99999999011105045_rp
!            vaulo_exm(25, 1:2) =    9.4095695945113766E-005_rp
!            vaulo_exm(26, 1:2) =   0.62481556900203306_rp
!            vaulo_exm(27, 1:2) =   0.34825026454480218_rp
!            vaulo_exm(28, 1:2) =    1.9934081572285325E-004_rp
!            vaulo_exm(29, 1:2) =   0.99683592429184664_rp
!            !vaulo_exm(30, 1:2) =    0.0000000000000000_rp                   
!         case default
!            call runend("EXM_INIT_VOLATGES: Undefined cell type.")
!         end select
!
!      case (EXM_CELL_NORMAL_70BPM)
!         select case (celltype)
!         case (EXM_CELLTYPE_ENDO)
!            voltage_local = -0.878052586716E+02_rp
!            viclo_exm(1:26) = 0.0_rp
!            vcolo_exm(1, 1:2) = 0.118326831303E-03_rp
!            vcolo_exm(2, 1:2) = 0.723782510580E+01_rp
!            vcolo_exm(3, 1:2) = 0.144381861856E+03_rp
!            vcolo_exm(4, 1:2) = 0.144381838468E+03_rp
!            vcolo_exm(5, 1:2) = 0.723768860481E+01_rp
!            vcolo_exm(6, 1:2) = 0.116484844848E-03_rp
!            vcolo_exm(7, 1:2) = 0.168920657095E+01_rp
!            vcolo_exm(8, 1:2) = 0.157494154916E+01_rp
!            vcolo_exm(9, 1:2) = 0.260097028716E-06_rp
!            vcolo_exm(10, 1:2) = 0.324744067913E-06_rp
!            vcolo_exm(11, 1:2) = 0.167951236385E-01_rp
!            vaulo_exm(1, 1:2) = 0.749079173748E-02_rp
!            vaulo_exm(2, 1:2) = 0.816932376673E+00_rp
!            vaulo_exm(3, 1:2) = 0.816941162831E+00_rp
!            vaulo_exm(4, 1:2) = 0.816438594186E+00_rp
!            vaulo_exm(5, 1:2) = 0.622090167494E+00_rp
!            vaulo_exm(6, 1:2) = 0.815305688267E+00_rp
!            vaulo_exm(7, 1:2) = 0.195426384241E-03_rp
!            vaulo_exm(8, 1:2) = 0.482779436873E+00_rp
!            vaulo_exm(9, 1:2) = 0.248901759996E+00_rp
!            vaulo_exm(10, 1:2) = 0.101446734364E-02_rp
!            vaulo_exm(11, 1:2) = 0.999538432028E+00_rp
!            vaulo_exm(12, 1:2) = 0.511668738557E+00_rp
!            vaulo_exm(13, 1:2) = 0.516901887213E-03_rp
!            vaulo_exm(14, 1:2) = 0.999538456645E+00_rp
!            vaulo_exm(15, 1:2) = 0.561580057347E+00_rp
!            vaulo_exm(16, 1:2) = 0.245226010924E-08_rp
!            vaulo_exm(17, 1:2) = 0.999999990372E+00_rp
!            vaulo_exm(18, 1:2) = 0.882561625277E+00_rp
!            vaulo_exm(19, 1:2) = 0.999999990372E+00_rp
!            vaulo_exm(20, 1:2) = 0.999211409367E+00_rp
!            vaulo_exm(21, 1:2) = 0.999831263545E+00_rp
!            vaulo_exm(22, 1:2) = 0.911276398365E-02_rp
!            vaulo_exm(23, 1:2) = 0.999999990362E+00_rp
!            vaulo_exm(24, 1:2) = 0.999999990363E+00_rp
!            vaulo_exm(25, 1:2) = 0.982010014179E-05_rp
!            vaulo_exm(26, 1:2) = 0.528283180130E+00_rp
!            vaulo_exm(27, 1:2) = 0.324804396838E+00_rp
!            vaulo_exm(28, 1:2) = 0.197273510439E-03_rp
!            vaulo_exm(29, 1:2) = 0.996811854061E+00_rp
!            outputs % S      = 0.877167079580E-01_rp
!            outputs % W      = 0.111808267657E+00_rp
!            outputs % CaTRPN = 0.293438024800E+00_rp
!            outputs % B      = 0.690381108301E+00_rp
!            outputs % zeta_s = 0.000000000000E+00_rp
!            outputs % zeta_w = 0.000000000000E+00_rp
!            outputs % Lambda = 0.100000000000E+01_rp
!            outputs % Ca50   = 0.184455240483E-03_rp
!            outputs % nbeats = 104_ip
!            outputs % toler  = 0.100000000000E-06_rp
!            outputs % toler  = 0.100000000000E-06_rp
!            outputs % rmse   = 0.944425989612E-07_rp
!         case (EXM_CELLTYPE_MID)
!            voltage_local = -0.875289884127E+02_rp
!            viclo_exm(1:26) = 0.0_rp
!            vcolo_exm(1, 1:2) = 0.112278054555E-03_rp
!            vcolo_exm(2, 1:2) = 0.877851764133E+01_rp
!            vcolo_exm(3, 1:2) = 0.142550239729E+03_rp
!            vcolo_exm(4, 1:2) = 0.142550207974E+03_rp
!            vcolo_exm(5, 1:2) = 0.877834713310E+01_rp
!            vcolo_exm(6, 1:2) = 0.111212639253E-03_rp
!            vcolo_exm(7, 1:2) = 0.178855084743E+01_rp
!            vcolo_exm(8, 1:2) = 0.158081620823E+01_rp
!            vcolo_exm(9, 1:2) = 0.847236257984E-06_rp
!            vcolo_exm(10, 1:2) = 0.105683695290E-05_rp
!            vcolo_exm(11, 1:2) = 0.275280380008E-01_rp
!            vaulo_exm(1, 1:2) = 0.770177275796E-02_rp
!            vaulo_exm(2, 1:2) = 0.810158121678E+00_rp
!            vaulo_exm(3, 1:2) = 0.810176393794E+00_rp
!            vaulo_exm(4, 1:2) = 0.809129729048E+00_rp
!            vaulo_exm(5, 1:2) = 0.611452878514E+00_rp
!            vaulo_exm(6, 1:2) = 0.806563972098E+00_rp
!            vaulo_exm(7, 1:2) = 0.205954902016E-03_rp
!            vaulo_exm(8, 1:2) = 0.463883172827E+00_rp
!            vaulo_exm(9, 1:2) = 0.227271991897E+00_rp
!            vaulo_exm(10, 1:2) = 0.103354801724E-02_rp
!            vaulo_exm(11, 1:2) = 0.999515422768E+00_rp
!            vaulo_exm(12, 1:2) = 0.479642076308E+00_rp
!            vaulo_exm(13, 1:2) = 0.526629002717E-03_rp
!            vaulo_exm(14, 1:2) = 0.999515466555E+00_rp
!            vaulo_exm(15, 1:2) = 0.523968911203E+00_rp
!            vaulo_exm(16, 1:2) = 0.261785295470E-08_rp
!            vaulo_exm(17, 1:2) = 0.999999989620E+00_rp
!            vaulo_exm(18, 1:2) = 0.851333896560E+00_rp
!            vaulo_exm(19, 1:2) = 0.999999989620E+00_rp
!            vaulo_exm(20, 1:2) = 0.998638318416E+00_rp
!            vaulo_exm(21, 1:2) = 0.999701917081E+00_rp
!            vaulo_exm(22, 1:2) = 0.767345137585E-02_rp
!            vaulo_exm(23, 1:2) = 0.999999989602E+00_rp
!            vaulo_exm(24, 1:2) = 0.999999989603E+00_rp
!            vaulo_exm(25, 1:2) = 0.142665291682E-04_rp
!            vaulo_exm(26, 1:2) = 0.560309894419E+00_rp
!            vaulo_exm(27, 1:2) = 0.372517016433E+00_rp
!            vaulo_exm(28, 1:2) = 0.203623514969E-03_rp
!            vaulo_exm(29, 1:2) = 0.996883357642E+00_rp
!            outputs % S      = 0.429941660108E-02_rp
!            outputs % W      = 0.406202412166E-02_rp
!            outputs % CaTRPN = 0.125236626983E+00_rp
!            outputs % B      = 0.987678042169E+00_rp
!            outputs % zeta_s = 0.000000000000E+00_rp
!            outputs % zeta_w = 0.000000000000E+00_rp
!            outputs % Lambda = 0.100000000000E+01_rp
!            outputs % Ca50   = 0.299581856665E-03_rp
!            outputs % nbeats = 221_ip
!            outputs % toler  = 0.100000000000E-06_rp
!            outputs % toler  = 0.100000000000E-06_rp
!            outputs % rmse   = 0.980941729834E-07_rp
!         case (EXM_CELLTYPE_EPI)
!            voltage_local = -0.877931532701E+02_rp
!            viclo_exm(1:26) = 0.0_rp
!            vcolo_exm(1, 1:2) = 0.994787619620E-04_rp
!            vcolo_exm(2, 1:2) = 0.769489351314E+01_rp
!            vcolo_exm(3, 1:2) = 0.143827951691E+03_rp
!            vcolo_exm(4, 1:2) = 0.143827923802E+03_rp
!            vcolo_exm(5, 1:2) = 0.769477670275E+01_rp
!            vcolo_exm(6, 1:2) = 0.992004473998E-04_rp
!            vcolo_exm(7, 1:2) = 0.207794427481E+01_rp
!            vcolo_exm(8, 1:2) = 0.192888143564E+01_rp
!            vcolo_exm(9, 1:2) = 0.468281397603E-06_rp
!            vcolo_exm(10, 1:2) = 0.585116054833E-06_rp
!            vcolo_exm(11, 1:2) = 0.206895410329E-01_rp
!            vaulo_exm(1, 1:2) = 0.749991478055E-02_rp
!            vaulo_exm(2, 1:2) = 0.816641316897E+00_rp
!            vaulo_exm(3, 1:2) = 0.816650190115E+00_rp
!            vaulo_exm(4, 1:2) = 0.816129466707E+00_rp
!            vaulo_exm(5, 1:2) = 0.621631279295E+00_rp
!            vaulo_exm(6, 1:2) = 0.815229260230E+00_rp
!            vaulo_exm(7, 1:2) = 0.195876224139E-03_rp
!            vaulo_exm(8, 1:2) = 0.486887454889E+00_rp
!            vaulo_exm(9, 1:2) = 0.257572541416E+00_rp
!            vaulo_exm(10, 1:2) = 0.101529529320E-02_rp
!            vaulo_exm(11, 1:2) = 0.999537631567E+00_rp
!            vaulo_exm(12, 1:2) = 0.999457318802E+00_rp
!            vaulo_exm(13, 1:2) = 0.517323962919E-03_rp
!            vaulo_exm(14, 1:2) = 0.999537633368E+00_rp
!            vaulo_exm(15, 1:2) = 0.999516007371E+00_rp
!            vaulo_exm(16, 1:2) = 0.245928662597E-08_rp
!            vaulo_exm(17, 1:2) = 0.999999990340E+00_rp
!            vaulo_exm(18, 1:2) = 0.899358308642E+00_rp
!            vaulo_exm(19, 1:2) = 0.999999990341E+00_rp
!            vaulo_exm(20, 1:2) = 0.999517320566E+00_rp
!            vaulo_exm(21, 1:2) = 0.999885908527E+00_rp
!            vaulo_exm(22, 1:2) = 0.497892588125E-02_rp
!            vaulo_exm(23, 1:2) = 0.999999990331E+00_rp
!            vaulo_exm(24, 1:2) = 0.999999990331E+00_rp
!            vaulo_exm(25, 1:2) = 0.897176485003E-05_rp
!            vaulo_exm(26, 1:2) = 0.497451830404E+00_rp
!            vaulo_exm(27, 1:2) = 0.291821510965E+00_rp
!            vaulo_exm(28, 1:2) = 0.197545150717E-03_rp
!            vaulo_exm(29, 1:2) = 0.996814989370E+00_rp
!            outputs % S      = 0.539101167324E-02_rp
!            outputs % W      = 0.536019939117E-02_rp
!            outputs % CaTRPN = 0.140557297625E+00_rp
!            outputs % B      = 0.983994547661E+00_rp
!            outputs % zeta_s = 0.000000000000E+00_rp
!            outputs % zeta_w = 0.000000000000E+00_rp
!            outputs % Lambda = 0.100000000000E+01_rp
!            outputs % Ca50   = 0.247894539227E-03_rp
!            outputs % nbeats = 120_ip
!            outputs % toler  = 0.100000000000E-06_rp
!            outputs % toler  = 0.100000000000E-06_rp
!            outputs % rmse   = 0.977785966109E-07_rp
!         case default
!            call runend("EXM_INIT_VOLATGES: Undefined cell type.")
!         end select
!
!      case (EXM_CELL_NORMAL_1000BCL)
!         select case (celltype)
!         case (EXM_CELLTYPE_ENDO)
!            !endocardium
!            vminimate_exm(celltype, mat) = -87.0_rp     ! INTRA (BI-DOMAIN) OR SINGLE (MONO) POTENTIAL
!            !voltage_local = -87.99_rp
!            viclo_exm(1:26) = 0.0_rp
!            vconc_initial(1, celltype, mat) = 0.0001_rp  !8.03149767106260e-05       !cai=6)
!            vconc_initial(2, celltype, mat) = 7.0_rp !8.17212071948942       !nass=3);
!            vconc_initial(3, celltype, mat) = 145.0_rp !143.675184333376       !ki=4);
!            vconc_initial(4, celltype, mat) = 145.0_rp  !143.675147146842       !kss=5);
!            vconc_initial(5, celltype, mat) = 7.0_rp  !8.17202753912090      ! nai=2
!            vconc_initial(6, celltype, mat) = 0.0001_rp  !7.96912949697674e-05       !cass=7);
!            vconc_initial(7, celltype, mat) = 1.2_rp  !2.14811455007091       !cansr=8);
!            vconc_initial(8, celltype, mat) = 1.2_rp  !2.03335899236798      !cajsr=9);
!            vauxi_exm_initial(1, celltype, mat) = 0.0_rp !0.00739719746920272       !m=10);
!            vauxi_exm_initial(2, celltype, mat) = 1.0_rp  !0.695621622011335       !hf=11);
!            vauxi_exm_initial(3, celltype, mat) = 1.0_rp  !0.695601842634086       !hs=12);
!            vauxi_exm_initial(4, celltype, mat) = 1.0_rp  !0.695486248719023       !j=13);
!            vauxi_exm_initial(5, celltype, mat) = 1.0_rp  !0.452023628358454       !hsp=14);
!            vauxi_exm_initial(6, celltype, mat) = 1.0_rp  !0.695403157533235       !jp=15);
!            vauxi_exm_initial(7, celltype, mat) = 0.0_rp !0.000190839777466418       !mL=16);
!            vauxi_exm_initial(8, celltype, mat) = 1.0_rp  !0.493606704642336       !hL=17);
!            vauxi_exm_initial(9, celltype, mat) = 1.0_rp  !0.264304293390731       !hLp=18);
!            vauxi_exm_initial(10, celltype, mat) = 0.0_rp !0.00100594231451985      !a=19);
!            vauxi_exm_initial(11, celltype, mat) = 1.0_rp  !0.999548606668578      !iF=20);
!            vauxi_exm_initial(12, celltype, mat) = 1.0_rp  !0.999488774162635      !iS=21);
!            vauxi_exm_initial(13, celltype, mat) = 0.0_rp !0.000512555980943569      !ap=22);
!            vauxi_exm_initial(14, celltype, mat) = 1.0_rp  !0.999548607287668      !iFp=23);
!            vauxi_exm_initial(15, celltype, mat) = 1.0_rp  !0.999488774162635      !iSp=24);
!            vauxi_exm_initial(16, celltype, mat) = 0.0_rp !2.38076098345898e-09      !d=25);
!            vauxi_exm_initial(17, celltype, mat) = 1.0_rp  !0.999999990696210      !ff=26);
!            vauxi_exm_initial(18, celltype, mat) = 1.0_rp  !0.904906458666787      !fs=27);
!            vauxi_exm_initial(19, celltype, mat) = 1.0_rp  !0.999999990696060      !fcaf=28);
!            vauxi_exm_initial(20, celltype, mat) = 1.0_rp  !0.999581201974281      !fcas=29);
!            vauxi_exm_initial(21, celltype, mat) = 1.0_rp  !0.999903346883777      !jca=30);
!            vauxi_exm_initial(22, celltype, mat) = 0.0_rp !0.00215555277945401      !nca=31);
!            vauxi_exm_initial(23, celltype, mat) = 1.0_rp  !0.999999990680285      !ffp=32);
!            vauxi_exm_initial(24, celltype, mat) = 1.0_rp  !0.999999990692529      !fcafp=33);
!            vauxi_exm_initial(25, celltype, mat) = 0.0_rp !8.64222375034682e-06      !xrf=34);
!            vauxi_exm_initial(26, celltype, mat) = 0.0_rp !0.487585264457487      !xrs=35);
!            vauxi_exm_initial(27, celltype, mat) = 0.0_rp !0.276203479404767      !xs1=36);
!            vauxi_exm_initial(28, celltype, mat) = 0.0_rp !0.000194412216700766      !xs2=37);
!            vauxi_exm_initial(29, celltype, mat) = 1.0_rp  !0.996778581263402      !xk1=38);
!            vconc_initial(9,  celltype, mat) = 0.0_rp !4.74946280300893e-07       !Jrelnp=39);
!            vconc_initial(10, celltype, mat) = 0.0_rp !5.93539009244893e-07      !Jrelp=40);
!            vconc_initial(11, celltype, mat) = 0.0_rp !0.0228529042639590      !CaMKt=41);
!         case (EXM_CELLTYPE_EPI)
!            vminimate_exm(celltype, mat) = -87.0_rp     ! INTRA (BI-DOMAIN) OR SINGLE (MONO) POTENTIAL
!            !voltage_local = -87.99_rp
!            viclo_exm(1:26) = 0.0_rp
!            vconc_initial(1, celltype, mat) = 0.0001_rp  !8.03149767106260e-05       !cai=6)
!            vconc_initial(2, celltype, mat) = 7.0_rp !8.17212071948942       !nass=3);
!            vconc_initial(3, celltype, mat) = 145.0_rp !143.675184333376       !ki=4);
!            vconc_initial(4, celltype, mat) = 145.0_rp  !143.675147146842       !kss=5);
!            vconc_initial(5, celltype, mat) = 7.0_rp  !8.17202753912090      ! nai=2
!            vconc_initial(6, celltype, mat) = 0.0001_rp  !7.96912949697674e-05       !cass=7);
!            vconc_initial(7, celltype, mat) = 1.2_rp  !2.14811455007091       !cansr=8);
!            vconc_initial(8, celltype, mat) = 1.2_rp  !2.03335899236798      !cajsr=9);
!            vauxi_exm_initial(1, celltype, mat) = 0.0_rp !0.00739719746920272       !m=10);
!            vauxi_exm_initial(2, celltype, mat) = 1.0_rp  !0.695621622011335       !hf=11);
!            vauxi_exm_initial(3, celltype, mat) = 1.0_rp  !0.695601842634086       !hs=12);
!            vauxi_exm_initial(4, celltype, mat) = 1.0_rp  !0.695486248719023       !j=13);
!            vauxi_exm_initial(5, celltype, mat) = 1.0_rp  !0.452023628358454       !hsp=14);
!            vauxi_exm_initial(6, celltype, mat) = 1.0_rp  !0.695403157533235       !jp=15);
!            vauxi_exm_initial(7, celltype, mat) = 0.0_rp !0.000190839777466418       !mL=16);
!            vauxi_exm_initial(8, celltype, mat) = 1.0_rp  !0.493606704642336       !hL=17);
!            vauxi_exm_initial(9, celltype, mat) = 1.0_rp  !0.264304293390731       !hLp=18);
!            vauxi_exm_initial(10, celltype, mat) = 0.0_rp !0.00100594231451985      !a=19);
!            vauxi_exm_initial(11, celltype, mat) = 1.0_rp  !0.999548606668578      !iF=20);
!            vauxi_exm_initial(12, celltype, mat) = 1.0_rp  !0.999488774162635      !iS=21);
!            vauxi_exm_initial(13, celltype, mat) = 0.0_rp !0.000512555980943569      !ap=22);
!            vauxi_exm_initial(14, celltype, mat) = 1.0_rp  !0.999548607287668      !iFp=23);
!            vauxi_exm_initial(15, celltype, mat) = 1.0_rp  !0.999488774162635      !iSp=24);
!            vauxi_exm_initial(16, celltype, mat) = 0.0_rp !2.38076098345898e-09      !d=25);
!            vauxi_exm_initial(17, celltype, mat) = 1.0_rp  !0.999999990696210      !ff=26);
!            vauxi_exm_initial(18, celltype, mat) = 1.0_rp  !0.904906458666787      !fs=27);
!            vauxi_exm_initial(19, celltype, mat) = 1.0_rp  !0.999999990696060      !fcaf=28);
!            vauxi_exm_initial(20, celltype, mat) = 1.0_rp  !0.999581201974281      !fcas=29);
!            vauxi_exm_initial(21, celltype, mat) = 1.0_rp  !0.999903346883777      !jca=30);
!            vauxi_exm_initial(22, celltype, mat) = 0.0_rp !0.00215555277945401      !nca=31);
!            vauxi_exm_initial(23, celltype, mat) = 1.0_rp  !0.999999990680285      !ffp=32);
!            vauxi_exm_initial(24, celltype, mat) = 1.0_rp  !0.999999990692529      !fcafp=33);
!            vauxi_exm_initial(25, celltype, mat) = 0.0_rp !8.64222375034682e-06      !xrf=34);
!            vauxi_exm_initial(26, celltype, mat) = 0.0_rp !0.487585264457487      !xrs=35);
!            vauxi_exm_initial(27, celltype, mat) = 0.0_rp !0.276203479404767      !xs1=36);
!            vauxi_exm_initial(28, celltype, mat) = 0.0_rp !0.000194412216700766      !xs2=37);
!            vauxi_exm_initial(29, celltype, mat) = 1.0_rp  !0.996778581263402      !xk1=38);
!            vconc_initial(9,  celltype, mat) = 0.0_rp !4.74946280300893e-07       !Jrelnp=39);
!            vconc_initial(10, celltype, mat) = 0.0_rp !5.93539009244893e-07      !Jrelp=40);
!            vconc_initial(11, celltype, mat) = 0.0_rp !0.0228529042639590      !CaMKt=41);
!         case (EXM_CELLTYPE_MID)
!            vminimate_exm(celltype, mat) = -87.0_rp     ! INTRA (BI-DOMAIN) OR SINGLE (MONO) POTENTIAL
!            !voltage_local = -87.99_rp
!            viclo_exm(1:26) = 0.0_rp
!            vconc_initial(1, celltype, mat) = 0.0001_rp  !8.03149767106260e-05       !cai=6)
!            vconc_initial(2, celltype, mat) = 7.0_rp !8.17212071948942       !nass=3);
!            vconc_initial(3, celltype, mat) = 145.0_rp !143.675184333376       !ki=4);
!            vconc_initial(4, celltype, mat) = 145.0_rp  !143.675147146842       !kss=5);
!            vconc_initial(5, celltype, mat) = 7.0_rp  !8.17202753912090      ! nai=2
!            vconc_initial(6, celltype, mat) = 0.0001_rp  !7.96912949697674e-05       !cass=7);
!            vconc_initial(7, celltype, mat) = 1.2_rp  !2.14811455007091       !cansr=8);
!            vconc_initial(8, celltype, mat) = 1.2_rp  !2.03335899236798      !cajsr=9);
!            vauxi_exm_initial(1, celltype, mat) = 0.0_rp !0.00739719746920272       !m=10);
!            vauxi_exm_initial(2, celltype, mat) = 1.0_rp  !0.695621622011335       !hf=11);
!            vauxi_exm_initial(3, celltype, mat) = 1.0_rp  !0.695601842634086       !hs=12);
!            vauxi_exm_initial(4, celltype, mat) = 1.0_rp  !0.695486248719023       !j=13);
!            vauxi_exm_initial(5, celltype, mat) = 1.0_rp  !0.452023628358454       !hsp=14);
!            vauxi_exm_initial(6, celltype, mat) = 1.0_rp  !0.695403157533235       !jp=15);
!            vauxi_exm_initial(7, celltype, mat) = 0.0_rp !0.000190839777466418       !mL=16);
!            vauxi_exm_initial(8, celltype, mat) = 1.0_rp  !0.493606704642336       !hL=17);
!            vauxi_exm_initial(9, celltype, mat) = 1.0_rp  !0.264304293390731       !hLp=18);
!            vauxi_exm_initial(10, celltype, mat) = 0.0_rp !0.00100594231451985      !a=19);
!            vauxi_exm_initial(11, celltype, mat) = 1.0_rp  !0.999548606668578      !iF=20);
!            vauxi_exm_initial(12, celltype, mat) = 1.0_rp  !0.999488774162635      !iS=21);
!            vauxi_exm_initial(13, celltype, mat) = 0.0_rp !0.000512555980943569      !ap=22);
!            vauxi_exm_initial(14, celltype, mat) = 1.0_rp  !0.999548607287668      !iFp=23);
!            vauxi_exm_initial(15, celltype, mat) = 1.0_rp  !0.999488774162635      !iSp=24);
!            vauxi_exm_initial(16, celltype, mat) = 0.0_rp !2.38076098345898e-09      !d=25);
!            vauxi_exm_initial(17, celltype, mat) = 1.0_rp  !0.999999990696210      !ff=26);
!            vauxi_exm_initial(18, celltype, mat) = 1.0_rp  !0.904906458666787      !fs=27);
!            vauxi_exm_initial(19, celltype, mat) = 1.0_rp  !0.999999990696060      !fcaf=28);
!            vauxi_exm_initial(20, celltype, mat) = 1.0_rp  !0.999581201974281      !fcas=29);
!            vauxi_exm_initial(21, celltype, mat) = 1.0_rp  !0.999903346883777      !jca=30);
!            vauxi_exm_initial(22, celltype, mat) = 0.0_rp !0.00215555277945401      !nca=31);
!            vauxi_exm_initial(23, celltype, mat) = 1.0_rp  !0.999999990680285      !ffp=32);
!            vauxi_exm_initial(24, celltype, mat) = 1.0_rp  !0.999999990692529      !fcafp=33);
!            vauxi_exm_initial(25, celltype, mat) = 0.0_rp !8.64222375034682e-06      !xrf=34);
!            vauxi_exm_initial(26, celltype, mat) = 0.0_rp !0.487585264457487      !xrs=35);
!            vauxi_exm_initial(27, celltype, mat) = 0.0_rp !0.276203479404767      !xs1=36);
!            vauxi_exm_initial(28, celltype, mat) = 0.0_rp !0.000194412216700766      !xs2=37);
!            vauxi_exm_initial(29, celltype, mat) = 1.0_rp  !0.996778581263402      !xk1=38);
!            vconc_initial(9,  celltype, mat) = 0.0_rp !4.74946280300893e-07       !Jrelnp=39);
!            vconc_initial(10, celltype, mat) = 0.0_rp !5.93539009244893e-07      !Jrelp=40);
!            vconc_initial(11, celltype, mat) = 0.0_rp !0.0228529042639590      !CaMKt=41);
!         case default
!            call runend("EXM_INIT_VOLTAGES: Undefined cell type.")
!         end select
!
      case (EXM_CELL_NORMAL_MODIFIED)
         select case (celltype)
         case (EXM_CELLTYPE_ENDO)
            !endocardium
            voltage_local = -87.5_rp     ! INTRA (BI-DOMAIN) OR SINGLE (MONO) POTENTIAL
            !voltage_local = -87.99_rp
            viclo_exm(1:26) = 0.0_rp
            vcolo_exm(1, 1:2) = 0.0001_rp  !8.03149767106260e-05       !cai=6)
            vcolo_exm(2, 1:2) = 7.0_rp !8.17212071948942       !nass=3);
            vcolo_exm(3, 1:2) = 145.0_rp !143.675184333376       !ki=4);
            vcolo_exm(4, 1:2) = 145.0_rp  !143.675147146842       !kss=5);
            vcolo_exm(5, 1:2) = 7.0_rp  !8.17202753912090      ! nai=2
            vcolo_exm(6, 1:2) = 0.0001_rp  !7.96912949697674e-05       !cass=7);
            vcolo_exm(7, 1:2) = 1.2_rp  !2.14811455007091       !cansr=8);
            vcolo_exm(8, 1:2) = 1.2_rp  !2.03335899236798      !cajsr=9);
            vcolo_exm(9, 1:2) = 0.0_rp !4.74946280300893e-07       !Jrelnp=39);
            vcolo_exm(10, 1:2) = 0.0_rp !5.93539009244893e-07      !Jrelp=40);
            vcolo_exm(11, 1:2) = 0.0_rp !0.0228529042639590      !CaMKt=41);
            vaulo_exm(1, 1:2) = 0.0_rp !0.00739719746920272       !m=10);
            vaulo_exm(2, 1:2) = 1.0_rp  !0.695621622011335       !hf=11);
            vaulo_exm(3, 1:2) = 1.0_rp  !0.695601842634086       !hs=12);
            vaulo_exm(4, 1:2) = 1.0_rp  !0.695486248719023       !j=13);
            vaulo_exm(5, 1:2) = 1.0_rp  !0.452023628358454       !hsp=14);
            vaulo_exm(6, 1:2) = 1.0_rp  !0.695403157533235       !jp=15);
            vaulo_exm(7, 1:2) = 0.0_rp !0.000190839777466418       !mL=16);
            vaulo_exm(8, 1:2) = 1.0_rp  !0.493606704642336       !hL=17);
            vaulo_exm(9, 1:2) = 1.0_rp  !0.264304293390731       !hLp=18);
            vaulo_exm(10, 1:2) = 0.0_rp !0.00100594231451985      !a=19);
            vaulo_exm(11, 1:2) = 1.0_rp  !0.999548606668578      !iF=20);
            vaulo_exm(12, 1:2) = 1.0_rp  !0.999488774162635      !iS=21);
            vaulo_exm(13, 1:2) = 0.0_rp !0.000512555980943569      !ap=22);
            vaulo_exm(14, 1:2) = 1.0_rp  !0.999548607287668      !iFp=23);
            vaulo_exm(15, 1:2) = 1.0_rp  !0.999488774162635      !iSp=24);
            vaulo_exm(16, 1:2) = 0.0_rp !2.38076098345898e-09      !d=25);
            vaulo_exm(17, 1:2) = 1.0_rp  !0.999999990696210      !ff=26);
            vaulo_exm(18, 1:2) = 1.0_rp  !0.904906458666787      !fs=27);
            vaulo_exm(19, 1:2) = 1.0_rp  !0.999999990696060      !fcaf=28);
            vaulo_exm(20, 1:2) = 1.0_rp  !0.999581201974281      !fcas=29);
            vaulo_exm(21, 1:2) = 1.0_rp  !0.999903346883777      !jca=30);
            vaulo_exm(22, 1:2) = 0.0_rp !0.00215555277945401      !nca=31);
            vaulo_exm(23, 1:2) = 1.0_rp  !0.999999990680285      !ffp=32);
            vaulo_exm(24, 1:2) = 1.0_rp  !0.999999990692529      !fcafp=33);
            vaulo_exm(25, 1:2) = 0.0_rp !8.64222375034682e-06      !xrf=34);
            vaulo_exm(26, 1:2) = 0.0_rp !0.487585264457487      !xrs=35);
            vaulo_exm(27, 1:2) = 0.0_rp !0.276203479404767      !xs1=36);
            vaulo_exm(28, 1:2) = 0.0_rp !0.000194412216700766      !xs2=37);
            vaulo_exm(29, 1:2) = 1.0_rp  !0.996778581263402      !xk1=38);
         case (EXM_CELLTYPE_EPI)
            !epicardial
            voltage_local = -87.5_rp     ! INTRA (BI-DOMAIN) OR SINGLE (MONO) POTENTIAL
            !voltage_local = -87.99_rp
            viclo_exm(1:26) = 0.0_rp
            vcolo_exm(1, 1:2) = 0.0001_rp  !8.03149767106260e-05       !cai=6)
            vcolo_exm(2, 1:2) = 7.0_rp !8.17212071948942       !nass=3);
            vcolo_exm(3, 1:2) = 145.0_rp !143.675184333376       !ki=4);
            vcolo_exm(4, 1:2) = 145.0_rp  !143.675147146842       !kss=5);
            vcolo_exm(5, 1:2) = 7.0_rp  !8.17202753912090      ! nai=2
            vcolo_exm(6, 1:2) = 0.0001_rp  !7.96912949697674e-05       !cass=7);
            vcolo_exm(7, 1:2) = 1.2_rp  !2.14811455007091       !cansr=8);
            vcolo_exm(8, 1:2) = 1.2_rp  !2.03335899236798      !cajsr=9);
            vcolo_exm(9, 1:2) = 0.0_rp !4.74946280300893e-07       !Jrelnp=39);
            vcolo_exm(10, 1:2) = 0.0_rp !5.93539009244893e-07      !Jrelp=40);
            vcolo_exm(11, 1:2) = 0.0_rp !0.0228529042639590      !CaMKt=41);
            vaulo_exm(1, 1:2) = 0.0_rp !0.00739719746920272       !m=10);
            vaulo_exm(2, 1:2) = 1.0_rp  !0.695621622011335       !hf=11);
            vaulo_exm(3, 1:2) = 1.0_rp  !0.695601842634086       !hs=12);
            vaulo_exm(4, 1:2) = 1.0_rp  !0.695486248719023       !j=13);
            vaulo_exm(5, 1:2) = 1.0_rp  !0.452023628358454       !hsp=14);
            vaulo_exm(6, 1:2) = 1.0_rp  !0.695403157533235       !jp=15);
            vaulo_exm(7, 1:2) = 0.0_rp !0.000190839777466418       !mL=16);
            vaulo_exm(8, 1:2) = 1.0_rp  !0.493606704642336       !hL=17);
            vaulo_exm(9, 1:2) = 1.0_rp  !0.264304293390731       !hLp=18);
            vaulo_exm(10, 1:2) = 0.0_rp !0.00100594231451985      !a=19);
            vaulo_exm(11, 1:2) = 1.0_rp  !0.999548606668578      !iF=20);
            vaulo_exm(12, 1:2) = 1.0_rp  !0.999488774162635      !iS=21);
            vaulo_exm(13, 1:2) = 0.0_rp !0.000512555980943569      !ap=22);
            vaulo_exm(14, 1:2) = 1.0_rp  !0.999548607287668      !iFp=23);
            vaulo_exm(15, 1:2) = 1.0_rp  !0.999488774162635      !iSp=24);
            vaulo_exm(16, 1:2) = 0.0_rp !2.38076098345898e-09      !d=25);
            vaulo_exm(17, 1:2) = 1.0_rp  !0.999999990696210      !ff=26);
            vaulo_exm(18, 1:2) = 1.0_rp  !0.904906458666787      !fs=27);
            vaulo_exm(19, 1:2) = 1.0_rp  !0.999999990696060      !fcaf=28);
            vaulo_exm(20, 1:2) = 1.0_rp  !0.999581201974281      !fcas=29);
            vaulo_exm(21, 1:2) = 1.0_rp  !0.999903346883777      !jca=30);
            vaulo_exm(22, 1:2) = 0.0_rp !0.00215555277945401      !nca=31);
            vaulo_exm(23, 1:2) = 1.0_rp  !0.999999990680285      !ffp=32);
            vaulo_exm(24, 1:2) = 1.0_rp  !0.999999990692529      !fcafp=33);
            vaulo_exm(25, 1:2) = 0.0_rp !8.64222375034682e-06      !xrf=34);
            vaulo_exm(26, 1:2) = 0.0_rp !0.487585264457487      !xrs=35);
            vaulo_exm(27, 1:2) = 0.0_rp !0.276203479404767      !xs1=36);
            vaulo_exm(28, 1:2) = 0.0_rp !0.000194412216700766      !xs2=37);
            vaulo_exm(29, 1:2) = 1.0_rp  !0.996778581263402      !xk1=38);
         case (EXM_CELLTYPE_MID)
            !MIDmyocardial
            voltage_local = -87.5_rp     ! INTRA (BI-DOMAIN) OR SINGLE (MONO) POTENTIAL
            !voltage_local = -87.99_rp
            viclo_exm(1:26) = 0.0_rp
            vcolo_exm(1, 1:2) = 0.0001_rp  !8.03149767106260e-05       !cai=6)
            vcolo_exm(2, 1:2) = 7.0_rp !8.17212071948942       !nass=3);
            vcolo_exm(3, 1:2) = 145.0_rp !143.675184333376       !ki=4);
            vcolo_exm(4, 1:2) = 145.0_rp  !143.675147146842       !kss=5);
            vcolo_exm(5, 1:2) = 7.0_rp  !8.17202753912090      ! nai=2
            vcolo_exm(6, 1:2) = 0.0001_rp  !7.96912949697674e-05       !cass=7);
            vcolo_exm(7, 1:2) = 1.2_rp  !2.14811455007091       !cansr=8);
            vcolo_exm(8, 1:2) = 1.2_rp  !2.03335899236798      !cajsr=9);
            vcolo_exm(9, 1:2) = 0.0_rp !4.74946280300893e-07       !Jrelnp=39);
            vcolo_exm(10, 1:2) = 0.0_rp !5.93539009244893e-07      !Jrelp=40);
            vcolo_exm(11, 1:2) = 0.0_rp !0.0228529042639590      !CaMKt=41);
            vaulo_exm(1, 1:2) = 0.0_rp !0.00739719746920272       !m=10);
            vaulo_exm(2, 1:2) = 1.0_rp  !0.695621622011335       !hf=11);
            vaulo_exm(3, 1:2) = 1.0_rp  !0.695601842634086       !hs=12);
            vaulo_exm(4, 1:2) = 1.0_rp  !0.695486248719023       !j=13);
            vaulo_exm(5, 1:2) = 1.0_rp  !0.452023628358454       !hsp=14);
            vaulo_exm(6, 1:2) = 1.0_rp  !0.695403157533235       !jp=15);
            vaulo_exm(7, 1:2) = 0.0_rp !0.000190839777466418       !mL=16);
            vaulo_exm(8, 1:2) = 1.0_rp  !0.493606704642336       !hL=17);
            vaulo_exm(9, 1:2) = 1.0_rp  !0.264304293390731       !hLp=18);
            vaulo_exm(10, 1:2) = 0.0_rp !0.00100594231451985      !a=19);
            vaulo_exm(11, 1:2) = 1.0_rp  !0.999548606668578      !iF=20);
            vaulo_exm(12, 1:2) = 1.0_rp  !0.999488774162635      !iS=21);
            vaulo_exm(13, 1:2) = 0.0_rp !0.000512555980943569      !ap=22);
            vaulo_exm(14, 1:2) = 1.0_rp  !0.999548607287668      !iFp=23);
            vaulo_exm(15, 1:2) = 1.0_rp  !0.999488774162635      !iSp=24);
            vaulo_exm(16, 1:2) = 0.0_rp !2.38076098345898e-09      !d=25);
            vaulo_exm(17, 1:2) = 1.0_rp  !0.999999990696210      !ff=26);
            vaulo_exm(18, 1:2) = 1.0_rp  !0.904906458666787      !fs=27);
            vaulo_exm(19, 1:2) = 1.0_rp  !0.999999990696060      !fcaf=28);
            vaulo_exm(20, 1:2) = 1.0_rp  !0.999581201974281      !fcas=29);
            vaulo_exm(21, 1:2) = 1.0_rp  !0.999903346883777      !jca=30);
            vaulo_exm(22, 1:2) = 0.0_rp !0.00215555277945401      !nca=31);
            vaulo_exm(23, 1:2) = 1.0_rp  !0.999999990680285      !ffp=32);
            vaulo_exm(24, 1:2) = 1.0_rp  !0.999999990692529      !fcafp=33);
            vaulo_exm(25, 1:2) = 0.0_rp !8.64222375034682e-06      !xrf=34);
            vaulo_exm(26, 1:2) = 0.0_rp !0.487585264457487      !xrs=35);
            vaulo_exm(27, 1:2) = 0.0_rp !0.276203479404767      !xs1=36);
            vaulo_exm(28, 1:2) = 0.0_rp !0.000194412216700766      !xs2=37);
            vaulo_exm(29, 1:2) = 1.0_rp  !0.996778581263402      !xk1=38);
         case default
            call runend("EXM_INIT_VOLATGES: Undefined cell type.")
         end select

      case default
         call runend("EXM_INIT_VOLATGES: Undefined condition type.")
      end select
   end subroutine exm_init_voltages


!------------------------------------------------------------------------
!> @addtogroup Exmedi
!> @{
!> @file    exm_ohara_updateaux.f90
!> @author  Mariano Vazquez
!> @brief   Updating cell model variables
!> @date   16/NOV/1966
!> @details Updating cell model variables
!!> @} 
!!!-----------------------------------------------------------------------
!subroutine exm_ohara_updateaux(ipoin,inf_results,tau_results,elmag_local,ko,conc6,dtimeEP,ituss_exm, tauHL)
!   use def_kintyp, only : ip,rp
!   use def_exmedi
!   use mod_exm_oharaprecalc
!   use mod_exm_cellmodel
!   implicit none
!
!   integer(ip), intent(in) :: ipoin, ituss_exm !< node
!   real(rp), intent(in)  :: inf_results(EXM_OHR_INF_EQUATIONS) !< inf
!   real(rp), intent(in)  :: tau_results(EXM_OHR_TAU_EQUATIONS) !< tau
!   real(rp), intent(in)  :: elmag_local !< elmag
!   real(rp), intent(in)  :: ko !< ko
!   real(rp), intent(in)  :: conc6 !< conc6
!   real(rp), intent(in)  :: dtimeEP !< dtimeEP
!   real(rp), intent(in)  :: tauHL ! tauHL
!
!   integer(ip) :: iauxi
!   real(rp) :: infs(nauxi_exm),taus(nauxi_exm)
!   real(rp) :: anca,dnca,km2n,vaux1,vaux4
!
!   real(rp), parameter   :: KMN = 0.002_rp
!   real(rp), parameter   :: K2N = 1000.0_rp
!
!   infs(1) = inf_results(OHR_M_INF)
!   taus(1) = tau_results(OHR_TAU_M)
!
!   infs(2) = inf_results(OHR_H_INF)
!   taus(2) = tau_results(OHR_TAU_H_FAST)
!
!   infs(3) = inf_results(OHR_H_INF)
!   taus(3) = tau_results(OHR_TAU_H_SLOW)
!
!   infs(4) = inf_results(OHR_H_INF)  ! j_inf = h_inf
!   taus(4) = tau_results(OHR_TAU_J)
!
!   infs(5) = inf_results(OHR_H_CAMK_INF)
!   taus(5) = 3.0_rp * tau_results(OHR_TAU_H_SLOW)
!
!   infs(6) = inf_results(OHR_H_INF)  ! j_camk_inf = j_inf = h_inf
!   taus(6) = 1.46_rp * tau_results(OHR_TAU_J)
!
!   infs(7) = inf_results(OHR_M_L_INF)
!   taus(7) = tau_results(OHR_TAU_M)  ! tau_m_l = tau_m
!
!   infs(8) = inf_results(OHR_H_L_INF)
!   taus(8) = tauHL
!
!   infs(9) = inf_results(OHR_H_L_CAMK_INF)
!   taus(9) = 3.0_rp * 200.0_rp
!
!   ! NOTE: OHR_INV_TAU_A_1/2 are infs, not taus!
!   infs(10) = inf_results(OHR_A_INF)
!   taus(10) = 1.0515_rp / (inf_results(OHR_INV_TAU_A_1) + inf_results(OHR_INV_TAU_A_2))
!
!   infs(11) = inf_results(OHR_I_INF)
!   taus(11) = tau_results(OHR_TAU_I_FAST)
!   if(ituss_exm == EXM_CELLTYPE_EPI) then
!      taus(11) = taus(11) * inf_results(OHR_DELTA_EPI)
!   end if
!
!   infs(12) = inf_results(OHR_I_INF)
!   taus(12) = tau_results(OHR_TAU_I_SLOW)
!   if(ituss_exm == EXM_CELLTYPE_EPI) then
!      taus(12) = taus(12) * inf_results(OHR_DELTA_EPI)
!   end if
!
!   infs(13) = inf_results(OHR_A_CAMK_INF)
!   taus(13) = taus(10) ! tau_a_camk_inf = tau_a
!
!   infs(14) = inf_results(OHR_I_INF) ! i_camk_inf = i_inf
!   taus(14) = tau_results(OHR_DELTA_CAMK_DEVELOP) * inf_results(OHR_DELTA_CAMK_RECOVER) * taus(11)
!
!   infs(15) = inf_results(OHR_I_INF) ! i_camk_inf = i_inf
!   taus(15) = tau_results(OHR_DELTA_CAMK_DEVELOP) * inf_results(OHR_DELTA_CAMK_RECOVER) * taus(12)
!
!   infs(16) = inf_results(OHR_D_INF)
!   taus(16) = tau_results(OHR_TAU_D)
!
!   infs(17) = inf_results(OHR_F_INF)
!   taus(17) = tau_results(OHR_TAU_F_FAST)
!
!   infs(18) = inf_results(OHR_F_INF)
!   taus(18) = tau_results(OHR_TAU_F_SLOW)
!
!   infs(19) = inf_results(OHR_F_INF) ! f_ca_inf = f_inf
!   taus(19) = tau_results(OHR_TAU_F_CA_FAST)
!
!   infs(20) = inf_results(OHR_F_INF) ! f_ca_inf = f_inf
!   taus(20) = tau_results(OHR_TAU_F_CA_SLOW)
!
!   infs(21) = inf_results(OHR_F_INF) ! j_ca_inf = f_ca_inf = f_inf
!   taus(21) = 75.0_rp
!
!   infs(23) = inf_results(OHR_F_INF) ! f_camk_inf = f_inf
!   taus(23) = 2.5_rp * tau_results(OHR_TAU_F_FAST)
!
!   infs(24) = inf_results(OHR_F_INF) ! f_ca_camk_inf = f_inf
!   taus(24) = 2.5_rp * tau_results(OHR_TAU_F_CA_FAST)
!
!   infs(25) = inf_results(OHR_X_R_INF)
!   taus(25) = tau_results(OHR_TAU_XR_FAST)
!
!   infs(26) = inf_results(OHR_X_R_INF)
!   taus(26) = tau_results(OHR_TAU_XR_SLOW)
!
!   infs(27) = inf_results(OHR_X_S1_INF)
!   taus(27) = tau_results(OHR_TAU_X_S1)
!
!   infs(28) = inf_results(OHR_X_S1_INF)
!   taus(28) = tau_results(OHR_TAU_X_S2)
!
!   infs(29) = 1.0_rp / (1.0_rp + safe_exp(-(elmag_local + 2.5538_rp * ko + 144.59_rp) / (1.5692_rp * ko + 3.8115_rp)))
!   taus(29) = tau_results(OHR_TAU_X_K1)
!
!   do iauxi=1,21
!      vauxi_exm(iauxi,ipoin,1) = infs(iauxi) - (infs(iauxi) - vauxi_exm(iauxi,ipoin,2)) * safe_exp(-dtimeEP/taus(iauxi))
!   end do
!   do iauxi=23,29
!      vauxi_exm(iauxi,ipoin,1) = infs(iauxi) - (infs(iauxi) - vauxi_exm(iauxi,ipoin,2)) * safe_exp(-dtimeEP/taus(iauxi))
!   end do
!        
!   !!!!! calculate icana and icak
!   km2n = vauxi_exm(21,ipoin,1)
!   vaux1 = KMN / conc6
!   vaux4 = (1.0_rp+vaux1) *(1.0_rp+vaux1) * (1.0_rp+vaux1) * (1.0_rp+vaux1)
!   anca = 1.0_rp / ((K2N / km2n) + vaux4)
!   dnca = (anca*k2n/km2n)
!   vauxi_exm(22,ipoin,1) = dnca - (dnca - vauxi_exm(22,ipoin,2)) * safe_exp(-dtimeEP*km2n)
!
!end subroutine exm_ohara_updateaux
!
!
!
!



end module mod_exm_ohararudy
