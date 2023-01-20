!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



module mod_chm_sectional_soot_model

  use def_kintyp, only   : ip,rp,r3p
  use def_domain, only   : ndime,mgaus,nelem

  implicit none

  !
  ! Soot constants
  !
  integer(ip), parameter                 ::       &
       n_soot_species = 8,                        & ! Number of species involved in soot calculations
       NCVmin = 32                                  ! Number of C atoms in Vmin

  real(rp), parameter                    ::       &
       pi = 3.141592653589793_rp,                 & ! pi
       MC = 12.0112_rp,                           & ! Atomic weight of Carbon (gm/mol)
       NAvo = 6.02214129e23_rp,                   & ! Avogadros number (1/mol)
       kBoltz = 1.38064852e-16_rp                   ! Boltzmann constsnt (cm2 g s-2 K-1)

  !
  ! Predetermined order of species involved in soot
  !
  integer(ip), parameter                 ::       &
       SSM_SPEC_POS_A4     = 1,                   &
       SSM_SPEC_POS_H      = 2,                   &
       SSM_SPEC_POS_H2     = 3,                   &
       SSM_SPEC_POS_OH     = 4,                   &
       SSM_SPEC_POS_O2     = 5,                   &
       SSM_SPEC_POS_H2O    = 6,                   &
       SSM_SPEC_POS_C2H2   = 7,                   &
       SSM_SPEC_POS_CO     = 8

  !
  ! Names associated with species involved in soot
  !
  character(Len=5)                       :: soot_species_names(n_soot_species)

  !
  ! Indeces associated with species involved in soot
  !
  integer(ip)                            :: soot_species_index(n_soot_species)

  !
  ! Molecular weights associated with soot species
  !
  real(rp)                               :: soot_species_W(n_soot_species)

  !
  ! Some precomputed constants
  !
  real(rp)                               :: &
      VC,                                   & ! volume associated with one carbon
      VPAH,                                 & ! volume of 1 PAH
      VC_NA,                                & ! volume associated with 1 Carbon times Avogadro number
      VPAH_NA                                 ! volume of 1 PAH times Avogadro number



  integer(ip)                           :: &
       gasCoupling_ssm,                    & ! Coupling with gas phase
       ID_NUCL,                            & ! Index nucleation process
       ID_COND,                            & ! Index condensation process
       ID_COAG,                            & ! Index coagulation process
       ID_SURF,                            & ! Index surface growth process
       idsm_0_ssm,                         & ! Index of last other unknown before sections
       idsm_first_ssm,                     & ! Index of first section in unknown vector
       idsm_last_ssm,                      & ! Index of last section in unknown vector
       nspec_ssm,                          & ! Number of species involved in surface growth
       nclas_ssm,                          & ! Total number of unknowns: Y_tot = Y_k + Y_s
       nsect_ssm,                          & ! Number of sections for soot sectional model
       kfl_spec_source_terms                 ! Wether or not account for species consumed by soot prcesses in the gas phase

  !
  ! Soot model variables
  !
  type(r3p),  pointer                   :: &
       Qsoot_gp_chm(:),                    & ! Soot volume fraction
       QsootDist_gp_chm(:),                & ! Soot volume fraction density
       NDsoot_gp_chm(:),                   & ! Number density
       Asoot_gp_chm(:),                    & ! Soot surface area
       Vsoot_gp_chm(:)                       ! Soot volume

  !
  ! Soot source terms
  !
  type(r3p),  pointer                   :: &
       Qnucl_gp_chm(:),                    & ! Nucleation
       Qcoag_gp_chm(:),                    & ! Coagulation
       Qcond_gp_chm(:),                    & ! Condensation
       Qsurf_gp_chm(:),                    & ! Surface growth
       Qtot_gp_chm(:)                        ! Total

  !
  ! Species source terms
  !
  type(r3p),  pointer                   :: &
       SSTnucl_gp_chm(:),                  & ! Nucleation
       SSTcond_gp_chm(:),                  & ! Condensation
       SSTsurf_gp_chm(:),                  & ! Surface growth and oxidation
       SSTtot_gp_chm(:)                      ! Total

  !
  ! Sectional soot model parameters
  !
  real(rp)                              :: &
       SootDensity_ssm,                    & ! Activation soot model
       RhoC_ssm,                           & ! Density of Carbon (gm/cm3)
       Vmax_ssm,                           & ! Index of soot sources models
       RadSoot_ssm                           ! Radiation parameter

  !
  ! Sectional variables
  !
  real(rp)                              :: &
       Vsec(100),                          & ! Sectional volume discretization
       Vpmean(100),                        & ! Mean particle volume in section
       Dpmean(100)                           ! Mean particle diameter
  !
  ! Unit conversion
  !
  real(rp), parameter                   :: &
       UnC = 1000.0_rp                       ! Unit conversion factor 1: CGS InOut, 1000: MKS InOut

  private

  !
  ! Public variables
  !
  public :: SootDensity_ssm
  public :: Vmax_ssm
  public :: RadSoot_ssm
  public :: idsm_0_ssm
  public :: idsm_first_ssm
  public :: idsm_last_ssm
  public :: nsect_ssm
  public :: nspec_ssm
  public :: nclas_ssm
  public :: gasCoupling_ssm
  public :: RhoC_ssm
  public :: ID_NUCL,ID_COND,ID_COAG,ID_SURF
  public :: Qnucl_gp_chm
  public :: Qcoag_gp_chm
  public :: Qcond_gp_chm
  public :: Qsurf_gp_chm
  public :: SSTtot_gp_chm
  public :: Qtot_gp_chm
  public :: Qsoot_gp_chm
  public :: NDsoot_gp_chm
  public :: QsootDist_gp_chm

  !
  ! Public subroutines
  !
  public :: assembly_soot_sources_ssm
  public :: initialize_constants_ssm
  public :: chm_initialization_ssm
  public :: chm_validation_test_ssm
  public :: chm_sectional_soot_memory_ssm
  public :: soot_volume_discretization_ssm
  public :: gather_soot_ssm
  public :: elmpre_soot_ssm

  contains

    subroutine chm_validation_test_ssm()

      use def_domain,   only : mgaus
      use mod_messages, only : messages_live

      implicit none

      integer(ip), parameter :: npts  = 200   ! Input from chem1d: number of points

      real(rp)               :: gpcon(mgaus,nclas_ssm)
      real(rp)               :: gptem(mgaus)
      real(rp)               :: gpden(mgaus)

      external               :: runend

      !
      ! Start soot sources calculation
      !
      call header()

      !
      ! Reading index input file
      !
      call initialize_constants_ssm()

      !
      ! Initialization soot source terms
      !
      call chm_initialization_ssm()

      !
      ! Reading chem1d input file
      !
      call chm_readChem1d_ssm(nclas_ssm,npts,gpden,gpcon,gptem)

      !
      ! Compute Volume discretization
      !    Vsec,Vpmean,Dpmean
      !
      call soot_volume_discretization_ssm()

      !
      ! Compute soot source terms
      !
      call chm_calculate_soot_sources_ssm(1_ip,npts,gpcon,gptem,gpden)

      !
      ! End soot source terms calculations
      !
      call footer()

      !
      ! Stop the calculation
      !
      call messages_live('SECTONAL METHOD VALIDATION TEST FINISHED.')
      call messages_live('INITIAL SOLUTION','END SECTION')
      call runend('O.K.!')

    end subroutine chm_validation_test_ssm

    subroutine gather_soot_ssm(pnode,lnods,elsoot)
      use def_master, only     :  conce
      implicit none
      integer(ip),           intent(in)  :: pnode
      integer(ip),           intent(in)  :: lnods(pnode)
      real(rp),              intent(out) :: elsoot(pnode,nsect_ssm)
      integer(ip)                        :: inode,ipoin

      elsoot   = 0.0_rp
      do inode=1,pnode
         ipoin=lnods(inode)
         elsoot(inode,1:nsect_ssm)   = conce(ipoin,idsm_first_ssm:idsm_last_ssm,1)
      enddo

    end subroutine gather_soot_ssm

    subroutine elmpre_soot_ssm(pnode,pgaus,gpsha,elsoot,gpsoot)
      implicit none
      integer(ip),          intent(in)  :: pnode
      integer(ip),          intent(in)  :: pgaus
      real(rp),             intent(in)  :: gpsha(pnode,pgaus)
      real(rp),             intent(in)  :: elsoot(pnode,nsect_ssm)
      real(rp),             intent(out) :: gpsoot(pgaus,nsect_ssm)

      integer(ip)              :: inode,igaus

      !
      ! Soot mass frctions
      !
      gpsoot = 0.0_rp
      do igaus = 1,pgaus
         do inode = 1,pnode
            gpsoot(igaus,1:nsect_ssm) = gpsoot(igaus,1:nsect_ssm) &
                                            + gpsha(inode,igaus) * elsoot(inode,1:nsect_ssm)
         end do
      end do

!!DMM      !
!!DMM      ! Divergence of thermophoretic velocity
!!DMM      !
!!DMM      do idsm=1,nsect_ssm
!!DMM         do igaus = 1,pgaus
!!DMM            do inode = 1,pnode
!!DMM               coeff(1:ndime) = elVtherm(1:ndime,inode) * elsoot(inode,idsm)
!!DMM               do idime = 1,ndime
!!DMM                  gpdivVt(igaus,idsm) = gpdivVt(igaus,idsm) + gpcar(idime,inode,igaus) * coeff(idime)
!!DMM               end do
!!DMM            end do
!!DMM         end do
!!DMM      end do

    end subroutine elmpre_soot_ssm


    subroutine assembly_soot_sources_ssm(ielem,pgaus,gpcon,gptem,gpden,gprhs,gpYk)
      implicit none

      integer(ip),        intent(in)     :: ielem
      integer(ip),        intent(in)     :: pgaus
      real(rp),           intent(in)     :: gpcon(pgaus,nclas_ssm,*)
      real(rp),           intent(in)     :: gptem(pgaus)
      real(rp),           intent(in)     :: gpden(pgaus)
      real(rp),           intent(inout)  :: gprhs(pgaus,nclas_ssm)
      real(rp), optional, intent(in)     :: gpYk(pgaus,n_soot_species) ! OPTIONAL INPUT: species mass fractions for soot processes

      integer(ip)                  :: idsm,iclas
      !
      ! Initialize soot source terms
      !
      Qnucl_gp_chm(ielem) % a = 0.0_rp
      Qcoag_gp_chm(ielem) % a = 0.0_rp
      Qcond_gp_chm(ielem) % a = 0.0_rp
      Qsurf_gp_chm(ielem) % a = 0.0_rp
      Qtot_gp_chm(ielem)  % a = 0.0_rp

      if (kfl_spec_source_terms /= 0) then
         SSTnucl_gp_chm(ielem) % a = 0.0_rp
         SSTcond_gp_chm(ielem) % a = 0.0_rp
         SSTsurf_gp_chm(ielem) % a = 0.0_rp
         SSTtot_gp_chm(ielem)  % a = 0.0_rp
      endif

      !
      ! Calculate soot source terms
      !
      if (present(gpYk)) then
         call chm_calculate_soot_sources_ssm(ielem,pgaus,gpcon,gptem,gpden,gpYk)
      else
         call chm_calculate_soot_sources_ssm(ielem,pgaus,gpcon,gptem,gpden)
      endif

      !
      ! Impose RHS sectional soot terms
      !
      do idsm=1,nsect_ssm
         iclas = idsm + idsm_0_ssm
         gprhs(1:pgaus,iclas) = Qtot_gp_chm(ielem) % a (1:pgaus,idsm,1)
      end do

      !
      ! Impose RHS species soot terms
      !
      if (kfl_spec_source_terms /= 0) then
         do iclas=1,nspec_ssm
            gprhs(1:pgaus,iclas) = SSTtot_gp_chm(ielem) % a (1:pgaus,iclas, 1)
         end do
      endif

    end subroutine assembly_soot_sources_ssm




   subroutine chm_sectional_soot_memory_ssm()
      use def_domain, only : ltype,ngaus
      use mod_memory, only : memory_alloca
      use def_master, only : mem_modul,modul
      implicit none

      integer(ip)     :: ielem,pgaus,pelty

      !
      ! Soot variables
      !
      call memory_alloca(mem_modul(1:2,modul),'QSOOT_GP_CHM','chm_solmem',    Qsoot_gp_chm,nelem)
      call memory_alloca(mem_modul(1:2,modul),'QSOOTDIST_GP_CHM','chm_solmem',QsootDist_gp_chm,nelem)
      call memory_alloca(mem_modul(1:2,modul),'NDSOOT_GP_CHM','chm_solmem',   NDsoot_gp_chm,nelem)
      call memory_alloca(mem_modul(1:2,modul),'ASOOT_GP_CHM','chm_solmem',    Asoot_gp_chm,nelem)
      call memory_alloca(mem_modul(1:2,modul),'VSOOT_GP_CHM','chm_solmem',    Vsoot_gp_chm,nelem)

      !
      ! Soot source terms
      !
      call memory_alloca(mem_modul(1:2,modul),'QNUCL_GP_CHM','chm_solmem',Qnucl_gp_chm,nelem)
      call memory_alloca(mem_modul(1:2,modul),'QCOAG_GP_CHM','chm_solmem',Qcoag_gp_chm,nelem)
      call memory_alloca(mem_modul(1:2,modul),'QCOND_GP_CHM','chm_solmem',Qcond_gp_chm,nelem)
      call memory_alloca(mem_modul(1:2,modul),'QSURF_GP_CHM','chm_solmem',Qsurf_gp_chm,nelem)
      call memory_alloca(mem_modul(1:2,modul),'QTOT_GP_CHM' ,'chm_solmem', Qtot_gp_chm,nelem)

      !
      ! Species source terms
      !
      call memory_alloca(mem_modul(1:2,modul),'SSTNUCL_GP_CHM','chm_solmem',SSTnucl_gp_chm,nelem)
      call memory_alloca(mem_modul(1:2,modul),'SSTCOND_GP_CHM','chm_solmem',SSTcond_gp_chm,nelem)
      call memory_alloca(mem_modul(1:2,modul),'SSTSURF_GP_CHM','chm_solmem',SSTsurf_gp_chm,nelem)
      call memory_alloca(mem_modul(1:2,modul),'SSTTOT_GP_CHM','chm_solmem' , SSTtot_gp_chm,nelem)

      do ielem = 1,nelem
         pelty = ltype(ielem)
         pgaus = ngaus(pelty)

         !
         ! Soot model variables
         !
         call memory_alloca(mem_modul(1:2,modul),'QSOOT_GP_CHM % A','chm_solmem',Qsoot_gp_chm(ielem)%a,pgaus,nsect_ssm,1_ip)
         call memory_alloca(mem_modul(1:2,modul),'QSOOTDIST_GP_CHM % A','chm_solmem',QsootDist_gp_chm(ielem)%a,pgaus,nsect_ssm,1_ip)
         call memory_alloca(mem_modul(1:2,modul),'NDSOOT_GP_CHM % A','chm_solmem',NDsoot_gp_chm(ielem)%a,pgaus,nsect_ssm,1_ip)
         call memory_alloca(mem_modul(1:2,modul),'ASOOT_GP_CHM % A','chm_solmem',Asoot_gp_chm(ielem)%a,pgaus,nsect_ssm,1_ip)
         call memory_alloca(mem_modul(1:2,modul),'VSOOT_GP_CHM % A','chm_solmem',Vsoot_gp_chm(ielem)%a,pgaus,nsect_ssm,1_ip)

         !
         ! Soot source terms
         !
         call memory_alloca(mem_modul(1:2,modul),'QNUCL_GP_CHM % A','chm_solmem',Qnucl_gp_chm(ielem)%a,pgaus,nsect_ssm,1_ip)
         call memory_alloca(mem_modul(1:2,modul),'QCOAG_GP_CHM % A','chm_solmem',Qcoag_gp_chm(ielem)%a,pgaus,nsect_ssm,1_ip)
         call memory_alloca(mem_modul(1:2,modul),'QCOND_GP_CHM % A','chm_solmem',Qcond_gp_chm(ielem)%a,pgaus,nsect_ssm,1_ip)
         call memory_alloca(mem_modul(1:2,modul),'QSURF_GP_CHM % A','chm_solmem',Qsurf_gp_chm(ielem)%a,pgaus,nsect_ssm,1_ip)
         call memory_alloca(mem_modul(1:2,modul),'QTOT_GP_CHM % A' ,'chm_solmem', Qtot_gp_chm(ielem)%a,pgaus,nsect_ssm,1_ip)

         !
         ! Species source terms
         !
         call memory_alloca(mem_modul(1:2,modul),'SSTNUCL_GP_CHM % A','chm_solmem',SSTnucl_gp_chm(ielem)%a,pgaus,nspec_ssm,1_ip)
         call memory_alloca(mem_modul(1:2,modul),'SSTCOND_GP_CHM % A','chm_solmem',SSTcond_gp_chm(ielem)%a,pgaus,nspec_ssm,1_ip)
         call memory_alloca(mem_modul(1:2,modul),'SSTSURF_GP_CHM % A','chm_solmem',SSTsurf_gp_chm(ielem)%a,pgaus,nspec_ssm,1_ip)
         call memory_alloca(mem_modul(1:2,modul),'SSTTOT_GP_CHM % A','chm_solmem' ,SSTtot_gp_chm(ielem)%a,pgaus,nspec_ssm,1_ip)

      end do

    end subroutine chm_sectional_soot_memory_ssm


    subroutine chm_calculate_soot_sources_ssm(ielem,pgaus,gpcon,gptem,gpden,gpYk)
      !-----------------------------------------------------------------------
      !****f* Chemic/
      ! NAME
      !
      !
      ! DESCRIPTION
      !
      ! USES
      ! USED BY
      !
      !***
      !-----------------------------------------------------------------------
      use def_chemic, only   : ittot_chm

      implicit none
      integer(ip),        intent(in)  :: ielem
      integer(ip),        intent(in)  :: pgaus
      real(rp),           intent(in)  :: gpcon(pgaus,nclas_ssm)     ! Mass fraction species Yk
      real(rp),           intent(in)  :: gptem(pgaus)               ! Temperature
      real(rp),           intent(in)  :: gpden(pgaus)               ! Density
      real(rp), optional, intent(in)  :: gpYk(pgaus,n_soot_species) ! OPTIONAL INPUT: species mass fractions for soot processes

      real(rp)                        :: gpsoot(pgaus,nsect_ssm)    ! Mass fraction soot species Ys
      real(rp)                        :: gpXk(pgaus,n_soot_species) ! Concentration
      !
      ! Compute sectional soot mass fractions
      !
      gpsoot(1:pgaus,1:nsect_ssm) = gpcon(1:pgaus,idsm_first_ssm:idsm_last_ssm)

      !
      ! Compute sectional soot model variables
      !
      Qsoot_gp_chm(ielem)%a(1:pgaus,1:nsect_ssm,1_ip)     = soot_vol_frac( pgaus,gpden(1:pgaus), &
                                                                           gpsoot(1:pgaus,1:nsect_ssm) )

      QsootDist_gp_chm(ielem)%a(1:pgaus,1:nsect_ssm,1_ip) = soot_vol_frac_dist( pgaus, &
                                                                          Qsoot_gp_chm(ielem)%a(1:pgaus,1:nsect_ssm,1_ip) )

      NDsoot_gp_chm(ielem)%a(1:pgaus,1:nsect_ssm,1_ip)    = soot_num_den( pgaus, &
                                                                          QsootDist_gp_chm(ielem)%a(1:pgaus,1:nsect_ssm,1_ip) )

      Asoot_gp_chm(ielem)%a(1:pgaus,1:nsect_ssm,1_ip)     = soot_surf_area( pgaus, &
                                                                            QsootDist_gp_chm(ielem)%a(1:pgaus,1:nsect_ssm,1_ip) )

      Vsoot_gp_chm(ielem)%a(1:pgaus,1:nsect_ssm,1_ip)     = soot_vol( pgaus, &
                                                                      QsootDist_gp_chm(ielem)%a(1:pgaus,1:nsect_ssm,1_ip), &
                                                                      NDsoot_gp_chm(ielem)%a(1:pgaus,1:nsect_ssm,1_ip) )

      !
      ! Compute mass concentration of ONLY species involved in soot processes
      !
      if (present(gpYk)) then
         gpXk = spec_con(pgaus,nclas_ssm,gpden,gpYk)
      else
         gpXk = spec_con(pgaus,nclas_ssm,gpden,gpcon)
      endif


      !
      ! Compute source term nucleation
      !
      if (ID_NUCL > 0) &
         call chm_source_nucleation_ssm(   pgaus,gpXk(1:pgaus,1:n_soot_species),gptem(1:pgaus), &
                                           Qnucl_gp_chm(ielem) % a(1:pgaus,1:nsect_ssm,1), &
                                           SSTnucl_gp_chm(ielem) % a(1:pgaus,1:nspec_ssm,1) )

      !
      ! Compute source term coagulation
      !
      if (ID_COAG > 0) &
         call chm_source_coagulation_ssm(  pgaus,gptem(1:pgaus), &
                                           NDsoot_gp_chm(ielem) % a(1:pgaus,1:nsect_ssm,1) , &
                                           Qcoag_gp_chm(ielem) % a(1:pgaus,1:nsect_ssm,1) )

      !
      ! Compute source term condensation
      !
      if (ID_COND > 0) &
         call chm_source_condensation_ssm( pgaus,gpXk(1:pgaus,1:n_soot_species),gptem(1:pgaus), &
                                           NDsoot_gp_chm(ielem) % a(1:pgaus,1:nsect_ssm,1) , &
                                           Qcond_gp_chm(ielem) % a(1:pgaus,1:nsect_ssm,1), &
                                           SSTcond_gp_chm(ielem) % a(1:pgaus,1:nspec_ssm,1) )

      !
      ! Compute source term surface growth
      !
      if (ID_SURF > 0) &
         call chm_source_surfgrowth_ssm(   pgaus,gpXk(1:pgaus,1:n_soot_species),gptem(1:pgaus), &
                                           NDsoot_gp_chm(ielem) % a(1:pgaus,1:nsect_ssm,1), &
                                           Asoot_gp_chm(ielem) % a(1:pgaus,1:nsect_ssm,1), &
                                           Vsoot_gp_chm(ielem) % a(1:pgaus,1:nsect_ssm,1), &
                                           Qsurf_gp_chm(ielem) % a(1:pgaus,1:nsect_ssm,1), &
                                           SSTsurf_gp_chm(ielem) % a(1:pgaus,1:nspec_ssm,1) )

      !
      ! Compute total source terms
      !
      ! Adjustment to avoid unphysical increment in coagulation source term at initial steps

      if (ittot_chm < 10_ip) then
         Qtot_gp_chm(ielem) % a (1:pgaus,1:nsect_ssm, 1) = RhoC_ssm * ( Qnucl_gp_chm(ielem) % a (1:pgaus,1:nsect_ssm, 1) &
                                                      + Qcond_gp_chm(ielem) % a (1:pgaus,1:nsect_ssm, 1) &
                                                      + Qsurf_gp_chm(ielem) % a (1:pgaus,1:nsect_ssm, 1) )
      else
         Qtot_gp_chm(ielem) % a (1:pgaus,1:nsect_ssm, 1) = RhoC_ssm * ( Qnucl_gp_chm(ielem) % a (1:pgaus,1:nsect_ssm, 1) &
                                                      + Qcond_gp_chm(ielem) % a (1:pgaus,1:nsect_ssm, 1) &
                                                      + Qsurf_gp_chm(ielem) % a (1:pgaus,1:nsect_ssm, 1) &
                                                      + Qcoag_gp_chm(ielem) % a (1:pgaus,1:nsect_ssm, 1) )
      end if

      SSTtot_gp_chm(ielem) % a (1:pgaus,1:nspec_ssm, 1) = UnC * ( SSTnucl_gp_chm(ielem) % a (1:pgaus,1:nspec_ssm, 1) &
                                                        + SSTcond_gp_chm(ielem) % a (1:pgaus,1:nspec_ssm, 1) &
                                                        + SSTsurf_gp_chm(ielem) % a (1:pgaus,1:nspec_ssm, 1) )

    end subroutine chm_calculate_soot_sources_ssm

    subroutine chm_initialization_ssm()
      !-----------------------------------------------------------------------
      !****f* Chemic/
      ! NAME
      !
      !
      ! DESCRIPTION
      !
      ! USES
      ! USED BY
      !
      !***
      !-----------------------------------------------------------------------

      use def_domain, only   : nelem
      implicit none
      integer(ip)     :: ielem                ! Elements of the mesh

      !
      ! Initialization soot and species source terms
      !
      do ielem = 1,nelem
         Qnucl_gp_chm(ielem) % a = 0.0_rp
         Qcoag_gp_chm(ielem) % a = 0.0_rp
         Qcond_gp_chm(ielem) % a = 0.0_rp
         Qsurf_gp_chm(ielem) % a = 0.0_rp
         Qtot_gp_chm(ielem)  % a = 0.0_rp

         SSTnucl_gp_chm(ielem) % a = 0.0_rp
         SSTcond_gp_chm(ielem) % a = 0.0_rp
         SSTsurf_gp_chm(ielem) % a = 0.0_rp
         SSTtot_gp_chm(ielem)  % a = 0.0_rp
      end do

    end subroutine chm_initialization_ssm



    subroutine initialize_constants_ssm()

#ifdef CANTERA
       use cantera,          only : speciesIndex
       use def_chemic,       only : gas_chm
#endif
       use def_chemic,       only : kfl_model_chm, kfl_solve_cond_CMC_chm
       use def_chemic,       only : kfl_yk_fw_ssm
       use def_kermod,       only : lookup_fw
       use def_chemic,       only : mixedEq_groups_chm,ngrou_chm
       use def_chemic,       only : W_k
       use mod_chm_mixedEq,  only : CHM_GR_SECTIONAL

       implicit none

       integer(ip)                                 :: ii,jj
       character(len=:), allocatable               :: soot_species_trimmed

       !
       ! Initialize names
       !
       soot_species_names(SSM_SPEC_POS_A4)   = 'A4'
       soot_species_names(SSM_SPEC_POS_H)    = 'H'
       soot_species_names(SSM_SPEC_POS_H2)   = 'H2'
       soot_species_names(SSM_SPEC_POS_OH)   = 'OH'
       soot_species_names(SSM_SPEC_POS_O2)   = 'O2'
       soot_species_names(SSM_SPEC_POS_H2O)  = 'H2O'
       soot_species_names(SSM_SPEC_POS_C2H2) = 'C2H2'
       soot_species_names(SSM_SPEC_POS_CO)   = 'CO'

       !
       ! Initialize which index to use
       !
       soot_species_index = 0_ip
       if (kfl_model_chm == 3 .or. (kfl_model_chm == 4 .and. kfl_solve_cond_CMC_chm == 1)) then
          !
          ! Finite rate model + CMC model
          !
          do ii = 1,n_soot_species
             soot_species_trimmed   = trim(adjustl(soot_species_names(ii)))
#ifdef CANTERA
             soot_species_index(ii) = speciesIndex(gas_chm, soot_species_trimmed)
#endif
          end do
       elseif(kfl_model_chm == 2 .and. kfl_yk_fw_ssm > 0) then
          !
          ! Mixed equation model with looked up species
          !
          do ii = 1,n_soot_species
             soot_species_trimmed   = trim(adjustl(soot_species_names(ii)))
             !
             ! Look through all columns of the table, to find species indeces
             !
             do jj = 1,lookup_fw(kfl_yk_fw_ssm)%main_table%nvar
                if (soot_species_trimmed == lookup_fw(kfl_yk_fw_ssm)%main_table%varname(jj)) then
                    soot_species_index(ii) = jj
                endif
             enddo
          end do
       endif

       !
       ! Initialize whether equations are present for soot precursor species
       !
       kfl_spec_source_terms = 1_ip
       if (kfl_yk_fw_ssm > 0) kfl_spec_source_terms = 0_ip

       !
       ! Initialize molecular weight of species involved with soot
       !
       !
       ! Defaults:
       !
       soot_species_W(SSM_SPEC_POS_A4)   = 0.202256_rp
       soot_species_W(SSM_SPEC_POS_H)    = 0.001008_rp
       soot_species_W(SSM_SPEC_POS_H2)   = 0.002016_rp
       soot_species_W(SSM_SPEC_POS_OH)   = 0.017007_rp
       soot_species_W(SSM_SPEC_POS_O2)   = 0.031998_rp
       soot_species_W(SSM_SPEC_POS_H2O)  = 0.018015_rp
       soot_species_W(SSM_SPEC_POS_C2H2) = 0.026038_rp
       soot_species_W(SSM_SPEC_POS_CO)   = 0.028010_rp

       !
       ! If Cantera is used, update, in case there are slight changes:
       !
       if (kfl_model_chm == 3) then
          soot_species_W(SSM_SPEC_POS_A4)  =W_k(soot_species_index(SSM_SPEC_POS_A4))
          soot_species_W(SSM_SPEC_POS_H)   =W_k(soot_species_index(SSM_SPEC_POS_H))
          soot_species_W(SSM_SPEC_POS_H2)  =W_k(soot_species_index(SSM_SPEC_POS_H2))
          soot_species_W(SSM_SPEC_POS_OH)  =W_k(soot_species_index(SSM_SPEC_POS_OH))
          soot_species_W(SSM_SPEC_POS_O2)  =W_k(soot_species_index(SSM_SPEC_POS_O2))
          soot_species_W(SSM_SPEC_POS_H2O) =W_k(soot_species_index(SSM_SPEC_POS_H2O))
          soot_species_W(SSM_SPEC_POS_C2H2)=W_k(soot_species_index(SSM_SPEC_POS_C2H2))
          soot_species_W(SSM_SPEC_POS_CO)  =W_k(soot_species_index(SSM_SPEC_POS_CO))
       endif


       !
       ! Initialize variables locating sectional equations among all unknowns
       !
       if (kfl_model_chm == 3 .or. (kfl_model_chm == 4 .and. kfl_solve_cond_CMC_chm == 1)) then
          !
          ! Finite rate model: hardcoded
          !
          idsm_0_ssm      = nclas_ssm - nsect_ssm
          idsm_first_ssm  = nclas_ssm - nsect_ssm + 1
          idsm_last_ssm   = nclas_ssm
       elseif (kfl_model_chm == 2) then
          !
          ! Mixed equation model: find the group with the correct type
          ! Assuming only one of these groups exist
          !
          do ii = 1,ngrou_chm
             if (mixedEq_groups_chm(ii) % kfl_grtype == CHM_GR_SECTIONAL) then
                idsm_0_ssm      = mixedEq_groups_chm(ii) % i_start-1
                idsm_first_ssm  = mixedEq_groups_chm(ii) % i_start
                idsm_last_ssm   = mixedEq_groups_chm(ii) % i_end
             endif
          enddo
       endif

       !
       ! Initialize constants related to converting soot volumes:
       !
       VC_NA   = MC/((RhoC_ssm/UnC))
       VPAH_NA = 0.5_rp * real(NCVmin,rp) * VC_NA
       VC      = VC_NA/NAvo
       VPAH    = VPAH_NA/NAvo

    end subroutine initialize_constants_ssm

    subroutine chm_readChem1d_ssm(nvar, npts, Rho, Yspec, Temp)
      !-----------------------------------------------------------------------
      !****f* Chemic/
      ! NAME
      !
      !
      ! DESCRIPTION
      !
      ! USES
      ! USED BY
      !
      !***
      !-----------------------------------------------------------------------

      use def_master,      only : intost
      use iso_fortran_env, only : iostat_eor,iostat_end

      implicit none

      integer(ip),  intent(in)            :: nvar
      integer(ip),  intent(in)            :: npts
      real(rp),     intent(out)           :: Rho(npts)
      real(rp),     intent(out)           :: Yspec(npts,nvar)
      real(rp),     intent(out)           :: Temp(npts)
      integer(ip)                         :: ii
      real(rp), dimension(npts)           :: X
      integer(ip)                         :: ioerror

      external                            :: runend

      !
      ! Read the variables from chem1d
      !
      open (111, file = 'yiend_chem1d.dat', action = 'read', status = 'old')
      do ii = 1,npts
          read(111,*, iostat=ioerror) X(ii),Temp(ii),Rho(ii),Yspec(ii,:)
          if ((ioerror == iostat_eor).or.(ioerror == iostat_end)) ioerror = 0
          if (ioerror /= 0 ) call runend('chm_readChem1d_ssm: cannot read chem1d solution. error code:'//trim(intost(ioerror)))
      end do
      close(111)

      print*,'    Inputs read !'

    end subroutine chm_readChem1d_ssm


    subroutine chm_source_nucleation_ssm(pgaus,gpXk,gptem,Qnucl,SSTnucl)
      !-----------------------------------------------------------------------
      !****f* Chemic/
      ! NAME
      !
      !
      ! DESCRIPTION
      !
      ! USES
      ! USED BY
      !
      !***
      !-----------------------------------------------------------------------
      use def_chemic , only : kfl_soot_chm
      use def_master,  only : momod,modul
      implicit none

      integer(ip),  intent(in)  :: pgaus                      ! Vector size
      real(rp),     intent(in)  :: gpXk(pgaus,n_soot_species) ! Concentration
      real(rp),     intent(in)  :: gptem(pgaus)               ! Temperature
      real(rp),     intent(out) :: Qnucl(pgaus,nsect_ssm)     ! Source term nucleation
      real(rp),     intent(out) :: SSTnucl(pgaus,nspec_ssm)   ! Species source term from nucleation

      integer(ip)               :: k, idsm
      real(rp)                  :: Betaij_Nuc(pgaus)          ! Collison frequancy factor
      real(rp)                  :: XPAH(pgaus)                ! Number density of PAH

      Qnucl   = 0.0_rp
      SSTnucl = 0.0_rp

      XPAH  = max(gpXk(:,SSM_SPEC_POS_A4), 0.0_rp)

      !
      ! Compute collision frequencies for nucleation
      !
      call soot_collision_freq(pgaus,VPAH,VPAH,gptem,Betaij_Nuc)

      idsm = 1
      do k = 1,pgaus
         Qnucl(k,idsm) = 2.0_rp*VPAH_NA*NAvo*Betaij_Nuc(k)*(XPAH(k)**2.0_rp)
      end do

      !
      ! Coupling to gas phase (gm/cm3)
      !
      if (kfl_spec_source_terms /= 0) then
         do k = 1,pgaus
            SSTnucl(k,soot_species_index(SSM_SPEC_POS_A4))  = -1.0_rp*(Qnucl(k,idsm)/(VPAH_NA))*(soot_species_W(SSM_SPEC_POS_A4)*UnC)
            SSTnucl(k,soot_species_index(SSM_SPEC_POS_H2))  =  5.0_rp*(Qnucl(k,idsm)/(VPAH_NA))*(soot_species_W(SSM_SPEC_POS_H2)*UnC)
         end do
      endif

      !
      ! Writing nucleation source
      !
      if (kfl_soot_chm == -1) then
         call soot_write_output('Qnucl.dat',Qnucl,pgaus,nsect_ssm)
         write(momod(modul) % lun_outpu,*)'# Nucleation module finished. '
      end if

    end subroutine chm_source_nucleation_ssm

    subroutine chm_source_condensation_ssm(pgaus,gpXk,gptem,NDsoot,Qcond,SSTcond)
      !-----------------------------------------------------------------------
      !****f* Chemic/
      ! NAME
      !
      !
      ! DESCRIPTION
      !
      ! USES
      ! USED BY
      !
      !***
      !-----------------------------------------------------------------------
      use def_chemic , only : kfl_soot_chm
      use def_master,  only : momod,modul

      implicit none

      integer(ip),  intent(in)  :: pgaus                                ! Vector size
      real(rp),     intent(in)  :: gpXk(pgaus,n_soot_species)           ! Concentration
      real(rp),     intent(in)  :: gptem(pgaus)                         ! Temperature
      real(rp),     intent(in)  :: NDsoot(pgaus,nsect_ssm)              ! NDSoot
      real(rp),     intent(out) :: Qcond(pgaus,nsect_ssm)               ! Source term condensation
      real(rp),     intent(out) :: SSTcond(pgaus,nspec_ssm)             ! Species source term from condensation

      integer(ip)               :: k, idsm
      integer(ip)               :: IST_type                             ! Intersectional dynamics type
      real(rp)                  :: Betaij_Cond(pgaus)                   ! Collison frequancy factor
      real(rp)                  :: XPAH(pgaus)                          ! Molar density of PAH
      real(rp)                  :: iVol                                 ! volume of section i
      real(rp)                  :: Deltarate(pgaus,nsect_ssm,3)         ! Intersectional matrix

      IST_type = 1

      Qcond   = 0.0_rp
      SSTcond = 0.0_rp
      Deltarate = 0.0_rp

      XPAH = max(gpXk(:,SSM_SPEC_POS_A4), 0.0_rp)

      do  idsm = 1,nsect_ssm

          iVol = Vpmean(idsm)
          call soot_collision_freq(pgaus,iVol,VPAH,gptem,Betaij_Cond)

          do k = 1,pgaus
             Deltarate(k,idsm,1) = VPAH_NA*Betaij_Cond(k)*(XPAH(k)*NDsoot(k,idsm))
             Deltarate(k,idsm,1) = max(Deltarate(k,idsm,1),1.0e-32_rp)

             !
             ! Contribution to gas phase  (gm/cm3)
             !
             if (kfl_spec_source_terms /= 0) then
                SSTcond(k,soot_species_index(SSM_SPEC_POS_A4)) = SSTcond(k,soot_species_index(SSM_SPEC_POS_A4)) &
                                  -1.0_rp*(Deltarate(k,idsm,1)/(VPAH_NA))*(soot_species_W(SSM_SPEC_POS_A4)*UnC)
                SSTcond(k,soot_species_index(SSM_SPEC_POS_H2)) = SSTcond(k,soot_species_index(SSM_SPEC_POS_H2)) &
                                  +5.0_rp*(Deltarate(k,idsm,1)/(VPAH_NA))*(soot_species_W(SSM_SPEC_POS_H2)*UnC)
             endif
          end do
      end do

      call soot_intersectdynamics(pgaus,IST_type,Deltarate)

      do k = 1,pgaus
           do idsm = 1,nsect_ssm
                if (idsm == 1) then
                   Qcond(k,idsm) = Deltarate(k,idsm,2)
                else
                   Qcond(k,idsm) = Deltarate(k,idsm-1,3) + &
                                    Deltarate(k,idsm,2)
                end if
           end do
      end do

      !
      ! Writing condensation source
      !
      if (kfl_soot_chm == -1) then
         call soot_write_output('Qcond.dat',Qcond,pgaus,nsect_ssm)
      end if

      if (kfl_soot_chm == -1) &
         write(momod(modul) % lun_outpu,*)'# Condensation module finished. '

    end subroutine chm_source_condensation_ssm

    subroutine chm_source_coagulation_ssm(pgaus,gptem,NDSoot,Qcoag)
      !-----------------------------------------------------------------------
      !****f* Chemic/
      ! NAME
      !
      !
      ! DESCRIPTION
      !
      ! USES
      ! USED BY
      !
      !***
      !-----------------------------------------------------------------------
      use def_chemic , only : kfl_soot_chm
      use def_master,  only : momod,modul

      implicit none

      integer(ip),  intent(in)  :: pgaus                      ! Vector size
      real(rp),     intent(in)  :: gptem(pgaus)               ! Temperature
      real(rp),     intent(in)  :: NDsoot(pgaus,nsect_ssm)    ! NDSoot
      real(rp),     intent(out) :: Qcoag(pgaus,nsect_ssm)     ! Source term condensation

      integer(ip)      :: i, j, k, idsm, jsect, ksect, ivol
      integer(ip)      :: CoagCutSec
      real(rp)         :: Vol_i, Vol_j, Vol_k, VTot            ! volume of section
      real(rp)         :: rcexp,ctexp
      real(rp)         :: Epsilonij, kfm, Betaij_fm
      real(rp)         :: Beta_ij(nsect_ssm,nsect_ssm)
      real(rp)         :: VkLow, VkUp
      real(rp)         :: split_fact, kron_del, gain, loss
      integer(ip)      :: iCoagOpt

      iCoagOpt = 2_ip           ! Optimization flag (1: Optimized, 2: Standard)

      Qcoag    = 0.0_rp
      Beta_ij  = 0.0_rp

      !
      ! Optimization of coagulation model
      !
      if (iCoagOpt == 1_ip) then
          !
          ! Cutoff section
          !
          if (mod(nsect_ssm,2_ip) == 0_ip) then
             CoagCutSec = nsect_ssm / 2_ip
          else
             CoagCutSec = (nsect_ssm-1_ip)/2_ip
          end if

      else
          CoagCutSec   = nsect_ssm
      end if

!--------------------------------------------------------------

      ctexp= 1.0_rp/3.0_rp
      rcexp= 1.0_rp/2.0_rp

      !
      ! Beta constants for free-molecular regime
      ! van der Wall collision efficiency
      !

      Epsilonij = 2.2_rp

      !
      ! Constant for Free-Molecular Regime
      !
      do k = 1,pgaus

         kfm = Epsilonij*((3.0_rp/(4.0_rp*pi))**(1.0_rp/6.0_rp)) &
                        *((6.0_rp*kBoltz*gptem(k)/(RhoC_ssm/UnC))**(rcexp))

         !
         ! Precalculation of all particle collision: Beta_ij
         !
         do idsm = 1, CoagCutSec
            Vol_i = Vpmean(idsm)
            do jsect = 1, CoagCutSec
                Vol_j = Vpmean(jsect)
                !
                ! For general collisons in free-molecular regime!
                !
                Betaij_fm = kfm  &
                            *( Vol_i**ctexp + Vol_j**ctexp )**2.0_rp &
                            *( (1.0_rp/Vol_i + 1.0_rp/ Vol_j)**rcexp )

                Beta_ij(idsm,jsect) = Betaij_fm
            end do
         end do

         !
         ! Calculation of collisions
         !
         do idsm = 1,nsect_ssm
            ksect = idsm
            Vol_k = Vpmean(ksect)
            !
            ! Increase Source due to particles coagulation: Vi + Vj = Vk
            ! Coag. of sections i, j <= k, volume -> k
            !
            gain = 0.0_rp

            do i = 1,ksect
               Vol_i = Vpmean(i)
               do j = i,ksect
                  Vol_j = Vpmean(j)
                  if (ksect == 1) then
                     VkLow = Vpmean(ksect)
                  else
                     VkLow = Vpmean(ksect-1)
                  end if

                  if (ksect == nsect_ssm) then
                     VkUp = Vpmean(ksect)
                  else
                     VkUp = Vpmean(ksect+1)  !  (i+1)th section
                  end if

                  VTot = Vol_i + Vol_j

                  !
                  ! Calculation of the spliting factor
                  !

                  if (Vol_k <= VTot .and. VTot <= VkUp) then
                     split_fact = (VkUp - VTot) / (VkUp - Vol_k)
                  elseif (VkLow <= VTot .and. VTot <= Vol_k) then
                     split_fact = (VkLow - VTot) / (VkLow - Vol_k)
                  else
                     split_fact = 0.0_rp
                  end if

                  if (j == i) then
                     kron_del = 1.0_rp      ! Kronecker delta
                  else
                     kron_del = 0.0_rp
                  end if

                  gain = gain + (1.0_rp-0.5_rp*kron_del)*split_fact*Beta_ij(i,j)&
                                 *NDsoot(k,i)*NDsoot(k,j)
               end do
            end do

            !
            ! Decrease Source due to particles coagulation: Vl + Vp -> Vl
            ! Coag. of sections k, j => k, volume -< k
            !
            loss = 0.0_rp

            do i= 1,nsect_ssm
               loss = loss + Beta_ij(idsm,i) &
                       * NDsoot(k,idsm)*NDsoot(k,i)
            end do


            !
            ! Coagulation source term
            !
            if (idsm < nsect_ssm) then
               ivol = idsm+1
               Qcoag(k,idsm) = (gain - loss) * &
                       ( Vsec(ivol)-Vsec(ivol-1) ) /log( Vsec(ivol)/Vsec(ivol-1) )
            endif
         end do
      end do

      Qcoag(1:pgaus,nsect_ssm) = 0.0_rp

      !
      ! Writing coagulation source
      !
      if (kfl_soot_chm == -1) then
         call soot_write_output('Qcoag.dat',Qcoag,pgaus,nsect_ssm)
      end if

      if (kfl_soot_chm == -1) &
         write(momod(modul) % lun_outpu,*)'# Coagulation module finished. '

    end subroutine chm_source_coagulation_ssm

    subroutine chm_source_surfgrowth_ssm(pgaus,gpXk,gptem,NDsoot,Asoot,Vsoot,Qsurf,SSTsurf)
      !-----------------------------------------------------------------------
      !****f* Chemic/
      ! NAME
      !
      !
      ! DESCRIPTION
      !
      ! USES
      ! USED BY
      !
      !***
      !-----------------------------------------------------------------------
      use def_chemic , only : kfl_soot_chm
      use def_master,  only : momod,modul

      implicit none

      integer(ip),  intent(in)   :: pgaus                       ! Vector size
      real(rp),     intent(in)   :: gpXk(pgaus,n_soot_species)  ! Concentration
      real(rp),     intent(in)   :: gptem(pgaus)                ! Temperature
      real(rp),     intent(in)   :: NDsoot(pgaus,nsect_ssm), &  ! Soot number density
                                    Asoot(pgaus,nsect_ssm),  &  ! Soot surface area density
                                    Vsoot(pgaus,nsect_ssm)      ! Soot volume density
      real(rp),     intent(out)  :: Qsurf(pgaus,nsect_ssm)      ! Source term condensation
      real(rp),     intent(out)  :: SSTsurf(pgaus,nspec_ssm)    ! Species source term from surfce growth

      integer(ip)                :: k, idsm
      real(rp), dimension(pgaus) :: conH, conH2, conOH, conO2, conH2O, conC2H2
      real(rp), dimension(pgaus) :: fR1, rR1, fR2, rR2, fR3, fR4, fR5, fR6
      real(rp), dimension(pgaus) :: RT, denomD, denomC, Ydep, Ycon, Kss
      real(rp), dimension(pgaus) :: sR1, sR2, sR3, sR4, sR5, sR6
      real(rp)                   :: VolAdd_R1, VolAdd_R2, VolAdd_R3, &
                                    VolAdd_R4, VolAdd_R5, VolAdd_R6
      integer(ip)                :: IST_type_R1, IST_type_R2, IST_type_R3, &
                                    IST_type_R4, IST_type_R5, IST_type_R6   ! Intersectional dynamics type
      real(rp)                   :: kP, kR
      real(rp), dimension(pgaus) :: ProdR, ProdP,  &
                                    ProdH, ProdOH,ProdH2, ProdO2, &
                                    ProdCO, ProdH2O,ProdC2H2, ProdHR3, ProdOHR6 ! Production of gas-phase species

      real(rp), dimension(pgaus,nsect_ssm,3) :: Deltarate_sg,   &  ! Intersectional matrix Surf
                                                Deltarate_oxO2, &  ! Intersectional matrix OxO2
                                                Deltarate_oxOH, &  ! Intersectional matrix OxOH
                                                Dummyrate
      real(rp), dimension(pgaus,nsect_ssm)   :: surfSource, oxiSourceO2, oxiSourceOH

      surfSource  = 0.0_rp
      oxiSourceO2 = 0.0_rp
      oxiSourceOH = 0.0_rp
      Qsurf       = 0.0_rp
      SSTsurf     = 0.0_rp

      ProdR    = 0.0_rp    ! Dummy
      ProdP    = 0.0_rp

      ProdH    = 0.0_rp
      ProdOH   = 0.0_rp
      ProdH2   = 0.0_rp
      ProdO2   = 0.0_rp
      ProdCO   = 0.0_rp
      ProdH2O  = 0.0_rp
      ProdC2H2 = 0.0_rp
      ProdHR3  = 0.0_rp     ! Prod of H due to SR3
      ProdOHR6 = 0.0_rp     ! Prod of OH due to SR6

      Deltarate_sg   = 0.0_rp
      Deltarate_oxO2 = 0.0_rp
      Deltarate_oxOH = 0.0_rp
      Dummyrate      = 0.0_rp

      RT = 1.987e-3_rp *gptem    ! R in kcal/mol

      conH    = max(gpXk(:,SSM_SPEC_POS_H),   0.0_rp)
      conH2   = max(gpXk(:,SSM_SPEC_POS_H2),  0.0_rp)
      conOH   = max(gpXk(:,SSM_SPEC_POS_OH),  0.0_rp)
      conO2   = max(gpXk(:,SSM_SPEC_POS_O2),  0.0_rp)
      conH2O  = max(gpXk(:,SSM_SPEC_POS_H2O), 0.0_rp)
      conC2H2 = max(gpXk(:,SSM_SPEC_POS_C2H2),0.0_rp)

      ! Surface growth is currently described by the following scheme
      ! (Appel et al., 2000, Combust. Flame 121:122-136):
      !       1.  Csoot-H + H = Csoot* + H2         (fR1, rR1)
      !       2.  Csoot-H + OH = Csoot* + H2O       (fR2, rR2)
      !       3.  Csoot*  + H -> Csoot              (fR3)
      !       4.  Csoot*  + C2H2 -> Csoot-H + H     (fR4)
      !       5.  Csoot*  + O2 -> products          (fR5)
      !       6.  Csoot   + OH -> products          (fR6)
      !
      !     Surface reactions coefficients
      !     The following rate coefficients are from Appel et al. (2000)
      !     Ea in (kcal k/mol)

      !
      ! Calculating rates of surface reactions
      !
      fR1 = 4.20e13_rp *exp(-13.0_rp/RT)*conH
      rR1 = 3.90e12_rp *exp(-11.0_rp/RT)*conH2
      fR2 = 1.00e10_rp *(gptem**0.734_rp)*exp(-1.43_rp/RT)*conOH
      rR2 = 3.68e8_rp  *(gptem**1.139_rp)*exp(-17.1_rp/RT)*conH2O
      fR3 = 2.00e13_rp *conH
      fR4 = 8.00e7_rp  *(gptem**1.56_rp)*exp(-3.8_rp/RT)*conC2H2
      fR5 = 2.20e12_rp *exp(-7.50_rp/RT)*conO2
      fR6 = 0.13_rp * conOH                       ! gamma = 0.13 from Neoh et al.

      !
      ! Fraction of radical sites via steady-state:
      !
      denomD = rR1 + rR2 + fR3 + fR4 + fR5     ! For radical depletion
      denomC = rR1 + rR2 + fR3 + fR5           ! For radical conservation

      denomD = max(denomD,1.0e-40_rp)
      denomC = max(denomC,1.0e-40_rp)

      Ydep = (fR1 + fR2)/denomD
      Ycon = (fR1 + fR2)/denomC

      !
      ! Steady-State ratio
      !
      Kss = RadSoot_ssm*(Ycon-Ydep)+Ydep      ! (size pgaus)

      !
      ! For surface reaction-1----------------
      !
      VolAdd_R1 = 0.0_rp
      sR1 = fR1 - rR1 * Kss
      kR = -1.0_rp
      kP = 1.0_rp
      IST_type_R1 = 0

      call soot_surf_reac_rates( pgaus,IST_type_R1,'R1',VolAdd_R1,gptem,sR1, &
                                 Dummyrate,kR,ProdH,kP,ProdH2,NDsoot,Asoot,Vsoot )

      ProdH2 = ProdH2/2.0_rp

      !
      ! For surface reaction-2----------------
      !
      VolAdd_R2 = 0.0_rp
      sR2 = fR2 - rR2 * Kss
      kR = -1.0_rp
      kP = 1.0_rp
      IST_type_R2 = 0

      call soot_surf_reac_rates( pgaus,IST_type_R2,'R2',VolAdd_R2,gptem,sR2, &
                                 Dummyrate,kR,ProdOH,kP,ProdH2O,NDsoot,Asoot,Vsoot )

      ProdH2 = ProdH2 - ProdH2O/2.0_rp

      !
      ! For surface reaction-3----------------
      !
      VolAdd_R3 = 0.0_rp
      sR3 = fR3 * Kss
      kR = -1.0_rp
      kP = 0.0_rp
      IST_type_R3 = 0

      call soot_surf_reac_rates( pgaus,IST_type_R3,'R3',VolAdd_R3,gptem,sR3, &
                                 Dummyrate,kR,ProdHR3,kP,ProdR,NDsoot,Asoot,Vsoot )

      ProdH  = ProdH  + ProdHR3
      ProdH2 = ProdH2 - ProdHR3/2.0_rp

      !
      ! For surface growth reaction ----------------
      !
      VolAdd_R4 = 2.0_rp*VC
      sR4 = fR4* Kss
      kR = -1.0_rp
      kP = 1.0_rp
      IST_type_R4 = 1

      call soot_surf_reac_rates( pgaus,IST_type_R4,'R4',VolAdd_R4,gptem,sR4, &
                                 Deltarate_sg,kR,ProdC2H2,kP,ProdH,NDsoot,Asoot,Vsoot )

      ProdH2 = ProdH2 - ProdC2H2/2.0_rp

      !
      ! For surface OxO2 reaction ----------------
      !
      VolAdd_R5 = -2.0_rp*VC
      sR5 = fR5* Kss
      kR = -1.0_rp
      kP = 2.0_rp
      IST_type_R5 = 2

      call soot_surf_reac_rates( pgaus,IST_type_R5,'R5',VolAdd_R5,gptem,sR5, &
                                 Deltarate_oxO2,kR,ProdO2,kP,ProdCO,NDsoot,Asoot,Vsoot )

      !
      ! For surface OxOH reaction ----------------
      !
      VolAdd_R6 = -1.0_rp*VC
      sR6 = fR6* NAvo
      kR = -1.0_rp
      kP = 1.0_rp
      IST_type_R6 = 2

      call soot_surf_reac_rates( pgaus,IST_type_R6,'R6',VolAdd_R6,gptem,sR6, &
                                 Deltarate_oxOH,kR,ProdOHR6,kP,ProdCO,NDsoot,Asoot,Vsoot )

      ProdOH = ProdOH + ProdOHR6
      ProdH2 = ProdH2 - ProdOHR6/2.0_rp


    do k = 1,pgaus
        do  idsm = 1,nsect_ssm
            !
            ! Calculation of surface growth source term
            !
            if (idsm == 1) then
               surfSource(k,idsm) = Deltarate_sg(k,idsm,2)
            else
               surfSource(k,idsm) = Deltarate_sg(k,idsm-1,3) + &
                                     Deltarate_sg(k,idsm,2)
            end if

            !
            ! Calculation of OH oxidation source term
            !
            if (idsm == nsect_ssm) then
               oxiSourceO2(k,idsm) = Deltarate_oxO2(k,idsm,2)
            else
               oxiSourceO2(k,idsm) = Deltarate_oxO2(k,idsm+1,3) + &
                                      Deltarate_oxO2(k,idsm,2)
            end if

            !
            ! Calculation of OH oxidation source term
            !
            if (idsm == nsect_ssm) then
               oxiSourceOH(k,idsm) = Deltarate_oxOH(k,idsm,2)
            else
               oxiSourceOH(k,idsm) = Deltarate_oxOH(k,idsm+1,3) + &
                                      Deltarate_oxOH(k,idsm,2)
            end if

           Qsurf(k,idsm) = surfSource(k,idsm) + oxiSourceO2(k,idsm) + &
                            oxiSourceOH(k,idsm)
        end do

        !
        ! Calculation of source terms for species due to surface growth and
        ! soot oxidation (gm/cm3)
        !
        if (kfl_spec_source_terms /= 0) then
           SSTsurf(k,soot_species_index(SSM_SPEC_POS_H))    = ProdH(k)   *(soot_species_W(SSM_SPEC_POS_H)   *UnC)
           SSTsurf(k,soot_species_index(SSM_SPEC_POS_OH))   = ProdOH(k)  *(soot_species_W(SSM_SPEC_POS_OH)  *UnC)
           SSTsurf(k,soot_species_index(SSM_SPEC_POS_H2))   = ProdH2(k)  *(soot_species_W(SSM_SPEC_POS_H2)  *UnC)
           SSTsurf(k,soot_species_index(SSM_SPEC_POS_O2))   = ProdO2(k)  *(soot_species_W(SSM_SPEC_POS_O2)  *UnC)
           SSTsurf(k,soot_species_index(SSM_SPEC_POS_CO))   = ProdCO(k)  *(soot_species_W(SSM_SPEC_POS_CO)  *UnC)
           SSTsurf(k,soot_species_index(SSM_SPEC_POS_H2O))  = ProdH2O(k) *(soot_species_W(SSM_SPEC_POS_H2O) *UnC)
           SSTsurf(k,soot_species_index(SSM_SPEC_POS_C2H2)) = ProdC2H2(k)*(soot_species_W(SSM_SPEC_POS_C2H2)*UnC)
        endif
    end do

    !
    ! Writing surface growth source
    !
    if (kfl_soot_chm == -1) then
       call soot_write_output('Qsurf.dat',Qsurf,pgaus,nsect_ssm)
    end if

    if (kfl_soot_chm == -1) &
       write(momod(modul) % lun_outpu,*)'# Surface growth module finished. '

    end subroutine chm_source_surfgrowth_ssm



    subroutine soot_volume_discretization_ssm()

      use def_master,  only : INOTSLAVE
      use def_master,  only : momod,modul
      !
      ! This subroutine discretizes the PSDF into number of sections
      !
      implicit none
      integer(ip)         :: ii, jj
      real(rp)            :: Vmin

      !
      ! Initialize
      !
      Vsec   = 0.0_rp
      Vpmean = 0.0_rp
      Dpmean = 0.0_rp

      Vmax_ssm = Vmax_ssm*(UnC**2.0_rp)      ! Vmax in (cm3)
      Vmin = real(NCVmin,rp)*VC

      !
      ! Volume of the smallest particle
      !
      Vsec(1) = Vmin

      !
      ! Mean particle size is linearily discretized in the log scale
      ! Profile of Roy 2014, PhD
      !
      do ii = 1,nsect_ssm
         jj = ii+1
         Vsec(jj) = Vmin*(Vmax_ssm / Vmin)**(real(ii,rp)/real(nsect_ssm,rp))
      end do

      !
      ! Calculation of the mean particle size in the log scale
      !
      do ii = 1,nsect_ssm
         jj = ii+1
         Vpmean(ii) = exp( ( log(Vsec(jj))+log(Vsec(jj-1)) ) /2.0_rp)
         Dpmean(ii) = (6.0_rp*Vpmean(ii)/pi)**(1.0_rp/3.0_rp)
      enddo

      !
      ! Write volume discretization
      !
      if (INOTSLAVE) then
         write(momod(modul) % lun_outpu,*)'----------------------------------- '
         write(momod(modul) % lun_outpu,*)' '
         write(momod(modul) % lun_outpu,*)'# Discretization of the PSDF: '
         do ii=1,nsect_ssm+1
            write(momod(modul) % lun_outpu,*)'  Section', ii-1,' Vsec =',Vsec(ii), ' [cm^3]'
         enddo
         write(momod(modul) % lun_outpu,*)' '
         write(momod(modul) % lun_outpu,*)'----------------------------------- '
      endif

    end subroutine soot_volume_discretization_ssm


    function soot_vol_frac(pgaus, mix_den, soot_massfrac)

      ! This function converts soot mass fraction to volume fraction

      implicit none

      integer(ip), intent(in) :: pgaus
      real(rp),    intent(in) :: mix_den(pgaus)
      real(rp),    intent(in) :: soot_massfrac(pgaus,nsect_ssm)
      integer(ip)             :: ii,jj
      real(rp)                :: soot_vol_frac(pgaus,nsect_ssm)

      do ii=1,pgaus
          do jj=1,nsect_ssm
             soot_vol_frac(ii,jj) = max(soot_massfrac(ii,jj),0.0_rp)*(mix_den(ii)/RhoC_ssm)
          end do
      end do

    end function soot_vol_frac

    function soot_vol_frac_dist(pgaus, sootvolfrac)

      !
      ! This function converts soot mass fraction to volume fraction density
      !
      implicit none

      integer(ip), intent(in) :: pgaus
      real(rp),    intent(in) :: sootvolfrac(pgaus,nsect_ssm)
      integer(ip)             :: ii, jj, ivol
      real(rp)                :: soot_vol_frac_dist(pgaus,nsect_ssm)

      do ii=1,pgaus
          do jj=1,nsect_ssm
             ivol = jj+1
             soot_vol_frac_dist(ii,jj) = sootvolfrac(ii,jj)/(Vsec(ivol) - Vsec(ivol-1))
          end do
      end do

    end function soot_vol_frac_dist

    function soot_num_den(pgaus, svfdist)

      !
      ! This function calculates soot number density in each section
      !
      implicit none

      integer(ip), intent(in) :: pgaus
      real(rp),    intent(in) :: svfdist(pgaus,nsect_ssm)
      integer(ip)             :: ii, jj, ivol
      real(rp)                :: soot_num_den(pgaus,nsect_ssm)

      do ii=1,pgaus
          do jj=1,nsect_ssm
             ivol = jj+1
             soot_num_den(ii,jj) = svfdist(ii,jj)*log(Vsec(ivol)/ Vsec(ivol-1))
          end do
      end do

    end function soot_num_den

    function soot_surf_area(pgaus, svfdist)

      !
      ! This function calculates soot surface area in each section
      !
      implicit none

      integer(ip), intent(in) :: pgaus
      real(rp),    intent(in) :: svfdist(pgaus,nsect_ssm)
      integer(ip)             :: ii, jj, ivol
      real(rp)                :: soot_surf_area(pgaus,nsect_ssm)
      real(rp)                :: ctexp

      ctexp = 2.0_rp/3.0_rp

      do ii=1,pgaus
          do jj=1,nsect_ssm
             ivol = jj+1
             soot_surf_area(ii,jj) = pi*(6.0_rp/pi)**ctexp &
                                  * svfdist(ii,jj) * (1.0_rp/ctexp) &
                                  * ( Vsec(ivol)**ctexp - Vsec(ivol-1)**ctexp )
          end do
      end do

    end function soot_surf_area

    function soot_vol(pgaus, svfdist, numden)

      !
      ! This function calculates soot volume in each section
      !
      implicit none

      integer(ip), intent(in) :: pgaus
      real(rp),    intent(in) :: svfdist(pgaus,nsect_ssm)
      real(rp),    intent(in) :: numden(pgaus,nsect_ssm)
      integer(ip)             :: ii, jj, ivol
      real(rp)                :: soot_vol(pgaus,nsect_ssm)

      do ii=1,pgaus
          do jj=1,nsect_ssm
             ivol = jj+1
             soot_vol(ii,jj) = svfdist(ii,jj) *( Vsec(ivol) - Vsec(ivol-1) ) &
                               / (numden(ii,jj)+1.0_rp)
          end do
      end do

    end function soot_vol

    function spec_con(pgaus, nsp, mix_den, spec_massfrac )

      !
      ! This function converts mass fraction to molar concentration
      !
      implicit none

      integer(ip), intent(in) :: pgaus
      integer(ip), intent(in) :: nsp
      real(rp),    intent(in) :: mix_den(pgaus)
      real(rp),    intent(in) :: spec_massfrac(pgaus,nsp)
      real(rp)                :: spec_con(pgaus,n_soot_species)
      integer(ip)             :: ii

      do ii = 1,n_soot_species
         ! To convert from (kg/kg) to (g/cm3)
         spec_con(:,ii) = (mix_den/UnC)*(spec_massfrac(:,soot_species_index(ii))/(soot_species_W(ii)*UnC))
      end do

    end function spec_con

    subroutine soot_collision_freq(pgaus,Vi,Vj,T,Betaij)
      !-----------------------------------------------------------------------
      !****f* Chemic/
      ! NAME
      !
      !
      ! DESCRIPTION
      !
      ! USES
      ! USED BY
      !
      !***
      !-----------------------------------------------------------------------

      ! This subroutine calculates collision coefficients

      implicit none

      integer(ip), intent(in)  :: pgaus
      real(rp),    intent(in)  :: Vi,Vj
      real(rp),    intent(in)  :: T(pgaus)
      real(rp),    intent(out) :: Betaij(pgaus)

      real(rp) :: Epsilonij
      real(rp) :: ctexp,rcexp
      real(rp) :: kfm(pgaus)

      ctexp = 1.0_rp/3.0_rp
      rcexp = 1.0_rp/2.0_rp

      !
      ! Beta constants for free-molecular regime
      ! van der Wall collision efficiency
      !
      Epsilonij = 2.2_rp

      !
      ! Constant for Free-Molecular Regime
      !
      kfm = Epsilonij *((3.0_rp/(4.0_rp*pi))**(1.0_rp/6.0_rp)) &
                      *( (6.0_rp*kBoltz*T/(RhoC_ssm/UnC))**(rcexp) )

      !
      ! For general collisons in free-molecular regime
      !

      Betaij = kfm * ( Vi**ctexp + Vj**ctexp )**2.0_rp &
                   * ( (1.0_rp/Vi + 1.0_rp/ Vj)**rcexp )

    end subroutine soot_collision_freq

    subroutine soot_intersectdynamics(pgaus,ST_type,DeltaQ)
      !-----------------------------------------------------------------------
      !****f* Chemic/
      ! NAME
      !
      !
      ! DESCRIPTION
      !
      ! USES
      ! USED BY
      !
      !***
      !-----------------------------------------------------------------------

      ! This subroutine distributes soot mass over sections

      implicit none

      integer(ip), intent(in)    :: pgaus
      integer(ip), intent(in)    :: ST_type
      real(rp),    intent(inout) :: DeltaQ(pgaus,nsect_ssm,3_ip)
      integer(ip)                :: isec, ivol, imeshp
      real(rp)                   :: IBout, IBin

      !
      ! Initialize output parts
      !
      DeltaQ(1:pgaus,1:nsect_ssm,2:3) = 0.0_rp


      do imeshp=1,pgaus
          do isec=1,nsect_ssm
             ivol=isec+1

             !
             ! Inter-sectional balance for Surface-Growth
             !
             if (ST_type==1) then

                 !
                 ! Sections i < nsect_ssm-1
                 !
                 if (isec < nsect_ssm) then
                    !
                    ! Inter-sectional Balance In : Dq_in
                    !
                    IBin = ((Vsec(ivol+1)-Vsec(ivol))/(Vsec(ivol)-Vsec(ivol-1)))* &
                             (log(Vsec(ivol)/Vsec(ivol-1))/log(Vsec(ivol+1)/Vsec(ivol)))

                    DeltaQ(imeshp,isec,2) = DeltaQ(imeshp,isec,1)/(1.0_rp-IBin)

                    !
                    ! Inter-sectional Balance Out : Dq_out
                    !
                    IBout = ((Vsec(ivol)-Vsec(ivol-1))/(Vsec(ivol+1)-Vsec(ivol)))*&
                            (log(Vsec(ivol+1)/Vsec(ivol))/log(Vsec(ivol)/Vsec(ivol-1)))

                    DeltaQ(imeshp,isec,3) = DeltaQ(imeshp,isec,1)/(1.0_rp-IBout)

                 !
                 ! Last section
                 !
                 elseif  (isec == nsect_ssm) then
                    DeltaQ(imeshp,isec,2) = DeltaQ(imeshp,isec,1)
                    DeltaQ(imeshp,isec,3) = 0.0_rp
                 endif

             !
             ! Inter-sectional balance for Oxidation
             !
             elseif (ST_type==2) then

                 if (isec > 1) then
                    !
                    ! Intersectional Balance In
                    !
                    IBIn  = ((Vsec(ivol-1)-Vsec(ivol-2))/(Vsec(ivol)-Vsec(ivol-1)))*&
                            (log(Vsec(ivol)/Vsec(ivol-1))/log(Vsec(ivol-1)/Vsec(ivol-2)))

                    DeltaQ(imeshp,isec,2) = DeltaQ(imeshp,isec,1)/(1.0_rp-IBIn)

                    !
                    ! Intersectional Balance Out
                    !
                    IBOut = ((Vsec(ivol)-Vsec(ivol-1))/(Vsec(ivol-1)-Vsec(ivol-2)))*&
                            (log(Vsec(ivol-1)/Vsec(ivol-2))/log(Vsec(ivol)/Vsec(ivol-1)))

                    DeltaQ(imeshp,isec,3) = DeltaQ(imeshp,isec,1)/(1.0_rp-IBOut)

                 elseif (isec == 1) then
                    DeltaQ(imeshp,isec,2) = DeltaQ(imeshp,isec,1)
                    DeltaQ(imeshp,isec,3) = 0.0_rp
                 end if
             end if
          end do
      end do

    end subroutine soot_intersectdynamics

    subroutine soot_surf_reac_rates(pgaus,ST_type,RName,VolAdd,Temp,Rcoef,DeltaQ,stoiMR,ProdR,stoiMP,ProdP,NDsoot,Asoot,Vsoot)
      !-----------------------------------------------------------------------
      !
      ! DESCRIPTION : This subroutine calculates the contribution of surface reaction
      !               in production rates of species and soot source
      !
      !
      ! INPUTS:
      !         pgaus   = Vector size
      !         ST_type = Type of intersectional dynamics
      !         RName   = Name of surface reaction
      !         VolAdd  = Volume of C added due to surface reaction
      !         Temp    = Gas temperature (K)
      !         Rcoef   = Surface rate of reaction
      !         DeltaQ  = Surface rate with inter-sectional distribution (cm3/cm3)
      !         stoiMR  = stoichimetric coefficient for reactant
      !         ProdR   = Production rate of reactant species (mol/cm3 s)
      !         stoiMP  = stoichimetric coefficient for product
      !         ProdP   = Production rate of product species (mol/cm3 s)
      !         NDsoot  = Number density of soot (# / cm3)
      !         Asoot   = Surface area density of soot (cm2 / cm3)
      !         Vsoot   = Volume of soot (cm3)
      !
      ! OUTPUT:
      !         DeltaQ  = Surface rate with inter-sectional distribution (cm3/cm3)
      !         ProdR   = Production rate of reactant species (mol/cm3 s)
      !         ProdP   = Production rate of product species (mol/cm3 s)
      !
      !-----------------------------------------------------------------------

      implicit none

      integer(ip), intent(in)        :: pgaus
      integer(ip), intent(in)        :: ST_Type
      character(len = 2), intent(in) :: RName
      real(rp),    intent(in)        :: VolAdd
      real(rp),    intent(in)        :: Temp(pgaus),Rcoef(pgaus)
      real(rp),    intent(out)       :: DeltaQ(pgaus,nsect_ssm,3)
      real(rp),    intent(in)        :: stoiMR,stoiMP
      real(rp),    intent(in)        :: NDsoot(pgaus,nsect_ssm), &
                                        Asoot(pgaus,nsect_ssm),  &
                                        Vsoot(pgaus,nsect_ssm)
      real(rp),    intent(inout)     :: ProdR(pgaus),ProdP(pgaus)

      integer(ip)                    :: idsm, k, iStericFac
      real(rp)                       :: Xsoot, VOH, Vol_i
      real(rp), dimension(pgaus)     :: TotalRate,par_a,par_b,nCarbon, &
                                        StericFact,Csoot,RateQ,AS,ModRate, &
                                        Betaij_OH

      iStericFac = 2       ! 1 = unity ; 2 = ABF2000

      !
      ! Nominal number of sites per cm^2
      !
      Xsoot = 2.3e+15_rp

      !
      ! Initilizing the rate summation for species production
      !
      TotalRate = 0.0_rp
      DeltaQ    = 0.0_rp

      do idsm = 1,nsect_ssm
          !
          ! Reaction name
          !
          if (RName == 'R6') then
             VOH = 1.2_rp/NAvo ! Molecular volume of OH
             Vol_i = Vpmean(idsm)
             call soot_collision_freq(pgaus,Vol_i,VOH,Temp,Betaij_OH)
             RateQ = Betaij_OH*NDsoot(1:pgaus,idsm) / NAvo

          !
          ! Steric factor  definition  (fraction of actie sites that will react)
          !
          else
             if (iStericFac == 1) then
                StericFact = 1.0_rp
             else if (iStericFac == 2) then
                par_a = 12.65_rp - 5.63e-03_rp * Temp
                par_b = -1.38_rp + 6.80e-04_rp * Temp
                nCarbon = max(Vsoot(1:pgaus,idsm),1.0e-40_rp)/VC          ! Average number of C atoms in a particles of section isec
                StericFact =  tanh((par_a/ log10(nCarbon)) + par_b)
                do k = 1,pgaus
                   if (StericFact(k) < 0.0_rp) then
                       StericFact(k) = 0.0_rp
                   end if
                end do
             end if
             AS = Asoot(1:pgaus,idsm)                     ! total soot surface area in section i
             Csoot = Xsoot * AS /NAvo                     ! concentration of sites available
             RateQ = Csoot*StericFact                     ! concentration of reactive sites
          end if

          TotalRate = TotalRate + RateQ                   ! summation of [S] over sections
          ModRate   = max(NAvo*(Rcoef*RateQ), 1.0e-32_rp)
          DeltaQ(1:pgaus,idsm,1) = VolAdd*ModRate         ! soot surface rate in section i
      end do

      !
      ! Inter-sectional distribution
      !
      if (ST_Type > 0) &
         call soot_intersectdynamics(pgaus,ST_Type,DeltaQ)

      !
      ! Total specie rate over all sections for Reation
      !
      ProdR = ProdR + stoiMR*(Rcoef*TotalRate)
      ProdP = ProdP + stoiMP*(Rcoef*TotalRate)

    end subroutine soot_surf_reac_rates

    subroutine soot_write_output(OutFname,OutVar,nrows,ncols)
      !-----------------------------------------------------------------------
      !****f* Chemic/
      ! NAME
      !
      !
      ! DESCRIPTION
      !
      ! USES
      ! USED BY
      !
      !***
      !-----------------------------------------------------------------------
      implicit none

      integer(ip), intent(in)  :: nrows,ncols
      character(9), intent(in) :: OutFname
      real(rp), intent(in)     :: OutVar(nrows,ncols)
      integer(ip)              :: i,k

      ! Writing nucleation source and contribution to gas phase

      open(1, file = OutFname, status='replace')

      do k=1,nrows
         do i=1,ncols
            write(1,*) OutVar(k,i)
         end do
      end do

      close(1)

      print*, '    Outputfile ', OutFname,' written!'

    end subroutine soot_write_output

    subroutine header()

        implicit none

        print*,' '
        print*,'---------------------------------------------------------------------'
        print*,'---------------------------------------------------------------------'
        print*,' '
        print*,'    SOOT SOURCE TERMS COMPUTATION MODULE '
        print*,' '
        print*,'---------------------------------------------------------------------'
        print*,'---------------------------------------------------------------------'
        print*,' '
        print*,'    Starting computations..'
        print*,' '

    end subroutine header

    subroutine footer()

        implicit none

        print*,' '
        print*,'    Computations done!'
        print*,' '
        print*,'---------------------------------------------------------------------'
        print*,'---------------------------------------------------------------------'
        print*,' '

    end subroutine footer

end module mod_chm_sectional_soot_model

