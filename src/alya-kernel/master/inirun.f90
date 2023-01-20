!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine inirun()
  !-----------------------------------------------------------------------
  !****f* master/inirun
  ! NAME
  !    inirun
  ! DESCRIPTION
  !    This subroutine initializes and defines some run parameters
  ! USES
  ! USED BY
  !    Alya
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_domain
  use def_meshin
  use def_master
  use def_solver
  use def_postpr
  use def_inpout
  use def_coupli
  use mod_opebcs
  use def_coupli
  use mod_memory,                   only : lun_memor
  use mod_memory,                   only : lun_varcount
  use mod_memory,                   only : memory_initialization
  use mod_postpr,                   only : postpr_initialization
  use mod_elsest,                   only : elsest_initialization
  use mod_couplings,                only : couplings_initialization
  use mod_htable,                   only : htable_initialization
  use mod_alya_direct_solver,       only : alya_direct_solver_initialization
  use mod_messages,                 only : messages_initialization
  use mod_alya2metis,               only : alya2metis_initialization
  use mod_elmgeo,                   only : elmgeo_element_type_initialization
  use mod_timings,                  only : timings_initialization
  use mod_mpio_seq_postpr,          only : mod_mpio_seq_postpr_initialization
  use mod_element_data_base,        only : element_data_base_initialization
  use mod_communications,           only : communications_initialization
  use mod_redistribution,           only : redistribution_initialization
  use mod_repartitioning,           only : repartitioning_initialization
  use mod_partition_sfc,            only : partition_sfc_initialization
  use mod_std,                      only : std_initialization
  use mod_output_postprocess,       only : output_postprocess_initialization
  use mod_output_postprocess,       only : output_postprocess_allocate
  use mod_arrays,                   only : arrays_initialization
  use mod_renumbering_nodes,        only : renumbering_nodes_initialization
  use mod_mpio_par_readom,          only : mpio_initialization
  use mod_alya2talp,                only : alya2talp_initialization
  use mod_alya2gmsh,                only : alya2gmsh_initialization
  use mod_AMR,                      only : AMR_initialization
  use mod_error_estimate,           only : error_estimate_initialization
  use mod_arrays,                   only : arrays_register
  use mod_optimum_partition,        only : optimum_partition_initialization
  use mod_tubes,                    only : tubes_initialization
  use mod_reset,                    only : reset_initialization
  use mod_ecoute,                   only : ecoute_set_read_unit
  use mod_ecoute,                   only : ecoute_set_write_unit
  use mod_one_sided_communications, only : par_one_sided_initialization
  use mod_ecoute,                   only : ecoute_initialization
  implicit none
  integer(ip) :: imodu,jmodu,ii
  external    :: runend

  new_periodicity = 1
  
  call cputim(cpu_initi)         ! Initial CPU time
  !
  ! Initialization of some modules
  !
  call memory_initialization()
  call output_postprocess_initialization()
  call postpr_initialization()
  call elsest_initialization(ielse)
  call couplings_initialization()
  
  call alya_direct_solver_initialization()
  call htable_initialization()
  call htable_initialization(htable_lninv_loc)
  call messages_initialization()
  call alya2metis_initialization()
  call elmgeo_element_type_initialization()
  call timings_initialization()
  call mod_mpio_seq_postpr_initialization()
  call mpio_initialization()
  call renumbering_nodes_initialization()
  call element_data_base_initialization()
  call communications_initialization()
  call redistribution_initialization()
  call repartitioning_initialization()
  call partition_sfc_initialization()
  call arrays_initialization()
  call std_initialization()
  call alya2talp_initialization()
  call alya2gmsh_initialization()
  call AMR_initialization()
  call error_estimate_initialization()
  call optimum_partition_initialization()
  call tubes_initialization()
  call reset_initialization()
  call par_one_sided_initialization()
  call ecoute_initialization(runend)

  ! File units moved to Initia()

  endst         = 1              ! Stop Alya of end of file found
  kfl_split_plus= 0              ! By default the symbol + is not a separator
  call ecoute_set_read_unit (lun_pdata) ! Reading file
  call ecoute_set_write_unit(lun_outpu) ! Writing file
  !
  ! Memory
  !
  iar3p        = 0               ! ger3p not allocated
  iasca        = 0               ! gesca not allocated
  iavec        = 0               ! gevec not allocated
  iaten        = 0               ! geten not allocated
  !
  ! Master
  !
  kfl_check_data_file = 0        ! Check data file
  cpu_outpu     = 0.0_rp         ! Output CPU
  cpu_other     = 0.0_rp         ! CPU's
  cpu_start     = 0.0_rp         ! CPU's
  cpu_domain    = 0.0_rp         ! CPU's
  ittim         = 0              ! First time step
  itti2         = 0              ! First time step
  ittim_ini     = 0              ! Initial time step (=1 of no restart)
  ittim_rst     = 0              ! Last restart time step
  kfl_goopt     = 1              ! Go in optimization
  kfl_gotim     = 1              ! Go in time
  kfl_stop      = 0              ! if stop/=0 the current time-step will be the last
  kfl_gocou     = 1              ! Go in coupling iterations
  kfl_reset     = -1             ! Reset step -1 is off, 0 is on but not required, 1 is do reset
  memke         = 0              ! Current and maximum memory
  dtold         = 0.0_rp         ! Old time step
  dtime         = 0.0_rp         ! Time step
  ittyp         = 0              ! We are in initial run
  isect         = 0              ! Live output sections
  inews         = 0              ! Live output New section
  file_opened   = .true.         ! File was openend successfully
  nspec         = 0              ! Species
  kfl_naked     = 1              ! File names old option
  !
  ! Parall service
  !
  kfl_paral     = -1             ! Alya not initiated by MPI
  kfl_ptask     =  1             ! 0=Master ONLY partition
  kfl_postp_par =  1             ! Postprocess in Master postprocess file
  icoml         =  1             ! Current level number
  npart         =  1             ! Sequential run
  npasi         =  0_ip
  npari         =  0_ip
  nparx         =  0_ip
  npasr         =  0_ip
  nparr         =  0_ip
  nparl         =  0_ip
  call vocabu(-1_ip,0_ip,0_ip)
  !
  ! Solvers
  !
  nusol         = 1              ! Number of solves
  mxdof         = 0              ! D.o.f. per node
  memit         = 0              ! Iterative solver memory
  smemo         = 0              ! Solver memory
  kfl_symgr     = 0              ! Symmetric graph not needed (define in ***_inivar)
  kfl_schur     = 0              ! No Schur complement solver exists
  kfl_aiipr     = 0              ! No Aii preconditioner exists

  nzmat         = 1              ! # Max matrix size over all modules
  nzmbt         = 1              ! # Max RHS eigen matrix size over all modules
  nzrhs         = 1              ! # Max RHS size over all modules
  nzpre         = 1              ! # Max Preconditioner size over all modules

  nzmax         = 1              ! = max(nzsol,nzsky,nzexp): Components of A
  nzrhx         = 1              ! = max(nzsol,nzsky,nzexp): Components of RHS
  nzprx         = 1              ! Components of Preconditioner

  neige         = 1              ! # Max eigenvalue vector size
  neiva         = 1              ! # Max eignvalues
  nzerr         = 1              ! # Max Error Estimator size over all modules
  !
  ! Domain
  !
  kfl_conma     = 0              ! Consistent mass not nedded
  nzsky         = 1              ! # of comp. in the skyline matrix of the graph
  nzsol         = 1              ! # of comp. in the CSR matrix = nzdom
  mpopo         = 0              ! # Max Node-element connectivity
  iffun         = 0              ! Read bc codes (no function)
  ifloc         = 0              ! Read bc codes (no local bc)
  ifbop         = 0              ! Read bc codes (not on boundary nodes)
  ifbes         = 1              ! Read boundary values
  !
  ! Variables read in readat
  !
  kfl_elses  = 0                 ! Elsest
  xscal(1)   = 1.0_rp            ! X scale factor
  xscal(2)   = 1.0_rp            ! Y scale factor
  xscal(3)   = 1.0_rp            ! Z scale factor
  ielse(1)   = 100               ! nx
  ielse(2)   = 100               ! ny
  ielse(3)   = 100               ! nz
  ielse(4)   = 0                 ! data format (0=type,1=list)
  ielse(5)   = 10                ! Maximum number of possible meshes
  ielse(6)   = 2                 ! Second try strategy: if box is not found
  ielse(7)   = 0                 ! Output unit
  ielse(8)   = 1                 ! Search strategy (0=bin,1=Quad)
  ielse(9)   = 100               ! Points per node for Quad/Octtree
  ielse(10)  = 0                 ! Result output frequency
  ielse(11)  = 0                 ! Neighboring-boxes-search radius, 0: search until all boxes finished
  ielse(12)  = 0                 ! Postprocess mesh
  ielse(13)  = 0                 ! Postprocess results
  relse      = 0.0_rp            ! Elsest
  relse(1)   = 0.01_rp           ! Tolerance for iteration
  !
  ! Global variables read by modules
  !
  kfl_coibm     = 0              ! Immbou coupling: Read by IMMBOU.
  kfl_advec     = 0              ! Mesh advection: Read by IMMBOU.
  kfl_async     = 1              ! Parall: asynchronous communications
  kfl_forca_res = 0              ! Forces in residual formulation for rigid body(alefor decides)
  kfl_htran     = 0              ! 0: compute enthalpy diffusive fluxes in a simplified way as div(rho*lambda*grad(h)); 1: compute detailed enthalpy diffusive fluxes from Fourier's law and species fluxes. Use chemic module to define the not default value.
  thicl         = 0.0_rp         ! Interface thicknes. Read by LEVELS.
  mcono         = 3              ! Max # codes per nodes
  !
  ! Required arrays
  !
  kfl_lface     = 0              ! List of global faces not required: LFACG
  kfl_lelp2     = 0              ! List of extended node-element graph: PELPO_2, LELPO_2
  kfl_lelbf     = 0              ! List of element boundary faces: LELBF
  kfl_symgr     = 0              ! If symmetric graph is needed
  kfl_conma     = 0              ! If consistent mass is needed
  kfl_element_bin = 0            ! If element bin is required
  !
  ! What is being done
  !
  read_restart      = .false.
  write_restart     = .false.
  do_repartitioning = .false.
  do_AMR            = .false.
  !
  ! Variables which can be read or given as a command option
  !
  kfl_rstar     = 0              ! Not a restart run
  mitim         = 0              ! Number of time steps
  !
  ! Meshing
  !
  kfl_ifbox     = 0              ! If bounding box is prescribed
  !
  ! Modules
  !
  do modul = 0,mmodu
     momod(modul) % cpu_modul = 0.0_rp
     momod(modul) % mem_modul = 0
     momod(modul) % kfl_delay = 0
     momod(modul) % kfl_solve = AT_EACH_TIME_STEP
     nullify(momod(modul) % solve)
     nullify(momod(modul) % solad)
     nullify(momod(modul) % postp)
     nullify(momod(modul) % eigen)
     nullify(momod(modul) % times)
  end do
  modul            = 0                       ! I am kernel
  kfl_timei        = 0                       ! Steady calculation
  cpu_modul        = 0.0_rp                  ! Module CPU
  cpu_modul(CPU_MINI_ASSEMBLY,:) = 1e6_rp    ! Min CPU for elements, to remove variablity effects
  cpu_modul(CPU_MINI_NODE,:)     = 1e6_rp    ! Min CPU for node, to remove variablity effects
  mem_modul        = 0                       ! Module memory to zero
  namod            = '      '                ! Module names
  exmod            = '   '                   ! Module extension
  kfl_delay        = 0                       ! Delay module
  kfl_conve        = 1                       ! Module convergence required
  kfl_modul        = 0                       ! Module not used
  kfl_solve        = AT_EACH_TIME_STEP       ! When module is solved
  ndela            = 0                       ! Steps to delay module
  itinn            = 0                       ! First internal iteration
  kfl_modul(0)     = 1                       ! Kernel is always On!
  kfl_modul(mmodu) = 1                       ! Kermod is always On!
  kfl_modul(mmodu) = 1                       ! Kermod is always On!
  do imodu = 0,mmodu
     lzone(imodu) = 1
     do jmodu = 0,mmodu
        kfl_coupl(jmodu,imodu) = 0
        kfl_cowhe(jmodu,imodu) = 0
     end do
     do ii = 1,20
        kfl_itask(ii,imodu) = 0
     end do
  end do
  namod(ID_KERMOD) = 'KERMOD' ; exmod(ID_KERMOD) = 'ker'
  namod(ID_KERNEL) = 'KERNEL' ; exmod(ID_KERNEL) = 'ke2'
  namod(ID_NASTIN) = 'NASTIN' ; exmod(ID_NASTIN) = 'nsi'
  namod(ID_TEMPER) = 'TEMPER' ; exmod(ID_TEMPER) = 'tem'
  namod(ID_CODIRE) = 'CODIRE' ; exmod(ID_CODIRE) = 'cdr'
  namod(ID_TURBUL) = 'TURBUL' ; exmod(ID_TURBUL) = 'tur'
  namod(ID_EXMEDI) = 'EXMEDI' ; exmod(ID_EXMEDI) = 'exm'
  namod(ID_NASTAL) = 'NASTAL' ; exmod(ID_NASTAL) = 'nsa'
  namod(ID_ALEFOR) = 'ALEFOR' ; exmod(ID_ALEFOR) = 'ale'
  namod(ID_LATBOL) = 'LATBOL' ; exmod(ID_LATBOL) = 'lat'
  namod(ID_GUSANO) = 'GUSANO' ; exmod(ID_GUSANO) = 'gus'
  namod(ID_SOLIDZ) = 'SOLIDZ' ; exmod(ID_SOLIDZ) = 'sld'
  namod(ID_GOTITA) = 'GOTITA' ; exmod(ID_GOTITA) = 'got'
  namod(ID_WAVEQU) = 'WAVEQU' ; exmod(ID_WAVEQU) = 'wav'
  namod(ID_LEVELS) = 'LEVELS' ; exmod(ID_LEVELS) = 'lev'
  namod(ID_QUANTY) = 'QUANTY' ; exmod(ID_QUANTY) = 'qua'
  namod(ID_MAGNET) = 'MAGNET' ; exmod(ID_MAGNET) = 'mag'
  namod(ID_PARTIS) = 'PARTIS' ; exmod(ID_PARTIS) = 'pts'
  namod(ID_NASEDG) = 'NASEDG' ; exmod(ID_NASEDG) = 'nsg'
  namod(ID_CHEMIC) = 'CHEMIC' ; exmod(ID_CHEMIC) = 'chm'
  namod(ID_HELMOZ) = 'HELMOZ' ; exmod(ID_HELMOZ) = 'hlm'
  namod(ID_IMMBOU) = 'IMMBOU' ; exmod(ID_IMMBOU) = 'ibm'
  namod(ID_RADIAT) = 'RADIAT' ; exmod(ID_RADIAT) = 'rad'
  namod(ID_CASIMI) = 'CASIMI' ; exmod(ID_CASIMI) = 'cas'
  namod(ID_POROUS) = 'POROUS' ; exmod(ID_POROUS) = 'por'
  namod(ID_XXXXXX) = 'XXXXXX' ; exmod(ID_XXXXXX) = 'xxx'
  namod(ID_NEUTRO) = 'NEUTRO' ; exmod(ID_NEUTRO) = 'neu'
  namod(ID_INSITU) = 'INSITU' ; exmod(ID_INSITU) = 'ins'
  namod(ID_SOLFE2) = 'SOLFE2' ; exmod(ID_SOLFE2) = 'fe2'
  !
  ! Postprocess of kernel
  !
  call output_postprocess_allocate(ID_KERNEL)
  kfl_reawr             = 0         ! Postprocess mode
  ivapo                 = 0
  kfl_ivari             = 0
  modul                 = ID_KERNEL ! Now we deal with kernel's postprocess
  call moddef(9_ip)
  postp(1) % kfl_oonce =  1         ! Kernel postprocess is only once by defalut
  postp(1) % npp_iniso =  1         ! Kernel postprocess always initial values
  !
  ! Geometry
  !
  call arrays_register( 15_ip,(/'LNODS','VECTO','NELEM','SECON'/),lnods)
  call arrays_register( 16_ip,(/'COORD','VECTO','NPOIN','SECON'/),coord)
  call arrays_register( 17_ip,(/'LTYPE','SCALA','NELEM','SECON'/),ltype)
  call arrays_register( 18_ip,(/'LNINV','SCALA','NPOIN','SECON'/),lninv_loc)
  call arrays_register( 23_ip,(/'LEINV','SCALA','NELEM','SECON'/),leinv_loc)

  postp(1) % wopos( 1,19)  = 'LELCH'
  postp(1) % wopos( 1,20)  = 'LNODB'
  postp(1) % wopos( 1,21)  = 'LTYPB'
  postp(1) % wopos( 1,22)  = 'LESUB'
  postp(1) % wopos( 1,24)  = 'LMATE'
  postp(1) % wopos( 1,25)  = 'LBINV'
  postp(1) % wopos( 1,26)  = 'LBOCH'
  postp(1) % wopos( 1,27)  = 'LNOCH'
  postp(1) % wopos( 1,33)  = 'LELBO'
  postp(1) % wopos( 1,75)  = 'LMAST'
  !
  ! Sets
  !
  postp(1) % wopos( 1,30)  = 'LESET'
  postp(1) % wopos( 1,31)  = 'LBSET'
  postp(1) % wopos( 1,34)  = 'LNSET'
  !
  ! Boundary conditions
  !
  postp(1) % wopos( 1,28)  = 'CODBO'
  postp(1) % wopos( 1,29)  = 'CODNO'
  !
  ! Fields
  !
  postp(1) % wopos( 1,32)  = 'XFIEL'


  postp(1) % wopos( 2,15)  = 'VECTO'
  postp(1) % wopos( 2,16)  = 'VECTO'
  postp(1) % wopos( 2,17)  = 'SCALA'
  postp(1) % wopos( 2,18)  = 'SCALA'
  postp(1) % wopos( 2,19)  = 'SCALA'
  postp(1) % wopos( 2,20)  = 'VECTO'
  postp(1) % wopos( 2,21)  = 'SCALA'
  postp(1) % wopos( 2,22)  = 'SCALA'
  postp(1) % wopos( 2,23)  = 'SCALA'
  postp(1) % wopos( 2,24)  = 'SCALA'
  postp(1) % wopos( 2,25)  = 'SCALA'
  postp(1) % wopos( 2,26)  = 'SCALA'
  postp(1) % wopos( 2,27)  = 'SCALA'

  postp(1) % wopos( 2,28)  = 'SCALA'
  postp(1) % wopos( 2,29)  = 'VECTO'

  postp(1) % wopos( 2,30)  = 'SCALA'
  postp(1) % wopos( 2,31)  = 'SCALA'
  postp(1) % wopos( 2,34)  = 'SCALA'

  postp(1) % wopos( 2,32)  = 'VECTO'
  postp(1) % wopos( 2,33)  = 'SCALA'
  postp(1) % wopos( 2,75)  = 'SCALA'

  postp(1) % wopos( 3,15)  = 'NELEM'
  postp(1) % wopos( 3,16)  = 'NPOIN'
  postp(1) % wopos( 3,17)  = 'NELEM'
  postp(1) % wopos( 3,18)  = 'NPOIN'
  postp(1) % wopos( 3,19)  = 'NELEM'
  postp(1) % wopos( 3,20)  = 'NBOUN'
  postp(1) % wopos( 3,21)  = 'NBOUN'
  postp(1) % wopos( 3,22)  = 'NELEM'
  postp(1) % wopos( 3,23)  = 'NELEM'
  postp(1) % wopos( 3,24)  = 'NELEM'
  postp(1) % wopos( 3,25)  = 'NBOUN'
  postp(1) % wopos( 3,26)  = 'NBOUN'
  postp(1) % wopos( 3,27)  = 'NPOIN'

  postp(1) % wopos( 3,28)  = 'NBOUN'
  postp(1) % wopos( 3,29)  = 'NPOIN'

  postp(1) % wopos( 3,30)  = 'NELEM'
  postp(1) % wopos( 3,31)  = 'NBOUN'
  postp(1) % wopos( 3,34)  = 'NPOIN'

  postp(1) % wopos( 3,32)  = 'UNKNO'
  postp(1) % wopos( 3,33)  = 'NBOUN'
  postp(1) % wopos( 3,75)  = 'NPOIN'
  !
  ! nullify pointers
  !
  nullify(enthalpy_transport)
  nullify(chemical_heat)
  nullify(radiative_heat)
  nullify(div_enthalpy_transport)
  nullify(gefil)
  nullify(veloc)
  nullify(press)
  nullify(tempe)
  nullify(densi)
  nullify(energ)
  nullify(visco)
  nullify(umome)
  nullify(untur)
  nullify(uncdr)
  nullify(elmag)
  nullify(dispm)
  nullify(velom)
  nullify(displ)
  nullify(spins)
  nullify(disfu)
  nullify(vdrop)
  nullify(cdrop)
  nullify(wavam)
  nullify(fleve)
  nullify(erres)
  nullify(vorti)
  nullify(conce)
  nullify(cdens)
  nullify(entha)
  nullify(therm)
  nullify(phion)
  nullify(fiber)
  nullify(fisoc)
  nullify(elmag_minmax)
  nullify(gpfib)
  nullify(radso)
  nullify(vconc)
  nullify(taulo)
  nullify(kfl_fixno_ale)
  nullify(bvess_ale)
  nullify(gisca)
  nullify(givec)
  nullify(giscp)
  nullify(givep)
  nullify(geten)
  nullify(gevec)
  nullify(gesca)
  nullify(getep)
  nullify(gevep)
  nullify(gescp)
  nullify(getex)
  nullify(gevex)
  nullify(gescx)
  nullify(ger3p)
  nullify(turmu)
  nullify(rhoon)
  nullify(forcf)
  nullify(forca)
  nullify(wmean)
  nullify(visck)
  nullify(massk)
  nullify(lescl)
  nullify(momentum_sink)
  nullify(heat_sink)
  nullify(mass_sink)

  nullify(condk)
  nullify(sphek)
  nullify(advec)
  nullify(sphec)
  nullify(vesgs)
  nullify(tesgs)
  nullify(veset)
  nullify(vbset)
  nullify(vnset)
  nullify(witne)
  nullify(unkno)
  nullify(eigen)
  nullify(eigva)
  nullify(amatr)
  nullify(bmatr)
  nullify(rhsid)
  nullify(pmatr)
  nullify(pschu)
  nullify(aii)
  nullify(aib)
  nullify(abi)
  nullify(abb)
  nullify(xxi)
  nullify(xxb)
  nullify(bbi)
  nullify(bbb)
  nullify(lumma)
  nullify(amatx)
  nullify(rhsix)
  nullify(unknx)
  nullify(pmatx)
  nullify(parin)
  nullify(pari1)
  nullify(pari2)
  nullify(pari3)
  nullify(paris)
  nullify(parig)
  nullify(lnbin)
  nullify(lgpar)
  nullify(lnwit)
  nullify(ledgg)
  nullify(lfacg)
  nullify(ledgc)
  nullify(lfacb)
  nullify(parlo)
  nullify(parhh)
  nullify(nelem_tot)
  nullify(npoin_tot)
  nullify(nboun_tot)
  nullify(npoia_tot)
  nullify(npoin_par)
  nullify(nelem_par)
  nullify(nboun_par)
  nullify(lninv_loc)
  nullify(leinv_loc)
  nullify(lbinv_loc)
  nullify(lginv_loc)
  nullify(lpoi4)
  nullify(recvbuf_ip)
  nullify(sendbuf_ip)
  nullify(recvbuf_rp)
  nullify(sendbuf_rp)
  nullify(displs)
  nullify(recvcounts)

  nullify(iaren_par)
  nullify(jaren_par)
  nullify(permr_par)
  nullify(invpr_par)
  nullify(lelbf)
  nullify(lelfa)
  nullify(parre)
  nullify(parr1)
  nullify(parrs)
  nullify(parr2)
  nullify(parr3)
  nullify(parcx)
  nullify(parx1)
  nullify(parx2)
  nullify(parx3)
  nullify(par3p)
  nullify(pai1p)
  nullify(kfl_fixno_dod)
  nullify(lnsub_dod)
  nullify(lbsub_dod)
  nullify(imbou)
  nullify(rbbou)
  nullify(lnint)
  nullify(lntib)
  nullify(lnti2)
  nullify(lndom)
  nullify(lntra)
  nullify(letib)
  nullify(massc)
  nullify(statt)
  nullify(vel1d)
  nullify(areas)
  nullify(flowr)
  nullify(gdepo)
  nullify(gdeinv)
  nullify(gdedet)
  
  nullify(gescar_mpio)
  nullify(gevecr_mpio)
  nullify(gescai_mpio)
  nullify(geveci_mpio)

  nullify(I_AM_IN_ZONE)
  nullify(I_AM_IN_SUBD)
  nullify(I_AM_IN_CODE)
  nullify(application_names)
  !
  ! Code/zone/subdomain (essentially used by service COUPLI)
  !
  !current_zone = lzone(ID_KERMOD)
  current_zone = 0
  current_code = 1
  current_subd = 1
  kfl_gozon    = 1
  !
  ! Others
  !
  IEMPTY    = .false.
  INOTEMPTY = .true.

end subroutine inirun
