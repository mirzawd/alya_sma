!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Exmedi
!> @{
!> @file    exm_iniunk.f90
!> @author  Mariano Vazquez
!> @date    26/12/2016
!> @brief   Starting: This routine sets up the initial conditions
!> @details Starting: This routine sets up the initial conditions
!> @}
!-----------------------------------------------------------------------
subroutine exm_iniunk
  use def_master
  use def_domain
  use mod_iofile
  use def_exmedi
  use mod_maths,               only : maths_distribute_iterations
  use mod_communications,      only : PAR_BROADCAST
  use mod_messages,            only : messages_live
  use mod_memory,              only : memory_initia, memory_alloca, memory_deallo
  use mod_exm_memory,          only : memory_alloca, memory_deallo
  use mod_exm_ohararudy
  use mod_exm_onetor
  use mod_exm_inicourte
  use mod_outfor
  use mod_eccoupling
  use mod_exm_fitzhugh_nagumo, only : exm_fitzhugh_nagumo_Initialize, exm_fitzhugh_nagumo_WriteInitialValuesToLog, exm_fitzhugh_nagumo_GetCa0
  use mod_exm_activation,      only : exm_stim_enumerate_activated_nodes
  use mod_exm_drugs,           only : exm_drugs_write_log

  use mod_exm_cellmodel

  implicit none



  integer(ip)                          :: kmodel_imate, kmodel_ipoin, imate
  integer(ip)                          :: any_active_material_id, ipair_active
  integer(ip)                          :: ipoin, icelltype, iconc, i
  integer(ip)                          :: onetor_status_tmp
  integer(ip)                          :: ipair, npairs,root_rank
  character(30)                        :: mesau(2)
  integer(ip), dimension(:,:), pointer :: mater_cell_pairs

  integer(ip), dimension(:),   pointer       :: list_pairs    
  real(rp),    dimension(nsdvs_ecc,mnode)    :: el_land
  real(rp),    dimension(nsdvs_ecc,mgaus)    :: gp_land
  TYPE(CELL_MODEL_OUTPUTS)                        :: cellm_out
  TYPE(CELL_MODEL_OUTPUTS), dimension(:), pointer :: cellm_out_allmaterials

  integer(ip)                          :: igaus, inode, pnode, pgaus, ielem, pelty

  character(15)                        :: steady_variablename, cellmodel_variablename !used for output

  !---------------------------------------------------------------------
  !
  !  Used only for TOR model until someone refactors is, to make old TOR code compatible with the new cell model usage
  !
  real(rp),    dimension(nsdvs_ecc)    :: tor_variables ! used temporarily
  integer(ip), parameter               :: torord_nstats = 3_ip
  real(rp)                             :: torord_stats(torord_nstats)  ! saves ohara stats: number of beats, tolerance, rmse
  !
  !----------------------------------------------------------------------------


  nullify(mater_cell_pairs,list_pairs)
  nullify(cellm_out_allmaterials)

  el_land = 0.0_rp

  !Find all the points that will be activated by the stimuli from the table
  !We are not saving these in the restart
  call exm_stim_enumerate_activated_nodes()


  if( kfl_rstar == 0_ip ) then
    !
    ! Initialize
    !
    call memory_initia(elmag)
    call memory_initia(vconc)
    !
    ! check if any TT-like model is present
    !
    vconc_initial    = 0.0_rp
    vauxi_exm_initial    = 0.0_rp
    kmodel_imate = 0_ip
    kmodel_ipoin = 0_ip

    do imate= 1,nmate
      !
      ! message: compute initial conditions for some cell models
      !
      kmodel_imate = kfl_cellmod(imate)

      if ( any( kmodel_imate == (/EXMSLD_CELL_TT2006, EXMSLD_CELL_OHARA, EXMSLD_CELL_OHARA_INAPA, EXMSLD_CELL_TORORD, EXMSLD_CELL_COURTE/) ) ) then
        call messages_live('  COMPUTING INITIAL PHYSIOLOGICAL CONDITIONS...')
        !if (kmodel_imate == 4 .or. kmodel_imate==5 .or. kmodel_imate==6 ) then
        call messages_live('    ITERATING ONE CELL MODEL...')
        mesau(1)= intost(moneclmate_exm(1,imate)) ! beats
        mesau(2)= intost(moneclmate_exm(2,imate)) ! cyclelength
        call messages_live('    BEATS=       '//trim(mesau(1)))
        call messages_live('    CYCLELENGTH= '//trim(mesau(2)))
      end if

      call messages_live('  INITIAL CONDITIONS FOR MATERIAL= '//trim(intost(imate)))
      select case(kmodel_imate)
      case (EXMSLD_CELL_FITZHUGH) 
        call messages_live('    SUB CELL MODEL: FITZHUGH NAGUMO')
      !case (EXMSLD_CELL_TT2006) 
      !   call messages_live('    SUB CELL MODEL: TT HETEROGENEOUS')
      case (EXMSLD_CELL_OHARA) 
        call messages_live('    SUB CELL MODEL: OHARA - RUDY')
      case (EXMSLD_CELL_OHARA_INAPA)
        call messages_live('    SUB CELL MODEL: OHARA - RUDY, INa Passini')
      case (EXMSLD_CELL_TORORD) 
        call messages_live('    SUB CELL MODEL: TOR-ORD')
      case (EXMSLD_CELL_COURTE) 
        call messages_live('    SUB CELL MODEL: COURTEMANCHE')
      end select

      if ( any( kmodel_imate == (/EXMSLD_CELL_OHARA, EXMSLD_CELL_OHARA_INAPA, EXMSLD_CELL_TORORD/) ) ) then
        select case ( kfl_steadystate_variable(imate) )
          case ( EXM_CELL_STEADY_VOLTAGE )
            call messages_live("    USING VOLTAGE TO TEST CELL CONVERGENCE FOR MATERIAL.")
          case ( EXM_CELL_STEADY_CALCIUM )
            call messages_live("    USING CALCIUM TO TEST CELL CONVERGENCE FOR MATERIAL.")
          case default
            call runend("EXM_INIUNK: Unknown name of the variable to determine the cell convergence.")
        end select
      end if

    end do


    npairs = nmate*ncelltypes_ecc

    call memory_alloca(mem_modul(1:2,modul),'MATER_CELL_PAIRS','exm_iniunk', mater_cell_pairs, npairs, 2_ip)
    !call memory_alloca(mem_modul(1:2,modul),'COURTE_STATS'    ,'exm_iniunk', courte_stats, npairs, 3_ip)  !EVA: REVISAR ESTO

    !save material and cell type to a linear array for paralelization
    ipair = 1_ip
    do imate= 1,nmate
      do icelltype = 1,ncelltypes_ecc
        mater_cell_pairs(ipair, 1_ip) = imate
        mater_cell_pairs(ipair, 2_ip) = icelltype
        ipair = ipair + 1_ip
      end do
    end do
    !
    ! LPAIRS(IPAIR) is the list of rank owner of pair IPAIR
    !
    call memory_alloca(mem_modul(1:2,modul),'LIST_PAIRS',                  'exm_iniunk', list_pairs, npairs)
    call memory_alloca(mem_modul(1:2,modul),'LAND_VARIABLES_ALLMATERIALS', 'exm_iniunk', cellm_out_allmaterials, npairs)


    call maths_distribute_iterations(npart,npairs,list_pairs,DEFAULT_VALUE=kfl_paral)

    ! Calculate initial conditions for the cell model
    do ipair= 1,npairs
      
      cellm_out_allmaterials(ipair) % S =         0.0_rp
      cellm_out_allmaterials(ipair) % W =         0.0_rp
      cellm_out_allmaterials(ipair) % CaTRPN =    0.0_rp
      cellm_out_allmaterials(ipair) % B =         0.0_rp
      cellm_out_allmaterials(ipair) % zeta_s =    0.0_rp
      cellm_out_allmaterials(ipair) % zeta_w =    0.0_rp
      cellm_out_allmaterials(ipair) % Ca50 =      0.805_rp/sms_conversion_currents_exm
      cellm_out_allmaterials(ipair) % Lambda =    1.0_rp
      cellm_out_allmaterials(ipair) % nbeats =    -1_ip
      cellm_out_allmaterials(ipair) % toler =     -1.0_rp
      cellm_out_allmaterials(ipair) % rmse =      -1.0_rp
      cellm_out_allmaterials(ipair) % success =   -1_ip
      cellm_out_allmaterials(ipair) % dt      =   -1.0_rp



      if( list_pairs(ipair) == kfl_paral ) then

        imate = mater_cell_pairs(ipair, 1_ip)
        icelltype = mater_cell_pairs(ipair, 2_ip)

        kmodel_imate = kfl_cellmod(imate)

        if ( any( kmodel_imate == (/EXMSLD_CELL_OHARA, EXMSLD_CELL_OHARA_INAPA/) ) ) then

          !
          ! Initialization of the cell model: done only by the master and then distributed to the slaves
          !
          !modifies Global vars: vminimate_exm(icelltype,imate), vauxi_exm_initial(:,icelltype,imate), vconc_initial(:,icelltype,imate)
          !Modifies Local: elmlo_exm, viclo_exm, vcolo_exm, vaulo_exm
          !Modifies: cellm_out only is land model is used and there is coupling with solidz
          if( icelltype <= ncelltypes_per_model(kmodel_imate) ) then 
              call exm_oneohr( imate, icelltype, cellm_out )
              cellm_out_allmaterials(ipair) = cellm_out
          end if

        !else if (kmodel_imate == EXMSLD_CELL_TT2006) then
        !   !this TT model needs deep revision, it may not be working correctly
        !   !maybe will be removed in the future
        !   !modifies: vminimate_exm(icelltype,imate), vauxi_exm_initial(:,icelltype,imate), vconc_initial(:,icelltype,imate)
        !   call exm_onecel(imate) ! these are precalculated, just sets up the variables runs extremeley fast. Needs to be refactored to pass also celltype


        else if (kmodel_imate == EXMSLD_CELL_TORORD) then

          !
          ! Initialization of the cell model: done only by the master then distributed to the slaves
          !
          if( icelltype <= ncelltypes_per_model(kmodel_imate) ) then 
              torord_stats = -1.0_rp
              onetor_status_tmp = -1_ip
              call exm_onetor(imate, icelltype, tor_variables, onetor_status_tmp, torord_stats )

              cellm_out % S =       tor_variables(1)
              cellm_out % W =       tor_variables(2)
              cellm_out % CaTRPN =  tor_variables(3)
              cellm_out % B =       tor_variables(4)
              cellm_out % zeta_s =  tor_variables(5)
              cellm_out % zeta_w =  tor_variables(6)
              cellm_out % Lambda =  tor_variables(7)
              cellm_out % Ca50   =  0.805_rp / sms_conversion_currents_exm
              cellm_out % nbeats =  torord_stats(1_ip)
              cellm_out % toler =   torord_stats(2_ip)
              cellm_out % rmse =    torord_stats(3_ip)
              cellm_out % success = onetor_status_tmp
              cellm_out % dt      = -1.0_rp

              cellm_out_allmaterials(ipair) = cellm_out
          end if


        else if (kmodel_imate == EXMSLD_CELL_COURTE) then
          !
          ! Initialization of the cell model: done only by the master and then distributed to the slaves
          !
          !modifies Global vars: vminimate_exm(icelltype,imate), vauxi_exm_initial(:,icelltype,imate), vconc_initial(:,icelltype,imate)
          !Modifies Local: elmlo_exm, viclo_exm, vcolo_exm, vaulo_exm
          !Modifies: cellm_out only is land model is used and there is coupling with solidz

          if( icelltype <= ncelltypes_per_model(kmodel_imate) ) then 
              call exm_onecourte(imate, icelltype, cellm_out )
              cellm_out_allmaterials(ipair) = cellm_out
          end if

        end if

      end if
    end do !ipair

    !Send oneohr_status from children to master

    !
    ! Owner of the pair broadcast the result to others
    !
    do ipair = 1,npairs
      root_rank = list_pairs(ipair)
      call PAR_BROADCAST(mater_cell_pairs(ipair,1), ROOT_RANK=root_rank)
      call PAR_BROADCAST(mater_cell_pairs(ipair,2), ROOT_RANK=root_rank)

      !do i = 1, courte_nstats
      !   call PAR_BROADCAST(courte_stats(ipair,i), ROOT_RANK=root_rank)
      !end do

      call PAR_BROADCAST(cellm_out_allmaterials(ipair) % S,       ROOT_RANK=root_rank)
      call PAR_BROADCAST(cellm_out_allmaterials(ipair) % W,       ROOT_RANK=root_rank)
      call PAR_BROADCAST(cellm_out_allmaterials(ipair) % CaTRPN,  ROOT_RANK=root_rank)
      call PAR_BROADCAST(cellm_out_allmaterials(ipair) % B,       ROOT_RANK=root_rank)
      call PAR_BROADCAST(cellm_out_allmaterials(ipair) % zeta_s,  ROOT_RANK=root_rank)
      call PAR_BROADCAST(cellm_out_allmaterials(ipair) % zeta_w,  ROOT_RANK=root_rank)
      call PAR_BROADCAST(cellm_out_allmaterials(ipair) % Ca50,    ROOT_RANK=root_rank)
      call PAR_BROADCAST(cellm_out_allmaterials(ipair) % Lambda,  ROOT_RANK=root_rank)
      call PAR_BROADCAST(cellm_out_allmaterials(ipair) % nbeats,  ROOT_RANK=root_rank)
      call PAR_BROADCAST(cellm_out_allmaterials(ipair) % toler,   ROOT_RANK=root_rank)
      call PAR_BROADCAST(cellm_out_allmaterials(ipair) % rmse,    ROOT_RANK=root_rank)
      call PAR_BROADCAST(cellm_out_allmaterials(ipair) % success, ROOT_RANK=root_rank)
      call PAR_BROADCAST(cellm_out_allmaterials(ipair) % dt,      ROOT_RANK=root_rank)
  
    end do

    do imate = 1,nmate
      do icelltype = 1,ncelltypes_ecc
        ipair     = (imate-1)*ncelltypes_ecc + icelltype
        root_rank = list_pairs(ipair)

        call PAR_BROADCAST(vminimate_exm(icelltype,imate), ROOT_RANK=root_rank)
        do iconc = 1,nconc_exm
          call PAR_BROADCAST(vconc_initial(iconc, icelltype, imate), ROOT_RANK=root_rank)
        end do
        do iconc = 1,nauxi_exm
          call PAR_BROADCAST(vauxi_exm_initial(iconc, icelltype, imate), ROOT_RANK=root_rank)
        end do

      end do
    end do


    !----------------------------------------
    !
    ! For inactive materials, set the values from one of the active materials that calculate initial voltage and calcium
    ! to reduce the gradients of ca0 and elmag across materials
    !
    any_active_material_id = -1_ip ! to find the first active material
    do imate = 1_ip,nmate
      if( all(kfl_cellmod(imate) .ne. (/EXMSLD_CELL_NOMODEL, EXMSLD_CELL_FITZHUGH/)) ) then
        any_active_material_id = imate
        exit
      end if
    end do

    if (any_active_material_id<0_ip) then
      call messages_live("EXMEDI: No acitve material is defined (no cell models)","WARNING")
    else
      ipair_active     = (any_active_material_id-1)*ncelltypes_ecc + 1_ip 

      do imate = 1,nmate
        if( any(kfl_cellmod(imate) == (/EXMSLD_CELL_NOMODEL, EXMSLD_CELL_FITZHUGH/)) ) then
          vconc_initial(1,1,imate) = vconc_initial(1, 1, any_active_material_id)
          vminimate_exm(1,imate) = vminimate_exm(1,any_active_material_id)
          exit
        end if
      end do
    end if
    !
    !
    !----------------------------------------


    if (INOTSLAVE) then
      call print_stats
    end if ! INOTSLAVE

    !--------------------------------------------------------------
    !
    ! If one of the oharas were executed
    ! Extract Ca0 and Ca50 to be saved in restart if needed
    ! This allows running only exmedi until electrphysiology converged, 
    ! save the restart and then run coupled problem
    !
    !--------------------------------------------------------------     
    if( kfl_exmsld_ecc  ) then
      do imate = 1,nmate
        kmodel_imate = kfl_cellmod(imate)

        select case( kmodel_imate )
        case (EXMSLD_CELL_TT2006, EXMSLD_CELL_OHARA, EXMSLD_CELL_OHARA_INAPA, EXMSLD_CELL_COURTE, EXMSLD_CELL_TORORD) 
          ! Save Ca0 values to automatically calibrate Hunter-McCulloch ECC model in sld_eccoup

          do icelltype = 1,ncelltypes_ecc
            call eccou_set_ca0( icelltype, imate, vconc_initial(1,icelltype,imate) * 1000.0_rp )

            ipair     = (imate-1)*ncelltypes_ecc + icelltype 
            call eccou_set_ca50( icelltype, imate, cellm_out_allmaterials(ipair) % Ca50 * 1000.0_rp )

            !sanity check. Ca0 is normalized in mod_exm_sld_eccoupling
            if ( vconc_initial(1,icelltype,imate) >= cellm_out_allmaterials(ipair) % Ca50 ) then
              if ( kmodel_imate == EXMSLD_CELL_TORORD ) then 
                call messages_live("EXM_INIUNK: Ca0 >= Ca50 for material "//trim(intost(imate))//" celltype "//trim(intost(icelltype)) ,"WARNING")
              else 
                !print *, "ca0,ca50:", vconc_initial(1,icelltype,imate), cellm_out_allmaterials(ipair) % Ca50
                call runend("EXM_INIUNK: Ca0 >= Ca50 for material "//trim(intost(imate))//" celltype "//trim(intost(icelltype)))
              end if
            end if
          end do
        case( EXMSLD_CELL_FITZHUGH ) 
          ! Save Ca0 values to automatically calibrate Hunter-McCulloch ECC model in sld_eccoup
          do icelltype = 1,ncelltypes_ecc
            call eccou_set_ca0 ( icelltype, imate, exm_fitzhugh_nagumo_GetCa0(imate) * 1000.0_rp )
            call eccou_set_ca50 ( icelltype, imate, 1.0_rp ) !the model does not need it, but to avoid a random value
          end do
        case (EXMSLD_CELL_NOMODEL) 
          do icelltype = 1,ncelltypes_ecc
            call eccou_set_ca0 ( icelltype, imate, vconc_initial(1,1,imate) * 1000.0_rp )
            call eccou_set_ca50 ( icelltype, imate, 1.0_rp ) !the model does not need it, but to avoid a random value
          end do
        case default
          call runend("exm_iniunk: unknown cell model in material "//trim(intost(imate)))
        end select
      end do
    end if

    !--------------------------------------------------------------
    !
    ! Distribute exmedi-solidz coupling parameters among cells and points
    !
    !--------------------------------------------------------------     
    if( kfl_exmsld_ecc ) then
      if( kfl_exmsld_lands_ecc )then

        call messages_live("EXM_INIUNK: PASSING LAND PARAMETERS TO SOLIDZ") 

        !pass calcium

        do ipoin = 1,npoin
          imate = nodemat(ipoin)

          if(kfl_celltype_fun .ne. 0_ip) then !if there is celltype field
            icelltype = int(celty_exm(ipoin), kind=ip)            
          else
            icelltype = 1_ip
          end if

          ipair     = (imate-1)*ncelltypes_ecc + icelltype            

          cellm_out = cellm_out_allmaterials(ipair)

          troponin_ecc(ipoin,1) = cellm_out % CaTRPN
          troponin_ecc(ipoin,2) = cellm_out % CaTRPN
        end do

        !interpolate land values for each gauss point from points
        do ielem = 1,nelem
          !
          ! Get element shape function
          !
          pelty = ltype(ielem)
          if( pelty > 0 ) then
            pgaus = ngaus(pelty)
            pnode = lnnod(ielem)
            
            do inode = 1,pnode
              ipoin          = lnods(inode,ielem)
              imate          = nodemat(ipoin)
              
              if(kfl_celltype_fun .ne. 0_ip) then !if there is celltype field
                icelltype      = int(celty_exm(ipoin), kind=ip)
              else 
                icelltype = 1_ip
              end if

              ipair             = (imate-1)*ncelltypes_ecc + icelltype            
              el_land(1, inode) = cellm_out_allmaterials(ipair) % S
              el_land(2, inode) = cellm_out_allmaterials(ipair) % W
              el_land(3, inode) = cellm_out_allmaterials(ipair) % CaTRPN
              el_land(4, inode) = cellm_out_allmaterials(ipair) % B
              el_land(5, inode) = cellm_out_allmaterials(ipair) % zeta_s
              el_land(6, inode) = cellm_out_allmaterials(ipair) % zeta_w
              el_land(7, inode) = cellm_out_allmaterials(ipair) % Lambda            
            end do

            gp_land = 0.0_rp
            do igaus = 1,pgaus
              do inode = 1,pnode
                do i = 1,nsdvs_ecc
                  gp_land(i,igaus) = gp_land(i,igaus) + el_land(i,inode) * elmar(pelty)%shape(inode,igaus)
                end do
              end do !inode
            end do !igauss

            call eccou_set_sdv_at_gp( ielem, gp_land )

          end if !pelty > 0 
        end do !ielem
      endif
    end if ! distribute parameters

    call memory_deallo(mem_modul(1:2,modul),'MATER_CELL_PAIRS'              ,'exm_iniunk', mater_cell_pairs)
    call memory_deallo(mem_modul(1:2,modul),'LAND_VARIABLES_ALLMATERIALS'   ,'exm_iniunk', cellm_out_allmaterials)

    call messages_live('  SETTING INITIAL CONDITIONS, NO RESTART...')
    !
    ! ELMAG
    !
    !TODO: check if this loop can be eliminated, seems to be repeated in the select below
    do ipoin = 1,npoin
      if (kfl_voini_exm(nodemat(ipoin)) == 1_ip) then
        elmag(ipoin,1) = vminimate_exm(1_ip, nodemat(ipoin))     ! Hoy restart
      end if
    end do
    !
    ! FISOC
    !
    do ipoin = 1,npoin
      fisoc(:, ipoin)   = -1.0_rp
      elmag_minmax(1, ipoin) = huge(0.1_rp)
      elmag_minmax(2, ipoin) = -huge(0.1_rp)
    end do
    !
    !  VAUXI_EXM and VCONC
    !
    do ipoin = 1,npoin
      imate = nodemat(ipoin)
      kmodel_ipoin=kfl_cellmod(imate)
      select case(kmodel_ipoin)
        !case (EXMSLD_CELL_NOMODEL) 
        ! nothing
        case (EXMSLD_CELL_FITZHUGH) 
          call exm_fitzhugh_nagumo_Initialize(ipoin, imate)
        !case (EXMSLD_CELL_TT2006)
        !   call exm_inihet(ipoin,imate)
        !   call runend('EXM_INIUNK: JAZMIN, REVISAR ESTO DEL INIHET!!!')
        case (EXMSLD_CELL_OHARA, EXMSLD_CELL_OHARA_INAPA)       !O'hara-Rudy
          call exm_iniohr(ipoin,imate)
        case (EXMSLD_CELL_TORORD)       ! ToR-ORd
          call exm_initor(ipoin,imate)
        case (EXMSLD_CELL_SCATRIA)      !Stem cell atria
          call exm_inisca(ipoin)
        case (EXMSLD_CELL_SCVENTRI)     !Stem cell ventricle
          call exm_iniscv(ipoin)
        case (EXMSLD_CELL_COURTE)   !Courtemanche atria
          call exm_inicourte(ipoin,imate)
        case (EXMSLD_CELL_NOMODEL)       !O'hara-Rudy
          vconc(1,ipoin,1) = vconc_initial(1,1,imate)
          elmag(ipoin,1:3) = vminimate_exm(1,imate)
      end select

      !Initialize min/max to the initial voltage
      elmag_minmax(:,ipoin) = elmag(ipoin,1)
    end do
    !
    ! Put all components of ELMAG to initial values
    !
    call exm_updunk(ITASK_INIUNK)

    ! Save coupling variables to exm.log
    if( kfl_exmsld_ecc ) then
      if (INOTMASTER) then
        call eccou_set_ca( vconc(1,:,1), 1000.0_rp ) ! pass calcium to eccoupling
      end if
      call eccou_finalize_initial_conditions()  
      call eccou_save_log()
    end if
   
  end if


contains

  subroutine print_stats
    implicit none

    character(len=200) stats_ohr


    !Master should print the information about the cell model convergence
    call outfor(25_ip, momod(modul) % lun_outpu,'OHARA CONVERGENCE')


    do ipair= 1,npairs
      imate = mater_cell_pairs(ipair, 1_ip)
      icelltype = mater_cell_pairs(ipair, 2_ip)
      kmodel_imate = kfl_cellmod(imate)

      if( any( kmodel_imate==(/EXMSLD_CELL_OHARA, EXMSLD_CELL_OHARA_INAPA, EXMSLD_CELL_TORORD, EXMSLD_CELL_COURTE/) ) ) then

        select case (kfl_steadystate_variable(imate))
        case (EXM_CELL_STEADY_VOLTAGE)
          steady_variablename = 'VOLTAGE'
        case (EXM_CELL_STEADY_CALCIUM)
          steady_variablename = 'CALCIUM'
        case default
          call runend("EXM_OCEOHR or EXM_OCETOR: Unknown name of the variable to determine the cell convergence.")
        end select


        select case ( kmodel_imate )
        case (EXMSLD_CELL_OHARA)
          cellmodel_variablename = 'OHARA'
        case (EXMSLD_CELL_OHARA_INAPA)
          cellmodel_variablename = 'OHARA-PASSINI'
        case (EXMSLD_CELL_TORORD)
          cellmodel_variablename = 'TOR'
        case (EXMSLD_CELL_COURTE)
          cellmodel_variablename = 'COURTE'
        case default
          call runend("EXM_INIUNK: Unknown cell type model.")
        end select


        if ( cellm_out_allmaterials(ipair) % success  >= 0_ip) then !if negative -- not executed
          stats_ohr = "MATERIAL "//trim(intost(imate))//&
                      ", CELLMODEL "//trim(cellmodel_variablename)//&
                      ", CELLTYPE "//trim(intost(icelltype))// &
                      ", BEATS=" //trim(intost(cellm_out_allmaterials(ipair) % nbeats))//&
                      ", TOL="//trim(retost(cellm_out_allmaterials(ipair) % toler))// &
                      ", RMSE="//trim(retost(cellm_out_allmaterials(ipair) % rmse ))

          if ( cellm_out_allmaterials(ipair) % success == EXM_CELLMODEL_NOTCONVERGED ) then
            call messages_live("CELL MODEL DID NOT CONVERGE. "//trim(stats_ohr) ,"WARNING")
            write(momod(modul)%lun_outpu,*) "CELL MODEL "//trim(steady_variablename)//" NOT CONVERGED. "//trim(stats_ohr)

            if( kfl_ignore_steadystate_celltypes_exm(imate, icelltype) == 0_ip ) then
              call runend("STEADY STATE NOT REACHED AFTER "//&
              trim(intost(cellm_out_allmaterials(ipair) % nbeats))//" BEATS. TRY INCREASING THE NUMBER OF BEATS. "//trim(stats_ohr))
            end if
          else if ( cellm_out_allmaterials(ipair) % success == EXM_CELLMODEL_CONVERGED ) then
            ! DO NOT ERASE THIS MESSAGE
            call messages_live("    CELL MODEL !!!DID!!! CONVERGE. "//trim(stats_ohr) )
            write(momod(modul)%lun_outpu,*) "CELL MODEL "//trim(steady_variablename)//" CONVERGED. "//trim(stats_ohr)
          else if ( cellm_out_allmaterials(ipair) % success == EXM_CELLMODEL_LOADED ) then
            ! DO NOT ERASE THIS MESSAGE
            call messages_live("    CELL MODEL LOADED. "//trim(stats_ohr) )
            write(momod(modul)%lun_outpu,*) "CELL MODEL "//trim(steady_variablename)//" LOADED. "//trim(stats_ohr)
          else if ( cellm_out_allmaterials(ipair) % success == EXM_CELLMODEL_NOTINITIALIZED ) then
            call messages_live("MATERIAL "//trim(intost(imate))//" CELLTYPE "//trim(intost(icelltype))//", VOLTAGES WERE NOT INITIALIZED. UNKNOWN MYOCYTE TYPE. MOST LIKELY REQUESTED NONEXISTING HARDCODED CELLTYPE.","WARNING")
            write(momod(modul)%lun_outpu,*) "OHARA NOT INITIALIZED. "//trim(stats_ohr) 
          else
            call runend("EXM_INIUNK: UNHANDLED RETURN VALUE OF EXM_ONEOHR - "//trim(intost(cellm_out_allmaterials(ipair) % success)))
          end if
        end if

                  
        !Save coupling related variables to the log
        write(momod(modul)%lun_outpu,          *) NEW_LINE('a')
        write(momod(modul)%lun_outpu,          *) " Coupling variables, MATERIAL "//trim(intost(imate))//" CELLTYPE "//trim(intost(icelltype))
        write(momod(modul)%lun_outpu,"(a,e16.8E3)") "   Ca0    = ", vconc_initial(1,icelltype,imate)
        write(momod(modul)%lun_outpu,"(a,e16.8E3)") "   Ca50   = ", cellm_out_allmaterials(ipair) % Ca50
        write(momod(modul)%lun_outpu,"(a,e16.8E3)") "   S      = ", cellm_out_allmaterials(ipair) % S
        write(momod(modul)%lun_outpu,"(a,e16.8E3)") "   W      = ", cellm_out_allmaterials(ipair) % W
        write(momod(modul)%lun_outpu,"(a,e16.8E3)") "   CaTRPN = ", cellm_out_allmaterials(ipair) % CaTRPN
        write(momod(modul)%lun_outpu,"(a,e16.8E3)") "   B      = ", cellm_out_allmaterials(ipair) % B
        write(momod(modul)%lun_outpu,"(a,e16.8E3)") "   zeta_s = ", cellm_out_allmaterials(ipair) % zeta_s
        write(momod(modul)%lun_outpu,"(a,e16.8E3)") "   zeta_w = ", cellm_out_allmaterials(ipair) % zeta_w
        write(momod(modul)%lun_outpu,"(a,e16.8E3)") "   Lambda = ", cellm_out_allmaterials(ipair) % Lambda
        write(momod(modul)%lun_outpu,"(a,e16.8E3)") "   dt(s)  = ", cellm_out_allmaterials(ipair) % dt
      end if

      write(momod(modul)%lun_outpu,*) NEW_LINE('a'), NEW_LINE('a')
    end do !ipair

    call exm_ohara_conductances_write_log()
    call exm_fitzhugh_nagumo_WriteInitialValuesToLog()
    call exm_drugs_write_log()

    flush(momod(modul)%lun_outpu) 
  end subroutine print_stats


end subroutine exm_iniunk
