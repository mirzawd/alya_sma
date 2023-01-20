!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine ker_parall(order)
  !-----------------------------------------------------------------------
  !****f* Parall/ker_parall
  ! NAME
  !    ker_parall
  ! DESCRIPTION
  !    This routine exchange data
  ! USES
  ! USED BY
  !    Reapro
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_kermod
  use def_domain
  use def_inpout
  use def_solver
  use mod_opebcs
  use mod_elsest,                  only : search_elsest_seq
  use mod_elsest,                  only : search_elsest_par
  use mod_witness,                 only : witness_parall
  use mod_spare_mesh,              only : spare_mesh_parall
  use mod_memory,                  only : memory_alloca
  use mod_memory,                  only : memory_deallo
  use mod_ker_regularization,      only : kfl_regularization, kfl_second, reg_type
  use mod_ker_subdomain,           only : ker_subdomain_parall
  use mod_ker_tendencies,          only : kfl_tendencies_ker
  use mod_domain,                  only : domain_memory_allocate
  use mod_domain,                  only : domain_memory_deallocate
  use mod_local_basis,             only : local_basis_parall
  use mod_windk,                   only : mod_windk_systems_par_exchange
  use mod_ker_discrete_function,   only : ker_discrete_function_parall
  use mod_ker_tendencies,          only : ker_tendencies_parall
  use mod_ker_space_time_function, only : space_time_function_parall
  use mod_AMR,                     only : AMR_parall
  use mod_tubes,                   only : tubes_parall
  use mod_eccoupling,              only : eccou_send_data
  use mod_biofibers,               only : biofibers
  use mod_interp_tab,              only : typ_lookup_table
  use mod_interp_tab,              only : typ_tab_coord
  use mod_interp_tab,              only : tab_par_exchange
  use mod_interp_tab,              only : fw_par_exchange
  use mod_interp_tab,              only : fw_allocate
  use mod_ker_proper,              only : ker_proper_parall 
  use mod_reset,                   only : reset_parall
  use mod_mesh_type_basic,         only : mesh_type_basic_broadcast
  use mod_output_postprocess,      only : output_postprocess_parall_old
  use mod_exchange,                only : exchange_init
  use mod_exchange,                only : exchange_add
  use mod_exchange,                only : exchange_end
  use mod_solver,                  only : solver_parall
  use mod_output_postprocess,      only : output_postprocess_parall
  implicit none
  integer(ip), intent(in)   :: order
  integer(ip)               :: ii,ji,ki,ifunc
  integer(ip)               :: idime,kdime,nexpr
  type(typ_lookup_table),     pointer :: ptr_lookup_tab                    
  type(typ_tab_coord),        pointer :: ptr_lookup_coords(:)  

  if( ISEQUEN ) return

  select case ( order )

  case( 1_ip )

     !-------------------------------------------------------------------
     !
     ! Broadcast data read in *.ker.dat file
     !
     !-------------------------------------------------------------------

     call exchange_init()
     !
     ! Physical problem
     !
     call exchange_add(kfl_vefun)
     call exchange_add(kfl_tefun)
     call exchange_add(kfl_cofun)
     call exchange_add(kfl_difun)
     call exchange_add(kfl_arfun)
     call exchange_add(kfl_fiber_long_fun)
     call exchange_add(kfl_fiber_tang_fun)
     call exchange_add(kfl_fiber_norm_fun)
     call exchange_add(kfl_celltype_fun)
     call exchange_add(random_from_file)
     call exchange_add(kfl_rough)
     call exchange_add(kfl_canhe)
     call exchange_add(kfl_heiov)
     call exchange_add(kfl_canla)
     call exchange_add(kfl_delta)
     call exchange_add(kfl_ustar)
     call exchange_add(kfl_walld)
     call exchange_add(kfl_walln)
     do ji = 1,2
        call exchange_add(kfl_walld_field(ji))
     end do
     call exchange_add(kfl_suppo)
     call exchange_add(kfl_extro)
     call exchange_add(kfl_prope)
     call exchange_add(kfl_kxmod_ker)
     call exchange_add(kfl_noslw_ker)
     call exchange_add(kfl_waexl_ker)
     call exchange_add(kfl_waexl_imp_ker)
     call exchange_add(kfl_twola_ker)
     call exchange_add(kfl_mlwm_ker)
     do ji = 1,9
        call exchange_add(krestr_2_codno(ji))
     end do
     call exchange_add(nrestr_2_codno)
     call exchange_add(kfl_wlaav_ker)
     call exchange_add(kfl_aveme_ker)
     call exchange_add(kfl_prtur_abl_ker)
     call exchange_add(kfl_algebra_operations)
     do ji = 1,mcodb+1
        call exchange_add(kfl_boexc_ker(ji))
     end do
     do ji = 1, nmate
        call exchange_add(kfl_dampi_ker(ji))
     end do

     call exchange_add(denme)
     call exchange_add(visme)
     call exchange_add(gasco)
     call exchange_add(conce_relhu)
     call exchange_add(difun_facto)
     call exchange_add(u_ref)
     call exchange_add(h_ref)
     call exchange_add(k_ref)
     call exchange_add(usref)

     do ji = 1,3
        call exchange_add(windg(ji))
     end do

     call exchange_add(delta_dom)
     call exchange_add(delmu_dom) 
     call exchange_add(delta_sc)
     call exchange_add(delmu_ml)
     call exchange_add(rough_dom)
     call exchange_add(grnor)
     call exchange_add(gravi(1))
     call exchange_add(gravi(2))
     call exchange_add(gravi(3))
     call exchange_add(thicl)
     call exchange_add(cmu_st)
     call exchange_add(dexlo_ker)
     call exchange_add(tpeav_ker)
     call exchange_add(kfl_tendencies_ker)
     call exchange_add(tlape_ker)
     call exchange_add(cpres_ker(1))
     call exchange_add(rpres_ker(1))
     call exchange_add(cpres_ker(2))
     call exchange_add(rpres_ker(2))
     call exchange_add(num_lobas)
     !
     ! Coupling between modules
     !
     do ji = 0,mmodu
        do ki = 0,mmodu
           call exchange_add(kfl_coupl(ki,ji))
           call exchange_add(kfl_cowhe(ki,ji))
        end do
     end do
     !
     ! Properties and subdomains
     !
     call ker_proper_parall()
     call ker_subdomain_parall()
     !
     ! Numerical problem 
     !
     call exchange_add(kfl_renumbering_npoin)
     call exchange_add(kfl_renumbering_nelem)
     call exchange_add(nsfc_renumbering_npoin)
     call exchange_add(ndivi)
     call exchange_add(multiply_with_curvature)
     call exchange_add(kfl_elndi)
     call exchange_add(kfl_mmpar)
     call exchange_add(kfl_edge_elements)
     call exchange_add(kfl_rotation_axe)
     call exchange_add(kfl_graph)
     call exchange_add(kfl_elm_graph)
     call exchange_add(kfl_savda)
     call exchange_add(kfl_data_base_array)
     call exchange_add(kfl_vector_size)
     call exchange_add(kfl_lface)
     call exchange_add(kfl_fv_data)
     call exchange_add(kfl_grpro)
     call exchange_add(kfl_fixsm)
     call exchange_add(kfl_matrix_grad)
     call exchange_add(kfl_conbc_ker)
     call exchange_add(kfl_element_bin)
     call exchange_add(kfl_elses)
     call exchange_add(kfl_temper_vect)
     call exchange_add(kfl_chemic_vect)
     call exchange_add(kfl_soot_vect)
     call exchange_add(kfl_conta)
     call exchange_add(contact_tol)

     call exchange_add(kfl_reset) 

     do ji=1,nelse
        call exchange_add(ielse(ji))
     end do
     do ji=1,nelse
        call exchange_add(relse(ji))
     end do
     call exchange_add(rotation_angle)
     do ji = 1,3
        call exchange_add(rotation_axis(ji))
     end do

     call exchange_add(kfl_cutel)
     call exchange_add(kfl_hangi)
     call exchange_add(kfl_lapla)
     call exchange_add(kfl_defor)
     call exchange_add(kfL_coo)
     call exchange_add(kfL_ell)
     call exchange_add(kfl_full_rows)
     call exchange_add(kfl_element_to_csr)
     call exchange_add(kfl_direct_solver)
     call exchange_add(npoin_mm)
     call exchange_add(nboun_mm)
     call exchange_add(mnodb_mm)

     call exchange_add(number_space_time_function)                ! Space/Time functions
     do ifunc = 1,max_space_time_function
        call exchange_add(space_time_function(ifunc) % ndime)
        call exchange_add(space_time_function(ifunc) % nexpr)
        call exchange_add(space_time_function(ifunc) % name)
     end do

     call exchange_add(number_time_function)                      ! Time functions
     do ifunc = 1,max_time_function
        call exchange_add(time_function(ifunc) % npara)
        call exchange_add(time_function(ifunc) % kfl_type)
        call exchange_add(time_function(ifunc) % name)
     end do
     call exchange_add(number_windk_systems)                      ! Windkessel systems
     do ifunc = 1,max_windk_systems
        call exchange_add(windk_systems(ifunc) % sysid)
        call exchange_add(windk_systems(ifunc) % name)
        call exchange_add(windk_systems(ifunc) % wdks_model)
        call exchange_add(windk_systems(ifunc) % nparam)
        call exchange_add(windk_systems(ifunc) % ID_IN)
        call exchange_add(windk_systems(ifunc) % ID_OUT)
        call exchange_add(windk_systems(ifunc) % tag_in)
        call exchange_add(windk_systems(ifunc) % tag_out)
        call exchange_add(windk_systems(ifunc) % ndxs)
        call exchange_add(windk_systems(ifunc) % iflow_nsi)

        call exchange_add( windk_systems(ifunc) % discrete_function )

     end do

     call exchange_add(number_pump_curve)                      ! Pumps system
     do ifunc = 1,max_pump_curve
        call exchange_add( pump_curve(ifunc) % nparam)
        call exchange_add( pump_curve(ifunc) % model)
        call exchange_add( pump_curve(ifunc) % vhvad)
     end do


     call exchange_add(number_lookup_tab)                          ! Lookup tables 
     call exchange_add(number_lookup_fw)                      

     call exchange_add(deformation_steps)
     call exchange_add(deformation_strategy)
     call exchange_add(kfl_duatss)
     call exchange_add(fact_duatss)
     call exchange_add(kfl_conma)
     call exchange_add(kfl_conma_weighted)
     call exchange_add(kfl_approx_inv_mass)
     call exchange_add(kfl_dmass)
     call exchange_add(kfl_logva)
     call exchange_add(reg_type)
     call exchange_add(num_spare_meshes)
     !
     ! Output
     !
     call exchange_add(nwitn_all)
     call exchange_add(nwitn)
     call exchange_add(mwitn)
     call exchange_add(nwitg)
     call exchange_add(nwith)
     call exchange_add(kfl_posdo)
     call exchange_add(kfl_posdi)
     call exchange_add(kfl_oumes)
     call exchange_add(kfl_wimes)
     call exchange_add(kfl_oustl)
     call exchange_add(npp_stepo)
     call exchange_add(nfilt)
     call exchange_add(kfl_abovx)
     call exchange_add(resvx)
     call exchange_add(travx)
     call exchange_add(kfl_vortx)
     call exchange_add(kfl_vortx_thres)
     call exchange_add(kfl_detection)
     call exchange_add(kfl_pixel)
     call exchange_add(kfl_livee)
     call exchange_add(plane_pixel)
     call exchange_add(variable_pixel)
     call exchange_add(number_pixel(1))
     call exchange_add(number_pixel(2))
     call exchange_add(nsteps_ensemble)

     call exchange_add(kfl_node_graph)
     call exchange_add(kfl_edge_graph)

     do ii = 1,2
        do ji = 1,3
           call exchange_add(bobvx(ii,ji))
        end do
     end do
     do ii = 1,3
        call exchange_add(travx(ii))
     end do
     call exchange_add(thr_veloc)
     call exchange_add(thr_qvort)
     call exchange_add(detection_length)
     call exchange_add(detection_velocity)
     call exchange_add(coord_plane_pixel)
     !
     ! adjoint and optimization
     !
     call exchange_add(kfl_cost_type)
     call exchange_add(kfl_dvar_type)
     call exchange_add(kfl_adj_prob)
     call exchange_add(kfl_cos_opt)
     call exchange_add(kfl_ndvars_opt)
     call exchange_add(kfl_nwall)
     !
     ! Broadcasting of logical data
     !
     call exchange_add(kfl_regularization)
     call exchange_add(kfl_second)
     !
     ! Solver and postprocess
     !
     call solver_parall()
     call output_postprocess_parall()

     call exchange_end()

     !-------------------------------------------------------------------
     !
     ! Local basis coordinate systems
     !
     !-------------------------------------------------------------------

     call local_basis_parall()

     !-------------------------------------------------------------------
     !
     ! Reset
     !
     !-------------------------------------------------------------------

     call reset_parall()

     !-------------------------------------------------------------------
     !
     ! Witness points, geometries and meshes 
     !
     !-------------------------------------------------------------------

     call witness_parall()

     !-------------------------------------------------------------------
     !
     ! Spare meshes
     !
     !-------------------------------------------------------------------

     call spare_mesh_parall()

     !-------------------------------------------------------------------
     !
     ! AMR
     !
     !-------------------------------------------------------------------

     call AMR_parall()

     !-------------------------------------------------------------------
     !
     ! Electro-mechanical coupling
     !
     !-------------------------------------------------------------------
     
     call eccou_send_data(0_ip)

     !-------------------------------------------------------------------
     !
     ! Bio-fibers 
     !
     !-------------------------------------------------------------------

     call biofibers % send_data()

     !-------------------------------------------------------------------
     !
     ! Boundary conditions
     !
     !-------------------------------------------------------------------

     call boundary_conditions_exchange(tncod_ker)
     call boundary_conditions_exchange(tgcod_ker)
     call boundary_conditions_exchange(tbcod_ker)

     !-------------------------------------------------------------------
     !
     ! Discrete functions
     !
     !-------------------------------------------------------------------

     call ker_discrete_function_parall()

     !-------------------------------------------------------------------
     !
     ! Tendencies
     !
     !-------------------------------------------------------------------

     call ker_tendencies_parall()

     !-------------------------------------------------------------------
     !
     ! Wall exchange
     !
     !-------------------------------------------------------------------

     call search_waexlo_seq % parall()
     call search_waexlo_par % parall()

     !-------------------------------------------------------------------
     !
     ! Element search strategy
     !
     !-------------------------------------------------------------------

     call search_elsest_seq % parall()
     call search_elsest_par % parall()

     !-------------------------------------------------------------------
     !
     ! Space/Time functions
     !
     !-------------------------------------------------------------------

     if( number_space_time_function > 0 ) then
        call exchange_init()
        if( ISLAVE ) then
           do ifunc = 1,number_space_time_function
              igene = ifunc
              call ker_memory(7_ip)
           end do
        end if
        do ifunc = 1,number_space_time_function
           nexpr = space_time_function(ifunc) % nexpr
           idime = space_time_function(ifunc) % ndime
           do kdime = 1,idime
              call exchange_add(space_time_function(ifunc) % expression(kdime))
           end do
        end do
        call exchange_end()
     end if

     !-------------------------------------------------------------------
     !
     ! Time functions
     !
     !-------------------------------------------------------------------

     if( number_time_function > 0 ) then
        call exchange_init()
        if( ISLAVE ) then
           do ifunc = 1,number_time_function
              igene = ifunc
              call ker_memory(8_ip)
           end do
        end if
        do ifunc = 1,number_time_function
           idime = time_function(ifunc) % npara
           do kdime = 1,idime
              call exchange_add(time_function(ifunc) % parameters(kdime))
           end do
        end do
        call exchange_end()
     end if
     
     !-------------------------------------------------------------------
     !
     ! Windkessel functions
     !
     !-------------------------------------------------------------------
     
     if( number_windk_systems > 0 ) then
        call exchange_init()
        if( ISLAVE ) then
           do ifunc = 1,number_windk_systems
              igene = ifunc
              call ker_memory(13_ip)
           end do
        end if
        do ifunc = 1,number_windk_systems
           idime = windk_systems(ifunc) % nparam
           do kdime = 1,idime
              call exchange_add(windk_systems(ifunc) % params(kdime))
           end do
        end do
        call exchange_end()
     end if

     if( number_pump_curve > 0 ) then
        call exchange_init()
        if( ISLAVE ) then
           do ifunc = 1,number_pump_curve
              igene = ifunc
              call ker_memory(14_ip)
           end do
        end if
        do ifunc = 1,number_pump_curve
           idime = pump_curve(ifunc) % nparam
           do kdime = 1,idime
              call exchange_add(pump_curve(ifunc) % params(kdime))
           end do
        end do
        call exchange_end()
     end if
     
     if( IMASTER ) then
        do ifunc = 1,number_windk_systems
           igene = ifunc
           call ker_memory(-14_ip)
        end do
     end if

     !
     ! Parallel exchange of lookup tables
     !
     do ifunc = 1,number_lookup_tab
        ptr_lookup_tab    => lookup_tab(ifunc)
        ptr_lookup_coords => lookup_coords(:,ifunc)
        call tab_par_exchange(ptr_lookup_coords,ptr_lookup_tab)
        lookup_tab(ifunc)      = ptr_lookup_tab
        lookup_coords(:,ifunc) = ptr_lookup_coords
     enddo

     !
     ! Parallel exchange of lookup frameworks
     !
     do ifunc = 1,number_lookup_fw
        call fw_par_exchange(1_ip,lookup_fw(ifunc))
        if (lookup_fw(ifunc) % kfl_tab_main /= 0) then
           if (ISLAVE) lookup_fw(ifunc) % main_table => lookup_tab(lookup_fw(ifunc) % kfl_tab_main)
           call fw_par_exchange(2_ip,lookup_fw(ifunc))
           do ii = 1, lookup_fw(ifunc) % main_table % ndim
              if (lookup_fw(ifunc) % kfl_scale(ii) == 1) then 
                 lookup_fw(ifunc) % scaling(ii) % tab => lookup_tab(lookup_fw(ifunc) % scaling(ii) % kfl_tab)
              endif
              if (ISLAVE) call fw_allocate(2_ip,lookup_fw(ifunc),ii)
           enddo
        endif
     enddo

     !
     ! Parallel exchange of artificial neural network frameworks
     !
     do ifunc = 1,max_ann_fw
        call ann_fw(ifunc) % par_comm ()
     enddo




  case( 2_ip )
     !
     ! Second round
     !

     call eccou_send_data(1_ip) !exchnge reallocated ca50, after reading ker.dat

     !-------------------------------------------------------------------
     !
     ! Exchange kernel data after reading module data
     !
     !-------------------------------------------------------------------

     call mod_windk_systems_par_exchange()

     !-------------------------------------------------------------------
     !
     ! Tubes functions
     !
     !-------------------------------------------------------------------
     
     call tubes_parall()

  end select

  npari = 0
  nparr = 0
  nparc = 0
  nullify(parin)
  nullify(parre)

end subroutine ker_parall

