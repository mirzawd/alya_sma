!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Domain
!> @{
!> @file    mod_AMR.f90
!> @author  houzeaux
!> @date    2020-03-07
!> @brief   Adaptive mesh refinement
!> @details All tools for adaptive mesh refinement
!-----------------------------------------------------------------------

module mod_AMR

  use def_master
  use def_domain
  use mod_memory
  use mod_communications
  use mod_maths
  use mod_parall
  use def_parall
  use def_coupli
  use mod_couplings
  use mod_mesh_type
  use mod_mesh_type_basic
  use def_inpout
  use mod_htable
  use def_AMR
  use mod_debug
  
  use mod_redistribution,                      only : commd_npoin
  use mod_redistribution,                      only : commd_nelem
  use mod_redistribution,                      only : commd_nboun
  use mod_output,                              only : output_open_files
  use mod_ecoute,                              only : ecoute
  use def_kermod,                              only : ndivi
  use def_kintyp_mesh_basic,                   only : mesh_type_basic
  use def_kintyp_mesh_basic,                   only : ELEMENT_TAG
  use def_kintyp_mesh_basic,                   only : NODE_TAG
  use mod_elmgeo,                              only : element_type
  use def_elmtyp,                              only : BOFEM
  use mod_meshes,                              only : meshes_submesh
  use mod_meshes,                              only : meshes_glue_two_meshes
  use mod_meshes,                              only : meshes_list_boundary_elements
  use mod_maths,                               only : maths_heap_sort
  use mod_graphs,                              only : graphs_number_to_linked_list
  use mod_graphs,                              only : graphs_elepoi
  use mod_graphs,                              only : graphs_elepoi_deallocate
  use mod_coupling_memory,                     only : cou_initialization
  use mod_coupling_memory,                     only : cou_deallocate
  use mod_couplings_communications,            only : COU_GENERATE_LOCAL_TRANSMISSION_MATRICES
  use mod_interpolation,                       only : COU_GET_INTERPOLATE_POINTS_VALUES
  use mod_messages,                            only : messages_live
  use mod_parall_destructor,                   only : parall_destructor
  use mod_par_tools,                           only : par_tools_ownership_interfaces
  use mod_par_additional_arrays,               only : par_global_variables_arrays
  use mod_domain,                              only : domain_memory_deallocate 
  use mod_par_global_numbering,                only : par_global_numbering_nodes
  use mod_par_global_numbering,                only : par_global_numbering_elements_boundaries
  use mod_kdtree,                              only : kdtree_construct
  use mod_mpio_par_postpr,                     only : posmpio_destructor
  use mod_mpio_par_postpr,                     only : posmpio_redistribute
  use mod_mesh_type_basic,                     only : mesh_type_basic_valid_mesh
  use mod_strings,                             only : integer_to_string
  use mod_renumbering,                         only : renumbering_node_arrays  
  use mod_AMR_interpolate,                     only : AMR_interpolate
  use mod_AMR_interpolate,                     only : current_interp
  use mod_alya2gmsh,                           only : alya2gmsh_remeshing
  use mod_strings,                             only : real_to_string
  use mod_domain,                              only : domain_memory_reallocate
  use mod_par_output_partition,                only : par_output_partition
  use mod_par_output_partition,                only : par_output_solvers
  use mod_communications_point_to_point_basic, only : PAR_INTERFACE_NODE_EXCHANGE
  use def_search_method,                       only : search_method_read_data
  use def_search_method,                       only : SEARCH_BIN    
  use def_search_method,                       only : SEARCH_OCTREE 
  use def_search_method,                       only : SEARCH_SKDTREE 
  use def_interpolation_method,                only : interpolation 
  use def_maths_bin,                           only : maths_bin
  use def_maths_tree,                          only : maths_octree
  use def_maths_tree,                          only : maths_skdtree
  use def_search_method,                       only : search_method
  use def_interpolation_method,                only : INT_ELEMENT_INTERPOLATION
  use def_interpolation_method,                only : INT_BOUNDARY_INTERPOLATION
  use def_interpolation_method,                only : INT_ELEMENT_VALUE
  use def_interpolation_method,                only : INT_BOUNDARY_VALUE
  use mod_par_additional_arrays,               only : par_bounding_box
  use mod_exchange,                            only : exchange_init
  use mod_exchange,                            only : exchange_add
  use mod_exchange,                            only : exchange_end
  use def_mpi
#include "def_mpi.inc"
  
  implicit none

  private
  !
  ! Generic
  !
  character(7)             :: vacal='mod_AMR'
  integer(ip)              :: num_passes
  real(rp),  pointer       :: solut(:)
  real(rp),  pointer       :: solut_multi(:,:)
  real(rp),  pointer       :: solut1(:)
  real(rp)                 :: time1, time2
  
  logical(lg), parameter   :: export_steps_amr = .false. ! for abel to debug
  logical(lg), parameter   :: export_solut_amr = .false. ! for abel to debug
  integer(ip)              :: iplot_adapt                ! for abel to debug
  integer(ip)              :: iplot_adapt_sol            ! for abel to debug

    
  public :: AMR
  public :: AMR_initialization
  public :: AMR_parall
  public :: AMR_readat
  
contains
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-02-29
  !> @brief   Initialization
  !> @details AMR initialization
  !> 
  !-----------------------------------------------------------------------

  subroutine AMR_variable(target_solut)

    real(rp), pointer, intent(in), optional :: target_solut(:)

    integer(ip) :: ipoin, ivaria, varia_id
    
    if(AMR_numVariables==1) then
      
      if(      (kfl_amr_varia.eq.AMR_GIVEN_VARIABLE)  .and.(.not.present(target_solut))) then
        call runend('AMR should be provided of the nodal field to adapt to.')
      else if( (kfl_amr_varia.eq.AMR_QUALITY_VARIABLE).and.(.not.present(target_solut))) then
        call runend('AMR should be provided of the nodal field to adapt to.')
      else if(((kfl_amr_varia.ne.AMR_GIVEN_VARIABLE).and.(kfl_amr_varia.ne.AMR_QUALITY_VARIABLE)) &
        .and.(present(target_solut))) then
        ! This could be accepted if we want to adapt in different parts of the code to different field, but just in case
        ! That is, if desired this runend can be removed
        call runend('A target field is provided to AMR while in the input file a different one for is prescrbed')
      end if
    
      if(present(target_solut)) then
        solut => target_solut
      else if( kfl_amr_varia == ID_TEMPE ) then
         solut => tempe(:,1)
      else if( kfl_amr_varia == ID_PRESS ) then
         solut => press(:,1)
      else if( kfl_amr_varia == ID_VELOC ) then
         call memory_alloca(memor_dom,'SOLUT1',vacal,solut1,npoin)
         do ipoin = 1,npoin
            solut1(ipoin) = sqrt(dot_product(veloc(1:ndime,ipoin,1),veloc(1:ndime,ipoin,1)))
         end do
         solut => solut1
      end if
      
    else
      !call runend('Not implemented yet')
      if(present(target_solut)) then
        call runend('Not compatible with mutiple solutions')
      end if
      
      if(npoin>0_ip) then
        
        call memory_alloca(memor_dom,'solut_multi',vacal,solut_multi,npoin,AMR_numVariables)
        do ivaria=1,AMR_numVariables
          varia_id = kfl_amr_varia_multi(ivaria) 
          if(     varia_id == ID_TEMPE ) then
            solut_multi(:,ivaria) = tempe(:,1)
         else if( varia_id == ID_PRESS ) then
            solut_multi(:,ivaria) = press(:,1)
         else if( varia_id == ID_VELOC ) then
            do ipoin = 1,npoin
               solut_multi(ipoin,ivaria) = sqrt(dot_product(veloc(1:ndime,ipoin,1),veloc(1:ndime,ipoin,1)))
            end do
         else
           call runend('Not implemented variable...any scalar variable can be added (ask Abel if necessary))')
         end if
        end do
        
      end if
      
    end if
    
  end subroutine AMR_variable
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-02-29
  !> @brief   Deallocate
  !> @details AMR deallocate
  !> 
  !-----------------------------------------------------------------------

  subroutine AMR_destructor()

    call memory_deallo(memor_dom,'SOLUT1',vacal,solut1)
    nullify(solut)

    call memory_deallo(memor_dom,'solut_multi',vacal,solut_multi)

  end subroutine AMR_destructor
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-02-29
  !> @brief   Initialization
  !> @details AMR initialization
  !> 
  !-----------------------------------------------------------------------

  subroutine AMR_initialization()

    nullify(solut)
    nullify(solut1)
    num_passes           =  0
    
    AMR_numVariables     =  1_ip
    nullify(solut_multi) ! to have more than one solution to adapt to
    
    num_amr              =  0
    kfl_amr              =  0
    kfl_amr_post         =  2
    kfl_amr_remesh       =  1
    kfl_amr_repart       =  0
    kfl_amr_freq         =  1e6
    nelem_amr            =  0
    npoin_amr            =  0
    maxit_amr            =  1
    kfl_amr_varia        =  0
    kfl_size_amr         =  0
    kfl_error_amr        =  AMR_LAPLACIAN
    min_size_amr         =  0.0_rp
    max_size_amr         =  0.0_rp
    kfl_mesh_size_amr    =  SEARCH_OCTREE ! SEARCH_BIN ! 
    mesh_size_param_amr  =  20.0_rp  
    
    AMR_qualityTarget    = 0.4_rp    ! Target quality to attain (sets the region to perform adaptivity and quality parameters)
    AMR_qualityActivation= 0.4_rp
    
    AMR_bouSmoothSize    = .false.
    AMR_anisotropic      = .false.
    AMR_LPnorm           = 1_ip      ! L2 norm by default
    AMR_iniTimeStep      = 0_ip
    kfl_amr_varia_multi  =-100_ip
    
    nullify(interp_AMR_npoin)
    nullify(interp_AMR_nelem)
    nullify(interp_AMR_nboun)

    call search_vol_seq % init()
    call search_vol_par % init()
    call search_bou_seq % init()
    call search_bou_par % init()

    search_vol_par % type       = SEARCH_BIN
    search_vol_par % param(1:3) = (/20.0_rp,20.0_rp,20.0_rp/)
    search_vol_seq % type       = SEARCH_OCTREE
    search_vol_seq % param(1:1) = (/20.0_rp/)

    search_bou_par % type       = SEARCH_SKDTREE
    search_bou_par % param(1:1) = (/1.0_rp/)
    search_bou_seq % type       = SEARCH_SKDTREE
    search_bou_seq % param(1:1) = (/1.0_rp/)

  end subroutine AMR_initialization
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-02-29
  !> @brief   Initialization
  !> @details AMR initialization
  !> 
  !-----------------------------------------------------------------------

  subroutine AMR_parall()

    integer(ip) :: ii
 
    call exchange_init() 
    call exchange_add(kfl_amr)
    call exchange_add(kfl_amr_post)
    call exchange_add(kfl_amr_remesh)
    call exchange_add(kfl_amr_repart)
    call exchange_add(kfl_amr_freq)
    call exchange_add(nelem_amr)
    call exchange_add(npoin_amr)
    call exchange_add(maxit_amr)
    call exchange_add(kfl_amr_varia)
    call exchange_add(kfl_size_amr)
    call exchange_add(kfl_error_amr)
    call exchange_add(min_size_amr)
    call exchange_add(max_size_amr)
    call exchange_add(kfl_mesh_size_amr)
    call exchange_add(mesh_size_param_amr)
    
    call exchange_add(AMR_qualityTarget)
    call exchange_add(AMR_qualityActivation)
    call exchange_add(AMR_bouSmoothSize)
    call exchange_add(AMR_anisotropic)
    call exchange_add(AMR_LPnorm)
    call exchange_add(AMR_iniTimeStep)
    call exchange_add(AMR_numVariables)
    call exchange_add(kfl_amr_varia_multi)
    !call exchange_add(is_alloca_hmesh)
    !call exchange_add(hmesh)
    
    
    call exchange_add(search_vol_seq % type)
    call exchange_add(search_vol_par % type)
    call exchange_add(search_bou_seq % type)
    call exchange_add(search_bou_par % type)
    call exchange_add(search_vol_seq % param)
    call exchange_add(search_vol_par % param)
    call exchange_add(search_bou_seq % param)
    call exchange_add(search_bou_par % param)

    call exchange_end()
    
  end subroutine AMR_parall
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-02-29
  !> @brief   Read data
  !> @details Read input data
  !> 
  !-----------------------------------------------------------------------

  subroutine AMR_readat()

    integer(ip) :: ivaria
    
    if( words(2) == 'OFF  ' ) then
       do while(words(1) /= 'ENDAD' )
          call ecoute('ker_readat')       
       end do
    else
       kfl_amr = 1
       call ecoute('ker_readat')
       do while(words(1) /= 'ENDAD' )
!          print*,words(1)
          if( words(1) == 'REPAR' ) then
             !
             ! Repartitioning method
             !
             if(      words(2) == 'SFC  ' ) then
                kfl_amr_repart = 0
             else if( words(2) == 'MAX  ' ) then
                kfl_amr_repart = 1
             else
                call runend('PAR_REAPRO: UNKNOWN REPART METHOD')
             end if
          else if( words(1) == 'REMES' ) then
             !
             ! Adaptation method
             !
             if(      words(2) == 'COPY ' ) then
                kfl_amr_remesh = AMR_COPY_MESH
             else if( words(2) == 'GMSH ' ) then
                kfl_amr_remesh = AMR_GMSH
             else if( words(2) == 'ALYA ' ) then
                kfl_amr_remesh = AMR_ALYA
             else
                call runend('PAR_REAPRO: UNKNOWN REMESH METHOD')
             end if

          else if( words(1) == 'POSTP' ) then
             !
             ! Postprocess file name strategy
             !
             if(      words(2) == 'TAG  ' ) then
                kfl_amr_post = 1
             else if( words(2) == 'DIREC' ) then
                kfl_amr_post = 2
             else if( words(2) == 'NONE ' ) then
                kfl_amr_post = 0
             else
                call runend('PAR_REAPRO: UNKNOWN REPARTITION POSTPROCESS STRATEGY')
             end if

          else if( words(1) == 'FREQU' ) then
             !
             ! Partition frequency
             !
             kfl_amr_freq = getint('FREQU',10_ip,'#FREQUENCY FOR REPARTITIONING')

          else if( words(1) == 'TARGE' ) then
             !
             ! Target elements
             !
             nelem_amr = getint('TARGE',0_ip,'#TARGET NUMBER OF ELEMS')
             npoin_amr = 0_ip

          else if( words(1) == 'NODES' ) then
             !
             ! Target elements
             !
             npoin_amr = getint('NODES',0_ip,'#TARGET NUMBER OF NODES')
             nelem_amr = 0_ip
            
          else if ( words(1) == 'ERROR' ) then
             !
             ! Error estimate
             !
             if( words(2) == 'SHARP' ) then
                kfl_error_amr = AMR_SHARP_LAPLACIAN
             else if( words(2) == 'LAPLA' ) then
                kfl_error_amr = AMR_LAPLACIAN
             else if( words(2) == 'DISTA' ) then
                kfl_error_amr = AMR_DISTANCE_POINT
                if (ndime == 2) then
                   amrp0(1) = getrea('amrp0_x',0.0_rp,'#amrp0_x')
                   amrp0(2) = getrea('amrp0_y',0.0_rp,'#amrp0_y')
                else if (ndime == 3) then
                   amrp0(1) = getrea('amrp0_x',0.0_rp,'#amrp0_x')
                   amrp0(2) = getrea('amrp0_y',0.0_rp,'#amrp0_y')
                   amrp0(3) = getrea('amrp0_z',0.0_rp,'#amrp0_z')
                end if
             else if(words(2) =='HESSI') then
               kfl_error_amr = AMR_HESSIAN
             end if

          else if( words(1) == 'VARIA' ) then
            !
            ! Variable to adapt to (or target of the adaptation -quality repair-)
            !
            ! nwopa or nnwor are the number of words (they seem to have the same value, check difference...)
            AMR_numVariables = nwopa-1
            if(AMR_numVariables==1) then
              kfl_amr_varia = getVaria(words(2))
!               if(      words(2) == 'TEMPE' ) then
!                   kfl_amr_varia = ID_TEMPE
!                else if( words(2) == 'VELOC' ) then
!                   kfl_amr_varia = ID_VELOC
!                else if( words(2) == 'PRESS' ) then
!                   kfl_amr_varia = ID_PRESS
!                else if( words(2) == 'GIVEN' ) then ! variable given directly to AMR from code
!                   kfl_amr_varia = AMR_GIVEN_VARIABLE
!                else if( words(2) == 'QUALI' ) then ! optimize without targetting a field, to repair mesh according to size and quality
!                   kfl_amr_varia = AMR_QUALITY_VARIABLE
!                else
!                   call runend('MOD_AMR: UNKNOWN VARIABLE')
!                end if
             else
               do ivaria=1,AMR_numVariables
                 kfl_amr_varia_multi(ivaria) = getVaria(words(ivaria+1))
               end do
             end if

          else if( words(1) == 'QUALT' ) then
            !
            ! Target quality
            !
            AMR_qualityTarget = getrea('QUALT',0.0_rp,'#TARGET QUALITY OF THE ADAPTATION PROCESS [0,1]')
            !print*,'   AMR_qualityTarget: ', AMR_qualityTarget, "  ", words(1),"  ", words(2)
            !
          else if( words(1) == 'QUALA' ) then
            !
            ! Target quality
            !
            AMR_qualityActivation = getrea('QUALA',0.0_rp,'#QUALITY TO ACTIVATE THE ADAPTATION PROCESS [0,1]')
            !print*,'   AMR_qualityActivation: ', AMR_qualityActivation, "  ", words(1),"  ", words(2)
            !
          else if( words(1) == 'BOUSM' ) then
            !
            ! Do smooth mesh field with current boundary or not
            !
            AMR_bouSmoothSize = words(2) == 'ON   '
            !
            !
            !print*,'   bound size smoothing: ', AMR_bouSmoothSize, "  ", words(2)

          else if( words(1) == 'TIMES' ) then
            !
            ! Time step to start adaptation
            !
            AMR_iniTimeStep = getint('TIMES',0_ip,'Time step to start adaptation')
            !
            !
            !print*,'   AMR_iniTimeStep: ', AMR_iniTimeStep

          else if( words(1) == 'LPNOR' ) then
             !
             ! Hp norm to define metric
             !
             AMR_LPnorm = getint('LPNOR',0_ip,'P TO SET HPNORM')
             
          else if( words(1) == 'ANISO' ) then
            !
            ! Choose if anisotropic adaptation or not
            !
            AMR_anisotropic = words(2) == 'ON   '
            
            if( AMR_anisotropic ) then
              call runend('AMR_anisotropic not allowed yet... on development')
            end if
            !
            ! 
          else if( words(1) == 'ITERA' ) then
             !
             ! Iterations
             !
             maxit_amr = getint('ITERA',1_ip,'#AMR NUMBER OF ITERATIONS')

          else if( words(1) == 'SIZES' ) then
             !
             ! Size strategy
             !
             if( words(2) == 'AUTOM' ) then
                kfl_size_amr = 0
             else if( words(2) == 'PRESC' ) then
                kfl_size_amr = 1
             end if
             if( exists('MINIM') ) min_size_amr = getrea('MINIM',0.0_rp,'#MINIMUM MESH SIZE')
             if( exists('MAXIM') ) max_size_amr = getrea('MAXIM',0.0_rp,'#MAXIMUM MESH SIZE')

          else if( words(1) == 'MESHS' ) then
             !
             ! MESH size interpolation strategy
             !
             call search_method_read_data(words,param,kfl_mesh_size_amr,mesh_size_param_amr)
             
          else if( words(1) == 'VOLSE' ) then
             !
             ! Volume interpolation strategy
             !
             call search_vol_seq % read_data(words,param)
             
          else if( words(1) == 'VOLPA' ) then
             !
             ! Volume interpolation strategy
             !
             call search_vol_par % read_data(words,param)
             
          else if( words(1) == 'BOUSE' ) then
             !
             ! Volume interpolation strategy
             !
             call search_bou_seq % read_data(words,param)
             
          else if( words(1) == 'BOUPA' ) then
             !
             ! Volume interpolation strategy
             !
             call search_bou_par % read_data(words,param)
                                   
          end if
          call ecoute('ker_readat')
       end do
    end if

  end subroutine AMR_readat
  !
  !
  !
  function getVaria(word) result(varia)
    character(5), intent(in) :: word
    integer(ip)              :: varia
    
    if(       word == 'TEMPE' ) then
        varia = ID_TEMPE
     else if( word == 'VELOC' ) then
        varia = ID_VELOC
     else if( word == 'PRESS' ) then
        varia = ID_PRESS
     else if( word == 'GIVEN' ) then ! variable given directly to AMR from code
        varia = AMR_GIVEN_VARIABLE
     else if( word == 'QUALI' ) then ! optimize without targetting a field, to repair mesh according to mesh size and quality
        varia = AMR_QUALITY_VARIABLE
     else
       varia = -100_ip
       call runend('MOD_AMR: UNKNOWN VARIABLE')
     end if
    
  end function getVaria
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-02-29
  !> @brief   Main driver
  !> @details Perform iterations
  !> 
  !-----------------------------------------------------------------------

  function AMR_tags(mesh) result(res)
    type(mesh_type_basic), intent(in) :: mesh  
    integer(ip)                       :: res
    real(rp),              pointer    :: tag(:)
    integer(ip)                       :: ipoin
    real(rp)                          :: resr
    
    nullify(tag)
    call memory_alloca(par_memor,'TAG',vacal,tag,mesh % npoin,INIT_VALUE=1.0_rp)
    call PAR_INTERFACE_NODE_EXCHANGE(tag,'SUM',mesh % comm)
    resr = 0.0_rp
    if( associated(mesh % tags) ) then
       do ipoin = 1,mesh % npoin
          resr = resr + mesh % tags(1) % values(ipoin) / tag(ipoin)
       end do
    end if
    call memory_deallo(par_memor,'TAG',vacal,tag)
    
    call PAR_SUM(resr)
    res = int(resr,ip)
    
  end function AMR_tags
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-02-29
  !> @brief   Main driver
  !> @details Perform iterations
  !> 
  !-----------------------------------------------------------------------

  subroutine AMR(target_solut_input,maxMeshSize_input,is_solut_mesh_size_input)
  
    use mod_out_paraview,             only: out_paraview_inp
    use mod_communications_global,    only: PAR_MAX
    use def_adapt,                    only: memor_adapt
    use mod_metricComputation,        only: compute_sizeField_mesh
    use mod_metricComputation,        only: compute_sizeField_smoothedBoundary
    use def_master,                   only: ittim

    use mod_error_estimate,           only : error_estimate_mesh_size
    use mod_error_estimate,           only : set_mesh_size
    use mod_error_estimate,           only : error_estimate_mesh_size_multipleSol
    
    real(rp), pointer, intent(in), optional :: target_solut_input(:)    ! scalar field to adapt the mesh 
    real(rp),          intent(in), optional :: maxMeshSize_input  ! max allowed mesh size (if not provided, will be computed inside)
    logical(lg),       intent(in), optional :: is_solut_mesh_size_input ! target_solut is nodal mesh size, not solution

    type(mesh_type_basic)        :: mesh_new         ! Adapted mesh
    type(mesh_type_basic)        :: mesh_cpy         ! Copy of the mesh
    type(mesh_type)              :: mesh_complete    ! Complete adapted mesh
    integer(ip)                  :: ipass
    integer(ip)                  :: ii,jj,kk
    integer(ip)                  :: num_tags,mowners
    integer(ip),         pointer :: lnod_owners(:,:)
    integer(ip),         pointer :: rank_nelem(:)
    real(rp),            pointer :: hh_opt(:)
    real(rp),            pointer :: hh(:)
    
    real(rp), pointer :: target_solut(:)
    real(rp)          :: maxMeshSize 
    logical(lg)       :: is_solut_mesh_size
    logical(lg)       :: hasBeenModified, is_modified
    logical(lg)       :: is_interp_allocated
    
    real(rp) :: hmin_debug, hmax_debug
    logical(lg), parameter :: hsize_debug = .false.
    real(rp) :: t_perf0,t_perf1
    logical(lg), parameter :: perform_debug = .false.
    
    real(rp), pointer :: nodeFields_export(:,:)
    
    
    ! To see if mes is invalid before starting.. check mesh
    !  - Necessary when using ALE.. since the boundary deformation/movement can invalidate the mesh...
    !call mescek(1_ip)
    
    hasBeenModified = .false.
    
    if(ittim<AMR_iniTimeStep) then
      return
    end if
    
    if(num_amr==0_ip) then
      iplot_adapt = 0_ip
      iplot_adapt_sol = 0_ip
    end if
    
    if(perform_debug) call cputim(t_perf0)
    
    num_passes = num_passes + 1

    if( num_passes == kfl_amr_freq .and. kfl_amr /= 0 ) then
       
       nullify(rank_nelem)
       nullify(hh_opt)
       nullify(lnod_owners)
       
       num_amr = num_amr + 1
       
       !-----------------------------------------------------------------
       !
       ! Copy original mesh MESH_NEW = MESHE(NDIVI)
       !
       !-----------------------------------------------------------------

       call mesh_complete % init      ('MESH_COMPLETE')
       
       call mesh_new      % init      ('MESH_NEW')
       call mesh_new      % copy      (meshe(ndivi))
       call mesh_new      % copy_comm (meshe(ndivi))
              
       call messages_live('---')
       call messages_live('ADAPTIVE MESH REFINEMENT '//integer_to_string(num_amr),'START SECTION')

       if( kfl_amr_remesh == AMR_COPY_MESH ) then
          call messages_live('LOOKS LIKE YOU ARE DEBUGGING, MESH WILL NOT CHANGE AS YOU ARE USING THE COPY MESH OPTION','WARNING')
       elseif( kfl_amr_remesh == AMR_ALYA) then
          call messages_live('NATIVE ALYA ADAPTIVE STRATEGY SELECTED')
       elseif( kfl_amr_remesh == AMR_GMSH) then
#ifndef ALYA_GMSH
          call messages_live('GMSH IS NOT ACTIVATED, COMPILE WITH -DALYA_GMSH','WARNING')
#endif
       else
         call runend('NOT KNOWN REMESHING STARTEGY')
       end if
       
       
       if( kfl_amr_remesh.eq.AMR_ALYA ) then
         ! Deal with different input option combinations (maybe simplify and deprecate some in the future)
         
         if(present(target_solut_input)) then
           nullify(target_solut)
           target_solut => target_solut_input
         end if
         
         if(present(is_solut_mesh_size_input)) then
           is_solut_mesh_size = is_solut_mesh_size_input
         else
           is_solut_mesh_size = .false.
         end if
         
         if( kfl_amr_varia == AMR_QUALITY_VARIABLE ) then

           if(mesh_new%npoin>0_ip) then
             call compute_sizeField_mesh(hh,mesh_new)
             maxMeshSize = maxval(hh)
             !call memory_deallo(memor_adapt,'hh','compute_sizeField_mesh',hh) ! LATER DEALLOCATED
           else
             maxMeshSize = 0.0_rp
           end if
           
           call PAR_MAX(maxMeshSize)
           maxMeshSize = AMR_factorMaxSize*maxMeshSize  ! max size after adaptivity bounded as XX times the current mesh size
           
           nullify(target_solut)
           target_solut => hh
           
           is_solut_mesh_size = .true.
           
         else
           ! set max size
           if(present(maxMeshSize_input)) then
             maxMeshSize = maxMeshSize_input
           else
             if(mesh_new%npoin>0_ip) then
               call compute_sizeField_mesh(hh,mesh_new)
               maxMeshSize = maxval(hh)
               call memory_deallo(memor_adapt,'hh','compute_sizeField_mesh',hh)
             else
               maxMeshSize = 0.0_rp
             end if
             call PAR_MAX(maxMeshSize)
             maxMeshSize = AMR_factorMaxSize*maxMeshSize  ! max size after adaptivity bounded as XX times the current mesh size
           end if
         
         end if
         
       else
         
         if(present(target_solut_input)) then
           nullify(target_solut)
           target_solut => target_solut_input
         end if
         
         if(present(is_solut_mesh_size_input)) then
           is_solut_mesh_size = is_solut_mesh_size_input
           if(is_solut_mesh_size) then
             call runend('Methods other than AMR_ALYA not implemented for this option')
           end if
         else
           is_solut_mesh_size = .false.
         end if
         
       end if
       
       !-----------------------------------------------------------------
       !
       ! Compute mesh size HH_OPT
       !
       !-----------------------------------------------------------------
       if( nelem_amr == 0 .and. npoin_amr==0) nelem_amr = nelem_total
       if(present(target_solut_input).or.is_solut_mesh_size) then ! improve this double call.. some chaos with input to AMR
         call AMR_variable(target_solut)
       else
         call AMR_variable()
       end if
       
       if(export_solut_amr.and.(npoin>0_ip)) then ! for debugging purposes
         iplot_adapt_sol = iplot_adapt_sol+1_ip
         !call out_paraview_inp(mesh_new,&
         !  filename=integer_to_string(kfl_paral)//"_"//integer_to_string(iplot_adapt_sol),nodeField=solut)
!          if(AMR_numVariables==1_ip) then
!            call out_paraview_inp(mesh_new,filename='/Users/abel/Desktop/Work/CasesAlya/paraview/'//&
!                integer_to_string(kfl_paral)//"meshSolAd_"//integer_to_string(iplot_adapt_sol),nodeField=solut)
!          else
!            call out_paraview_inp(mesh_new,filename='/Users/abel/Desktop/Work/CasesAlya/paraview/'//&
!                integer_to_string(kfl_paral)//"meshSolAd_"//integer_to_string(iplot_adapt_sol),nodeFields=solut_multi)
!          end if
         ! To force plot ||veloc||
         !print*,meshe(ndivi)%kfl_codbo
         !print*,tgcod(1) % l(1) % kfl_value
         !print*,tgcod(1) % l(2) % kfl_value
!          print*,'size(meshe(ndivi)%kfl_codno ',size(meshe(ndivi)%kfl_codno,1),size(meshe(ndivi)%kfl_codno,2)
!          print*,"veloc: ",size(veloc,1),size(veloc,2)
!          print*,"codbo: ",size(meshe(ndivi)%kfl_codbo(:))
!          !print*,"size(lnset): ",size(meshe(ndivi)%lnset)
!          !print*,"size(lbset): ",size(meshe(ndivi)%lbset)
!          !print*,size(meshe(ndivi)%lnset)
!          print*,"size(meshe(ndivi)%tags): ",size(meshe(ndivi)%tags)
!          !print*,meshe(ndivi)%tags(1)%name
!          !print*,associated(meshe(ndivi)%tags(1)%values)
!          !print*,"size(meshe(ndivi)%tags(1)%values): ",size(meshe(ndivi)%tags(1)%values)
!          print*,'size(mesh_new%tags): ',size(mesh_new%tags)
!
!          print*,' '
!          print*,'size(lnodb)',size(meshe(ndivi)%lnodb,1),' ',size(meshe(ndivi)%lnodb,2)
!          print*,'associated(meshe(ndivi)%kfl_codbo): ',associated(meshe(ndivi)%kfl_codbo)
!          print*,'associated(meshe(ndivi)%lbset): ',associated(meshe(ndivi)%lbset)
         nullify(nodeFields_export)
         call memory_alloca(memor_dom,'nodeFields_export','mod_AMR',nodeFields_export,INT(size(veloc,2),ip),2_ip)
         !nodeFields_export(:,1) = sqrt(SUM(veloc(1:ndime,:,1)**2,1))
         !nodeFields_export(:,2) = kfl_codno(1,:)!press(:,1)!1.0_rp*meshe(ndivi)%kfl_codno(:) !kfl_codbo
         nodeFields_export(:,1) = sqrt(SUM(veloc(1:ndime,:,1)**2,1))
         nodeFields_export(:,2) = press(:,1)!kfl_codno(1,:)!press(:,1)!1.0_rp*meshe(ndivi)%kfl_codno(:) !kfl_codbo
         !nodeFields_export(:,2) = meshe(ndivi)%lnset !press(:,1)!1.0_rp!*mesh_complete%kfl_codno(:)!meshe(ndivi)%kfl_codno !kfl_codbo
         call out_paraview_inp(mesh_new,&
           filename='./'//integer_to_string(kfl_paral)//"meshSolAd_"//integer_to_string(iplot_adapt_sol),nodeFields=nodeFields_export)
         call memory_deallo(memor_dom,'nodeFields_export','mod_AMR',nodeFields_export)
!          call out_paraview_inp(mesh_new,&
!            filename=&'./'//&!'/Users/abel/Desktop/Work/CasesAlya/paraview/'//&
!              integer_to_string(kfl_paral)//"meshSolAd_"//integer_to_string(iplot_adapt_sol),nodeField=sqrt(SUM(veloc(1:ndime,:,1)**2,1)))
             
       end if
       
       ! To see if mes is invalid before starting.. check mesh
       !  - Necessary when using ALE.. since the boundary deformation/movement can invalidate the mesh...
       ! - i put it here to make it break after exporting, but should be at the beginning
       call mescek(1_ip)
       
       if( kfl_amr_remesh.eq.AMR_ALYA ) then
         
         if(is_solut_mesh_size) then
           !
           if(npoin_amr>0_ip) then
             call set_mesh_size(solut,hh_opt,&
              NPOIN_OPT=npoin_amr,maxMeshSize=maxMeshSize,isNodalField=.true.)
            else
             call set_mesh_size(solut,hh_opt,&
              NELEM_OPT=nelem_amr,maxMeshSize=maxMeshSize,isNodalField=.true.)
           end if
           !
         else
           
           if(AMR_numVariables==1_ip) then
             if(npoin_amr>0_ip) then
               call error_estimate_mesh_size(solut,hh_opt,&
                NPOIN_OPT=npoin_amr,METHOD=kfl_error_amr,maxMeshSize=maxMeshSize,isNodalField=.true.)
             else
               call error_estimate_mesh_size(solut,hh_opt,&
                NELEM_OPT=nelem_amr,METHOD=kfl_error_amr,maxMeshSize=maxMeshSize,isNodalField=.true.)
             end if
           else
             if(npoin_amr>0_ip) then
               call error_estimate_mesh_size_multipleSol(solut_multi,hh_opt,&
                NPOIN_OPT=npoin_amr,METHOD=kfl_error_amr,maxMeshSize=maxMeshSize)
             else
               call error_estimate_mesh_size_multipleSol(solut_multi,hh_opt,&
                NELEM_OPT=nelem_amr,METHOD=kfl_error_amr,maxMeshSize=maxMeshSize)
             end if
           end if
           
         end if
         
         if( kfl_amr_varia == AMR_QUALITY_VARIABLE ) then
           call memory_deallo(memor_adapt,'hh','compute_sizeField_mesh',hh)
         end if
         
       else
         call error_estimate_mesh_size(solut,hh_opt,NELEM_OPT=nelem_amr,METHOD=kfl_error_amr,maxMeshSize=maxMeshSize,isNodalField=.false.)
       end if
     
       ! Here we should add the possibility to gradate the size field according to the boundary mesh size
       ! 1- compute mesh size hh
       ! 2- update hh_opt with respect to hh on the boundary
       ! 3- gradate hh_opt again
       if(AMR_bouSmoothSize) then
         if(.not.is_alloca_hmesh) then
           is_alloca_hmesh=.true.  ! use this to not re-compute size field if it has been computed already, for instance to set maxSize TODO
           call compute_sizeField_mesh(hmesh,mesh_new)
           
           call compute_sizeField_smoothedBoundary(mesh_new,hmesh,hh_opt)

           is_alloca_hmesh=.false.
           call memory_deallo(memor_adapt,'hh','compute_sizeField_mesh',hmesh)
         end if
       end if
       !-----------------------------------------------------------------
       !
       ! Set tags for frozen nodes, mesh size and global numbering
       !
       !-----------------------------------------------------------------

       call AMR_set_tags(mesh_new,hh_opt)
       call memory_deallo(par_memor,'HH_OPT',vacal,hh_opt)
       
       !-----------------------------------------------------------------
       !
       ! Adaptivity iterations
       !
       !-----------------------------------------------------------------
       call PAR_BARRIER()
       call cputim(time1)
       call messages_live('PARALLEL ITERATIONS','START SECTION')
       
       ipass    = 0
       num_tags = 1
       do while( ipass < maxit_amr .and. num_tags /= 0 )
          ipass = ipass + 1
          call messages_live('ITERATION '//integer_to_string(ipass),'START SECTION')
          !
          ! Adaptation, migration strategy and repartitioning
          !
          !block
          !  if( mesh_new % nelem > 0 ) then
          !     call mesh_new % output  (filename='mesh-before-'//integer_to_string(kfl_paral))
          !     call mesh_new % results (xx=mesh_new % tags(2) % values,names='TAGSS',filename='mesh-before-'//integer_to_string(kfl_paral),where='ON ELEMENTS')
          !  end if
          !end block
          
          if( kfl_amr_remesh.eq.AMR_ALYA ) then
            if(export_steps_amr.and.(npoin>0_ip)) then ! for debugging purposes
              iplot_adapt = iplot_adapt+1_ip
              call out_paraview_inp(mesh_new,&
                filename=integer_to_string(kfl_paral)//"_"//integer_to_string(iplot_adapt),nodeField=mesh_new%tags(2)%values)
            end if
            
            if(hsize_debug) then
              if(mesh_new%npoin>0_ip) then
                hmin_debug = minval(mesh_new % tags(2) % values)
                hmax_debug = maxval(mesh_new % tags(2) % values)
              else
                hmin_debug = huge(hmin_debug)
                hmax_debug = tiny(hmax_debug)
              end if
              call PAR_BARRIER()
              !print*,kfl_paral,' --->  hmin_debug: ',hmin_debug,'   hmax_debug: ',hmax_debug
              call PAR_MIN(hmin_debug)
              call messages_live(' ====> hmin = '//real_to_string(hmin_debug))
              call PAR_MAX(hmax_debug)
              call messages_live(' ====> hmax = '//real_to_string(hmax_debug))
            end if
            
            call times(2) % ini() ; call AMR_adapt_mesh_ALYA(mesh_new,mesh_cpy,is_modified) ; call times(2) % add()
            
            !is_modified = .true.
!             if(.not.is_modified) then
!               ipass = maxit_amr !dont iterate more
!             end if
            
!             !is_modified = .true.
!             hasBeenModified = hasBeenModified.or.is_modified
!             !if(.not.hasBeenModified) then
!             if( (ipass==1_ip).and.(.not.hasBeenModified) ) then
!              call messages_live('MESH QUALITY IS OK: MESH NOT ADAPTED')
!              call messages_live('ITERATION','END SECTION')
!              call messages_live('PARALLEL ITERATIONS','END SECTION')
!              go to 666
!             end if
            hasBeenModified = hasBeenModified.or.is_modified
            if(.not.hasBeenModified) then
              is_interp_allocated = .false.
              call messages_live('MESH QUALITY IS OK: MESH NOT ADAPTED')
              call messages_live('ITERATION','END SECTION')
              call messages_live('PARALLEL ITERATIONS','END SECTION')
              go to 666
            else if(hasBeenModified.and.(.not.is_modified)) then
              call messages_live('MESH QUALITY IS OK: MESH NOT ADAPTED')
              ipass = maxit_amr !dont iterate more
              go to 777
            else
              is_interp_allocated = .true.
            end if
            
             
            if(hsize_debug) then !mesh_cpy (remove this if when solved issue)
              if(mesh_cpy%npoin>0_ip) then
                hmin_debug = minval(mesh_cpy % tags(2) % values)
                hmax_debug = maxval(mesh_cpy % tags(2) % values)
              else
                hmin_debug = huge(hmin_debug)
                hmax_debug = tiny(hmax_debug)
              end if
              call PAR_BARRIER()
              !print*,kfl_paral,' --->  hmin_debug: ',hmin_debug,'   hmax_debug: ',hmax_debug
              call PAR_MIN(hmin_debug)
              call messages_live(' ====> previous hmin = '//real_to_string(hmin_debug))
              call PAR_MAX(hmax_debug)
              call messages_live(' ====> previous hmax = '//real_to_string(hmax_debug))
            end if
            
            !!!***the size interpolation should not be necessary.. I already have the sizes on my mesh..
            !print*,'I did reactivate the node sizes to see if I find the negative'
            if(hsize_debug) then
              print*,'This is my result finding the ids in my mesh'
              if(mesh_new%npoin>0_ip) then
                hmin_debug = minval(mesh_new % tags(2) % values)
                hmax_debug = maxval(mesh_new % tags(2) % values)
              else
                hmin_debug = huge(hmin_debug)
                hmax_debug = tiny(hmax_debug)
              end if
              call PAR_BARRIER()
              !print*,kfl_paral,' --->  hmin_debug: ',hmin_debug,'   hmax_debug: ',hmax_debug
              call PAR_MIN(hmin_debug)
              call messages_live(' ====> hmin = '//real_to_string(hmin_debug))
              call PAR_MAX(hmax_debug)
              call messages_live(' ====> hmax = '//real_to_string(hmax_debug))
              
              if(hmin_debug<1.0E-10) then
                call runend('Not ok... mesh size should not go to zero... Interpolation or invalid mesh?')
              end if
            end if
            !call times(8) % ini() ; call AMR_interpolate_size_nodes (mesh_new,mesh_cpy)                ; call times(8) % add()
 
            if(export_steps_amr.and.(npoin>0_ip)) then ! for debugging purposes
              iplot_adapt = iplot_adapt+1_ip
              call out_paraview_inp(mesh_new,&
                filename=integer_to_string(kfl_paral)//"_"//integer_to_string(iplot_adapt),nodeField=mesh_new%tags(2)%values)
            end if
            
            if(hsize_debug) then
              print*,'This is after AMR_interpolate_size_nodes'
              if(mesh_new%npoin>0_ip) then
                hmin_debug = minval(mesh_new % tags(2) % values)
                hmax_debug = maxval(mesh_new % tags(2) % values)
              else
                hmin_debug = huge(hmin_debug)
                hmax_debug = tiny(hmax_debug)
              end if
              call PAR_BARRIER()
              !print*,kfl_paral,' --->  hmin_debug: ',hmin_debug,'   hmax_debug: ',hmax_debug
              call PAR_MIN(hmin_debug)
              call messages_live(' ====> hmin = '//real_to_string(hmin_debug))
              call PAR_MAX(hmax_debug)
              call messages_live(' ====> hmax = '//real_to_string(hmax_debug))
              
              if(hmin_debug<1.0E-10) then
                call runend('Not ok... mesh size should not go to zero... Interpolation or invalid mesh?')
              end if
            end if
            
          else
            
            hasBeenModified     = .true. ! assume mesh is modified 
            is_interp_allocated = .true.
            
            call times(2) % ini() ; call AMR_adapt_mesh        (mesh_new,mesh_cpy)                     ; call times(2) % add()
            call times(8) % ini() ; call AMR_interpolate_size  (mesh_new,mesh_cpy)                     ; call times(8) % add()
            
          end if
          
          777 continue
          
          call times(6) % ini() ; call AMR_global_numbering  (mesh_new)                                ; call times(6) % add() 
          call times(3) % ini() ; call AMR_repartitioning    (mesh_new,rank_nelem,ipass)               ; call times(3) % add()
          call times(5) % ini() ; call PAR_communicator_pre  (mesh_new,rank_nelem,lnod_owners,mowners) ; call times(5) % add() 
          call times(4) % ini() ; call AMR_new_mesh          (rank_nelem,mesh_new)                     ; call times(4) % add()
          call times(5) % ini() ; call PAR_communicator_post (mesh_new,lnod_owners,mowners)            ; call times(5) % add() 

          if(export_steps_amr.and.(npoin>0_ip)) then ! for debugging purposes
            iplot_adapt = iplot_adapt+1_ip
            call out_paraview_inp(mesh_new,&
              filename=integer_to_string(kfl_paral)//"_"//integer_to_string(iplot_adapt),nodeField=mesh_new%tags(2)%values)
          end if

          !if( mesh_new % nelem > 0 ) then
          !   block
          !     call mesh_new % output  (filename='mesh-after-'//integer_to_string(kfl_paral))
          !     call mesh_new % results (xx=mesh_new % tags(2) % values,names='TAGSS',filename='mesh-after-'//integer_to_string(kfl_paral),where='ON ELEMENTS')
          !   end block
          !end if
          !call runend('O.K.!')
          !
          ! Tags
          !
          num_tags = AMR_tags(mesh_new)
          call messages_live('NUMBER OF FIXED NODES= '//integer_to_string(num_tags))
          call messages_live('ITERATION','END SECTION')
          call memory_deallo(memor_dom,'RANK_NELEM',vacal,rank_nelem)

          !if( islave ) then
          !   call mesh_new % output (filename='tags-'//integer_to_string(kfl_paral)//'-'//integer_to_string(ipass))
          !   call mesh_new % results(xx=mesh_new % tags(1) % values,names='TAGSS',filename='tags-'//integer_to_string(kfl_paral)//'-'//integer_to_string(ipass),where='ON NODES')
          !end if
          
       end do

       call PAR_BARRIER()
       call cputim(time2)
       call messages_live('PARALLEL ITERATIONS'//', TIME = '//real_to_string(time2-time1),'END SECTION')
       
       !-----------------------------------------------------------------
       !
       ! Create new mesh
       !
       !-----------------------------------------------------------------
       !
       ! Output
       !
       call messages_live('NEW MESH HAS BEEN GENERATED:')
       if( kfl_amr_remesh.eq.AMR_ALYA ) then
         ii = mesh_new % npoin
         call PAR_SUM(ii)
         jj = mesh_new % npoin
         call PAR_MIN(jj)
         kk = mesh_new % npoin
         call PAR_MAX(kk)
         call messages_live('  - TOTAL NUMBER OF NODES= '//integer_to_string(ii))
         call messages_live('  - MIN   NUMBER OF NODES= '//integer_to_string(jj))
         call messages_live('  - MAX   NUMBER OF NODES= '//integer_to_string(kk))       
       else
         ii = mesh_new % nelem
         call PAR_SUM(ii)
         jj = mesh_new % nelem
         call PAR_MIN(jj)
         kk = mesh_new % nelem
         call PAR_MAX(kk)
         call messages_live('  - TOTAL NUMBER OF ELEMENTS= '//integer_to_string(ii))
         call messages_live('  - MIN   NUMBER OF ELEMENTS= '//integer_to_string(jj))
         call messages_live('  - MAX   NUMBER OF ELEMENTS= '//integer_to_string(kk))       
       end if
       !
       ! Point complete mesh to basic one 
       !
       call mesh_type_basic_to_complete(mesh_new,mesh_complete)    ! MESH_COMPLETE => MESH_NEW       
       !
       ! Create boundary mesh
       !       
       call AMR_boundary_mesh(mesh_complete)

       !-----------------------------------------------------------------
       !
       ! Destruct and recreate domain 
       !
       !-----------------------------------------------------------------
       !
       ! Partitiong redistribution communicators
       !
       call commd_npoin % deallo(par_memor,COMM_NAME='COMMD_NPOIN')
       call commd_nelem % deallo(par_memor,COMM_NAME='COMMD_NPOIN')
       call commd_nboun % deallo(par_memor,COMM_NAME='COMMD_NPOIN')
       !
       ! Output parallel info
       !
       call par_output_partition(2_ip)
       call par_output_solvers()       
       !
       ! Update domain
       !
       call times(7) % ini()
       call AMR_update_communication_arrays(mesh_complete)
       call AMR_domain                     (mesh_complete)         ! Missing arrays: LNOCH, KFL_CODNO, etc.
       call times(7) % add()
       call AMR_update_original_mesh       (mesh_complete)         ! COORD, LNODS, etc. = MESH % COMPLETE
       call AMR_own_nodes() 
       call AMR_destructor()
       call mesh_type_update_last_mesh()
       call mesh_new      % deallo()
       call mesh_complete % null_arrays()
       call mesh_complete % deallo()
       ! 
       ! Destroy
       !
       call memory_allocation_mode('DEALLOCATE BEFORE ALLOCATING')    
       call messages_live('DESTRUCT SOLVERS')
       call solver_destructor()
       call messages_live('DESTRUCT DOMAIN')
       call coupli_destructor()
       call domain_destructor()
       call parall_destructor(COMM_MY_CODE_ARRAY=.false.)
       call posmpio_destructor()
       call output_open_files(AMR=.true.)

       call posmpio_redistribute(FORCE_GENERATE_COMMUNICATOR=.true.)

       call par_global_variables_arrays()
       !
       ! Reconstruct domain and restart run
       !   
       call Domain()
       !
       ! Interpolate
       !
       call PAR_BARRIER()
       call messages_live('INTERPOLATION MODULE ARRAYS','START SECTION')
       call cputim(time1)          
       call Interp()    
       call cputim(time2);time2 = time2-time1; call PAR_MAX(time2) ; call cputim(time1)
       call messages_live('INTERPOLATION MODULE ARRAYS'//', TIME = '//real_to_string(time2),'END SECTION')
       !
       ! Beginning of the run
       !
       call Solmem()  
       call Begrun()
       
       666 continue
       
       call memory_allocation_mode('DO NOT DEALLOCATE BEFORE ALLOCATING')

       call messages_live('ADAPTIVE MESH REFINEMENT','END SECTION')
       call messages_live('---')
       
       !--------------------------------------------------------------------
       !
       ! Deallocate interpolation
       !
       !--------------------------------------------------------------------
       
       call search_vol_seq % deallo()
       call search_vol_par % deallo()
       call search_bou_seq % deallo()
       call search_bou_par % deallo()
       if(is_interp_allocated) then
         deallocate(interp_AMR_npoin)
         deallocate(interp_AMR_nelem)
         deallocate(interp_AMR_nboun)
       end if
       if(hasBeenModified) then
         if(associated(velom)) velom= 0.0_rp
       else
         ! mesh not modified, skipped everything in 666
         ! then, AMR module has to be destroyed (but not dmain, solver, etc)
         call AMR_destructor()
       end if

       num_passes = 0
       
    end if

    
    if(perform_debug) then
      call cputim(t_perf1)
      print*,'****Total time: ',t_perf1-t_perf0
    end if

    if(export_solut_amr.and.(npoin>0_ip)) then ! for debugging purposes
      !print*,'===> NOT COMMIT THIS in merge with master.. will fail for not simplicial meshes'
      iplot_adapt_sol = iplot_adapt_sol+1_ip
!       if(AMR_numVariables==1_ip) then
!         call out_paraview_inp(mesh_new,&
!           filename=&
!             '/Users/abel/Desktop/Work/CasesAlya/paraview/'//&
!             integer_to_string(kfl_paral)//"meshSolAd_"//integer_to_string(iplot_adapt_sol),nodeField=solut)
!       else
!         call out_paraview_inp(mesh_new,&
!           filename=&
!             '/Users/abel/Desktop/Work/CasesAlya/paraview/'//&
!             integer_to_string(kfl_paral)//"meshSolAd_"//integer_to_string(iplot_adapt_sol),nodeFields=solut_multi)
!       end if
      ! To force plot ||veloc||
      call mesh_new      % init      ('MESH_NEW')
      call mesh_new      % copy      (meshe(ndivi))
!       call out_paraview_inp(mesh_new,&
!         filename=&
!           './'//&!'/Users/abel/Desktop/Work/CasesAlya/paraview/'//&
!           integer_to_string(kfl_paral)//"meshSolAd_"//integer_to_string(iplot_adapt_sol),nodeField=sqrt(SUM(veloc(1:ndime,:,1)**2,1)))
      !
!       print*,meshe(ndivi)%kfl_codbo
!      print*,size(kfl_codno,1),size(kfl_codno,1)
!       print*,mesh_complete%lnset
!       print*,'size(mesh_complete%lnset): ',size(mesh_complete%lnset)
!       print*,'size(meshe(ndivi)%lnset): ',size(meshe(ndivi)%lnset)
!       print*,"veloc: ",size(veloc,1),size(veloc,2)
!       print*,"codbo: ",size(mesh_complete%kfl_codbo(:))
!       print*,"codbo: ",size(meshe(ndivi)%kfl_codbo(:))
! print*,'meshe(ndivi)%for'

      nullify(nodeFields_export)
      call memory_alloca(memor_dom,'nodeFields_export','mod_AMR',nodeFields_export,INT(size(veloc,2),ip),2_ip)
      !nodeFields_export(:,1) = sqrt(SUM(veloc(1:ndime,:,1)**2,1))
      !nodeFields_export(:,2) = kfl_codno(1,:)!meshe(ndivi)%lnset !press(:,1)!1.0_rp!*mesh_complete%kfl_codno(:)!meshe(ndivi)%kfl_codno !kfl_codbo
      nodeFields_export(:,1) = sqrt(SUM(veloc(1:ndime,:,1)**2,1))
      nodeFields_export(:,2) = press(:,1)!kfl_codno(1,:)!meshe(ndivi)%lnset !press(:,1)!1.0_rp!*mesh_complete%kfl_codno(:)!meshe(ndivi)%kfl_codno !kfl_codbo
      call out_paraview_inp(mesh_new,&
        filename='./'//integer_to_string(kfl_paral)//"meshSolAd_"//integer_to_string(iplot_adapt_sol),nodeFields=nodeFields_export)
      call memory_deallo(memor_dom,'nodeFields_export','mod_AMR',nodeFields_export)
      !
      call mesh_new      % deallo()
    end if

  end subroutine AMR

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-02-29
  !> @brief   Set tags
  !> @details Set tags
  !>          when using copy mesh option, tags are used to
  !>          renumber the mesh right after the adaptivity step. In the other
  !>          case, tags will be used to control the interface motion
  !>          TAG = 1 ... interface nodes frozen
  !>              = 2 ... MESH SIZE
  !>              = 3 ... LNINV_LOC
  !>              = 4 ... LEINV_LOC
  !> 
  !-----------------------------------------------------------------------

  subroutine AMR_set_tags(mesh_new,hh_opt)
    
    type(mesh_type_basic), intent(inout) :: mesh_new
    real(rp), pointer,     intent(in)    :: hh_opt(:)
    integer(ip)                          :: ielem,ipoin,ii
    
    if(kfl_amr_remesh == AMR_ALYA) then
      
       mesh_new % ntags = 2
       call mesh_new % alloca_tags()
       mesh_new % tags(1) % type = NODE_TAG
       mesh_new % tags(2) % type = NODE_TAG  ! the mesh sizes are a nodal field
       call mesh_new % alloca_tag()
       do ii = 1,mesh_new % comm % bound_dim
          ipoin = mesh_new % comm % bound_perm(ii)
          mesh_new % tags(1) % values(ipoin) = 1.0_rp
       end do
       do ipoin = 1,mesh_new % npoin
          mesh_new % tags(2) % values(ipoin) = hh_opt(ipoin)
       end do
       
    else
      
      if( kfl_amr_remesh == AMR_COPY_MESH ) then
         mesh_new % ntags = 4
         call mesh_new % alloca_tags()
         mesh_new % tags(1) % type = NODE_TAG
         mesh_new % tags(2) % type = ELEMENT_TAG
         mesh_new % tags(3) % type = NODE_TAG
         mesh_new % tags(4) % type = ELEMENT_TAG
         call mesh_new % alloca_tag()
         do ipoin = 1,mesh_new % npoin
            mesh_new % tags(3) % values(ipoin) = real(lninv_loc(ipoin),rp)
         end do
         do ielem = 1,mesh_new % nelem
            mesh_new % tags(4) % values(ielem) = real(leinv_loc(ielem),rp)
         end do
      else       
         mesh_new % ntags = 2
         call mesh_new % alloca_tags()
         mesh_new % tags(1) % type = NODE_TAG
         mesh_new % tags(2) % type = ELEMENT_TAG
         call mesh_new % alloca_tag()
      end if
      do ii = 1,mesh_new % comm % bound_dim
         ipoin = mesh_new % comm % bound_perm(ii)
         mesh_new % tags(1) % values(ipoin) = 1.0_rp
      end do
      do ielem = 1,mesh_new % nelem
         mesh_new % tags(2) % values(ielem) = hh_opt(ielem)
      end do
      
    end if
    
  end subroutine AMR_set_tags
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-02-29
  !> @brief   Global numbering
  !> @details Define a global numbering
  !> 
  !-----------------------------------------------------------------------

  subroutine AMR_global_numbering(mesh_new)
    
    type(mesh_type_basic), intent(inout) :: mesh_new
    integer(ip)                          :: ielem,ipoin
    
    if( kfl_amr_remesh == AMR_COPY_MESH ) then
       do ipoin = 1 , mesh_new % npoin 
          mesh_new % lninv_loc(ipoin) = int(mesh_new % tags(3) % values(ipoin),ip)
       end do
       do ielem = 1 , mesh_new % nelem
          mesh_new % leinv_loc(ielem) = int(mesh_new % tags(4) % values(ielem),ip)
       end do
    else
       call mesh_type_basic_global_numbering(mesh_new) 
    end if

  end subroutine AMR_global_numbering
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-02-29
  !> @brief   Determine own nodes
  !> @details Determine own nodes and renumber mesh
  !> 
  !-----------------------------------------------------------------------

  subroutine AMR_own_nodes()
    
    integer(ip), pointer :: permr(:)
    
    nullify(permr)
    call par_tools_ownership_interfaces(npoin,commd,permr)
    call renumbering_node_arrays(permr)
    !
    ! Own nodes
    !
    npoi1     = commd % npoi1
    npoi2     = commd % npoi2
    npoi3     = commd % npoi3
    npoin_own = commd % npoi3
    call memory_deallo(par_memor,'PERMR',vacal,permr)

  end subroutine AMR_own_nodes

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-02-29
  !> @brief   Repartition mesh
  !> @details Move interface elements to create a new partition
  !> 
  !-----------------------------------------------------------------------

  subroutine AMR_repartitioning(mesh_new,rank_nelem,ipass)

    use mod_maths,    only :  maths_sfc_1d_to2d3d_tab
    use def_master,   only :  npart, kfl_paral

    type(mesh_type_basic),    intent(inout) :: mesh_new
    integer(ip),     pointer, intent(inout) :: rank_nelem(:)
    integer(ip),              intent(in)    :: ipass

    integer(ip),     pointer                :: rank_npoin(:)
    integer(ip),     pointer                :: rank_npoin_new(:)
    integer(ip),     pointer                :: ord(:)
    integer(ip),     pointer                :: inv(:)
    integer(ip),     pointer                :: npoin_gat(:)

    integer(ip)                             :: ipoin,ielem,inode
    integer(ip)                             :: min_rank,kpass
    integer(ip)                             :: x,y,z, root
    integer(ip)                             :: ipart, sfc_dim
    integer(ip)                             :: iaux
    real(rp)                                :: raux
    integer(2)                              :: typ, typ_aux
    integer(ip)                             :: npoin_own, ii

    call messages_live('MOVING INTERFACES')

    if( IPARALL ) then
       !
       ! Copy communicator 
       !
       nullify(rank_npoin)
       nullify(rank_npoin_new)
       nullify(ord)
       nullify(inv)
       nullify(npoin_gat)
       
       call memory_alloca(memor_dom,'RANK_NPOIN',    vacal,rank_npoin    ,mesh_new % npoin)
       call memory_alloca(memor_dom,'RANK_NPOIN_NEW',vacal,rank_npoin_new,mesh_new % npoin)
       call memory_alloca(memor_dom,'RANK_NELEM',    vacal,rank_nelem    ,mesh_new % nelem)
       !
       ! Define SFC orientation and reordering arrays
       !
       if(kfl_amr_repart == 0) then 

          iaux    = mod(7_ip*ipass+num_passes*11_ip,24_ip)
          typ     = int(iaux,kind=2)
          typ_aux = typ
          raux    = real(npart,rp)
          root    = 2**(ceiling(log(raux)/log(8.0)))
          sfc_dim = root**3_ip
          call memory_alloca(memor_dom,'ORD', vacal,ord       ,sfc_dim)
          call memory_alloca(memor_dom,'INV', vacal,inv       ,sfc_dim)    
          do ipart = 1, sfc_dim
             call  maths_sfc_1d_to2d3d_tab(root,ipart,x,y,z,typ)
             typ = typ_aux
             ord(ipart)      = (x-1)*root*root+(y-1)*root+z
             inv(ord(ipart)) = ipart 
          enddo

       else if(kfl_amr_repart == 1) then

          npoin_own = mesh_new % npoin
          call memory_alloca(par_memor,'NPOIN_GAT',vacal,npoin_gat,npart+1_ip,LBOUN=0_ip)
          call memory_alloca(par_memor,'ORD'      ,vacal,ord      ,npart+1_ip,LBOUN=0_ip)
          call memory_alloca(par_memor,'INV'      ,vacal,inv      ,npart+1_ip,LBOUN=0_ip)
          call PAR_ALLGATHER(npoin_own,npoin_gat,1_4)

          do ii=1,npart
             inv(ii)=ii
          enddo
          call heapsorti2(2_ip,npart,npoin_gat(1:),inv(1:))
          do ii=1,npart
             ord(inv(ii))=ii
          enddo

       else

          call runend('PAR_REAPRO: UNKNOWN REPART METHOD')

       end if
       !
       ! Evaluate rank_npoin depending on the SFC ordering
       !
       !if( kfl_amr_remesh == 0 ) then
       !   do ielem = 1,mesh_new % nelem
       !      rank_nelem(ielem) = kfl_paral
       !   end do
       !else
       do ipoin = 1,mesh_new % npoin
          rank_npoin(ipoin)     = ord(kfl_paral)
          rank_npoin_new(ipoin) = huge(1_ip)
       end do
       !
       ! Move interface 
       !
       do kpass = 1,2
          call PAR_INTERFACE_NODE_EXCHANGE(rank_npoin,'MIN',COMM = mesh_new % comm)
          do ielem = 1,mesh_new % nelem
             min_rank = huge(1_ip)
             do inode = 1,nnode(mesh_new % ltype(ielem))
                ipoin    = mesh_new % lnods(inode,ielem)
                min_rank = min(min_rank,rank_npoin(ipoin))
             end do
             do inode = 1,nnode(mesh_new % ltype(ielem))
                ipoin                 = mesh_new % lnods(inode,ielem)
                rank_npoin_new(ipoin) = min(rank_npoin_new(ipoin),min_rank)
             end do
             rank_nelem(ielem) = min_rank
          end do
          do ipoin = 1,mesh_new % npoin
             rank_npoin(ipoin) = rank_npoin_new(ipoin)
          end do
       end do
       !
       ! Revert the SFC reordering
       !
       do ielem = 1,mesh_new % nelem
          rank_nelem(ielem) = inv(rank_nelem(ielem))
       enddo

       !end if

       call memory_deallo(memor_dom,'RANK_NPOIN'    ,vacal,rank_npoin)
       call memory_deallo(memor_dom,'RANK_NPOIN_NEW',vacal,rank_npoin_new)
       call memory_deallo(memor_dom,'ORD'           ,vacal,ord)
       call memory_deallo(memor_dom,'INV'           ,vacal,inv)
       call memory_deallo(par_memor,'NPOIN_GAT'     ,vacal,npoin_gat)

    else

       call memory_alloca(memor_dom,'RANK_NELEM',    vacal,rank_nelem    ,mesh_new % nelem)       
       do ielem = 1,mesh_new % nelem
          rank_nelem(ielem) = 0 !int(mesh_new % comm % RANK4,ip)
       end do
       
    end if

  end subroutine AMR_repartitioning

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-02
  !> @brief   Adapt mesh
  !> @details Main Adaptation subroutine
  !> 
  !-----------------------------------------------------------------------

  subroutine AMR_mesh_size(mesh,mesh_size,mesh_bak)

    type(mesh_type_basic),           intent(in)    :: mesh
    real(rp),              pointer,  intent(inout) :: mesh_size(:)
    type(mesh_type_basic), optional, intent(in)    :: mesh_bak
    integer(ip)                                    :: ielem,ipoin,idime,pnode
    integer(ip)                                    :: inode,pelty
    real(rp)                                       :: xx(3),r

    if( present(mesh_bak) ) then

       nullify(mesh_size)
       allocate(mesh_size(mesh_bak % nelem))
       do ielem = 1,mesh_bak % nelem
          pelty = mesh_bak % ltype(ielem)
          pnode = element_type(pelty) % number_nodes       
          xx = 0.0_rp
          do inode = 1,pnode
             ipoin = mesh_bak % lnods(inode,ielem)
             xx(1:ndime) = xx(1:ndime) + mesh_bak % coord(1:ndime,ipoin)
          end do
          xx = xx / real(pnode,rp)
          mesh_size(ielem) = 0.0_rp
          do idime = 1,ndime
             mesh_size(ielem) = mesh_size(ielem) + sqrt((xx(idime) - 0.5_rp)*(xx(idime) - 0.5_rp))
          end do

          r                = sqrt( (xx(1)-0.5_rp)**2 + (xx(2)-0.5_rp)**2 )
          mesh_size(ielem) = 0.1_rp-0.09_rp*exp(-100.0_rp * (r-0.3_rp)**2)
          !mesh_size(ielem) = 0.05_rp * sqrt(mesh_size(ielem))+0.005_rp
          !mesh_size(ielem) = 0.05_rp
       end do
    else

       nullify(mesh_size)
       allocate(mesh_size(mesh % nelem))
       do ielem = 1,mesh % nelem
          pelty = mesh % ltype(ielem)
          pnode = element_type(pelty) % number_nodes       
          xx = 0.0_rp
          do inode = 1,pnode
             ipoin = mesh % lnods(inode,ielem)
             xx(1:ndime) = xx(1:ndime) + mesh % coord(1:ndime,ipoin)
          end do
          xx = xx / real(pnode,rp)
          mesh_size(ielem) = 0.0_rp
          do idime = 1,ndime
             mesh_size(ielem) = mesh_size(ielem) + sqrt((xx(idime) - 0.5_rp)*(xx(idime) - 0.5_rp))
          end do

          r                = sqrt( (xx(1)-0.5_rp)**2 + (xx(2)-0.5_rp)**2 )
          mesh_size(ielem) = 0.1_rp-0.09_rp*exp(-100.0_rp * (r-0.3_rp)**2)
          !mesh_size(ielem) = 0.05_rp * sqrt(mesh_size(ielem))+0.005_rp
          !mesh_size(ielem) = 0.05_rp
       end do
    end if

  end subroutine AMR_mesh_size

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-02
  !> @brief   Adapt mesh
  !> @details Main Adaptation subroutine
  !>
  !>          MESH_CPY: Original mesh
  !>          MESH_EXT: Extracted mesh with mask
  !>          MESH_NOT: Extracted not valid mesh
  !>          MESH_OFF: Mesh not remeshed (other than triangle or tetrahedra)
  !>          MESH_BOU: Boundary of MESH_EXT
  !>          MESH_NEW: New mesh
  !>
  !>          +--------------------------------+
  !>          |                                |
  !>          |                                |
  !>          |            MESH_CPY            |
  !>          |                                |
  !>          |                                |
  !>          +--------------------------------+
  !>                          ||
  !>                          || Extract non-tringale elements from mesh 
  !>                          \/
  !>          +---------------------+----------+
  !>          |                     |          |
  !>          |                     |          |
  !>          |      MESH_EXT       | MESH_OFF |
  !>          |                     |          |
  !>          |                     |          |
  !>          +---------------------+----------+
  !>                     ||
  !>                     || Extract valid mesh to get manifold mesh (take biggest cluster)
  !>                     \/
  !>          +----------+----------+
  !>          |          |          |
  !>          |          |          |
  !>          | MESH_EXT | MESH_NOT |
  !>          |          |          |
  !>          |          |          |
  !>          +----------+----------+
  !>               ||
  !>               || Extract and append boundary mesh
  !>               \/
  !>          +----------+----------+
  !>          |          |          |
  !>          |          |          |
  !>          | MESH_EXT | MESH_BOU |
  !>          |          |          |
  !>          |          |          |
  !>          +----------+----------+
  !>                    ||
  !>                    ||
  !>                    || REMESHING
  !>                    ||
  !>                    \/
  !>               +----------+
  !>               |          |
  !>               |          |
  !>               | MESH_NEW |
  !>               |          |
  !>               |          |
  !>               +----------+
  !>                    ||
  !>                    || Reconstruct global mesh
  !>                    ||
  !>                    \/
  !>          +----------+----------+----------+
  !>          |          |          |          |
  !>          |          |          |          |
  !>          | MESH_NEW | MESH_OFF | MESH_NOT |
  !>          |          |          |          |
  !>          |          |          |          |
  !>          +----------+----------+----------+
  !>                          ||
  !>                          || Append to new mesh
  !>                          ||
  !>                          \/
  !>          +--------------------------------+
  !>          |                                |
  !>          |                                |
  !>          |            MESH_NEW            |
  !>          |                                |
  !>          |                                |
  !>          +--------------------------------+
  !>
  !> 
  !-----------------------------------------------------------------------

  subroutine AMR_adapt_mesh(mesh_new,mesh_cpy)

    type(mesh_type_basic),            intent(inout) :: mesh_new
    type(mesh_type_basic),            intent(inout) :: mesh_cpy
    type(mesh_type_basic)                           :: mesh_ext
    type(mesh_type_basic)                           :: mesh_off
    type(mesh_type_basic)                           :: mesh_bou
    type(mesh_type_basic)                           :: mesh_not
    integer(ip)                                     :: ipoin
    integer(ip)                                     :: ii,itag,jpoin
    integer(ip),            pointer                 :: invpn(:)
    real(rp),               pointer                 :: mesh_size(:)
    logical(lg),            pointer                 :: lmask_nelem(:)
    integer(ip),            pointer                 :: res(:)
    integer(ip),            pointer                 :: permu(:)
    type(maths_bin)                                 :: bin
    !
    ! Initialize meshes
    !    
    nullify(lmask_nelem,invpn,mesh_size,permu,res)

    call mesh_cpy % init('MESH_COPY')
    call mesh_ext % init('MESH_EXTRACTED')
    call mesh_off % init('MESH_MASK')
    call mesh_bou % init('MESH_BOUNDARY')
    call mesh_not % init('MESH_NON_VALID')

    call messages_live('REMESHING')

    if( mesh_new % nelem > 0 ) then
       !
       ! Copy mesh
       !
       call mesh_cpy % copy      (mesh_new)
       call mesh_cpy % copy_comm (mesh_new)
       call mesh_new % deallo    ()
       call mesh_new % init      ('MESH_NEW')
       !
       ! Mask elements that are triangles or tetrahedra
       !
       call memory_alloca(memor_dom,'LMASK_NELEM',vacal,lmask_nelem,mesh_cpy % nelem)
       where( element_type(abs(mesh_cpy % ltype(:))) % topology == 1 ) 
          lmask_nelem = .true.
       elsewhere
          lmask_nelem = .false.
       end where
       !
       ! Extract mesh: MESH_EXT = MESH_CPY - MESH_OFF
       !
       call mesh_ext % extract(mesh_cpy,lmask_nelem,mesh_cmp=mesh_off)
       !
       ! Create a valid mesh: MESH_EXT = MESH_EXT - MESH_NOT
       !
       call mesh_type_basic_valid_mesh(mesh_ext,mesh_not)
       !
       ! Mesh size 
       !
       mesh_size => mesh_cpy % tags(2) % values
       !
       ! Extract boundary mesh MESH_BOU of mesh MESH_EXT 
       !
       call mesh_bou % boundary_mesh(mesh_ext) 
       !
       ! Mesh size
       !
       if( kfl_size_amr == 0 .and. associated(mesh_size) ) then
          min_size_amr = minval(mesh_size)
          max_size_amr = maxval(mesh_size)
       end if
       !
       ! Remesh. The mesh size MESH_SIZE is given on original mesh MESH_CPY 
       !       
       if( kfl_amr_remesh == AMR_GMSH ) then
          if( mesh_bou % npoin > 0 ) call alya2gmsh_remeshing(mesh_new,mesh_cpy,mesh_bou,mesh_size,min_size_amr,max_size_amr)
          if( mesh_new % npoin > 0 ) mesh_new % lninv_loc = 0
          if( mesh_not % npoin > 0 ) mesh_not % lninv_loc = 0
       else if( kfl_amr_remesh == AMR_COPY_MESH ) then
          mesh_new = mesh_ext
        else
          call runend('Check AMR strategy: should bot been here (AMR_GMSH or AMR_COPY)')
       end if
       !
       ! Append extracted meshes and copy comm
       !
       call mesh_off % append    (mesh_not)
       call mesh_new % append    (mesh_off)
       call mesh_new % collapse  ()   
       call mesh_new % copy_comm (mesh_cpy)
       !
       ! Tags. When using the copy mesh option, tags are already copied
       !
       if( kfl_amr_remesh /= 0 ) then
          mesh_new % ntags = mesh_cpy % ntags
          call mesh_new % alloca_tags()
          do itag = 1,mesh_new % ntags
             mesh_new % tags(itag) % type = mesh_cpy % tags(itag) % type
             call mesh_new % alloca_tag(itag)
          end do
       else
          do ii = 1,mesh_new % comm % bound_dim
             ipoin = mesh_new % comm % bound_perm(ii)
             if( int(mesh_new % tags(1) % values(ipoin),ip) == 1 ) mesh_new % tags(1) % values(ipoin) = 2.0_rp
          end do
          do ipoin = 1,mesh_new % npoin
             mesh_new % tags(1) % values(ipoin) = max(mesh_new % tags(1) % values(ipoin)-1.0_rp,0.0_rp)
          end do
       end if
       !
       ! Recovering original numbering in permutation arrays
       !
       call memory_alloca(memor_dom,'PERMU',vacal,PERMU,mesh_cpy % npoin)       
       call bin % init  ()
       call bin % input (LIMIT=20_ip)
       call bin % fill  (mesh_new % coord)       
       do ii = 1,mesh_cpy % comm % bound_dim
          ipoin =  mesh_cpy % comm % bound_perm(ii)
          res   => mesh_new % identical_coord(mesh_cpy % coord(:,ipoin),bin)
          if( associated(res) ) then
             permu(ipoin) = res(1)
             deallocate(res)
          end if
       end do
       do ii = 1,mesh_cpy % comm % bound_dim
          ipoin =  mesh_cpy % comm % bound_perm(ii)
          if( permu(ipoin) == 0 ) then
             write(6,*)'a=',kfl_paral,mesh_cpy % coord(:,ipoin)
             call runend('MOD_AMR: PERMUTATION MOT FOUND')
          end if
       end do
       call bin      % deallo()
       call mesh_new % comm % permute(permu)
       do ii = 1,mesh_cpy % comm % bound_dim
          ipoin = mesh_cpy % comm % bound_perm(ii)
          jpoin = permu(ipoin)
          mesh_new % tags(1) % values(jpoin) = mesh_cpy % tags(1) % values(ipoin)
       end do
       call memory_deallo(memor_dom,'PERMU',vacal,PERMU)
       !block
       !  call mesh_new % output  (filename='mesh-before-'//integer_to_string(kfl_paral))
       !  call mesh_new % results (xx=mesh_new % tags(1) % values,names='TAGSS',filename='mesh-before-'//integer_to_string(kfl_paral),where='ON NODES')
       !end block
    end if
    !
    ! Deallocate temporary meshes
    !
    ii = mesh_new % nelem
    call PAR_SUM(ii)
    call messages_live('   NUMBER OF ELEMENTS= '//integer_to_string(ii))

    !call mesh_cpy % deallo()
    call mesh_bou % deallo()
    call mesh_not % deallo()
    call mesh_off % deallo()
    call mesh_ext % deallo()
    call memory_deallo(memor_dom,'LMASK_NELEM',vacal,lmask_nelem)

  end subroutine AMR_adapt_mesh
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  abel gargallo (from houzeaux AMR_adapt_mesh)
  !> @date    2020-06-02
  !> @brief   Adapt mesh
  !> @details Main Adaptation subroutine
  !>
  !>          MESH_CPY: Original mesh
  !>          MESH_EXT: Extracted mesh with mask
  !>          MESH_NOT: Extracted not valid mesh
  !>          MESH_OFF: Mesh not remeshed (other than triangle or tetrahedra)
  !>          MESH_BOU: Boundary of MESH_EXT
  !>          MESH_NEW: New mesh
  !>
  !>          +--------------------------------+
  !>          |                                |
  !>          |                                |
  !>          |            MESH_CPY            |
  !>          |                                |
  !>          |                                |
  !>          +--------------------------------+
  !>                    ||
  !>                    || Extract interface from mesh 
  !>                    \/
  !>          +---------------------+----------+
  !>          |                     |          |
  !>          |                     |          |
  !>          |      MESH_EXT       | MESH_OFF |
  !>          |                     |          |
  !>          |                     |          |
  !>          +---------------------+----------+
  !>                    ||
  !>                    ||
  !>                    || REMESHING
  !>                    ||
  !>                    \/
  !>               +----------+
  !>               |          |
  !>               |          |
  !>               | MESH_NEW |
  !>               |          |
  !>               |          |
  !>               +----------+
  !>                    ||
  !>                    || Reconstruct global mesh
  !>                    ||
  !>                    \/
  !>          +----------+----------+
  !>          |          |          |
  !>          |          |          |
  !>          | MESH_NEW | MESH_OFF |
  !>          |          |          |
  !>          |          |          |
  !>          +----------+----------+
  !>                     ||
  !>                     || Append to new mesh
  !>                     ||
  !>                     \/
  !>          +--------------------------------+
  !>          |                                |
  !>          |                                |
  !>          |            MESH_NEW            |
  !>          |                                |
  !>          |                                |
  !>          +--------------------------------+
  !>
  !> 
  !-----------------------------------------------------------------------

  subroutine AMR_adapt_mesh_ALYA(mesh_new,mesh_cpy,is_modified)
    
    use mod_adapt, only: adapt_mesh_to_metric
    use mod_metric, only: mesh_metric_type
    use mod_metric, only: compute_metric_fromElemSizeField
    use mod_metric, only: compute_metric_fromNodalSizeField
    use mod_communications_global, only: PAR_OR

    type(mesh_type_basic),            intent(inout) :: mesh_new
    type(mesh_type_basic),            intent(inout) :: mesh_cpy
    logical(lg),            optional, intent(inout) :: is_modified
    type(mesh_type_basic)                           :: mesh_ext
    type(mesh_type_basic)                           :: mesh_off
!    type(mesh_type_basic)                           :: mesh_not
    integer(ip)                                     :: ipoin
    integer(ip)                                     :: ii,itag,jpoin
    integer(ip),            pointer                 :: invpn(:)
    real(rp),               pointer                 :: mesh_size(:)
    logical(lg),            pointer                 :: lmask_nelem(:)
    integer(ip),            pointer                 :: res(:)
    integer(ip),            pointer                 :: permu(:)
    type(maths_bin)                                 :: bin
    type(mesh_metric_type)                          :: metric
    !
    logical(lg) :: compute_permu_adapti = .false. !permu_adapti is not always equal to permu below... some error? TODO and then use .true.
    integer(ip), pointer :: permu_adapti(:)
    
    logical(lg) :: isModifiedMesh 
    real(rp) :: t_perf0,t_perf1
    logical(lg), parameter :: perform_debug = .false.
    !
    if(perform_debug) call cputim(t_perf0)
    
    !
    ! Initialize meshes
    !    
    nullify(lmask_nelem,invpn,mesh_size,permu,res)

    call mesh_cpy % init('MESH_COPY')
    call mesh_ext % init('MESH_EXTRACTED')
    call mesh_off % init('MESH_MASK')
    !call mesh_not % init('MESH_NON_VALID')

    call messages_live('REMESHING')

    isModifiedMesh = .false.
    
    if( mesh_new % nelem > 0 ) then
       !
       ! Copy mesh
       !
       call mesh_cpy % copy      (mesh_new)
       call mesh_cpy % copy_comm (mesh_new)
       call mesh_new % deallo    ()
       call mesh_new % init      ('MESH_NEW')
       !
       ! Mask elements that are triangles or tetrahedra
       !
       call memory_alloca(memor_dom,'LMASK_NELEM',vacal,lmask_nelem,mesh_cpy % nelem)
       where( element_type(abs(mesh_cpy % ltype(:))) % topology == 1 ) 
          lmask_nelem = .true.
       elsewhere
          lmask_nelem = .false.
       end where
       !
       ! Extract mesh: MESH_EXT = MESH_CPY - MESH_OFF
       !
       call mesh_ext % extract(mesh_cpy,lmask_nelem,mesh_cmp=mesh_off)
       !
       ! Create a valid mesh: MESH_EXT = MESH_EXT - MESH_NOT
       !
       !call mesh_type_basic_valid_mesh(mesh_ext,mesh_not) !ABEL
       !
       ! Mesh size 
       !
       mesh_size => mesh_ext % tags(2) % values
       !
       ! Mesh size
       !
       if( kfl_size_amr == 0 .and. associated(mesh_size) ) then
          min_size_amr = minval(mesh_size)
          max_size_amr = maxval(mesh_size)
       end if
       !
       ! Remesh. The mesh size MESH_SIZE is given on original mesh MESH_CPY 
       !
       call compute_metric_fromNodalSizeField(mesh_size,mesh_ext,metric)
       
       if(perform_debug) then
         call cputim(t_perf1)
         print*,'   AMR_ADAPT up to adaptMesh: ',t_perf1-t_perf0
       end if
       if(perform_debug) call cputim(t_perf0)
       if(compute_permu_adapti) then
         call adapt_mesh_to_metric(&
          mesh_ext,&
          metric,&
          lock_valid_elems =.true.,&
          isModifiedMesh   =isModifiedMesh,&
          thresholdQuality =AMR_qualityTarget,&
          activationQuality=AMR_qualityActivation,&
          mapNodes_input_to_adapted_mesh=permu_adapti)
       else
         call adapt_mesh_to_metric(&
          mesh_ext,&
          metric,&
          lock_valid_elems =.true.,&
          isModifiedMesh   =isModifiedMesh,&
          thresholdQuality =AMR_qualityTarget,&
          activationQuality=AMR_qualityActivation)
       end if
       if(perform_debug) then
         call cputim(t_perf1)
         print*,'   AMR_ADAPT cost adaptMesh: ',t_perf1-t_perf0
       end if
       if(perform_debug) call cputim(t_perf0)
       !call metric%deallo()
       
       !
       ! Do the following when mesh is not modified.. right now commented.. should be checked
!        if(.not.isModifiedMesh) then
!          call mesh_new % copy      (mesh_cpy)
!          call mesh_new % copy_comm (mesh_cpy)
!          return
!        end if
       !
       call mesh_new % copy      (mesh_ext)
       !
       ! Append extracted meshes and copy comm
       !
       !call mesh_off % append    (mesh_not)
       call mesh_new % append    (mesh_off)
       call mesh_new % collapse  ()   
       call mesh_new % copy_comm (mesh_cpy)
       !
       ! Tags. 
       !
       mesh_new % ntags = mesh_cpy % ntags
       call mesh_new % alloca_tags()
       do itag = 1,mesh_new % ntags
          mesh_new % tags(itag) % type = mesh_cpy % tags(itag) % type
          call mesh_new % alloca_tag(itag)
       end do
       ! Lo seguent no puc pel collapse, si tingues la permutacio si que podria
       !mesh_new % tags(2) % values(1:mesh_new%npoin) = mesh_cpy % tags(2) % values(1:mesh_new%npoin) 
       !ABEL: I did this:
!        if( kfl_amr_remesh /= 0 ) then
          ! mesh_new % ntags = mesh_cpy % ntags
!           call mesh_new % alloca_tags()
!           do itag = 1,mesh_new % ntags
!              mesh_new % tags(itag) % type = mesh_cpy % tags(itag) % type
!              call mesh_new % alloca_tag(itag)
!           end do
!        else
!           do ii = 1,mesh_new % comm % bound_dim
!              ipoin = mesh_new % comm % bound_perm(ii)
!              if( int(mesh_new % tags(1) % values(ipoin),ip) == 1 ) mesh_new % tags(1) % values(ipoin) = 2.0_rp
!           end do
!           do ipoin = 1,mesh_new % npoin
!              mesh_new % tags(1) % values(ipoin) = max(mesh_new % tags(1) % values(ipoin)-1.0_rp,0.0_rp)
!           end do
!        end if
       !
       ! Recovering original numbering in permutation arrays
       ! ====> using this to compute size field from adapted mesh (to avoid interpolation)
       !
       if(.not.compute_permu_adapti) then
         call memory_alloca(memor_dom,'PERMU',vacal,PERMU,mesh_cpy % npoin)
         call bin % init  ()
         call bin % input (LIMIT=20_ip)
         call bin % fill  (mesh_new % coord)
         do ii = 1,mesh_cpy % comm % bound_dim
            ipoin =  mesh_cpy % comm % bound_perm(ii)
            res   => mesh_new % identical_coord(mesh_cpy % coord(:,ipoin),bin)
            if( associated(res) ) then
               permu(ipoin) = res(1)
               !write(6,*)permu(ipoin), " vs ",permu_adapti(ipoin)
               deallocate(res)
            end if
         end do
         ! 
         ! To compute mesh size
         !
         if(.not.metric%isSizeField) then
           call runend('Extend mesh tags to work with metrics...')
         end if
         mesh_new % tags(2) % values = huge(1.0_rp)
         do ipoin = 1,mesh_ext % npoin
           res   => mesh_new % identical_coord(mesh_ext % coord(:,ipoin),bin)
           if( associated(res) ) then
             ii = res(1)
             if( ii == 0 ) then
                write(6,*)'a=',kfl_paral,mesh_ext % coord(:,ipoin)
                call runend('MOD_AMR: point in adapted mesh not found')
             end if
             deallocate(res)
             mesh_new % tags(2) % values(ii) = metric%M(1_ip,1_ip,ipoin)
           else
             write(6,*)'a=',kfl_paral,mesh_ext % coord(:,ipoin)
             call runend('MOD_AMR: point in adapted mesh not found')
           end if
         end do
       
         call metric%deallo()
         !
       end if
       
       do ii = 1,mesh_cpy % comm % bound_dim
          ipoin =  mesh_cpy % comm % bound_perm(ii)
          if( permu(ipoin) == 0 ) then
             write(6,*)'a=',kfl_paral,mesh_cpy % coord(:,ipoin)
             call runend('MOD_AMR: PERMUTATION MOT FOUND')
          end if
       end do
       if(.not.compute_permu_adapti) then
         call bin      % deallo()
       end if
       call mesh_new % comm % permute(permu)
       do ii = 1,mesh_cpy % comm % bound_dim
          ipoin = mesh_cpy % comm % bound_perm(ii)
          jpoin = permu(ipoin)
          mesh_new % tags(1) % values(jpoin) = mesh_cpy % tags(1) % values(ipoin)
       end do
       call memory_deallo(memor_dom,'PERMU',vacal,PERMU)
       !block
       !  call mesh_new % output  (filename='mesh-before-'//integer_to_string(kfl_paral))
       !  call mesh_new % results (xx=mesh_new % tags(1) % values,names='TAGSS',filename='mesh-before-'//integer_to_string(kfl_paral),where='ON NODES')
       !end block
    end if
    !
    ! Deallocate temporary meshes
    !
    call PAR_OR(isModifiedMesh)
    
    if(isModifiedMesh) then
      call messages_live('ACTIVATED -> MIN_Q < Q_ACTI' )
      !ii = mesh_new % nelem
      !call PAR_SUM(ii)
      !call messages_live('   NUMBER OF ELEMENTS= '//integer_to_string(ii))
      ii = mesh_new % npoin
      call PAR_SUM(ii)
      call messages_live('   NUMBER OF NODES= '//integer_to_string(ii))
    else
      call messages_live('NOT ACTIVATED -> MIN_Q > Q_ACTI' )
    end if
    
    if(present(is_modified)) then      
      is_modified = isModifiedMesh
    end if

    !call mesh_cpy % deallo()
    !call mesh_bou % deallo() !ABEL
    !call mesh_not % deallo() !ABEL
    call mesh_off % deallo()
    call mesh_ext % deallo()
    call memory_deallo(memor_dom,'LMASK_NELEM',vacal,lmask_nelem)
    
    if(perform_debug) then
      call cputim(t_perf1)
      print*,'   AMR_ADAPT from adaptMesh: ',t_perf1-t_perf0
    end if
    
  end subroutine AMR_adapt_mesh_ALYA

  !-----------------------------------------------------------------------
  !> 
  !> @author  borrell
  !> @date    2021-10-06
  !> @brief   Interface nodes
  !> @details Determine interface nodes
  !> 
  !-----------------------------------------------------------------------
  
  subroutine PAR_communicator_post(mesh,lnod_owners,mowners)

     type(mesh_type_basic),        intent(inout) :: mesh
     integer(ip),          pointer,intent(inout) :: lnod_owners(:,:)
     integer(ip),                  intent(inout) :: mowners
     integer(ip)                                 :: bound_dim,ii,ipoin
     integer(ip)                                 :: my_rank,ipart,ipoin_global
     integer(ip)                                 :: my_size
     integer(ip)                                 :: jj,kk,nrows
     integer(ip)                                 :: lid
     type(hash_t)                                :: ht     
     integer(ip),          pointer               :: list_ranks(:)
     type(i1p),            pointer               :: list_nodes(:)
     
     my_rank = int(mesh % comm % RANK4,ip)
     my_size = int(mesh % comm % SIZE4,ip)

     nullify(list_ranks)
     nullify(list_nodes)
     call memory_alloca(memor_dom,'LIST_RANKS',vacal,list_ranks,my_size)
     call memory_alloca(memor_dom,'LIST_NODES',vacal,list_nodes,my_size)

     call ht % init()
     call htaini(ht,mowners*10_ip,lidson=.true.,REPEATED_ELEMENTS=.true.)
     
     do ipoin = 1,memory_size(lnod_owners,2_ip)
        ipoin_global = lnod_owners(0,ipoin)
        if( lnod_owners(2,ipoin) > 0 ) then
           loop_ii: do ii = 1,mowners
              ipart = lnod_owners(ii,ipoin)
              if(ipart > 0 .and. ipart /= my_rank ) then
                 call htaadd(ht,ipart,lid=lid,icont=ipoin_global)
                 list_ranks(lid) = ipart
              else if ( ipart < 0 ) then
                 exit loop_ii
              end if
           end do loop_ii
        end if
     end do
     !
     ! Order neighbors
     !
     if( ht % nelem > 0 ) call maths_heap_sort(2_ip,ht % nelem,list_ranks)
     bound_dim = 0
     do ii = 1,ht % nelem
        ipart = list_ranks(ii)
        list_nodes(ii) % l => htalid_list(ht,ipart)
        if( associated(list_nodes(ii) % l) ) then
           nrows     = memory_size(list_nodes(ii) % l)
           bound_dim = bound_dim + nrows
           call maths_heap_sort(2_ip,nrows,list_nodes(ii) % l)
        end if
     end do
     call htades(ht)
     !
     ! Allocate new communicator
     !
     call mesh % comm % deallo(INITIALIZE=.false.)
     mesh % comm % nneig     = ht % nelem
     mesh % comm % bound_dim = bound_dim
     call mesh % comm % alloca()
     !
     ! Reconstruct communicator
     !     
     call ht % init()
     call htaini(ht,mesh % npoin*10_ip,lidson=.true.)
     call htaadd(ht,mesh % lninv_loc)
     jj = 0
     mesh % comm % bound_size(1) = 1
     do ii = 1,mesh % comm % nneig
        nrows                          = memory_size(list_nodes(ii) % l)
        mesh % comm % neights(ii)      = list_ranks(ii)
        mesh % comm % bound_size(ii+1) = mesh % comm % bound_size(ii) + nrows
        do kk = 1,nrows
           jj           = jj + 1
           ipoin_global = list_nodes(ii) % l(kk)
           ipoin        = htalid(ht,ipoin_global)
           mesh % comm % bound_perm(jj) = ipoin
        end do
     end do
     !
     ! Deallocate
     !
     call htades(ht)
     call memory_deallo(memor_dom,'LIST_NODES' ,vacal,list_nodes)
     call memory_deallo(memor_dom,'LIST_RANKS' ,vacal,list_ranks)
     call memory_deallo(par_memor,'LNOD_OWNERS',vacal,lnod_owners)

  end subroutine PAR_communicator_post
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  borrell
  !> @date    2021-10-06
  !> @brief   Interface nodes
  !> @details Determine interface nodes
  !> 
  !-----------------------------------------------------------------------
  
  subroutine PAR_communicator_pre(mesh,lrank_elem,lnod_owners,mowners)

    type(mesh_type_basic),           intent(in)    :: mesh
    integer(ip),            pointer, intent(in)    :: lrank_elem(:)
    integer(ip),            pointer, intent(inout) :: lnod_owners(:,:)
    integer(ip),                     intent(inout) :: mowners

    integer(ip)                                    :: ielem, ipoin,  nrank
    integer(ip)                                    :: irank, ii, jj, nsend, nrecv
    integer(ip)                                    :: ninter, stride, ibuff, offset,nsize
    integer(ip)                                    :: lid, nnse, nnre, ipart,ielpo,ineig
    integer(ip)                                    :: dummi(1), mepoi_loc,naux,kk,nn
    type(hash_t)                                   :: ht     
    integer(ip),            pointer                :: pelpo_loc(:), lelpo_loc(:)
    integer(ip),            pointer                :: lnown(:),lown(:,:)
    integer(ip),            pointer                :: lnsend(:),lnrecv(:)
    integer(ip),            pointer                :: bsend(:), brecv(:)
    integer(ip),            pointer                :: loffse(:),loffre(:)
    integer(ip),            pointer                :: contre(:)

    logical(lg)                                    :: isin
    MY_MPI_COMM                                    :: PAR_COMM4

    if( IPARALL ) then

       nullify(pelpo_loc,lelpo_loc)
       nullify(lnown,lown)
       nullify(lnsend,lnrecv)
       nullify(loffse,loffre)
       nullify(bsend,brecv)
       nullify(contre)
       nullify(lnod_owners)

       PAR_COMM4 = mesh % comm % PAR_COMM_WORLD
       nsize     = int(mesh % comm % SIZE4,ip)
       !
       ! How many ranks exist in my array lrank_elem array
       ! Multiply by two just in case some of my neighbors have extra ranks
       !
       call htaini(ht,nsize*10_ip)
       do ielem = 1,mesh % nelem
          call htaadd(ht,lrank_elem(ielem))
       end do
       nrank = ht % nelem
       call PAR_MAX(nrank)
       nrank = 2 * nrank
       call htades(ht)
       !
       ! Generate the graph of neighbour elements per point
       !
       call graphs_elepoi(mesh,mepoi_loc,pelpo_loc,lelpo_loc,memor=memor_dom)
       !
       ! Evaluate lnown: number of owners per node 
       !          lown : list owners per node
       !
       call memory_alloca(memor_dom,'LNOWN',vacal,lnown,mesh % npoin)
       call memory_alloca(memor_dom,'LOWN' ,vacal,lown , nrank, mesh % npoin)
       call htaini(ht,nrank*10,lidson = .true.)

       do ipoin = 1,mesh % npoin
          call htares(ht)
          do ielpo = pelpo_loc(ipoin),pelpo_loc(ipoin+1)-1
             ielem = lelpo_loc(ielpo)
             call htaadd(ht,lrank_elem(ielem),lid)
             lown(lid,ipoin) = lrank_elem(ielem)
          enddo
          lnown(ipoin) = ht % nelem
       end do
       call htades(ht)

       naux = 0
       do ii = 1,mesh  % comm % bound_dim
          ipoin = mesh % comm % bound_perm(ii)
          naux  = max(naux,lnown(ipoin))
       end do
       call PAR_MAX(naux)
       !
       ! Add ranks to which interface nodes belong (on the other side of the interface)
       !
       call memory_alloca(memor_dom,'BSEND' ,vacal,bsend ,mesh % comm % bound_dim * naux,INIT_VALUE=-1_ip)
       call memory_alloca(memor_dom,'BRECV' ,vacal,brecv ,mesh % comm % bound_dim * naux)
       call memory_alloca(memor_dom,'LOFFSE',vacal,loffse,mesh % comm % nneig)

       nsend = 0
       do ineig = 1,mesh % comm % nneig
          loffse(ineig) = nsend * naux  
          nsend         = nsend + mesh % comm % bound_size(ineig+1)-mesh % comm % bound_size(ineig)
       end do

       do ineig = 1,mesh % comm % nneig
          do ii = mesh % comm % bound_size(ineig),mesh % comm % bound_size(ineig+1)-1
             ipoin = mesh % comm % bound_perm(ii)
             do kk = 1,lnown(ipoin)
                bsend(loffse(ineig)+kk) = lown(kk,ipoin)
             end do
             loffse(ineig) = loffse(ineig) + naux
          end do
       end do

       nsend = 0
       do ineig = 1,mesh % comm % nneig
          loffse(ineig) = nsend * naux  
          nsend         = nsend + mesh % comm % bound_size(ineig+1)-mesh % comm % bound_size(ineig)
       end do

       do ineig = 1,mesh % comm % nneig
          ipart = mesh % comm % neights(ineig)
          nsend = (mesh % comm % bound_size(ineig+1)-mesh % comm % bound_size(ineig)) * naux
          ii    = loffse(ineig) + 1
          call PAR_SEND_RECEIVE(nsend,nsend,bsend(ii:ii+nsend-1),brecv(ii:ii+nsend-1),'IN THE WORLD',ipart)
       end do

       do ineig = 1,mesh % comm % nneig
          do ii = mesh % comm % bound_size(ineig),mesh % comm % bound_size(ineig+1)-1
             ipoin = mesh % comm % bound_perm(ii)
             loop_kk: do kk = 1,naux
                irank = brecv(loffse(ineig)+kk)
                if( irank >= 0 ) then
                   isin = .false.
                   do lid = 1,lnown(ipoin)
                      if( irank == lown(lid,ipoin) ) isin = .true.
                   end do
                   if( .not. isin ) then
                      lnown(ipoin)             = lnown(ipoin) + 1
                      lown(lnown(ipoin),ipoin) = irank
                   end if
                else
                   exit loop_kk
                end if
             end do loop_kk
             loffse(ineig) = loffse(ineig) + naux
          end do
       end do

       call memory_deallo(memor_dom,'BSEND' ,vacal,bsend )
       call memory_deallo(memor_dom,'BRECV' ,vacal,brecv )
       call memory_deallo(memor_dom,'LOFFSE',vacal,loffse)
       !
       ! Evaluate lnsend: number of points to send per rank
       !          lnrecv: number of points to recv per rank
       !          loffse: offsets in the send buffer
       !          loffre: offsets in teh recv buffer
       !           nsend: total number of elements to send
       !           nrecv: total number of elements to recv
       !           nrank: max number of ranks per point (gloval)
       !
       call memory_alloca(memor_dom,'LNSEND',vacal,lnsend,nsize+1,LBOUN=0_ip)
       call memory_alloca(memor_dom,'LNRECV',vacal,lnrecv,nsize+1,LBOUN=0_ip)

       nrank = 0
       do ipoin = 1, mesh % npoin
          nrank = max(nrank,lnown(ipoin))
          if( lnown(ipoin) > 1 ) then
             do ii = 1,lnown(ipoin)
                irank         = lown(ii,ipoin)
                lnsend(irank) = lnsend(irank)+1_ip
             end do
          end if
       end do

       call PAR_MAX(nrank)
       call PAR_ALLTOALL(1_ip,1_ip,lnsend,lnrecv,PAR_COMM_IN=PAR_COMM4)

       call memory_alloca(memor_dom,'LOFFSE',vacal,loffse,nsize+2_ip,LBOUN=0_ip)
       call memory_alloca(memor_dom,'LOFFRE',vacal,loffre,nsize+2_ip,LBOUN=0_ip)

       nn    = nrank+1
       nrecv = 0
       nsend = 0
       do irank = 0,nsize
          loffse(irank) = nsend * nn  !note: for each point are received nrank+1 integers
          loffre(irank) = nrecv * nn
          nrecv         = nrecv + lnrecv(irank)
          nsend         = nsend + lnsend(irank)
       end do
       loffse(nsize+1) = nsend * nn
       loffre(nsize+1) = nrecv * nn
       !
       ! bsend and brecv are the buffers to perform communications
       !
       call memory_alloca(memor_dom,'BSEND',vacal,bsend,nsend * nn,INIT_VALUE=-1_ip)
       call memory_alloca(memor_dom,'BRECV',vacal,brecv,nrecv * nn)
       !
       ! Fill buffer send
       !
       do ipoin = 1, mesh % npoin
          if( lnown(ipoin) > 1 ) then
             do ii = 1, lnown(ipoin)
                irank = lown(ii,ipoin)
                bsend(loffse(irank)+1) = mesh % lninv_loc(ipoin)
                do jj=2,lnown(ipoin)+1
                   bsend(loffse(irank)+jj) = lown(jj-1,ipoin)
                enddo
                loffse(irank) = loffse(irank)+nn !the slot of one point used
             end do
          end if
       end do
       !
       ! Reevaluate loffse that has been changed
       !
       nsend = 0
       do irank = 0, nsize
          loffse(irank) = nsend * nn
          nsend         = nsend + lnsend(irank)
       enddo
       loffse(nsize+1) = nsend * nn

       !write(6,*)'recv=',kfl_paral,loffre(:)
       !write(6,*)'send=',kfl_paral,loffse(:)
       !write(6,*)'recv=',kfl_paral,lnsend(1:)
       !write(6,*)'send=',kfl_paral,lnrecv(1:)
       !
       ! Send and recv
       !
       call PAR_START_NON_BLOCKING_COMM(1_ip,nsize+1_ip)
       call PAR_SET_NON_BLOCKING_COMM_NUMBER(1_ip)
       do ipart = 1,nsize
          nnse = lnsend(ipart) * nn
          nnre = lnrecv(ipart) * nn
          if( nnse > 0 .and. nnre > 0 ) then
             call PAR_SEND_RECEIVE(nnse,nnre,bsend(loffse(ipart)+1:loffse(ipart+1)),&
                  brecv(loffre(ipart)+1:loffre(ipart+1)),'IN THE WORLD',ipart,'NON BLOCKING')
          else
             if(nnse>0) then
                call PAR_SEND_RECEIVE(nnse,0_ip,bsend(loffse(ipart)+1:loffse(ipart+1)),dummi(1:1),&
                     &'IN THE WORLD',ipart,'NON BLOCKING')
             endif
             if(nnre>0) then
                call PAR_SEND_RECEIVE(0_ip,nnre,dummi(1:1),brecv(loffre(ipart)+1:loffre(ipart+1)),&
                     &'IN THE WORLD',ipart,'NON BLOCKING')
             end if
          end if
       end do
       call PAR_END_NON_BLOCKING_COMM(1_ip)
       !call PAR_BARRIER()

       !if( inotmaster ) then
       !   write(6,*)'a=',kfl_paral,': ',bsend
       !   write(6,*)'b=',kfl_paral,': ',brecv
       !end if
       !
       ! Evaluate htnodes
       !
       call memory_alloca(memor_dom,'CONTRE',vacal,contre,max(nrecv,1_ip),'INITIALIZE')
       stride = nrank+1
       call htaini(ht,10*nrecv, lidson=.true.)
       mowners=0
       do ibuff = 0, nrecv-1
          call htaadd(ht,brecv(ibuff*stride+1),lid)
          contre(lid) = contre(lid)+1
          mowners     = max(mowners,contre(lid))
       enddo
       ninter = ht % nelem
       !
       ! Evaluate lnod_owners
       !
       call memory_alloca(memor_dom,'LNOD_OWNERS',vacal,lnod_owners,mowners*nrank+1_ip,ninter,INIT_VALUE=-1_ip,LBOUN1=0_ip,LBOUN2=1_ip)

       if( associated(contre) ) contre(:)=0
       mowners=0
       do ibuff = 0, nrecv-1
          offset = ibuff*stride
          lid = htalid(ht, brecv(offset+1))

          !if( brecv(offset+1)==26 )  write(6,*)'WWWWWW=',kfl_paral,': ',brecv(offset+1:(ibuff+1)*stride-1)
          do ii = 1, nrank
             if(brecv(offset+1+ii) /= -1_ip) then
                isin = .false.
                do jj=1,contre(lid) 
                   if(brecv(offset+1+ii)==lnod_owners(jj,lid)) isin = .true.
                enddo
                if(.not. isin) then
                   if( contre(lid) == 0 ) then
                      lnod_owners(contre(lid),lid)=brecv(offset+1)
                   end if
                   contre(lid) = contre(lid) + 1_ip
                   if(contre(lid) > mowners ) mowners = contre(lid)
                   lnod_owners(contre(lid),lid)=brecv(offset+1+ii)
                end if
             end if
          end do
       end do
       call htades(ht)

       !if( kfl_paral==1) then
       !   do ii = 1,ninter
       !      if( lnod_owners(2,ii)/=-1 ) write(6,*)'ll=',lnod_owners(0,ii),': ',lnod_owners(1:,ii)
       !   end do
       !end if
       !call runend('O.K.!')
       !
       ! Deallocate pointers
       !
       call memory_deallo(memor_dom,'LNOWN'  ,vacal,lnown)
       call memory_deallo(memor_dom,'LOWN'   ,vacal,lown )
       call memory_deallo(memor_dom,'LNSEND' ,vacal,lnsend)
       call memory_deallo(memor_dom,'LNRECV' ,vacal,lnrecv)
       call memory_deallo(memor_dom,'LOFFSE' ,vacal,loffse)
       call memory_deallo(memor_dom,'LOFFRE' ,vacal,loffre)
       call memory_deallo(memor_dom,'BSEND'  ,vacal,bsend)
       call memory_deallo(memor_dom,'BRECV'  ,vacal,brecv)
       call memory_deallo(memor_dom,'CONTRE' ,vacal,contre)
       call graphs_elepoi_deallocate(pelpo_loc,lelpo_loc,memor=memor_dom)

    end if
    
  end subroutine PAR_communicator_pre
   
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-02-29
  !> @brief   Interface nodes
  !> @details Determine interface nodes
  !> 
  !-----------------------------------------------------------------------

  subroutine AMR_new_mesh(rank_nelem,mesh_new)

    integer(ip),           pointer, intent(in)    :: rank_nelem(:)
    type(mesh_type_basic),          intent(inout) :: mesh_new
    integer(ip)                                   :: ielem,ipoin,nsize
    integer(ip)                                   :: rank,ipart,ii
    integer(ip)                                   :: mepoi_loc,ielpo,itag
    integer(ip)                                   :: mnode_max
    MY_MPI_COMM                                   :: PAR_COMM4
    integer(4)                                    :: SIZE4
    integer(4)                                    :: RANK4
    logical(lg),           pointer                :: list_neighbors(:)
    integer(ip),           pointer                :: nelem_send(:)
    integer(ip),           pointer                :: npoin_send(:)
    integer(ip),           pointer                :: npoin_recv(:)
    integer(ip),           pointer                :: nelem_recv(:)
    type(i1p),             pointer                :: inte_npoin(:)
    type(mesh_type_basic), pointer                :: geome_send(:)
    type(mesh_type_basic), pointer                :: geome_recv(:)
    logical(lg),           pointer                :: lmask(:)
    integer(ip),           pointer                :: lelpo_loc(:)
    integer(ip),           pointer                :: pelpo_loc(:)
    integer(ip),           pointer                :: rank_npoin(:)

    if( ISEQUEN ) return
    
    call messages_live('MIGRATING MESH')
    
    nullify(inte_npoin)
    nullify(list_neighbors)
    nullify(nelem_send)
    nullify(nelem_recv)
    nullify(npoin_send)
    nullify(npoin_recv)
    nullify(geome_send)
    nullify(geome_recv)
    nullify(lmask)
    nullify(pelpo_loc)
    nullify(lelpo_loc)
    nullify(rank_npoin)
    !
    ! Here we assume that rank 0 cannot receive elements
    !
    PAR_COMM4 = mesh_new % comm % PAR_COMM_WORLD
    SIZE4     = mesh_new % comm % SIZE4
    RANK4     = mesh_new % comm % RANK4
    nsize     = int(SIZE4,ip)-1_ip

    call memory_alloca(memor_dom,'NELEM_SEND'    ,vacal,nelem_send    ,nsize+1_ip,LBOUN=0_ip)
    call memory_alloca(memor_dom,'NPOIN_SEND'    ,vacal,npoin_send    ,nsize+1_ip,LBOUN=0_ip)
    call memory_alloca(memor_dom,'NELEM_RECV'    ,vacal,nelem_recv    ,nsize+1_ip,LBOUN=0_ip)
    call memory_alloca(memor_dom,'NPOIN_RECV'    ,vacal,npoin_recv    ,nsize+1_ip,LBOUN=0_ip)
    call memory_alloca(memor_dom,'LIST_NEIGHBORS',vacal,list_neighbors,nsize)

    !--------------------------------------------------------------------
    !
    ! Determine migrating nodes INTE_NPOIN(IPOIN) % L(:) = RANK
    !
    !--------------------------------------------------------------------

    call graphs_elepoi(mesh_new,mepoi_loc,pelpo_loc,lelpo_loc,memor=memor_dom)
    call memory_alloca(memor_dom,'RANK_NPOIN',vacal,rank_npoin,mepoi_loc)
    call memory_alloca(memor_dom,'INTE_NPOIN',vacal,inte_npoin,mesh_new % npoin)

    do ipoin = 1,mesh_new % npoin
       ipart      = 0
       rank_npoin = 0
       do ielpo = pelpo_loc(ipoin),pelpo_loc(ipoin+1)-1
          ielem = lelpo_loc(ielpo)
          rank  = rank_nelem(ielem)
          if( count( rank_npoin == rank ) == 0 ) then
             ipart                = ipart + 1
             rank_npoin(ipart)    = rank
             list_neighbors(rank) = .true.
          end if
       end do
       call memory_alloca(memor_dom,'INTE_NPOIN % L',vacal,inte_npoin(ipoin) % l,ipart)
       inte_npoin(ipoin) % l(1:ipart) = rank_npoin(1:ipart)
    end do
    call memory_deallo(memor_dom,'RANK_NPOIN',vacal,rank_npoin)
    call graphs_elepoi_deallocate(pelpo_loc,lelpo_loc,memor=memor_dom)
    
    !-----------------------------------------------------------------
    !
    ! Send buffer sizes
    !    
    !-----------------------------------------------------------------
    
    do ipoin = 1,mesh_new % npoin
       do ii = 1,memory_size(inte_npoin(ipoin) % l)
          ipart = inte_npoin(ipoin) % l(ii)
          npoin_send(ipart) = npoin_send(ipart) + 1
       end do
    end do
    do ielem = 1,mesh_new % nelem
       ipart = rank_nelem(ielem)
       nelem_send(ipart) = nelem_send(ipart) + 1
    end do    
    call PAR_ALLTOALL(1_ip,1_ip,npoin_send,npoin_recv,PAR_COMM_IN=PAR_COMM4)
    call PAR_ALLTOALL(1_ip,1_ip,nelem_send,nelem_recv,PAR_COMM_IN=PAR_COMM4)
    
    !-----------------------------------------------------------------
    !
    ! Extract submesh GEOME_SEND from MESH_NEW to be sent
    !
    !-----------------------------------------------------------------
    
    allocate(geome_send(nsize))
    allocate(geome_recv(nsize))
    do ipart = 1,nsize
       call geome_send(ipart) % init()
       call geome_recv(ipart) % init()
    end do

    mnode_max = mesh_new % mnode
    call PAR_MAX(mnode_max)
    
    geome_send % mnode = mnode_max
    geome_send % ndime = mesh_new % ndime
    geome_recv % mnode = mnode_max
    geome_recv % ndime = mesh_new % ndime

    call memory_alloca(memor_dom,'LMASK',vacal,lmask,mesh_new % nelem)
    do ipart = 1,nsize              
       if( nelem_send(ipart) > 0 ) then
          where( rank_nelem == ipart ) lmask = .true.
          call geome_send(ipart) % extract(mesh_new,lmask)          
          lmask = .false.
       end if
    end do
    call memory_deallo(memor_dom,'LMASK',vacal,lmask)

    !-----------------------------------------------------------------
    !
    ! Send/Recv submeshes GEOME_SEND/GEOME_RECV
    !    
    !-----------------------------------------------------------------
    !
    ! If we have tags
    !
    do ipart = 1,nsize
       geome_recv(ipart) % ntags = mesh_new % ntags      
    end do
    !
    ! Allocate basic arrays
    !
    do ipart = 1,nsize
       call geome_recv(ipart) % alloca(NELEM=nelem_recv(ipart),NPOIN=npoin_recv(ipart))
       geome_recv(ipart) % nelem = nelem_recv(ipart)
       geome_recv(ipart) % npoin = npoin_recv(ipart)
    end do
    !
    ! Assume tags exist everywhere
    !
    if( mesh_new % ntags > 0 ) then
       do ipart = 1,nsize
          call geome_recv(ipart) % alloca(NTAGS=mesh_new % ntags)
          do itag = 1,mesh_new % ntags
             geome_recv(ipart) % tags(itag) % type = mesh_new % tags(itag) % type
             call geome_recv(ipart) % alloca_tag(itag)
          end do
       end do
    end if
    !
    ! Although there is no tags, the tag type should be defined to enable the send/recv
    !
    do ipart = 1,nsize
       if( .not. associated(geome_send(ipart) % tags) ) call geome_send(ipart) % alloca(NTAGS=mesh_new % ntags)
    end do
    !
    ! Send receive mesh and tags
    !
    do ipart = 1,nsize
       call mesh_type_basic_send_recv(geome_send(ipart),geome_recv(ipart),ipart,PAR_COMM4)
    end do

    !-----------------------------------------------------------------
    !
    ! Create new mesh MESH_NEW from received meshes GEOME_RECV
    !    
    !-----------------------------------------------------------------
    
    call mesh_new % deallo()
    call mesh_new % init  ('MESH_NEW')    

    do ipart = 1,nsize
       call mesh_new % merge(geome_recv(ipart))
    end do

    !   if( inotmaster ) then
    !      call mesh_new % output (filename='tags-'//integer_to_string(kfl_paral))
    !      call mesh_new % results(xx=mesh_new % tags(1) % values,names='TAGSS',filename='tags-'//integer_to_string(kfl_paral),where='ON NODES')
    !   end if

    !--------------------------------------------------------------------
    !
    ! Deallocate memory
    !
    !--------------------------------------------------------------------

    call memory_deallo(memor_dom,'INTE_NPOIN'    ,vacal,inte_npoin)
    call memory_deallo(memor_dom,'LIST_NEIGHBORS',vacal,list_neighbors)
    call memory_deallo(memor_dom,'NELEM_SEND'    ,vacal,nelem_send)
    call memory_deallo(memor_dom,'NPOIN_SEND'    ,vacal,npoin_send)
    call memory_deallo(memor_dom,'NPOIN_RECV'    ,vacal,npoin_recv)
    call memory_deallo(memor_dom,'NELEM_RECV'    ,vacal,nelem_recv)

    do ipart = 1,nsize 
       call geome_send(ipart) % deallo()
       call geome_recv(ipart) % deallo()
    end do
    
    mesh_new % comm % SIZE4          = SIZE4  
    mesh_new % comm % RANK4          = RANK4  
    mesh_new % comm % PAR_COMM_WORLD = PAR_COMM4

  end subroutine AMR_new_mesh
 
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-03-02
  !> @brief   Domain
  !> @details Recompute the domain
  !> 
  !-----------------------------------------------------------------------

  subroutine AMR_update_original_mesh(mesh_complete)
    
    type(mesh_type), intent(inout) :: mesh_complete

    call messages_live('UPDATE DOMAIN')
    !
    ! Copy new mesh
    !
    nelem = mesh_complete % nelem
    npoin = mesh_complete % npoin
    nboun = mesh_complete % nboun
#ifndef PNODE_VALUE
    mnode = mesh_complete % mnode
    call PAR_MAX(mnode)
#endif
    !
    ! Deallocate old mesh
    !
    call memory_deallo(memor_dom,'LNODS'    ,'mesh_type_save_original_mesh',lnods    )
    call memory_deallo(memor_dom,'LTYPE'    ,'mesh_type_save_original_mesh',ltype    )
    call memory_deallo(memor_dom,'LEINV_LOC','mesh_type_save_original_mesh',leinv_loc)

    call memory_deallo(memor_dom,'COORD'    ,'mesh_type_save_original_mesh',coord    )
    call memory_deallo(memor_dom,'LNINV_LOC','mesh_type_save_original_mesh',lninv_loc)
    
    call memory_deallo(memor_dom,'LNODB'    ,'mesh_type_save_original_mesh',lnodb    )
    call memory_deallo(memor_dom,'LTYPB'    ,'mesh_type_save_original_mesh',ltypb    )
    call memory_deallo(memor_dom,'LELBO'    ,'mesh_type_save_original_mesh',lelbo    )
    call memory_deallo(memor_dom,'LBOCH'    ,'mesh_type_save_original_mesh',lboch    )
    call memory_deallo(memor_dom,'LBINV_LOC','mesh_type_save_original_mesh',lbinv_loc)
    !
    ! Copy new mesh from MESH_COMPLETE
    !
    call memory_copy(memor_dom,'MESH_COMPLETE % LNODS'    ,'mesh_type_save_original_mesh',mesh_complete % lnods    ,lnods    ,'DO_NOT_DEALLOCATE',COPY_NAME='LNODS')
    call memory_copy(memor_dom,'MESH_COMPLETE % LTYPE'    ,'mesh_type_save_original_mesh',mesh_complete % ltype    ,ltype    ,'DO_NOT_DEALLOCATE',COPY_NAME='LTYPE')
    call memory_copy(memor_dom,'MESH_COMPLETE % LEINV_LOC','mesh_type_save_original_mesh',mesh_complete % leinv_loc,leinv_loc,'DO_NOT_DEALLOCATE',COPY_NAME='LEINV_LOC')

    call memory_copy(memor_dom,'MESH_COMPLETE % COORD'    ,'mesh_type_save_original_mesh',mesh_complete % coord    ,coord    ,'DO_NOT_DEALLOCATE',COPY_NAME='COORD')
    call memory_copy(memor_dom,'MESH_COMPLETE % LNINV_LOC','mesh_type_save_original_mesh',mesh_complete % lninv_loc,lninv_loc,'DO_NOT_DEALLOCATE',COPY_NAME='LNINV_LOC')

    call memory_copy(memor_dom,'MESH_COMPLETE % LNODB'    ,'mesh_type_save_original_mesh',mesh_complete % lnodb    ,lnodb    ,'DO_NOT_DEALLOCATE',COPY_NAME='LNODB')
    call memory_copy(memor_dom,'MESH_COMPLETE % LTYPB'    ,'mesh_type_save_original_mesh',mesh_complete % ltypb    ,ltypb    ,'DO_NOT_DEALLOCATE',COPY_NAME='LTYPB')
    call memory_copy(memor_dom,'MESH_COMPLETE % LELBO'    ,'mesh_type_save_original_mesh',mesh_complete % lelbo    ,lelbo    ,'DO_NOT_DEALLOCATE',COPY_NAME='LELBO')
    call memory_copy(memor_dom,'MESH_COMPLETE % LBOCH'    ,'mesh_type_save_original_mesh',mesh_complete % lboch    ,lboch    ,'DO_NOT_DEALLOCATE',COPY_NAME='LBOCH')
    call memory_copy(memor_dom,'MESH_COMPLETE % LBINV_LOC','mesh_type_save_original_mesh',mesh_complete % lbinv_loc,lbinv_loc,'DO_NOT_DEALLOCATE',COPY_NAME='LBINV_LOC')
    !
    ! Deallocate others
    !
    call domain_memory_deallocate('LNNOD')
    call domain_memory_deallocate('LNNOB')
    call domain_memory_deallocate('LGAUS')
    call domain_memory_deallocate('LPOTY')
    call domain_memory_deallocate('LNLEV')

  end subroutine AMR_update_original_mesh

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-03-02
  !> @brief   Domain
  !> @details Create boundary mesh
  !> 
  !-----------------------------------------------------------------------

  subroutine AMR_boundary_mesh(mesh_complete)

    type(mesh_type),       intent(inout) :: mesh_complete
    integer(ip)                          :: iboun

    call messages_live('CREATE BOUNDARY MESH')

    if( mesh_complete % nelem > 0 ) then
       !
       ! Creat local boundary meshes
       !
       call meshes_list_boundary_elements(&
            mesh_complete % nelem,mesh_complete % npoin,mesh_complete % lnods,mesh_complete % ltype,&
            mesh_complete % nboun,mesh_complete % lelbo,mesh_complete % lnodb,mesh_complete % ltypb,&
            memor_dom,&
            LELBO_NAME=trim(mesh_complete % name)//' % LELBO',&
            LNODB_NAME=trim(mesh_complete % name)//' % LNODB',&
            LTYPB_NAME=trim(mesh_complete % name)//' % LTYPB')
       !
       ! Remove internal faces
       !
       call AMR_remove_internal_boundaries(mesh_complete)
       !
       ! Compute LBOCH and LBINV_LOC
       !
       call memory_alloca(memor_dom,trim(mesh_complete % name)//' % LBOCH'    ,'mesh_type_save_original_mesh',mesh_complete % lboch    ,mesh_complete % nboun)
       call memory_alloca(memor_dom,trim(mesh_complete % name)//' % LBINV_LOC','mesh_type_save_original_mesh',mesh_complete % lbinv_loc,mesh_complete % nboun)
       do iboun = 1,mesh_complete % nboun
          mesh_complete % lboch(iboun) = BOFEM
       end do
       
    end if

   call par_global_numbering_elements_boundaries(mesh_complete % nboun,mesh_complete % lbinv_loc)

  end subroutine AMR_boundary_mesh

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-03-02
  !> @brief   Communicator
  !> @details Recompute communication arrays
  !> 
  !-----------------------------------------------------------------------

  subroutine AMR_update_communication_arrays(mesh_complete)

    type(mesh_type), intent(inout) :: mesh_complete
    
    call messages_live('UPDATE COMMUNICATION ARRAYS')

    call PAR_COMM_MY_CODE_ARRAY(1) % deallo(par_memor)
    call PAR_COMM_MY_CODE_ARRAY(1) % copy  (mesh_complete % comm,par_memor)
    
    PAR_COMM_MY_CODE_ARRAY(1) % PAR_COMM_WORLD =  PAR_COMM_MY_CODE
    PAR_COMM_MY_CODE_ARRAY(1) % RANK4          =  int(kfl_paral,4)
    PAR_COMM_MY_CODE_ARRAY(1) % SIZE4          =  int(npart+1,4)
    commd                                      => PAR_COMM_MY_CODE_ARRAY(1)
    !
    ! Determine own nodes
    !
    call memory_copy(memor_dom,'COMM % BOUND_PERM',vacal,&
         PAR_COMM_MY_CODE_ARRAY(1) % bound_perm,&
         PAR_COMM_MY_CODE_ARRAY(1) % bound_invp,'DO_NOT_DEALLOCATE',COPY_NAME='COMMD % BOUND_INVP')
    
  end subroutine AMR_update_communication_arrays

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-03-02
  !> @brief   Domain
  !> @details Recompute the domain
  !> 
  !-----------------------------------------------------------------------

  subroutine AMR_domain(mesh_new)

    type(mesh_type),         intent(inout) :: mesh_new 
    integer(ip)                            :: ifiel
    real(rp),    pointer                   :: bobox(:,:,:)
    real(rp),    pointer                   :: bbbox(:,:,:)
    real(rp),    pointer                   :: subox(:,:,:)
    real(rp),    pointer                   :: sbbox(:,:,:)
    real(rp),    pointer                   :: centr(:,:)
    real(rp),    pointer                   :: cbntr(:,:)
    
    integer(ip) :: ii,count_not_found
    
    integer(ip) :: num_to_find
    integer(ip) :: num_found

    call messages_live('ARRAY INTERPOLATIONS')
    !
    ! Nullify pointers
    !
    nullify(bobox)
    nullify(bbbox)
    nullify(subox)
    nullify(sbbox)
    nullify(centr)
    nullify(cbntr)

    !--------------------------------------------------------------------
    !
    ! Initialize search strategies
    !
    !--------------------------------------------------------------------

    call PAR_BARRIER()
    call cputim(time1)
    
    !--------------------------------------------------------------------
    !
    ! Interpolation strategies
    !
    !--------------------------------------------------------------------

    call messages_live('INTERPOLATION STRATEGIES','START SECTION')

    allocate(interp_AMR_npoin)
    allocate(interp_AMR_nelem)
    allocate(interp_AMR_nboun)

    call interp_AMR_npoin % init ()
    call interp_AMR_nelem % init ()
    call interp_AMR_nboun % init ()

    !--------------------------------------------------------------------
    !
    ! Search strategies, BOBOX, BBBOX and SUBOX
    !
    !--------------------------------------------------------------------

    call messages_live('SET SEARCH METHODS')

    call meshe(ndivi) % element_bb  (bobox,MEMORY_COUNTER=memor_dom)
    call meshe(ndivi) % boundary_bb (bbbox,MEMORY_COUNTER=memor_dom)
    
    call par_bounding_box(mesh_new % comm % PAR_COMM_WORLD,bobox,subox,MEMORY_COUNTER=memor_dom)
    call par_bounding_box(mesh_new % comm % PAR_COMM_WORLD,bbbox,sbbox,MEMORY_COUNTER=memor_dom)
    
    call search_vol_seq % set    (BOBOX=bobox,NAME='VOL_SEQ')
    call search_vol_par % set    (BOBOX=subox,NAME='VOL_PAR')
    call search_bou_seq % set    (BOBOX=bbbox,NAME='BOU_SEQ')
    call search_bou_par % set    (BOBOX=sbbox,NAME='BOU_PAR')
    
    call memory_deallo(memor_dom,'BOBOX','par_bounding_box',bobox)
    call memory_deallo(memor_dom,'BBBOX','par_bounding_box',bbbox) 
    call memory_deallo(memor_dom,'SUBOX','par_bounding_box',subox)
    call memory_deallo(memor_dom,'SBBOX','par_bounding_box',sbbox)
    
    !--------------------------------------------------------------------
    !
    ! Create element interpolation. Interpolation on host element.
    !
    !--------------------------------------------------------------------

    call messages_live('SET NODE INTERPOLATION METHOD')
    call interp_AMR_npoin % input      (search_vol_seq % method,search_vol_par % method,&
         &                             INTERPOLATION_METHOD=INT_ELEMENT_INTERPOLATION,&
         &                             COMM=mesh_new % comm % PAR_COMM_WORLD,&
         &                             FORCE_FIND   =.false.,&       ! XXX ABEL: check if this helps always finding nodes in prev mesh
         &                             MAX_IT_FORCE = 6_ip )
    call interp_AMR_npoin % preprocess (mesh_new % coord,meshe(ndivi))
   
    if(kfl_paral==-1_ip) then ! SEQUENTIAL CHECK
      count_not_found = 0_ip
      do ii = 1_ip,size(interp_AMR_npoin % found)
         if( .not. interp_AMR_npoin % found(ii) ) then
            count_not_found = count_not_found + 1_ip
            print*,'point not found= ',ii,lpoty(ii)
         end if
      end do
      if(count_not_found>0_ip) then
        print*,'Total number of points=',size(interp_AMR_npoin % found)
        print*,'Total points not found=',count_not_found
        call runend(' Point not found, cannot continue procedure')
      end if
    else ! PARALLEL CHECK
      
!       print*,'AMR_domain: CHECK IF ALL POINTS FOUND'
      
!       num_to_find =  size(mesh_new % coord,2)
!       call PAR_SUM(num_to_find)
!
!       num_found   = count(interp_AMR_npoin % found)
!       call PAR_SUM(num_found  )
!
!       if(kfl_paral==0) then
!         if(num_to_find==num_found) then
!           print*,'AMR_domain: ALL POINTS FOUND'
!         else if(num_to_find<num_found) then
!           print*,'AMR_domain: some points missing'
!           print*,'  Num points to find: ',num_to_find
!           print*,'  Num points found  : ',num_found
!           call runend('AMR_domain: some points missing')
!         else if(num_to_find<num_found) then
!           print*,'AMR_domain: some points found in more than one processor...!'
!           print*,'  Num points to find: ',num_to_find
!           print*,'  Num points found  : ',num_found
!           call runend('AMR_domain: some points double-found')
!         end if
!       end if
    end if
!      call runend('ENSURE ALL POINTS ARE FOUND: NOW THIS IS NOT WORKING AND HAS A SEGMENTATION FAULT')
        
    !--------------------------------------------------------------------
    !
    ! Create element value. One value per element on centroid.
    !
    !--------------------------------------------------------------------

    call messages_live('SET ELEMENT INTERPOLATION METHOD')
    call mesh_new         % centroid   (centr)
    call interp_AMR_nelem % input      (search_vol_seq % method,search_vol_par % method,&
         &                             INTERPOLATION_METHOD=INT_ELEMENT_VALUE,&
         &                             COMM=mesh_new % comm % PAR_COMM_WORLD)
    call interp_AMR_nelem % preprocess (centr,meshe(ndivi))
    call mesh_new         % centroid   (centr,ONLY_DEALLOCATE=.true.)
   
    if(kfl_paral==-1_ip) then ! SEQUENTIAL CHECK
      count_not_found = 0_ip
      do ii = 1_ip,size(interp_AMR_nelem % found)
         if( .not. interp_AMR_nelem % found(ii) ) then
            count_not_found = count_not_found + 1_ip
            print*,'element not found= ',ii,lpoty(ii)
         end if
      end do
      if(count_not_found>0_ip) then
        print*,'Total number of elements=',size(interp_AMR_nelem % found)
        print*,'Total elements not found=',count_not_found
        call runend(' Element not found, cannot continue procedure')
      end if
    end if

    !--------------------------------------------------------------------
    !
    ! Create boundary value using kdtree
    !
    !--------------------------------------------------------------------
    call messages_live('SET BOUNDARY INTERPOLATION METHOD')

    call mesh_new            % associate_boundary ()
    call mesh_new % boundary % centroid           (cbntr)
    call interp_AMR_nboun    % input              (search_bou_seq % method,search_bou_par % method,&
         &                                        INTERPOLATION_METHOD=INT_BOUNDARY_VALUE,&
         &                                        COMM=mesh_new % comm % PAR_COMM_WORLD)
    call interp_AMR_nboun    % preprocess         (cbntr,meshe(ndivi))
    call mesh_new % boundary % centroid           (cbntr,ONLY_DEALLOCATE=.true.)
    call mesh_new            % associate_boundary (ONLY_DEALLOCATE=.true.)
    
    call messages_live('INTERPOLATION STRATEGIES','END SECTION')
   
    if(kfl_paral==-1_ip) then ! SEQUENTIAL CHECK
      count_not_found = 0_ip
      do ii = 1_ip,size(interp_AMR_nboun % found)
         if( .not. interp_AMR_nboun % found(ii) ) then
            count_not_found = count_not_found + 1_ip
            print*,'bound not found= ',ii,lpoty(ii)
         end if
      end do
      if(count_not_found>0_ip) then
        print*,'Total number of bounds=',size(interp_AMR_nboun % found)
        print*,'Total bounds not found=',count_not_found
        call runend(' Bound not found, cannot continue procedure')
      end if
    end if

    !--------------------------------------------------------------------
    !
    ! Interpolate domain arrays 
    !
    !--------------------------------------------------------------------
 
    call messages_live('INTERPOLATION MESH ARRAYS','START SECTION')

    npoin_new = mesh_new % npoin
    nboun_new = mesh_new % nboun
    nelem_new = mesh_new % nelem
    
    !call AMR_interpolate(lgrou_dom,'NPOIN',1_ip,memor_dom,'LGROU_DOM')
    
    do ifiel = 1,nfiel
       if( kfl_field(2,ifiel) == NPOIN_TYPE ) then
          if( kfl_field(4,ifiel) == 1 ) then
             call runend('AMR: NODAL FIELD WITH AMR NOT CODED')
             call AMR_interpolate(xfiel(ifiel) % a,'NPOIN',2_ip,memor_dom,VARIABLE_NAME='XFIEL % A')
          end if
       end if
    end do

    call AMR_interpolate(lelch,    'NELEM',1_ip,memor_dom,'LELCH')
    call AMR_interpolate(lesub,    'NELEM',1_ip,memor_dom,'LESUB')
    call AMR_interpolate(lmate,    'NELEM',1_ip,memor_dom,'LMATE')
    call AMR_interpolate(leset,    'NELEM',1_ip,memor_dom,'LESET')
    do ifiel = 1,nfiel
       if( kfl_field(2,ifiel) == NELEM_TYPE ) then
          if( kfl_field(4,ifiel) == 1 ) then
             call AMR_interpolate(xfiel(ifiel) % a,'NELEM',2_ip,memor_dom,VARIABLE_NAME='XFIEL % A')
          end if
       end if
    end do

    call AMR_interpolate(kfl_codbo,'NBOUN',2_ip,memor_dom,'KFL_CODBO')    
    call AMR_interpolate(lbset,    'NBOUN',1_ip,memor_dom,'LBSET')

    do ifiel = 1,nfiel
       if( kfl_field(2,ifiel) == NBOUN_TYPE ) then
          if( kfl_field(4,ifiel) == 1 ) then
             call AMR_interpolate(xfiel(ifiel) % a,'NBOUN',2_ip,memor_dom,VARIABLE_NAME='XFIEL % A')
          end if
       end if
    end do

    call messages_live('INTERPOLATION MESH ARRAYS','END SECTION')

    nelem = mesh_new % nelem
    npoin = mesh_new % npoin
    nboun = mesh_new % nboun

    !--------------------------------------------------------------------
    !
    ! Arrays which do not make sense to be interpolated
    !
    !--------------------------------------------------------------------
     
    call domain_memory_reallocate('LGROU_DOM')
    call domain_memory_reallocate('LNOCH')
    call domain_memory_reallocate('LMAST')
    call domain_memory_reallocate('LNSET')
    call domain_memory_reallocate('KFL_CODNO')

    call par_global_variables_arrays()
    call par_mesh_dimensions()

  end subroutine AMR_domain

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-03-09
  !> @brief   Remove internal boundaries
  !> @details Identify common faces with neighbors and remove them
  !> 
  !-----------------------------------------------------------------------

  subroutine AMR_remove_internal_boundaries(mesh)

    type(mesh_type),     intent(inout) :: mesh
    integer(ip),         pointer       :: nboun_send(:)
    type(i2p),           pointer       :: lboun_send(:)
    integer(ip),         pointer       :: nboun_recv(:)
    type(i2p),           pointer       :: lboun_recv(:)
    integer(ip),         pointer       :: nsubd_npoin(:)
    type(i1p),           pointer       :: lsubd_npoin(:)
    integer(ip),         pointer       :: lboun_tot(:,:)
    integer(ip),         pointer       :: knodb(:)
    integer(ip)                        :: ii,ipoin,iboun,inodb,pnodb
    integer(ip)                        :: ineig,dom_i,kneig,kboun
    integer(ip)                        :: lnodb_loc(mnodb)
    MY_MPI_COMM                        :: PAR_COMM_AMR4
    logical(lg)                        :: same_boundary
    integer(ip),         pointer       :: lnodb_cpy(:,:)   
    integer(ip),         pointer       :: ltypb_cpy(:)    
    integer(ip),         pointer       :: lelbo_cpy(:)  

    if( mesh % comm % nneig > 0 ) then

       nullify(nboun_send)
       nullify(lboun_send)
       nullify(nboun_recv)
       nullify(lboun_recv)
       nullify(nsubd_npoin)
       nullify(lsubd_npoin)
       nullify(lboun_tot)
       nullify(knodb)

       nullify(lnodb_cpy)
       nullify(ltypb_cpy)
       nullify(lelbo_cpy)

       PAR_COMM_AMR4 = mesh % comm % PAR_COMM_WORLD

       call memory_alloca(memor_dom,'NBOUN_SEND' ,vacal,nboun_send ,mesh % comm % nneig)
       call memory_alloca(memor_dom,'LBOUN_SEND' ,vacal,lboun_send ,mesh % comm % nneig)
       call memory_alloca(memor_dom,'LBOUN_TOT'  ,vacal,lboun_tot  ,mesh % comm % nneig,mesh % nboun)
       call memory_alloca(memor_dom,'NSUBD_NPOIN',vacal,nsubd_npoin,mesh % npoin)
       call memory_alloca(memor_dom,'LSUBD_NPOIN',vacal,lsubd_npoin,mesh % npoin)
       do ipoin = 1,mesh % npoin
          call memory_alloca(memor_dom,'LSUBD_NPOIN % L',vacal,lsubd_npoin(ipoin)%l,mesh % comm % nneig)
       end do
       call memory_alloca(memor_dom,'KNODB' ,vacal,knodb,mesh % comm % nneig)
       !
       ! List of node neighbors LSUBD_NPOIN(1:NPOIN) % L(:) 
       !
       do ineig = 1,mesh % comm % nneig
          do ii = mesh % comm % bound_size(ineig),mesh % comm % bound_size(ineig+1)-1
             ipoin                                      = mesh % comm % bound_perm(ii)
             nsubd_npoin(ipoin)                         = nsubd_npoin(ipoin) + 1
             lsubd_npoin(ipoin) % l(nsubd_npoin(ipoin)) = ineig
          end do
       end do
       !
       ! NBOUN_SEND(INEIG): Number of boundaries to send
       !
       do iboun = 1,mesh % nboun
          knodb = 0
          pnodb = nnode(abs(mesh % ltypb(iboun)))
          do inodb = 1,pnodb
             ipoin = mesh % lnodb(inodb,iboun)
             do ii = 1,nsubd_npoin(ipoin)
                ineig        = lsubd_npoin(ipoin) % l(ii)
                knodb(ineig) = knodb(ineig) + 1
             end do
          end do
          kneig = 0
          do ineig = 1,mesh % comm % nneig
             if( knodb(ineig) == pnodb ) then
                kneig                  = kneig + 1
                nboun_send(ineig)      = nboun_send(ineig) + 1
                lboun_tot(kneig,iboun) = ineig
             end if
          end do
       end do
       call memory_deallo(memor_dom,'NSUBD_NPOIN',vacal,nsubd_npoin)
       call memory_deallo(memor_dom,'LSUBD_NPOIN',vacal,lsubd_npoin)
       !
       ! LBOUN_SEND(INEIG) % L(:,:): List of boundaries to send
       !
       do ineig = 1,mesh % comm % nneig
          call memory_alloca(memor_dom,'LBOUN_SEND % L',vacal,lboun_send(ineig) % l,mesh % mnodb,nboun_send(ineig))
          nboun_send(ineig) = 0
       end do

       do iboun = 1,mesh % nboun
          loop_kneig: do kneig = 1,mesh % comm % nneig
             ineig = lboun_tot(kneig,iboun)
             if( ineig == 0 ) then
                exit loop_kneig
             else
                pnodb             = nnode(abs(mesh % ltypb(iboun)))
                nboun_send(ineig) = nboun_send(ineig) + 1
                do inodb = 1,pnodb
                   ipoin = mesh % lnodb(inodb,iboun)
                   lboun_send(ineig) % l(inodb,nboun_send(ineig)) = mesh % lninv_loc(ipoin)
                end do
                call maths_heap_sort(2_ip,pnodb,lboun_send(ineig) % l(:,nboun_send(ineig)))
             end if
          end do loop_kneig
       end do
       !
       ! Send/Receive list
       !
       call memory_alloca(memor_dom,'NBOUN_RECV' ,vacal,nboun_recv ,mesh % comm % nneig)
       call memory_alloca(memor_dom,'LBOUN_RECV' ,vacal,lboun_recv ,mesh % comm % nneig)

       do ineig = 1,mesh % comm % nneig
          dom_i = mesh % comm % neights(ineig)
          call PAR_SEND_RECEIVE(nboun_send(ineig),nboun_recv(ineig),'IN MY CODE',dom_i,PAR_COMM_IN=PAR_COMM_AMR4)
          call memory_alloca(memor_dom,'LBOUN_RECV % L',vacal,lboun_recv(ineig) % l,mesh % mnodb,nboun_recv(ineig))
          call PAR_SEND_RECEIVE(lboun_send(ineig) %l,lboun_recv(ineig) %l,'IN MY CODE',dom_i,PAR_COMM_IN=PAR_COMM_AMR4)
       end do
       !
       ! Identify internal boundaries
       !
       do iboun = 1,mesh % nboun
          pnodb              = nnode(abs(mesh % ltypb(iboun)))
          lnodb_loc(1:pnodb) = mesh % lninv_loc(mesh % lnodb(1:pnodb,iboun))
          lboun_tot(1,iboun) = 1
          call maths_heap_sort(2_ip,pnodb,lnodb_loc)
          loop_ineig: do ineig = 1,mesh % comm % nneig
             do kboun = 1,nboun_recv(ineig)
                same_boundary = .true.
                loop_inodb: do inodb = 1,pnodb
                   if( lnodb_loc(inodb) /= lboun_recv(ineig) %l(inodb,kboun) ) then
                      same_boundary = .false.
                      exit loop_inodb
                   end if
                end do loop_inodb
                if( same_boundary ) then
                   lboun_tot(1,iboun) = 0
                   exit loop_ineig
                end if
             end do
          end do loop_ineig
       end do
       !
       ! Deallocate
       !
       call memory_deallo(memor_dom,'LBOUN_RECV',vacal,lboun_recv )
       call memory_deallo(memor_dom,'NBOUN_RECV',vacal,nboun_recv )
       call memory_deallo(memor_dom,'LBOUN_SEND',vacal,lboun_send )
       call memory_deallo(memor_dom,'NBOUN_SEND',vacal,nboun_send )
       call memory_deallo(memor_dom,'KNODB'     ,vacal,knodb      )
       !
       ! Resize boundary arrays
       !
       kboun = 0
       do iboun = 1,mesh % nboun
          kboun = kboun + lboun_tot(1,iboun)       
       end do
       call memory_copy  (memor_dom,trim(mesh % name)//' % LNODB',vacal,mesh % lnodb,lnodb_cpy,COPY_NAME='LNODB_CPY')
       call memory_copy  (memor_dom,trim(mesh % name)//' % LTYPB',vacal,mesh % ltypb,ltypb_cpy,COPY_NAME='LTYPB_CPY')
       call memory_copy  (memor_dom,trim(mesh % name)//' % LELBO',vacal,mesh % lelbo,lelbo_cpy,COPY_NAME='LELBO_CPY')
       call memory_alloca(memor_dom,trim(mesh % name)//' % LNODB',vacal,mesh % lnodb,mnodb,kboun,'DO_NOT_INITIALIZE')
       call memory_alloca(memor_dom,trim(mesh % name)//' % LTYPB',vacal,mesh % ltypb,kboun      ,'DO_NOT_INITIALIZE')
       call memory_alloca(memor_dom,trim(mesh % name)//' % LELBO',vacal,mesh % lelbo,kboun      ,'DO_NOT_INITIALIZE')

       kboun = 0
       do iboun = 1,mesh % nboun
          if( lboun_tot(1,iboun) == 1 ) then
             kboun                           = kboun + 1
             pnodb                           = nnode(abs(ltypb_cpy(iboun)))
             mesh % lnodb(1:pnodb,kboun) = lnodb_cpy(1:pnodb,iboun)
             mesh % ltypb(kboun)         = ltypb_cpy(iboun)
             mesh % lelbo(kboun)         = lelbo_cpy(iboun)
          end if
       end do
       mesh % nboun = kboun

       call memory_deallo(memor_dom,'LBOUN_TOT' ,vacal,lboun_tot)
       call memory_deallo(memor_dom,'LNODB_CPY' ,vacal,lnodb_cpy)
       call memory_deallo(memor_dom,'LTYPB_CPY' ,vacal,ltypb_cpy)
       call memory_deallo(memor_dom,'LELBO_CPY' ,vacal,lelbo_cpy)

    end if

  end subroutine AMR_remove_internal_boundaries

  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-10-06
  !> @brief   Interpolate mesh size
  !> @details Interpolate mesh size form mesh_in to mesh_out
  !> 
  !-----------------------------------------------------------------------
  
  subroutine AMR_interpolate_size(mesh_out,mesh_in)

    type(mesh_type_basic), intent(inout) :: mesh_in
    type(mesh_type_basic), intent(inout) :: mesh_out
    real(rp),              pointer       :: bobox(:,:,:)
    real(rp),              pointer       :: centroid(:,:)
    class(interpolation),  pointer       :: interp
!    real(rp)                             :: xx(3)

    call messages_live('INTERPOLATE MESH SIZE')

    !write(6,*)'A=',kfl_paral,mesh_out % nelem
    
    if( associated(mesh_in % tags) .and. mesh_out % nelem > 0 ) then
       
       nullify(bobox)
       nullify(centroid)
       nullify(search_amr)
       nullify(interp)
       !
       ! Select search method
       !
       select case ( kfl_mesh_size_amr )
       case        ( SEARCH_BIN        ) ; allocate( maths_bin    :: search_amr)
       case        ( SEARCH_OCTREE     ) ; allocate( maths_octree :: search_amr)
       case default                      ; call runend('MOD_AMR: WRONG SEARCH STRATEGY')
       end select
       call search_amr % init ()
       !
       ! Set search strategy
       !
       select type ( search_amr   )
       class is    ( maths_octree ) ; call search_amr % input (PARAM=mesh_size_param_amr) 
       class is    ( maths_bin    ) ; call search_amr % input (PARAM=mesh_size_param_amr) 
       class default                ; call runend('MOD_AMR: WRONG SEARCH STRATEGY')
       end select
        
       call mesh_in    % element_bb (bobox)
       call search_amr % fill       (BOBOX=bobox)
       call mesh_in    % element_bb (bobox,ONLY_DEALLOCATE=.true.)
       !
       ! Compute centroids
       !
       call mesh_out % centroid(centroid)
       !
       ! Interpolation strategy
       !
       allocate(interp)
       call interp     % init       ()
       call interp     % input      (SEARCH_METHOD_SEQ=search_amr,&
            &                        INTERPOLATION_METHOD=INT_ELEMENT_VALUE,NAME='AMR')
       call interp     % preprocess (centroid,mesh_in)
       call interp     % values     (mesh_in % tags(2) % values,mesh_out % tags(2) % values)
       !
       ! Deallocate
       !
       call interp     % deallo     ()
       call search_amr % deallo     ()
       call mesh_out   % centroid   (centroid,ONLY_DEALLOCATE=.true.)
       deallocate(interp)
       if( associated(search_amr) ) deallocate(search_amr)
       !
       ! Postprocess
       !
       !call mesh_out % output (filename='mesh_out-'//integer_to_string(kfl_paral))
       !call mesh_out % results(xx=mesh_out % tags(2) % values,names='SIZES',filename='mesh_out-'//integer_to_string(kfl_paral),where='ON ELEMENTS')
       !call mesh_in  % output (filename='mesh_in-'//integer_to_string(kfl_paral))

       !call mesh_in  % results(xx=mesh_in % tags(2) % values,names='SIZES',filename='mesh_in-'//integer_to_string(kfl_paral),where='ON ELEMENTS')

    end if

    call mesh_in % deallo()

  end subroutine AMR_interpolate_size
  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume+abel
  !> @date    2021-10-06
  !> @brief   Interpolate mesh size on the nodes of the mesh
  !> @details Interpolate mesh size form mesh_in to mesh_out
  !> 
  !-----------------------------------------------------------------------
  
  subroutine AMR_interpolate_size_nodes(mesh_out,mesh_in)

    use mod_quality, only: isPointInTet
    
    type(mesh_type_basic), intent(inout) :: mesh_in
    type(mesh_type_basic), intent(inout) :: mesh_out
    type(maths_octree),    target        :: oct
    type(maths_bin),       target        :: bin
    real(rp),              pointer       :: bobox(:,:,:)
    real(rp),              pointer       :: points_to_interpolate(:,:)
    type(interpolation)                  :: interp
!    real(rp)                             :: xx(3)

    call messages_live('INTERPOLATE MESH SIZE ON NODES')

    if( associated(mesh_in % tags) .and. mesh_out % nelem > 0 ) then
       
       nullify(bobox)
       nullify(points_to_interpolate)                                                      
       !
       ! Select search method
       !
       select case ( kfl_mesh_size_amr )
       case        ( SEARCH_BIN        ) ; search_amr => bin
       case        ( SEARCH_OCTREE     ) ; search_amr => oct
       end select
       call search_amr % init ()
       !
       ! Set in which points we want to interpolate
       !
       points_to_interpolate => mesh_out%coord
       !
       ! Set search strategy
       !
       select type ( search_amr   )
       class is    ( maths_octree ) ; call search_amr % input (PARAM=mesh_size_param_amr) 
       class is    ( maths_bin    ) ; call search_amr % input (PARAM=mesh_size_param_amr) 
       end select

       call mesh_in    % element_bb (bobox)
       call search_amr % fill       (BOBOX=bobox)
       call mesh_in    % element_bb (bobox,ONLY_DEALLOCATE=.true.)
       !
       ! Interpolation strategy
       !
       call interp     % init       ()
       call interp     % input      (search_amr,INTERPOLATION_METHOD=INT_ELEMENT_INTERPOLATION) 
       call interp     % preprocess (points_to_interpolate,mesh_in)      
       call interp     % values     (mesh_in % tags(2) % values,mesh_out % tags(2) % values)
       !
       ! Postprocess
       !
       ! for nodal values
       !
!        call mesh_out % output (filename='mesh_out-'//integer_to_string(kfl_paral))
!        call mesh_out % results(xx=mesh_out % tags(2) % values,names='SIZES',filename='mesh_out-'//integer_to_string(kfl_paral),where='ON NODES')
!        call mesh_in  % output (filename='mesh_in-'//integer_to_string(kfl_paral))
!        call mesh_in  % results(xx=mesh_in % tags(2) % values,names='SIZES',filename='mesh_in-'//integer_to_string(kfl_paral),where='ON NODES')
       !
       ! for elem values
       !
!        call mesh_out % output (filename='mesh_out-'//integer_to_string(kfl_paral))
!        call mesh_out % results(xx=mesh_out % tags(2) % values,names='SIZES',filename='mesh_out-'//integer_to_string(kfl_paral),where='ON ELEMENTS')
!        call mesh_in  % output (filename='mesh_in-'//integer_to_string(kfl_paral))
!        call mesh_in  % results(xx=mesh_in % tags(2) % values,names='SIZES',filename='mesh_in-'//integer_to_string(kfl_paral),where='ON ELEMENTS')
       !
       !
       ! Deallocate
       !
       call interp     % deallo     ()       
       call search_amr % deallo     ()
       !call mesh_out   % centroid   (centroid,ONLY_DEALLOCATE=.true.)
    end if
  
    call mesh_in % deallo()
    
  end subroutine AMR_interpolate_size_nodes
  
end module mod_AMR
!> @}
