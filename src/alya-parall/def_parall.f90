!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Parall
!> @{
!> @file    def_parall.f90
!> @author  houzeaux
!> @date    2020-05-09
!> @brief   Variables
!> @details Variables of parall 
!-----------------------------------------------------------------------

module def_parall

  use def_kintyp_basic,  only : ip,rp,lg
  use def_kintyp_comm,   only : comm_data_par
  use def_kintyp_comm,   only : typ_optimum_npart
  use def_kintyp_domain, only : nelty

  !------------------------------------------------------------------------
  ! Parameters
  !------------------------------------------------------------------------

  integer(ip), parameter   :: &
       lun_domai_par = 5503,  &    ! Graph partition
       lun_trace_par = 5504,  &    ! Trace
       lun_conve_par = 5506,  &    ! Time statistics
       lun_rstar_par = 5516,  &    ! Restart
       lun_autom_par = 5517,  &    ! Automatic repartitioning
       lun_autor_par = 5518,  &    ! Automatic repartitioning report
       lun_aonlp_par = 10000       ! additional number for ONLY_PREPROCESS
  integer(ip),   parameter :: &
       nvarp_par=20                ! # postprocess variables

  integer(ip),   parameter      :: &
       OPENACC_POLICY_1PPDEV  = 0, & ! One process per device
       OPENACC_POLICY_NPPDEV  = 1, & ! Several processes per device
       OPENACC_POLICY_ALL2DEV = 2    ! All the processes use a device

  !------------------------------------------------------------------------
  ! File names
  !------------------------------------------------------------------------

  character(150)           :: &
       fil_rstar_par,         &    ! Restart file
       fil_conne_par               ! MPI conenctivity

  !------------------------------------------------------------------------
  ! Reapro
  !------------------------------------------------------------------------

  integer(ip)                            :: &
       npart_par,                           &    ! Number of subdomains
       npart_empty_par,                     &    ! Number of empty partitions
       kfl_ascii_par,                       &    ! Restart SCII format(=1)
       kfl_bytes_par,                       &    ! Integer bytes for files
       kfl_parti_par,                       &    ! Partition type
       kfl_fileh_par,                       &    ! File hierarchy file
       kfl_filio_par,                       &    ! Open and close files in preprocess
       kfl_weigh_par,                       &    ! Weighting of Metis graph
       kfl_repart_par,                      &    ! Repartitioning
       kfl_repart_module_par,               &    ! Module for repartitioning
       kfl_repart_criterion_par,            &    ! Criterion for repartitioning
       kfl_repart_method_par,               &    ! Method for repartitioning
       kfl_repart_freq_par,                 &    ! Partition frequency
       kfl_repart_windo_par,                &    ! Linear regression window
       kfl_repart_post_par,                 &    ! Repartitioning postprocess strategy
       kfl_virfi_par,                       &    ! Virtual files
       kfl_global_numbering_par,            &    ! Global numbering strategy
       kfl_order_exchange_par,              &    ! Order exchange
       kfl_connectivity_par,                &    ! Output statistics
       nsire_par,                           &    ! Number of simultaneous reading (restart)
       kfl_matri_par,                       &    ! Global matrix postprocess
       kfl_partition_par,                   &    ! Partition method
       kfl_interface_parti_par,             &    ! Interface partitioning strategy
       kfl_parseq_par,                      &    ! Parallel or sequential partitioning
       kfl_openacc_policy,                  &    ! Mapping rank-device policy
       kfl_openacc_ranks_per_dev,           &    ! Number of processes per device
       kfl_openacc_contiguous_ranks_per_dev,&    ! Number of consecutive processes per device
       kfl_openacc_streams_per_dev,         &    ! Streams per device
       boxes_coarse_par(3),                 &    ! Number of boxes coarse bin
       boxes_fine_par(3),                   &    ! Number of boxes fine bin
       sfc_criteria,                        &
       sfc_check,                           &
       sfc_dim_bin_core,                    &
       kfl_automatic_npart,                 &    ! Automatic partitioning
       kfl_output_partition,                &    ! Partition output
       kfl_output_node_comm_arrays,         &    ! Comm arrays output
       kfl_output_edge_comm_arrays               ! Comm arrays output
        
  real(rp)                               :: &
       rmbyt_par,                           &    ! Max number of Gb for virtual files
       vect_partition_par(3),               &    ! Direction of partition for oriented bin
       weights_elements_par(nelty),         &    ! Weights for element types
       weights_materials_par(10_ip),        &    ! Weights for materials 
       repart_threshold_lb,                 &    ! Threshold for repartitioning
       repart_toler,                        &    ! Tolerance for repartitioning
       repart_min_weight,                   &    ! Minimum weight per process
       repart_max_weight,                   &    ! Maximum weight per process
       repart_wfact,                        &    ! WLR weigh ingrease factor
       min_efficiency_par                        ! Minimum efficiency
  
  logical(lg)                               &
       kfl_openacc_multithreading                ! run parallel on CPU (threads)

  character(20)                          :: &
       method_redistribution_par                 ! Redistribution strategy
  
  type(typ_optimum_npart)                :: &
       optimum_partition                         ! Automatic partitioning
  
  !------------------------------------------------------------------------
  ! Others
  !------------------------------------------------------------------------

  integer(8)               :: &
       par_memor(2)                ! Memory counter  
  !
  ! Partition scalar data
  !
  integer(ip)              :: &
       nproc_par,             &    ! Number of processes
       iproc_par,             &    ! My process
       nbcol,                 &    ! Number of colours
       nneig,                 &    ! Local number of neighbour domains
       slfbo,                 &    ! Init of self boundary
       iproc_part_par              ! My order into the partition slaves

  real(rp) :: &
       cpu_paral(50)
  !
  ! Mesh graph
  !
  integer(ip), pointer     :: &
       padja_par(:),          &
       ladja_par(:)

  integer(ip), pointer     :: &
       lepar_par(:),          & ! Domain of every element
       lnpar_par(:),          & ! Domain of every node
       lbpar_par(:),          & ! Domain of every boundary element
       leper_par(:),          &
       lbper_par(:),          & ! Boundary permutation
       lneig_par(:),          & ! number of neightbours of every domain
       ginde_par(:,:),        &
       lcomm_par(:,:),        &
       slfbo_par(:),          &
       leind_par(:),          &
       lbind_par(:),          &
       lnods_par(:,:),        &
       nhang_par(:),          &
       xlnin_loc(:)           

  integer(ip), pointer     :: &
       neighDom(:),           &
       xadjDom(:),            &
       adjDom(:),             &
       translDual(:),         &
       iaDual(:),             &
       jaDual(:),             &
       colours(:)

  integer(ip), pointer     :: &
       badj(:),               &
       bdom(:),               &
       bpoin(:)

  integer(ip), pointer     :: &
       leinv_par(:),          &  ! inverse perm of elements
       lbinv_par(:),          &  ! inverse perm of boundaries
       lnper_par(:),          &  ! perm of nodes
       lninv_par(:)              ! inverse perm of nodes (1..npoin)
  !
  ! Asynchronous communications
  !
  integer(4),  pointer     :: &
       ireq4(:)                  ! Request for asynchronous communication
  integer(ip)              :: &
       ipass_par                 ! Request for asynchronous communication
  !
  ! Mesh multiplication
  !
  integer(ip), pointer     :: &
       lowns_par(:),          &  ! List of my own nodes
       lownr_par(:)              ! List of my neighbor's own nodes
  !
  ! Sets
  !
  integer(ip), pointer     :: &
       lnsec_par(:,:)            ! Node set
  integer(ip), pointer     :: &
       nnset_par(:),          &  ! Node set per subdomain
       nwitn_par(:)              ! Witness points per subdomain
  !
  ! Repartitioning
  !
  integer(ip)              :: &
       num_repart_par            ! Number of repartitioning

end module def_parall
