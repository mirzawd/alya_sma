!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine par_reapro()
  !------------------------------------------------------------------------
  !****f* Parall/par_reapro
  ! NAME
  !    par_reapro
  ! DESCRIPTION
  !    This routine reads service data
  ! OUTPUT
  ! USED BY
  !    Parall
  !***
  !------------------------------------------------------------------------
  
  use def_parame
  use def_master
  use def_parall
  use def_parall,                    only : kfl_openacc_policy
  use def_parall,                    only : kfl_openacc_ranks_per_dev
  use def_parall,                    only : kfl_openacc_contiguous_ranks_per_dev
  use def_parall,                    only : kfl_openacc_streams_per_dev
  use def_parall,                    only : kfl_openacc_multithreading
  use def_inpout
  use mod_elmgeo,                    only : element_type
  use def_domain,                    only : nelty
  use mod_parall,                    only : par_omp_granularity
  use mod_parall,                    only : par_omp_granularity
  use mod_parall,                    only : par_omp_coloring_alg
  use mod_parall,                    only : par_topo_num_cores_per_node
  use mod_parall,                    only : par_topo_num_nodes
  use mod_parall,                    only : par_hybrid
  use mod_parall,                    only : par_one_sided
  use mod_parall,                    only : PAR_OPENMP_COLORING
  use mod_parall,                    only : PAR_OPENMP_NO_COLORING
  use mod_parall,                    only : PAR_OMPSS
  use mod_parall,                    only : PAR_HYBRID_OFF
  use mod_parall,                    only : par_omp_nelem_chunk
  use mod_parall,                    only : par_omp_nboun_chunk
  use mod_parall,                    only : par_omp_npoin_chunk
  use mod_parall,                    only : par_omp_partition_alg
  use mod_parall,                    only : PAR_METIS4
  use mod_parall,                    only : PAR_SFC
  use mod_parall,                    only : PAR_ZOLTAN
  use mod_parall,                    only : PAR_ORIENTED_BIN
  use mod_parall,                    only : PAR_NUMBERING
  use mod_parall,                    only : PAR_USING_RANK 
  use mod_parall,                    only : PAR_RANDOM
  use mod_parall,                    only : PAR_SEQUENTIAL_PARTITION
  use mod_parall,                    only : PAR_PARALLEL_PARTITION
  use mod_parall,                    only : PAR_WEIGHT_GAUSS     
  use mod_parall,                    only : PAR_WEIGHT_OFF       
  use mod_parall,                    only : PAR_WEIGHT_ELEMENT   
  use mod_parall,                    only : PAR_WEIGHT_SQUARE    
  use mod_parall,                    only : PAR_WEIGHT_MATERIAL
  use mod_parall,                    only : PAR_WEIGHT_STUPID
  use mod_ecoute,                    only : ecoute
  use mod_messages,                  only : messages_live
  use mod_optimum_partition,         only : optimum_partition_read
  implicit none
  integer(ip)  :: ielty,imate,imodu,kfl_servi
  character(5) :: welem

  kfl_servi = 0
  
  if( IMASTER .or. ISEQUEN ) then

     npart_par                 = nproc_par-1                 ! Number of subdomains
     npart                     = npart_par                   ! IDEM
     npart_empty_par           = 0                           ! Empty partitions
     kfl_ascii_par             = 0                           ! Restart file format
     kfl_bytes_par             = ip                          ! Integer bytes for restart files
     kfl_parti_par             = 2                           ! Partition type (nodes/faces)
     kfl_postp_par             = 1                           ! Postprocess type (master  =1, slave  =0)
     kfl_fileh_par             = 0                           ! No file hierarchy
     kfl_filio_par             = 0                           ! Do not open and close files in preprocess
     kfl_virfi_par             = 0                           ! No virtual files
     kfl_global_numbering_par  = 0                           ! Global numbering strategy
     kfl_order_exchange_par    = 0                           ! Order exchange update
     kfl_connectivity_par      = 0                           ! MPI connectivity
     kfl_async                 = 1                           ! Asynchronous (by default)
     nsire_par                 = nproc_par+1                 ! Number of simultaneous readings
     rmbyt_par                 = 1                           ! Max number of Gb for virtual files
     par_omp_granularity       = 10                          ! Granularity
     par_omp_coloring_alg      = 1                           ! Default openmp coloring algorithm
     par_omp_nelem_chunk       = 0                           ! Element chunk size
     par_omp_nboun_chunk       = 0                           ! Node chunk size
     par_omp_npoin_chunk       = 0                           ! Boundary chunk size
     par_omp_partition_alg     = PAR_METIS4                  ! Partition method for OmpSS
     kfl_matri_par             = 0                           ! Global matrix postprocess
     par_one_sided             = 0                           ! No one-sided
     
     kfl_partition_par         = PAR_METIS4                  ! Partition method
     kfl_parseq_par            = PAR_SEQUENTIAL_PARTITION    ! Sequential partitioning
     kfl_interface_parti_par   = 1                           ! Interface partitioning strategy
     kfl_weigh_par             = PAR_WEIGHT_GAUSS            ! Weight on Gauss points
     kfl_repart_par            = 0                           ! No repartitioning
     kfl_repart_freq_par       = 10                          ! Repartioting frequency
     kfl_repart_windo_par      = 10                          ! Window for linear regression
     kfl_repart_module_par     = ID_KERMOD                   ! Module to be used for repartitioning
     kfl_repart_criterion_par  = CPU_ASSEMBLY                ! Assemlby is the criterion for load balance
     kfl_repart_method_par     = CPU_REPART_LOAD             ! LOAD reescaling is the method for load blance
     kfl_repart_post_par       = 0                           ! Repartition is indicated with a tag
     weights_elements_par      = 1.0_rp                      ! Relative weights of elements
     weights_materials_par     = 1.0_rp                      ! Relative weights of elements
     repart_threshold_lb       = 0.9_rp                      ! threshold load balance for repartitioning
     repart_toler              = 0.00_rp                     ! tolerance used for repartitioning
     repart_min_weight         = 0.01_rp                     ! Minimum normalized weight per mpi-process for redistribuiton
     repart_max_weight         = 100.0_rp                    ! Maximum normalized weight per mpi-process for redistribuiton
     repart_wfact              = 1.0_rp                      ! Increase weight factor for WLR
     boxes_coarse_par          = 0                           ! Number of boxes coarse bin
     boxes_fine_par            = 0                           ! Number of boxes fine bin
     vect_partition_par        = (/ 1.0_rp,0.0_rp,0.0_rp /)  ! Direction of partition for oriented bin
     vect_partition_par        = 0.0_rp
     
     kfl_openacc_policy                   = 0_ip             ! Rank-device mapping policy
     kfl_openacc_ranks_per_dev            = 1_ip             ! Number of processes per device
     kfl_openacc_contiguous_ranks_per_dev = 1_ip             ! Number of contiguous processes per device
     kfl_openacc_streams_per_dev          = 1_ip             ! Number of streams per device
     kfl_openacc_multithreading           = .False.          ! Run parallel on CPU (threads)

     sfc_check                 = 0_ip
     sfc_criteria              = 2_ip
     sfc_dim_bin_core          = 256_ip
 
     method_redistribution_par = 'SYNCHRONOUS'               ! Redistribution strategy

     kfl_output_partition        = 1                         ! Partition output
     kfl_output_node_comm_arrays = 0                         ! Comm arrays output
     kfl_output_edge_comm_arrays = 0                         ! Comm arrays output
     !
     ! Reach the section
     !
     !rewind(lisda)
     !do while(words(1)/='RUNDA')
     !   call ecoute('PAR_REAPRO')
     !end do
     !do while(words(1)/='ENDRU')
     !   call ecoute('PAR_REAPRO')
     !end do
     !do while(words(1)/='PROBL')
     !   call ecoute('PAR_REAPRO')
     !end do
     !
     ! Read data
     !
     !do while(words(1)/='ENDPR')
        !call ecoute('PAR_REAPRO')
        !if(words(1)=='PARAL') then
           if(exists('ON   ')) then
              kfl_servi =1
              do while(words(1)/='ENDPA')
                if( words(1) == 'IO' .and. trim(words(2)) == 'ON   ' ) then
                   !
                   ! Partitioning parameters
                   !
                   call runend('PAR_REAPRO: MPI IO SECTION IS NOW APART FROM PARALL SECTION')

                else if( words(1) == 'OPTIM' ) then
                   !
                   ! Automatic partitioning
                   !
                   call optimum_partition_read()
                   
                else if( words(1) == 'REDIS' ) then
                   !
                   ! Redistribution strategy, sendrecv or alltoallv
                   !
                   if( words(2) == 'SYNCH' ) then
                      method_redistribution_par = 'SYNCHRONOUS'
                   else if( words(2) == 'ASYNC' ) then
                      method_redistribution_par = 'ASYNCHRONOUS'
                   else if( words(2) == 'ALLTO' ) then
                      method_redistribution_par = 'ALLTOALLV'
                   end if
                   
                else if( words(1) == 'OUTPU' ) then

                   if( trim(words(2)) == '' ) then
                      do while( words(1) /= 'ENDOU' )
                         if( words(1) == 'PARTI' ) then
                            if( option('PARTI') ) then                    ! Partition output
                               kfl_output_partition = 1       
                            else
                               kfl_output_partition = 0
                            end if
                         else if( words(1) == 'NODEC' ) then
                            if( option('NODEC') ) then                    ! Comm arrays output
                               kfl_output_node_comm_arrays = 1
                            else
                               kfl_output_node_comm_arrays = 0
                            end if
                         else if( words(1) == 'EDGEC' ) then
                            if( option('EDGEC') ) then                    ! Comm arrays output
                               kfl_output_edge_comm_arrays = 1
                            else
                               kfl_output_edge_comm_arrays = 0
                            end if
                         end if
                         call ecoute('par_reapro')
                      end do
                   else
                      !
                      ! Output done by master or slave
                      !
                      if( exists('CONNE') ) kfl_connectivity_par = 1
                   end if
                   
                 else if( words(1) == 'POSTP' ) then
                    !
                    ! Postprocess
                    !
                    if(words(2)=='MASTE' ) then
                       kfl_postp_par=1
                    else
                       kfl_postp_par=0
                    end if

                 else if( words(1) == 'MATRI' ) then
                    !
                    ! Global matrix postprocess
                    !
                    if( option('MATRI') ) kfl_matri_par = 1

                 else if( words(1) == 'TASK ' ) then
                    !
                    ! Task: all/pre or post
                    !
                    if( words(2)=='ONLYP' .or. words(2)=='PREPR' .or. words(2)=='WRITE' ) then

                       kfl_ptask=0
                       call vocabu(-1_ip,0_ip,0_ip)
                       if( nproc_par > 1 ) then
                          npart_par = nproc_par-1
                          npart     = nproc_par-1
                       else
                          npart_par = getint('SUBDO',1_ip,'#Number of subdomains')
                          npart     = npart_par
                       end if
                       if( exists('ASCII') .or. exists('FORMA') ) kfl_ascii_par = 1
                       if( exists('BINAR') .or. exists('UNFOR') ) kfl_ascii_par = 0
                       if( exists('FOURB') .or. exists('4BYTE') ) kfl_bytes_par = 4
                       if( exists('EIGHT') .or. exists('8BYTE') ) kfl_bytes_par = 8

                    else if(words(2)=='PARTI' ) then

                       kfl_ptask=1
                       call vocabu(-1_ip,0_ip,0_ip)

                    else if(words(2)=='READP' ) then

                       kfl_ptask=2
                       call vocabu(-1_ip,0_ip,0_ip)
                       if(exists('SIMUL')) then
                          nsire_par=getint('SIMUL',nproc_par+1_ip,'#Number of simultaneous readings')
                       end if
                       if( exists('ASCII') .or. exists('FORMA') ) kfl_ascii_par = 1
                       if( exists('BINAR') .or. exists('UNFOR') ) kfl_ascii_par = 0

                    end if
                    if(npart_par<2.and.nproc_par==1) then
                       print*,'npart_par',npart_par
                       call runend('PAR_REAPRO:WRONG NUMBER OF SUBDOMAINS')
                    end if

                 else if( words(1) == 'ORDER') then

                    if( option('ORDER') ) then
                       kfl_order_exchange_par = 1
                    else
                       kfl_order_exchange_par = 0
                    end if
                    
                 else if( words(1) == 'GLOBA' ) then
                    !
                    ! Global numbering
                    !
                    if( words(2) == 'LEXIC' ) then
                       kfl_global_numbering_par = 1
                    else
                       kfl_global_numbering_par = 0
                    end if
                    
                 else if( words(1) == 'PARTI' .and. trim(words(2)) == '' ) then
                    !
                    ! Partitioning parameters
                    !
                    call ecoute('par_reapro')
                    do while( words(1) /= 'ENDPA' )

                       if( words(1) == 'EXECU' ) then
                          !
                          ! Execution mode: Sequential or parallel partitioning
                          !
                          if( words(2) == 'SEQUE' ) then
                             kfl_parseq_par = PAR_SEQUENTIAL_PARTITION
                          else
                             kfl_parseq_par = PAR_PARALLEL_PARTITION
                          end if

                       else if( words(1) == 'METHO' ) then
                          !
                          ! Partitioning method
                          !
                          if( words(2) == 'METIS' ) then
                             kfl_partition_par    = PAR_METIS4
                          else if( words(2) == 'SFC  ' ) then
                             kfl_partition_par    = PAR_SFC
                          else if( words(2) == 'ZOLTA' ) then
                             kfl_partition_par    = PAR_ZOLTAN
                          else if( words(2) == 'ORIEN' ) then
                             kfl_partition_par    = PAR_ORIENTED_BIN
                          else if( words(2) == 'FIELD' ) then
                             kfl_partition_par    = -getint('FIELD',-1_ip,'#Partition is taken from an element field')
                          else if( words(2) == 'NUMBE' ) then
                             kfl_partition_par    = PAR_NUMBERING
                          else if( words(2) == 'RANK ' ) then
                             kfl_partition_par    = PAR_USING_RANK 
                          else if( words(2) == 'RANDO' ) then
                             kfl_partition_par    = PAR_RANDOM 
                          else
                             call runend('PAR_REAPRO: UNKNOWN PARTITIONING METHOD')
                          end if
                          
                       else if( words(1) == 'ELEME' ) then
                          !
                          ! Element graph
                          !
                          if( words(2) == 'NODES' ) then
                             kfl_parti_par = 1
                          else if( words(2) == 'FACES' ) then
                             kfl_parti_par = 2
                          end if

                      else if( words(1) == 'INTER' ) then
                          !
                          ! Interface partitioning strategy
                          !
                          if(      words(2) == 'METIS' ) then
                             kfl_interface_parti_par = 1
                          else if( words(2) == 'CHUNK' ) then
                             kfl_interface_parti_par = 2
                          end if

                       else if( words(1) == 'BOXES' ) then
                          !
                          ! Number of boxes
                          !
                          if( exists('COARS') ) then
                             boxes_coarse_par(1:3) = int(param(1:3)) ! Coarse bin
                          else
                             boxes_fine_par(1:3)   = int(param(1:3)) ! Fine bin
                          end if

                       else if( words(1) == 'CRITE' ) then
                          !
                          sfc_criteria= int(param(1),ip)

                       else if( words(1) == 'DIMBI' ) then
                          !
                          sfc_dim_bin_core= int(param(1),ip)

                       else if( words(1) == 'CHECK' ) then
                          !
                          if(words(2)=='ON   ') then
                                sfc_check = 1_ip
                          end if

                       else if( words(1) == 'DIREC' ) then
                          !
                          ! Direction of oriented bin
                          !
                          if(      words(2) == 'X    ' ) then
                             vect_partition_par(1:3)   = (/ 1.0_rp,0.0_rp,0.0_rp /)
                          else if( words(2) == 'Y    ' ) then
                             vect_partition_par(1:3)   = (/ 0.0_rp,1.0_rp,0.0_rp /)
                          else if( words(2) == 'Z    ' ) then
                             vect_partition_par(1:3)   = (/ 0.0_rp,0.0_rp,1.0_rp /)
                          else
                             vect_partition_par(1:3) = param(1:3)
                          end if

                       else if( words(1) == 'WEIGH' ) then
                          !
                          ! Weight for elements
                          !
                          if(      words(2) == 'GAUSS' ) then
                             kfl_weigh_par = PAR_WEIGHT_GAUSS 
                          else if( words(2) == 'OFF  ' .or. words(2) == 'NONE ' ) then
                             kfl_weigh_par = PAR_WEIGHT_OFF                        
                          else if( words(2) == 'SQUAR' ) then
                             kfl_weigh_par = PAR_WEIGHT_SQUARE
                          else if( words(2) == 'FIELD' ) then
                             kfl_weigh_par =  getint('FIELD',1_ip,'#Field for partition weights')
                          else if( words(2) == 'ELEME' ) then
                             kfl_weigh_par = PAR_WEIGHT_ELEMENT
                             do ielty = 1,nelty
                               welem =  element_type(ielty) % name 
                               if( exists(welem) ) then
                                  weights_elements_par(ielty) = getrea(welem,1.0_rp,'#Element weight for partition')
                                end if
                             end do
                          else if( words(2) == 'MATER' ) then
                             kfl_weigh_par = PAR_WEIGHT_MATERIAL
                             do imate = 1,min(nnpar,size(weights_materials_par,KIND=ip))
                                weights_materials_par(imate) = param(imate+1)
                             end do
                          else if( words(2) == 'STUPI' ) then
                             kfl_weigh_par = PAR_WEIGHT_STUPID
                          end if
                          
                       else if( words(1) == 'EMPTY' ) then
                          !
                          ! Number of empty partitions
                          !
                          npart_empty_par = getint('EMPTY',0_ip,'#NUMBER OF EMPTY PARTITIONS')

                       end if

                       call ecoute('par_reapro')
                    end do

                 else if( words(1) == 'REPAR' ) then
                    !
                    ! Repartitioning: REPARTITIONING: ON, MODULE=NASTIN, CRITERION=ASSEMBLY
                    !
                    if( words(2) == 'OFF  ' ) then
                       do while( words(1) /= 'ENDRE' )
                          call ecoute('par_reapro')
                       end do
                    end if
                    do while( words(1) /= 'ENDRE' )
                       kfl_repart_par = 1
                       if( words(1) == 'MODUL' ) then
                          !
                          ! Leading module
                          !
                          imodu = idmod(words(2))
                          kfl_repart_module_par = imodu
                          if( imodu == -1 ) call runend('PAR_REAPRO: WRONG MODULE')

                       else if( words(1) == 'CRITE' ) then
                          !
                          ! Criterion
                          !
                          if(      words(2) == 'ASSEM' ) then
                             kfl_repart_criterion_par = CPU_ASSEMBLY
                          else if( words(2) == 'MINAS' ) then
                             kfl_repart_criterion_par = CPU_MINI_ASSEMBLY
                          else if( words(2) == 'ALWAY' ) then
                             kfl_repart_criterion_par  = -1
                          else
                             call runend('PAR_REAPRO: UNKNOWN CRITERION')
                          end if
                       
                       else if( words(1) == 'THRES' ) then
                          !
                          ! Threshold load balance
                          !
                          repart_threshold_lb = getrea('THRES',1.0_rp,'#THRESHOLD FOR LOAD BALANCE')
                       
                       else if( words(1) == 'TOLER' ) then
                          !
                          ! Threshold load balance
                          !
                          repart_toler = getrea('TOLER',0.01_rp,'#TOLERANCE FOR LOAD BALANCE')
                       
                       else if( words(1) == 'MINIM' ) then
                          !
                          ! Minimum weight per process
                          !
                          repart_min_weight = getrea('MINIM',0.01_rp,'#MINIM WEIGHT PER PROCESS (NORMALIZED)')
                       
                       else if( words(1) == 'MAXIM' ) then
                          !
                          ! Minimum weight per process
                          !
                          repart_max_weight = getrea('MAXIM',100.0_rp,'#MAXIM WEIGHT PER PROCESS (NORMALIZED)')

                       else if( words(1) == 'FREQU' ) then
                          !
                          ! Partition frequency
                          !
                          kfl_repart_freq_par = getint('FREQU',10_ip,'#FREQUENCY FOR REPARTITIONING')
                       
                       else if( words(1) == 'WINDO' ) then
                          !
                          ! Window of the linear regression
                          !
                          kfl_repart_windo_par = getint('WINDO',10_ip,'#WINDOW FOR LINAR REGRESSION')
                       
                       else if( words(1) == 'POSTP' ) then
                          !
                          ! Postprocess file name strategy
                          !
                          if(      words(2) == 'TAG  ' ) then
                             kfl_repart_post_par = 1
                          else if( words(2) == 'DIREC' ) then
                             kfl_repart_post_par = 2
                          else if( words(2) == 'NONE ' ) then
                             kfl_repart_post_par = 0
                          else
                             call runend('PAR_REAPRO: UNKNOWN REPARTITION POSTPROCESS STRATEGY')
                          end if                          
                       
                       else if( words(1) == 'METHO' ) then
                          !
                          ! Criterion
                          !
                          if(      words(2) == 'WLR  ' ) then
                             kfl_repart_method_par = CPU_REPART_WLR
                             if ( words(3) == 'EXP  ') then
                                repart_wfact = getrea('EXP  ',1.5_rp,'WEIGHT ENCREASE FACTOR')
                             else if ( words(3) == '     ') then
                                repart_wfact = 1.5_rp
                             else
                                call runend('PAR_REAPRO: UNKNOWN WLR STRATEGY')
                             endif
                          else if( words(2) == 'OLR  ' ) then
                             kfl_repart_method_par = CPU_REPART_OLR
                          else if( words(2) == 'LOAD ' ) then
                             kfl_repart_method_par  = CPU_REPART_LOAD
                          else
                             call runend('PAR_REAPRO: UNKNOWN REPARTITION METHOD')
                          end if
                       end if
                       call ecoute('par_reapro')
                       
                    end do                   

                  else if( words(1) == 'PARTI' ) then
                    !
                    ! Partitioning type
                    !
                    if(words(2)=='NODES' ) then
                       kfl_parti_par=1
                    else if(words(2)=='FACES' ) then
                       kfl_parti_par=2
                    end if

                 else if( words(1) == 'FILEH' ) then
                    !
                    ! File hierarchy
                    !
                    if(words(2)=='NONE ' ) then
                       kfl_fileh_par=0
                    else if(words(2)=='YES  '.or.words(2)=='ONELE'.or.words(2)=='ON   ' ) then
                       kfl_fileh_par=1
                    else if(words(2)=='TWOLE' ) then
                       kfl_fileh_par=2
                    end if

                 else if( words(1) == 'FILEO' ) then
                    !
                    ! Open/Close files
                    !
                    if(words(2)=='YES  '.or.words(2)=='ON   ' ) then
                       kfl_filio_par=1
                    else
                       kfl_filio_par=0
                    end if

                 else if( words(1) == 'WEIGH' ) then
                    !
                    ! Weight for graph
                    !
                    if(      words(2) == 'GAUSS' ) then
                       kfl_weigh_par =  0
                    else if( words(2) == 'OFF  ' ) then
                       kfl_weigh_par = -1                            
                    else if( words(2) == 'FIELD' ) then
                       kfl_weigh_par =  getint('FIELD',1_ip,'#Field for partition weights')
                    end if

                 else if( words(1) == 'VIRTU' ) then
                    !
                    ! Virtual file
                    !
                    if(words(2)=='ON   ' ) then
                       kfl_virfi_par = 1
                       kfl_ascii_par = 0
                       rmbyt_par     = getrea('MAXIM',1.0_rp,'#MAXIMUM NUMBER OF GB FOR VIRTUAL FILES')
                    else
                       kfl_virfi_par = 0
                    end if

                 else if( words(1) == 'COMMU' ) then
                    if( words(2) == 'SYNCH' ) then
                       kfl_async = 0
                    else if( words(2) == 'ASYNC' ) then
                       kfl_async = 1
                    end if

                 else if( words(1) == 'TOPOL' ) then
                    !
                    ! Topology of the supercomputer
                    !
                    call ecoute('PAR_REAPRO')
                    do while( words(1) /= 'ENDTO' )
                       if(      words(1) == 'NODES' ) then
                          par_topo_num_nodes          = getint('NODES',1_ip,'#Number of computing nodes')
                       else if( words(1) == 'CORES' ) then
                          par_topo_num_cores_per_node = getint('CORES',1_ip,'#Number of core per node')
                       end if
                       call ecoute('PAR_REAPRO')
                    end do

                 else if( words(1) == 'ONESI' ) then
                    !
                    ! One sided
                    !
                    if( option('ONESI') ) par_one_sided = 1
                    
                 else if( words(1) == 'HYBRI' ) then
                    !
                    ! Hybrid parallelization
                    !
                    if( words(2) == 'OPENM' ) then
                       if( exists('NOCOL') ) then
                          par_hybrid = PAR_OPENMP_NO_COLORING
                       else if( exists('COLOR') ) then
                          par_hybrid = PAR_OPENMP_COLORING
                       else
                          call runend('PAR_REAPRO: UNKNOWN HYBRID PARALLELIZATION')
                       end if
                    else if( words(2) == 'OMPSS' ) then
                       par_hybrid = PAR_OMPSS
                    else if( words(2) == 'OFF  ' ) then
                       par_hybrid = PAR_HYBRID_OFF
                    else
                       call runend('PAR_REAPRO: UNKNOWN HYBRID PARALLELIZATION')
                    end if

                 else if( words(1) == 'OPENM' ) then
                    !
                    ! OpenMP stuffs
                    !
                    call ecoute('PAR_REAPRO')
                    do while( words(1) /= 'ENDOP' )
                       if( words(1) == 'GRANU' ) then
                          !
                          ! Granularity
                          !
                          par_omp_granularity = getint('GRANU',10_ip,'#Granularity of OpenMP blocks')

                       else if( words(1) == 'COLOR' ) then
                          !
                          ! Coloring technique
                          !
                          if( words(2) == 'MINCO' ) then
                            par_omp_coloring_alg = 1
                          else
                            par_omp_coloring_alg = 0
                          end if

                       else if( words(1) == 'ELEME' ) then
                          !
                          ! Element chunk size
                          !
                          par_omp_nelem_chunk = getint('ELEME',100_ip,'#Element chunk size')

                       else if( words(1) == 'NODEC' ) then
                          !
                          ! Node chunk size
                          !
                          par_omp_npoin_chunk = getint('NODAC',100_ip,'#Node chunk size')

                       else if( words(1) == 'BOUND' ) then
                          !
                          ! Node chunk size
                          !
                          par_omp_nboun_chunk = getint('BOUND',100_ip,'#Boundary chunk size')

                       else if( words(1) == 'PARTI' ) then
                          !
                          ! Partitioning method for OmpSS
                          !
                          if( words(2) == 'METIS' ) then
                             par_omp_partition_alg = PAR_METIS4
                          else if( words(2) == 'SFC  ' ) then
                             par_omp_partition_alg = PAR_SFC
                          else
                             call runend('PAR_REAPRO: WRONG PARTITIONING FOR OMPSS') 
                          end if
                          
                       end if
                       call ecoute('PAR_REAPRO')
                    end do
                 else if( words(1) == 'OPENA' ) then
                    !
                    ! OpenACC
                    !
                    call ecoute('PAR_REAPRO')
                    do while( words(1) /= 'ENDOP' )

                      if( words(1) == 'POLIC') then
                        kfl_openacc_policy = getint('POLIC',0_ip,'#OpenACC Rank-device mapping policy')
                      end if
                      if( words(1) == 'RANKS') then
                        kfl_openacc_ranks_per_dev = getint('RANKS',1_ip,'#OpenACC # of processes per device')
                      end if
                      if( words(1) == 'CONTI') then
                        kfl_openacc_contiguous_ranks_per_dev = getint('CONTI',1_ip,'#OpenACC # of contiguous processes per device')
                      end if
                      if( words(1) == 'STREA') then
                        kfl_openacc_streams_per_dev = getint('STREA',1_ip,'#OpenACC # of streams per device')
                      end if
                      if( words(1) == 'MULTI') then
                        kfl_openacc_multithreading = (exists('ON   ') .or. exists('YES  '))
                      end if

                      call ecoute('PAR_REAPRO')
                    end do

                 end if
                 call ecoute('PAR_REAPRO')
              end do
           end if
        !end if
     !end do
     !
     ! Check errors and warnings
     !
     if(IMASTER.and.kfl_servi/=1) &
          call runend('MPI WAS INITIATED AND PARALL SERVICE IS OFF')
     if(IMASTER.and.npart_par<=1) &
          call runend('WRONG NUMBER OF SUBDOMAINS TO USE PARALL SERVICE')
     if(kfl_ptask==0 .and. kfl_filio_par==0 .and. npart_par>1000) then
        call messages_live('OPEN AND CLOSE ONCE RESTART FILES MAY LEAD TO PROBLEMS','WARNING')
     end if
     if( kfl_ptask == 2 ) kfl_virfi_par = 0

  end if

end subroutine par_reapro

