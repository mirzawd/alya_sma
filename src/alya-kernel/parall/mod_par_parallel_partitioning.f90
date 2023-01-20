!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Parall
!> @{
!> @file    mod_par_parallel_partitioning.f90
!> @author  houzeaux
!> @date    2018-10-09
!> @brief   Parallel partitioning
!> @details Tools for Parallel partitioning
!>          WARNING: if you introduce a new way of imposing wiehgts on
!>          the elements, you have to edit 3 subroutines
!>          par_partition_total_weight
!>          par_partition_weights_rp
!>          par_partition_weights_ip
!>
!>          kfl_weigh_par =  0 ...  Weight on Gauss points
!>                        = -1 ...  No weight    
!>                        = -2 ...  Weight on element types
!>                        = -3 ...  Weight on square of Gauss points
!>                        = -4 ...  Weight on material
!>                        >  0 ...  Weights read from a field
!>
!-----------------------------------------------------------------------

module mod_par_parallel_partitioning

  use def_master
  use def_elmtyp
  use def_kintyp_comm,            only : comm_data_par
  use def_kintyp_basic,           only : ip,rp
  use def_master,                 only : kfl_paral
  use def_domain,                 only : ndime
  use def_domain,                 only : nelem
  use def_domain,                 only : ltype
  use def_domain,                 only : ngaus
  use def_domain,                 only : xfiel
  use def_domain,                 only : lmate
  use def_domain,                 only : npoin
  use def_domain,                 only : lnods
  use def_domain,                 only : coord
  use def_domain,                 only : nnode
  use def_domain,                 only : memor_dom
  use def_parall,                 only : boxes_coarse_par
  use def_parall,                 only : boxes_fine_par
  use def_parall,                 only : kfl_partition_par
  use def_parall,                 only : lepar_par
  use def_parall,                 only : kfl_weigh_par
  use def_parall,                 only : vect_partition_par
  use def_parall,                 only : weights_elements_par
  use def_parall,                 only : weights_materials_par
  use def_parall,                 only : npart_empty_par
  use mod_parall,                 only : PAR_SFC
  use mod_parall,                 only : PAR_ORIENTED_BIN
  use mod_parall,                 only : PAR_ZOLTAN
  use mod_parall,                 only : PAR_USING_RANK
  use mod_parall,                 only : PAR_RANDOM
  use mod_parall,                 only : commd
  use mod_parall,                 only : par_memor
  use mod_parall,                 only : PAR_COMM_MY_CODE_ARRAY
  use mod_parall,                 only : PAR_COMM_MY_CODE_WM
  use mod_parall,                 only : PAR_COMM_MY_CODE
  use mod_parall,                 only : PAR_CODE_SIZE
  use mod_parall,                 only : PAR_COMM_WORLD
  use mod_parall,                 only : PAR_COMM_MY_CODE_ARRAY
  use mod_parall,                 only : PAR_WEIGHT_GAUSS     
  use mod_parall,                 only : PAR_WEIGHT_OFF       
  use mod_parall,                 only : PAR_WEIGHT_ELEMENT   
  use mod_parall,                 only : PAR_WEIGHT_SQUARE    
  use mod_parall,                 only : PAR_WEIGHT_MATERIAL  
  use mod_parall,                 only : PAR_WEIGHT_STUPID  
  use mod_communications,         only : PAR_COMM_RANK_AND_SIZE
  use mod_communications,         only : PAR_BARRIER
  use mod_communications,         only : PAR_BROADCAST
  use mod_communications,         only : PAR_MAX
  use mod_communications,         only : PAR_SUM
  use mod_par_interface_exchange, only : par_interface_exchange
  use mod_partition_sfc,          only : partition_sfc
  use mod_memory,                 only : memory_alloca
  use mod_memory,                 only : memory_deallo
  use mod_memory,                 only : memory_size
  use mod_memory,                 only : memory_resize
  use mod_memory,                 only : memory_copy
  use mod_redistribution,         only : redistribution_comm_initialization
  use mod_redistribution,         only : redistribution_domain
  use mod_messages,               only : messages_live
  use mod_partitioning,           only : partitioning
  use mod_domain,                 only : domain_number_gauss_points

  implicit none

  private
  
  real(rp), pointer :: coord_sav(:,:)

  interface par_partition_weights
     module procedure &
          par_partition_weights_rp,&
          par_partition_weights_ip
  end interface par_partition_weights

  interface par_partition_total_weight
     module procedure &
          par_partition_total_weight_ip
  end interface par_partition_total_weight

  public :: par_parallel_partitioning
  public :: par_partition_weights
  public :: par_partition_total_weight

  character(100), PARAMETER :: vacal = "mod_par_parallel_partitioning"
  
contains

  !-----------------------------------------------------------------------
  !>
  !> @author  houzeaux
  !> @date    2018-10-09
  !> @brief   Main subroutine
  !> @details Perform parallel partitioning
  !>
  !-----------------------------------------------------------------------

  subroutine par_parallel_partitioning(READ_MESH,LCORRECT,TIMINGS)

    use def_domain
    use mod_memory
    use def_coupli
    use mod_mpio_config, only: mpio_config
    
    logical(lg),       intent(in),  optional   :: READ_MESH
    real(rp), pointer, intent(in),  optional   :: LCORRECT(:)
    real(rp),          intent(out), optional   :: TIMINGS(:)
    integer(ip)                                :: pnode,npart_sav,ifiel
!    integer(ip)                                :: ipart
    integer(ip)                                :: idime,ielem,ielty,max_subd
    integer(ip)                                :: npart_empty,npart_eff
    integer(ip), pointer                       :: ia(:),ja(:)       ! We had to put this for PGI
    logical(lg), pointer                       :: lmask(:)          ! We had to put this for PGI
    integer(ip), pointer                       :: lweig(:)
    real(rp),    pointer                       :: xcoor(:,:)
    real(rp),    pointer                       :: rweig(:)
    integer(ip), pointer                       :: lsubd(:)
    logical(lg)                                :: if_read_mesh
    real(rp)                                   :: time1,time2

    call messages_live('PARALLEL PARTITIONING','START SECTION')

    nullify(xcoor,rweig,lweig,ia,ja,lmask)
    nullify(ia,ja,lmask)
    nullify(lsubd)
    
    if( present(READ_MESH) ) then
       if_read_mesh = READ_MESH
    else
       if_read_mesh = .true.       
    end if

    call cputim(time1)
    if( if_read_mesh ) then

       !--------------------------------------------------------------------
       !
       ! Broadcast data from reaset and reabcs
       !
       !--------------------------------------------------------------------

       if( IMASTER ) then
          do kfl_desti_par = 1,npart
             call par_sendat(2_ip)
          end do
       else
          kfl_desti_par = 0
          call par_sendat(2_ip)
       end if

       !--------------------------------------------------------------------
       !
       ! Sequential reading: distribute element evenly
       !
       !--------------------------------------------------------------------

       if( .not. mpio_config%input%enabled .or. .not. mpio_config%input%parallel ) then
          !
          ! Normal sequential reading, data not scattered yet
          !
          call par_partition_scatter_from_master()

       else if( mpio_config%input%parallel ) then
          !
          ! Parallel reading, data already scattered
          !
          call par_partition_scatter_from_mpio()

       end if
       
    end if

    if( present(TIMINGS) ) then
       call PAR_BARRIER() 
       call cputim(time2)
       TIMINGS(1) = time2-time1
       time1 = time2
    end if

    !--------------------------------------------------------------------
    !
    ! Weights and center of gravity
    !
    !--------------------------------------------------------------------

    if( INOTMASTER .and. ( &
         &      kfl_partition_par == PAR_SFC          &
         & .or. kfl_partition_par == PAR_ORIENTED_BIN & 
         & .or. kfl_partition_par == PAR_ZOLTAN       ) ) then
       call memory_alloca(par_memor,'XCOOR','par_parallel_partitioning',xcoor,ndime,nelem)

       if(      kfl_partition_par == PAR_ZOLTAN  &
         & .or. kfl_partition_par == PAR_SFC     ) then
          call memory_alloca(par_memor,'LWEIG','par_parallel_partitioning',lweig,nelem)
          call par_partition_weights(lweig)
       else if( kfl_partition_par == PAR_ORIENTED_BIN ) then
          call memory_alloca(par_memor,'LWEIG','par_parallel_partitioning',lweig,nelem)
          call par_partition_weights(lweig)
       end if
       do ielem = 1,nelem
          ielty = ltype(ielem)
          pnode = nnode(ielty)
          do idime = 1,ndime
             xcoor(idime,ielem) = sum(coord(idime,lnods(1:pnode,ielem)))/real(pnode,rp)
          end do
       end do
    end if

    if( present(TIMINGS) ) then
       call PAR_BARRIER() 
       call cputim(time2)
       TIMINGS(2) = time2-time1
       time1 = time2
    end if

    !--------------------------------------------------------------------
    !
    ! Partitioning: LEPAR_PAR(IELEM)=subdomain of element IELEM
    !
    !--------------------------------------------------------------------

    if(      kfl_partition_par == PAR_SFC ) then
       call messages_live('PARTITIONING MESH WITH SFC')
    else if( kfl_partition_par == PAR_ORIENTED_BIN ) then
       call messages_live('PARTITIONING MESH WITH ORIENTED BIN')
    else if( kfl_partition_par == PAR_ZOLTAN ) then
       call messages_live('PARTITIONING MESH WITH ZOLTAN HSFC')
    else if( kfl_partition_par == PAR_USING_RANK ) then
       call messages_live('NO PARTITIONING: USING CURRENT MPI RANK AS SUBDOMAIN NUMBER')
    else if( kfl_partition_par == PAR_RANDOM ) then
       call messages_live('NO PARTITIONING: ASSIGNING SUBDOMAINS RANDOMLY')
    else if( kfl_partition_par < 0 ) then
       call messages_live('PARTITIONING TAKEN FROM A FIELD')
    end if

    if( INOTMASTER ) then
       call memory_alloca(par_memor,'LEPAR_PAR','partitioning_oriented_bin',lepar_par,nelem)
       if( kfl_partition_par < 0 ) then
          !
          ! Read from a field
          !
          ifiel = -kfl_partition_par
          if( kfl_field(2,ifiel) /= NELEM_TYPE           ) call runend('PAR_PARALLEL_PARTITIONING: FIELD SHOULD BE OF ELEMENT TYPE')
          if( memory_size(xfiel(ifiel) % a,2_ip) < nelem ) call runend('PAR_PARALLEL_PARTITIONING: WRONG FIELD DIMENSIONS')
          do ielem = 1,nelem
             lepar_par(ielem) = int(xfiel(ifiel) % a(1,ielem,1),ip)
          end do
          if( nelem > 0 ) then
             max_subd = maxval(lepar_par)
          else
             max_subd = 0
          end if
          call PAR_MAX(max_subd,PAR_COMM_MY_CODE_WM,INCLUDE_ROOT=.true.)
       else
          !
          ! Normal partitioning
          !
          npart_sav = npart
          npart_eff = npart-npart_empty_par
          call partitioning(&
               kfl_partition_par,nelem,npart_eff,lepar_par,lweig=lweig,&
               rweig=rweig,COMM4=PAR_COMM_MY_CODE_WM,ia=ia,ja=ja,lmask=lmask,&
               xcoor=xcoor,boxes_coarse=boxes_coarse_par,boxes_fine=boxes_fine_par,&
               direction=vect_partition_par,ndime_=ndime,lcorr_=LCORRECT)
          if( npart_sav /= npart ) &
               call runend('PAR_PARALLEL_PARTITIONING: CANNOT CHANGE NUMBER OF PARTITIONS!')
       end if
    end if

    if( present(TIMINGS) ) then
       call PAR_BARRIER() 
       call cputim(time2)
       TIMINGS(3) = time2-time1
       time1 = time2
    end if
    !
    ! Check partitioning
    !
    !call memory_alloca(par_memor,'LSUBD','partitioning_oriented_bin',lsubd,npart)
    !do ielem = 1,nelem
    !   ipart = lepar_par(ielem)
    !   lsubd(ipart) = lsubd(ipart) + 1
    !end do
    !call PAR_SUM(npart,lsubd)
    !if( IMASTER .and. minval(lsubd) == 0 ) call runend('PAR_PARALLEL_PARTITIONING: SOME PARTITIONS ARE EMPTY!')
    !call memory_deallo(par_memor,'LSUBD','partitioning_oriented_bin',lsubd)

    !--------------------------------------------------------------------
    !
    ! Redistribute
    !
    !--------------------------------------------------------------------

    if( kfl_partition_par /= PAR_USING_RANK ) then
       call messages_live('REDISTRIBUTE GEOMETRY AMONG SLAVES','START SECTION')
       call redistribution_comm_initialization(PAR_COMM_MY_CODE)
       call redistribution_domain(lepar_par,VERBOSE=.false.)
       call messages_live('REDISTRIBUTE GEOMETRY AMONG SLAVES','END SECTION')
    end if  

    if( present(TIMINGS) ) then
       call PAR_BARRIER() 
       call cputim(time2)
       TIMINGS(4) = time2-time1
       time1 = time2
    end if

    !--------------------------------------------------------------------
    !
    ! Construct communication arrays
    !
    !--------------------------------------------------------------------

    call messages_live('CONSTRUCT COMMUNICATION ARRAYS','START SECTION')
    call par_interface_exchange(PAR_COMM_MY_CODE_ARRAY(1))!,COMM_NAME='COMMD')!,TYPE_OF_COUPLING=SAME_COORDINATE)
    call messages_live('CONSTRUCT COMMUNICATION ARRAYS','END SECTION')
    !
    ! Some warnings...
    !
    npart_empty = 0
    if( nelem == 0 ) npart_empty = 1
    call PAR_SUM(npart_empty)
    if( npart_empty > 0 ) call messages_live('EMPTY SUBDOMAINS HAVE BEEN FOUND= '//trim(intost(npart_empty)),'WARNING')    
    !
    ! Deallocate memory
    !
    call memory_deallo(par_memor,'RWEIG'    ,'par_parallel_partitioning',rweig)
    call memory_deallo(par_memor,'LWEIG'    ,'par_parallel_partitioning',lweig)
    call memory_deallo(par_memor,'XCOOR'    ,'par_parallel_partitioning',xcoor)
    call memory_deallo(par_memor,'LEPAR_PAR','par_parallel_partitioning',lepar_par)

    call messages_live('PARALLEL PARTITIONING','END SECTION')

    if( present(TIMINGS) ) then
       call PAR_BARRIER()
       call cputim(time2)
       TIMINGS(5) = time2-time1
       time1 = time2
       call PAR_MAX(5_ip,TIMINGS)
    end if

  end subroutine par_parallel_partitioning

  !-----------------------------------------------------------------------
  !>
  !> @author  houzeaux
  !> @date    2018-10-13
  !> @brief   Compute weight
  !> @details Weights given to the partitioner
  !>
  !-----------------------------------------------------------------------

  subroutine par_partition_weights_rp(rweig,TOTAL_WEIGHT)

    real(rp),    pointer, optional, intent(inout) :: rweig(:)     !< Element weight
    real(rp),             optional, intent(out)   :: TOTAL_WEIGHT !< Total weight
    integer(ip)                                   :: ielem,imate
    integer(ip)                                   :: ielty,pgaus
    real(rp)                                      :: rr
    
    if( nelem <= 0 ) return
    if( present(TOTAL_WEIGHT) ) TOTAL_WEIGHT = 0.0_rp
    
    if( kfl_weigh_par == PAR_WEIGHT_GAUSS  ) then
       !
       ! Weight on Gauss points
       !
       do ielem = 1,nelem
          ielty = ltype(ielem)
          pgaus = domain_number_gauss_points(ielty)
          if( present(rweig) )        rweig(ielem) = real(pgaus,rp)
          if( present(TOTAL_WEIGHT) ) TOTAL_WEIGHT = TOTAL_WEIGHT + real(pgaus,rp)
       end do

    else if( kfl_weigh_par == PAR_WEIGHT_OFF ) then
       !
       ! No weight
       !
       do ielem = 1,nelem
          if( present(rweig) ) rweig(ielem) = 1.0_rp
       end do
       if( present(TOTAL_WEIGHT) ) TOTAL_WEIGHT = real(nelem,rp)

    else if( kfl_weigh_par == PAR_WEIGHT_ELEMENT ) then
       !
       ! Weight given by user
       !
       do ielem = 1,nelem
          ielty = abs(ltype(ielem))
          if( present(rweig) )        rweig(ielem) = weights_elements_par(ielty)
          if( present(TOTAL_WEIGHT) ) TOTAL_WEIGHT = TOTAL_WEIGHT + weights_elements_par(ielty)
       end do

    else if( kfl_weigh_par == PAR_WEIGHT_SQUARE ) then
       !
       ! Weight on square of Gauss points
       !
       do ielem = 1,nelem
          ielty = ltype(ielem)
          pgaus = domain_number_gauss_points(ielty)
          if( present(rweig) )        rweig(ielem) = real(pgaus,rp)**2
          if( present(TOTAL_WEIGHT) ) TOTAL_WEIGHT = TOTAL_WEIGHT + real(pgaus,rp)**2
       end do

    else if( kfl_weigh_par == PAR_WEIGHT_MATERIAL ) then
       !
       ! Weight given by material
       !
       if( maxval(lmate) > size(weights_materials_par,KIND=ip) ) then
          call runend('YOU MAVE TOO MANY MATERIALS TO ASSIGN WEIGHTS FOR PARTITION')
       end if
       do ielem = 1,nelem
          imate = abs(lmate(ielem))
          if( present(rweig) )        rweig(ielem) = weights_materials_par(imate)
          if( present(TOTAL_WEIGHT) ) TOTAL_WEIGHT = TOTAL_WEIGHT + weights_materials_par(imate)
       end do

    else if( kfl_weigh_par == PAR_WEIGHT_STUPID ) then
       !
       ! Weight given by material
       !
       do ielem = 1,nelem
          CALL RANDOM_NUMBER(rr)
          if( present(rweig) )        rweig(ielem) = rr
          if( present(TOTAL_WEIGHT) ) TOTAL_WEIGHT = TOTAL_WEIGHT + rr
       end do
       
    else if( kfl_weigh_par > 0 ) then
       !
       ! Weights read from a field
       !
       do ielem = 1,nelem
          if( present(rweig) )        rweig(ielem) = xfiel(kfl_weigh_par) % a(1,ielem,1)
          if( present(TOTAL_WEIGHT) ) TOTAL_WEIGHT = TOTAL_WEIGHT + xfiel(kfl_weigh_par) % a(1,ielem,1)
       end do
       
    end if

  end subroutine par_partition_weights_rp

  subroutine par_partition_weights_ip(lweig)

    integer(ip), pointer, intent(inout) :: lweig(:)     !< Weights
    integer(ip)                         :: ielem,ielty
    integer(ip)                         :: imate,pgaus

    if( nelem <= 0 ) return
    if( kfl_weigh_par == PAR_WEIGHT_GAUSS ) then
       !
       ! Weight on Gauss points
       !
       do ielem = 1,nelem
          ielty        = abs(ltype(ielem))
          pgaus        = domain_number_gauss_points(ielty)
          lweig(ielem) = pgaus
       end do

    else if( kfl_weigh_par == PAR_WEIGHT_OFF ) then
       !
       ! No weight
       !
       do ielem = 1,nelem
          ielty        = abs(ltype(ielem))
          lweig(ielem) = 1_ip
       end do
       
    else if( kfl_weigh_par == PAR_WEIGHT_ELEMENT ) then
       !
       ! Weight given by user
       !
       do ielem = 1,nelem
          ielty        = abs(ltype(ielem))
          lweig(ielem) = int(weights_elements_par(ielty),KIND=ip)
       end do
       
    else if( kfl_weigh_par == PAR_WEIGHT_SQUARE ) then
       !
       ! Weight on square of Gauss points
       !
       do ielem = 1,nelem
          ielty        = ltype(ielem)
          pgaus        = domain_number_gauss_points(ielty)
          lweig(ielem) = pgaus*pgaus
       end do

    else if( kfl_weigh_par == PAR_WEIGHT_MATERIAL ) then
       !
       ! Weight given by material
       !
       if( maxval(lmate) > size(weights_materials_par,KIND=ip) ) then
          call runend('YOU MAVE TOO MANY MATERIALS TO ASSIGN WEIGHTS FOR PARTITION')
       end if
       do ielem = 1,nelem
          imate = abs(lmate(ielem))
          lweig(ielem) = int(weights_materials_par(imate),KIND=ip)
       end do

    else if( kfl_weigh_par > 0 ) then
       !
       ! Weights read from a field
       !
       do ielem = 1,nelem
          ielty = abs(ltype(ielem))
          lweig(ielem) = int(xfiel(kfl_weigh_par) % a(1,ielem,1),KIND=ip)
       end do

    end if

  end subroutine par_partition_weights_ip

  subroutine par_partition_total_weight_ip(TOTAL_WEIGHT)
    integer(ip),          intent(out) :: TOTAL_WEIGHT
    integer(ip)                       :: ielem,pgaus
    integer(ip)                       :: ielty
    
    TOTAL_WEIGHT = 0

    if( kfl_weigh_par == 0 ) then
       !
       ! Weight on Gauss points
       !
       do ielem = 1,nelem
          ielty = ltype(ielem)
          TOTAL_WEIGHT = TOTAL_WEIGHT + domain_number_gauss_points(ielty)
       end do

    else if( kfl_weigh_par == -1 ) then
       !
       ! Weights on nodes
       !
       do ielem = 1,nelem
          ielty = abs(ltype(ielem))
          TOTAL_WEIGHT = TOTAL_WEIGHT + 1_ip
       end do

    else if( kfl_weigh_par == -2 ) then
       !
       ! Weight given by user
       !
       do ielem = 1,nelem
          ielty = abs(ltype(ielem))
          TOTAL_WEIGHT = TOTAL_WEIGHT + int(weights_elements_par(ielty),ip)
       end do

    else if( kfl_weigh_par == -3 ) then
       !
       ! Weight on Gauss points
       !
       do ielem = 1,nelem
          ielty        = ltype(ielem)
          pgaus        = domain_number_gauss_points(ielty)
          TOTAL_WEIGHT = TOTAL_WEIGHT + pgaus*pgaus
       end do

    else if( kfl_weigh_par > 0 ) then
       !
       ! Weights read from a field
       !
       do ielem = 1,nelem
          ielty = abs(ltype(ielem))
          TOTAL_WEIGHT = TOTAL_WEIGHT + int(xfiel(kfl_weigh_par) % a(1,ielem,1),ip)
       end do

    end if

  end subroutine par_partition_total_weight_ip

  !-----------------------------------------------------------------------
  !>
  !> @author  houzeaux
  !> @date    2018-10-10
  !> @brief   Scatter geometry from master
  !> @details When master reads the geometry sequentially, scatter it the
  !>          to the slaves by distributing evenly the elements
  !>
  !-----------------------------------------------------------------------

  subroutine par_partition_scatter_from_master()

    use def_domain, only : nboun
    use def_domain, only : npoin

    integer(ip), pointer    :: lpart(:)
    integer(ip)             :: ii,nsize

    nullify(lpart)


    call messages_live('SCATTER GEOMETRY FROM MASTER TO SLAVES','START SECTION')

    call redistribution_comm_initialization(PAR_COMM_WORLD)

    if( IMASTER ) then
       call memory_alloca(par_memor,'LNINV_LOC','savmsh',lninv_loc,npoin,'IDENTITY')
       call memory_alloca(par_memor,'LEINV_LOC','savmsh',leinv_loc,nelem,'IDENTITY')
       call memory_alloca(par_memor,'LBINV_LOC','savmsh',lbinv_loc,nboun,'IDENTITY')
       call memory_alloca(par_memor,'LPART','par_parallel_partitioning',lpart,nelem)
       nsize = nelem / npart
       do ii = 1,npart-1
          lpart((ii-1)*nsize+1:ii*nsize) = ii
       end do
       lpart((npart-1)*nsize+1:) = npart
    end if

    call redistribution_domain(lpart)
    call memory_deallo(par_memor,'LPART','par_parallel_partitioning',lpart)
    call messages_live('SCATTER GEOMETRY FROM MASTER TO SLAVES','END SECTION')

  end subroutine par_partition_scatter_from_master

  !-----------------------------------------------------------------------
  !>
  !> @author  houzeaux
  !> @date    2018-11-16
  !> @brief   Scatter from MPIO reading
  !> @details Paritions of elements, nodes and boudnaries have been
  !>          read using MPIO (NELEM_PART_P, NPOIN_PART_P, NPOIN_PART_P).
  !>          The blocks are the fixed size chunks reads.
  !>          These blocks NELEM_BLOK_P, NPOIN_BLOK_P, NPOIN_BLOK_P may
  !>          be different from what was actually read because NELEM/NPART
  !>          is not necessarily an integer!
  !>
  !>          \verbatim
  !>                                        not read
  !>                                           |
  !>          +------+------+------+------+--+---+
  !>          |      |      |      |      |  |000|
  !>          +------+------+------+------+--+---+
  !>          <------>                    <-->
  !>        NELEM_BLOCK_P              NELEM_PART_P
  !>
  !>
  !>          NELEM       NPOIN       NBOUN
  !>          ----------------------------------
  !>
  !>               <----- COMM_NPOIN
  !>               <----------------- COMM_NBOUN
  !>          ----------> COMM_NELEM (merge)
  !>                      ----------> COMM_NELBO
  !>
  !>          \endverbatim
  !>
  !-----------------------------------------------------------------------

  subroutine par_partition_scatter_from_mpio()

    use def_kintyp
    use def_elmtyp
    use def_master
    use def_domain
    use mod_parall
    use mod_redistribution
    use mod_communications
    use mod_memory
    use mod_maths
    use mod_mpio_par_readom, only : nelem_blok_p
    use mod_mpio_par_readom, only : nboun_blok_p
    use mod_mpio_par_readom, only : npoin_blok_p
    use mod_mpio_par_readom, only : nelem_part_p
    use mod_mpio_par_readom, only : nboun_part_p
    use mod_mpio_par_readom, only : npoin_part_p
    use mod_htable,          only : HtableMaxPrimeNumber
    use mod_htable,          only : hash_t
    use mod_htable,          only : htaini
    use mod_htable,          only : htares
    use mod_htable,          only : htaadd
    use mod_htable,          only : htades
    
    use def_parall
    integer(ip)          :: chunk_npoin
    integer(ip)          :: chunk_nelem
    integer(ip)          :: chunk_nboun
    type(comm_data_par)  :: comm_npoin
    type(comm_data_par)  :: comm_nboun
    integer(ip), pointer :: list_npoin(:)
    integer(ip), pointer :: list_nelem(:)
    integer(ip)          :: inode,ipoin,ielem
    integer(ip)          :: kelem,kpoin,kboun
    integer(ip)          :: iboun,inodb,rank
    integer(ip)          :: ifiel,max_npoin
    type(hash_chunkdist) :: htable_npoin
    type(hash_t)         :: ht

    call comm_npoin % init(COMM_NAME='COMM_NPOIN')
    call comm_nboun % init(COMM_NAME='COMM_NBOUN')

    nullify(list_npoin)
    nullify(list_nelem)
    !
    ! Dimensions
    !
    chunk_npoin = npoin_blok_p   ! Size of blocks of nodes
    chunk_nelem = nelem_blok_p   ! Size of blocks of elements
    chunk_nboun = nboun_blok_p   ! Size of blocks of boundaries

    nelem       = nelem_part_p   ! Number of elements read; final number of elements

    call messages_live('REDISTRIBUTE GEOMETRY AMONG SLAVES','START SECTION')

    if( kfl_bouel == 0 )  then
       call messages_live('COMPUTE LELBO AS IT WAS MISSING IN THE DOMAIN FILES')
    end if

    call PAR_BARRIER()
    call messages_live('NODAL ARRAYS')
    
    !-------------------------------------------------------------------
    !
    ! Get node communicator COMM_NPOIN
    !
    !-------------------------------------------------------------------

    if( INOTMASTER ) then
       max_npoin = nelem*mnode
       call memory_alloca(memor_dom,'LIST_NPOIN','memgeo',list_npoin,max_npoin)
       call htaini( ht, max_npoin, lidson=.false., AUTOMATIC_SIZE=.true.)
       call htares( ht, list_npoin )
       npoin = 0
       do ielem = 1,nelem
          do inode = 1,nnode(ltype(ielem))
             ipoin = lnods(inode,ielem)
             call htaadd( ht, ipoin)
             npoin = npoin + 1_ip
          end do
       end do
       npoin = ht % nelem
       if( npoin > 0 ) then
          call maths_heap_sort(2_ip,npoin ,list_npoin)
          call memory_resize(memor_dom,'LIST_NPOIN','memgeo',list_npoin,npoin)
       end if
       call htades( ht )

       call generate_comm_chunkdist(list_npoin,chunk_npoin,PAR_COMM_MY_CODE_WM,comm_npoin,dist_dom='SOURCE')
       call redistribution_generate_hash_chunkdist(list_npoin,chunk_npoin,comm_npoin,htable_npoin)
       call memory_deallo(memor_dom,'LIST_NPOIN','memgeo',list_npoin)
       !
       ! LNINV_LOC
       !
       if( .not. associated(lninv_loc) ) then
          call memory_alloca(memor_dom,'LNINV_LOC','memgeo',lninv_loc,npoin_part_p,'IDENTITY')
          rank = int(comm_npoin % RANK4,ip)
          kpoin = rank * chunk_npoin
          do ipoin = 1,npoin_part_p
             lninv_loc(ipoin) = lninv_loc(ipoin) + kpoin
          end do
       end if
       !
       ! Coordinates
       !
       if( kfl_bouel == 0 )  then
          !nullify(coord_sav)
          !call memory_copy(memor_dom,'COORD_SAVE','memgeo',coord,coord_sav,'DO_NOT_DEALLOCATE')
          call par_partition_lelbo2(chunk_nelem,chunk_npoin,comm_npoin) 
          kfl_bouel = 1
       end if
    end if
        
    !-------------------------------------------------------------------
    !
    ! NPOIN arrays redistribution
    !
    !-------------------------------------------------------------------
    
    if( INOTMASTER ) then       
       call redistribution_array(coord,    'NPOIN',MEMOR=memor_dom,VARIABLE_NAME='COORD'    ,COMM=comm_npoin,PERMUTE_RECV=.true.)
       call redistribution_array(lninv_loc,'NPOIN',MEMOR=memor_dom,VARIABLE_NAME='LNINV_LOC',COMM=comm_npoin,PERMUTE_RECV=.true.)
       call redistribution_array(lnoch,    'NPOIN',MEMOR=memor_dom,VARIABLE_NAME='LNOCH'    ,COMM=comm_npoin,PERMUTE_RECV=.true.)
       call redistribution_array(lmast,    'NPOIN',MEMOR=memor_dom,VARIABLE_NAME='LMAST'    ,COMM=comm_npoin,PERMUTE_RECV=.true.)
       call redistribution_array(kfl_codno,'NPOIN',MEMOR=memor_dom,VARIABLE_NAME='KFL_CODNO',COMM=comm_npoin,PERMUTE_RECV=.true.)
       call redistribution_array(lnset,    'NPOIN',MEMOR=memor_dom,VARIABLE_NAME='LNSET'    ,COMM=comm_npoin,PERMUTE_RECV=.true.)
       call redistribution_array(lgrou_dom,'NPOIN',MEMOR=memor_dom,VARIABLE_NAME='LGROU_DOM',COMM=comm_npoin,PERMUTE_RECV=.true.)
       do ifiel = 1,nfiel
          if( kfl_field(2,ifiel) == NPOIN_TYPE ) then
             if( kfl_field(4,ifiel) == 1 ) then
                call redistribution_array_RP_3(xfiel(ifiel) % a,'NPOIN',POSIT=2_ip,MEMOR=memor_dom,VARIABLE_NAME='XFIEL % A',COMM=comm_npoin,PERMUTE_RECV=.true.)
             else
                call redistribution_array_RP_3(xfiel(ifiel) % a,'NPOIN',POSIT=2_ip,MEMOR=memor_dom,VARIABLE_NAME='XFIEL % A',COMM=comm_npoin,PERMUTE_RECV=.true.,EXPLODE_DIM=3_ip)
             end if
          end if
       end do

       npoin = comm_npoin % lrecv_dim

    end if

    call PAR_BARRIER()
    call messages_live('BOUNDARY ARRAYS')

    !-------------------------------------------------------------------
    !
    ! Get node communicator COMM_NBOUN
    !
    !-------------------------------------------------------------------
    
    if( INOTMASTER ) then
       !
       ! Get boundary communicator COMM_NBOUN
       !
       call memory_alloca(memor_dom,'LIST_NELEM','memgeo',list_nelem,nboun_part_p)
       !
       ! Compute LELBO if it's not available...
       !
       !if( kfl_bouel == 0 )  then
       !   call messages_live('COMPUTE BOUNDARY/ELEMENT CONNECTIVITY AS IT WAS MISSING IN THE DOMAIN FILES')
       !   rank = int(PAR_MY_CODE_RANK_WM,ip)
       !   call par_partition_lelbo(chunk_nelem,chunk_npoin,htable_npoin)
       !   !call par_partition_lelbo_old(chunk_nelem)
       !   call memory_deallo(memor_dom,'COORD_SAVE','memgeo',coord_sav)
       !   kfl_bouel = 1
       !end if
       do iboun = 1,nboun_part_p
          ielem = lelbo(iboun)
          list_nelem(iboun) = ielem
       end do
       call generate_comm_chunkdist(list_nelem,chunk_nelem,PAR_COMM_MY_CODE_WM,comm_nboun,dist_dom='TARGET')
       call memory_deallo(memor_dom,'LIST_NELEM','memgeo',list_nelem)
       !
       ! LBINV_LOC
       !
       if( .not. associated(lbinv_loc) ) then
          call memory_alloca(memor_dom,'LBINV_LOC','memgeo',lbinv_loc,nboun_part_p,'IDENTITY')
          rank = int(comm_nboun % RANK4,ip)
          kboun = rank * chunk_nboun
          do iboun = 1,nboun_part_p
             lbinv_loc(iboun) = lbinv_loc(iboun) + kboun
          end do
       end if
    end if

    !-------------------------------------------------------------------
    !
    ! NBOUN arrays redistribution
    !
    !-------------------------------------------------------------------

    if( INOTMASTER ) then
       call redistribution_array(ltypb,    'NBOUN',MEMOR=memor_dom,VARIABLE_NAME='LTYPB',    COMM=comm_nboun,PERMUTE_RECV=.false.)
       call redistribution_array(lboch,    'NBOUN',MEMOR=memor_dom,VARIABLE_NAME='LBOCH',    COMM=comm_nboun,PERMUTE_RECV=.false.)
       call redistribution_array(lbset,    'NBOUN',MEMOR=memor_dom,VARIABLE_NAME='LBSET',    COMM=comm_nboun,PERMUTE_RECV=.false.)
       call redistribution_array(kfl_codbo,'NBOUN',MEMOR=memor_dom,VARIABLE_NAME='KFL_CODBO',COMM=comm_nboun,PERMUTE_RECV=.false.)
       call redistribution_array(lbinv_loc,'NBOUN',MEMOR=memor_dom,VARIABLE_NAME='LBINV_LOC',COMM=comm_nboun,PERMUTE_RECV=.false.)
       call redistribution_array(lelbo,    'NBOUN',MEMOR=memor_dom,VARIABLE_NAME='LELBO'    ,COMM=comm_nboun,PERMUTE_RECV=.false.)
       call redistribution_array(lnodb,    'NBOUN',MEMOR=memor_dom,VARIABLE_NAME='LNODB'    ,COMM=comm_nboun,PERMUTE_RECV=.false.)
       do ifiel = 1,nfiel
          if( kfl_field(2,ifiel) == NBOUN_TYPE ) then
             call redistribution_array_RP_3(xfiel(ifiel) % a,'NBOUN',POSIT=2_ip,MEMOR=memor_dom,VARIABLE_NAME='XFIEL % A',COMM=comm_nboun,PERMUTE_RECV=.false.)
          end if
       end do
    end if
    
    !-------------------------------------------------------------------
    !
    ! Permute LNODS, LNODB to local numbering and recompute LEINV_LOC
    !
    !-------------------------------------------------------------------
    
    if( INOTMASTER ) then

       do ielem = 1,nelem
          do inode = 1,nnode(ltype(ielem))
             ipoin = lnods(inode,ielem)
             lnods(inode,ielem) = redistribution_lid_chunkdist(ipoin,htable_npoin)
          end do
       end do
       !
       ! LEINV_LOC
       !
       if( .not. associated(leinv_loc) ) then
          call memory_alloca(memor_dom,'LNINV_LOC','memgeo',leinv_loc,chunk_nelem,'IDENTITY')
          rank = int(PAR_MY_CODE_RANK_WM,ip)
          kelem = rank * chunk_nelem
          do ielem = 1,chunk_nelem
             leinv_loc(ielem) = leinv_loc(ielem) + kelem
          end do
       end if
       !
       ! Permute LNODB
       !
       nboun = comm_nboun % lrecv_dim
       do iboun = 1,nboun
          ielem = comm_nboun % lrecv_perm(iboun)
          lelbo(iboun) = ielem
          do inodb = 1,nnode(ltypb(iboun))
             ipoin = lnodb(inodb,iboun)
             lnodb(inodb,iboun) = redistribution_lid_chunkdist(ipoin,htable_npoin)
          end do
       end do

    end if

    !-------------------------------------------------------------------
    !
    ! Deallocate
    !
    !-------------------------------------------------------------------

    call commd_npoin_from_mpio % init(COMM_NAME='COMMD_NPOIN_FROM_MPIO')
    call commd_nelem_from_mpio % init(COMM_NAME='COMMD_NELEM_FROM_MPIO')
    call commd_nboun_from_mpio % init(COMM_NAME='COMMD_NBOUN_FROM_MPIO')

    call commd_npoin_from_mpio % copy(comm_npoin,par_memor)
    call commd_nboun_from_mpio % copy(comm_nboun,par_memor)
    
    if( INOTMASTER ) then
       call comm_npoin % deallo(par_memor,INITIALIZE=.true.)
       call comm_nboun % deallo(par_memor,INITIALIZE=.true.)
       call redistribution_deallocate_hash(htable_npoin)
    end if

    call messages_live('REDISTRIBUTE GEOMETRY AMONG SLAVES','END SECTION')

  end subroutine par_partition_scatter_from_mpio

  subroutine par_partition_lelbo_old(chunk_nelem)

    use mod_communications,  only : PAR_MIN
    use mod_communications,  only : PAR_MAX
    use mod_communications,  only : PAR_ALLGATHER
    use mod_communications,  only : PAR_SEND_RECEIVE
    use mod_communications,  only : PAR_SEND_RECEIVE_IP
    use mod_communications,  only : PAR_ALLTOALL
    use mod_redistribution,  only : hash_chunkdist
    use mod_mpio_par_readom, only : nboun_part_p
    use mod_meshes,          only : meshes_list_boundary_elements
    use mod_htable,          only : hash_t
    use mod_htable,          only : htaini
    use mod_htable,          only : htares
    use mod_htable,          only : htaadd
    use mod_htable,          only : htades
    use mod_htable,          only : htalid
    use mod_memory,          only : memory_copy
    use def_domain

    integer(ip),          intent(in) :: chunk_nelem

    integer(ip)                      :: ipart,isend,irecv
    integer(ip)                      :: ipoin_min_lnods
    integer(ip)                      :: ipoin_max_lnods
    integer(ip)                      :: ipoin_min_lnodb
    integer(ip)                      :: ipoin_max_lnodb
    integer(ip), pointer             :: ipoin_min_lnods_gat(:)
    integer(ip), pointer             :: ipoin_max_lnods_gat(:)
    integer(ip), pointer             :: nelem_send(:)
    integer(ip), pointer             :: nelem_recv(:)
    type(i2p),   pointer             :: lnods_typ(:)
    type(i1p),   pointer             :: ltype_typ(:)
    type(i1p),   pointer             :: lperm_typ(:)

    integer(ip)                      :: ielem_boun,ielem,ipoin_loc
    integer(ip)                      :: ipoin,inode,npoin_loc,max_npoin
    integer(ip), pointer             :: lnods_boun(:,:)
    integer(ip), pointer             :: lnods_copy(:,:)
    integer(ip), pointer             :: ltype_boun(:)
    integer(ip)                      :: nelem_boun
    integer(ip), pointer             :: lperm(:)
    integer(ip)                      :: number_boundary_elements
    integer(ip), pointer             :: list_npoin(:)

    integer(ip)                      :: nelem_save
    integer(ip), pointer             :: lnods_save(:,:)
    integer(ip), pointer             :: ltype_save(:)
    type(hash_t)                     :: ht

    nullify(ipoin_min_lnods_gat)
    nullify(ipoin_max_lnods_gat)
    nullify(lnods_typ)
    nullify(ltype_typ)
    nullify(lperm_typ)
    nullify(nelem_send)
    nullify(nelem_recv)
    nullify(lnods_copy)
    nullify(lnods_boun)
    nullify(ltype_boun)
    nullify(lperm)
    nullify(list_npoin)

    !--------------------------------------------------------------------
    !
    ! Eliminate interior elements
    !
    !--------------------------------------------------------------------

    nelem_save =  nelem
    lnods_save => lnods
    ltype_save => ltype

    if( nelem > 0 ) then
       !
       ! Renumber mesh as LNODS is in global numbering
       !
       max_npoin = nelem * mnode
       call memory_copy  (par_memor,'LNODS'     ,'par_partition_lelbo',lnods,lnods_copy,'DO_NOT_DEALLOCATE')
       call memory_alloca(par_memor,'LIST_NPOIN','par_partition_lelbo',list_npoin,max_npoin)
       call htaini( ht, max_npoin, lidson=.true., AUTOMATIC_SIZE=.true.)
       call htares( ht, list_npoin )
       do ielem = 1,nelem
          do inode = 1,nnode(ltype(ielem))
             ipoin = lnods(inode,ielem)
             call htaadd( ht, ipoin)
          end do
       end do
       do ielem = 1,nelem
          do inode = 1,nnode(ltype(ielem))
             ipoin = lnods(inode,ielem)
             ipoin_loc = htalid(ht, ipoin)
             lnods_copy(inode,ielem) = ipoin_loc
          end do
       end do
       npoin_loc = ht % nelem
       call htades( ht )
       call memory_deallo(par_memor,'LIST_NPOIN','par_partition_lelbo',list_npoin)
       !
       ! Identify interior elements
       !
       call meshes_list_boundary_elements(&
            nelem,npoin,lnods_copy,ltype,&
            number_boundary_elements,&
            lperm,MEMORY_COUNTER=par_memor)
       nelem_boun = number_boundary_elements
       call memory_alloca(par_memor,'LNODS_BOUN','par_partition_lelbo',lnods_boun,mnode,nelem_boun)
       call memory_alloca(par_memor,'LTYPE_BOUN','par_partition_lelbo',ltype_boun,nelem_boun)
       do ielem_boun = 1,nelem_boun
          ielem                    = lperm(ielem_boun)
          lnods_boun(:,ielem_boun) = lnods(:,ielem) 
          ltype_boun(ielem_boun)   = ltype(ielem) 
       end do
       call memory_deallo(par_memor,'LNODS_COPY','par_partition_lelbo',lnods_copy)
       nelem =  nelem_boun
       lnods => lnods_boun
       ltype => ltype_boun
    end if

    !--------------------------------------------------------------------
    !
    ! Gather minimum node number in LNODS fo all partitions in IPOIN_MIN_LNODS_GAT and IPOIN_MAX_LNODS_GAT
    !
    !--------------------------------------------------------------------

    if( nelem > 0 ) then

       ipoin_min_lnods = minval(lnods,lnods/=0)
       ipoin_max_lnods = maxval(lnods)

    else

       ipoin_min_lnods =  huge(1_ip)
       ipoin_max_lnods = -huge(1_ip)
       
    end if

    call memory_alloca(par_memor,'IPOIN_MIN_GAT','par_partition_lelbo',ipoin_min_lnods_gat,npart,lboun=0_ip)
    call memory_alloca(par_memor,'IPOIN_MAX_GAT','par_partition_lelbo',ipoin_max_lnods_gat,npart,lboun=0_ip)

    call PAR_ALLGATHER(ipoin_min_lnods,ipoin_min_lnods_gat,PAR_COMM_IN=PAR_COMM_MY_CODE_WM)
    call PAR_ALLGATHER(ipoin_max_lnods,ipoin_max_lnods_gat,PAR_COMM_IN=PAR_COMM_MY_CODE_WM)

    call memory_alloca(par_memor,'LNODS_TYP','par_partition_lelbo',lnods_typ,npart,lboun=0_ip)
    call memory_alloca(par_memor,'LTYPE_TYP','par_partition_lelbo',ltype_typ,npart,lboun=0_ip)
    call memory_alloca(par_memor,'LTYPE_TYP','par_partition_lelbo',lperm_typ,npart,lboun=0_ip)

    call memory_alloca(par_memor,'NELEM_SEND','par_partition_lelbo',nelem_send,npart,lboun=0_ip)
    call memory_alloca(par_memor,'NELEM_RECV','par_partition_lelbo',nelem_recv,npart,lboun=0_ip)

    !--------------------------------------------------------------------
    !
    ! NELEM_SEND and NELEM_RECV
    !
    !--------------------------------------------------------------------

    do ipart = 0,npart-1
       isend = 0
       if( nboun_part_p > 0 ) then
          ipoin_min_lnodb = minval(lnodb,lnodb/=0)
          ipoin_max_lnodb = maxval(lnodb)
          if( ipoin_min_lnodb <= ipoin_max_lnods_gat(ipart) .and. ipoin_max_lnodb >= ipoin_min_lnods_gat(ipart) ) then
             isend = 1              ! Ask IPART to send his elements
          end if
       end if
       call PAR_SEND_RECEIVE(isend,irecv,dom_i=ipart,PAR_COMM_IN=PAR_COMM_MY_CODE_WM)
       if( irecv == 1 ) then
          nelem_send(ipart) = nelem ! I have to send my LNODS
       else
          nelem_send(ipart) = 0
       end if
    end do
    call PAR_ALLTOALL(1_ip,1_ip,nelem_send,nelem_recv,PAR_COMM_IN=PAR_COMM_MY_CODE_WM)

    !--------------------------------------------------------------------
    !
    ! Decide who will send LNODS
    !
    !--------------------------------------------------------------------

    do ipart = 0,npart-1
       if( nelem_send(ipart) /= 0 .or. nelem_recv(ipart) /= 0 ) then
          call memory_alloca(par_memor,'LNODS_TYP','par_partition_lelbo',lnods_typ(ipart) % l,mnode,max(1_ip,nelem_recv(ipart)))
          call memory_alloca(par_memor,'LTYPE_TYP','par_partition_lelbo',ltype_typ(ipart) % l,      max(1_ip,nelem_recv(ipart)))
          call memory_alloca(par_memor,'LPERM_TYP','par_partition_lelbo',lperm_typ(ipart) % l,      max(1_ip,nelem_recv(ipart)))
          call PAR_SEND_RECEIVE_IP(nelem_send(ipart)      ,nelem_recv(ipart)      ,ltype,ltype_typ(ipart) % l,dom_i=ipart,PAR_COMM_IN=PAR_COMM_MY_CODE_WM)
          call PAR_SEND_RECEIVE_IP(nelem_send(ipart)*mnode,nelem_recv(ipart)*mnode,lnods,lnods_typ(ipart) % l,dom_i=ipart,PAR_COMM_IN=PAR_COMM_MY_CODE_WM)
          call PAR_SEND_RECEIVE_IP(nelem_send(ipart)      ,nelem_recv(ipart)      ,lperm,lperm_typ(ipart) % l,dom_i=ipart,PAR_COMM_IN=PAR_COMM_MY_CODE_WM)
       end if
    end do

    !--------------------------------------------------------------------
    !
    ! Generate partial LELBO using global and local numbering
    !
    !--------------------------------------------------------------------

    if( nboun_part_p > 0 ) then
       do ipart = 0,npart-1
          if( nelem_recv(ipart) > 0 ) then
             call par_parallel_partitioning_boundary_elements(&
                  nelem_recv(ipart),lnods_typ(ipart) % l,ltype_typ(ipart)%l,&
                  nboun_part_p,lnodb,ltypb,&
                  lelbo,chunk_nelem,ipart,PERMUTATION=lperm_typ(ipart)%l)
          end if
       end do
    end if

    call memory_deallo(par_memor,'LNODS_TYP'     ,'par_partition_lelbo',lnods_typ)
    call memory_deallo(par_memor,'LTYPE_TYP'     ,'par_partition_lelbo',ltype_typ)
    call memory_deallo(par_memor,'LTYPE_TYP'     ,'par_partition_lelbo',lperm_typ)
    call memory_deallo(par_memor,'IPOIN_MIN_GAT' ,'par_partition_lelbo',ipoin_min_lnods_gat)
    call memory_deallo(par_memor,'IPOIN_MAX_GAT' ,'par_partition_lelbo',ipoin_max_lnods_gat)
    call memory_deallo(par_memor,'NELEM_SEND'    ,'par_partition_lelbo',nelem_send)
    call memory_deallo(par_memor,'NELEM_RECV'    ,'par_partition_lelbo',nelem_recv)
    call memory_deallo(par_memor,'LNODS_BOUN'    ,'par_partition_lelbo',lnods_boun)
    call memory_deallo(par_memor,'LTYPE_BOUN'    ,'par_partition_lelbo',ltype_boun)
    call memory_deallo(par_memor,'LPERM'         ,'par_partition_lelbo',lperm)
    !
    ! Recover original mesh
    !
    nelem =  nelem_save
    lnods => lnods_save
    ltype => ltype_save

  end subroutine par_partition_lelbo_old

  subroutine par_partition_lelbo(chunk_nelem,chunk_npoin,htable_npoin)

    use mod_communications,  only : PAR_MIN
    use mod_communications,  only : PAR_MAX
    use mod_communications,  only : PAR_ALLGATHER
    use mod_communications,  only : PAR_SEND_RECEIVE
    use mod_communications,  only : PAR_SEND_RECEIVE_IP
    use mod_communications,  only : PAR_ALLTOALL
    use mod_redistribution,  only : hash_chunkdist
    use mod_redistribution,  only : redistribution_array
    use mod_redistribution,  only : generate_comm_chunkdist
    use mod_redistribution,  only : redistribution_generate_hash_chunkdist
    use mod_redistribution,  only : redistribution_lid_chunkdist
    use mod_mpio_par_readom, only : nboun_part_p
    use mod_meshes,          only : meshes_list_boundary_elements
    use mod_htable,          only : hash_t
    use mod_htable,          only : htaini
    use mod_htable,          only : htares
    use mod_htable,          only : htaadd
    use mod_htable,          only : htades
    use mod_htable,          only : htalid
    use mod_memory,          only : memory_copy
    use mod_maths,           only : maths_heap_sort
    use def_domain

    integer(ip),          intent(in) :: chunk_nelem
    integer(ip),          intent(in) :: chunk_npoin
    type(hash_chunkdist), intent(in) :: htable_npoin

    integer(ip)                      :: ipart,isend,irecv
    integer(ip)                      :: ipoin_min_lnods
    integer(ip)                      :: ipoin_max_lnods
    integer(ip)                      :: ipoin_min_lnodb
    integer(ip)                      :: ipoin_max_lnodb
    integer(ip)                      :: ielem_boun,ielem
    integer(ip)                      :: iboun
!    integer(ip)                      :: icoun1,icoun2,icoun3
    integer(ip)                      :: ipoin,inode,npoin_loc,max_npoin
    integer(ip)                      :: nelem_boun,idime
    integer(ip)                      :: number_boundary_elements
 
    integer(ip)                      :: nelem_save
    type(comm_data_par)              :: comm_npoin

    real(rp)                         :: comin(3),comax(3)
    real(rp)                         :: comin_loc(3),comax_loc(3)
    type(hash_t)                     :: ht
    type(hash_chunkdist)             :: htable_npoin_loc

    integer(ip), pointer             :: ipoin_min_lnods_gat(:)
    integer(ip), pointer             :: ipoin_max_lnods_gat(:)
    integer(ip), pointer             :: lnods_rcv(:,:)
    integer(ip), pointer             :: ltype_rcv(:)
    integer(ip), pointer             :: lperm_rcv(:)
    integer(ip), pointer             :: nelem_send(:)
    integer(ip), pointer             :: nelem_recv(:)
    integer(ip), pointer             :: lperm(:)
    integer(ip), pointer             :: list_npoin(:)
    integer(ip), pointer             :: lnods_boun(:,:)
    integer(ip), pointer             :: lnods_copy(:,:)
    integer(ip), pointer             :: ltype_boun(:)
    integer(ip), pointer             :: lnods_save(:,:)
    integer(ip), pointer             :: ltype_save(:)
    logical(lg), pointer             :: intersection(:)
    logical(lg), pointer             :: consider_node(:)
    real(rp),    pointer             :: comin_gat(:,:)
    real(rp),    pointer             :: comax_gat(:,:)
    
    nullify(ipoin_min_lnods_gat)
    nullify(ipoin_max_lnods_gat)
    nullify(lnods_rcv)
    nullify(ltype_rcv)
    nullify(lperm_rcv)
    nullify(nelem_send)
    nullify(nelem_recv)
    nullify(lperm)
    nullify(list_npoin)
    nullify(lnods_boun)
    nullify(lnods_copy)
    nullify(ltype_boun)
    nullify(lnods_save)
    nullify(ltype_save)
    nullify(intersection)
    nullify(consider_node)
    nullify(comin_gat)
    nullify(comax_gat)

    !--------------------------------------------------------------------
    !
    ! Eliminate interior elements
    !
    !--------------------------------------------------------------------

    nelem_save =  nelem
    lnods_save => lnods
    ltype_save => ltype
    comin      =  huge(1.0_rp) * 0.1_rp
    comax      = -huge(1.0_rp) * 0.1_rp
    comin_loc  =  huge(1.0_rp) * 0.1_rp
    comax_loc  = -huge(1.0_rp) * 0.1_rp       

    call memory_alloca(par_memor,'COMIN_GAT'    ,'par_partition_lelbo',comin_gat,3_ip,npart,lboun1=1_ip,lboun2=0_ip)
    call memory_alloca(par_memor,'COMAX_GAT'    ,'par_partition_lelbo',comax_gat,3_ip,npart,lboun1=1_ip,lboun2=0_ip)
    call memory_alloca(par_memor,'CONSIDER_NODE','par_partition_lelbo',consider_node,npoin)

    if( nelem > 0 ) then
       !
       ! Renumber mesh as LNODS is in global numbering
       !
       call memory_copy(par_memor,'LNODS_COPY','par_partition_lelbo',lnods,lnods_copy,'DO_NOT_DEALLOCATE')
       do ielem = 1,nelem
          do inode = 1,nnode(ltype(ielem))
             ipoin = lnods(inode,ielem)
             lnods_copy(inode,ielem) = redistribution_lid_chunkdist(ipoin,htable_npoin)
          end do
       end do
       !
       ! Identify interior elements and mark nodes
       !
       call meshes_list_boundary_elements(&
            nelem,npoin,lnods_copy,ltype,&
            number_boundary_elements,&
            lperm,MEMORY_COUNTER=par_memor)
       nelem_boun = number_boundary_elements
       call memory_alloca(par_memor,'LNODS_BOUN','par_partition_lelbo',lnods_boun,mnode,nelem_boun)
       call memory_alloca(par_memor,'LTYPE_BOUN','par_partition_lelbo',ltype_boun,nelem_boun)
       do ielem_boun = 1,nelem_boun
          ielem                    = lperm(ielem_boun)
          lnods_boun(:,ielem_boun) = lnods(:,ielem) 
          ltype_boun(ielem_boun)   = ltype(ielem) 
          do inode = 1,nnode(ltype(ielem))
             ipoin = lnods_copy(inode,ielem)
             consider_node(ipoin) = .true.
          end do
       end do    
       call memory_deallo(par_memor,'LNODS_COPY','par_partition_lelbo',lnods_copy)
       nelem =  nelem_boun
       lnods => lnods_boun
       ltype => ltype_boun
       !
       ! Bounding box of 
       !
       if( npoin > 0 ) then
          do idime = 1,ndime
             comin(idime) = minval(coord(idime,:),consider_node)
             comax(idime) = maxval(coord(idime,:),consider_node)
          end do
          comin = comin - 0.01_rp*(comax-comin)
          comax = comax + 0.01_rp*(comax-comin)
       end if
       
    end if

    call PAR_ALLGATHER(comin,comin_gat,3_4,PAR_COMM_IN=PAR_COMM_MY_CODE_WM)
    call PAR_ALLGATHER(comax,comax_gat,3_4,PAR_COMM_IN=PAR_COMM_MY_CODE_WM)
    call memory_deallo(par_memor,'CONSIDER_NODE','par_partition_lelbo',consider_node)

    !--------------------------------------------------------------------
    !
    ! Get NPOIN communicator
    !
    !--------------------------------------------------------------------
    !
    ! Mark required nodes
    !
    max_npoin = nboun_part_p*mnodb
    call memory_alloca(par_memor,'LIST_NPOIN','memgeo',list_npoin,max_npoin)
    call htaini( ht, 2_ip*max_npoin, lidson=.false., AUTOMATIC_SIZE=.true.)
    call htares( ht, list_npoin )
    do iboun = 1,nboun_part_p
       do inode = 1,nnode(ltypb(iboun))
          ipoin = lnodb(inode,iboun)
          call htaadd( ht, ipoin)
       end do
    end do
    npoin_loc = ht % nelem
    !
    ! List of required nodes
    !
    if( npoin_loc > 0 ) then
       call maths_heap_sort(2_ip,npoin_loc,list_npoin)
       call memory_resize(par_memor,'LIST_NPOIN','memgeo',list_npoin,npoin_loc)
    else
       call memory_deallo(par_memor,'LIST_NPOIN','memgeo',list_npoin)          
    end if
    call htades(ht)
    call generate_comm_chunkdist(list_npoin,chunk_npoin,PAR_COMM_MY_CODE_WM,comm_npoin,dist_dom='SOURCE')
    call redistribution_generate_hash_chunkdist(list_npoin,chunk_npoin,comm_npoin,htable_npoin_loc)
    call memory_deallo(par_memor,'LIST_NPOIN','memgeo',list_npoin)
    call redistribution_array(coord_sav,'NPOIN',MEMOR=par_memor,VARIABLE_NAME='COORD_SAV',COMM=comm_npoin,PERMUTE_RECV=.true.)
    call comm_npoin % deallo(par_memor,INITIALIZE=.true.)
    if( npoin_loc > 0 ) then
       do idime = 1,ndime
          comin_loc(idime) = minval(coord_sav(idime,:))
          comax_loc(idime) = maxval(coord_sav(idime,:)) 
       end do
    end if
    
    !--------------------------------------------------------------------
    !
    ! Intersection boundary bounding box with LNODS bounding box
    !
    !--------------------------------------------------------------------

    call memory_alloca(par_memor,'INTERSECTION','par_partition_lelbo',intersection,npart,lboun=0_ip)
    intersection = .true.
    do ipart = 0,npart-1
       do idime = 1,ndime
          if(    comin_loc(idime) > comax_gat(idime,ipart) .or. &
               & comax_loc(idime) < comin_gat(idime,ipart)  ) then
             intersection(ipart) = .false.
          end if
       end do
    end do
    call memory_deallo(par_memor,'COMIN_GAT','par_partition_lelbo',comin_gat)
    call memory_deallo(par_memor,'COMAX_GAT','par_partition_lelbo',comax_gat)

    !--------------------------------------------------------------------
    !
    ! Gather minimum node number in LNODS fo all partitions in IPOIN_MIN_LNODS_GAT and IPOIN_MAX_LNODS_GAT
    !
    !--------------------------------------------------------------------

    if( nelem > 0 ) then

       ipoin_min_lnods = minval(lnods,lnods/=0)
       ipoin_max_lnods = maxval(lnods)

    else

       ipoin_min_lnods =  huge(1_ip)
       ipoin_max_lnods = -huge(1_ip)

    end if

    call memory_alloca(par_memor,'IPOIN_MIN_GAT','par_partition_lelbo',ipoin_min_lnods_gat,npart,lboun=0_ip)
    call memory_alloca(par_memor,'IPOIN_MAX_GAT','par_partition_lelbo',ipoin_max_lnods_gat,npart,lboun=0_ip)

    call PAR_ALLGATHER(ipoin_min_lnods,ipoin_min_lnods_gat,PAR_COMM_IN=PAR_COMM_MY_CODE_WM)
    call PAR_ALLGATHER(ipoin_max_lnods,ipoin_max_lnods_gat,PAR_COMM_IN=PAR_COMM_MY_CODE_WM)

    !--------------------------------------------------------------------
    !
    ! NELEM_SEND and NELEM_RECV
    !
    !--------------------------------------------------------------------

    call memory_alloca(par_memor,'NELEM_SEND','par_partition_lelbo',nelem_send,npart,lboun=0_ip)
    call memory_alloca(par_memor,'NELEM_RECV','par_partition_lelbo',nelem_recv,npart,lboun=0_ip)

    do ipart = 0,npart-1
       isend = 0
       if( nboun_part_p > 0 ) then
          ipoin_min_lnodb = minval(lnodb,lnodb/=0)
          ipoin_max_lnodb = maxval(lnodb)
          if( ipoin_min_lnodb <= ipoin_max_lnods_gat(ipart) .and. ipoin_max_lnodb >= ipoin_min_lnods_gat(ipart) ) then
             isend = 1              ! Ask IPART to send his elements
          end if
          if( .not. intersection(ipart) ) isend = 0
       end if
       call PAR_SEND_RECEIVE(isend,irecv,dom_i=ipart,PAR_COMM_IN=PAR_COMM_MY_CODE_WM)
       if( irecv == 1 ) then
          nelem_send(ipart) = nelem ! I have to send my LNODS
       else
          nelem_send(ipart) = 0
       end if
    end do
    call memory_deallo(par_memor,'IPOIN_MIN_GAT' ,'par_partition_lelbo',ipoin_min_lnods_gat)
    call memory_deallo(par_memor,'IPOIN_MAX_GAT' ,'par_partition_lelbo',ipoin_max_lnods_gat)
    call PAR_ALLTOALL(1_ip,1_ip,nelem_send,nelem_recv,PAR_COMM_IN=PAR_COMM_MY_CODE_WM)

    !--------------------------------------------------------------------
    !
    ! Decide who will send LNODS
    !
    !--------------------------------------------------------------------

    !icoun1=0
    !icoun2=0
    !icoun3=0
    !do ipart = 0,npart-1
    !   if( intersection(ipart) ) then
    !      if( nelem_send(ipart) /= 0 .or. nelem_recv(ipart) /= 0 ) then
    !         icoun1 = icoun1 + 1
    !         if( nboun_part_p > 0 .and. nelem_recv(ipart) > 0 ) then
    !            icoun2 = icoun2 + 1
    !         end if
    !      end if
    !   end if
    !   if( intersection(ipart) ) icoun3 = icoun3 + 1
    !end do

    do ipart = 0,npart-1
       if( nelem_send(ipart) /= 0 .or. nelem_recv(ipart) /= 0 ) then
          call memory_alloca(par_memor,'LNODS_RCV','par_partition_lelbo',lnods_rcv,mnode,max(1_ip,nelem_recv(ipart)))
          call memory_alloca(par_memor,'LTYPE_RCV','par_partition_lelbo',ltype_rcv,      max(1_ip,nelem_recv(ipart)))
          call memory_alloca(par_memor,'LPERM_RCV','par_partition_lelbo',lperm_rcv,      max(1_ip,nelem_recv(ipart)))
          call PAR_SEND_RECEIVE_IP(nelem_send(ipart)      ,nelem_recv(ipart)      ,ltype,ltype_rcv,dom_i=ipart,PAR_COMM_IN=PAR_COMM_MY_CODE_WM)
          call PAR_SEND_RECEIVE_IP(nelem_send(ipart)*mnode,nelem_recv(ipart)*mnode,lnods,lnods_rcv,dom_i=ipart,PAR_COMM_IN=PAR_COMM_MY_CODE_WM)
          call PAR_SEND_RECEIVE_IP(nelem_send(ipart)      ,nelem_recv(ipart)      ,lperm,lperm_rcv,dom_i=ipart,PAR_COMM_IN=PAR_COMM_MY_CODE_WM)
          if( nboun_part_p > 0 .and. nelem_recv(ipart) > 0 ) then
             call par_parallel_partitioning_boundary_elements(&
                  nelem_recv(ipart),lnods_rcv,ltype_rcv,&
                  nboun_part_p,lnodb,ltypb,&
                  lelbo,chunk_nelem,ipart,PERMUTATION=lperm_rcv)
          end if
          call memory_deallo(par_memor,'LNODS_RCV','par_partition_lelbo',lnods_rcv)
          call memory_deallo(par_memor,'LTYPE_RCV','par_partition_lelbo',ltype_rcv)
          call memory_deallo(par_memor,'LPERM_RCV','par_partition_lelbo',lperm_rcv)             
       end if
    end do

    call memory_deallo(par_memor,'INTERSECTION'  ,'par_partition_lelbo',intersection)
    call memory_deallo(par_memor,'NELEM_SEND'    ,'par_partition_lelbo',nelem_send)
    call memory_deallo(par_memor,'NELEM_RECV'    ,'par_partition_lelbo',nelem_recv)
    call memory_deallo(par_memor,'LNODS_BOUN'    ,'par_partition_lelbo',lnods_boun)
    call memory_deallo(par_memor,'LTYPE_BOUN'    ,'par_partition_lelbo',ltype_boun)
    call memory_deallo(par_memor,'LPERM'         ,'par_partition_lelbo',lperm)
    !
    ! Recover original mesh
    !
    nelem =  nelem_save
    lnods => lnods_save
    ltype => ltype_save

  end subroutine par_partition_lelbo

  subroutine par_parallel_partitioning_boundary_elements(&
       nelem_loc,lnods_loc,ltype_loc,&
       nboun_loc,lnodb_loc,ltypb_loc,&
       lelbo_loc,chunk_nelem,rank,PERMUTATION)

    use def_kintyp,         only : ip,i1p
    use mod_memory,         only : memory_alloca
    use mod_memory,         only : memory_deallo
    use def_domain,         only : mnode
    use mod_parall,         only : par_memor
    use mod_maths,          only : maths_heap_sort
    use mod_messages,       only : messages_live
    use mod_redistribution, only : redistribution_lid_chunkdist
    use mod_htable,         only : hash_t,htaini,htaadd
    use mod_htable,         only : htalid,htades
    use mod_htable,         only : HTableMaxPrimeNumber

    implicit none

    integer(ip),          intent(in)             :: nelem_loc
    integer(ip), pointer, intent(in)             :: lnods_loc(:,:)
    integer(ip), pointer, intent(in)             :: ltype_loc(:)
    integer(ip),          intent(in)             :: nboun_loc
    integer(ip), pointer, intent(in)             :: lnodb_loc(:,:)
    integer(ip), pointer, intent(in)             :: ltypb_loc(:)
    integer(ip), pointer, intent(inout)          :: lelbo_loc(:)
    integer(ip), pointer, intent(in),   optional :: PERMUTATION(:)

    integer(ip),          intent(in)             :: chunk_nelem
    integer(ip),          intent(in)             :: rank

    integer(ip)                                  :: ipoin,inodb,pnode,pnodb
    integer(ip)                                  :: jpoin,inode,ielem,ii
    integer(ip)                                  :: iboun,kelem,nlist
    integer(ip)                                  :: npoin_loc
    integer(ip)                                  :: ielem_offset
    type(i1p),   pointer                         :: lpoel(:)
    integer(ip), pointer                         :: lnoel(:)
    integer(ip), pointer                         :: liste(:)
    integer(ip), pointer                         :: listp(:)
    type(hash_t)                                 :: htable_npoin

    nullify(lpoel)
    nullify(lnoel)
    nullify(liste)
    nullify(listp)

    ielem_offset    = rank * chunk_nelem
    !
    ! Allocate memory
    !
    call messages_live('COMPUTE BOUNDARY/ELEMENT CONNECTIVITY')
    !
    ! Generate htable
    !
    npoin_loc = nelem_loc*mnode

    call memory_alloca(par_memor,'LNOEL','lbouel',lnoel,npoin_loc)

    call htaini(htable_npoin,nelem_loc*mnode,lidson=.true.,AUTOMATIC_SIZE=.true.)
    !
    ! Compute LNOEL and local ID
    !
    do ielem = 1,nelem_loc
       do inode = 1,nnode(ltype_loc(ielem))
          ipoin = lnods_loc(inode,ielem)
          call htaadd(htable_npoin,ipoin,jpoin)
          lnoel(jpoin) = lnoel(jpoin) + 1
       end do
    end do
    !
    ! Maxmimum number of elements connected to nodes
    !
    npoin_loc = htable_npoin % nelem
    call memory_alloca(par_memor,'LPOEL','lbouel',lpoel,npoin_loc)
    !
    ! Construct the list of elements connected to nodes LPOEL(:) % L
    !
    do ipoin = 1,npoin_loc
       call memory_alloca(par_memor,'LPOEL % L','lbouel',lpoel(ipoin) % l,lnoel(ipoin))
       lnoel(ipoin) = 0
    end do

    do ielem = 1,nelem_loc
       pnode = nnode(ltype_loc(ielem))
       do inode = 1,pnode
          ipoin = lnods_loc(inode,ielem)
          jpoin                          = htalid(htable_npoin,ipoin)
          lnoel(jpoin)                   = lnoel(jpoin) + 1
          lpoel(jpoin) % l(lnoel(jpoin)) = ielem
       end do
    end do
    call memory_alloca(par_memor,'LISTE','lbouel',liste,nelem_loc)
    call memory_alloca(par_memor,'LISTP','lbouel',listp,nelem_loc)
    nlist = 0

    !$*OMP PARALLEL   DO                                           &
    !$*OMP SCHEDULE ( STATIC )                                     &
    !$*OMP DEFAULT  ( NONE )                                       &
    !$*OMP SHARED   ( nboun_loc,ltypb_loc,lnodb_loc,lnoel,lpoel,   &
    !$*OMP            lelbo_loc                                  ) &
    !$*OMP PRIVATE  ( iboun,inodb,ipoin,kelem,jpoin,pnodb,ielem  )

    do iboun = 1,nboun_loc

       do ii = 1,nlist
          liste(listp(ii)) = 0_ip
       end do
       nlist = 0

       if( lelbo_loc(iboun) == 0 ) then
          pnodb = nnode(ltypb_loc(iboun))
          do inodb = 1,pnodb
             ipoin = lnodb_loc(inodb,iboun)
             jpoin = htalid(htable_npoin,ipoin)
             if( jpoin /= 0 ) then
                do kelem = 1,lnoel(jpoin)
                   ielem        = lpoel(jpoin) % l(kelem)
                   liste(ielem) = liste(ielem) + 1
                   nlist         = nlist + 1
                   listp(nlist)  = ielem
                   if( liste(ielem) == pnodb ) then
                      if( present(PERMUTATION) ) then
                         lelbo_loc(iboun) = PERMUTATION(ielem) + ielem_offset
                      else
                         lelbo_loc(iboun) = ielem + ielem_offset
                      end if
                      go to 1
                   end if
                end do
             else
                goto 1
             end if
          end do
       end if

1      continue

    end do

    !$*OMP END PARALLEL DO
    !
    ! Deallocate memory
    !
    call memory_deallo(par_memor,'LNOEL','lbouel',lnoel)
    call memory_deallo(par_memor,'LPOEL','lbouel',lpoel)
    call memory_deallo(par_memor,'LISTE','lbouel',liste)
    call memory_deallo(par_memor,'LISTP','lbouel',liste)
    call htades(htable_npoin)

  end subroutine par_parallel_partitioning_boundary_elements

  
subroutine par_partition_lelbo2(chunk_nelem,chunk_npoin,comm_npoin)

   use def_domain,          only  : mnodb 
   use def_domain,          only  : lnodb
   use def_domain,          only  : ltypb
   use mod_htable,          only  : hash_t
   use mod_htable,          only  : htaini
   use mod_htable,          only  : htares
   use mod_htable,          only  : htaadd
   use mod_htable,          only  : htades
   use mod_htable,          only  : htalid
   use mod_redistribution,  only  : generate_comm_chunkdist
   use mod_redistribution,  only  : generate_comm_sizes
   use mod_communications,  only  : PAR_SEND_RECEIVE_TO_ALL
   use mod_parall,          only  : PAR_TRANSPOSE_COMMUNICATION_ARRAY
   use mod_mpio_par_readom, only  : nboun_part_p
   use mod_mpio_par_readom, only  : npoin_part_p
   use def_domain,          only  : lelbo ! <- output

   implicit none

   integer(ip),          intent(in) :: chunk_nelem
   integer(ip),          intent(in) :: chunk_npoin
   type(comm_data_par)              :: comm_npoin

   integer(ip),          pointer    :: list_bpoin(:)     ! boundary points
   integer(ip),          pointer    :: list_nboup(:)     ! #(boundary elements containing point)
   integer(ip),          pointer    :: list_boupo(:,:)   ! list of boundary elements containing point
   integer(ip),          pointer    :: lnrank(:)         ! #(ranks to which each point is sent by comm_npoin)
   integer(ip),          pointer    :: lranks(:,:)       ! list of ranks to which each point is sent
   integer(ip),          pointer    :: buff0(:,:)        
   integer(ip),          pointer    :: lnodb_buf(:,:)
   integer(ip),          pointer    :: ltypb_buf(:)
   integer(ip),          pointer    :: lelbo_val(:)
   integer(ip),          pointer    :: lelbo_buf(:)
   integer(4),           pointer    :: lnsend(:)
   integer(ip),          pointer    :: all2send(:)
   integer(ip),          pointer    :: nperm(:)
   type(i1p),            pointer    :: send_perm(:)
   integer(ip)                      :: mranks, isend, irecv
   integer(ip)                      :: iboun,  inodb, ipoin
   integer(ip)                      :: icont,  ineig, pos
   integer(ip)                      :: nsendp, iperm, ibuff
   integer(ip)                      :: irank, nboun_buf, ipart
   integer(4)                       :: PAR_CODE_SIZE_WM4, ipart4
   type(hash_t)                     :: ht
   type(comm_data_par)              :: comm0           ! communicator from nodes to boundaries
   type(comm_data_par)              :: comm1           ! communicator from boundaries to elements
   type(comm_data_par)              :: comm2           ! transpose of comm1
   !
   ! Initializations
   !
   nullify(list_bpoin,list_nboup,list_boupo,lelbo_buf)
   nullify(lnsend,all2send,nperm,send_perm,lnrank)
   nullify(lranks,lnodb_buf,ltypb_buf,buff0,lelbo_val)
  
   call comm0 % init(COMM_NAME='COMM0')
   call comm1 % init(COMM_NAME='COMM1')
   call comm2 % init(COMM_NAME='COMM2')
   call PAR_COMM_RANK_AND_SIZE(PAR_COMM_MY_CODE_WM,ipart4,PAR_CODE_SIZE_WM4)
   ipart = int(ipart4,ip)
   !
   ! Generate COMM0 and LIST_BOUPO
   !
   call memory_alloca(par_memor,'LIST_BPOIN',vacal,list_bpoin,nboun_part_p*mnodb)
   call memory_alloca(par_memor,'LIST_NBOUP',vacal,list_nboup,nboun_part_p*mnodb)

   call htaini( ht, 3*chunk_npoin, lidson=.true., AUTOMATIC_SIZE=.true.)
   call htares( ht, list_bpoin )
   do iboun = 1,nboun_part_p
      do inodb = 1,nnode(ltypb(iboun))
         ipoin = lnodb(inodb,iboun)
         call htaadd( ht, ipoin,lid=pos)
         list_nboup(pos)                 = list_nboup(pos) + 1_ip 
      end do
   end do

   if( nboun_part_p> 0 ) then
      pos = maxval(list_nboup)
      list_nboup = 0
   end if
   call memory_alloca(par_memor,'LIST_BOUPO',vacal,list_boupo,pos,ht % nelem)
   do iboun = 1,nboun_part_p
      do inodb = 1,nnode(ltypb(iboun))
         ipoin                           = lnodb(inodb,iboun)
         pos                             = htalid(ht,ipoin)
         list_nboup(pos)                 = list_nboup(pos) + 1_ip 
         list_boupo(list_nboup(pos),pos) = iboun
      end do
   end do
   
   call memory_resize(par_memor,'LIST_BPOIN',vacal,list_bpoin,ht % nelem)
   call generate_comm_chunkdist(list_bpoin,chunk_npoin,PAR_COMM_MY_CODE_WM,comm0,dist_dom='SOURCE')
   call memory_deallo(par_memor,'LIST_BPOIN',vacal,list_bpoin)

   !
   ! Evaluate field LRANKS
   !
   call memory_alloca(par_memor,'LNRANK',vacal,lnrank,npoin_part_p)
   do isend = 1,comm_npoin % lsend_dim
      ipoin = comm_npoin % lsend_perm(isend)
      lnrank( ipoin ) = lnrank( ipoin ) + 1_ip
   enddo
   mranks = maxval(lnrank)
   
   call PAR_MAX(mranks,PAR_COMM_MY_CODE_WM,INCLUDE_ROOT=.true.)

   call memory_alloca(par_memor,'LRANKS',vacal,lranks,mranks,npoin_part_p)
   if( npoin_part_p > 0 ) then
      lnrank=0_ip
      lranks=-1_ip
   end if
   icont = 1_ip 
   do ineig = 1, comm_npoin % nneig
      do isend = comm_npoin % lsend_size(ineig), comm_npoin % lsend_size(ineig+1)-1
         ipoin = comm_npoin % lsend_perm(icont)
         lnrank( ipoin ) = lnrank( ipoin ) + 1_ip
         lranks(lnrank(ipoin),ipoin) = comm_npoin % neights(ineig)
         icont = icont + 1_ip
      enddo
   enddo
   call memory_deallo(par_memor,'LNRANK',vacal,lnrank)
   !
   ! Redistribute LRANKS
   !
   call memory_alloca(par_memor,'BUFF0',vacal,buff0,mranks,comm0 % lrecv_dim)
   call PAR_SEND_RECEIVE_TO_ALL(lranks,buff0,comm0,'ASYNCHRONOUS',PERMUTE_SEND=.true.,PERMUTE_RECV=.true.)
   call memory_deallo(par_memor,'LRANKS',vacal,lranks)
   !
   ! Generate COMM1
   !

   !  I. Evaluate lnsend
   call memory_alloca(par_memor,'LNSEND',vacal,lnsend,PAR_CODE_SIZE_WM4,lboun=0_4)
   do irecv = 1, comm0 % lrecv_dim
      do ibuff= 1, mranks
         if(buff0(ibuff,irecv) /= -1) then
            irank         = buff0(ibuff,irecv)
            lnsend(irank) = lnsend(irank) + int(list_nboup(irecv),4) 
         endif
      enddo
   enddo

   ! II. Evaluate communication sizes, dims etc
   call comm1 % deallo(par_memor,INITIALIZE=.true.)
   call comm1 % init  (COMM_NAME='COMM1')
   call generate_comm_sizes(lnsend,nsendp,comm1,PAR_COMM_MY_CODE_WM)
   
   !III. Evaluate lsend_perm
   call memory_alloca(par_memor,'all2send', vacal,all2send,int(PAR_CODE_SIZE_WM4,ip),lboun=0_ip)
   call memory_alloca(par_memor,'nperm',    vacal,nperm,    max(1_ip,nsendp))
   call memory_alloca(par_memor,'send_perm',vacal,send_perm,max(1_ip,nsendp))
   
   isend = 1_ip
   do irank = 0_ip,PAR_CODE_SIZE_WM4-1
      if(lnsend(irank) > 0_4) then
         all2send(irank) = isend
         call memory_alloca(par_memor,'send_perm % l',vacal,send_perm(isend) % l,int(lnsend(irank),ip)) 
         isend = isend + 1_ip
      endif
   end do

   do irecv = 1, comm0 % lrecv_dim
      do ibuff= 1, mranks
         if(buff0(ibuff,irecv) /= -1) then
            irank = buff0(ibuff,irecv)
            isend =  all2send(irank)
            do iboun = 1, list_nboup(irecv)
               nperm(isend) = nperm(isend) + 1_ip
               send_perm(isend) % l(nperm(isend)) = list_boupo(iboun,irecv)
            enddo
         endif
      enddo
   enddo

   iperm=1_ip
   do isend=1,nsendp
      if(nperm(isend) /= int(size(send_perm(isend) % l),ip)) then
         call runend("mod_par_parallel_partitioning: wrong evaluation of cont_perm")
      endif
      do ibuff=1,nperm(isend)
         comm1 % lsend_perm(iperm) = send_perm(isend) % l(ibuff)
         iperm = iperm + 1_ip
      enddo
   enddo
   call memory_deallo(par_memor,'lnsend',   vacal, lnsend)
   call memory_deallo(par_memor,'send_perm',vacal, send_perm)
   call memory_deallo(par_memor,'all2send', vacal, all2send)
   call memory_deallo(par_memor,'nperm',    vacal, nperm)
   call memory_deallo(par_memor,'buff0',    vacal, buff0)
   !
   ! Transpose communicator comm2
   !
   call comm2 % copy(comm1)
   call PAR_TRANSPOSE_COMMUNICATION_ARRAY(comm2,par_memor)
 
   !
   ! Send candidate boundaries 
   !
   call memory_alloca(par_memor,'LNODB_BUF',vacal,lnodb_buf,mnodb,comm1 % lrecv_dim)
   call memory_alloca(par_memor,'LTYPB_BUF',vacal,ltypb_buf,comm1 % lrecv_dim)
   call PAR_SEND_RECEIVE_TO_ALL(lnodb,lnodb_buf,comm1,'ASYNCHRONOUS',PERMUTE_SEND=.true.)
   call PAR_SEND_RECEIVE_TO_ALL(ltypb,ltypb_buf,comm1,'ASYNCHRONOUS',PERMUTE_SEND=.true.)

   !
   ! Evaluate if I've the boundary
   !
   call memory_alloca(par_memor,'LELBO_VAL',vacal,lelbo_val,comm1 % lrecv_dim,wzero='INITIALIZE')
   call PAR_MAX(pos,PAR_COMM_MY_CODE_WM,rank_max_owner=iboun,INCLUDE_ROOT=.true.)

   nboun_buf = comm1 % lrecv_dim

   if( nboun_buf > 0 ) then
      call par_parallel_partitioning_boundary_elements(&
           nelem,lnods,ltype,&
           nboun_buf,lnodb_buf,ltypb_buf,&
           lelbo_val,chunk_nelem,ipart)
   end if
   !
   ! Receive result for each candidate sent
   !   
   call memory_alloca(par_memor,'LELBO_BUF',vacal,lelbo_buf,comm2 % lrecv_dim)   
   call PAR_SEND_RECEIVE_TO_ALL(lelbo_val,lelbo_buf,comm2,'ASYNCHRONOUS') 
   do ibuff = 1,comm2 % lrecv_dim
      if( lelbo_buf(ibuff) /= 0 ) then
         iboun = comm2 % lrecv_perm(ibuff)
         lelbo(iboun) = lelbo_buf(ibuff)
      end if
   end do
   !
   ! Deallocations
   !
   call htades( ht )
   call memory_deallo(par_memor,'lelbo_buf',vacal,lelbo_buf)
   call memory_deallo(par_memor,'lelbo_val',vacal,lelbo_val)
   call memory_deallo(par_memor,'list_nboup',vacal,list_nboup)
   call memory_deallo(par_memor,'list_boupo',vacal,list_boupo)
   call memory_deallo(par_memor,'lnodb_buf',vacal,lnodb_buf)
   call memory_deallo(par_memor,'ltypb_buf',vacal,ltypb_buf)
   call comm0 % deallo(par_memor,INITIALIZE=.true.)
   call comm1 % deallo(par_memor,INITIALIZE=.true.)
   call comm2 % deallo(par_memor,INITIALIZE=.true.)
   
 end subroutine par_partition_lelbo2

end module mod_par_parallel_partitioning
!> @}

