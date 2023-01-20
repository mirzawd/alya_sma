!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> 
!> @defgroup Partitioning_Toolbox
!> @{
!> @name    Partitioning toolbox
!> @file    mod_paritioning.f90
!> @author  Guillaume Houzeaux
!> @date    24/11/2012
!> @brief   Partitioning
!> @details Partitioning toolbox
!> 
!------------------------------------------------------------------------

module mod_partitioning

  use def_kintyp,               only : ip,rp,lg,i1p
  use def_master,               only : kfl_paral
  use def_master,               only : lninv_loc
  use def_domain,               only : memor_dom
  use def_domain,               only : lmate
  use mod_maths,                only : maths_mapping_coord_to_3d
  use mod_maths,                only : maths_mapping_1d_to_3d
  use mod_maths,                only : maths_local_orthonormal_basis
  use mod_memory,               only : memory_alloca
  use mod_memory,               only : memory_deallo
  use mod_memory,               only : memory_size
  use mod_elmgeo,               only : elmgeo_number_nodes
  use mod_communications,       only : PAR_COMM_RANK_AND_SIZE
  use mod_communications,       only : PAR_ALLGATHER
  use mod_communications,       only : PAR_SUM
  use mod_communications,       only : PAR_MIN
  use mod_communications,       only : PAR_MAX
  use mod_communications,       only : PAR_INTERFACE_NODE_EXCHANGE
  use mod_communications_tools, only : PAR_COMM_TO_INT
  use mod_partition_sfc,        only : partition_sfc
  use mod_partition_sfc,        only : partition_sfc_statistics
  use mod_alya2metis,           only : alya2metis_METIS_PartGraph
  use mod_alya2zoltan,          only : alya2zoltan_Zoltan_LB_Partition
  use mod_partition_sfc_2,      only : partition_sfc_2
  use mod_optional_argument,    only : optional_argument
  use def_mpi
#include "def_mpi.inc"
  implicit none
  private

  public :: partitioning
  public :: partitioning_frontal
  public :: partitioning_frontal_materials
  public :: partitioning_metis
  public :: partitioning_zoltan
  public :: partitioning_oriented_bin

contains

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    24/11/2016
  !> @brief   Partitioning
  !> @details Partition a list into NPART partitions. NPART can be 
  !>          changed.
  !>          Input can be a graph or coordinates
  !>          METIS4           = 0 (topological)
  !>          SFC              = 1 (geometrical)
  !>          Oriented bin     = 2 (geometrical)
  !>          Frontal approach = 3 (topological)
  !>
  !----------------------------------------------------------------------

  subroutine partitioning(imeth,nn,npart,lpart,&
       lweig,rweig,COMM4,ia,ja,lmask,xcoor,&
       boxes_coarse,boxes_fine,direction,ndime_,lcorr_)

    integer(ip),          intent(in)            :: imeth           !< Strategy
    integer(ip),          intent(in)            :: nn              !< Number of nodes
    integer(ip),          intent(inout)         :: npart           !< Number of partitions
    integer(ip), pointer, intent(inout)         :: lpart(:)        !< List of domain
    integer(ip), pointer, intent(in),  optional :: lweig(:)        !< Weights (int)
    real(rp),    pointer, intent(in),  optional :: rweig(:)        !< Weights (real)
    MY_MPI_COMM   ,       intent(in),  optional :: COMM4           !< Communicator (-1 for sequential)
    integer(ip), pointer, intent(in),  optional :: ia(:)           !< Graph IA for topological methods
    integer(ip), pointer, intent(in),  optional :: ja(:)           !< Graph JA
    logical(lg), pointer, intent(in),  optional :: lmask(:)        !< Mask
    real(rp),    pointer, intent(in),  optional :: xcoor(:,:)      !< Coordinates for geometrical methods
    integer(ip),          intent(in),  optional :: boxes_coarse(3) !< Number of boxes (coarse)
    integer(ip),          intent(in),  optional :: boxes_fine(3)   !< Number of boxes (fine)
    real(rp),             intent(in),  optional :: direction(3)    !< Direction for oriented bin
    integer(ip),          intent(in),  optional :: ndime_          !< Problem dimension
    real(rp),    pointer, intent(in),  optional :: lcorr_(:)       !< Partition correction coeficients
    integer(ip)                                 :: kdime

    !
    ! Special cases
    !

       select case ( imeth )

       case ( 0_ip )
          !
          ! METIS
          !
          call partitioning_metis(npart,nn,ia,ja,lpart,lweig,COMM4)

       case ( 1_ip )
          !
          ! Space filling curve (Hilbert)
          !
          if(  present(xcoor)        .and. present(lweig) .and. present(COMM4) .and. &
               present(boxes_coarse) .and. present(boxes_fine) ) then
             if( present(ndime_) ) then
                kdime = ndime_
             else
                kdime = size(xcoor,1,KIND=ip)
                call PAR_MAX(kdime,COMM4,INCLUDE_ROOT=.true.)
             end if
             call partition_sfc(lpart,npart,xcoor,nn,kdime,lweig,boxes_coarse,boxes_fine,lcorr_,COMM4) 
             !call partition_sfc_2(lpart,npart,xcoor,nn,kdime,lweig,lcorr_,COMM4)
          else
             call runend('WRONG NUMBER OF ARGUMENTS FOR SFC PARTITIONING')
          end if

       case ( 2_ip )
          !
          ! Oriented bin
          !
          if(  present(xcoor)      .and. present(lweig)     .and. &
               present(boxes_fine) .and. present(direction) ) then
             if( present(ndime_) ) then
                kdime = ndime_
             else
                kdime = size(xcoor,1,KIND=ip)
                call PAR_MAX(kdime,COMM4,INCLUDE_ROOT=.true.)
             end if
             call partitioning_oriented_bin(&
                  kdime,npart,nn,boxes_fine,direction, &
                  xcoor,lpart,'LOAD BALANCE WITH OPENMP',lweig,COMM4)
          else
             call runend('WRONG NUMBER OF ARGUMENTS FOR ORTIENTED BIN PARTITIONING')
          end if

       case ( 3_ip )
          !
          ! Frontal
          !
          call partitioning_frontal(&
               npart,nn,ia,ja,lpart,lweig,COMM4)

       case ( 4_ip )
          !
          ! Partitioning based on numbering          
          !
          call partitioning_numbering(&
               npart,nn,lpart,lweig,COMM4)
          
       case ( 5_ip )
          !
          ! Partitioning using rank        
          !
          call partitioning_using_rank(&
               npart,nn,lpart,COMM4)

       case ( 6_ip )
          !
          ! Partitioning using rank        
          !
          call partitioning_random(&
               npart,nn,lpart,COMM4)
 
       case ( 7_ip )
          !
          ! ZOLTAN
          !
          if( present(ndime_) ) then
             kdime = ndime_
          else
             kdime = size(xcoor,1,KIND=ip)
             call PAR_MAX(kdime,COMM4,INCLUDE_ROOT=.true.)
          end if
          call partitioning_zoltan(npart,nn,kdime,xcoor,lpart,lweig,PAR_COMM_TO_INT(COMM4))

         
       case default
          !
          ! Unknown
          !
          call runend('PARTITIONING: UNKNOWN METHOD')

       end select

    !end if

  end subroutine partitioning

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    24/11/2016
  !> @brief   Partitioning
  !> @details Partition a list using oriented bin. Layers are created
  !>          along the prescribed direction PVECT_IN. Examples:
  !>
  !>          PVECT_IN=(1,0)      PVECT_IN=(0,1)      PVECT_IN=(0,-1)
  !>          +---+---+---+---+   +---------------+   +---------------+
  !>          |   |   |   |   |   |       3       |   |       1       |
  !>          |   |   |   |   |   +---------------+   +---------------+
  !>          | 1 | 2 | 3 | 4 |   |       2       |   |       2       |
  !>          |   |   |   |   |   +---------------+   +---------------+
  !>          |   |   |   |   |   |       1       |   |       3       |
  !>          +---+---+---+---+   +---------------+   +---------------+
  !>
  !----------------------------------------------------------------------

  subroutine partitioning_oriented_bin(ndime,npart,nn,boxes,pvect_in,xx,lpart,message,lweight,COMM4)

    integer(ip),                    intent(in)    :: ndime        !< Dimension
    integer(ip),                    intent(inout) :: npart        !< Number partitions
    integer(ip),                    intent(in)    :: nn           !< # of entities (node, elements, etc.)
    integer(ip),                    intent(in)    :: boxes(*)     !< # of boxes in each direction
    real(rp),                       intent(in)    :: pvect_in(*)  !< Orientation of bin
    real(rp),    pointer,           intent(in)    :: xx(:,:)      !< Coordinates of entities
    integer(ip), pointer,           intent(inout) :: lpart(:)     !< partition result
    character(*),         optional, intent(in)    :: message      !< Message
    integer(ip), pointer, optional, intent(in)    :: lweight(:)   !< Weight of entities
    MY_MPI_COMM   ,       optional, intent(in)    :: COMM4        !< Communicator
    real(rp)                                      :: comin(3)
    real(rp)                                      :: comax(3)
    real(rp)                                      :: pvect(3)
    real(rp)                                      :: rotma(ndime,ndime)
    real(rp)                                      :: e1(ndime,ndime)
    real(rp)                                      :: e2(ndime,ndime)
    real(rp)                                      :: pnorm
    integer(ip)                                   :: ii,jj,kk,idime
    integer(ip)                                   :: nx,ny,nz,ix,jdime
    integer(ip)                                   :: ipart,kpart
    integer(ip)                                   :: iweight,nweight,nt
    integer(ip)                                   :: weight_target
    logical(lg)                                   :: iflayers
    integer(ip), allocatable                      :: box_weight(:,:,:)
    integer(ip), allocatable                      :: box_num(:,:,:)
    integer(ip), pointer                          :: part_num(:)
    type(i1p),   pointer                          :: box_list(:,:,:) 
    real(rp),    pointer                          :: xx_new(:,:) 
    logical(lg)                                   :: use_openmp
    logical(lg)                                   :: use_mpi

    nullify(box_list)
    nullify(xx_new)
    nullify(part_num)
    use_openmp = .false.
    use_mpi    = .false.
    !
    ! Message
    !
    iflayers = .false.
    if( present(message) ) then
       if( index(message,'LAYER') /= 0 ) then
          iflayers = .true.
       else if( index(message,'LOAD BALANCE') /= 0 ) then 
          iflayers = .false.
       end if
       if( index(message,'WITH OPENMP') /= 0 ) use_openmp = .true.
    end if
    if( present(COMM4) ) then
       if( PAR_COMM_TO_INT(COMM4) /= -1 ) use_mpi = .true.
    end if
    !
    ! Allocate
    !
    if( .not. associated(lpart) ) &
         call memory_alloca(memor_dom,'LPART','partitioning_oriented_bin',lpart,nn)
    !
    ! Compute rotation matrix
    !
    e1              = 0.0_rp
    e1(1,1)         = 1.0_rp
    e1(2,2)         = 1.0_rp
    e1(ndime,ndime) = 1.0_rp
    pvect           = 0.0_rp

    pvect(1:ndime)  = pvect_in(1:ndime)
    pnorm           = sqrt(dot_product(pvect,pvect))
    if( pnorm == 0.0_rp ) call runend('ORIENTED BIN REQUIRES A DIRECTION')
    pvect           = pvect / pnorm
    e2(1:ndime,1)   = pvect(1:ndime)
    call maths_local_orthonormal_basis(ndime,e2)

    do idime = 1,ndime
       do jdime = 1,ndime
          rotma(idime,jdime) = dot_product(e1(1:ndime,jdime),e2(1:ndime,idime))
       end do
    end do
    !
    ! Transform coordinates
    !
    call memory_alloca(memor_dom,'XX_NEW','partitioning_oriented_bin',xx_new,ndime,nn)

    if( use_openmp ) then
       !$OMP PARALLEL DO SCHEDULE ( STATIC )       & 
       !$OMP DEFAULT  ( NONE )                     &
       !$OMP PRIVATE  ( ix,idime,jdime )           &
       !$OMP SHARED   ( ndime,nn,rotma,xx,xx_new )
       do ix = 1,nn
          xx_new(1:ndime,ix) = 0.0_rp
          do idime = 1,ndime
             do jdime = 1,ndime
                xx_new(idime,ix) = xx_new(idime,ix) + rotma(idime,jdime) * xx(jdime,ix)
             end do
          end do
       end do
       !$OMP END PARALLEL DO       
    else
       do ix = 1,nn
          xx_new(1:ndime,ix) = 0.0_rp
          do idime = 1,ndime
             do jdime = 1,ndime
                xx_new(idime,ix) = xx_new(idime,ix) + rotma(idime,jdime) * xx(jdime,ix)
             end do
          end do
       end do
    end if
    !
    ! Min and max coordinates
    !
    do idime = 1,ndime
       comin(idime) = minval(xx_new(idime,:))
       comax(idime) = maxval(xx_new(idime,:))
    end do
    if( use_mpi ) then
       call PAR_MIN(ndime,comin,COMM4,INCLUDE_ROOT=.true.)
       call PAR_MAX(ndime,comax,COMM4,INCLUDE_ROOT=.true.)
    end if
    !
    ! Number of boxes
    !
    if( boxes(1) == 0 .or. boxes(2) == 0 ) call runend('ORIENTED BIEN REQUIRE NUMBER OF BOXES')
    nx = boxes(1)
    ny = boxes(2)
    if( ndime == 3 ) then
       if( boxes(3) == 0 ) call runend('ORIENTED BIEN REQUIRE NUMBER OF BOXES')
       nz = boxes(3)
    else
       nz = 1
    end if
    nt = nx*ny*nz
    !
    ! WEIGHT_TARGET: target weight by subdomain
    ! 
    if( present(lweight) ) then
       nweight = sum(lweight(1:nn))
    else
       nweight = nn
    end if
    if( use_mpi ) then
       call PAR_SUM(nweight,COMM4,INCLUDE_ROOT=.true.)
    end if
    weight_target = max(nweight/npart,1_ip)
    ! 
    ! Number of nodes per box and weight per box
    !
    allocate(box_num   (nx,ny,nz))
    allocate(box_weight(nx,ny,nz))
    call memory_alloca(memor_dom,'BOX_LIST','partitioning_oriented_bin',box_list,nx,ny,nz)
    call memory_alloca(memor_dom,'PART_NUM','partitioning_oriented_bin',part_num,npart)

    box_num    = 0_ip
    box_weight = 0_ip

    if( use_openmp ) then
       !
       ! Reduction variables cannot be alloctable with OpenMP
       !
       !$*OMP  PARALLEL  DO SCHEDULE ( STATIC )                         & 
       !$*OMP  DEFAULT   ( NONE )                                       &
       !$*OMP  PRIVATE   ( ix,ii,jj,kk )                                &
       !$*OMP  SHARED    ( xx_new,nn,ndime,boxes,comin,comax,lweight )  & 
       !$*OMP  REDUCTION ( +:box_num,box_weight )       
       do ix = 1,nn
          call maths_mapping_coord_to_3d(ndime,boxes,comin,comax,xx_new(:,ix),ii,jj,kk)
          box_num(ii,jj,kk) = box_num(ii,jj,kk) + 1
          if( present(lweight) ) then
             box_weight(ii,jj,kk) = box_weight(ii,jj,kk) + lweight(ix)
          else
             box_weight(ii,jj,kk) = box_weight(ii,jj,kk) + 1
          end if
       end do
       !$*OMP END PARALLEL DO       
    else
       do ix = 1,nn
          call maths_mapping_coord_to_3d(ndime,boxes,comin,comax,xx_new(:,ix),ii,jj,kk)
          box_num(ii,jj,kk) = box_num(ii,jj,kk) + 1
          if( present(lweight) ) then
             box_weight(ii,jj,kk) = box_weight(ii,jj,kk) + lweight(ix)
          else
             box_weight(ii,jj,kk) = box_weight(ii,jj,kk) + 1
          end if
       end do
    end if

    if( use_mpi ) call PAR_SUM(nx,ny,nz,box_weight,COMM4,INCLUDE_ROOT=.true.)
    !
    ! Fill boxes woth nodes
    !
    do kk = 1,nz
       do jj = 1,ny
          do ii = 1,nx             
             call memory_alloca(memor_dom,'BOX_LIST % L','partitioning_oriented_bin',box_list(ii,jj,kk) % l,box_num(ii,jj,kk))
          end do
       end do
    end do
    box_num = 0_ip
    if( use_openmp ) then
       do ix = 1,nn
          call maths_mapping_coord_to_3d(ndime,boxes,comin,comax,xx_new(:,ix),ii,jj,kk)
          box_num(ii,jj,kk) = box_num(ii,jj,kk) + 1
          box_list(ii,jj,kk) % l(box_num(ii,jj,kk)) = ix
       end do
    else
       do ix = 1,nn
          call maths_mapping_coord_to_3d(ndime,boxes,comin,comax,xx_new(:,ix),ii,jj,kk)
          box_num(ii,jj,kk) = box_num(ii,jj,kk) + 1
          box_list(ii,jj,kk) % l(box_num(ii,jj,kk)) = ix
       end do
    end if
    !
    ! Partitioning
    !
    if( iflayers ) then
       !
       ! Fill in layers entirely
       !
       ii = 1
       do ipart = 1,npart
          iweight = 0
          do while( iweight < weight_target .and. ii <= nx )             
             do jj = 1,ny
                do kk = 1,nz
                   do ix = 1,box_num(ii,jj,kk)
                      lpart( box_list(ii,jj,kk) % l(ix) ) = ipart
                   end do
                end do
             end do
             iweight = iweight + sum(box_weight(ii,:,:))
             ii = ii + 1
          end do
          part_num(ipart) = iweight
       end do
    else
       !
       ! Fill by layers until target weight is reached
       !
       ii = 1 ; jj = 1 ; kk = 1
       do ipart = 1,npart
          iweight = 0
          do while( iweight < weight_target .and. ii <= nx )             
             do while( iweight < weight_target .and. jj <= ny )        
                do while( iweight < weight_target .and. kk <= nz )        
                   do ix = 1,box_num(ii,jj,kk)
                      lpart( box_list(ii,jj,kk) % l(ix) ) = ipart
                   end do
                   iweight = iweight + box_weight(ii,jj,kk)
                   kk = kk + 1
                end do
                if( kk > nz ) then
                   kk = 1 ; jj = jj + 1
                end if
             end do
             if( jj > ny ) then
                jj = 1 ; ii = ii + 1
             end if
          end do
          part_num(ipart) = iweight
       end do

    end if
    !
    ! Count number of partitions
    !
    kpart = 0
    do ipart = 1,npart
       if( part_num(ipart) /= 0 ) then
          kpart = kpart + 1
          part_num(ipart) = kpart
       end if
    end do
    npart = kpart
    !
    ! Renumber partition number
    ! Assign last partition number to non-marked nodes
    !
    do ix = 1,nn
       if(     lpart(ix) == 0 ) then
          lpart(ix) = kpart
       else if( lpart(ix) > 0 ) then
          lpart(ix) = part_num(lpart(ix))
       end if
    end do
    !
    ! Deallocate memory
    !
    deallocate(box_weight)
    deallocate(box_num)
    call memory_deallo(memor_dom,'BOX_LIST'  ,'partitioning_oriented_bin',box_list  )
    call memory_deallo(memor_dom,'XX_NEW'    ,'partitioning_oriented_bin',xx_new    )
    call memory_deallo(memor_dom,'PART_NUM'  ,'partitioning_oriented_bin',part_num  )

  end subroutine partitioning_oriented_bin

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    24/11/2016
  !> @brief   Partitioning
  !> @details Partition a graph using METIS
  !>
  !----------------------------------------------------------------------

  subroutine partitioning_metis(npart,nn,ia,ja,lpart,lweig,COMM4,onwhat)

    integer(ip),                    intent(inout) :: npart
    integer(ip),                    intent(in)    :: nn
    integer(ip), pointer,           intent(in)    :: ia(:)
    integer(ip), pointer,           intent(in)    :: ja(:)
    integer(ip), pointer,           intent(inout) :: lpart(:)
    integer(ip), pointer, optional, intent(in)    :: lweig(:)
    MY_MPI_COMM   ,       optional, intent(in)    :: COMM4
    character(*),         optional                :: onwhat

!    integer(ip)                                   :: kfl_fortr_par
    integer(ip)                                   :: kpart
!    integer(ip)                                   :: optio_par(8)
    integer(ip)                                   :: kfl_weigh_par
!    integer(ip)                                   :: edgecut,kfl_weigh_par
!    integer(ip)                                   :: wedge_par(2)
    integer(ip), pointer                          :: wvert_par(:)
    integer(ip), pointer                          :: ia_new(:)
    integer(ip), pointer                          :: ja_new(:)
    integer(ip), pointer                          :: lpart_new(:)
    integer(ip), pointer                          :: invpr(:)
    integer(ip), pointer                          :: npart_gat(:) 
    integer(ip)                                   :: nn_new,npart_loc
    integer(ip)                                   :: my_rank,comm_size
    MY_MPI_COMM                                   :: COMM
    integer(4)                                    :: my_rank4,comm_size4
    logical(lg)                                   :: new_graph

    nullify(ia_new)
    nullify(ja_new)
    nullify(lpart_new)
    nullify(invpr)
    !
    ! Parallel partitioning
    !
    if( present(COMM4) ) then 
       call PAR_COMM_RANK_AND_SIZE(COMM4,my_rank4,comm_size4)
       my_rank   = int(my_rank4,ip)
       comm_size = int(comm_size4,ip)
       COMM      = COMM4
       npart_loc = max( npart / comm_size , 1_ip )
    else
       npart_loc = npart
    end if
    !
    ! Check if we have prescribed nodes. In this case, compute a new graph
    ! by eliminating prescribed nodes 
    !
    if( minval(lpart) < 0 ) new_graph = .true.
    if( new_graph ) then
       call runend('PARTITIONING_METIS: CANNOT DEPEND ON GRAPHS TO AVOID CIRCULAR DEPENDENCE')
       !call graphs_compress(nn,ia,ja,lpart,nn_new,ia_new,ja_new,invpr,memor_dom)
       !call memory_alloca(memor_dom,'LPART_NEW','partitioning_metis',lpart_new,nn_new) 
    else
       nn_new    =  nn
       ia_new    => ia
       ja_new    => ja
       lpart_new => lpart
    end if
    !
    ! Weight is given by user
    !
    kfl_weigh_par =  2_ip ! 0: No weights, 1: edges, 2: vertices, 3: both
    if( present(lweig) ) then
       if( new_graph ) then
          call memory_alloca(memor_dom,'WVERT_PAR','partitioning_metis',wvert_par,nn_new)  
          wvert_par(1:nn_new) = lweig(invpr(1:nn_new))
       else
          wvert_par => lweig
       end if
    else
       call memory_alloca(memor_dom,'WVERT_PAR','partitioning_metis',wvert_par,nn_new)
       wvert_par = 1_ip
    end if
    !
    ! Partition with METIS
    !
    call alya2metis_METIS_PartGraph(npart_loc,nn_new,ia_new,ja_new,wvert_par,lpart_new)

    !kfl_fortr_par = 1_ip
    !optio_par(1)  = 0_ip 
    !optio_par(2)  = 3_ip
    !optio_par(3)  = 1_ip
    !optio_par(4)  = 1_ip
    !optio_par(5)  = 0_ip      
    !if( npart_loc <= 8 ) then
    !   call metis_partgraphrecursive(                      &
    !        nn_new, ia_new, ja_new, wvert_par, wedge_par,  &
    !        kfl_weigh_par, kfl_fortr_par, npart_loc,       &
    !        optio_par, edgecut, lpart_new )
    !else
    !   call metis_partgraphkway(                           &
    !        nn_new, ia_new, ja_new, wvert_par, wedge_par,  &
    !        kfl_weigh_par, kfl_fortr_par, npart_loc,       &
    !        optio_par, edgecut, lpart_new )
    !end if
    !
    ! Deallocate memory
    !
    if( present(lweig) ) then
       if( new_graph ) then
          call memory_deallo(memor_dom,'WVERT_PAR','partitioning_metis',wvert_par)  
       end if
    else
       call memory_deallo(memor_dom,'WVERT_PAR','partitioning_metis',wvert_par)  
    end if
    !
    ! Permute if a new graph was created
    !
    if( new_graph ) then
       lpart(invpr(1:nn_new)) = lpart_new(1:nn_new)
       call memory_deallo(memor_dom,'INVPR'    ,'partitioning_metis',invpr)       
       call memory_deallo(memor_dom,'LPART_NEW','partitioning_metis',lpart_new)       
       !call graphs_dealep(ia_new,ja_new)
    end if
    !
    ! Renumber and recompute the partitions in parallel
    !
    if( present(COMM4) ) then
       call memory_alloca(memor_dom,'NPART_GAT','partitioning_frontal',npart_gat,comm_size)       
       call PAR_ALLGATHER(npart_loc,npart_gat,1_4,'CURRENT COMMUNICATOR',COMM4)
       kpart = sum(npart_gat(1:my_rank))
       where ( lpart(1:nn) > 0 ) lpart(1:nn) = lpart(1:nn) + kpart
       call PAR_SUM(npart_loc,COMM4,INCLUDE_ROOT=.true.)
    end if
    npart = npart_loc
    !
    ! If partition is on nodes, distribute values on interface
    !
    if( present(COMM4) .and. present(onwhat) ) then
       if( trim(onwhat) == 'NODE' ) then
          call PAR_INTERFACE_NODE_EXCHANGE(lpart,'DISTRIBUTE')
       end if
    end if

  end subroutine partitioning_metis
  !----------------------------------------------------------------------
  !>
  !> @author  Guillermo Oyarzun
  !> @date    11/11/2019
  !> @brief   Partitioning
  !> @details Partition a graph using ZOLTAN
  !>
  !----------------------------------------------------------------------

  subroutine partitioning_zoltan(npart,nn,ndime,xcoor,lpart,lweig,COMM4,onwhat)

    integer(ip),                    intent(inout) :: npart
    integer(ip),                    intent(in)    :: nn
    integer(ip),                    intent(in)    :: ndime
    real(rp),        intent(in),    pointer       :: xcoor(:,:)
    integer(ip), pointer,           intent(inout) :: lpart(:)
    integer(ip), pointer, optional, intent(in)    :: lweig(:)
    integer(4),           optional, intent(in)    :: COMM4
    character(*),         optional                :: onwhat
    !
    ! Partition with ZOLTAN
    !
    call alya2zoltan_Zoltan_LB_Partition(npart,nn,ndime,lninv_loc,xcoor,lpart,COMM4)
    !
    ! If partition is on nodes, distribute values on interface
    !
    if( present(COMM4) .and. present(onwhat) ) then
       if( trim(onwhat) == 'NODE' ) then
          call PAR_INTERFACE_NODE_EXCHANGE(lpart,'DISTRIBUTE')
       end if
    end if

  end subroutine partitioning_zoltan

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    24/11/2016
  !> @brief   Partitioning
  !> @details Partition a graph into NPART partitions using a frontal 
  !>          technique. NPART can be changed.
  !>          Nodes of the graph with LPART(II)=-1 or with 
  !>          MASK(II)=.TRUE. are excluded from the partition.
  !>
  !>    The main variables are 4:\n
  !>    
  !>     lqueu(:):   array containing the points that a group has.\n
  !>     lgrou(:):   array containg the number of gruops (max. number of groups = max. number of points)\n 
  !>     lfron(:):   contains the points that are in the front.\n
  !>     lmark(:):   contains the points which do not belong to a group.\n
  !>
  !>     Description of other variables:\n
  !>
  !>      ipoin:  point number (it is related with a node)\n
  !>      npogr:  number of points per group. \n
  !>      npomi:  minum number of points per group. \n
  !>      nmark:  number of marked points. \n
  !>      nmarkt: total number of marked points. \n
  !>      kpoin:  point number in lqueu. \n 
  !>      kgrou:  group number. \n
  !>      nfron, nfnew: maximum number of points in lfron(:) in each loop (maxnfron = npoin)\n
  !>      jpoin, iqueu, jqueu, igrou: counting variables.\n
  !> 
  !>      Description of the main loop:\n
  !>
  !>      open-loop:\n
  !>          neigh (loop):...................Loop on the queue. Computes de nonzero ipoin, looks if they belong to a group,\n
  !>                                          if the group is full and if they are marked, allocating them in lqueu.\n
  !>          Adding, cleaning and\n
  !>          compressing nodes of the front..The last point added is then used to make the new front. This is done in three steps:\n
  !>                                               1. First looks for the nonzero points of the mesh that follow this last point \n
  !>                                                  (izdom = r_dom(ipoin), r_dom(ipoin+1)-1)\n
  !>                                               2. Then looks if they belong to a group and if they are marked (HOLE point),\n
  !>                                                  if these two are satisfied, then allocates the point in lgrou, computing the new value\n
  !>                                                  for nfron. \n
  !>                                               3. Cleans lgrou (unmarking the points)\n
  !>                                               4. Computes the new fornt nfnew \n
  !>                                          using the same array in lfron in each loop to allocate the point in the front\n
  !>          Not enough nodes................Looks for nodes that belong to other groups using the last points added.\n
  !>          Final check.....................Checks if all the points are marked anf if not looks for the ones that are not marked.\n
  !>                                  
  !>           !>
  !----------------------------------------------------------------------

  subroutine partitioning_frontal_materials(npart,nn,ia,ja,lpart,lweig,COMM4,onwhat,ja_type)

    use def_domain
    use mod_communications
    use mod_graphs
    
    integer(ip),                    intent(inout) :: npart
    integer(ip),                    intent(in)    :: nn
    integer(ip), pointer, optional, intent(in)    :: ia(:)
    integer(ip), pointer, optional, intent(in)    :: ja(:)
    integer(ip), pointer,           intent(inout) :: lpart(:)
    integer(ip), pointer, optional, intent(in)    :: lweig(:)
    integer(4),           optional, intent(in)    :: COMM4
    character(*),         optional, intent(in)    :: onwhat
    type(i1p),   pointer, optional, intent(in)    :: ja_type(:)
    integer(ip), pointer                          :: lmate_nodes_min(:)
    integer(ip), pointer                          :: lmate_nodes_max(:)
    integer(ip)                                   :: pnode,ielem,inode,ipoin
    integer(ip)                                   :: nn_loc,ii_loc,npart_loc
    integer(ip)                                   :: ii,jj,iz,imate,iz_loc
    integer(ip), pointer                          :: ia_loc(:),ja_loc(:)
    integer(ip), pointer                          :: lpart_loc(:)
    integer(ip), pointer                          :: permn(:)
    integer(ip), pointer                          :: invpn(:)
    integer(ip)                                   :: npart_sav

    nullify(ia_loc)
    nullify(ja_loc)
    nullify(lpart_loc)
    nullify(permn)
    nullify(invpn)
    nullify(lmate_nodes_min)
    nullify(lmate_nodes_max)

    call memory_alloca(memor_dom,'LMATE_NODES_MIN','partitioning_frontal_materials',lmate_nodes_min,nn)
    call memory_alloca(memor_dom,'LMATE_NODES_MAX','partitioning_frontal_materials',lmate_nodes_max,nn)
    do ielem = 1,nelem
       pnode = elmgeo_number_nodes(ltype(ielem),lnods(:,ielem))
       do inode = 1,pnode
          ipoin = lnods(inode,ielem)         
          if( lmate_nodes_min(ipoin) == 0 ) then
             lmate_nodes_min(ipoin) = lmate(ielem)   ! Only one materials
             lmate_nodes_max(ipoin) = lmate(ielem)   ! Only one materials
          else if( lmate_nodes_min(ipoin) /= lmate(ielem) ) then
             lmate_nodes_min(ipoin) = -1             ! More than 1 material
             lmate_nodes_max(ipoin) = -1             ! More than 1 material
          end if
       end do
    end do
    call PAR_INTERFACE_NODE_EXCHANGE(lmate_nodes_min,'MIN')
    call PAR_INTERFACE_NODE_EXCHANGE(lmate_nodes_max,'MAX')
    do ii = 1,nn
       if( lmate_nodes_min(ii) /= lmate_nodes_max(ii) ) lmate_nodes_min(ii) = -1 
    end do
    
    call memory_alloca(memor_dom,'PERMN'    ,'partitioning_frontal_materials',permn,nn)
    call memory_alloca(memor_dom,'INVPN'    ,'partitioning_frontal_materials',invpn,nn)
    call memory_alloca(memor_dom,'IA_LOC'   ,'partitioning_frontal_materials',ia_loc,nn+1)
    call memory_alloca(memor_dom,'JA_LOC'   ,'partitioning_frontal_materials',ja_loc,size(ja,KIND=ip))
    call memory_alloca(memor_dom,'LPART_LOC','partitioning_frontal_materials',lpart_loc,nn)

    npart_sav = npart
    npart = 0
    do imate = 1,nmate
       nn_loc = 0
       do ii = 1,nn
          if( lmate_nodes_min(ii) == imate ) then             
             nn_loc = nn_loc + 1
             permn(nn_loc) = ii
             invpn(ii)     = nn_loc
          end if
       end do
       iz_loc = 0
       do ii_loc = 1,nn_loc
          ii             = permn(ii_loc)
          ia_loc(ii_loc) = 0
          do iz = ia(ii),ia(ii+1)-1
             jj = ja(iz)
             if( lmate_nodes_min(jj) == imate ) then
                iz_loc         = iz_loc + 1
                ia_loc(ii_loc) = ia_loc(ii_loc) + 1
                ja_loc(iz_loc) = invpn(jj)
             end if
          end do
          lpart_loc(ii_loc) = 0
       end do
       call graphs_number_to_linked_list(nn_loc,ia_loc)
       npart_loc = int(real(nn_loc,rp)/real(nn,rp)*real(npart_sav,rp),ip)
       call partitioning_frontal(npart_loc,nn_loc,ia_loc,ja_loc,lpart_loc)
       
       do ii_loc = 1,nn_loc
          ii        = permn(ii_loc)
          lpart(ii) = lpart_loc(ii_loc) + npart
       end do
       npart = npart + npart_loc
    end do
    do ii = 1,nn
       if( lmate_nodes_min(ii) == -1 ) lpart(ii) = -1
    end do

    call memory_deallo(memor_dom,'LMATE_NODES_MIN','partitioning_frontal_materials',lmate_nodes_min)
    call memory_deallo(memor_dom,'LMATE_NODES_MAX','partitioning_frontal_materials',lmate_nodes_max)
    call memory_deallo(memor_dom,'PERMN'          ,'partitioning_frontal_materials',permn)
    call memory_deallo(memor_dom,'INVPN'          ,'partitioning_frontal_materials',invpn)
    call memory_deallo(memor_dom,'IA_LOC'         ,'partitioning_frontal_materials',ia_loc)
    call memory_deallo(memor_dom,'JA_LOC'         ,'partitioning_frontal_materials',ja_loc)
    call memory_deallo(memor_dom,'LPART_LOC'      ,'partitioning_frontal_materials',lpart_loc)

  end subroutine partitioning_frontal_materials
  
  subroutine partitioning_frontal(npart,nn,ia,ja,lpart,lweig,COMM4,onwhat,ja_type)

    integer(ip),                    intent(inout) :: npart
    integer(ip),                    intent(in)    :: nn
    integer(ip), pointer, optional, intent(in)    :: ia(:)
    integer(ip), pointer, optional, intent(in)    :: ja(:)
    integer(ip), pointer,           intent(inout) :: lpart(:)
    integer(ip), pointer, optional, intent(in)    :: lweig(:)
    MY_MPI_COMM   ,       optional, intent(in)    :: COMM4
    character(*),         optional, intent(in)    :: onwhat
    type(i1p),   pointer, optional, intent(in)    :: ja_type(:)
    integer(ip)                                   :: ii,nqueu,npogr,ipart,jpart,iqueu
    integer(ip)                                   :: kpart,kk,izdom,jpoin,npomi,npart_loc
    integer(ip)                                   :: nfron,nfold,jqueu,ifront,nfnew,icheck
    integer(ip)                                   :: nmarkt,nmark,wqueu
    integer(ip)                                   :: my_rank,comm_size
    MY_MPI_COMM                                   :: COMM
    integer(4)                                    :: my_rank4,comm_size4
    integer(ip), pointer                          :: lqueu(:) 
    integer(ip), pointer                          :: lfron(:) 
    integer(ip), pointer                          :: npart_gat(:) 
    logical(lg), pointer                          :: lmark(:) 
    logical(lg)                                   :: if_mpi
    logical(lg), pointer                          :: lpart_check(:)

    nullify(lqueu)
    nullify(lfron)
    nullify(lmark)
    nullify(lpart_check)
    !
    ! Compute the number of partitions per worker in the communicator
    !
    if_mpi    = .false.
    npart_loc = npart

    if( present(COMM4) ) then
       if( PAR_COMM_TO_INT(COMM4) /= -1_4 ) then
          call PAR_COMM_RANK_AND_SIZE(COMM4,my_rank4,comm_size4)
          if_mpi    = .true.
          my_rank   = int(my_rank4,ip)
          comm_size = int(comm_size4,ip)
          COMM      = COMM4
          npart_loc = max( npart / comm_size , 1_ip )
       end if
    end if
    !
    ! Allocate LPART if necessary
    !
    if( .not. associated(lpart) ) then
       call memory_alloca(memor_dom,'LPART','partitioning_frontal',lpart,nn)       
    end if
    !
    ! Allocate memory for LQUEU, LMARK, LFRON
    ! 
    call memory_alloca(memor_dom,'LQUEU','partitioning_frontal',lqueu,nn)
    call memory_alloca(memor_dom,'LMARK','partitioning_frontal',lmark,nn)
    call memory_alloca(memor_dom,'LFRON','partitioning_frontal',lfron,nn)
    !
    ! Count imposed nodes. NMARKT nodes are free
    !
    nmarkt = nn 
    do ii = 1,nn
       if( lpart(ii) == -1 ) then
          nmarkt = nmarkt - 1
          lmark(ii) = .true.
       end if
    end do
    if( nmarkt <= 0 ) call runend('PARTITIONING_FRONTAL: ALL NODES ARE MARKED')
    !
    ! Only one partition!
    !
    if( npart_loc == 1 ) then
       do ii = 1,nn
          if( .not. lmark(ii) ) lpart(ii) = 1
       end do
       goto 1
    end if
    !
    ! Limit number of groups to number of nodes
    !
    if( npart_loc >= nmarkt ) then
       ipart = 0
       do ii = 1,nn
          if( .not. lmark(ii) ) then
             ipart = ipart + 1
             lpart(ii) = ipart
          end if
       end do
       npart_loc = nmarkt
       goto 1
    end if
    !
    ! Compute points per group and minimal threshold
    !
    !rpogr = max(real(nmarkt,rp)/real(npart_loc,rp),1.0_rp)
    !rpomi = rpogr/3.0_rp

    npogr = max(nmarkt/npart_loc,1_ip)
    npomi = npogr/3
    !
    ! Find first non marked point as a starting node
    !
    kk = 1
    do while( lmark(kk) ) 
       kk = kk + 1
    end do
    !
    ! Initialize main loop
    !
    nmark        = 0
    ipart        = 1
    lqueu(1)     = kk
    lpart(kk)    = ipart
    nmark        = nmark + 1
    nqueu        = 1
    iqueu        = 1
    wqueu        = 1
    nfron        = 0

    if( present(ja_type) ) then
       
       open_loop_type: do
          !
          ! Loop on queue
          ! 
          neigh_type: do 
             ii = lqueu(iqueu)
             !
             ! Loop on neighbors
             ! 
             do izdom = 1,memory_size(ja_type(ii)%l)
                jpoin = ja_type(ii) % l(izdom)
                !
                ! Has the group of this point been assigned before ?
                !
                if( lpart(jpoin) == 0 ) then
                   !
                   ! Did we reach the maximum number of points per group
                   !
                   if( present(lweig) ) then
                      if( wqueu >= npogr ) exit neigh_type
                   else
                      if( nqueu == npogr ) exit neigh_type
                   end if
                   lpart(jpoin) = ipart
                   nmark        = nmark + 1    
                   nqueu        = nqueu + 1
                   if( present(lweig) ) wqueu = wqueu + lweig(jpoin)
                   lqueu(nqueu) = jpoin
                   !
                   ! Does this point belong to the current front
                   !
                   if( .not. lmark(jpoin) ) lmark(jpoin) = .true.
                endif
             end do
             !
             ! Did we exhaust the queue?
             !        
             if( iqueu == nqueu ) exit neigh_type
             iqueu = iqueu + 1

          end do neigh_type
          !
          ! Add nodes to the front
          !
          nfold = nfron + 1

          do jqueu = iqueu,nqueu
             ii = lqueu(jqueu)            
             do izdom = 1,memory_size(ja_type(ii)%l)
                jpoin = ja_type(ii) % l(izdom)
                if( lpart(jpoin) == 0 ) then
                   if( .not. lmark(jpoin) ) then
                      nfron = nfron + 1
                      if( nfron > nn ) nfron = nn
                      lfron(nfron) = jpoin
                      lmark(jpoin) = .true.
                   end if
                end if
             end do
          end do
          !
          ! Clean up the last points
          !
          do ifront = nfold,nfron
             ii = lfron(ifront)
             lmark(ii) = .false.
          end do
          !
          ! Compress the front
          !
          nfnew = 0
          do ifront = 1,nfron
             ii = lfron(ifront)
             if( .not. lmark(ii) ) then
                nfnew = nfnew + 1
                lfron(nfnew) = lfron(ifront)
             end if
          end do
          nfron = nfnew
          !
          ! Special case: Do we have enough nodes   
          !
          if( nqueu < npomi ) then
             !
             ! Find a neighboring group
             !
             jpart = 0
             glue_type: do iqueu = 1,nqueu
                ii = lqueu(iqueu)
                do izdom = 1,memory_size(ja_type(ii)%l)
                   jpoin = ja_type(ii) % l(izdom)
                   jpart = lpart(jpoin)
                   if( jpart /= ipart .and. jpart /= -1 ) then
                      exit glue_type
                   endif
                enddo
             end do glue_type

             if( jpart == 0 ) then
                lpart(ii) = -1
             else
                do ii = 1,nn
                   kpart = lpart(ii)
                   if( kpart == ipart ) then
                      lpart(ii) = jpart 
                   end if
                end do
             end if

             nqueu = 0
             if( present(lweig) ) wqueu = 0
             ipart = ipart-1

          end if
          !
          ! Is the front empty --> go home
          !
          if( nfron == 0 ) then
             !
             ! Did we mark all the points
             !
             if( nmark == nmarkt ) then
                !
                ! Is the front empty --> go home
                !
                exit open_loop_type

             else
                !
                ! Find a non marked point
                !
                icheck = 0_ip
                do ii = 1,nn
                   if( lpart(ii) == 0 )then
                      lfron(1) = ii
                      icheck = 1_ip
                      exit
                   end if
                end do
                if( icheck == 0 ) then
                   call runend('GRODOM: nmark/=nmarkt and point not found')
                end if
             end if

          end if
          !
          ! Find new seed
          !
          jpoin        = lfron(1)
          !
          ! Initialize new round
          !
          ipart        = ipart+1
          lqueu(1)     = jpoin
          lpart(jpoin) = ipart
          nmark        = nmark + 1    
          nqueu        = 1
          wqueu        = 1
          iqueu        = 1 
          lmark(jpoin) = .true.

       end do open_loop_type

    else
       
       open_loop: do
          !
          ! Loop on queue
          ! 
          neigh: do 
             ii = lqueu(iqueu)
             !
             ! Loop on neighbors
             ! 
             do izdom = ia(ii),ia(ii+1)-1
                jpoin = ja(izdom)
                !
                ! Has the group of this point been assigned before ?
                !
                if( lpart(jpoin) == 0 ) then
                   !
                   ! Did we reach the maximum number of points per group
                   !
                   if( present(lweig) ) then
                      if( wqueu >= npogr ) exit neigh
                   else
                      if( nqueu == npogr ) exit neigh
                   end if
                   lpart(jpoin) = ipart
                   nmark        = nmark + 1    
                   nqueu        = nqueu + 1
                   if( present(lweig) ) wqueu = wqueu + lweig(jpoin)
                   lqueu(nqueu) = jpoin
                   !
                   ! Does this point belong to the current front
                   !
                   if( .not. lmark(jpoin) ) lmark(jpoin) = .true.
                endif
             end do
             !
             ! Did we exhaust the queue?
             !        
             if( iqueu == nqueu ) exit neigh
             iqueu = iqueu + 1

          end do neigh
          !
          ! Add nodes to the front
          !
          nfold = nfron + 1

          do jqueu = iqueu,nqueu
             ii = lqueu(jqueu)
             do izdom = ia(ii),ia(ii+1)-1
                jpoin = ja(izdom)
                if( lpart(jpoin) == 0 ) then
                   if( .not. lmark(jpoin) ) then
                      nfron = nfron + 1
                      if( nfron > nn ) nfron = nn
                      lfron(nfron) = jpoin
                      lmark(jpoin) = .true.
                   end if
                end if
             end do
          end do
          !
          ! Clean up the last points
          !
          do ifront = nfold,nfron
             ii = lfron(ifront)
             lmark(ii) = .false.
          end do
          !
          ! Compress the front
          !
          nfnew = 0
          do ifront = 1,nfron
             ii = lfron(ifront)
             if( .not. lmark(ii) ) then
                nfnew = nfnew + 1
                lfron(nfnew) = lfron(ifront)
             end if
          end do
          nfron = nfnew
          !
          ! Special case: Do we have enough nodes   
          !
          if( nqueu < npomi ) then
             !
             ! Find a neighboring group
             !
             jpart = 0
             glue: do iqueu = 1,nqueu
                ii = lqueu(iqueu)
                do izdom = ia(ii),ia(ii+1)-1
                   jpoin = ja(izdom)
                   jpart = lpart(jpoin)
                   if( jpart /= ipart .and. jpart /= -1 ) then
                      exit glue
                   endif
                enddo
             end do glue

             if( jpart == 0 ) then
                lpart(ii) = -1
             else
                do ii = 1,nn
                   kpart = lpart(ii)
                   if( kpart == ipart ) then
                      lpart(ii) = jpart 
                   end if
                end do
             end if

             nqueu = 0
             if( present(lweig) ) wqueu = 0
             ipart = ipart-1

          end if
          !
          ! Is the front empty --> go home
          !
          if( nfron == 0 ) then
             !
             ! Did we mark all the points
             !
             if( nmark == nmarkt ) then
                !
                ! Is the front empty --> go home
                !
                exit open_loop

             else
                !
                ! Find a non marked point
                !
                icheck = 0_ip
                do ii = 1,nn
                   if( lpart(ii) == 0 )then
                      lfron(1) = ii
                      icheck = 1_ip
                      exit
                   end if
                end do
                if( icheck == 0 ) then
                   call runend('GRODOM: nmark/=nmarkt and point not found')
                end if
             end if

          end if
          !
          ! Find new seed
          !
          jpoin        = lfron(1)
          !
          ! Initialize new round
          !
          ipart        = ipart+1
          lqueu(1)     = jpoin
          lpart(jpoin) = ipart
          nmark        = nmark + 1    
          nqueu        = 1
          wqueu        = 1
          iqueu        = 1 
          lmark(jpoin) = .true.

       end do open_loop
    end if
    !
    ! Recompute ipart: do not know why but we can end up with
    ! LPART(II) > IPART BECAUSE OF IPART=IPART-1
    !
    ipart = maxval(lpart(1:nn))
    call memory_alloca(memor_dom,'LPART_CHECK','partitioning_frontal',lpart_check,ipart)      
    do ii = 1,nn
       ipart              = lpart(ii)
       if( ipart > 0 ) lpart_check(ipart) = .true.
    end do
    npart_loc = count(lpart_check) 
    call memory_deallo(memor_dom,'LPART_CHECK','partitioning_frontal',lpart_check)       

    !ipart = max(ipart,maxval(lpart(1:nn)))
    !!
    !! Redefine ngrou
    !!
    !npart_loc = max(ipart,npart_loc)

1   continue
    !
    ! Renumber and recompute the partitions in parallel
    !
    if( if_mpi ) then
       call memory_alloca(memor_dom,'NPART_GAT','partitioning_frontal',npart_gat,comm_size)       
       call PAR_ALLGATHER(npart_loc,npart_gat,1_4,'CURRENT COMMUNICATOR',COMM4)
       kpart = sum(npart_gat(1:my_rank))
       do ii = 1,nn
          if( lpart(ii) > 0 ) lpart(ii) = kpart + lpart(ii)
       end do
       call PAR_SUM(npart_loc,COMM4,INCLUDE_ROOT=.true.)
    end if
    !
    ! If partition is on nodes, distribute values on interface
    !
    if( if_mpi .and. present(onwhat) ) then
       if( trim(onwhat) == 'NODE' ) then
          call PAR_INTERFACE_NODE_EXCHANGE(lpart,'DISTRIBUTE')
       end if
    end if
    !
    ! Update number of partitions
    !
    npart = npart_loc
    !
    ! Deallocate memory
    !
    call memory_deallo(memor_dom,'LQUEU','partitioning_frontal',lqueu)
    call memory_deallo(memor_dom,'LMARK','partitioning_frontal',lmark)
    call memory_deallo(memor_dom,'LFRON','partitioning_frontal',lfron)

  end subroutine partitioning_frontal

  logical(lg) function consider_frontal(ipoin,imate,lmate_nodes)

    integer(ip),                    intent(in) :: ipoin
    integer(ip),                    intent(in) :: imate
    integer(ip), optional, pointer, intent(in) :: lmate_nodes(:)

    if( present(lmate_nodes) ) then
       if( lmate_nodes(ipoin) == imate ) then
          consider_frontal = .true.
       else
          consider_frontal = .false.
       end if
    else       
       consider_frontal = .true.
    end if
    
  end function consider_frontal
    
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-05-03
  !> @brief   Partitioning based on numbering
  !> @details Partition the mesh using the element numbering
  !> 
  !-----------------------------------------------------------------------

  subroutine partitioning_numbering(npart,nn,lpart,lweig,COMM4)

    integer(ip),                    intent(inout) :: npart
    integer(ip),                    intent(in)    :: nn
    integer(ip), pointer,           intent(inout) :: lpart(:)
    integer(ip), pointer, optional, intent(in)    :: lweig(:)
    MY_MPI_COMM   ,       optional, intent(in)    :: COMM4
!    integer(ip)                                   :: ii,nn_part,nn_all
!    integer(ip), pointer                          :: nn_gat(:)
!    integer(ip), pointer                          :: start_part(:)
!    integer(ip)                                   :: nn_start,ipart
!    integer(4)                                    :: my_rank4,comm_size4

    call runend('PARTITION_NUMBERING: NOT CODED')
!!$    call PAR_COMM_RANK_AND_SIZE(COMM4,my_rank4,comm_size4)
!!$    nullify(nelem_gat)
!!$    nullify(start_part)
!!$    call memory_alloca(memor_dom,'NN_GAT','partitioning_numbering',nn_gat,PAR_WORLD_SIZE,LBOUN=0_ip)
!!$    nn_part = nn
!!$    if( IMASTER ) nn_part = 0
!!$    call PAR_ALLGATHER(nn_part,nn_gat)
!!$    nn_all  = sum(nn_gat)
!!$    nn_part = nn / npart
!!$    !
!!$    !
!!$    !
!!$    call memory_alloca(memor_dom,'START_PART','partitioning_numbering',start_part,npart)
!!$    start_part(1) = 1
!!$    do ipart = 2,npart
!!$       start_part(ipart) = start_part(ipart-1) + nn_part
!!$    end do
!!$    nn_start = sum(nn_gat(1:my_rank4-1))
!!$    !
!!$    !
!!$    !
!!$    if( .not. associated(lpart) ) &
!!$         call memory_alloca(memor_dom,'LPART','partitioning_numbering',lpart,nn)
!!$    do ii = 1,nn
!!$       ii_global = ii + nn_start
!!$       lpart(ii) = 
!!$    end do
!!$    
  end subroutine partitioning_numbering
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-05-03
  !> @brief   Partitioning based on numbering
  !> @details Partition the mesh using the element numbering
  !> 
  !-----------------------------------------------------------------------

  subroutine partitioning_using_rank(npart,nn,lpart,COMM4)

    integer(ip),                    intent(inout) :: npart
    integer(ip),                    intent(in)    :: nn
    integer(ip), pointer,           intent(inout) :: lpart(:)
    MY_MPI_COMM   ,       optional, intent(in)    :: COMM4
    integer(ip)                                   :: ii
    integer(4)                                    :: my_rank4,comm_size4

    if( .not. associated(lpart) ) &
         call memory_alloca(memor_dom,'LPART','partitioning_numbering',lpart,nn)

    if( present(COMM4) ) then
       call PAR_COMM_RANK_AND_SIZE(COMM4,my_rank4,comm_size4)
       if( npart /= comm_size4 ) call runend('PARTITIONING_USING_RANK: NOT POSSIBLE')       
       if( associated(lpart) ) then
          do ii = 1,nn
             lpart(ii) = int(my_rank4,KIND=ip) + 1_ip
          end do
       end if
    else
       do ii = 1,nn
          lpart(ii) = 1_ip
       end do
    end if
    
  end subroutine partitioning_using_rank

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-05-03
  !> @brief   Partitioning based on... nothing!
  !> @details Random partitioning
  !> 
  !-----------------------------------------------------------------------

  subroutine partitioning_random(npart,nn,lpart,COMM4)

    use mod_random, only : random_initialization
    use mod_random, only : random_generate_number
    
    integer(ip),                    intent(inout) :: npart
    integer(ip),                    intent(in)    :: nn
    integer(ip), pointer,           intent(inout) :: lpart(:)
    MY_MPI_COMM,          optional, intent(in)    :: COMM4
    integer(ip)                                   :: ii,ierro
    real(rp)                                      :: rr,rpart
    integer(ip), pointer                          :: num_part(:)
    
    nullify(num_part)

    call random_initialization(.false.)
    
    if( .not. associated(lpart) ) &
         call memory_alloca(memor_dom,'LPART','partitioning_numbering',lpart,nn)

10  continue
    ierro = 1

    do while( ierro == 1 )
       rpart = real(npart,rp)
       do ii = 1,nn
          rr = random_generate_number(.false.)      
          lpart(ii) = min(npart,int(rr*rpart+1.0_rp,ip))
       end do

       ierro = 0
       
       if( 1 == 1 ) then
          call memory_alloca(memor_dom,'NUM_PART','partitioning_numbering',num_part,npart)
          do ii = 1,nn
             num_part(lpart(ii)) = 1
          end do
          call PAR_SUM(num_part,'IN MY CODE WITHOUT MASTER',INCLUDE_ROOT=.true.)
          ierro = 0
          do ii = 1,npart
             if( num_part(ii) == 0 ) ierro = 1
          end do
          call memory_deallo(memor_dom,'NUM_PART','partitioning_numbering',num_part)
       end if
    end do
    
  end subroutine partitioning_random

end module mod_partitioning
!> @} 
