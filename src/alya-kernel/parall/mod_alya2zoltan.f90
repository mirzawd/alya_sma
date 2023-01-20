!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Parall
!> Brige to ZOLTAN 4 and 5
!> @{
!> @file    mod_alya2zoltan.f90
!> @author  houzeaux
!> @date    2019-01-07
!> @brief   Bridge to ZOLTAN
!> @details Bridge to usefull functions of ZOLTAN 4, 5.0.2 and 5.1.0
!>          WANRING: Integer types of Alya and ZOLTAN must correspond.
!>          The following subroutines are interfaced:
!>          ZOLTAN_SetDefaultOptions (only ZOLTAN5)
!>          ZOLTAN_NodeND
!>          ZOLTAN_PartGraphRecursive
!>          ZOLTAN_PartGraphKway
!>          Manual available at:
!>          http://glaros.dtc.umn.edu/gkhome/fetch/sw/zoltan/manual.pdf
!-----------------------------------------------------------------------

module mod_alya2zoltan

  use def_kintyp,         only : ip,rp
#ifdef I_AM_NOT_ALYA 
  use mod_communications
#else
  use mod_communications, only : PAR_COMM_RANK_AND_SIZE
  use mod_communications, only : PAR_MAX
  use mod_communications, only : PAR_SUM
  use mod_parall,         only : par_memor
  use mod_memory,         only : memory_alloca
  use mod_memory,         only : memory_deallo
  use def_mpi
#include "def_mpi.inc"
#endif

#ifdef ZOLTAN  
  use  zoltan
#endif

#ifdef ZOLTAN  
  !
  ! Mesh data for RCB
  !
  !  integer :: numGlobObjs, numLocObjs
  !  integer(ZOLTAN_INT), dimension(:), allocatable :: GIDs
  !  real,                dimension(:), allocatable :: xcoords, ycoords, zcoords
  !
  ! Zoltan data to store in module
  !
  LOGICAL                                    :: changes 
  INTEGER(Zoltan_INT)                        :: numGidEntries, numLidEntries
  INTEGER(Zoltan_INT)                        :: numImport, numExport
  INTEGER(Zoltan_INT), POINTER, DIMENSION(:) :: importGlobalGids, exportGlobalGids 
  INTEGER(Zoltan_INT), POINTER, DIMENSION(:) :: importLocalGids, exportLocalGids
  INTEGER(Zoltan_INT), POINTER, DIMENSION(:) :: importProcs, exportProcs
  INTEGER(Zoltan_INT), POINTER, DIMENSION(:) :: importToPart, exportToPart

  integer(ip)                                :: ndime
  integer(ip)                                :: numLocObjs
  integer(ip),         pointer               :: GIDs(:)
  real(rp),            pointer               :: xcoords(:)
  real(rp),            pointer               :: ycoords(:)
  real(rp),            pointer               :: zcoords(:)
#endif
  real(rp)                                   :: time_partition
  real(rp)                                   :: time_comm
 
  private

  public :: alya2zoltan_initialization      ! Module initialization
  public :: alya2zoltan_Zoltan_LB_Partition ! Graph partitioning
  public :: alya2zoltan_results             ! Timings

contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-01-07
  !> @brief   ZOLTAN Initialization
  !> @details ZOLTAN Initialization
  !> 
  !-----------------------------------------------------------------------

  subroutine alya2zoltan_initialization()

#if defined  ZOLTAN 
    integer(Zoltan_INT) :: ierr
    real(Zoltan_FLOAT)  :: version       
    ierr = Zoltan_Initialize(version)
#endif

  end subroutine alya2zoltan_initialization

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-01-07
  !> @brief   ZOLTAN part graph
  !> @details ZOLTAN graph partititioning
  !> 
  !-----------------------------------------------------------------------
!  subroutine alya2zoltan_Zoltan_LB_Partition(nparts,nn,coord_z,lninv_z,weight,lpart)


  subroutine alya2zoltan_Zoltan_LB_Partition(nparts,nn,sdime,lninv_z,coord_z,lpart_z,PAR_COMM_)

    integer(ip),                    intent(in)    :: nparts       !< Number of partitions
    integer(ip),                    intent(in)    :: nn           !< Size of graph
    integer(ip),                    intent(in)    :: sdime        !< Size of graph
    integer(ip),           pointer, intent(in)    :: lninv_z(:)
    real(rp),              pointer, intent(in)    :: coord_z(:,:)
    integer(ip),           pointer, intent(inout) :: lpart_z(:)
    integer(4),  optional,          intent(in)    :: PAR_COMM_
#ifdef ZOLTAN  
    integer(ip)                                   :: ii
    type(Zoltan_Struct),   pointer                :: zz_obj
    integer(ZOLTAN_INT)                           :: ierr
    !!#endif
    integer(4)                                    :: my_rank4,comm_size4
    real(rp)                                      :: time1,time2

    call cputim(time1)      

    call PAR_COMM_RANK_AND_SIZE(PAR_COMM_,my_rank4,comm_size4)
    ndime = sdime

    nullify(GIDs)
    nullify(xcoords)
    nullify(ycoords)
    nullify(zcoords)
    nullify(zz_obj)

    if( nn <  1 ) then
       return
    else if( nn <= 1 .or. nparts <= 1 ) then
       lpart_z = 1_ip
       return
    else if( nparts > nn ) then
       do ii = 1,nn
          lpart_z(ii) = ii 
       end do
       return
    else
       !
       ! Initialize
       !
       call alya2zoltan_initialization()
       !
       ! Variables
       !
       GIDs       => lninv_z
       numLocObjs =  nn
       allocate(xcoords(nn))
       allocate(ycoords(nn))
       do ii = 1,nn
          xcoords(ii) = coord_z(1,ii)
          ycoords(ii) = coord_z(2,ii)
       end do
       if( ndime == 3 ) then
          allocate(zcoords(nn))
          do ii = 1,nn
             zcoords(ii) = coord_z(3,ii)
          end do
       end if
       !
       ! ZOLTAN
       !
       zz_obj => Zoltan_Create(PAR_COMM_)

       !       
       ! General Zoltan Parameters
       !
       ierr = Zoltan_Set_Param(zz_obj, "LB_METHOD", "HSFC")

       !
       ! Register query functions
       !
       ierr = Zoltan_Set_Fn(zz_obj, ZOLTAN_NUM_OBJ_FN_TYPE, zoltNumObjs)
       ierr = Zoltan_Set_Fn(zz_obj, ZOLTAN_OBJ_LIST_FN_TYPE,zoltGetObjs)
       ierr = Zoltan_Set_Fn(zz_obj, ZOLTAN_NUM_GEOM_FN_TYPE,zoltNumGeom)
       ierr = Zoltan_Set_Fn(zz_obj, ZOLTAN_GEOM_FN_TYPE,    zoltGeom   )
       !
       ! Use Zoltan to partition the vertices in the simple mesh.
       !
       ! Params:
       !     zz_obj           -- input (all remaining fields are output)
       !     changes          -- 1 if partition was changed, 0 otherwise 
       !     numGidEntries    -- Number of integers used for a global ID 
       !     numLidEntries    -- Number of integers used for a local ID 
       !     numImport        -- Number of vertices to be sent to me 
       !     importGlobalGids -- Global IDs of vertices to be sent to me 
       !     importLocalGids  -- Local IDs of vertices to be sent to me 
       !     importProcs      -- Process rank for source of each incoming vertex 
       !     importToPart     -- New part for each incoming vertex 
       !     numExport        -- Number of vertices I must send to other processes
       !     exportGlobalGids -- Global IDs of the vertices I must send 
       !     exportLocalGids  -- Local IDs of the vertices I must send 
       !     exportProcs      -- Process to which I send each of the vertices 
       !     exportToPart     -- Part to which each vertex will belong 
       !
       ierr = Zoltan_LB_Partition(zz_obj, changes, numGidEntries, numLidEntries, &
            &                     numImport, importGlobalGids, importLocalGids, importProcs, importToPart, &
            &                     numExport, exportGlobalGids, exportLocalGids, exportProcs, exportToPart)


       !Copying back the new partition
       do ii = 1,nn
          lpart_z(ii)= my_rank4+1
       end do
       do ii = 1,numExport
          lpart_z(exportLocalGids(ii)) = exportToPart(ii) +1
       end do
       !
       ! Destroy and clean
       !
       call Zoltan_Destroy(zz_obj)
       call zoltanCleanUp()

    end if

    call cputim(time2)      
    time_partition = time2-time1
    time_comm      = 0.0_rp

#else
    !
    ! Nothing
    !
    call runend('TO USE ZOLTAN, PLEASE COMPILE WITH PROPER MACRO -DZOLTAN AND LIKE WITH IT!')

#endif

  end subroutine alya2zoltan_Zoltan_LB_Partition
  
#ifdef ZOLTAN  
  !
  ! Frees arrays allocated by Zoltan_LB_Partition
  !
  subroutine zoltanCleanup()
    use zoltan
    implicit none

    integer :: error

    error = Zoltan_LB_Free_Part(importGlobalGids, importLocalGids, importProcs, importToPart)
    error = Zoltan_LB_Free_Part(exportGlobalGids, exportLocalGids, exportProcs, exportToPart)

  end subroutine zoltanCleanup
  !
  ! User defined query function to register with Zoltan
  !
  integer function zoltNumObjs(data, ierr)
    use zoltan 
    implicit none

    ! Local declarations
    INTEGER(Zoltan_INT), intent(in)  :: data(*)
    integer(ZOLTAN_INT), intent(out) :: ierr

    zoltNumObjs = numLocObjs
    ierr = ZOLTAN_OK

  end function zoltNumObjs
  !
  ! User defined query function to register with Zoltan
  !
  subroutine zoltGetObjs (data, num_gid_entries, num_lid_entries, global_ids, & 
       local_ids, wgt_dim, obj_wgts, ierr)
    use zoltan
    implicit none

    integer(ZOLTAN_INT), intent(in)  :: data(*)
    integer(ZOLTAN_INT), intent(in)  :: num_gid_entries 
    integer(ZOLTAN_INT), intent(in)  :: num_lid_entries
    integer(ZOLTAN_INT), intent(out) :: global_ids(*)
    integer(ZOLTAN_INT), intent(out) :: local_ids(*)
    integer(ZOLTAN_INT), intent(in)  :: wgt_dim 
    real(ZOLTAN_FLOAT), intent(out)  :: obj_wgts(*)
    integer(ZOLTAN_INT), intent(out) :: ierr

    integer :: i

    do i= 1, numLocObjs
       global_ids(i) = GIDs(i)
       local_ids(i) = i
    end do

    ierr = ZOLTAN_OK

  end subroutine zoltGetObjs
  !
  ! User defined query function to register with Zoltan
  !
  integer function zoltNumGeom(data, ierr)
    use zoltan 
    implicit none
    integer(ZOLTAN_INT), intent(in) :: data(*)
    integer(ZOLTAN_INT)             :: ierr

    zoltNumGeom = ndime
    ierr        = ZOLTAN_OK

  end function zoltNumGeom
  !
  ! User defined query function to register with Zoltan
  !
  subroutine zoltGeom(data, num_gid_entries, num_lid_entries, global_id, &
       local_id, geom_vec, ierr)
    use zoltan
    implicit none

    integer(ZOLTAN_INT), intent(in)  :: data(*)
    integer(ZOLTAN_INT), intent(in)  :: num_gid_entries 
    integer(ZOLTAN_INT), intent(in)  :: num_lid_entries
    integer(ZOLTAN_INT), intent(in)  :: global_id
    integer(ZOLTAN_INT), intent(in)  :: local_id
    real(ZOLTAN_DOUBLE), intent(out) :: geom_vec(*)
    integer(ZOLTAN_INT), intent(out) :: ierr

    geom_vec(1) =  xcoords(local_id)
    geom_vec(2) =  ycoords(local_id)
    if( ndime == 3 ) geom_vec(3) =  zcoords(local_id)

    ierr = ZOLTAN_OK
   
  end subroutine zoltGeom
#endif
  
  subroutine alya2zoltan_results(&
       npart,lepar,lweig,min_weight,max_weight,&
       ave_weight,load_balance,time_total,&
       time_mpi,memor_max,COMM4)

    integer(ip),           intent(in)  :: npart
    integer(ip), pointer,  intent(in)  :: lepar(:)
    integer(ip), pointer,  intent(in)  :: lweig(:)
    real(rp),              intent(out) :: min_weight
    real(rp),              intent(out) :: max_weight
    real(rp),              intent(out) :: ave_weight
    real(rp),              intent(out) :: load_balance
    real(rp),              intent(out) :: time_total
    real(rp),              intent(out) :: time_mpi
    real(rp),              intent(out) :: memor_max
    MY_MPI_COMM   , optional, intent(in)  :: COMM4
    real(rp),    pointer               :: weights(:)
    integer(ip)                        :: ii,ipart
    character(100), parameter          :: vacal = "partition_sfc_results"

    nullify(weights)
    call memory_alloca(par_memor,'WEIGHTS',vacal,weights,npart)
    if( associated(lepar) ) then
       do ii = 1,size(lepar)
          ipart = lepar(ii)
          weights(ipart) = weights(ipart) + real(lweig(ii),rp)
       end do
    end if
    memor_max = real(par_memor(2),rp)
    
    call PAR_SUM(weights,  COMM4)
    call PAR_MAX(time_comm,COMM4)
    call PAR_MAX(memor_max,COMM4)

    min_weight   = minval(weights)
    max_weight   = maxval(weights)
    ave_weight   = sum(weights)/real(npart,rp)
    load_balance = ave_weight/max_weight
    time_total   = time_partition
    time_mpi     = time_comm

    call memory_deallo(par_memor,'WEIGHTS',vacal,weights)

  end subroutine alya2zoltan_results
  
end module mod_alya2zoltan
!> @}
