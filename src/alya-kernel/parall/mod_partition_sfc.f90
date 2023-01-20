!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Parall
!> @{
!> @file    mod_partition_sfc.f90
!> @author  Ricard Borrell and Guillaume Houzeaux
!> @date    2019-09-29
!> @brief   SFC partitioning
!> @details Paralle partitioning using SFC
!-----------------------------------------------------------------------

module mod_partition_sfc

  use def_kintyp_basic,    only : ip,rp,lg,r1p,i1p
  use def_kintyp_comm,     only : comm_data_par
  use mod_maths,           only : maths_sfc_1d_to2d3d_tab
  use mod_maths,           only : maths_sfc_par_part_2
  
#ifdef I_AM_NOT_ALYA 
  use mod_communications
#else
  use mod_memory,          only : memory_alloca
  use mod_memory,          only : memory_deallo
  use mod_parall,          only : par_memor
  use mod_communications,  only : PAR_SUM
  use mod_communications,  only : PAR_MAX
  use mod_communications,  only : PAR_BARRIER
  use mod_communications,  only : PAR_ALLTOALL
  use mod_communications,  only : PAR_SEND_RECEIVE_TO_ALL
  use mod_communications,  only : PAR_COMM_RANK_AND_SIZE
  use def_mpi
#endif
#include "def_mpi.inc"
  
  implicit none

  private
  
  integer(ip), parameter  :: DIM_BIN_CORE        = 128_ip          ! Number of boxes per direction per partitioning process  
  integer(ip), parameter  :: DEFAULT_COARSE_BIN  = 1_ip            ! <1> The minimum larger than the number of procs.
  real(rp),    parameter  :: zeror               = epsilon(1.0_rp) ! Almost zero...
  real(rp),    parameter  :: epsil               = epsilon(1.0_rp)

  integer(ip)             :: dboxf(3)                              ! #bin boxes per direction in the fine bin
  integer(ip)             :: dboxc(3)                              ! #bin boxes per direction in the coarse bin
  integer(ip)             :: dboxl(3)                              ! #bin boxes per direction in the local bin
  integer(ip)             :: nboxc                                 ! #boxes of the coarse bin
  integer(ip)             :: nboxl                                 ! #boxes of the local bin

  real(rp)                :: min_coord(3)                          ! bounding box min
  real(rp)                :: max_coord(3)                          ! bounding bos max

  integer(ip), pointer    :: weigc(:)                              ! Weight of the boxes in the coarse bin
  integer(ip), pointer    :: weigc_loc(:)                          ! Weight of the boxes in the coarse bin within one phase
  integer(ip), pointer    :: weigl(:)                              ! Weight of the boxes in the local bin
  
  integer(ip), pointer    :: partl(:)                              ! Partition of the local bin

  integer(ip), pointer    :: send_buff(:)                          ! Send commuincation buffer
  integer(ip), pointer    :: recv_buff(:)                          ! Receive communication buffer

  integer(ip), pointer    :: lboxes(:)                             ! local box assigned to each element/node (if negative is in the owned coarse box) !rick:better explain
  real(rp),    pointer    :: lcorr(:)                              ! dystribution conrrection

  integer(ip), pointer    :: lepar(:)     
  real(rp),    pointer    :: lenti(:,:)
  integer(ip), pointer    :: lweig(:)   
  integer(ip)             :: nenti
  integer(ip)             :: sdim
  integer(ip)             :: npart_sfc

#ifndef I_AM_NOT_ALYA
  MY_MPI_COMM             :: PAR_COMM
  integer(4)              :: PAR_MY_RANK
  integer(4)              :: PAR_MY_SIZE
#endif
  
  integer(ip)             :: boxes_coarse(3)
  integer(ip)             :: boxes_fine(3)

  real(rp)                :: time_ini
  real(rp)                :: time_partition
  real(rp)                :: time_comm

  integer(ip)             :: nphase   
  integer(ip)             :: iphase   
  integer(ip)             :: iboxc0,iboxc1,iboxc1_
  integer(ip)             :: gbuff0

  real(rp)                :: totalw

  integer(ip), pointer    :: lboxc(:)
  integer(ip), pointer    :: lboxl(:)
  integer(ip), pointer    :: neighs(:)     ! 1: I have elements/points in the corresponding coarse box, 0: none
  integer(4),  pointer    :: lnsend(:)     !number of "local" boxes send/received
  integer(ip), pointer    :: bufwei(:)    ! OPT ! Local buffer to accumulate weights
  type(i1p),   pointer    :: bufse2(:)    ! OPT ! Buffer to send "local" boxes info to the partitioning processes
  integer(ip)             :: nneigs
  integer(ip)             :: ineis

  integer(ip), pointer    :: reord(:)
  integer(ip), pointer    :: invor(:) 
  type(comm_data_par)     :: comm1           ! First point-to-point communication structure
  type(comm_data_par)     :: comm2           ! Second point-to-point communication structure


  logical(lg)             :: if_lcorr
  logical(lg)             :: if_lweig
  !
  ! Output statistics
  !
  real(rp)                :: number_touched_bin
  real(rp)                :: percentage_touched_bin
  real(rp)                :: max_weight_bin
  real(rp)                :: ave_weight_bin
  !
  ! Public functions
  !
  public :: partition_sfc
  public :: partition_sfc_statistics
  public :: partition_sfc_initialization
  public :: partition_sfc_results

contains
  
  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    29/09/2019
  !> @brief   Initialization
  !> @details Initialization of moduel variables
  !
  !----------------------------------------------------------------------
  
  subroutine partition_sfc_initialization()
    !
    ! Initialization and nullify
    !
    number_touched_bin     = 0.0_rp
    percentage_touched_bin = 0.0_rp
    max_weight_bin         = 0.0_rp
    ave_weight_bin         = 0.0_rp
    time_ini               = 0.0_rp
    time_partition         = 0.0_rp
    time_comm              = 0.0_rp
    
    nullify(weigc)                  
    nullify(weigc_loc)            
    nullify(weigl)                
    nullify(partl)                
    nullify(send_buff)            
    nullify(recv_buff)            
    nullify(lboxes)               
    nullify(lcorr)             
    nullify(lepar)     
    nullify(lenti)
    nullify(lweig)
    nullify(lboxc)
    nullify(lboxl)
    nullify(neighs)     
    nullify(lnsend)
    nullify(bufwei)     
    nullify(bufse2)     
    nullify(reord)
    nullify(invor)    

    comm1 % nneig = 0
    comm2 % nneig = 0
    nullify(comm1 % neights)
    nullify(comm2 % neights)
    nullify(comm1 % lsend_size)
    nullify(comm2 % lsend_size)
    nullify(comm1 % lrecv_size)
    nullify(comm2 % lrecv_size)
    
#ifdef I_AM_NOT_ALYA
    par_memor = 0_8
#endif

  end subroutine partition_sfc_initialization
  
  !----------------------------------------------------------------------
  !
  !> @author  Ricard Borrell & Juan Cajas
  !> @date    29/11/2016
  !> @brief   Partition of the mesh using Space Filling Curves (SFC)
  !> @details Main subroutine
  !
  !----------------------------------------------------------------------
  
  subroutine partition_sfc(lepar_,npart_sfc_,lenti_,nenti_,sdim_,lweig_,boxes_coarse_,boxes_fine_,lcorr_,PAR_COMM_)

    use def_master
    integer(ip), pointer,           intent(inout) :: lepar_(:)         !< partition rank of each entity 
    integer(ip),                    intent(inout) :: npart_sfc_        !< number of partitions
    real(rp),    pointer,           intent(in)    :: lenti_(:,:)       !< coordinates of the entities LENTI_(SDIM_,NENTI_)
    integer(ip),                    intent(in)    :: nenti_            !< number of entities 
    integer(ip),                    intent(in)    :: sdim_             !< space dimension (2D/3D)
    integer(ip), pointer, optional, intent(in)    :: lweig_(:)         !< weight of each entitiy LWEIG_(NENTI_)
    integer(ip),          optional, intent(in)    :: boxes_coarse_(3)  !< coarse boxed (to distribute bounding box)
    integer(ip),          optional, intent(in)    :: boxes_fine_(3)    !< find boxes (to define SFC)            
    real(rp),    pointer, optional, intent(in)    :: lcorr_(:)         !< partitition correction coeficients           
    MY_MPI_COMM   ,       optional, intent(in)    :: PAR_COMM_         !< communicator (empty if sequential)

    call partition_sfc_constructor(lepar_,npart_sfc_,lenti_,nenti_,sdim_,lweig_,boxes_coarse_,boxes_fine_,lcorr_,PAR_COMM_)

    call partition_sfc_setup()
    call partition_sfc_define_distribution(PAR_COMM_)
 
    do iphase = 1_ip,nphase

       call partition_sfc_redistribute_weights()
       call partition_sfc_local_partition()
       call partition_sfc_redistribute_result()

    end do

    call partition_sfc_destructor()

  end subroutine partition_sfc

  !-----------------------------------------------------------------------
  !> 
  !> @author  Ricard Borrell
  !> @date    2019-09-29
  !> @brief   Constructor
  !> @details SFC constructor
  !> 
  !-----------------------------------------------------------------------

  subroutine partition_sfc_constructor(&
       lepar_,npart_sfc_,lenti_,nenti_,sdim_,lweig_,boxes_coarse_,boxes_fine_,lcorr_,PAR_COMM_)

    integer(ip),  pointer,           intent(inout) :: lepar_(:)      !< List of elements partition indices !rick: pq ha de ser in tambe?
    integer(ip),                     intent(inout) :: npart_sfc_     !< #partitions
    real(rp),     pointer,           intent(in)    :: lenti_(:,:) 
    integer(ip),                     intent(in)    :: nenti_
    integer(ip),                     intent(in)    :: sdim_
    integer(ip),  pointer, optional, intent(in)    :: lweig_(:) 
    integer(ip),           optional, intent(in)    :: boxes_coarse_(3)
    integer(ip),           optional, intent(in)    :: boxes_fine_(3)                
    real(rp),     pointer, optional, intent(in)    :: lcorr_(:)   
    MY_MPI_COMM   ,        optional, intent(in)    :: PAR_COMM_
    integer(ip)                                    :: ii
    character(100), parameter                      :: vacal = "partition_sfc_constructor"
    !
    ! Deal with optional arguments
    !
    if( present(lweig_) ) then
       if( associated(lweig_) ) then
          if_lweig = .true.
       else
          if_lweig = .false.
       end if
    else
       if_lweig = .false.
    end if       
    if( present(lcorr_) ) then
       if( associated(lcorr_) ) then
          if_lcorr = .true.
       else
          if_lcorr = .false.
       end if
    else
       if_lcorr = .false.
    end if       
    !
    ! Take module arguments
    !
    if( .not. associated(lepar_) ) then
       nullify(lepar_)
       call memory_alloca(par_memor,'lepar_',vacal,lepar_,nenti_)
    end if
    lepar     => lepar_
    npart_sfc =  npart_sfc_
    lenti     => lenti_
    nenti     =  nenti_
    sdim      =  sdim_
    
    if( if_lweig ) then 
       lweig => lweig_
    else
       nullify(lweig)
       call memory_alloca(par_memor,'lweig',vacal,lweig,max(1_ip,nenti_))
       lweig = 1_ip
    end if
    
    boxes_coarse(1:3) = 0_ip
    boxes_fine(1:3)   = 0_ip
    if(present(boxes_coarse_)) boxes_coarse(1:3) = boxes_coarse_(1:3)
    if(present(boxes_fine_))   boxes_fine(1:3)   = boxes_fine_(1:3)
    if(present(PAR_COMM_)) then
       PAR_COMM=PAR_COMM_
       call PAR_COMM_RANK_AND_SIZE(PAR_COMM,PAR_MY_RANK,PAR_MY_SIZE)
    else   
       PAR_MY_SIZE=1_4
    endif
    !
    !  Check boxes arguments
    !
    if(         present(boxes_coarse_)  .and. ( .not. present(boxes_fine_)) .or. &
         (.not. present(boxes_coarse_)) .and.         present(boxes_fine_)) then
       call runend("ERROR mod_partition_sfc: both coarse boxes and fine boxes must be present or absent in the arguments")
    endif
    if( present(boxes_coarse_) .and. present(boxes_fine_) ) then
       if(      (boxes_coarse_(1)/=0_ip .and. boxes_fine_(1)==0_ip) .or. &
            &   (boxes_coarse_(1)==0_ip .and. boxes_fine_(1)/=0_ip) )then
          call runend("ERROR mod_partition_sfc: both coarse boxes and fine boxes must be defined or not defined")
       endif
    end if
    do ii = 2,sdim
       if( boxes_coarse(1) /= boxes_coarse(ii) ) then
          call runend("ERROR mod_partition_sfc: the boxes grid must be isotropic")
       endif
       if( boxes_fine(1) /= boxes_fine(ii) ) then
          call runend("ERROR mod_partition_sfc: the boxes grid must be isotropic")
       endif
    enddo
    !
    ! Get correction coeficients from file
    !
    if( if_lcorr  ) then 
       lcorr => lcorr_
    else
       nullify(lcorr)
       call memory_alloca(par_memor,'lcorr',vacal,lcorr,max(1_ip,npart_sfc))
       lcorr(1:npart_sfc) = 1.0_rp
    endif
    
    ineis = 1_ip
     
    call PAR_BARRIER(COMM=PAR_COMM_)
    call cputim(time_ini)

  end subroutine partition_sfc_constructor
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  bsc21943
  !> @date    2019-02-05
  !> @brief   SFC statistics
  !> @details x1 = Ratio of touched bin (give an idea on the 
  !>               anisotropy of the mesh)
  !>          x2 = Max weight
  !>          x3 = Average weight
  !>          x4 = Granularity = max/sum (realtive weight of the 
  !>               heaviest bin)
  !> 
  !-----------------------------------------------------------------------

  subroutine partition_sfc_results(&
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
    
  end subroutine partition_sfc_results
  
  subroutine partition_sfc_statistics(x1,x2,x3,x4,x5,i1,i2,COMM4)

    real(rp),       intent(out)          :: x1 
    real(rp),       intent(out)          :: x2
    real(rp),       intent(out)          :: x3
    real(rp),       intent(out)          :: x4
    real(rp),       intent(out)          :: x5
    integer(ip),    intent(out)          :: i1(:)
    integer(ip),    intent(out)          :: i2(:)
    MY_MPI_COMM   , intent(in), optional :: COMM4
    integer(ip)                          :: boxes(2)

    x1 = number_touched_bin
    x2 = max_weight_bin
    x3 = ave_weight_bin
    x5 = real(nboxc,rp)*real(nboxl,rp)

    if( present(COMM4) ) then
       call PAR_SUM(x1,COMM4)
       call PAR_MAX(x2,COMM4)
       call PAR_SUM(x3,COMM4)
       call PAR_MAX(x5,COMM4)
    else
       call PAR_SUM(x1,PAR_COMM)
       call PAR_MAX(x2,PAR_COMM)
       call PAR_SUM(x3,PAR_COMM)
       call PAR_MAX(x5,PAR_COMM)
    end if

    if( x1 /= 0.0_rp ) then
       x4 = x2 / (x3+zeror)
       x3 = x3 / (x1+zeror)
       x1 = x1 / (x5+zeror)
    else
       x1 = 0.0_rp
       x2 = 0.0_rp
       x3 = 0.0_rp
       x4 = 0.0_rp
       x5 = 0.0_rp
    end if

    boxes(1) = dboxf(1)
    boxes(2) = dboxc(1)
    if( present(COMM4) ) then
       call PAR_MAX(2_ip,boxes,COMM=COMM4)
    else
       call PAR_MAX(2_ip,boxes,COMM=PAR_COMM)
    end if
    i1(:) = boxes(1)
    i2(:) = boxes(2)
   
  end subroutine partition_sfc_statistics

  !-----------------------------------------------------------------------
  !> 
  !> @author  Ricard Borrell
  !> @date    2019-09-29
  !> @brief   SFC setup
  !> @details Define bins
  !> 
  !-----------------------------------------------------------------------

  subroutine partition_sfc_setup()

    implicit none
    real(rp),       pointer   :: aux_coord(:)
    integer(ip)               :: idime, lev,ii
    real(rp)                  :: r1
    integer(ip)               :: iproc,x,y,z,coo1d
    character(100), parameter :: vacal = "partition_sfc_setup"
    !
    ! Set number of slaves according to hilbert restriction
    !
    if(DEFAULT_COARSE_BIN==0_ip) then
       lev = int(log(real(PAR_MY_SIZE,rp))/log((2.0_rp)**(sdim)),ip)
    elseif(DEFAULT_COARSE_BIN==1_ip) then
       r1  = log(real(PAR_MY_SIZE,rp))/log((2.0_rp)**(sdim))
       lev = int(r1,ip)
       if(r1 > real(lev,rp)) then
          lev = lev + 1_ip
       endif
    else
       call runend("ERROR mod_partition_sfc: unexpected DEFAULT_COARSE_BIN")
    endif
    dboxc(1:3) = 2_ip**lev
    dboxf(1:3) = DIM_BIN_CORE*dboxc(1:3)
    if( boxes_fine(1) /= 0 .and. boxes_coarse(1) /= 0 ) then
       dboxc(1:3) = boxes_coarse(1)
       dboxf(1:3) = boxes_fine(1)               
    endif
    dboxl(1:3) = dboxf(1:3) / dboxc(1:3)
    !
    ! Bounding box
    !
    nullify(aux_coord)
    allocate(aux_coord(6))
    aux_coord(1:3)=-huge(1.0_rp)*0.1_rp
    aux_coord(4:6)=-huge(1.0_rp)*0.1_rp
    if(nenti>0) then
       do idime = 1_ip,sdim
          aux_coord(idime)        =  maxval(lenti(idime,1:nenti))
          aux_coord(3_ip+idime)   = -minval(lenti(idime,1:nenti))
       end do
    endif
    call PAR_MAX(aux_coord,PAR_COMM,INCLUDE_ROOT=.true.)

    max_coord(1:3) =  aux_coord(1:3)
    min_coord(1:3) = -aux_coord(4:6)
    max_coord(1:3) = max_coord(1:3) + (max_coord(1:3)-min_coord(1:3))*0.001_rp + zeror
    min_coord(1:3) = min_coord(1:3) - (max_coord(1:3)-min_coord(1:3))*0.001_rp - zeror
    deallocate(aux_coord)
    !
    ! Eval bin sizes
    !
    nboxc = dboxc(1)**sdim
    nboxl = dboxl(1)**sdim
    !
    ! Phases
    !
    nphase = int(nboxc/PAR_MY_SIZE,ip)
    if(mod(nboxc,int(PAR_MY_SIZE,ip)) > 0_ip) nphase = nphase + 1_ip;
    !
    ! totalw
    !
    totalw=0.0_rp
    do ii = 1_ip,nenti
       totalw = totalw + real(lweig(ii),rp) 
    enddo
    call PAR_SUM(totalw, PAR_COMM,INCLUDE_ROOT=.true.)
    !
    ! Ordering of coarse bins
    !
    gbuff0=0_ip
    nullify(reord,invor)
    call memory_alloca(par_memor,'reord'    ,vacal,reord,nboxc)
    call memory_alloca(par_memor,'invor'    ,vacal,invor,nboxc)
    do iproc=1,nboxc
       if(sdim==2_ip) then
          call maths_sfc_1d_to2d3d_tab(dboxc(1),iproc,x,y)
          coo1d=(y-1)*dboxc(1)+x
       else
          call maths_sfc_1d_to2d3d_tab(dboxc(1),iproc,x,y,z)
          coo1d=(z-1)*dboxc(1)*dboxc(1)+(y-1)*dboxc(1)+x
       endif
       reord(iproc) = coo1d
       invor(coo1d) = iproc
    enddo
    !
    ! SFC statistics
    !
    number_touched_bin     =  0.0_rp
    percentage_touched_bin =  0.0_rp
    max_weight_bin         = -1.0_rp
    ave_weight_bin         =  0.0_rp

  end subroutine partition_sfc_setup
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  Ricard Borrell
  !> @date    2019-09-29
  !> @brief   Define distribution
  !> @details Define distribution
  !> 
  !-----------------------------------------------------------------------

  subroutine partition_sfc_define_distribution(PAR_COMM_)

    MY_MPI_COMM               :: PAR_COMM_
    integer(ip)               :: bcoof(3),bcooc(3),bcool(3) 
    integer(ip)               :: ielem, gbuff, ibuff
    integer(ip)               :: iboxc, iboxl, ineig, kbuff
    type(i1p),      pointer   :: lenbc(:)
    integer(ip),    pointer   :: auxbf(:)
    character(100), parameter :: vacal = "partition_sfc_define_distribution"
    real(rp)                  :: dboxf_rp(3)
    real(rp)                  :: dboxc_rp(3)
    real(rp)                  :: dx(3)
    !
    ! Allocate pointers and initialize
    !
    nullify(weigc,weigc_loc,lboxes,weigl,lnsend)
    nullify(lboxc,lboxl,neighs,bufwei,bufse2,lenbc,auxbf)

    call memory_alloca(par_memor,'weigc_loc' ,vacal,weigc_loc, max(1_ip,int(PAR_MY_SIZE,ip)))
    call memory_alloca(par_memor,'weigc'     ,vacal,weigc,     max(1_ip,nboxc))
    call memory_alloca(par_memor,'weigl'     ,vacal,weigl,     max(1_ip,nboxl))
    call memory_alloca(par_memor,'neighs'    ,vacal,neighs,    max(1_ip,nboxc))
    call memory_alloca(par_memor,'lboxc'     ,vacal,lboxc,     max(1_ip,nenti))
    call memory_alloca(par_memor,'lboxl'     ,vacal,lboxl,     max(1_ip,nenti))

    neighs(:) = 0_ip
    nneigs    = 0_ip
    bcool     = 1_ip
    lboxc     = 0_ip
    lboxl     = 0_ip
    dboxf_rp  = real(dboxf,rp)
    dboxc_rp  = real(dboxc,rp)
    dx        = (max_coord-min_coord)

    if( sdim == 2 ) then
       do ielem = 1_ip, nenti
          !
          ! Eval. index of the mc box (box coordinates) at the different bins: bcoof,bcooc,bcool
          !
          bcoof(1:2) = int( ( (lenti(1:2,ielem)-min_coord(1:2)-epsil) / dx(1:2) )*dboxf_rp(1:2), ip ) + 1_ip
          bcooc(1:2) = int( ( (lenti(1:2,ielem)-min_coord(1:2)-epsil) / dx(1:2) )*dboxc_rp(1:2), ip ) + 1_ip
          bcoof(1:2) = min(bcoof(1:2),dboxf(1:2))
          bcooc(1:2) = min(bcooc(1:2),dboxc(1:2))
          bcool(1:2) = bcoof(1:2)-(bcooc(1:2)-1)*dboxl(1:2)
          bcoof(3)   = 1
          bcooc(3)   = 1
          !
          ! Eval. position of the element in the boxes grid according to lexicografic order
          !
          lboxc(ielem) = invor(dboxc(2_ip)*dboxc(1_ip)*(bcooc(3)-1_ip) + dboxc(1_ip)*(bcooc(2)-1_ip) + bcooc(1))
          lboxl(ielem) = dboxl(2_ip)*dboxl(1_ip)*(bcool(3)-1_ip) + dboxl(1_ip)*(bcool(2)-1_ip) + bcool(1)

          if(neighs(lboxc(ielem)) == 0_ip) then
             nneigs = nneigs + 1_ip
          endif
          neighs(lboxc(ielem)) = neighs(lboxc(ielem)) + 1_ip
       end do
    else
       do ielem = 1_ip, nenti
          !
          ! Eval. index of the mc box (box coordinates) at the different bins: bcoof,bcooc,bcool
          !
          bcoof(1:3) = int( ( (lenti(1:3,ielem)-min_coord(1:3)-epsil) / dx(1:3) )*dboxf_rp(1:3), ip ) + 1_ip
          bcooc(1:3) = int( ( (lenti(1:3,ielem)-min_coord(1:3)-epsil) / dx(1:3) )*dboxc_rp(1:3), ip ) + 1_ip
          bcoof(1:3) = min(bcoof(1:3),dboxf(1:3))
          bcooc(1:3) = min(bcooc(1:3),dboxc(1:3))
          bcool(1:3) = bcoof(1:3)-(bcooc(1:3)-1)*dboxl(1:3)
          !
          ! Eval. position of the element in the boxes grid according to lexicografic order
          !
          lboxc(ielem) = invor(dboxc(2_ip)*dboxc(1_ip)*(bcooc(3)-1_ip) + dboxc(1_ip)*(bcooc(2)-1_ip) + bcooc(1))
          lboxl(ielem) = dboxl(2_ip)*dboxl(1_ip)*(bcool(3)-1_ip) + dboxl(1_ip)*(bcool(2)-1_ip) + bcool(1)

          if(neighs(lboxc(ielem)) == 0_ip) then
             nneigs = nneigs + 1_ip
          endif
          neighs(lboxc(ielem)) = neighs(lboxc(ielem)) + 1_ip
       end do
    end if
    !
    !  Allocate and fill lenbc
    !
    call memory_alloca(par_memor,'lenbc'     ,vacal,lenbc,max(1_ip,nboxc))
    do iboxc = 1_ip, nboxc
       if(neighs(iboxc) > 0_ip) then
          call memory_alloca(par_memor,'lenbc % l',vacal,lenbc(iboxc) % l,neighs(iboxc))
       endif
    enddo
    neighs(:)=0_ip
    do ielem = 1_ip, nenti
       iboxc = lboxc(ielem)
       neighs(iboxc)  = neighs(iboxc) + 1_ip
       lenbc(iboxc) % l(neighs(iboxc)) = ielem 
    end do

    !
    ! Store information being sent in temporal buffers
    !
    call memory_alloca(par_memor,'lboxes',vacal,lboxes,max(1_ip,nenti))
    call memory_alloca(par_memor,'bufwei',vacal,bufwei,max(1_ip,2*nboxl))
    call memory_alloca(par_memor,'auxbf' ,vacal,auxbf, max(1_ip,nenti))

    ineig = 0_ip
    gbuff = 1_ip !accumulates total of "local" boxes sent

    call memory_alloca(par_memor,'bufse2',vacal,bufse2,max(1_ip,nneigs))
    call memory_alloca(par_memor,'lnsend',vacal,lnsend,int(max(1_ip,nboxc),4),lboun=0_4)
    
    do iboxc = 1_ip, nboxc

       if( neighs(iboxc) > 0_ip) then
          ineig = ineig + 1_ip
          !
          ! Count boxes to be sent and its weight
          !
          do ibuff = 1_ip, neighs(iboxc) 
             ielem = lenbc(iboxc) % l(ibuff)
             iboxl = lboxl(ielem)
             !
             ! In bufwei(iboxl*2) we store store iboxf
             ! In bufwei((iboxl-1)*2+1) we accumulate the weight of the bi
             !
             if( bufwei(iboxl*2) == 0_ip ) then
                lnsend(iboxc-1_ip)     = lnsend(iboxc-1_ip) + 1_4
                auxbf(lnsend(iboxc-1)) = iboxl
                bufwei( iboxl*2 )      = 1_ip
             endif
             bufwei((iboxl-1)*2+1) = bufwei((iboxl-1)*2+1) + lweig(ielem)
             lboxes(ielem)         = iboxl
          end do
          !
          ! Save info to be sent in temporal buffer
          !
          call memory_alloca(par_memor,'bufse2 % l',vacal,bufse2(ineig) % l,max(1_ip,2_ip*lnsend(iboxc-1)))
          ibuff = 1_ip
          !
          ! we send:
          !  1. local index of the box
          !  2. asssociated weight
          !
          do kbuff=1,lnsend(iboxc-1_ip)
             iboxl                            = auxbf(kbuff)
             bufse2(ineig) % l((ibuff-1)*2+1) = iboxl                  !iboxl
             bufse2(ineig) % l(ibuff*2)       = bufwei((iboxl-1)*2+1)  !weight
             bufwei(iboxl*2)                  = gbuff*2
             ibuff                            = ibuff + 1_ip
             gbuff                            = gbuff + 1_ip
          end do
          !
          ! Store the place where the result will be obtained in the buffer
          !
          do ibuff = 1_ip, neighs(iboxc) 
             ielem         = lenbc(iboxc) % l(ibuff)
             lboxes(ielem) = int(bufwei(lboxes(ielem)*2),ip)
          end do
          !
          ! Initialize to 0.0 the components of bufwei used
          !
          do kbuff = 1,lnsend(iboxc-1_ip)
             iboxl                 = auxbf(kbuff)
             bufwei((iboxl-1)*2+1) = 0_ip
             bufwei(iboxl*2)       = 0_ip
          end do
       end if
    end do
    call memory_deallo(par_memor,'lenbc',vacal,lenbc)
    call memory_deallo(par_memor,'auxbf',vacal,auxbf)

  end subroutine partition_sfc_define_distribution

  !-----------------------------------------------------------------------
  !> 
  !> @author  Ricard Borrell
  !> @date    2019-09-29
  !> @brief   Box storing
  !> @details Create a bin structure of the geometrical domain storing the
  !>          elements in the boxes
  !> 
  !-----------------------------------------------------------------------

  subroutine partition_sfc_redistribute_weights()

    integer(ip)               :: ineig, neigh, icont
    integer(ip)               :: iboxc,iboxl, ibuff
    character(100), parameter :: vacal = "partition_sfc_redistribute_weights"
    real(rp)                  :: time1,time2
    integer(4),     pointer   :: lnsend_loc(:) 
    integer(4),     pointer   :: lnrecv_loc(:) 
    
    nullify(send_buff,recv_buff,lnsend_loc,lnrecv_loc)
    !
    ! iboxc0 and iboxc1 limits of the coarse bins treated in this phase
    !
    iboxc0  = (iphase-1_ip)*PAR_MY_SIZE
    iboxc1  = iboxc0+PAR_MY_SIZE-1_ip
    iboxc1_ = min(iboxc1,nboxc-1)
    weigl   = 0_ip

    call memory_alloca(par_memor,'lnsend_loc',vacal,lnsend_loc,int(max(1_ip,int(PAR_MY_SIZE,ip)),4),lboun=0_4)
    call memory_alloca(par_memor,'lnrecv_loc',vacal,lnrecv_loc,int(max(1_ip,int(PAR_MY_SIZE,ip)),4),lboun=0_4)
    
    do iboxc = iboxc0,iboxc1_
       lnsend_loc(iboxc-iboxc0) = lnsend(iboxc)
    end do
    !
    ! Share communication requierements with others slaves
    !
    call cputim(time1)
    call PAR_ALLTOALL(lnsend_loc,lnrecv_loc,PAR_COMM_IN=PAR_COMM)
    call cputim(time2)
    time_comm = time_comm + (time2-time1)
    !
    ! Allocate par_memory for communication structure
    !
    comm1 % PAR_COMM_WORLD = PAR_COMM
    comm2 % PAR_COMM_WORLD = PAR_COMM
    comm1 % nneig          = 0_ip
    
    do iboxc = iboxc0, iboxc1 
       if(lnrecv_loc(iboxc-iboxc0) > 0_ip .or. lnsend_loc(iboxc-iboxc0) > 0_ip) then
          comm1 % nneig = comm1 % nneig + 1
       endif
    end do
    comm2 % nneig = comm1 % nneig
    call memory_alloca(par_memor,'COMM1 % NEIGHTS'   ,vacal,comm1 % neights,   max(1_ip, comm1 % nneig   ))
    call memory_alloca(par_memor,'COMM1 % LSEND_SIZE',vacal,comm1 % lsend_size,max(1_ip, comm1 % nneig+1 ))
    call memory_alloca(par_memor,'COMM1 % LRECV_SIZE',vacal,comm1 % lrecv_size,max(1_ip, comm1 % nneig+1 ))
    call memory_alloca(par_memor,'COMM2 % NEIGHTS'   ,vacal,comm2 % neights,   max(1_ip, comm2 % nneig   ))
    call memory_alloca(par_memor,'COMM2 % LSEND_SIZE',vacal,comm2 % lsend_size,max(1_ip, comm2 % nneig+1 ))
    call memory_alloca(par_memor,'COMM2 % LRECV_SIZE',vacal,comm2 % lrecv_size,max(1_ip, comm2 % nneig+1 ))
    
    comm1 % lsend_dim     = 0_ip
    comm1 % lrecv_dim     = 0_ip
    comm1 % lsend_size(1) = 1_ip
    comm1 % lrecv_size(1) = 1_ip
    ineig                 = 1_ip

    do iboxc = iboxc0, iboxc1
       if(lnrecv_loc(iboxc-iboxc0) > 0_ip .or. lnsend_loc(iboxc-iboxc0) > 0_ip) then
          comm1 % neights(ineig)      = iboxc-iboxc0
          comm1 % lsend_dim           = comm1 % lsend_dim + 2*lnsend_loc(iboxc-iboxc0)
          comm1 % lsend_size(ineig+1) = comm1 % lsend_size(ineig) + 2_ip*lnsend_loc(iboxc-iboxc0)
          comm1 % lrecv_dim           = comm1 % lrecv_dim + 2*lnrecv_loc(iboxc-iboxc0)
          comm1 % lrecv_size(ineig+1) = comm1 % lrecv_size(ineig) + 2_ip*lnrecv_loc(iboxc-iboxc0)
          ineig                       = ineig + 1_ip
       endif
    end do
    comm2 % lsend_dim     = comm1 % lrecv_dim
    comm2 % lrecv_dim     = comm1 % lsend_dim
    comm2 % neights(:)    = comm1 % neights(:)
    comm2 % lsend_size(:) = comm1 % lrecv_size(:)
    comm2 % lrecv_size(:) = comm1 % lsend_size(:)
    !
    ! Store data in unified buffer, perform boxes transfer, save local boxes infor in weiglh/sub_bins
    !
    call memory_alloca( par_memor,'SEND_BUFF',vacal,send_buff,max(1_ip,2*comm1 % lsend_dim))
    call memory_alloca( par_memor,'RECV_BUFF',vacal,recv_buff,max(1_ip,2*comm1 % lrecv_dim))
    icont = 1_ip
    do ineig = 1_ip, comm1 % nneig
       neigh = comm1 % neights(ineig)
       if(lnsend_loc(neigh) > 0) then
          do ibuff = 1_ip, lnsend_loc(neigh)
             send_buff(((icont-1)*2)+1)  = bufse2(ineis) % l((ibuff-1)*2+1) !iboxl
             send_buff(icont*2)          = bufse2(ineis) % l(ibuff*2)       !weight
             icont                       = icont + 1_ip
          end do
          ineis = ineis + 1_ip
       end if
    end do
    call cputim(time1)
    call PAR_SEND_RECEIVE_TO_ALL(send_buff,recv_buff,comm1,'ASYNCHRONOUS')
    call cputim(time2)
    time_comm = time_comm + (time2-time1)    
    icont = 1_ip
    do ineig = 1_ip, comm1 % nneig
       neigh = comm1 % neights(ineig)
       do ibuff = 1_ip, lnrecv_loc(neigh)
          iboxl = recv_buff((icont-1)*2+1)
          if( weigl(iboxl) == 0_ip ) then
             number_touched_bin = number_touched_bin + 1.0_rp
          end if
          ave_weight_bin = ave_weight_bin + real(recv_buff(icont*2),rp)
          weigl(iboxl)   = weigl(iboxl) + recv_buff(icont*2)
          icont          = icont + 1             
          max_weight_bin = max(max_weight_bin,real(weigl(iboxl),rp))
       end do
    end do
    !
    ! Eval weigh distribution in coarse bin
    !   
    weigc_loc = 0
    do icont = 1_ip, nboxl
       weigc_loc(PAR_MY_RANK + 1_ip) = weigc_loc(PAR_MY_RANK + 1_ip) + weigl(icont)
    end do
    call cputim(time1)
    call PAR_SUM(weigc_loc, PAR_COMM, INCLUDE_ROOT=.true.)

    call memory_deallo(par_memor,'LNSEND_LOC',vacal,lnsend_loc)
    call memory_deallo(par_memor,'LNRECV_LOC',vacal,lnrecv_loc)
    
    call cputim(time2)
    time_comm = time_comm + (time2-time1)    
  
  end subroutine partition_sfc_redistribute_weights

  !-----------------------------------------------------------------------
  !> 
  !> @author  Ricard Borrell
  !> @date    2019-09-29
  !> @brief   Local bin partitioning
  !> @details Partition of the local bins by means of SFC
  !> 
  !-----------------------------------------------------------------------

  subroutine partition_sfc_local_partition

    weigc(iboxc0+1:iboxc1_+1) = weigc_loc(1:iboxc1_-iboxc0+1)

    if( weigc_loc(PAR_MY_RANK+1_ip) > 0 ) then
       call maths_sfc_par_part_2(sdim,dboxl,weigl,npart_sfc,dboxc,iboxc0+PAR_MY_RANK+1_ip,weigc,partl,lcorr,totalw,PART_NAME='PARTL')
    end if

  end subroutine partition_sfc_local_partition

  !-----------------------------------------------------------------------
  !> 
  !> @author  Ricard Borrell
  !> @date    2019-09-29
  !> @brief   Distribute result
  !> @details Distribute partitioning result
  !> 
  !-----------------------------------------------------------------------
  
  subroutine partition_sfc_redistribute_result()

    implicit none
    integer(ip)               :: ielem, ibuff,gbuffb
    character(100), parameter :: vacal = "partition_sfc_redistribute_result"
    real(rp)                  :: time1,time2
    !
    ! Here we use the receive buffer to send and the sent buffer to receive...
    !
    gbuffb = gbuff0
    do ibuff = 1_ip,comm1 % lrecv_dim/2_ip
       recv_buff(ibuff*2) = partl(recv_buff((ibuff-1)*2+1))
    end do
    call cputim(time1)
    call PAR_SEND_RECEIVE_TO_ALL(recv_buff,send_buff,comm2,'ASYNCHRONOUS')
    call cputim(time2)
    time_comm = time_comm + (time2-time1)
    !
    ! Assign elements to slave
    !
    do ielem = 1_ip, nenti
       if(lboxc(ielem) >= iboxc0+1_ip .and. lboxc(ielem) <= iboxc1_+1_ip) then
          lepar(ielem) = int(send_buff(lboxes(ielem)-gbuffb))
          gbuff0       = max(gbuff0,lboxes(ielem))
       endif
    end do
    !
    ! Deallocate communication pointers
    !
    call memory_deallo(par_memor,'SEND_BUFF'         ,vacal,send_buff         )
    call memory_deallo(par_memor,'RECV_BUFF'         ,vacal,recv_buff         )
    call memory_deallo(par_memor,'PARTL'             ,vacal,partl             )
    call memory_deallo(par_memor,'COMM1 % NEIGHTS'   ,vacal,comm1 % neights   )
    call memory_deallo(par_memor,'COMM2 % NEIGHTS'   ,vacal,comm2 % neights   )
    call memory_deallo(par_memor,'COMM1 % LSEND_SIZE',vacal,comm1 % lsend_size)
    call memory_deallo(par_memor,'COMM2 % LSEND_SIZE',vacal,comm2 % lsend_size)
    call memory_deallo(par_memor,'COMM1 % LRECV_SIZE',vacal,comm1 % lrecv_size)
    call memory_deallo(par_memor,'COMM2 % LRECV_SIZE',vacal,comm2 % lrecv_size)

  end subroutine partition_sfc_redistribute_result

  !-----------------------------------------------------------------------
  !> 
  !> @author  Ricard Borrell
  !> @date    2019-09-29
  !> @brief   Destructor
  !> @details Deallocate par_memory
  !> 
  !-----------------------------------------------------------------------

  subroutine partition_sfc_destructor()

    character(100), PARAMETER :: vacal = "par_dealocate_sfc"

    call memory_deallo(par_memor,'LBOXES'    ,vacal,lboxes)
    call memory_deallo(par_memor,'WEIGL'     ,vacal,weigl)
    call memory_deallo(par_memor,'WEIGC'     ,vacal,weigc)
    call memory_deallo(par_memor,'WEIGC_LOC' ,vacal,weigc_loc)
    call memory_deallo(par_memor,'BUFWEI'    ,vacal,bufwei)
    call memory_deallo(par_memor,'LNSEND'    ,vacal,lnsend)
    call memory_deallo(par_memor,'LBOXC'     ,vacal,lboxc)
    call memory_deallo(par_memor,'LBOXL'     ,vacal,lboxl)
    call memory_deallo(par_memor,'NEIGHS'    ,vacal,neighs)
    call memory_deallo(par_memor,'BUFSE2'    ,vacal,bufse2)
    if(.not. if_lcorr) then
       call memory_deallo(par_memor,'LCORR'     ,vacal,lcorr)
    else
       nullify(lcorr)
    endif
    if(.not. if_lweig) then
       call memory_deallo(par_memor,'LWEIG'     ,vacal,lweig)
    else
       nullify(lweig)
    endif
    call memory_deallo(par_memor,'REORD'     ,vacal,reord)
    call memory_deallo(par_memor,'INVOR'     ,vacal,invor)

    call PAR_BARRIER(COMM=PAR_COMM)
    call cputim(time_partition)
    time_partition = time_partition-time_ini

  end subroutine partition_sfc_destructor

end module mod_partition_sfc
!> @}
