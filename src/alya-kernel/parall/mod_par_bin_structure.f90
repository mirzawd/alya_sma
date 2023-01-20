!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Parall
!> @{
!> @file    mod_par_bin_structure.f90
!> @author  houzeaux
!> @date    2018-04-13
!> @brief   Subdomain bin structure
!> @details ???
!>          ???
!-----------------------------------------------------------------------

module mod_par_bin_structure

  use def_kintyp,            only :  ip,rp,i1p,lg
  use def_master,            only :  IPARALL,ISEQUEN,IMASTER
  use def_master,            only :  ioutp,lun_outpu,kfl_paral
  use def_master,            only :  npart
  use def_domain,            only :  ndime
  use mod_optional_argument, only :  optional_argument
  use mod_communications,    only :  PAR_MIN
  use mod_communications,    only :  PAR_MAX
  use mod_communications,    only :  PAR_SUM
  use mod_communications,    only :  PAR_COMM_RANK_AND_SIZE
  use mod_communications,    only :  PAR_ALLGATHERV
  use mod_communications,    only :  PAR_ALLGATHER
  use mod_communications,    only :  PAR_DEFINE_COMMUNICATOR
  use mod_parall,            only :  par_memor
  use mod_parall,            only :  PAR_WORLD_SIZE
  use mod_parall,            only :  PAR_MY_CODE_RANK
  use mod_parall,            only :  PAR_COMM_WORLD
  use mod_memory,            only :  memory_alloca
  use mod_memory,            only :  memory_deallo
  use mod_memory,            only :  memory_size
  use mod_maths,             only :  maths_mapping_3d_to_1d
  use mod_maths,             only :  maths_mapping_1d_to_3d_x
  use mod_maths,             only :  maths_mapping_1d_to_3d_y
  use mod_maths,             only :  maths_mapping_1d_to_3d_z
  use mod_maths,             only :  maths_mapping_coord_to_3d
  use mod_maths,             only :  maths_in_box
  use mod_parall,            only :  typ_bin_structure
  use mod_outfor,            only :  outfor
  use mod_messages,          only :  messages_live
  use def_mpi
#include "def_mpi.inc"
  implicit none

  private

  integer(ip),            parameter :: number_bins=20_ip
  real(rp),               parameter :: eps=epsilon(1.0_rp)
  type(typ_bin_structure)           :: bin_structure_par
  
  interface par_bin_structure
     module procedure &
          &    par_bin_structure_output,&
          &    par_bin_structure_empty,&
          &    par_bin_structure_type
  end interface par_bin_structure

  interface par_bin_structure_deallocate 
     module procedure &
          &    par_bin_structure_deallocate_type,&
          &    par_bin_structure_deallocate_empty
  end interface par_bin_structure_deallocate

  interface par_bin_structure_initialization 
     module procedure &
          &    par_bin_structure_initialization_type,&
          &    par_bin_structure_initialization_empty
  end interface par_bin_structure_initialization

  public :: par_bin_structure
  public :: par_bin_structure_initialization
  public :: par_bin_structure_deallocate
  public :: par_bin_structure_partition_bounding_box
  public :: par_bin_structure_points_in_partitions
  
contains

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    01/02/2014
  !> @brief   Initialize bin structure
  !> @details Initialize bin structure
  !>
  !----------------------------------------------------------------------

  subroutine par_bin_structure_deallocate_type(bin_structure)

    type(typ_bin_structure), intent(inout) :: bin_structure

    call memory_deallo(par_memor,'PAR_BIN_SIZE'  ,'par_bin_structure_deallocate',bin_structure % size)
    call memory_deallo(par_memor,'PAR_BIN_PART'  ,'par_bin_structure_deallocate',bin_structure % part)
    call memory_deallo(par_memor,'PAR_PART_COMIN','par_bin_structure_deallocate',bin_structure % part_comin)
    call memory_deallo(par_memor,'PAR_PART_COMAX','par_bin_structure_deallocate',bin_structure % part_comax)

  end subroutine par_bin_structure_deallocate_type

  subroutine par_bin_structure_deallocate_empty()

    use mod_parall,         only :  par_bin_part
    use mod_parall,         only :  par_bin_boxes
    use mod_parall,         only :  par_bin_comin
    use mod_parall,         only :  par_bin_comax
    use mod_parall,         only :  par_bin_size
    use mod_parall,         only :  par_part_comin
    use mod_parall,         only :  par_part_comax

    call memory_deallo(par_memor,'PAR_BIN_SIZE'  ,'par_bin_structure_deallocate',par_bin_size)
    call memory_deallo(par_memor,'PAR_BIN_PART'  ,'par_bin_structure_deallocate',par_bin_part)
    call memory_deallo(par_memor,'PAR_PART_COMIN','par_bin_structure_deallocate',par_part_comin)
    call memory_deallo(par_memor,'PAR_PART_COMAX','par_bin_structure_deallocate',par_part_comax)

    par_bin_boxes =  0_ip                       
    par_bin_comin = -huge(1.0_rp)*0.1_rp ! To avoid overflow when doing comax-comin              
    par_bin_comax =  huge(1.0_rp)*0.1_rp        

  end subroutine par_bin_structure_deallocate_empty

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    01/02/2014
  !> @brief   Initialize bin structure
  !> @details Initialize bin structure
  !>
  !----------------------------------------------------------------------

  subroutine par_bin_structure_initialization_type(bin_structure)

    type(typ_bin_structure), intent(inout) :: bin_structure

    bin_structure % boxes =  0_ip                       
    bin_structure % comin = -huge(1.0_rp)*0.1_rp ! To avoid overflow when doing comax-comin              
    bin_structure % comax =  huge(1.0_rp)*0.1_rp        
    nullify(bin_structure % size)            
    nullify(bin_structure % part)                      
    nullify(bin_structure % part_comin)                
    nullify(bin_structure % part_comax)                

  end subroutine par_bin_structure_initialization_type

  subroutine par_bin_structure_initialization_empty()

    use mod_parall,         only :  par_bin_part
    use mod_parall,         only :  par_bin_boxes
    use mod_parall,         only :  par_bin_comin
    use mod_parall,         only :  par_bin_comax
    use mod_parall,         only :  par_bin_size
    use mod_parall,         only :  par_part_comin
    use mod_parall,         only :  par_part_comax

    par_bin_boxes =  1
    par_bin_comin = -huge(1.0_rp)*0.1_rp
    par_bin_comax =  huge(1.0_rp)*0.1_rp       
    nullify(par_bin_size)
    nullify(par_bin_part)
    nullify(par_part_comin)
    nullify(par_part_comax)

    call par_bin_structure_initialization_type(bin_structure_par)
 
  end subroutine par_bin_structure_initialization_empty

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    01/02/2014
  !> @brief   Subdomain bin structure
  !> @details Subdomain bin structure with type as argument
  !>
  !----------------------------------------------------------------------

  subroutine par_bin_structure_type(bin_structure,comin,comax,ABSOLUTE_TOLERANCE,COMM,VERBOSE)

    type(typ_bin_structure), intent(inout)          :: bin_structure      !< Bin structure
    real(rp),                intent(in),   optional :: comin(3)           !< Coordinates
    real(rp),                intent(in),   optional :: comax(3)           !< Coordinates
    real(rp),                intent(in),   optional :: ABSOLUTE_TOLERANCE !< Tolerance (same units as coord)
    MY_MPI_COMM   ,          intent(in),   optional :: COMM               !< Communicator
    logical(lg),             intent(in),   optional :: VERBOSE            !< Output information

    call par_bin_structure_output(&
         bin_structure % boxes,bin_structure % comin,bin_structure % comax,&
         bin_structure % size,bin_structure % part,bin_structure % part_comin,&
         bin_structure % part_comax,comin,comax,ABSOLUTE_TOLERANCE,&
         COMM,VERBOSE)

  end subroutine par_bin_structure_type

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    01/02/2014
  !> @brief   Subdomain bin structure
  !> @details Subdomain bin structure with no arguments: thus
  !>          choose Alya's global subdomains
  !>
  !----------------------------------------------------------------------

  subroutine par_bin_structure_empty()

    use mod_parall,         only :  par_bin_part
    use mod_parall,         only :  par_bin_boxes
    use mod_parall,         only :  par_bin_comin
    use mod_parall,         only :  par_bin_comax
    use mod_parall,         only :  par_bin_size
    use mod_parall,         only :  par_part_comin
    use mod_parall,         only :  par_part_comax

    call par_bin_structure_output(&
         par_bin_boxes,par_bin_comin,par_bin_comax,&
         par_bin_size,par_bin_part,par_part_comin,&
         par_part_comax,VERBOSE=.true.)

    bin_structure_par % boxes      =   par_bin_boxes     
    bin_structure_par % comin      =   par_bin_comin 
    bin_structure_par % comax      =   par_bin_comax
    bin_structure_par % size       =>  par_bin_size   
    bin_structure_par % part       =>  par_bin_part              
    bin_structure_par % part_comin =>  par_part_comin              
    bin_structure_par % part_comax =>  par_part_comax               
    
  end subroutine par_bin_structure_empty

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    01/02/2014
  !> @brief   Subdomain bin structure
  !> @details Compute the parallel bin structure
  !>          Limit the bin size to 100 in each dimension
  !>
  !----------------------------------------------------------------------

  subroutine par_bin_structure_output(&
       par_bin_boxes,par_bin_comin,par_bin_comax,&
       par_bin_size,par_bin_part,par_part_comin,&
       par_part_comax,comin_opt,comax_opt,&
       ABSOLUTE_TOLERANCE,COMM,VERBOSE)

    integer(ip), intent(out)                     :: par_bin_boxes(3)     !< # boxes in each direction
    real(rp),    intent(out)                     :: par_bin_comin(3)     !< Minimum box coordinates
    real(rp),    intent(out)                     :: par_bin_comax(3)     !< Maximum box coordinates
    integer(ip), intent(inout),pointer           :: par_bin_size(:,:,:)  !< # partitions per box
    type(i1p),   intent(inout),pointer           :: par_bin_part(:,:,:)  !< Bin structure of world partition
    real(rp),    intent(inout),pointer           :: par_part_comin(:,:)  !< Subdomain minimum coordinates
    real(rp),    intent(inout),pointer           :: par_part_comax(:,:)  !< Subdomain maximum coordinates
    real(rp),    intent(in),   optional          :: comin_opt(3)         !< Subdomain minimum coordinates
    real(rp),    intent(in),   optional          :: comax_opt(3)         !< Subdomain maximum coordinates
    real(rp),    intent(in),   optional          :: ABSOLUTE_TOLERANCE   !< Tolerance (same units as coord)
    MY_MPI_COMM   , intent(in),   optional          :: COMM                 !< Communicator
    logical(lg), intent(in),   optional          :: VERBOSE              !< Output information

    integer(ip)                                  :: PAR_CURRENT_SIZE     ! Size of the communicator
    integer(ip)                                  :: PAR_CURRENT_RANK
    MY_MPI_COMM                                  :: PAR_CURRENT_COMM4
    integer(ip)                                  :: idime
    integer(ip)                                  :: nx,ny,nz             ! Number if bins in each directions
    integer(ip)                                  :: imin,imax
    integer(ip)                                  :: jmin,jmax
    integer(ip)                                  :: kmin,kmax
    integer(ip)                                  :: ii,jj,kk,iboxe
    integer(ip)                                  :: iboxe1,iboxe2
    integer(ip)                                  :: nxny
    integer(ip)                                  :: nnbox,nzbox,kboxe
    integer(ip)                                  :: ipart
    integer(ip)                                  :: mdime
    real(rp)                                     :: dnx,dny,dnz
    real(rp)                                     :: delta(3)             ! Bin size
    real(rp)                                     :: comin(3),comax(3)    ! Bounding box
    real(rp)                                     :: xmibo(3)             ! Minimum box size
    real(rp)                                     :: deltx,delty,deltz
    integer(ip), pointer                         :: npart_per_box(:)
    integer(ip), pointer                         :: list_boxes(:)
    integer(ip), pointer                         :: my_list_boxes(:)
    integer(ip), pointer                         :: number_boxes(:)
    integer(4),  pointer                         :: recvcounts4(:)
    real(rp)                                     :: toler
    logical(lg)                                  :: ifverbose
    !
    ! VERY IMPORTANT PARAMETER: tolerance added to the bounding boxes of the subdomains
    !
    toler = 0.01_rp
    deltx = 0.0_rp
    delty = 0.0_rp
    deltz = 0.0_rp
    !
    ! Initialize
    !
    nullify( npart_per_box )
    nullify( list_boxes    )
    nullify( my_list_boxes )
    nullify( number_boxes  )
    nullify( recvcounts4   )    
    !
    ! Communicator
    !
    if( present(COMM) ) then
       call PAR_COMM_RANK_AND_SIZE(COMM,PAR_CURRENT_RANK,PAR_CURRENT_SIZE)
       PAR_CURRENT_COMM4 = COMM
    else 
       PAR_CURRENT_COMM4 = PAR_COMM_WORLD
       PAR_CURRENT_SIZE  = PAR_WORLD_SIZE
    end if

    mdime = ndime
    call PAR_MAX(mdime,PAR_CURRENT_COMM4)

    ifverbose = optional_argument(.true.,VERBOSE)
    if( ifverbose ) call messages_live('COMPUTE BIN STRUCTURE OF PARTITIONS')
    !
    ! COMIN and COMAX: My bounding box
    ! XMIBO: minimum partition size in each dimension
    !
    if( present(comin_opt) .and. present(comax_opt) ) then
       comin          = 0.0_rp
       comax          = 0.0_rp
       comin(1:mdime) = comin_opt(1:mdime)
       comax(1:mdime) = comax_opt(1:mdime)       
       xmibo(1:mdime) = abs(comax(1:mdime) - comin(1:mdime))
    else
       call par_bin_structure_partition_bounding_box(comin,comax,xmibo,toler,ABSOLUTE_TOLERANCE)
    end if
    !
    ! All gather the bounding boxes
    !
    call memory_alloca(par_memor,'PAR_PART_COMIN','par_bin_structure',par_part_comin,3_ip,PAR_CURRENT_SIZE,'INITIALIZE',1_ip,0_ip)   
    call memory_alloca(par_memor,'PAR_PART_COMAX','par_bin_structure',par_part_comax,3_ip,PAR_CURRENT_SIZE,'INITIALIZE',1_ip,0_ip)
    call PAR_ALLGATHER(comin,par_part_comin,3_4,PAR_COMM_IN=PAR_CURRENT_COMM4)
    call PAR_ALLGATHER(comax,par_part_comax,3_4,PAR_COMM_IN=PAR_CURRENT_COMM4)
    !
    ! Minimum subdomain size
    !
    ! XMIBO(1:MDIME)         = minimum box size in each direction 
    ! PAR_BIN_COMIN(1:MDIME) = minimum position of all-domain bounding box
    ! PAR_BIN_COMAX(1:MDIME) = maximum position of all-domain bounding box
    !
    call PAR_MIN(mdime,xmibo(1:mdime),        PAR_CURRENT_COMM4)
    par_bin_comin(1:mdime) = comin(1:mdime)
    call PAR_MIN(mdime,par_bin_comin(1:mdime),PAR_CURRENT_COMM4)
    par_bin_comax(1:mdime) = comax(1:mdime)
    call PAR_MAX(mdime,par_bin_comax(1:mdime),PAR_CURRENT_COMM4)
    !
    ! Size of bin
    ! DELTA(1:MDIME)         = size of boxes in each dimension
    ! PAR_BIN_BOXES(1:MDIME) = nb of boxes in each dimension
    !
    par_bin_comin = par_bin_comin - eps
    par_bin_comax = par_bin_comax + eps
    par_bin_boxes = 10

    do idime = 1,mdime
       delta(idime)         = (par_bin_comax(idime)-par_bin_comin(idime) ) / max(xmibo(idime),eps)
       if( delta(idime) >= real(huge(1_ip),rp) ) then 
          par_bin_boxes(idime) = number_bins
       else
          par_bin_boxes(idime) = min(number_bins,max(1_ip,int(delta(idime),ip)))
       end if
       delta(idime)        = (par_bin_comax(idime)-par_bin_comin(idime))/real(par_bin_boxes(idime),rp)
    end do
    !
    ! Fill in bin structure
    !
    nx    = par_bin_boxes(1)
    ny    = par_bin_boxes(2)
    nz    = par_bin_boxes(3)
    dnx   = real(nx,rp)
    dny   = real(ny,rp)
    dnz   = real(nz,rp)
    nnbox = nx*ny*nz

    allocate( npart_per_box(nnbox) )
    do iboxe = 1,nnbox
       npart_per_box(iboxe) = 0
    end do

    !----------------------------------------------------------------------
    !
    ! Fill in boxes with my subdomain
    ! The second conditions tests the case with empty subdomain (comin>comax)
    !
    !----------------------------------------------------------------------      

    if( PAR_MY_CODE_RANK /= 0  .and. comin(1) <= comax(1) ) then

       deltx = 1.0_rp / ( par_bin_comax(1)-par_bin_comin(1) )
       if( mdime >= 2 ) then
          delty = 1.0_rp / ( par_bin_comax(2)-par_bin_comin(2) )
       end if
       if( mdime == 3 ) then
          deltz = 1.0_rp / ( par_bin_comax(3)-par_bin_comin(3) )
       end if
       if( mdime >= 2 ) nxny  = par_bin_boxes(1) * par_bin_boxes(2)

       if( mdime == 1 ) then

          imin = int(( ( comin(1) - par_bin_comin(1) ) * deltx ) * dnx , ip ) + 1
          imax = int(( ( comax(1) - par_bin_comin(1) ) * deltx ) * dnx , ip ) + 1  

          imin = max(imin,1_ip)
          imax = min(imax,par_bin_boxes(1))

          do ii = imin,imax
             iboxe = ii
             npart_per_box(iboxe) = 1
          end do

       else if( mdime == 2 ) then

          imin = int(( ( comin(1) - par_bin_comin(1) ) * deltx ) * dnx , ip ) + 1
          imax = int(( ( comax(1) - par_bin_comin(1) ) * deltx ) * dnx , ip ) + 1      

          jmin = int(( ( comin(2) - par_bin_comin(2) ) * delty ) * dny , ip ) + 1
          jmax = int(( ( comax(2) - par_bin_comin(2) ) * delty ) * dny , ip ) + 1      

          imin = max(imin,1_ip)
          imax = min(imax,par_bin_boxes(1))
          jmin = max(jmin,1_ip)
          jmax = min(jmax,par_bin_boxes(2))

          do ii = imin,imax
             do jj = jmin,jmax
                iboxe = (jj-1_ip) * nx + ii
                npart_per_box(iboxe) = 1
             end do
          end do

       else if( mdime == 3 ) then

          imin = int(( ( comin(1) - par_bin_comin(1) ) * deltx ) * dnx , ip ) + 1
          imax = int(( ( comax(1) - par_bin_comin(1) ) * deltx ) * dnx , ip ) + 1      

          jmin = int(( ( comin(2) - par_bin_comin(2) ) * delty ) * dny , ip ) + 1
          jmax = int(( ( comax(2) - par_bin_comin(2) ) * delty ) * dny , ip ) + 1      

          kmin = int(( ( comin(3) - par_bin_comin(3) ) * deltz ) * dnz , ip ) + 1
          kmax = int(( ( comax(3) - par_bin_comin(3) ) * deltz ) * dnz , ip ) + 1      

          imin = max(imin,1_ip)
          imax = min(imax,par_bin_boxes(1))
          jmin = max(jmin,1_ip)
          jmax = min(jmax,par_bin_boxes(2))
          kmin = max(kmin,1_ip)
          kmax = min(kmax,par_bin_boxes(3))

          do kk = kmin,kmax
             iboxe2 = nxny * (kk-1)
             do jj = jmin,jmax
                iboxe1 = iboxe2 + par_bin_boxes(1) * (jj-1)
                do ii = imin,imax
                   iboxe = iboxe1 + ii
                   npart_per_box(iboxe) = 1
                end do
             end do
          end do

       end if
    end if
    !
    ! All gather NUMBER_BOXES 
    ! 
    nzbox = 0  
    if( PAR_MY_CODE_RANK /= 0 ) then
       do iboxe = 1,nnbox
          nzbox = nzbox + npart_per_box(iboxe)
       end do
    end if
    if( nzbox > 0 ) then
       allocate( my_list_boxes(nzbox) )
       nzbox = 0  
       do iboxe = 1,nnbox
          if( npart_per_box(iboxe) == 1 ) then
             nzbox = nzbox + 1
             my_list_boxes(nzbox) = iboxe
          end if
       end do
    else
       nullify( my_list_boxes )
    end if
    allocate ( number_boxes(0:PAR_CURRENT_SIZE-1) )
    call PAR_ALLGATHER(nzbox,number_boxes,1_4,PAR_COMM_IN=PAR_CURRENT_COMM4)
    !
    ! NZBOX = TOTAL NUMBER OF BOXES
    !
    nzbox = 0
    do ipart = 0,PAR_CURRENT_SIZE-1
       nzbox = nzbox + number_boxes(ipart)
    end do
    !
    ! LIST_BOXES = All gather list of boxes
    !
    if( nzbox > 0 ) allocate( list_boxes(nzbox) )
    do iboxe = 1,nzbox
       list_boxes(iboxe) = 0
    end do
    allocate( recvcounts4(0:PAR_CURRENT_SIZE-1) )  
    do ipart = 0,PAR_CURRENT_SIZE-1
       recvcounts4(ipart) = int(number_boxes(ipart),4)
    end do
    call PAR_ALLGATHERV(my_list_boxes,list_boxes,recvcounts4,PAR_COMM_IN=PAR_CURRENT_COMM4)
    deallocate( recvcounts4 )

    !----------------------------------------------------------------------
    !
    ! Fill in bin
    ! PAR_BIN_SIZE(I,J,K) = Number of partitions
    ! PAR_BIN_PART(I,J,K) % L = List of partitions in MPI_COMM_WORLD
    !
    !----------------------------------------------------------------------

    if( .not. associated(par_bin_part) ) &
         call memory_alloca(par_memor,'PAR_BIN_PART','par_bin_structure',par_bin_part,nx,ny,nz)
    if( .not. associated(par_bin_size) ) &
         call memory_alloca(par_memor,'PAR_BIN_SIZE','par_bin_structure',par_bin_size,nx,ny,nz)

    par_bin_size = 0
    nzbox = 0
    do ipart = 0,PAR_CURRENT_SIZE-1
       do iboxe = 1,number_boxes(ipart)
          nzbox = nzbox + 1
          kboxe = list_boxes(nzbox)
          ii    = min(maths_mapping_1d_to_3d_x(nx,ny,nz,kboxe),nx)
          jj    = min(maths_mapping_1d_to_3d_y(nx,ny,nz,kboxe),ny)
          kk    = min(maths_mapping_1d_to_3d_z(nx,ny,nz,kboxe),nz)
          par_bin_size(ii,jj,kk) = par_bin_size(ii,jj,kk) + 1
       end do
    end do

    do kk = 1,nz
       do jj = 1,ny
          do ii = 1,nx
             if( par_bin_size(ii,jj,kk) > 0 ) then
                if( associated(par_bin_part(ii,jj,kk) % l) ) then
                   call memory_deallo(par_memor,'PAR_BIN_PART % L','par_bin_structure',par_bin_part(ii,jj,kk) % l)
                end if
                call memory_alloca(par_memor,'PAR_BIN_PART % L','par_bin_structure',par_bin_part(ii,jj,kk) % l,par_bin_size(ii,jj,kk) )
                par_bin_size(ii,jj,kk) = 0
             else
                nullify( par_bin_part(ii,jj,kk) % l )
             end if
          end do
       end do
    end do

    nzbox = 0
    do ipart = 0,PAR_CURRENT_SIZE-1
       do iboxe = 1,number_boxes(ipart)
          nzbox = nzbox + 1
          kboxe = list_boxes(nzbox)
          ii    = min(maths_mapping_1d_to_3d_x(nx,ny,nz,kboxe),nx)
          jj    = min(maths_mapping_1d_to_3d_y(nx,ny,nz,kboxe),ny)
          kk    = min(maths_mapping_1d_to_3d_z(nx,ny,nz,kboxe),nz)
          par_bin_size(ii,jj,kk) = par_bin_size(ii,jj,kk) + 1
          par_bin_part(ii,jj,kk) % l(par_bin_size(ii,jj,kk)) = ipart
       end do
    end do
    !
    ! Deallocate memory
    !
    if( associated(npart_per_box) ) deallocate( npart_per_box )
    if( associated(list_boxes)    ) deallocate( list_boxes    )
    if( associated(my_list_boxes) ) deallocate( my_list_boxes )
    if( associated(number_boxes)  ) deallocate( number_boxes  )
    !
    ! Output bin
    !   
    if( ifverbose ) then
       ioutp(1) = par_bin_boxes(1) 
       ioutp(2) = par_bin_boxes(2) 
       ioutp(3) = par_bin_boxes(3) 
       call outfor(68_ip,lun_outpu,'')
    end if
    !
    ! Code/rank in each bin (II,JJ,KK)
    !
    !if( PAR_MY_WORLD_RANK == 0 ) then
    !   do ii=1,nx;do jj=1,ny ; do kk=1,nz
    !      if( par_bin_size(ii,jj,kk) > 0 ) then
    !         do iboxe = 1,par_bin_size(ii,jj,kk)
    !            ipart = par_bin_part(ii,jj,kk) % l(iboxe)
    !            icode = PAR_COMM_WORLD_TO_CODE_PERM(1,ipart)
    !            ipart = PAR_COMM_WORLD_TO_CODE_PERM(2,ipart)
    !         end do
    !      end if
    !   end do; end do; end do
    !end if
    !call runend('O.K.!')
    !end if ! GUIGUI

  end subroutine par_bin_structure_output

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-04-13
  !> @brief   Compute partition bounding box
  !> @details Compute partition bounding box up to a tolerance
  !>          An absolute tolerance can be given. Then the bounding
  !>          box is increased on both side by a factor of:
  !>          ABSOLUTE_TOLERANCE <  0 : -ABSOLUTE_TOLERANCE * Vol_ave
  !>                             >= 0 : ABSOLUTE_TOLERANCE
  !>          where Vol_ave is the average element volume in the mesh
  !> 
  !-----------------------------------------------------------------------

  subroutine par_bin_structure_partition_bounding_box(&
       comin,comax,delta,RELATIVE_TOLERANCE,ABSOLUTE_TOLERANCE)

    use def_domain, only : npoin,coord,voave

    real(rp),    intent(out)          :: comin(3)           !< Minimum coordinates of bounding box
    real(rp),    intent(out)          :: comax(3)           !< Maximum coordinates of bounding box
    real(rp),    intent(out)          :: delta(3)           !< Size in each direction
    real(rp),    intent(in), optional :: RELATIVE_TOLERANCE !< Relative tolerance
    real(rp),    intent(in), optional :: ABSOLUTE_TOLERANCE !< Tolerance (same units as coord)
    integer(ip)                       :: idime
    real(rp)                          :: toler_abs
    real(rp)                          :: toler_rel
    !
    ! Tolerances
    !
    toler_abs = 0.0_rp
    toler_rel = 0.01_rp    
    if( present(RELATIVE_TOLERANCE) ) toler_rel = RELATIVE_TOLERANCE
    if( present(ABSOLUTE_TOLERANCE) ) then
       if( ABSOLUTE_TOLERANCE < 0.0_rp ) then
          if( voave > 0.0_rp ) then
             toler_abs = -ABSOLUTE_TOLERANCE * voave**(1.0_rp/real(ndime,rp))
          else
             toler_abs = 0.0_rp
          end if
       else
          toler_abs =  ABSOLUTE_TOLERANCE
       end if
    end if
    !
    ! COMIN and COMAX: My bounding box
    !    
    comin(1:3) =  huge(1.0_rp)*0.1_rp
    comax(1:3) = -huge(1.0_rp)*0.1_rp    
    if( PAR_MY_CODE_RANK /= 0 ) then
       if( npoin > 0 ) then
          do idime = 1,ndime
             comin(idime) = minval(coord(idime,1:npoin))
             comax(idime) = maxval(coord(idime,1:npoin))
          end do
       end if
       delta(1:ndime) = abs(comax(1:ndime) - comin(1:ndime))
       comin(1:ndime) = comin(1:ndime) - toler_abs - toler_rel * delta(1:ndime)
       comax(1:ndime) = comax(1:ndime) + toler_abs + toler_rel * delta(1:ndime)
    end if

  end subroutine par_bin_structure_partition_bounding_box

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-07-07
  !> @brief   Points in partitions
  !> @details Return the list of points partitions
  !> 
  !-----------------------------------------------------------------------
  
  subroutine par_bin_structure_points_in_partitions(xx,npoin_in_part,lpoin_in_part,NUMBER_POINTS,INCLUDE_MYSELF,bin_structure)

    real(rp),                          pointer, intent(in)    :: xx(:,:)
    integer(ip),                       pointer, intent(inout) :: npoin_in_part(:)
    type(i1p),                         pointer, intent(inout) :: lpoin_in_part(:)
    logical(lg),             optional,          intent(in)    :: INCLUDE_MYSELF
    integer(ip),             optional,          intent(in)    :: NUMBER_POINTS
    type(typ_bin_structure), optional,          intent(in)    :: bin_structure
    integer(ip)                                               :: kpart,cpart,pp
    integer(ip)                                               :: ii,jj,kk,nn,ksize
    integer(ip),                       pointer                :: aux(:)
    type(i1p),                         pointer                :: lpart_in_poin(:)
    logical(lg)                                               :: if_include_myself

    nullify(aux,lpart_in_poin)

    nn                = memory_size(xx,2_ip)
    nn                = optional_argument(nn,NUMBER_POINTS)    
    if_include_myself = optional_argument(.false.,INCLUDE_MYSELF)    
    !
    ! Allocate
    !
    if( .not. associated(npoin_in_part) )&
         call memory_alloca(par_memor,'NPOIN_IN_PART','par_bin_structure_points_in_partitions',npoin_in_part,npart+1,lboun=0_ip)
    if( .not. associated(lpoin_in_part) )&
         call memory_alloca(par_memor,'LPOIN_IN_PART','par_bin_structure_points_in_partitions',lpoin_in_part,npart+1,lboun=0_ip)

    call memory_alloca(par_memor,'LPART_IN_POIN','par_bin_structure_points_in_partitions',lpart_in_poin,nn)
    call memory_alloca(par_memor,'AUX','par_bin_structure_points_in_partitions',aux,npart+1)
    !
    ! NPOIN_IN_PART(CPART) = number of points in CPART
    ! LPART_IN_POIN(PP)    = list of partitions for point PP
    !
    if( present(bin_structure) ) then
       do pp = 1,nn
          call maths_mapping_coord_to_3d(ndime,bin_structure % boxes,bin_structure % comin,bin_structure % comax,xx(1:ndime,pp),ii,jj,kk)
          if( ii /= 0 ) then
             ksize = 0
             do kpart = 1,bin_structure % size(ii,jj,kk)
                cpart = bin_structure % part(ii,jj,kk) % l(kpart)
                if( cpart /= kfl_paral .or. if_include_myself ) then
                   if( maths_in_box(ndime,xx(1:ndime,pp),bin_structure % part_comin(1:ndime,cpart),bin_structure % part_comax(1:ndime,cpart)) ) then
                      npoin_in_part(cpart) = npoin_in_part(cpart) + 1
                      ksize                = ksize + 1                   
                      aux(ksize)           = cpart
                   end if
                end if
             end do
             if( ksize > 0 ) then
                call memory_alloca(par_memor,'LPART_IN_POIN % L','par_bin_structure_points_in_partitions',lpart_in_poin(pp) % l,ksize)
                lpart_in_poin(pp) % l(1:ksize) = aux(1:ksize)
             end if
          end if
       end do
    else
       do pp = 1,nn
          call maths_mapping_coord_to_3d(ndime,bin_structure_par % boxes,bin_structure_par % comin,bin_structure_par % comax,xx(1:ndime,pp),ii,jj,kk)
          if( ii /= 0 ) then
             ksize = 0
             do kpart = 1,bin_structure_par % size(ii,jj,kk)
                cpart = bin_structure_par % part(ii,jj,kk) % l(kpart)
                if( cpart /= kfl_paral .or. if_include_myself ) then
                   if( maths_in_box(ndime,xx(1:ndime,pp),bin_structure_par % part_comin(1:ndime,cpart),bin_structure_par % part_comax(1:ndime,cpart)) ) then
                      npoin_in_part(cpart) = npoin_in_part(cpart) + 1
                      ksize                = ksize + 1                   
                      aux(ksize)           = cpart
                   end if
                end if
             end do
             if( ksize > 0 ) then
                call memory_alloca(par_memor,'LPART_IN_POIN % L','par_bin_structure_points_in_partitions',lpart_in_poin(pp) % l,ksize)
                lpart_in_poin(pp) % l(1:ksize) = aux(1:ksize)
             end if
          end if
       end do
    end if
    !
    ! LPOIN_IN_PART(CPART) = list of points of partition CPART
    !
    aux = 0_ip
    do pp = 1,nn
       do ii = 1,memory_size(lpart_in_poin(pp) % l)
          cpart      = lpart_in_poin(pp) % l(ii)
          aux(cpart) = aux(cpart) + 1
       end do
    end do
    do cpart = 1,npart
       call memory_alloca(par_memor,'LPOIN_IN_PART % L','par_bin_structure_points_in_partitions',lpoin_in_part(cpart) % l,aux(cpart))      
       aux(cpart) = 0
    end do
    do pp = 1,nn
       do ii = 1,memory_size(lpart_in_poin(pp) % l)
          cpart                                = lpart_in_poin(pp) % l(ii)
          aux(cpart)                           = aux(cpart) + 1
          lpoin_in_part(cpart) % l(aux(cpart)) = pp
       end do
    end do
    call memory_deallo(par_memor,'AUX','par_bin_structure_points_in_partitions',aux)
    call memory_deallo(par_memor,'LPART_IN_POIN','par_bin_structure_points_in_partitions',lpart_in_poin)

  end subroutine par_bin_structure_points_in_partitions

end module mod_par_bin_structure

  
