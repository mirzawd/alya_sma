!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Adaptivity
!> @{
!> @file    mod_adapt.f90
!> @author  abel.gargallo
!> @date    2021-03-31
!> @brief   mod_adapt
!> @details mod_adapt
!>
!>          To add further details
!>
!>
!-----------------------------------------------------------------------

MODULE mod_adapt
!***************************************************************
!*
!*  Module for performing adaptive topological operation
!*
!***************************************************************
use def_kintyp_basic,       only: ip,rp,lg
use def_kintyp_mesh_basic,  only: mesh_type_basic
use mod_memory,             only: memory_alloca, memory_deallo
!
use def_adapt,              only: memor_adapt 
use mod_metric,             only: mesh_metric_type
use mod_debugTools,         only: out_debug_text, out_debug_paraview
use mod_debugTools,         only: out_performance, time_total, time_split, time_coll, time_swap, time_smooth
use mod_debugTools,         only: out_performance_profile, out_performance_profile_global
use mod_meshTopology,       only: edgeData_type
!
implicit none
!
integer(ip), parameter :: adapt_strat_split    = 1
integer(ip), parameter :: adapt_strat_collapse = 2
integer(ip), parameter :: adapt_strat_swap     = 3
integer(ip), parameter :: adapt_strat_smooth   = 4
!
real(rp),    parameter :: minQ_validity   = 1.0e-6! 0.00001_rp!1.0e-5
logical(lg), parameter :: do_warnValidity = .true.
!
logical(lg), parameter :: output_steps_paraview       = .true.
logical(lg), parameter :: output_steps_swap_paraview  = .false.
logical(lg), parameter :: output_steps_edgeRep_paraview  = .true.
!
logical(lg), parameter :: out_verbose = .true.
logical(lg), parameter :: out_time    = .true.

character(len=200), parameter :: file_out =TRIM('mesh_adapt')

! ---- ADAPTIVITY PARAMETERS ----------------------------
!
! In adapt_mesh_to_metric: quality to set the subregion in where the mesh is adapted (below threshold)
!
real(rp)  :: threshold_quality_to_repair 
integer(ip), parameter:: numNieghLevels = 2_ip
!
! In adapt_mesh_to_metric_private:
!
integer(ip),parameter :: numRepairs = 4
  !real(rp)   ,parameter:: tols_coll(numRepairs)   = (/    0.5,   0.5, 1.0/3.0, 1.0/3.0 /)
!   logical(lg),parameter:: doCheckColl(numRepairs) = (/.false.,.false.,  .false.,  .false./) ! check improvement from prev configuration
  !logical(lg),parameter:: doCheckColl(numRepairs) = (/.true.,.true.,  .true.,  .true./) ! check improvement from prev configuration
real(rp)   ,parameter:: tols_coll(numRepairs)   = (/     0.1,     0.2,     0.2,     0.4 /)
! real(rp)   ,parameter:: tols_coll(numRepairs)   = (/     0.00001,     0.00002,     0.00002,     0.00004 /)
logical(lg),parameter:: doCheckColl(numRepairs) = (/  .true.,  .true.,  .true.,  .true. /) ! check improvement from prev configuration
! logical(lg),parameter:: doCheckColl(numRepairs) = (/  .false.,  .false.,  .true.,  .true. /) ! check improvement from prev configuration
  
real(rp)   ,parameter:: tols_split(numRepairs)  = (/ 10.0,    3.0,     2.0,     1.5/) !sqrt(2.0)
! logical(lg),parameter:: doCheckSplit(numRepairs)= (/.false.,.false.,  .true.,  .true./) ! check improvement from prev configuration
logical(lg),parameter:: doCheckSplit(numRepairs)= (/.true.,.true.,  .true.,  .true./) ! check improvement from prev configuration
! logical(lg),parameter:: doCheckSplit(numRepairs)= (/.false.,.true.,  .true.,  .true./) ! check improvement from prev configuration
  
integer(ip) :: num_swaps(       numRepairs  )
real(rp)    :: tols_swap_ini(   numRepairs  )
real(rp)    :: tols_swap_step(  numRepairs  )

real(rp)    :: tols_smooth(     numRepairs  ) != (/ 0.1, 0.2, 0.3, 0.4/) ! XXYYZZ aqui canviat per abelcruz

!
! Marks
!
integer(ip), parameter :: mark_isNodeDeleted        = -1000_ip
integer(ip), parameter :: mark_isNodeNotFromIniMesh =     0_ip
!
! Global module variables
!
integer(ip) :: numNodes_current = 0
integer(ip) :: numNodes_alloca  = 0
integer(ip), pointer :: isModifiedNode(:)
integer(ip), parameter :: mark_isMod = 1_ip
integer(ip), parameter :: mark_noMod = 0_ip
type(edgeData_type)     :: attempt_edgeSplit
type(edgeData_type)     :: attempt_edgeColl
type(edgeData_type)     :: attempt_edgeSwap
integer(ip), pointer :: attemp_smoothNode(:)
!
!
!
private
!
public :: adapt_mesh_to_metric
!
!
!
CONTAINS
!
!
!
subroutine setAdaptParameters(thresholdQuality)
  implicit none
  
  real(rp), optional, intent(in) ::thresholdQuality
  
  integer(ip), parameter :: n_1 = 0_ip
  integer(ip), parameter :: n_2 = 0_ip
  
  if(.not.present(thresholdQuality)) then
    !call setAdaptParameters(0.4_rp) ! -> removed call to itself to avoid declaring it as a recursion...
    threshold_quality_to_repair = 0.4_rp ! abelcruz !2.0_rp!0.5_rp!0.3_rp

    num_swaps     = (/    n_1,    n_2,   1_ip,   1_ip/)
    tols_swap_ini = (/ 0.1_rp, 0.1_rp, 0.2_rp, 0.3_rp/)
    tols_swap_step= (/ 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp/)
    tols_smooth   = (/ 0.1_rp, 0.1_rp, 0.2_rp, 0.3_rp/)
    
  else
    threshold_quality_to_repair = thresholdQuality+0.001_rp
    
    if(     thresholdQuality<0.1001_rp) then
!       num_swaps     = (/   1_ip,   1_ip,   1_ip,   1_ip/)
!       tols_swap_ini = (/ 0.1_rp, 0.1_rp, 0.1_rp, 0.1_rp/)
!       tols_swap_step= (/ 0.1_rp, 0.1_rp, 0.1_rp, 0.1_rp/)
!       tols_smooth(:)   = 0.1_rp
      num_swaps     = (/     n_1,     n_2,   1_ip,   1_ip/)
      tols_swap_ini = (/ 0.05_rp, 0.05_rp, 0.1_rp, 0.1_rp/)
      tols_swap_step= (/  0.0_rp,  0.0_rp, 0.0_rp, 0.0_rp/)
      tols_smooth(:)   = 0.1_rp
    !
    else if(thresholdQuality<0.2001_rp) then
!       num_swaps     = (/   1_ip,   1_ip,   2_ip,   2_ip/)
!       tols_swap_ini = (/ 0.1_rp, 0.1_rp, 0.1_rp, 0.1_rp/)
!       tols_swap_step= (/ 0.1_rp, 0.1_rp, 0.1_rp, 0.1_rp/)
!       tols_smooth   = (/ 0.1_rp, 0.1_rp, 0.2_rp, 0.2_rp/)
      num_swaps     = (/    n_1,    n_2,   1_ip,   1_ip/)
      tols_swap_ini = (/ 0.1_rp, 0.1_rp, 0.2_rp, 0.2_rp/)
      tols_swap_step= (/ 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp/)
      tols_smooth   = (/ 0.1_rp, 0.1_rp, 0.1_rp, 0.2_rp/)
    !
    else if(thresholdQuality<0.3001_rp) then
!       num_swaps     = (/   1_ip,   1_ip,   2_ip,   3_ip/)
!       tols_swap_ini = (/ 0.1_rp, 0.1_rp, 0.1_rp, 0.1_rp/)
!       tols_swap_step= (/ 0.1_rp, 0.1_rp, 0.1_rp, 0.1_rp/)
!       tols_smooth   = (/ 0.1_rp, 0.1_rp, 0.2_rp, 0.3_rp/)

      num_swaps     = (/    n_1,    n_2,   1_ip,   1_ip/)
      tols_swap_ini = (/ 0.1_rp, 0.1_rp, 0.2_rp, 0.2_rp/)
      tols_swap_step= (/ 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp/)
      tols_smooth   = (/ 0.1_rp, 0.1_rp, 0.2_rp, 0.2_rp/)
    !
    else if(thresholdQuality<0.4001_rp) then
!       num_swaps     = (/   1_ip,   2_ip,   2_ip,   2_ip/)
!       tols_swap_ini = (/ 0.1_rp, 0.1_rp, 0.1_rp, 0.2_rp/)
!       tols_swap_step= (/ 0.1_rp, 0.1_rp, 0.2_rp, 0.2_rp/)
!       tols_smooth   = (/ 0.1_rp, 0.2_rp, 0.3_rp, 0.4_rp/)
      
      num_swaps     = (/    n_1,    n_2,   1_ip,   1_ip/)
      tols_swap_ini = (/ 0.1_rp, 0.1_rp, 0.2_rp, 0.3_rp/)
      tols_swap_step= (/ 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp/)
      tols_smooth   = (/ 0.1_rp, 0.1_rp, 0.2_rp, 0.3_rp/)
    !
    else if(thresholdQuality<0.5001_rp) then
!       num_swaps     = (/   1_ip,   2_ip,   2_ip,   3_ip/)
!       tols_swap_ini = (/ 0.1_rp, 0.1_rp, 0.1_rp, 0.1_rp/)
!       tols_swap_step= (/ 0.1_rp, 0.1_rp, 0.2_rp, 0.2_rp/)
!       tols_smooth   = (/ 0.1_rp, 0.2_rp, 0.3_rp, 0.5_rp/)
      num_swaps     = (/    n_1,    n_2,   1_ip,   1_ip/)
      tols_swap_ini = (/ 0.1_rp, 0.2_rp, 0.2_rp, 0.4_rp/)
      tols_swap_step= (/ 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp/)
      tols_smooth   = (/ 0.1_rp, 0.2_rp, 0.3_rp, 0.4_rp/)
    !
    else if(thresholdQuality<0.6001_rp) then
!       num_swaps     = (/   1_ip,   2_ip,   2_ip,   3_ip/)
!       tols_swap_ini = (/ 0.1_rp, 0.1_rp, 0.1_rp, 0.1_rp/)
!       tols_swap_step= (/ 0.1_rp, 0.1_rp, 0.3_rp, 0.25_rp/)
!       tols_smooth   = (/ 0.1_rp, 0.2_rp, 0.3_rp, 0.6_rp/)
      num_swaps     = (/    n_1,    n_2,   1_ip,   1_ip/)
      tols_swap_ini = (/ 0.1_rp, 0.2_rp, 0.3_rp, 0.5_rp/)
      tols_swap_step= (/ 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp/)
      tols_smooth   = (/ 0.1_rp, 0.2_rp, 0.3_rp, 0.5_rp/)
    !
    else if(thresholdQuality<0.7001_rp) then
!       num_swaps     = (/   1_ip,   2_ip,   2_ip,   3_ip/)
!       tols_swap_ini = (/ 0.1_rp, 0.1_rp, 0.1_rp, 0.1_rp/)
!       tols_swap_step= (/ 0.1_rp, 0.1_rp, 0.3_rp, 0.3_rp/)
!       tols_smooth   = (/ 0.1_rp, 0.2_rp, 0.3_rp, 0.7_rp/)
      num_swaps     = (/    n_1,    n_2,   1_ip,   1_ip/)
      tols_swap_ini = (/ 0.1_rp, 0.2_rp, 0.3_rp, 0.6_rp/)
      tols_swap_step= (/ 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp/)
      tols_smooth   = (/ 0.1_rp, 0.2_rp, 0.3_rp, 0.6_rp/)
    !
    else if(thresholdQuality<0.8001_rp) then
!       num_swaps     = (/   1_ip,   2_ip,   2_ip,   3_ip/)
!       tols_swap_ini = (/ 0.1_rp, 0.1_rp, 0.1_rp, 0.1_rp/)
!       tols_swap_step= (/ 0.1_rp, 0.1_rp, 0.3_rp, 0.35_rp/)
!       tols_smooth   = (/ 0.1_rp, 0.2_rp, 0.4_rp, 0.8_rp/)
      num_swaps     = (/    n_1,    n_2,   1_ip,   1_ip/)
      tols_swap_ini = (/ 0.1_rp, 0.2_rp, 0.4_rp, 0.6_rp/)
      tols_swap_step= (/ 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp/)
      tols_smooth   = (/ 0.1_rp, 0.3_rp, 0.4_rp, 0.6_rp/)
    !
    else if(thresholdQuality<0.9001_rp) then
      num_swaps     = (/   1_ip,   2_ip,   2_ip,   3_ip/)
      tols_swap_ini = (/ 0.1_rp, 0.1_rp, 0.1_rp, 0.1_rp/)
      tols_swap_step= (/ 0.1_rp, 0.1_rp, 0.3_rp, 0.4_rp/)
      tols_smooth   = (/ 0.1_rp, 0.2_rp, 0.4_rp, 0.8_rp/)
    !
    else if(thresholdQuality<1.0001_rp) then
      num_swaps     = (/   1_ip,   2_ip,   2_ip,   3_ip/)
      tols_swap_ini = (/ 0.1_rp, 0.1_rp, 0.1_rp, 0.1_rp/)
      tols_swap_step= (/ 0.1_rp, 0.1_rp, 0.3_rp, 0.4_rp/)
      tols_smooth   = (/ 0.1_rp, 0.2_rp, 0.4_rp, 0.9_rp/)
    else
      call runend('Quality should be in range [0,1]')
    end if
    
  end if
  
end subroutine setAdaptParameters
!
!
!
subroutine adapt_mesh_to_metric(&
  mesh,&
  metric,&
  mapNodes_input_to_adapted_mesh,&
  lock_valid_elems,&
  isModifiedMesh,&
  thresholdQuality,&
  activationQuality)
  !
  use mod_quality,      only: compute_mesh_quality_sizeShape, quality_deallo
  use mod_quality,      only: print_q_stats, do_print_q_stats
  use mod_meshTopology, only: node_to_elems_type
!  use def_master,       only: kfl_paral
  use mod_strings,      only: real_to_string
  use mod_communications_global,    only: PAR_MIN
  use mod_messages,     only : messages_live
  !
  type(mesh_type_basic),          intent(inout) :: mesh
  type(mesh_metric_type),         intent(inout) :: metric
  integer(ip), pointer, optional, intent(inout) :: mapNodes_input_to_adapted_mesh(:)
  logical(lg)         , optional, intent(in)    :: lock_valid_elems
  logical(lg)         , optional, intent(inout) :: isModifiedMesh
  real(rp)            , optional, intent(inout) :: thresholdQuality
  real(rp)            , optional, intent(inout) :: activationQuality
  !
  real(rp), pointer    :: q(:)
  integer(ip)          :: iele_good!, ielem, iele_bad
  integer(ip), pointer :: lnods_good(:,:)
  logical(lg), pointer :: is_bad_ele(:)
  integer(ip), pointer :: mapNodes_input_to_adapted_mesh_private(:)
  !
  real(rp)    :: minQ
  logical(lg) ::isq_allo 
  
  real(rp)    :: t0,t1
  
  integer(ip) :: ii
  !
  isq_allo = .false.
  time_total = 0.0_rp
  if(out_performance) then
    call cpu_time(t0)
    call compute_mesh_quality_sizeShape(mesh,metric,q)
    isq_allo = .true.
    call do_print_q_stats(q)
  end if
  
  call setAdaptParameters(thresholdQuality)
  
  call out_firstLast(mesh,metric,file_out)

  if(do_warnValidity) then ! check mesh validity.. just in case
    if(.not.isq_allo) call compute_mesh_quality_sizeShape(mesh,metric,q)
    call print_q_stats(q)
    if(minval(q)<minQ_validity) then
      call messages_live('Invalid mesh BEFORE mesh adaption','WARNING')
    end if
    isq_allo = .true.
  end if

  if(present(activationQuality)) then
    if(.not.isq_allo) call compute_mesh_quality_sizeShape(mesh,metric,q)
    isq_allo = .true.
  
    minQ = minval(q) 
    !call PAR_MIN(minQ)
    if(minQ>activationQuality) then
      !call messages_live('NOT ACTIVATED -> MIN QUAL: '//real_to_string(minQ)//' > '//real_to_string(activationQuality) )
      if(present(isModifiedMesh)) isModifiedMesh = .false. 
      return
    else
      ! There are bad elements, so we mark it as potentially modified
      ! Another thing is that we cannot modify it.. but may be due to parallel constrains,
      ! so we risk entering again with a different partition... safe decicions right now
      !call messages_live('ACTIVATED -> MIN QUAL: '//real_to_string(minQ)//' < '//real_to_string(activationQuality) )
      if(present(isModifiedMesh)) isModifiedMesh = .true. 
    end if
  end if

  if(present(lock_valid_elems).and.(.not.lock_valid_elems)) then
    if(present(isModifiedMesh)) then
      ! First I thought to update this in adapt_mesh_to_metric_private is mesh is actually modified
      ! Then I decided that not... becaus in parallel iterations mesh could not be modified due
      ! to inner region boundaries but being bad
      ! Another thing is that we cannot modify it.. but may be due to parallel constrains,
      ! so we risk entering again with a different partition... safe decicions right now
      isModifiedMesh = .true. 
    end if
    !
    call adapt_mesh_to_metric_private(mesh,metric,mapNodes_input_to_adapted_mesh_private)
    !
  else
    
    if(.not.isq_allo) call compute_mesh_quality_sizeShape(mesh,metric,q)
    
    nullify(is_bad_ele)
    call memory_alloca(memor_adapt,'is_bad_ele','mod_adapt',is_bad_ele,mesh%nelem)
    is_bad_ele = q(:)<threshold_quality_to_repair
    
    if(count(is_bad_ele).eq.0_ip) then
      if(present(mapNodes_input_to_adapted_mesh).and.mesh%npoin>0_ip) then
        nullify(mapNodes_input_to_adapted_mesh)
        call memory_alloca(memor_adapt,'mapNodes_input_to_adapted_mesh','mod_adapt',mapNodes_input_to_adapted_mesh,mesh%npoin)
        !mapNodes_input_to_adapted_mesh = (/(ielem, ielem=1_ip,mesh%npoin, 1_ip)/)
        do ii=1_ip,mesh%npoin
          mapNodes_input_to_adapted_mesh(ii) = ii
        end do
      end if
      if(present(isModifiedMesh)) isModifiedMesh = .false.
      return
    else
      if(present(isModifiedMesh)) then
        isModifiedMesh = .true.
      end if
    end if
    !
    !
    !
    call  set_submesh(mesh,numNieghLevels,is_bad_ele,iele_good,lnods_good)
    !
    call memory_deallo(memor_adapt,'is_bad_ele','mod_adapt',is_bad_ele)
    !
    !
    !
    call adapt_mesh_to_metric_private(mesh,metric,mapNodes_input_to_adapted_mesh_private)    
    !
    !
    !
    call rebuild_mesh_fromSubmesh(mesh,iele_good,lnods_good,mapNodes_input_to_adapted_mesh_private)
    !
    !
    !
  end if
  
  if(present(mapNodes_input_to_adapted_mesh)) then ! if present, allocate it
    nullify(mapNodes_input_to_adapted_mesh)
    call memory_alloca(memor_adapt,'mapNodes_input_to_adapted','mod_adapt',mapNodes_input_to_adapted_mesh,size(mapNodes_input_to_adapted_mesh_private))
    mapNodes_input_to_adapted_mesh = mapNodes_input_to_adapted_mesh_private
    call memory_deallo(memor_adapt,'mapNodes_input_to_adapted_mesh','mod_adapt',mapNodes_input_to_adapted_mesh_private)
  end if
  
  call out_firstLast(mesh,metric,file_out)
  
  if(do_warnValidity) then
    if(isq_allo) call quality_deallo(q)
    call compute_mesh_quality_sizeShape(mesh,metric,q)
    call print_q_stats(q)
    if(minval(q)<minQ_validity) then
      call messages_live('Invalid mesh AFTER mesh adaption','WARNING')
    end if
  end if
  
  if(isq_allo) call quality_deallo(q)
  
  if(out_performance) then
    call cpu_time(t1)
    time_total = t1-t0

    call compute_mesh_quality_sizeShape(mesh,metric,q)
    call do_print_q_stats(q)
    
    call out_performance_profile()
    call out_performance_profile_global()
  end if
  
end subroutine adapt_mesh_to_metric
!
!
!
!-----------------------------------------------------------------------
!
!> @author  Abel Gargallo
!> @date    28/04/2021
!> @brief   set_submesh
!> @details set_submesh 
!>          set_submesh from mesh, extracting the subelems of mesh
!>          lnods_good is the complementary elements of sub_mesh
!>          once it works replace the variable names
!
!-----------------------------------------------------------------------
subroutine set_submesh(mesh,numNieghLevels,is_bad_ele,iele_good,lnods_good,q,q_fixed) !lnods_good is the complementary of the submesh in mesh
  use mod_meshTopology, only: node_to_elems_type
  use mod_quality, only: quality_deallo
  implicit none
  type(mesh_type_basic),      intent(inout) :: mesh
  integer(ip),                intent(in)    :: numNieghLevels
  logical(lg),                intent(inout) :: is_bad_ele(:)!mesh%nelem)
  integer(ip),                intent(out)   :: iele_good
  integer(ip), pointer,       intent(inout) :: lnods_good(:,:)
  real(rp), pointer, optional,intent(inout) :: q(:)
  real(rp), pointer, optional,intent(inout) :: q_fixed(:)
  !
  type(node_to_elems_type)  :: node_to_elems
  integer(ip)               :: ielem, iele_bad, inode
  integer(ip), pointer      :: list_good(:)
  integer(ip), pointer      :: list_bad(:)
  integer(ip), pointer      :: lnods_bad(:,:)
!   integer(ip), pointer      :: lnods_aux(:,:)
  logical(lg), pointer      :: is_bad_ele_prev(:)
  integer(ip)               :: ltype_save
  integer(ip)               :: ilevel
  integer(ip), pointer      :: elems_adjToNode(:)
  real(rp)   , pointer      :: q_copy(:)
  integer(ip)               :: numElemNodes
  integer(ip)               :: ii,jj
  !
  if(numNieghLevels>0_ip) then
    call node_to_elems%set(mesh)
    nullify(is_bad_ele_prev)
    call memory_alloca(memor_adapt,'is_bad_ele_prev','mod_adapt',is_bad_ele_prev,mesh%nelem)
    do ilevel = 1,numNieghLevels
      is_bad_ele_prev = is_bad_ele
      do ielem=1,mesh%nelem
        if(is_bad_ele_prev(ielem)) then
          do inode=1,size(mesh%lnods,1)
            call node_to_elems%get_elems(mesh%lnods(inode,ielem),elems_adjToNode)
            is_bad_ele(elems_adjToNode) = .true.
            call node_to_elems%delete_elems(elems_adjToNode)
          end do
        end if
      end do
    end do
    call memory_deallo(memor_adapt,'is_bad_ele_prev','mod_adapt',is_bad_ele_prev)
    call node_to_elems%deallo()
  end if
  
  nullify(list_good)
  call memory_alloca(memor_adapt,'list_good','mod_adapt',list_good,mesh%nelem)
  nullify(list_bad)
  call memory_alloca(memor_adapt,'list_bad' ,'mod_adapt',list_bad ,mesh%nelem)
  iele_good = 0_ip
  iele_bad  = 0_ip
  do ielem=1,mesh%nelem
    if(is_bad_ele(ielem)) then
      iele_bad = iele_bad+1_ip
      list_bad(iele_bad) = ielem
    else
      iele_good = iele_good+1_ip
      list_good(iele_good) = ielem
    end if
  end do
  
  if(iele_good>0_ip) then
    numElemNodes = INT(size(mesh%lnods,1),ip)
    
    nullify(lnods_good)
    call memory_alloca(memor_adapt,'lnods_good','mod_adapt',lnods_good,numElemNodes,iele_good)
    do jj=1,iele_good
      lnods_good(:,jj)=mesh%lnods(:,list_good(jj)) !lnods_good = mesh%lnods(:,list_good(1:iele_good))
    end do

    nullify(lnods_bad)
    call memory_alloca(memor_adapt,'lnods_bad' ,'mod_adapt',lnods_bad ,numElemNodes,iele_bad)
    !
    do jj=1,iele_bad
        lnods_bad(:,jj)=mesh%lnods(:,list_bad(jj)) !lnods_bad = mesh%lnods(:,list_bad(1:iele_bad))
    end do

    mesh%nelem = iele_bad
    call memory_deallo(memor_adapt,'mesh%lnods' ,'mod_adapt',mesh%lnods )
    nullify(mesh%lnods)
    call memory_alloca(memor_adapt,'lnods','mod_adapt',mesh%lnods,numElemNodes,iele_bad)
    do jj=1,iele_bad
      mesh%lnods(:,jj) = lnods_bad(:,jj) !mesh%lnods = lnods_bad!(:,1:iele_bad)
    end do
    call memory_deallo(memor_adapt,'lnods_bad' ,'mod_adapt',lnods_bad )

    ltype_save = mesh%ltype(1_ip)
    call memory_deallo(memor_adapt,'ltype'      ,'mod_adapt',mesh%ltype )
    call memory_deallo(memor_adapt,'leinv_loc'  ,'mod_adapt',mesh%leinv_loc )
    call memory_deallo(memor_adapt,'perme'      ,'mod_adapt',mesh%perme )
    nullify(mesh%ltype)
    call memory_alloca(memor_adapt,'ltype'      ,'mod_adapt',mesh%ltype,mesh%nelem)
    nullify(mesh%leinv_loc)
    call memory_alloca(memor_adapt,'leinv_loc'  ,'mod_adapt',mesh%leinv_loc,mesh % nelem)
    nullify(mesh%perme)
    call memory_alloca(memor_adapt,'perme'      ,'mod_adapt',mesh%perme,mesh % nelem)
    do ielem = 1,mesh % nelem
       mesh %     ltype(ielem) = ltype_save
       mesh % leinv_loc(ielem) = ielem
       mesh %     perme(ielem) = ielem
    end do
    
    if(present(q)) then
      nullify(q_fixed)
      call memory_alloca(memor_adapt,'q_fixed','mod_adapt',q_fixed,iele_good)
      do ii=1,iele_good
        q_fixed(ii) = q(list_good(ii)) !q_fixed = q(list_good(1:iele_good))
      end do
    
      nullify(q_copy)
      call memory_alloca(memor_adapt,'q_copy','mod_adapt',q_copy,iele_bad)
      do ii=1,iele_bad
        q_copy(ii) = q(list_bad(ii)) !q_copy = q(list_bad(1:iele_bad))
      end do

      call quality_deallo(q)
      nullify(q)
      call memory_alloca(memor_adapt,'q','compute_mesh_quality',q,iele_bad)
      do ii=1,iele_bad
        q(ii) = q_copy(ii) !q = q_copy
      end do
      
      call memory_deallo(memor_adapt,'q_copy','mod_adapt',q_copy)
    end if
  end if
  !
  call memory_deallo(memor_adapt,'list_good','mod_adapt',list_good)
  call memory_deallo(memor_adapt,'list_bad' ,'mod_adapt',list_bad )
  !
end subroutine set_submesh
!
!
!
subroutine rebuild_mesh_fromSubmesh(mesh,iele_good,lnods_good,mapNodes_ini_to_final,q,q_fixed,isNewElem)
  use mod_quality, only: quality_deallo
  implicit none
  type(mesh_type_basic),      intent(inout) :: mesh
  integer(ip),                intent(in)    :: iele_good
  integer(ip),pointer,        intent(inout) :: lnods_good(:,:)
  integer(ip),optional,       intent(in)    :: mapNodes_ini_to_final(:)
  real(rp), pointer, optional,intent(inout) :: q(:)
  real(rp), pointer, optional,intent(inout) :: q_fixed(:)
  logical(lg),pointer,optional,intent(inout):: isNewElem(:)
  !
  integer(ip)               :: ielem, iele_aux, inode, ltype_save
  integer(ip), pointer      :: lnods_aux(:,:)
  real(rp),    pointer      :: q_copy(:)
  logical(lg), pointer      :: isNewElem_copy(:)
  integer(ip)               :: ii
  !
  if(iele_good>0_ip) then
    iele_aux = mesh%nelem
    nullify(lnods_aux)
    call memory_alloca(memor_adapt,'lnods_aux' ,'mod_adapt',lnods_aux ,INT(size(mesh%lnods,1),ip),iele_aux)
    lnods_aux = mesh%lnods
    do ii=1,iele_aux
      lnods_aux(:,ii)=mesh%lnods(:,ii)
    end do

    mesh%nelem = iele_good + iele_aux
    call memory_deallo(memor_adapt,'mesh%lnods' ,'mod_adapt',mesh%lnods )
    nullify(mesh%lnods)
    call memory_alloca(memor_adapt,'mesh%lnods','mod_adapt',mesh%lnods,INT(size(lnods_good,1),ip),mesh%nelem)
    if(present(mapNodes_ini_to_final)) then
      do inode=1,size(mesh%lnods,1)
        mesh%lnods(inode,1:iele_good) = mapNodes_ini_to_final(lnods_good(inode,:))
      end do
    else
      do ii=1,iele_good
        mesh%lnods(:,ii)=lnods_good(:,ii) !mesh%lnods(:,1:iele_good) = lnods_good(:,:)
      end do
    end if
    mesh%lnods(:,(iele_good+1):mesh%nelem) = lnods_aux
    call memory_deallo(memor_adapt,'lnods_good' ,'mod_adapt',lnods_good)
    call memory_deallo(memor_adapt,'lnods_aux' ,'mod_adapt',lnods_aux)

    ltype_save = mesh%ltype(1_ip)
    call memory_deallo(memor_adapt,'ltype'      ,'mod_adapt',mesh%ltype )
    call memory_deallo(memor_adapt,'leinv_loc'  ,'mod_adapt',mesh%leinv_loc )
    call memory_deallo(memor_adapt,'perme'      ,'mod_adapt',mesh%perme )
    nullify(mesh%ltype)
    call memory_alloca(memor_adapt,'ltype'      ,'mod_adapt',mesh%ltype,mesh%nelem)
    nullify(mesh%leinv_loc)
    call memory_alloca(memor_adapt,'leinv_loc'  ,'mod_adapt',mesh%leinv_loc,mesh % nelem)
    nullify(mesh%perme)
    call memory_alloca(memor_adapt,'perme'      ,'mod_adapt',mesh%perme,mesh % nelem)
    do ielem = 1,mesh % nelem
       mesh %     ltype(ielem) = ltype_save
       mesh % leinv_loc(ielem) = ielem
       mesh %     perme(ielem) = ielem
    end do
    if(present(q)) then
      nullify(q_copy)
      call memory_alloca(memor_adapt,'q_copy' ,'mod_adapt',q_copy,size(q))
      do ii=1,size(q)
        q_copy(ii) = q(ii) !q_copy = q
      end do
      call quality_deallo(q)
      nullify(q)
      call memory_alloca(memor_adapt,'q' ,'compute_mesh_quality',q,mesh%nelem)
      do ii=1,iele_good
        q(ii) = q_fixed(ii) !q(1:iele_good) = q_fixed
      end do
      call memory_deallo(memor_adapt,'q_fixed','mod_adapt',q_fixed)
      do ii=1,size(q_copy)
        q(iele_good+ii) = q_copy(ii) !q((iele_good+1):mesh%nelem) = q_copy
      end do
      call memory_deallo(memor_adapt,'q_copy' ,'mod_adapt',q_copy)
    end if
    if(present(isNewElem)) then
      nullify(isNewElem_copy)
      call memory_alloca(memor_adapt,'isNewElem_copy' ,'mod_adapt',isNewElem_copy,INT(size(isNewElem),ip))
      !isNewElem_copy = isNewElem
      do ii=1,size(isNewElem_copy)
        isNewElem_copy(ii) = isNewElem(ii)
      end do
      call memory_deallo(memor_adapt,'isNewElem' ,'mod_adapt',isNewElem)
      nullify(isNewElem)
      call memory_alloca(memor_adapt,'isNewElem' ,'mod_adapt',isNewElem,mesh%nelem)
      do ii=1,iele_good
        isNewElem(ii) = .false. !isNewElem(1:iele_good) =.false.
      end do
      do ii=1,size(isNewElem_copy)
        isNewElem(iele_good+ii) = isNewElem_copy(ii) !isNewElem((iele_good+1):mesh%nelem) = isNewElem_copy
      end do
      call memory_deallo(memor_adapt,'isNewElem_copy' ,'mod_adapt',isNewElem_copy)
    end if
  end if
  !
end subroutine rebuild_mesh_fromSubmesh
!
!
!
subroutine adapt_mesh_to_metric_private(mesh,metric,mapNodes_input_to_adapted_mesh)
  use mod_smoothing,    only: smooth_mesh
  use mod_meshTopology, only: setBoundaryEdgeData
  use mod_meshTopology, only: getIsBoundaryNode
  use mod_quality,      only: compute_mesh_quality_sizeShape
  !use def_master, only: kfl_paral
  use mod_debugTools,   only: deb_collTimeTotal,  deb_splitTimeTotal, deb_swapTimeTotal, deb_optiTimeTotal, deb_timeTotal
  !
  implicit none
  !
  type(mesh_type_basic),            intent(inout) :: mesh
  type(mesh_metric_type),           intent(inout) :: metric
  integer(ip), pointer, optional,   intent(inout) :: mapNodes_input_to_adapted_mesh(:)
  
  integer(ip) :: i_repair
  
  real(rp)    :: tol_length_collapse, tol_length_split, tol_qual_swap, tol_qual_smooth
  logical(lg) :: do_check_qual_coll , do_check_qual_split
  integer(ip) :: iswap
  integer(ip) :: iout
  real(rp) :: t0,t1

  integer(ip), pointer  :: mapNodes_new_to_old(:)
  integer(ip) :: iaux, node_old, node_new
  
  type(edgeData_type) :: edgeData_boun
  logical(lg), pointer:: isBouNode_iniMesh(:)
!  real(rp), pointer :: q(:)
  
  logical(lg), pointer :: isRemovedNode(:)
  real(rp) :: t0priv,t1priv
  !
  if(out_performance) call cpu_time(t0priv)
  
  time_split  = 0.0_rp
  time_coll   = 0.0_rp
  time_swap   = 0.0_rp
  time_smooth = 0.0_rp
  !
  nullify(mapNodes_new_to_old)
  !call memory_alloca(memor_adapt,'mapNodes_new_to_old','mod_adapt',mapNodes_new_to_old,size(mesh%coord,2))
  call memory_alloca(memor_adapt,'mapNodes_new_to_old','mod_adapt',mapNodes_new_to_old,mesh%npoin)
  do iaux = 1_ip,mesh%npoin
    mapNodes_new_to_old(iaux)=iaux ! identity mapping to start
  end do
  
  if(present(mapNodes_input_to_adapted_mesh)) then ! if present, allocate it
    nullify(mapNodes_input_to_adapted_mesh)
    call memory_alloca(memor_adapt,'mapNodes_input_to_adapted_mesh','mod_adapt',mapNodes_input_to_adapted_mesh,mesh%npoin)! size(mesh%coord,2)
  end if
  
  edgeData_boun = setBoundaryEdgeData(mesh%lnods)
  call getIsBoundaryNode(mesh,isBouNode_iniMesh)
  
  if(out_debug_text) then
    print*,'Num points and elems before repair'
    print*,size(mesh%coord,2)
    print*,size(mesh%lnods,2)
  end if

  iout = -1_ip
  call out_with_q(mesh,metric,file_out,iout)
  
  do i_repair = 1,numRepairs
    
    if(out_debug_text.and.out_verbose) then
      print*,"i_repair: ",i_repair
      print*,'Num points:',size(mesh%coord,2)
      print*,'Num elems: ',size(mesh%lnods,2)
    end if

    if(out_debug_text.and.out_verbose) print*,"_____COLLAPSE___________________________________________________________________________"
    tol_length_collapse = tols_coll(  i_repair)
    do_check_qual_coll  = doCheckColl(i_repair)
    if(out_debug_text.and.out_verbose) print '(" -> collapse: ",f5.3,1x , L1)',tol_length_collapse,do_check_qual_coll
    if((out_debug_text.and.out_verbose.and.out_time).or.out_performance) call cpu_time(t0)
    !
    call repair_mesh_edgeLengths(mesh,metric,tol_length_collapse,do_check_qual_coll,&
      mapNodes_new_to_old, edgeData_boun, isBouNode_iniMesh) ! collapse
    !
    if((out_debug_text.and.out_verbose.and.out_time).or.out_performance) call cpu_time(t1)
    if(out_debug_text.and.out_verbose.and.out_time) print '("    Time = ",f6.0," sec")',(t1-t0)
    if(out_performance) time_coll = time_coll+(t1-t0)
    if(out_performance) deb_collTimeTotal = deb_collTimeTotal + (t1-t0)
    call out_with_q(mesh,metric,file_out,iout)
    !
    if(out_debug_text.and.out_verbose) print*,"_____SPLIT______________________________________________________________________________"
    tol_length_split    = tols_split(  i_repair)
    do_check_qual_split = doCheckSplit(i_repair)
    if(out_debug_text.and.out_verbose) print '(" -> split: ",f5.3,1x , L1)',tol_length_split,do_check_qual_split
    if((out_debug_text.and.out_verbose.and.out_time).or.out_performance) call cpu_time(t0)
    !
    call repair_mesh_edgeLengths(mesh,metric,tol_length_split   ,do_check_qual_split,&
      mapNodes_new_to_old, edgeData_boun, isBouNode_iniMesh) !split
    !
    if((out_debug_text.and.out_verbose.and.out_time).or.out_performance) call cpu_time(t1)
    if(out_debug_text.and.out_verbose.and.out_time) print '("    Time = ",f6.0," sec")',(t1-t0)
    if(out_performance) time_split = time_split+(t1-t0)
    if(out_performance) deb_splitTimeTotal = deb_splitTimeTotal + (t1-t0)
    call out_with_q(mesh,metric,file_out,iout)
    !
    if(out_debug_text.and.out_verbose) print*,"_____SWAP______________________________________________________________________________"
    do iswap=1,num_swaps(i_repair)
      tol_qual_swap       = tols_swap_ini(i_repair) + real(iswap-1,rp)*tols_swap_step(i_repair)
      if(out_debug_text.and.out_verbose) print '(" -> swap: target q ",f3.1)',tol_qual_swap
      if((out_debug_text.and.out_verbose.and.out_time).or.out_performance) call cpu_time(t0)
      !
      call swap_mesh_edges(mesh,metric,tol_qual_swap, mapNodes_new_to_old, edgeData_boun, isBouNode_iniMesh)
      !
      if((out_debug_text.and.out_verbose.and.out_time).or.out_performance) call cpu_time(t1)
      if(out_debug_text.and.out_verbose.and.out_time) print '("    Time = ",f6.0," sec")',(t1-t0)
      if(out_performance) time_swap = time_swap+(t1-t0)
      if(out_performance) deb_swapTimeTotal = deb_swapTimeTotal + (t1-t0)
      call out_with_q(mesh,metric,file_out,iout)
    end do
    !
    if(out_debug_text.and.out_verbose) print*,"______SMOOTH____________________________________________________________________________"
    tol_qual_smooth     = tols_smooth( i_repair)
    if(out_debug_text.and.out_verbose) print '(" -> smoothing: target q ",f3.1)',tol_qual_smooth
    if((out_debug_text.and.out_verbose.and.out_time).or.out_performance) call cpu_time(t0)
    !
    call smooth_mesh(mesh,metric,tol_qual_smooth)
    !
    if((out_debug_text.and.out_verbose.and.out_time).or.out_performance) call cpu_time(t1)
    if(out_debug_text.and.out_verbose.and.out_time) print '("    Time = ",f6.0," sec")',(t1-t0)
    if(out_performance) time_smooth = time_smooth+(t1-t0)
    if(out_performance) deb_optiTimeTotal = deb_optiTimeTotal + (t1-t0)
    call out_with_q(mesh,metric,file_out,iout)
    !
  end do
  
  ! Remove unused nodes
  nullify(isRemovedNode)
  call memory_alloca(memor_adapt,'isRemovedNode','mod_adapt',isRemovedNode,mesh%npoin)
  !isRemovedNode(:) = mapNodes_new_to_old(:)==mark_isNodeDeleted
  do iaux = 1_ip,mesh%npoin
    isRemovedNode(iaux) = mapNodes_new_to_old(iaux)==mark_isNodeDeleted
  end do
  call removeDeletedNodes(mesh,metric,isRemovedNode,mapNodes_new_to_old)
  call memory_deallo(memor_adapt,'isRemovedNode','mod_adapt',isRemovedNode)
  !
  ! Finalazing adaptation:
  !
  if(out_debug_text.and.out_verbose) print*,'Mesh adaptation finalized'
  if(out_debug_text.and.out_verbose) print*,'Num points:',size(mesh%coord,2)
  if(out_debug_text.and.out_verbose) print*,'Num elems: ',size(mesh%lnods,2)
  
  if(present(mapNodes_input_to_adapted_mesh)) then
    mapNodes_input_to_adapted_mesh(:) = 0_ip
    do node_new=1,size(mapNodes_new_to_old)
      node_old = mapNodes_new_to_old(node_new)
      if( node_old>0_ip) then
        mapNodes_input_to_adapted_mesh(node_old) = node_new
      end if
    end do
  end if
  call memory_deallo(memor_adapt,'mapNodes_new_to_old','mod_adapt',mapNodes_new_to_old)
  call memory_deallo(memor_adapt,'isBouNode_iniMesh','mod_adapt',isBouNode_iniMesh)
  
  call edgeData_boun    %deallo()
  call attempt_edgeColl %deallo()
  call attempt_edgeSplit%deallo()
  call attempt_edgeSwap %deallo()
  
  if(out_performance) then
    call cpu_time(t1priv)
    deb_timeTotal = deb_timeTotal +(t1priv-t0priv)
  end if
  
  return 
end subroutine adapt_mesh_to_metric_private
!
!
!
!--------- SUBMODULE 1: _adapt_lengthRepair ---------------------------------------------------------------------------
! include "_adapt_lengthRepair.f90"
!
!
!
subroutine repair_mesh_edgeLengths(mesh,metric,tol_length,do_checkQual,mapNodes_new_to_old,edgeData_boun,isBouNode_iniMesh)
  use mod_meshTopology,       only: setIsBoundaryNodeFromIniMesh, deleteIsBouNode

!  use def_master,    only: kfl_paral
  use mod_out_paraview, only: out_paraview_inp
  
  implicit none
  
  type(mesh_type_basic),   intent(inout) :: mesh
  type(mesh_metric_type),  intent(inout) :: metric
  real(rp),                intent(in)    :: tol_length
  logical(lg),             intent(in)    :: do_checkQual
  integer(ip), pointer,    intent(inout) :: mapNodes_new_to_old(:)
  type(edgeData_type),     intent(in)    :: edgeData_boun
  logical(lg), pointer,    intent(in)    :: isBouNode_iniMesh(:)
  !
  integer(ip), parameter :: max_iter_frozen = 1000
  integer(ip), parameter :: max_iter_global = 100
  integer(ip), parameter :: numNieghLevels = 1
  !
  logical(lg) :: isModifiedMesh_frozen, hasBeenModified
  integer(ip) :: iter
  integer(ip) :: counter_global
  logical(lg) :: isModified_global

  integer(ip) :: current_strategy
  integer(ip) :: numTaggedEdges
  integer(ip), pointer :: taggedEdges_to_node(:,:)
  logical(lg), pointer :: isPerformedEdge(:)
  logical(lg), pointer :: isRemovedNode(:)

  logical(lg), pointer   :: isBoundaryNode(:)
  integer(ip) :: numPerfEdges
  integer(ip) :: iout

  logical(lg), pointer :: isNewElem(:)
  integer(ip) :: numElems_fixed
  integer(ip), pointer :: lnods_fixed(:,:)
  logical(lg) :: do_useSubmesh
  integer(ip) :: ii
    
  iout= -1_ip
  
  if ( tol_length.lt.(1.0_rp) ) then ! COLLAPSE
    current_strategy = adapt_strat_collapse
    if(.not.attempt_edgeColl%isAlloca())  then
      call attempt_edgeColl%alloca( INT(size(mesh%lnods,1)*size(mesh%lnods,2),ip) )
    end if
  else !SPLIT 
    current_strategy = adapt_strat_split
    if(.not.attempt_edgeSplit%isAlloca()) then
      call attempt_edgeSplit%alloca( INT(size(mesh%lnods,1)*size(mesh%lnods,2),ip) )
    end if
  end if
  
  counter_global    = 0_ip
  isModified_global = .true.
  numPerfEdges      = 0_ip
  do while(isModified_global)
    !
    do_useSubmesh = counter_global>0 ! in the first iteration, check all mesh
    if(do_useSubmesh) then
      call set_submesh(mesh,numNieghLevels,isNewElem(:),numElems_fixed,lnods_fixed)!,q=q,q_fixed=q_fixed)
      call memory_deallo(memor_adapt,'isNewElem','mod_adapt',isNewElem)
    end if
    nullify(isNewElem)
    call memory_alloca(memor_adapt,'isNewElem','mod_adapt',isNewElem,mesh%nelem)
    do ii=1,mesh%nelem
      isNewElem(ii) = .false. !isNewElem(:) = .false.
    end do
    !
    call set_tagged_edges_byLength(current_strategy,tol_length,mesh,metric,numTaggedEdges,taggedEdges_to_node)
    
    if(numTaggedEdges>0) then
      
      call setIsBoundaryNodeFromIniMesh(isBoundaryNode,isBouNode_iniMesh,mapNodes_new_to_old)
    
      nullify(isPerformedEdge)
      call memory_alloca(memor_adapt,'isPerformedEdge','repair_mesh_edgeLengths',isPerformedEdge,int(size(taggedEdges_to_node,2),ip))
      do ii=1,size(isPerformedEdge)
        isPerformedEdge(ii) = .false. !isPerformedEdge = .false.
      end do
      
      !! The next is if I delete deleted nodes each iter, if not, it should be updated (althought it does not affect anywhere)
      nullify(isRemovedNode)
      call memory_alloca(memor_adapt,'isRemovedNode','repair_mesh_edgeLengths',isRemovedNode,mesh%npoin)
      do ii=1,size(isRemovedNode)
        isRemovedNode(ii) = .false. !isRemovedNode = .false.
      end do
    
      hasBeenModified = .false.
      isModifiedMesh_frozen = .true.
      iter = 0
      do while ( isModifiedMesh_frozen.and.(iter<max_iter_frozen) )
        iter = iter+1
        call repair_mesh_edgeLengths_frozen(mesh,metric,&
          current_strategy,do_checkQual,taggedEdges_to_node,&
          isBoundaryNode,isPerformedEdge,isRemovedNode,isModifiedMesh_frozen,numPerfEdges,&
          mapNodes_new_to_old,edgeData_boun,isNewElem)

        hasBeenModified = isModifiedMesh_frozen.or.hasBeenModified
      end do
      
      call updateNodeData(mesh,metric,isRemovedNode,mapNodes_new_to_old) !removeDeletedNodes
      
      call memory_deallo(memor_adapt,'taggedEdges_to_node','mod_adapt',taggedEdges_to_node)
      call memory_deallo(memor_adapt,'isPerformedEdge',    'repair_mesh_edgeLengths',isPerformedEdge    )
      call memory_deallo(memor_adapt,'isRemovedNode',      'repair_mesh_edgeLengths',isRemovedNode      )
      call deleteIsBouNode(isBoundaryNode)
      !
    else
      !
      hasBeenModified = .false.
      !
    end if
    !
    if(do_useSubmesh) then
      call rebuild_mesh_fromSubmesh(mesh,numElems_fixed,lnods_fixed,isNewElem=isNewElem)!,q=q,q_fixed=q_fixed)
    end if
    !
    if(output_steps_edgeRep_paraview)  call out_with_q(mesh,metric,file_out//"_edgeRep",iout)
    !
    counter_global = counter_global+1
    isModified_global = hasBeenModified .and.(counter_global<max_iter_global)
    !
  end do
  !
  return 
end subroutine repair_mesh_edgeLengths
!
!
!
subroutine repair_mesh_edgeLengths_frozen(mesh,metric,&
  current_strategy,do_checkQual,taggedEdges_to_node,&
  isBoundaryNode,isPerformedEdge,isRemovedNode,isModifiedMesh,numPerfEdges,&
  mapNodes_new_to_old,edgeData_boun,&
  isNewElem)
  use mod_cavity, only: cavity_insert_point, cavity_reinsert, delete_Tcav_remeshed
  use mod_quality,            only: compute_minQuality_sizeShape_cavity
  use mod_meshTopology,       only: node_to_elems_type!, delete_arrayElems
  use def_master,    only: kfl_paral
  use mod_debugTools, only: out_performance, deb_splitLoop, deb_collLoop

  implicit none
  
  type(mesh_type_basic),    intent(inout) :: mesh
  type(mesh_metric_type),   intent(inout) :: metric
  integer(ip),              intent(in)    :: current_strategy
  logical(lg),              intent(in)    :: do_checkQual
  integer(ip), pointer,     intent(in)    :: taggedEdges_to_node(:,:)
  logical(lg), pointer,     intent(in)    :: isBoundaryNode(:)
  logical(lg), pointer,     intent(inout) :: isPerformedEdge(:)
  logical(lg), pointer,     intent(inout) :: isRemovedNode(:)
  logical(lg),              intent(out)   :: isModifiedMesh
  integer(ip),              intent(inout) :: numPerfEdges
  integer(ip), pointer,     intent(in)    :: mapNodes_new_to_old(:)
  type(edgeData_type),      intent(in)    :: edgeData_boun
  logical(lg), pointer,     intent(inout) :: isNewElem(:)
  !
  integer(ip) :: numTaggedEdges, numTaggedEdges_remaining
  logical(lg), pointer :: isFrozenElem(:)
  integer(ip) :: iedge
  integer(ip) :: n1,n2
  integer(ip) :: n1_initialMesh,n2_initialMesh
  integer(ip), pointer :: elems_cavity(:)
  logical(lg) :: isValidRemesh, isRemovedEdge,isEdgeAllowed
  logical(lg) :: isBoundEdge, isPotentialBoundEdge
  integer(ip) :: last_point_new, last_point_global, last_elem_new, ini_elem
  real(rp)    :: p_insert(mesh%ndime)
  type(mesh_metric_type):: metric_new
  integer(ip) :: max_newPoints, max_newElems
  real(rp),    pointer :: coord_new(:,:)
  integer(ip), pointer :: lnods_new(:,:)
  integer(ip), pointer :: Tcav_remeshed(:,:)
  integer(ip) :: num_new_elems
  integer(ip) :: nodeId_toInsert
  logical(lg) :: isNodeId_existent
  logical(lg) :: isExceededAllocatedSize
  real(rp)    :: q_prev, q_new
  type(node_to_elems_type) :: node_to_elems
  real(rp) :: M_p_insert(1,1)
  real(rp) :: t0,t1
  !
  if(mesh%ndime==0_ip) call runend('possible memory problem creating p_insert(ndime) of ndime=0')
  
  call node_to_elems%set(mesh)
  
  if(out_performance) call cpu_time(t0)
  
  numTaggedEdges = size(taggedEdges_to_node,2)
  numTaggedEdges_remaining = numTaggedEdges- count(isPerformedEdge)
  max_newPoints = numTaggedEdges_remaining 
  nullify(coord_new)
  call memory_alloca(memor_adapt,'coord_new','repair_mesh_edgeLengths_frozen',coord_new,mesh%ndime,max_newPoints)
  
  select case (current_strategy)
    case  (adapt_strat_collapse)
      if(mesh%ndime==2) then ! each node has an adjacency approx of 6 elems, 2 nodes per edge
        max_newElems = numTaggedEdges_remaining*6_ip*2_ip
      else ! each node has an adjacency approx of 16 elems, 2 nodes per edge
        max_newElems = numTaggedEdges_remaining*16_ip*2_ip
      end if
    case  (adapt_strat_split)
      if(mesh%ndime==2) then ! and edge is shared by 2 triangles (each one split in 2)
        max_newElems = numTaggedEdges_remaining*4_ip
      else ! edge is contained at most by 10 elems (in average), each element is split in 2
        max_newElems = numTaggedEdges_remaining*10_ip*2_ip
      end if
    case default
      call runend('Not implemented strategy in repair_mesh_edges_frozen')
  end select
  nullify(lnods_new)
  call memory_alloca(memor_adapt,'lnods_new','repair_mesh_edgeLengths_frozen',lnods_new,int(size(mesh%lnods,1),ip),max_newElems)  
  call metric_new%alloca(metric%dim,max_newPoints,metric%isSizeField)
  
  last_point_new    = 0_ip
  last_point_global = size(mesh%coord,2)
  last_elem_new     = 0_ip !size(mesh%lnods,2)
  
  nullify(isFrozenElem)
  call memory_alloca(memor_adapt,'isFrozenElem'    ,'repair_mesh_edgeLengths_frozen',isFrozenElem,mesh%nelem)
  isFrozenElem(:) = .false.
  isModifiedMesh  = .false.
  do iedge = 1,numTaggedEdges
    if(isPerformedEdge(iedge)) cycle
    
    n1 = taggedEdges_to_node(1,iedge)
    n2 = taggedEdges_to_node(2,iedge)
    isRemovedEdge        = isRemovedNode(n1).or.isRemovedNode(n2)
    if(isRemovedEdge) cycle
    
    isPotentialBoundEdge = isBoundaryNode(n1).and.isBoundaryNode(n2)
    if(isPotentialBoundEdge) then
      if(current_strategy.eq.adapt_strat_collapse) then !cannot collapse two boundary nodes, even if edge is interior
        cycle
      else
        n1_initialMesh = mapNodes_new_to_old(n1)
        n2_initialMesh = mapNodes_new_to_old(n2)
        isBoundEdge = edgeData_boun%isEdge(n1_initialMesh,n2_initialMesh)
        if(isBoundEdge) cycle
      end if 
    end if
    
    isEdgeAllowed = .not.isAttemptedEdgeLengthRepair(n1,n2,current_strategy)
    if(isEdgeAllowed) then
      select case (current_strategy)
        case  (adapt_strat_collapse)
          call node_to_elems%get_edgeAdjacency(n1,n2,elems_cavity)
        case  (adapt_strat_split)
          call node_to_elems%get_edgeContainers(n1,n2,elems_cavity)
      end select
      
      isEdgeAllowed = .not.(any(isFrozenElem(elems_cavity)))
      if( isEdgeAllowed ) then
        call setAttemptedEdgeLengthRepair(n1,n2,current_strategy)
        
        if((current_strategy.eq.adapt_strat_collapse).and.(isBoundaryNode(n1).or.isBoundaryNode(n2))) then
          isNodeId_existent = .true. ! if one of the two nodes is boundary, the collapse keeps its id
          if(    isBoundaryNode(n1)) then
            nodeId_toInsert = n1
          elseif(isBoundaryNode(n2)) then
            nodeId_toInsert = n2
          else
            call runend('Not possible to collapse and edge with two boundary nodes...')
          end if
        else
          last_point_global = last_point_global + 1_ip
          last_point_new    = last_point_new    + 1_ip
          call metric%generate_edge_point(n1,n2,mesh%coord(:,n1),mesh%coord(:,n2),p_insert,M_p_insert)
          nodeId_toInsert = last_point_global
          isNodeId_existent = .false.
        end if
        
        if( do_checkQual ) q_prev = compute_minQuality_sizeShape_cavity(mesh%coord,mesh%lnods(:,elems_cavity),metric)
        
        if(isNodeId_existent) then
          call cavity_reinsert(    nodeId_toInsert ,mesh%lnods(:,elems_cavity), Tcav_remeshed)
        else
          call cavity_insert_point(nodeId_toInsert, mesh%lnods(:,elems_cavity), Tcav_remeshed)
        end if
        
        num_new_elems = size(Tcav_remeshed,2)
        isExceededAllocatedSize = (last_elem_new + num_new_elems)>max_newElems ! no need to check the nodes (exactly calculated)
        isValidRemesh = .not.isExceededAllocatedSize
      
        if(isValidRemesh) then ! check quality
          
          if(isNodeId_existent) then
            q_new = compute_minQuality_sizeShape_cavity(mesh%coord,Tcav_remeshed(:,:),metric)
          else
            coord_new(:,last_point_new) = p_insert
            call metric_new%set_metric_node( last_point_new , M_p_insert )
            q_new = compute_minQuality_sizeShape_cavity( mesh%coord,Tcav_remeshed(:,:),metric,&
              p_insert, metric_new%M(:,:,last_point_new) )
          end if
          
          isValidRemesh = (q_new>minQ_validity) !.and.(q_new>q_prev) !(q_new>q_prev/2.0_rp)
          if( do_checkQual ) isValidRemesh = (q_new>q_prev).and.isValidRemesh
          
        end if
      
        if(isValidRemesh) then
          
          isModifiedMesh = .true.
          numPerfEdges = numPerfEdges+1_ip
          isPerformedEdge(iedge)     = .true.
          isFrozenElem(elems_cavity) = .true.
          
          if(current_strategy==adapt_strat_collapse) then
            isRemovedNode(n1) = .true.
            isRemovedNode(n2) = .true.
            if(isNodeId_existent) isRemovedNode(nodeId_toInsert) = .false.
          end if
          
          ini_elem      = last_elem_new+1
          last_elem_new = last_elem_new+num_new_elems
          lnods_new(:,ini_elem:last_elem_new) = Tcav_remeshed
          call delete_Tcav_remeshed(Tcav_remeshed)
        else
          if(.not.isNodeId_existent) then
            last_point_global = last_point_global-1_ip
            last_point_new    = last_point_new   -1_ip
          end if
          if(isExceededAllocatedSize) then
            call node_to_elems%delete_elems(elems_cavity)
            go to 666
          end if
        end if
        
      end if
      
      call node_to_elems%delete_elems(elems_cavity)
    end if
  end do
  
  666 continue
  
  if(out_performance) then
    call cpu_time(t1)
    select case (current_strategy)
      case  (adapt_strat_collapse)
        deb_collLoop = deb_collLoop + (t1-t0)
      case  (adapt_strat_split)
        deb_splitLoop = deb_splitLoop + (t1-t0)
    end select
  end if
  
  call node_to_elems%deallo()
  if (isModifiedMesh) then
!     print*,kfl_paral,' edge_length -> outside mesh%nelem: ',mesh%nelem
!     print*,kfl_paral,' edge_length-> outside mesh%name: ',mesh%name
    call set_mesh_new_modify(mesh,coord_new,lnods_new,last_point_new,last_elem_new,&
      isDeletedElem=isFrozenElem,isNewElem=isNewElem)!,isRemovedNode)
    call metric%merge(metric_new,last_point_new)!,isRemovedNode-> do not remove here)
  else
    call metric_new%deallo() ! performed in metric%merge
    call memory_deallo(memor_adapt,'coord_new','set_tagged_edges',coord_new)! performed in set_mesh_new
    call memory_deallo(memor_adapt,'lnods_new','set_tagged_edges',lnods_new)! performed in set_mesh_new
  end if
  !
  call memory_deallo(memor_adapt,'isFrozenElem'    ,'repair_mesh_edgeLengths_frozen',isFrozenElem)
  return 
end subroutine repair_mesh_edgeLengths_frozen
!
!
!
function isAttemptedEdgeLengthRepair(n1,n2,current_strategy) result(isAttempted)
  implicit none
  
  integer(ip), intent(in) :: n1,n2,current_strategy
  logical(lg) :: isAttempted
  
  ! IMPLEMENT:
  ! check if edge swap exists
  ! TODO: if exists check if any of the surrounding nodes has been modified at some point -> is this easy?
  !    I need a track of the nodes of the mesh to know if at this current iteration the state of the node has changed
  
  !isAttempted = attempt_edgeSwap%isEdge(n1,n2)
  if(current_strategy.eq.adapt_strat_collapse) then
    isAttempted = attempt_edgeColl%isEdge(n1,n2)
  else
    isAttempted = attempt_edgeSplit%isEdge(n1,n2)
  end if
  
end function isAttemptedEdgeLengthRepair
!
!
!
subroutine setAttemptedEdgeLengthRepair(n1,n2,current_strategy)
  implicit none
  integer(ip), intent(in) :: n1,n2,current_strategy
  !
  if(current_strategy.eq.adapt_strat_collapse) then
    call attempt_edgeColl%setEdge(n1,n2)
  else
    call attempt_edgeSplit%setEdge(n1,n2)
  end if
  !
end subroutine setAttemptedEdgeLengthRepair
!
!
!
subroutine set_tagged_edges_byLength(current_strategy,tol_length,mesh,metric,numTaggedEdges,taggedEdges_to_node)
  use mod_meshEdges,  only: compute_mesh_edgeLengths
  use mod_debugTools, only: out_performance,deb_collTagEdges, deb_splitTagEdges, deb_splitEdgeLen, deb_collEdgeLen
  implicit none
  
  integer(ip),            intent(in)      :: current_strategy
  real(rp),               intent(in)      :: tol_length
  type(mesh_type_basic),  intent(in)      :: mesh
  type(mesh_metric_type), intent(in)      :: metric

  integer(ip), pointer,   intent(inout)   :: taggedEdges_to_node(:,:)
  integer(ip),            intent(out)     :: numTaggedEdges
  !
  real(rp),    pointer  :: edge_lengths(:)
  integer(ip), pointer  :: edge_to_node(:,:)
  integer(ip), pointer  :: perm_edges_byLength(:)
  integer(ip)           :: iedge, numEdges, theEdge
  real(rp)              :: l_edge
  !
  real(rp) :: t0,t1
  
  if(out_performance) call cpu_time(t0)
  
  call compute_mesh_edgeLengths(mesh,metric,edge_to_node,edge_lengths)
  numEdges = size(edge_to_node,2)
  
  if(out_performance) then
    call cpu_time(t1)
    select case (current_strategy)
      case  (adapt_strat_collapse)
        deb_collEdgeLen  = deb_collEdgeLen  + (t1-t0)
      case  (adapt_strat_split)
        deb_splitEdgeLen = deb_splitEdgeLen + (t1-t0)
    end select
  end if
  
  !
  ! SORTING
  !
  select case (current_strategy)
    case  (adapt_strat_collapse)
      call sort_perm_reals_increasing(edge_lengths,perm_edges_byLength)
    case  (adapt_strat_split)
      call sort_perm_reals_decreasing(edge_lengths,perm_edges_byLength)
    case default
      call runend('Nos implemented strategy in set_tagged_edges')
  end select

  loopEdges: do iedge=1,numEdges
    theEdge = perm_edges_byLength(iedge)
    l_edge = edge_lengths(theEdge)
    select case (current_strategy)
      case  (adapt_strat_collapse)
        if(l_edge>tol_length) exit loopEdges
      case  (adapt_strat_split)
        if(l_edge<tol_length) exit loopEdges
    end select
  end do loopEdges
  numTaggedEdges = iedge-1

  nullify(taggedEdges_to_node)
  call memory_alloca(memor_adapt,'taggedEdges_to_node','mod_adapt',taggedEdges_to_node,2_ip,numTaggedEdges)
  if(numTaggedEdges>0) then
    taggedEdges_to_node(:,:) = edge_to_node( : , perm_edges_byLength(1:numTaggedEdges) )
  end if
  
  call memory_deallo(memor_adapt,'perm_edges_byLength','mod_adapt',perm_edges_byLength)
  !
  ! WITHOUT SORTING
  !
!   nullify(list_edges)
!   call memory_alloca(memor_adapt,'list_edges','set_tagged_edges',list_edges,numEdges)
!   select case (current_strategy)
!     case  (adapt_strat_collapse)
!       !
!       numTaggedEdges = 0
!       do iedge=1,numEdges
!         if( edge_lengths(iedge) < tol_length ) then
!           numTaggedEdges = numTaggedEdges+1
!           list_edges(numTaggedEdges) = iedge
!         end if
!       end do
!       !
!     case  (adapt_strat_split)
!       !
!       numTaggedEdges = 0
!       do iedge=1,numEdges
!         if( edge_lengths(iedge) > tol_length ) then
!           numTaggedEdges = numTaggedEdges+1
!           list_edges(numTaggedEdges) = iedge
!         end if
!       end do
!       !
!   end select
!
!   nullify(taggedEdges_to_node)
!   if(numTaggedEdges>0) then
!     call memory_alloca(memor_adapt,'tagged_edges','set_tagged_edges',taggedEdges_to_node,2_ip,numTaggedEdges)
!     taggedEdges_to_node(:,:) = edge_to_node( : , list_edges(1:numTaggedEdges) )
!   end if
!   call memory_deallo(memor_adapt,'list_edges'  ,'set_tagged_edges',list_edges)
  !
  ! Deallocate and nullify
  call memory_deallo(memor_adapt,'edge_to_node','set_tagged_edges',edge_to_node)
  call memory_deallo(memor_adapt,'edge_lengths','set_tagged_edges',edge_lengths)
  
  if(out_performance) then
    call cpu_time(t1)
    select case (current_strategy)
      case  (adapt_strat_collapse)
        deb_collTagEdges = deb_collTagEdges + (t1-t0)
      case  (adapt_strat_split)
        deb_splitTagEdges = deb_splitTagEdges + (t1-t0)
    end select
  end if
  
end subroutine set_tagged_edges_byLength
!
!
!
!--------- SUBMODULE 2: _adapt_swap ---------------------------------------------------------------------------
! include "_adapt_swap.f90"
!
!
!
!-----------------------------------------------------------------------
!
!> @author  Abel Gargallo
!> @date    28/04/2021
!> @brief   swap_mesh_edges
!> @details swap_mesh_edges 
!>          Process:
!>          - while the mesh is modified
!>          - try to swap mesh edges using swap_mesh_edges_frozenTaggedEdges
!
!-----------------------------------------------------------------------
subroutine swap_mesh_edges(mesh,metric,tol_qual,mapNodes_new_to_old,edgeData_boun,isBouNode_iniMesh)
  use mod_meshTopology,       only: setIsBoundaryNodeFromIniMesh, deleteIsBouNode
  use mod_quality,    only: compute_mesh_quality_sizeShape, quality_deallo
  use mod_debugTools,   only: out_performance, deb_swapLoop, deb_swapRest, deb_swapIni, deb_tagEdgesByQ,&
    deb_swapPreTagQ,deb_swapPreTagNE
  implicit none
  !
  type(mesh_type_basic),  intent(inout)   :: mesh
  type(mesh_metric_type), intent(in)      :: metric
  real(rp),               intent(in)      :: tol_qual 
  integer(ip), pointer,   intent(inout)   :: mapNodes_new_to_old(:)
  type(edgeData_type),    intent(in)      :: edgeData_boun
  logical(lg), pointer                    :: isBouNode_iniMesh(:)
  !
  integer(ip), parameter :: max_iter_frozen = 1000
  integer(ip), parameter :: numNieghLevels = 1
  !
  logical(lg) :: isModifiedMesh
  integer(ip) :: iter
  !
  logical(lg), pointer :: isBoundaryNode(:)
  real(rp),    pointer :: q(:)
  !
  logical(lg), pointer :: isNewElem(:)
  integer(ip) :: numElems_fixed
  integer(ip), pointer :: lnods_fixed(:,:)
  real(rp)   , pointer :: q_fixed(:)
  logical(lg)          :: do_useSubmesh
  !
  integer(ip) :: iout = -1_ip
  !
  real(rp) :: t0,t1
  !
  if(output_steps_swap_paraview) call out_with_q(mesh,metric,file_out//"_swap",iout)
  !
  if(.not.attempt_edgeSwap%isAlloca()) call attempt_edgeSwap%alloca( INT(size(mesh%lnods,1)*size(mesh%lnods,2),ip) )

  call setIsBoundaryNodeFromIniMesh(isBoundaryNode,isBouNode_iniMesh,mapNodes_new_to_old)
  !
  if(out_performance) call cpu_time(t0)
  call compute_mesh_quality_sizeShape(mesh,metric,q)
  if(out_performance) then
    call cpu_time(t1)
    deb_swapPreTagQ = deb_swapPreTagQ + (t1-t0)
  end if
  !
  isModifiedMesh = .true.
  iter = 0
  do while ( isModifiedMesh.and.(iter<max_iter_frozen) )
    !
    iter = iter+1
    do_useSubmesh = iter>1 ! in the first iteration, check all mesh
    !
    if(do_useSubmesh) then
      call set_submesh(mesh,numNieghLevels,isNewElem(:),numElems_fixed,lnods_fixed,q=q,q_fixed=q_fixed)
      call memory_deallo(memor_adapt,'isNewElem','mod_adapt',isNewElem)
    end if
    !
    call swap_mesh_edges_frozenTaggedEdges(mesh,metric,isBoundaryNode,tol_qual,isModifiedMesh,q,isNewElem)
    !
    if(do_useSubmesh) then
      call rebuild_mesh_fromSubmesh(mesh,numElems_fixed,lnods_fixed,q=q,q_fixed=q_fixed,isNewElem=isNewElem)
    end if
    !
    if(output_steps_swap_paraview) call out_with_q(mesh,metric,file_out//"_swap",iout)
  end do
  !
  call quality_deallo(q)
  !
  call deleteIsBouNode(isBoundaryNode)
  !
end subroutine swap_mesh_edges
!
!
!
!-----------------------------------------------------------------------
!
!> @author  Abel Gargallo
!> @date    28/04/2021
!> @brief   swap_mesh_edges_frozen
!> @details swap_mesh_edges_frozen 
!>          (frozen means checks if that elements exist to reuse topological computed date before recomputing it)
!>          Process:
!>          - sort edges according to lowest adjacent elem quality
!>          - swap the selected mesh edges using swap_mesh_edges_frozenTopo
!
!-----------------------------------------------------------------------
subroutine swap_mesh_edges_frozenTaggedEdges(mesh,metric,isBoundaryNode,tol_qual,isModifiedMesh,q,isNewElem)
  use mod_meshTopology, only: node_to_elems_type
  use mod_quality,      only: compute_mesh_quality_sizeShape, quality_deallo
  use mod_debugTools,   only: out_performance, deb_swapLoop, deb_swapRest, deb_swapIni, deb_tagEdgesByQ,&
    deb_swapPreTagQ,deb_swapPreTagNE
  !
  implicit none
  !
  type(mesh_type_basic),  intent(inout) :: mesh
  type(mesh_metric_type), intent(in)    :: metric
  real(rp),               intent(in)    :: tol_qual 
  logical(lg), pointer,   intent(in)    :: isBoundaryNode(:)
  logical(lg),            intent(out)   :: isModifiedMesh
  real(rp),    pointer,   intent(inout) :: q(:)
  logical(lg), pointer,   intent(inout) :: isNewElem(:)
  !
  integer(ip), parameter   :: maxIter_fixTopo = 100_ip
  
  type(node_to_elems_type) :: node_to_elems
  integer(ip)              :: numTaggedEdges
  integer(ip), pointer     :: taggedEdges_to_node(:,:)
  real(rp)                 :: t0,t1,t0bis,t1bis
  logical(lg)              :: isModifiedMesh_inner, isFrozenLock, doRetrySwapping
  
  integer(ip)              :: iter_fixTopo
  integer(ip)              :: ii
  !
  if(out_performance) call cpu_time(t0bis)
  if(out_performance) call cpu_time(t0)
  
  call node_to_elems%set(mesh) 
  
  if(out_performance) then
    call cpu_time(t1)
    deb_swapPreTagNE = deb_swapPreTagNE + (t1-t0)
  end if
  !
  if(out_performance) call cpu_time(t0)
  call set_tagged_edges_byQ(tol_qual,mesh,node_to_elems,q,numTaggedEdges,taggedEdges_to_node)
  if(out_performance) then
    call cpu_time(t1)
    deb_tagEdgesByQ = deb_tagEdgesByQ + (t1-t0   )
    deb_swapIni     = deb_swapIni     + (t1-t0bis)
  end if
  !
  nullify(isNewElem)
  call memory_alloca(memor_adapt,'isNewElem','mod_adapt',isNewElem,mesh%nelem)
  do ii=1,mesh%nelem
    isNewElem(ii) = .false. !isNewElem(:) = .false.
  end do
  
  isModifiedMesh        = .false.
  doRetrySwapping       = .true.
  iter_fixTopo          =  0_ip
  do while(doRetrySwapping.and.(iter_fixTopo<maxIter_fixTopo))
    iter_fixTopo = iter_fixTopo+1
    !
    call swap_mesh_edges_frozenTopo(mesh,metric,isBoundaryNode,tol_qual,&
      numTaggedEdges,taggedEdges_to_node,node_to_elems,&
      isModifiedMesh_inner,q,isFrozenLock,isNewElem)
    !
    if(isModifiedMesh_inner) isModifiedMesh = .true.
    doRetrySwapping = isModifiedMesh_inner.and.isFrozenLock ! if is modified but not frozen, it can stop looping
    !
  end do
  !
  call node_to_elems%deallo()
  call memory_deallo(memor_adapt,'taggedEdges_to_node','mod_adapt',taggedEdges_to_node)
  !
end subroutine swap_mesh_edges_frozenTaggedEdges
!-----------------------------------------------------------------------
!
!> @author  Abel Gargallo
!> @date    28/04/2021
!> @brief   swap_mesh_edges_frozen
!> @details swap_mesh_edges_frozen 
!>          (frozen means checks if that elements exist to reuse topological computed date before recomputing it)
!>          Process:
!>          - fixed the topological information between elements (frozing deleted elements),
!>          - for each edge
!>          -   compute cavity of the edge
!>          -   for each boundary node of the cavity
!>          -     reinsert the node
!>          -     compute quality of the new cavity
!>          -     save the swap if it improves current quality
!>          - that is perform the reinsertion of the _edge and cavity _node with highest _q  (_q,_edge,_node)
!
!-----------------------------------------------------------------------
subroutine swap_mesh_edges_frozenTopo(mesh,metric,isBoundaryNode,tol_qual,&
      numTaggedEdges,taggedEdges_to_node,node_to_elems,&
      isModifiedMesh,q,isFrozenLock,isNewElem)
  use mod_meshTopology, only: node_to_elems_type
  use mod_quality,      only: compute_mesh_quality_sizeShape, quality_deallo
  use mod_debugTools,   only: out_performance, deb_swapLoop, deb_swapRest, deb_swapIni, deb_tagEdgesByQ,&
    deb_swapPreTagQ,deb_swapPreTagNE,deb_swapSetMesh
  use def_master,    only: kfl_paral
  !
  implicit none
  !
  type(mesh_type_basic),   intent(inout) :: mesh
  type(mesh_metric_type),  intent(in)    :: metric
  real(rp),                intent(in)    :: tol_qual 
  logical(lg), pointer,    intent(in)    :: isBoundaryNode(:)
  integer(ip)         ,    intent(in)    :: numTaggedEdges
  integer(ip), pointer,    intent(in)    :: taggedEdges_to_node(:,:)
  type(node_to_elems_type),intent(inout) :: node_to_elems
  logical(lg),             intent(out)   :: isModifiedMesh
  real(rp),    pointer,    intent(inout) :: q(:)
  logical(lg),             intent(out)   :: isFrozenLock
  logical(lg), pointer,    intent(inout) :: isNewElem(:)
  !
  integer(ip)              :: max_newElems, count_newElems, count_newPoints
  integer(ip), pointer     :: lnods_new(:,:)
  real(rp)   , pointer     :: q_new(:)
  logical(lg), pointer     :: isFrozenElem(:)
  integer(ip)              :: iedge, edgeNodes(2)
  real(rp),    pointer     :: coords_new_fake(:,:)
  real(rp)                 :: t0,t1
  !
  if(mesh%ndime==2) then ! each node has an adjacency approx of 6 elems, 2 nodes per edge
    max_newElems = numTaggedEdges*6_ip*2_ip
  else ! each node has an adjacency approx of 16 elems, 2 nodes per edge
    max_newElems = numTaggedEdges*16_ip*2_ip
  end if
  !
  count_newElems        =  0_ip
  nullify(lnods_new)
  call memory_alloca(memor_adapt,'lnods_new'   ,'swap_edges',lnods_new,int(size(mesh%lnods,1),ip),max_newElems)
  nullify(q_new)
  call memory_alloca(memor_adapt,'q_new'       ,'swap_edges',q_new,max_newElems) 
  nullify(isFrozenElem)
  call memory_alloca(memor_adapt,'isFrozenElem','swap_edges',isFrozenElem,mesh%nelem)
  
  
  if(out_performance) call cpu_time(t0)
  isModifiedMesh          = .false.
  isFrozenElem(:)         = .false.
  isFrozenLock            = .false.
  do iedge=1,numTaggedEdges
    !
    edgeNodes = taggedEdges_to_node(:,iedge)
    !
    call swap_edge(edgeNodes,                       &
      mesh,metric,node_to_elems,isBoundaryNode,q,   &
      isModifiedMesh,isFrozenElem,                  &
      count_newElems,lnods_new,q_new,               &
      isFrozenLock)
    !
  end do
  !
  if(out_performance) then
    call cpu_time(t1)
    deb_swapLoop = deb_swapLoop + (t1-t0)
  end if
  !
  if( isModifiedMesh ) then
    if(count_newElems>0_ip) then
      !
      call cpu_time(t0)
      count_newPoints = 0_ip
      nullify(coords_new_fake)
      call memory_alloca(memor_adapt,'coords_new_fake','swap_edges',coords_new_fake,mesh%ndime,2_ip)
  
!       print*,kfl_paral,' edge_swap -> outside size(isFrozenElem)',size(isFrozenElem)
!       print*,kfl_paral,' edge_swap -> outside count(isFrozenElem)',count(isFrozenElem)
!       print*,kfl_paral,' edge_swap -> outside mesh%nelem: ',mesh%nelem
!       print*,kfl_paral,' edge_swap -> outside mesh%name:  ',mesh%name
!       print*,kfl_paral,' count_newPoints,count_newElems :', count_newPoints,'  ',count_newElems
!       print*,kfl_paral,' size(coords_new_fake): ',size(coords_new_fake)
!       print*,kfl_paral,' size(lnods_new):       ',size(lnods_new)
!       print*,kfl_paral,' size(q_new):           ',size(q_new)
!       print*,kfl_paral,' size(q):               ',size(q)
!       print*,kfl_paral,' size(isNewElem):       ',size(isNewElem)
!       print*,kfl_paral,' -> going inside set_mesh_new_modify!'
      call set_mesh_new_modify(mesh,coords_new_fake,lnods_new,count_newPoints,count_newElems,&
        isDeletedElem=isFrozenElem,q=q,q_new=q_new,isNewElem=isNewElem)
  
      call memory_deallo(memor_adapt,'coords_new_fake','swap_edges',coords_new_fake)
      if(out_performance) then
        call cpu_time(t1)
        deb_swapSetMesh = deb_swapSetMesh + (t1-t0)
      end if
      !
      call cpu_time(t0)
      call node_to_elems%deallo()
      call node_to_elems%set(mesh)
      if(out_performance) then
        call cpu_time(t1)
        deb_swapPreTagNE = deb_swapPreTagNE + (t1-t0)
      end if
      !
    end if
  else
    call memory_deallo(memor_adapt,'lnods_new'      ,'swap_edges',lnods_new)
    call memory_deallo(memor_adapt,'q_new'          ,'swap_edges',    q_new)
  end if
  
  call memory_deallo(memor_adapt,'isFrozenElem'    ,'swap_edges',isFrozenElem)
  !
end subroutine swap_mesh_edges_frozenTopo
!
!
!
subroutine swap_edge(edge,mesh,metric,node_to_elems,isBoundaryNode,q,&
  isModifiedMesh,isFrozenElem,count_newElems,lnods_new,q_new,isFrozenLock)
  use mod_meshTopology, only: node_to_elems_type, getNodesCav, deleteNodesCav
  use mod_meshTopology, only: get_boundary_faces, delete_boundary_faces
  use mod_quality,      only: compute_minQuality_sizeShape_cavity, compute_quality_sizeShape_cavity
  use mod_cavity,       only: cavity_reinsert, cavity_reinsert_fromCavBoun
  implicit none 
  !
  integer(ip),              intent(in)    :: edge(2)
  type(mesh_type_basic),    intent(in)    :: mesh
  type(mesh_metric_type),   intent(in)    :: metric
  type(node_to_elems_type), intent(in)    :: node_to_elems
  logical(lg), pointer,     intent(in)    :: isBoundaryNode(:)
  real(rp),                 intent(in)    :: q(:)
  logical(lg),              intent(inout) :: isModifiedMesh
  logical(lg),              intent(inout) :: isFrozenElem(:)
  integer(ip),              intent(inout) :: count_newElems
  integer(ip), pointer,     intent(inout) :: lnods_new(:,:)
  real(rp)   , pointer,     intent(inout) ::     q_new(:)
  logical(lg),              intent(inout) :: isFrozenLock
  !
  integer(ip)          :: n1, n2
  integer(ip), pointer :: elems_cavity(:), nodesCav(:)
  integer(ip)          :: numNodCav
  real(rp)             :: q_bestSwap, q_current
  real(rp)   , pointer :: qvec_bestSwap(:),qcav_remeshed(:)
  integer(ip), pointer :: T_bestSwap(:,:), Tcav_remeshed(:,:), boundary_faces(:,:)
  integer(ip)          :: knode, pointToReinsert
  integer(ip)          :: newCavElems, newElemsBestSwap
  logical(lg)          :: isImproved, isBoundEdge, isFrozenCavity
  logical(lg)          :: isExceededAllocatedSize
  !
  n1 = edge(1)
  n2 = edge(2)
  !
  isBoundEdge = isBoundaryNode(n1).and.isBoundaryNode(n2)
  if(isBoundEdge           ) return !cycle
  !
  if(isAttemptedSwap(n1,n2)) return !cycle
  !
  call node_to_elems%get_edgeContainers(n1,n2,elems_cavity)
  isFrozenCavity = any(isFrozenElem(elems_cavity))
  if( .not.isFrozenCavity ) then
    !
    call setAttemptedSwap(n1,n2)
    !
    call getNodesCav(       mesh%lnods(:,elems_cavity),numNodCav,nodesCav)
    call get_boundary_faces(mesh%lnods(:,elems_cavity),boundary_faces    )
    !
    nullify(Tcav_remeshed)
    call memory_alloca(memor_adapt,'Tcav_remeshed','swap_edge',Tcav_remeshed,&
      int(size(boundary_faces,1)+1,ip), int(size(boundary_faces,2),ip) )
    nullify(qcav_remeshed)
    call memory_alloca(memor_adapt,'qcav_remeshed','swap_edge',qcav_remeshed,int(size(boundary_faces,2),ip) )
    nullify(T_bestSwap)
    call memory_alloca(memor_adapt,'T_bestSwap','swap_edge',T_bestSwap,&
      int(size(Tcav_remeshed,1)   ,ip), int(size(Tcav_remeshed,2) ,ip) )
    nullify(qvec_bestSwap)
    call memory_alloca(memor_adapt,'qvec_bestSwap','swap_edge',qvec_bestSwap,int(size(Tcav_remeshed,2) ,ip) )
    !
    q_bestSwap = minval(q(elems_cavity)) !-> we know they are not frozen (these elements are unmodified)
    isImproved = .false.
    do knode =1,numNodCav
      !
      pointToReinsert = nodesCav(knode)
      !
      if((pointToReinsert.ne.n1).and.(pointToReinsert.ne.n2)) then
        !
        call cavity_reinsert_fromCavBoun(pointToReinsert,boundary_faces,Tcav_remeshed,newCavElems)
        !
        ! aquesta lhem dusar per guardar la nova qualitat i estalviar de recalcularla TODO XXX
        !q_current = compute_minQuality_sizeShape_cavity(mesh%coord,Tcav_remeshed(:,1:newCavElems),metric)
        call compute_quality_sizeShape_cavity(qcav_remeshed(1:newCavElems),q_current,&
          mesh%coord,Tcav_remeshed(:,1:newCavElems),metric)
        if( q_current > q_bestSwap + EPSILON(q_bestSwap) ) then
          isImproved       = .true.
          q_bestSwap       = q_current
          newElemsBestSwap = newCavElems
          T_bestSwap(   :,1:newElemsBestSwap) = Tcav_remeshed(:,1:newCavElems)
          qvec_bestSwap(  1:newElemsBestSwap) = qcav_remeshed(  1:newCavElems)
        end if
        !
        if(mesh%ndime==2) EXIT ! only one node must be reinserted in 2D
        !
      end if
    end do
    !
    if( isImproved ) then
      !
      isModifiedMesh          = .true.
      !
      isExceededAllocatedSize = (count_newElems+newElemsBestSwap) > size(lnods_new,2)
      if( isExceededAllocatedSize ) call amplify_lnods(lnods_new,count_newElems,q_new)
      !
      isFrozenElem(elems_cavity) = .true.
      lnods_new(:,(count_newElems+1):(count_newElems+newElemsBestSwap)) =    T_bestSwap(:,1:newElemsBestSwap)
      q_new    (  (count_newElems+1):(count_newElems+newElemsBestSwap)) = qvec_bestSwap(  1:newElemsBestSwap)
      count_newElems = count_newElems+newElemsBestSwap
      !
    end if
    !
    call memory_deallo(memor_adapt,'Tcav_remeshed','swap_edge',Tcav_remeshed)
    call memory_deallo(memor_adapt,'qcav_remeshed','swap_edge',qcav_remeshed)
    call memory_deallo(memor_adapt,'T_bestSwap'   ,'swap_edge',   T_bestSwap)
    call memory_deallo(memor_adapt,'qvec_bestSwap','swap_edge',qvec_bestSwap)
    call deleteNodesCav(nodesCav)
    call delete_boundary_faces(boundary_faces)
    !
  else
    isFrozenLock = .true.
  end if
  !
  call node_to_elems%delete_elems(elems_cavity)
  !
end subroutine swap_edge
!
!
!
!-----------------------------------------------------------------------
!
!> @author  Abel Gargallo
!> @date    28/04/2021
!> @brief   swap_mesh_edges_frozen
!> @details swap_mesh_edges_frozen 
!>          (frozen means checks if that elements exist to reuse topological computed date before recomputing it)
!>          Process:
!>          - compute low quality elems
!>          - sort edges according to lowest adjacent elem quality
!>          - for each edge
!>          -   compute cavity of the edge
!>          -   for each boundary node of the cavity
!>          -     reinsert the node
!>          -     compute quality of the new cavity
!>          -     save the swap if it improves current quality
!>          - that is perform the reinsertion of the _edge and cavity _node with highest _q  (_q,_edge,_node)
!
!-----------------------------------------------------------------------
subroutine swap_mesh_edges_frozen_separedSwapEdge(mesh,metric,isBoundaryNode,tol_qual,isModifiedMesh,q)
  use mod_meshTopology, only: node_to_elems_type
  use mod_quality,      only: compute_mesh_quality_sizeShape, quality_deallo
  use mod_debugTools,   only: out_performance, deb_swapLoop, deb_swapRest, deb_swapIni, deb_tagEdgesByQ,&
    deb_swapPreTagQ,deb_swapPreTagNE
  !
  implicit none
  !
  type(mesh_type_basic),  intent(inout) :: mesh
  type(mesh_metric_type), intent(in)    :: metric
  real(rp),               intent(in)    :: tol_qual 
  logical(lg), pointer,   intent(in)    :: isBoundaryNode(:)
  logical(lg),            intent(out)   :: isModifiedMesh
  real(rp),    pointer,   intent(inout) :: q(:)
  !
  type(node_to_elems_type) :: node_to_elems
  logical(lg), pointer     :: isFrozenElem(:)
  integer(ip)              :: numTaggedEdges
  integer(ip), pointer     :: taggedEdges_to_node(:,:)
  integer(ip), pointer     :: lnods_new(:,:)
  real(rp)   , pointer     :: q_new(:)
  integer(ip)              :: max_newElems, count_newElems, count_newPoints
  integer(ip)              :: iedge, edgeNodes(2)
  real(rp),    pointer     :: coords_new_fake(:,:)
  real(rp)                 :: t0,t1,t0bis,t1bis
  logical(lg)              :: isFrozenLock
  !
  if(out_performance) call cpu_time(t0bis)
  !
  if(out_performance) call cpu_time(t0)
  call node_to_elems%set(mesh)  ! NO COST
  if(out_performance) then
    call cpu_time(t1)
    deb_swapPreTagNE = deb_swapPreTagNE + (t1-t0)
  end if
  !
  if(out_performance) call cpu_time(t0)
  call set_tagged_edges_byQ(tol_qual,mesh,node_to_elems,q,numTaggedEdges,taggedEdges_to_node)
  if(out_performance) then
    call cpu_time(t1)
    deb_tagEdgesByQ = deb_tagEdgesByQ + (t1-t0)
  end if
  !
  if(mesh%ndime==2) then ! each node has an adjacency approx of 6 elems, 2 nodes per edge
    max_newElems = numTaggedEdges*6_ip*2_ip
  else ! each node has an adjacency approx of 16 elems, 2 nodes per edge
    max_newElems = numTaggedEdges*16_ip*2_ip
  end if
  nullify(lnods_new)
  call memory_alloca(memor_adapt,'lnods_new','swap_edges',lnods_new,&
    int(size(mesh%lnods,1),ip),max_newElems)  
  nullify(q_new)
  call memory_alloca(memor_adapt,'q_new'    ,'swap_edges',q_new,max_newElems)  
  !
  if(out_performance) then
    call cpu_time(t1bis)
    deb_swapIni = deb_swapIni + (t1bis-t0bis)
  end if
  !
  if(out_performance) call cpu_time(t0)
  count_newElems          =  0_ip
  isModifiedMesh          = .false.
  nullify(isFrozenElem)
  call memory_alloca(memor_adapt,'isFrozenElem'    ,'swap_edges',isFrozenElem,mesh%nelem)  
  isFrozenElem(:)         = .false.
  do iedge=1,numTaggedEdges
    !
    edgeNodes = taggedEdges_to_node(:,iedge)
    !
    call swap_edge(edgeNodes,                       &
      mesh,metric,node_to_elems,isBoundaryNode,q,   &
      isModifiedMesh,isFrozenElem,                  &
      count_newElems,lnods_new,q_new,               &
      isFrozenLock)
    !
  end do
  !
  if(out_performance) then
    call cpu_time(t1)
    deb_swapLoop = deb_swapLoop + (t1-t0)
  end if
  !
  if( isModifiedMesh ) then
    if(count_newElems>0_ip) then
      count_newPoints = 0_ip
      nullify(coords_new_fake)
      call memory_alloca(memor_adapt,'coords_new_fake','swap_edges',coords_new_fake,mesh%ndime,2_ip)
      call set_mesh_new_modify(mesh,coords_new_fake,lnods_new,count_newPoints,count_newElems,&
        isDeletedElem=isFrozenElem,q=q,q_new=q_new)
      call memory_deallo(memor_adapt,'coords_new_fake','swap_edges',coords_new_fake)
    end if
  else
    call memory_deallo(memor_adapt,'lnods_new'      ,'swap_edges',lnods_new)
    call memory_deallo(memor_adapt,'q_new'          ,'swap_edges',    q_new)
  end if
  !
  call memory_deallo(memor_adapt,'isFrozenElem'    ,'swap_edges',isFrozenElem)
  call memory_deallo(memor_adapt,'taggedEdges_to_node','mod_adapt',taggedEdges_to_node)
  !
  if(out_performance) then
    call cpu_time(t1bis)
    deb_swapRest = deb_swapRest + (t1bis-t0bis) - (t1-t0)
  end if
end subroutine swap_mesh_edges_frozen_separedSwapEdge
!
!
!
subroutine amplify_lnods(lnods_new,count_newElems,q_new)
  implicit none

  integer(ip), pointer,     intent(inout) :: lnods_new(:,:)
  integer(ip),              intent(in   ) :: count_newElems
  real(rp)   , pointer,     intent(inout) ::     q_new(:)

  integer(ip), pointer :: lnods_new_copy(:,:)
  real(rp)   , pointer :: q_new_copy(:)
  !
  nullify(lnods_new_copy)
  call memory_alloca(memor_adapt,'lnods_new_copy','amplify_lnods',lnods_new_copy,&
    int(size(lnods_new,1),ip),count_newElems)
  lnods_new_copy = lnods_new(:,1:count_newElems)
  call memory_deallo(memor_adapt,'lnods_new','swap_edges',lnods_new)
  nullify(lnods_new)
  call memory_alloca(memor_adapt,'lnods_new','swap_edges',lnods_new,&
    int(size(lnods_new_copy,1),ip),count_newElems*2)
  lnods_new(:,1:count_newElems) = lnods_new_copy
  call memory_deallo(memor_adapt,'lnods_new_copy','amplify_lnods',lnods_new_copy)
  !
  nullify(q_new_copy)
  call memory_alloca(memor_adapt,'q_new_copy','amplify_lnods',q_new_copy,count_newElems)
  q_new_copy = q_new(1:count_newElems)
  call memory_deallo(memor_adapt,'q_new','swap_edges',q_new)
  nullify(q_new)
  call memory_alloca(memor_adapt,'q_new','swap_edges',q_new,count_newElems*2)
  q_new(1:count_newElems) = q_new_copy
  call memory_deallo(memor_adapt,'q_new_copy','amplify_lnods',q_new_copy)
  !
end subroutine amplify_lnods
!
!
!
!-----------------------------------------------------------------------
!
!> @author  Abel Gargallo
!> @date    28/04/2021
!> @brief   swap_mesh_edges_frozen_allInOne
!> @details swap_mesh_edges_frozen_allInOne 
!>          (frozen means checks if that elements exist to reuse topological computed date before recomputing it)
!>          Process:
!>          - compute low quality elems
!>          - sort edges according to lowest adjacent elem quality
!>          - for each edge
!>          -   compute cavity of the edge
!>          -   for each boundary node of the cavity
!>          -     reinsert the node
!>          -     compute quality of the new cavity
!>          -     save the swap if it improves current quality
!>          - that is perform the reinsertion of the _edge and cavity _node with highest _q  (_q,_edge,_node)
!
!-----------------------------------------------------------------------
subroutine swap_mesh_edges_frozen_allInOne(mesh,metric,isBoundaryNode,tol_qual,isModifiedMesh)
  use mod_meshTopology, only: node_to_elems_type, getNodesCav, deleteNodesCav, delete_arrayElems
  use mod_quality,      only: compute_mesh_quality_sizeShape, quality_deallo
  use mod_quality,      only: compute_minQuality_sizeShape_cavity, compute_minQuality_sizeshape_subset
  use mod_cavity,       only: cavity_reinsert, cavity_reinsert_fromCavBoun
  use mod_meshTopology, only: get_boundary_faces, delete_boundary_faces
  use mod_debugTools,   only: out_performance, deb_swapLoop, deb_swapRest, deb_swapIni, deb_tagEdgesByQ,&
    deb_swapPreTagQ,deb_swapPreTagNE
  
  implicit none
  type(mesh_type_basic),  intent(inout) :: mesh
  type(mesh_metric_type), intent(in)    :: metric
  real(rp),               intent(in)    :: tol_qual 
  logical(lg), pointer,   intent(in)    :: isBoundaryNode(:)
  logical(lg),            intent(out)   :: isModifiedMesh

  type(node_to_elems_type) :: node_to_elems
  integer(ip) :: iedge, n1, n2
  logical(lg), pointer :: isFrozenElem(:)
  integer(ip), pointer :: elems_cavity(:), nodesCav(:)
  real(rp),    pointer :: q(:)
  integer(ip) :: numNodCav
  
  real(rp) :: q_bestSwap, q_current
  integer(ip), pointer :: T_bestSwap(:,:), Tcav_remeshed(:,:), lnods_new(:,:)
  integer(ip) :: knode, pointToReinsert
  logical(lg) :: isImproved, isBoundEdge
  integer(ip)         :: max_newElems, count_newElems, count_newPoints
  real(rp), pointer   :: coords_new_fake(:,:)
  
  integer(ip) :: numTaggedEdges
  integer(ip), pointer :: taggedEdges_to_node(:,:)
  logical(lg) :: isFrozenCavity
  real(rp)    :: eps_dim
  
  integer(ip), pointer :: boundary_faces(:,:)
  integer(ip) :: newCavElems
  integer(ip) :: newElemsBestSwap
  real(rp) :: t0,t1,t0bis,t1bis
  
  eps_dim = 0.1_rp**(14_rp/real(mesh%ndime,rp)) !1.0/( huge(1.0_rp)**(1.0/mesh%ndime) ) 
  !
  if(out_performance) call cpu_time(t0bis)
  !
  if(out_performance) call cpu_time(t0)
  call node_to_elems%set(mesh)  ! NO COST
  if(out_performance) then
    call cpu_time(t1)
    deb_swapPreTagNE = deb_swapPreTagNE + (t1-t0)
  end if
  if(out_performance) call cpu_time(t0)
  call compute_mesh_quality_sizeShape(mesh,metric,q)! aquesta ens la podem estalviar TODO XXX
  if(out_performance) then
    call cpu_time(t1)
    deb_swapPreTagQ = deb_swapPreTagQ + (t1-t0)
  end if

  if(out_performance) call cpu_time(t0)
  call set_tagged_edges_byQ(tol_qual,mesh,node_to_elems,q,numTaggedEdges,taggedEdges_to_node)
  if(out_performance) then
    call cpu_time(t1)
    deb_tagEdgesByQ = deb_tagEdgesByQ + (t1-t0)
  end if
  
  if(mesh%ndime==2) then ! each node has an adjacency approx of 6 elems, 2 nodes per edge
    max_newElems = numTaggedEdges*6_ip*2_ip
  else ! each node has an adjacency approx of 16 elems, 2 nodes per edge
    max_newElems = numTaggedEdges*16_ip*2_ip
  end if
  nullify(lnods_new)
  call memory_alloca(memor_adapt,'lnods_new','repair_mesh_edgeLengths_frozen',lnods_new,int(size(mesh%lnods,1),ip),max_newElems)  

  if(out_performance) then
    call cpu_time(t1bis)
    deb_swapIni = deb_swapIni + (t1bis-t0bis)
  end if
  
  if(out_performance) call cpu_time(t0)
  isModifiedMesh  = .false.
  nullify(isFrozenElem)
  call memory_alloca(memor_adapt,'isFrozenElem','repair_mesh_edgeLengths_frozen',isFrozenElem,mesh%nelem)
  isFrozenElem(:) = .false.
  count_newElems  = 0_ip
  loopEdges: do iedge=1,numTaggedEdges
    !
    n1 = taggedEdges_to_node(1,iedge)
    n2 = taggedEdges_to_node(2,iedge)
    !
    isBoundEdge = isBoundaryNode(n1).and.isBoundaryNode(n2)
    if(isBoundEdge) cycle
    if(isAttemptedSwap(n1,n2)) cycle
    !
    call node_to_elems%get_edgeContainers(n1,n2,elems_cavity)
    isFrozenCavity = any(isFrozenElem(elems_cavity))
    if( .not.isFrozenCavity ) then
      
      call setAttemptedSwap(n1,n2)
      
      call getNodesCav(       mesh%lnods(:,elems_cavity),numNodCav,nodesCav)
      call get_boundary_faces(mesh%lnods(:,elems_cavity),boundary_faces)
      
      nullify(Tcav_remeshed)
      call memory_alloca(memor_adapt,'Tcav_remeshed','cavity_insert_point',Tcav_remeshed,&
        int(size(boundary_faces,1)+1,ip), int(size(boundary_faces,2),ip) )

      nullify(T_bestSwap)
      call memory_alloca(memor_adapt,'T_bestSwap','swap_mesh_edges_frozen',T_bestSwap,int(size(Tcav_remeshed,1),ip),int(size(Tcav_remeshed,2),ip))
      
      !q_bestSwap = compute_minQuality_sizeshape_subset(mesh,metric,elems_cavity) !current configuration
      q_bestSwap = minval(q(elems_cavity)) !-> we know they are not frozen (these elements are unmodified)
      isImproved = .false.
      do knode =1,numNodCav
        
        pointToReinsert = nodesCav(knode)

        if((pointToReinsert.ne.n1).and.(pointToReinsert.ne.n2)) then
          
          call cavity_reinsert_fromCavBoun(pointToReinsert,boundary_faces,Tcav_remeshed,newCavElems)
          
          ! aquesta lhem dusar per guardar la nova qualitat i estalviar de recalcularla TODO XXX
          q_current = compute_minQuality_sizeShape_cavity(mesh%coord,Tcav_remeshed(:,1:newCavElems),metric)
          if( q_current > q_bestSwap + eps_dim) then
            isImproved = .true.
            q_bestSwap = q_current
            newElemsBestSwap = newCavElems
            T_bestSwap(:,1:newElemsBestSwap) = Tcav_remeshed(:,1:newCavElems)
          end if
        end if
      end do
      
      if( isImproved ) then
        isModifiedMesh = .true.
        if( (count_newElems+size(T_bestSwap,2)) > max_newElems ) go to 6666 ! exceeded allocated size
        
        isFrozenElem(elems_cavity) = .true.
        lnods_new(:,(count_newElems+1):(count_newElems+newElemsBestSwap)) = T_bestSwap(:,1:newElemsBestSwap)
        count_newElems = count_newElems+newElemsBestSwap
      end if
      
      call memory_deallo(memor_adapt,'Tcav_remeshed','swap_mesh_edges_frozen',Tcav_remeshed)
      call memory_deallo(memor_adapt,'T_bestSwap'   ,'swap_mesh_edges_frozen',   T_bestSwap)
      call deleteNodesCav(nodesCav)
      call delete_boundary_faces(boundary_faces)
      
    end if !isFrozen
    
    call delete_arrayElems(elems_cavity)
    
  end do loopEdges

  if(out_performance) then
    call cpu_time(t1)
    deb_swapLoop = deb_swapLoop + (t1-t0)
  end if
  
  6666 continue
  if( isModifiedMesh ) then
    count_newPoints = 0_ip
    nullify(coords_new_fake)
    call memory_alloca(memor_adapt,'coords_new_fake','swap_mesh_edges_frozen',coords_new_fake,mesh%ndime,2_ip)!mesh%ndime,1_ip)
    call set_mesh_new_modify(mesh,coords_new_fake,lnods_new,count_newPoints,count_newElems,&
      isDeletedElem=isFrozenElem)
    call memory_deallo(memor_adapt,'coords_new_fake','swap_mesh_edges_frozen',coords_new_fake)
  else
    call memory_deallo(memor_adapt,'lnods_new'      ,'set_tagged_edges',lnods_new)
  end if
  
  call quality_deallo(q)

  call memory_deallo(memor_adapt,'isFrozenElem','repair_mesh_edgeLengths_frozen',isFrozenElem)
  call memory_deallo(memor_adapt,'taggedEdges_to_node','mod_adapt',taggedEdges_to_node)
  
  if(out_performance) then
    call cpu_time(t1bis)
    deb_swapRest = deb_swapRest + (t1bis-t0bis) - (t1-t0)
  end if
  
  return 
end subroutine swap_mesh_edges_frozen_allInOne
!
!
!
function isAttemptedSwap(n1,n2) result(isAttempted)
  implicit none
  
  integer(ip), intent(in) :: n1,n2
  logical(lg) :: isAttempted
  
  !isAttempted = .false.
  
  ! IMPLEMENT:
  ! check if edge swap exists
  ! if exists check if any of the surrounding nodes has been modified at some point -> is this easy?
  !    I need a track of the nodes of the mesh to know if at this current iteration the state of the node has changed
  isAttempted = attempt_edgeSwap%isEdge(n1,n2)
  
end function isAttemptedSwap
!
!
!
subroutine setAttemptedSwap(n1,n2)
  implicit none
  integer(ip), intent(in) :: n1,n2

  call attempt_edgeSwap%setEdge(n1,n2)
  
end subroutine setAttemptedSwap
!
!
!
subroutine set_tagged_edges_byQ(tol_q,mesh,node_to_elems,q,numTaggedEdges,taggedEdges_to_node)
  use mod_meshEdges,    only: compute_mesh_edgeQuality
  use mod_meshTopology, only: node_to_elems_type
  implicit none
  
  real(rp),                intent(in) :: tol_q
  type(mesh_type_basic),   intent(in) :: mesh
  type(node_to_elems_type),intent(in) :: node_to_elems
  real(rp), pointer,       intent(in) :: q(:)
  !
  integer(ip), pointer, intent(inout) :: taggedEdges_to_node(:,:)
  integer(ip),          intent(out)   :: numTaggedEdges
  !
  real(rp),    pointer  :: edge_q(:)
  integer(ip), pointer  :: edge_to_node(:,:)
  integer(ip), pointer  :: perm_edges_byQ(:)
  integer(ip)           :: iedge, numEdges, theEdge
  real(rp)              :: q_edge
  !
  call compute_mesh_edgeQuality(mesh,node_to_elems,q,edge_to_node,edge_q)
  call sort_perm_reals_increasing(edge_q,perm_edges_byQ)
  
  numEdges = size(edge_to_node,2)
  loopEdges: do iedge=1,numEdges
    theEdge = perm_edges_byQ(iedge)
    q_edge  = edge_q(theEdge)
    if(q_edge>tol_q) exit loopEdges
  end do loopEdges
  numTaggedEdges = iedge-1
  
  nullify(taggedEdges_to_node)
  call memory_alloca(memor_adapt,'taggedEdges_to_node','mod_adapt',taggedEdges_to_node,2_ip,numTaggedEdges)
  if(numTaggedEdges>0) then
    taggedEdges_to_node(:,:) = edge_to_node( : , perm_edges_byQ(1:numTaggedEdges) )
  end if
  
  ! Deallocate and nullify
  call memory_deallo(memor_adapt,'edge_to_node',  'set_tagged_edges_byQ',edge_to_node)
  call memory_deallo(memor_adapt,'edge_q',  'set_tagged_edges_byQ',edge_q)
  call memory_deallo(memor_adapt,'perm_edges_byQ','set_tagged_edges_byQ',perm_edges_byQ)
  !
end subroutine set_tagged_edges_byQ
!
!
!
!--------- SUBMODULE 3: _adapt_auxMesh ---------------------------------------------------------------------------
! include "_adapt_auxMesh.f90"
!
!
!
subroutine out_with_q(mesh,metric,file_out,iout)
  use mod_out_paraview, only: out_paraview_inp
  use mod_quality,      only: compute_mesh_quality_sizeShape
  
  use def_master, only: kfl_paral
  use mod_strings,               only : integer_to_string

  use def_master, only: kfl_paral
  
  implicit none
  !
  type(mesh_type_basic)  , intent(in):: mesh
  type(mesh_metric_type) , intent(in):: metric
  character(len=*), intent(in) :: file_out
  integer(ip) , intent(inout) :: iout
  !
  character(20) :: iout_char
  real(rp), pointer :: q(:)
  !
  if(out_debug_paraview.and.output_steps_paraview) then
    iout = iout+1
    iout_char=int2str(iout)
    !print*,"iout: ",iout
    call compute_mesh_quality_sizeShape(mesh,metric,q)
    !print*,"kfl_paral: ",kfl_paral,"  -> minval(q): ",minval(q),"- maxval(q): ",maxval(q)
    call out_paraview_inp(mesh,filename=TRIM(file_out)//'_'//integer_to_string(kfl_paral)//'_'//iout_char,nodeField=metric%M(1,1,:),elemField=q)
    
    call memory_deallo(memor_adapt,'QUALITY','compute_mesh_quality_sizeShape',q)
    
  end if
  !
end subroutine out_with_q

  !
  !
  !
subroutine out_firstLast(mesh,metric,file_out)

  use def_master,     only: kfl_paral
  use mod_strings,    only: integer_to_string
  use mod_quality,    only: compute_mesh_quality_sizeShape
  use mod_debugTools, only: out_debug_paraview_firstLast, out_id_firstLast
  use mod_out_paraview, only: out_paraview_inp

  implicit none
  !
  type(mesh_type_basic)  , intent(in):: mesh
  type(mesh_metric_type) , intent(in):: metric
  character(len=*), intent(in) :: file_out
  !
  character(20) :: iout_char
  real(rp), pointer :: q(:)
  real(rp) :: t0,t1
  !    
  if(out_debug_paraview_firstLast) then
    call cpu_time(t0)
    
    !if(kfl_paral.eq.1) then
    out_id_firstLast = out_id_firstLast+1
    !end if
    !call PAR_BARRIER()
  
    iout_char=int2str(out_id_firstLast)
    call compute_mesh_quality_sizeShape(mesh,metric,q)
    call out_paraview_inp(mesh,&
      filename='./'//&!'/Users/abel/Desktop/Work/CasesAlya/paraview/'//&
      TRIM(file_out)//'_FL_'//integer_to_string(kfl_paral)//'_'//iout_char,&
      nodeField=metric%M(1,1,:),elemField=q)

    call memory_deallo(memor_adapt,'QUALITY','compute_mesh_quality_sizeShape',q)
  
    call cpu_time(t1)
    print*,'    Time exporting: ',t1-t0
    
!     if(kfl_paral<0_ip) then !ONLY SEQUENTIAL
!     end if
  end if
  !
end subroutine out_firstLast
!
!
!
function int2str(num) result(str)
  implicit none
  integer(ip), intent(in):: num
  character(20) :: str
  ! convert integer to string using formatted write
  write(str, '(i20)') num
  str = adjustl(str)
end function int2str
!
!
!
subroutine removeDeletedNodes(mesh_new,metric_new,isDeletedPoint,mapNodes_new_to_old)
  use mod_out_paraview
  implicit none
  
  type(mesh_type_basic),  intent(inout) :: mesh_new
  type(mesh_metric_type), intent(inout) :: metric_new
  !logical(lg),            intent(inout) :: isDeletedPoint(:)  
  logical(lg), pointer,   intent(inout) :: isDeletedPoint(:)  
  integer(ip), pointer,   intent(inout) :: mapNodes_new_to_old(:)
  !
  integer(ip) :: map_old_to_new_point(mesh_new%npoin)
  integer(ip) :: ipoin,ielem,numPoints,numDelPoints,numPoints_preAdapt,countPoints
  integer(ip), pointer :: mapNodes_new_to_old_aux(:)
  type(mesh_type_basic)  :: mesh ! to copy the old mesh (change names to improve routine)
  !
  numDelPoints = count(isDeletedPoint)
  numPoints = mesh_new%npoin - numDelPoints

  numPoints_preAdapt = size(isDeletedPoint)
  
  nullify(mapNodes_new_to_old_aux)
  call memory_alloca(memor_adapt,'mapNodes_new_to_old_aux','mod_adapt',mapNodes_new_to_old_aux,numPoints_preAdapt)
  mapNodes_new_to_old_aux(:) = mapNodes_new_to_old
  
  call memory_deallo(memor_adapt,'mapNodes_new_to_old','mod_adapt',mapNodes_new_to_old)
  nullify(mapNodes_new_to_old)
  call memory_alloca(memor_adapt,'mapNodes_new_to_old','mod_adapt',mapNodes_new_to_old,numPoints)
  mapNodes_new_to_old(:) = mark_isNodeNotFromIniMesh !(/ ( 0_ip*iaux, iaux = 1_ip, numPoints) /) ! identity mapping to start
  
  !
  ! *** COPY MESH DELETING REMOVED NODES ****************************************
  call mesh%init(mesh_new%name)
  call mesh%copy(mesh_new)
  call mesh_new%deallo()
  
  call mesh_new%init(mesh%name)
  call mesh_new%alloca(mesh%ndime,mesh%mnode,mesh%nelem,numPoints)
  mesh_new%ltype(:) = mesh%ltype(1)
  
  countPoints = 0_ip
  do ipoin=1,numPoints_preAdapt!size(isDeletedPoint)!mesh%npoin
    if(isDeletedPoint(ipoin)) then
      map_old_to_new_point(ipoin) = mark_isNodeDeleted!-1000_ip
    else
      countPoints = countPoints + 1_ip
      mesh_new%coord(:,countPoints) = mesh%coord(:,ipoin)
      map_old_to_new_point(ipoin) = countPoints
      !
      mapNodes_new_to_old(countPoints) = mapNodes_new_to_old_aux(ipoin)
    end if
  end do
  do ipoin = (numPoints_preAdapt+1),mesh%npoin
    countPoints = countPoints+1
    mesh_new%coord(:,countPoints) = mesh%coord(:,ipoin)
    map_old_to_new_point(ipoin) = countPoints
  end do
  if(countPoints.ne.numPoints) then
    print*,'error in num points'
    print*,countPoints
    print*,numPoints
    call runend('error in removeDeletedNodes')
  end if

  do ielem=1,mesh%nelem
!     do ipoin=1,3
!       if(map_old_to_new_point(mesh%lnods(ipoin,ielem))<0_ip) then
!         print*,mesh%lnods(ipoin,ielem),'   ',map_old_to_new_point(mesh%lnods(ipoin,ielem))
!         print*,mesh%lnods(:,ielem)
!         call out_paraview_inp(mesh,'./unitt/unitt_adapti_hessian_alya/mesh_kk',nodeField=map_old_to_new_point*1.0_rp)
!       end if
!     end do
    mesh_new%lnods(:,ielem) = map_old_to_new_point(mesh%lnods(:,ielem))
  end do
  
  call mesh%deallo()
  
  ! *** DELETING REMOVED NODES OF THE METIRC ****************************************
  call metric_new%remove_nodes(isDeletedPoint)
  
  call memory_deallo(memor_adapt,'mapNodes_new_to_old_aux','mod_adapt',mapNodes_new_to_old_aux)
  
end subroutine removeDeletedNodes
!
!
!
subroutine updateNodeData(mesh_new,metric_new,isDeletedPoint,mapNodes_new_to_old)
  use mod_out_paraview
  implicit none
  
  type(mesh_type_basic),  intent(inout) :: mesh_new
  type(mesh_metric_type), intent(inout) :: metric_new
  !logical(lg),            intent(inout) :: isDeletedPoint(:)  
  logical(lg), pointer,   intent(inout) :: isDeletedPoint(:)  
  integer(ip), pointer,   intent(inout) :: mapNodes_new_to_old(:)
  !
  integer(ip) :: ipoin,ielem,numPoints,numDelPoints,numPoints_preAdapt,countPoints
  
  integer(ip), pointer :: mapNodes_new_to_old_aux(:)
  integer(ip) :: iaux
  
  type(mesh_type_basic)  :: mesh ! to copy the old mesh (change names to improve routine)
  !
  numDelPoints = 0_ip
  numPoints = mesh_new%npoin - numDelPoints
  numPoints_preAdapt = size(isDeletedPoint)
  
  nullify(mapNodes_new_to_old_aux)
  call memory_alloca(memor_adapt,'mapNodes_new_to_old_aux','mod_adapt',mapNodes_new_to_old_aux,numPoints_preAdapt)
  mapNodes_new_to_old_aux(:) = mapNodes_new_to_old
  
  call memory_deallo(memor_adapt,'mapNodes_new_to_old','mod_adapt',mapNodes_new_to_old)
  nullify(mapNodes_new_to_old)
  call memory_alloca(memor_adapt,'mapNodes_new_to_old','mod_adapt',mapNodes_new_to_old,numPoints)
  mapNodes_new_to_old(:) = mark_isNodeNotFromIniMesh !(/ ( 0_ip*iaux, iaux = 1_ip, numPoints) /) ! identity mapping to start
  
  !
  ! We just need to increase size of the node mapping..
  ! ... and we mark deleted nodes...
  !
  mapNodes_new_to_old(1:numPoints_preAdapt) = mapNodes_new_to_old_aux
  
  do ipoin=1,numPoints_preAdapt!size(isDeletedPoint)!mesh%npoin
    if(isDeletedPoint(ipoin)) then
      mapNodes_new_to_old(ipoin) = mark_isNodeDeleted !-1000_ip ! old node deleted
    end if
  end do
  
  call memory_deallo(memor_adapt,'mapNodes_new_to_old_aux','mod_adapt',mapNodes_new_to_old_aux)
  
end subroutine updateNodeData
! 
!
!
function generate_mid_point(x1,x2,M1,M2) result(x_new)
  implicit none
  real(rp), intent(in) :: x1(:), x2(:)
  real(rp), optional, intent(in) :: M1(:,:), M2(:,:)
  real(rp) :: x_new(size(x1))
  
  !real(rp) :: alpha
  
  call runend('should not be used any more.. moved to metric')
  
  if(present(M1)) then
    !print*,"implement generate_mid_point with metric"
    !call runend("epa")
    !x_new = (x1+x2)/2.0
    
   ! alpha = 
    
  else
    x_new = (x1+x2)/2.0
  end if
  
end function generate_mid_point
!
!
!
subroutine set_mesh_new_modify(&
  mesh_new,coord_new,lnods_new,num_new_points,num_new_elems,isDeletedElem,&
  q,q_new,isNewElem)!,isDeletedPoint)
  !
  use mod_quality,    only: quality_deallo
  use mod_debugTools, only: deb_setNewMeshModify
  use def_master,    only: kfl_paral
  
  implicit none
  !
  type(mesh_type_basic),           intent(inout) :: mesh_new
  real(rp),              pointer,  intent(inout) :: coord_new(:,:)
  integer(ip),           pointer,  intent(inout) :: lnods_new(:,:)
  integer(ip),                     intent(in)    :: num_new_points
  integer(ip),                     intent(in)    :: num_new_elems
  logical(lg),           pointer,  intent(in)    :: isDeletedElem(:)
  real(rp),    optional, pointer,  intent(inout) :: q(:)
  real(rp),    optional, pointer,  intent(inout) :: q_new(:)
  logical(lg), optional, pointer,  intent(inout) :: isNewElem(:)
!   logical(lg), optional,           intent(in)    :: isDeletedPoint(:)
  
  integer(ip) :: numElems, numDelElems, numPoints, ielem, ipoin
  integer(ip) :: numDelPoints
  integer(ip) :: ini_elems, fin_elems, ini_points, fin_points
  
!   integer(ip) :: map_old_to_new_point(mesh_new%npoin)
  type(mesh_type_basic)                :: mesh ! to copy the old mesh (change names to improve routine)
  !
  real(rp)    :: q_copy(mesh_new%nelem)
  logical(lg) :: isNewElem_copy(mesh_new%nelem)
  !
  real(rp) :: t0,t1
  !
  integer(ip) :: ii
  !
  if(out_performance) call cpu_time(t0)
  
!   if(mesh%nelem==0_ip) then
!     return
!   end if
  
!   if(present(isDeletedPoint)) then
!     !print*,'Implement removing points from mesh%coord'
!     numDelPoints = count(isDeletedPoint)
!     call runend('do not use inside the frozen loop, since tagged edge nodes are not updated')
!   else
!     !print*,'Only points are added, not removed'
!     numDelPoints = 0_ip
!   end if
  numDelPoints = 0_ip
  
  call mesh%init(mesh_new%name)
  call mesh%copy(mesh_new)
  numDelElems = count(isDeletedElem)
  numElems  = mesh%nelem - numDelElems  + num_new_elems
  numPoints = mesh%npoin - numDelPoints + num_new_points
  call mesh_new%deallo()
  
  call mesh_new%init(mesh%name)
  call mesh_new%alloca(mesh%ndime,mesh%mnode,numElems,numPoints)
  
  if(present(q)) then
    do ii=1,size(q)
      q_copy(ii) = q(ii)!q_copy = q
    end do
    call quality_deallo(q)
    nullify(q)
    call memory_alloca(memor_adapt,'q','compute_mesh_quality',q,numElems)
  end if
  if(present(isNewElem)) then
    if(.not.associated(isNewElem)) call runend('isNewElem should be associated..')
    do ii=1,size(isNewElem)
      isNewElem_copy(ii) = isNewElem(ii) !isNewElem_copy = isNewElem
    end do
    call memory_deallo(memor_adapt,'isNewElem','mod_adapt',isNewElem)
    nullify(isNewElem)
    call memory_alloca(memor_adapt,'isNewElem','mod_adapt',isNewElem,numElems)
  end if

!   if(present(isDeletedPoint)) then
!     numPoints = 0_ip
!     do ipoin=1,size(isDeletedPoint)!mesh%npoin
!       if(.not.isDeletedPoint(ipoin)) then
!         numPoints = numPoints + 1_ip
!         mesh_new%coord(:,numPoints) = mesh%coord(:,ipoin)
!         map_old_to_new_point(ipoin) = numPoints
!       end if
!     end do
!   else
!     mesh_new%coord(:,1:mesh%npoin) = mesh%coord
!     numPoints = mesh%npoin
!   end if
  !
  do ii=1,mesh%npoin
    mesh_new%coord(:,ii) = mesh%coord(:,ii)
  end do
  numPoints = mesh%npoin
  
  if( num_new_points > 0 ) then
    ini_points = numPoints+1_ip
    fin_points = numPoints+num_new_points
    !mesh_new%coord(:,ini_points:fin_points) = coord_new(:,1:num_new_points)
    do ii=1,num_new_points
      mesh_new%coord(:,ini_points+ii-1) = coord_new(:,ii)
    end do
  end if
  
  !mesh_new%ltype(:) = mesh%ltype(1)
  numElems = 0_ip
  do ielem=1,mesh%nelem
    if(.not.isDeletedElem(ielem)) then
      numElems = numElems + 1_ip
      mesh_new%ltype(numElems) = mesh%ltype(1)
!       if(present(isDeletedPoint)) then
!         mesh_new%lnods(:,numElems) = map_old_to_new_point(mesh%lnods(:,ielem))
!       else
!         mesh_new%lnods(:,numElems) = mesh%lnods(:,ielem)
!       end if
      mesh_new%lnods(:,numElems) = mesh%lnods(:,ielem)
      if(present(q)        )         q(numElems) =          q_copy(ielem)
      if(present(isNewElem)) isNewElem(numElems)  = isNewElem_copy(ielem)
    end if
  end do
  ini_elems = numElems+1_ip
  fin_elems = numElems+num_new_elems
  !mesh_new%lnods(:,ini_elems:fin_elems) = lnods_new(:,1:num_new_elems)
  do ii=1,num_new_elems
    numElems = numElems+1_ip!ini_elems+ii-1
    mesh_new%lnods(:,numElems) = lnods_new(:,ii)
    mesh_new%ltype(  numElems) = mesh%ltype(1)
    if(present(q)) q(numElems) = q_new(ii)
    if(present(isNewElem)) isNewElem(numElems) = .true.
  end do
  
!   if(present(q)) then
!     !q(ini_elems:fin_elems) = q_new(1:num_new_elems)
!     do ii=1,num_new_elems
!       q(ini_elems+ii-1) = q_new(ii)
!     end do
!   end if
!
!   if(present(isNewElem)) then
!     !!isNewElem(1:) = isNewElem_copy
!     !!isNewElem = .false.
!
!     !isNewElem(ini_elems:fin_elems) = .true.
!     do ii=ini_elems,fin_elems
!       isNewElem(ii) = .true.
!     end do
!   end if
  
  call memory_deallo(memor_adapt,'coord_new','set_mesh_new_modify',coord_new)
  call memory_deallo(memor_adapt,'lnods_new','set_mesh_new_modify',lnods_new) 
  call mesh%deallo()
  if(present(q)) call memory_deallo(memor_adapt,'q_new','set_mesh_new_modify',q_new)

  if(out_performance) then
    call cpu_time(t1)
    deb_setNewMeshModify = deb_setNewMeshModify + (t1-t0)
  end if
  
  return
end subroutine set_mesh_new_modify
!
!
!
subroutine sort_perm_reals_decreasing(f,sort_perm)
  use mod_maths_sort, only: maths_heap_sort !maths_quick_sort
  implicit none
  !
  real(rp)   , pointer, intent(in)    :: f(:)
  integer(ip), pointer, intent(inout) :: sort_perm(:)
  !
  real(rp), pointer :: f_toSort(:)
  integer(ip) :: n_f
  integer(ip) ::task_sort
  integer(ip) :: iaux
  !
  n_f = size(f)
  nullify(sort_perm)
  call memory_alloca(memor_adapt,'sort_perm','mod_adapt',sort_perm,n_f)
  nullify(f_toSort)
  call memory_alloca(memor_adapt,'f_toSort' ,'sort_perm_reals_decreasing',f_toSort ,n_f)
  do iaux = 1_ip,n_f
    f_toSort (iaux) = f(iaux)
    sort_perm(iaux) = iaux
  end do
  !
  task_sort = 1_ip  ! 1->Decreasing value, 2->Increasing value
  call maths_heap_sort(itask=task_sort, nrows=n_f, rvin=f_toSort, ivo1=sort_perm)
  !
  call memory_deallo(memor_adapt,'f_toSort' ,'sort_perm_reals_decreasing',f_toSort)
  !
!   f_toSort = (/1,2,3,4,5/)
!   sort_perm = (/ (iaux, iaux = 1, 5) /)
!   call maths_heap_sort(itask=task_sort, nrows=n_f, rvin=f_toSort, ivo1=sort_perm)!, ELIMINATE_REPLICATES=.false.)
!   print*,'f_toSort: ',f_toSort
!   print*,'sort_perm:',sort_perm
!   print*,'---'
!
!   f_toSort = (/1,3,2,4,5/)
!   sort_perm = (/ (iaux, iaux = 1, 5) /)
!   call maths_heap_sort(itask=task_sort, nrows=n_f, rvin=f_toSort, ivo1=sort_perm)!, ELIMINATE_REPLICATES=.false.)
!   print*,'f_toSort: ',f_toSort
!   print*,'sort_perm:',sort_perm
!   print*,'---'
!
!   f_toSort = (/5,4,3,2,1/)
!   sort_perm = (/ (iaux, iaux = 1, 5) /)
!   call maths_heap_sort(itask=task_sort, nrows=n_f, rvin=f_toSort, ivo1=sort_perm)!, ELIMINATE_REPLICATES=.false.)
!   print*,'f_toSort: ',f_toSort
!   print*,'sort_perm:',sort_perm
!   print*,'---'
  
end subroutine sort_perm_reals_decreasing
!
!
!
subroutine sort_perm_reals_increasing(f,sort_perm)
  use mod_maths_sort, only: maths_heap_sort !maths_quick_sort
  implicit none
  !
  real(rp)   , pointer, intent(in)    :: f(:)
  integer(ip), pointer, intent(inout) :: sort_perm(:)
  !
  real(rp), pointer :: f_toSort(:)
  integer(ip) :: n_f
  integer(ip) ::task_sort
  integer(ip) :: iaux
  !
  n_f = size(f)
  nullify(sort_perm)
  call memory_alloca(memor_adapt,'sort_perm','mod_adapt',sort_perm,n_f)
  nullify(f_toSort)
  call memory_alloca(memor_adapt,'f_toSort' ,'sort_perm_reals_increasing',f_toSort ,n_f)
  do iaux = 1_ip,n_f
    f_toSort (iaux) = f(iaux)
    sort_perm(iaux) = iaux
  end do
  !
  task_sort = 2_ip  ! 1->Decreasing value, 2->Increasing value
  call maths_heap_sort(itask=task_sort, nrows=n_f, rvin=f_toSort, ivo1=sort_perm)
  !
  call memory_deallo(memor_adapt,'f_toSort' ,'sort_perm_reals_increasing',f_toSort)
  
end subroutine sort_perm_reals_increasing
!
!
!
!
!
!
END MODULE mod_adapt

!> @}

!
!
!
! function set_mesh_new_deprecated(&
!   mesh,coord_new,lnods_new,num_new_points,num_new_elems,isDeletedElem,isDeletedPoint) result(mesh_new)
!   type(mesh_type_basic)                :: mesh_new
!   type(mesh_type_basic), intent(inout) :: mesh
!   real(rp),    pointer,  intent(inout) :: coord_new(:,:)
!   integer(ip), pointer,  intent(inout) :: lnods_new(:,:)
!   integer(ip),           intent(in)    :: num_new_points
!   integer(ip),           intent(in)    :: num_new_elems
!   logical(lg),           intent(in)    :: isDeletedElem(:)
!   logical(lg), optional, intent(in)    :: isDeletedPoint(:)
!
!   integer(ip) :: numElems, numDelElems, numPoints, ielem, numDelPoints, ipoin
!   integer(ip) :: ini_elems, fin_elems, ini_points, fin_points
!
!   integer(ip) :: map_old_to_new_point(mesh%npoin)
!
!   if(present(isDeletedPoint)) then
!     !print*,'Implement removing points from mesh%coord'
!     numDelPoints = count(isDeletedPoint)
!     call runend('do not use inside the frozen loop, since tagged edge nodes are not updated')
!   else
!     !print*,'Only points are added, not removed'
!     numDelPoints = 0_ip
!   end if
!
!   numDelElems = count(isDeletedElem)
!   numElems  = mesh%nelem - numDelElems  + num_new_elems
!   numPoints = mesh%npoin - numDelPoints + num_new_points
!
!   call mesh_new%init(mesh%name)
!   call mesh_new%alloca(mesh%ndime,mesh%mnode,numElems,numPoints)
!
!   if(present(isDeletedPoint)) then
!     numPoints = 0_ip
!     do ipoin=1,mesh%npoin
!       if(.not.isDeletedPoint(ipoin)) then
!         numPoints = numPoints + 1_ip
!         mesh_new%coord(:,numPoints) = mesh%coord(:,ipoin)
!         map_old_to_new_point(ipoin) = numPoints
!       end if
!     end do
!   else
!     mesh_new%coord(:,1:mesh%npoin) = mesh%coord
!     numPoints = mesh%npoin
!   end if
!   ini_points = numPoints+1_ip
!   fin_points = numPoints+num_new_points
!   mesh_new%coord(:,ini_points:fin_points) = coord_new(:,1:num_new_points)
!
!   mesh_new%ltype(:) = mesh%ltype(1)
!
!   numElems = 0_ip
!   do ielem=1,mesh%nelem
!     if(.not.isDeletedElem(ielem)) then
!       numElems = numElems + 1_ip
!       if(present(isDeletedPoint)) then
!         mesh_new%lnods(:,numElems) = map_old_to_new_point(mesh%lnods(:,ielem))
!       else
!         mesh_new%lnods(:,numElems) = mesh%lnods(:,ielem)
!       end if
!     end if
!   end do
!   ini_elems = numElems+1_ip
!   fin_elems = numElems+num_new_elems
!   mesh_new%lnods(:,ini_elems:fin_elems) = lnods_new(:,1:num_new_elems)
!
!   call memory_deallo(memor_adapt,'coord_new','set_tagged_edges',coord_new)
!   call memory_deallo(memor_adapt,'lnods_new','set_tagged_edges',lnods_new)
!
!   call mesh%deallo()
!   return
! end function set_mesh_new_deprecated
