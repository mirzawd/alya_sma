!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Adaptivity
!> @{
!> @file    mod_meshSmoothing.f90
!> @author  abel.gargallo
!> @date    2021-03-31
!> @brief   mod_metric
!> @details mod_metric
!>
!>          To add further details
!>
!>
!-----------------------------------------------------------------------

MODULE mod_smoothing
!***************************************************************
!*
!*  Module for metric computation and operations
!*
!***************************************************************
use def_kintyp_basic, only : ip,rp,lg

use mod_memory,             only: memory_alloca, memory_deallo
use def_adapt,              only: memor_adapt

use mod_debugTools,         only: out_debug_text, out_debug_paraview

implicit none

private

public :: smooth_mesh


real(rp) , parameter      :: tol_x = 0.05_rp
real(rp) , parameter      :: tol_f = 0.05_rp

logical(lg), parameter    :: doFunOptization = .true.
logical(lg), parameter    :: doLapSmoothing  = .false.




logical(lg), parameter    :: out_time = .false.

!
!
!
CONTAINS
!
!
!
function my_norm2(v) result(v_norm2)
  implicit none
  
  real(rp), intent(in) :: v(:)
  real(rp) :: v_norm2
  
  v_norm2 = sqrt(sum(v**2_ip))
  
  return 
end function my_norm2
!
!
!
subroutine smooth_mesh(mesh,metric,min_q_smooth)
  use def_kintyp_mesh_basic,  only : mesh_type_basic
  use mod_metric,             only : mesh_metric_type
  
  use mod_meshTopology,       only : node_to_elems_type
  use mod_meshTopology,       only : node_to_nodes_type
  use mod_meshTopology,       only : getNodesCav
  use mod_meshTopology,       only : getIsBoundaryNode, deleteIsBouNode
  
  use mod_quality,            only : compute_mesh_quality_sizeShape
  use mod_quality,            only : compute_distortion_sizeshape_subset
  use mod_quality,            only : compute_minQuality_sizeshape_subset
  use mod_quality,            only : compute_quality_sizeshape_subset
  use mod_quality,            only : quality_deallo

  use mod_objectiveFun,       only : set_fun_meshData, delete_fun_meshData
  use mod_objectiveFun,       only : set_fun_nodeData, delete_fun_nodeData
  use mod_objectiveFun,       only : f_node_distortion
  
  use mod_optimization,       only : optimizeFunction, NEWTON_RAPHSON, STEEPEST_DESCENT
  use mod_messages,           only : messages_live

  implicit none
  !
  type(mesh_type_basic) ,   intent(inout) :: mesh
  type(mesh_metric_type),   intent(inout) :: metric
  real(rp),                 intent(in)    :: min_q_smooth
  !
  type(node_to_elems_type)  :: node_to_elems
  logical(lg),  pointer     :: isBoundaryNode(:)
  real(rp),     pointer     :: q(:)
  logical(lg)               :: isNotConverged
  real(rp)                  :: minQ_node
  integer(ip)               :: inode
  integer(ip), pointer      :: elems_adjToNode(:)
  integer(ip)               :: numAdjNodes
  integer(ip), pointer      :: nodes_adjToNode(:)
  real(rp)                  :: x_cm(mesh%ndime),x_prev(mesh%ndime),x_opt(mesh%ndime)
  real(rp)                  :: distortion_prev, distortion_cm
  real(rp)                  :: eps_dim, L_char
  real(rp)                  :: v_char(mesh%ndime)
  !real(rp)                  :: distortion_node
  real(rp)                  :: f0,f1
  real(rp)                  :: e_rel_x, e_rel_f
  integer(ip)               :: iter
  integer(ip), parameter    :: max_iter = 100_ip
  logical(lg)               :: isConvergedNode(mesh%npoin)
  !integer(ip), allocatable  :: targetNodes(:)
  integer(ip), pointer  :: targetNodes(:)
  integer(ip)               :: numTargetNodes,i_targetNode
  real(rp) :: t0,t1
  
  x_prev(:) = 0.0_rp

  eps_dim = 0.1_rp**(14_rp/mesh%ndime)
  !
  if(out_time) call cpu_time(t0)
  call node_to_elems%set(mesh)
  call set_fun_meshData(mesh,metric)
  if(out_debug_text.and.out_time) call cpu_time(t1)
  if(out_debug_text.and.out_time) print '("   Time n-t-e + fundata = ",f6.0," sec")',(t1-t0)

  if(out_debug_text.and.out_time) call cpu_time(t0)
  call compute_mesh_quality_sizeShape(mesh,metric,q)

  call getIsBoundaryNode(mesh,isBoundaryNode)
  if(out_debug_text.and.out_time) call cpu_time(t1)
  if(out_debug_text.and.out_time) print '("   Time q + bound = ",f6.0," sec")',(t1-t0)

  isConvergedNode = isBoundaryNode
  !print*,'listNodes = (/(inode, inode=1,mesh%npoin)/)'
!   listNodes = (/(inode, inode=1,mesh%npoin)/)
  !print*,'targetNodes     = pack(listNodes,mask=.not.isConvergedNode)'
!   targetNodes     = pack(listNodes,mask=.not.isConvergedNode)
!   !targetNodes     = pack(listNodes,mask=.not.isConvergedNode(listNodes))
!   print*,'numTargetNodes  = size(targetNodes)'
!   numTargetNodes  = size(targetNodes)
  ! canviar aixo per comptar not is converged i allocatar target
  numTargetNodes = count(.not.isConvergedNode)
  nullify(targetNodes)
  call memory_alloca(memor_adapt,'targetNodes','smooth_mesh',targetNodes,numTargetNodes)
  numTargetNodes = 0_ip
  do inode=1_ip,mesh%npoin
    if(.not.isConvergedNode(inode)) then
      numTargetNodes=numTargetNodes+1_ip
      targetNodes(numTargetNodes)=inode
    end if
  end do 
  
  isNotConverged     = numTargetNodes>0_ip
  iter = 0_ip
  if(out_debug_text.and.out_time) call cpu_time(t0)
  do while(isNotConverged.and.(iter<max_iter)) 
    iter = iter+1_ip
    
    do i_targetNode=1,numTargetNodes
      inode = targetNodes(i_targetNode)
      
      if( .not.isConvergedNode(inode) ) then !.and.(.not.isBoundaryNode(inode))

        if(out_debug_text.and.out_time) call cpu_time(t0)
        call node_to_elems%get_elems(inode,elems_adjToNode)
        if(associated(elems_adjToNode)) then
          minQ_node = minval(  q(elems_adjToNode)  )
        else
          isConvergedNode(inode) = .true.
          minQ_node = 1.0_rp
        end if
        
        if( minQ_node < min_q_smooth ) then
          call set_fun_nodeData(inode,elems_adjToNode)
          call getNodesCav(mesh%lnods(:,elems_adjToNode),numAdjNodes,nodes_adjToNode) ! create structure node_to_node
          
          x_prev = mesh%coord(:,inode)
          
          v_char = mesh%coord(:,nodes_adjToNode(1))-mesh%coord(:,inode)
          L_char = my_norm2(v_char) !norm2(v_char) does not compile in p9 ...
          v_char = mesh%coord(:,nodes_adjToNode(numAdjNodes))-mesh%coord(:,inode)
          L_char = L_char + my_norm2(v_char)
          
          if(doLapSmoothing) then
            if(out_debug_text.and.out_time) call cpu_time(t0)
            distortion_prev = f_node_distortion( mesh%ndime ,x_prev)

            x_cm = sum(mesh%coord(:,nodes_adjToNode),2)/numAdjNodes
            distortion_cm = f_node_distortion( mesh%ndime ,x_cm)

            if(distortion_cm<distortion_prev-eps_dim) then ! undo
              mesh%coord(:,inode) = x_cm
              q(elems_adjToNode) = compute_quality_sizeshape_subset(mesh,metric,elems_adjToNode)
            
              e_rel_f = abs((distortion_cm-distortion_prev)/distortion_prev)
              e_rel_x = my_norm2(x_prev-x_cm)/L_char
              if( (e_rel_x>tol_x).or.(e_rel_f>tol_f) ) then
                isConvergedNode(inode) = .false. ! for efficicency (maybe should put it as converged if some iterations have been conv )
              end if
            else
              mesh%coord(:,inode) = x_prev
            end if
            if(out_debug_text.and.out_time) call cpu_time(t1)
            if(out_debug_text.and.out_time) print '("   Time FIRST = ",f6.0," sec")',(t1-t0)
          end if
          
          if(doFunOptization) then

            minQ_node = compute_minQuality_sizeshape_subset(mesh,metric,elems_adjToNode)
            if( minQ_node < min_q_smooth ) then
              !
              if(out_debug_text.and.out_time) call cpu_time(t0)
              !
              !x_opt = optimizeFunction(f_node_distortion,mesh%ndime,mesh%coord(:,inode),f0,f1,NEWTON_RAPHSON)
              !x_opt = optimizeFunction(f_node_distortion,mesh%ndime,mesh%coord(:,inode),f0,f1,STEPEST_DESCENT)
              x_opt = optimizeFunction(f_node_distortion,mesh%ndime,mesh%coord(:,inode),f0,f1) ! default method in mod_optimization
              !
              if(out_debug_text.and.out_time) call cpu_time(t1)
              if(out_debug_text.and.out_time) print '("   Time optimizeFunction = ",f6.0," sec")',(t1-t0)
              if(out_debug_text.and.out_time) call cpu_time(t0)
              !
              if(f1<f0) then
                mesh%coord(:,inode) = x_opt
                q(elems_adjToNode) = compute_quality_sizeshape_subset(mesh,metric,elems_adjToNode)
                e_rel_f = abs((f0-f1)/f1)               
                e_rel_x = my_norm2(x_prev-x_opt)/L_char
                if( (e_rel_x>tol_x).or.(e_rel_f>tol_f) ) then
                  isConvergedNode(inode) = .false.
                else
                  isConvergedNode(inode) = .true.
                end if
              else
                isConvergedNode(inode) = .true.
              end if
              if(out_debug_text.and.out_time) call cpu_time(t1)
              if(out_debug_text.and.out_time) print '("   Time rest = ",f6.0," sec")',(t1-t0)
            end if
          end if
          
          call delete_fun_nodeData()
        else
          isConvergedNode(inode) = .true.
        end if
        
        call memory_deallo(memor_adapt,'elems_adjToNode','smooth_mesh',elems_adjToNode)
        
      end if
      
    end do
    
    !deallocate(targetNodes)
    !targetNodes     = pack(listNodes,mask=.not.isConvergedNode)
    !numTargetNodes  = size(targetNodes)
    call memory_deallo(memor_adapt,'targetNodes','smooth_mesh',targetNodes)
    numTargetNodes = count(.not.isConvergedNode)
    nullify(targetNodes)
    call memory_alloca(memor_adapt,'targetNodes','smooth_mesh',targetNodes,numTargetNodes)
    numTargetNodes = 0_ip
    do inode=1_ip,mesh%npoin
      if(.not.isConvergedNode(inode)) then
        numTargetNodes=numTargetNodes+1_ip
        targetNodes(numTargetNodes)=inode
      end if
    end do 
    isNotConverged     = numTargetNodes>0_ip
  end do
  
  !call memory_deallo(memor_adapt,'q','smooth_mesh',q)
  call quality_deallo(q)
  !call memory_deallo(memor_adapt,'isBoundaryNode','smooth_mesh',isBoundaryNode)
  call deleteIsBouNode(isBoundaryNode)
  call memory_deallo(memor_adapt,'targetNodes','smooth_mesh',targetNodes)
  call delete_fun_meshData()
  
end subroutine smooth_mesh
!
!
!
! subroutine smoothNode()
!
! end subroutine smoothNode
!
!
!
! fer un mod_objectiveFunction_mesh
!    guardar punter a malla
!    guardar node inode que es toca
!     per tant fer funcio objectiu nomes amb x (de linode per la malla apuntada)
!    fer funcio set_mesh
!    fer funcio set_node_id
!    fer funcio delete_mesh
!    fer funcio delete_node_id
!
!
! subroutine f()
! end subroutine f
!
!
!
END MODULE mod_smoothing

!> @}
