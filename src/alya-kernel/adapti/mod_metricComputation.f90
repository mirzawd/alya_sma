!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Adaptivity
!> @{
!> @file    mod_metricComputation.f90
!> @author  abel.gargallo
!> @date    2021-03-31
!> @brief   mod_metricComputation to compute metric from solution following the ideas from Alauzet&Loiselle
!> @details mod_metricComputation:
!>
!>          Note: this module needs the master (for the mesh topology to set metric complexity and gradation)
!>
!>
!-----------------------------------------------------------------------

MODULE mod_metricComputation
!***************************************************************
!*
!*  Module for metric computation and operations
!*
!***************************************************************
use def_kintyp_basic, only : ip,rp,lg
use mod_metric,       only : mesh_metric_type

implicit none

private

interface compute_sizeField_mesh
   module procedure &
        compute_sizeField_meshTypeBasic,&
        compute_sizeField_meshType
end interface compute_sizeField_mesh



public :: compute_metric_from_sol_viaHessian
public :: compute_metric_from_hessian
public :: compute_sizeField_from_sol_viaHessian
public :: compute_sizeField_from_multiSol_viaHessian
public :: compute_sizeField_mesh
public :: compute_smoothSizeField
public :: compute_sizeField_smoothedBoundary

logical(lg), parameter :: do_gradate_metrics = .true.!.false.!

real(rp),    parameter :: hmin_sec = TINY(1.0_rp)**0.33! to allow it squaring


!
!
!
CONTAINS
!
!
!
subroutine compute_sizeField_smoothedBoundary(mesh,hmesh,htarget) !and optional mesh size field
  use def_kintyp_mesh_basic, only: mesh_type_basic
  use mod_metric, only: compute_metric_fromNodalSizeField
  !
  type(mesh_type_basic), intent(in)     :: mesh
  real(rp), pointer,     intent(in)     :: hmesh(:)
  real(rp), pointer,     intent(inout)  :: htarget(:)
  
  type(mesh_metric_type) :: metric
  type(mesh_type_basic)  :: mesh_bou
  
  !
  if(mesh%npoin==0_ip) then
    return
  end if
  
  call mesh_bou % init('BOUNDARY')
  call mesh_bou % boundary_mesh(mesh)
  
  htarget(mesh_bou%permn) = hmesh(mesh_bou%permn)
  !htarget(mesh_bou%permn) = ( hmesh(mesh_bou%permn) + htarget(mesh_bou%permn)) /2.0
  !htarget(mesh_bou%permn) = sqrt( hmesh(mesh_bou%permn) * htarget(mesh_bou%permn))
  
  call compute_metric_fromNodalSizeField(htarget,mesh%ndime,metric)
  
  call compute_size_gradation(metric)
  
  htarget = metric%M(1,1,:)
  
  htarget(mesh_bou%permn) = hmesh(mesh_bou%permn) ! enforce boundary size
  
  call mesh_bou%deallo()
  call metric%deallo()
  
  !
end subroutine compute_sizeField_smoothedBoundary
!
!
!
subroutine compute_sizeField_meshTypeBasic(hh,mesh)
  !
  use def_kintyp_mesh_basic, only: mesh_type_basic
  use mod_meshTopology, only: compute_mesh_edges
  !
  type(mesh_type_basic), intent(in) :: mesh
  real(rp), pointer, intent(inout) :: hh(:)
  !
  integer(ip), pointer :: edge_to_node(:,:)
  real(rp),    pointer :: coord(:,:)
  !
  if(mesh%npoin==0) then
    return
  end if
  
  call compute_mesh_edges(mesh,edge_to_node)
  
  nullify(coord)
  coord        => mesh % coord(:,1:mesh%npoin)

  call compute_sizeField_meshData(hh,coord,edge_to_node)
  
end subroutine compute_sizeField_meshTypeBasic
!
!
!
subroutine compute_sizeField_meshType(hh,mesh)

  use def_kintyp_mesh, only: mesh_type
  use mod_graphs, only: graphs_edges
  use def_adapt,              only: memor_adapt
  
  type(mesh_type), intent(inout) :: mesh
  !
  real(rp), pointer, intent(inout) :: hh(:)
  
  integer(ip), pointer :: edge_to_node(:,:)
  real(rp),    pointer :: coord(:,:)
  !
  if(mesh % npoin==0_ip) then
    nullify(hh)
    return
  end if

  nullify( mesh % edge_to_node )
  nullify( mesh % ledgs        )
  nullify( mesh % lnned        )
  nullify( mesh % ledgb        )
  nullify( mesh % lnneb        )
  nullify( mesh % r_edg        )
  nullify( mesh % c_edg        )
  
  call graphs_edges(mesh,memor_adapt)

  if(mesh % nedge==0_ip) then
    nullify(hh)
    return
  end if

  nullify(edge_to_node)
  edge_to_node => mesh % edge_to_node
  nullify(coord)
  coord        => mesh % coord(:,1:mesh%npoin)
  
  call compute_sizeField_meshData(hh,coord,edge_to_node)
    
end subroutine compute_sizeField_meshType
!
!
!
subroutine compute_sizeField_meshData(hh,coord,edge_to_node)
  !
  use mod_memory,             only: memory_alloca, memory_deallo
  use def_adapt,              only: memor_adapt
  use mod_metric, only: interpolation_law, law_geometric_h, law_linear_h
  
!  use def_master,             only: kfl_paral
  !
  real(rp),    pointer,  intent(in) :: coord(:,:)
  integer(ip), pointer,  intent(in) :: edge_to_node(:,:)
  
  real(rp),    pointer,  intent(inout):: hh(:)
  !
  integer(ip) :: node_adj(size(coord,2))
  
  integer(ip) :: iedge, n1, n2, npoin
  integer(ip) :: nedge
  real(rp)    :: len_eucl
  !
  logical(lg), parameter :: do_minSize = .false.
  !
  npoin = size(coord,2)
  
  !print*,kfl_paral,' -> size coord: ',size(coord,1),"  ",size(coord,2),' -> size edge_to_node: ',size(edge_to_node,1),"  ",size(edge_to_node,2)
  
  nullify(hh)
  if(npoin>0_ip) then
    call memory_alloca(memor_adapt,'hh','compute_sizeField_from_sol_viaHessian',hh,npoin)
  else
    return
  end if
  
  if(.not.do_minSize) then
    !
    ! ****Compute node size as the interpolated size of the adjacent edges
    !
    node_adj = 0_ip
    if(interpolation_law.eq.law_geometric_h) then
      hh       = 1.0_rp
    elseif(interpolation_law.eq.law_linear_h) then
      hh       = 0.0_rp
    else
      call runend('implement law..')
    end if
  
    nedge = size(edge_to_node,2)
    do iedge=1,nedge
  
      !print*,kfl_paral,"  ",iedge,size(edge_to_node)
      n1 = edge_to_node(1,iedge)
      n2 = edge_to_node(2,iedge)
    
      !len_eucl = norm2( coord(:,n1) - coord(:,n2) )
      len_eucl = sqrt( sum( (coord(:,n1) - coord(:,n2))**2 ) )
    
      node_adj(n1) = node_adj(n1)+1_ip
      node_adj(n2) = node_adj(n2)+1_ip
    
      if(interpolation_law.eq.law_geometric_h) then
        hh(n1) = hh(n1)*len_eucl
        hh(n2) = hh(n2)*len_eucl
      elseif(interpolation_law.eq.law_linear_h) then
        hh(n1) = hh(n1) + len_eucl
        hh(n2) = hh(n2) + len_eucl
      end if
    
    end do
  
    if(interpolation_law.eq.law_geometric_h) then
      hh = hh**(1.0_rp/node_adj)
    elseif(interpolation_law.eq.law_linear_h) then
      hh = hh/node_adj
    end if
  
  else
    !
    ! ****Compute node size as the minimum of the adjacent edges
    !
    hh = huge(0.0_rp)
    nedge = size(edge_to_node,2)
    do iedge=1,nedge
      n1 = edge_to_node(1,iedge)
      n2 = edge_to_node(2,iedge)
      len_eucl = sqrt( sum( (coord(:,n1) - coord(:,n2))**2 ) )
      hh(n1) = min(hh(n1),len_eucl)
      hh(n2) = min(hh(n2),len_eucl)
    end do
    !
  end if

    
end subroutine compute_sizeField_meshData
!
!
!
subroutine compute_smoothSizeField(hh,hh_in,ndime,maxMeshSize,numTargetNodes)
  use mod_metric,             only: compute_metric_fromNodalSizeField
  use def_adapt,              only: memor_adapt
  use mod_memory,             only: memory_alloca, memory_deallo

  real(rp), pointer, intent(inout) :: hh(:)
  real(rp), pointer, intent(in)    :: hh_in(:)
  integer(ip),       intent(in)    :: ndime
  real(rp)   ,       intent(in)    :: maxMeshSize
  integer(ip), optional, intent(in):: numTargetNodes

  type(mesh_metric_type) :: metric
  logical(lg) :: isIsotropic
  
  isIsotropic = .true.
  
  call compute_metric_fromNodalSizeField(hh_in,ndime,metric)
  
  call compute_size_gradation(metric)  ! gradate should be an atribute of metric metric%gradate(optional type?, optional parameter?)
  
  if(present(numTargetNodes)) then
    if(numTargetNodes>0_ip) then
      call metric%setMetricComplexity(numTargetNodes)
    end if
  end if

  nullify(hh)
  if(associated(hh_in)) then
    call memory_alloca(memor_adapt,'hh','compute_sizeField_from_sol_viaHessian',hh,size(hh_in))
    hh(:) = metric%M(1,1,:)
  end if

  call metric%deallo()

end subroutine compute_smoothSizeField
!
!
!
subroutine compute_sizeField_from_sol_viaHessian(hh,u,ndime,maxMeshSize,numTargetNodes)
  use mod_hessianReconstruction, only: compute_hessian_from_scalar
  use def_adapt,              only: memor_adapt
  use mod_memory,             only: memory_alloca, memory_deallo

  real(rp), pointer, intent(inout) :: hh(:)
  real(rp), pointer, intent(in)    :: u(:)
  integer(ip),       intent(in)    :: ndime
  real(rp)   ,       intent(in)    :: maxMeshSize
  integer(ip), optional, intent(in):: numTargetNodes

  type(mesh_metric_type) :: metric
  logical(lg) :: isIsotropic
  
  isIsotropic = .true.

  call compute_metric_from_sol_viaHessian(metric,u,ndime,maxMeshSize,isIsotropic)
  
  if(present(numTargetNodes)) then
    if(numTargetNodes>0_ip) then
      call metric%setMetricComplexity(numTargetNodes)
    end if
  end if

  nullify(hh)
  if(associated(u)) then
    call memory_alloca(memor_adapt,'hh','compute_sizeField_from_sol_viaHessian',hh,size(u))
    hh(:) = metric%M(1,1,:)
  end if

  call metric%deallo()

end subroutine compute_sizeField_from_sol_viaHessian
!
!
!
subroutine compute_sizeField_from_multiSol_viaHessian(hh,u,ndime,maxMeshSize,numTargetNodes)
  use mod_hessianReconstruction, only: compute_hessian_from_scalar
  use def_adapt,              only: memor_adapt
  use mod_memory,             only: memory_alloca, memory_deallo
!  use def_master, only: kfl_paral

  real(rp), pointer, intent(inout) :: hh(:)
  real(rp), pointer, intent(in)    :: u(:,:)
  integer(ip),       intent(in)    :: ndime
  real(rp)   ,       intent(in)    :: maxMeshSize
  integer(ip), optional, intent(in):: numTargetNodes

  type(mesh_metric_type) :: metric
  logical(lg) :: isIsotropic
  
  isIsotropic = .true.

  call compute_metric_from_multiSol_viaHessian(metric,u,ndime,maxMeshSize,isIsotropic)
  
  !print*,kfl_paral," --metric%isSizeField: ",metric%isSizeField
  
  if(present(numTargetNodes)) then
    if(numTargetNodes>0_ip) then
      call metric%setMetricComplexity(numTargetNodes)
    end if
  end if

  nullify(hh)
  if(associated(u)) then
    call memory_alloca(memor_adapt,'hh','compute_sizeField_from_sol_viaHessian',hh,size(u,1))
    hh(:) = metric%M(1,1,:)
  end if

  call metric%deallo()

end subroutine compute_sizeField_from_multiSol_viaHessian
!
!
!
subroutine compute_metric_from_sol_viaHessian(metric,u,ndime,hmax,isIsotropic)
  use mod_hessianReconstruction, only: compute_hessian_from_scalar,deallo_hessian_from_scalar
!  use def_adapt,              only: memor_adapt
  use mod_memory,             only: memory_alloca, memory_deallo
!  use def_domain, only: memor_dom

  real(rp), pointer, intent(in)    :: u(:)
  integer(ip),       intent(in)    :: ndime
  real(rp)   ,       intent(in)    :: hmax
  logical(lg),       intent(in)    :: isIsotropic
  type(mesh_metric_type),intent(inout) :: metric
  
  real(rp), pointer :: hessu(:,:,:)
  integer(ip) :: npoin
  
  if(.not.isIsotropic) then
    call runend('compute_metric_from_sol_viaHessian: Not implemented many things yet for anisotropic')
  end if
  !
  if(associated(u)) then
    npoin = size(u)
  else
    npoin = 0_ip
  end if
  
  call compute_hessian_from_scalar(u,ndime,hessu)
  
  call compute_metric_from_hessian(metric,hessu,ndime,hmax,isIsotropic)
  
  call compute_size_gradation(metric)  ! gradate should be an atribute of metric metric%gradate(optional type?, optional parameter?)
  
  if(associated(hessu)) then
    !call memory_deallo(memor_adapt,'HESSU','mod_metricComputation',hessu)
    !call memory_deallo(memor_dom,'HESSI','mod_projec',hessu)
    call deallo_hessian_from_scalar(hessu)
  end if
  
end subroutine compute_metric_from_sol_viaHessian
!
!
!
subroutine compute_metric_from_multiSol_viaHessian(metric,u,ndime,hmax,isIsotropic)
  use mod_hessianReconstruction, only: compute_hessian_from_scalar,deallo_hessian_from_scalar
  use mod_memory,             only: memory_alloca, memory_deallo
!  use def_domain, only: memor_dom

  real(rp), pointer, intent(in)    :: u(:,:)
  integer(ip),       intent(in)    :: ndime
  real(rp)   ,       intent(in)    :: hmax
  logical(lg),       intent(in)    :: isIsotropic
  type(mesh_metric_type),intent(inout) :: metric
  type(mesh_metric_type) :: metric_isol
  
  real(rp), pointer :: hessu(:,:,:), u_i(:)
  integer(ip) :: npoin
  
  integer(ip) :: isol, numSols
  
  if(.not.isIsotropic) then
    call runend('compute_metric_from_sol_viaHessian: Not implemented many things yet for anisotropic')
  end if
  !
  if(associated(u)) then
    npoin = size(u)
  else
    npoin = 0_ip
    call metric%alloca(ndime,npoin,isIsotropic)
    return
  end if
  
  numSols = size(u,2)
  do isol=1_ip,numSols
    
    nullify(u_i)
    u_i => u(:,isol)
    call compute_hessian_from_scalar(u_i,ndime,hessu)  
    
    if(isol==1) then
      call compute_metric_from_hessian(metric,hessu,ndime,hmax,isIsotropic)
    else
      call compute_metric_from_hessian(metric_isol,hessu,ndime,hmax,isIsotropic)
      call metric%intersect(metric_isol)
    end if
    
    !call memory_deallo(memor_dom,'HESSI','mod_projec',hessu)
  end do
  
  call deallo_hessian_from_scalar(hessu)
  
  !metric = intersectMetrics(metric_multi)
  call compute_size_gradation(metric)  ! gradate should be an atribute of metric metric%gradate(optional type?, optional parameter?)
  
end subroutine compute_metric_from_multiSol_viaHessian
!
!
!
subroutine compute_metric_from_hessian(metric,hessu,ndime,hmax,isIsotropic)
  use mod_hessianReconstruction, only: compute_hessian_from_scalar
  use mod_memory,             only: memory_alloca, memory_deallo

  real(rp), pointer, intent(in)    :: hessu(:,:,:)
  integer(ip),       intent(in)    :: ndime
  real(rp),          intent(in)    :: hmax
  logical(lg),       intent(in)    :: isIsotropic
  type(mesh_metric_type),intent(inout) :: metric
  
  integer(ip) :: inode, idim,npoin
  real(rp)    :: eig_val(ndime), eig_vec(ndime,ndime)
  real(rp)    :: h_node(1,1), eig_max, M_node(ndime,ndime)
  
  real(rp)    :: eig_bound, fact_maxh
!  real(rp)    :: cd, my_eps
  
  if(.not.isIsotropic) then
    call runend('compute_metric_from_hessian: Not implemented many things yet for anisotropic')
  end if
  !
  !
  if(associated(hessu)) then
    npoin = size(hessu,3)
  else
    npoin = 0_ip
  end if
  call metric%alloca(ndime,npoin,isIsotropic)
  
  !---- I started from this from the literature
!   cd        = 0.5_rp*(ndime/(1.0_rp + ndime))**2
!   my_eps    = 1.0_rp !0.1! desired tolerance (maybe I can improve it..or take it out..)
!   eig_bound = (my_eps/hmax**2)/cd ! from "Mesh-based anisotropic mesh adaptation", Alauzet
  !---- And switched to this.. directly imposing maximum mesh size in flat areas
  fact_maxh = 10.0_rp
  eig_bound = 1.0_rp/(fact_maxh*hmax)**2  ! here a not-infinite large mesh size that will be truncated after the complexity
  
  do inode=1,npoin
  
    call compute_eigen(ndime,hessu(:,:,inode),eig_val,eig_vec)
  
    if(isIsotropic) then
    
      eig_max = maxval(abs(eig_val))
      if(eig_max<eig_bound) then
       eig_max = eig_bound
      end if
      h_node  = 1.0_rp/sqrt(eig_max)
      
      if(h_node(1,1)<hmin_sec) then
        h_node = hmin_sec
      end if
      
      call metric%set_metric_node(inode,h_node)
    
    else
    
      eig_val = abs(eig_val)
      M_node = 0.0_rp
      !M_node(start\end\skip) = eig_val
      do idim=1,ndime
        M_node(idim,idim) = eig_val(idim)
      end do
      M_node = matmul(M_node,eig_vec)
      M_node = matmul(transpose(eig_vec),M_node)
    
      print*,'double check'
      call runend('double check anisotropic metric definition')
    
      call metric%set_metric_node(inode,M_node)
    
    end if
  
  end do
  
  ! compute eigenvalues (Alauzet): DONT USE THIS ---> use LOSEILLE VERSION to set in terms of dofs once all metric is defined
  !M=RΛ ̃tR where Λ ̃=diag(λ ̃), i
  !λ ̃=min(max(􏰉c|λ_i| ,h_max^-2) /eps,h_min^-2) 
  !c = 0.5*(d/(d+1))^2
  ! eps -> error -> can be changed to dofs
  !print*,'go from hess to  metric -> incorporate here notion of dofs? or do version via error/dofs'
  
end subroutine compute_metric_from_hessian
!
!
!
!-----------------------------------------------------------------------
!
!> @author  Abel Gargallo
!> @brief   compute_size_gradation
!> @details  see Borouchaci et al, MESH GRADATION CONTROL
!>   Two options:(1) Hv-correction and (2) Hc-correction
!>    
!>   (1) Hv-correction
!>  • Loop while the metric at a point is modi􏰁ed 
!>    – Loop over the edges of T(􏰀)
!>    ∗ Let PQ be the current edge
!>    ∗ Compute lP(PQ) (resp. lQ(PQ)) the length of PQ in the metric Md(P) (resp. Md(Q)) ∗
!>    ∗ Compute 􏰂P = (1 + 􏰅lP (PQ))−2 and 􏰂Q = (1 + 􏰅lQ (PQ))−2
!>    ∗ Replace Md(P) (resp. Md(Q)) by (Md(P);􏰂PMd(Q)) (resp. (Md(Q);􏰂QMd(P))).
!>    
!>   (2) Hc-correction
!>   • While the H-shock along an edge is ¿􏰁 – Loop over the edges of T(􏰀)
!>   * Let PQ be the current edge
!>   * Compute h(P) (resp. h(Q)) the unit length in the vicinity of P (resp. Q) in the direction −→
!>   PQ (suppose h(Q)¿h(P))
!>   * Compute l(PQ), the length of PQ
!>   * Compute c(PQ), the H-shock on PQ
!>   2
!>   * If c(PQ)¿􏰁 then replace Md(Q) by Md(Q)=􏰈 , where 􏰈=
!>   􏰊 􏰁 􏰋l(PQ) c(PQ)
!
!-----------------------------------------------------------------------
!
!
subroutine compute_size_gradation(metric)
  !
  use def_adapt,              only: memor_adapt
  use mod_memory,             only: memory_alloca, memory_deallo
  
  use mod_graphs, only: graphs_edges
  use def_domain, only: meshe, npoin
  use def_kermod, only: ndivi

  !
  implicit none
  !
  type(mesh_metric_type),intent(inout) :: metric
  !
  integer(ip), pointer :: edge_to_node(:,:)
  real(rp),    pointer :: coord(:,:)
  integer(ip) :: iedge, n1, n2
  integer(ip) :: nedge
  
  integer(ip) :: maxIter, iter
  
  logical(lg), pointer :: isNodeModified(:)
  logical(lg), pointer :: isNodeModified_prev(:)
  logical(lg) :: isNode1_mod, isNode2_mod
  
  real(rp) :: len_eucl
  
  integer(ip), parameter :: correction_hv = 1_ip
  integer(ip), parameter :: correction_hc = 2_ip ! USE THISSSS
  integer(ip), parameter :: type_correction  = correction_hc
  
  !print*,'compute_size_gradation.... ',do_gradate_metrics," ",type_correction," ",npoin
  
  if(.not.do_gradate_metrics) then
    return
  end if
  
  if(type_correction == correction_hv) then
    call runend('Better use HC correction than hv, looks better.')
  end if
  
  if(npoin.eq.0) then
    return
  end if
  
  if(.not.metric%isSizeField) then
    call runend('Imiplement compute_size_gradation for anisotropic meshes')
  end if
  
  
  if(.not.associated(meshe(ndivi) % edge_to_node)) then
    nullify( meshe(ndivi) % edge_to_node )
    nullify( meshe(ndivi) % ledgs        )
    nullify( meshe(ndivi) % lnned        )
    nullify( meshe(ndivi) % ledgb        )
    nullify( meshe(ndivi) % lnneb        )
    nullify( meshe(ndivi) % r_edg        )
    nullify( meshe(ndivi) % c_edg        )
  
    call graphs_edges(meshe(ndivi),memor_adapt)
  end if
  
  nedge        =  meshe(ndivi) % nedge
  nullify(edge_to_node)
  edge_to_node => meshe(ndivi) % edge_to_node
  nullify(coord)
  coord        => meshe(ndivi) % coord
  
  nullify(isNodeModified_prev)
  call memory_alloca(memor_adapt,'isNodeModified_prev',',mod_metricComputation',isNodeModified_prev,metric%npoin)
  nullify(isNodeModified)
  call memory_alloca(memor_adapt,'isNodeModified'     ,',mod_metricComputation',isNodeModified     ,metric%npoin)
  isNodeModified_prev(:) = .true.
  
  maxIter = 100_ip
  iter    = 0_ip
  do while( any(isNodeModified_prev).and.(iter<maxIter) )
    
    iter = iter+1_ip
    
    do iedge=1,nedge
  
      n1 = edge_to_node(1,iedge)
      n2 = edge_to_node(2,iedge)
      
      if( isNodeModified_prev(n1).or.isNodeModified_prev(n2) ) then
        
        !len_eucl = norm2( coord(:,n1) - coord(:,n2) )
        len_eucl = sqrt( sum( (coord(:,n1) - coord(:,n2))**2 ) )
    
        if(    type_correction.eq.correction_hv) then
          
          call compute_Hv_correction(len_eucl, metric%M(1,1,n1), metric%M(1,1,n2), isNode1_mod, isNode2_mod)
          
        elseif(type_correction.eq.correction_hc) then
          
          call compute_Hc_correction(len_eucl, metric%M(1,1,n1), metric%M(1,1,n2), isNode1_mod, isNode2_mod)
          
        else
          call runend('compute_size_gradation: implement other gradation types')
        end if
        
        isNodeModified(n1) = isNodeModified(n1).or.isNode1_mod
        isNodeModified(n2) = isNodeModified(n2).or.isNode2_mod
        
      end if
      
    end do
    
    isNodeModified_prev(:) = isNodeModified(:)
    isNodeModified = .false.
    
  end do
  
  call memory_deallo(memor_adapt,'isNodeModified_prev',',mod_metricComputation',isNodeModified_prev)
  call memory_deallo(memor_adapt,'isNodeModified',     ',mod_metricComputation',isNodeModified)
  
end subroutine compute_size_gradation
!
!
!
subroutine compute_Hv_correction(len_eucl,hp,hq,isNode_p, isNode_q)
  
  implicit none

  real(rp),    intent(in)    :: len_eucl
  real(rp),    intent(inout) :: hp,hq
  logical(lg), intent(inout)   :: isNode_p,isNode_q
  
  real(rp) :: lpq_p, lpq_q
  real(rp) :: gamma_p, gamma_q
  
  real(rp), parameter :: alpha = 1.3_rp ! 1.3
  
  call runend('DOES NOT WORK AS WELL AS THE HC')
  
  lpq_p = len_eucl / hp
  lpq_q = len_eucl / hq
  
  gamma_p = 1.0_rp + alpha*lpq_p
  gamma_q = 1.0_rp + alpha*lpq_q
  
  isNode_p = hp > hq*gamma_p
  if( isNode_p ) then ! metric intersection
    hp =   hq*gamma_p
  end if
  
  isNode_q = hq > hp*gamma_q
  if( isNode_q ) then ! metric intersection
    hq =   hp*gamma_q
  end if
  
end subroutine compute_Hv_correction
!
!
!
subroutine compute_Hc_correction(len_eucl,h1,h2,isNode1, isNode2)
  
  implicit none

  real(rp),    intent(in)    :: len_eucl
  real(rp),    intent(inout) :: h1,h2
  logical(lg), intent(inout)   :: isNode1, isNode2
  
  real(rp) :: hq_corr, cpq, lpq, hp, hq
  real(rp) :: eta
  
  real(rp), parameter :: beta = 4.0_rp ! 4.0_rp ! 2.0_rp !
  real(rp), parameter :: tol_h = 1.0e-12_rp
  
  isNode1 = .false.
  isNode2 = .false.
  
  if( abs(h1-h2)/max(h1,h2) < tol_h ) then
    !lpq = len_eucl / hp
    return ! if they are equal.. ther is no shock...
  elseif( h2>h1 ) then
    hp = h1
    hq = h2
  else
    hp = h2
    hq = h1
  end if
  
  lpq = len_eucl *( hq-hp ) / (hp*hq*log10(hq/hp))
  
  cpq = hshock(hp,hq,lpq)

  if( cpq > beta ) then
    
    eta = (beta/cpq)**lpq
    hq_corr = hq * eta 
    
    if( h2>h1 ) then
      h2 = hq_corr
      isNode1 = .false.
      isNode2 = .true.
    else
      h1 = hq_corr
      isNode1 = .true.
      isNode2 = .false.
    end if
    
  end if
  
end subroutine compute_Hc_correction
!
!
!
function hshock(hp,hq,lpq) result(c)
  
  implicit none
  real(rp), intent(in) :: hp,hq,lpq
  real(rp) :: c
  
  ! lpq = ||pq||/hp
  
  c = max(hp/hq,hq/hp)**(1.0_rp/lpq)
  
end function hshock
!
!
!
! subroutine compute_elemSizeField_from_sol_viaHessian(hh_elem,u,ndime,nelem,isIsotropic_inp)
!   use mod_hessianReconstruction, only: compute_hessian_from_scalar
!   use def_adapt,              only: memor_adapt
!   use mod_memory,             only: memory_alloca, memory_deallo
!
!   real(rp), pointer, intent(in)    :: u(:)
!   integer(ip),       intent(in)    :: ndime,nelem
!   logical(lg), optional, intent(in):: isIsotropic_inp
!   real(rp), pointer, intent(inout) :: hh_elem(:)
!   real(rp), pointer :: hh_node(:)
!
!   call compute_sizeField_from_sol_viaHessian(hh_node,u,ndime,isIsotropic_inp)
!
!   nullify(hh_elem)
!   call memory_alloca(memor_adapt,'hh_elem','compute_sizeField_from_sol_viaHessian',hh_elem,nelem)
!
!   call memory_deallo(memor_adapt,'hh','compute_elemSizeField_from_sol_viaHessian',hh)
!
! end subroutine compute_elemSizeField_from_sol_viaHessian
!
!
!
subroutine compute_eigen(ndime,M,values,vectors)
  use mod_maths_solver, only: maths_eigen_symmetric_matrix
  
  implicit none
  
  integer(ip), intent(in)  :: ndime
  real(rp),    intent(in)  :: M(ndime,ndime)
  real(rp),    intent(out) :: values(ndime), vectors(ndime,ndime)
  
  call maths_eigen_symmetric_matrix(ndime,M,values,vectors)
  
end subroutine compute_eigen

!
!
!
END MODULE mod_metricComputation

!> @}




  
!   ! compute eigenvalues (Alauzet): DONT USE THIS ---> use LOSEILLE VERSION to set in terms of dofs once all metric is defined
!   !M=RΛ ̃tR where Λ ̃=diag(λ ̃), i
!   !λ ̃=min(max(􏰉c|λ_i| ,h_max^-2) /eps,h_min^-2)
!   !c = 0.5*(d/(d+1))^2
!   ! eps -> error -> can be changed to dofs
!   !print*,'go from hess to  metric -> incorporate here notion of dofs? or do version via error/dofs'
!
!   do inode=1,npoin
!
!     print*,'compute_eigen'
!     call compute_eigen(ndime,hessu(:,:,inode),eig_val,eig_vec)
!     print*,'compute_eigen done'
!
!     ! compute the proper value (without correcting it yet)
!     if(isIsotropic) then
!
!       eig_max = maxval(abs(eig_vec))
!       h_node  = 1.0_rp/sqrt(eig_max)
!       call metric%set_metric_node(inode,h_node)
!
!     else
!
!       eig_val = abs(eig_val)
!       M_node = 0.0_rp
!       !M_node(start\end\skip) = eig_val
!       do idim=1,ndime
!         M_node(idim,idim) = eig_val(idim)
!       end do
!       M_node = matmul(M_node,eig_vec)
!       M_node = matmul(transpose(eig_vec),M_node)
!
!       print*,'double check'
!       call runend('double check anisotropic metric definition')
!
!       call metric%set_metric_node(inode,M_node)
!
!     end if
!
!   end do
!
!   call memory_deallo(memor_adapt,'HESSU','compute_metric_from_sol_viaHessian',hessu)
