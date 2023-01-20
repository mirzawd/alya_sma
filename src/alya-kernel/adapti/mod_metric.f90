!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Adaptivity
!> @{
!> @file    mod_metric.f90
!> @author  abel.gargallo
!> @date    2021-03-31
!> @brief   mod_metric
!> @details mod_metric
!>
!>          To add further details
!>
!>
!-----------------------------------------------------------------------

MODULE mod_metric
!***************************************************************
!*
!*  Module for metric computation and operations
!*
!***************************************************************
use def_kintyp_basic, only : ip,rp,lg

implicit none

private
!
!
!save
!
! type metric
!   real(rp), pointer :: M(:,:)
! !    character(20)                      :: dim
! !  contains
! !    !procedure,                 pass    :: init                   ! Initialize all
! end type metric
!
!
!
type mesh_metric_type
   integer(ip) :: npoin          ! number of metric points stored
   integer(ip) :: dim            ! dim of the space to construct a dim x dim metric
   logical(lg) :: isSizeField    ! if it is a size field then true, general metric (dimxdim) then false
   real(rp), pointer :: M(:,:,:) ! dim x dim x npoin for general metrics, 1 x 1 x npoin for sizeField
   !integer(ip)                      :: dim_allocated?
   integer(ip) :: numTargetNodes ! numTargetNodes to reproduce metric
 contains
   !procedure,                 pass    :: init                   ! 
   procedure, pass :: alloca                  =>   allocate_mesh_metric
   procedure, pass :: deallo                  => deallocate_mesh_metric
   procedure, pass :: copy
   procedure, pass :: intersect               => intersectMetrics
   procedure, pass :: get_metric_node
   procedure, pass :: get_size_node
   procedure, pass :: get_metric_elem
   procedure, pass :: get_metric_elem_newCav
   procedure, pass :: interpolate             => interpolate_metric_midpoint
   procedure, pass :: generate_edge_point
   procedure, pass :: set_metric_node
   procedure, pass :: merge                   => merge_metrics
   procedure, pass :: remove_nodes            => remove_deleted_nodes
   !generic,   public                  :: assignment(=) => equal
   procedure, pass :: setMetricComplexity
end type mesh_metric_type


integer(ip), parameter :: law_linear_h                     = 0
integer(ip), parameter :: law_linear_eigen                 = 1
integer(ip), parameter :: lwa_linear_hinv                  = 2
integer(ip), parameter :: law_geometric_h                  = 3
integer(ip), parameter :: interpolation_law = law_geometric_h

integer(ip), parameter :: p_norm_metric_default = 4_ip ! default: 2


interface compute_metric_fromNodalSizeField
   module procedure compute_metric_fromNodalSizeField_mesh,&
        &           compute_metric_fromNodalSizeField_ndime
end interface compute_metric_fromNodalSizeField


public :: mesh_metric_type
public :: compute_metric_fromElemSizeField
public :: compute_metric_fromNodalSizeField

public :: law_linear_h, law_linear_eigen, lwa_linear_hinv, law_geometric_h ! From Alauzet, Size gradation
public :: interpolation_law
!
!
!
CONTAINS
!
!
!
subroutine copy(metric,metric_copy)
  !
  implicit none
  !
  class(mesh_metric_type), intent(in   ) :: metric
  type(mesh_metric_type),  intent(inout) :: metric_copy
  !
  call metric_copy%alloca(metric%dim,metric%npoin,metric%isSizeField)
  !
  metric_copy%M = metric%M
  !
end subroutine
!
!
!
subroutine compute_metric_fromElemSizeField(elemSizeField,mesh,metric)
  !
  use def_kintyp_mesh_basic,  only: mesh_type_basic
  !
  implicit none
  !
  real(rp),               pointer, intent(in)    :: elemSizeField(:)
  type(mesh_type_basic)          , intent(in)    :: mesh
  class(mesh_metric_type),         intent(inout) :: metric
  
  integer(ip) :: adjacencyNode(mesh%npoin)
  real(rp) :: nodal_size(mesh%npoin)
  
  logical(lg) :: isSizeField = .true.
  
  integer(ip) :: inode,ielem,theNode
  
  call metric%alloca(mesh%ndime,mesh%npoin,isSizeField)
  
  adjacencyNode = 0_ip
  nodal_size = 0.0_rp
  do ielem=1,mesh%nelem
    do inode=1,(mesh%ndime+1)! only for simplices
      
      theNode = mesh%lnods(inode,ielem)
      
      adjacencyNode(theNode) = adjacencyNode(theNode) + 1_ip
      nodal_size(theNode) = nodal_size(theNode) + elemSizeField(ielem)
    end do
  end do
  nodal_size = nodal_size/adjacencyNode
  
  metric%M(1,1,:) = nodal_size
  
  return
end subroutine compute_metric_fromElemSizeField
!
!
!
subroutine compute_metric_fromNodalSizeField_mesh(nodalSizeField,mesh,metric)

  use def_kintyp_mesh_basic,  only: mesh_type_basic
  !
  implicit none
  ! 
  real(rp),               pointer, intent(in)    :: nodalSizeField(:)
  type(mesh_type_basic)          , intent(in)    :: mesh
  class(mesh_metric_type),         intent(inout) :: metric
  
  logical(lg) :: isSizeField = .true.
  
  call metric%alloca(mesh%ndime,mesh%npoin,isSizeField)
  
  metric%M(1,1,:) = nodalSizeField
  
  return
end subroutine compute_metric_fromNodalSizeField_mesh
!
!
!
subroutine compute_metric_fromNodalSizeField_ndime(nodalSizeField,ndime,metric)

  use def_kintyp_mesh_basic,  only: mesh_type_basic
  !
  implicit none
  !
  
  real(rp),               pointer, intent(in)    :: nodalSizeField(:)
  integer(ip)                    , intent(in)    :: ndime
  class(mesh_metric_type),         intent(inout) :: metric
  
  logical(lg) :: isSizeField = .true.
  integer(ip) :: npoin
  
  if(.not.associated(nodalSizeField)) then
    call metric%alloca(ndime,0_ip,isSizeField)
    return
  end if
  
  npoin = size(nodalSizeField)
  call metric%alloca(ndime,npoin,isSizeField)
  
  metric%M(1,1,:) = nodalSizeField
  
  return
end subroutine compute_metric_fromNodalSizeField_ndime
!
!
!
subroutine allocate_mesh_metric(metric,dim,npoin,isSizeField)
  use mod_memory,             only: memory_alloca
  use def_adapt,              only: memor_adapt
  implicit none
  
  class(mesh_metric_type),  intent(inout) :: metric
  integer(ip), intent(in) :: npoin, dim
  logical(lg), intent(in) :: isSizeField
  
  metric%dim         = dim
  metric%npoin       = npoin
  metric%isSizeField = isSizeField
  
  nullify(metric%M)
  if(isSizeField) then
    call memory_alloca(memor_adapt,'metric%M','mod_metric',metric%M,1_ip,1_ip,npoin)
  else !general anisotropic metric
    call memory_alloca(memor_adapt,'metric%M','mod_metric',metric%M,dim,dim,npoin)
  end if
  
  metric%numTargetNodes = 0_ip
  
  return
end subroutine allocate_mesh_metric
!
!
!
subroutine deallocate_mesh_metric(metric)
  use mod_memory,             only: memory_deallo
  use def_adapt,              only: memor_adapt
  implicit none
  class(mesh_metric_type),  intent(inout) :: metric
  
  metric%dim      = 0_ip
  metric%npoin    = 0_ip
  
  call memory_deallo(memor_adapt,'metric%M','mod_metric',metric%M)
  
  return
end subroutine deallocate_mesh_metric
!
!
! 
subroutine merge_metrics(metric1,metric2,numPoints2,isDeletedPoint1)
  use mod_debugTools, only: out_performance, deb_metricMerge
  
  implicit none
  
  class(mesh_metric_type),  intent(inout) :: metric1
  type(mesh_metric_type),   intent(inout) :: metric2
  integer(ip), intent(in) :: numPoints2
  logical(lg), intent(in), optional    :: isDeletedPoint1(:)
  
  integer(ip) :: numPoints, numPointsDelete1,ipoin,countPoints
  real(rp)    :: M1(size(metric1%M,1),size(metric1%M,2),metric1%npoin)
  
  real(rp) :: t0,t1
  !
  if(out_performance) call cpu_time(t0)
  
  M1 = metric1%M
  
  if(present(isDeletedPoint1)) then
    print*,'error-> warning warning, cannot delete from emtric if mesh preserves them'
    call runend("check if used")
    numPointsDelete1 = count(isDeletedPoint1)
  else
    numPointsDelete1 = 0
  end if
  numPoints = metric1%npoin - numPointsDelete1  + numPoints2

  !call memory_deallo(memor_adapt,'metric1%M','merge_metrics',metric1%M)
  !nullify(metric1%M)
  !call memory_alloca(memor_adapt,'metric1%M','merge_metrics',metric1%M,metric1%dim,metric1%dim,numPoints)
  call metric1%deallo()
  call metric1%alloca(metric2%dim,numPoints,metric2%isSizeField)
  
  if(present(isDeletedPoint1)) then
    countPoints = 0
    do ipoin=1,metric1%npoin
      if(.not.isDeletedPoint1(ipoin)) then
        countPoints = countPoints+1
        metric1%M(:,:,countPoints) = M1(:,:,ipoin)
      end if
    end do
  else
    countPoints = size(M1,3)
    metric1%M(:,:,1:countPoints) = M1
  end if
  
  countPoints = countPoints+1
  !metric1%M(:,:,countPoints:numPoints) = metric2%M(:,:,1:numPoints2)
  do ipoin=1,numPoints2
    metric1%M(:,:,countPoints+ipoin-1) = metric2%M(:,:,ipoin)
  end do 
    
  metric1%npoin = numPoints
  
  call metric2%deallo()
  
  if(out_performance) then
    call cpu_time(t1)
    deb_metricMerge = deb_metricMerge + (t1-t0)
  end if
  
  return
end subroutine merge_metrics
!
!
! 
subroutine remove_deleted_nodes(metric,isDeletedPoint)
  use mod_memory,             only: memory_alloca, memory_deallo
!  use def_adapt,              only: memor_adapt
  implicit none
  
  class(mesh_metric_type),  intent(inout) :: metric
  !logical(lg), intent(in), optional    :: isDeletedPoint(:)
  logical(lg), pointer, intent(in), optional    :: isDeletedPoint(:)
  
  integer(ip) :: numPoints, numPointsDelete,ipoin,countPoints,numPoints_preAdapt
  real(rp)    :: M(size(metric%M,1),size(metric%M,2),size(metric%M,3))
  integer(ip) :: dim,npoin
  logical(lg) :: isSizeField
  
  numPointsDelete = count(isDeletedPoint)
  
  if(numPointsDelete>0) then
    dim=metric%dim
    npoin=metric%npoin
    isSizeField=metric%isSizeField
    M(:,:,:) = metric%M
  
    numPoints = metric%npoin - numPointsDelete 
    
    call metric%deallo()
    call metric%alloca(dim,numPoints,isSizeField)

    countPoints = 0
    numPoints_preAdapt = size(isDeletedPoint)
    do ipoin=1,numPoints_preAdapt
      if(.not.isDeletedPoint(ipoin)) then
        countPoints = countPoints+1
        metric%M(:,:,countPoints) = M(:,:,ipoin)
      end if
    end do
    !metric%M(:,:,(countPoints+1):numPoints) = M(:,:,(numPoints_preAdapt+1):npoin)
    do ipoin=1,(numPoints-countPoints)
      metric%M(:,:,countPoints+ipoin) = M(:,:,(numPoints_preAdapt+ipoin))
    end do 

    metric%npoin = numPoints
  end if
  
  return
end subroutine remove_deleted_nodes
!
!
!
! subroutine set_metric(metric,dim,Mvalues)
!   implicit none
!   class(mesh_metric_type),  intent(inout) :: metric
!   integer(ip),              intent(in)    :: dim
!   real(rp),                 intent(in)    :: Mvalues(:,:,:)
!
!   call metric%alloca(dim,npoin,isSizeField)
!
!   metric%M(:,:,inode) = Mnode
!
! end subroutine set_metric
!
!
!
subroutine set_metric_node(metric,inode,Mnode)
  implicit none
  class(mesh_metric_type),  intent(inout) :: metric 
  integer(ip),              intent(in)    :: inode
  real(rp),                 intent(in)    :: Mnode(:,:)
  !
  metric%M(:,:,inode) = Mnode
  !
end subroutine set_metric_node
!
!
!
function get_metric_node(M,inode) result(Mnode)
  implicit none
  class(mesh_metric_type),  intent(in)    :: M 
  integer(ip),              intent(in)    :: inode
  
  real(rp)                          :: Mnode(M%dim,M%dim)
!  integer(ip) :: idim
  real(rp) :: h_node
!  real(rp) :: m_const
  
  if(M%isSizeField) then
    h_node = M%M(1,1,inode)
    Mnode = computeMetric_fromSizeField(h_node,M%dim)
  else
    Mnode = M%M(:,:,inode)  
  end if
  
end function get_metric_node
!
!
!
function get_size_node(M,inode) result(h_node)
  implicit none
  class(mesh_metric_type),  intent(in)    :: M 
  integer(ip),              intent(in)    :: inode
  
  real(rp) :: h_node
  
  if(M%isSizeField) then
    h_node = M%M(1,1,inode)
  else
    call runend('Anisotropic metric, cannot recall get_size_node')
  end if
  
end function get_size_node
!
!
!
function computeMetric_fromSizeField(h,dim) result(M)
  implicit none
  
  real(rp),     intent(in) :: h
  integer(ip),  intent(in) :: dim
  real(rp) :: M(dim,dim)
  
  real(rp) :: m_diag
  integer(ip) :: idim
  
  m_diag = (1.0_rp/h)**2
  
  M(:,:) = 0.0_rp
  do idim=1,dim
    M(idim,idim) = m_diag
  end do
  
end function computeMetric_fromSizeField
!
!
!
function get_metric_elem(M,Telem) result(Melem)
  implicit none
  
  class(mesh_metric_type),  intent(in)  :: M 
  integer(ip),              intent(in)  :: Telem(:)
  real(rp)                              :: Melem(M%dim,M%dim)
  
  real(rp) :: h_elem
  
  if(M%isSizeField) then
    if(interpolation_law.eq.law_geometric_h) then

      h_elem = PRODUCT( M%M( 1,1,Telem ) )
      h_elem = h_elem**(1.0_rp/size(Telem))

    else
      !h_elem = sum( M%M(1,1,Telem) )/size(Telem)
      call runend("implement interpolation law")
    end if
    !
    !h_elem = minval( M%M( 1,1,Telem ) )
    !
    Melem = computeMetric_fromSizeField(h_elem,M%dim)
    !
  else
    Melem(:,:) = sum( M%M(:,:,Telem) , 3 )/size(Telem)
    call runend("Improve interpolation metric element for non-size field")
  end if
    
end function get_metric_elem
!
!
!
function get_metric_elem_newCav(M,Telem,m_newNode) result(Melem)
  implicit none
  class(mesh_metric_type),  intent(in)    :: M 
  integer(ip),              intent(in)    :: Telem(:)
  real(rp),                 intent(in)    :: m_newNode(:,:)
  
  real(rp)                          :: Melem(M%dim,M%dim)
  real(rp) :: h_elem
  
  if(M%isSizeField) then

    if(interpolation_law.eq.law_geometric_h) then
      h_elem = PRODUCT( M%M( 1,1,Telem(1:(size(Telem)-1)) ) )* m_newNode(1,1)
      h_elem = h_elem**(1.0_rp/size(Telem))
      !h_elem = min( minval( M%M( 1,1,Telem(1:(size(Telem)-1)) ) ) , m_newNode(1,1) )
    else
      call runend("implement interpolation law")
    end if
    !h_elem = (sum( M%M( 1,1,Telem(1:(size(Telem)-1)) ) ) + m_newNode(1,1) )/size(Telem)
    Melem = computeMetric_fromSizeField(h_elem,M%dim)
  else
    call runend("Improve interpolation metric element for non-size field")
  end if
    
end function get_metric_elem_newCav
!
!
!
function interpolate_metric_midpoint(M,n1,n2) result(M_mid)
  
  implicit none

  class(mesh_metric_type), intent(in)    :: M 
  integer(ip),             intent(in)    :: n1,n2
  
  real(rp)                          :: M_mid(size(M%M,1),size(M%M,2))!M%dim,M%dim)
  
  call runend('should not be used any more...')
  
  if(m%isSizeField) then
    !!! Linear on h:
!     M_mid(:,:) = ( M%M(:,:,n1)+M%M(:,:,n2) )/2.0_rp
    !!! Geometric on h
    if(interpolation_law.eq.law_geometric_h) then
      M_mid(:,:) = sqrt( M%M(:,:,n1)*M%M(:,:,n2) )
    else
      call runend('interpolate_metric_midpoint: defaulted geometric law on h')
    end if
    !!! Linear interp on h^-1
    
    !!! Linear interp on eigen
    !
  else
    M_mid(:,:) = ( M%M(:,:,n1)+M%M(:,:,n2) )/2.0_rp
    call runend("Improve interpolate_metric_midpoint for non-size field")
  end if
  
  return 
end function interpolate_metric_midpoint
!
!
!
subroutine generate_edge_point(M,n1,n2,x1,x2,xnew,Mnew)
  
  implicit none

  class(mesh_metric_type), intent(in)    :: M 
  integer(ip),             intent(in)    :: n1,n2
  real(rp),                intent(in)    :: x1(:), x2(:)
  real(rp),                intent(out)   :: xnew(:)
  real(rp),                intent(out)   :: Mnew(size(M%M,1),size(M%M,2))!M%dim,M%dim)
  
  real(rp)  :: t, hp, hq, xp(M%dim), xq(M%dim)
  
  if(m%isSizeField) then

    if(interpolation_law.eq.law_geometric_h) then
      
      if( M%M(1,1,n1) < M%M(1,1,n2) ) then
        hp = M%M(1,1,n1)
        hq = M%M(1,1,n2)
        xp = x1
        xq = x2
      else
        hp = M%M(1,1,n2)
        hq = M%M(1,1,n1)
        xp = x2
        xq = x1
      end if
      
      t = 1.0_rp -  (1.0_rp/hp) / ( (1.0_rp/hp) + (1.0_rp/hq)  ) 
      
      xnew = xp + t*(xq-xp)
      
      Mnew(1,1) =  (  hp**(1.0_rp-t)  )*(   hq**t   )
      
    else
      call runend('interpolate_metric_midpoint: defaulted geometric law on h')
    end if
    
  else
    !M_mid(:,:) = ( M%M(:,:,n1)+M%M(:,:,n2) )/2.0_rp
    call runend("Improve interpolate_metric_midpoint for non-size field")
  end if
  
  return 
end subroutine generate_edge_point

!
!
!
subroutine setMetricComplexity(metric,numTargetNodes,pnorm_inp)
  !
  ! see CONTINUOUS MESH FRAMEWORK II, from Adrien Loiselle and Frederic Alauzet
  !
  use mod_memory,    only: memory_alloca, memory_deallo
  use def_adapt,     only: memor_adapt
!  use def_master,    only: kfl_paral
  use def_domain,    only: npoin
  
  implicit none
  
  class(mesh_metric_type), intent(inout) :: metric
  integer(ip), intent(in) :: numTargetNodes
  integer(ip), optional, intent(in) :: pnorm_inp
  
  real(rp)                :: scalar_factor
  real(rp), pointer       :: detM(:)
  real(rp), pointer       :: integrand(:)
  integer(ip)             :: pnorm
!  integer(ip)             :: inode
  
  integer(ip) :: d, dim
  
  if(present(pnorm_inp)) then
    pnorm = pnorm_inp
  else
    pnorm = p_norm_metric_default
  end if
  
  d   = metric%dim ! this is my claim, I replaced the 3's in the paper by dimension.. recheck if necessary
  dim = metric%dim 
  
  metric%numTargetNodes = numTargetNodes
  
  if(metric%isSizeField) then
    ! I am writting it in terms of a size field instead of a metric
    nullify(detM)
    call memory_alloca(memor_adapt,'detM',     'mod_metric',detM,     npoin)
    nullify(integrand)
    call memory_alloca(memor_adapt,'integrand','mod_metric',integrand,npoin)
    !call memory_alloca(memor_adapt,'detp','mod_metric',detp,dimensionCambiante)
    !allocate(detp(metric%npoin)) ! using memory_alloca theres some error in master when dimension is 0 !XXYYZZ
    if(metric%npoin>0_ip) then
      !detp(:) = metric%M(1,1,:)**(-1.0_rp*pnorm*d/((2.0_rp*pnorm)+d))
      detM      = 1.0_rp/metric%M(1,1,:)**(2_ip*dim) 
      integrand = detM**(pnorm/( (2.0_rp*pnorm)+d ) )
      
    end if
    
    !print*,'integrate_scalar... kfl_paral:',kfl_paral
    !print*,integrand
    scalar_factor = integrate_scalar(integrand)
    !print*,"scalar_factor: ",scalar_factor,'    kfl_paral:',kfl_paral
    !print*,'...integrate_scalar kfl_paral:',kfl_paral
    !call runend('integro un campo>1 y me da <1... siendo la medida del dominio 1..?')
    
    if(metric%npoin>0_ip) then
      !scalar_factor = (scalar_factor/numTargetNodes)**(2.0_rp/d)
      !metric%M(1,1,:) = scalar_factor*(detM(:)**(-1.0_rp/pnorm))*metric%M(1,1,:)
      scalar_factor   = (numTargetNodes/scalar_factor)**(2.0_rp/d) 
      !detM            = detM*(scalar_factor*(detM**(-1.0_rp/(2.0_rp*pnorm + d))))**dim
      !metric%M(1,1,:) = 1.0_rp/detM**( 1.0_rp/(2_ip*dim) )
      ! alternatively:
      detM            = ( scalar_factor*(detM**(-1.0_rp/(2.0_rp*pnorm + d))) )*( detM**(1.0_rp/dim) )
      metric%M(1,1,:) = 1.0_rp/sqrt(detM)
    end if

    call memory_deallo(memor_adapt,'detM',     'mod_metric',     detM)
    call memory_deallo(memor_adapt,'integrand','mod_metric',integrand)
    !deallocate(detp)
    
  else
    
    print*,'IMPLEMENT...'
    call runend('implement setMetricComplexity for anisotropic metrics (easy, just do it)')
    ! its the same but with detM, instead of 1/h^dim
!     do inode=1,metric%npoin
!       nodal_field(inode) =
!     end do
    
  end if
  
end subroutine setMetricComplexity
!
!
!
function integrate_scalar(field) result(integral)
  use mod_func,      only: func_ptr, func_initialization, func_identity
  use mod_integrals, only: field_arrays, integrals_volume
  use def_domain,    only: nelem, npoin
  use mod_memory,    only: memory_alloca, memory_deallo
  use def_adapt,     only: memor_adapt
!  use def_master,    only: kfl_paral
  
  implicit none
  
  real(rp), pointer, intent(in) :: field(:)
  real(rp) :: integral

  type(func_ptr)                                      :: array_func(1,1)
  type(field_arrays)                                  :: array_fields(1)
  real(rp)                                            :: integral_clusters(1,1) !(1,nclus) see mod_clusters
  integer(ip),           pointer   :: legro(:)
  
  nullify(legro)
  call memory_alloca(memor_adapt,'legro','mod_metric',legro,nelem)
  if(nelem>0_ip) then
    legro(:) = 1_ip
  end if

  call func_initialization(array_func)
  array_func(1,1) % f => func_identity

  nullify(array_fields(1)%a)
  call memory_alloca(memor_adapt,'array_fields','mod_metric',array_fields(1)%a,1_ip,npoin)
  if(npoin>0_ip) then
    array_fields(1) % a(1,:) = field(:)
  end if
  
  !print*,'integrals_volume...: kfl_paral:',kfl_paral
  call integrals_volume(array_fields,array_func,legro,integral_clusters)
  !print*,'...integrals_volume: kfl_paral:',kfl_paral

  call memory_deallo(memor_adapt,'array_fields','mod_metric',array_fields(1)%a)
  call memory_deallo(memor_adapt,'legro','mod_metric',legro)

  integral = integral_clusters(1,1)
  
end function integrate_scalar
!
!
!
! function intersectMetricList(metrics) result(metric_intersect)
!   implicit none
!   !
!   type(mesh_metric_type), intent(in) :: metrics(:)
!   type(mesh_metric_type) :: metric_intersect
!   !
!   integer(ip) :: imetric
!   !
!   call metrics(1_ip)%copy(metric_intersect)
!
!   do imetric = 2_ip,size(metrics)
!     call intersectMetrics(metric_intersect,metrics(imetric))
!   end do
!
! end function intersectMetricList
!
!
!
subroutine intersectMetrics(metric1,metric2)
  implicit none 
  !
  class(mesh_metric_type), intent(inout) :: metric1
  type(mesh_metric_type) , intent(in)    :: metric2
  !
  logical(lg) :: isSizeField
  
  isSizeField = metric1%isSizeField
  
  if(.not.isSizeField) then
    call runend('Anisotropic metric intersection not computed yet -> TO CODE...')
  end if
  
  metric1%M(:,1,1)=min(metric1%M(:,1,1),metric2%M(:,1,1))
  
end subroutine intersectMetrics
!
!
!
END MODULE mod_metric

!> @}
