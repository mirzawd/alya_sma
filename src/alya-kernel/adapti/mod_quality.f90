!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Adaptivity
!> @{
!> @file    mod_quality.f90
!> @author  abel.gargallo
!> @date    2021-03-31
!> @brief   mod_quality
!> @details mod_quality
!>
!>          To add further details
!>
!>
!-----------------------------------------------------------------------

MODULE mod_quality
!***************************************************************
!*
!*  Module for computing mesh quality
!*
!***************************************************************
use def_kintyp_basic,       only : ip,rp,lg
use def_kintyp_mesh_basic,  only: mesh_type_basic
use mod_metric,             only: mesh_metric_type

use def_adapt,              only : memor_adapt 
use mod_memory,             only: memory_alloca, memory_deallo

implicit none

!
! Method parameters
!
! integer(ip), parameter :: exp_size  = 2_ip
! integer(ip), parameter :: exp_shape = 3_ip
integer(ip), parameter :: exp_size  = 1_ip
integer(ip), parameter :: exp_shape = 1_ip
! real(rp), parameter :: exp_size  = 1.0_rp
! real(rp), parameter :: exp_shape = 1.1_rp

!
! Private variables
!
!real(rp), parameter :: delta_reg_zero = 0.000001_rp
real(rp), parameter :: delta_reg_zero = 0.0001_rp

real(rp), private, target :: Winv_tri(2,2) = &
  reshape((/ &
               1.0,             0.0,&
    -1.0/sqrt(3.0),    2.0/sqrt(3.0) &
    /),  (/2,2/)  )
    
real(rp), private, target :: Winv_tet(3,3) = &
  reshape((/ &
                         1.0,                      0.0,                     0.0, &
        -(1.0/3.0)*sqrt(3.0),      (2.0/3.0)*sqrt(3.0),                     0.0, &
    -(1.0/3.0)*sqrt(6.0)/2.0, -(1.0/3.0)*sqrt(6.0)/2.0,            sqrt(6.0)/2.0 &
    /),  (/3,3/)  )!/(sqrt(2.0))**(1.0/3.0)
! real(rp), private, target :: Winv_tri(2,2) = &
!   transpose(reshape((/ &
!     1.0, -1.0/sqrt(3.0),&
!     0.0,  2.0/sqrt(3.0) &
!     /),  (/2,2/)  ))!*(sqrt(3.0)/2.0)**(1.0/2.0)
!
! real(rp), private, target :: Winv_tet(3,3) = &
!   transpose(reshape((/ &
!     1.0, -(1.0/3.0)*sqrt(3.0), -(1.0/3.0)*sqrt(6.0)/2.0, &
!     0.0,  (2.0/3.0)*sqrt(3.0), -(1.0/3.0)*sqrt(6.0)/2.0, &
!     0.0,     0.0             ,            sqrt(6.0)/2.0 &
!     /),  (/3,3/)  ))!/(sqrt(2.0))**(1.0/3.0)

integer(ip), private :: dimW_saved = 0
real(rp), private, pointer :: Winv(:,:)

private

public :: compute_mesh_quality_sizeShape      ,& ! used in adapt
          compute_minQuality_sizeShape_subset ,& ! used in adapt
          compute_minQuality_sizeShape_cavity ,& ! used in adapt
          compute_quality_sizeShape_cavity    ,& ! used in adapt
          compute_quality_sizeShape_subset    ,& ! used in smoothing
          compute_distortion_sizeShape_subset ,& ! used in smoothing
          compute_mesh_quality_shape          ,&         
          quality_deallo                      ,&
          do_print_q_stats                    ,&
          print_q_stats                       ,&
          isPointInTet


!save

CONTAINS
!
!
!
subroutine quality_deallo(q)
  real(rp),   pointer ,intent(inout)  :: q(:)
  
  call memory_deallo(memor_adapt,'q','compute_mesh_quality',q)
  
end subroutine quality_deallo
!
!
!
subroutine compute_mesh_quality_sizeShape(mesh,metric,q)
  !***************************************************
  !*
  !*   Computes the quality of the 3D mesh
  !*   Output: quality_sizeShape 
  !*
  !***************************************************
  use mod_debugTools, only: out_performance, deb_qualityMesh
  implicit none
  
  type(mesh_type_basic),   intent(in) :: mesh
  type(mesh_metric_type),  intent(in) :: metric
  real(rp),   pointer ,intent(inout)  :: q(:)
  
  real (rp)   :: Xelem(mesh%ndime,mesh%mnode)
  integer(ip) :: ielem
  real(rp)    :: M(mesh%ndime,mesh%ndime)
  
  real(rp) :: t0,t1
  !
  if(out_performance) call cpu_time(t0)
  !
  call setWinv(mesh%ndime)
  
  nullify(q)
  call memory_alloca(memor_adapt,'q','compute_mesh_quality',q,mesh%nelem)
  
  do ielem=1,mesh%nelem
    Xelem    = mesh%coord( : , mesh%lnods(:,ielem) )
    M        = metric%get_metric_elem(mesh%lnods(:,ielem))
    q(ielem) = quality_sizeShape(Xelem,M)
  end do
  
  if(out_performance) then
    call cpu_time(t1)
    deb_qualityMesh = deb_qualityMesh + (t1-t0)
  end if
  
end subroutine compute_mesh_quality_sizeShape
!
!
!
subroutine compute_mesh_quality_shape(mesh,q)

  !***************************************************
  !*
  !*   Computes the quality of the 3D mesh
  !*   Output: quality3D 
  !*
  !***************************************************
  
  implicit none
  
  type(mesh_type_basic), intent(in)   :: mesh
  real(rp),   pointer ,intent(inout)  :: q(:)
  
  real (rp)   :: Xelem(mesh%ndime,mesh%mnode)
  integer(ip) :: ielem
  real(rp)    :: M(mesh%ndime,mesh%ndime)
  integer(ip) :: idime
  
  call setWinv(mesh%ndime)
  
  nullify(q)
  call memory_alloca(memor_adapt,'q','compute_mesh_quality',q,mesh%nelem)
  
  M = 0.0
  do idime=1,mesh%ndime
    M(idime,idime)=1.0_rp
  end do
  
  do ielem=1,mesh%nelem
    Xelem = mesh%coord( : , mesh%lnods(:,ielem) )
    q(ielem) = quality_shape(Xelem,M)
  end do

end subroutine compute_mesh_quality_shape
!
!
!
function compute_minQuality_sizeShape_cavity(X,T,metric,p_new,m_new) result(q_min)

  !***************************************************
  !*
  !*   Computes the quality of the 3D mesh
  !*   Output: quality3D 
  !*
  !***************************************************
  use mod_debugTools, only: out_performance, deb_minQCavity
  implicit none
  
  real(rp),                intent(in) :: X(:,:) ! coords
  integer(ip),             intent(in) :: T(:,:) ! connectivity
  type(mesh_metric_type),  intent(in) :: metric
  real(rp),    optional,   intent(in) :: p_new(:)
  real(rp),    optional,   intent(in) :: m_new(:,:)
  
  real(rp) :: q_min
  
  real (rp)   :: elemQual
  integer(ip) :: ielem
  real(rp)    :: Xelem(size(X,1),size(T,1))
  real(rp)    :: Melem(size(X,1),size(X,1))
  
  real(rp) :: t0,t1
  !
  if(out_performance) call cpu_time(t0)
  !
!   if(present(p_new).and.(T(size(T,1),1)<=size(X,2))) then
!     print*,"Last node should be new, that is  out of X, check mod_cavity for sincronization with this function.."
!     call runend("Last node should be new, that is  out of X, check mod_cavity for sincronization with this function..")
!   end if

  call setWinv( int(size(X,1),ip) )
  q_min = 1.0

  if(present(p_new)) then
    do ielem=1,size(T,2)
      Xelem(:,1:(size(T,1)-1)) = X( : , T( 1:(size(T,1)-1),ielem ) )
      Xelem(:,size(T,1))       = p_new
      Melem = metric%get_metric_elem_newCav(T(:,ielem),m_new)
      elemQual = quality_sizeShape(Xelem,Melem)
      q_min = min(elemQual,q_min)
    end do
  else
    do ielem=1,size(T,2)
      Xelem = X( : , T(:,ielem) )
      Melem = metric%get_metric_elem(T(:,ielem))
      elemQual = quality_sizeShape(Xelem,Melem)
      q_min = min(elemQual,q_min)
    end do
  end if

  if(out_performance) then
    call cpu_time(t1)
    deb_minQCavity = deb_minQCavity + (t1-t0)
  end if
!
!   loop_elements: do ielem=1,size(T,2)
!
!     if(present(p_new)) then
!       Xelem(:,1:(size(T,1)-1)) = X( : , T( 1:(size(T,1)-1),ielem ) )
!       Xelem(:,size(T,1))       = p_new
!       Melem = metric%get_metric_elem_newCav(T(:,ielem),m_new)
!     else
!       Xelem = X( : , T(:,ielem) )
!       Melem = metric%get_metric_elem(T(:,ielem))
!     end if
!
!     elemQual = quality_sizeShape(Xelem,Melem)
!     q_min = min(elemQual,q_min)
!
!   end do loop_elements
end function compute_minQuality_sizeShape_cavity
!
!
!
subroutine compute_quality_sizeShape_cavity(q,q_min,X,T,metric,p_new,m_new)

  !***************************************************
  !*
  !*   Computes the quality of the 3D mesh
  !*   Output: quality3D 
  !*
  !***************************************************
  use mod_debugTools, only: out_performance, deb_minQCavity
  implicit none
  
  real(rp),                intent(in) :: X(:,:) ! coords
  integer(ip),             intent(in) :: T(:,:) ! connectivity
  type(mesh_metric_type),  intent(in) :: metric
  real(rp),    optional,   intent(in) :: p_new(:)
  real(rp),    optional,   intent(in) :: m_new(:,:)
  
  real(rp), intent(out) :: q(size(T,2)), q_min
  
  integer(ip) :: ielem
  real(rp)    :: Xelem(size(X,1),size(T,1))
  real(rp)    :: Melem(size(X,1),size(X,1))
  
  real(rp) :: t0,t1
  !
  if(out_performance) call cpu_time(t0)
  !
  call setWinv( int(size(X,1),ip) )
  q_min = 1.0

  if(present(p_new)) then
    do ielem=1,size(T,2)
      Xelem(:,1:(size(T,1)-1)) = X( : , T( 1:(size(T,1)-1),ielem ) )
      Xelem(:,size(T,1))       = p_new
      Melem = metric%get_metric_elem_newCav(T(:,ielem),m_new)
      q(ielem) = quality_sizeShape(Xelem,Melem)
      q_min = min(q(ielem),q_min)
    end do
  else
    do ielem=1,size(T,2)
      Xelem = X( : , T(:,ielem) )
      Melem = metric%get_metric_elem(T(:,ielem))
      q(ielem) = quality_sizeShape(Xelem,Melem)
      q_min = min(q(ielem),q_min)
    end do
  end if

  if(out_performance) then
    call cpu_time(t1)
    deb_minQCavity = deb_minQCavity + (t1-t0)
  end if
end subroutine compute_quality_sizeShape_cavity
!
!
!
function compute_minQuality_sizeShape_subset(mesh,metric,list_elems) result(q_min)

  !***************************************************
  !*
  !*   Computes the quality of the 3D mesh
  !*   Output: quality3D 
  !*
  !***************************************************
  
  use mod_debugTools, only: out_performance, deb_minQSubset
  implicit none
  
  type(mesh_type_basic),   intent(in) :: mesh
  type(mesh_metric_type),  intent(in) :: metric
  integer(ip), intent(in) :: list_elems(:)
  
  real(rp) :: q_min
  
  real (rp)   :: elemQual
  integer(ip) :: ielem, theElem
  real (rp)   :: Xelem(mesh%ndime,mesh%mnode)
  real(rp)    :: Melem(mesh%ndime,mesh%ndime)
  
  real(rp) :: t0,t1
  !
  if(out_performance) call cpu_time(t0)
  !
  call setWinv( mesh%ndime )

  q_min = 1.0

  loop_elements: do ielem=1,size(list_elems)

    theElem = list_elems(ielem)
    Xelem = mesh%coord( : , mesh%lnods(:,theElem) )
    Melem = metric%get_metric_elem(mesh%lnods(:,theElem))
    elemQual = quality_sizeShape(Xelem,Melem)

    q_min = min(elemQual,q_min)

  end do loop_elements
  
  if(out_performance) then
    call cpu_time(t1)
    deb_minQSubset = deb_minQSubset + (t1-t0)
  end if
  
  return
end function compute_minQuality_sizeShape_subset
!
!
!
function compute_quality_sizeShape_subset(mesh,metric,list_elems) result(q_subset)

  !***************************************************
  !*
  !*   Computes the quality of the 3D mesh
  !*   Output: quality3D 
  !*
  !***************************************************
  
  use mod_debugTools, only: out_performance, deb_qualSubset
  implicit none
  
  type(mesh_type_basic),   intent(in) :: mesh
  type(mesh_metric_type),  intent(in) :: metric
  integer(ip), intent(in) :: list_elems(:)
  
  real(rp) :: q_subset(size(list_elems))
  
  integer(ip) :: ielem, theElem
  real (rp)   :: Xelem(mesh%ndime,mesh%mnode)
  real(rp)    :: Melem(mesh%ndime,mesh%ndime)
  
  real(rp) :: t0,t1
  !
  if(out_performance) call cpu_time(t0)
  !
  call setWinv(mesh%ndime)

  loop_elements: do ielem=1,size(list_elems)

    theElem = list_elems(ielem)
    Xelem = mesh%coord( : , mesh%lnods(:,theElem) )
    Melem = metric%get_metric_elem(mesh%lnods(:,theElem))
    q_subset(ielem) = quality_sizeShape(Xelem,Melem)

  end do loop_elements
  
  if(out_performance) then
    call cpu_time(t1)
    deb_qualSubset = deb_qualSubset + (t1-t0)
  end if
  
  return
end function compute_quality_sizeShape_subset
!
!
!
function compute_distortion_sizeShape_subset(mesh,metric,list_elems,delta_inp) result(eta)

  !***************************************************
  !*
  !*   Computes the quality of the 3D mesh
  !*   Output: quality3D 
  !*
  !***************************************************
  
  use mod_debugTools, only: out_performance, deb_distSubset
  implicit none
  
  type(mesh_type_basic),   intent(in) :: mesh
  type(mesh_metric_type),  intent(in) :: metric
  integer(ip),             intent(in) :: list_elems(:)
  real(rp), optional,      intent(in) :: delta_inp
  
  real(rp)    :: eta
  
  real (rp)   :: eta_elem
  integer(ip) :: ielem, theElem, pnorm
  real (rp)   :: Xelem(mesh%ndime,mesh%mnode)
  real(rp)    :: Melem(mesh%ndime,mesh%ndime)
  real(rp) :: delta
  
  real(rp) :: t0,t1
  !
  if(out_performance) call cpu_time(t0)
  !
  
  if(present(delta_inp)) then
    delta = delta_inp
  else
    delta = delta_reg_zero
  end if
  
  call setWinv(mesh%ndime)

  pnorm = 2
  
  eta = 0.0_rp

  loop_elements: do ielem=1,size(list_elems)

    theElem  = list_elems(ielem)
    Xelem    = mesh%coord( : , mesh%lnods(:,theElem) )
    Melem    = metric%get_metric_elem(mesh%lnods(:,theElem))
    eta_elem = distortion_sizeShape(Xelem,Melem)

    eta = eta + eta_elem**pnorm
    
  end do loop_elements
  
  eta = eta**(1.0_rp/pnorm)
  
  if(out_performance) then
    call cpu_time(t1)
    deb_distSubset = deb_distSubset + (t1-t0)
  end if
  
  return
end function compute_distortion_sizeShape_subset
!
!
!------- Private functions
!
subroutine setWinv(dim)
  implicit none
  integer(ip), intent(in) :: dim
  
  if(dim==dimW_saved) return
  
  if(dim==2_ip) then
    Winv => Winv_tri
  else
    Winv => Winv_tet
  end if
  
end subroutine setWinv
!
!
!
subroutine deleteWinv(dim)
  implicit none
  integer(ip), intent(in) :: dim
  
  dimW_saved = 0
  nullify(Winv)
  
  return 
end subroutine deleteWinv
!
!
!
function quality_shape(Xelem,M) result(qual)
!***
!*  Shape quality measure
!*** 
implicit none
real(rp), intent(in) :: Xelem(:,:)
real(rp), intent(in) :: M(:,:)
real(rp)   :: qual

real(rp) :: distortion 

distortion = distortion_shape(Xelem,M)
qual = 1.0/distortion

end function quality_shape
!
!
!
function quality_sizeShape(Xelem,M) result(qual)
!***
!*  Shape quality measure
!*** 
implicit none
real(rp), intent(in) :: Xelem(:,:)
real(rp), intent(in) :: M(:,:)
real(rp)   :: qual

real(rp) :: distortion 

distortion = distortion_sizeShape(Xelem,M)
qual = 1.0/distortion

end function quality_sizeShape
!
!
!
function distortion_shape(Xelem,M,delta_inp) result(distortion)
  !***
  !*  Shape distortion measure
  !*** 
  implicit none
  real(rp), intent(in)  :: Xelem(:,:)
  real(rp), intent(in)  :: M(:,:)
  real(rp), intent(in), optional :: delta_inp
  real(rp)              :: distortion

  integer(ip) :: i,dim
  real(rp)    :: A(size(Xelem,1),size(Xelem,1))
  real(rp)    :: S(size(Xelem,1),size(Xelem,1))
  real(rp)    :: nFrobSquared,det_reg
  real(rp)    :: detSign, detS, delta, detdet
  real(rp)    :: StMS(size(Xelem,1),size(Xelem,1))

  if(present(delta_inp)) then
    delta = delta_inp
  else
    delta = delta_reg_zero
  end if
  
  dim = size(A,2)
  do i=1,dim
    A(:,i) = Xelem(:,i+1) - Xelem(:,1)
  end do
  
  S = matmul(A,Winv)
  StMS = matmul(transpose(S),matmul(M,S))
  nFrobSquared = compute_trace(StMS)

  detS = determinant(S) ! just for the sign
  detSign = sign( 1.0_rp, detS)
  
  detdet = determinant(StMS)
  det_reg = 0.5_rp*(detSign*sqrt(abs( detdet )) + sqrt(detdet + 4.0_rp*delta**2_ip ) ) ! + 1.0e-15
  
  distortion = nFrobSquared /(  dim*( det_reg**(2.0_rp/dim) )  )
  
!   delta = delta_reg_zero
!
!   distortion = 0.0
!
!   dim = size(A,2)
!   do i=1,dim
!     A(:,i) = Xelem(:,i+1) - Xelem(:,1)
!   end do
!
!   S = matmul(A,Winv)
!   StMS = matmul(transpose(S),matmul(M,S))
!   nFrobSquared = compute_trace(StMS)
!
!   detS = determinant(S) ! just for the sign
!   detSign = sign( 1.0_rp, detS)
!
!   detdet = determinant(StMS)
!   det_reg = 0.5_rp*(detSign*sqrt(abs( detdet )) + sqrt(detdet + 4.0_rp*delta**2_ip ) ) ! + 1.0e-15
!
!   distortion = nFrobSquared /(  dim*( det_reg**(2.0_rp/dim) )  )

end function distortion_shape
!
!
!
function distortion_sizeShape(Xelem,M,delta_inp) result(distortion)
  !***
  !*  SizeShape distortion measure
  !*** 
  use def_master, only: kfl_paral
  implicit none
  real(rp), intent(in)  :: Xelem(:,:)
  real(rp), intent(in)  :: M(:,:)
  real(rp)              :: distortion
  real(rp), intent(in), optional :: delta_inp

  integer(ip) :: i,dim
  real(rp)    :: A(size(Xelem,1),size(Xelem,1))
  real(rp)    :: S(size(Xelem,1),size(Xelem,1))
  real(rp)    :: distortion_size, distortion_shape
  real(rp)    :: StMS(size(Xelem,1),size(Xelem,1))
  

  real(rp) :: nFrobSquared
  real(rp) :: det_reg
  real(rp) :: detdet
  real(rp) :: detSign
  real(rp) :: detS
  real(rp) :: delta 
  real(rp) :: mu_size
  real(rp) :: det_aux
  real(rp) :: detWithSign
  real(rp) :: detWithAbs
  real(rp) :: detWithDelta
  real(rp) :: mu1, mu2, det_reg_inv

  if(present(delta_inp)) then
    delta = delta_inp
  else
    delta = delta_reg_zero
  end if
  
  dim = size(A,2)
  do i=1,dim
    A(:,i) = Xelem(:,i+1) - Xelem(:,1)
  end do

  S = matmul(A,Winv)
  StMS = matmul(transpose(S),matmul(M,S))
  nFrobSquared = compute_trace(StMS)

  detS = determinant(S) ! just for the sign
  detSign = sign( 1.0_rp, detS)
  
  detdet = determinant(StMS)
  !!det_reg = 0.5_rp*(detSign*sqrt(abs( detdet )) + sqrt(detdet + 4.0_rp*delta**2_ip ) ) ! + 1.0e-15
  !det_reg = 0.5_rp*(detS + sqrt(detS*detS + 4.0_rp*delta*delta ) )
  detWithAbs   = sqrt(abs( detdet ))
  !detWithSign  = detSign*detWithAbs
  !detWithDelta = sqrt(detdet + 4.0_rp*delta*delta ) 
  !det_reg = 0.5_rp*( detWithSign + detWithDelta) 
  !det_reg = (0.5_rp*( detWithSign + detWithAbs )) + delta
  det_reg = ((detSign+1_rp)/2_rp)*detWithAbs + delta
  
  det_reg = det_reg + epsilon(1.0_rp)
  
  if(dim==0) then
    print*,kfl_paral,' ',size(Xelem,1)
    print*,kfl_paral,' ',size(Xelem,2)
    print*,kfl_paral,' ',size(A,1)
    print*,kfl_paral,' ',size(A,2)
    print*,kfl_paral,' ','det_reg: ',det_reg
    print*,kfl_paral,' ','dim: ',dim
    print*,kfl_paral,' ','op x:   ',(2.0_rp/real(dim,rp))
    print*,kfl_paral,' ','op xx:  ',( det_reg**(2.0_rp/real(dim,rp)) )
    print*,kfl_paral,' ','op xxx: ',dim*( det_reg**(2.0_rp/real(dim,rp)) )
    
    call runend('eee dimension zero!')
  end if
  !det_aux = real(dim,rp)*(  det_reg**(2.0_rp/real(dim,rp))  ) 
  det_aux =  det_reg*det_reg
  det_aux =  det_aux**(1.0_rp/real(dim,rp))
  det_aux = real(dim,rp)*det_aux
  distortion_shape = nFrobSquared / det_aux
  
  mu_size = min(det_reg, 1.0_rp/det_reg)
  distortion_size = 1.0_rp/mu_size
  !mu_size =  exp(1.0_rp)*( det_reg*exp(-det_reg) + exp(-1.0_rp/det_reg)/det_reg )/2.0_rp
!   mu1 = det_reg*exp(-det_reg)
!   if(det_reg<epsilon(1.0_rp)*1000) then
!     print*,'mu1: ', mu1
!     print*,'det_reg: ', det_reg
!     print*,'1.0_rp/det_reg: ',1.0_rp/det_reg
!   end if
!   det_reg_inv = 1.0_rp/det_reg
!   mu2 = det_reg_inv/exp(det_reg_inv)
!   mu_size = exp(1.0_rp)*( mu1 + mu2 )/2.0_rp
!   distortion_size = 1.0_rp/mu_size

  !distortion = sqrt(distortion_size*distortion_shape)
  distortion = (distortion_size**exp_size)*(distortion_shape**exp_shape)
  distortion = distortion**(1.0_rp/(exp_size+exp_shape))
!   distortion = ( ((distortion_size**1)+(distortion_shape**1))/2.0 )
!   distortion = sqrt( ((distortion_size**2)+(distortion_shape**2))/2.0 )

  !if( isnan(distortion) ) then
!   if( distortion/=distortion ) then ! a nan is not equal to a nan
!     print*,"A: ",A
!     print*,"determinant(StMS): ",determinant(StMS)
!     print*,"sqrt(determinant(StMS)): ",sqrt(determinant(StMS))
!     print*,"det: ",det
!     print*,"det_reg: ",det_reg
!     print*,"distortion_shape: ",distortion_shape
!     print*,"distortion_size: ",distortion_size
!     print*,"distortion: ",distortion
!     print*,'in distortion size-shape'
!     call runend('Inifinite distortion value... nan nan')
!   end if

end function distortion_sizeShape
!
!
!
function compute_nFrobSquared(S) result(nFrobS) 
  implicit none 
  
  real(rp),intent(in) :: S(:,:)
  real(rp) :: nFrobS
  integer(ip) :: i

  nFrobS = 0.0
  do i=1,size(S,2)
    nFrobS = nFrobS + sum((S(:,i)**2))
  end do
end function compute_nFrobSquared
!
!
!
function compute_trace(S) result(trace) 
  implicit none 
  
  real(rp),intent(in) :: S(:,:)
  real(rp) :: trace
  integer(ip) :: i

  trace = 0.0
  do i=1,size(S,2)
    trace = trace + S(i,i)
  end do
end function compute_trace
!
!
!
subroutine print_q_stats(quality)
  use mod_debugTools, only: out_debug_text
  implicit none
  real(rp),intent(in)    :: quality(:)

  real(rp)::meanQ,stDev,maxQ,minQ

  meanQ = (sum(quality)/float(size(quality)) )
  stDev = sqrt( (sum((quality-meanQ)**2))/(1.0_rp*size(quality)) )
  maxQ = maxval(quality)
  minQ = minval(quality)

  if(out_debug_text) then
   write(*,1) maxQ,minQ,meanQ,stDev
   1  format(/,&
   5x,'|-------Quality statistics--------|',/, &
   5x,'| Maximum quality value: ',f8.2,' |',/, &
   5x,'| Minimum quality value: ',f8.2,' |',/, &
   5x,'| Mean    quality value: ',f8.2,' |',/, &
   5x,'| Standard  deviation  : ',f8.2,' |',/, &
   /) !5x,'| Num inverted elements: ',i8,' |',
  end if

return
end subroutine print_q_stats
!
!
!
subroutine do_print_q_stats(quality)
  use mod_debugTools, only: out_debug_text
  implicit none
  real(rp),intent(in)    :: quality(:)

  real(rp)::meanQ,stDev,maxQ,minQ
  integer(ip) :: numInv

  meanQ = (sum(quality)/float(size(quality)) )
  stDev = sqrt( (sum((quality-meanQ)**2))/(1.0_rp*size(quality)) )
  maxQ = maxval(quality)
  minQ = minval(quality)
  
  numInv = count(quality<1.0E-6)

   write(*,1) maxQ,minQ,meanQ,stDev,numInv
   1  format(/,&
   5x,'|-------Quality statistics--------|',/, &
   5x,'| Maximum quality value: ',f8.2,' |',/, &
   5x,'| Minimum quality value: ',f8.2,' |',/, &
   5x,'| Mean    quality value: ',f8.2,' |',/, &
   5x,'| Standard  deviation  : ',f8.2,' |',/, &
   5x,'| Num inverted elements: ',i8,' |', &
   /) 

return
end subroutine do_print_q_stats
!
!
!
function determinant(S) result(detS)
  implicit none 

  real(rp), intent(in) :: S(:,:)
  real(rp)  :: detS

  if (size(S,1)==3) then
    detS = &
      S(1,1)*S(2,2)*S(3,3)  &
    - S(1,1)*S(2,3)*S(3,2)  &
    - S(1,2)*S(2,1)*S(3,3)  &
    + S(1,2)*S(2,3)*S(3,1)  &
    + S(1,3)*S(2,1)*S(3,2)  &
    - S(1,3)*S(2,2)*S(3,1)
  elseif (size(S,1)==2) then
    detS = S(1,1)*S(2,2) - S(1,2)*S(2,1)
  else
    call runend('Not implemented determinant')
  end if
end function determinant
!
!
!
function isPointInTet(p,coordTet,out_values) result(isInside)
  implicit none 

  real(rp), intent(in) :: p(:), coordTet(:,:)
  logical(lg), intent(in), optional :: out_values
  integer(ip)          :: isInside
  real(rp)    :: detS, S(3,3)
  integer(ip) :: iface,inode
  
  logical(lg) :: isInside_clear, isInside_doubt
  
  integer(ip), parameter :: list_faces_TET04(3,4)  = reshape ( (/ 1,3,2,   2,3,4,   1,2,4,   3,1,4 /), (/3,4 /) )
  
  isInside_clear = .true.
  isInside_doubt = .true.
  do iface=1,4
    
    !S = coordTet(:,list_faces_TET04(:,iface))-p
    do inode=1,3
      S(:,inode) = coordTet(:,list_faces_TET04(inode,iface))-p
    end do
    detS = determinant(S)

    isInside_clear = isInside_clear.and.(detS>=0.0_rp)
    isInside_doubt = isInside_doubt.and.(detS>=-1.0E-12_rp)

  end do

  isInside = 0_ip
  if(isInside_doubt) isInside=2_ip
  if(isInside_clear) isInside=1_ip
  
  if(present(out_values)) then
    if(out_values.and.(isInside>0)) then
      do iface=1,4
        do inode=1,3
          S(:,inode) = coordTet(:,list_faces_TET04(inode,iface))-p
        end do
        detS = determinant(S)
        print*,'  Sign face: ',iface,': ',detS
      end do
    end if
  end if
end function isPointInTet
!
!
!
pure function matinv2(A) result(B)
  !! Performs a direct calculation of the inverse of a 2×2 matrix.
  real(rp), intent(in) :: A(2,2)   !! Matrix
  real(rp)             :: B(2,2)   !! Inverse matrix
  real(rp)             :: detinv

  ! Calculate the inverse determinant of the matrix
  detinv = 1/(A(1,1)*A(2,2) - A(1,2)*A(2,1))

  ! Calculate the inverse of the matrix
  B(1,1) = +detinv * A(2,2)
  B(2,1) = -detinv * A(2,1)
  B(1,2) = -detinv * A(1,2)
  B(2,2) = +detinv * A(1,1)
end function

pure function matinv3(A) result(B)
  !! Performs a direct calculation of the inverse of a 3×3 matrix.
  real(rp), intent(in) :: A(3,3)   !! Matrix
  real(rp)             :: B(3,3)   !! Inverse matrix
  real(rp)             :: detinv

  ! Calculate the inverse determinant of the matrix
  detinv = 1/(A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2)&
            - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1)&
            + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1))

  ! Calculate the inverse of the matrix
  B(1,1) = +detinv * (A(2,2)*A(3,3) - A(2,3)*A(3,2))
  B(2,1) = -detinv * (A(2,1)*A(3,3) - A(2,3)*A(3,1))
  B(3,1) = +detinv * (A(2,1)*A(3,2) - A(2,2)*A(3,1))
  B(1,2) = -detinv * (A(1,2)*A(3,3) - A(1,3)*A(3,2))
  B(2,2) = +detinv * (A(1,1)*A(3,3) - A(1,3)*A(3,1))
  B(3,2) = -detinv * (A(1,1)*A(3,2) - A(1,2)*A(3,1))
  B(1,3) = +detinv * (A(1,2)*A(2,3) - A(1,3)*A(2,2))
  B(2,3) = -detinv * (A(1,1)*A(2,3) - A(1,3)*A(2,1))
  B(3,3) = +detinv * (A(1,1)*A(2,2) - A(1,2)*A(2,1))
end function
!
!
!
END MODULE mod_quality
!
!
!> @}


  
!   if(allocated(Winv)) then
!     if ( size(Winv,1).ne.mesh%ndime ) then
!       deallocate(Winv)
!     end if
!   end if
!   if(.not.allocated(Winv)) then
!     allocate(Winv(mesh%ndime,mesh%ndime))
!     if(mesh%ndime==2) then
!       Winv = Winv_tri
!     else
!       Winv = Winv_tet
!     end if
!   end if
