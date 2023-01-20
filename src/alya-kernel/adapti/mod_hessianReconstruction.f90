!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Adaptivity
!> @{
!> @file    mod_hessianReconstruction.f90
!> @author  abel.gargallo
!> @date    2021-03-31
!> @brief   mod_hessianReconstruction to compute metric from solution following the ideas from Alauzet&Loiselle
!> @details mod_hessianReconstruction
!>
!>          To add further details
!>
!>
!-----------------------------------------------------------------------

MODULE mod_hessianReconstruction
!***************************************************************
!*
!*  Module for metric computation and operations
!*
!***************************************************************
use def_kintyp_basic,       only: ip,rp,lg
use def_adapt,              only: memor_adapt
use mod_memory,             only: memory_alloca, memory_deallo


implicit none

private

public :: compute_hessian_from_scalar,deallo_hessian_from_scalar
!
!
!
CONTAINS
!
!
!
subroutine compute_hessian_from_scalar(u,ndime,hessu)! result(hessu)
  use mod_gradie, only: gradie, projec_hessian
  implicit none
  !
  real(rp), pointer, intent(in)    :: u(:)
  integer(ip),       intent(in)    :: ndime
  real(rp), pointer, intent(inout) :: hessu(:,:,:)
  !
  nullify(hessu)
  !
  if(associated(u)) then
    call projec_hessian(u,hessu)
  end if
  
end subroutine compute_hessian_from_scalar
!
!
!
subroutine deallo_hessian_from_scalar(hessu)

  use def_domain, only: memor_dom
  implicit none
  
  real(rp), pointer, intent(inout) :: hessu(:,:,:)
  
  call memory_deallo(memor_dom,'HESSI','mod_projec',hessu)
  
end subroutine deallo_hessian_from_scalar

!
!
!
END MODULE mod_hessianReconstruction

!> @}


!   real(rp), pointer :: gradu(:,:)
!   integer(ip)       :: idime, inode
!   integer(ip)       :: npoin
!
!   !- reconstruct smooth gradient  (higher-order approx)
!   !- reconstruct smooth hessian   (higher-order approx)

!   npoin = size(u)
!
!   nullify(gradu)
!   call memory_alloca(memor_adapt,'GRADU','compute_metric_from_sol_viaHessian',gradu,ndime,npoin)
!   call gradie(u,gradu)
!
!   print*,'hessu'
!   nullify(hessu)
!   call memory_alloca(memor_adapt,'HESSU','compute_metric_from_sol_viaHessian',hessu,ndime,ndime,npoin)
!   print*,'gradie per dim'
!   do idime=1,ndime
!     call gradie(gradu(idime,:),hessu(idime,:,:))
!   end do
!   call memory_deallo(memor_adapt,'GRADU','compute_metric_from_sol_viaHessian',gradu)
!
!   !print*,'SHOULD I SYMMETRIZE HESSIAN HERE????...YES... AT LEAST UNTIL I USE ANOTHER HESSIAN COMPUTATION'
!   !hessu = (hessu + TRANSPOSE(hessu))/2.0 ! ensure hessian is symmetric
!   print*,'symetrize?'
!   do inode=1,npoin
!     !hessu(:,:,inode) = (hessu(:,:,inode) + transpose(hessu(:,:,inode)))/2.0_rp
!     hessu(:,:,inode) = 0.0_rp
!     hessu(1,1,inode) = 1.0_rp
!     hessu(2,2,inode) = 1.0_rp
!   end do
!
!   print*,'end loop and deallo gradu'
!
! 
