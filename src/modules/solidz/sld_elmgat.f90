!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    sld_elmgat.f90
!> @author  Guillaume Houzeaux
!> @date    07/05/2013
!> @brief   Element gather
!> @details Gather some global arrays to elemental arrays
!> @}
!-----------------------------------------------------------------------

subroutine sld_elmgat(pnode,lnods,eldis,elddi,eldip,elvel,elacc,elcod, &
     eldix,elvex,elmof,elrst)

  use def_kintyp, only     :  ip,rp
  use def_master, only     :  displ, TIME_N, ITER_K
  use def_domain, only     :  ndime, coord
  use def_solidz, only     :  ncomp_sld, nvoig_sld
  use def_solidz, only     :  veloc_sld, accel_sld, ddisp_sld
  use def_solidz, only     :  kfl_xfeme_sld, dxfem_sld, vxfem_sld
  use def_solidz, only     :  fiemo_sld, kfl_moduf_sld
  use def_solidz, only     :  kfl_restr_sld, restr_sld
  use def_solidz, only     :  kfl_ninex_sld, disep_sld

  implicit none

  integer(ip), intent(in)  :: pnode                         !< Number of element nodes
  integer(ip), intent(in)  :: lnods(pnode)                  !< Element node connectivity
  real(rp),    intent(out) :: elcod(ndime,pnode)            !< Element node coordinates
  real(rp),    intent(out) :: eldis(ndime,pnode,ncomp_sld)  !< Element node displacement
  real(rp),    intent(out) :: elddi(ndime,pnode,ncomp_sld)  !< Element node increment of displacement
  real(rp),    intent(out) :: eldip(ndime,pnode)            !< Element node displacement perturbation
  real(rp),    intent(out) :: elvel(ndime,pnode,ncomp_sld)  !< Element node velocity
  real(rp),    intent(out) :: elacc(ndime,pnode,ncomp_sld)  !< Element node acceleration
  real(rp),    intent(out) :: eldix(ndime,pnode,3)          !< Element XFEM node displacement
  real(rp),    intent(out) :: elvex(ndime,pnode,3)          !< Element XFEM node velocity
  real(rp),    intent(out) :: elmof(pnode)                  !< Element node modulator field
  real(rp)                 :: elrst(6,pnode)                !< Element node residual stresses

  integer(ip)              :: inode,ipoin,ivoig
  !
  ! Current global values
  !
  do inode = 1,pnode
     ipoin = lnods(inode)
     elcod(1:ndime,inode)        = coord(1:ndime,ipoin)
     eldis(1:ndime,inode,ITER_K) = displ(1:ndime,ipoin,ITER_K)
     eldis(1:ndime,inode,TIME_N) = displ(1:ndime,ipoin,TIME_N)
     elddi(1:ndime,inode,ITER_K) = ddisp_sld(1:ndime,ipoin,ITER_K)
     elddi(1:ndime,inode,TIME_N) = ddisp_sld(1:ndime,ipoin,TIME_N)
     elvel(1:ndime,inode,ITER_K) = veloc_sld(1:ndime,ipoin,ITER_K)
     elacc(1:ndime,inode,ITER_K) = accel_sld(1:ndime,ipoin,ITER_K)
     elacc(1:ndime,inode,TIME_N) = accel_sld(1:ndime,ipoin,TIME_N)
  end do
  !
  ! Inexact Newton
  !
  eldip = 0.0_rp
  if ( kfl_ninex_sld == 1 ) then
     do inode = 1,pnode
        ipoin = lnods(inode)
        eldip(1:ndime,inode) = disep_sld(1:ndime,ipoin)
     end do
  end if
  !
  ! XFEM
  !
  eldix = 0.0_rp
  elvex = 0.0_rp
  if ( kfl_xfeme_sld == 1 ) then
     do inode = 1,pnode
        ipoin = lnods(inode)
        eldix(1:ndime,inode,1) = dxfem_sld(1:ndime,ipoin,1)
        eldix(1:ndime,inode,2) = dxfem_sld(1:ndime,ipoin,3)
        eldix(1:ndime,inode,3) = dxfem_sld(1:ndime,ipoin,4)
        elvex(1:ndime,inode,1) = vxfem_sld(1:ndime,ipoin,1)
     end do
  end if
  !
  ! Rigidity modulator field
  !
  elmof = 0.0_rp
  if ( kfl_moduf_sld(1) > 0 ) then
     do inode = 1, pnode
        ipoin = lnods(inode)
        elmof(inode) = fiemo_sld(1,ipoin)
     end do
  end if
  !
  ! Residual stresses
  !
  elrst = 0.0_rp
  if ( kfl_restr_sld < 0 ) then
     do ivoig= 1, nvoig_sld
        do inode = 1, pnode
           ipoin = lnods(inode)
           elrst(ivoig,inode) = restr_sld(ivoig,ipoin)
        end do
     end do
  end if

end subroutine sld_elmgat
