!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup Partis
!> @{
!> @file    pts_stk_local.f90
!> @author  Laura Nicolaou
!> @brief   This routine computes the instantaneous and effective Stokes
!>          numbers based on local properties of the flow field
!> @details
!>          Characteristic fluid timescale given by max. eigenvalue
!>          of velocity gradient tensor (Nicolaou & Zaki, JAS, 2016)
!>
!> @}
!------------------------------------------------------------------------

subroutine pts_stk_local(ielem, tau, t_inject, Stk_inst, Stk_eff)

  use def_kintyp, only          : ip,rp
  use def_master, only          : dtime, cutim, advec
  use def_domain, only          : ndime, nnode, mnode, elmar, lnods, ltype, coord
  !use def_partis
  implicit none

  integer(ip), intent(in)    :: ielem
  real(rp),    intent(in)    :: tau, t_inject
  real(rp),    intent(out)   :: Stk_inst
  real(rp),    intent(inout) :: Stk_eff
  integer(ip)                :: idime, jdime, inode, jnode      !counters
  integer(ip)                :: pelty,pnode,ipoin
  real(rp)                   :: elvel(ndime,mnode),elcod(ndime,mnode)
  real(rp)                   :: detjm,xjaci(9),xjacm(9),cartc(ndime,mnode)
  real(rp)                   :: gvelo(ndime,ndime),lambda_max

  pelty = ltype(ielem)
  pnode = nnode(pelty)
  !
  ! Gather
  !
  do inode = 1,pnode 
     ipoin = lnods(inode,ielem)
     elcod(1:ndime,inode) = coord(1:ndime,ipoin)     
     elvel(1:ndime,inode) = advec(1:ndime,ipoin,1)
  end do
  !
  ! Cartesian deriative at center of gravity
  ! 
  call elmder(&
       pnode,ndime,elmar(pelty)%dercg,&
       elcod,cartc,detjm,xjacm,xjaci)
  !
  ! Velocity gradient tensor GVELO(i,j) = duj/dxi
  !
  gvelo = 0.0_rp
  do jnode = 1,pnode 
     do idime = 1,ndime 
        do jdime = 1,ndime     
           gvelo(idime,jdime) = gvelo(idime,jdime) & 
                + cartc(idime,jnode)*elvel(jdime,jnode) 
        end do
     end do
  end do

  if( ndime == 3 ) then
     !
     ! Instantaneous Stokes number
     !
     call veloc_grad_tensor_lambda_max(gvelo, lambda_max)
     Stk_inst = tau*lambda_max
     !
     ! Effective Stokes number (running average)
     !
     if( cutim-t_inject /= 0.0_rp )then
        Stk_eff = (Stk_eff*(cutim-dtime-t_inject) + Stk_inst*dtime)/(cutim-t_inject)
     else
        Stk_eff = 0.0_rp
     end if
  else
     Stk_inst = 0.0_rp
     Stk_eff  = 0.0_rp
  end if

end subroutine pts_stk_local

