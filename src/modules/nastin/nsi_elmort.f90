!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> 
!> @author  guillaume
!> @date    2022-09-14
!> @brief   This subroutine assemble the projected residuals
!> @details This subroutine assemble the projected residuals
!> 
!-----------------------------------------------------------------------

subroutine nsi_elmort(&
     ielem,pgaus,pnode,ndofn,elvel,elpre,rmom1,rmom2,gprhs,&
     gpgpr,gpsha,gpvol,gpden,gpadv,gpcar,gpsp1,gpsp2,gpst1)

  use def_kintyp, only       :  ip,rp
  use def_elmtyp, only       :  ELEXT
  use def_domain, only       :  ndime,mnode,lnods,lelch
  use def_nastin, only       :  kfl_stabi_nsi,vepr2_nsi,prpr2_nsi, kfl_regim_nsi
  use def_nastin, only       :  grpr2_nsi,kfl_rmom2_nsi
  implicit none
  integer(ip), intent(in)    :: ielem,pgaus,pnode,ndofn
  real(rp),    intent(in)    :: elvel(ndime,pnode)
  real(rp),    intent(in)    :: elpre(pnode)
  real(rp),    intent(inout) :: rmom1(pnode,pgaus)
  real(rp),    intent(in)    :: rmom2(ndime,ndime,pnode,pgaus)
  real(rp),    intent(in)    :: gprhs(ndofn,pgaus)
  real(rp),    intent(in)    :: gpgpr(ndime,pgaus)
  real(rp),    intent(in)    :: gpsha(pnode,pgaus)
  real(rp),    intent(in)    :: gpvol(pgaus)
  real(rp),    intent(in)    :: gpden(pgaus)
  real(rp),    intent(in)    :: gpadv(ndime,pgaus)
  real(rp),    intent(in)    :: gpcar(ndime,mnode,pgaus)
  real(rp),    intent(in)    :: gpsp1(pgaus)
  real(rp),    intent(in)    :: gpsp2(pgaus)
  real(rp),    intent(in)    :: gpst1(pgaus)
  integer(ip)                :: idime,igaus,inode,jdime
  real(rp)                   :: gpvep(ndime,pgaus)
  real(rp)                   :: gpprp(pgaus)
  real(rp)                   :: gpgrp(ndime,pgaus)
  real(rp)                   :: elvep(ndime,pnode)
  real(rp)                   :: elprp(pnode)
  real(rp)                   :: elgrp(ndime,pnode)
  real(rp)                   :: fact1,dummr, taupr(pgaus)

  if( kfl_stabi_nsi == 0 ) then

     return

  else if( kfl_stabi_nsi == 1 ) then

     !-------------------------------------------------------------------
     !
     ! Full OSS
     !
     ! GPVEP: Projection of momentum residual
     ! = - tau1 * R(u)
     ! =   tau1 * [ du/dt + rho*(a.grad) u + grad(p) - div[2*mu*eps(u)] - f ]
     !
     !-------------------------------------------------------------------
     
     if (kfl_regim_nsi==3) then
        do igaus = 1,pgaus
           taupr(igaus) = 1.0_rp
        end do
     else
        do igaus = 1,pgaus
          taupr(igaus)=  gpst1(igaus)
        end do
     end if
     do igaus = 1,pgaus
        do idime = 1,ndime
           gpvep(idime,igaus) = gprhs(idime,igaus) - gpgpr(idime,igaus)
           do inode=1,pnode
              gpvep(idime,igaus) = gpvep(idime,igaus)&
                   - rmom1(inode,igaus) * elvel(idime,inode)
           end do                      
           gpvep(idime,igaus) = - taupr(igaus) * gpvep(idime,igaus) * gpvol(igaus)
        end do
     end do

     if( kfl_rmom2_nsi /= 0 ) then
        do igaus = 1,pgaus
           fact1 =  taupr(igaus)* gpvol(igaus)
           do idime = 1,ndime
              do inode = 1,pnode
                 do jdime = 1,ndime
                    gpvep(idime,igaus) = gpvep(idime,igaus) & 
                         + rmom2(idime,jdime,inode,igaus) * elvel(jdime,inode)*fact1
                 end do
              end do
           end do
        end do
     end if

  else if( kfl_stabi_nsi == 2 ) then

     !-------------------------------------------------------------------
     !
     ! Split OSS
     !
     ! GPVEP: Projection of tau1' * rho * (a.grad) u
     ! GPGRP: Projection of tau1' * [ grad(p) - f ]
     !
     !-------------------------------------------------------------------

     do igaus = 1,pgaus

        do idime = 1,ndime
           gpvep(idime,igaus) = 0.0_rp
           do inode = 1,pnode
              gpvep(idime,igaus) = gpvep(idime,igaus)&
                   + rmom1(inode,igaus) * elvel(idime,inode)
           end do
           fact1              = gpsp1(igaus) * gpvol(igaus)
           gpvep(idime,igaus) = fact1 * gpvep(idime,igaus) 
           gpgrp(idime,igaus) = fact1 * ( gpgpr(idime,igaus) - gprhs(idime,igaus) ) 
        end do

     end do

  end if

  !----------------------------------------------------------------------
  !
  ! Projection of tau2 * div(u)
  !
  !----------------------------------------------------------------------

  do igaus = 1,pgaus
     gpprp(igaus) = 0.0_rp
     do idime = 1,ndime
        do inode=1,pnode
           gpprp(igaus) = gpprp(igaus) &
                + gpcar(idime,inode,igaus) * elvel(idime,inode)
        end do
     end do
     gpprp(igaus) = gpprp(igaus) * gpsp2(igaus) * gpvol(igaus)
  end do

  !----------------------------------------------------------------------
  !
  ! Assemble RHS
  !
  !----------------------------------------------------------------------

  if( ndime == 2 ) then

     do inode = 1,pnode
        elvep(1,inode) = 0.0_rp
        elvep(2,inode) = 0.0_rp
        elprp(inode)   = 0.0_rp
     end do
     do igaus = 1,pgaus
        do inode = 1,pnode
           elvep(1,inode) = elvep(1,inode) + gpvep(1,igaus) * gpsha(inode,igaus)
           elvep(2,inode) = elvep(2,inode) + gpvep(2,igaus) * gpsha(inode,igaus)
           elprp(inode)   = elprp(inode)   + gpprp(igaus)   * gpsha(inode,igaus)
        end do
     end do
     if( kfl_stabi_nsi == 2 ) then
        do inode = 1,pnode
           elgrp(1,inode) = 0.0_rp
           elgrp(2,inode) = 0.0_rp
        end do
        do igaus = 1,pgaus
           do inode = 1,pnode
              elgrp(1,inode) = elgrp(1,inode) + gpgrp(1,igaus) * gpsha(inode,igaus)
              elgrp(2,inode) = elgrp(2,inode) + gpgrp(2,igaus) * gpsha(inode,igaus)
           end do
        end do
     end if

  else

     do inode = 1,pnode
        elvep(1,inode) = 0.0_rp
        elvep(2,inode) = 0.0_rp
        elvep(3,inode) = 0.0_rp
        elprp(inode)   = 0.0_rp
     end do
     do igaus = 1,pgaus
        do inode = 1,pnode
           elvep(1,inode) = elvep(1,inode) + gpvep(1,igaus) * gpsha(inode,igaus)
           elvep(2,inode) = elvep(2,inode) + gpvep(2,igaus) * gpsha(inode,igaus)
           elvep(3,inode) = elvep(3,inode) + gpvep(3,igaus) * gpsha(inode,igaus)
           elprp(inode)   = elprp(inode)   + gpprp(igaus)   * gpsha(inode,igaus)
        end do
     end do
     if( kfl_stabi_nsi == 2 ) then
        do inode = 1,pnode
           elgrp(1,inode) = 0.0_rp
           elgrp(2,inode) = 0.0_rp
           elgrp(3,inode) = 0.0_rp
        end do
        do igaus = 1,pgaus
           do inode = 1,pnode
              elgrp(1,inode) = elgrp(1,inode) + gpgrp(1,igaus) * gpsha(inode,igaus)
              elgrp(2,inode) = elgrp(2,inode) + gpgrp(2,igaus) * gpsha(inode,igaus)
              elgrp(3,inode) = elgrp(3,inode) + gpgrp(3,igaus) * gpsha(inode,igaus)
           end do
        end do
     end if
  end if
  !
  ! Extension elements (Dodeme)
  !
  if( lelch(ielem) == ELEXT ) then
     call elmext(&
          3_ip,ndime,pnode,dummr,dummr,dummr,dummr,dummr,&
          dummr,elvep,dummr)
     call elmext(&
          3_ip, 1_ip,pnode,dummr,dummr,dummr,dummr,dummr,&
          dummr,elprp,dummr)
     if( kfl_stabi_nsi == 2 ) then
        call elmext(&
             3_ip,ndime,pnode,dummr,dummr,dummr,dummr,dummr,&
             dummr,elgrp,dummr)        
     end if
  end if

  call assrhs(ndime,pnode,lnods(1,ielem),elvep,vepr2_nsi)
  call assrhs(1_ip, pnode,lnods(1,ielem),elprp,prpr2_nsi)
  if( kfl_stabi_nsi == 2 ) then
     call assrhs(ndime,pnode,lnods(1,ielem),elgrp,grpr2_nsi) 
  end if

end subroutine nsi_elmort
