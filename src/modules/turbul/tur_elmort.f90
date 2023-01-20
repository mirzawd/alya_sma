!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tur_elmort(&
     gpres, gpcon, gprec,grtup,gprhs, ielem,pgaus,pnode, gpsha,gpvol)
  !------------------------------------------------------------------------
  !****f* Turbul/tur_elmort
  ! NAME 
  !    tur_elmort
  ! DESCRIPTION
  !    This subroutine assemble the projected residuals
  ! USES
  ! USED BY
  !    tur_elmop2
  !***
  !------------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp
  use def_domain, only       :  ndime,lnods
  use def_turbul, only       :  kfl_ortho_tur, tupr2_tur, tuprr_tur, kfl_shock_tur, tupgr_tur
  implicit none
  integer(ip), intent(in)    :: ielem,pgaus,pnode
  real(rp),    intent(in)    :: gpres(pgaus), gpcon(pgaus), gprec(pgaus), gprhs(pgaus)
  real(rp),    intent(in)    :: gpsha(pnode,pgaus)
  real(rp),    intent(in)    :: gpvol(pgaus), grtup(ndime,pgaus)
  real(rp)                   :: gptup(pgaus), gprep(pgaus) ! projection of turbul
  real(rp)                   :: eltup(pnode), eltu2(pnode), grpro(ndime, pgaus), elgrp(ndime, pnode)
  integer(ip)                :: igaus,inode

  if( kfl_ortho_tur <= 0 .and.kfl_shock_tur.eq.0) then ! no orthogonal projection

     return

  else if( kfl_ortho_tur == 1 ) then ! FULL OSS

     !-------------------------------------------------------------------
     !
     ! GPTUP: Projection of turbulence residual
     !
     !-------------------------------------------------------------------

     do igaus =1, pgaus
        gptup(igaus) =  - gpres(igaus) * gpvol(igaus)
     end do   

  else if( kfl_ortho_tur == 2 ) then ! split oss

     !-------------------------------------------------------------------
     !
     ! GPTUP: Projection of rho*(a.grad)u
     !
     !-------------------------------------------------------------------

     do igaus =1, pgaus
        gptup(igaus) =  gpcon(igaus) * gpvol(igaus)
        gprep(igaus) =  gprec(igaus) * gpvol(igaus)
     end do
     ! Assemble RHS
     do inode = 1,pnode
        eltu2(inode) = 0.0_rp
     end do
     do igaus = 1,pgaus
        do inode = 1,pnode
           eltu2(inode) = eltu2(inode) &
                + gprep(igaus) * gpsha(inode,igaus) 
        end do
     end do
     call assrhs(1_ip,pnode,lnods(1,ielem),eltu2,tuprr_tur)
  end if
  if (kfl_shock_tur.ne.0) then
     !
     ! GPTUP : Projection of turbulence gradient
     !
     do igaus =1, pgaus
        grpro(1:ndime,igaus) =  grtup(1:ndime, igaus) * gpvol(igaus)
     end do

     elgrp(1:ndime, 1:pnode) =0.0_rp
     do igaus =1, pgaus
        do inode=1, pnode
           elgrp(1:ndime, inode)=  elgrp(1:ndime, inode) &
                + grpro(1:ndime, igaus)*gpsha(inode, igaus)
        end do

     end do
     ! Assemble RHS
     call assrhs(ndime, pnode, lnods(1, ielem), elgrp, tupgr_tur)
  end if

  !----------------------------------------------------------------------
  !
  ! Assemble RHS
  !
  !----------------------------------------------------------------------
  if (kfl_shock_tur==0) then
     do inode = 1,pnode
        eltup(inode) = 0.0_rp
     end do
     do igaus = 1,pgaus
        do inode = 1,pnode
           eltup(inode) = eltup(inode) &
                + gptup(igaus) * gpsha(inode,igaus) 
        end do
     end do     
     call assrhs(1_ip,pnode,lnods(1,ielem),eltup,tupr2_tur)
  end if

end subroutine tur_elmort
