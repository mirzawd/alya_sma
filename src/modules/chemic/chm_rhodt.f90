!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine chm_rhodt(pnode,pgaus,porde,lnods_loc,gpsha,gpvol,gpden,dt_rho_chm_loc)
  !------------------------------------------------------------------------
  ! NAME
  !    chm_rhodt
  ! DESCRIPTION
  !    Projection of rho/dt for explicit time step
  ! USES
  ! USED BY
  !    chm_elmcfi
  !***
  !------------------------------------------------------------------------

  use def_master,      only : dtinv
  use def_kintyp,      only : ip, rp
  implicit none

  integer(ip), intent(in)  :: pnode
  integer(ip), intent(in)  :: pgaus
  integer(ip), intent(in)  :: porde
  integer(ip), intent(in)  :: lnods_loc(pnode)
  real(rp),    intent(in)  :: gpsha(pnode,pgaus)
  real(rp),    intent(in)  :: gpvol(pgaus)
  real(rp),    intent(in)  :: gpden(pgaus)
  real(rp),    intent(out) :: dt_rho_chm_loc(*)
  integer(ip) :: inode,jnode,ipoin,igaus
  real(rp)     :: eldtrho(pnode), fact
  real(rp)                :: elmat(pnode,pnode)
  real(rp)                :: trace, elmass


  if( porde == 1 ) then
     !
     ! Element assembly
     !
     eldtrho = 0.0_rp
     do igaus = 1,pgaus
        fact = gpvol(igaus) / ( gpden(igaus)  * dtinv )
        do inode = 1,pnode
           eldtrho(inode) = eldtrho(inode) + gpsha(inode,igaus) * fact
        end do
     end do
     !
     ! Nodal projection
     !
     do inode = 1,pnode
        ipoin = lnods_loc(inode)
        dt_rho_chm_loc(ipoin) = dt_rho_chm_loc(ipoin) + eldtrho(inode)
     end do
  else
     !
     ! Element assembly
     !
     eldtrho = 0.0_rp


    ! T = 0.0_rp
    ! d = 0.0_rp
    ! do igaus = 1,pgaus
    !    fact = gpvol(igaus) / ( gpden(igaus) * dtinv )
    !    T = T + fact
    !    do inode = 1,pnode
    !       eldtrho(inode) = eldtrho(inode) + gpsha(inode,igaus)**2 * fact
    !    end do
    ! end do
    ! do inode = 1,pnode
    !    d = d + eldtrho(inode)
    ! end do
    ! do inode = 1,pnode
    !    eldtrho(inode) = eldtrho(inode) * T / d
    ! end do

     do inode=1,pnode
        do jnode=1,pnode
           elmat(inode,jnode)=0.0_rp
        end do
     end do

     do igaus=1,pgaus
        do inode=1,pnode
           fact=gpvol(igaus)*gpsha(inode,igaus)/(gpden(igaus)*dtinv)
           do jnode=1,pnode
              elmat(inode,jnode)=elmat(inode,jnode) +fact*gpsha(jnode,igaus)
           end do
        end do
     end do

     trace  = 0.0_rp
     elmass = 0.0_rp
     do inode = 1,pnode
        trace = trace + elmat(inode,inode)
        do jnode = 1,pnode
           elmass = elmass + elmat(inode,jnode)
        end do
     end do

     !
     ! Nodal projection
     !
     do inode = 1,pnode
        ipoin = lnods_loc(inode)
        dt_rho_chm_loc(ipoin) = dt_rho_chm_loc(ipoin) + elmat(inode,inode)*(elmass/trace)
     end do
  end if


end subroutine chm_rhodt
