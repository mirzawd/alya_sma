!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tem_rhocpdt(ielem,pnode,pgaus,porde,lnods_loc,gpsha,gpvol,gpden,dt_rho_cp_tem_loc)
  use def_parame
  use def_elmtyp
  use def_master
  use def_domain
  use def_kermod
  use mod_ker_proper 
  use def_temper
  implicit none

  integer(ip), intent(in)  :: ielem
  integer(ip), intent(in)  :: pnode
  integer(ip), intent(in)  :: pgaus
  integer(ip), intent(in)  :: porde
  integer(ip), intent(in)  :: lnods_loc(pnode)
  real(rp),    intent(in)  :: gpsha(pnode,pgaus)
  real(rp),    intent(in)  :: gpvol(pgaus)
  real(rp),    intent(in)  :: gpden(pgaus)
  real(rp),    intent(out) :: dt_rho_cp_tem_loc(*)
  integer(ip) :: inode,jnode,ipoin,igaus
  real(rp) :: eldtcprho(pnode), fact
! real(rp)                 :: d,T
  real(rp)                :: elmat(pnode,pnode)
  real(rp)                :: trace, elmass       

  if (kfl_rhs_scal_tem == 0 ) then

     if( porde == 1 ) then
        !
        ! Element assembly
        !
        eldtcprho = 0.0_rp
        do igaus = 1,pgaus
           fact = gpvol(igaus) / ( gpden(igaus)  * dtinv )
           do inode = 1,pnode
              eldtcprho(inode) = eldtcprho(inode) + gpsha(inode,igaus) * fact
           end do
        end do
        !
        ! Nodal projection
        !
        do inode = 1,pnode
           ipoin = lnods_loc(inode)
           dt_rho_cp_tem_loc(ipoin) = dt_rho_cp_tem_loc(ipoin) + eldtcprho(inode) 
        end do
     else
        !
        ! Element assembly
        !
        eldtcprho = 0.0_rp

        !  T = 0.0_rp
        !  d = 0.0_rp
        !  do igaus = 1,pgaus
        !     fact = gpvol(igaus) / ( gpden(igaus) * dtinv )
        !     T = T + fact
        !     do inode = 1,pnode
        !        eldtcprho(inode) = eldtcprho(inode) + gpsha(inode,igaus)**2 * fact
        !     end do
        !  end do
        !  do inode = 1,pnode
        !     d = d + eldtcprho(inode)
        !  end do
        !  do inode = 1,pnode
        !     eldtcprho(inode) = eldtcprho(inode) * T / d
        !  end do

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
           dt_rho_cp_tem_loc(ipoin) = dt_rho_cp_tem_loc(ipoin) + elmat(inode,inode)*(elmass/trace) 
        end do
     end if
  endif


end subroutine tem_rhocpdt
