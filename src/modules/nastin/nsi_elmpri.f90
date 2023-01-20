!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_elmpri(&
     pnode,pgaus,lnods,gpcar,gpvol,gpden,gpvis,&
     gppor,gpsha,chale,dtinv_loc,elmap)
  !----------------------------------------------------------------------
  !****f* Nastin/nsi_elmsch
  ! NAME 
  !    nsi_elmsch
  ! DESCRIPTION
  !    Compute the Schur complement preconditioner
  ! USES 
  ! USED BY
  !***
  !----------------------------------------------------------------------
  use def_kintyp, only     :  ip,rp 
  use def_nastin, only     :  staco_nsi,kfl_taust_nsi
  use def_domain, only     :  mnode,ndime
  use mod_tauadr, only     :  tauadr
  implicit none
  integer(ip), intent(in)  :: pnode,pgaus
  integer(ip), intent(in)  :: lnods(pnode)
  real(rp),    intent(in)  :: gpcar(ndime,mnode,pgaus),gpvol(pgaus)
  real(rp),    intent(in)  :: gpden(pgaus),gpvis(pgaus),gppor(pgaus)
  real(rp),    intent(in)  :: gpsha(pnode,pgaus)
  real(rp),    intent(in)  :: chale(ndime)
  real(rp),    intent(in)  :: dtinv_loc
  real(rp),    intent(out) :: elmap(pnode,pnode)
  integer(ip)              :: inode,jnode,kdime,igaus
  real(rp)                 :: fact1,fact3
  real(rp)                 :: tau,adv,dif,rea

  !----------------------------------------------------------------------
  !
  ! Elemental matrix: ( tau*grad p , grad q )
  !
  !----------------------------------------------------------------------

  do jnode = 1,pnode
     do inode = 1,pnode
        elmap(inode,jnode) = 0.0_rp
     end do
  end do

  do igaus = 1,pgaus
     
     adv = 0.0_rp                          ! Convective term
     dif = gpvis(igaus)                    ! Viscous term
     rea = gppor(igaus)                    ! Reaction term
     call tauadr(&
          kfl_taust_nsi,staco_nsi,adv,dif,rea,&
          chale(1),chale(2),tau,gpden(igaus)*dtinv_loc)
     fact3 = tau * gpvol(igaus)
     !fact3 = gpvol(igaus)
     do inode = 1,pnode
        do jnode = inode+1,pnode
           fact1 = 0.0_rp
           do kdime = 1,ndime
              fact1 = fact1 + gpcar(kdime,inode,igaus)&
                   &        * gpcar(kdime,jnode,igaus)
           end do
           fact1 = fact1 * fact3
           elmap(inode,jnode) = elmap(inode,jnode) + fact1
           elmap(jnode,inode) = elmap(jnode,inode) + fact1
        end do
        fact1 = 0.0_rp
        do kdime = 1,ndime
           fact1 = fact1 + gpcar(kdime,inode,igaus)&
                &        * gpcar(kdime,inode,igaus)
        end do
        fact1 = fact1 * fact3
        elmap(inode,inode) = elmap(inode,inode) + fact1
     end do
  end do

end subroutine nsi_elmpri
