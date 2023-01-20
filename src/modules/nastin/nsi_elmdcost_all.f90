!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!




subroutine nsi_elmdcost_all(&
		    elvel,pnode,pgaus,gpsha,gpvol,lnods,elrbu,elrbp)
  !------------------------------------------------------------------------
  !****f* Nastin/nsi_elmdcost_all
  ! NAME 
  !    nsi_elmdcost_all
  ! DESCRIPTION
  !    This subroutine calculates elemental contributions of the dcost .. d(f)/d(U)
  ! USES
  ! USED BY
  !------------------------------------------------------------------------

  use def_kintyp, only       :  ip,rp
  use def_domain, only       :  ndime
  use def_kermod, only       :  kfl_cost_type
  use def_master
  use mod_matrix
  
  implicit none

  integer(ip), intent(in)    :: pnode,pgaus
  
  
  real(rp)    :: gpunk(pgaus)
  integer(ip) :: inode,igaus,idime
  
  integer(ip),    intent(in)        :: lnods(pnode)
  real(rp),       intent(in)        :: gpvol(pgaus)                          ! |J|*w
  real(rp),       intent(in)        :: elvel(ndime,pnode),gpsha(pnode,pgaus)
  real(rp),       intent(inout)     :: elrbu(ndime,pnode)
  real(rp),       intent(inout)     :: elrbp(pnode)
  
  real(rp)                          :: elmdcost(ndime,pnode)
  
  if (kfl_cost_type == 2) then
    !
    ! initializations
    !
    do inode = 1,pnode
      do idime = 1,ndime
	  elmdcost(idime,inode) = 0.0_rp
      enddo
    end do
    do igaus=1,pgaus
	  gpunk(igaus) = 0.0_rp
    end do
    !
    ! calculate elmdcost
    !
    if (kfl_coupl(ID_CHEMIC,ID_NASTIN) == 0) then
      do igaus=1,pgaus
	  do inode =1,pnode
	    gpunk(igaus) = gpunk(igaus) + gpsha(inode,igaus)*elvel(1,inode)
	  end do
	  do inode =1,pnode
	    elmdcost(1,inode) = elmdcost(1,inode) + gpvol(igaus)*2*gpsha(inode,igaus)*gpunk(igaus)
	  end do
      end do 
    endif
    !
    ! elrbu = elrbu + elmdcost
    !  
    do inode = 1,pnode
      do idime = 1,ndime
	  elrbu(idime,inode) = elrbu(idime,inode) - elmdcost(idime,inode)
      enddo
    end do
  
  endif
  
 
end subroutine nsi_elmdcost_all
