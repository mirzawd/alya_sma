!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!




subroutine tur_elmdcost_all(&
		    eltur,pnode,pgaus,gpsha,gpvol,lnods,elrhs)
  !------------------------------------------------------------------------
  !****f* Nastin/tur_elmdcost_all
  ! NAME 
  !    tur_elmdcost_all
  ! DESCRIPTION
  !    This subroutine calculates elemental contributions of the dcost .. d(f)/d(U)
  ! USES
  ! USED BY
  !    tur_elmope
  !------------------------------------------------------------------------

  use def_kintyp, only       :  ip,rp
  use def_kermod, only       :  kfl_cost_type
  use def_master
  use mod_matrix
  use def_turbul, only       :  iunkn_tur,nturb_tur
  
  implicit none

  integer(ip), intent(in)    :: pnode,pgaus
  
  
  real(rp)    :: gpunk(pgaus)
  integer(ip) :: inode,igaus
  
  integer(ip),    intent(in)        :: lnods(pnode)
  real(rp),       intent(in)        :: gpvol(pgaus)                          ! |J|*w
  real(rp),       intent(in)        :: eltur(nturb_tur,pnode,3),gpsha(pnode,pgaus)
  real(rp),       intent(inout)     :: elrhs(pnode)
  
  real(rp)                          :: elmdcost(pnode)
  
  
  if (kfl_cost_type == 6 .and. iunkn_tur == 1) then
    !
    ! initializations
    !
    elmdcost = 0.0_rp
    gpunk = 0.0_rp
    !
    ! calculate elmdcost
    !
    do igaus=1,pgaus
      do inode =1,pnode
        gpunk(igaus) = gpunk(igaus) + gpsha(inode,igaus)*eltur(iunkn_tur,inode,1)
      end do
      do inode =1,pnode
!         elmdcost(inode) = elmdcost(inode) + gpvol(igaus)*2*gpsha(inode,igaus)*gpunk(igaus)
!         elmdcost(inode) = elmdcost(inode) + 2*gpsha(inode,igaus)*gpunk(igaus)
        elmdcost(inode) = elmdcost(inode) + gpsha(inode,igaus)
      end do
    end do 
    !
    ! elrhs = elrhs + elmdcost
    !  
    do inode = 1,pnode
      elrhs(inode) = elrhs(inode) - elmdcost(inode)
    end do
  endif
  
 
end subroutine tur_elmdcost_all
