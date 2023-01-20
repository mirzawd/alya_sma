!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!




subroutine tur_elmcost_all(&
		      eltur,pnode,pgaus,gpsha,gpvol,lnods,costf)

  !------------------------------------------------------------------------
  !****f* Nastin/tur_elmcost_all
  ! NAME 
  !    tur_elmcost_all
  ! DESCRIPTION
  !    This subroutine calculates elemental contributions of the cost .. f
  ! USES
  ! USED BY
  !    tur_elmope
  !------------------------------------------------------------------------

  use def_kintyp, only       :  ip,rp
  use def_kermod, only       :  kfl_cost_type
  use def_turbul, only       :  iunkn_tur,nturb_tur

  implicit none
  
  integer(ip), intent(in)    :: pnode,pgaus
  
  real(rp)    :: gpunk(pgaus)
  integer(ip) :: inode,igaus
  
  real(rp),       intent(in)        :: gpvol(pgaus)                      ! |J|*w
  real(rp),       intent(in)        :: eltur(nturb_tur,pnode,3),gpsha(pnode,pgaus)
  real(rp),       intent(inout)     :: costf
  integer(ip),    intent(in)        :: lnods(pnode)

  if (kfl_cost_type == 6 .and. iunkn_tur == 1) then  
    !
    !  Initialization
    !
    do igaus=1,pgaus
	  gpunk(igaus) = 0.0_rp
    end do
    !
    ! read the veloc from the forward values
    !    
    do igaus=1,pgaus
      do inode =1,pnode
        gpunk(igaus) = gpunk(igaus) + gpsha(inode,igaus)*eltur(iunkn_tur,inode,1)
      end do
      !$OMP ATOMIC
!       costf = costf + gpvol(igaus)*gpunk(igaus)**2
!       costf = costf + gpunk(igaus)**2
      costf = costf + gpunk(igaus)
    end do
  endif

end subroutine tur_elmcost_all
