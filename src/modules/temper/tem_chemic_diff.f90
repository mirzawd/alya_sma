!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tem_chemic_diff(&
     itask,ielem,pgaus,pnode,elmdiff)
  !-----------------------------------------------------------------------
  !****f* Temper/tem_chemic_diff
  ! NAME
  !   tem_radiat
  ! DESCRIPTION
  !    Couple to the heat source from chemical processes
  ! USES
  ! USED BY
  !    tem_elmop2 
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp
  use def_master, only       :  chemicalHeatdiff
  use def_kermod, only       :  kfl_ndvars_opt

  implicit none 
  integer(ip), intent(in)    :: itask,ielem,pgaus,pnode
  real(rp),    intent(out)   :: elmdiff(pgaus,kfl_ndvars_opt)
  integer(ip)                :: igaus,ides_var
  
  select case (itask)
  case(1)
      ! chemicalHeatdiff
      do ides_var=1,kfl_ndvars_opt
	  do igaus=1,pgaus
	      elmdiff(igaus,ides_var) = chemicalHeatdiff(ielem)%a(igaus,ides_var)
	  end do
      enddo
      
  end select
  
  
end subroutine tem_chemic_diff
