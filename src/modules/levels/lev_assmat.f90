!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine lev_assmat(&
     elmat,amatr,pnode,mnode,nequa,lnods,lplev,algso,kfl_timet_lev)
  !-----------------------------------------------------------------------
  !****f* Levels/lev_assmat
  ! NAME 
  !    lev_assmat
  ! DESCRIPTION
  !    This routine performs the assembly of the big matrix
  ! USES
  !    skyase
  !    csrase
  ! USED BY
  !    lev_elmope
  !***
  !-----------------------------------------------------------------------
  use def_kintyp
  use def_solver
  implicit none
  integer(ip), intent(in)    :: pnode,nequa,mnode,kfl_timet_lev
  integer(ip), intent(in)    :: algso
  integer(ip), intent(in)    :: lplev(nequa)
  integer(ip), intent(in)    :: lnods(pnode)
  real(rp),    intent(in)    :: elmat(mnode,mnode)
  real(rp),    intent(inout) :: amatr(*)

  if(algso==0.and.kfl_timet_lev==2) then
     call runend('LEV_ASSMAT: OPTION NO LONGER EXISTS')
  else
     call csrase(elmat,amatr,1_ip,pnode,mnode,lnods,2_ip)
  end if
  
end subroutine lev_assmat
