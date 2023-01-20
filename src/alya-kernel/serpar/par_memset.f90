!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine par_memset(itask)
  !-----------------------------------------------------------------------
  !****f* Parall/par_openfi
  ! NAME
  !    par_openfi
  ! DESCRIPTION
  !    Allocate memory for node set and witness points treatment
  ! USED BY
  !    par_senset
  !    par_outprt
  !***
  !-----------------------------------------------------------------------
  use  def_parame
  use  def_master
  use  def_domain
  use  def_parall
  use  mod_memchk        
  use mod_parall, only : par_memor
  implicit none
  integer(ip), intent(in) :: itask
  integer(4)              :: istat

  select case(itask)

  case(1)
     !
     ! Node sets
     !
     allocate(nnset_par(0:npart_par),stat=istat)
     call memchk(zero,istat,par_memor,'NNSET_PAR','par_senset',nnset_par)
     if(nnset>0) then
        allocate(lnsec_par(nnset,2),stat=istat)
        call memchk(zero,istat,par_memor,'LNSEC_PAR','par_senset',lnsec_par)
     end if

  case(2)
     !
     ! Witness points
     !
     allocate(nwitn_par(npart_par),stat=istat)
     call memchk(zero,istat,par_memor,'NWITN_PAR','par_senset',nwitn_par)
     
  end select

end subroutine par_memset
