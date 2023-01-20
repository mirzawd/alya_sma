!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine par_deallo(itask)
  !------------------------------------------------------------------------
  !****f* Parall/par_deallo
  ! NAME
  !    par_deallo
  ! DESCRIPTION
  !    Deallocate memory for master arrays used by Parall 
  ! OUTPUT
  ! USED BY
  !    Parall
  !***
  !------------------------------------------------------------------------
  use def_parame
  use def_master
  use def_parall
  use mod_memchk
  use mod_parall, only : par_memor
  implicit none
  integer(ip), intent(in) :: itask
  integer(4)              :: istat

  select case (itask)

  case(1_ip)
     !
     ! PARR1
     !
     call memchk(two,istat,par_memor,'PARR1','par_deallo',parr1)
     deallocate(parr1,stat=istat)
     if(istat/=0) call memerr(two,'par_deallo','PARR1',0_ip) 

  case(2_ip)
     !
     ! PARR2
     !
     call memchk(two,istat,par_memor,'PARR2','par_deallo',parr2)
     deallocate(parr2,stat=istat)
     if(istat/=0) call memerr(two,'par_deallo','PARR2',0_ip) 

  case(3)
     !
     ! PARR3
     !
     call memchk(two,istat,par_memor,'PARR3','par_deallo',parr3)
     deallocate(parr3,stat=istat)
     if(istat/=0) call memerr(two,'par_deallo','PARR3',0_ip) 

  case(4_ip)
     !
     ! PARI1
     !
     call memchk(two,istat,par_memor,'PARI1','par_deallo',pari1)
     deallocate(pari1,stat=istat)
     if(istat/=0) call memerr(two,'par_deallo','PARI1',0_ip) 

  end select

end subroutine par_deallo
