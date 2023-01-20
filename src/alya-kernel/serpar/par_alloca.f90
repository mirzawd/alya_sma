!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine par_alloca(itask)
  !------------------------------------------------------------------------
  !****f* Parall/par_alloca
  ! NAME
  !    par_alloca
  ! DESCRIPTION
  !    Allocates memory for master arrays used by Parall 
  ! OUTPUT
  ! USED BY
  !    Parall
  !***
  !------------------------------------------------------------------------
  use      def_parame
  use      def_master
  use      def_parall
  use      def_domain
  use      mod_memchk
  use mod_parall, only : par_memor
  implicit none
  integer(ip), intent(in) :: itask
  integer(4)              :: istat

  select case (itask)

  case(1_ip)
     !
     ! PARR1(NPOIN)
     !
     allocate(parr1(npoin),stat=istat)     
     call memchk(zero,istat,par_memor,'PARR1','par_alloca',parr1)

  case(2_ip)
     !
     ! PARR2(NDIME,NPOIN)
     !
     allocate(parr2(ndime,npoin),stat=istat)     
     call memchk(zero,istat,par_memor,'PARR2','par_alloca',parr2)

  case(3_ip)
     !
     ! PARI1(NPOIN)
     !
     allocate(pari1(npoin),stat=istat)     
     call memchk(zero,istat,par_memor,'PARI1','par_alloca',pari1)

  case(4_ip)
     !
     ! PARR1(PARII)
     !
     allocate(parr1(parii),stat=istat)     
     call memchk(zero,istat,par_memor,'PARR1','par_alloca',parr1)

  case(5_ip)
     !
     ! PARI1(NELEM)
     !
     allocate(pari1(nelem),stat=istat)     
     call memchk(zero,istat,par_memor,'PARI1','par_alloca',pari1)


  end select

end subroutine par_alloca
