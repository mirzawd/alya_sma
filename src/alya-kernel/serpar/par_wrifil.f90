!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine par_wrifil(iunit)
  !------------------------------------------------------------------------
  !****f* Parall/par_send
  ! NAME
  !    par_send
  ! DESCRIPTION
  !    This routine Send all buffers to process 'kfl_desti'
  ! OUTPUT
  !   
  ! USED BY
  !    Parall
  !***
  !------------------------------------------------------------------------
  use def_parall
  use def_master
  use mod_iofile
  use mod_memchk
  use mod_parall, only : par_memor
  implicit none
  integer(ip), intent(in) :: iunit
  integer(ip)             :: ipari
  integer(4)              :: iunit4,istat
  integer(4),  pointer    :: parin4(:)

  iunit4 = int(iunit,4_4)

  if( kfl_bytes_par == 4 .and. ip /= 4 ) then
     allocate(parin4(npari),stat=istat) 
     call memchk(zero,istat,par_memor,'parin4','par_wrifil',parin4)
     do ipari = 1,npari
        parin4(ipari) = int(parin(ipari),4)
     end do
  end if


  if( kfl_ascii_par == 0 ) then

     !-------------------------------------------------------------------
     !
     ! Binary file
     !
     !-------------------------------------------------------------------

     if( kfl_bytes_par == 4 .and. ip /= 4 ) then
        write(iunit4) int(npari,4),int(nparr,4),int(nparc,4)
        if( npari > 0 )   write(iunit4) parin4(1:npari) 
     else
        write(iunit4)  npari,nparr,nparc
        if( npari > 0 )   write(iunit4) parin(1:npari) 
     end if
     if( nparr > 0 )   write(iunit4) parre(1:nparr)
     if( nparc > 0 )   write(iunit4) parch(1:nparc)

  else

     !-------------------------------------------------------------------
     !
     ! ASCII file
     !
     !-------------------------------------------------------------------

     if( kfl_bytes_par == 4 .and. ip /= 4 ) then
        write(iunit,*)  int(npari,4),int(nparr,4),int(nparc,4)
        write(iunit,*)  trim(strin),' ',trim(strre),' ',trim(strch)
        if( npari > 0 ) write(iunit,*) parin4(1:npari)  
     else
        write(iunit,*)  npari,nparr,nparc
        write(iunit,*)  trim(strin),' ',trim(strre),' ',trim(strch)
        if( npari > 0 ) write(iunit,*) parin(1:npari)  
     end if
     if( nparr > 0 ) write(iunit,*) parre(1:nparr) 
     if( nparc > 0 ) write(iunit,*) parch(1:nparc)
     strin='NULL'
     strre='NULL'
     strch='NULL'

  end if

  if( kfl_bytes_par == 4 .and. ip /= 4 ) then
     call memchk(two,istat,par_memor,'parin4','par_wrifil',parin4)
     deallocate(parin4,stat=istat) 
     if(istat/=0) call memerr(two,'PARIN4','par_wrifil',0_ip)
  end if

end subroutine par_wrifil
