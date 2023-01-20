!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine opeclo(itask)
  !-----------------------------------------------------------------------
  !****f* domain/domain
  ! NAME
  !    domain
  ! DESCRIPTION
  !    This is the main routine of the domain. It performs the operations 
  !    needed to build up the domain data for the run.
  ! USED BY
  !    Turnon
  !-----------------------------------------------------------------------
  use def_kintyp, only : ip,elm
  use def_elmtyp, only : PYR05
  use def_domain, only : nnode,memor_dom,lexis,nelty
  use def_domain, only : ngaus,elmar,mgaus,mnode
  use mod_memchk, only : memchk
  implicit none
  integer(ip), intent(in)    :: itask
  type(elm),   pointer, save :: elmar_tmp(:) => null()
  integer(ip), pointer, save :: ngaus_tmp(:) => null()
  integer(ip),          save :: mgaus_tmp
  integer(ip)                :: ielty
  integer(4)                 :: istat

  if( itask == 1 ) then

     allocate(elmar_tmp(nelty),stat=istat)
     call memchk(0_ip,istat,memor_dom,'ELMAR_tmp','opeclo',elmar_tmp)
     allocate(ngaus_tmp(nelty),stat=istat)
     call memchk(0_ip,istat,memor_dom,'NGAUS_tmp','opeclo',ngaus_tmp)

     mgaus_tmp = mgaus
     do ielty = 1,nelty
        if( lexis(ielty) == 1 ) then 
           elmar_tmp(ielty) = elmar(ielty)
           ngaus_tmp(ielty) = ngaus(ielty)
        end if
     end do
     
     mgaus = max(mnode,mgaus)
     do ielty = 1,nelty
        if( lexis(ielty) == 1 .and. ielty /= PYR05 ) then 
           ngaus(ielty)         =  nnode(ielty)
           elmar(ielty) % shape => elmar(ielty) % shapc
           elmar(ielty) % deriv => elmar(ielty) % deric
           elmar(ielty) % heslo => elmar(ielty) % heslc
           elmar(ielty) % weigp => elmar(ielty) % weigc
        end if
     end do
 
  else if( itask == 2 ) then

     mgaus = mgaus_tmp
     do ielty = 1,nelty
        if( lexis(ielty) == 1 ) then 
           elmar(ielty) = elmar_tmp(ielty)
           ngaus(ielty) = ngaus_tmp(ielty)
        end if
     end do

     call memchk(2_ip,istat,memor_dom,'NGAUS_tmp','opeclo',ngaus_tmp)
     deallocate(ngaus_tmp,stat=istat)
     if(istat/=0) call memerr(2_ip,'NGAUS_TMP','opeclo',0_ip)

     call memchk(2_ip,istat,memor_dom,'ELMAR_tmp','opeclo',elmar_tmp)
     deallocate(elmar_tmp,stat=istat)
     if(istat/=0) call memerr(2_ip,'ELMAR_TMP','opeclo',0_ip)

  end if

end subroutine opeclo
  
