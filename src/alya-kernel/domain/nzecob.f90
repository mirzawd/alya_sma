!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nzecob(&
     nlist,ncoef,nzbou,liste,touch)
  !-----------------------------------------------------------------------
  !
  ! This routine computes the number of non-zero coefficients of a
  ! mesh graph stored in compressed sparse row (CSR) format 
  !                 
  !-----------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp,lg
  use def_domain, only       :  nnode,mnodb,lnodb,ltypb
  implicit none
  integer(ip), intent(in)    :: nlist,ncoef
  integer(ip), intent(in)    :: liste(nlist)
  integer(ip), intent(inout) :: nzbou
  logical(lg), intent(inout) :: touch(ncoef)
  integer(ip)                :: jboun,jnode,jpoin,nnodj,jposi,jlist
  integer(ip)                :: kboun,knode,kpoin,nnodk,kposi,klist

  do jlist=1,nlist                                      ! Loop over those elements 
     jboun=liste(jlist)                                 ! where the point is
     nnodj=nnode(ltypb(jboun))
     do jnode=1,nnodj
        jpoin=lnodb(jnode,jboun)
        jposi=(jlist-1)*mnodb+jnode
        if(.not.touch(jposi)) then                      ! Position not touched           
           do klist=1,nlist                             ! Search other elements 
              kboun=liste(klist)                        ! where JPOIN is and 
              nnodk=nnode(ltypb(kboun))
              do knode=1,nnodk                          ! touch their position
                 kpoin=lnodb(knode,kboun)
                 if(kpoin==jpoin) then
                    kposi=(klist-1)*mnodb+knode
                    touch(kposi)=.true.
                 end if
              end do
           end do
           nzbou = nzbou+1
        end if
     end do
  end do

end subroutine nzecob
