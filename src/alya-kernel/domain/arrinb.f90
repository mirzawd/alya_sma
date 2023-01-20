!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine arrinb(&
     nlist,ncoef,liste,touch,nzbou,ipoin)
  !-----------------------------------------------------------------------
  !                 
  ! This routine constructs the arrays of indexes for a mesh graph.
  ! These are organized as follows (CSR format):
  !   
  ! R_BOU(IPOIN) = coefficient of the graph where IPOIN starts,
  ! C_BOU(NZBOU) = column of the NZBOU coefficient of the graph.
  !              
  !-----------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp,lg
  use def_domain, only       :  mnodb,nnode,lnodb,ltypb,r_bou,c_bou
  implicit none
  integer(ip), intent(in)    :: nlist
  integer(ip), intent(in)    :: ncoef,ipoin
  integer(ip), intent(in)    :: liste(nlist)
  integer(ip), intent(inout) :: nzbou
  logical(lg), intent(inout) :: touch(ncoef)
  integer(ip)                :: jboun,jnode,jpoin,nnodj,jposi,jlist
  integer(ip)                :: kboun,knode,kpoin,nnodk,kposi,klist

  r_bou(ipoin) = nzbou

  do jlist = 1,nlist
     jboun = liste(jlist)
     nnodj = nnode(ltypb(jboun))
     do jnode = 1,nnodj
        jpoin = lnodb(jnode,jboun)
        jposi = (jlist-1)*mnodb+jnode
        if(.not.touch(jposi)) then
           do klist = 1,nlist
              kboun = liste(klist)
              nnodk = nnode(ltypb(kboun))
              do knode = 1,nnodk
                 kpoin = lnodb(knode,kboun)
                 if(kpoin==jpoin) then
                    kposi = (klist-1)*mnodb+knode
                    touch(kposi) = .true.
                 end if
              end do
           end do
           c_bou(nzbou) = jpoin
           nzbou = nzbou+1
        end if
     end do
  end do

end subroutine arrinb
