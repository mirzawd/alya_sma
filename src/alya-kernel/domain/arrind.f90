!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine arrind(&
     nlist,ncoef,liste,touch,nzdom,ipoin)
  !-----------------------------------------------------------------------
  !                 
  ! This routine constructs the arrays of indexes for a mesh graph.
  ! These are organized as follows (CSR format):
  !   
  ! R_DOM(IPOIN) = coefficient of the graph where IPOIN starts,
  ! C_DOM(NZDOM) = column of the NZDOM coefficient of the graph.
  !              
  !-----------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp,lg
  use def_domain, only       :  mnode,lnods,r_dom,c_dom,lnnod
  implicit none
  integer(ip), intent(in)    :: nlist,ncoef,ipoin
  integer(ip), intent(in)    :: liste(nlist)
  integer(ip), intent(inout) :: nzdom
  logical(lg), intent(inout) :: touch(ncoef)
  integer(ip)                :: jelem,jnode,jpoin,nnodj,jposi,jlist
  integer(ip)                :: kelem,knode,kpoin,nnodk,kposi,klist

  r_dom(ipoin) = nzdom

  do jlist = 1,nlist
     jelem = liste(jlist)
     nnodj = lnnod(jelem)
     !nnodj = nnode(ltype(jelem))
     do jnode = 1,nnodj
        jpoin = lnods(jnode,jelem)
        jposi = (jlist-1)*mnode+jnode
        if(.not.touch(jposi)) then
           do klist = 1,nlist
              kelem = liste(klist)
              nnodk = lnnod(kelem)
              !nnodk = nnode(ltype(kelem))
              do knode = 1,nnodk
                 kpoin = lnods(knode,kelem)
                 if(kpoin==jpoin) then
                    kposi = (klist-1)*mnode+knode
                    touch(kposi) = .true.
                 end if
              end do
           end do
           c_dom(nzdom) = jpoin
           nzdom = nzdom+1
        end if
     end do
  end do

end subroutine arrind
