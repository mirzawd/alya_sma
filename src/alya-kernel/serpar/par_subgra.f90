!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine par_subgra( &
     nbnodes, ia, ja, dom, part,    &
     nodeI, nodeB, xadj, adj,       &
     permI, invI,permB, invB ,      &
     adjsize ,dsiz1, dsiz2, dsize,  &
     ierro )
  !-------------------------------------------------------------------------------
  !****f* parall/par_subgra
  ! NAME
  !    par_subgra
  ! DESCRIPTION
  !    Create a subgraph for a subdomain using LNPAR_PAR:
  !    NBNODINTER ..... # of internal nodes
  !    NBNODBOUND ..... # of boundary nodes
  !    ADJ, XADJ ...... Adjacancy for internal nodes (subgraph)
  !    PERMI, INVPI ... Permutation/inverse for interior nodes
  !    PERMB, INVPB ... Permutation/inverse for boundary nodes
  ! INPUT
  !    nelem
  !    ia
  !    ja
  !    dom
  !    part
  ! OUTPUT
  !    xadj
  !    adj
  !    permGD
  !    invGD
  ! USED BY
  !***
  !-------------------------------------------------------------------------------
  use def_kintyp, only     :  ip,rp
  implicit none
  integer(ip), intent(in)  :: nbnodes,dom
  integer(ip), intent(in)  :: ia(*), ja(*), part(*)
  integer(ip), intent(out) :: nodeI, nodeB, adjsize
  integer(ip), intent(out) :: xadj(*), adj(*)
  integer(ip), intent(out) :: permI(*), invI(*)
  integer(ip), intent(out) :: permB(*), invB(*)
  integer(ip), intent(in)  :: dsiz1,dsiz2,dsize
  integer(ip), intent(out) :: ierro
  integer(ip)              :: next, vv, ww, ii

  nodeI   = 1
  nodeB   = 1
  xadj(1) = 1
  next    = 1
  ierro   = 0
  
  do vv = 1, nbnodes

     if(part(vv) == dom) then
        permI(vv)   = nodeI

        if( nodeI > dsize ) then
           ierro = 3
           write(*,*)'nodeI > dsize',nodeI,dsize
           return
        end if
        
        invI(nodeI) = vv

        do ii = ia(vv), ia(vv+1)-1
           ww = ja(ii)
           if (ww/=vv) then
              if (part(ww) == dom) then
                 if( next > dsiz2 ) then
                    ierro = 1
                    write(*,*)'next > dsiz2',next,dsiz2
                    return
                    !call runend('PAR_SUBGRA: RESIZE XADJ - next > dsiz2')
                 end if
                 adj(next) = ww
                 next      = next + 1
              end if
           end if
        end do

        nodeI       = nodeI + 1
        if( nodeI > dsiz1 ) then
           ierro = 2
           write(*,*)'nodeI > dsiz1',nodeI,dsiz1
           return
           !call runend('PAR_SUBGRA: RESIZE XADJ - nodeI > dsiz1')
        end if
        xadj(nodeI) = next

     else if(part(vv) == -dom) then
        permB(vv)   = nodeB
        if( nodeB > dsize ) then
           ierro = 3
           write(*,*)'nodeB > dsize',nodeB,dsize
           return
        end if
        invB(nodeB) = vv
        nodeB       = nodeB + 1

     end if

  end do

  if( next-1 > dsiz2 ) call runend('PAR_SUBGRA: RESIZE ADJ')
  do ii = 1, next-1
     adj(ii) = permI(adj(ii))
  end do

  adjsize = next
  nodeI   = nodeI - 1
  nodeB   = nodeB - 1

end subroutine par_subgra


