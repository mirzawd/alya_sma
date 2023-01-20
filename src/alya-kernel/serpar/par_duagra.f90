!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup Parall
!> @{
!> @file    par_duagra.f90
!> @author  Guillaume Houzeaux
!> @brief   Parallel subdomain dual graph
!> @details Compute the subdomain dual graph.
!!
!!          \verbatim
!!
!!            A simple example to illustrate arrays:
!!
!!            +----+------------+
!!            |    |      2     |
!!            | 1  +------+-----+
!!            |    |      |     |
!!            +----+      |  3  |
!!            |           |     |
!!            |     4     |     |
!!            |           |     |
!!            +-----------+-----+
!!
!!            - IADUAL and JADUAL: Dual graph
!!            - NBDUAL: number of nodes of the dual graph
!!            - TRANSDUAL: from adjacancies to edges of the dual graph
!!
!!                   1
!!            (1)---------(2)
!!              \         /|
!!                \   3 /  |
!!               2  \ /    | 4            
!!                  / \    |              
!!                /     \  |
!!              /         \|
!!             (3)--------(4)         
!!                    5
!!
!!            NBDUAL = 5
!!            I = 1 ... JADUAL(IADUAL(I):IADUAL(ii+1)-1) = 3,4,2
!!            I = 2 ... JADUAL(IADUAL(I):IADUAL(ii+1)-1) = 4,5,1
!!            I = 3 ... JADUAL(IADUAL(I):IADUAL(ii+1)-1) = 5,1,4
!!            I = 4 ... JADUAL(IADUAL(I):IADUAL(ii+1)-1) = 2,5,1,3
!!            I = 5 ... JADUAL(IADUAL(I):IADUAL(ii+1)-1) = 2,4,3
!!
!!            I = 1 ... ADJDOM(XADJDOM(I+0)) = 2   =>   TRANSDUAL(XADJDOM(I+0)) = 1 : subd1 and subd2: edge 1
!!            I = 1 ... ADJDOM(XADJDOM(I+1)) = 4   =>   TRANSDUAL(XADJDOM(I+1)) = 2 : subd1 and subd2: edge 2
!!            I = 2 ... ADJDOM(XADJDOM(I+0)) = 1   =>   TRANSDUAL(XADJDOM(I+0)) = 1 : subd2 and subd1: edge 1
!!            I = 2 ... ADJDOM(XADJDOM(I+1)) = 3   =>   TRANSDUAL(XADJDOM(I+1)) = 3 : subd2 and subd3: edge 3
!!            I = 2 ... ADJDOM(XADJDOM(I+2)) = 4   =>   TRANSDUAL(XADJDOM(I+2)) = 4 : subd2 and subd4: edge 4
!!            I = 3 ... ADJDOM(XADJDOM(I+0)) = 2   =>   TRANSDUAL(XADJDOM(I+0)) = 3 : subd3 and subd2: edge 3
!!            I = 3 ... ADJDOM(XADJDOM(I+1)) = 4   =>   TRANSDUAL(XADJDOM(I+1)) = 5 : subd3 and subd4: edge 5
!!            I = 4 ... ADJDOM(XADJDOM(I+0)) = 1   =>   TRANSDUAL(XADJDOM(I+0)) = 2 : subd4 and subd2: edge 2
!!            I = 4 ... ADJDOM(XADJDOM(I+1)) = 2   =>   TRANSDUAL(XADJDOM(I+1)) = 4 : subd4 and subd2: edge 4
!!            I = 4 ... ADJDOM(XADJDOM(I+2)) = 3   =>   TRANSDUAL(XADJDOM(I+2)) = 5 : subd4 and subd3: edge 5
!!
!!          \endverbatim
!> @} 
!------------------------------------------------------------------------
subroutine par_duagra(nbDual,memor)
  use def_parame
  use def_parall
  use mod_memchk
  use mod_parall, only : par_memor
  use mod_memory, only : memory_alloca
  implicit none
  integer(ip), intent(out) :: nbDual
  integer(8),  intent(out) :: memor(2)
  integer(ip)              :: ii, jj, vv, ww, t1, t2
  integer(4)               :: istat
  logical(lg)              :: fin
  integer(ip), allocatable :: iwa(:)
  !
  ! Construir vector para asignar un id a cada arista y asi traducir las 
  ! aristas como nodos del nuevo grafo 
  !
  call memory_alloca(par_memor ,'TRANSLDUAL','par_memory',translDual,xadjDom(npart_par+1)-1)
  
  nbDual = 0
  do vv = 1, npart_par
     do ii = xadjDom(vv), xadjDom(vv+1)-1
        ww = adjDom(ii)
        if (vv < ww) then
           nbDual = nbDual + 1
           translDual(ii) = nbDual
        else
           fin = .false.
           jj  = xadjDom(ww)
           do while ((jj < xadjDom(ww+1)) .and. (.not. fin))
              if (adjDom(jj) == vv) then
                 translDual(ii) = translDual(jj)
                 fin = .true.             
              endif
              jj = jj + 1
           enddo
        endif
     enddo
  enddo
  !
  ! Compute size of jaDual (nodes)
  !
  allocate(iwa(nbDual),stat=istat)
  call memchk(zero,istat,memor,'iwa','par_metis',iwa)

  do ii = 1, nbDual
     iwa(ii) = 0
  enddo

  do ii = 1, xadjDom(npart_par+1)-1
     vv = adjDom(ii)
     t1 = translDual(ii)
     do jj = xadjDom(vv), xadjDom(vv+1)-1
        if (translDual(jj) /= t1) then
           iwa(t1) = iwa(t1) + 1
        endif
     enddo
  enddo
  !
  ! Construct iaDual (nodes)
  !
  call memory_alloca(par_memor ,'IADUAL','par_memory',iaDual,nbDual+1)

  iaDual(1) = 1
  do ii = 1, nbDual
     iaDual(ii+1)   = iaDual(ii) + iwa(ii)
     iwa(ii) = iaDual(ii)
  enddo
  !
  ! Construct jaDual (not ordered)
  !
  call memory_alloca(par_memor ,'JADUAL','par_memory',jaDual,iaDual(nbDual+1)-1)

  do ii = 1,xadjDom(npart_par+1)-1
     vv = adjDom(ii)
     t1 = translDual(ii)
     do jj = xadjDom(vv),xadjDom(vv+1)-1
        t2 = translDual(jj)
        if (t1 /= t2) then
           jaDual(iwa(t1)) = t2
           iwa(t1)         = iwa(t1) + 1
        endif
     enddo
  enddo

  call memchk(two,istat,memor,'IWA','par_duagra',iwa)
  deallocate(iwa,stat=istat)
  if(istat/=0) call memerr(two,'IWA','par_duagra',0_ip)

end subroutine par_duagra
