!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup Parall
!> @{
!> @file    par_colgra.f90
!> @author  Guillaume Houzeaux
!> @brief   Color dual graph
!> @details Color the subdomain dual graph. The strategy consists in
!!          using the minimum number of colors so that each node of the
!!          dual graph does not have the same color as its neighbors.
!!          Remember that a node of the dual graph is an edge connecting
!!          two adjacent subdomains. 
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
!!            NBCOLORS: Number of colors
!!            COLOURS: Colors of the dual graph
!!            Colors are in brackets near the node of the dual graph:
!!
!!                   1[1]
!!            (1)---------(2)
!!              \    3[2] /|
!!                \     /  |
!!            2[2]  \ /    | 4[3]            
!!                  / \    |              
!!                /     \  |
!!              /         \|
!!            (3)---------(4)        
!!                  5[1]
!!
!!            NBCOLORS = 3
!!            I = 1 ... COLORS(I) = 1
!!            I = 2 ... COLORS(I) = 2
!!            I = 3 ... COLORS(I) = 2
!!            I = 4 ... COLORS(I) = 3
!!            I = 5 ... COLORS(I) = 1
!!
!!          \endverbatim
!!
!> @} 
!------------------------------------------------------------------------
subroutine par_colgra(nbDual,nbColours,memor)
  use def_parame
  use def_parall
  use mod_memchk
  use mod_memory, only : memory_alloca
  use mod_parall, only : par_memor
  implicit none
  integer(ip),     intent(in)  :: nbDual             !< Number nodes of the dual graph (number of subdomain edges)	
  integer(ip),     intent(out) :: nbColours          !< Number of colors
  integer(8),      intent(out) :: memor(2)           !< Memory counter
  integer(ip),     parameter   :: maxColors = 500_ip
  integer(ip)                  :: colNode, ii, vv, ww, node
  integer(4)                   :: istat
  logical(lg)                  :: colourFound
  type :: tSortVect 
     integer(ip)               :: node
     integer(ip)               :: nbAdj
  end type tSortVect
  type(tSortVect)              :: aux
  type(tSortVect), allocatable :: sortVec(:)
  integer(ip),     allocatable :: mask(:)

  if( nbDual == 0 ) return

  call memory_alloca(par_memor ,'COLOURS','par_memory',colours,nbDual)

  allocate(sortVec(nbDual))
  !call memchk(zero,istat,memor,'sortVec','par_metis',sortVec)

  allocate(mask(0:maxColors),stat=istat)
  call memchk(zero,istat,memor,'mask','par_metis',mask)

  do vv = 1, nbDual
     sortVec(vv)%node  = vv
     sortVec(vv)%nbAdj = iaDual(vv+1) - iaDual(vv)
  end do
  !
  ! Order sortVec, de momento burbuja
  !
  do vv = nbDual, 1, -1
     do ww = 1, vv-1
        if( sortVec(ww)%nbAdj > sortVec(ww+1)%nbAdj ) then
           aux           = sortVec(ww) 
           sortVec(ww)   = sortVec(ww+1)
           sortVec(ww+1) = aux
        end if
     end do
  end do
  !
  ! Initialize color for each dual graph node
  !
  do vv = 1,nbDual
     colours(vv) = 0
  end do
  nbColours = 0
  !
  ! Use the minimum number of colors so that each node of the graph does not have the same
  ! color as its neighbors
  !
  do ii= 0,maxColors
     mask(ii) = 0
  end do
  !
  ! Subdomain 1 has color 1... just to start with something!
  !
  nbColours                = 1
  colours(sortVec(1)%node) = 1

  do vv = 2,nbDual

     node = sortVec(vv)%node

     do ii = iaDual(node),iaDual(node+1)-1
        ! Should check that colours(jaDual(ii)) < maxColors: very unlikely to happen but we never know!...
        if( colours(jaDual(ii)) > maxColors ) write(*,*) "colours(jaDual(ii)) ",colours(jaDual(ii))
        mask(colours(jaDual(ii))) = vv
     end do

     colourFound = .false.
     ii          = 1
     do while( ii <= nbColours .and. .not. colourFound )
        if( mask(ii) /= vv ) then
           colNode     = ii
           colourFound = .true.
        end if
        ii = ii + 1
     end do

     if ( .not. colourFound ) then
        nbColours = nbColours + 1
        colNode   = nbColours
     end if

     colours(sortVec(vv)%node) = colNode

  end do

  call memchk(two,istat,memor,'mask','par_colgra',mask)
  deallocate(mask,stat=istat)
  if(istat/=0) call memerr(two,'mask','par_colgra',0_ip)
  !call memchk(two,istat,memor,'sortVec','par_colgra',sortVec)
  deallocate(sortVec,stat=istat)
  !if(istat/=0) call memerr(two,'sortVec','par_colgra',0_ip)

end subroutine par_colgra
