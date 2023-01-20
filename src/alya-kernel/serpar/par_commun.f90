!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup Parall
!> @{
!> @file    par_commun.f90
!> @author  Guillaume Houzeaux
!> @brief   Compute the communication strategy
!> @details Compute the communication strategy. Some variables
!!          have been computed in previous subroutines:
!!          - 1. par_domgra: subdomain adjacancy graph: ADJDOM and XADJDOM
!!          - 2. par_duagra: dual graph: IADUAL, JADUAL and TRANSLDUAL
!!          - 3. par_colgra: color graph: NBCOLOR, COLORS
!!          - 4. par_commun: communcation scheduling: LCOMM_PAR
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
!!            1. PAR_DOMGRA
!!            -------------
!!
!!            XADJDOM and ADJDOM: Subdomain adjancancy graph
!!
!!            I = 1 ... ADJDOM(XADJDOM(I):XADJDOM(I+1)-1) = 2,4
!!            I = 2 ... ADJDOM(XADJDOM(I):XADJDOM(I+1)-1) = 1,3,4
!!            I = 3 ... ADJDOM(XADJDOM(I):XADJDOM(I+1)-1) = 2,4
!!            I = 3 ... ADJDOM(XADJDOM(I):XADJDOM(I+1)-1) = 1,2,3
!!
!!            2. PAR_DUALGRA
!!            --------------
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
!!            3. PAR_COLGRA
!!            --------------
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
!!            4. PAR_COMMUN
!!            -------------
!!
!!            LCOMM_PAR: Communication scheduling, organized by color
!!            L_COMM_PAR(COLOR,SUBDOMAIN)= ADJACENT SUBDOMAIN
!!
!!            L_COMM_PAR( 1, 1 ) = 2  <= During the first communication stage, color 1,
!!            L_COMM_PAR( 1, 2 ) = 1     we have: 1 <=> 2 and 3 <=> 4
!!            L_COMM_PAR( 1, 3 ) = 4
!!            L_COMM_PAR( 1, 4 ) = 3
!!            L_COMM_PAR( 2, 1 ) = 4  <= During the 2nd communcation statge, color 3,
!!            L_COMM_PAR( 2, 2 ) = 3     we have: 1 <=> 4 and 2 <=> 3
!!            L_COMM_PAR( 2, 3 ) = 2
!!            L_COMM_PAR( 2, 4 ) = 1
!!            L_COMM_PAR( 3, 1 ) = 0  <= subdomain 1 is not involved in color 3
!!            L_COMM_PAR( 3, 2 ) = 4  <= tage 4, only 2 and 4 communicate
!!            L_COMM_PAR( 3, 3 ) = 0  <= subdomain 3 is not involved in color 3
!!            L_COMM_PAR( 3, 4 ) = 2
!!
!!          \endverbatim
!!
!> @} 
!------------------------------------------------------------------------
subroutine par_commun(nbcolours,nbdual)
  use def_parame
  use def_parall
  use def_master
  use mod_memory
  use def_kintyp_comm
  use mod_parall, only : par_memor
  implicit none
  integer(ip),    intent(in)  :: nbcolours            !< Number of colors
  integer(ip),    intent(in)  :: nbdual               !< Number of nodes of the dual graph
  integer(ip)                 :: adjncy,ii,jj,vv,ww
!  integer(4)                  :: istat
  type(tAdj_par), pointer     :: invTransl(:)

#ifdef EVENT
  call mpitrace_user_function(1)
#endif

  nullify(invTransl)
  call memory_alloca(par_memor,'LCOMM_PAR','par_commun',lcomm_par,nbcolours,npart_par)
  call memory_alloca(par_memor,'INVTRANSL','par_commun',invTransl,nbdual,'DO_NOT_INITIALIZE')

  !allocate(invTransl(nbdual),stat=istat)
  !call par_memchk(zero,istat,par_memor,'INVTRANSL','par_commun',invTransl)
  !
  ! Construct inverse of TRANSL
  ! TRANSDUAL gives the dual graph node JJ (edge connecting adjacent edges)
  ! INVTRANSL(JJ) % NODE1/NODE2 gives for a dual graph node JJ the two adjacent 
  ! subdomains NODE1 and NODE2 forming this communication edge
  !
  do vv = 1, npart_par
     do ii = xadjDom(vv),xadjDom(vv+1)-1
        ww = adjDom(ii)
        if( vv < ww ) then
           jj                    = translDual(ii)
           invTransl(jj) % node1 = vv
           invTransl(jj) % node2 = ww
        endif
     enddo
  enddo
  !
  ! Initialize communication array LCOMM_PAR
  !
  do ii = 1, nbcolours
     do jj = 1, npart_par
        lcomm_par(ii,jj) = 0
     enddo
  enddo
  !
  ! Construct communication array LCOMM_PAR
  ! Nodes of the dual graph of the same color can communicate at the same time
  ! so the strategy is the following:
  ! for each color, make adjacent subdomains of this color communicate. Each 
  ! color consists of a stage of the communcation scheduling.
  ! Communications are automatically ordered.
  !
  do ii = 1,nbcolours
     do adjncy = 1,nbdual
        if( colours(adjncy) == ii ) then
           vv               = invTransl(adjncy) % node1  
           ww               = invTransl(adjncy) % node2
           lcomm_par(ii,vv) = ww
           lcomm_par(ii,ww) = vv
        endif
     end do
  end do

  call memory_deallo(par_memor,'INVTRANSL','par_commun',invTransl)


#ifdef EVENT
  call mpitrace_user_function(0)
#endif

end subroutine par_commun

