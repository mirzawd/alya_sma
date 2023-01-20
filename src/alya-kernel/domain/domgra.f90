!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Domain
!> @{
!> @file    domgra.f90
!> @author  houzeaux
!> @date    2018-11-22
!> @brief   Compute graph
!> @details The variables computed are:
!>          NZDOM ... Number of nonzero coefficients of the graph
!>          R_DOM ... Pointer to the array of rows r_dom(npoin+1) (r_dom(ipoin) =
!>                    coefficient of the graph where row ipoin starts)
!>          C_DOM ... Pointer to the array of columns c_dom(nzdom) (c_dom (izdom)
!>                    = column of the izdom coefficient of mesh graph)
!>          MEPOI ... Max number of elements connected to nodes
!>          PELPO ... Node to element connectivity linked list
!>          LELPO ... Node to element connectivity linked list
!>          
!> @} 
!-----------------------------------------------------------------------

subroutine domgra(itask)

  use def_kintyp,   only : ip
  use def_master,   only : IMASTER
  use def_master,   only : INOTMASTER
  use def_domain,   only : npoin,nelem,mnode
  use def_domain,   only : lnods,lnnod,ltype
  use def_domain,   only : r_dom,c_dom
  use def_domain,   only : r_sol,c_sol,nzsym
  use def_domain,   only : nzsol
  use def_domain,   only : bandw_dom,profi_dom
  use def_domain,   only : pelpo,lelpo,mepoi
  use def_domain,   only : nzdom,memor_dom
  use mod_graphs,   only : graphs_poipoi
  use mod_memory,   only : memory_alloca
  use mod_memory,   only : memory_deallo
  use mod_elmgeo,   only : elmgeo_number_nodes
  use mod_messages, only : messages_live

  implicit none

  integer(ip), intent(in) :: itask
  integer(ip)             :: ielem,kfl_lnnod
  integer(ip), pointer    :: lnnod_loc(:)

  if(      itask == 1 ) then
     call messages_live('MASTER COMPUTES GRAPH TO CARRY OUT PARTITONING AND GATHERING OF ARRAYS')
  else if( itask == 2 ) then
     call messages_live('COMPUTE GRAPH')
  end if

  if( ( itask == 1 .and. IMASTER ) .or. ( itask == 2 .and. INOTMASTER ) ) then
     !
     ! Error checking... for example, if we have read with MPIO and 
     ! paritioning is parallel... Master does not have LTYPE and LNODS.
     ! We could do a gather just for this special case. Shall we?
     !
     nullify(lnnod_loc)
     kfl_lnnod = 0
     if( itask == 1 ) then
        if( .not. associated(ltype) .and. IMASTER ) call runend('YOU SHOULD USE THE PARALLEL FRONTAL APPROACH FOR GROUPS')     
     end if
     !
     ! Compute LNNOD_LOC if not already computed
     !
     if( .not. associated(lnnod) ) then
        kfl_lnnod = 1
        call memory_alloca(memor_dom,'LNNOD_LOC','memgeo',lnnod_loc,nelem)
        do ielem = 1,nelem
           lnnod_loc(ielem) = elmgeo_number_nodes(ltype(ielem),lnods(:,ielem))
        end do
     else
        lnnod_loc => lnnod
     end if
     !
     ! Compute graph
     !
     if( npoin > 0 ) then
        call graphs_poipoi(&
             npoin,nelem,mnode,lnods,lnnod_loc,ltype,r_dom,c_dom,&
             bandw_dom,profi_dom,pelpo,lelpo,mepoi,&
             IA_NAME='R_DOM',JA_NAME='C_DOM',memor=memor_dom)
        nzdom = r_dom(npoin+1)-1
     else
        call memory_alloca(memor_dom,'R_DOM','memgeo',r_dom,1_ip)
        nzdom = 0
     end if
     !
     ! Deallocate temporary LNNOD_LOC
     !
     if( kfl_lnnod == 1 ) then
        call memory_deallo(memor_dom,'LNNOD_LOC','memgeo',lnnod_loc)
     end if
     !
     ! Solver graph
     !
     c_sol => c_dom
     r_sol => r_dom
     nzsol =  nzdom
     nzsym =  (nzsol-npoin)/2 + npoin
 
  end if

end subroutine domgra
