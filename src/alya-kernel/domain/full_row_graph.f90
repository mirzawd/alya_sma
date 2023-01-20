!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Domain
!> @{
!> @file    full_row_graph.f90
!> @author  Guillaume Houzeaux
!> @date    27/03/2017
!> @brief   Full row graph
!> @details Full row graph for own nodes
!>          NZDOM_OWN .................... Own graph size
!>          ( R_DOM_OWN , C_DOM_OWN ) .... Own graph
!>          R_DOM_INI .................... Own graph (split between own and others)
!>          R_DOM_END .................... Own graph (split between own and others)
!>
!>          R_DOM_MIN is computed to loop on own or other's boundary
!>          node positions in the matrix. For example:
!>          R_DOM_OWN(II) => R_DOM_END(II+1)-1 : own nodes
!>          R_DOM_INI(II) => R_DOM_OWN(II+1)-1 : other nodes
!> @}
!-----------------------------------------------------------------------

subroutine full_row_graph(meshe)

  use def_kintyp,   only     :  ip,rp
  use def_master,   only     :  IMASTER,ISEQUEN,IPARALL
  use def_domain,   only     :  r_dom_own,c_dom_own
  use def_domain,   only     :  nzdom,memor_dom
  use def_domain,   only     :  mesh_type
  use def_domain,   only     :  npoin_halo
  use def_domain,   only     :  nzdom_own
  use def_kermod,   only     :  kfl_full_rows
  use mod_graphs,   only     :  graphs_elepoi
  use mod_graphs,   only     :  graphs_poipoi
  use mod_memory,   only     :  memory_alloca_min
  use mod_memory,   only     :  memory_alloca
  use mod_messages, only     : messages_live
  implicit none

  type(mesh_type), intent(inout) :: meshe    !< Mesh type
  integer(ip)                    :: bandw
  integer(ip)                    :: ii,iz,jz
  real(rp)                       :: profi
  
  if( kfl_full_rows /= 0 ) then

     call messages_live('PARALL: COMPUTE GRAPH FOR OWN NODES (FOR FULL ROW MATRIX FORMAT)')

     if( IMASTER ) then
        !
        ! Master allocates minimum memory
        !
        call memory_alloca_min(memor_dom,'MESHE % R_DOM_OWN','full_row_graph',meshe % r_dom_own)
        call memory_alloca_min(memor_dom,'MESHE % C_DOM_OWN','full_row_graph',meshe % c_dom_own)

     else

        if( ISEQUEN ) then
           !
           ! Sequential: full row graph is original graph
           !
           meshe % nzdom_own =  meshe % nzdom
           meshe % r_dom_own => meshe % r_dom
           meshe % c_dom_own => meshe % c_dom
           nzdom_own         =  nzdom
           
        else
           !
           ! Parallel
           !
           call graphs_poipoi(&
                meshe % npoin_own,meshe % nelem_2,meshe % mnode,meshe % lnods,&
                meshe % lnnod,meshe % ltype,meshe % r_dom_own,meshe % c_dom_own,&
                bandw,profi,IA_NAME='MESHE % R_DOM_OWN',JA_NAME='MESHE % C_DOM_OWN',&
                memor=memor_dom)
        end if

        r_dom_own         => meshe % r_dom_own
        c_dom_own         => meshe % c_dom_own
        nzdom_own         =  r_dom_own(meshe % npoin_own+1)-1
        meshe % nzdom_own =  nzdom_own
        !
        ! Some check
        !
        if( maxval(c_dom_own(1:nzdom_own)) > npoin_halo ) &
             call runend('FULL_ROW_GRAPH: WRONG FULL ROW GRAPH')
        !
        ! Compute R_DOM_MIN
        !
        if( IPARALL ) then
           call memory_alloca(memor_dom,'MESHE % R_DOM_END','memgeo',meshe % r_dom_end,meshe % npoin_own+1_ip)
           call memory_alloca(memor_dom,'MESHE % R_DOM_INI','memgeo',meshe % r_dom_ini,meshe % npoin_own+1_ip)
           meshe % r_dom_end(1) = 1
           meshe % r_dom_ini(1) = 1
           do ii = 1,meshe % npoin_own
              iz = meshe % r_dom_own(ii)
              jz = 1
              iz_loop: do while( meshe % c_dom_own(iz) <= meshe % npoin_own .and. iz <= meshe % r_dom_own(ii+1)-1 )
                 iz = iz + 1
                 jz = jz + 1
                 if(  iz <= meshe % r_dom_own(ii+1)-1 ) exit
              end do iz_loop
              meshe % r_dom_end(ii+1) = iz
              meshe % r_dom_ini(ii)   = iz
           end do
        end if

     end if

  end if

end subroutine full_row_graph
