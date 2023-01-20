!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Parall
!> @{                                                                   
!> @file    par_exgraf.f90
!> @author  Guillaume Houzeaux
!> @date    19/11/2013
!> @brief   Exchange a matrix on interface                          
!> @details        
!>             subdomain 2  
!>
!>              37     8 
!>               o------o
!>           2  / \  2 / \
!>             / 1 \  / 1 \
!>            /     \/     \
!>         2 o      o 9     o 7
!>           o------o-------o
!>           |      |       |   
!>           |   1  |   1   |
!>           |      |       |   subdomain 1
!>           o------o-------o
!>           6      5       4
!>
!>           subdomain 1 and subdomain 2 share nodes 2 and 13 but
!>           subdomain 2 could not have 9-2 in its graph
!>               
!> @}                                                                   
!-----------------------------------------------------------------------

subroutine par_exgraf()

  use def_kintyp
  use def_master
  use def_domain
  use mod_memory
  use mod_parall,         only : commd
  use mod_communications, only : PAR_SEND_RECEIVE
  use mod_parall,         only : NODE_IN_NEIGHBOR
  use mod_parall,         only : par_memor
  use def_solver,         only : solve_sol
  implicit none
  integer(ip)          :: ipoin,ii,jj,bsize,ll,idofn,jdofn
  integer(ip)          :: kpoin,kk,jpoin,dom_i,ineig,izdom
  real(rp),    pointer :: loc_sparr1(:),loc_rparr1(:)
  type(i1p),   pointer :: listn(:)
  integer(ip), pointer :: list_nodes(:)
  integer(ip), pointer :: list_nodes_loc(:)

  !-------------------------------------------------------------
  !
  ! Count information to send and receive from my neighbor
  !
  !-------------------------------------------------------------

  nullify(loc_sparr1)
  nullify(loc_rparr1)
  nullify(listn)
  nullify(list_nodes)
  nullify(list_nodes_loc)  

  do ineig = 1,commd % nneig
     !
     ! Loop over subdomains INEIG
     !
     dom_i = commd % neights(ineig)
     bsize = commd % bound_size(ineig+1) - commd % bound_size(ineig)

     call memory_alloca(par_memor,'LIST_NODEA',    'par_exgraf',list_nodes,bsize)
     call memory_alloca(par_memor,'LOST_NODES_LOC','par_exgraf',list_nodes_loc,bsize)
     call memory_alloca(par_memor,'LISTN',         'par_exgraf',listn,bsize)

     jj = 0
     do ii = commd % bound_size(ineig),commd % bound_size(ineig+1)-1        
        !
        ! Loop over the nodes IPOIN shared with INEIG
        !        
        ipoin = commd % bound_perm(ii)
        jj    = jj + 1
        kk    = 0
        do izdom = r_dom(ipoin),r_dom(ipoin+1)-1
           jpoin = c_dom(izdom)
           !
           ! JPOIN is also in INEIG?
           ! 
           if(  NODE_IN_NEIGHBOR(jpoin,ineig,commd) ) then 
              kk = kk + 1
              list_nodes(kk) = jpoin
           end if
        end do
        !
        ! Fill in neighbors of IPOIN in INEIG: KK nodes common
        !
        if( kk > 0 ) then
           allocate( listn(jj) % l(kk) )
           do ll = 1,kk
              jpoin              = list_nodes(ll)
              listn(jj) % l(ll)  = jpoin
              list_nodes_loc(ll) = lninv_loc(jpoin)
           end do
           call heapsorti2(1_ip,kk,list_nodes_loc,listn(jj) % l)
        else
           nullify(listn(jj) % l)
        end if
     end do
     !
     ! Count number of coefficients to send
     !
     bsize = 0
     do ii = 1,jj
        if( associated(listn(ii) % l) ) then
           bsize = bsize + size(listn(ii) % l)
        end if
     end do

if((kfl_paral==3.and.dom_i==4).or.(kfl_paral==4.and.dom_i==3)) print*,'BSIZE=',bsize 
     bsize = pard1 * pard1 * bsize 
     call memory_alloca(par_memor,'LOC_SPARR1','solope',loc_sparr1,bsize)
     call memory_alloca(par_memor,'LOC_RPARR1','solope',loc_rparr1,bsize)
     !
     ! Fill in what to send
     !  
     jj = 0
     kk = 0
     do ii = commd % bound_size(ineig),commd % bound_size(ineig+1)-1  
        ipoin = commd % bound_perm(ii)
        jj = jj + 1
        do ll = 1,size(listn(jj) % l)
           jpoin = listn(jj) % l(ll)
           izdom = r_dom(ipoin)
           do while( izdom <= r_dom(ipoin+1)-1 )
              kpoin = c_dom(izdom)
              if( kpoin == jpoin ) then
                 do idofn = 1,pard1
                    do jdofn = 1,pard1
                       kk = kk + 1
                       loc_sparr1(kk) = geten(idofn,jdofn,izdom)
                    end do
                 end do
                 izdom = r_dom(ipoin+1)+2
              end if
              izdom = izdom + 1
           end do
           if( izdom /= r_dom(ipoin+1)+3 ) call runend('PAR_EXGRAF: COULD NOT FIND NODE')
        end do
     end do
     !
     ! Exchange values
     !
     if( bsize > 0 ) call PAR_SEND_RECEIVE(bsize,bsize,loc_sparr1,loc_rparr1,'IN MY CODE',dom_i)
     !
     ! Assemble the contribution of my neighbor
     !
     jj = 0
     kk = 0
     do ii = commd % bound_size(ineig),commd % bound_size(ineig+1)-1  
        ipoin = commd % bound_perm(ii)
        jj = jj + 1
        do ll = 1,size(listn(jj) % l)
           jpoin = listn(jj) % l(ll)
           izdom = r_dom(ipoin)
           do while( izdom <= r_dom(ipoin+1)-1 )
              kpoin = c_dom(izdom)
              if( kpoin == jpoin ) then
                 do idofn = 1,pard1
                    do jdofn = 1,pard1
                       kk = kk + 1 
                       geten(idofn,jdofn,izdom) = geten(idofn,jdofn,izdom) + loc_rparr1(kk) 
                    end do
                 end do
                 izdom = r_dom(ipoin+1)
              end if
              izdom = izdom + 1
           end do
        end do
     end do

     call memory_deallo( par_memor,'LOC_SPARR1','solope',loc_sparr1)
     call memory_deallo( par_memor,'LOC_RPARR1','solope',loc_rparr1) 

     call memory_deallo( par_memor,'LIST_NODEA',    'par_exgraf', list_nodes     )
     call memory_deallo( par_memor,'LOST_NODES_LOC','par_exgraf', list_nodes_loc )
     call memory_deallo( par_memor,'LISTN',         'par_exgraf', listn          )

  end do

end subroutine par_exgraf

  
