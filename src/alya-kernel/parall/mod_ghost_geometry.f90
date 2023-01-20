!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



module mod_ghost_geometry
  use def_parame
  use def_master
  use def_domain
  use mod_parall,              only : commd,par_memor
  use mod_communications,      only : PAR_SEND_RECEIVE
  use mod_communications,      only : PAR_GHOST_ELEMENT_EXCHANGE
  use mod_communications,      only : PAR_GHOST_NODE_EXCHANGE
  use mod_communications,      only : PAR_GHOST_BOUNDARY_EXCHANGE
  use mod_communications,      only : PAR_SUM
  use mod_communications,      only : PAR_BARRIER
  use mod_renumbering,         only : renumbering_node_arrays
  use mod_memory,              only : memory_alloca
  use mod_memory,              only : memory_deallo
  use mod_memory,              only : memory_resize
  use mod_mesh_type,           only : mesh_type_update_last_mesh
  use mod_messages,            only : livinf
  use mod_domain,              only : domain_memory_allocate
  use mod_domain,              only : domain_memory_deallocate
  use mod_elmgeo,              only : elmgeo_number_nodes
  use mod_maths_arrays
  implicit none
  private

  public :: par_ghost_geometry

contains

  subroutine par_ghost_geometry()


    use def_domain,           only  : htable_lninv_loc  
    use mod_htable,           only  : HtableMaxPrimeNumber
    use mod_htable,           only  : hash_t
    use mod_htable,           only  : htaini
    use mod_htable,           only  : htaadd
    use mod_htable,           only  : htades
    use mod_htable,           only  : htalid

    !-----------------------------------------------------------------------
    !****f* Parall/par_fringe
    ! NAME
    !    par_fringe
    ! DESCRIPTION
    !    This routine increases geometrical arrays by one layer of elements
    !    given by the nieghbors.
    !    TO DO: identify repeated nodes using lninv_loc
    ! OUTPUT
    ! USED BY
    !    nsi_autint
    !***
    !-----------------------------------------------------------------------
    type neigeo
       integer(ip)          :: nelem_4
       integer(ip)          :: npoin_4
       integer(ip)          :: nboun_4
       integer(ip), pointer :: lnods_4(:,:) ! NELEM
       integer(ip), pointer :: ltype_4(:)   ! NELEM
       integer(ip), pointer :: lnnod_4(:)   ! NELEM
       integer(ip), pointer :: lelch_4(:)   ! NELEM
       integer(ip), pointer :: lesub_4(:)   ! NELEM
       integer(ip), pointer :: lmate_4(:)   ! NELEM
       integer(ip), pointer :: leinv_4(:)   ! NELEM
       integer(ip), pointer :: leldo_4(:)   ! NELEM
       integer(ip), pointer :: lbodo_4(:)   ! NBOUN
       integer(ip), pointer :: lnodb_4(:,:) ! NBOUN
       integer(ip), pointer :: ltypb_4(:)   ! NBOUN
       integer(ip), pointer :: lboch_4(:)   ! NBOUN
       integer(ip), pointer :: lboel_4(:,:) ! NBOUN
       integer(ip), pointer :: lelbo_4(:)   ! NBOUN
       integer(ip), pointer :: lninv_4(:)   ! NPOIN
       integer(ip), pointer :: lnoch_4(:)   ! NPOIN
       integer(ip), pointer :: lmast_4(:)   ! NPOIN
       real(rp),    pointer :: coord_4(:,:) ! NPOIN 
       integer(ip), pointer :: lslav_4(:,:) 
    end type neigeo
    type myygeo
       integer(ip), pointer :: leldo_3(:)   ! NELEM
       integer(ip), pointer :: lbodo_3(:)   ! NBOUN       
    end type myygeo
    integer(ip)          :: ipoin,inode,ielem,iboun,inodb,pnodb
    integer(ip)          :: pelty,jnode,jelem,jboun,idime,ii,kk
    integer(ip)          :: dom_i,jj,ineig,kneig,ielpo,jpoin,ll
    integer(ip)          :: pblty,pnode,jneig,nneig
    integer(ip)          :: kelem,kboun,kpoin
    integer(ip)          :: kelem_int,kboun_int,kpoin_int
    integer(ip)          :: kelem_rea,kboun_rea,kpoin_rea
    integer(ip), target  :: dumis(3),dumir(3)
    integer(ip), pointer :: lnewe(:)
    integer(ip), pointer :: lnewb(:)
    integer(ip), pointer :: lnewn(:)  
    integer(ip), pointer :: lnndo(:)  
    type(i1p),   pointer :: lnnei(:) 
    !
    ! My fringe geometry
    !
    integer(ip)          :: nelem_3
    integer(ip)          :: nboun_3
    integer(ip)          :: npoin_3
    integer(ip), pointer :: lnods_3(:,:) ! NELEM
    integer(ip), pointer :: ltype_3(:)   ! NELEM
    integer(ip), pointer :: lnnod_3(:)   ! NELEM
    integer(ip), pointer :: lelch_3(:)   ! NELEM
    integer(ip), pointer :: lesub_3(:)   ! NELEM
    integer(ip), pointer :: lmate_3(:)   ! NELEM
    integer(ip), pointer :: leinv_3(:)   ! NELEM
    integer(ip), pointer :: leldo_3(:)   ! NELEM
    integer(ip), pointer :: lnodb_3(:,:) ! NBOUN
    integer(ip), pointer :: ltypb_3(:)   ! NBOUN
    integer(ip), pointer :: lboch_3(:)   ! NBOUN
    integer(ip), pointer :: lboel_3(:,:) ! NBOUN
    integer(ip), pointer :: lelbo_3(:)   ! NBOUN
    integer(ip), pointer :: lninv_3(:)   ! NPOIN
    integer(ip), pointer :: lnoch_3(:)   ! NPOIN
    integer(ip), pointer :: lmast_3(:)   ! NPOIN
    real(rp),    pointer :: coord_3(:,:) ! NPOIN 
    !
    ! My neighboring subdomain geometry
    !
    type(neigeo),pointer :: negeo(:)
    integer(ip)          :: nelem_4
    integer(ip)          :: nboun_4
    integer(ip)          :: npoin_4
    integer(ip)          :: max_npoin_4
    !
    ! The geometry I send
    !
    type(myygeo),pointer :: mygeo(:)
    !
    ! Ghost element and boundary communicator
    !
    integer(ip), pointer :: list_real_halo(:)
    integer(ip), pointer :: numb_real_halo(:)
    integer(ip)          :: ii_send_boun_dim 
    integer(ip)          :: ii_recv_boun_dim 
    integer(ip)          :: ii_send_elem_dim 
    integer(ip)          :: ii_recv_elem_dim 
    integer(ip)          :: jj_send_boun_dim 
    integer(ip)          :: jj_recv_boun_dim 
    integer(ip)          :: jj_send_elem_dim 
    integer(ip)          :: jj_recv_elem_dim 
    !
    ! New merged mesh
    !
    integer(ip), pointer :: lninv_2(:)   ! NPOIN 
    integer(ip), pointer :: lnoch_2(:)   ! NPOIN 
    integer(ip), pointer :: lmast_2(:)   ! NPOIN 
    real(rp),    pointer :: coord_2(:,:) ! NPOIN 
    !
    ! Communication
    !
    integer(ip), pointer :: fring_perm(:)
    !
    ! Mesh renumbering
    !
    integer(ip), pointer :: permr(:)
    !
    ! Hash table
    !
    type(hash_t)         :: ht
    integer(ip)          :: lid
    integer(ip), pointer :: lpoins(:)
    integer(ip), pointer :: lneigs(:)
    logical(lg)          :: isin

    real(rp)             :: time1

    if( ISEQUEN ) return

    call livinf(0_ip,'PARALL: GET HALO GEOMETRY',0_ip)

    call cputim(time1)

    if( INOTMASTER .and. INOTEMPTY ) then
       !
       ! Nullify pointer
       !
       nullify(lnewe)
       nullify(lnewb)
       nullify(lnewn)  
       nullify(lnndo)  
       nullify(lnnei) 
       nullify(lnods_3) 
       nullify(ltype_3)   
       nullify(lnnod_3)   
       nullify(lelch_3)   
       nullify(lesub_3)   
       nullify(lmate_3)   
       nullify(leinv_3)   
       nullify(leldo_3)   
       nullify(lnodb_3) 
       nullify(ltypb_3)   
       nullify(lboch_3)   
       nullify(lboel_3) 
       nullify(lelbo_3) 
       nullify(lninv_3)   
       nullify(lnoch_3)   
       nullify(lmast_3)   
       nullify(coord_3) 
       nullify(negeo)
       nullify(mygeo)
       nullify(lninv_2)   
       nullify(lnoch_2)   
       nullify(lmast_2)   
       nullify(coord_2) 
       nullify(fring_perm)
       nullify(permr)

       nullify(list_real_halo)
       nullify(numb_real_halo)

       nullify(lpoins)
       nullify(lneigs)

       !----------------------------------------------------------------------
       !
       ! NEGEO(NNEIG): Structure for my neighboring subdomain geometry
       ! FRCOM(NNEIG): Communication structure
       !
       !----------------------------------------------------------------------

       nneig = commd % nneig    
       allocate( negeo(nneig) )
       do ineig = 1,nneig
          nullify( negeo(ineig) % lnods_4 )
          nullify( negeo(ineig) % ltype_4 )
          nullify( negeo(ineig) % lnnod_4 )
          nullify( negeo(ineig) % lelch_4 )
          nullify( negeo(ineig) % lesub_4 )
          nullify( negeo(ineig) % lmate_4 )
          nullify( negeo(ineig) % leinv_4 )
          nullify( negeo(ineig) % leldo_4 )
          nullify( negeo(ineig) % lnodb_4 )
          nullify( negeo(ineig) % ltypb_4 )
          nullify( negeo(ineig) % lboch_4 )
          nullify( negeo(ineig) % lboel_4 )
          nullify( negeo(ineig) % lelbo_4 )
          nullify( negeo(ineig) % lbodo_4 )
          nullify( negeo(ineig) % lninv_4 )
          nullify( negeo(ineig) % lnoch_4 )
          nullify( negeo(ineig) % lmast_4 )
          nullify( negeo(ineig) % coord_4 )
          nullify( negeo(ineig) % lslav_4 )
          negeo(ineig) % npoin_4 = 0
          negeo(ineig) % nelem_4 = 0
          negeo(ineig) % nboun_4 = 0
       end do
       allocate( mygeo(nneig) )
       do ineig = 1,nneig
          nullify( mygeo(ineig) % leldo_3 )
          nullify( mygeo(ineig) % lbodo_3 )       
       end do

       !----------------------------------------------------------------------
       !
       ! Neighboring subdomains for each node:
       ! LNNDO(IPOIN) ............ Number of neighboring subdomains
       ! LNNEIG(IPOIN) % L (:) ... List of neighboring subdomains
       ! 
       !----------------------------------------------------------------------

       call memory_alloca(par_memor,'LNNEI','mod_ghost_geometry',lnnei,npoin)
       call memory_alloca(par_memor,'LNNDO','mod_ghost_geometry',lnndo,npoin)

       do ipoin = 1,npoin
          lnndo(ipoin) = 0
       end do
       do ineig = 1,commd % nneig
          dom_i = commd % neights(ineig)
          if( dom_i /= kfl_paral ) then
             do jj = commd % bound_size(ineig),commd % bound_size(ineig+1)-1
                ipoin        = commd % bound_perm(jj)
                lnndo(ipoin) = lnndo(ipoin) + 1
             end do
          end if
       end do
       do ipoin = 1,npoin
          if( lnndo(ipoin) > 0 ) then
             call memory_alloca(par_memor,'LNNEI % L','mod_ghost_geometry',lnnei(ipoin) % l,lnndo(ipoin))
             lnndo(ipoin) = 0
          end if
       end do
       do ineig = 1,commd % nneig
          dom_i = commd % neights(ineig)
          if( dom_i /= kfl_paral ) then
             do jj = commd % bound_size(ineig),commd % bound_size(ineig+1)-1
                ipoin        = commd % bound_perm(jj)
                lnndo(ipoin) = lnndo(ipoin) + 1
                lnnei(ipoin) % l ( lnndo(ipoin) ) = ineig
             end do
          end if
       end do

       !----------------------------------------------------------------------
       !
       ! Create my fringe geometry with my neighbor INEIG
       ! 
       !----------------------------------------------------------------------
       !
       ! Allocate memory
       !
       call memory_alloca(par_memor,'lnewe','mod_ghost_geometry',lnewe,nelem)
       call memory_alloca(par_memor,'lnewb','mod_ghost_geometry',lnewb,nboun)
       call memory_alloca(par_memor,'lnewn','mod_ghost_geometry',lnewn,npoin) 
       !
       ! Get my fringe elements, boundary and nodes
       !
       do ineig = 1,commd % nneig

          dom_i         = commd % neights(ineig)

          if( dom_i /= kfl_paral ) then

             kfl_desti_par = dom_i

             do ielem = 1,nelem
                lnewe(ielem) = 0
             end do
             do iboun = 1,nboun
                lnewb(iboun) = 0
             end do
             do ipoin = 1,npoin
                lnewn(ipoin) = 0 
             end do

             do ipoin = npoi1+1,npoin

                kneig = 0
                do while( kneig < lnndo(ipoin) )
                   kneig = kneig + 1
                   if( lnnei(ipoin) % l(kneig) == ineig ) then
                      kneig = lnndo(ipoin) + 1
                   end if
                end do

                if( kneig == lnndo(ipoin) + 1 ) then
                   do ielpo = pelpo(ipoin),pelpo(ipoin+1)-1
                      jelem = lelpo(ielpo)
                      if( lnewe(jelem) == 0 ) then
                         lnewe(jelem) = 1
                         pelty = abs(ltype(jelem))
                         do jnode = 1,lnnod(jelem) !nnode(pelty)
                            jpoin = lnods(jnode,jelem) 
                            lnewn(jpoin) = 1
                         end do
                      end if
                   end do
                end if
             end do

             do iboun = 1,nboun
                pnodb = nnode(abs(ltypb(iboun)))
                ielem = lelbo(iboun)
                if( lnewe(ielem) > 0 ) then
                   lnewb(iboun) = 1
                end if
             end do
             !
             ! NELEM_3, NBOUN_3, NPOIN_3
             !
             nelem_3 = 0
             do ielem = 1,nelem
                if( lnewe(ielem) > 0 ) then
                   nelem_3      = nelem_3 + 1
                   lnewe(ielem) = nelem_3
                end if
             end do
             npoin_3 = 0
             do ipoin = 1,npoin
                if( lnewn(ipoin) > 0 ) then
                   npoin_3      = npoin_3 + 1
                   lnewn(ipoin) = npoin_3
                end if
             end do
             nboun_3 = 0
             do iboun = 1,nboun
                if( lnewb(iboun) > 0 ) then
                   nboun_3      = nboun_3 + 1
                   lnewb(iboun) = nboun_3
                end if
             end do
             !
             ! Allocate memory
             !
             call memory_alloca(par_memor,'lnods_3','mod_ghost_geometry',lnods_3,mnode,nelem_3) 
             call memory_alloca(par_memor,'ltype_3','mod_ghost_geometry',ltype_3,nelem_3)
             call memory_alloca(par_memor,'lnnod_3','mod_ghost_geometry',lnnod_3,nelem_3) 
             call memory_alloca(par_memor,'lelch_3','mod_ghost_geometry',lelch_3,nelem_3) 
             call memory_alloca(par_memor,'lesub_3','mod_ghost_geometry',lesub_3,nelem_3) 
             call memory_alloca(par_memor,'lmate_3','mod_ghost_geometry',lmate_3,nelem_3) 
             call memory_alloca(par_memor,'leinv_3','mod_ghost_geometry',leinv_3,nelem_3) 
             call memory_alloca(par_memor,'leldo_3','mod_ghost_geometry',leldo_3,nelem_3)

             call memory_alloca(par_memor,'lnodb_3','mod_ghost_geometry',lnodb_3,mnodb,nboun_3)             
             call memory_alloca(par_memor,'ltypb_3','mod_ghost_geometry',ltypb_3,nboun_3)                   
             call memory_alloca(par_memor,'lboch_3','mod_ghost_geometry',lboch_3,nboun_3)                   
             call memory_alloca(par_memor,'lboel_3','mod_ghost_geometry',lboel_3,mnodb,nboun_3)             
             call memory_alloca(par_memor,'lelbo_3','mod_ghost_geometry',lelbo_3,nboun_3)                 

             call memory_alloca(par_memor,'lninv_3','mod_ghost_geometry',lninv_3,npoin_3)      
             call memory_alloca(par_memor,'lnoch_3','mod_ghost_geometry',lnoch_3,npoin_3)      
             call memory_alloca(par_memor,'lmast_3','mod_ghost_geometry',lmast_3,npoin_3)      
             call memory_alloca(par_memor,'coord_3','mod_ghost_geometry',coord_3,ndime,npoin_3)
             !
             ! My send geometry
             !
             call memory_alloca(par_memor,'MYGEO % LELDO_3','mod_ghost_geometry',mygeo(ineig) % leldo_3,nelem_3)
             call memory_alloca(par_memor,'MYGEO % LBODO_3','mod_ghost_geometry',mygeo(ineig) % lbodo_3,nboun_3)
             !
             ! COORD_3 and LNINV_3
             !
             do ipoin = 1,npoin
                if( lnewn(ipoin) > 0 ) then
                   jpoin                  = lnewn(ipoin)
                   lninv_3(jpoin)         = lninv_loc(ipoin)
                   lnoch_3(jpoin)         = lnoch(ipoin)
                   lmast_3(jpoin)         = lmast(ipoin)
                   coord_3(1:ndime,jpoin) = coord(1:ndime,ipoin)
                end if
             end do
             !
             ! LNODS_3, LTYPE_3, LNNOD_3, LELCH_3, LESUB_3, LMATE_3, LEINV_3
             !
             do ielem = 1,nelem
                if( lnewe(ielem) > 0 ) then
                   jelem          = lnewe(ielem)
                   pelty          = abs(ltype(ielem))
                   ltype_3(jelem) = ltype(ielem)
                   lnnod_3(jelem) = lnnod(ielem)
                   lelch_3(jelem) = lelch(ielem)
                   lesub_3(jelem) = lesub(ielem)
                   lmate_3(jelem) = lmate(ielem)
                   leinv_3(jelem) = leinv_loc(ielem)
                   leldo_3(jelem) = ielem
                   mygeo(ineig) % leldo_3(jelem) = ielem
                   do inode = 1,lnnod(ielem)
                      lnods_3(inode,jelem) = lnewn(lnods(inode,ielem)) 
                   end do
                end if
             end do
             !
             ! LNODB_3, LTYPB_3, LBOCH_3, LBOEL_3, LELBO_3
             !
             do iboun = 1,nboun
                if( lnewb(iboun) > 0 ) then
                   jboun          = lnewb(iboun)
                   pblty          = abs(ltypb(iboun))
                   pnodb          = nnode(pblty)
                   ltypb_3(jboun) = ltypb(iboun)
                   lboch_3(jboun) = lboch(iboun)
                   mygeo(ineig) % lbodo_3(jboun) = iboun
                   do inodb = 1,pnodb
                      lnodb_3(inodb,jboun) = lnewn(lnodb(inodb,iboun)) 
                      lboel_3(inodb,jboun) = lboel(inodb,iboun) 
                   end do
                   lelbo_3(jboun) = lnewe(lelbo(iboun))
                end if
             end do
             !
             ! Send/recv geometry size to/from my neighbor INEIG
             ! Myself (NELEM_3,NBOUN_3,NPOIN_3) => 
             !                                  <= Neighbor (NELEM_4,NBOUN_4,NPOIN_4)
             ! 
             dumis(1)               =  nelem_3
             dumis(2)               =  nboun_3
             dumis(3)               =  npoin_3

             call PAR_SEND_RECEIVE(3_ip,3_ip,dumis,dumir,'IN MY CODE',dom_i)

             negeo(ineig) % nelem_4 = dumir(1)
             negeo(ineig) % nboun_4 = dumir(2)
             negeo(ineig) % npoin_4 = dumir(3)
             nelem_4                = negeo(ineig) % nelem_4
             nboun_4                = negeo(ineig) % nboun_4
             npoin_4                = negeo(ineig) % npoin_4
             !
             ! Allocate memory for my neighbor INEIG's geometry
             !
             call memory_alloca(par_memor,'negeo % lnods_4','mod_ghost_geometry',negeo(ineig) % lnods_4 ,mnode,nelem_4)
             call memory_alloca(par_memor,'negeo % ltype_4','mod_ghost_geometry',negeo(ineig) % ltype_4 ,nelem_4)      
             call memory_alloca(par_memor,'negeo % lnnod_4','mod_ghost_geometry',negeo(ineig) % lnnod_4 ,nelem_4)      
             call memory_alloca(par_memor,'negeo % lelch_4','mod_ghost_geometry',negeo(ineig) % lelch_4 ,nelem_4)      
             call memory_alloca(par_memor,'negeo % lesub_4','mod_ghost_geometry',negeo(ineig) % lesub_4 ,nelem_4)      
             call memory_alloca(par_memor,'negeo % lmate_4','mod_ghost_geometry',negeo(ineig) % lmate_4 ,nelem_4)      
             call memory_alloca(par_memor,'negeo % leinv_4','mod_ghost_geometry',negeo(ineig) % leinv_4 ,nelem_4)      
             call memory_alloca(par_memor,'negeo % leldo_4','mod_ghost_geometry',negeo(ineig) % leldo_4 ,nelem_4)      
             call memory_alloca(par_memor,'negeo % lnodb_4','mod_ghost_geometry',negeo(ineig) % lnodb_4 ,mnodb,nboun_4)
             call memory_alloca(par_memor,'negeo % ltypb_4','mod_ghost_geometry',negeo(ineig) % ltypb_4 ,nboun_4)      
             call memory_alloca(par_memor,'negeo % lboch_4','mod_ghost_geometry',negeo(ineig) % lboch_4 ,nboun_4)      
             call memory_alloca(par_memor,'negeo % lboel_4','mod_ghost_geometry',negeo(ineig) % lboel_4 ,mnodb,nboun_4)
             call memory_alloca(par_memor,'negeo % lbodo_4','mod_ghost_geometry',negeo(ineig) % lbodo_4 ,nboun_4)      
             call memory_alloca(par_memor,'negeo % lelbo_4','mod_ghost_geometry',negeo(ineig) % lelbo_4 ,nboun_4)      
             call memory_alloca(par_memor,'negeo % lninv_4','mod_ghost_geometry',negeo(ineig) % lninv_4 ,npoin_4)      
             call memory_alloca(par_memor,'negeo % lnoch_4','mod_ghost_geometry',negeo(ineig) % lnoch_4 ,npoin_4)      
             call memory_alloca(par_memor,'negeo % lmast_4','mod_ghost_geometry',negeo(ineig) % lmast_4 ,npoin_4)      
             call memory_alloca(par_memor,'negeo % coord_4','mod_ghost_geometry',negeo(ineig) % coord_4 ,ndime,npoin_4)
             !
             ! Send/recv geometry size to/from my neighbor INEIG
             !
             call PAR_SEND_RECEIVE( lnods_3 , negeo(ineig) % lnods_4 , 'IN MY CODE' , dom_i)
             call PAR_SEND_RECEIVE( ltype_3 , negeo(ineig) % ltype_4 , 'IN MY CODE' , dom_i)
             call PAR_SEND_RECEIVE( lnnod_3 , negeo(ineig) % lnnod_4 , 'IN MY CODE' , dom_i)
             call PAR_SEND_RECEIVE( lelch_3 , negeo(ineig) % lelch_4 , 'IN MY CODE' , dom_i)
             call PAR_SEND_RECEIVE( lesub_3 , negeo(ineig) % lesub_4 , 'IN MY CODE' , dom_i)
             call PAR_SEND_RECEIVE( lmate_3 , negeo(ineig) % lmate_4 , 'IN MY CODE' , dom_i)
             call PAR_SEND_RECEIVE( leinv_3 , negeo(ineig) % leinv_4 , 'IN MY CODE' , dom_i)
             call PAR_SEND_RECEIVE( leldo_3 , negeo(ineig) % leldo_4 , 'IN MY CODE' , dom_i) ! SACAR
             call PAR_SEND_RECEIVE( lnodb_3 , negeo(ineig) % lnodb_4 , 'IN MY CODE' , dom_i)
             call PAR_SEND_RECEIVE( ltypb_3 , negeo(ineig) % ltypb_4 , 'IN MY CODE' , dom_i)
             call PAR_SEND_RECEIVE( lboch_3 , negeo(ineig) % lboch_4 , 'IN MY CODE' , dom_i)
             call PAR_SEND_RECEIVE( lboel_3 , negeo(ineig) % lboel_4 , 'IN MY CODE' , dom_i)
             call PAR_SEND_RECEIVE( lelbo_3 , negeo(ineig) % lelbo_4 , 'IN MY CODE' , dom_i)
             call PAR_SEND_RECEIVE( lninv_3 , negeo(ineig) % lninv_4 , 'IN MY CODE' , dom_i)
             call PAR_SEND_RECEIVE( lnoch_3 , negeo(ineig) % lnoch_4 , 'IN MY CODE' , dom_i)
             call PAR_SEND_RECEIVE( lmast_3 , negeo(ineig) % lmast_4 , 'IN MY CODE' , dom_i)
             call PAR_SEND_RECEIVE( coord_3 , negeo(ineig) % coord_4 , 'IN MY CODE' , dom_i)

             call PAR_SEND_RECEIVE( mygeo(ineig) % leldo_3 , negeo(ineig) % leldo_4 , 'IN MY CODE' , dom_i)
             call PAR_SEND_RECEIVE( mygeo(ineig) % lbodo_3 , negeo(ineig) % lbodo_4 , 'IN MY CODE' , dom_i)
             !
             ! Deallocate memory
             !
             call memory_deallo(par_memor,'lnods_3','mod_ghost_geometry',lnods_3) 
             call memory_deallo(par_memor,'ltype_3','mod_ghost_geometry',ltype_3)
             call memory_deallo(par_memor,'lnnod_3','mod_ghost_geometry',lnnod_3) 
             call memory_deallo(par_memor,'lelch_3','mod_ghost_geometry',lelch_3) 
             call memory_deallo(par_memor,'lesub_3','mod_ghost_geometry',lesub_3) 
             call memory_deallo(par_memor,'lmate_3','mod_ghost_geometry',lmate_3) 
             call memory_deallo(par_memor,'leinv_3','mod_ghost_geometry',leinv_3) 
             call memory_deallo(par_memor,'leldo_3','mod_ghost_geometry',leldo_3)

             call memory_deallo(par_memor,'lnodb_3','mod_ghost_geometry',lnodb_3)             
             call memory_deallo(par_memor,'ltypb_3','mod_ghost_geometry',ltypb_3)                   
             call memory_deallo(par_memor,'lboch_3','mod_ghost_geometry',lboch_3)                   
             call memory_deallo(par_memor,'lboel_3','mod_ghost_geometry',lboel_3)             
             call memory_deallo(par_memor,'lelbo_3','mod_ghost_geometry',lelbo_3)   

             call memory_deallo(par_memor,'lninv_3','mod_ghost_geometry',lninv_3)      
             call memory_deallo(par_memor,'lnoch_3','mod_ghost_geometry',lnoch_3)      
             call memory_deallo(par_memor,'lmast_3','mod_ghost_geometry',lmast_3)      
             call memory_deallo(par_memor,'coord_3','mod_ghost_geometry',coord_3)

          end if

       end do

       call memory_deallo(par_memor,'LNNEI','mod_ghost_geometry',lnnei)
       call memory_deallo(par_memor,'LNNDO','mod_ghost_geometry',lnndo)

       call memory_deallo(par_memor,'par_fringe','COMMD % GHOST_SEND_BOUN_SIZE',lnnei)
       call memory_deallo(par_memor,'par_fringe','COMMD % GHOST_SEND_BOUN_SIZE',lnndo)

       !----------------------------------------------------------------------
       !
       ! Communication arrays: exchange global boundary number
       ! 
       !----------------------------------------------------------------------

       call memory_alloca(par_memor,'COMMD % GHOST_SEND_BOUN_SIZE','mod_ghost_geometry',commd % ghost_send_boun_size,commd % nneig+1_ip)
       call memory_alloca(par_memor,'COMMD % GHOST_RECV_BOUN_SIZE','mod_ghost_geometry',commd % ghost_recv_boun_size,commd % nneig+1_ip)
       call memory_alloca(par_memor,'COMMD % GHOST_SEND_ELEM_SIZE','mod_ghost_geometry',commd % ghost_send_elem_size,commd % nneig+1_ip)
       call memory_alloca(par_memor,'COMMD % GHOST_RECV_ELEM_SIZE','mod_ghost_geometry',commd % ghost_recv_elem_size,commd % nneig+1_ip)
       commd % ghost_send_boun_size(1) = 1
       commd % ghost_recv_boun_size(1) = 1
       commd % ghost_send_elem_size(1) = 1
       commd % ghost_recv_elem_size(1) = 1
       commd % ghost_send_boun_dim     = 0 
       commd % ghost_recv_boun_dim     = 0 
       commd % ghost_send_elem_dim     = 0 
       commd % ghost_recv_elem_dim     = 0 

       do ineig = 1,commd % nneig

          if( associated(mygeo(ineig) % lbodo_3) ) then
             ii_send_boun_dim = size( mygeo(ineig) % lbodo_3 , kind=ip )
          else
             ii_send_boun_dim = 0
          end if
          if( associated(negeo(ineig) % lbodo_4) ) then
             ii_recv_boun_dim = size( negeo(ineig) % lbodo_4 , kind=ip )
          else
             ii_recv_boun_dim = 0
          end if
          if( associated(mygeo(ineig) % leldo_3) ) then
             ii_send_elem_dim = size( mygeo(ineig) % leldo_3 , kind=ip )
          else
             ii_send_elem_dim = 0
          end if
          if( associated(negeo(ineig) % leldo_4) ) then
             ii_recv_elem_dim = size( negeo(ineig) % leldo_4 , kind=ip )
          else
             ii_recv_elem_dim = 0
          end if

          commd % ghost_send_boun_size(ineig+1) = commd % ghost_send_boun_size(ineig) + ii_send_boun_dim
          commd % ghost_recv_boun_size(ineig+1) = commd % ghost_recv_boun_size(ineig) + ii_recv_boun_dim
          commd % ghost_send_elem_size(ineig+1) = commd % ghost_send_elem_size(ineig) + ii_send_elem_dim
          commd % ghost_recv_elem_size(ineig+1) = commd % ghost_recv_elem_size(ineig) + ii_recv_elem_dim

       end do

       commd % ghost_send_boun_dim = commd % ghost_send_boun_size(commd % nneig + 1) - 1
       commd % ghost_recv_boun_dim = commd % ghost_recv_boun_size(commd % nneig + 1) - 1
       commd % ghost_send_elem_dim = commd % ghost_send_elem_size(commd % nneig + 1) - 1
       commd % ghost_recv_elem_dim = commd % ghost_recv_elem_size(commd % nneig + 1) - 1

       call memory_alloca(par_memor,'COMMD % GHOST_SEND_BOUN_PERM','mod_ghost_geometry',commd % ghost_send_boun_perm,commd % ghost_send_boun_dim)
       call memory_alloca(par_memor,'COMMD % GHOST_RECV_BOUN_PERM','mod_ghost_geometry',commd % ghost_recv_boun_perm,commd % ghost_recv_boun_dim)
       call memory_alloca(par_memor,'COMMD % GHOST_SEND_ELEM_PERM','mod_ghost_geometry',commd % ghost_send_elem_perm,commd % ghost_send_elem_dim)
       call memory_alloca(par_memor,'COMMD % GHOST_RECV_ELEM_PERM','mod_ghost_geometry',commd % ghost_recv_elem_perm,commd % ghost_recv_elem_dim)

       ii_send_boun_dim = 0
       ii_recv_boun_dim = 0
       ii_send_elem_dim = 0
       ii_recv_elem_dim = 0

       do ineig = 1,commd % nneig

          dom_i = commd % neights(ineig)

          if( dom_i /= kfl_paral ) then

             do jj_send_boun_dim = 1,commd % ghost_send_boun_size(ineig+1)-commd % ghost_send_boun_size(ineig)
                ii_send_boun_dim = ii_send_boun_dim + 1
                commd % ghost_send_boun_perm(ii_send_boun_dim) = mygeo(ineig) % lbodo_3(jj_send_boun_dim)
             end do

             do jj_recv_boun_dim = 1,commd % ghost_recv_boun_size(ineig+1)-commd % ghost_recv_boun_size(ineig)
                ii_recv_boun_dim = ii_recv_boun_dim + 1
                commd % ghost_recv_boun_perm(ii_recv_boun_dim) = nboun + ii_recv_boun_dim 
             end do

             do jj_send_elem_dim = 1,commd % ghost_send_elem_size(ineig+1)-commd % ghost_send_elem_size(ineig)
                ii_send_elem_dim = ii_send_elem_dim + 1
                commd % ghost_send_elem_perm(ii_send_elem_dim) = mygeo(ineig) % leldo_3(jj_send_elem_dim)
             end do

             do jj_recv_elem_dim = 1,commd % ghost_recv_elem_size(ineig+1)-commd % ghost_recv_elem_size(ineig)
                ii_recv_elem_dim = ii_recv_elem_dim + 1
                commd % ghost_recv_elem_perm(ii_recv_elem_dim) = nelem + ii_recv_elem_dim
             end do

          end if

       end do

       do ineig = 1,nneig
          call memory_deallo(par_memor,'MYGEO % LELDO_3','mod_ghost_geometry',mygeo(ineig) % leldo_3)
          call memory_deallo(par_memor,'MYGEO % LBODO_3','mod_ghost_geometry',mygeo(ineig) % lbodo_3)
       end do
       deallocate(mygeo)

       !----------------------------------------------------------------------
       !
       ! Merge my neighboring subdomain geometries with mine and themselves
       ! 
       !----------------------------------------------------------------------
       !
       ! Communication arrays: Fringe nodes communication arrays
       !
       call memory_alloca(par_memor,'COMMD % GHOST_RECV_NODE_SIZE','mod_ghost_geometry',commd % ghost_recv_node_size,commd % nneig+1)

       commd % ghost_recv_node_dim = 0
       do ineig = 1,commd % nneig
          npoin_4                             = negeo(ineig) % npoin_4
          commd % ghost_recv_node_size(ineig) = npoin_4
          commd % ghost_recv_node_dim         = commd % ghost_recv_node_dim + npoin_4
       end do
       commd % ghost_recv_node_size(commd%nneig+1) = 0

       call memory_alloca(par_memor,'FRING_PERM','mod_ghost_geometry',fring_perm,max(1_ip,commd % ghost_recv_node_dim))

       ii = 0
       do ineig = 1,commd % nneig
          
          do ipoin = 1,negeo(ineig) % npoin_4
             ii             = ii + 1
             jpoin          = negeo(ineig) % lninv_4(ipoin)
             fring_perm(ii) = jpoin
          end do
       end do
       !
       ! Merge my neighbors nodes
       !
       do ineig = 1,commd % nneig
          npoin_4 = negeo(ineig) % npoin_4
          call memory_alloca(par_memor,'negeo % lslav_4','mod_ghost_geometry',negeo(ineig) % lslav_4,2_ip,npoin_4) 
          do ipoin = 1,npoin_4
             negeo(ineig) % lslav_4(1,ipoin) = 0
             negeo(ineig) % lslav_4(2,ipoin) = 0
          end do
       end do
       !
       ! Check if JPOIN in JNEIG is IPOIN in INEIG 
       !
!!$       do ineig = 1,commd % nneig
!!$          npoin_4 = negeo(ineig) % npoin_4
!!$          do ipoin = 1,npoin_4
!!$             if( negeo(ineig) % lslav_4(1,ipoin) == 0 ) then
!!$                do jneig = ineig+1,commd%nneig
!!$                   do jpoin = 1,negeo(jneig) % npoin_4
!!$                      if( negeo(ineig) % lninv_4(ipoin) == negeo(jneig) % lninv_4(jpoin) ) then
!!$                         negeo(jneig) % lslav_4(1,jpoin) =  ineig 
!!$                         negeo(jneig) % lslav_4(2,jpoin) = -ipoin
!!$                      end if
!!$                   end do
!!$                end do
!!$             end if
!!$          end do
!!$       end do
       if( commd % nneig > 0 ) then
          max_npoin_4 = maxval(negeo(:) % npoin_4)
          call memory_alloca(par_memor,'LPOINS','par_ghost_geometry',lpoins,commd % nneig*max_npoin_4)
          call memory_alloca(par_memor,'LNEIGS','par_ghost_geometry',lneigs,commd % nneig*max_npoin_4)
          call htaini( ht,commd % nneig*max_npoin_4, lidson=.true., AUTOMATIC_SIZE=.true.,memor_opt=par_memor)
          do ineig = 1,commd % nneig
             do ipoin = 1,negeo(ineig) % npoin_4
                call htaadd(ht,int(negeo(ineig) % lninv_4(ipoin),ip), lid, isin)
                if(isin) then
                   negeo(ineig) % lslav_4(1,ipoin) = lneigs(lid)
                   negeo(ineig) % lslav_4(2,ipoin) = lpoins(lid)
                else
                   lneigs(lid) =  ineig
                   lpoins(lid) = -ipoin
                endif
             end do
          end do
          call htades(ht,memor_opt=par_memor)
          call memory_deallo(par_memor,'LPOINS','par_ghost_geometry',lpoins)
          call memory_deallo(par_memor,'LNEIGS','par_ghost_geometry',lneigs)
       end if
       !
       ! Merge my neighbors' nodes with mine
       !
!!$       do ineig = 1,commd % nneig
!!$          npoin_4 = negeo(ineig) % npoin_4
!!$          do ipoin = 1,npoin_4   
!!$             do jj = commd % bound_size(ineig),commd % bound_size(ineig+1)-1
!!$                jpoin = commd % bound_perm(jj)
!!$                if( negeo(ineig) % lninv_4(ipoin) == lninv_loc(jpoin) ) then
!!$                   negeo(ineig) % lslav_4(1,ipoin) = -1
!!$                   negeo(ineig) % lslav_4(2,ipoin) = -jpoin
!!$                   lid = htalid(htable_lninv_loc,int(negeo(ineig) % lninv_4(ipoin),ip))                   
!!$                   !if(kfl_paral==1) print*,'LID=',jpoin,lninv_loc(jpoin)      
!!$                end if
!!$             end do
!!$          end do
!!$       end do
       do ineig = 1,commd % nneig
          do ipoin = 1,negeo(ineig) % npoin_4
             lid = htalid(htable_lninv_loc,int(negeo(ineig) % lninv_4(ipoin),ip))
             if(lid>0) then
                negeo(ineig) % lslav_4(1,ipoin) = -1
                negeo(ineig) % lslav_4(2,ipoin) = -lid
             endif
          end do
       end do
       !
       ! Communication arrays: Definitive communication strategy: compact FRING_PERM to COMMD % FRING_PERM
       !
       commd % ghost_recv_node_dim = 0
       do ineig = 1,commd % nneig
          do ipoin = 1,negeo(ineig) % npoin_4
             if( negeo(ineig) % lslav_4(1,ipoin) == 0 ) then
                commd % ghost_recv_node_dim = commd % ghost_recv_node_dim + 1           
             end if
          end do
       end do

       call memory_alloca(par_memor,'COMMD % GHOST_RECV_NODE_PERM','mod_ghost_geometry',commd % ghost_recv_node_perm,commd % ghost_recv_node_dim)

       ii = 0
       do ineig = 1,commd%nneig
          commd % ghost_recv_node_size(ineig) = 0
          do ipoin = 1,negeo(ineig) % npoin_4
             if( negeo(ineig) % lslav_4(2,ipoin) == 0 ) then
                ii                                  = ii + 1
                commd % ghost_recv_node_size(ineig) = commd % ghost_recv_node_size(ineig) + 1
                fring_perm(ii)                      = negeo(ineig) % lninv_4(ipoin)
                commd % ghost_recv_node_perm(ii)    = npoin + ii
             end if
          end do
       end do
       !
       ! Communication arrays: Send my to my slave what I want from them
       !
       npasr = 0
       nparr = 0 
       call memory_alloca(par_memor,'COMMD % GHOST_SEND_NODE_SIZE','mod_ghost_geometry',commd % ghost_send_node_size,commd % nneig+1_ip)

       commd % ghost_send_node_dim = 0
       do ineig = 1,commd % nneig
          dom_i =  commd % neights(ineig)
          call PAR_SEND_RECEIVE(1_ip,1_ip,commd % ghost_recv_node_size(ineig:ineig),commd % ghost_send_node_size(ineig:ineig),'IN MY CODE',dom_i)
          commd % ghost_send_node_dim = commd % ghost_send_node_dim + commd % ghost_send_node_size(ineig)
       end do

       call memory_alloca(par_memor,'COMMD % GHOST_SEND_NODE_PERM','mod_ghost_geometry',commd % ghost_send_node_perm,max(1_ip,commd % ghost_send_node_dim))

       ii = 1
       kk = 1
       do ineig = 1,nneig
          dom_i =  commd % neights(ineig)
          if( kk > size(fring_perm,kind=ip) ) print*,'OUPS'
          if(      commd % ghost_recv_node_size(ineig) == 0 .and. commd % ghost_send_node_size(ineig) /= 0 ) then
             call PAR_SEND_RECEIVE(commd % ghost_recv_node_size(ineig),commd % ghost_send_node_size(ineig),dumir,commd % ghost_send_node_perm(ii:),'IN MY CODE',dom_i)
          else if( commd % ghost_recv_node_size(ineig) /= 0 .and. commd % ghost_send_node_size(ineig) == 0 ) then
             call PAR_SEND_RECEIVE(commd % ghost_recv_node_size(ineig),commd % ghost_send_node_size(ineig),fring_perm(kk:),dumir,'IN MY CODE',dom_i)
          else if( commd % ghost_recv_node_size(ineig) /= 0 .and. commd % ghost_send_node_size(ineig) /= 0 ) then
             call PAR_SEND_RECEIVE(commd % ghost_recv_node_size(ineig),commd % ghost_send_node_size(ineig),fring_perm(kk:),commd % ghost_send_node_perm(ii:),'IN MY CODE',dom_i)
          end if
          !call PAR_SEND_RECEIVE(commd % ghost_recv_node_size(ineig),commd % ghost_send_node_size(ineig),fring_perm(kk:),commd % ghost_send_node_perm(ii:),'IN MY CODE',dom_i)
          ii    = ii + commd % ghost_send_node_size(ineig)
          kk    = kk + commd % ghost_recv_node_size(ineig)
       end do

       call memory_deallo(par_memor,'FRING_PERM','mod_ghost_geometry',fring_perm)
       !
       ! Communication arrays: convert size COMMD % GHOST_RECV_NODE_SIZE to graph
       !
       ii = commd % ghost_recv_node_size(1)
       jj = commd % ghost_send_node_size(1)
       commd % ghost_recv_node_size(1) = 1
       commd % ghost_send_node_size(1) = 1
       do ineig = 1,commd % nneig
          kk = commd % ghost_recv_node_size(ineig+1)
          ll = commd % ghost_send_node_size(ineig+1)
          commd % ghost_recv_node_size(ineig+1) = commd % ghost_recv_node_size(ineig) + ii
          commd % ghost_send_node_size(ineig+1) = commd % ghost_send_node_size(ineig) + jj
          ii = kk
          jj = ll
       end do
       !
       ! Communication arrays: transform to local numbering
       !
       do ii = 1,commd % ghost_send_node_dim
          jpoin = commd % ghost_send_node_perm(ii)
          !ipoin = 0
          !do while( ipoin < npoin )
          !   ipoin = ipoin + 1
          !   if( lninv_loc(ipoin) == jpoin ) then
          !      commd % ghost_send_node_perm(ii) = ipoin
          !      ipoin = npoin + 2
          !   end if
          !end do
          ipoin = htalid(htable_lninv_loc,jpoin)
          if( ipoin > 0 )then
             commd % ghost_send_node_perm(ii) = ipoin
          else
             call runend('PAR_FRINGE: NODE NOT FOUND')
          end if
       end do
       !
       ! Renumber master nodes: they are NPOIN_2 new nodes
       !
       npoin_2 = 0
       do ineig = 1,commd%nneig
          npoin_4 = negeo(ineig) % npoin_4
          do ipoin = 1,npoin_4
             if( negeo(ineig) % lslav_4(1,ipoin) == 0 ) then
                npoin_2 = npoin_2 + 1
                negeo(ineig) % lslav_4(2,ipoin) = npoin_2 + npoin
             end if
          end do
       end do
       !
       ! Coordinates of newly created nodes
       !
       call memory_alloca(par_memor,'COORD_2','mod_ghost_geometry',coord_2,ndime,npoin_2) 
       call memory_alloca(par_memor,'LNINV_2','mod_ghost_geometry',lninv_2,npoin_2)       
       call memory_alloca(par_memor,'LNOCH_2','mod_ghost_geometry',lnoch_2,npoin_2)       
       call memory_alloca(par_memor,'LMAST_2','mod_ghost_geometry',lmast_2,npoin_2)       
       jpoin = 0
       do ineig = 1,commd%nneig
          npoin_4 = negeo(ineig) % npoin_4
          do ipoin = 1,npoin_4
             if( negeo(ineig) % lslav_4(1,ipoin) == 0 ) then
                jpoin = jpoin + 1
                lninv_2(jpoin) = negeo(ineig) % lninv_4(ipoin)
                lnoch_2(jpoin) = negeo(ineig) % lnoch_4(ipoin)
                lmast_2(jpoin) = negeo(ineig) % lmast_4(ipoin)
                do idime = 1,ndime
                   coord_2(idime,jpoin) = negeo(ineig) % coord_4(idime,ipoin)
                end do
             end if
          end do
       end do
       !
       ! Renumber slave nodes according to their slave
       !
       do jneig = 1,commd%nneig
          npoin_4 = negeo(jneig) % npoin_4
          do jpoin = 1,npoin_4

             if( negeo(jneig) % lslav_4(1,jpoin) > 0 ) then       ! JPOIN is a slave of another neighbor

                ineig =  negeo(jneig) % lslav_4(1,jpoin)
                ipoin = -negeo(jneig) % lslav_4(2,jpoin)
                kpoin =  negeo(ineig) % lslav_4(2,ipoin)          ! KPOIN global numbering of my master
                negeo(jneig) % lslav_4(2,jpoin) = kpoin

             else if( negeo(jneig) % lslav_4(1,jpoin) < 0 ) then  ! JPOIN is a slave of my own mesh

                ipoin = -negeo(jneig) % lslav_4(2,jpoin)
                negeo(jneig) % lslav_4(2,jpoin) = ipoin

             end if
          end do
       end do
       !
       ! Renumber element/boundary nodes
       !
       nelem_2 = 0
       nboun_2 = 0

       do ineig = 1,commd % nneig

          nelem_4 = negeo(ineig) % nelem_4
          nelem_2 = nelem_2 + nelem_4

          do ielem = 1,nelem_4

             pelty = negeo(ineig) % ltype_4(ielem)
             pnode = elmgeo_number_nodes(pelty,negeo(ineig) % lnods_4(:,ielem)) !nnode(pelty)
             do inode = 1,pnode
                ipoin = negeo(ineig) % lnods_4(inode,ielem)
                jpoin = negeo(ineig) % lslav_4(2,ipoin)
                negeo(ineig) % lnods_4(inode,ielem) = jpoin 
             end do

          end do

          nboun_4 = negeo(ineig) % nboun_4
          nboun_2 = nboun_2 + nboun_4

          do iboun = 1,nboun_4

             pblty = abs(negeo(ineig) % ltypb_4(iboun))
             pnodb = nnode(pblty)
             do inodb = 1,pnodb
                ipoin = negeo(ineig) % lnodb_4(inodb,iboun)
                jpoin = negeo(ineig) % lslav_4(2,ipoin)
                negeo(ineig) % lnodb_4(inodb,iboun) = jpoin 
             end do

          end do

       end do
       !
       ! Copy old mesh and reallocate
       !
       call memory_alloca(par_memor,'lnods_3','mod_ghost_geometry', lnods_3,mnode,nelem) 
       call memory_alloca(par_memor,'ltype_3','mod_ghost_geometry', ltype_3,nelem)       
       call memory_alloca(par_memor,'lnnod_3','mod_ghost_geometry', lnnod_3,nelem)       
       call memory_alloca(par_memor,'lelch_3','mod_ghost_geometry', lelch_3,nelem)       
       call memory_alloca(par_memor,'lesub_3','mod_ghost_geometry', lesub_3,nelem)       
       call memory_alloca(par_memor,'lmate_3','mod_ghost_geometry', lmate_3,nelem)       
       call memory_alloca(par_memor,'leinv_3','mod_ghost_geometry', leinv_3,nelem)       

       call memory_alloca(par_memor,'lnodb_3','mod_ghost_geometry', lnodb_3,mnodb,nboun) 
       call memory_alloca(par_memor,'ltypb_3','mod_ghost_geometry', ltypb_3,nboun)       
       call memory_alloca(par_memor,'lboch_3','mod_ghost_geometry', lboch_3,nboun)       
       call memory_alloca(par_memor,'lboel_3','mod_ghost_geometry', lboel_3,mnodb,nboun) 
       call memory_alloca(par_memor,'lelbo_3','mod_ghost_geometry', lelbo_3,nboun)       

       call memory_alloca(par_memor,'coord_3','mod_ghost_geometry', coord_3,ndime,npoin) 
       call memory_alloca(par_memor,'lninv_3','mod_ghost_geometry', lninv_3,npoin)       
       call memory_alloca(par_memor,'lnoch_3','mod_ghost_geometry', lnoch_3,npoin)       
       call memory_alloca(par_memor,'lmast_3','mod_ghost_geometry', lmast_3,npoin)       

       nelem_3 = nelem
       nboun_3 = nboun
       npoin_3 = npoin
       kelem   = 0
       kboun   = 0
       kpoin   = 0  
       call par_ghost_mshcpy(&
            1_ip      , &
            kelem     , kboun   , kpoin   , nelem_3 , nboun_3 , npoin_3 , &
            lnods_3   , ltype_3 , lnnod_3 , lelch_3 , lesub_3 , lmate_3 , leinv_3 , &
            lnodb_3   , ltypb_3 , lboch_3 , lboel_3 , lelbo_3 , coord_3 , lninv_3 , &
            lnoch_3   , lmast_3 , lnods   , ltype   , lnnod   , lelch   , lesub   , lmate   , &
            leinv_loc , lnodb   , ltypb   , lboch   , lboel   , lelbo   , coord   , &
            lninv_loc , lnoch   , lmast   )

       call domain_memory_deallocate('GEOMETRY')
       call domain_memory_deallocate('LBOEL')
       call domain_memory_deallocate('LNNOD')

       call memory_deallo(par_memor,'LNINV_LOC','mod_ghost_geometry',lninv_loc)
       call memory_deallo(par_memor,'LEINV_LOC','mod_ghost_geometry',leinv_loc)

       nelem = nelem + nelem_2
       nboun = nboun + nboun_2
       npoin = npoin + npoin_2

       call domain_memory_allocate('GEOMETRY') ! Geometrical arrays read                                 
       call domain_memory_allocate('LBOEL')    ! Geometrical arrays read                                 
       call domain_memory_allocate('LNNOD')    ! Geometrical arrays computed from read geometrical arrays

       call memory_alloca(par_memor,'LNINV_LOC','mod_ghost_geometry',lninv_loc,npoin)
       call memory_alloca(par_memor,'LEINV_LOC','mod_ghost_geometry',leinv_loc,nelem)

       kelem = 0
       kboun = 0
       kpoin = 0

       call par_ghost_mshcpy(&
            1_ip    , &
            kelem   , kboun   , kpoin   , nelem_3 , nboun_3 , npoin_3   , &
            lnods   , ltype   , lnnod   , lelch   , lesub   , lmate     , leinv_loc , &
            lnodb   , ltypb   , lboch   , lboel   , lelbo   , coord   , lninv_loc , &
            lnoch   , lmast   , lnods_3 , ltype_3 , lnnod_3 , lelch_3 , lesub_3   , lmate_3   , &
            leinv_3 , lnodb_3 , ltypb_3 , lboch_3 , lboel_3 , lelbo_3 , coord_3   , &
            lninv_3 , lnoch_3 , lmast_3 )

       call memory_deallo(par_memor,'lnods_3','mod_ghost_geometry', lnods_3) 
       call memory_deallo(par_memor,'ltype_3','mod_ghost_geometry', ltype_3)       
       call memory_deallo(par_memor,'lnnod_3','mod_ghost_geometry', lnnod_3)       
       call memory_deallo(par_memor,'lelch_3','mod_ghost_geometry', lelch_3)       
       call memory_deallo(par_memor,'lesub_3','mod_ghost_geometry', lesub_3)       
       call memory_deallo(par_memor,'lmate_3','mod_ghost_geometry', lmate_3)       
       call memory_deallo(par_memor,'leinv_3','mod_ghost_geometry', leinv_3)       

       call memory_deallo(par_memor,'lnodb_3','mod_ghost_geometry', lnodb_3) 
       call memory_deallo(par_memor,'ltypb_3','mod_ghost_geometry', ltypb_3)       
       call memory_deallo(par_memor,'lboch_3','mod_ghost_geometry', lboch_3)       
       call memory_deallo(par_memor,'lboel_3','mod_ghost_geometry', lboel_3) 
       call memory_deallo(par_memor,'lelbo_3','mod_ghost_geometry', lelbo_3)       

       call memory_deallo(par_memor,'coord_3','mod_ghost_geometry', coord_3) 
       call memory_deallo(par_memor,'lninv_3','mod_ghost_geometry', lninv_3)       
       call memory_deallo(par_memor,'lnoch_3','mod_ghost_geometry', lnoch_3)       
       call memory_deallo(par_memor,'lmast_3','mod_ghost_geometry', lmast_3)       
       !
       ! LELDO
       !
       call memory_alloca(par_memor,'LELDO','mod_ghost_geometry',leldo,2_ip,nelem_2) ! NELEM
       jelem = 0 
       do ineig = 1,commd%nneig
          do ielem = 1,negeo(ineig) % nelem_4
             jelem = jelem + 1
             leldo(1,jelem) = ineig
             leldo(2,jelem) = negeo(ineig) % leldo_4(ielem)
          end do
       end do
       !
       ! Merge my neighbors' mesh
       !
       !
       ! LNINV, LNOCH, LMAST, COORD
       !  
       kpoin = npoin_3
       do ipoin = 1,npoin_2
          kpoin = kpoin + 1
          lninv_loc(kpoin) = lninv_2(ipoin)
          lnoch(kpoin)     = lnoch_2(ipoin)
          lmast(kpoin)     = lmast_2(ipoin)
          do idime = 1,ndime
             coord(idime,kpoin) = coord_2(idime,ipoin)
          end do
       end do

       do ineig = 1,commd%nneig
          lnods_3 => negeo(ineig) % lnods_4
          ltype_3 => negeo(ineig) % ltype_4
          lnnod_3 => negeo(ineig) % lnnod_4
          lelch_3 => negeo(ineig) % lelch_4
          lesub_3 => negeo(ineig) % lesub_4
          lmate_3 => negeo(ineig) % lmate_4
          leinv_3 => negeo(ineig) % leinv_4
          lnodb_3 => negeo(ineig) % lnodb_4
          ltypb_3 => negeo(ineig) % ltypb_4
          lboch_3 => negeo(ineig) % lboch_4
          lboel_3 => negeo(ineig) % lboel_4
          lelbo_3 => negeo(ineig) % lelbo_4
          coord_3 => negeo(ineig) % coord_4
          lninv_3 => negeo(ineig) % lninv_4
          lnoch_3 => negeo(ineig) % lnoch_4
          lmast_3 => negeo(ineig) % lmast_4
          nelem_3 =  negeo(ineig) % nelem_4
          nboun_3 =  negeo(ineig) % nboun_4
          npoin_3 =  negeo(ineig) % npoin_4 

          call par_ghost_mshcpy(&
               0_ip    , &
               kelem   , kboun   , kpoin   , nelem_3 , nboun_3 , npoin_3   , &
               lnods   , ltype   , lnnod   , lelch   , lesub   , lmate     , leinv_loc , &
               lnodb   , ltypb   , lboch   , lboel   , lelbo   , coord   , lninv_loc , &
               lnoch   , lmast   , lnods_3 , ltype_3 , lnnod_3 , lelch_3 , lesub_3   , lmate_3   , &
               leinv_3 , lnodb_3 , ltypb_3 , lboch_3 , lboel_3 , lelbo_3 , coord_3   , &
               lninv_3 , lnoch_3 , lmast_3 )

       end do
       !
       ! Deallocate memory
       !
       call memory_deallo(par_memor,'lnewe','mod_ghost_geometry',lnewe)
       call memory_deallo(par_memor,'lnewb','mod_ghost_geometry',lnewb)
       call memory_deallo(par_memor,'lnewn','mod_ghost_geometry',lnewn) 

       do ineig = 1,commd % nneig

          call memory_deallo(par_memor,'negeo % lnods_4','mod_ghost_geometry',negeo(ineig) % lnods_4 )
          call memory_deallo(par_memor,'negeo % ltype_4','mod_ghost_geometry',negeo(ineig) % ltype_4 )      
          call memory_deallo(par_memor,'negeo % lnnod_4','mod_ghost_geometry',negeo(ineig) % lnnod_4 )      
          call memory_deallo(par_memor,'negeo % lelch_4','mod_ghost_geometry',negeo(ineig) % lelch_4 )      
          call memory_deallo(par_memor,'negeo % lesub_4','mod_ghost_geometry',negeo(ineig) % lesub_4 )      
          call memory_deallo(par_memor,'negeo % lmate_4','mod_ghost_geometry',negeo(ineig) % lmate_4 )      
          call memory_deallo(par_memor,'negeo % leinv_4','mod_ghost_geometry',negeo(ineig) % leinv_4 )      
          call memory_deallo(par_memor,'negeo % leldo_4','mod_ghost_geometry',negeo(ineig) % leldo_4 )

          call memory_deallo(par_memor,'negeo % lnodb_4','mod_ghost_geometry',negeo(ineig) % lnodb_4 )
          call memory_deallo(par_memor,'negeo % ltypb_4','mod_ghost_geometry',negeo(ineig) % ltypb_4 )      
          call memory_deallo(par_memor,'negeo % lboch_4','mod_ghost_geometry',negeo(ineig) % lboch_4 )      
          call memory_deallo(par_memor,'negeo % lboel_4','mod_ghost_geometry',negeo(ineig) % lboel_4 )
          call memory_deallo(par_memor,'negeo % lbodo_4','mod_ghost_geometry',negeo(ineig) % lbodo_4 )      
          call memory_deallo(par_memor,'negeo % lelbo_4','mod_ghost_geometry',negeo(ineig) % lelbo_4 )

          call memory_deallo(par_memor,'negeo % lninv_4','mod_ghost_geometry',negeo(ineig) % lninv_4 )      
          call memory_deallo(par_memor,'negeo % lnoch_4','mod_ghost_geometry',negeo(ineig) % lnoch_4 )      
          call memory_deallo(par_memor,'negeo % lmast_4','mod_ghost_geometry',negeo(ineig) % lmast_4 )      
          call memory_deallo(par_memor,'negeo % coord_4','mod_ghost_geometry',negeo(ineig) % coord_4 )
          call memory_deallo(par_memor,'negeo % lslav_4','mod_ghost_geometry',negeo(ineig) % lslav_4 )

       end do
       if( associated( negeo   ) ) deallocate( negeo   )

       call memory_deallo(par_memor,'COORD_2','mod_ghost_geometry',coord_2) 
       call memory_deallo(par_memor,'LNINV_2','mod_ghost_geometry',lninv_2)       
       call memory_deallo(par_memor,'LNOCH_2','mod_ghost_geometry',lnoch_2)       
       call memory_deallo(par_memor,'LMAST_2','mod_ghost_geometry',lmast_2)       

       nelem   = nelem - nelem_2
       nboun   = nboun - nboun_2
       npoin   = npoin - npoin_2
       nelem_2 = nelem + nelem_2
       nboun_2 = nboun + nboun_2
       npoin_2 = npoin + npoin_2

       !----------------------------------------------------------------------
       !
       ! Test of the ghost element and node exchange
       ! 
       !----------------------------------------------------------------------

       kelem_int = 0
       kboun_int = 0
       kpoin_int = 0
       kelem_rea = 0
       kboun_rea = 0
       kpoin_rea = 0
       !
       ! Integer on elements
       !
       call memgen(1_ip,nelem_2,0_ip)
       do ielem = 1,nelem
          gisca(ielem) = kfl_paral
       end do
       call PAR_GHOST_ELEMENT_EXCHANGE(gisca,'SUBSTITUTE','IN MY CODE')
       do ineig = 1,commd % nneig
          dom_i = commd % neights(ineig)
          do ii = commd % ghost_recv_elem_size(ineig),commd % ghost_recv_elem_size(ineig+1)-1
             jelem = commd % ghost_recv_elem_perm(ii)
             if( gisca(jelem) /= dom_i ) kelem_int = kelem_int + 1
          end do
       end do
       call memgen(3_ip,nelem_2,0_ip)
       !
       ! Integer on boundaries
       !
       call memgen(1_ip,nboun_2,0_ip)
       do iboun = 1,nboun_2
          gisca(iboun) = kfl_paral
       end do
       call PAR_GHOST_BOUNDARY_EXCHANGE(gisca,'SUBSTITUTE','IN MY CODE')
       do ineig = 1,commd % nneig
          dom_i = commd % neights(ineig)
          do ii = commd % ghost_recv_boun_size(ineig),commd % ghost_recv_boun_size(ineig+1)-1
             jboun = commd % ghost_recv_boun_perm(ii)
             if( gisca(jboun) /= dom_i ) kboun_int = kboun_int + 1
          end do
       end do
       call memgen(3_ip,nboun_2,0_ip)
       !
       ! Integer on nodes
       !
       call memgen(1_ip,npoin_2,0_ip)
       do ipoin = 1,npoin_2
          gisca(ipoin) = kfl_paral
       end do
       call PAR_GHOST_NODE_EXCHANGE(gisca,'SUBSTITUTE','IN MY CODE')
       do ineig = 1,commd % nneig
          dom_i = commd % neights(ineig)
          do ii = commd % ghost_recv_node_size(ineig),commd % ghost_recv_node_size(ineig+1)-1
             jpoin = commd % ghost_recv_node_perm(ii)
             if( gisca(jpoin) /= dom_i ) kboun_int = kboun_int + 1
          end do
       end do
       call memgen(3_ip,npoin_2,0_ip)
       !
       ! Real on elements
       !
       call memgen(0_ip,nelem_2,0_ip)
       do ielem = 1,nelem
          gesca(ielem) = real(kfl_paral,rp)
       end do
       call PAR_GHOST_ELEMENT_EXCHANGE(gesca,'SUBSTITUTE','IN MY CODE')
       do ineig = 1,commd % nneig
          dom_i = commd % neights(ineig)
          do ii = commd % ghost_recv_elem_size(ineig),commd % ghost_recv_elem_size(ineig+1)-1
             jelem = commd % ghost_recv_elem_perm(ii)
             if( int(gesca(jelem),ip) /= dom_i ) kelem_rea = kelem_rea + 1
          end do
       end do
       call memgen(2_ip,nelem_2,0_ip)
       !
       ! Real on boundaries
       !
       call memgen(0_ip,nboun_2,0_ip)
       do iboun = 1,nboun_2
          gesca(iboun) = real(kfl_paral,rp)
       end do
       call PAR_GHOST_BOUNDARY_EXCHANGE(gesca,'SUBSTITUTE','IN MY CODE')
       do ineig = 1,commd % nneig
          dom_i = commd % neights(ineig)
          do ii = commd % ghost_recv_boun_size(ineig),commd % ghost_recv_boun_size(ineig+1)-1
             jboun = commd % ghost_recv_boun_perm(ii)
             if( int(gesca(jboun),ip) /= dom_i ) kboun_rea = kboun_rea + 1
          end do
       end do
       call memgen(2_ip,nboun_2,0_ip)
       !
       ! Real on nodes
       !
       call memgen(0_ip,npoin_2,0_ip)
       do ipoin = 1,npoin_2
          gesca(ipoin) = real(kfl_paral,rp)
       end do
       call PAR_GHOST_NODE_EXCHANGE(gesca,'SUBSTITUTE','IN MY CODE')
       do ineig = 1,commd % nneig
          dom_i = commd % neights(ineig)
          do ii = commd % ghost_recv_node_size(ineig),commd % ghost_recv_node_size(ineig+1)-1
             jpoin = commd % ghost_recv_node_perm(ii)
             if( int(gesca(jpoin),ip) /= dom_i ) kboun_rea = kboun_rea + 1
          end do
       end do
       call memgen(2_ip,npoin_2,0_ip)

    else

       nelem_2   = 0
       nboun_2   = 0
       npoin_2   = 0
       kelem_int = 0
       kboun_int = 0
       kpoin_int = 0
       kelem_rea = 0
       kboun_rea = 0
       kpoin_rea = 0

    end if
    !
    ! Check errors
    !
    call PAR_SUM(kelem_int,'IN MY CODE')
    call PAR_SUM(kboun_int,'IN MY CODE')
    call PAR_SUM(kpoin_int,'IN MY CODE')
    call PAR_SUM(kelem_rea,'IN MY CODE')
    call PAR_SUM(kboun_rea,'IN MY CODE')
    call PAR_SUM(kpoin_rea,'IN MY CODE')
    if( kelem_int /= 0 ) call runend('PAR_GHOST_GEOMETRY: WRONG GHOST ELEMENT COMMUNICATOR FOR INT:   '//trim(intost(kelem_int)))
    if( kboun_int /= 0 ) call runend('PAR_GHOST_GEOMETRY: WRONG GHOST BOUNDARY COMMUNICATOR FOR INT:  '//trim(intost(kboun_int)))
    if( kpoin_int /= 0 ) call runend('PAR_GHOST_GEOMETRY: WRONG GHOST NODE COMMUNICATOR FOR INT:      '//trim(intost(kpoin_int)))
    if( kelem_rea /= 0 ) call runend('PAR_GHOST_GEOMETRY: WRONG GHOST ELEMENT COMMUNICATOR FOR REAL:  '//trim(intost(kelem_rea)))
    if( kboun_rea /= 0 ) call runend('PAR_GHOST_GEOMETRY: WRONG GHOST BOUNDARY COMMUNICATOR FOR REAL: '//trim(intost(kboun_rea)))
    if( kpoin_rea /= 0 ) call runend('PAR_GHOST_GEOMETRY: WRONG GHOST NODE COMMUNICATOR FOR REAL:     '//trim(intost(kpoin_rea)))
    !
    ! Exchange additional arrays 
    !
    if( INOTEMPTY ) then

       call memory_resize(memor_dom,'LNNOB','mod_ghost_geometry',lnnob,nboun_2)
       call PAR_GHOST_BOUNDARY_EXCHANGE(lnnob,'SUBSTITUTE','IN MY CODE')
       call memory_resize(memor_dom,'LBINV_LOC','mod_ghost_geometry',lbinv_loc,nboun_2)
       call PAR_GHOST_BOUNDARY_EXCHANGE(lbinv_loc,'SUBSTITUTE','IN MY CODE')
       if( neset > 0 ) then
          call memory_resize(memor_dom,'LESET','mod_ghost_geometry',leset,nelem_2)
          call PAR_GHOST_ELEMENT_EXCHANGE(leset,'SUBSTITUTE','IN MY CODE')
       end if
       if( nbset > 0 ) then
          call memory_resize(memor_dom,'LBSET','mod_ghost_geometry',lbset,nboun_2)
          call PAR_GHOST_BOUNDARY_EXCHANGE(lbset,'SUBSTITUTE','IN MY CODE')
       end if

    end if
    !
    ! Renumber halo nodes so that halo nodes connected to my own nodes (interior+own boundary) are numbered first
    !
    if( INOTEMPTY ) then
       call memory_alloca(par_memor,'PERMR_GHOST','mod_ghost_geometry',permr,npoin_2)
       do ielem = nelem+1,nelem_2
          loop_inode: do inode = 1,lnnod(ielem)
             ipoin = lnods(inode,ielem)
             if( ipoin <= npoi1 .or. ( ipoin >= npoi2 .and. ipoin <= npoi3 ) ) then
                do jnode = 1,lnnod(ielem)
                   jpoin = lnods(jnode,ielem)
                   if( jpoin > npoin ) permr(jpoin) = 1
                end do
                exit loop_inode                
             end if
          end do loop_inode
       end do
       do ipoin = 1,npoin
          permr(ipoin) = ipoin 
       end do
       jpoin = npoin
       do ipoin = npoin+1,npoin_2
          if( permr(ipoin) == 1 ) then
             jpoin = jpoin + 1
             permr(ipoin) = jpoin
          end if
       end do
       npoin_halo = jpoin
       do ipoin = npoin+1,npoin_2
          if( permr(ipoin) == 0 ) then
             jpoin = jpoin + 1
             permr(ipoin) = jpoin
          end if
       end do
       call renumbering_node_arrays(permr,npoin_2,nelem_2,nboun_2)
       call memory_deallo(par_memor,'PERMR_GHOST','mod_ghost_geometry',permr)

!!$       if( kfl_paral == 1 ) then
!!$          open(unit=90,file='subdo.post.msh',status='unknown')
!!$          open(unit=91,file='subdo.post.res',status='unknown')
!!$          write(90,*) 'MESH SUBDOMAIN dimension 2 Elemtype Triangle Nnode 3'
!!$          write(90,*) 'coordinates'
!!$          do ipoin = 1,npoin_2
!!$             write(90,'(i7,1x,3(1x,e13.6))') ipoin,coord(1:ndime,ipoin)
!!$          end do
!!$          write(90,*) 'end coordinates'
!!$          write(90,*) 'elements'
!!$          do ielem = 1,nelem
!!$             write(90,'(10(1x,i7))') ielem,lnods(:,ielem),1
!!$          end do
!!$          do ielem = nelem+1,nelem_2
!!$             write(90,'(10(1x,i7))') ielem,lnods(:,ielem),2
!!$          end do
!!$          write(90,*) 'end elements'
!!$          write(91,*) 'GiD Post Results File 1.0'
!!$          write(91,*) ' '
!!$          write(91,*) 'Result TYPE ALYA 0.0 Scalar OnNodes'
!!$          write(91,*) 'ComponentNames TYPE'
!!$          write(91,*) 'values'
!!$           do ipoin = 1,npoi1
!!$             write(91,'(i7,1x,3(1x,e13.6))') ipoin,1.0_rp
!!$          end do
!!$           do ipoin = npoi2,npoi3
!!$             write(91,'(i7,1x,3(1x,e13.6))') ipoin,2.0_rp
!!$          end do
!!$          do ipoin = npoi3+1,npoin
!!$             write(91,'(i7,1x,3(1x,e13.6))') ipoin,3.0_rp
!!$          end do
!!$          do ipoin = npoin+1,npoin_2
!!$             write(91,'(i7,1x,3(1x,e13.6))') ipoin,4.0_rp
!!$          end do         
!!$          write(91,*) 'end values'
!!$          close(90)
!!$          close(91)
!!$       end if
    end if
    !
    ! Increase hash table to take halo nodes into account
    !
    if( INOTEMPTY ) then       
       if(npoin_2 > npoin) then
          call htaadd( htable_lninv_loc,  npoin_2-npoin, lninv_loc(npoin+1:npoin_2), memor_opt=memor_dom)
       end if
       if(htable_lninv_loc % nelem /= npoin_2 ) then 
          call runend("domain : repeated elements in lninv_loc")
       end if
    end if
    !
    ! Repoint postprocess mesh, as variables have been deallocated
    !
    call mesh_type_update_last_mesh()

  end subroutine par_ghost_geometry

  subroutine par_ghost_mshcpy( &
       itask     ,                                                   &
       kelem     , kboun   , kpoin   , nelem   , nboun   , npoin   , &
       lnods_3   , ltype_3 , lnnod_3 , lelch_3 , lesub_3 , lmate_3 , leinv_3 , &
       lnodb_3   , ltypb_3 , lboch_3 , lboel_3 , lelbo_3 , coord_3 , lninv_3 , &
       lnoch_3   , lmast_3 , lnods   , ltype   , lnnod   , lelch   , lesub   , lmate   , &
       leinv_loc , lnodb   , ltypb   , lboch   , lboel   , lelbo   , coord   , &
       lninv     , lnoch   , lmast   )
    !-----------------------------------------------------------------------
    !****f* Parall/mshcpy
    ! NAME
    !    mshcpy
    ! DESCRIPTION
    !    This routine copies a mesh over another
    ! OUTPUT
    ! USED BY
    !    nsi_autint
    !***
    !-----------------------------------------------------------------------
    
    integer(ip), intent(in)             :: itask
    integer(ip), intent(inout)          :: kelem,kboun,kpoin
    integer(ip), intent(in)             :: nelem,nboun,npoin
    integer(ip), intent(inout), pointer :: lnods_3(:,:)
    integer(ip), intent(inout), pointer :: ltype_3(:)
    integer(ip), intent(inout), pointer :: lnnod_3(:)
    integer(ip), intent(inout), pointer :: lelch_3(:)
    integer(ip), intent(inout), pointer :: lesub_3(:)
    integer(ip), intent(inout), pointer :: lmate_3(:)
    integer(ip), intent(inout), pointer :: leinv_3(:)
    integer(ip), intent(inout), pointer :: lnodb_3(:,:)
    integer(ip), intent(inout), pointer :: ltypb_3(:) 
    integer(ip), intent(inout), pointer :: lboch_3(:)
    integer(ip), intent(inout), pointer :: lboel_3(:,:)
    integer(ip), intent(inout), pointer :: lelbo_3(:)
    real(rp),    intent(inout), pointer :: coord_3(:,:)
    integer(ip), intent(inout), pointer :: lninv_3(:)
    integer(ip), intent(inout), pointer :: lnoch_3(:)
    integer(ip), intent(inout), pointer :: lmast_3(:)
    integer(ip), intent(in),    pointer :: lnods(:,:)
    integer(ip), intent(in),    pointer :: ltype(:)
    integer(ip), intent(in),    pointer :: lnnod(:)
    integer(ip), intent(in),    pointer :: lelch(:)
    integer(ip), intent(in),    pointer :: lesub(:)
    integer(ip), intent(in),    pointer :: lmate(:)
    integer(ip), intent(in),    pointer :: leinv_loc(:)
    integer(ip), intent(in),    pointer :: lnodb(:,:)
    integer(ip), intent(in),    pointer :: ltypb(:)
    integer(ip), intent(in),    pointer :: lboch(:)
    integer(ip), intent(in),    pointer :: lboel(:,:)
    integer(ip), intent(in),    pointer :: lelbo(:)
    real(rp),    intent(in),    pointer :: coord(:,:)
    integer(ip), intent(in),    pointer :: lninv(:)
    integer(ip), intent(in),    pointer :: lnoch(:)
    integer(ip), intent(in),    pointer :: lmast(:)
    integer(ip)                         :: inode,ipoin,ielem,iboun,pelty,pblty
    integer(ip)                         :: idime,inodb,pnodb
    !
    ! LNODS_3, LTYPE_3, LNNOD_3, LELCH_3, LESUB_3, LMATE_3, LEINV_3
    !
    do ielem = 1,nelem
       kelem = kelem + 1
       pelty = abs(ltype(ielem))
       ltype_3(kelem) = ltype(ielem)
       lnnod_3(kelem) = lnnod(ielem)
       lelch_3(kelem) = lelch(ielem)
       lesub_3(kelem) = lesub(ielem)
       lmate_3(kelem) = lmate(ielem)
       leinv_3(kelem) = leinv_loc(ielem)
       do inode = 1,lnnod(ielem) !nnode(pelty)
          lnods_3(inode,kelem) = lnods(inode,ielem)
       end do
    end do
    !
    ! LNODB_3, LTYPB_3, LBOCH_3, LBOEL_3, LELBO_3
    !
    do iboun = 1,nboun
       kboun = kboun + 1
       pblty = abs(ltypb(iboun))
       pnodb = nnode(pblty)
       ltypb_3(kboun) = ltypb(iboun)
       lboch_3(kboun) = lboch(iboun)
       do inodb = 1,pnodb
          lnodb_3(inodb,kboun) = lnodb(inodb,iboun) 
          lboel_3(inodb,kboun) = lboel(inodb,iboun) 
       end do
       lelbo_3(kboun) = lelbo(iboun)
       if( itask == 0 ) lelbo_3(kboun) = lelbo(iboun) + kelem - nelem
    end do
    !
    ! COORD_3, LNINV_3, LNOCH_3, LMAST_3
    !
    if( itask == 1 ) then
       do ipoin = 1,npoin
          kpoin = kpoin + 1
          lninv_3(kpoin) = lninv(ipoin)
          lnoch_3(kpoin) = lnoch(ipoin)
          lmast_3(kpoin) = lmast(ipoin)
          do idime = 1,ndime
             coord_3(idime,kpoin) = coord(idime,ipoin)
          end do
       end do
    end if

  end subroutine par_ghost_mshcpy

end module mod_ghost_geometry
