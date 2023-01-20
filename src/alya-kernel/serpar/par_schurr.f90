!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine par_schurr()
  !------------------------------------------------------------------------
  !****f* Parall/par_schurr
  ! NAME 
  !    par_schurr
  ! DESCRIPTION
  !    This routine constructs the IB graphs for Schur complement 
  !    type solvers
  ! USES
  ! USED BY 
  !***
  !------------------------------------------------------------------------
  use def_kintyp
  use def_master
  use def_domain
  use def_solver
  use def_parall
  use mod_memchk
  use mod_graphs
  use mod_parall, only : PAR_INTEGER
  use mod_parall, only : PAR_COMM_MY_CODE,commd
  use mod_parall, only : par_memor
  use def_mpi
  implicit none
  
  integer(ip)             :: ipoin,ii,kpoin,izdom,jpoin,kbb
  integer(ip)             :: ibb,kk,ll,jj,ks,kr,qpoin,isize
  integer(ip)             :: ifoun,spoin,kzdom,pzdom,jzdom
  integer(ip)             :: nsolv,isolv,imodu,dummi
  integer                 :: istat,bsizs4,bsizr4,dom_i
  integer(ip),  pointer   :: gisnd(:,:) 
  integer(ip),  pointer   :: gircv(:,:) 
  integer(ip),  pointer   :: lcorn(:) 
  type(i1p),    pointer   :: type_prec(:)
  integer(ip),  pointer   :: permr_tmp(:),invpr_tmp(:)

  if( ( kfl_schur == 2 .or. kfl_schur == 3 ) .and. INOTMASTER ) then

     !----------------------------------------------------------------------
     !
     ! LCORN: list of corner nodes
     !
     !----------------------------------------------------------------------

     allocate( lcorn(npoin_bb) , stat = istat )
     do kpoin = 1,npoin_bb
        lcorn(kpoin) = 0
     end do
     do ii = 1,nneig
        do jj = commd % bound_size(ii),commd % bound_size(ii+1)-1
           ipoin = commd % bound_perm(jj) ! Local numbering
           kpoin = ipoin - npoi1          ! Boundary numbering
           ifoun = 0
           kkneig: do kk = 1,nneig
              if( kk /= ii ) then
                 call par_insubd(kk,ipoin,ifoun)
                 if( ifoun == 1 ) exit kkneig 
              end if
           end do kkneig
           if( ifoun == 0 ) lcorn(kpoin) = ii
        end do
     end do

     !----------------------------------------------------------------------
     !
     ! Preconditioner graph: R_DOM_PREC, C_DOM_PREC
     !
     !----------------------------------------------------------------------

     allocate( type_prec(npoin_bb) , stat = istat )
     do ipoin = 1,npoin_bb
        kk =  r_dom_abb(ipoin+1)-r_dom_abb(ipoin) 
        allocate( type_prec(ipoin) % l( kk ) , stat = istat )
        ll = 0
        do jj = r_dom_abb(ipoin),r_dom_abb(ipoin+1)-1
           ll = ll + 1
           type_prec(ipoin) % l(ll) = c_dom_abb(jj)
        end do
     end do
     !
     ! Add nodes to graph 
     ! 
     if( kfl_schur == 3 ) then
        do kpoin = 1,npoin_bb
           do pzdom = 1,r_dom_abb(kpoin),r_dom_abb(kpoin+1)-1
              ipoin = c_dom_abb(pzdom)
              if( lcorn(ipoin) == lcorn(kpoin) .and. lcorn(ipoin) /= 0 ) then
                 do jzdom = r_dom_abi(kpoin),r_dom_abi(kpoin+1)-1      
                    qpoin = c_dom_abi(jzdom)  

                    ifoun = 0
                    kzdom = r_dom_aib(qpoin)
                    kzdom2: do while( kzdom <= r_dom_aib(qpoin+1)-1 )
                       if( c_dom_aib(kzdom) == ipoin ) then
                          ifoun = 1
                          exit kzdom2
                       end if
                       kzdom = kzdom + 1
                    end do kzdom2
                    !
                    ! Add node IPOIN to graph of KPOIN
                    !
                    if( ifoun == 0 ) then
                       isize = size( type_prec(kpoin) % l )
                       ll    = 1
                       ifoun = 0
                       popo: do while( ll <= isize )
                          if( type_prec(kpoin) % l(ll) == ipoin ) then
                             ifoun = 1
                             exit popo
                          end if
                          ll = ll + 1
                       end do popo
                       if( ifoun == 0 ) then
                          call memgen(1_ip,isize,0_ip)
                          do ll = 1,isize
                             gisca(ll) = type_prec(kpoin) % l(ll)
                          end do
                          deallocate( type_prec(kpoin) % l , stat = istat )
                          isize = isize + 1
                          allocate( type_prec(kpoin) % l(isize) , stat = istat )
                          do ll = 1,isize-1
                             type_prec(kpoin) % l(ll) = gisca(ll)
                          end do
                          type_prec(kpoin) % l(isize) = ipoin
                          call memgen(3_ip,isize,0_ip)
                       end if
                       !
                       ! Add node IPOIN to graph of KPOIN
                       !
                       isize = size( type_prec(ipoin) % l )
                       ll    = 1
                       ifoun = 0
                       popo2: do while( ll <= isize )
                          if( type_prec(ipoin) % l(ll) == kpoin ) then
                             ifoun = 1
                             exit popo2
                          end if
                          ll = ll + 1
                       end do popo2
                       if( ifoun == 0 ) then
                          call memgen(1_ip,isize,0_ip)
                          do ll = 1,isize
                             gisca(ll) = type_prec(ipoin) % l(ll)
                          end do
                          deallocate( type_prec(ipoin) % l , stat = istat )
                          isize = isize + 1
                          allocate( type_prec(ipoin) % l(isize) , stat = istat )
                          do ll = 1,isize-1
                             type_prec(ipoin) % l(ll) = gisca(ll)
                          end do
                          type_prec(ipoin) % l(isize) = kpoin
                          call memgen(3_ip,isize,0_ip)
                       end if
                    end if
                 end do
              end if
           end do
        end do
     end if

     allocate( gisnd(2,r_dom_abb(npoin_bb+1)) , stat = istat )
     allocate( gircv(2,r_dom_abb(npoin_bb+1)) , stat = istat )

     do ii = 1,nneig

        ks    = 0
        dom_i = commd % neights(ii)
        do jj = commd % bound_size(ii),commd % bound_size(ii+1)-1
           !
           ! IPOIN-JPOIN: local numbering
           !     <=>
           ! KPOIN-QPOIN: boundary numbering
           !
           ipoin = commd % bound_perm(jj)
           kpoin = ipoin - npoi1
           do izdom = r_dom_abb(kpoin),r_dom_abb(kpoin+1)-1
              qpoin = c_dom_abb(izdom)
              jpoin = qpoin + npoi1
              if( kpoin >= qpoin ) then
                 call par_insubd(ii,jpoin,ifoun)
                 if( ifoun == 1 ) then
                    ks          = ks + 1
                    gisnd(1,ks) = lninv_loc(ipoin) ! Global numbering
                    gisnd(2,ks) = lninv_loc(jpoin) ! Global numbering 
                 end if
              end if
           end do

        end do
        !
        ! Send/Recv size of graphs
        !
        bsizs4 = int(1_ip,4)
        bsizr4 = int(1_ip,4)
#ifdef MPI_OFF
#else
        call MPI_Sendrecv( ks, bsizs4,         &
             PAR_INTEGER,  dom_i, 0,         &
             kr,  bsizr4,                      &
             PAR_INTEGER,  dom_i, 0,         &
             PAR_COMM_MY_CODE, status, istat     )
#endif
        bsizs4 = int(2*ks,4)
        bsizr4 = int(2*kr,4)
#ifdef MPI_OFF
#else
        call MPI_Sendrecv( gisnd, bsizs4,      &
             PAR_INTEGER,  dom_i, 0,         &
             gircv,  bsizr4,                   &
             PAR_INTEGER,  dom_i, 0,         &
             PAR_COMM_MY_CODE, status, istat     )
#endif
        !
        ! Modify my graph
        !
        do jj = 1,kr
           kpoin = gircv(1,jj)
           qpoin = gircv(2,jj)
           ipoin = 1
           do while( lninv_loc(ipoin) /= kpoin )
              ipoin = ipoin + 1
           end do
           jpoin = 1
           do while( lninv_loc(jpoin) /= qpoin )
              jpoin = jpoin + 1
           end do

           if( ipoin /= jpoin .and. ( lcorn(ipoin-npoi1) == ii .and. lcorn(jpoin-npoi1) == ii ) ) then
              !
              ! Check if I already have JPOIN in IPOIN line
              !
              ll    = 1
              ifoun = 0
              isize = size( type_prec(ipoin-npoi1) % l )
              ipoin1: do while( ll <= isize )
                 spoin = type_prec(ipoin-npoi1) % l(ll)
                 if( spoin + npoi1 == jpoin ) then
                    ifoun = 1
                    exit ipoin1
                 end if
                 ll = ll + 1
              end do ipoin1

              if( ifoun == 0 ) then

                 isize = size( type_prec(ipoin-npoi1) % l )
                 call memgen(1_ip,isize,0_ip)
                 do ll = 1,isize
                    gisca(ll) = type_prec(ipoin-npoi1) % l(ll)
                 end do
                 deallocate( type_prec(ipoin-npoi1) % l , stat = istat )
                 isize = isize + 1
                 allocate( type_prec(ipoin-npoi1) % l(isize) , stat = istat )
                 do ll = 1,isize-1
                    type_prec(ipoin-npoi1) % l(ll) = gisca(ll)
                 end do
                 type_prec(ipoin-npoi1) % l(isize) = jpoin-npoi1
                 call memgen(3_ip,isize,0_ip)

                 isize = size( type_prec(jpoin-npoi1) % l )
                 call memgen(1_ip,isize,0_ip)
                 do ll = 1,isize
                    gisca(ll) = type_prec(jpoin-npoi1) % l(ll)
                 end do
                 deallocate( type_prec(jpoin-npoi1) % l , stat = istat )
                 isize = isize + 1
                 allocate( type_prec(jpoin-npoi1) % l(isize) , stat = istat )
                 do ll = 1,isize-1
                    type_prec(jpoin-npoi1) % l(ll) = gisca(ll)
                 end do
                 type_prec(jpoin-npoi1) % l(isize) = ipoin-npoi1
                 call memgen(3_ip,isize,0_ip)

              end if

           end if
        end do
     end do
     !
     ! Deallocate
     !
     deallocate( gisnd , stat = istat )
     deallocate( gircv , stat = istat )
     deallocate( lcorn , stat = istat )
     !
     ! R_DOM_PREC, BZDOM_PREC
     !
     allocate( r_dom_prec(npoin_bb+1) , stat = istat )
     call memchk(zero,istat,par_memor,'R_DOM_PREC','par_schurr',r_dom_prec)     

     do ipoin = 1,npoin_bb
        r_dom_prec(ipoin) = size( type_prec(ipoin) % l )
     end do
     kbb           = r_dom_prec(1)
     r_dom_prec(1) = 1
     do kpoin = 1,npoin_bb
        ibb                 = r_dom_prec(kpoin+1)
        r_dom_prec(kpoin+1) = r_dom_prec(kpoin) + kbb
        kbb                 = ibb
     end do
     nzdom_prec = r_dom_prec(npoin_bb+1) - 1

     !-------------------------------------------------------------------
     !
     ! Copy graph: C_DOM_PREC
     !
     !-------------------------------------------------------------------

     allocate( c_dom_prec(nzdom_prec) , stat = istat )
     call memchk(zero,istat,par_memor,'C_DOM_PREC','par_schurr',c_dom_prec)     

     do kpoin = 1,npoin_bb
        kk = 1
        do jj = r_dom_prec(kpoin),r_dom_prec(kpoin+1)-1
           c_dom_prec(jj) = type_prec(kpoin) % l(kk)
           kk             = kk + 1
        end do
     end do

     !-------------------------------------------------------------------
     !
     ! Order (and) renumber graph if there is a sparse direct solver
     ! 
     !-------------------------------------------------------------------

     dummi = 0
     do imodu = 1,mmodu
        solve_sol => momod(imodu) % solve
        nsolv     =  size(solve_sol)        
        do isolv = 1,nsolv
           if( solve_sol(isolv) % kfl_schur == 1 ) then
              if( solve_sol(isolv) % kfl_scpre == 1 ) then
                 dummi = 1
              end if
           end if
        end do
     end do

     if( dummi == 1 ) then
        !
        ! Renumber graph R_DOM_PREC, C_DOM_PREC for sparse direct solver
        ! and compute permutation arrays PERMR_PREC, INVPR_PREC
        ! Uses internally METIS_NodeND
        ! 
        call graphs_permut(&
             npoin_bb,nzdom_prec,r_dom_prec,c_dom_prec,permr_prec,invpr_prec,memor=memor_dom)  
        !
        ! 1. Copy second permutation using postordering: PERMR_TMP, INVPR_TMP
        !    R_DOM_PREC AND C_DOM_PREC were already been reordered
        ! 2. Renumber the graph R_DOM_PREC, C_DOM_PRED
        !
        call graphs_iniper(&
             npoin_bb,permr_tmp,invpr_tmp,memor=memor_dom)
        call graphs_postorder(&
             npoin_bb,r_dom_prec,c_dom_prec,permr_tmp,invpr_tmp,memor=memor_dom)        
        call graphs_rengra(&
             npoin_bb,nzdom_prec,permr_tmp,invpr_tmp,r_dom_prec,c_dom_prec,memor=memor_dom)
        !
        ! 1. Compose permutation PERMR_PREC, INVPR_PREC and PERMR_TMP, INVPR_TMP
        ! 2. Deallocate memory
        !
        call graphs_comper(&
             npoin_bb,permr_prec,invpr_prec,permr_tmp,invpr_tmp,memor=memor_dom)
        call graphs_deaper(&
             permr_tmp,invpr_tmp,memor=memor_dom)
     else
        !
        ! 1. Create PERMR_PREC = Id and INVPR_PREC = Id as they are used anyway
        ! 2. Order graph
        !
        call graphs_iniper(npoin_bb,permr_prec,invpr_prec,memor=memor_dom)
        do kpoin = 1,npoin_bb
           kk = r_dom_prec(kpoin+1)-r_dom_prec(kpoin)
           call heapsorti1(2_ip,kk,c_dom_prec(r_dom_prec(kpoin)))
        end do

     end if

     !-------------------------------------------------------------------
     !
     ! Deallocate temporary graph TYPE_PREC
     !
     !-------------------------------------------------------------------

     do ipoin = 1,npoin_bb
        deallocate( type_prec(ipoin) % l , stat = istat )
     end do
     deallocate( type_prec , stat = istat )

  end if

end subroutine par_schurr
