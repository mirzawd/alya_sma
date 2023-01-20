!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine solpre()
  !------------------------------------------------------------------------
  !****f* master/solpre
  ! NAME 
  !    solpre
  ! DESCRIPTION
  !    This routine constructs the IB graphs for Schur complement 
  !    type solvers and AII type preconditioners
  ! USES
  ! USED BY 
  !***
  !------------------------------------------------------------------------
  use def_kintyp
  use def_master
  use def_domain
  use def_solver
  use mod_memchk
  use mod_graphs
  use mod_memory, only : memory_alloca
  use mod_parall, only : par_memor
  implicit none 
  integer(ip)             :: ipoin,ii,ib,kpoin,izdom,jpoin,kbb,kii,kib,kbi
  integer(ip)             :: ibb,iii,iib,ibi,kk
  integer(4)              :: istat
  logical(lg)             :: onlyAii

  if( kfl_schur >= 1 .or. kfl_aiipr >= 1 ) then 

     nzdom_aii = 0
     nzdom_aib = 0
     nzdom_abi = 0
     nzdom_abb = 0
     !
     ! Check if only Aii graph is needed
     !
     if( kfl_aiipr >= 1 .and. kfl_schur == 0 ) then
        onlyAii = .true.
     else
        onlyAii = .false.
     end if
     
     if( IMASTER ) then

        npoin_ii =  0
        npoin_bb =  0

        allocate( r_dom_aii(1) , stat = istat )
        allocate( r_dom_aib(1) , stat = istat )
        allocate( r_dom_abi(1) , stat = istat )
        allocate( r_dom_abb(1) , stat = istat )
        allocate( c_dom_aii(1) , stat = istat )
        allocate( c_dom_aib(1) , stat = istat ) 
        allocate( c_dom_abi(1) , stat = istat )
        allocate( c_dom_abb(1) , stat = istat )

     else

        npoin_ii  = npoi1
        npoin_bb  = npoin-npoi1

        !----------------------------------------------------------------------
        !
        ! Memory for Graphs and matrices
        !
        !----------------------------------------------------------------------

        if( onlyAii ) then

           allocate( r_dom_aii(npoin_ii+1) , stat = istat )
           call memchk(zero,istat,par_memor,'R_DOM_AII','par_schurr',r_dom_aii)      
 
           do ipoin = 1,npoi1
              ii = 0
              do izdom = r_dom(ipoin),r_dom(ipoin+1)-1
                 jpoin = c_dom(izdom)
                 if( jpoin <= npoi1 ) then
                    nzdom_aii = nzdom_aii + 1
                    ii        = ii + 1
                 end if
              end do
              r_dom_aii(ipoin) = ii
           end do

           allocate( c_dom_aii(nzdom_aii) , stat = istat )
           call memchk(zero,istat,par_memor,'C_DOM_AII','par_schurr',c_dom_aii)      

           !----------------------------------------------------------------------
           !
           ! Fill in graphs: C_DOM_A**
           !
           !----------------------------------------------------------------------

           nzdom_aii = 0

           do ipoin = 1,npoi1
              do izdom = r_dom(ipoin),r_dom(ipoin+1)-1
                 jpoin = c_dom(izdom)
                 if( jpoin <= npoi1 ) then
                    nzdom_aii            = nzdom_aii + 1
                    c_dom_aii(nzdom_aii) = jpoin
                 end if
              end do
           end do

           !----------------------------------------------------------------------
           !
           ! Reconstruct R_DOM_A**
           !
           !----------------------------------------------------------------------

           kii          = r_dom_aii(1)
           r_dom_aii(1) = 1
           do ipoin = 1,npoi1
              iii                = r_dom_aii(ipoin+1)
              r_dom_aii(ipoin+1) = r_dom_aii(ipoin) + kii
              kii                = iii
           end do

        else

           allocate( r_dom_aii(npoin_ii+1) , stat = istat )
           call memchk(zero,istat,par_memor,'R_DOM_AII','par_schurr',r_dom_aii)     
           allocate( r_dom_aib(npoin_ii+1) , stat = istat )
           call memchk(zero,istat,par_memor,'R_DOM_AIB','par_schurr',r_dom_aib)     
           allocate( r_dom_abi(npoin_bb+1) , stat = istat )
           call memchk(zero,istat,par_memor,'R_DOM_ABI','par_schurr',r_dom_abi)     
           allocate( r_dom_abb(npoin_bb+1) , stat = istat )
           call memchk(zero,istat,par_memor,'R_DOM_ABB','par_schurr',r_dom_abb)     

           do ipoin = 1,npoi1
              ii = 0
              ib = 0
              do izdom = r_dom(ipoin),r_dom(ipoin+1)-1
                 jpoin = c_dom(izdom)
                 if( jpoin <= npoi1 ) then
                    nzdom_aii = nzdom_aii + 1
                    ii        = ii + 1
                 else
                    nzdom_aib = nzdom_aib + 1           
                    ib        = ib + 1
                 end if
              end do
              r_dom_aii(ipoin) = ii
              r_dom_aib(ipoin) = ib
           end do

           kpoin = 0
           do ipoin = npoi1+1,npoin
              kpoin = kpoin + 1
              ib = 0
              ii = 0
              do izdom = r_dom(ipoin),r_dom(ipoin+1)-1
                 jpoin = c_dom(izdom)
                 if( jpoin <= npoi1 ) then
                    nzdom_abi = nzdom_abi + 1
                    ib        = ib + 1
                 else
                    nzdom_abb = nzdom_abb + 1           
                    ii        = ii + 1
                 end if
              end do
              r_dom_abi(kpoin) = ib
              r_dom_abb(kpoin) = ii
           end do

           allocate( c_dom_aii(nzdom_aii) , stat = istat )
           call memchk(zero,istat,par_memor,'C_DOM_AII','par_schurr',c_dom_aii)     
           allocate( c_dom_aib(nzdom_aib) , stat = istat )
           call memchk(zero,istat,par_memor,'C_DOM_AIB','par_schurr',c_dom_aib)     
           allocate( c_dom_abi(nzdom_abi) , stat = istat )
           call memchk(zero,istat,par_memor,'C_DOM_ABI','par_schurr',c_dom_abi)     
           allocate( c_dom_abb(nzdom_abb) , stat = istat )
           call memchk(zero,istat,par_memor,'C_DOM_ABB','par_schurr',c_dom_abb)     

           !----------------------------------------------------------------------
           !
           ! Fill in graphs: C_DOM_A**
           !
           !----------------------------------------------------------------------

           nzdom_aii = 0
           nzdom_aib = 0
           nzdom_abi = 0
           nzdom_abb = 0

           do ipoin = 1,npoi1
              do izdom = r_dom(ipoin),r_dom(ipoin+1)-1
                 jpoin = c_dom(izdom)
                 if( jpoin <= npoi1 ) then
                    nzdom_aii            = nzdom_aii + 1
                    c_dom_aii(nzdom_aii) = jpoin
                 else
                    nzdom_aib            = nzdom_aib + 1           
                    c_dom_aib(nzdom_aib) = jpoin - npoi1
                 end if
              end do
           end do

           kpoin = 0
           do ipoin = npoi1+1,npoin
              kpoin = kpoin + 1
              do izdom = r_dom(ipoin),r_dom(ipoin+1)-1
                 jpoin = c_dom(izdom)
                 if( jpoin <= npoi1 ) then
                    nzdom_abi            = nzdom_abi + 1
                    c_dom_abi(nzdom_abi) = jpoin
                 else
                    nzdom_abb            = nzdom_abb + 1   
                    c_dom_abb(nzdom_abb) = jpoin - npoi1
                 end if
              end do
           end do

           !----------------------------------------------------------------------
           !
           ! Reconstruct R_DOM_A**
           !
           !----------------------------------------------------------------------

           kii           = r_dom_aii(1)
           kib           = r_dom_aib(1)
           r_dom_aii(1) = 1
           r_dom_aib(1) = 1
           do ipoin = 1,npoi1
              iii                = r_dom_aii(ipoin+1)
              iib                = r_dom_aib(ipoin+1)
              r_dom_aii(ipoin+1) = r_dom_aii(ipoin) + kii
              r_dom_aib(ipoin+1) = r_dom_aib(ipoin) + kib
              kii                = iii
              kib                = iib
           end do

           kbb           = r_dom_abb(1)
           kbi           = r_dom_abi(1)
           r_dom_abb(1) = 1
           r_dom_abi(1) = 1
           kpoin        = 0
           do ipoin = npoi1+1,npoin
              kpoin              = kpoin + 1
              ibb                = r_dom_abb(kpoin+1)
              ibi                = r_dom_abi(kpoin+1)
              r_dom_abb(kpoin+1) = r_dom_abb(kpoin) + kbb
              r_dom_abi(kpoin+1) = r_dom_abi(kpoin) + kbi
              kbb                = ibb
              kbi                = ibi
           end do

        end if

        !-------------------------------------------------------------------
        !
        ! Renumber Aii:
        !
        !             permr(jpoin)
        ! ipoin (new) <----------- jpoin (old)
        ! 
        !-------------------------------------------------------------------
        !
        ! Decomment following line to permute and change csrshu
        !
        if( .not. onlyAii ) then
           nullify(permr_aii)
           nullify(invpr_aii)
           !call memory_alloca(memor_dom,'PERMR_AII','solpre',permr_aii,npoin_ii,'IDENTITY')
           !call memory_alloca(memor_dom,'INVPR_AII','solpre',invpr_aii,npoin_ii,'IDENTITY')
           !do kpoin = 1,npoin_ii
           !   invpr_aii(kpoin) = kpoin
           !   permr_aii(kpoin) = kpoin
           !end do
          ! call graphs_permut_metis_postordering(npoin_ii,nzdom_aii,r_dom_aii,c_dom_aii,permr_aii,invpr_aii)
        end if
        
        !-------------------------------------------------------------------
        !
        ! Order graph: Aii and Abb
        ! 
        !-------------------------------------------------------------------

        do kpoin = 1,npoin_ii
           kk = r_dom_aii(kpoin+1)-r_dom_aii(kpoin)
           if( kk > 0 ) call heapsorti1(2_ip,kk,c_dom_aii(r_dom_aii(kpoin)))
        end do
        if( .not. onlyAii ) then
           do kpoin = 1,npoin_bb
              kk = r_dom_abb(kpoin+1)-r_dom_abb(kpoin)
              if( kk > 0 ) call heapsorti1(2_ip,kk,c_dom_abb(r_dom_abb(kpoin)))
           end do
        end if

     end if

  end if

  !-------------------------------------------------------------------
  !
  ! Parall service
  !
  !-------------------------------------------------------------------

  call par_schurr()

end subroutine solpre
