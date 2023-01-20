!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine par_cgschu(nbvar,bi,bb,xi,xb,aii,aib,abi,abb)
  !------------------------------------------------------------------------
  !****f* Parall/par_cgschu
  ! NAME 
  !    par_cgschu 
  ! DESCRIPTION
  !    This routine drives the library to solve a linear system of equations
  !    by iterative methods.
  ! USES
  ! USED BY 
  !***
  !------------------------------------------------------------------------
  use def_kintyp, only    :  ip,rp,i1p
  use def_master, only    :  npoi1,npoi2,npoi3
!  use def_master, only    :  kfl_paral
  use def_master, only    :  lninv_loc,parr1,intost
  use def_master, only    :  IMASTER,INOTMASTER,NPOIN_TYPE
  use def_solver, only    :  resi1,iters,solve_sol,memit,resi2
  use def_solver, only    :  SOL_DIAGONAL,SOL_APPROXIMATE_SCHUR,SOL_ABB
  use def_solver, only    :  SOL_MOD_DIAGONAL
  use def_solver, only    :  nzdom_prec
  use def_domain, only    :  c_dom_aii,r_dom_aii,c_dom_aib,r_dom_aib
  use def_domain, only    :  c_dom_abi,r_dom_abi,c_dom_abb,r_dom_abb
  use def_domain, only    :  c_dom_prec,r_dom_prec
  use def_domain, only    :  invpr_aii
  use def_domain, only    :  npoin
  use def_domain, only    :  npoin_ii,npoin_bb
  use def_domain, only    :  permr_prec,invpr_prec
  use mod_parall, only    :  PAR_REAL
  use mod_skyline, only   :  chofac,chosol
  use mod_alya_direct_solver, only : alya_Symbolical_CSR_LU
  use mod_alya_direct_solver, only : alya_Numerical_CSR_LU
  use def_parall
  use mod_memchk
  use mod_csrdir
  use mod_graphs
  use mod_parall, only    :  PAR_COMM_MY_CODE,PAR_INTEGER,commd
  use mod_communications_global, only : PAR_SUM
  use def_mpi
  implicit none


  integer(ip), intent(in)   :: nbvar
  !
  ! => MATRICES ---------------------------------------------------------
  real(rp),     intent(in)    :: aii(nbvar,nbvar,*)
  real(rp),     intent(in)    :: aib(nbvar,nbvar,*)
  real(rp),     intent(in)    :: abi(nbvar,nbvar,*)
  real(rp),     intent(in)    :: abb(nbvar,nbvar,*)
  real(rp),     intent(in)    :: bi(*)
  real(rp),     intent(in)    :: bb(*)
  real(rp),     intent(inout) :: xb(*)
  real(rp),     intent(out)   :: xi(*)
  ! <= MATRICES ---------------------------------------------------------
  !
  ! => SKYLINE FOR AII --------------------------------------------------
  integer(ip)               :: nskyl
  integer(ip),  pointer     :: iskyl(:)
  real(rp),     pointer     :: askyl(:)
  ! <= SKYLINE OR AII  --------------------------------------------------
  !
  ! => SPARSE DIRECT SOLVER FOR AII -------------------------------------
  integer(ip), pointer       :: iL(:),jL(:)
  real(rp),    pointer       :: Ln(:)
  integer(ip), pointer       :: iU(:),jU(:)
  real(rp),    pointer       :: Un(:)
  integer(ip), pointer       :: invpR(:),invpC(:)
  ! <= SPARSE DIRECT SOLVER FOR AII -------------------------------------
  !
  ! => CONJUGATE GRADIENT  ----------------------------------------------
  integer(ip)               :: nrows,ierro,maxiter
  integer(ip)               :: npoii2,npoii3,npoii1,npoiin
  integer(ip)               :: idprecon,kfl_cvgso,kk
  real(rp)                  :: alpha,beta,rho,stopcri,resid,eps
  real(rp)                  :: invnb,newrho,dummr,bnorm,resin
  real(rp)                  :: invb2,rnorm,xnorm,resfi
  real(rp),     pointer     :: rr(:),pp(:),zz(:),invdiag(:)
  real(rp),     pointer     :: ww(:),zi(:),sb(:),zirhs(:)
  ! <= CONJUGATE GRADIENT -----------------------------------------------
  !
  ! => SKYLINE FOR PRECONDITIONER P -------------------------------------
  integer(4)                :: dom_i
  integer(4)                :: bsizs4,bsizr4
  integer(ip),  pointer     :: gisnd(:,:) 
  integer(ip),  pointer     :: gircv(:,:) 
  real(rp),     pointer     :: gesnd(:,:) 
  real(rp),     pointer     :: gercv(:,:) 
  real(rp),     pointer     :: prec(:,:,:)
  real(rp),     pointer     :: pre2(:,:,:)
  integer(ip)               :: ifoun
  integer(ip)               :: jzdom,kzdom,lzdom,pzdom,ks,kr
  integer(ip)               :: nskyl_p
  integer(ip),  pointer     :: iskyl_p(:)
  real(rp),     pointer     :: askyl_p(:)
  real(rp)                  :: xdiag
  ! <= SKYLINE FOR PRECONDITIONER P -------------------------------------
  !
  ! => SPARSE DIRECT SOLVER FOR PRECONDITIONER P ------------------------
  integer(ip), pointer      :: iLpre(:),jLpre(:)
  real(rp),    pointer      :: Lnpre(:)
  integer(ip), pointer      :: iUpre(:),jUpre(:)
  real(rp),    pointer      :: Unpre(:)
  integer(ip), pointer      :: lcorn(:)
  ! <= SPARSE DIRECT SOLVER FOR PRECONDITIONER P ------------------------
  !
  integer(ip)               :: ipoin,izdom,jpoin,kpoin,lpoin
  integer(ip)               :: ii,jj,ll,kskyl,info,ipoin_bb,jpoin_bb
  integer(ip)               :: ipoin_new,kpoin_new,jpoin_new
  integer(ip)               :: ipoin_bb_new,jpoin_bb_new
  integer(4)                :: istat
  real(rp)                  :: time1,time2,time3,time4,time5
  real(rp)                  :: time6,time7,time8,time9,time10
  real(rp)                  :: ctim1,ctim2,ctim3
  logical(lg)               :: Checkmode    
  !
  ! Variables
  !
  maxiter   = solve_sol(1) % miter
  eps       = solve_sol(1) % solco
  idprecon  = solve_sol(1) % kfl_preco
  kfl_cvgso = solve_sol(1) % kfl_cvgso
  Checkmode = .false.

  ctim1     = 0.0_rp
  ctim2     = 0.0_rp
  ctim3     = 0.0_rp
  time1     = 0.0_rp
  time2     = 0.0_rp
  time3     = 0.0_rp
  time4     = 0.0_rp
  time5     = 0.0_rp
  time10    = 0.0_rp

  call cputim(time1)
  
  if( IMASTER ) then

     npoin_ii =  0
     npoin_bb =  0
     nrows    =  0
     npoiin   =  0
     npoii1   =  0
     npoii2   =  0
     npoii3   = -1

  else

     nrows    = nbvar * npoin_bb
     npoiin   = npoin
     npoii1   = npoi1
     npoii2   = npoi2
     npoii3   = npoi3

     if( solve_sol(1) % kfl_scaii == 0 ) then

        !----------------------------------------------------------------------
        !
        ! Skyline for Aii
        !
        !----------------------------------------------------------------------

        if( nbvar == 1 ) then
           !
           ! NDOFN = 1
           !
           allocate(iskyl(npoin_ii+1),stat=istat)
           call memchk(0_ip,istat,memit,'ISKYL','par_cgschu',iskyl)

           do ipoin = 1,npoin_ii+1
              iskyl(ipoin) = npoin_ii
           end do

           do ipoin = 1,npoin_ii
              do izdom = r_dom_aii(ipoin),r_dom_aii(ipoin+1)-1
                 jpoin = c_dom_aii(izdom)     
                 if( ipoin >= jpoin .and. jpoin < iskyl(ipoin+1) ) iskyl(ipoin+1) = jpoin
              end do
           end do
        else
           call runend('SCHUR: NOT CODED')
        end if

        nskyl    =  1
        iskyl(1) =  1       
        do ipoin = 1,npoin_ii
           kskyl = ipoin - iskyl(ipoin+1) + 1
           nskyl = nskyl + kskyl
           iskyl(ipoin+1) = nskyl
        end do
        nskyl = nskyl - 1 
        !
        ! Assemble Aii
        !
        allocate(askyl(nskyl),stat=istat)
        call memchk(0_ip,istat,memit,'ASKYL','par_cgschu',askyl)
        do ipoin = 1,npoin_ii
           do izdom = r_dom_aii(ipoin),r_dom_aii(ipoin+1)-1
              jpoin = c_dom_aii(izdom)
              if( jpoin < ipoin ) then
                 kskyl        = iskyl(ipoin+1) - 1 - (ipoin-jpoin)
                 askyl(kskyl) = askyl(kskyl) + aii(1,1,izdom)
              else if( ipoin == jpoin ) then
                 kskyl        = iskyl(ipoin+1) - 1
                 askyl(kskyl) = askyl(kskyl) + aii(1,1,izdom)  
              end if
           end do
        end do
        !
        ! Factorize Aii
        !
        call chofac(npoin_ii*nbvar,nskyl,iskyl,askyl,info) 
        if( info /= 0 ) call runend('SCHUR: ERROR WHILE DOING CHOLESKY FACTORIZATION')

     else if( solve_sol(1) % kfl_scaii == 1 ) then

        !----------------------------------------------------------------------
        !
        ! Sparse direct solver for Aii
        !
        !----------------------------------------------------------------------

        nullify(Il,Jl,Ln,iU,jU,Un)                            ! Permutation arrays
        nullify(invpR,invpC)
     
        call alya_Symbolical_CSR_LU(npoin_ii,r_dom_aii,c_dom_aii,iL,jL,iU,jU)
        call alya_Numerical_CSR_LU (npoin_ii,nbvar,r_dom_aii,c_dom_aii,aii,iL,jL,Ln,iU,jU,Un,info)

        !call CSR_LU_Factorization(&                           ! CSR LU Factorization for Aii
        !     npoin_ii,nbvar,r_dom_aii,c_dom_aii,aii,IL,JL,&
        !     LN,IU,JU,UN,info)
        if( info /= 0 ) call runend('PAR_CGSCHU: SINGULAR MATRIX')

     end if

  end if

  call cputim(time2)

  !----------------------------------------------------------------------
  !
  ! Preconditioner
  !
  !----------------------------------------------------------------------

  if( ( idprecon == SOL_APPROXIMATE_SCHUR .or. idprecon == SOL_ABB ) .and. INOTMASTER ) then
     !
     ! LCORN: list of corner nodes
     !
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
     !
     ! Allocate memory: PREC
     !
     allocate( prec(1,1,nzdom_prec)   , stat = istat )
     allocate( pre2(1,1,nzdom_prec)   , stat = istat )
     do izdom = 1,nzdom_prec
        prec(1,1,izdom) = 0.0_rp
     end do
     !
     ! P = Abb. Their graph can be different.
     !
     do kpoin = 1,npoin_bb
        kpoin_new = permr_prec(kpoin)
        do kzdom = r_dom_abb(kpoin),r_dom_abb(kpoin+1)-1
           jpoin     = c_dom_abb(kzdom)
           jpoin_new = permr_prec(jpoin)
           pzdom     = r_dom_prec(kpoin_new)
           do while( pzdom <= r_dom_prec(kpoin_new+1)-1 )
              if( c_dom_prec(pzdom) == jpoin_new ) then
                 prec(1,1,pzdom) = abb(1,1,kzdom)
                 pzdom = r_dom_prec(kpoin_new+1)
              end if
              pzdom = pzdom + 1
           end do
        end do
     end do
     !
     ! P <= Abb - Abi diag(Aii)^{-1} Aib
     !                                            ipoin
     !                                              |
     !            ipoin               lpoin   +-      -+
     !             |                    |     |     .  | / diag(.)
     !          +-       -+    +-         -+  |     .  | / diag(.)
     !          |         |    |           |  |     .  | / diag(.)
     ! kpoin => |  x      | =  | ......... |  |     .  | / diag(.)
     !          +-       -+    +-          +  |     .  | / diag(.)
     !             |                          |     .  | / diag(.)
     !           pzdom                        +-      -+
     !  
     !          -----------    -------------  ----------
     !             P               Abi            Aib     diag(Aii^-1) 
     !
     if( idprecon == SOL_APPROXIMATE_SCHUR ) then
        do kpoin = 1,npoin_bb
           kpoin_new = permr_prec(kpoin)
           do pzdom = r_dom_prec(kpoin_new),r_dom_prec(kpoin_new+1)-1
              ipoin_new = c_dom_prec(pzdom) 
              ipoin     = invpr_prec(ipoin_new) 
              if( lcorn(ipoin) == lcorn(kpoin) .and. lcorn(ipoin) /= 0 ) then
                 do jzdom = r_dom_abi(kpoin),r_dom_abi(kpoin+1)-1      
                    lpoin = c_dom_abi(jzdom)                           ! LPOIN in local

                    ifoun = 0
                    lzdom = r_dom_aii(lpoin)
                    lzdom1: do while( lzdom <= r_dom_aii(lpoin+1)-1 )
                       if( c_dom_aii(lzdom) == lpoin ) then
                          ifoun = 1
                          exit lzdom1
                       end if
                       lzdom = lzdom + 1
                    end do lzdom1

                    xdiag = 1.0_rp / Aii(1,1,lzdom)

                    ifoun = 0
                    kzdom = r_dom_aib(lpoin)
                    kzdom2: do while( kzdom <= r_dom_aib(lpoin+1)-1 )
                       if( c_dom_aib(kzdom) == ipoin ) then
                          prec(1,1,pzdom) = prec(1,1,pzdom) - Abi(1,1,jzdom) * xdiag * Aib(1,1,kzdom)
                          ifoun = 1
                          exit kzdom2
                       end if
                       kzdom = kzdom + 1
                    end do kzdom2
                 end do
              end if
           end do
        end do
     end if
     !
     ! Temporal preconditioner
     !
     do izdom = 1,nzdom_prec
        pre2(1,1,izdom) = prec(1,1,izdom)
     end do
     !
     ! Exchange matrix
     !
     allocate( gisnd(2,r_dom_prec(npoin_bb+1)) , stat = istat )
     allocate( gircv(2,r_dom_prec(npoin_bb+1)) , stat = istat )
     allocate( gesnd(1,r_dom_prec(npoin_bb+1)) , stat = istat )
     allocate( gercv(1,r_dom_prec(npoin_bb+1)) , stat = istat )

     do ii = 1,nneig

        ks    = 0
        dom_i = commd % neights(ii)
        do jj = commd % bound_size(ii),commd % bound_size(ii+1)-1
           !
           ! IPOIN-JPOIN: local numbering
           !     <=>
           ! KPOIN-LPOIN: boundary numbering
           !
           ipoin    = commd % bound_perm(jj)
           ipoin_bb = ipoin - npoi1           
           if( lcorn(ipoin_bb) == ii ) then
              ipoin_bb_new = permr_prec(ipoin_bb)
              do izdom = r_dom_prec(ipoin_bb_new),r_dom_prec(ipoin_bb_new+1)-1
                 jpoin_bb_new = c_dom_prec(izdom)
                 jpoin_bb     = invpr_prec(jpoin_bb_new) 
                 jpoin        = jpoin_bb + npoi1
                 if( ipoin_bb >= jpoin_bb .and. lcorn(jpoin_bb) == ii ) then
                    ks          = ks + 1
                    gisnd(1,ks) = lninv_loc(ipoin) ! Global numbering ahora
                    gisnd(2,ks) = lninv_loc(jpoin) ! Global numbering ahora
                    gesnd(1,ks) = pre2(1,1,izdom)  ! Value
                 end if
              end do
           end if
        end do

        bsizs4           = int(1_ip,4)
        bsizr4           = int(1_ip,4)
#ifdef MPI_OFF
#else
        call MPI_Sendrecv( ks, bsizs4,         &
             PAR_INTEGER,  dom_i, 0_4,         &
             kr,  bsizr4,                      &
             PAR_INTEGER,  dom_i, 0_4,         &
             PAR_COMM_MY_CODE, status, istat     )
#endif
        bsizs4 = int(2*ks,4)
        bsizr4 = int(2*kr,4)
#ifdef MPI_OFF
#else
        call MPI_Sendrecv( gisnd, bsizs4,      &
             PAR_INTEGER,  dom_i, 0_4,         &
             gircv,  bsizr4,                   &
             PAR_INTEGER,  dom_i, 0_4,         &
             PAR_COMM_MY_CODE, status, istat     )
#endif

        bsizs4 = int(ks,4)
        bsizr4 = int(kr,4)
#ifdef MPI_OFF
#else
        call MPI_Sendrecv( gesnd, bsizs4,      &
             PAR_REAL, dom_i, 0_4, &
             gercv, bsizr4,                    &
             PAR_REAL, dom_i, 0_4, &
             PAR_COMM_MY_CODE, status, istat     )
#endif
        !
        ! Fill in my matrix
        !
        do jj = 1,kr
           !
           ! I receive KPOIN-LPOIN in global numbering
           ! They correspond to IPOIN-JPOIN
           !
           kpoin = gircv(1,jj)
           lpoin = gircv(2,jj)
           ipoin = 1
           do while( lninv_loc(ipoin) /= kpoin )
              ipoin = ipoin + 1
           end do
           jpoin = 1
           do while( lninv_loc(jpoin) /= lpoin )
              jpoin = jpoin + 1
           end do
           ipoin_bb = ipoin - npoi1
           jpoin_bb = jpoin - npoi1
           !
           ! Then look for corresponding IZDOM
           ! KPOIN is in boundary numbering
           !
           if( lcorn(ipoin_bb) == ii .and. lcorn(jpoin_bb) == ii ) then
              !
              ! IZDOM: IPOIN-JPOIN. Only for non-corner nodes.
              !
              ipoin_bb_new = permr_prec(ipoin_bb)
              jpoin_bb_new = permr_prec(jpoin_bb)
              izdom        = r_dom_prec(ipoin_bb_new)
              izdom1: do while( izdom <= r_dom_prec(ipoin_bb_new+1)-1 )
                 if( c_dom_prec(izdom) == jpoin_bb_new ) exit izdom1
                 izdom = izdom + 1
              end do izdom1
              if( izdom /= r_dom_prec(ipoin_bb_new+1) ) then
                 prec(1,1,izdom) = prec(1,1,izdom) + gercv(1,jj)
              end if
              !
              ! IZDOM: JPOIN-IPOIN
              !
              izdom = r_dom_prec(jpoin_bb_new)
              izdom2: do while( izdom <= r_dom_prec(jpoin_bb_new+1)-1 )
                 if( c_dom_prec(izdom) == ipoin_bb_new ) exit izdom2
                 izdom = izdom + 1
              end do izdom2
              if( izdom /= r_dom_prec(jpoin_bb_new+1) .and. ipoin /= jpoin ) then
                 prec(1,1,izdom) = prec(1,1,izdom) + gercv(1,jj)
              end if
           end if

        end do

     end do
     !
     deallocate( pre2  , stat = istat )
     deallocate( gisnd , stat = istat )
     deallocate( gesnd , stat = istat )
     deallocate( gircv , stat = istat )
     deallocate( gercv , stat = istat )

  end if

  call cputim(time3)
 
  !----------------------------------------------------------------------
  !
  ! End preconditioner
  !
  !----------------------------------------------------------------------

  !----------------------------------------------------------------------
  !
  ! Conjugate gradient
  !
  !----------------------------------------------------------------------

  if( INOTMASTER ) then

     allocate(rr(nrows),stat=istat) 
     call memchk(0_ip,istat,memit,'RR','par_cgschu',rr)

     allocate(pp(nrows),stat=istat)  
     call memchk(0_ip,istat,memit,'PP','par_cgschu',pp)

     allocate(zz(nrows),stat=istat) 
     call memchk(0_ip,istat,memit,'ZZ','par_cgschu',zz)

     allocate(ww(nrows),stat=istat) 
     call memchk(0_ip,istat,memit,'WW','par_cgschu',ww)

     allocate(invdiag(nrows),stat=istat)
     call memchk(0_ip,istat,memit,'INVDIAG','par_cgschu',invdiag)

     allocate(sb(nrows),stat=istat)
     call memchk(0_ip,istat,memit,'SB','par_cgschu',sb)

     allocate(zi(npoi1*nbvar),stat=istat) 
     call memchk(0_ip,istat,memit,'ZI','par_cgschu',zi)

     allocate(zirhs(npoi1*nbvar),stat=istat) 
     call memchk(0_ip,istat,memit,'ZIRHS','par_cgschu',zirhs)

  end if
  !
  ! Schur system RHS: bb <= bb - Abi Aii^{-1} bi
  ! Original bb and bi were exchanged previously
  !
  if( INOTMASTER ) then
     do ipoin = 1,npoin_ii
        zi(ipoin) = bi(ipoin)
     end do
     if( solve_sol(1) % kfl_scaii == 0 ) then
        call chosol(&
             npoin_ii*nbvar,nskyl,iskyl,1_ip,askyl,zi,&
             npoin_ii*nbvar,info)
     else
        call CSR_LUsol(&
             npoin_ii,nbvar,invpr_aii,invpr_aii,IL,JL,& ! Solve Aii^{-1}
             LN,IU,JU,UN,bi,zi) 
     end if
     do kpoin = 1,npoin_bb
        ww(kpoin) = 0.0_rp
        do izdom = r_dom_abi(kpoin),r_dom_abi(kpoin+1)-1
           jpoin = c_dom_abi(izdom)
           ww(kpoin) = ww(kpoin) - abi(1,1,izdom) * zi(jpoin)
        end do
     end do
     parr1 => ww
     call par_slesch()    
     do kpoin = 1,npoin_bb
        sb(kpoin) = bb(kpoin) + ww(kpoin)
     end do
  end if
  !
  ! Initial solution
  !
  ierro = 0
  resi1 = 1.0_rp
  resi2 = 1.0_rp
  resin = 1.0_rp
  solve_sol(1) % iters = 0
  call outcso(dummr,dummr,dummr)
  !
  ! INVDIAG
  !
  ipoin = npoi1
  do kpoin = 1,npoin_bb
     ipoin = ipoin + 1
     jj    =  r_dom_abb(kpoin)
     ll    = -1
     do while( jj < r_dom_abb(kpoin+1) .and. ll == -1 )
        if( c_dom_abb(jj) == kpoin ) then
           ll = jj
        end if
        jj = jj + 1
     end do
     if( ll /= -1 ) then
        invdiag(kpoin) = abb(1,1,ll)
     end if

  end do
  !
  ! P <= diag(Abb) - diag( - Abi diag(Aii^{-1}) Aib )
  !
  if( idprecon == SOL_APPROXIMATE_SCHUR .or. idprecon == SOL_MOD_DIAGONAL ) then
     do kpoin = 1,npoin_bb
        do izdom = r_dom_abb(kpoin),r_dom_abb(kpoin+1)-1
           if( c_dom_abb(izdom) == kpoin ) then
              do jzdom = r_dom_abi(kpoin),r_dom_abi(kpoin+1)-1
                 lpoin = c_dom_abi(jzdom)                           ! LPOIN in local

                 lzdom = r_dom_aii(lpoin)
                 lzdom3: do while( lzdom <= r_dom_aii(lpoin+1)-1 )
                    if( c_dom_aii(lzdom) == lpoin ) exit lzdom3
                    lzdom = lzdom + 1
                 end do lzdom3
                 xdiag = 1.0_rp / Aii(1,1,lzdom)

                 kzdom = r_dom_aib(lpoin)
                 kzdom4: do while( kzdom <= r_dom_aib(lpoin+1)-1 )
                    if( c_dom_aib(kzdom) == kpoin ) then
                       invdiag(kpoin) = invdiag(kpoin) - Abi(1,1,jzdom) * xdiag * Aib(1,1,kzdom)
                       exit kzdom4
                    end if
                    kzdom = kzdom + 1
                 end do kzdom4
              end do
           end if
        end do
     end do
  end if
  !
  ! Parallel exchange
  !
  if( INOTMASTER ) then
     parr1 => invdiag
     call par_slesch()    
  end if
  do kpoin = 1,npoin_bb
     invdiag(kpoin) = 1.0_rp / invdiag(kpoin)
  end do

  if( ( idprecon == SOL_APPROXIMATE_SCHUR .or. idprecon == SOL_ABB ) .and. INOTMASTER ) then

     !-------------------------------------------------------------------
     !
     ! Put diagonal on corner nodes
     !
     !-------------------------------------------------------------------

     do kpoin = 1,npoin_bb
        !
        ! This node is repeated: this is a corner node
        !
        kpoin_new = permr_prec(kpoin)
        if( lcorn(kpoin) == 0 ) then
           do jpoin = 1,npoin_bb
              jpoin_new = permr_prec(jpoin)
              do izdom = r_dom_prec(jpoin_new),r_dom_prec(jpoin_new+1)-1
                 if( c_dom_prec(izdom) == kpoin_new ) then
                    prec(1,1,izdom) = 0.0_rp
                 end if
              end do
           end do
           do izdom = r_dom_prec(kpoin_new),r_dom_prec(kpoin_new+1)-1
              jpoin_new = c_dom_prec(izdom)
              jpoin     = invpr_prec(jpoin_new)
              if( jpoin == kpoin ) then
                 prec(1,1,izdom) = 1.0_rp / invdiag(kpoin)
              else
                 prec(1,1,izdom) = 0.0_rp
              end if
           end do
        else
           do izdom = r_dom_prec(kpoin_new),r_dom_prec(kpoin_new+1)-1
              jpoin_new = c_dom_prec(izdom)
              jpoin     = invpr_prec(jpoin_new)
              if( lcorn(jpoin) /= lcorn(kpoin) ) then
                 prec(1,1,izdom) = 0.0_rp
              end if
           end do
        end if
     end do

     deallocate( lcorn , stat = istat )

     if( solve_sol(1) % kfl_scpre == 0 ) then

        !-------------------------------------------------------------------
        !
        ! Skyline for P
        !
        !-------------------------------------------------------------------

        allocate(iskyl_p(npoin_bb+1),stat=istat)
        call memchk(0_ip,istat,memit,'ISKYL_P','par_cgschu',iskyl_p)

        do ipoin = 1,npoin_bb+1
           iskyl_p(ipoin) = npoin_bb
        end do
        do ipoin = 1,npoin_bb
           do izdom = r_dom_prec(ipoin),r_dom_prec(ipoin+1)-1
              jpoin = c_dom_prec(izdom)     
              if( ipoin >= jpoin .and. jpoin < iskyl_p(ipoin+1) ) iskyl_p(ipoin+1) = jpoin
           end do
        end do
        nskyl_p    =  1
        iskyl_p(1) =  1
        do ipoin = 1,npoin_bb
           kskyl            = ipoin   - iskyl_p(ipoin+1) + 1
           nskyl_p          = nskyl_p + kskyl
           iskyl_p(ipoin+1) = nskyl_p
        end do
        nskyl_p = nskyl_p - 1 
        !
        ! Assemble P
        !
        allocate(askyl_p(nskyl_p),stat=istat)
        call memchk(0_ip,istat,memit,'ASKYL_P','par_cgschu',askyl_p)

        do ipoin = 1,npoin_bb
           do izdom = r_dom_prec(ipoin),r_dom_prec(ipoin+1)-1
              jpoin = c_dom_prec(izdom)
              if( jpoin < ipoin ) then
                 kskyl          = iskyl_p(ipoin+1) - 1 - (ipoin-jpoin)
                 askyl_p(kskyl) = askyl_p(kskyl) + prec(1,1,izdom)
              else if( ipoin == jpoin ) then
                 kskyl          = iskyl_p(ipoin+1) - 1
                 askyl_p(kskyl) = askyl_p(kskyl) + prec(1,1,izdom)  
              end if
           end do
        end do
        !
        ! Factorize P
        !
        call chofac(npoin_bb*nbvar,nskyl_p,iskyl_p,askyl_p,info) 
        if( info /= 0 ) call runend('SCHUR: ERROR WHILE DOING CHOLESKY FACTORIZATION')

     else if( solve_sol(1) % kfl_scpre == 1 ) then

        !-------------------------------------------------------------------
        !
        ! Sparse for P
        !
        !-------------------------------------------------------------------

        !if( kfl_paral== 1 ) then
        !   open(unit=100+kfl_paral,file='matrice-'//trim(intost(kfl_paral))//'.ps',status='unknown')
        !   call pspltm(&
        !        npoin_bb,npoin_bb,1_ip,0_ip,c_dom_prec,r_dom_prec,prec,&
        !        'matrice',0_ip,18.0_rp,'cm',&
        !        0_ip,0_ip,2_ip,100+kfl_paral)
        !   close(unit=100+kfl_paral)
        !end if

        nullify(Ilpre,Jlpre,Lnpre,iUpre,jUpre,Unpre)                  ! Permutation arrays
        call CSR_LU_Factorization(&                                   ! CSR LU Factorization  
             npoin_bb,nbvar,r_dom_prec,c_dom_prec,prec,ILpre,JLpre,&
             LNpre,IUpre,JUpre,UNpre,info)
        if( info /= 0 ) call runend('PAR_CGSCHU: SINGULAR MATRIX')

        !if( kfl_paral== 1 ) then
        !   open(unit=200+kfl_paral,file='matriceL-'//trim(intost(kfl_paral))//'.ps',status='unknown')
        !   call pspltm(&
        !        npoin_bb,npoin_bb,1_ip,0_ip,jLpre,iLpre,prec,&
        !        'matrice',0_ip,18.0_rp,'cm',&
        !        0_ip,0_ip,2_ip,200+kfl_paral)
        !   close(unit=200+kfl_paral)
        !end if

     end if

  end if
  !
  ! Check symmetry
  !
  if( Checkmode .and. (idprecon == SOL_APPROXIMATE_SCHUR .or. idprecon == SOL_ABB) ) then
     if( INOTMASTER ) then 
        do kpoin = 1,npoin_bb
           do izdom = r_dom_prec(kpoin),r_dom_prec(kpoin+1)-1
              jpoin = c_dom_prec(izdom)
              ii    = r_dom_prec(jpoin)
              ipoin = 0
              caca2: do while( ii <= r_dom_prec(jpoin+1)-1)
                 if( c_dom_prec(ii) == kpoin ) then
                    ipoin = kpoin
                    exit caca2
                 end if
                 ii = ii + 1
              end do caca2
              if( ipoin == 0 ) print*,'PROBLEMMMMMMMMMMMMMMMMMMMMMMMMMM-par_cgschu'
              if( abs(prec(1,1,ii)-prec(1,1,izdom)) > 1.0e-6_rp ) then
                 print*,'NON-SYMMETRIC: IPOIN, JPOIN=',lninv_loc(ipoin+npoi1),lninv_loc(jpoin+npoi1)
                 print*,'NON-SYMMETRIC: VALUES=      ',prec(1,1,ii),prec(1,1,izdom)
              end if
           end do
        end do
     end if
  end if
  !
  ! w = L^-1 b, BNORM = || L^-1 w ||
  !
  call cputim(time4)

  if( idprecon == SOL_DIAGONAL .or. idprecon == SOL_MOD_DIAGONAL ) then

     do kpoin = 1,npoin_bb
        ww(kpoin) = invdiag(kpoin) * sb(kpoin) 
     end do

  else if( idprecon == SOL_APPROXIMATE_SCHUR .or. idprecon == SOL_ABB ) then

     if( INOTMASTER ) then
        do kpoin = 1,npoin_bb
           ww(kpoin) = sb(kpoin) 
        end do
        if( solve_sol(1) % kfl_scpre == 0 ) then
           call chosol(&
                npoin_bb*nbvar,nskyl_p,iskyl_p,1_ip,askyl_p,ww,&
                npoin_bb*nbvar,info)
        else
           call CSR_LUsol(&
                npoin_bb,nbvar,invpr_prec,invpr_prec,ILpre,JLpre,&
                LNpre,IUpre,JUpre,UNpre,sb,ww)
        end if
     end if

  end if

  bnorm = 0.0_rp
  do ipoin = npoii2,npoii3
     kpoin = ipoin - npoi1
     bnorm = bnorm + ww(kpoin) * ww(kpoin)
  end do
  call PAR_SUM(bnorm)
  bnorm = sqrt(bnorm)
  !
  ! || L^-1 b|| = 0: trivial solution x = 0
  !    
  if( bnorm <= 1.0e-12_rp ) then
     ipoin = 0
     do kpoin = 1,nrows
        xb(kpoin) = 0.0_rp
     end do
     resid                =  0.0_rp
     resin                =  0.0_rp
     resi1                =  0.0_rp
     invnb                =  0.0_rp
     solve_sol(1) % resin =  0.0_rp 
     solve_sol(1) % resi2 =  0.0_rp
     ierro                = -1
     goto 10
  end if
  invnb = 1.0_rp / bnorm

  !----------------------------------------------------------------------
  !
  ! Initial residual: r = b - Ax
  !
  !----------------------------------------------------------------------
  !
  ! XNORM = ||x||
  !
  xnorm = 0.0_rp
  do ipoin = npoii2,npoii3
     kpoin = ipoin - npoi1
     xnorm = xnorm + abs(xb(kpoin))
  end do
  call PAR_SUM(xnorm)

  if( xnorm == 0.0_rp ) then
     !
     ! Initial x is zero: r = b
     !
     do ii= 1,nrows
        rr(ii) = sb(ii)
     end do

  else 
     !
     ! r = b - A x
     !
     if( INOTMASTER ) then
        if( solve_sol(1) % kfl_scaii == 0 ) then
           call par_schuax(&
                nbvar,zi,xb,aib,abi,abb,r_dom_aib,c_dom_aib,r_dom_abb,&
                c_dom_abb,r_dom_abi,c_dom_abi,askyl,iskyl,nskyl,rr)
        else
           call par_schua2(&
                1_ip,nbvar,zi,xb,aib,abi,abb,r_dom_aib,c_dom_aib,r_dom_abb,&
                c_dom_abb,r_dom_abi,c_dom_abi,rr)        
           do ii = 1,npoi1
              zirhs(ii) = zi(ii)
           end do
           call CSR_LUsol(&
                npoi1,nbvar,invpr_aii,invpr_aii,IL,JL,&  ! Solve Aii^{-1}
                LN,IU,JU,UN,zirhs,zi)
           call par_schua2(&
                2_ip,nbvar,zi,xb,aib,abi,abb,r_dom_aib,c_dom_aib,r_dom_abb,&
                c_dom_abb,r_dom_abi,c_dom_abi,rr)        
        end if
     end if
     do kpoin = 1,npoin_bb
        rr(kpoin) = sb(kpoin) - rr(kpoin) 
     end do
  end if

  rnorm = 0.0_rp
  invb2 = 0.0_rp
  do ipoin = npoii2,npoii3
     kpoin = ipoin - npoi1
     rnorm = rnorm + rr(kpoin) * rr(kpoin)
     invb2 = invb2 + sb(kpoin) * sb(kpoin)
  end do
  call PAR_SUM(rnorm)
  call PAR_SUM(invb2)

  rnorm = sqrt(rnorm)
  invb2 = sqrt(invb2)

  if( invb2 /= 0.0_rp ) invb2 = 1.0_rp / invb2 
  solve_sol(1) % resi2 = rnorm * invb2
  !
  ! Preconditioned residual: z = L^{-1} r
  !
  if( idprecon == SOL_DIAGONAL .or. idprecon == SOL_MOD_DIAGONAL ) then

     do ipoin = 1,nrows
        zz(ipoin) = rr(ipoin) * invdiag(ipoin)
     end do

  else if( idprecon == SOL_APPROXIMATE_SCHUR .or. idprecon == SOL_ABB ) then

     if( INOTMASTER ) then
        do ipoin = 1,nrows
           zz(ipoin) = rr(ipoin) 
        end do
        if( solve_sol(1) % kfl_scpre == 0 ) then
           call chosol(&
                npoin_bb*nbvar,nskyl_p,iskyl_p,1_ip,askyl_p,zz,&
                npoin_bb*nbvar,info)
        else       
           call CSR_LUsol(&
                npoin_bb,nbvar,invpr_prec,invpr_prec,ILpre,JLpre,&
                LNpre,IUpre,JUpre,UNpre,rr,zz)           
        end if
     end if

  end if
  !
  ! Initial rho = < L^-1 r , r > 
  !
  newrho = 0.0_rp
  do ipoin = npoii2,npoii3
     kpoin = ipoin - npoi1
     newrho = newrho + rr(kpoin) * zz(kpoin)
  end do
  call PAR_SUM(newrho)
  resid = sqrt(newrho)     
  resin = resid * invnb
  resi1 = resin
  solve_sol(1)%resin = resin
  !
  ! Adaptive residual
  !
  if( solve_sol(1)%kfl_adres == 1 ) then
     eps = max( solve_sol(1) % resin*solve_sol(1) % adres , solve_sol(1) % solmi )
  end if
  !
  ! stop criterion = ||b|| * eps
  !
  stopcri = eps * bnorm
  !
  ! Test if the initial guess is the solution
  !
  if ( resid <= stopcri ) then
     ierro = -2 
     resfi =  resin
     goto 10
  end if
  !
  ! Initial direction: not needed by GMRES nor RICHARDSON
  !
  do ii = 1, nrows
     pp(ii) = zz(ii)
  end do

  call cputim(time5)

  !-----------------------------------------------------------------
  !
  !  MAIN LOOP
  !
  !-----------------------------------------------------------------

  do while( solve_sol(1) % iters < maxiter .and. resid > stopcri )
     !
     ! q^{k+1} =  A R^-1 p^{k+1}
     !
     ! A p = Abb p - Abi Aii^{-1} Aib p
     !
     call cputim(time6)
     if( INOTMASTER ) then
        if( solve_sol(1) % kfl_scaii == 0 ) then
           call par_schuax(&
                nbvar,zi,pp,aib,abi,abb,r_dom_aib,c_dom_aib,r_dom_abb,&
                c_dom_abb,r_dom_abi,c_dom_abi,askyl,iskyl,nskyl,zz)
        else
           call par_schua2(&
                1_ip,nbvar,zi,pp,aib,abi,abb,r_dom_aib,c_dom_aib,r_dom_abb,&
                c_dom_abb,r_dom_abi,c_dom_abi,zz)
           do ii = 1,npoi1
              zirhs(ii) = zi(ii)
           end do
           call CSR_LUsol(&
                npoi1,nbvar,invpr_aii,invpr_aii,IL,JL,&  ! Solve Aii^{-1}
                LN,IU,JU,UN,zirhs,zi)
           call par_schua2(&
                2_ip,nbvar,zi,pp,aib,abi,abb,r_dom_aib,c_dom_aib,r_dom_abb,&
                c_dom_abb,r_dom_abi,c_dom_abi,zz)
        end if
     end if
     call cputim(time7)
     ctim1 = ctim1 + time7 - time6
     !
     ! alpha = rho^k / <p^{k+1},q^{k+1}>
     !
     alpha = 0.0_rp
     do ipoin = npoii2,npoii3
        kpoin = ipoin - npoi1
        alpha = alpha + pp(kpoin) * zz(kpoin)
     end do
     call PAR_SUM(alpha)

     if( alpha == 0.0_rp ) then
        ierro = 2
        goto 10
     end if
     rho   = newrho
     alpha = newrho / alpha
     !
     ! x^{k+1} = x^k + alpha*p^{k+1}
     ! r^{k+1} = r^k - alpha*q^{k+1}
     !
     do ii = 1,nrows
        xb(ii) = xb(ii) + alpha * pp(ii)
        rr(ii) = rr(ii) - alpha * zz(ii)
     end do
     !
     !  L z^{k+1} = r^{k+1} 
     !
     call cputim(time8)
     ctim2 = ctim2 + time8 - time7
     if( idprecon == SOL_DIAGONAL .or. idprecon == SOL_MOD_DIAGONAL ) then

        do ii = 1,nrows
           zz(ii) = invdiag(ii) * rr(ii)
        end do

     else if( idprecon == SOL_APPROXIMATE_SCHUR .or. idprecon == SOL_ABB ) then

        if( INOTMASTER ) then
           do ii = 1,nrows
              zz(ii) = rr(ii) 
           end do
           if( solve_sol(1) % kfl_scpre == 0 ) then
              call chosol(&
                   npoin_bb*nbvar,nskyl_p,iskyl_p,1_ip,askyl_p,zz,&
                   npoin_bb*nbvar,info)
           else          
              call CSR_LUsol(&
                   npoin_bb,nbvar,invpr_prec,invpr_prec,ILpre,JLpre,&
                   LNpre,IUpre,JUpre,UNpre,rr,zz)              
           end if
        end if

     end if

     call cputim(time9)
     ctim3 = ctim3 + time9 - time8

     newrho = 0.0_rp
     do ipoin = npoii2,npoii3
        kpoin = ipoin - npoi1
        newrho = newrho + zz(kpoin) * rr(kpoin)
     end do
     call PAR_SUM(newrho)
     !
     ! beta  = rho^k / rho^{k-1}  
     !
     beta = newrho / rho
     !
     ! p^{k+1} = z^k + beta*p^k - W.mu^k
     !
     do ii = 1,nrows
        pp(ii) = zz(ii) + beta * pp(ii)
     end do

     resid  = sqrt(newrho)
     resi2  = resi1
     resi1  = resid * invnb
     solve_sol(1) % iters = solve_sol(1) % iters + 1
     !
     ! Convergence post process:
     ! kk    = iteration number
     ! resi1 = preconditioned residual
     !
     if( kfl_cvgso == 1 ) &
          call outcso(dummr,dummr,dummr)

  end do
  
  solve_sol(1) % resfi = resi1

  !-----------------------------------------------------------------
  !
  !  END MAIN LOOP
  !
  !-----------------------------------------------------------------


  !-----------------------------------------------------------------
  !
  ! Solve interior nodes: Aii ui = bi - Aib ub
  !
  !-----------------------------------------------------------------

  do ipoin = 1,npoii1
     xi(ipoin) = bi(ipoin)
     do izdom = r_dom_aib(ipoin),r_dom_aib(ipoin+1)-1
        jpoin = c_dom_aib(izdom)
        xi(ipoin) = xi(ipoin) - aib(1,1,izdom) * xb(jpoin)
     end do
  end do
  if( INOTMASTER ) then
     if( solve_sol(1) % kfl_scaii == 0 ) then
        call chosol(&
             npoin_ii*nbvar,nskyl,iskyl,1_ip,askyl,xi,&
             npoin_ii*nbvar,info) 
     else
        do ipoin = 1,npoii1
           zi(ipoin) = xi(ipoin)
        end do
        call CSR_LUsol(&
             npoin_ii,nbvar,invpr_aii,invpr_aii,IL,JL,&  ! Solve Aii^{-1}
             LN,IU,JU,UN,zi,xi) 
     end if
  end if

  call cputim(time10)

  !-----------------------------------------------------------------
  !
  ! Final residual
  !
  !-----------------------------------------------------------------

  resfi = resi1
  if( INOTMASTER ) then
     if( solve_sol(1) % kfl_scaii == 0 ) then
        call par_schuax(&
             nbvar,zi,xb,aib,abi,abb,r_dom_aib,c_dom_aib,r_dom_abb,&
             c_dom_abb,r_dom_abi,c_dom_abi,askyl,iskyl,nskyl,rr)          
     else
        call par_schua2(&
             1_ip,nbvar,zi,xb,aib,abi,abb,r_dom_aib,c_dom_aib,r_dom_abb,&
             c_dom_abb,r_dom_abi,c_dom_abi,rr)
        do ii = 1,npoi1
           zirhs(ii) = zi(ii)
        end do
        call CSR_LUsol(&
             npoi1,nbvar,invpr_aii,invpr_aii,IL,JL,& ! Solve Aii^{-1}
             LN,IU,JU,UN,zirhs,zi)
        call par_schua2(&
             2_ip,nbvar,zi,xb,aib,abi,abb,r_dom_aib,c_dom_aib,r_dom_abb,&
             c_dom_abb,r_dom_abi,c_dom_abi,rr)
     end if
  end if
  do ii = 1, nrows
     rr(ii) = sb(ii) - rr(ii) 
  end do
  rnorm = 0.0_rp
  do ipoin = npoii2,npoii3
     kpoin = ipoin - npoi1
     rnorm = rnorm + rr(kpoin) * rr(kpoin)
  end do
  call PAR_SUM(rnorm)
  solve_sol(1) % resf2 = rnorm * invb2

10 continue

  if( solve_sol(1) % kfl_solve == 1 ) then
     if( ierro > 0 ) write(solve_sol(1) % lun_solve,120) iters
  end if

  !-----------------------------------------------------------------
  !
  ! Deallocate memory
  !
  !-----------------------------------------------------------------

  if( INOTMASTER ) then 

     call memchk(2_ip,istat,memit,'ZI','par_cgschu',zi)
     deallocate(zi,stat=istat)
     if(istat/=0) call memerr(2_ip,'ZI','par_cgschu',0_ip)

     call memchk(2_ip,istat,memit,'SB','par_cgschu',sb)
     deallocate(sb,stat=istat)
     if(istat/=0) call memerr(2_ip,'SB','par_cgschu',0_ip)

     call memchk(2_ip,istat,memit,'INVDIAG','par_cgschu',invdiag)
     deallocate(invdiag,stat=istat)
     if(istat/=0) call memerr(2_ip,'INVDIAG','par_cgschu',0_ip)

     call memchk(2_ip,istat,memit,'WW','par_cgschu',ww)
     deallocate(ww,stat=istat)
     if(istat/=0) call memerr(2_ip,'WW','par_cgschu',0_ip)

     call memchk(2_ip,istat,memit,'ZZ','par_cgschu',zz)
     deallocate(zz,stat=istat)
     if(istat/=0) call memerr(2_ip,'ZZ','par_cgschu',0_ip)

     call memchk(2_ip,istat,memit,'PP','par_cgschu',pp)
     deallocate(pp,stat=istat)
     if(istat/=0) call memerr(2_ip,'PP','par_cgschu',0_ip)

     call memchk(2_ip,istat,memit,'RR','par_cgschu',rr)
     deallocate(rr,stat=istat)
     if(istat/=0) call memerr(2_ip,'RR','par_cgschu',0_ip)

     if( solve_sol(1) % kfl_scaii == 0 ) then
        call memchk(2_ip,istat,memit,'ISKYL','par_cgschu',iskyl)
        deallocate(iskyl,stat=istat)
        if(istat/=0) call memerr(2_ip,'ISKYL','par_cgschu',0_ip)

        call memchk(2_ip,istat,memit,'ASKYL','par_cgschu',askyl)
        deallocate(askyl,stat=istat)
        if(istat/=0) call memerr(2_ip,'ASKYL','par_cgschu',0_ip)
     end if
     if( solve_sol(1) % kfl_scaii == 1 ) then
        call memchk(2_ip,istat,memit,'ZIRHS','par_cgschu',zirhs)
        deallocate(zirhs,stat=istat)
        if(istat/=0) call memerr(2_ip,'ZIRHS','par_cgschu',0_ip)
        if( INOTMASTER ) then         
           call CSR_LUfin(&
                IL,JL,LN,IU,JU,UN)      
        end if
     end if

     if( idprecon == SOL_APPROXIMATE_SCHUR .or. idprecon == SOL_ABB ) then

        call memchk(2_ip,istat,memit,'PREC','par_cgschu',prec)
        deallocate(prec,stat=istat)
        if(istat/=0) call memerr(2_ip,'PREC','par_cgschu',0_ip)

        if( solve_sol(1) % kfl_scpre == 0 ) then

           call memchk(2_ip,istat,memit,'ISKYL_P','par_cgschu',iskyl_p)
           deallocate(iskyl_p,stat=istat)
           if(istat/=0) call memerr(2_ip,'ISKYL_P','par_cgschu',0_ip)

           call memchk(2_ip,istat,memit,'ASKYL_P','par_cgschu',askyl_p)
           deallocate(askyl_p,stat=istat)
           if(istat/=0) call memerr(2_ip,'ASKYL_P','par_cgschu',0_ip)

        else

           if( INOTMASTER ) then         
              call CSR_LUfin(&
                   ILpre,JLpre,LNpre,IUpre,JUpre,UNpre)      
           end if

        end if

     end if

  end if
  !
  ! CPU time
  !
  solve_sol(1) % cpu_schur(1) = solve_sol(1) % cpu_schur(1) + time2  - time1 ! Invert Aii
  solve_sol(1) % cpu_schur(2) = solve_sol(1) % cpu_schur(2) + time3  - time2 ! Compute preconditioner P
  solve_sol(1) % cpu_schur(3) = solve_sol(1) % cpu_schur(3) + time4  - time3 ! Factorize preconditioner P
  solve_sol(1) % cpu_schur(4) = solve_sol(1) % cpu_schur(4) + time5  - time4 ! Solver initialization
  solve_sol(1) % cpu_schur(5) = solve_sol(1) % cpu_schur(5) + time10 - time5 ! Main loop
  solve_sol(1) % cpu_schur(6) = solve_sol(1) % cpu_schur(6) + ctim1          ! Main loop: q = A p = Abb p - Abi Aii^{-1} Aib p
  solve_sol(1) % cpu_schur(7) = solve_sol(1) % cpu_schur(7) + ctim2          ! Main loop: x = x + alpha*p
  solve_sol(1) % cpu_schur(8) = solve_sol(1) % cpu_schur(8) + ctim3          ! Main loop: L z = q

110 format(i5,18(2x,e12.6))
120 format(&
       & '# Error at iteration ',i6,&
       & 'Dividing by zero: alpha = rho^k / <p^{k+1},q^{k+1}>')

end subroutine par_cgschu
