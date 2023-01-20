!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



! 1)borro todo lo de own first
! 2)en unkwsmp2alya   habia una parte que decia
    ! Exchange - This was needed en AGMG because each subdomain had its own nodes ok (1,npoi1  & npoin2,npoin3)
    ! and it needed to pass teh values to the neighbors - I guess teh same will happen in WSMP
! esto lo he borrado pero puede que lo necesite.revisar manula de psblas o consultar Salvatore
! 3)  ojo hay muchas partes que se repiten con lo de agmg ver de unir todo aun mÃs
!  por ejemplo   call mat2mathal  &call PAR_INTERFACE_MATRIX_W_HALOS_EXCHANGE
! falta copiar el valor de salida a unkno
!
!

#ifdef INC_PSBLAS



! 2check
! 1) ./mod_alya2wsmp.f90:274:       allocate(aa_wsmp(nd_aux,nd_aux,pt % nz_2)) 
! 2) unknowns vectorial case - how are they ordered  ?   !1x,1y,2x,2y ?  Altrenativelly it could be 1x,2x....1y,2y
!  Morever it could be in a vector rhs(ndime*npoin) or rhs(ndime,npoin)  or rhs(npoin,ndime)  ! Chech
! I will suppose rhs(ndime*npoin)  ordered   !1x,1y,2x,2y   -- see unkalya2agmg
! 3) I trust that your solver does not change miter , ia, ja  ! AGMG did alter them -- in that case you have to make an auxiliary copy!!
! 4) for WSMP from what I saw in the manual, I belive that the global numbering should be changed so that all all nodes from teh fisrt subdomain are first, tehn those from teh second one ...
!   In AMGX I do not know if this is necesary. It is something I have not done.


!>   
!------------------------------------------------------------------------
!> @addtogroup Algebraic_Solver
!> @{
!> @name    alya2psblas
!> @file    mod_alya2psblas.f90
!> @author  Herbert Owen
!> @date    12/04/2016
!> @brief   Bridge to PSBLAS
!> @details Module to call PSBLAS iterative solver
!------------------------------------------------------------------------
module mod_alya2psblas

  use def_kintyp,         only : ip,rp,lg   !To revise what is actually needed
  use def_kintyp_solvers, only : soltyp
  use mod_memory,         only : memory_alloca
  use mod_memory,         only : memory_deallo ! ver que este deallocando todo lo que no necesito
  use def_solver,         only : memit
  use def_master,         only : comm_data_par 
  use def_domain,         only : mesh_type
  use def_master,         only : kfl_paral
  use def_master,         only : lninv_loc
  use mod_communications, only : PAR_MAX
  use mod_alya2agmg,      only : agmg_pt_nd,mathal2matfinal,calc_lexplode_ndofn

  use psb_const_mod,      only : itwo,psb_root_
  use psb_base_mod
  use psb_prec_mod
  use psb_krylov_mod
  use psb_util_mod
  use mld_prec_mod
  
  implicit none
  private

  public :: alya2psblas_initialization
  public :: alya2psblas_solution
 
  
  type psblas_pt
     integer(ip)          :: npoin    ! added to avoid needing to pass meshe to alya2psblas_solution
     integer(ip)          :: npoin_own
     integer(ip)          :: npoin_2  ! added to avoid needing to pass meshe to alya2psblas_solution
     integer(ip)          :: npoi1    ! added to avoid needing to pass meshe to alya2psblas_solution
     integer(ip)          :: npoi2    ! added to avoid needing to pass meshe to alya2psblas_solution
     integer(ip)          :: npoi3    ! added to avoid needing to pass meshe to alya2psblas_solution
     integer(ip)          :: nz_own
     integer(ip)          :: nz_2
     integer(ip), pointer :: lwhalo_inv(:)
  end type psblas_pt

  integer(ip),parameter   :: kfl_debug=0
  
  type(psblas_pt)         :: pt

  ! parallel environment
  integer(psb_ipk_) :: ictxt, iam, np
  !psblas
  integer(psb_ipk_) :: info
  ! descriptor
  type(psb_desc_type)   :: desc_a

  ! sparse matrix and preconditioner
  type(psb_dspmat_type) :: apsblas
  type(psb_dprec_type)  :: prec
  type(mld_dprec_type)  :: mlprec
  ! dense vectors
  type(psb_d_vect_type) :: xpsblas,bpsblas
  ! input parameters
  character(len=20) :: kmethd, ptype

  ! solver parameters
  integer(psb_ipk_) :: iter, itmax,itrace, istopc, irst
  !    integer(psb_long_int_k_) :: amatsize, precsize, descsize, d2size
  real(psb_dpk_)    :: err, eps

  integer(psb_ipk_) :: n_nz,n_np
  
  integer(psb_ipk_),allocatable          :: jrow(:)
  integer(psb_ipk_),allocatable          :: irow(:)
  integer(psb_ipk_),allocatable          :: irw(:)
  integer(psb_ipk_),allocatable          :: vl3d(:)
  integer(ip),allocatable          :: vl3d_aux(:)

  type(agmg_pt_nd), pointer   :: ptnd(:)

  logical :: init_desc
  logical, parameter :: dump_data=.false.

contains 

  !----------------------------------------------------------------------- 
  !
  !> @author  Herbert Owen
  !> @date    12/04/1016
  !> @brief   Initialize PSBLAS 
  !> @details Initialize PSBLAS for a given problem
  !>          The ia and ja I need here are from the graph with halos only when I receive the matrix and wanto to build las
  !
  !-----------------------------------------------------------------------

  subroutine alya2psblas_initialization(meshe,solve)        ! (meshe,COMM,solve)

    use def_master, only               : INOTMASTER
    use def_master, only               : ISEQUEN
    use mod_alya2agmg, only            : calc_lmat2mathal,max_ndofn_agmg_sol
    use mod_parall, only               : PAR_COMM_MY_CODE_WM
    use mod_communications_tools, only : PAR_COMM_TO_INT
    ! en maphys eran de intent inout y ademas estab COMM de intent in
    type(mesh_type),     intent(in) :: meshe            !< Mesh type     ! en maphys era inout
    type(soltyp),        intent(in) :: solve            !< Solver type   ! en maphys era inout
    
    integer(psb_ipk_)          :: nz2_aux,nd_aux,iz,ipoin,idime,i
    integer(ip)                :: nz2_aux_ip,nd_aux_ip
    integer(ip), save          :: ipass = 0

    character(len=20)          :: name

    integer(4)                 :: PAR_COMM_TO_USE


    if ( ISEQUEN ) return

    if( INOTMASTER ) then
       !
       ! mod_par_additional_arrays: se construyen los arrays de intercambio: ja, nzdom,nzdom_ii para cada neighbor   
       ! ESTO SE LLAMA DESDE  Domain.f90 L180   por ahora le meto con #ifdef PARAL_AGMG
       ! buscar bound_mat_
       ! tambien hay unos nullify en modu_parall   -- estos los dejo que se hagan siempre
       !
       nd_aux = solve % ndofn 
       nd_aux_ip = solve % ndofn 

       if( ipass == 0 ) then
          ipass=1
          nullify(pt % lwhalo_inv)

          nullify(ptnd)
          allocate(ptnd(max_ndofn_agmg_sol))  ! por ahoar usare el valor de agmg -- en algun momento hay que limpiar bien y no repetir
          do i=1,max_ndofn_agmg_sol
             nullify(ptnd(i) % ia_ndofn)
             nullify(ptnd(i) % ja_ndofn)
             nullify(ptnd(i) % lndofn_id)
             nullify(ptnd(i) % lndofn_jd)
             nullify(ptnd(i) % lndofn_iz)
             nullify(ptnd(i) % lfinal_rhs)
             nullify(ptnd(i) % lfinal_rhs_inv)
          end do
       end if
       !
       ! Obtain things that do not depend on ndofn
       !
       if(.not.associated(pt % lwhalo_inv)) then

          call memory_alloca(memit,'PT % LWHALO_INV','alya2agmg_initialization',pt % lwhalo_inv,solve % nnz)
          !
          ! obtain lmat2mathal
          ! meshe % r_dom_2,meshe % c_dom_2  have been obtained in par_matrix_w_halos_exchange_on_interface_nodes graphs_poipoi  
          !
          nz2_aux = meshe % r_dom_2(meshe % npoin_2 + 1) - 1
          pt % nz_2 = nz2_aux   ! stored here to have it for alya2agmg_solution
          pt % npoin = meshe % npoin
          pt % npoin_2 = meshe % npoin_2
          pt % npoi1 = meshe % npoi1
          pt % npoi2 = meshe % npoi2
          pt % npoi3 = meshe % npoi3
          !
          ! Premutation array lwhalo_inv of matrix without to matrix with halo
          !
          nz2_aux_ip = nz2_aux
          call calc_lmat2mathal(meshe % npoin,meshe % nzdom ,meshe % npoin_2, & 
               nz2_aux_ip,meshe % r_dom,meshe % c_dom, &
               meshe % r_dom_2,meshe % c_dom_2,pt % lwhalo_inv)
       
          pt % nz_own = meshe % r_dom_2(meshe % npoi3 +1) -1
          pt % npoin_own = meshe % npoi3
       end if
       !
       ! from ppde3d
       !
       PAR_COMM_TO_USE = PAR_COMM_TO_INT(PAR_COMM_MY_CODE_WM)   ! communicator without master
       call psb_init(ictxt, basectxt=PAR_COMM_TO_USE)

       call psb_info(ictxt,iam,np)

       if (iam < 0) then 
          ! This should not happen, but just in case
          call psb_exit(ictxt)
          stop
       endif
       if(psb_get_errstatus() /= 0) stop ' error in mod_alya2psblas 1'
       name='alya'
       call psb_set_errverbosity(itwo)
       call psb_cd_set_large_threshold(itwo)

       !
       if (iam == psb_root_) then 
          write(*,*) 'Welcome to PSBLAS version: ',psb_version_string_
          write(*,*) 'This is the ',trim(name),' sample program',iam,np
       end if
       !
       !  get parameters  obtenerlos yo
       !
       !       call get_parms(ictxt,kmethd,ptype,afmt,idim,istopc,itmax,itrace,irst)
       ! afmt e idim son ppde3d
!       ptype='DIAG'
       ptype='NONE'
       if (nd_aux == 1) then   ! solucion temporal
          kmethd = 'CG'
       else
          kmethd = 'GMRES'   ! me dijo que no conoce GMRES y puso BICGStab ya me va bien despeus revisar bien que hay
       end if
       istopc = 2   ! 1: use the normwise backward error, 2: use the scaled 2-norm of the residual. Default: 2
       itrace=1     ! info every itrace iter - o no info
       itmax=1000
       irst=50      ! restart -  BiCGSTABL or RGMRES methods, otherwise it is ignored.

       ! Los punto 1 y 2 se hacen una sola vez  por eso pase a init --si se resolvieran varios problemas con psblas habría que repensar
       ! Esto lo saque de lo que mando Salvatore por email
       !
       !1. Initialize the communication descriptor
       !
       n_nz = pt % nz_own * nd_aux * nd_aux
       n_np = pt % npoin_own * nd_aux
       allocate(jrow(n_nz))
       allocate(irow(n_nz))
       allocate(vl3d(n_np))  ! perhaps use directly lninv_loc(1:n_np)
       allocate(irw(n_np))

       do ipoin = 1, n_np/nd_aux
          do idime = 1,nd_aux
             irw  ( idime + (ipoin-1)*nd_aux ) = idime + (ipoin-1)*nd_aux
             vl3d ( idime + (ipoin-1)*nd_aux ) = idime + (lninv_loc(ipoin)-1)*nd_aux
          end do
       end do

       !
       ! Obtain things that depend on ndofn - borrowed from agmg
       !
       if(.not.associated(ptnd(nd_aux) % lndofn_id)) then

          call memory_alloca(memit,'PTND(SOLVE%NDOF) % IA_NDOFN','alya2agmg_initialization',ptnd(nd_aux) % ia_ndofn, &
               pt % npoin_own * nd_aux_ip + 1)
          call memory_alloca(memit,'PTND(SOLVE%NDOF) % JA_NDOFN','alya2agmg_initialization',ptnd(nd_aux) % ja_ndofn, &
               pt % nz_own * nd_aux_ip**2_ip)
          call memory_alloca(memit,'PTND(SOLVE%NDOF) % LNDOFN_ID','alya2agmg_initialization',ptnd(nd_aux) % lndofn_id, &
               pt % nz_own * nd_aux_ip**2_ip)
          call memory_alloca(memit,'PTND(SOLVE%NDOF) % LNDOFN_JD','alya2agmg_initialization',ptnd(nd_aux) % lndofn_jd, &
               pt % nz_own * nd_aux_ip**2_ip)
          call memory_alloca(memit,'PTND(SOLVE%NDOF) % LNDOFN_IZ','alya2agmg_initialization',ptnd(nd_aux) % lndofn_iz, &
               pt % nz_own * nd_aux_ip**2_ip)

!          call calc_lexplode_ndofn(nd_aux, pt % npoin_own, pt % nz_own, meshe % r_dom_2, pt % ja_yvan, &    ! cambio ja_yvan a c_dom_2  revsar si ok
          call calc_lexplode_ndofn(nd_aux_ip, pt % npoin_own, pt % nz_own, meshe % r_dom_2, meshe % c_dom_2, &
               ptnd(nd_aux) % ia_ndofn,   ptnd(nd_aux) % ja_ndofn, ptnd(nd_aux) % lndofn_id, &
               ptnd(nd_aux) % lndofn_jd,  ptnd(nd_aux) % lndofn_iz)

       end if

       do ipoin = 1, n_np
          do iz = ptnd(nd_aux) % ia_ndofn(ipoin), ptnd(nd_aux) % ia_ndofn(ipoin+1)-1
             irow(iz) = ipoin
          end do
       end do
       
       jrow(1:n_nz) = ptnd(nd_aux) % ja_ndofn(1:n_nz)
       !jrow(1:n_nz) = meshe % c_dom_2(1:n_nz)


       block
         integer(psb_ipk_)              :: iam, np, nr, nc
         real(rp)                       :: bnrm2, xnrm2, anrmi
         real(rp), allocatable          :: vv(:)
!         integer(psb_ipk_), allocatable :: icol(:)  ! estaba al pedo
         character(80)                  :: fname
         if (dump_data) then 
           call psb_info(ictxt,iam,np)
           write(*,*) iam,np,' Before CDALL'
           write(*,*) iam,np, 'Checks: ',size(jrow), size(vl3d),&
                & maxval(jrow(1:n_nz)),size(pt%lwhalo_inv), size(lninv_loc)
           call psb_barrier(ictxt)
           nr = size(vl3d)
           nc = maxval(jrow(1:n_nz))
           
           write(fname,'(a,i2.2,a)') 'glb-j-',iam,'.mtx' 
           call mm_array_write(int(lninv_loc(1:nc),psb_ipk_),'Alya global col indices',info,filename=fname)
         end if
       end block

       !
       ! An explanation is in order here.
       ! A descriptor can be built in two ways.
       ! One is implicitly through the use of psb_spins.
       ! One is explicitly through the use of psb_cdins.
       ! When (later in the program) you call psb_spins with local=.true.
       ! you are assuming that the descriptor has *already*
       ! been fully built and  assembled; in other words, the implicit
       ! build of the descriptor through psb_spins does not work if you
       ! call psb_spins with local=.true.
       ! The solution as implemented below is to call psb_cdins
       ! and psb_cdall here, since you already have global and local indices
       ! sorted out. 
       ! The second version is even better, because we do away with IA and only use
       ! column indices IGCOL and LIDX, with much smaller vectors; waiting for tests.
       ! After cdasb, you can then go on to call spins with local=.true.
       !
       init_desc = .true.
       call psb_cdall(ictxt,desc_a,info,vl=vl3d(:))
       block
         integer(psb_ipk_) :: nr, nc, i, ic, ic_aux
         integer(psb_ipk_), allocatable :: igrow(:), igcol(:), lidx(:)
         if (.false.) then
            if(nd_aux/=1) print*,'alya2psblas: this part is not ready for nd-aux/=0'
            if(nd_aux/=1) call runend('alya2psblas: this part is not ready for nd-aux/=0')
            allocate(igrow(n_nz),igcol(n_nz))
            igrow(1:n_nz) = lninv_loc(irow(1:n_nz))
            igcol(1:n_nz) = lninv_loc(jrow(1:n_nz))
            call psb_cdins(n_nz,igrow,igcol,desc_a,info)
         else
            nr = size(vl3d)
            nc = maxval(jrow(1:n_nz))
            allocate(vl3d_aux(nc))
            do ipoin = 1, nc/nd_aux
               do idime = 1,nd_aux
                  vl3d_aux ( idime + (ipoin-1)*nd_aux ) = idime + (lninv_loc(ipoin)-1)*nd_aux
               end do
            end do
            igcol = vl3d_aux(1:nc)             
            lidx  = (/(i,i=1,nc)/)
!            write(kfl_paral+1100,*)lidx
!            write(kfl_paral+1200,*)jrow
!            write(kfl_paral+1300,*)irow
!            write(kfl_paral+1400,*)irw
!            write(kfl_paral+1500,*)vl3d
            call psb_cdins(nc,igcol,desc_a,info,lidx=lidx)
         end if

       end block
       call psb_cdasb(desc_a,info)
       
     
    end if   ! inotmaster

  end subroutine alya2psblas_initialization

  
  !-----------------------------------------------------------------------
  !
  !> @author  Herbert Owen
  !> @date    12/04/1016
  !> @brief   Solve Ax=b
  !> @details Solve algebraic system using PSBLAS.
  !>          
  !
  !-----------------------------------------------------------------------
  subroutine alya2psblas_solution(solve,amatr,rhsid,unkno)
    use mod_communications, only : PAR_INTERFACE_MATRIX_W_HALOS_EXCHANGE
    use mod_parall,         only : PAR_COMM_MY_CODE_ARRAY
    use def_master,         only : INOTMASTER,ISEQUEN,ISLAVE,IMASTER

    use def_domain,         only :  meshe  ! for residual with alya matrix & also for ia & ja in serial version
    use def_kermod,         only :  ndivi

    use def_solver,         only :  SOL_SOLVER_AGMG, SOL_GAMGX

    use mod_alya2agmg,      only :  mat2mathal
    use mod_communications, only :  PAR_INTERFACE_NODE_EXCHANGE

    type(soltyp),                intent(inout) :: solve                !< Solver type

    !  OJO inivar  L119  :  solve_sol(1) % nnz  =  meshe % nzdom
    real(rp),                    intent(in)    :: amatr(solve % ndofn,solve % ndofn,solve % nnz) !< Matrix
    real(rp),                    intent(inout) :: rhsid(solve % ndofn*solve % nequa) !< RHS
    real(rp),                    intent(inout) :: unkno(solve % ndofn*solve % nequa) !< Solution

    real(psb_dpk_),allocatable    :: rhsid_auxt(:),unkno_auxt(:)  ! Missin see if thsi is is actually necesarry
    real(psb_dpk_),allocatable    :: temp(:)
    real(rp),pointer              :: temp2(:)
    real(rp),pointer              :: temp3(:,:)

    real(rp),pointer              :: aa_hal(:,:,:)   ! no se si habra problema que sea pointer estando dentro de un modulo
    real(psb_dpk_),allocatable    :: aa_psblas(:)   ! For the moment I will leave it without pointer  ! ojo puse rp pero puede que sea psb_dpk_  no dio problema con dpk
    

    integer(psb_ipk_)  :: ipoin,hh_aux,idime
    integer(psb_ipk_)  :: nd_aux,niter_aux
    integer(ip)        :: nd_aux_ip
    real(rp)     :: auxtol
    ! Timers
    real(psb_dpk_) :: t1, t2, tprec, tslv 

    real(rp),allocatable    :: rr(:)
    real(rp)                :: resif,rhsnorm,resii
    integer(ip), parameter  :: kfl_verbor = 0_ip      ! 10 print all   , 0 nothing 
    logical, parameter      :: use_mldprec=.true.

    nullify(temp2)
    nullify(temp3)

    niter_aux = solve % miter ! I leave it in case AMGX-Vishal needs it
    !
    ! rhsnorm
    !
    nd_aux = solve % ndofn  ! guess the master must have some correct value
    nd_aux_ip = nd_aux
    call norm2x(nd_aux_ip,rhsid,rhsnorm)
    !
    ! initial residual & auxtol
    !
    auxtol = solve % solco

    allocate(rr(nd_aux*max(1,meshe(ndivi) % npoin)))
    rr = 0.0_rp
    ! rr = amatr * x
    call bcsrax( 1_ip, meshe(ndivi) % npoin, nd_aux, amatr, meshe(ndivi) % c_dom, meshe(ndivi) % r_dom, unkno, rr )
    if( INOTMASTER ) rr = rr - rhsid
    


    call norm2x(nd_aux_ip,rr,resii)
    deallocate(rr)
    if (solve % kfl_adres == 1_ip) auxtol = max(auxtol,solve % adres*resii/rhsnorm) ! adaptive

    if( ISLAVE ) then

      ! so that there is no problem when passing to agmg - This will imply using palin allocate and deallocate
      ! Guillueme said that if there is no interface there should be no problem in defining it as pointer and recieving
      ! it as a normal array , but just in case
      !
      ! I will first transform amatr to amatr_hal  and then add the contribution fron neigbours.
      ! then I will transform al that in to final form that AGMG needs
      ! it can be done all in one step but I believe it is more complicated.
      ! I wil try to have it working this way an then do the other option.
      !
      ! Allocate matrices - later I can allocate them somewere else and only once for the whole progarm
      !
      ! Ok, so me adding init_desc was not required. However the code below shows that:
      ! 1. If you have a matrix which stays the same, you do not need to do anything,
      ! 2. If you have a matrix which changes *values* but *not* the pattern, then you can
      !    use a%reinit as below
      ! 3. If it is a different pattern, or a different size, you have to restart the descriptor.    
      !
      if (init_desc) then
        call psb_spall(apsblas,desc_a,info)
        call psb_geall(bpsblas,desc_a,info)
        call psb_geall(xpsblas,desc_a,info)
        init_desc = .false.
      else
        call apsblas%reinit(clear=.true.)
        call xpsblas%zero()
        call bpsblas%zero()
      endif

      nullify(aa_hal)
      call memory_alloca(memit,'AA_HAL','alya2psblas_solution',aa_hal,nd_aux_ip,nd_aux_ip,pt % nz_2)    ! ojo aca habia puesto nz_own en lugar de nz2 y la cagaba
      aa_hal=0.0_rp
      allocate(aa_psblas(nd_aux * nd_aux * pt % nz_own))
      !
      ! Pass the matrix aa to a matrix aa_hal  - In alya2psblas_initialization
      !            I create a pointer that makes this automatic ej lmat2mathal 
      !
      call mat2mathal(solve % nnz,pt % nz_2,nd_aux_ip,pt % lwhalo_inv,amatr,aa_hal)
      !
      ! Add the contibution lo the matrix with halos aa_hal from the neighbors
      ! there is a call to call PAR_INTERFACE_MATRIX_EXCHANGE in solope tu use as example
      !
      call PAR_INTERFACE_MATRIX_W_HALOS_EXCHANGE(nd_aux_ip,amatr,aa_hal,PAR_COMM_MY_CODE_ARRAY(1))
      !
      ! Now convert to the final format also reorder rhs and initial guess - borrowed from agmg
      !
      call mathal2matfinal  (pt % nz_own,pt % nz_2,nd_aux_ip,ptnd(nd_aux) % lndofn_id,ptnd(nd_aux) % lndofn_jd, &   !volver a ponerlo
           ptnd(nd_aux) % lndofn_iz,aa_hal,aa_psblas)
      !
      ! Solve - call PSBLAS
      !
      !
      ! 3. On each process, read from file into memory triples
      !IROW(1:nz), ICOL(1:nz), VAL(1:nz)
      ! here I need t create irow fron ia   easy  perhaps do it in the inic phase does not change
      !

      ! local Whether the entries in the indices vectors ia, ja are already in local
      !       call psb_spins(nz,irow,icol,val,a,desc_a,info)
      ! irow icol val son de entrada integer y real
      allocate(rhsid_auxt(solve % ndofn*solve % nequa))
      rhsid_auxt(1:solve % ndofn*solve % nequa) = rhsid(1:solve % ndofn*solve % nequa)
      call psb_spins(n_nz,irow,jrow,aa_psblas,apsblas,desc_a,info,local=.true.)   !PTND(SOLVE%NDOF) % IA_NDOFN'
      deallocate(aa_psblas)
      !
      ! 4. Read the RHS b(1:nl) on process 0 (1:n1 on 1 etc.) and the global row indices IR(1:nl)
      ! Here since I will give all the matrix at once IR = 3DVL
      !

      call psb_geins(n_np,irw(1:n_np),rhsid_auxt(1:n_np),bpsblas,desc_a,info,local=.true.)
      deallocate(rhsid_auxt)

      !
      ! perhaps create unkno_aux -not sure 
      !
      allocate(unkno_auxt(solve % ndofn*solve % nequa))
      unkno_auxt(1:solve % ndofn*solve % nequa) = unkno(1:solve % ndofn*solve % nequa)
      call psb_geins(n_np,irw(1:n_np),unkno_auxt(1:n_np),     xpsblas,desc_a,info,local=.true.)
      deallocate(unkno_auxt)
      !
      ! para ir tengo 2 opciones o uso el mismo 3DVL  que en psb_cdall   con lo cual no pongo local
      ! o pongo local e ir es la identidad
      ! from manual call psb_geins(m, irw, val, x, desc_a, info [,dupl,local])   ! es como que le falta x a lo que em dijo Salvatore x es el unico intent inout
      ! el resto son de intent in
      ! lo mismo pasa con psb_spins donde solo a,desc_a son inout 
      ! con lo cual le puedo pasar aa_hal  y rhs que no lso av atocar  
      !
      !5 . Assembly
      !
      call psb_spasb(apsblas,desc_a,info)
      call psb_geasb(bpsblas,desc_a,info)
      call psb_geasb(xpsblas,desc_a,info)
      block
        integer(psb_ipk_) :: iam, np, ictxt, nr
        real(rp)    :: bnrm2, xnrm2, anrmi
        real(rp), allocatable :: vv(:)
        character(80) :: fname
        if (dump_data) then 
          ictxt = desc_a%get_ctxt()
          call psb_info(ictxt,iam,np)
          write(*,*) iam,np,' Assembly before call to Krylov method'
          call psb_barrier(ictxt)
          write(*,*) iam,np,' Sanity tests'
          write(*,*) iam,np,' Global/local rows/cols',&
               & desc_a%get_global_rows(),desc_a%get_global_cols(),&
               & desc_a%get_local_rows(),desc_a%get_local_cols()
          write(*,*) iam,np,' Matrix sizes',apsblas%get_nrows(),apsblas%get_ncols()
          write(*,*) iam,np,' Vector sizes',xpsblas%get_nrows(),bpsblas%get_nrows()
          call psb_barrier(ictxt)                  
          bnrm2 = psb_genrm2(bpsblas,desc_a,info)
          xnrm2 = psb_genrm2(xpsblas,desc_a,info)
          anrmi = psb_spnrmi(apsblas,desc_a,info)
          write(*,*) iam,np,' Norms: ',bnrm2, xnrm2, anrmi
          nr = desc_a%get_local_rows()
          write(fname,'(a,i2.2,a)') 'amat-',iam,'.mtx'
          call apsblas%print(fname,head='Test mat Alya')
          write(fname,'(a,i2.2,a)') 'xv-',iam,'.mtx'
          vv = xpsblas%get_vect()
          call mm_array_write(vv(1:nr),'Alya X Vector',info,filename=fname)
          write(fname,'(a,i2.2,a)') 'bv-',iam,'.mtx'
          vv = bpsblas%get_vect()
          call mm_array_write(vv(1:nr),'Alya B Vector',info,filename=fname)
          write(fname,'(a,i2.2,a)') 'irow-',iam,'.mtx'
          call mm_array_write(irow,'Alya IROW indices',info,filename=fname)
          write(fname,'(a,i2.2,a)') 'jrow-',iam,'.mtx'
          call mm_array_write(jrow,'Alya JROW indices',info,filename=fname)
          write(fname,'(a,i2.2,a)') 'irw-',iam,'.mtx'
          call mm_array_write(irw,'Alya IRW indices',info,filename=fname)
          write(fname,'(a,i2.2,a)') 'vl3d-',iam,'.mtx'
          call mm_array_write(vl3d,'Alya VL3D indices',info,filename=fname)
          write(*,*) iam,np,'Checks on IROW/JORW:',&
               & minval(irow),maxval(irow),minval(jrow),maxval(jrow)
          call psb_barrier(ictxt)
          write(*,*) iam,np,' Sanity tests done'
        end if
      end block


      !6. Choose the preconditioner to be used with psb_precset and build it with psb_precbld
      ! no encontre ningun ejemplo de psb_precset - en spde3d.f90 se usa psb_precinit - lo copio
      if (use_mldprec) then
        if(kfl_verbor > 0_ip .and. iam == psb_root_) write(psb_out_unit,'("Setting preconditioner to : ",a)') 'ML'
        call mlprec%init(ictxt,'ML',info)   ! new libraries from git
        !call mlprec%init('ML',info)        !  salvatore compiled lib
        call mlprec%set('ml_cycle', 'WCYCLE',    info)
        call mlprec%set('outer_sweeps',  1,info)
        call mlprec%set('smoother_type',   'FBGS',     info)
        call mlprec%set('smoother_sweeps', 2,    info)
        call mlprec%set('coarse_solve',  'BJAC',    info)
        call mlprec%set('coarse_subsolve',  'ILU',    info)
        call mlprec%set('coarse_mat',  'DIST',    info)
        call mlprec%set('coarse_fillin', 1,    info)
        call mlprec%set('coarse_sweeps', 4,    info)
        call psb_barrier(ictxt)
        t1 = psb_wtime()
        call mlprec%hierarchy_build(apsblas,desc_a,info)
        if (info == psb_success_) call mlprec%smoothers_build(apsblas,desc_a,info)
        !all psb_barrier(ictxt)
        t2 = psb_wtime()
        tprec = t2 - t1
        if(info /= psb_success_) stop 'error mld_precbld'
        if(kfl_verbor > 0_ip .and. iam == psb_root_) write(psb_out_unit,'("Preconditioner build time : ",f20.4)') tprec
        call mlprec%descr(iout=psb_out_unit)
      else
        if(kfl_verbor > 0_ip .and. iam == psb_root_) write(psb_out_unit,'("Setting preconditioner to : ",a)')ptype
        !call psb_precinit(prec,ptype,info) ! salvatore compiled lib
        call prec%init(ictxt,ptype,info)    ! new libraries from git
        if(info /= psb_success_) stop 'error psb_precinit'
        
        call psb_barrier(ictxt)
        t1 = psb_wtime()
        call psb_precbld(apsblas,desc_a,prec,info)
        t2 = psb_wtime()
        tprec = t2 - t1
        if(info /= psb_success_) stop 'error psb_precbld'
        if(kfl_verbor > 0_ip .and. iam == psb_root_) write(psb_out_unit,'("Preconditioner build time : ",f20.4)') tprec       
      end if

      !7. Call the iterative method of choice, e.g. psb_bicgstab 
      !call psb_set_debug_level(9999)
      eps = 1.d-9
      call psb_barrier(ictxt)
      t1 = psb_wtime()
      if (use_mldprec) then
         call psb_krylov(kmethd,apsblas,mlprec,bpsblas,xpsblas,eps,desc_a,info,&
              & itmax=itmax,iter=iter,err=err,itrace=itrace,istop=istopc,irst=irst)
      else
         call psb_krylov(kmethd,apsblas,prec,bpsblas,xpsblas,eps,desc_a,info,&
              & itmax=itmax,iter=iter,err=err,itrace=itrace,istop=istopc,irst=irst)
      end if     
      t2 = psb_wtime()
      tslv = t2 - t1
      if(info /= psb_success_) then 
        print*,info
        stop 'error psb_krylov'  
      end if
      if(kfl_verbor > 0_ip .and. iam == psb_root_) write(psb_out_unit,'("Krylov solve         time : ",f20.4)') tslv
      if(kfl_verbor > 0_ip .and. iam == psb_root_) write(psb_out_unit,'("  linear system size      : ",I20)') desc_a%get_global_rows()
      if(kfl_verbor > 0_ip .and. iam == psb_root_) write(psb_out_unit,'("  estimator on exit       : ",g21.14)') err
      !       if(iam == psb_root_) print*,'xpsblas',xpsblas
      if(kfl_verbor > 0_ip .and. iam == psb_root_) print*,'info,iter,err',info,iter,err    !,cond
      call memory_deallo(memit,'AA_HAL'     ,'alya2agmg_solution',aa_hal) 
      !
      ! NEW
      !
      temp = xpsblas%get_vect()
      allocate(temp2(nd_aux * meshe(ndivi) % npoin))
      allocate(temp3(nd_aux,meshe(ndivi) % npoin))
      temp2 = 0.0_rp
      temp2(1:n_np) = temp(1:n_np)
      hh_aux = n_np/nd_aux
      !      temp3(nd_aux,1:hh_aux) = reshape(temp2(1:n_np), (/nd_aux,1:hh_aux/))   ! no  se porque no funciona lo hago amano
      do ipoin = 1,meshe(ndivi) % npoin
         do idime = 1,nd_aux
            temp3(idime,ipoin) = temp2(idime + (ipoin-1) * nd_aux)
         end do
      end do        
      call PAR_INTERFACE_NODE_EXCHANGE(temp3,'SUM')
      !      temp2(1:n_np) = reshape(temp3(nd_aux,1:hh_aux), (/n_np/) )
      do ipoin = 1,meshe(ndivi) % npoin
         do idime = 1,nd_aux
            temp2(idime + (ipoin-1) * nd_aux) = temp3(idime,ipoin)
         end do
      end do     
      unkno(1:meshe(ndivi) % npoin * nd_aux) = temp2(1:meshe(ndivi) % npoin * nd_aux)
      deallocate(temp2)
      deallocate(temp3)     
      ! Should I deallocate temp here?

    else if ( ISEQUEN ) then


    end if

    !
    ! Calculate final residual with alya matrix - This serves as a check
    !
    allocate(rr(max(1,nd_aux*meshe(ndivi) % npoin)))
    rr = 0.0_rp
    ! rr = amatr * x
    call bcsrax( 1_ip, meshe(ndivi) % npoin, nd_aux, amatr, meshe(ndivi) % c_dom, meshe(ndivi) % r_dom, unkno, rr )
    if( INOTMASTER ) rr = rr - rhsid
    call norm2x(nd_aux_ip,rr,resif)
    deallocate(rr)

    solve % iters = int(iter,ip)
    if( IMASTER ) err = 0.0_rp
    call PAR_MAX(err)
    solve % xorth = err ! beware I will output here the value of teh error given by psblas it should coincide with the one calculated by alya
    solve % resin = 0.0_rp ! I do not have a preconditioned residual
    solve % resi2 = resii/rhsnorm    ! this is the one that typically appears in Alya - Column 4 - Non preconditioned initial residual 
    solve % resfi = 0.0_rp ! resif
    solve % resf2 = resif/rhsnorm
    call PAR_MAX(solve % iters)

    if( kfl_verbor > 0_ip .and. INOTMASTER .and. iam == psb_root_) print*,'solve % resin,solve % resi2,solve % resfi,solve % resf2,solve % iters',&
         solve % resin,solve % resi2,solve % resfi,solve % resf2,solve % iters
  end subroutine alya2psblas_solution


end module mod_alya2psblas

#endif

