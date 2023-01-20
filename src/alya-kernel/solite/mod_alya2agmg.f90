!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!
! ja_yvan luego entra en call calc_lexplode_ndofn
! listrank is used by  calc_listrank_ndofn to obtain listrank_ndofn
!
! como ya no se necesita mas lo de own first  al   llamar  calc_lmatyvan  en lugar de pasarle ia_ownfirst & ja_ownfirst
!se pasa ia2 y ja2  o sea   meshe % r_dom_2,meshe % c_dom_2
!
! Eliminate l_rhs_ownfirst & l_rhs_inv_ownfirst  they are just teh identity now ! ELIMINATED
!
! Notice that in line 1022 I have        !ELIMINATED
!        lndofn_rhs(jj) = jj 
!        lndofn_rhs_inv(jj) = jj
! This means a lot of usesless code - before those lines it was beeing calculated in someother way - erase
! Then it appears in calc_lfinal_rhs just making thing confusing- it needs to be elimitaded from all of this file!!!
! 
!



!
! Al hacer svn commit se queja - pero es porque yo lo hice allocatable no pointer
! si es allocatable no tiene sentido nullify
!
!----WARNING!! 2[Priority]: The pointer aa_agmg shoud be nullified: nullify(aa_agmg);
!----WARNING!! 3[Priority]: The pointer rhs_agmg shoud be nullified: nullify(rhs_agmg);
!----WARNING!! 4[Priority]: The pointer unk_agmg shoud be nullified: nullify(unk_agmg);
!----WARNING!! 5[Priority]: The pointer ia_copy shoud be nullified: nullify(ia_copy);
!----WARNING!! 6[Priority]: The pointer ja_copy shoud be nullified: nullify(ja_copy);
!----WARNING!! 7[Priority]: The pointer listrank_copy shoud be nullified: nullify(listrank_copy);
!----WARNING!! 12[Priority]: The pointer kowner shoud be nullified: nullify(kowner);
!----WARNING!! 2[Priority]: The pointer aa_wsmp shoud be nullified: nullify(aa_wsmp);
!----WARNING!! 3[Priority]: The pointer rhs_wsmp shoud be nullified: nullify(rhs_wsmp);
!----WARNING!! 4[Priority]: The pointer unk_wsmp shoud be nullified: nullify(unk_wsmp);
!----WARNING!! 5[Priority]: The pointer ia_copy shoud be nullified: nullify(ia_copy);
!----WARNING!! 6[Priority]: The pointer ja_copy shoud be nullified: nullify(ja_copy);
!----WARNING!! 7[Priority]: The pointer listrank_copy shoud be nullified: nullify(listrank_copy);  



! a pensar -- en el RHS hago algo ?? debería traer contrib de vecinos?? NO ya esta confirme con JC

! aclarar bien cuando recibo un pointer si le tengo que poner que es pointer , o target, o nada
! hable con JC - si le voya  dar el tamaño y ya viene alocatado el no le pondría pointer.
! Vi que es una opción valida y lo haré así. En alya2maphys te pone target pero ahí es necesario
! porque lo apunta. 
!
!>    Si haces una búsqueda de bound_matrix, lo encontrarás en varios sitios. Lo tuyo
!podría ser bound_matrix_halo

!1. mod_parall: nulificación en varios sitios   --- ENTIENDO QUE ESTO SE LLAMA DE ALGUN LADO YO NO necesito volver a llamar
!2. mod_par_additional_arrays: se construyen los arrays de intercambio: ja, nzdom,nzdom_ii para cada neighbor   
!--- ESTO SE LLAMA DESDE  Domain.f90 L180   por ahoar le meto con #ifdef PARAL_AGMG - decidir luego con guillaume si lo dejo para todo el mundo o pongo un flag
!La diferencia es que para ti no se recibe lo mismo que se manda

!3. mod_communications: es donde hago el intercambio. Lo uso para el RAS desde solope.f90
!(hh PAR_INTERFACE_MATRIX_EXCHANGE_RP  l9983 mod_communications)
!./kernel/solite/solope.f90:334:        call PAR_INTERFACE_MATRIX_EXCHANGE(nbvar,solve_sol(1) % direct_solver_RAS % aa,PAR_COMM_MY_CODE_ARRAY(1))

! alguna snotas viejas
!
! 1) hacer ejemplo 3 elem alineados 3 subdom y ver si funciona
!
!
! 2) sobre si tener los nodos bien ordenados o sea que npoi2=npoi1+1   -- o sea intrior, bound_own, bound_other   y que no haya intrior,bound_other, bound_own, bound_other
!  par_array L180   lo hace el master  --- Jodido para meterme.
!
!
! Nuevos temas
!
! 1) Mirando el ejemplo de la pag 19 para TASK 1   me di cuenta que no necesito mover los componentes de la matriz (ver ja=[7 1 5
!  A mi me va a quedar 1 7 5 ; no está mal pero no es necesario. Solo es una cuestion de numeracion.
!
! 2) me doy cuenta que tener las unknowns 1x,2x,...npx,1y,2y,...npy ; no va a ser valido  si quiero tener own first!!!!   tenia anotado que esto era
!prefered by Notay pero no estoy seguro!!!
! voy a tener que cambiar toto lo de calc_lexplode  y ponerlo 1x,1y,2x,2y...  con lo cual mucho de lo que hecho es al pedo.
!
! Entre esto y lo anterior parece que mucho de lo que he hecho es al pedo.
!
! To commit
! 1) mod_alya2agmg: Changed to integer(4)   ia_copy,ja_copy,listrank_copy
!

!>   
!------------------------------------------------------------------------
!> @addtogroup Algebraic_Solver
!> @{
!> @name    alya2agmg
!> @file    mod_alya2agmg.f90
!> @author  Herbert Owen
!> @date    12/04/2016
!> @brief   Bridge to AGMG
!> @details Module to call AGMG iterative solver
!------------------------------------------------------------------------
module mod_alya2agmg

  use def_kintyp,         only : ip,rp,lg   !To revise what is actually needed
  use def_kintyp_solvers,         only : soltyp
  use mod_memory,         only : memory_alloca
  use mod_memory,         only : memory_deallo ! ver que este deallocando todo lo que no necesito
  use def_solver,         only : memit
  use def_kintyp_comm,    only : comm_data_par 
  use def_domain,         only : mesh_type
  use def_master,         only : kfl_paral
  use def_master,         only : lninv_loc
  use mod_communications_tools, only : PAR_COMM_TO_INT
  use mod_communications_global, only : PAR_MAX
  implicit none
  private

  public :: alya2agmg_initialization
  public :: alya2agmg_solution

  public :: mat2mathal          ! for  alya2wsmp_solution
  public :: calc_lmat2mathal    ! for  alya2wsmp_initialization
  !  public :: calc_lown_first     ! for  alya2wsmp_initialization

  public :: agmg_pt_nd             ! for  alya2psblas
  public :: mathal2matfinal        ! for  alya2psblas
  public :: calc_lexplode_ndofn    ! for  alya2psblas
  public :: max_ndofn_agmg_sol    ! for  alya2psblas

  type agmg_pt
     integer(ip)          :: npoin    ! added to avoid needing to pass meshe to alya2agmg_solution
     integer(ip)          :: npoin_own
     integer(ip)          :: npoin_2  ! added to avoid needing to pass meshe to alya2agmg_solution
     integer(ip)          :: npoi1    ! added to avoid needing to pass meshe to alya2agmg_solution
     integer(ip)          :: npoi2    ! added to avoid needing to pass meshe to alya2agmg_solution
     integer(ip)          :: npoi3    ! added to avoid needing to pass meshe to alya2agmg_solution
     integer(ip)          :: nz_own
     integer(ip)          :: nz_2
     integer(ip), pointer :: lwhalo_inv(:)
     integer(ip), pointer :: ja_yvan(:)
     integer(ip), allocatable :: listrank(:)   ! will use a normal allocate because I need to give it a lower and upper bound

  end type agmg_pt

  type agmg_pt_nd
     integer(ip), pointer :: ia_ndofn(:)
     integer(ip), pointer :: ja_ndofn(:)
     integer(ip), pointer :: lndofn_id(:)
     integer(ip), pointer :: lndofn_jd(:)
     integer(ip), pointer :: lndofn_iz(:)
     integer(ip), pointer :: lfinal_rhs(:)
     integer(ip), pointer :: lfinal_rhs_inv(:)
     integer(ip), allocatable :: listrank_ndofn(:)   ! will use a normal allocate because I need to give it a lower and upper bound

  end type agmg_pt_nd

  integer(ip),parameter       :: max_ndofn_agmg_sol=10
  integer(ip),parameter       :: kfl_debug=0

  type(agmg_pt)               :: pt
  type(agmg_pt_nd), pointer   :: ptnd(:)

contains 

  !-----------------------------------------------------------------------
  !
  !> @author  Herbert Owen
  !> @date    12/04/1016
  !> @brief   Initialize AGMG 
  !> @details Initialize AGMG for a given problem
  !>          This part I will think it all as if it were a scalar problem the fact that it can be
  !>          a vector problem such as the velocity I belive I can treat it later.
  !>          The ia and ja I need here are from the graph with halos only when I receive the matrix and wanto to build las
  !
  !-----------------------------------------------------------------------

  subroutine alya2agmg_initialization(meshe,solve)        ! (meshe,COMM,solve)

    use def_master, only       :  INOTMASTER
    use def_master, only       :  ISEQUEN


    ! en maphys eran de intent inout y ademas estab COMM de intent in
    type(mesh_type),     intent(in) :: meshe            !< Mesh type     ! en maphys era inout
    type(soltyp),        intent(in) :: solve            !< Solver type   ! en maphys era inout

    integer(ip)                :: i,nz2_aux,nd_aux
    integer(ip), save          :: ipass = 0

    if ( ISEQUEN ) return

    if( INOTMASTER ) then
       !
       ! mod_par_additional_arrays: se construyen los arrays de intercambio: ja, nzdom,nzdom_ii para cada neighbor   
       ! ESTO SE LLAMA DESDE  Domain.f90 L180   por ahora le meto con #ifdef PARAL_AGMG
       ! buscar bound_mat_
       ! tambien hay unos nullify en mod_parall   -- estos los dejo que se hagan siempre
       !
       nd_aux = solve % ndofn 

       if( ipass == 0 ) then
          ipass=1
          nullify(pt % lwhalo_inv)
          nullify(pt % ja_yvan)

          nullify(ptnd)
          allocate(ptnd(max_ndofn_agmg_sol))  
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

          call calc_lmat2mathal(meshe % npoin,meshe % nzdom ,meshe % npoin_2, & 
               nz2_aux,meshe % r_dom,meshe % c_dom, &
               meshe % r_dom_2,meshe % c_dom_2,pt % lwhalo_inv)
          if(kfl_debug==1) then
             if ( meshe%npoi2 - meshe%npoi1 /= 1 ) write(kfl_paral+500,*)'meshe%npoi1,meshe%npoi2',meshe%npoi1,meshe%npoi2
             write(*,*)'meshe%npoi1',meshe%npoi1
             write(*,*)'meshe%npoi2',meshe%npoi2
             write(*,*)'meshe%npoi3',meshe%npoi3
             write(*,*)'meshe%npoin',meshe%npoin
          end if


          pt % nz_own = meshe % r_dom_2(meshe % npoi3 + 1) - 1  ! Ahora guillauem hace que npoi2=npoi1+1  y npoi3 es el npoin_own
          pt % npoin_own = meshe % npoi3

          call memory_alloca(memit,'PT % JA_YVAN','alya2agmg_initialization',pt % ja_yvan,pt % nz_own)
          allocate(pt % listrank(pt % npoin_own+1:meshe % npoin_2)) ! it is actually a bit more than what I actually need
          !
          ! obtain xxx_yvan so that the numbering for non-own nodes in ja is as needed by agmg
          !
          call calc_lmatyvan(pt % npoin_own, pt % nz_own, meshe % npoi1, meshe % npoin_2, meshe % r_dom_2, meshe %c_dom_2, &   
               pt % ja_yvan, pt % listrank)

          if(1==2) then
             write(kfl_paral+500,*)'pt % ja_yvan(1:15)',pt % ja_yvan(1:15)
          end if
       end if

       !
       ! Obtain things that depend on ndofn
       ! 
       if(.not.associated(ptnd(nd_aux) % lndofn_id)) then

          call memory_alloca(memit,'PTND(SOLVE%NDOF) % IA_NDOFN','alya2agmg_initialization',ptnd(nd_aux) % ia_ndofn, &
               pt % npoin_own * nd_aux + 1)
          call memory_alloca(memit,'PTND(SOLVE%NDOF) % JA_NDOFN','alya2agmg_initialization',ptnd(nd_aux) % ja_ndofn, &
               pt % nz_own * nd_aux**2_ip)
          call memory_alloca(memit,'PTND(SOLVE%NDOF) % LNDOFN_ID','alya2agmg_initialization',ptnd(nd_aux) % lndofn_id, &
               pt % nz_own * nd_aux**2_ip)
          call memory_alloca(memit,'PTND(SOLVE%NDOF) % LNDOFN_JD','alya2agmg_initialization',ptnd(nd_aux) % lndofn_jd, &
               pt % nz_own * nd_aux**2_ip)
          call memory_alloca(memit,'PTND(SOLVE%NDOF) % LNDOFN_IZ','alya2agmg_initialization',ptnd(nd_aux) % lndofn_iz, &
               pt % nz_own * nd_aux**2_ip)

          call calc_lexplode_ndofn(nd_aux, pt % npoin_own, pt % nz_own, meshe % r_dom_2, pt % ja_yvan, &
               ptnd(nd_aux) % ia_ndofn,   ptnd(nd_aux) % ja_ndofn, ptnd(nd_aux) % lndofn_id, &
               ptnd(nd_aux) % lndofn_jd,  ptnd(nd_aux) % lndofn_iz)

          call memory_alloca(memit,'PTND(SOLVE%NDOF) % LFINAL_RHS','alya2agmg_initialization',ptnd(nd_aux) % lfinal_rhs, &
               pt % npoin_own * nd_aux)
          call memory_alloca(memit,'PTND(SOLVE%NDOF) % LFINAL_RHS_INV','alya2agmg_initialization',ptnd(nd_aux) % lfinal_rhs_inv, &
               meshe % npoin_2 * nd_aux)

          call calc_lfinal_rhs(nd_aux, pt % npoin_own, meshe % npoin_2, &
               ptnd(nd_aux) % lfinal_rhs, ptnd(nd_aux) % lfinal_rhs_inv)

          ! it is actually a bit more than what I actually need
          allocate(ptnd(nd_aux) % listrank_ndofn(pt % npoin_own * nd_aux+1 : meshe % npoin_2 * nd_aux)) 

          call calc_listrank_ndofn( pt % npoin_own, meshe % npoin_2, nd_aux , pt % listrank, ptnd(nd_aux) % listrank_ndofn)

          ! deallocate all except those that will be used in alya2agmg_solution
          ! ptnd() % lndofn_id,ptnd() % lndofn_jd,ptnd() % lndofn_iz, ptnd() % lfinal_rhs_inv, ptnd() % lfinal_rhs_inv
          ! ptnd() % ia_ndofn, ptnd() % ja_ndofn
          !          call memory_deallo(memit,'PTND(ND_AUX) % LNDOFN_IZ','alya2agmg_initialization',ptnd(nd_aux) % lndofn_iz)  ! now it used in place of l_last3
       end if

    end if   ! inotmaster

  end subroutine alya2agmg_initialization


  !-----------------------------------------------------------------------
  !
  !> @author  Herbert Owen
  !> @date    12/04/1016
  !> @brief   Solve Ax=b
  !> @details Solve algebraic system using AGMG.
  !>          
  !
  !-----------------------------------------------------------------------
  subroutine alya2agmg_solution(solve,amatr,rhsid,unkno)
    use mod_communications, only : PAR_INTERFACE_MATRIX_W_HALOS_EXCHANGE
    use mod_parall,         only : PAR_COMM_MY_CODE_ARRAY
    use mod_parall,         only : PAR_COMM_MY_CODE_WM
   
    use def_master,         only : INOTMASTER,ISEQUEN,ISLAVE

    use def_domain,     only :  meshe  ! for residual with alya matrix & also for ia & ja in serial version
    use def_kermod,     only :  ndivi

    use def_solver,     only :  SOL_SOLVER_AGMG, SOL_GAMGX

    type(soltyp),                intent(inout) :: solve                !< Solver type

    !  OJO inivar  L119  :  solve_sol(1) % nnz  =  meshe % nzdom
    real(rp),                    intent(in)    :: amatr(solve % ndofn,solve % ndofn,solve % nnz) !< Matrix
    real(rp),                    intent(inout) :: rhsid(solve % ndofn*solve % nequa) !< RHS
    real(rp),                    intent(inout) :: unkno(solve % ndofn*solve % nequa) !< Solution

    real(rp),pointer              :: aa_hal(:,:,:)
    real(rp),allocatable          :: aa_agmg(:)   ! For the moment I will leave it without pointer
    real(rp),allocatable          :: rhs_agmg(:)  ! For the moment I will leave it without pointer
    real(rp),allocatable          :: unk_agmg(:)  ! For the moment I will leave it without pointer

    !
    ! Changed to kind=4 ia_copy,ja_copy,listrank_copy,ifirstlistrank_aux,niter_aux,iprint,nz_aux, & created ndtnp_aux,nkryl_aux
    ! the conversions are added for all of them
    !
    integer(4),allocatable       :: ia_copy(:)  ! For the moment I will leave it without pointer
    integer(4),allocatable       :: ja_copy(:)  ! For the moment I will leave it without pointer
    integer(4)                   :: niter_aux,iprint,nz_aux
#ifdef AMGX
    integer(4),allocatable       :: listrank_copy(:)  ! For the moment I will leave it without pointer
    integer(4)                   :: ifirstlistrank_aux
    integer(ip)                  :: lowerb_lr(1),upperb_lr(1)
#endif



    integer(ip)  :: nd_aux,i
    integer(4)   :: PAR_COMM_TO_USE
    real(rp)     :: auxtol

    real(rp),allocatable    :: rr(:)
    real(rp)                :: resif,rhsnorm,resii
    integer(ip)             :: nbrows,korder ! for sequential
    integer(ip)             :: j,idof,irow,neqn
    integer(ip)             :: iblock(4)


#ifdef PARAL_AGMG
    integer(4)              :: ndtnp_aux,nkryl_aux
    integer(ip),save        :: kfl_already_set = 0
#endif

    korder = 1 !1x,1y,2x,2y see matreorder - for sequential
    ! Only the first processor prints output
    ! if 'OUTPUT: CONVERGENCE' option is enabled
    if ((solve % kfl_cvgso == 1) .and. (kfl_paral==1)) then
       iprint = 6  ! output unit for agmg
    else
       iprint = -1 ! 
    end if
    niter_aux = int(solve % miter,4) ! on output it will have the number of iteration done - I can put it to the sol file
    ! for tol the same thing does not happen it does not give you the residual
    ! perhaps it thinks we should always converge and thus it does not make sense to output it
    ! well it is in the agmg output file
    !
    ! rhsnorm
    !
    nd_aux = solve % ndofn  ! guess the master must have some correct value
    call norm2x(nd_aux,rhsid,rhsnorm)
    !
    ! initial residual & auxtol
    !
    auxtol = solve % solco

    allocate(rr(nd_aux*max(1_ip,meshe(ndivi) % npoin)))
    rr = 0.0_rp
    ! rr = amatr * x
    call bcsrax( 1_ip, meshe(ndivi) % npoin, nd_aux, amatr, meshe(ndivi) % c_dom, meshe(ndivi) % r_dom, unkno, rr )
    if( INOTMASTER ) rr = rr - rhsid
    call norm2x(nd_aux,rr,resii)
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
       nullify(aa_hal)
       call memory_alloca(memit,'AA_HAL','alya2agmg_solution',aa_hal,nd_aux,nd_aux,pt % nz_2)
       aa_hal=0.0_rp
       allocate(aa_agmg(nd_aux * nd_aux * pt % nz_2))
       !
       ! Pass the matrix aa to a matrix aa_hal  - In alya2agmg_initialization I create a pointer that makes this automatic ej lmat2mathal 
       !
       call mat2mathal(solve % nnz,pt % nz_2,nd_aux,pt % lwhalo_inv,amatr,aa_hal)

       !
       ! Add the contibution lo the matrix with halos aa_hal from the neighbors
       ! there is a call to call PAR_INTERFACE_MATRIX_EXCHANGE in solope tu use as example
       !
       call PAR_INTERFACE_MATRIX_W_HALOS_EXCHANGE(nd_aux,amatr,aa_hal,PAR_COMM_MY_CODE_ARRAY(1))
       !    call PAR_INTERFACE_MATRIX_EXCHANGE(nbvar,solve_sol(1) % direct_solver_RAS % aa,PAR_COMM_MY_CODE_ARRAY(1))  ! de solope guillaume

       !
       ! Now convert to the final format also reorder rhs and initial guess  
       !
       call mathal2matfinal  (pt % nz_own,pt % nz_2,nd_aux,ptnd(nd_aux) % lndofn_id,ptnd(nd_aux) % lndofn_jd, &
            ptnd(nd_aux) % lndofn_iz,aa_hal,aa_agmg)       
       !
       ! Convert the initial guess and rhs to AGMG ordering
       !
       allocate(unk_agmg(nd_aux * pt % npoin_own))
       allocate(rhs_agmg(nd_aux * pt % npoin_own))

       call unkalya2agmg(pt % npoin_own, pt % npoin_2, nd_aux, ptnd(nd_aux) % lfinal_rhs, rhsid, rhs_agmg)
       call unkalya2agmg(pt % npoin_own, pt % npoin_2, nd_aux, ptnd(nd_aux) % lfinal_rhs, unkno, unk_agmg)    
       !
       ! AGMG Solve
       !
       ! I will directly control cg or GCR with solve % nkryd 

       !       call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)   ! JC told me to use it in this way
       PAR_COMM_TO_USE = PAR_COMM_TO_INT(PAR_COMM_MY_CODE_WM)   ! communicator without master

#ifdef AMGX

       allocate( ia_copy((pt % npoin_own * nd_aux) + 1))
       ia_copy( 1: (pt % npoin_own * nd_aux) + 1) = int(ptnd(nd_aux) % ia_ndofn ( 1: (pt % npoin_own * nd_aux) + 1),4)

       nz_aux = ia_copy ( ( pt % npoin_own * nd_aux ) + 1 )  - 1
       allocate( ja_copy( nz_aux )  )
       ja_copy ( 1:nz_aux ) = int(ptnd(nd_aux) % ja_ndofn ( 1:nz_aux ),4)

       ! listrank and ifirstlistrank may be modified on output - see dagmg_par
       lowerb_lr = lbound(ptnd(nd_aux) % listrank_ndofn,kind=ip)
       upperb_lr = ubound(ptnd(nd_aux) % listrank_ndofn,kind=ip)
       allocate(listrank_copy(lowerb_lr(1) : upperb_lr(1)))
       listrank_copy(lowerb_lr(1) : upperb_lr(1)) = int(ptnd(nd_aux) % listrank_ndofn(lowerb_lr(1) : upperb_lr(1)),4)

       ifirstlistrank_aux = int(pt % npoin_own*nd_aux+1,4)

       call amgx_wrap(solve % nequa, solve % ndofn, solve % miter,amatr,rhsid,unkno, &
            ia_copy, ja_copy, listrank_copy, lowerb_lr(1), upperb_lr(1) , pt % npoin_own, &
            pt % npoin_2, nd_aux, PAR_COMM_MY_CODE_ARRAY(1))


       deallocate(ia_copy)
       deallocate(ja_copy)
       deallocate(listrank_copy)

#endif


#ifdef PARAL_AGMG
       !    call dagmgpar(N,a,ja,ia,f,x,ijob,iprint,nrest,iter,tol & MPI_COMM , listrank , ifirstlistrank )
       !              A,IA,JA are "output" parameters because on exit the
       !              entries of each row may occur in a different order (The
       !              matrix is mathematically the same, but stored in
       !              different way).
       ! ver de pasarle copia no sea que los cambie y luego yo me jodo
       !
       ! Create copy ia and ja. With aa_agmg there is no problem because it creates it every time it enters here
       ! this I believe is not the most efficient thing to do but I am not sure how to do it otherwise
       !
       allocate( ia_copy((pt % npoin_own * nd_aux) + 1))
       ia_copy( 1: (pt % npoin_own * nd_aux) + 1) = int(ptnd(nd_aux) % ia_ndofn ( 1: (pt % npoin_own * nd_aux) + 1),4)

       nz_aux = ia_copy ( ( pt % npoin_own * nd_aux ) + 1 )  - 1
       allocate( ja_copy( nz_aux )  )
       ja_copy ( 1:nz_aux ) = int(ptnd(nd_aux) % ja_ndofn ( 1:nz_aux ),4)

       ! listrank and ifirstlistrank may be modified on output - see dagmg_par
       lowerb_lr = lbound(ptnd(nd_aux) % listrank_ndofn,kind=ip)
       upperb_lr = ubound(ptnd(nd_aux) % listrank_ndofn,kind=ip)
       allocate(listrank_copy(lowerb_lr(1) : upperb_lr(1)))
       listrank_copy(lowerb_lr(1) : upperb_lr(1)) = int(ptnd(nd_aux) % listrank_ndofn(lowerb_lr(1) : upperb_lr(1)),4)

       if(1==2) then
          write(kfl_paral+500,*)'ptnd(nd_aux) % ja_ndofn(1:13)',ptnd(nd_aux) % ja_ndofn(1:13)
          write(kfl_paral+500,*)'ptnd(nd_aux) % listrank_ndofn',ptnd(nd_aux) % listrank_ndofn
       end if

       !
       ! test Ricard
       !
       !       if(kfl_debug==1) call test_ricard(solve,amatr,aa_agmg,ia_copy,nd_aux, pt % nz_2, pt % npoin_own)
       !    
       ifirstlistrank_aux = int(pt % npoin_own*nd_aux+1,4)

       if(nd_aux==1.and.kfl_debug==1) then
          write(6000+kfl_paral,'(300(e10.3,1x))')aa_agmg
          write(6000+kfl_paral,'(300(i5,1x))')ja_copy
          write(6000+kfl_paral,'(300(i5,1x))')ia_copy
          write(6000+kfl_paral,'(300(e10.3,1x))')rhs_agmg
          write(6000+kfl_paral,'(300(e10.3,1x))')unk_agmg
          write(6000+kfl_paral,'(300(i5,1x))') meshe(ndivi) % npoi1
          do i=1,meshe(ndivi) % npoi1
             write(6000+kfl_paral,'(300(i5,1x))')i,lninv_loc(i)
          end do
          write(6000+kfl_paral,*)'OWNBOUN', meshe(ndivi) % npoi2, meshe(ndivi) % npoi3
          do i=meshe(ndivi) % npoi2, meshe(ndivi) % npoi3
             write(6000+kfl_paral,'(300(i5,1x))')i,lninv_loc(i)
          end do
       end if

       if( solve % kfl_clean_precond == 0  )   then          ! NEVER CHANGE
          if( solve % kfl_algso == SOL_SOLVER_AGMG ) then
             if (kfl_already_set == 0 ) then
                ndtnp_aux = int(pt % npoin_own*nd_aux,4)
                nkryl_aux = int(solve % nkryd,4)
                call dagmgpar(ndtnp_aux, aa_agmg, ja_copy, ia_copy, rhs_agmg, unk_agmg,  &
                     1, iprint, nkryl_aux, niter_aux, auxtol , PAR_COMM_TO_USE , &     ! 1: performs setup only
                     listrank_copy , ifirstlistrank_aux )
                kfl_already_set = 1
             end if
             ndtnp_aux = int(pt % npoin_own*nd_aux,4)
             nkryl_aux = int(solve % nkryd,4)
             call dagmgpar(ndtnp_aux, aa_agmg, ja_copy, ia_copy, rhs_agmg, unk_agmg,  &
                  12, iprint, nkryl_aux, niter_aux, auxtol , PAR_COMM_TO_USE , &        ! 12:  solves , initial guess in x(1:n)
                  listrank_copy , ifirstlistrank_aux )

          end if
       else
          ndtnp_aux = int(pt % npoin_own*nd_aux,4)
          nkryl_aux = int(solve % nkryd,4)
          if( solve % kfl_algso == SOL_SOLVER_AGMG ) &
               call dagmgpar(ndtnp_aux, aa_agmg, ja_copy, ia_copy, rhs_agmg, unk_agmg,  &
               10, iprint, nkryl_aux, niter_aux, auxtol , PAR_COMM_TO_USE , &         ! 10: performs setup + solve + memory release, initial guess in x(1:n)
               listrank_copy , ifirstlistrank_aux )
       end if

       if(nd_aux==1.and.kfl_debug==1) then
          write(6000+kfl_paral,'(300(e10.3,1x))')unk_agmg
       end if

       flush(iprint)

       solve % iters = int(niter_aux,ip)

       deallocate(ia_copy)
       deallocate(ja_copy)
       deallocate(listrank_copy)

#endif

       !
       ! put solution to alyas ordering
       !
       call unkagmg2alya(pt % npoin_own, pt % npoin, pt % npoi1, pt % npoi2, pt % npoi3, nd_aux,  &
            ptnd(nd_aux) % lfinal_rhs_inv, unk_agmg, unkno)

       call memory_deallo(memit,'AA_HAL'     ,'alya2agmg_solution',aa_hal)
       deallocate(aa_agmg)
       deallocate(rhs_agmg)
       deallocate(unk_agmg)     

    else if ( ISEQUEN ) then

       nbrows = meshe(ndivi) % npoin
       neqn = nbrows*nd_aux
       nz_aux  = meshe(ndivi) % r_dom(nbrows+1)-1 !TODO: Possible change of value in conversion from INTEGER(8) to INTEGER(4)

       allocate(aa_agmg(nd_aux*nd_aux*nz_aux))
       allocate(ja_copy(nd_aux*nd_aux*nz_aux))
       allocate(ia_copy(nd_aux*nbrows+1))
       allocate(rhs_agmg(nd_aux*nbrows))
       allocate(unk_agmg(nd_aux*nbrows))

       call matreorder(rhsid, unkno, amatr, meshe(ndivi) % r_dom, meshe(ndivi) % c_dom, nd_aux, nbrows, korder, aa_agmg, &
            ja_copy, ia_copy, rhs_agmg, unk_agmg)       !
       !
       ! SOLVE
       !
       !
       ! Missing see how to add the CPU Times for agmg directly into solve_sol(1) % cputi(??)
       !
       if (nd_aux >4 ) call runend('nsi_agmgsol: nd_aux>4')
       iblock(1) = 1
       iblock(2) = 1+nbrows
       iblock(3) = 1+2*nbrows
       iblock(4) = 1+3*nbrows
       !
       !#ifdef solve_w_agmg now I will use PARAL_AGMG - actually later I should call it directly AGMG
#ifdef PARAL_AGMG       
       call dagmg(nd_aux*nbrows, aa_agmg, ja_copy, ia_copy, rhs_agmg, unk_agmg, 10, iprint, solve % nkryd, &
            niter_aux, auxtol) !non_block
       ! call dagmg(nd_aux*nbrows,aa_agmg,ja_copy, ia_copy, rhs_agmg, unk_agmg, 10, iprint, solve % nkryd, &
       !    niter_aux, auxtol, nd_aux, iblock) ! block
#endif  
       !
       ! now pass unk_agmg to unkno
       !
       if (korder==1) then
          unkno(1:nd_aux*nbrows) = unk_agmg(1:nd_aux*nbrows)
       else
          do idof = 1,nd_aux
             do irow = 1,nbrows
                i = idof+(irow-1)*nd_aux
                j = irow+(idof-1)*nbrows
                unkno(i) = unk_agmg(j)
             end do
          end do
       end if

       solve % iters = int(niter_aux,ip)

       deallocate(aa_agmg)
       deallocate(ja_copy)
       deallocate(ia_copy)
       deallocate(rhs_agmg)
       deallocate(unk_agmg)

    end if

    !
    ! Calculate final residual with alya matrix - This serves as a check
    !
    allocate(rr(max(1_ip,nd_aux*meshe(ndivi) % npoin)))
    rr = 0.0_rp
    ! rr = amatr * x
    call bcsrax( 1_ip, meshe(ndivi) % npoin, nd_aux, amatr, meshe(ndivi) % c_dom, meshe(ndivi) % r_dom, unkno, rr )
    if( INOTMASTER ) rr = rr - rhsid
    call norm2x(nd_aux,rr,resif)
    deallocate(rr)

    solve % xorth = rhsnorm   ! beware I will output it in the colum that corresponds xorth that I do not have it for AGMG
    solve % resin = resii     ! I do not have a preconditioned residual I will output this
    solve % resi2 = resii/rhsnorm    ! this is the one that typically appears in Alya - Column 4 - Non preconditioned initial residual 
    solve % resfi = resif
    solve % resf2 = resif/rhsnorm
    call PAR_MAX(solve % iters)


  end subroutine alya2agmg_solution


  subroutine calc_lmat2mathal(npoin,nz,npoin_2,nz2,ia,ja,ia2,ja2,lwhalo_inv)
    !
    ! Obtains lwhalo such that amatr_w_halo(iz) = amatr(lwhalo(iz))
    ! Later I saw that what I really needed was  lwhalo_inv such that amatr_w_halo(lwhalo_inv(iz)) = amatr(iz)
    ! ia2 and ja2 are already known - they are intent in
    ! the unknowns do not change order - then rhs does not change
    !
    integer(ip),intent(in)  :: npoin,nz,npoin_2,nz2
    integer(ip),intent(in)  :: ia(npoin+1),ja(nz)
    integer(ip),intent(in)  :: ia2(npoin_2+1),ja2(nz2)
    !    integer(ip),intent(out) :: lwhalo(nz2)
    integer(ip),intent(out) :: lwhalo_inv(nz)

    integer(ip)  :: ipoin,iz,jj,jp,kp

    do ipoin=1,npoin
       !
       ! 
       !       
       do jj = ia(ipoin), ia(ipoin+1)-1
          iz = ia2(ipoin)    ! I recalculate it even if ipoin has not changed because I alter it below
          jp = ja(jj)
          do while( iz <= ia2(ipoin+1)-1 )
             kp = ja2(iz)
             if( kp == jp ) then
                !               lwhalo(iz) = jj
                lwhalo_inv(jj) = iz
                iz = ia2(ipoin+1)
             end if
             iz = iz + 1
          end do
          if( iz /= ia2(ipoin+1)+1 ) then
             write(kfl_paral+500,*)'COULD NOT FIND NODE: ipoin,iz',ipoin,iz
             call runend('calc_lmat2mathal: COULD NOT FIND NODE')
          end if
       end do
    end do

  end subroutine calc_lmat2mathal

  subroutine mat2mathal(nz,nz2,ndofn,lwhalo_inv,amatr,amatr_w_halo)

    !
    ! amatr_w_halo(lwhalo_inv(iz)) = amatr(iz)
    ! 

    integer(ip),intent(in)          :: nz,nz2,ndofn
    integer(ip),intent(in)          :: lwhalo_inv(nz)
    real(rp),intent(in), target     :: amatr(ndofn,ndofn,nz)         !I put target following what is done in alya2maphys
    real(rp),intent(out), target    :: amatr_w_halo(ndofn,ndofn,nz2)

    integer(ip)  :: iz,jz

    do iz=1,nz
       jz=lwhalo_inv(iz)
       if (jz<1) then
          write(kfl_paral+500,*)'iz,jz',iz,jz
          call runend('mat2mathal:jz<1')
       end if
       amatr_w_halo(:,:,jz) = amatr(:,:,iz)
    end do

  end subroutine mat2mathal




  subroutine calc_lmatyvan(npoin_own,nz_own,npoi1,npoin_2,ia2,ja2,ja_yvan,listrank)

    ! This subroutine fulfills the second part of the following requirement 
    ! Consistency of local orderings. Besides the condition that they are larger than N, indexes of nonlocal variable may be chosen arbitrarily,
    ! !!!!!! providing that their ordering is consistent with the local ordering on their “home” task!!!!!
    !
    ! Obtains ja_yvan & listrank ;  l_mat_yvan & ia_yvan should later be eliminated
    ! Just obtain a new ja in the part corresponding to not-own nodes
    !
    ! it seems I do not need to do anything to the matrix - just obtain a new ja in the part corresponding to not-own nodes
    ! In order to having working quickly I will leave l_mat_yvan so that it does nothing - later I will no longer calculate it  ! I have eliminated it now!!!
    ! and eliminate it from l_last3
    !
    use mod_parall,    only : PAR_COMM_MY_CODE_ARRAY
    use mod_maths,     only : maths_heap_sort

    integer(ip),intent(in)  :: npoin_own
    integer(ip),intent(in)  :: nz_own
    integer(ip),intent(in)  :: npoi1
    integer(ip),intent(in)  :: npoin_2
    integer(ip),intent(in)  :: ia2(npoin_own+1)
    integer(ip),intent(in)  :: ja2(nz_own)
    integer(ip),intent(out) :: ja_yvan(nz_own)  ! beware I have obtained it only for the row of own nodes
    integer(ip),intent(out) :: listrank(npoin_own+1:npoin_2)


    integer(ip)                         :: nreorder
    integer(ip),allocatable             :: kowner(:)   ! i will use normal allocate not to make them pointers
    integer(ip),allocatable             :: num_in_owner(:)
    integer(ip),allocatable             :: index_table(:)
    integer(ip),allocatable             :: aux_copy(:)
    integer(ip),allocatable             :: index_table_2(:)
    integer(ip),allocatable             :: jauxi(:)
    integer(ip),allocatable             :: ordered_list_inv(:)

    integer(ip)  :: ipoin,jj,jp
    integer(ip)  :: jpoin,nzeros,kdomain,istar,iend,nordered,iauxi

    allocate(jauxi(npoin_own+1:npoin_2))
    allocate(kowner(npoin_own+1:npoin_2))
    allocate(num_in_owner(npoin_own+1:npoin_2))
    allocate(index_table(npoin_2-npoin_own))
    allocate(aux_copy(npoin_2-npoin_own))
    allocate(index_table_2(npoin_2-npoin_own))  ! this is actually much more than I need
    allocate(ordered_list_inv(npoin_2-npoin_own))

    listrank = -1_ip !initialization

    jauxi = 0      
    num_in_owner=0  
    kowner=0       
    index_table=0   
    index_table_2=0 
    aux_copy=0       
    ordered_list_inv=0 
    !
    ! par_node_number_in_owner that obtains PAR_COMM_MY_CODE_ARRAY(1) % node_number_in_owner has been called in domain
    !
    ! I have to order all non-own nodes conected to own nodes of this subdomain according to the numbering in their master subdomain.
    ! Moreover as two non-own nodes could have the same number in their master subdomain the right thing to do seems to order all
    ! non-own nodes corresponding to one master subdomain, then those corresponding to another subdomain and so on. 
    !
    !
    !
    ! feb 17: For all points jp conected to own nodes obtain their OWNER subdomain (kowner(jp))
    ! and the local numbering in the owner subdomain  (num_in_owner(jp)).
    ! I have to do it this way because in Alya halos are all nodes conected to a node in my boundary.
    ! While for AGMG only those connected to own nodes are of interest (which I belive is a more natural definition of halo).
    !
    nreorder = 0
    do ipoin=1,npoin_own
       do jj = ia2(ipoin), ia2(ipoin+1)-1
          jp = ja2(jj)
          if ( jp <= npoin_own ) then
             ja_yvan(jj) = ja2(jj) ! just leave the same value in own nodes
          else
             if(jauxi(jp)==0) then
                nreorder = nreorder + 1
                kowner(jp) = PAR_COMM_MY_CODE_ARRAY(1) % bound_owner_rank(jp-npoi1)
                num_in_owner(jp) = PAR_COMM_MY_CODE_ARRAY(1) % node_number_in_owner(jp-npoi1)
                jauxi(jp)=1
             end if
          end if
       end do
    end do

    if ( nreorder > 0 ) then          ! Feb 17: What happens if nreorder=0. Should we have an else? 
       do ipoin=1,npoin_2-npoin_own      ! Order by Owner
          index_table(ipoin) = ipoin
          aux_copy(ipoin) = kowner(ipoin+npoin_own)
       end do
       iauxi=npoin_2-npoin_own   ! because it is intent inout in maths_heap_sort
       call maths_heap_sort(2_ip,iauxi,aux_copy,' ',index_table)  ! 2 increasing
       !
       nzeros=0
       zeros: do ipoin=1,npoin_2-npoin_own
          if (kowner(index_table(ipoin)+npoin_own) == 0) then
             nzeros = nzeros + 1
          else
             exit zeros
          end if
       end do zeros

       ! first ordered according to the subdomain to which they belong and then
       ! and then according to the numbering in their owner
       ! understand well why there are zeros - guess it is alya understands halos
       ! different from agmg  !!!!!! esto tendre que arreglarlo

       kdomain = kowner(index_table(nzeros+1)+npoin_own)
       !       kount = 1  ! related to the fact that I use nzeros+2  2 lines later
       istar = nzeros+1   ! The position where a group of nodes of the same subdomain starts. 
       nordered=0
       do ipoin=nzeros+2,npoin_2-npoin_own    ! For ech owner order by number in owner  
          if (kowner(index_table(ipoin)+npoin_own) /= kdomain .or. ipoin == npoin_2-npoin_own) then
             ! new domain
             ! process the old one
             iend = ipoin-1
             if (ipoin == npoin_2-npoin_own) iend = ipoin
             do jpoin=1,iend-istar+1
                aux_copy(jpoin)=num_in_owner(index_table(jpoin+istar-1)+npoin_own)
                index_table_2(jpoin) = jpoin
             end do

             iauxi=iend-istar+1   ! because it is intent in out in maths_heap_sort
             call maths_heap_sort(2_ip,iauxi,aux_copy,' ',index_table_2)
             !
             ! Obtain ordered list
             !
             do jpoin=1,iend-istar+1
                listrank(jpoin+nordered+npoin_own) = kowner(npoin_own+index_table(index_table_2(jpoin)+istar-1))
                ordered_list_inv(index_table(index_table_2(jpoin)+istar-1)) = jpoin + nordered + npoin_own !ojo probando npoin_own
             end do
             ! now prepare to continue
             nordered = nordered + (iend-istar+1)
             if (ipoin==npoin_2-npoin_own) then   ! small check
                if (nordered+nzeros/=npoin_2-npoin_own) then
                   print*,'nordered,nzeros,npoin_2,npoin_own',nordered,nzeros,npoin_2,npoin_own
                   call runend('nordered+nzeros/=npoin_2-npoin_own')
                end if
             else
                istar = ipoin
                kdomain = kowner(index_table(ipoin)+npoin_own)
             end if
          end if
       end do
       !
       ! now I can complete ja
       !
       do ipoin=1,npoin_own

          do jj = ia2(ipoin), ia2(ipoin+1)-1
             jp = ja2(jj)
             if ( jp > npoin_own ) then
                ja_yvan(jj) = ordered_list_inv(jp-npoin_own) 
             end if
          end do
       end do
    end if

    deallocate(jauxi)
    deallocate(kowner)
    deallocate(num_in_owner)
    deallocate(index_table)
    deallocate(aux_copy)
    deallocate(index_table_2)
    deallocate(ordered_list_inv)    

  end subroutine calc_lmatyvan


  subroutine calc_lexplode_ndofn(ndofn,npoin_own,nz_own,ia_in,ja_in,ia_ndofn,ja_ndofn, &
       lndofn_id,lndofn_jd,lndofn_iz)
    !
    ! This subroutine creates the arrays to transform a matrix for a vector unknown as it is stored in Alya
    ! to a general matrix where each component of the vector is treated as a scalar unknown. This is need for solvers that have no specific
    ! format for vector unknowns such as AGMG.
    !
    ! I will use what I called korder==1 in matreorder. Initially I had used 2 but realized that it would not be compatible
    ! with the fact that own unknowns must go first.
    ! so I am now using 1x,1y,1x,2y.....nx,ny
    ! Somewhere I had witten that Yvan prefered   1x,2x,...npx,1y,2y,...npy ;  
    ! Re check this with Yvan my idea was that having them closer would favour grouping in agmg but actually there should be no reason
    ! The numbering sould not affect the way they are agregated. Re think. Ask Yvan if AGMG has been tested in non scalar problems. 
    !

    !
    ! Obtains lndofn_id,lndofn_jd,lndofn_iz such that  amatr_ndofn(iz) = amatr_yvan_num(lndof_id(iz),lndofn_jd(iz),lndofn_iz(iz))
    ! Rhs - see details just before its calculation.
    !
    integer(ip),intent(in)  :: ndofn,npoin_own,nz_own
    integer(ip),intent(in)  :: ia_in(npoin_own+1)
    integer(ip),intent(in)  :: ja_in(nz_own)

    integer(ip),intent(out) :: ia_ndofn(npoin_own*ndofn+1)
    integer(ip),intent(out) :: ja_ndofn(nz_own*ndofn*ndofn)

    integer(ip),intent(out) :: lndofn_id(nz_own*ndofn*ndofn)
    integer(ip),intent(out) :: lndofn_jd(nz_own*ndofn*ndofn)
    integer(ip),intent(out) :: lndofn_iz(nz_own*ndofn*ndofn)

    integer(ip) :: nz,ii_ndofn,iz_ndofn,kk1,kk2,ii,jj,iz

    ii_ndofn = 0
    iz_ndofn = 0
    !
    ! matrix - I will have to go back to what I had called korder=1 in matreorder
    !  ->  1x,1y,1x,2y......npx,npy   - Th other option is incompatible with the fact that own have to go first.
    !  strangelly when I was working in serial I thought the other option was the one preffered by yvan - check
    !

    do ii = 1,npoin_own
       do kk2 = 1,ndofn
          do iz = ia_in(ii),ia_in(ii+1)-1
             jj = ja_in(iz)
             do kk1 = 1,ndofn
                iz_ndofn = iz_ndofn + 1
                !amatr_ndofn(iz_ndofn) = amatr(kk1,kk2,iz)
                lndofn_iz(iz_ndofn) = iz
                lndofn_id(iz_ndofn) = kk1
                lndofn_jd(iz_ndofn) = kk2
                ja_ndofn(iz_ndofn)    = (jj-1)*ndofn + kk1
             end do
          end do
       end do
       nz = ia_in(ii+1)-ia_in(ii)
       do kk1 = 1,ndofn
          ii_ndofn = ii_ndofn + 1
          ia_ndofn(ii_ndofn) = nz*ndofn
       end do
    end do

    kk1        = ia_ndofn(1)
    ia_ndofn(1) = 1 
    do ii = 2,ndofn*npoin_own+1
       kk2         = ia_ndofn(ii)
       ia_ndofn(ii) = ia_ndofn(ii-1) + kk1
       kk1         = kk2
    end do

  end subroutine calc_lexplode_ndofn

  subroutine mathal2matfinal  (nz_own,nz2,ndofn,lndofn_id,lndofn_jd,lndofn_iz,amatr_w_halo,amatr_final)
    !
    ! amatr_ndofn(iz) =  amatr_w_halo(lndofn_id(iz),lndofn_jd(iz),lndofn_iz(iz))
    ! 
    integer(ip),intent(in)          :: nz_own,ndofn,nz2
    integer(ip),intent(in)          :: lndofn_id(nz_own*ndofn*ndofn)
    integer(ip),intent(in)          :: lndofn_jd(nz_own*ndofn*ndofn)
    integer(ip),intent(in)          :: lndofn_iz(nz_own*ndofn*ndofn)
    real(rp),intent(in), target     :: amatr_w_halo(ndofn,ndofn,nz2)         !I put target folling what is done in alya2maphys
    real(rp),intent(out)            :: amatr_final(ndofn*ndofn*nz_own)

    integer(ip)  :: iz,jz,id,jd

    do iz=1,nz_own*ndofn*ndofn
       id = lndofn_id(iz)
       jd = lndofn_jd(iz)
       jz = lndofn_iz(iz)
       if (id<1) call runend('mathal2matfinal: id<1')
       if (jd<1) call runend('mathal2matfinal: jd<1')
       if (jz<1) call runend('mathal2matfinal: jz<1')
       amatr_final(iz) =  amatr_w_halo(id,jd,jz)
    end do

  end subroutine mathal2matfinal

  subroutine calc_lfinal_rhs(ndofn,npoin_own,npoin_2, &
       lfinal_rhs,lfinal_rhs_inv)

    ! The explanation is the following:  
    ! o   rhs_exploded(i)  =  rhs_w_halo(j) 
    !                      =  rhs_w_halo(lfinal_rhs(i))
    !
    !
    integer(ip),intent(in)          :: npoin_own,npoin_2,ndofn  !,npoi1,npoi2,npoi3
    integer(ip),intent(out)         :: lfinal_rhs(npoin_own*ndofn)
    integer(ip),intent(out)         :: lfinal_rhs_inv(npoin_2*ndofn)
    !    
    integer(ip)          :: i,idofn,ipoin
    !
    ! Obtain lfinal_rhs
    !
    do idofn=1,ndofn
       do ipoin=1,npoin_own
          i = idofn + (ipoin-1) * ndofn
          lfinal_rhs(i) = idofn+(ipoin-1)*ndofn
       end do
    end do
    !
    ! Obtain lfinal_rhs_inv
    !
    do i=1,npoin_own*ndofn
       lfinal_rhs_inv(lfinal_rhs(i)) = i
    end do

  end subroutine calc_lfinal_rhs

  subroutine unkalya2agmg(npoin_own,npoin_2,ndofn,lfinal_rhs,rhs_w_halo,rhs_exploded)
    !
    ! puts unknown to AGMG's ordering
    ! some unification of names is still missing
    ! here I refer to the final agmg ordering perhaps at some point
    ! At some other points I acll this exploded
    ! I should unify names at some point
    ! 
    ! rhs_exploded(i) = rhs_w_halo(lfinal_rhs(i))
    ! or I could have written it for unkno instead of rhs
    !
    ! At some point I will have to define some unified naming
    !
    integer(ip),intent(in)          :: npoin_own,npoin_2,ndofn
    integer(ip),intent(in)          :: lfinal_rhs(npoin_own*ndofn)
    real(rp),intent(in)             :: rhs_w_halo(npoin_2*ndofn)
    real(rp),intent(out)            :: rhs_exploded(npoin_own*ndofn)

    integer(ip)          :: i

    do i=1,npoin_own*ndofn
       rhs_exploded(i) = rhs_w_halo(lfinal_rhs(i))
    end do

  end subroutine unkalya2agmg


  subroutine unkagmg2alya(npoin_own,npoin,npoi1,npoi2,npoi3,ndofn,lfinal_rhs_inv,unkno_exploded,unkno_w_halo)
    !
    ! puts unknown back to Alya's ordering
    !
    ! unkno_w_halo(i) = unkno_exploded(lfinal_rhs_inv(i))
    !
    use mod_communications, only : PAR_INTERFACE_NODE_EXCHANGE

    integer(ip),intent(in)          :: npoin_own,npoin,ndofn
    integer(ip),intent(in)          :: npoi1,npoi2,npoi3
    integer(ip),intent(in)          :: lfinal_rhs_inv(npoin*ndofn)
    real(rp),intent(in)             :: unkno_exploded(npoin_own*ndofn)
    real(rp),intent(out)            :: unkno_w_halo(npoin*ndofn)

    real(rp),pointer                :: unkno_aux(:,:)
    integer(ip)                     :: i,ipoin,idofn

    nullify(unkno_aux)
    call memory_alloca(memit,'UNKNO_AUX','unkagmg2alya',unkno_aux,ndofn,npoin)
    unkno_aux = 0.0_rp

    !    write(*,*)'unkno_exploded',unkno_exploded
    do i=1,npoi1*ndofn
       unkno_w_halo(i) = unkno_exploded(lfinal_rhs_inv(i))
       !       write(*,*)'i,lfinal_rhs_inv(i)',i,lfinal_rhs_inv(i)
       !       flush(6)
    end do
    do i=(npoi2-1)*ndofn+1, npoi3*ndofn
       unkno_w_halo(i) = unkno_exploded(lfinal_rhs_inv(i))
       !       write(*,*)'i,lfinal_rhs_inv(i)',i,lfinal_rhs_inv(i)
       !       flush(6)
    end do
    !    write(*,*)'unkno_w_halo',unkno_w_halo


    !
    ! Exchange - otherwise in the boundaries it did not come out well
    !   
    do i=1, npoi1*ndofn
       ipoin = (i+ndofn-1)/ndofn
       idofn = i-(ipoin-1)*ndofn
       unkno_aux(idofn,ipoin) = unkno_w_halo(i)
    end do
    do i = (npoi2-1)*ndofn+1, npoi3*ndofn !revisar
       ipoin = (i+ndofn-1)/ndofn
       idofn = i-(ipoin-1)*ndofn
       unkno_aux(idofn,ipoin) = unkno_w_halo(i)
    end do
    !    write(*,*)'unkno_aux',unkno_aux

    call PAR_INTERFACE_NODE_EXCHANGE(unkno_aux,'SUM'       ,'IN MY CODE')
    !    write(*,*)'unkno_aux',unkno_aux

    do i=1, npoin*ndofn
       ipoin = (i+ndofn-1)/ndofn
       idofn = i-(ipoin-1)*ndofn
       unkno_w_halo(i) = unkno_aux(idofn,ipoin)
    end do


    call memory_deallo(memit,'UNKNO_AUX','unkagmg2alya',unkno_aux)

  end subroutine unkagmg2alya


  subroutine calc_listrank_ndofn(npoin_own,npoin_2,ndofn,listrank,listrank_ndofn)

    ! 
    ! Obtains listrank_ndofn from listrank
    !
    integer(ip),intent(in)  :: npoin_own
    integer(ip),intent(in)  :: npoin_2
    integer(ip),intent(in)  :: ndofn 
    integer(ip),intent(in)  :: listrank(npoin_own+1:npoin_2)   ! actually more than needed
    integer(ip),intent(out) :: listrank_ndofn(npoin_own*ndofn+1:npoin_2*ndofn)

    integer(ip) :: ipoin,idofn,iauxi

    do ipoin = npoin_own+1,npoin_2  ! more than actually needed some will be 0
       iauxi = listrank(ipoin)
       do idofn = 1,ndofn
          !
          ! added -1 because I was getting Subscript #1 of the array LPROC has value 4 which is greater than the upper bound of 3
          ! in 5160  dagmg_par.f90, I relate this with the fact that the master is not included
          !
          listrank_ndofn((ipoin-1)*ndofn+idofn) = iauxi - 1   
       end do
    end do

  end subroutine calc_listrank_ndofn



  subroutine test_ricard(solve,amatr,aa_agmg,ia_agmg,nd_aux, nz_2, npoin_own)

    use def_domain,     only :  meshe ! to avoid passing meshe that I do not have in alya2agmg-solution
    use def_kermod,     only :  ndivi
    use def_master,     only :  INOTMASTER


    type(soltyp),      intent(in)    :: solve                ! Actually only for nnz

    integer(ip),       intent(in)    :: nd_aux, nz_2, npoin_own

    real(rp),          intent(in)    :: amatr(nd_aux,nd_aux,solve % nnz) !< Matrix
    real(rp),          intent(in)    :: aa_agmg(nd_aux * nd_aux * nz_2)   
    integer(ip),       intent(in)    :: ia_agmg( (npoin_own * nd_aux) + 1 )  

    real(rp),allocatable    :: xx(:)
    real(rp),allocatable    :: bb(:)

    real(rp)                :: suma,max_agmg_row
    integer(ip)             :: j,ipoin,npdiff

    integer(ip),parameter   :: nconn=7
    integer(ip)             :: jpoin(nconn),kpoin,l,ixy,jxy



    allocate(xx(nd_aux*max(1_ip,meshe(ndivi) % npoin)))
    allocate(bb(nd_aux*max(1_ip,meshe(ndivi) % npoin)))

    xx = 1.0_rp   
    bb = 0.0_rp

    ! bb = amatr * x
    call bcsrax( 1_ip, meshe(ndivi) % npoin, nd_aux, amatr, meshe(ndivi) % c_dom, meshe(ndivi) % r_dom, xx, bb ) 
    if( INOTMASTER ) then
       !
       ! check
       !   
       npdiff = (meshe(ndivi) % npoi2-meshe(ndivi) % npoi1-1)*nd_aux
       do ipoin=1,npoin_own*nd_aux
          suma = 0.0_rp
          do j = ia_agmg(ipoin),ia_agmg(ipoin+1)-1
             suma = suma + aa_agmg(j)
          end do
          max_agmg_row = maxval(abs(aa_agmg(ia_agmg(ipoin):ia_agmg(ipoin+1)-1)))
          if (max_agmg_row == 0.0_rp) then
             write(kfl_paral+5000,'(a,1x,i8,200(1x,e14.7))') 'ipoin,bb(ipoin),aa_agmg(ia_agmg(ipoin):ia_agmg(ipoin+1)-1)',ipoin,bb(ipoin),aa_agmg(ia_agmg(ipoin):ia_agmg(ipoin+1)-1)
             !             call runend ('mod_alya2agmg:test_ricard--max_agmg_row == 0.0_rp')
          end if
          if (ipoin <= meshe(ndivi) % npoi1 * nd_aux) then
             if (max_agmg_row /= 0.0_rp) then
                if ( (bb(ipoin)- suma)/max_agmg_row >1.0e-6_rp)  write(kfl_paral+3000,'(a,1x,i8,200(1x,e14.7))') &
                     'ipoin,bb(ipoin),suma',ipoin,bb(ipoin),suma,aa_agmg(ia_agmg(ipoin):ia_agmg(ipoin+1)-1)
             else
                if ( abs(bb(ipoin)- suma) >1.0e-8_rp)  write(kfl_paral+3000,'(a,1x,i8,200(1x,e14.7))') &
                     'ipoin,bb(ipoin),suma',ipoin,bb(ipoin),suma,aa_agmg(ia_agmg(ipoin):ia_agmg(ipoin+1)-1)
             end if

          else
             if (max_agmg_row /= 0.0_rp) then
                if ( (bb(ipoin+npdiff) - suma)/max_agmg_row >1.0e-6_rp) write(kfl_paral+3000,'(a,1x,i8,200(1x,e14.7))') &
                     'ipoin,bb(ipoin+npdiff),suma',ipoin,bb(ipoin+npdiff),suma,aa_agmg(ia_agmg(ipoin):ia_agmg(ipoin+1)-1)
             else
                if ( abs(bb(ipoin+npdiff) - suma) >1.0e-8_rp) write(kfl_paral+3000,'(a,1x,i8,200(1x,e14.7))') &
                     'ipoin,bb(ipoin+npdiff),suma',ipoin,bb(ipoin+npdiff),suma,aa_agmg(ia_agmg(ipoin):ia_agmg(ipoin+1)-1)             
             end if
          end if

       end do
    end if

    if (1==2) then  ! This part can help for looking column by column once 1 row were they differ has been identified
       jpoin(1)=89  ! all nodes connected to the one that gives wrong
       jpoin(2)=99
       jpoin(3)=108
       jpoin(4)=103
       jpoin(5)=77
       jpoin(6)=68
       jpoin(7)=80


       !    ixy =1 !x
       ixy =2 !y

       do jxy=1,nd_aux
          do l=1,nconn
             xx = 0.0_rp
             bb = 0.0_rp
             do ipoin=1,meshe(ndivi) % npoin
                if (lninv_loc(ipoin)==jpoin(l)) then
                   xx((ipoin-1)*nd_aux+jxy) = 1.0_rp
                end if
             end do
             call bcsrax( 1_ip, meshe(ndivi) % npoin, nd_aux, amatr, meshe(ndivi) % c_dom, meshe(ndivi) % r_dom, xx, bb ) 
             if( INOTMASTER ) then
                do ipoin=1,meshe(ndivi) % npoin 
                   if (lninv_loc(ipoin)==89) then
                      kpoin=(ipoin-1)*nd_aux+ixy
                      write(kfl_paral+5000,'(a,1x,i8,200(1x,e14.7))') '!!!!!kpoin,bb(kpoin+npdiff)',kpoin,bb(kpoin+npdiff)
                   end if
                end do
             end if
          end do
       end do

    end if


    deallocate(xx)
    deallocate(bb)


  end subroutine test_ricard



end module mod_alya2agmg
!> @}
