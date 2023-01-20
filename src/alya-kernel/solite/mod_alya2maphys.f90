!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup Algebraic_Solver
!> @{
!> @name    Bridge to MaPhys
!> @file    mod_alya2maphys.f90
!> @author  Guillaume Houzeaux
!> @date    07/03/2016
!> @brief   Bridge to MAPHYS
!> @details Module to call MAPHYS iterative solver
!------------------------------------------------------------------------
module mod_alya2maphys

  use def_kintyp,               only : ip,rp,lg
  use def_kintyp_solvers,       only : soltyp
  use mod_memory,               only : memory_alloca
  use mod_memory,               only : memory_deallo
  use def_solver,               only : memit
  use def_solver,               only : SOL_SOLVER_MAPHYS_UNSYMMETRIC
  use def_solver,               only : SOL_SOLVER_MAPHYS_SYMMETRIC
  use def_solver,               only : SOL_SOLVER_MAPHYS_SPD
  use def_solver,               only : SOL_MATRIX_HAS_CHANGED
  use mod_matrix,               only : matrix_full_to_half
  use def_kintyp_comm,          only : comm_data_par 
  use def_domain,               only : mesh_type
  use def_master,               only : INOTMASTER,IMASTER,kfl_paral
  use mod_iofile,               only : iofile
  use mod_parall,               only : PAR_MY_CODE_RANK_WM
  use mod_parall,               only : PAR_COMM_MY_CODE_WM
  use mod_communications,       only : PAR_SUM,PAR_MAX
  use mod_communications,       only : PAR_COMM_RANK_AND_SIZE
  use mod_communications,       only : PAR_ALLGATHER
  use mod_communications,       only : PAR_INTERFACE_NODE_EXCHANGE
  use mod_communications,       only : PAR_BROADCAST
  use mod_par_output_partition, only : par_output_maphys
  use mod_par_output_partition, only : par_output_global_matrix
  use mod_graphs,               only : graphs_csr_to_coo
  use mod_matrix,               only : matrix_CSR_SpMV
  use mod_matrix,               only : matrix_csr_to_coo
  use mod_messages,             only : livinf
  use mod_communications_tools, only : PAR_COMM_TO_INT

#ifdef MAPHYS
  use dmph_maphys_mod
  use dmph_toolkit_mod
#endif
  implicit none
#ifdef MAPHYS
#include "mph_defs_f.h"
#endif
  private

  integer(ip) :: mem_maphys_facto     ! IINFO(9)
  integer(ip) :: mem_maphys_schur     ! IINFO(12)
  integer(ip) :: mem_maphys_precond   ! IINFO(17)
  integer(ip) :: mem_maphys_solver    ! IINFO(37)

  real(rp)    :: time_maphys_facto    ! RINFO(5)
  real(rp)    :: time_maphys_analysis ! RINFO(4)
  real(rp)    :: time_maphys_precond  ! RINFO(6)
  real(rp)    :: time_maphys_solver   ! RINFO(7)

  public :: alya2maphys_initialization
  public :: alya2maphys_solution
  public :: mem_maphys_facto    
  public :: mem_maphys_schur    
  public :: mem_maphys_precond  
  public :: mem_maphys_solver   
  public :: time_maphys_facto    
  public :: time_maphys_analysis    
  public :: time_maphys_precond  
  public :: time_maphys_solver   
   
contains 

  !-----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    14/03/1016
  !> @brief   Initialize MAPHYS
  !> @details Initialize MAPHYS for a given problem
  ! 
  !-----------------------------------------------------------------------

  subroutine alya2maphys_initialization(meshe,COMM,solve)

    type(mesh_type),     intent(inout) :: meshe                !< Mesh type
    type(comm_data_par), intent(in)    :: COMM                 !< Communication array
    type(soltyp),        intent(inout) :: solve                !< Solver type
    integer(ip)                        :: info,ipoin,kpoin
    integer(ip)                        :: ii,jj,ineig
    integer(ip)                        :: jdofn
    integer(ip)                        :: jpoin
    integer(ip)                        :: nnz
    integer(ip)                        :: ndofn
    integer(ip),         pointer       :: num_interface_nodes(:)
    integer(ip),         pointer       :: interface_nodes(:)
    integer(ip),         pointer       :: rows(:),cols(:)

    ! MaPHyS specific, I32 only
    integer(4)                         :: job,sym    
    integer(4)                         :: myndof,mysizeIntrf,mynbvi,sizeIntrf
    integer(4),          pointer       :: myinterface(:)
    integer(4),          pointer       :: myindexVi(:)
    integer(4),          pointer       :: myptrindexVi(:)
    integer(4),          pointer       :: myindexintrf(:)

#ifdef MAPHYS
    character(MAPHYS_STRL)             :: matrixfile,rhsfile,initguessfile
    character(MAPHYS_STRL)             :: outrhsfile,outsolfile
#endif

#ifndef MAPHYS

    call runend('ALYA2MAPHYS_INITIALIZATION: LINK WITH MAPHYS AND USE -DMAPHYS')

#else

    call livinf(0_ip,'INITIALIZE MAPHYS FOR SOLVER: '//trim(solve % wprob),0_ip)

    if( INOTMASTER ) then
       !
       ! Nullify MAPHYS communication data structure
       !
       ndofn       = solve % ndofn
       myndof      = 0
       mynbvi      = 0
       mysizeIntrf = 0
       nullify(myinterface )
       nullify(myindexVi   )
       nullify(myptrindexVi) 
       nullify(myindexintrf)
       nullify(num_interface_nodes)
       nullify(interface_nodes)
       !
       ! Initialize MAPHYS instance
       !
       solve % mphs % comm = PAR_COMM_TO_INT(PAR_COMM_MY_CODE_WM)   ! communicator without master
       solve % mphs % job  = -1
       call DMPH_maphys_driver(solve % mphs)
       info = solve % mphs % iinfog(1)

       if( info < 0 ) call runend('ALYA2MAPHYS: INITIALIZATION FAILED')

       !-----------------------------------------------------------------
       !
       ! Default values
       !
       !-----------------------------------------------------------------
       !
       ! Output
       !
       solve % mphs % ICNTL( 1) = solve % lun_exsol   ! Output unit for error messages
       solve % mphs % ICNTL( 2) = solve % lun_exsol   ! Output unit for warning messages
       solve % mphs % ICNTL( 3) = solve % lun_exsol   ! Output unit for statistics messages
       solve % mphs % ICNTL( 4) = 3                   ! Print errors
       solve % mphs % ICNTL( 5) = 1                   ! Print ICNTL
       solve % mphs % ICNTL( 6) = 1                   ! Output output solver info (0=nothing, 1=end of the solve, 2=all)
       solve % mphs % ICNTL(14) = 1                   ! Output in org format
       !
       ! Solver
       !
       solve % mphs % ICNTL(20) = 3                   ! Let MAPHYS decide on the solver (CG if SPD, GMRES otherwise)
       !
       ! GMRES
       !
       if(      solve % kfl_ortho == 0 ) then
          solve % mphs % ICNTL(22) = 2                ! Classical Gram-Schmidt
       else if( solve % kfl_ortho == 1 ) then
          solve % mphs % ICNTL(22) = 0                ! Modified  Gram-Schmidt
       else if( solve % kfl_ortho == 3 ) then
          solve % mphs % ICNTL(22) = 3                ! Iterative classical Gram-Schmidt
       end if
       solve % mphs % ICNTL(22) = 3                   ! Orthogonolization
       solve % mphs % ICNTL(26) = solve % nkryd       ! Max Krylov dimension
       !
       ! Convergence
       !
       solve % mphs % ICNTL(24) = solve % miter       ! Maximum number of total number of iteration
       solve % mphs % ICNTL(25) = 1                   ! Residual explicitly computed using a SpMV
       solve % mphs % ICNTL(27) = 1                   ! Schur complement assembled in a dense matrix (1=yes, 2=no: better memory usage and perf in 3D)
       solve % mphs % RCNTL(21) = solve % solco       ! Stopping criterion based on Schur. system backward error || b_s - S x_s || / || b_s ||
       !
       ! Preconditioner
       !
       solve % mphs % ICNTL(21) = 1                   ! Preconditioner (1=dense,2=sparse with threashold)
       solve % mphs % RCNTL(11) = 1.0e-6_rp           ! Threshold to sparsify Schur complement (if ICNTL(21)=2)
       !
       ! Others
       !
       solve % mphs % ICNTL(23) = 1                   ! (0=initial guess, 1=initial guess provided by user)
       solve % mphs % ICNTL(43) = 2                   ! Distributed interface

       !-----------------------------------------------------------------
       !
       ! Read configuration file
       !
       !-----------------------------------------------------------------

       if( PAR_MY_CODE_RANK_WM == 0 ) then
          if( trim(solve % conf_file) == "" ) then
             call runend('ALYA2MAPHYS_INITIALIZATION: NO CONFIGURATION FILE AVAILABLE. PUT CONFIG_FILE: conf_maphys.txt')
          else
             call iofile(0_ip,90_ip,trim(solve % conf_file),'MAPHYS CONFIG FILE','old')
          end if
          matrixfile    = ''
          rhsfile       = ''
          initguessfile = ''
          outrhsfile    = ''
          outsolfile    = ''
          call DMPH_read_param_freeformat(                              &
               90_ip,solve % mphs % icntl,solve % mphs % rcntl,sym,job, &
               matrixfile,rhsfile,initguessfile,                        &
               outrhsfile,outsolfile)          
          call iofile(2_ip,90_ip,trim(solve % conf_file),'MAPHYS CONFIG FILE')
       end if

       call PAR_BROADCAST(size(solve % mphs % ICNTL),             solve % mphs % ICNTL,'IN MY CODE WITHOUT MASTER')
       call PAR_BROADCAST(int(size(solve % mphs % RCNTL),kind=ip),solve % mphs % RCNTL,'IN MY CODE WITHOUT MASTER')
       !
       ! Data structures for local subdomains (cannot be deallocated)
       !
       myndof      = meshe % npoin * ndofn
       mysizeIntrf = ( meshe % npoin - meshe % npoi1 ) * ndofn
       mynbvi      = COMM % nneig 
       sizeIntrf   = ( COMM % bound_size(COMM % nneig+1) - 1_ip ) * ndofn

       call memory_alloca(memit,'myinterface' ,'alya2maphys_initialization',myinterface ,mysizeIntrf)
       call memory_alloca(memit,'myindexVi'   ,'alya2maphys_initialization',myindexVi   ,mynbvi)
       call memory_alloca(memit,'myptrindexVi','alya2maphys_initialization',myptrindexVi,mynbvi+1_ip)
       call memory_alloca(memit,'myindexintrf','alya2maphys_initialization',myindexintrf,sizeIntrf) 
       !
       ! List of interface nodes
       !
       kpoin = meshe % npoi1
       do ipoin = 1,mysizeIntrf
          kpoin = kpoin + 1
          myinterface(ipoin) = meshe % lninv_loc(kpoin) 
       end do
       !
       ! Neighbors: remove master and shift ranks
       !
       myindexVi(1:mynbvi) = COMM % neights(1:mynbvi) - 1

       if( ndofn == 1 ) then
          do ineig = 1,COMM % nneig 
             myptrindexVi(ineig) = COMM % bound_size(ineig)
             do ii = COMM % bound_size(ineig),COMM % bound_size(ineig+1)-1
                kpoin = meshe % lninv_loc(COMM % bound_perm(ii))
                jj    = 1
                do while( myinterface(jj) /= kpoin )
                   jj = jj + 1
                end do
                myindexintrf(ii) = jj
             end do
          end do
       else
          do ineig = 1,COMM % nneig 
             myptrindexVi(ineig) = ( COMM % bound_size(ineig) - 1 ) * ndofn + 1
             do ii = COMM % bound_size(ineig),COMM % bound_size(ineig+1)-1
                kpoin = meshe % lninv_loc(COMM % bound_perm(ii))
                jj    = 1
                do while( myinterface(jj) /= kpoin )
                   jj = jj + 1
                end do
                do jdofn = 1,ndofn
                   myindexintrf( (ii-1)*ndofn + jdofn) = (jj-1)*ndofn + jdofn
                end do
             end do
          end do
          kpoin = meshe % npoi1
          jpoin = 0
          do ipoin = 1,mysizeIntrf / ndofn
             kpoin = kpoin + 1
             do jdofn = 1,ndofn
                jpoin = jpoin + 1
                myinterface(jpoin) = (meshe % lninv_loc(kpoin)-1)*ndofn + jdofn
             end do
          end do

       end if
       print*,'b=',kfl_paral

       myptrindexVi(COMM % nneig+1) = ( COMM % bound_size(COMM % nneig+1) - 1 ) * ndofn + 1
       !
       ! Row-columm graph
       ! 
       nullify(solve % mphs % rows,solve % mphs % cols)     
       nullify(rows,cols)     
       if( solve % kfl_symeq == 1 ) then
          call graphs_csr_to_coo(solve % nequa,ndofn,solve % ia,solve % ja,nnz,rows,cols,'SYMMETRIC',memor=memor_dom)
       else
          call graphs_csr_to_coo(solve % nequa,ndofn,solve % ia,solve % ja,nnz,rows,cols,memor=memor_dom)
       end if

       solve % mphs % nnz  = nnz
       allocate(solve % mphs % rows(nnz))
       allocate(solve % mphs % cols(nnz))
       solve % mphs % rows = rows 
       solve % mphs % cols = cols
       !
       ! MAPHYS: Create domain 
       !
       call DMPH_distributed_interface(&
            solve % mphs , myndof       , &
            mysizeIntrf  , myinterface  , &
            mynbvi       , myindexVi    , &
            myptrindexVi , myindexintrf   )
       !
       ! Deallocate memory
       !
       call memory_deallo(memit,'myinterface' ,'alya2maphys_initialization',myinterface )
       call memory_deallo(memit,'myindexVi'   ,'alya2maphys_initialization',myindexVi   )
       call memory_deallo(memit,'myptrindexVi','alya2maphys_initialization',myptrindexVi)
       call memory_deallo(memit,'myindexintrf','alya2maphys_initialization',myindexintrf) 
       call memory_deallo(memit,'rows'        ,'alya2maphys_initialization',rows        ) 
       call memory_deallo(memit,'cols'        ,'alya2maphys_initialization',cols        ) 
       print*,'c=',kfl_paral

    end if
#endif

  end subroutine alya2maphys_initialization

  !-----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    14/03/1016
  !> @brief   Solve Ax=b
  !> @details Solve algebraic system using MAPHYS.
  !>          MAPHYS requires non-assemble matrix and RHS
  !
  !-----------------------------------------------------------------------

  subroutine alya2maphys_solution(solve,amatr,rhsid,unkno)

#ifdef MAPHYS
    use DMPH_schur_mod
#endif

    type(soltyp),                intent(inout) :: solve                !< Solver type
    real(rp),            target, intent(in)    :: amatr(solve % nzmat) !< Matrix
    real(rp),            target, intent(inout) :: rhsid(solve % nunkn) !< RHS
    real(rp),            target, intent(inout) :: unkno(solve % nunkn) !< Solution

#ifdef MAPHYS
    integer(ip)                                :: iallo_values
    integer(ip)                                :: ndofn,nzmat
    real(rp)                                   :: rnorm,bnorm
    real(rp),            pointer               :: resid(:),rhsid_cpy(:)


    !--------------------------------------------------------------------
    !
    ! Initial residual: residual must be exchanged
    !
    !--------------------------------------------------------------------

    nullify(resid,rhsid_cpy)
    ndofn = solve % ndofn

    if( INOTMASTER ) then
       allocate(resid    (solve % nunkn)) 
       allocate(rhsid_cpy(solve % nunkn))
       call matrix_CSR_SpMV(1_ip,solve % nequa,ndofn,ndofn,ndofn,ndofn,solve % ia,solve % ja,amatr,unkno,resid)
       rhsid_cpy(1:solve % nunkn) = rhsid(1:solve % nunkn)
       call PAR_INTERFACE_NODE_EXCHANGE(ndofn,resid    ,'SUM','IN MY CODE')
       call PAR_INTERFACE_NODE_EXCHANGE(ndofn,rhsid_cpy,'SUM','IN MY CODE')
       resid(1:solve % nunkn) = rhsid_cpy(1:solve % nunkn) - resid(1:solve % nunkn) 
    else
       allocate(resid    (1)) 
       allocate(rhsid_cpy(1))
    end if

    call norm2x(ndofn,resid,    rnorm)
    call norm2x(ndofn,rhsid_cpy,bnorm)   
    solve % bnorm = bnorm
    
    deallocate(rhsid_cpy)
    deallocate(resid)

    if( bnorm <= solve % bnorm_min ) then

       !--------------------------------------------------------------------
       !
       ! RHS norm is zero
       !
       !--------------------------------------------------------------------

       if( INOTMASTER ) unkno(1:solve % nunkn) =  0.0_rp
       solve % resin   =  0.0_rp 
       solve % resi2   =  0.0_rp
       solve % resfi   =  0.0_rp 
       solve % resf2   =  0.0_rp

    else

       solve % resin = rnorm / bnorm 
       solve % resi2 = solve % resin
       !
       ! Adaptive residual
       !
       if( solve % kfl_adres == 1 ) then
          solve % solco = max( solve % resin * solve % adres , solve % solmi )
          solve % mphs % RCNTL(21) = solve % solco  
       end if

       !--------------------------------------------------------------------
       !
       ! Solve system using MAPHYS
       !
       !--------------------------------------------------------------------

       if( INOTMASTER ) then
          !
          ! Matrix. Copy half of matrix if problem is symmetric
          !
          iallo_values = 0
          if( ndofn == 1 ) then
             if( solve % kfl_symeq == 1 ) then
                !
                ! dof=1, symmetric matrix
                !
                iallo_values = 1 
                allocate(solve % mphs % values(solve % mphs % nnz))
                call matrix_full_to_half(solve % nequa,solve % ndofn,solve % ia,solve % ja,amatr,solve % mphs % values)
             else
                !
                ! dof=1, unsymmetric matrix
                !
                solve % mphs % values => amatr
             end if
          else         
             if( solve % kfl_symeq == 1 ) then
                !
                ! dof>1, symmetric matrix
                !
                iallo_values = 1
                nzmat        = solve % nzmat
                allocate( solve % mphs % values(nzmat) )
                call matrix_csr_to_coo(solve % nequa,solve % ndofn,solve % ia,solve % ja,amatr,solve % mphs % values,'SYMMETRIC',memor=memor_dom)
             else
                !
                ! dof>1, unsymmetric matrix
                !
                iallo_values = 1
                nzmat        = solve % nzmat
                allocate( solve % mphs % values(nzmat) )
                call matrix_csr_to_coo(solve % nequa,ndofn,solve % ia,solve % ja,amatr,solve % mphs % values,memor=memor_dom)             
             end if

          end if

          solve % mphs % n   =  solve % nequa * ndofn
          solve % mphs % rhs => rhsid
          solve % mphs % sol => unkno
          solve % resfi      =  solve % resin
          !
          ! If matrix has just been assembled, recomptue everything
          ! If not, matrix has not changed
          !
          if( solve % kfl_assem == SOL_MATRIX_HAS_CHANGED ) then
             solve % mphs % job =  6 
          else
             solve % mphs % job =  4
          end if
          !
          ! Solver: Unsymmetric (0), SPD (1), symmetric(2)
          !
          if(      solve % kfl_algso == SOL_SOLVER_MAPHYS_UNSYMMETRIC ) then
             solve % mphs % sym =  0   
          else if( solve % kfl_algso == SOL_SOLVER_MAPHYS_SPD ) then
             solve % mphs % sym =  1
          else if( solve % kfl_algso == SOL_SOLVER_MAPHYS_SYMMETRIC ) then
             solve % mphs % sym =  2 
          end if
          !
          ! Call MAPHYS
          !
          call DMPH_maphys_driver(solve % mphs)
          !
          ! Deallocate matrix in the symmetric case
          !
          if( iallo_values == 1 ) then
             deallocate(solve % mphs % values)        
          end if
          !
          ! Get solver info. Final eesidual in MAPHYS is computed as ||b-Ax||/||b|
          ! But inside MAPHYS as ||b_s-Sx_s||/||b|| = ||b-Ax||/||b| because Schur
          ! is explicitly computed. Final Schur error is stored in RINFOG(3)
          !
          solve % resfi = solve % mphs % RINFOG(4)

          !call matrix_CSR_SpMV(solve % nequa,solve % ndofn,solve % ia,solve % ja,amatr,unkno,resid)
          !resid(1:solve % nunkn) = rhsid(1:solve % nunkn) - resid(1:solve % nunkn) 
          !call PAR_INTERFACE_NODE_EXCHANGE(solve % ndofn,resid,'SUM','IN MY CODE')

       end if

       !--------------------------------------------------------------------
       !
       ! Final residual
       !
       !--------------------------------------------------------------------

       call PAR_MAX(solve % resfi,'IN MY CODE')
       solve % resf2 = solve % resfi 
       !
       ! Timings and iterations
       !
       if( INOTMASTER ) then
          solve % cputi(3) = solve % cputi(3) + solve % mphs % RINFO(30)
       else
          solve % cputi(3) = solve % cputi(3) + 0.0_rp
       end if
       solve % iters    = solve % mphs % IINFOG(5)
       call PAR_MAX(solve % iters,'IN MY CODE')

    end if

    !--------------------------------------------------------------------
    !
    ! Alya requires RHS RHSID to be assembled
    !
    !--------------------------------------------------------------------
 
    call PAR_INTERFACE_NODE_EXCHANGE(ndofn,rhsid,'SUM','IN MY CODE')

    !--------------------------------------------------------------------
    !
    ! Statistics
    !
    !--------------------------------------------------------------------

    mem_maphys_facto     = int(solve % mphs % IINFO( 9),ip)   
    mem_maphys_schur     = int(solve % mphs % IINFO(12),ip)  
    mem_maphys_precond   = int(solve % mphs % IINFO(17),ip)
    mem_maphys_solver    = int(solve % mphs % IINFO(37),ip)

    time_maphys_facto    = real(solve % mphs % RINFO( 5),rp)
    time_maphys_analysis = real(solve % mphs % RINFO( 4),rp) 
    time_maphys_precond  = real(solve % mphs % RINFO( 6),rp)  
    time_maphys_solver   = real(solve % mphs % RINFO( 7),rp) 

    call par_output_maphys(&
         mem_maphys_facto, mem_maphys_schur, mem_maphys_precond ,mem_maphys_solver,&
         time_maphys_facto,time_maphys_analysis,time_maphys_precond,time_maphys_solver)

#endif

  end subroutine alya2maphys_solution

end module mod_alya2maphys
!> @} 
