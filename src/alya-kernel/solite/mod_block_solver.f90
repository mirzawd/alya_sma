!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!




!-----------------------------------------------------------------------
!> @addtogroup Solver
!> @{
!> @file    mod_orthomin.f90
!> @author  guillaume
!> @date    2021-03-29
!> @brief   Orthomin
!> @details Block Orthomin
!-----------------------------------------------------------------------

module mod_block_solver

  use def_kintyp_basic,      only : ip,rp
  use mod_matrix,            only : matrix_CSR_SpMV
  use mod_solver,            only : solver_parallel_SpMV => solver_parallel_SpMV_0
  use mod_solver,            only : solver_vector_difference_norm
  use mod_memory_basic,      only : memory_alloca
  use mod_memory_basic,      only : memory_deallo
  use mod_solver,            only : solver_old_to_new
  use mod_solver,            only : solver_old_to_blk
  use def_mat_blk,           only : mat_blk
  use def_mat_csr,           only : mat_csr
  use def_mat_coo,           only : mat_coo
  use def_mat_ell,           only : mat_ell
  use def_mat_den,           only : mat_den
  use mod_solver_operations, only : parallel_mv
  use mod_solver_operations, only : parallel_exchange
  use mod_solver_operations, only : parallel_L2norm
  use mod_solver_operations, only : parallel_double_dot_product
  use def_solver
  implicit none

  private
  
  character(20), parameter :: vacal = 'mod_orthomin'
  real(rp),      parameter :: epsil = epsilon(1.0_rp)
  integer(ip)              :: n1
  integer(ip)              :: n2
  integer(ip)              :: ndof1
  integer(ip)              :: ndof2
  type(soltyp),  pointer   :: solve1(:)
  type(soltyp),  pointer   :: solve2(:)
  real(rp),      pointer   :: u1(:)
  real(rp),      pointer   :: u2(:)
  real(rp),      pointer   :: b1(:)
  real(rp),      pointer   :: b2(:)

  public :: block_solver
  
contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-04-06
  !> @brief   Initialization
  !> @details Initialization of block solver
  !> 
  !-----------------------------------------------------------------------

  subroutine block_init(solve,amatr,rhsid,unkno,qmatr)

    type(soltyp),                     intent(inout) :: solve
    real(rp),               pointer,  intent(in)    :: amatr(:)
    real(rp),               pointer,  intent(in)    :: rhsid(:)
    real(rp),               pointer,  intent(in)    :: unkno(:)
    real(rp),     optional, pointer,  intent(in)    :: qmatr(:)
    
    nullify(solve1)
    nullify(solve2)
    nullify(u1)
    nullify(u2)

    n1     =  solve % block_solve(1) % nequa
    n2     =  solve % block_solve(2) % nequa
    ndof1  =  solve % block_solve(1) % ndofn
    ndof2  =  solve % block_solve(2) % ndofn
    solve1 => solve % block_solve(1:)
    solve2 => solve % block_solve(2:)

    if( n1 > 0 ) then
       u1     => unkno(solve % block_rhs(1):solve % block_rhs(2)-1)
       u2     => unkno(solve % block_rhs(2):solve % block_rhs(2)+ndof2*n2)
       b1     => rhsid(solve % block_rhs(1):)
       b2     => rhsid(solve % block_rhs(2):)
    else
       allocate(u1(1),u2(1),b1(1),b2(1))
    end if
    
  end subroutine block_init

  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-04-06
  !> @brief   End
  !> @details End of block solver
  !> 
  !-----------------------------------------------------------------------

  subroutine block_end(solve,amatr,rhsid,unkno,qmatr)

    type(soltyp),                     intent(inout) :: solve
    real(rp),               pointer,  intent(in)    :: amatr(:)
    real(rp),               pointer,  intent(in)    :: rhsid(:)
    real(rp),               pointer,  intent(in)    :: unkno(:)
    real(rp),     optional, pointer,  intent(in)    :: qmatr(:)

    if( n1 > 0 ) then
       continue
    else
       deallocate(u1,u2,b1,b2)
    end if
    
  end subroutine block_end

  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-04-06
  !> @brief   Block solver driver
  !> @details IBlock solver driver
  !> 
  !-----------------------------------------------------------------------
  
  subroutine block_solver(solve,amatr,rhsid,unkno,qmatr)

    type(soltyp),                                    intent(inout) :: solve
    real(rp),                  contiguous, pointer,  intent(in)    :: amatr(:)
    real(rp),                  contiguous, pointer,  intent(in)    :: rhsid(:)
    real(rp),                  contiguous, pointer,  intent(in)    :: unkno(:)
    real(rp),      optional,   contiguous, pointer,  intent(in)    :: qmatr(:)
    type(mat_blk)                                                  :: A
    class(*),                              pointer                 :: Q
    
    call solver_old_to_blk(solve,amatr,A)
    !
    ! Initialize
    !
    call block_init(solve,amatr,rhsid,unkno,qmatr)
    !
    ! Call block solver
    !
    select case ( solve % kfl_algso )
       
    case ( SOL_SOLVER_GAUSS_SEIDEL )
       
       call block_gauss_seidel(solve,A)
       
    case ( SOL_SOLVER_ORTHOMIN )
       
       if( present(qmatr) ) then
          nullify(Q)
          call solver_old_to_new(solve2(1),qmatr,Q)
          call block_orthomin   (solve,A,Q)
          if( associated(Q) )   deallocate(Q)         
       end if
       
    case default
       call runend('MOD_BLOCK_SOLVER: SOLVER DOES NOT EXIST') 
    end select
    !
    ! End
    !
    call block_end(solve,amatr,rhsid,unkno,qmatr)
    
 end subroutine block_solver

 !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-03-29
  !> @brief   Block Gauss-Seidel
  !> @details Block Gauss-Seidel
  !> 
  !-----------------------------------------------------------------------

  subroutine block_gauss_seidel(solve,A)

    type(soltyp),           intent(inout) :: solve
    type(mat_blk)                         :: A
    real(rp),     pointer                 :: u1_new(:)
    real(rp),     pointer                 :: u2_new(:)
    real(rp),     pointer                 :: rhs1(:)
    real(rp),     pointer                 :: rhs2(:)
    integer(ip)                           :: i1,i2
    real(rp)                              :: relax,resid
    real(rp)                              :: rela1
        
    nullify(u1_new,u2_new,rhs1,rhs2)    

    call memory_alloca(memit,'U1_NEW',vacal,u1_new,n1*ndof1)
    call memory_alloca(memit,'U2_NEW',vacal,u2_new,n2*ndof2)
    call memory_alloca(memit,'RHS1'  ,vacal,rhs1  ,n1*ndof1)
    call memory_alloca(memit,'RHS2'  ,vacal,rhs2  ,n2*ndof2)
    do i1 = 1,n1*ndof1
       u1_new(i1) = u1(i1)
    end do
    do i2 = 1,n2*ndof2
       u2_new(i2) = u2(i2)
    end do
    
    relax         = solve % relax
    rela1         = 1.0_rp - relax 
    solve % iters = 0
    resid         = 1.0_rp
    solve % resin = resid
    solve % resi2 = resid
    
    do while( solve % iters < solve % miter .and. resid > solve % solco )

       resid         = 0.0_rp
       solve % iters = solve % iters + 1 
       !
       ! Solve first block: A11 u1 = b1 - A12 u2
       !
       call parallel_mv(solve1(1),a % vA(1,2),u2,rhs1,MPI=.false.)       
       do i1 = 1,n1*ndof1
          rhs1(i1) = b1(i1) - rhs1(i1)
       end do
       call block_solve(solve1,rhs1,u1_new,a % vA(1,1))
       resid = max(resid,solver_vector_difference_norm(solve1(1),2_ip,u1,u1_new,relax))
       do i1 = 1,n1*ndof1
          u1(i1)     = relax * u1_new(i1) + rela1 * u1(i1)
          u1_new(i1) = u1(i1)
       end do
       !
       ! Solve second block: A22 u2 = b2 - A21 u1
       !
       call parallel_mv(solve2(1),a % vA(2,1),u1,rhs2,MPI=.false.)
       resid = max(resid,solver_vector_difference_norm(solve2(1),2_ip,u2,u2_new,relax))
       do i2 = 1,n2*ndof2
          rhs2(i2) = b2(i2) - rhs2(i2)
       end do
       call block_solve(solve2,rhs2,u2_new,a % vA(2,2))
       do i2 = 1,n2*ndof2
          u2(i2)     = relax * u2_new(i2) + rela1 * u2(i2)
          u2_new(i2) = u2(i2)
       end do       
       
    end do
    
    solve % resfi = resid
    
    call memory_deallo(memit,'U1_NEW',vacal,u1_new)
    call memory_deallo(memit,'U2_NEW',vacal,u2_new)
    call memory_deallo(memit,'RHS1'  ,vacal,rhs1  )
    call memory_deallo(memit,'RHS2'  ,vacal,rhs2  )
           
  end subroutine block_gauss_seidel

  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-04-12
  !> @brief   Block solve
  !> @details Block solve
  !> 
  !-----------------------------------------------------------------------

  subroutine block_solve(solve,b,x,A)
    
    type(soltyp), pointer, intent(inout) :: solve(:)
    real(rp),     pointer, intent(inout) :: b(:)
    real(rp),     pointer, intent(inout) :: x(:)
    class(*),              intent(inout) :: A
    real(rp),     target                 :: y(2)
    real(rp),     pointer                :: bb(:)
    real(rp),     pointer                :: xx(:)

    solve_sol => solve

    if( associated(x) ) then
       xx => x
    else
       xx => y
    end if
    if( associated(b) ) then
       bb => b
    else
       bb => y
    end if

    select type ( A )
       
    class is ( mat_csr )
       
       if( associated(A % vA) ) then
          call solver(bb,xx,A % vA)
       else
          call solver(bb,xx,y)
       end if
       
    class is ( mat_coo )
       
       if( associated(A % vA) ) then
          call solver(bb,xx,A % vA)
       else
          call solver(bb,xx,y)
       end if
       
    class is ( mat_ell )
       
       if( associated(A % vA) ) then
          call solver(bb,xx,A % vA)
       else
          call solver(bb,xx,y)
       end if
       
    class is ( mat_den )
       
       if( associated(A % vA) ) then
          call solver(bb,xx,A % vA)
       else
          call solver(bb,xx,y)
       end if
       
    end select
    
    
  end subroutine block_solve
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-05-10
  !> @brief   Schur complement system with a preconditioned Orthomin(1) iteration
  !> @details Schur complement system with a preconditioned Orthomin(1) iteration
  !>          INPUT  u1^{k+1}
  !>          OUTPUT u1^{k+2}, u2^{k+1}, r2^{k+1}
  !>
  !>          A11 u1^{k+1} = b1 - A12 u2^k
  !>          r2^k         = b2 - A21 u1^{k+1} - A22 u2^{k}
  !>          Q du2        = r2^k
  !>          A11 du1      = A12 du2
  !>          x2           = A22 du2 - A21 du1
  !>          alpha        = (r2^k,x2) / (x2,x2)
  !>
  !>          u2^{k+1}     = u2^k     + alpha * du2
  !>          u1^{k+2}     = u1^{k+1} - alpha * du1
  !> 
  !-----------------------------------------------------------------------
  
  subroutine block_orthomin(solve,A,Q)

    type(soltyp),           intent(inout) :: solve
    type(mat_blk),          intent(inout) :: A
    class(*),               intent(inout) :: Q
    integer(ip)                           :: i1,i2
    real(rp)                              :: ak,denom,numer,resid,invb2
    real(rp),       pointer               :: x1(:)
    real(rp),       pointer               :: x2(:)
    real(rp),       pointer               :: r2(:)
    real(rp),       pointer               :: du1(:)
    real(rp),       pointer               :: du2(:)

    nullify(x1,x2,r2,du1,du2)

    call memory_alloca(memit,'DU1',vacal,du1,n1*ndof1)
    call memory_alloca(memit,'DU2',vacal,du2,n2*ndof2)
    call memory_alloca(memit,'R2' ,vacal,r2, n2*ndof2)
    call memory_alloca(memit,'X1' ,vacal,x1, n1*ndof1)
    call memory_alloca(memit,'X2' ,vacal,x2, n2*ndof2)

    solve % iters = 0
    resid         = 1.0_rp
    solve % resin = resid
    solve % resi2 = resid
    !
    ! Initial residual and ||b2||
    !
    do i2 = 1,n2*ndof2
       r2(i2) = b2(i2)
    end do
    call parallel_exchange(solve2(1),r2)
    call parallel_L2norm  (solve2(1),r2,invb2)
    if( invb2 > epsil ) invb2 = 1.0_rp / sqrt(invb2)
    
    do while( solve % iters < solve % miter .and. resid > solve % solco )

       solve % iters = solve % iters + 1 
       !
       ! Solve first block: A11 u1^{k+1} = b1 - A12 u2^k
       !
       call parallel_mv(solve1(1),a % vA(1,2),u2,du1,MPI=.false.)

       do i1 = 1,n1*ndof1                                                               ! du1 <= b1 - A12 u2^k
          du1(i1) = b1(i1) - du1(i1)
       end do
       call block_solve(solve1,du1,u1,a % vA(1,1))
       !
       ! R2 = b2 - A21 u1^{k+1} - A22 u2^{k}
       !
       call parallel_mv(solve1(1),a % vA(2,1),u1,r2,MPI=.false.)                        ! r2 =  A21 u1^{k+1} 
       call parallel_mv(solve1(1),a % vA(2,2),u2,r2,MPI=.false.,INITIALIZATION=.false.) ! r2 <= r2 + A22 u2^k        
       do i2 = 1,n2*ndof2                                                               ! r2 <= b2 - r2
          r2(i2) = b2(i2) - r2(i2)
       end do
       !
       ! Solve system Q du2 = r2
       !
       do i2 = 1,n2*ndof2
          du2(i2) = 0.0_rp                                                              ! Du2=0
       end do
       call block_solve(solve2,r2,du2,Q)                                                ! Solve system
       !
       ! Solve A11 du1 = A12 du2
       !
       call parallel_mv(solve1(1),a % vA(1,2),du2,x1,MPI=.false.)                       ! x1 = A12 du2
       do i1 = 1,n1*ndof1
          du1(i1) = 0.0_rp
       end do
       call block_solve(solve1,x1,du1,a % vA(1,1))                                      ! Solve system
       !
       ! Compute x1 = A22 du2 - A21 du1
       !
       call parallel_mv(solve2(1),a % vA(2,2),du2,x2,MPI=.false.)                       ! x2 = A22 du2
       call parallel_mv(solve2(1),a % vA(2,1),du1,x1,MPI=.false.)                       ! x1 = A21 du1^{k+1}
       do i2 = 1,n2*ndof2
          x2(i2) = x2(i2) - x1(i2)
       end do
       call parallel_exchange(solve2(1),x2) 
       !
       ! Compute ak = (rr,x2) / (x2,x2)
       !
       call parallel_double_dot_product(solve2(1),x2,r2,x2,numer,denom)!,MY_TIMING,OPENMP)
       call parallel_L2norm            (solve2(1),r2,resid) 

       if( denom == 0.0_rp ) then
          ak = 1.0_rp
       else
          ak = numer / denom
       end if
       resid = resid * invb2
       !
       ! Actualize u2^{k+1} = u2^k     + ak * z
       !           u1^{k+2} = u1^{k+1} - ak * v
       !
       do i2 = 1,n2*ndof2
          u2(i2) = u2(i2) + ak * du2(i2)
       end do
       do i1 = 1,n1*ndof1
          u1(i1) = u1(i1) - ak * du1(i1)
       end do

    end do

    solve % resfi = resid

    call memory_deallo(memit,'DU1',vacal,du1)
    call memory_deallo(memit,'DU2',vacal,du2)
    call memory_deallo(memit,'R2' ,vacal,r2 )
    call memory_deallo(memit,'X1' ,vacal,x1 )
    call memory_deallo(memit,'X2' ,vacal,x2 )

  end subroutine block_orthomin
  
end module mod_block_solver
!> @}

