!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!




!-----------------------------------------------------------------------
!> @addtogroup Solvers
!> @{
!> @file    def_dcg.f90
!> @author  houzeaux
!> @date    2022-03-01
!> @brief   Deflated conjugate gradient
!> @details Deflated conjugate gradient
!>          Input:
!>          ------
!>          ngrou ... number of groups
!>          lgrou ... list of groups
!>          Intermediate input:
!>          -------------------
!>          acoarse % iA
!>          acoarse % jA
!>          Output:
!>          -------
!>          acoarse % vA
!>          Factorization
!>         
!-----------------------------------------------------------------------

module def_dcg

  use def_kintyp,             only : ip,rp
  use def_mat,                only : mat
  use def_mat_csr,            only : mat_csr
  use def_mat_sky,            only : mat_sky
  use def_mat_den,            only : mat_den
  use def_iterative_solvers,  only : krylov_symmetric
  use def_direct_solvers,     only : direct_factorization
  use mod_memory_basic,       only : memory_alloca
  use mod_memory_basic,       only : memory_deallo
  use mod_memory_basic,       only : memory_size
  use def_solver,             only : SOL_SKYLINE
  use def_solver,             only : SOL_DENSE
  use def_solver,             only : SOL_CSR
  use def_lu_factorization,   only : lu_factorization
  use def_cholesky,           only : cholesky
  implicit none
  private
  
  type, extends(krylov_symmetric) :: dcg
     class(mat),                  pointer :: acoarse
     class(direct_factorization), pointer :: dir
     integer(ip)                          :: ngrou
     integer(ip),                 pointer :: lgrou(:)
     real(rp),                    pointer :: va(:,:,:)
     integer(ip),                 pointer :: ia(:)
     integer(ip),                 pointer :: ja(:)
   contains
     procedure,                   pass    :: init 
     procedure,                   pass    :: alloca
     procedure,                   pass    :: deallo
     procedure,                   pass    :: solve
     procedure,                   pass    :: setup
     procedure,                   pass    :: set
     ! Coarse solver
     procedure,                   pass    :: coarse_correction
     procedure,                   pass    :: set_coarse
     procedure,                   pass    :: fine_to_coarse 
     procedure,                   pass    :: coarse_to_fine
  end type dcg

  character(3), parameter :: vacal = 'dcg'

  public :: dcg
  
contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-03-10
  !> @brief   Initialize
  !> @details Initialize
  !> 
  !-----------------------------------------------------------------------
  
  subroutine init(self)
    
    class(dcg), intent(inout) :: self

    self % ngrou = 0
    nullify(self % acoarse)
    nullify(self % lgrou)
    nullify(self % va)
    nullify(self % ia)
    nullify(self % ja)
    nullify(self % dir)
    
  end subroutine init  
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-03-10
  !> @brief   Initialize
  !> @details Initialize
  !> 
  !-----------------------------------------------------------------------
  
  subroutine set(self,ngrou,lgrou,ia,ja)
    
    class(dcg),           intent(inout) :: self
    integer(ip),          intent(in)    :: ngrou
    integer(ip), pointer, intent(in)    :: lgrou(:)
    integer(ip), pointer, intent(in)    :: ia(:)
    integer(ip), pointer, intent(in)    :: ja(:)
    !
    ! Allocate coarse matrix
    !
    select case ( self % input % kfl_defas )
    case ( SOL_CSR     ) ; allocate( mat_csr :: self % acoarse )
    case ( SOL_DENSE   ) ; allocate( mat_den :: self % acoarse )
    case ( SOL_SKYLINE ) ; allocate( mat_sky :: self % acoarse )
    end select
    
    call self % acoarse % init  ()    
    if( ngrou > 0 ) then
       self % acoarse % nrows = ngrou
       self % acoarse % ncols = ngrou
    end if
    !
    ! Choose sparse solver
    !
    select case ( self % input % kfl_defas )
    case ( SOL_CSR     ) ; allocate( lu_factorization :: self % dir )
    case ( SOL_SKYLINE ) ; allocate( cholesky         :: self % dir )
    end select
    call self % dir % init()
    ! 
    ! Graph for coarse matrix
    !
    if( associated(self % acoarse) .and. ngrou > 0 ) then
       self % ngrou =  ngrou
       self % lgrou => lgrou
       select type ( ac => self % acoarse )
       class is ( mat_csr ) 
          ac   % nz    =  memory_size(ja)
          ac   % ia    => ia
          ac   % ja    => ja
          !select type ( dir => self % dir )
          !class is ( lu_factorization ) ; call dir % symbolical(self % acoarse)
          !end select
       class is ( mat_sky )
          call ac % csr2sky(ia=ia,ja=ja)
       end select
    end if

  end subroutine set
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-03-10
  !> @brief   Deallocate
  !> @details Deallocate
  !> 
  !-----------------------------------------------------------------------
  
  subroutine deallo(self)
     
    class(dcg), intent(inout) :: self

    if( associated(self % acoarse) ) then
       select type ( ac => self % acoarse )
       class is ( mat_csr ); call ac % deallo_matrix()
       class is ( mat_sky ); call ac % deallo_matrix()
       end select
       deallocate(self % acoarse)
    end if
    
  end subroutine deallo
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-03-10
  !> @brief   Allocate
  !> @details Allocate
  !> 
  !-----------------------------------------------------------------------
  
  subroutine alloca(self)

    class(dcg), intent(inout) :: self          !< Solver

  end subroutine alloca
 
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-03-10
  !> @brief   Setup
  !> @details Setup
  !> 
  !-----------------------------------------------------------------------
  
 subroutine setup(self,a)

    class(dcg),                         intent(inout) :: self          !< Solver
    class(mat),                         intent(inout) :: a             !< Matrix A

    self % acoarse % ndof1 = self % ndof
    self % acoarse % ndof2 = self % ndof

    if( associated(self % acoarse) ) then
       select type ( ac => self % acoarse )
       class is ( mat_csr ) ; call ac % alloca_matrix()
       class is ( mat_sky ) ; call ac % alloca_matrix()
       end select
    else
       call runend('DEF_DCG: COARSE SOLVER NOT ALLOCATED')
    end if
    
    !call self % dir % setup(self % acoarse) !hoy
    select type ( v => self % dir )
    class is ( direct_factorization )
       call v % set(self % acoarse)
    end select
    
  end subroutine setup
 
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-03-10
  !> @brief   Iteration
  !> @details Iteration
  !> 
  !-----------------------------------------------------------------------
  
  subroutine solve(self,b,x,a)
    
    class(dcg),               intent(inout) :: self              !< Solver
    real(rp),        pointer, intent(in)    :: b(:)              !< RHS b
    real(rp),        pointer, intent(inout) :: x(:)              !< Solve x
    class(mat),               intent(in)    :: a                 !< Matrix A
    integer(ip)                             :: n,nn
    integer(ip)                             :: ndof,ngrou
    integer(ip)                             :: ierro,i
    real(rp)                                :: alpha,beta,rho
    real(rp)                                :: newrho,resid
    real(rp)                                :: denom

    real(rp),        pointer                :: r(:)
    real(rp),        pointer                :: p(:)
    real(rp),        pointer                :: q(:)
    real(rp),        pointer                :: w(:)
    real(rp),        pointer                :: mu(:)
    real(rp),        pointer                :: murhs(:)
    
    nullify(r,p,q,w,mu,murhs)
    
    ndof  = self % ndof
    n     = self % n
    nn    = ndof * n
    ngrou = self % ngrou
    
    call memory_alloca(self % memor,'R'    ,vacal,r    ,nn)
    call memory_alloca(self % memor,'P'    ,vacal,p    ,nn)
    call memory_alloca(self % memor,'Q'    ,vacal,q    ,nn)
    call memory_alloca(self % memor,'W'    ,vacal,w    ,nn)
    call memory_alloca(self % memor,'MU'   ,vacal,mu   ,ndof*ngrou)
    call memory_alloca(self % memor,'MURHS',vacal,murhs,ndof*ngrou+1_ip)
    !
    ! Deflation  for initial solution
    !
    call self % alloca                 ()                        ! Allocate coarse matrix
    call self % set_coarse             (a)                       ! Set coarse matrix
    call self % dir % numerical        (self % acoarse)          ! A' = LU
    !
    ! Solve r <= ( W A'^-1 W^T ) r
    !
    call self % parallel_residual      (a,x,b,r)                 ! r = b - A x
    call self % coarse_correction      (mu,murhs,r)              ! r <= ( W A'^-1 W^T ) r
    do i = 1,nn                                                  ! x = x + ( W A'^-1 W^T ) r
       x(i) = x(i) + r(i)                                       
    end do

    call self % start(b,x,a,r,q,ierro)                           ! Starting operations

    if( ierro /= 0 ) goto 1   
    newrho = self % parallel_dot(r,q)
    if( newrho <= 0.0_rp ) goto 1
    resid = sqrt(newrho)    

    call self % init_iterations(resid)                           ! Residuals and tolerance
    do i = 1,nn
       p(i) = q(i)
    end do
    !
    ! p <= p - ( W A'^-1 W^T ) A q
    !
    call self % parallel_mv      (a,q,w)
    call self % coarse_correction(mu,murhs,w)
    do i = 1,nn
       p(i) = p(i) - w(i)
    end do

    do while( self % output % iters < self % input % miter .and. resid > self % output % toler )
       !
       ! q^{k+1} =  A p^{k+1}
       !
       call self % parallel_mv(a,p,q)
       !
       ! alpha = rho^k / <p^{k+1},q^{k+1}>
       !       
       denom = self % parallel_dot(p,q)
       if( denom <= 0.0_rp ) then
          ierro = 2
          goto 1
       end if
       rho   = newrho
       alpha = newrho / denom
       !
       ! x^{k+1} = x^k + alpha * p^{k+1}
       ! r^{k+1} = r^k - alpha * q^{k+1}
       !
       do i = 1,nn
          x(i) = x(i) + alpha * p(i)
          r(i) = r(i) - alpha * q(i)
       end do
       !
       !  L z^{k+1} = r^{k+1}
       !
       call self % preconditioning(r,q)
       newrho = self % parallel_dot(r,q)
       !
       ! Solve w = ( W A'^-1 W^T ) A q^k
       !
       call self % parallel_mv      (a,q,w)     
       call self % coarse_correction(mu,murhs,w)
       !
       ! beta  = rho^k / rho^{k-1}
       !
       beta = newrho / rho
       !
       ! p^{k+1} = z^k + beta*p^k
       !
       do i = 1,nn
          p(i) = q(i) + beta * p(i) - w(i)
       end do
       !
       ! Update and output residuals
       !
       resid                       = sqrt(newrho)
       self % output % resip_old   = self % output % resip_final
       self % output % resip_final = resid * self % output % bnorp_inv
       self % output % iters       = self % output % iters + 1
       call self % outcvg(b,x,a,r)
       
    end do

1   continue

    call self % end(b,x,a,r)             ! Ending operations
    call self % deallo()                 ! Deallocate coarse matrix
    call self % dir % partial_cleaning()
    
    call memory_deallo(self % memor,'R'    ,vacal,r)
    call memory_deallo(self % memor,'P'    ,vacal,p)
    call memory_deallo(self % memor,'Q'    ,vacal,q)
    call memory_deallo(self % memor,'W'    ,vacal,w)
    call memory_deallo(self % memor,'MU'   ,vacal,mu)
    call memory_deallo(self % memor,'MURHS',vacal,murhs)

    if( self % input % lun_cvgso /= 0 .and. self % iwrite .and. ierro > 0 ) then
       if( ierro == 2 ) write(self % input % lun_cvgso,2) self % output % iters
    end if
    
2   format(&
       & '# Error at iteration ',i6,&
       & 'Dividing by zero: alpha = rho^k / <p^{k+1},q^{k+1}>')
    
  end subroutine solve

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-03-10
  !> @brief   Set coarse matrix
  !> @details Set coarse matrix
  !> 
  !-----------------------------------------------------------------------
  
  subroutine set_coarse(self,a)

    class(dcg),               intent(inout) :: self              !< Solver
    class(mat),               intent(in)    :: a                 !< Matrix A
    integer(ip)                             :: iz,i,j
    integer(ip)                             :: igrou,jgrou
    integer(ip)                             :: izgro,ndof,kskyl
    
    if( self % ngrou > 0 ) then
       select type ( ac => self % acoarse )
          
       class is ( mat_csr )
          !
          ! Matrix in CSR format
          !
          select type ( a )
          class is ( mat_csr )
             ndof = a % ndof1
             do i = 1,a % nrows
                if( self % lgrou(i) > 0 ) then
                   igrou = self % lgrou(i)
                   do iz = a % ia(i),a % ia(i+1)-1
                      j = a % ja(iz)
                      if( self % lgrou(j) > 0 ) then
                         jgrou = self % lgrou(j)
                         izgro = ac % iA(igrou)
                         iifzgro1: do while( ac% jA(izgro) /= jgrou )
                            izgro = izgro + 1
                         end do iifzgro1
                         ac % vA(1:ndof,1:ndof,izgro) = ac % vA(1:ndof,1:ndof,izgro) &
                              + a % vA(1:ndof,1:ndof,iz)
                      end if
                   end do
                end if
             end do
          end select
          call self % parallel_sum(ac % vA)
          
       class is ( mat_sky )
          !
          ! Matrix in skyline format
          !
          select type ( a )
          class is ( mat_csr )
             do i = 1,a % nrows
                if( self % lgrou(i) > 0 ) then
                   igrou = self % lgrou(i)
                   do iz = a % ia(i),a % ia(i+1)-1
                      j = a % ja(iz)
                      if( j < i ) then
                         if( self % lgrou(j) > 0 ) then
                            jgrou = self % lgrou(j)
                            if( igrou > jgrou ) then
                               kskyl          = ac % iskyl(igrou+1) - 1 - (igrou-jgrou)
                               ac % vA(kskyl) = ac % vA(kskyl) + a % vA(1,1,iz)
                            else if( igrou < jgrou ) then
                               kskyl          = ac % iskyl(jgrou+1) - 1 - (jgrou-igrou)
                               ac % vA(kskyl) = ac % vA(kskyl) + a % vA(1,1,iz)
                            else
                               kskyl          = ac % iskyl(igrou+1) - 1
                               ac % vA(kskyl) = ac % vA(kskyl) + 2.0_rp * a % vA(1,1,iz)
                            end if
                         end if
                      else if( i == j ) then
                         kskyl          = ac % iskyl(igrou+1) - 1
                         ac % vA(kskyl) = ac % vA(kskyl) + a % vA(1,1,iz)  
                      end if

                   end do
                end if
             end do
          end select
          call self % parallel_sum(ac % vA)
       end select
    end if

  end subroutine set_coarse

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-03-10
  !> @brief   Pass a vector from fine to coarse
  !> @details Pass a vector from fine to coarse
  !>          XC = W^T.XF 
  !> 
  !-----------------------------------------------------------------------

  subroutine fine_to_coarse(self,xf,xc)

    class(dcg),          intent(inout) :: self              !< Solver
    real(rp),   pointer, intent(inout) :: xf(:)
    real(rp),   pointer, intent(inout) :: xc(:)
    integer(ip)                        :: i,igrou,kgrou,k,idof

    if( self % ngrou > 0 ) then

       do igrou = 1,self % ngrou
          xc(igrou) = 0.0_rp
       end do

       if( self % ndof == 1 ) then
          do i = 1,self % nown
             igrou = self % lgrou(i)
             if( igrou > 0 ) xc(igrou) = xc(igrou)+xf(i)
          end do
       else
          do i = 1,self % nown
             igrou  = self % lgrou(i)
             if( igrou > 0 ) then 
                do idof = 1,self % ndof
                   kgrou     = (igrou-1) * self % ndof + idof
                   k         = (i-1)     * self % ndof + idof
                   xc(kgrou) = xc(kgrou) + xf(k)
                end do
             end if
          end do
       end if

       call self % parallel_sum(xc)

    end if
    
  end subroutine fine_to_coarse
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-03-10
  !> @brief   Pass a vector from coarse to fine
  !> @details Pass a vector from coarse to fine
  !>          XF = W.XC 
  !> 
  !-----------------------------------------------------------------------

  subroutine coarse_to_fine(self,xc,xf)

    class(dcg),          intent(inout) :: self              !< Solver
    real(rp),   pointer, intent(inout) :: xc(:)
    real(rp),   pointer, intent(inout) :: xf(:)
    integer(ip)                        :: i,igrou,kgrou,k,idof

    if( self % n > 0 .and. self % ngrou > 0 ) then
       
       if( self % ndof == 1 ) then

          do i = 1,self % n
             igrou = self % lgrou(i)
             if( igrou > 0 ) then
                xf(i) = xc(igrou)
             else
                xf(i) = 0.0_rp
             end if
          end do

       else

          do i = 1,self % n
             igrou = self % lgrou(i)
             if( igrou > 0 ) then
                do idof = 1,self % ndof 
                   kgrou = (igrou-1) * self % ndof + idof
                   k     = (i-1)     * self % ndof + idof
                   xf(k) = xc(kgrou)
                end do
             else
                do idof = 1,self % ndof 
                   k     = (i-1)     * self % ndof + idof
                   xf(k) = 0.0_rp
                end do
             end if
          end do

       end if
       
    end if

  end subroutine coarse_to_fine

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-03-10
  !> @brief   Coarse correction
  !> @details Compute a coarse correction
  !>          w <= A' mu = W A'^-1 W^T w
  !> 
  !-----------------------------------------------------------------------
  
  subroutine coarse_correction(self,mu,murhs,w)
    
    class(dcg),                    intent(inout) :: self              !< Solver
    real(rp),             pointer, intent(inout) :: mu(:)
    real(rp),             pointer, intent(inout) :: murhs(:)
    real(rp),             pointer, intent(inout) :: w(:)
    integer(ip)                                  :: i
    
    if( self % ngrou > 0 ) then
       call self % fine_to_coarse (w,murhs)                 ! W^T w
       call self % dir % solve    (murhs,mu,self % acoarse) ! A' mu = murhs
       call self % coarse_to_fine (mu,w)                    ! w <= mu
    else
       do i = 1,self % nn
          w(i) = 0.0_rp
       end do
    end if

  end subroutine coarse_correction
  
end module def_dcg
!> @}
