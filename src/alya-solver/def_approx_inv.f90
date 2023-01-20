!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!




!-----------------------------------------------------------------------
!> @addtogroup Solvers
!> @{
!> @file    def_approx_inv.f90
!> @author  houzeaux
!> @date    2022-03-01
!> @brief   Approx LU given a sparsity pattern
!> @details A = L U or A = L L^t
!>         
!-----------------------------------------------------------------------

module def_approx_inv

  use def_kintyp,            only : ip,rp,lg,i1p
  use def_mat,               only : mat
  use def_mat_den,           only : mat_den
  use def_mat_csr,           only : mat_csr
  use def_mat_dia,           only : mat_dia
  use def_direct_solvers,    only : direct_solver
  use mod_memory_basic,      only : memory_alloca
  use mod_memory_basic,      only : memory_deallo
  use mod_memory_basic,      only : memory_size
  use mod_maths_arrays,      only : maths_findloc
  use mod_strings
  implicit none
  private
  
  type, extends(direct_solver) :: approx_inv
     type(mat_csr)           :: a
     type(mat_dia)           :: dia
   contains
     procedure,         pass :: init 
     procedure,         pass :: alloca
     procedure,         pass :: deallo
     procedure,         pass :: solve
     procedure,         pass :: setup
     procedure,         pass :: insert_ordered
     procedure,         pass :: identity
  end type approx_inv

  character(10), parameter :: vacal = 'approx_inv'

  public :: approx_inv
  
contains
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-03-12
  !> @brief   Initialization
  !> @details Initialization
  !> 
  !-----------------------------------------------------------------------

  subroutine init(self)
    
    class(approx_inv), intent(inout) :: self
    
    call self % a   % init()
    call self % dia % init()

  end subroutine init  
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-03-12
  !> @brief   Deallocate
  !> @details Deallocate
  !> 
  !-----------------------------------------------------------------------

  subroutine deallo(self)
    
    class(approx_inv), intent(inout) :: self

    call self % a   % deallo()
    call self % dia % deallo()
    
  end subroutine deallo
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-03-12
  !> @brief   Alllocate
  !> @details Alllocate
  !> 
  !-----------------------------------------------------------------------

  subroutine alloca(self)

    class(approx_inv),                intent(inout) :: self          !< Solver
    
  end subroutine alloca
 
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-03-12
  !> @brief   Set up solver
  !> @details Set up solver. Graph is assmued to be ordered
  !>          in increasing order
  !> 
  !-----------------------------------------------------------------------

  subroutine setup(self,a)

    class(approx_inv),                  intent(inout) :: self          !< Solver
    class(mat),                         intent(inout) :: a             !< Matrix A
    type(mat_den)                                     :: matden
    type(mat_dia)                                     :: dial
    integer(ip)                                       :: ii,iz,jj,ir,nz,istat
    integer(ip)                                       :: n,nn,nr,jz,maxnr
    integer(ip)                                       :: kk,li,i1,i2,pos,ll,kz
    integer(ip)                                       :: xx,yy,ndof,nend
    integer(ip)                                       :: idof,jdof
    type(i1p),    pointer                             :: il(:)
    type(i1p),    pointer                             :: il_tmp(:)
    real(rp),     pointer                             :: rhs(:)
    real(rp),     pointer                             :: sol(:)

    if( self % n > 0 ) then

       nullify(il,il_tmp,rhs,sol)
       
       select type ( a )
       class is ( mat_csr )

          select case   ( self % input % kfl_full_rows )
          case ( 0_ip ) ; nend = self % nint
          case default  ; nend = self % n
          end select
          ndof = a % ndof1 
          n    = nend
          nn   = nend * ndof
          
          !----------------------------------------------------------------
          !
          ! Diagonal matrix
          !
          !----------------------------------------------------------------

          call a    % diag             (self % dia)
          call self % parallel_exchange(self % dia)
          call self % dia % inverse    ()
          
          !----------------------------------------------------------------
          !
          ! Extended graph IL. For starting graph, consider only
          ! interior nodes when using partial row
          !
          !----------------------------------------------------------------

          call memory_alloca(self % memor,'IL'    ,vacal,il    ,nend)
          call memory_alloca(self % memor,'IL_TMP',vacal,il_tmp,nend)

          do ii = 1,nend
             nr = count(a % ja(a % ia(ii):a % ia(ii+1)-1)<=nend)
             call memory_alloca(self % memor,'IL % L',vacal,il(ii) % l,nr,'DO_NOT_INITIALIZE')
             nr = 0
             do iz = a % ia(ii),a % ia(ii+1)-1
                if( a % ja(iz) <= nend ) then
                   nr = nr + 1 ; il(ii) % l(nr) = a % ja(iz)
                end if
             end do
          end do
          
          do li = 1,self % input % levels_sparsity-1
             do ii = 1,nend
                nr = memory_size(il(ii) % l)
                call memory_alloca(self % memor,'IL_TMP % Ls',vacal,il_tmp(ii) % l,nr,'DO_NOT_INITIALIZE')
                il_tmp(ii) % l(1:nr) = il(ii) % l(1:nr) 
             end do
             do ii = 1,nend
                do i1 = 1,memory_size(il_tmp(ii)%l)
                   kk = il_tmp(ii) % l(i1)          
                   do i2 = 1,memory_size(il_tmp(kk)%l)
                      jj = il_tmp(kk) % l(i2)
                      call self % insert_ordered(jj,il(ii)%l,pos)
                   end do
                end do
             end do
             do ii = 1,nend
                call memory_deallo(self % memor,'IL_TMP % L',vacal,il_tmp(ii) % l)
             end do
          end do
          call memory_deallo(self % memor,'IL_TMP',vacal,il_tmp)
          !
          ! Allocate matrix and RHS
          !
          maxnr = 0
          do ii = 1,nend
             maxnr = max(maxnr,memory_size(il(ii)%l))
          end do

          !----------------------------------------------------------------
          !
          ! Copy graph A % IA, A % JA <= IL
          !
          !----------------------------------------------------------------

          call self % a % init()    
          nz = 0
          if( self % input % kfl_symm == 1 ) then
             do ii = 1,nend
                iz0: do kk = 1,memory_size(il(ii)%l)
                   if( il(ii)%l(kk) == ii ) exit iz0
                end do iz0
                nz = nz + kk
             end do
          else
             do ii = 1,nend
                nz = nz + memory_size(il(ii)%l)
             end do
          end if

          self % a % nz    = nz
          self % a % nrows = nend
          self % a % ncols = nend
          self % a % ndof1 = ndof
          self % a % ndof2 = ndof
          call self % a % alloca()
          self % a % ia(1) = 1

          if( self % input % kfl_symm == 1 ) then
             do ii = 1,nend
                kk = 0
                iz1: do kk = 1,memory_size(il(ii)%l)
                   if( il(ii)%l(kk) == ii ) exit iz1     
                end do iz1
                self % a % ia(ii+1) = self % a % ia(ii) + kk
                kk = 0
                iz2: do iz = self % a % ia(ii),self % a % ia(ii+1)-1
                   kk = kk + 1
                   self % a % ja(iz) = il(ii)%l(kk)
                   if( self % a % ja(iz) == ii ) exit iz2
                end do iz2
             end do
          else
             do ii = 1,nend
                self % a % ia(ii+1) = self % a % ia(ii) + memory_size(il(ii)%l)
                kk = 0
                do iz = self % a % ia(ii),self % a % ia(ii+1)-1
                   kk = kk + 1
                   self % a % ja(iz) = il(ii)%l(kk)
                end do
             end do
          end if
          !
          ! Allocate and deallocate
          !
          call memory_alloca(self % memor,'RHS'   ,vacal,rhs,maxnr)
          call memory_alloca(self % memor,'SOL'   ,vacal,sol,maxnr)
          call memory_deallo(self % memor,'IL'    ,vacal,il)

          !----------------------------------------------------------------
          !
          ! Fill in preconditioner matrix SELF % A
          !
          !----------------------------------------------------------------

          do ii = 1,nend
             do idof = 1,ndof
                nr = (self % a % ia(ii+1) - self % a % ia(ii)) * ndof
                call matden % init()
                matden % nrows = nr * ndof
                matden % ncols = nr * ndof
                matden % ndof1 = 1
                matden % ndof2 = 1
                call matden % alloca()
                !
                ! Densify A => SELF % A % VA
                !
                xx = 0
                do iz = self % a % ia(ii),self % a % ia(ii+1)-1
                   kk = self % a % ja(iz)
                   xx = xx + 1
                   yy = 0
                   do jz = self % a % ia(ii),self % a % ia(ii+1)-1
                      ll = self % a % ja(jz)
                      yy = yy + 1
                      do kz = a % ia(kk),a % ia(kk+1)-1
                         if( a % ja(kz) == ll ) then
                            do jdof = 1,ndof
                               matden % vA(1,1,(xx-1)*ndof+idof,(yy-1)*ndof+idof) &
                                    = a % va(jdof,idof,kz)
                               if( ll == ii .and. jdof == idof ) ir = (yy-1)*ndof+idof
                            end do
                         end if
                      end do
                   end do
                end do
                !
                ! Output dense matrices
                !
                !call matden % output('DENSE','mat'//integer_to_string(ii)//'.txt')
                !call a   % output('DENSE','original.txt')
                !
                ! RHS
                !
                rhs(1:nr) = 0.0_rp
                select case   ( self % input % kfl_symm )
                case ( 1_ip ) ; rhs(nr) = 1.0_rp
                case default  ; rhs(ir) = 1.0_rp
                end select
                !
                ! Solve dense system MAT SOL = RHS
                !
                call matden % solve(rhs,sol,istat) !,PIVOTING=.true.)

                if( istat == 1 ) then
                   !
                   ! If not invertible, take diagonal
                   !
                   ir = 0
                   do iz = self % a % ia(ii),self % a % ia(ii+1)-1
                      do jdof = 1,ndof
                         ir = ir + 1
                         if( self % a % ja(iz) == ii ) self % a % vA(idof,idof,iz) &
                              = self % dia % vA(idof,ii)
                      end do
                   end do
                else
                   !
                   ! Substitute row of II by solution SOL of system
                   !
                   ir = 0
                   do iz = self % a % ia(ii),self % a % ia(ii+1)-1
                      do jdof = 1,ndof
                         ir = ir + 1
                         self % a % vA(jdof,idof,iz) = sol(ir)
                      end do
                   end do
                end if

                call matden % deallo()

             end do

          end do
          
          !----------------------------------------------------------------
          !
          ! Scale for symmetric case A <= A / sqrt(diag(A))
          !
          !----------------------------------------------------------------

          if( self % input % kfl_symm == 1 ) then
             call dial % init      ()             
             call self % a % diag  (dial)
             call dial % inverse   ()
             call self % a % scale (dial,SQUARE_ROOT=.true.)
              call dial % deallo    ()
          end if

          !----------------------------------------------------------------
          !
          ! Deallocate
          !
          !----------------------------------------------------------------

          call memory_deallo(self % memor,'RHS',vacal,rhs)
          call memory_deallo(self % memor,'SOL',vacal,sol)

          !----------------------------------------------------------------
          !
          ! Check identity
          !
          !----------------------------------------------------------------

          !call identity(self,a)

       end select

    end if

  end subroutine setup
 
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-03-12
  !> @brief   Solve
  !> @details Solve x = A'^-1 b ... unsymmetric case
  !>                x = L L^t b ... symmetric case
  !> 
  !-----------------------------------------------------------------------

  subroutine solve(self,b,x,a)
    
    class(approx_inv),        intent(inout) :: self       !< Solver
    real(rp),        pointer, intent(in)    :: b(:)       !< RHS b
    real(rp),        pointer, intent(inout) :: x(:)       !< Solve x
    class(mat),               intent(in)    :: a          !< Matrix A
    real(rp),        pointer                :: w(:)
    integer(ip)                             :: i,idof,nend
    !
    ! Interface nodes, use diagonal
    !
    if( self % input % kfl_full_rows == 0 ) then
       do i = self % nint+1,self % n
          do idof = 1,self % ndof
             x((i-1)*self % ndof+idof) =  b((i-1)*self % ndof+idof) * self % dia % vA(idof,i)
          end do
       end do
    end if

    block
      real(rp), pointer :: y(:)
      nullify(y)
      nend = self % nint
      call memory_alloca(self % memor,'W',vacal,y,nend)
      if(self % nint > 0 ) y = 1.0_rp
      if( self % input % kfl_symm == 1 ) then
         nullify(w)
         call memory_alloca(self % memor,'W',vacal,w,nend)
         call self % a % mv(y,w,n1=1_ip,n2=nend)
         call self % a % mv(w,x,n1=1_ip,n2=nend,TRANSPOSE=.true.)
         call memory_deallo(self % memor,'W',vacal,w)
      else
         call self % a % mv(y,x,n1=1_ip,n2=nend)
      end if
    end block
    !
    ! Use mv for interior nodes only
    !
    select case   ( self % input % kfl_full_rows )
    case ( 0_ip ) ; nend = self % nint
    case default  ; nend = self % n
    end select
    
    if( self % input % kfl_symm == 1 ) then
       nullify(w)
       call memory_alloca(self % memor,'W',vacal,w,nend)
       call self % a % mv(b,w,n1=1_ip,n2=nend)
       call self % a % mv(w,x,n1=1_ip,n2=nend,TRANSPOSE=.true.)
       call memory_deallo(self % memor,'W',vacal,w)
    else
       call self % a % mv(b,x,n1=1_ip,n2=nend)
    end if
    
  end subroutine solve
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-03-12
  !> @brief   Insert
  !> @details Insert a new point in increasing order
  !>
  !-----------------------------------------------------------------------

  subroutine insert_ordered(self,i,l,pos)
    
    class(approx_inv),               intent(inout) :: self       !< Solver
    integer(ip),                    intent(in)    :: i
    integer(ip),           pointer, intent(inout) :: l(:)
    integer(ip), optional,          intent(out)   :: pos
    integer(ip),           pointer                :: l_tmp(:)
    integer(ip)                                   :: j,n,iloc
    logical(lg)                                   :: found

    n  = memory_size(l)
    if( present(pos) ) pos = 0
    
    if( n > 0 ) then
       iloc = maths_findloc(l,i) 
       if( iloc == 0 ) then
          j     = 0
          n     = size(l)
          found = .false.
          do while( j < n )
             j = j + 1
             if( l(j) > i ) then
                found = .true.
                exit
             end if
          end do
          if( .not. found  ) j = n+1
          if( present(pos) ) pos = j
          allocate(l_tmp(n+1))
          l_tmp(1:j-1)   = l(1:j-1)
          l_tmp(j)       = i
          l_tmp(j+1:n+1) = l(j:n)
          deallocate(l)
          l => l_tmp
       end if
    else
       call memory_alloca(self % memor,'IL % L',vacal,l,1_ip)
       l(1) = i
    end if

  end subroutine insert_ordered

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-03-12
  !> @brief   Check identity
  !> @details Check identity
  !> 
  !-----------------------------------------------------------------------

  subroutine identity(self,a)

    class(approx_inv), intent(inout) :: self          !< Solver
    class(mat),       intent(inout) :: a             !< Matrix A
    type(mat_den)                   :: ad
    type(mat_den)                   :: ainv
    type(mat_den)                   :: id
    integer(ip)                     :: ii,jj,nn

    nn = a % nrows * a % ndof1
    select type ( a )
    class is    ( mat_csr )
       call ad   % init()
       call ainv % init()
       call id   % init()
       call ad   % csr2den(a)
       call ainv % csr2den(self % a)
       id % nrows = nn
       id % ncols = nn
       id % ndof1 = 1
       id % ndof2 = 1
       call id % alloca()
       call ad % output('DENSE','a.txt')
       do ii = 1,nn
          do jj = 1,nn
             id % vA(1,1,ii,jj) = dot_product(ad%va(1,1,ii,1:nn),ainv%va(1,1,1:nn,jj))
             if( ii == jj .and. abs(id % vA(1,1,ii,jj)-1.0_rp) > 1.0e-12_rp ) stop 1
             if( ii /= jj .and. abs(id % vA(1,1,ii,jj)-0.0_rp) > 1.0e-12_rp ) stop 1
          end do
       end do
       call id % output('DENSE','id.txt')
    end select

  end subroutine identity
  
end module def_approx_inv
!> @}

