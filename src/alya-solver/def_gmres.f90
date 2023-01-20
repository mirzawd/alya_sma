!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!




!-----------------------------------------------------------------------
!> @addtogroup Solvers
!> @{
!> @file    def_gmres.f90
!> @author  houzeaux
!> @date    2022-03-01
!> @brief   Conjugate gradient
!> @details Conjugate gradient
!>         
!-----------------------------------------------------------------------

module def_gmres

  use def_kintyp,            only : ip,rp,lg
  use def_mat,               only : mat
  use def_iterative_solvers, only : krylov_unsymmetric
  use mod_memory_basic,      only : memory_alloca
  use mod_memory_basic,      only : memory_deallo
  implicit none
  private
  
  type, extends(krylov_unsymmetric) :: gmres
   contains
     procedure,         pass :: init 
     procedure,         pass :: alloca
     procedure,         pass :: deallo
     procedure,         pass :: solve
     procedure,         pass :: setup
  end type gmres

  character(8), parameter :: vacal = 'gmres'
  real(rp),     parameter :: epsmac=1.0e-16_rp

  public :: gmres
  
contains

  subroutine init(self)
    
    class(gmres), intent(inout) :: self

  end subroutine init  
  
  subroutine deallo(self)
    
    class(gmres), intent(inout) :: self
        
  end subroutine deallo
  
  subroutine alloca(self)

    class(gmres),             intent(inout) :: self          !< Solver
    
  end subroutine alloca
 
  subroutine setup(self,a)

    class(gmres),                       intent(inout) :: self          !< Solver
    class(mat),                         intent(inout) :: a             !< Matrix A

  end subroutine setup
 
  subroutine solve(self,b,x,a)
    class(gmres),             intent(inout) :: self              !< Solver
    real(rp),        pointer, intent(in)    :: b(:)              !< RHS b
    real(rp),        pointer, intent(inout) :: x(:)              !< Solve x
    class(mat),               intent(in)    :: a                 !< Matrix A
    integer(ip)                             :: n,nn
    integer(ip)                             :: ndof,kryldim
    integer(ip)                             :: ierro,i
    integer(ip)                             :: ii,jj,kk,jj1,idx
    real(rp)                                :: resid,gama
    real(rp)                                :: time1,time2
    real(rp),        pointer                :: bp(:)
    real(rp),        pointer                :: kryl(:,:)
    real(rp),        pointer                :: cc(:)
    real(rp),        pointer                :: ss(:)
    real(rp),        pointer                :: rs(:)
    real(rp),        pointer                :: hh(:)
    real(rp),        pointer                :: r(:)
    real(rp),        pointer                :: kp(:)
    real(rp),        pointer                :: kq(:)
    logical(lg)                             :: convergence,fin2

    nullify(bp,kryl,cc,ss,rs,hh,r)

    ndof    =  self % ndof
    n       =  self % n
    nn      =  self % nn
    kryldim =  self % input % nkryd  

    call memory_alloca(self % memor,'CC'  ,'gmres',cc  ,kryldim)
    call memory_alloca(self % memor,'SS'  ,'gmres',ss  ,kryldim)
    call memory_alloca(self % memor,'RS'  ,'gmres',rs  ,kryldim+1_ip)
    call memory_alloca(self % memor,'HH'  ,'gmres',hh  ,kryldim+((kryldim*(kryldim+1_ip))/2_ip),'DO_NOT_INITIALIZE')
    call memory_alloca(self % memor,'BP'  ,'gmres',bp  ,nn)
    call memory_alloca(self % memor,'R'   ,'gmres',r   ,nn)
    call memory_alloca(self % memor,'KRYL','gmres',kryl,max(1_ip,nn),kryldim+1_ip,'DO_NOT_INITIALIZE')
        
    kp => kryl(:,1)
    call self % start (b,x,a,r,kp,ierro)    
    if( ierro /= 0 ) goto 1

    call self % preconditioning(b,bp)                            ! bp = L^-1 b

    resid = self % parallel_dot(kp,kp)

    if( resid <= 0.0_rp ) goto 1
    
    resid = sqrt(resid)    
    call self % init_iterations(resid)                           ! Residuals and tolerance

    convergence                 = .false.
    fin2                        = .false.

    do while( .not. convergence ) 
       !
       ! Initial residual: kryl(1) = L^-1 ( b - A x ) = bp - L^-1 A x
       ! L kryl(1) = A R^-1 x
       !
       kq => kryl(:,1)
       call self % preconditioning(x,kq,r,a)

       do i = 1,nn
          kryl(i,1) = bp(i) - kryl(i,1)         
       end do       
       kp                          => kryl(:,1)
       resid                       =  self % parallel_L2norm(kp)
       self % output % resip_old   =  self % output % resip_final
       self % output % resip_final =  resid * self % output % bnorp_inv

       call self % outcvg(b,x,a)

       if( resid <= self % output % toler ) then
          !
          ! The initial guess is the solution
          !
          convergence = .true.
       else
          !
          ! Initialize 1-st term of the rhs of hessenberg system
          !
          rs(1) = resid
          !
          ! Ortonormalize kryl(*,1)
          !
          resid = 1.0_rp / resid
          do i = 1,nn
             kryl(i,1) = kryl(i,1) * resid
          end do
       end if

       jj   = 0
       idx  = -1
       fin2 = convergence
       !
       ! Inner loop. Restarted each kryldim iterations
       !
       do while( .not. fin2 ) 

          self % output % iters = self % output % iters + 1
          jj                    = jj + 1
          jj1                   = jj + 1
          idx                   = idx + jj
          !
          ! L kryl(jj1) = A R^-1 kryl(jj)
          !
          kp => kryl(:,jj)
          kq => kryl(:,jj1)
          call self % preconditioning(kp,kq,r,a)
          !
          ! Orthogonalization: Modified Gram-Schmidt
          ! For i= 1, j
          !     H(i,j) = <v_i, v_j1>
          !       v_j1 = v_j1 - H(i,j) * v_i
          !
          call cputim(time1)
          do ii = 1,jj
             kp    => kryl(:,ii)
             kq    => kryl(:,jj1)
             resid = self % parallel_dot(kp,kq)
             hh(idx+ii) = resid
             do kk = 1,nn
                kryl(kk,jj1) = kryl(kk,jj1) - resid * kryl(kk,ii)
             end do
          end do
          call cputim(time2)
          self % output % cputi(5) = self % output % cputi(5) + time2 - time1
          !
          ! H(jj1,jj) = ||kryl(*,jj1)||
          !
          kp          => kryl(:,jj1)
          resid       = self % parallel_L2norm(kp)
          hh(idx+jj1) = resid

          if( resid == 0.0_rp ) then
             fin2                  = .true.
             convergence           = .true.
             idx                   = idx - jj
             jj                    = jj - 1
             self % output % iters = self % output % iters - 1
             if(self % input % lun_solve /= 0 ) &
                  write(self % input % lun_solve,*) '||kryl(*,jj1)|| = 0 ',jj1,self % output % iters
             goto 2
          else
             !
             ! Ortonormalize kryl(*,jj1)
             !
             resid = 1.0_rp / resid
             do kk = 1,nn
                kryl(kk,jj1) = kryl(kk,jj1) * resid
             end do
          end if
          !
          ! Update factorization of H. Perform previous
          ! transformations on jj-th column of H
          !
          do ii = 1,jj-1
             kk         =   ii + 1
             resid      =   hh(idx+ii)
             hh(idx+ii) =   cc(ii) * resid + ss(ii) * hh(idx+kk)
             hh(idx+kk) = - ss(ii) * resid + cc(ii) * hh(idx+kk)
          end do
          gama = hh(idx+jj)*hh(idx+jj) + hh(idx+jj1)*hh(idx+jj1)
          gama = sqrt(gama)
          !
          ! if gama is zero then take any small
          ! value will affect only residual estimate
          !
          if( gama == 0.0_rp ) then
             gama = epsmac
             if( self % input % lun_solve /= 0 ) &
                  write(self % input % lun_solve,*) 'gama==0.0 ',self % output % iters
          end if
          !
          ! Get next plane rotation
          !
          gama       =   1.0_rp / gama
          cc(jj)     =   hh(idx+jj)  * gama
          ss(jj)     =   hh(idx+jj1) * gama
          hh(idx+jj) =   cc(jj) * hh(idx+jj) + ss(jj) * hh(idx+jj1)
          !
          ! Update the rhs of the LS problem
          !
          rs(jj1)    = - ss(jj) * rs(jj)
          rs(jj)     =   cc(jj) * rs(jj)
          !
          ! Convergence Test
          !
          resid      =   abs( rs(jj1) )

          self % output % resip_old   = self % output % resip_final
          self % output % resip_final = resid * self % output % bnorp_inv

          call self % outcvg(b,x,a)

          if( resid <= self % output % toler ) then
             convergence = .true.
             fin2        = .true.
          else
             if( self % output % iters >= self % input % miter ) then
                convergence = .true.
                fin2        = .true.
             else
                if( jj >= kryldim ) then
                   fin2 = .true.
                end if
             end if
          end if

2        continue
       end do

       !----------------------------------------------------------------------
       !
       ! END INNER LOOP
       !
       !----------------------------------------------------------------------
       !
       ! Compute y => Solve upper triangular system
       !
       do ii= jj,2,-1
          rs(ii) = rs(ii) / hh(idx+ii)
          resid  = rs(ii)

          do kk = 1,ii-1
             rs(kk) = rs(kk) - hh(idx+kk) * resid
          end do

          idx = idx - ii
       end do

       if( hh(1) /= 0.0_rp ) rs(1) = rs(1) / hh(1)
       !
       ! Linear combination of kryl(*,jj)'s to get the solution.
       !
       do ii = 1,jj
          resid = rs(ii)
          do kk = 1,nn
             x(kk) = x(kk) + resid * kryl(kk,ii)
          end do
       end do

    end do

1   continue
    
    call self % end(b,x,a)

    call memory_deallo(self % memor,'CC'  ,'gmres',cc  )
    call memory_deallo(self % memor,'SS'  ,'gmres',ss  )
    call memory_deallo(self % memor,'RS'  ,'gmres',rs  )
    call memory_deallo(self % memor,'HH'  ,'gmres',hh  )
    call memory_deallo(self % memor,'BP'  ,'gmres',bp  )
    call memory_deallo(self % memor,'R'   ,'gmres',r   )
    call memory_deallo(self % memor,'KRYL','gmres',kryl)

  end subroutine solve

end module def_gmres
!> @}

