!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!
!> @defgroup Maths_Toolbox
!> Toolbox for mathematical solvers
!> @{
!> @name    ToolBox for mathematics operations
!> @file    mod_maths.f90
!> @author  Guillaume Houzeaux
!> @brief   ToolBox for algebraic solvers
!> @details ToolBox for algebraic solvers
!
!-----------------------------------------------------------------------

module mod_maths_solver

  use def_maths
  use mod_maths_basic
  implicit none 

  private
  
  public :: maths_conjugate_gradient            ! CG solver
  public :: maths_direct                        ! Direct solver
  public :: maths_invert_matrix                 ! Matrix inversion
  public :: maths_schur_complement              ! Schur complement
  public :: maths_eigen_symmetric_matrix        ! Bridgr to 2D and 3D
  public :: maths_eigen_3x3_symmetric_matrix    ! Eigenvalues and vectors for real symmetric 3x3 system
  public :: maths_eigen_2x2_symmetric_matrix    ! Eigenvalues and vectors for real symmetric 2x2 system
  public :: maths_solve_overdetermined_system   ! Overdetermined solver
  public :: maths_eigen_jacobi                  ! Jacobi method to solve eigenvalue problem                  

contains
  
  !-----------------------------------------------------------------------
  !
  !> @date    27/06/2014
  !> @author  Guillaume Houzeaux
  !> @brief   Inverse a matrix
  !> @details Invert
  !
  !-----------------------------------------------------------------------

  pure subroutine maths_invert_matrix(nsize,a,deter,inva,ACCURACY)
    
    integer(ip), intent(in)             :: nsize
    real(rp),    intent(inout)          :: a(nsize,nsize)
    real(rp),    intent(out),  optional :: inva(nsize,nsize)
    real(rp),    intent(out),  optional :: deter
    logical(lg), intent(in),   optional :: accuracy
    logical(lg)                         :: accurate_method
    real(rp)                            :: denom,t1,t2,t3,t4,det
    real(rp)                            :: tt(6),tp,tn
    real(rp),    allocatable            :: b(:,:)
    integer(ip)                         :: ii,jj,n

    accurate_method = .false.
    if( present(ACCURACY) ) then
       accurate_method = ACCURACY
    end if

    if( present(inva) ) then

       select case ( nsize )

       case ( 1_ip )

          det = a(1,1)
          if( abs(det) == 0.0_rp ) goto 10
          inva(1,1) = 1.0_rp/a(1,1)

       case( 2_ip )

          det = a(1,1)*a(2,2)-a(2,1)*a(1,2)
          if( abs(det) == 0.0_rp ) goto 10
          denom  = 1.0_rp/det
          inva(1,1) = a(2,2)*denom
          inva(2,2) = a(1,1)*denom
          inva(2,1) =-a(2,1)*denom
          inva(1,2) =-a(1,2)*denom

       case ( 3_ip )

          if( accurate_method ) then
             tt(1) =  a(1,1)*a(2,2)*a(3,3)
             tt(2) = -a(1,1)*a(3,2)*a(2,3)
             tt(3) = -a(1,2)*a(2,1)*a(3,3)
             tt(4) =  a(1,2)*a(3,1)*a(2,3)
             tt(5) =  a(1,3)*a(2,1)*a(3,2)
             tt(6) = -a(1,3)*a(3,1)*a(2,2)
             tp    =  sum(tt,mask=tt >= 0.0_rp)
             tn    =  sum(tt,mask=tt <  0.0_rp)
             det   =  tp-tn
          else
             t1    = a(2,2)*a(3,3) - a(3,2)*a(2,3)
             t2    =-a(2,1)*a(3,3) + a(3,1)*a(2,3)
             t3    = a(2,1)*a(3,2) - a(3,1)*a(2,2)
             det   = a(1,1)*t1 + a(1,2)*t2 + a(1,3)*t3
          end if

          if( abs(det) == 0.0_rp ) goto 10
          denom     = 1.0_rp/det
          inva(1,1) = t1*denom
          inva(2,1) = t2*denom
          inva(3,1) = t3*denom
          inva(2,2) = ( a(1,1)*a(3,3) - a(3,1)*a(1,3))*denom
          inva(3,2) = (-a(1,1)*a(3,2) + a(1,2)*a(3,1))*denom
          inva(3,3) = ( a(1,1)*a(2,2) - a(2,1)*a(1,2))*denom
          inva(1,2) = (-a(1,2)*a(3,3) + a(3,2)*a(1,3))*denom
          inva(1,3) = ( a(1,2)*a(2,3) - a(2,2)*a(1,3))*denom
          inva(2,3) = (-a(1,1)*a(2,3) + a(2,1)*a(1,3))*denom

       case ( 4_ip )

          t1=    a(2,2)*a(3,3)*a(4,4) + a(2,3)*a(3,4)*a(4,2)&
               + a(2,4)*a(3,2)*a(4,3) - a(2,3)*a(3,2)*a(4,4)&
               - a(2,2)*a(3,4)*a(4,3) - a(2,4)*a(3,3)*a(4,2)
          t2=  - a(2,1)*a(3,3)*a(4,4) - a(2,3)*a(3,4)*a(4,1)&
               - a(2,4)*a(3,1)*a(4,3) + a(2,4)*a(3,3)*a(4,1)&
               + a(2,3)*a(3,1)*a(4,4) + a(2,1)*a(3,4)*a(4,3)
          t3=    a(2,1)*a(3,2)*a(4,4) + a(2,2)*a(3,4)*a(4,1)&
               + a(2,4)*a(3,1)*a(4,2) - a(2,4)*a(3,2)*a(4,1)&
               - a(2,2)*a(3,1)*a(4,4) - a(2,1)*a(3,4)*a(4,2)
          t4=  - a(2,1)*a(3,2)*a(4,3) - a(2,2)*a(3,3)*a(4,1)&
               - a(2,3)*a(3,1)*a(4,2) + a(2,3)*a(3,2)*a(4,1)&
               + a(2,2)*a(3,1)*a(4,3) + a(2,1)*a(3,3)*a(4,2)
          det= a(1,1)*t1 + a(1,2)*t2 + a(1,3)*t3 + a(1,4)*t4
          if( abs(det) == 0.0_rp ) goto 10
          denom     = 1.0_rp/det
          inva(1,1) = t1*denom
          inva(2,1) = t2*denom
          inva(3,1) = t3*denom
          inva(4,1) = t4*denom
          inva(1,2) =(- a(1,2)*a(3,3)*a(4,4) - a(1,3)*a(3,4)*a(4,2)&
               &      - a(1,4)*a(3,2)*a(4,3) + a(1,3)*a(3,2)*a(4,4)&
               &      + a(1,2)*a(3,4)*a(4,3) + a(1,4)*a(3,3)*a(4,2))*denom
          inva(2,2) =(  a(1,1)*a(3,3)*a(4,4) + a(1,3)*a(3,4)*a(4,1)&
               &      + a(1,4)*a(3,1)*a(4,3) - a(1,4)*a(3,3)*a(4,1)&
               &      - a(1,3)*a(3,1)*a(4,4) - a(1,1)*a(3,4)*a(4,3))*denom
          inva(3,2) =(- a(1,1)*a(3,2)*a(4,4) - a(1,2)*a(3,4)*a(4,1)&
               &      - a(1,4)*a(3,1)*a(4,2) + a(1,4)*a(3,2)*a(4,1)&
               &      + a(1,2)*a(3,1)*a(4,4) + a(1,1)*a(3,4)*a(4,2))*denom
          inva(4,2) =(  a(1,1)*a(3,2)*a(4,3) + a(1,2)*a(3,3)*a(4,1)&
               &      + a(1,3)*a(3,1)*a(4,2) - a(1,3)*a(3,2)*a(4,1)&
               &      - a(1,2)*a(3,1)*a(4,3) - a(1,1)*a(3,3)*a(4,2))*denom
          inva(1,3) =(  a(1,2)*a(2,3)*a(4,4) + a(1,3)*a(2,4)*a(4,2)&
               &      + a(1,4)*a(2,2)*a(4,3) - a(1,3)*a(2,2)*a(4,4)&
               &      - a(1,2)*a(2,4)*a(4,3) - a(1,4)*a(2,3)*a(4,2))*denom
          inva(2,3) =(- a(1,1)*a(2,3)*a(4,4) - a(1,3)*a(2,4)*a(4,1)&
               &      - a(1,4)*a(2,1)*a(4,3) + a(1,4)*a(2,3)*a(4,1)&
               &      + a(1,3)*a(2,1)*a(4,4) + a(1,1)*a(2,4)*a(4,3))*denom
          inva(3,3) =(  a(1,1)*a(2,2)*a(4,4) + a(1,2)*a(2,4)*a(4,1)&
               &      + a(1,4)*a(2,1)*a(4,2) - a(1,4)*a(2,2)*a(4,1)&
               &      - a(1,2)*a(2,1)*a(4,4) - a(1,1)*a(2,4)*a(4,2))*denom
          inva(4,3) =(- a(1,1)*a(2,2)*a(4,3) - a(1,2)*a(2,3)*a(4,1)&
               &      - a(1,3)*a(2,1)*a(4,2) + a(1,3)*a(2,2)*a(4,1)&
               &      + a(1,2)*a(2,1)*a(4,3) + a(1,1)*a(2,3)*a(4,2))*denom
          inva(1,4) =(- a(1,2)*a(2,3)*a(3,4) - a(1,3)*a(2,4)*a(3,2)&
               &      - a(1,4)*a(2,2)*a(3,3) + a(1,4)*a(2,3)*a(3,2)&
               &      + a(1,3)*a(2,2)*a(3,4) + a(1,2)*a(2,4)*a(3,3))*denom
          inva(2,4) =(  a(1,1)*a(2,3)*a(3,4) + a(1,3)*a(2,4)*a(3,1)&
               &      + a(1,4)*a(2,1)*a(3,3) - a(1,4)*a(2,3)*a(3,1)&
               &      - a(1,3)*a(2,1)*a(3,4) - a(1,1)*a(2,4)*a(3,3))*denom
          inva(3,4) =(- a(1,1)*a(2,2)*a(3,4) - a(1,2)*a(2,4)*a(3,1)&
               &      - a(1,4)*a(2,1)*a(3,2) + a(1,4)*a(2,2)*a(3,1)&
               &      + a(1,2)*a(2,1)*a(3,4) + a(1,1)*a(2,4)*a(3,2))*denom
          inva(4,4) =(  a(1,1)*a(2,2)*a(3,3) + a(1,2)*a(2,3)*a(3,1)&
               &      + a(1,3)*a(2,1)*a(3,2) - a(1,3)*a(2,2)*a(3,1)&
               &      - a(1,2)*a(2,1)*a(3,3) - a(1,1)*a(2,3)*a(3,2))*denom

       case default

          inva(1:nsize,1:nsize) = a(1:nsize,1:nsize)
          det = 1.0_rp

          do n = 1,nsize
             t1 = inva(n,n)
             if( t1 == 0.0_rp ) then
                det = 0.0_rp
                goto 10
             end if
             do jj = 1,nsize
                inva(n,jj) = -inva(n,jj)/t1
             end do
             do ii = 1,nsize
                if ( n /= ii ) then
                   do jj = 1,nsize
                      if( n /= jj ) inva(ii,jj) = inva(ii,jj) + inva(ii,n) * inva(n,jj)
                   end do
                end if
                inva(ii,n) = inva(ii,n)/t1
             end do
             inva(n,n) = 1.0_rp/t1
          end do

       end select

    else

       if( nsize <= 4 ) then
          allocate( b(nsize,nsize) )
          b(1:nsize,1:nsize) = a(1:nsize,1:nsize)
       end if

       select case ( nsize )

       case ( 1_ip )

          det = b(1,1)
          if( abs(det) == 0.0_rp ) goto 20
          a(1,1) = 1.0_rp/b(1,1)

       case( 2_ip )

          det = b(1,1)*b(2,2)-b(2,1)*b(1,2)
          if( abs(det) == 0.0_rp ) goto 20
          denom  = 1.0_rp/det
          a(1,1) = b(2,2)*denom
          a(2,2) = b(1,1)*denom
          a(2,1) =-b(2,1)*denom
          a(1,2) =-b(1,2)*denom

       case ( 3_ip )

          t1  = b(2,2)*b(3,3) - b(3,2)*b(2,3)
          t2  =-b(2,1)*b(3,3) + b(3,1)*b(2,3)
          t3  = b(2,1)*b(3,2) - b(3,1)*b(2,2)
          det = b(1,1)*t1 + b(1,2)*t2 + b(1,3)*t3
          if( abs(det) == 0.0_rp ) goto 20
          denom  = 1.0_rp/det
          a(1,1) = t1*denom
          a(2,1) = t2*denom
          a(3,1) = t3*denom
          a(2,2) = ( b(1,1)*b(3,3) - b(3,1)*b(1,3))*denom
          a(3,2) = (-b(1,1)*b(3,2) + b(1,2)*b(3,1))*denom
          a(3,3) = ( b(1,1)*b(2,2) - b(2,1)*b(1,2))*denom
          a(1,2) = (-b(1,2)*b(3,3) + b(3,2)*b(1,3))*denom
          a(1,3) = ( b(1,2)*b(2,3) - b(2,2)*b(1,3))*denom
          a(2,3) = (-b(1,1)*b(2,3) + b(2,1)*b(1,3))*denom

       case ( 4_ip )

          t1=    b(2,2)*b(3,3)*b(4,4) + b(2,3)*b(3,4)*b(4,2)&
               + b(2,4)*b(3,2)*b(4,3) - b(2,3)*b(3,2)*b(4,4)&
               - b(2,2)*b(3,4)*b(4,3) - b(2,4)*b(3,3)*b(4,2)
          t2=  - b(2,1)*b(3,3)*b(4,4) - b(2,3)*b(3,4)*b(4,1)&
               - b(2,4)*b(3,1)*b(4,3) + b(2,4)*b(3,3)*b(4,1)&
               + b(2,3)*b(3,1)*b(4,4) + b(2,1)*b(3,4)*b(4,3)
          t3=    b(2,1)*b(3,2)*b(4,4) + b(2,2)*b(3,4)*b(4,1)&
               + b(2,4)*b(3,1)*b(4,2) - b(2,4)*b(3,2)*b(4,1)&
               - b(2,2)*b(3,1)*b(4,4) - b(2,1)*b(3,4)*b(4,2)
          t4=  - b(2,1)*b(3,2)*b(4,3) - b(2,2)*b(3,3)*b(4,1)&
               - b(2,3)*b(3,1)*b(4,2) + b(2,3)*b(3,2)*b(4,1)&
               + b(2,2)*b(3,1)*b(4,3) + b(2,1)*b(3,3)*b(4,2)
          det= b(1,1)*t1 + b(1,2)*t2 + b(1,3)*t3 + b(1,4)*t4
          if( abs(det) == 0.0_rp ) goto 20
          denom  = 1.0_rp/det
          a(1,1) = t1*denom
          a(2,1) = t2*denom
          a(3,1) = t3*denom
          a(4,1) = t4*denom
          a(1,2) =(- b(1,2)*b(3,3)*b(4,4) - b(1,3)*b(3,4)*b(4,2)&
               &   - b(1,4)*b(3,2)*b(4,3) + b(1,3)*b(3,2)*b(4,4)&
               &   + b(1,2)*b(3,4)*b(4,3) + b(1,4)*b(3,3)*b(4,2))*denom
          a(2,2) =(  b(1,1)*b(3,3)*b(4,4) + b(1,3)*b(3,4)*b(4,1)&
               &   + b(1,4)*b(3,1)*b(4,3) - b(1,4)*b(3,3)*b(4,1)&
               &   - b(1,3)*b(3,1)*b(4,4) - b(1,1)*b(3,4)*b(4,3))*denom
          a(3,2) =(- b(1,1)*b(3,2)*b(4,4) - b(1,2)*b(3,4)*b(4,1)&
               &   - b(1,4)*b(3,1)*b(4,2) + b(1,4)*b(3,2)*b(4,1)&
               &   + b(1,2)*b(3,1)*b(4,4) + b(1,1)*b(3,4)*b(4,2))*denom
          a(4,2) =(  b(1,1)*b(3,2)*b(4,3) + b(1,2)*b(3,3)*b(4,1)&
               &   + b(1,3)*b(3,1)*b(4,2) - b(1,3)*b(3,2)*b(4,1)&
               &   - b(1,2)*b(3,1)*b(4,3) - b(1,1)*b(3,3)*b(4,2))*denom
          a(1,3) =(  b(1,2)*b(2,3)*b(4,4) + b(1,3)*b(2,4)*b(4,2)&
               &   + b(1,4)*b(2,2)*b(4,3) - b(1,3)*b(2,2)*b(4,4)&
               &   - b(1,2)*b(2,4)*b(4,3) - b(1,4)*b(2,3)*b(4,2))*denom
          a(2,3) =(- b(1,1)*b(2,3)*b(4,4) - b(1,3)*b(2,4)*b(4,1)&
               &   - b(1,4)*b(2,1)*b(4,3) + b(1,4)*b(2,3)*b(4,1)&
               &   + b(1,3)*b(2,1)*b(4,4) + b(1,1)*b(2,4)*b(4,3))*denom
          a(3,3) =(  b(1,1)*b(2,2)*b(4,4) + b(1,2)*b(2,4)*b(4,1)&
               &   + b(1,4)*b(2,1)*b(4,2) - b(1,4)*b(2,2)*b(4,1)&
               &   - b(1,2)*b(2,1)*b(4,4) - b(1,1)*b(2,4)*b(4,2))*denom
          a(4,3) =(- b(1,1)*b(2,2)*b(4,3) - b(1,2)*b(2,3)*b(4,1)&
               &   - b(1,3)*b(2,1)*b(4,2) + b(1,3)*b(2,2)*b(4,1)&
               &   + b(1,2)*b(2,1)*b(4,3) + b(1,1)*b(2,3)*b(4,2))*denom
          a(1,4) =(- b(1,2)*b(2,3)*b(3,4) - b(1,3)*b(2,4)*b(3,2)&
               &   - b(1,4)*b(2,2)*b(3,3) + b(1,4)*b(2,3)*b(3,2)&
               &   + b(1,3)*b(2,2)*b(3,4) + b(1,2)*b(2,4)*b(3,3))*denom
          a(2,4) =(  b(1,1)*b(2,3)*b(3,4) + b(1,3)*b(2,4)*b(3,1)&
               &   + b(1,4)*b(2,1)*b(3,3) - b(1,4)*b(2,3)*b(3,1)&
               &   - b(1,3)*b(2,1)*b(3,4) - b(1,1)*b(2,4)*b(3,3))*denom
          a(3,4) =(- b(1,1)*b(2,2)*b(3,4) - b(1,2)*b(2,4)*b(3,1)&
               &   - b(1,4)*b(2,1)*b(3,2) + b(1,4)*b(2,2)*b(3,1)&
               &   + b(1,2)*b(2,1)*b(3,4) + b(1,1)*b(2,4)*b(3,2))*denom
          a(4,4) =(  b(1,1)*b(2,2)*b(3,3) + b(1,2)*b(2,3)*b(3,1)&
               &   + b(1,3)*b(2,1)*b(3,2) - b(1,3)*b(2,2)*b(3,1)&
               &   - b(1,2)*b(2,1)*b(3,3) - b(1,1)*b(2,3)*b(3,2))*denom

       case default

          det = 1.0_rp
          do n = 1,nsize
             t1 = a(n,n)
             if( t1 == 0.0_rp ) then
                det = 0.0_rp
                goto 20
             end if
             do jj = 1,nsize
                a(n,jj) = -a(n,jj)/t1
             end do
             do ii = 1,nsize
                if ( n /= ii ) then
                   do jj = 1,nsize
                      if( n /= jj ) a(ii,jj) = a(ii,jj) + a(ii,n) * a(n,jj)
                   end do
                end if
                a(ii,n) = a(ii,n)/t1
             end do
             a(n,n) = 1.0_rp/t1
          end do

       end select

20     if( nsize <= 4 ) deallocate(b)

    end if

10  if( present(deter) ) deter = det

  end subroutine maths_invert_matrix

  !-----------------------------------------------------------------------
  !> 
  !> @author  Guillaume
  !> @date    2020-05-20
  !> @brief   Conjugate gradient
  !> @details Very simple sequential conjugate gradient
  !> 
  !-----------------------------------------------------------------------

  pure subroutine maths_conjugate_gradient(n,A,b,x,TOLERANCE,ITERATIONS,RELATIVE_TOLERANCE,ierr)

    integer(ip), intent(in)             :: n
    real(rp),    intent(in)             :: A(n,n)
    real(rp),    intent(in)             :: b(n)
    real(rp),    intent(inout)          :: x(n)
    real(rp),    intent(in),   optional :: TOLERANCE
    integer(ip), intent(in),   optional :: ITERATIONS
    logical(lg), intent(in),   optional :: RELATIVE_TOLERANCE
    integer(ip), intent(out),  optional :: ierr
    integer(ip)                         :: maxiter,iiter,ii,ierr_loc
    real(rp),    allocatable            :: rr(:),pp(:),zz(:),qq(:),invdiag(:)
    real(rp)                            :: alpha,beta,rho,newrho,toler,resid,bnorm
    
    allocate(rr(n))
    allocate(pp(n))
    allocate(zz(n))
    allocate(qq(n))
    allocate(invdiag(n))
    
    ierr_loc = 0
    !
    ! diag(A)
    !
    do ii = 1,n
       if( abs(A(ii,ii)) /= 0.0_rp ) then
          invdiag(ii) = 1.0_rp / A(ii,ii)
       else
          invdiag(ii) = 1.0_rp
       end if
    end do
    !
    ! Initial residuals
    !
    rr       = b - matmul(A,x)
    zz       = invdiag * rr
    pp       = zz
    qq       = invdiag * b
    newrho   = dot_product(rr,zz)       
    bnorm    = sqrt(dot_product(b(1:n),qq(1:n)))    
    if( bnorm < 1.0e-16_rp ) bnorm = 1.0_rp
    resid    = sqrt(newrho)/bnorm
    !
    ! Options
    !
    if( present(TOLERANCE) ) then
       toler = TOLERANCE
    else
       toler = 1.0e-6_rp
    end if
    if( present(ITERATIONS)) then
       maxiter = ITERATIONS
    else
       maxiter = 100
    end if
    
    if( present(RELATIVE_TOLERANCE) ) then
       if( RELATIVE_TOLERANCE ) toler = TOLERANCE * resid
    end if
    !
    ! Main loop
    !
    iiter    = 0
    do while( iiter < maxiter .and. resid > toler )
       iiter = iiter + 1
       qq     = matmul(A,pp)
       alpha  = dot_product(pp,qq)
       if( alpha <= 0.0_rp ) then
          ierr_loc = 2
          goto 10
       end if
       alpha  = newrho / alpha
       x      = x + alpha * pp
       rr     = rr - alpha * qq
       zz     = invdiag * rr
       rho    = newrho
       newrho = dot_product(rr,zz)
       if( newrho <= 0.0_rp ) then
          ierr_loc = 3
          goto 10
       end if
       beta   = newrho / rho
       resid  = sqrt(newrho)/bnorm
       pp     = zz + beta * pp
    end do
    
    if( iiter >= maxiter .or. resid > toler ) ierr_loc = 1

10  continue

    if( present(ierr) ) ierr = ierr_loc
    
    deallocate(rr)
    deallocate(pp)
    deallocate(zz)
    deallocate(qq)
    deallocate(invdiag)

  end subroutine maths_conjugate_gradient
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  Guillaume
  !> @date    2020-05-20
  !> @brief   Conjugate gradient
  !> @details Very simple sequential conjugate gradient
  !> 
  !-----------------------------------------------------------------------

  subroutine maths_direct(n,A,b,x,ierr)

    integer(ip), intent(in)             :: n
    real(rp),    intent(in)             :: A(n,n)
    real(rp),    intent(in)             :: b(n)
    real(rp),    intent(out)            :: x(n)
    integer(ip), intent(out),  optional :: ierr
    real(rp)                            :: C(n,n)
    real(rp)                            :: deter
    
    C = A
    call maths_invert_matrix(n,C,deter)
    x = matmul(C,b)

    if( present(ierr) ) then
       ierr = 0
       if( deter <= 0.0_rp ) ierr = 1
    end if
    
  end subroutine maths_direct

  !-----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @brief   Schur complement
  !> @details Compute the Schur complement of dense matrices
  !>          A11 = A11 - A12.A22^-1.A21
  !>          A11(n1,n2), A12(n1,n3), A22(n3,n3), A21(n3,n2)
  !>
  !-----------------------------------------------------------------------

  subroutine maths_schur_complement(n1,n2,n3,A11,A12,A22,A21)

    integer(ip), intent(in)    :: n1
    integer(ip), intent(in)    :: n2
    integer(ip), intent(in)    :: n3
    real(rp),    intent(inout) :: A11(n1,n2)
    real(rp),    intent(in)    :: A12(n1,n3)
    real(rp),    intent(in)    :: A22(n3,n3)
    real(rp),    intent(in)    :: A21(n3,n2)
    real(rp)                   :: B21(n3,n2)
    integer(ip)                :: i,j,k
    real(rp)                   :: deter
    real(rp)                   :: A22inv(n3,n3)
    !
    ! Using the fact that: A_{n x m} X B_{m x p} => (AB)_ij = sum_k=1^m Aik*Bkj
    !
    if( n3 == 1 ) then
       !
       ! B21_{n3 x n2} = A22^-1_{n3 x n3} X A21_{n3 X n2}
       !
       B21 = 0.0_rp
       do j = 1,n2
          B21(1,j) = B21(1,j) + A21(1,j) / A22(1,1)
       end do
       !
       ! A11_{n1 x n2} = A11_{n1 x n2} - A12_{n1 x n3} X B21_{n3 x n2}
       !
       do i = 1,n1
          do j = 1,n2
             A11(i,j) = A11(i,j) - A12(i,1) * B21(1,j)
          end do
       end do

    else
       !
       ! Invert A22
       !
       A22inv = A22
       call maths_invert_matrix(n3,A22inv,deter)
       !
       ! B21_{n3 x n2} = A22^-1_{n3 x n3} X A21_{n3 X n2}
       !
       B21 = 0.0_rp
       do i = 1,n3
          do j = 1,n2
             do k = 1,n3
                B21(i,j) = B21(i,j) + A21(k,j) * A22inv(i,k)
             end do
          end do
       end do
       !
       ! A11_{n1 x n2} = A11_{n1 x n2} - A12_{n1 x n3} X B21_{n3 x n2}
       !
       do i = 1,n1
          do j = 1,n2
             do k = 1,n3
                A11(i,j) = A11(i,j) - A12(i,k) * B21(k,j)
             end do
          end do
       end do
    end if

  end subroutine maths_schur_complement

  !-----------------------------------------------------------------------
  !> 
  !> @author  lehmkuhl and houzeaux
  !> @date    2020-03-20
  !> @brief   Eigenvalues of a symetric real 3x3 matrix
  !> @details based on:
  !>          M.J. Kronenburg, A Method for Fast Diagonalization of a 2x2 or
  !>          3x3 Real Symmetric Matrix, arXiv:1306.6291, 2015
  !>          https://arxiv.org/abs/1306.6291
  !>
  !>          See also: Smith, Oliver K. (April 1961), "Eigenvalues of a
  !>          symmetric 3x3 matrix.", Communications of the ACM, 4 (4): 168
  !>          For the genral case, see:
  !>          https://www.geometrictools.com/Documentation/RobustEigenSymmetric3x3.pdf
  !>   
  !-----------------------------------------------------------------------

  subroutine maths_eigen_3x3_symmetric_matrix(Ain,eigVal,eigVec,SOLUTION_METHOD,QUAD_REAL,maxit,toler,iterations)

    real(rp),               intent(in)  :: Ain(3,3)
    real(rp),               intent(out) :: eigVal(3)
    real(rp),     optional, intent(out) :: eigVec(3,3)
    character(*), optional, intent(in)  :: SOLUTION_METHOD
    logical(lg),  optional, intent(in)  :: QUAD_REAL
    integer(ip),  optional, intent(in)  :: maxit
    real(rp),     optional, intent(in)  :: toler
    integer(ip),  optional, intent(out) :: iterations
    integer(ip)                         :: ii,jj,ierr
    character(10)                       :: method
    logical(lg)                         :: if_quad_real

    if( present(QUAD_REAL) ) then
       if_quad_real = QUAD_REAL
    else
       if_quad_real = .false.
    end if

    if( present(SOLUTION_METHOD) ) then
       method = trim(SOLUTION_METHOD)
    else
       method = 'direct'
    end if

    if( .not. if_quad_real ) then

       block

         real(rp) :: A(3,3),D(3,3)
         real(rp) :: delta
         real(rp) :: bb,ee,alpha,p1,pp,qq
         real(rp) :: traceA,extraA,strangA

         alpha = maxval(abs(Ain))

         if( alpha == 0.0_rp ) then

            eigVal(:)   = 0.0_rp
            if( present(eigVec) ) then
               eigVec(:,1) = (/ 1.0_rp , 0.0_rp , 0.0_rp /)
               eigVec(:,2) = (/ 0.0_rp , 1.0_rp , 0.0_rp /)
               eigVec(:,3) = (/ 0.0_rp , 0.0_rp , 1.0_rp /)
            end if

         else

            alpha = 1.0_rp / alpha  
            A     = Ain * alpha

            p1 = A(1,2)**2 + A(1,3)**2 + A(2,3)**2

            if( p1 < epsilon(1.0_rp)*alpha ) then 
               !
               ! A is diagonal
               !
               eigVal(1)   = A(1,1) 
               eigVal(2)   = A(2,2)
               eigVal(3)   = A(3,3)
               if( present(eigVec) ) then
                  eigVec(:,1) = (/ 1.0_rp , 0.0_rp , 0.0_rp /)
                  eigVec(:,2) = (/ 0.0_rp , 1.0_rp , 0.0_rp /)
                  eigVec(:,3) = (/ 0.0_rp , 0.0_rp , 1.0_rp /)
               end if

            else

               traceA    = A(1,1) + A(2,2) + A(3,3)
               extraA    = A(1,2)**2 + A(1,3)**2 + A(2,3)**2
               strangA   = A(1,1)*A(2,3)**2 + A(2,2)*A(1,3)**2 + A(3,3)*A(1,2)**2
               bb        = traceA
               pp        = 0.5_rp*( (A(1,1)-A(2,2))**2 + (A(1,1)-A(3,3))**2 + (A(2,2)-A(3,3))**2 ) + 3.0_rp * extraA
               qq        =    18.0_rp * ( A(1,1)*A(2,2)*A(3,3) + 3.0_rp*A(1,2)*A(1,3)*A(2,3) ) &
                    &      +   2.0_rp * ( A(1,1)**3 + A(2,2)**3 + A(3,3)**3 ) & 
                    &      +   9.0_rp * traceA * extraA &
                    &      -   3.0_rp * ( A(1,1) + A(2,2) ) * ( A(1,1) + A(3,3) ) * ( A(2,2) + A(3,3) ) &
                    &      -  27.0_rp * strangA
               ee        = qq/(2.0_rp*sqrt(pp*pp*pp))
               ee        = min(1.0_rp,max(-1.0_rp,ee))
               Delta     = acos(ee)
               eigVal(1) = 1.0_rp/3.0_rp * ( bb + 2.0_rp * sqrt(pp) * cos(delta/3.0_rp))
               eigVal(2) = 1.0_rp/3.0_rp * ( bb + 2.0_rp * sqrt(pp) * cos((delta+2.0_rp*pi)/3.0_rp))
               eigVal(3) = 1.0_rp/3.0_rp * ( bb + 2.0_rp * sqrt(pp) * cos((delta-2.0_rp*pi)/3.0_rp))

               if( present(eigVec) ) then

                  select case ( trim(method) )

                  case ( 'cg' ) 

                     do ii = 1,3
                        D = A
                        do jj = 1,3
                           D(jj,jj) = D(jj,jj) - eigVal(ii)
                        end do
                        if( D(1,1) < 0.0_rp ) then
                           do jj = 1,3
                              D(jj,jj) = -D(jj,jj)
                           end do
                        end if
                        eigVec(:,ii)  = 0.0_rp
                        eigVec(ii,ii) = 1.0_rp
                        call maths_conjugate_gradient(3_ip,D,(/0.0_rp,0.0_rp,0.0_rp/),eigVec(:,ii),&
                             TOLERANCE=1.0e-8_rp,ITERATIONS=100_ip,IERR=ierr)
                     end do

                  case ( 'direct' ) 

                     if( abs(qq**2-4.0_rp*pp**3) < epsil ) then
                        !
                        ! Only two eigenvalues
                        !
                        call maths_eigen_3x3_symmetric_matrix_vect(1_ip,A,eigVal,eigVec)

                     else
                        !
                        ! General case
                        !
                        call maths_eigen_3x3_symmetric_matrix_vect(2_ip,A,eigVal,eigVec)

                     end if

                  end select
                  !
                  ! Renormalize
                  !
                  do ii = 1,3
                     eigVec(:,ii) = eigVec(:,ii) / sqrt(dot_product(eigVec(:,ii),eigVec(:,ii)))
                  end do
               end if
            end if
            eigVal(:) = eigVal(:) / alpha

         end if

       end block

    else
       !
       ! QUAD PRECISION
       !
       block

         real(qp), parameter :: epsil_qp = epsilon(1.0_qp)
         real(qp), parameter :: pi_qp    = 3.141592653589793238462643383279502884197_qp
         real(qp)            :: eigVec_qp(3,3)
         real(qp)            :: eigVal_qp(3)
         real(qp)            :: A_qp(3,3),D_qp(3,3),R_qp(2,2),Ain_qp(3,3)
         real(qp)            :: delta_qp
         real(qp)            :: bb_qp,ee_qp,vv_qp,ww_qp,alpha_qp,p1_qp,pp_qp,qq_qp
         real(qp)            :: R1_qp(3,3),R2_qp(3,3),R3_qp(3,3),phi2_qp,phi3_qp
         real(qp)            :: f1_qp(2),f2_qp(2),g1_qp(2),g2_qp(2),Rg_qp(2)
         real(qp)            :: g1x_qp,g1y_qp,g2x_qp,g2y_qp,psi1_qp,psi2_qp,phi11_qp,phi12_qp
         real(qp)            :: f1n_qp,f2n_qp,l12_qp,l23_qp,l13_qp,l31_qp,lambda_qp,lambda3_qp
         real(qp)            :: diff_min_qp,phi1_min_qp
         real(qp)            :: phi2_min_qp,phi3_min_qp,traceA_qp,extraA_qp,strangA_qp
         real(qp)            :: phi_sign_qp(4,2)

         Ain_qp   = real(Ain,qp)
         alpha_qp = maxval(abs(Ain_qp))

         if( alpha_qp == 0.0_qp ) then

            eigVal_qp(:)   = 0.0_qp
            if( present(eigVec) ) then
               eigVec_qp(:,1) = (/ 1.0_qp , 0.0_qp , 0.0_qp /)
               eigVec_qp(:,2) = (/ 0.0_qp , 1.0_qp , 0.0_qp /)
               eigVec_qp(:,3) = (/ 0.0_qp , 0.0_qp , 1.0_qp /)
            end if

         else

            alpha_qp = 1.0_qp / alpha_qp
            A_qp     = Ain_qp * alpha_qp

            p1_qp = A_qp(1,2)**2 + A_qp(1,3)**2 + A_qp(2,3)**2

            if( p1_qp == 0.0_qp ) then 
               !
               ! A is diagonal
               !
               eigVal_qp(1)   = A_qp(1,1) 
               eigVal_qp(2)   = A_qp(2,2)
               eigVal_qp(3)   = A_qp(3,3)
               if( present(eigVec) ) then
                  eigVec_qp(:,1) = (/ 1.0_qp , 0.0_qp , 0.0_qp /)
                  eigVec_qp(:,2) = (/ 0.0_qp , 1.0_qp , 0.0_qp /)
                  eigVec_qp(:,3) = (/ 0.0_qp , 0.0_qp , 1.0_qp /)
               end if

            else

               traceA_qp    = A_qp(1,1) + A_qp(2,2) + A_qp(3,3)
               extraA_qp    = A_qp(1,2)**2 + A_qp(1,3)**2 + A_qp(2,3)**2
               strangA_qp   = A_qp(1,1)*A_qp(2,3)**2 + A_qp(2,2)*A_qp(1,3)**2 + A_qp(3,3)*A_qp(1,2)**2
               bb_qp        = traceA_qp
               pp_qp        =     0.5_qp * ( (A_qp(1,1)-A_qp(2,2))**2 + (A_qp(1,1)-A_qp(3,3))**2 + (A_qp(2,2)-A_qp(3,3))**2 ) + 3.0_qp * extraA_qp
               qq_qp        =    18.0_qp * ( A_qp(1,1)*A_qp(2,2)*A_qp(3,3) + 3.0_qp*A_qp(1,2)*A_qp(1,3)*A_qp(2,3) ) &
                    &         +   2.0_qp * ( A_qp(1,1)**3 + A_qp(2,2)**3 + A_qp(3,3)**3 ) & 
                    &         +   9.0_qp * traceA_qp * extraA_qp &
                    &         -   3.0_qp * ( A_qp(1,1) + A_qp(2,2) ) * ( A_qp(1,1) + A_qp(3,3) ) * ( A_qp(2,2) + A_qp(3,3) ) &
                    &         -  27.0_qp * strangA_qp
               ee_qp        = qq_qp/(2.0_qp*abs(pp_qp)*sqrt(pp_qp))
               ee_qp        = min(1.0_qp,max(-1.0_qp,ee_qp))
               Delta_qp     = acos(ee_qp)
               eigVal_qp(1) = 1.0_qp/3.0_qp * ( bb_qp + 2.0_qp * sqrt(pp_qp) * cos( delta_qp              /3.0_qp))
               eigVal_qp(2) = 1.0_qp/3.0_qp * ( bb_qp + 2.0_qp * sqrt(pp_qp) * cos((delta_qp+2.0_qp*pi_qp)/3.0_qp))
               eigVal_qp(3) = 1.0_qp/3.0_qp * ( bb_qp + 2.0_qp * sqrt(pp_qp) * cos((delta_qp-2.0_qp*pi_qp)/3.0_qp))

               if( present(eigVec) ) then

                  select case ( trim(method) )

                  case ( 'cg' ) 

                     !write(*,*) 'MOD_MATHS_SOLVER: CG ALGORITHM NOT AVAILABLE FOR QP'

                  case ( 'direct' ) 

                     diff_min_qp = huge(1.0_qp)
                     f1_qp   = (/ A_qp(1,2),-A_qp(1,3) /)
                     f2_qp   = (/ A_qp(2,2)-A_qp(3,3),-2.0_qp*A_qp(2,3) /)
                     f1n_qp  = sqrt(dot_product(f1_qp,f1_qp))  
                     f2n_qp  = sqrt(dot_product(f2_qp,f2_qp))
                     psi1_qp = maths_angle_vector_2D_qp(f1_qp)
                     psi2_qp = maths_angle_vector_2D_qp(f2_qp)

                     if( abs(qq_qp**2-4.0_qp*pp_qp**3) < epsil_qp ) then
                        !
                        ! Only two eigenvalues
                        !
                        if(      abs(eigVal_qp(1)-eigVal_qp(2)) < epsil_qp ) then
                           lambda3_qp = eigVal_qp(3)
                           lambda_qp  = eigVal_qp(1)
                        else if( abs(eigVal_qp(1)-eigVal_qp(3)) < epsil_qp ) then
                           lambda3_qp = eigVal_qp(2)
                           lambda_qp  = eigVal_qp(1)
                        else 
                           lambda3_qp = eigVal_qp(1)
                           lambda_qp  = eigVal_qp(2)
                        end if
                        eigVal_qp(1) = lambda_qp
                        eigVal_qp(2) = lambda_qp
                        eigVal_qp(3) = lambda3_qp
                        phi_sign_qp  = reshape((/ 1.0_qp , -1.0_qp ,  1.0_qp , -1.0_qp,    &
                             &                 0.0_qp ,  0.0_qp , -1.0_qp , -1.0_qp /), (/4,2/))

                        l13_qp  = lambda_qp-lambda3_qp
                        pp_qp   = A_qp(1,1)-lambda3_qp 
                        if( l13_qp == 0.0_qp ) then
                           vv_qp= 1.0_qp
                        else
                           vv_qp = pp_qp / l13_qp
                        end if
                        ww_qp   = 1.0_qp
                        vv_qp   = abs(min(1.0_qp,max(0.0_qp,vv_qp)))

                        do ii = 1,2

                           if( ii == 1 ) then
                              phi2_qp =  acos(sqrt(vv_qp))
                           else
                              phi2_qp = -acos(sqrt(vv_qp))
                           end if
                           phi3_qp = 0.0_qp

                           g1_qp   = (/ 0.0_qp , 0.5_qp * l13_qp * sin(2.0_qp*phi2_qp) /)
                           g2_qp   = (/ l13_qp * vv_qp , 0.0_qp /)


                           R_qp     = maths_rotation_matrix_2D_qp(-psi1_qp)
                           Rg_qp    = matmul(R_qp,g1_qp)
                           phi11_qp = maths_angle_vector_2D_qp(Rg_qp)

                           R_qp     = maths_rotation_matrix_2D_qp(-psi2_qp)
                           Rg_qp    = matmul(R_qp,g2_qp)
                           phi12_qp = 0.5_qp * maths_angle_vector_2D_qp(Rg_qp)

                           if( abs(phi11_qp-phi12_qp) < diff_min_qp ) then
                              diff_min_qp = abs(phi11_qp-phi12_qp)
                              phi2_min_qp = phi2_qp
                              phi3_min_qp = phi3_qp
                              if( f1n_qp >= f2n_qp ) then
                                 phi1_min_qp = phi11_qp
                              else
                                 phi1_min_qp = phi12_qp
                              end if
                           end if
                        end do

                     else
                        !
                        ! General case
                        !
                        phi_sign_qp  = reshape( (/ 1.0_qp , -1.0_qp ,  1.0_qp , -1.0_qp,    &
                             &                  1.0_qp ,  1.0_qp , -1.0_qp , -1.0_qp /), (/4,2/))

                        l12_qp       = eigVal_qp(1)-eigVal_qp(2)
                        l23_qp       = eigVal_qp(2)-eigVal_qp(3)
                        l31_qp       = eigVal_qp(3)-eigVal_qp(1)

                        pp_qp =  ( A_qp(1,2)**2 + A_qp(1,3)**2 + (A_qp(1,1)-eigVal_qp(3)) * (A_qp(1,1)+l31_qp-eigVal_qp(2)) )
                        if( l23_qp * l31_qp == 0.0_qp ) then
                           vv_qp = 1.0_qp
                        else 
                           vv_qp = pp_qp / ( l23_qp * l31_qp )
                        end if
                        vv_qp = min(1.0_qp,max(0.0_qp,vv_qp))

                        qq_qp = A_qp(1,1) - eigVal_qp(3) - l23_qp * vv_qp
                        if( l12_qp * vv_qp == 0.0_qp ) then
                           ww_qp = 1.0_qp
                        else
                           ww_qp = qq_qp / ( l12_qp * vv_qp )
                        end if
                        ww_qp   = min(1.0_qp,max(0.0_qp,ww_qp))

                        do ii = 1,4

                           phi2_qp  = phi_sign_qp(ii,1) * acos(sqrt(vv_qp))
                           phi3_qp  = phi_sign_qp(ii,2) * acos(sqrt(ww_qp))

                           g1x_qp   = 0.5_qp * l12_qp * cos(phi2_qp) * sin(2.0_qp*phi3_qp)
                           g1y_qp   = 0.5_qp * ( l12_qp * ww_qp + l23_qp ) * sin(2.0_qp*phi2_qp)
                           g2x_qp   = l12_qp * ( 1.0_qp + (vv_qp-2.0_qp)*ww_qp) + l23_qp * vv_qp
                           g2y_qp   = l12_qp * sin(phi2_qp) * sin(2.0_qp*phi3_qp)

                           g1_qp    = (/ g1x_qp , g1y_qp /)
                           g2_qp    = (/ g2x_qp , g2y_qp /)

                           R_qp     = maths_rotation_matrix_2D_qp(-psi1_qp)
                           Rg_qp    = matmul(R_qp,g1_qp)
                           phi11_qp = maths_angle_vector_2D_qp(Rg_qp)

                           R_qp     = maths_rotation_matrix_2D_qp(-psi2_qp)
                           Rg_qp    = matmul(R_qp,g2_qp)
                           phi12_qp = 0.5_qp * maths_angle_vector_2D_qp(Rg_qp)

                           if( abs(phi11_qp-phi12_qp) < diff_min_qp ) then
                              diff_min_qp = abs(phi11_qp-phi12_qp)
                              phi2_min_qp = phi2_qp
                              phi3_min_qp = phi3_qp
                              if( f1n_qp >= f2n_qp ) then
                                 phi1_min_qp = phi11_qp
                              else
                                 phi1_min_qp = phi12_qp
                              end if
                           end if

                        end do

                     end if
                     R1_qp          = maths_rotation_matrix_3D_qp(phi1_min_qp,1_ip)
                     R2_qp          = maths_rotation_matrix_3D_qp(phi2_min_qp,2_ip)
                     R3_qp          = maths_rotation_matrix_3D_qp(phi3_min_qp,3_ip)
                     D_qp           = matmul(R2_qp,R3_qp)
                     D_qp           = matmul(R1_qp,D_qp)
                     eigVec_qp(:,1) = D_qp(:,1)
                     eigVec_qp(:,2) = D_qp(:,2)
                     eigVec_qp(:,3) = D_qp(:,3)
                  end select
                  !
                  ! Renormalize
                  !
                  do ii = 1,3
                     eigVec_qp(:,ii) = eigVec_qp(:,ii) / sqrt(dot_product(eigVec_qp(:,ii),eigVec_qp(:,ii)))
                  end do
                  eigVec = real(eigVec_qp,rp)
               end if
            end if

            eigVal_qp(:) = eigVal_qp(:) / alpha_qp

         end if

         eigVal = real(eigVal_qp,rp)

       end block

    end if
    !
    ! Use iterative method if not converged
    !
    if( present(iterations) ) iterations = 0
    if( present(toler) .and. present(maxit) ) then
       block
         real(rp) :: err,maxerr,dd1,dd2,den,vec1(3),vec2(3)
         maxerr = 0.0_rp
         do ii = 1,3
            vec1   = eigVal(ii)*eigVec(:,ii)
            vec2   = matmul(Ain,eigVec(:,ii))
            err    = dot_product(vec1-vec2,vec1-vec2)
            dd1    = dot_product(vec1,vec1)
            dd2    = dot_product(vec2,vec2)
            den    = sqrt(max(dd1,dd2))
            err    = sqrt(err)
            maxerr = max(maxerr,err/max(den,epsilon(1.0_rp)))
         end do
         if( maxerr > toler ) then
            call maths_eigen_jacobi(3_ip,Ain,toler,maxit,eigval,eigvec,ITERATIONS=iterations)
         end if
       end block
    end if
    !
    ! Sort decreasing order
    !
    if( present(eigVec) ) then
       if( eigVal(1) < eigVal(2) )then
          call maths_swap_rp(eigVal(1),eigVal(2))
          call maths_swap_rp(eigVec(:,1),eigVec(:,2))
       endif
       if( eigval(1) < eigVal(3) )then
          call maths_swap_rp(eigVal(1),eigVal(3))
          call maths_swap_rp(eigVec(:,1),eigVec(:,3)) 
       endif
       if( eigval(2) < eigVal(3) )then
          call maths_swap_rp(eigVal(2),eigVal(3))
          call maths_swap_rp(eigVec(:,2),eigVec(:,3)) 
       endif
    else
       if( eigVal(1) < eigVal(2) )then
          call maths_swap_rp(eigVal(1),eigVal(2))
       endif
       if( eigval(1) < eigVal(3) )then
          call maths_swap_rp(eigVal(1),eigVal(3))
       endif
       if( eigval(2) < eigVal(3) )then
          call maths_swap_rp(eigVal(2),eigVal(3))
       endif
    endif

  end subroutine maths_eigen_3x3_symmetric_matrix
   
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-03-20
  !> @brief   Eigenvalues of a symmetric real 2x2 matrix
  !> @details Eigenvalues of a symmetric real 2x2 matrix
  !>   
  !-----------------------------------------------------------------------
  
  pure subroutine maths_eigen_2x2_symmetric_matrix(A,eigVal,eigVec)
     
     real(rp),           intent(in)  :: A(2,2)
     real(rp),           intent(out) :: eigVal(2)
     real(rp), optional, intent(out) :: eigVec(2,2)
     real(rp)                        :: delta,traceA

     if( abs(A(1,2)) == 0.0_rp ) then
        eigVal(1)   = A(1,1)
        eigVal(2)   = A(2,2)
        if( present(eigVec) ) then
           eigVec(:,1) = (/ 1.0_rp , 0.0_rp /)
           eigVec(:,2) = (/ 0.0_rp , 1.0_rp /)
        end if
     else
        traceA    = A(1,1)+A(2,2)
        delta     = (A(1,1)-A(2,2))**2 + 4.0_rp*A(1,2)**2
        eigVal(1) = 0.5_rp*(traceA+sqrt(delta))
        eigVal(2) = 0.5_rp*(traceA-sqrt(delta))
        if( present(eigVec) ) then
           eigVec(1,1) = 1.0_rp
           eigVec(2,1) = (eigVal(1)-A(1,1))/A(1,2)
           eigVec(:,1) = eigVec(:,1) / sqrt(dot_product(eigVec(:,1),eigVec(:,1)))           
           eigVec(2,2) = (eigVal(2)-A(1,1))/A(1,2)
           eigVec(1,2) = 1.0_rp
           eigVec(:,2) = eigVec(:,2) / sqrt(dot_product(eigVec(:,2),eigVec(:,2)))
        end if
     end if
     ! Sort decreasing order
     if( present(eigVec) ) then
        if( eigVal(1) < eigVal(2) )then
           call maths_swap_rp(eigVal(1),eigVal(2))
           call maths_swap_rp(eigVec(:,1),eigVec(:,2))
        endif
     else
        if( eigVal(1) < eigVal(2) )then
           call maths_swap_rp(eigVal(1),eigVal(2))
        endif
     endif
     
  end subroutine maths_eigen_2x2_symmetric_matrix

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-03-20
  !> @brief   Eigenvalues of a symmetric real nxn matrix
  !> @details Eigenvalues of a symmetric real nxn matrix
  !>   
  !-----------------------------------------------------------------------
  
  subroutine maths_eigen_symmetric_matrix(n,A,eigVal,eigVec,maxit,toler,iterations)

    integer(ip),           intent(in)    :: n
    real(rp),              intent(in)    :: A(n,n)
    real(rp),              intent(out)   :: eigVal(n)
    real(rp),    optional, intent(out)   :: eigVec(n,n)
    integer(ip), optional, intent(in)    :: maxit
    real(rp),    optional, intent(in)    :: toler
    integer(ip), optional, intent(out)   :: iterations

    if( n == 2 ) then
       call maths_eigen_2x2_symmetric_matrix(A,eigVal,eigVec)
    else if( n == 3 ) then
       call maths_eigen_3x3_symmetric_matrix(A,eigVal,eigVec,maxit=maxit,toler=toler,iterations=iterations)
    else
       call runend('MOD_MATHS: EIGEN PROBLEM ONLY FOR 2X2 AND 3X3 MATRIX')
    end if

  end subroutine maths_eigen_symmetric_matrix

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-10-13
  !> @brief   Jacobi iterative solver for eigenvalue problems
  !> @details Jacobi iterative solver for eigenvalue problems:
  !>          https://en.wikipedia.org/wiki/Jacobi_eigenvalue_algorithm
  !> 
  !-----------------------------------------------------------------------

  subroutine maths_eigen_jacobi(n,A,toler,maxit,eval,evec,ITERATIONS)

    integer(ip),           intent(in)    :: n
    real(rp),              intent(in)    :: A(n,n)
    real(rp),              intent(in)    :: toler
    integer(ip),           intent(in)    :: maxit
    real(rp),              intent(inout) :: eval(n)
    real(rp),              intent(inout) :: evec(n,n)
    integer(ip), optional, intent(out)   :: ITERATIONS
    integer(ip)                          :: i,k,l,m,state,iiter
    real(rp)                             :: s,c,t,p,y,d,r
    real(rp)                             :: Eik,Eil
    integer(ip)                          :: ind(n)
    real(rp)                             :: B(n,n)
    logical(lg)                          :: changed(n)
    
    B     = A
    state = n
    !
    ! init e, E
    !
    evec = 0.0_rp
    do i = 1,n
       evec(i,i) = 1.0_rp
       eval(i)   = B(i,i)
    end do
    !
    ! init arrays ind, changed
    !
    do k = 1,n
       ind(k)     = maths_eigen_maxind(n,B,k)
       changed(k) = .true.       
    end do

    iiter = 0
    do while( state /= 0 .and. iiter < maxit )
       iiter = iiter + 1
       m = 1
       do k = 2,n-1
          if( abs(B(k,ind(k))) > abs(B(m,ind(m))) ) m=k
       end do
       k = m
       l = ind(m)
       p = B(k,l)
       !
       ! calculate c = cos(phi), s = sin(phi)
       !
       y = 0.5_rp*(eval(l)-eval(k))
       d = abs(y) + sqrt(p*p+y*y)
       r = sqrt(p*p+d*d)
       c = d/r
       s = p/r
       t = p*p/d
       if( y < 0.0_rp ) then
          s = -s
          t = -t
       end if
       B(k,l) = 0.0_rp
       call maths_eigen_update(n,A,k,-t,toler,eval,evec,state,changed)
       call maths_eigen_update(n,A,l, t,toler,eval,evec,state,changed)
       !
       ! rotate rows and columns k and l
       !
       do i = 1,k-1
          call maths_eigen_rotate(n,B,i,k,i,l,c,s)
       end do
       do i = k+1,l-1
          call maths_eigen_rotate(n,B,k,i,i,l,c,s)
       end do
       do i = l+1,n
          call maths_eigen_rotate(n,B,k,i,l,i,c,s)
       end do
       !
       ! rows k, l have changed, update rows indk, indl
       !
       do i = 1,n
          Eik       = evec(i,k)
          Eil       = evec(i,l)
          evec(i,k) = c*Eik -s*Eil
          evec(i,l) = s*Eik +c*Eil
       end do
       ind(k)    = maths_eigen_maxind(n,B,k)
       ind(l)    = maths_eigen_maxind(n,B,l)
    end do
    
    if( present(ITERATIONS) ) ITERATIONS = iiter
    
  end subroutine maths_eigen_jacobi

  pure function maths_eigen_maxind(n,S,k) result(m)
    integer(ip), intent(in)  :: n ! dimension    
    real(rp),    intent(in)  :: S(n,n)
    integer(ip), intent(in)  :: k
    integer(ip)              :: m
    integer(ip)              :: i

    m=k+1
    do i = k+2,n       
       if( abs(S(k,i)) > abs(S(k,m)) ) m=i
    end do
    
  end function maths_eigen_maxind

  pure subroutine maths_eigen_update(n,A,k,t,toler,eval,evec,state,changed)

    integer(ip), intent(in)    :: n
    integer(ip), intent(in)    :: k
    real(rp),    intent(in)    :: t
    real(rp),    intent(in)    :: A(n,n)    
    real(rp),    intent(in)    :: toler
    real(rp),    intent(inout) :: eval(n)
    real(rp),    intent(in)    :: evec(n,n)
    integer(ip), intent(inout) :: state
    logical(lg), intent(inout) :: changed(n)
    real(rp)                   :: y,err,den,dd1,dd2
    real(rp)                   :: xnorm,vec1(n),vec2(n)
    
    if( 1 == 2 ) then
       y       = eval(k)
       eval(k) = y + t
       xnorm   = max(eval(k)*eval(k),epsil)
       if( changed(k) .and. abs(t)/xnorm<toler ) then
          changed(k) = .false.
          state      = state - 1
       else if ( .not. changed(k) .and. abs(t)/xnorm>toler ) then
          changed(k) = .true.
          state      = state + 1
       end if
    else
       y       = eval(k)
       eval(k) = y + t
       vec1    = eval(k)*evec(:,k)
       vec2    = matmul(A,evec(:,k))
       err     = dot_product(vec1-vec2,vec1-vec2)
       dd1     = dot_product(vec1,vec1)
       dd2     = dot_product(vec2,vec2)
       den     = max(dd1,dd2,epsilon(1.0_rp))
       err     = sqrt(err/den)
       if( changed(k) .and. err<toler ) then
          changed(k) = .false.
          state      = state - 1
       else if ( .not. changed(k) .and. err>toler ) then
          changed(k) = .true.
          state      = state + 1
       end if
    end if
    
  end subroutine maths_eigen_update

  pure subroutine maths_eigen_rotate(n,A,k,l,i,j,c,s)
    integer(ip), intent(in)     :: n ! dimension    
    real(rp),    intent(inout)  :: A(n,n)
    integer(ip), intent(in)     :: k,l,i,j
    real(rp),    intent(in)     :: c,s
    real(rp)                    :: Skl,Sij
    
    Skl    = A(k,l)
    Sij    = A(i,j)
    A(k,l) = c*Skl -s*Sij
    A(i,j) = s*Skl +c*Sij
    
  end subroutine maths_eigen_rotate
  
  !-----------------------------------------------------------------------
  !>
  !> @author  houzeaux
  !> @date    2019-02-25
  !> @brief   Solve an overdetermined system
  !> @details Solve Ax=b with A(n,m) n > m
  !>
  !-----------------------------------------------------------------------

  subroutine maths_solve_overdetermined_system(n,m,A,b,x,error)

    integer(ip), intent(in)  :: n
    integer(ip), intent(in)  :: m
    real(rp),    intent(in)  :: A(n,m)
    real(rp),    intent(in)  :: b(n,1)
    real(rp),    intent(out) :: x(m,1)
    real(rp),    intent(out) :: error
    integer(ip)              :: ii
    real(rp),    allocatable :: At(:,:)
    real(rp),    allocatable :: At_A(:,:)
    real(rp),    allocatable :: At_A_inv(:,:)
    real(rp),    allocatable :: At_b(:,:)
    real(rp)                 :: deter,numer
    real(rp)                 :: dummr,denom

    allocate(At(m,n))
    allocate(At_A(m,m))
    allocate(At_A_inv(m,m))
    allocate(At_b(m,1))
    !
    ! x = (A^t A)^-1 (A^t b)
    !
    At   = transpose(A)
    At_A = matmul(At,A)
    call maths_invert_matrix(m,At_A,deter,At_A_inv)
    if( abs(deter) < epsil ) then
       error = -1.0_rp
    else
       At_b = matmul(At,b)
       x    = matmul(At_A_inv,At_b)
       !
       ! Error
       !
       error = 0.0_rp
       denom = 0.0_rp
       numer = 0.0_rp
       do ii = 1,n
          dummr = dot_product(x(:,1),A(ii,:))
          numer = numer + (dummr-b(ii,1))**2
          denom = denom + (b(ii,1))**2
       end do
       error = 100.0_rp*sqrt(numer/(denom+epsilon(1.0_rp)))
    end if

    deallocate(At,At_A,At_A_inv,At_b)

  end subroutine maths_solve_overdetermined_system

  pure subroutine maths_eigen_3x3_symmetric_matrix_vect(itask,A,eigVal,eigVec)

    integer(ip),            intent(in)    :: itask
    real(rp),               intent(in)    :: A(3,3)
    real(rp),               intent(inout) :: eigVal(3)
    real(rp),     optional, intent(out)   :: eigVec(3,3)
    integer(ip)                           :: ii,nsign


    real(rp) :: D(3,3)

    real(rp) :: vv,ww,pp,qq
    real(rp) :: R1(3,3),R2(3,3),R3(3,3)

    real(rp) :: f1(2),f2(2),psi1,psi2,Rpsi1(2,2),Rpsi2(2,2)
    real(rp) :: f1n,f2n,l12,l23,l13,l31,lambda,lambda3
    real(rp) :: diff_min,phi1_min
    real(rp) :: phi2_min,phi3_min
    real(rp) :: phi_sign(4,2)

    real(rp) :: phi2(4),phi3(4)
    real(rp) :: g2x
    real(rp) :: g1(4,2),g2(4,2)
    real(rp) :: phi11(4),phi12(4)
    real(rp) :: Rg(4,2)

    real(rp) :: g1s(2),g2s(2),Rgs(2)
    
    f1        = (/ A(1,2),-A(1,3) /)
    f2        = (/ A(2,2)-A(3,3),-2.0_rp*A(2,3) /)
    f1n       = sqrt(dot_product(f1,f1))  
    f2n       = sqrt(dot_product(f2,f2))
    psi1      = maths_angle_vector_2D(f1)
    psi2      = maths_angle_vector_2D(f2)
    Rpsi1     = maths_rotation_matrix_2D(-psi1)
    Rpsi2     = maths_rotation_matrix_2D(-psi2)

    select case ( itask )

    case ( 1_ip )
       !
       ! General case: abs(qq**2-4.0_rp*pp**3) < epsil
       !
       nsign = 2
       
       if(      abs(eigVal(1)-eigVal(2)) < epsil ) then
          lambda3 = eigVal(3)
          lambda  = eigVal(1)
       else if( abs(eigVal(1)-eigVal(3)) < epsil ) then
          lambda3 = eigVal(2)
          lambda  = eigVal(1)
       else 
          lambda3 = eigVal(1)
          lambda  = eigVal(2)
       end if
       eigVal(1) = lambda
       eigVal(2) = lambda
       eigVal(3) = lambda3
       phi_sign  = reshape((/ 1.0_rp , -1.0_rp ,  1.0_rp , -1.0_rp,    &
            &                 0.0_rp ,  0.0_rp , -1.0_rp , -1.0_rp /), (/4,2/))
       
       l13  = lambda-lambda3
       pp   = A(1,1)-lambda3 
       if( l13 == 0.0_rp ) then
          vv= 1.0_rp
       else
          vv = pp / l13 
       end if
       ww   = 1.0_rp
       vv   = abs(min(1.0_rp,max(0.0_rp,vv)))

       do ii = 1,2

          if( ii == 1 ) then
             phi2(ii) =  acos(sqrt(vv))
          else
             phi2(ii) = -acos(sqrt(vv))
          end if
          phi3(ii)  = 0.0_rp
          g1s(:)    = (/ 0.0_rp , 0.5_rp * l13 * sin(2.0_rp*phi2(ii)) /)
          g2s(:)    = (/ l13 * vv , 0.0_rp /)
          
          Rgs       = matmul(Rpsi1,g1s)
          phi11(ii) = maths_angle_vector_2D(Rgs)
          
          Rgs       = matmul(Rpsi2,g2s)
          phi12(ii) = 0.5_rp * maths_angle_vector_2D(Rgs)
          
       end do


    case ( 2_ip )
       !
       ! General case
       !
       nsign     = 4
       phi_sign  = reshape( (/ 1.0_rp , -1.0_rp ,  1.0_rp , -1.0_rp,    &
            &                  1.0_rp ,  1.0_rp , -1.0_rp , -1.0_rp /), (/4,2/))

       l12       = eigVal(1)-eigVal(2)
       l23       = eigVal(2)-eigVal(3)
       l31       = eigVal(3)-eigVal(1)

       pp =  ( A(1,2)**2 + A(1,3)**2 + (A(1,1)-eigVal(3)) * (A(1,1)+l31-eigVal(2)) )
       if( l23 * l31 == 0.0_rp ) then
          vv = 1.0_rp
       else 
          vv = pp / ( l23 * l31 )
       end if
       vv = min(1.0_rp,max(0.0_rp,vv))

       qq = A(1,1) - eigVal(3) - l23 * vv
       if( l12 * vv == 0.0_rp ) then
          ww = 1.0_rp
       else
          ww = qq / ( l12 * vv )
       end if
       ww       = min(1.0_rp,max(0.0_rp,ww))       
       g2x      = l12 * ( 1.0_rp + (vv-2.0_rp)*ww) + l23 * vv
       !
       ! Vectorized operations
       !
       phi2(:)  = phi_sign(:,1) * acos(sqrt(vv))
       phi3(:)  = phi_sign(:,2) * acos(sqrt(ww))

       g1(:,1)  = 0.5_rp * l12 * cos(phi2(:)) * sin(2.0_rp*phi3(:))
       g1(:,2)  = 0.5_rp * ( l12 * ww + l23 ) * sin(2.0_rp*phi2(:))
       g2(:,1)  = g2x
       g2(:,2)  = l12 * sin(phi2(:)) * sin(2.0_rp*phi3(:))

       Rg(:,:)  = maths_matmul(Rpsi1,g1)
       phi11(:) = maths_angle_vector_2D(Rg)

       Rg(:,:)  = maths_matmul(Rpsi2,g2)
       phi12(:) = 0.5_rp * maths_angle_vector_2D(Rg)

    end select
    
    diff_min = huge(1.0_rp)
    
    do ii = 1,nsign
       if( abs(phi11(ii)-phi12(ii)) < diff_min ) then
          diff_min = abs(phi11(ii)-phi12(ii))
          phi2_min = phi2(ii)
          phi3_min = phi3(ii)
          if( f1n >= f2n ) then
             phi1_min = phi11(ii)
          else
             phi1_min = phi12(ii)
          end if
       end if
    end do

    R1          = maths_rotation_matrix_3D(phi1_min,1_ip)
    R2          = maths_rotation_matrix_3D(phi2_min,2_ip)
    R3          = maths_rotation_matrix_3D(phi3_min,3_ip)
    D           = matmul(R2,R3)
    D           = matmul(R1,D)
    eigVec(:,1) = D(:,1)
    eigVec(:,2) = D(:,2)
    eigVec(:,3) = D(:,3)
    
  contains

    pure function maths_matmul(A,B) result(C)

      real(rp),  intent(in)  :: A(2,2)
      real(rp),  intent(in)  :: B(:,:)
      real(rp)               :: C(size(B,1),2)
      integer(ip)            :: i,j

      C = 0.0_rp
      
      if( size(B,1) == 4 ) then
         do j = 1,2
            do i = 1,2
               C(1:4,i) = C(1:4,i) + A(i,j) * B(1:4,j)
            end do
         end do
      else if( size(B,1) == 2 ) then
         do j = 1,2
            do i = 1,2
               C(1:2,i) = C(1:2,i) + A(i,j) * B(1:2,j)
            end do
         end do         
      end if

    end function maths_matmul

  end subroutine maths_eigen_3x3_symmetric_matrix_vect
  
end module mod_maths_solver
!> @}
