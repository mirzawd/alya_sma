!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!***************************************************************
!*
!*              Module for optimization
!* 
!***************************************************************
MODULE mod_optimization

  use def_kintyp_basic, only: ip, rp, lg
  use mod_debugTools,         only: out_debug_text
  
  IMPLICIT NONE
  SAVE
  
  integer(ip), parameter ::   NEWTON_RAPHSON = 0_ip
  integer(ip), parameter :: STEEPEST_DESCENT = 1_ip  ! DONT USE THIS -> SLOWERRR
  integer(ip), parameter :: typeFunMin_default  = NEWTON_RAPHSON
  
  real(rp), private :: tolOptimization    = 0.01_rp
  
  
  integer(ip), parameter :: DER_CENTERED   = 0_ip
  integer(ip), parameter :: DER_FORWARD    = 1_ip
  integer(ip), parameter :: typeNumerDer  = DER_FORWARD
    
  public :: optimizeFunction
  public :: setToleranceOptimization
  public :: NEWTON_RAPHSON, STEEPEST_DESCENT

  private
  
CONTAINS
  !
  !
  !
  function optimizeFunction(objectiveFunction,dim,x0,f0,f1,algorithm_input) result(x)
    use mod_debugTools, only: out_performance, deb_optiFun 
    implicit none
    !
    interface
      function objectiveFunction(dim,x) result(f)
        use def_kintyp_basic,       only: ip,rp
        integer(ip), intent(in)  :: dim
        real(rp),    intent(in)  :: x(dim)
        real(rp) :: f
     end function objectiveFunction
    end interface
    !
    integer(ip), intent(in)  :: dim
    real(rp),    intent(in)  :: x0(dim)
    real(rp),    intent(out) :: f0,f1
    integer(ip), intent(in), optional :: algorithm_input
    !
    real(rp)    :: x(dim)
    integer(ip) :: algorithm
    !
    real(rp) :: t0,t1
    !
    if(out_performance) call cpu_time(t0)
    !
    x = x0
    
    if(present(algorithm_input)) then
      algorithm = algorithm_input
    else
      algorithm = typeFunMin_default
    end if
    
    select case (algorithm)
    case (NEWTON_RAPHSON)
      call   newtonRaphson_numerical(objectiveFunction,dim,x,f0,f1)
    case (STEEPEST_DESCENT)
      call steepestDescent_numerical(objectiveFunction,dim,x,f0,f1)
    case default
      call runend('not implemented optimization procedure')
    end select
    
    if(out_performance) then
      call cpu_time(t1)
      deb_optiFun = deb_optiFun + (t1-t0)
    end if
    
    return
  end function optimizeFunction
  !
  !
  !
  subroutine setToleranceOptimization(tol_inp)
    implicit none
    
    real(rp), intent(in) :: tol_inp
    
    tolOptimization = tol_inp
    
    return
  end subroutine setToleranceOptimization
  !
  !
  !
  function optimizeFunction_analytical(OF,DOF,HOF,dim,x0,f0,f1,algorithm_input) result(x)
    implicit none
    !
    interface
      function OF(dim,x)  result(f) 
        use def_kintyp_basic, only: ip,rp
        integer(ip), intent(in)  :: dim
        real(rp),    intent(in)  :: x(dim)
        real(rp) :: f
     end function OF
    end interface
    !
    interface
      function DOF(dim,x) result(df)
        use def_kintyp_basic, only: ip,rp
        integer(ip), intent(in)  :: dim
        real(rp),    intent(in)  :: x(dim)
        real(rp) :: df(dim)
     end function DOF
    end interface
    !
    interface
      function HOF(dim,x) result(hf)
        use def_kintyp_basic, only: ip,rp
        integer(ip), intent(in)  :: dim
        real(rp),    intent(in)  :: x(dim)
        real(rp) :: hf(dim,dim)
     end function HOF
    end interface
    !
    integer(ip), intent(in)  :: dim
    real(rp),    intent(in)  :: x0(dim)
    real(rp),    intent(out) :: f0,f1
    integer(ip), intent(in), optional :: algorithm_input
    !
    real(rp)    :: x(dim)
    !
    integer(ip) :: algorithm
    !
    
    x = x0
    
    if(present(algorithm_input)) then
      algorithm = algorithm_input
    else
      algorithm = NEWTON_RAPHSON
    end if
    
    select case (algorithm)
    case (NEWTON_RAPHSON)
      call newtonRaphson_analytical(OF,DOF,HOF,dim,x,f0,f1)
    case default
      call runend('not implemented optimization procedure')
    end select
    
    
    return
  end function optimizeFunction_analytical
  !
  !
  !
  subroutine newtonRaphson_numerical(objectiveFunction,dim,x,f0,f1)
    !************************************************************************
    !**** Newton-Raphson minimization procedure
    !************************************************************************
    use mod_numDer,     only: numericalDerivative_centered, numericalHessian_centered, numericalHessian_forward
    use mod_debugTools, only: out_performance, deb_FunMin
    implicit none
    !*** input-output variables
    interface
      function objectiveFunction(dim,x) result(f)
        use def_kintyp_basic,       only: ip,rp
        integer(ip), intent(in)  :: dim
        real(rp),    intent(in)  :: x(dim)
        real(rp) :: f
     end function objectiveFunction
    end interface
    
    integer(ip), intent(in)  :: dim
    real(rp), intent(inout) :: x(dim)
    real(rp), intent(out)  :: f0,f1
    !*** inner variables
    integer(ip) :: iter, iterMaxNR
    real(rp)  :: f, fnew, fold, relErr
    real(rp)  :: H(dim,dim), Hinv(dim,dim), df(dim), pk(dim)
    logical(lg)  :: isNotFinished, isConverged, isIterViolated
    logical(lg)  :: errorflag
    real(rp) :: t0,t1
    !
    if(out_performance) call cpu_time(t0)
  
    iterMaxNR   = 1000_ip
    isNotFinished  = .true.
  
    iter = 0_ip
    do while( isNotFinished )
      
      if(typeNumerDer==DER_CENTERED) then
        f =objectiveFunction(                             dim,x)
        df=numericalDerivative_centered(objectiveFunction,dim,x)
        H =numericalHessian_centered(   objectiveFunction,dim,x)
      else
        call numericalHessian_forward(  objectiveFunction,dim,x,f,df,H)
      end if
      
      fold =f
      if(iter.eq.0)  f0 = f
      
      Hinv = invertMatrix(dim,H,errorflag)
      
      if(errorflag) then
        f1 = f
        return
      else
        pk = -matmul(Hinv,df)
        
        if( checkNansInVector(dim,pk) ) then
          f1 = f
          return
        end if

        call backtrackingLineSearch(objectiveFunction,dim,x,f,df,pk) !x updated
                
        fnew = f ! backtracking returns f as the new value
        relErr = abs(fold-fnew)/max(fold,fnew)!1e-6)
        isConverged     = relErr<tolOptimization
        isIterViolated  = iter>iterMaxNR
        isNotFinished   = (.not.isConverged).and.(.not.isIterViolated)
        iter = iter+1
        
        f1 = fnew
      end if
    end do
  
    if(out_performance) then
      call cpu_time(t1)
      deb_FunMin = deb_FunMin + (t1-t0)
    end if
    
    return
  end subroutine newtonRaphson_numerical
  !
  !
  !
  subroutine newtonRaphson_analytical(objectiveFunction,DOF,HOF,dim,x,f0,f1)
    !************************************************************************
    !**** Newton-Raphson minimization procedure
    !************************************************************************
    implicit none
    !*** input-output variables
    !
    interface
      function objectiveFunction(dim,x) result(f)
        use def_kintyp_basic, only: ip,rp
        integer(ip), intent(in)  :: dim
        real(rp),    intent(in)  :: x(dim)
        real(rp) :: f
     end function objectiveFunction
    end interface
    !
    interface
      function DOF(dim,x) result(df)
        use def_kintyp_basic, only: ip, rp
        integer(ip), intent(in)  :: dim
        real(rp),    intent(in)  :: x(dim)
        real(rp) :: df(dim)
     end function DOF
    end interface
    !
    interface
      function HOF(dim,x) result(hf)
        use def_kintyp_basic, only: ip, rp
        integer(ip), intent(in)  :: dim
        real(rp),    intent(in)  :: x(dim)
        real(rp) :: hf(dim,dim)
     end function HOF
    end interface
    !
    integer(ip), intent(in)  :: dim
    real(rp), intent(inout) :: x(dim)
    real(rp), intent(out)  :: f0,f1
    !*** inner variables
    integer(ip) :: iter, iterMaxNR
    real(rp)  :: f, fnew, fold, relErr
    real(rp)  :: H(dim,dim), Hinv(dim,dim), df(dim), pk(dim)
    logical  :: isNotFinished, errorflag
  
    iterMaxNR   = 1000
    isNotFinished  = .true.
  
    iter = 0 
    do while( isNotFinished )
      
      f  = objectiveFunction(dim,x)
      df = DOF(dim,x)
      H  = HOF(dim,x)
      
      fold =f
      if(iter.eq.0)  f0 = f
      
      Hinv = invertMatrix(dim,H,errorflag)
      
      if(errorflag) then
        call runend('OptimizationModule: not possible to invert matrix in analytical mode')
      end if
      
      pk = -matmul(Hinv,df)
      
      if( checkNansInVector(dim,pk) ) then
        f1 = f
        return
      end if
      
      call backtrackingLineSearch(objectiveFunction,dim,x,f,df,pk) !x updated
      
      fnew = f ! backtracking returns f as the new value
      relErr = abs(fold-fnew)/max(fold,fnew)!1e-6)
      isNotFinished = (relErr>tolOptimization).and.(iter.le.iterMaxNR)
      iter = iter+1
      
      f1 = fnew
    end do

    return
  end subroutine newtonRaphson_analytical
  !
  !
  !
  subroutine steepestDescent_numerical(objectiveFunction,dim,x,f0,f1)
    !************************************************************************
    !**** Newton-Raphson minimization procedure
    !************************************************************************
    use mod_numDer, only: numericalDerivative_centered, numericalDerivative_forward
    use mod_debugTools, only: out_performance, deb_FunMin
    implicit none
    !*** input-output variables
    interface
      function objectiveFunction(dim,x) result(f)
        use def_kintyp_basic,       only: ip,rp
        integer(ip), intent(in)  :: dim
        real(rp),    intent(in)  :: x(dim)
        real(rp) :: f
     end function objectiveFunction
    end interface
    
    integer(ip), intent(in)  :: dim
    real(rp), intent(inout) :: x(dim)
    real(rp), intent(out)  :: f0,f1
    !*** inner variables
    integer(ip) :: iter, iterMaxNR
    real(rp)  :: detM, detInv, f, fnew, fold, relErr, auxDelta
    real(rp)  :: H(dim,dim), Hinv(dim,dim), df(dim), pk(dim)
    logical(lg)  :: isNotFinished, isConverged, isIterViolated
    !logical(lg)  :: errorflag
    real(rp) :: t0,t1
    !
    if(out_performance) call cpu_time(t0)
  
    iterMaxNR   = 1000_ip
    isNotFinished  = .true.
  
    iter = 0_ip
    do while( isNotFinished )
      
      if(typeNumerDer==DER_CENTERED) then
        f =objectiveFunction(                             dim,x)
        df=numericalDerivative_centered(objectiveFunction,dim,x)
      else
        call numericalDerivative_forward(  objectiveFunction,dim,x,f,df)
      end if
      if(iter.eq.0)  f0 = f
      fold =f
            
      pk = -df
      call backtrackingLineSearch(objectiveFunction,dim,x,f,df,pk) !x updated
              
      fnew = f ! backtracking returns f as the new value
      relErr = abs(fold-fnew)/max(fold,fnew)!1e-6)
      isConverged     = relErr<tolOptimization
      isIterViolated  = iter>iterMaxNR
      isNotFinished   = (.not.isConverged).and.(.not.isIterViolated)
      iter = iter+1
      
      f1 = fnew
    end do
  
    if(out_performance) then
      call cpu_time(t1)
      deb_FunMin = deb_FunMin + (t1-t0)
    end if
    
    return
  end subroutine steepestDescent_numerical
  !
  !
  !
  function checkNansInVector(dim,v) result(areThereNans)
    implicit none
    integer(ip), intent(in) :: dim
    real(rp), intent(in) :: v(dim)
    
    integer(ip) :: i
    logical(lg) :: areThereNans
    
    areThereNans = .false.
    
    loopDim: do i = 1,dim
    
      !if(isnan(v(i))) then
      if(v(i).ne.v(i)) then
        areThereNans = .true.
        exit loopDim
      end if
      
    end do loopDim
    
    
    return
  end function checkNansInVector
  !
  !
  !
  subroutine backtrackingLineSearch(objectiveFunction,dim,x,f,df,pk)
    !************************************************************************
    !**** Backtracking line search for  distortion minimization
    !************************************************************************
    use mod_debugTools, only: out_performance, deb_BLS
    implicit none
    !*** input-output variables
    interface
      function objectiveFunction(dim,x) result(f)
        use def_kintyp_basic, only: ip,rp
        integer(ip), intent(in)  :: dim
        real(rp),    intent(in)  :: x(dim)
        real(rp) :: f
     end function objectiveFunction
    end interface
    !
    integer(ip), intent(in)     :: dim
    real(rp),    intent(inout)  :: x(dim)
    real(rp),    intent(inout)  :: f
    real(rp),    intent(in)     :: df(dim)
    real(rp),    intent(in)     :: pk(dim)
    !*** inner variables
    integer(ip) :: iter, iterMaxBLS
    real(rp)   :: factu,x0(dim),alpha,constant
    logical(lg) :: isNotFinished
    logical(lg) :: isIterViolated, isConverged
    !
    real(rp) :: t0,t1
    ! method parameters
    real(rp), parameter :: rho  = 0.5_rp
    real(rp), parameter :: c   = 1.0e-4_rp
    constant=c*sum(pk*df)
    alpha = 1.0_rp
    !
    if(out_performance) call cpu_time(t0)
    
    !  loop variables
    iterMaxBLS = 50_ip
    isNotFinished = .true.
  
    x0 = x
    iter = 0 
    do while( isNotFinished )
      
      x=x0+pk*alpha
      factu = objectiveFunction(dim,x)

      isIterViolated = iter>=iterMaxBLS
      isConverged    = factu <= f + alpha*constant
      isNotFinished = ( .not.isIterViolated ).and.( .not.isConverged )
      alpha = alpha*rho
      iter = iter+1
    end do

    if(isIterViolated) then
      x = x0
      if(out_debug_text) then
        print*,"     iter backtrack: ",iter," from max: ",iterMaxBLS
        print*,"     fpre : ",f
        print*,"     factu: ",factu
        print*,"     constant: ",constant
        print*,"     alpha: ",alpha
        print*,"     x    : ",x
        print*,"     x0   : ",x0
      end if
    else
      f = factu
    end if
  
    if(out_performance) then
      call cpu_time(t1)
      deb_BLS = deb_BLS + (t1-t0)
    end if
    
  end subroutine backtrackingLineSearch
  !
  !
  !
  function determinant3D(M) result(det)
    implicit none
    real(rp), intent(in) :: M(3,3)
    real(rp) :: a,b,c,d,e,f,g,h,i,det

    a = M(1,1)
    b = M(1,2)
    c = M(1,3)
    d = M(2,1)
    e = M(2,2)
    f = M(2,3)
    g = M(3,1)
    h = M(3,2)
    i = M(3,3)

    det = (a*e*i - a*f*h - b*d*i + b*f*g + c*d*h - c*e*g)

    return
  end function determinant3D
  !
  !
  !
  function invertMatrix(dim,M,errorFlag) result(Minv)
    implicit none
    integer(ip), intent(in) :: dim
    real(rp), intent(in) :: M(dim,dim)
    real(rp) :: Minv(dim,dim)
    logical(lg), INTENT(OUT) :: errorflag !Return error status. -1 for error, 0 for normal
    
    if(size(M,1)==3_ip) then
      Minv = invMatrix3(M,errorFlag)
    else if(size(M,1)==2_ip) then
      Minv = invMatrix2(M,errorFlag)
    else
      print*,"Not implemented matrix inversion for this dimension"
      call runend("Not implemented matrix inversion for this dimension")
    end if
    
    return
  end function
  !
  !
  !
  function invMatrix3(M,errorFlag) result(Minv)
    implicit none
    real(rp), intent(in) :: M(3,3)
    real(rp) :: Minv(3,3)
    real(rp) :: a,b,c,d,e,f,g,h,i,det
    logical(lg), intent(out) :: errorFlag

    a = M(1,1)
    b = M(1,2)
    c = M(1,3)
    d = M(2,1)
    e = M(2,2)
    f = M(2,3)
    g = M(3,1)
    h = M(3,2)
    i = M(3,3)

    det = (a*e*i - a*f*h - b*d*i + b*f*g + c*d*h - c*e*g)

    if(abs(det)<1e-14_rp) then
      Minv = 0.0_rp
      errorFlag = .true.
    else
      errorFlag = .false.
    end if

    Minv(1,1) =  (e*i - f*h)/det
    Minv(1,2) = -(b*i - c*h)/det
    Minv(1,3) =  (b*f - c*e)/det
    Minv(2,1) = -(d*i - f*g)/det
    Minv(2,2) =  (a*i - c*g)/det
    Minv(2,3) = -(a*f - c*d)/det
    Minv(3,1) =  (d*h - e*g)/det
    Minv(3,2) = -(a*h - b*g)/det
    Minv(3,3) =  (a*e - b*d)/det

    return
  end function invMatrix3
  !
  !
  !
  function invMatrix2(M,errorFlag) result(Minv)
    implicit none
    real(rp),    intent(in)  :: M(2,2)
    logical(lg), intent(out) :: errorFlag
    real(rp) :: Minv(2,2)

    real(rp) :: det

    det = M(1,1)*M(2,2) - M(1,2)*M(2,1)

    if(abs(det)<1e-14_rp) then
      Minv = 0.0
      errorFlag = .true.
    else
      errorFlag = .false.
    end if

    Minv(1,1) =  M(2,2)/det
    Minv(1,2) = -M(1,2)/det
    Minv(2,1) = -M(2,1)/det
    Minv(2,2) =  M(1,1)/det
    
    return
  end function invMatrix2  
  !
  !
  !
END MODULE mod_optimization
