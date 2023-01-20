!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Kermod
!> @{
!> @file    mod_ker_regularization.f90
!> @author  Miguel Zavala
!> @date    2018-12-29
!> @brief   Do not know
!> @details No idea...
!>         
!-----------------------------------------------------------------------


module mod_ker_regularization
  use def_kintyp, only : rp, ip

  implicit none
  save 
  real(rp),    private,parameter  :: delta = 1.0_rp
  real(rp),    private,parameter  :: trasl = 9.0_rp
  real(rp),    private,parameter  :: compr = 20.0_rp
  integer(ip), public, parameter  :: reg_identity = -1
  integer(ip), public, parameter  :: reg_garantzha = 0
  integer(ip), public, parameter  :: reg_exp = 1
  integer(ip), public, parameter  :: reg_exp_linear = 2
  integer(ip), public, parameter  :: reg_exp_quadratic = 3
  integer(ip), public, parameter  :: reg_exp_garantzha = 4
  integer(ip), public, parameter  :: reg_exp_traslation = 5
  integer(ip), public, parameter  :: reg_exp_compression = 6
  integer(ip), public, parameter  :: reg_exp_quadratic_comp   = 7
  integer(ip), public, parameter  :: reg_exp_log   = 8 !exponential-cubic-log
  integer(ip), public, parameter  :: reg_exp_cub   = 9 !exponential-cubic-linear
  integer(ip), public             :: reg_type  ! 0: Garantzha, 1: exp, 2: exp and x +1, 3: exp and 1+x^2/2
  logical    , public             :: kfl_regularization
  logical    , public             :: kfl_second
  ! cubic function parameters
  ! functions f1= exp(x), f2= 1+x+0.5*x^2+d_cub*x^3, f3 = g_cub + h_cub*x
  real(rp),    private,parameter  :: x_cub = 2.0_rp  
  real(rp),    private            :: g_cub, h_cub, d_cub
!  real(rp),    private            :: min_val=1.0e-11_rp
  real(rp),    private            :: min_k =0.0_rp
  real(rp),    private            :: min_e =0.0_rp

! functions  
  public  ::  regul_k,  regul_e, inv_regul_k, inv_regul_e
  public  ::  dregularization, d2regularization
  public  ::  regularization, inv_regularization
  private ::  fgarantzha, dgarantzha, d2garantzha  
!  private ::  fexp, dexp, d2exp 
  private ::  fexp_linear, dexp_linear,d2exp_linear 
  private ::  fexp_quadratic, dexp_quadratic, d2exp_quadratic
  private ::  fexp_log, dexp_log, d2exp_log
  private ::  fexp_cub, dexp_cub, d2exp_cub
  private ::  inv_garantzha, inv_exp_linear, inv_exp_quadratic, inv_exp_log, inv_exp_cub
  public  ::  set_regularization_type
contains


  subroutine  set_regularization_type(type) 
    implicit none
    integer(ip), intent(in)  :: type
    
    reg_type = type
    
    print *, 'regul type=', type
    return
  end subroutine  set_regularization_type
  subroutine set_regularization_pars()
    implicit none

    select case (reg_type)
    case(reg_exp_cub)
       ! functions f1= exp(x), f2= 1+x+0.5*x^2+d_cub*x^3, f3 = g_cub + h_cub*x
       d_cub = -1.0_rp/(6.0_rp*x_cub) !continuity of second derivatives at x_cub
       h_cub = 1.0_rp + 0.5_rp*x_cub  !continuity of first derivatives at x_cub
       g_cub = 1.0_rp + (1.0_rp - h_cub)*x_cub + 0.5_rp*x_cub*x_cub &
            + d_cub*(x_cub**3.0_rp)   !continuity of functions at x_cub
  !    print *,  d_cub, h_cub, g_cub
    end select
  end subroutine set_regularization_pars
  
  function regul_k(x) result(value)
    implicit none
    real(rp), intent(in) :: x
    real(rp)             :: value

    value =regularization(x) + min_k
    return

  end function regul_k

  function regul_e(x) result(value)
    implicit none
    real(rp), intent(in) :: x
    real(rp)             :: value
    value =regularization(x) + min_e
    return
  end function regul_e

  function inv_regul_k(f) result(value)
    implicit none
    real(rp), intent(in) :: f
    real(rp)             :: value, f_in

    f_in = max(10.0_rp*epsilon(f), f-min_k)
    
    value = inv_regularization(f_in)
    
    return
  end function inv_regul_k

  function inv_regul_e(f) result(value)
    implicit none
    real(rp), intent(in) :: f
    real(rp)             :: value, f_in

    f_in = max(10.0_rp*epsilon(f), f-min_e)
    
    value = inv_regularization(f_in)
    
    return
  end function inv_regul_e
  
  function regularization(x)  result( value)

    implicit none
    real(rp), intent(in) :: x
    real(rp)             :: value

    select case (reg_type)
    case (reg_garantzha)
       value = fgarantzha(x)  
    case(reg_exp)
       value = exp(x)
    case(reg_exp_linear)
       value = fexp_linear(x)
    case(reg_exp_quadratic)
       value = fexp_quadratic(x)  ! fexp_quadratic
    case(reg_exp_garantzha) 
       value = fexp_garantzha(x)
    case(reg_exp_traslation) 
       value = fexp_traslation(x)
    case(reg_exp_compression) 
       value = fexp_compression(x)
    case(reg_exp_quadratic_comp)
       value = fexp_quadratic(compr*x)  ! fexp_quadratic
    case(reg_exp_log)
       value = fexp_log(x)  
    case(reg_exp_cub)
       value = fexp_cub(x)  
    case(reg_identity)
       value = x
  !  case  default
  !     call runend('')
    end select
    
    return
    
  end function regularization

  function inv_regularization(x)  result( value)

    implicit none
    real(rp), intent(in) :: x
    real(rp)             :: value

    select case (reg_type)
    case (reg_garantzha)
       value = inv_garantzha(x)  
    case(reg_exp)
       value = log(x)
    case(reg_exp_linear)
       value = inv_exp_linear(x)
    case(reg_exp_quadratic)
       value = inv_exp_quadratic(x)
    case (reg_exp_garantzha) 
       value = inv_exp_garantzha(x)
    case(reg_exp_traslation) 
       value = inv_exp_traslation(x)
    case(reg_exp_compression) 
       value = inv_exp_compression(x)
    case(reg_exp_quadratic_comp)
       value = inv_exp_quadratic(x)/compr
    case(reg_exp_log)
       value = inv_exp_log(x)
    case(reg_exp_cub)
       value = inv_exp_cub(x)     
    case(reg_identity)
       value = x
  !  case  default
  !     call runend('')
    end select
    
    return
    
  end function inv_regularization
  
function dregularization(x)  result( value)

    implicit none
    real(rp), intent(in) :: x
    real(rp)             :: value

    select case (reg_type)
    case (reg_garantzha)
       value = dgarantzha(x)  
    case(reg_exp)
       value = exp(x)
    case(reg_exp_linear)
       value = dexp_linear(x)
    case(reg_exp_quadratic)
       value = dexp_quadratic(x)
    case (reg_exp_garantzha) 
       value = dexp_garantzha(x)
    case(reg_exp_traslation) 
       value = dexp_traslation(x)
    case(reg_exp_compression) 
       value = dexp_compression(x)
    case(reg_exp_quadratic_comp)
       value = compr*dexp_quadratic(compr*x)  
    case(reg_exp_log)
       value = dexp_log(x)
    case(reg_exp_cub)
       value = dexp_cub(x)  
    case(reg_identity)
       value = 1.0_rp
  !  case  default
  !     call runend('')
    end select
    
    return
    
  end function dregularization
  
  function d2regularization(x)  result( value)

    implicit none
    real(rp), intent(in) :: x
    real(rp)             :: value
    
    select case (reg_type)
    case (reg_garantzha)
       value = d2garantzha(x)  
    case(reg_exp)
       value = exp(x)
    case(reg_exp_linear)
       value = d2exp_linear(x)
    case(reg_exp_quadratic)
       value = d2exp_quadratic(x)
    case (reg_exp_garantzha) 
       value = d2exp_garantzha(x)
    case(reg_exp_traslation) 
       value = d2exp_traslation(x)
    case(reg_exp_compression) 
       value = d2exp_compression(x)
    case(reg_exp_quadratic_comp)
       value = compr*compr*d2exp_quadratic(compr*x)
    case(reg_exp_log)
       value = d2exp_log(x)
    case(reg_exp_cub)
       value = d2exp_cub(x)  
    case(reg_identity)
       value = 0.0_rp
  !  case  default
  !     call runend('')
    end select
    
    return
    
  end function d2regularization

  function fgarantzha(x) result (func)
    implicit none
    real(rp), intent(in) :: x
    real(rp)             :: func
    
    func = 0.5_rp*(x + sqrt(x*x +4.0_rp*delta*delta))
    
  end function fgarantzha

  function inv_garantzha(x) result (func)
    implicit none
    real(rp), intent(in) :: x
    real(rp)             :: func

    func = x - delta*delta/x
  end function inv_garantzha
  function dgarantzha(x) result (dfunc)
    implicit none
    real(rp), intent(in) :: x
    real(rp)             :: dfunc
    
    dfunc = 0.5_rp*(1.0_rp  + x/sqrt(x*x +4.0_rp*delta*delta))
    
  end function dgarantzha
  function d2garantzha(x) result (d2func)
    implicit none
    real(rp), intent(in) :: x
    real(rp)           :: d2func

    d2func = 2.0_rp*delta*delta/(x*x +4.0_rp*delta*delta)**1.5_rp
    
  end function d2garantzha
  
  function fexp_linear(x) result (func)
    implicit none
    real(rp), intent(in) :: x
    real(rp)             :: func
    if (x<0.0_rp) then
       func = exp(x)
    else
       func = x+1.0_rp
    end if
    
  end function fexp_linear

  function inv_exp_linear(x) result (func)
    implicit none
    real(rp), intent(in) :: x
    real(rp)             :: func

    if (x<1.0_rp) then
       func = log(x)
    else
       func = x - 1.0_rp
    end if
    
  end function inv_exp_linear
  
  function dexp_linear(x) result (dfunc)
    implicit none
    real(rp), intent(in) :: x
    real(rp)             :: dfunc

    if (x<0.0_rp) then
       dfunc = exp(x)
    else
       dfunc = 1.0_rp
    end if
        
  end function dexp_linear

  function d2exp_linear(x) result (d2func)
    implicit none
    real(rp), intent(in) :: x
    real(rp)             :: d2func

    if (x<0.0_rp) then
       d2func = exp(x)
    else
       d2func = 0.0_rp
    end if   
    
  end function d2exp_linear
  
   function fexp_quadratic(x) result (func)
    implicit none
    real(rp), intent(in) :: x
    real(rp)             :: func
    if (x<0.0_rp) then
       func = exp(x)
    else
       func = 0.5_rp*x*x+ x + 1.0_rp
    end if
    
  end function fexp_quadratic

  function inv_exp_quadratic(x) result (func)
    implicit none
    real(rp), intent(in) :: x
    real(rp)             :: func

    if (x<1.0_rp) then
       func = log(x)
    else
       func = - 1.0_rp + sqrt(2.0_rp*x - 1.0_rp )
    end if
    
  end function inv_exp_quadratic
  
  function dexp_quadratic(x) result (dfunc)
    implicit none
    real(rp), intent(in) :: x
    real(rp)             :: dfunc

    if (x<0.0_rp) then
       dfunc = exp(x)
    else
       dfunc = x + 1.0_rp
    end if
        
  end function dexp_quadratic

  function d2exp_quadratic(x) result (d2func)
    implicit none
    real(rp), intent(in) :: x
    real(rp)             :: d2func

    if (x<0.0_rp) then
       d2func = exp(x)
    else
       d2func = 1.0_rp
    end if   
    
  end function d2exp_quadratic
  
  function fexp_garantzha(x) result (func)
    implicit none
    real(rp), intent(in) :: x
    real(rp)             :: func
    if (x<0.0_rp) then
       func = exp(x)
    else
       func =x + sqrt(x*x +1.0_rp)
    end if
    
  end function fexp_garantzha

  function inv_exp_garantzha(x) result (func)
    implicit none
    real(rp), intent(in) :: x
    real(rp)             :: func

    if (x<1.0_rp) then
       func = log(x)
    else
       func =0.5_rp*(x*x -1.0_rp)/x
    end if
    
  end function inv_exp_garantzha
  
  function dexp_garantzha(x) result (dfunc)
    implicit none
    real(rp), intent(in) :: x
    real(rp)             :: dfunc

    if (x<0.0_rp) then
       dfunc = exp(x)
    else
       dfunc = 1.0_rp + x/sqrt(x*x + 1.0_rp)
    end if
        
  end function dexp_garantzha

  function d2exp_garantzha(x) result (d2func)
    implicit none
    real(rp), intent(in) :: x
    real(rp)             :: d2func

    if (x<0.0_rp) then
       d2func = exp(x)
    else
       d2func = 1.0_rp/((x*x +1.0_rp)**1.5_rp)
    end if   
    
  end function d2exp_garantzha

 function fexp_traslation(x) result (func)
    implicit none
    real(rp), intent(in) :: x
    real(rp)             :: func
   
    func = exp(x-trasl)
    
  end function fexp_traslation

  function inv_exp_traslation(x) result (func)
    implicit none
    real(rp), intent(in) :: x
    real(rp)             :: func

    
    func = trasl + log(x)
    
  end function inv_exp_traslation
  
  function dexp_traslation(x) result (dfunc)
    implicit none
    real(rp), intent(in) :: x
    real(rp)             :: dfunc

    dfunc = exp(x-trasl)
       
  end function dexp_traslation

  function d2exp_traslation(x) result (d2func)
    implicit none
    real(rp), intent(in) :: x
    real(rp)             :: d2func

    d2func = exp(x - trasl)
    
  end function d2exp_traslation

 function fexp_compression(x) result (func)
    implicit none
    real(rp), intent(in) :: x
    real(rp)             :: func
    
    func = exp(compr*x)
    
  end function fexp_compression

  function inv_exp_compression(x) result (func)
    implicit none
    real(rp), intent(in) :: x
    real(rp)             :: func

    func = log(x)/compr
    
  end function inv_exp_compression
  
  function dexp_compression(x) result (dfunc)
    implicit none
    real(rp), intent(in) :: x
    real(rp)             :: dfunc

    dfunc = compr*exp(compr*x)
        
  end function dexp_compression

  function d2exp_compression(x) result (d2func)
    implicit none
    real(rp), intent(in) :: x
    real(rp)             :: d2func

    d2func = compr*compr*exp(compr*x)
    
  end function d2exp_compression
  
 function fexp_log(x) result (func)
    implicit none
    real(rp), intent(in) :: x
    real(rp)             :: func
    if (x.lt.0.0_rp) then
       func = exp(x)
    else if (x.lt.1.0_rp) then
       func = 0.5_rp*x*x+ x + 1.0_rp - x*x*x/3.0_rp
    else
       func = 13.0_rp/6.0_rp + log(x)
    end if
    
  end function fexp_log

  function inv_exp_log(f) result (x)
    implicit none
    real(rp), intent(in) :: f
    real(rp)             :: x
    real(rp)             :: x_i, err, toler
    integer (ip)         :: iiter 

    if (f.lt.1.0_rp) then
       x = log(f)
    else if(f.lt.13.0_rp/6.0_rp) then
       ! look for a cero (x) of function  = 1- f + x - 0.5*x^2 -x^3/3
       x = 0.5_rp ! initial guess (0<x<1)
       x_i = 0.0_rp
       iiter = 0
       err = 10.0_rp
       toler = 10.0_rp*epsilon(x)
       do while(err.gt.toler.and.iiter.lt.100)
          x_i = x
          x = x_i - (fexp_log(x_i)-f)/dexp_log(x_i)
          iiter= iiter +1   
          err = max(abs(x_i-x), abs(fexp_log(x)-f))          
       end do
       if (iiter==100) call runend('mod_ker_regularization:inv function did not converge')
    else
       x = exp(f-13.0_rp/6.0_rp)
    end if
    
  end function inv_exp_log
  
  function dexp_log(x) result (dfunc)
    implicit none
    real(rp), intent(in) :: x
    real(rp)             :: dfunc

    if (x.lt.0.0_rp) then
       dfunc = exp(x)
    else if (x.lt.1.0_rp) then
       dfunc = x + 1.0_rp - x*x ! maximum = 1.25@x=0.5
    else 
       dfunc = 1.0_rp/x
    end if
        
  end function dexp_log

  function d2exp_log(x) result (d2func)
    implicit none
    real(rp), intent(in) :: x
    real(rp)             :: d2func

    if (x.lt.0.0_rp) then
       d2func = exp(x)
    else if(x.lt.1.0_rp) then
       d2func = 1.0_rp -2.0_rp*x
    else
       d2func = -1.0_rp/(x*x)
    end if   
    
  end function d2exp_log

 function fexp_cub(x) result (func)
    implicit none
    real(rp), intent(in) :: x
    real(rp)             :: func

    if (x.lt.0.0_rp) then
       func = exp(x)
    else if (x.lt.x_cub) then
       func = 0.5_rp*x*x+ x + 1.0_rp +d_cub* x*x*x 
    else
       func =g_cub + h_cub*x 
    end if
    
!    func = func + min_val
  end function fexp_cub

  function inv_exp_cub(f) result (x)
    implicit none
    real(rp), intent(in) :: f
    real(rp)             :: x
    real(rp)             :: x_i, err, toler
    integer (ip)         :: iiter 

!    f_in =max(min_val +10.0*epsilon(f), f)
    
    if (f.lt.1.0_rp) then
       x = log(f)
    else if (f.lt.(g_cub + h_cub*x_cub)) then
       ! look for a cero (x) of function  = 1- f + x - 0.5*x^2 -x^3/3
       x = x_cub*0.5_rp ! initial guess (0<x<x_c)
       x_i = 0.0_rp
       iiter = 0
       err = 10.0_rp
       toler = 10.0_rp*epsilon(x)
       do while(err.gt.toler.and.iiter.lt.100)
          x_i = x
          x = x_i - (fexp_cub(x_i)-f)/dexp_cub(x_i)
          iiter= iiter +1   
          err = max(abs(x_i-x), abs(fexp_cub(x)-f))          
!          print *, 'iiter, err', iiter, err
       end do
       if (iiter==100) call runend('mod_ker_regularization:inv function did not converge')
    else
       x = (f-g_cub)/h_cub
    end if
    
  end function inv_exp_cub
  
  function dexp_cub(x) result (dfunc)
    implicit none
    real(rp), intent(in) :: x
    real(rp)             :: dfunc

    if (x.lt.0.0_rp) then
       dfunc = exp(x)
    else if (x.lt.x_cub) then
       dfunc = x + 1.0_rp + 3.0_rp*d_cub*x*x 
    else 
       dfunc = h_cub
    end if
    
  end function dexp_cub

  function d2exp_cub(x) result (d2func)
    implicit none
    real(rp), intent(in) :: x
    real(rp)             :: d2func

    if (x.lt.0.0_rp) then
       d2func = exp(x)
    else if(x.lt.x_cub) then
       d2func = 1.0_rp + 6.0_rp*d_cub*x
    else
       d2func = 0.0_rp
    end if   
    
  end function d2exp_cub




end module mod_ker_regularization
!> @}
