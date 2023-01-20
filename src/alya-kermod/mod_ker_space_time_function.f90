!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup Function
!> @{
!> @name    ToolBox for space/time functions
!> @file    mod_ker_space_time_function.f90
!> @author  Guillaume Houzeaux
!> @date    22/02/2013
!> @brief   ToolBox for space/time functions
!> @details Allocate memory, parse formulas, etc.
!>          This module is based in fparser module from Roland Schmehl.
!>
!>          \verbatim
!>          !------- -------- --------- --------- --------- --------- --------- --------- -------
!>          ! Fortran 90 function parser v1.0
!>          !------- -------- --------- --------- --------- --------- --------- --------- -------
!>          !
!>          ! 
!>          !
!>          ! The source code is inspired by an old version of the code available here:
!>          ! https://github.com/jacobwilliams/fortran_function_parser
!>          !
!>          ! Please send comments, corrections or questions to the author:
!>          ! Roland Schmehl <Roland.Schmehl@mach.uni-karlsruhe.de>
!>          !
!>          !------- -------- --------- --------- --------- --------- --------- --------- -------
!>          ! The function parser concept is based on a C++ class library written by Warp
!>          ! <warp@iki.fi> available from:
!>          ! http://www.students.tut.fi/~warp/FunctionParser/fparser.zip
!>          !------- -------- --------- --------- --------- --------- --------- --------- -------
!>          \endverbatim
!>
!>          To add a function xxx, modify the following:
!>          1. Define variable number Cxxx
!>          2. Add function name in Funcs 'xxx  '
!>          3. Code your function case (Cxxx)
!> 
!> @{
!------------------------------------------------------------------------

module mod_ker_space_time_function

  use def_kintyp,         only : ip,rp,lg
  use def_master,         only : npari,nparr,nparc 
  use def_master,         only : parin,parre
  use def_master,         only : ISLAVE
  use def_master,         only : IMASTER
  use def_master,         only : igene
  use def_master,         only : mem_modul
  use def_master,         only : modul
  use def_kermod,         only : number_space_time_function,space_time_function
  use mod_memory,         only : memory_alloca
  use mod_memory,         only : memory_deallo
  use mod_communications, only : PAR_BROADCAST
  use fparser
  
  implicit none
  private
  save
  !
  ! Variable names and values
  !
  integer(ip),    parameter                  :: numbervariables=4
  character(10)                              :: variablenames(numbervariables)      
  real(rp)                                   :: variablesvalues(numbervariables)
  integer(ip),    pointer                    :: function_position(:)
  integer(ip)                                :: ker_n_functions 
  !
  ! Interface for space/time functions evaluation
  !
  interface ker_space_time_function
     module procedure ker_space_time_function_scalar,   &
          &           ker_space_time_function_scalar_1, &
          &           ker_space_time_function_vector,   &
          &           ker_space_time_function_vector_1
  end interface

  public :: ker_n_functions 
  public :: ker_init_space_time_function
  public :: ker_space_time_function
  public :: space_time_function_number
  public :: space_time_function_parall
  
contains

  
  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    25/02/2013
  !> @brief   Initialize space/time functions
  !> @details Allocate memory for FUNCTION_NUMBER functions
  !
  !----------------------------------------------------------------------

  subroutine ker_space_time_function_destructor()

    if( associated(function_position) ) then
       deallocate(function_position)
    end if
    
  end subroutine ker_space_time_function_destructor
  
  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    25/02/2013
  !> @brief   Initialize space/time functions
  !> @details Allocate memory for FUNCTION_NUMBER functions
  !
  !----------------------------------------------------------------------

  subroutine ker_init_space_time_function()
    integer(ip) :: ifunc
    integer(ip) :: idime
    integer(ip) :: function_number
    !
    ! Variables
    !
    variablenames(1) = 'x'
    variablenames(2) = 'y'
    variablenames(3) = 'z'
    variablenames(4) = 't'
    !
    ! Allocate 
    !
    if( number_space_time_function > 0 ) &
         allocate(function_position(number_space_time_function))
    !
    ! FUNCTION_NUMBER= number of functions
    !
    function_number = 0
    do ifunc = 1,number_space_time_function
       function_position(ifunc) = function_number
       function_number = function_number + size(space_time_function(ifunc) % expression)
    end do
    !
    ! Initialize function parser for function_numberfunctions
    !
    if( function_number > 0 ) call initf(function_number) 
    !
    ! Parse and bytecompile ifunc-th function string
    !
    function_number = 0 
    do ifunc = 1,number_space_time_function
       do idime = 1,size(space_time_function(ifunc) % expression)
          function_number = function_number + 1
          call parsef(function_number,space_time_function(ifunc) % expression(idime), variablenames)
       end do
    end do

   ker_n_functions = function_number  

  end subroutine ker_init_space_time_function

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    25/02/2013
  !> @brief   Get a function number
  !> @details Get a space and time function number given the function
  !>          name
  !
  !----------------------------------------------------------------------

  function space_time_function_number(wfname)
    integer(ip)                :: space_time_function_number
    character(*),  intent(in)  :: wfname
    integer(ip)                :: ifunc

    space_time_function_number = 0
    do ifunc = 1,number_space_time_function
       if( trim(wfname) == trim(space_time_function(ifunc) % name) ) then
          space_time_function_number = ifunc
       end if
    end do
    if( space_time_function_number == 0 ) &
         call runend('SPACE TIME FUNCTION '//trim(wfname)//' DOES NOT EXIST')

  end function space_time_function_number

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    25/02/2013
  !> @brief   Evaluate a space/time function
  !> @details Evaluate multidimensional arrays using space/time functions
  !
  !----------------------------------------------------------------------

  subroutine ker_space_time_function_vector(ifunc,x,y,z,t,xvalu)
    integer(ip),  intent(in)  :: ifunc
    real(rp),     intent(in)  :: x
    real(rp),     intent(in)  :: y
    real(rp),     intent(in)  :: z
    real(rp),     intent(in)  :: t
    real(rp),     intent(out) :: xvalu(:)
    integer(ip)               :: idime,nsize,nexpr
    integer(ip)               :: current_function

    variablesvalues(1) = x
    variablesvalues(2) = y
    variablesvalues(3) = z
    variablesvalues(4) = t
    !
    ! Look for function position 
    !
    current_function = function_position(ifunc)
    nsize = size(xvalu)
    nexpr = size(space_time_function(ifunc) % expression)
   
    if( nexpr == 1 ) then
       !
       ! Only one function for all degrees of freedom
       !
       current_function = current_function + 1
       do idime = 1,nsize
          xvalu(idime) = evalf(current_function,variablesvalues)
          if( EvalErrType > 0 ) write(*,*) 'Error was found while evaluating a space/time function: ',EvalErrMsg()
       end do
    else if( nexpr == nsize ) then
       !
       ! One function per degree of freedom
       !
       do idime = 1,nsize
          current_function = current_function + 1
          xvalu(idime) = evalf(current_function,variablesvalues)
          if( EvalErrType > 0 ) write(*,*) 'Error was found while evaluating a space/time function: ',EvalErrMsg()
       end do
    else
       !
       ! Wrong combination
       !
       write(*,*) 'Error space/time function dimension'
    end if

  end subroutine ker_space_time_function_vector

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    25/02/2013
  !> @brief   Evaluate a space/time function
  !> @details Evaluate multidimensional arrays using space/time functions
  !
  !----------------------------------------------------------------------

  subroutine ker_space_time_function_vector_1(ifunc,x,t,xvalu)
    integer(ip),  intent(in)  :: ifunc
    real(rp),     intent(in)  :: x(:)
    real(rp),     intent(in)  :: t
    real(rp),     intent(out) :: xvalu(:)
    integer(ip)               :: idime,nsize,nexpr
    integer(ip)               :: current_function,ndime,ndim2,ndim3

    ndime = size(x)
    ndim2 = min(2_ip,ndime)
    ndim3 = min(3_ip,ndime)
    variablesvalues(1) = x(1)
    variablesvalues(2) = x(ndim2)
    variablesvalues(3) = x(ndim3)
    variablesvalues(4) = t
    !
    ! Look for function position 
    !
    current_function = function_position(ifunc)
    nsize = size(xvalu)
    nexpr = size(space_time_function(ifunc) % expression)
   
    if( nexpr == 1 ) then
       !
       ! Only one function for all degrees of freedom
       !
       current_function = current_function + 1
       do idime = 1,nsize
          xvalu(idime) = evalf(current_function,variablesvalues)
          if( EvalErrType > 0 ) write(*,*) 'Error was found while evaluating a space/time function: ',EvalErrMsg()
       end do
    else if( nexpr == nsize ) then
       !
       ! One function per degree of freedom
       !
       do idime = 1,nsize
          current_function = current_function + 1
          xvalu(idime) = evalf(current_function,variablesvalues)
          if( EvalErrType > 0 ) write(*,*) 'Error was found while evaluating a space/time function: ',EvalErrMsg()
       end do
    else
       !
       ! Wrong combination
       !
       write(*,*) 'Error space/time function dimension'
    end if

  end subroutine ker_space_time_function_vector_1

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    25/02/2013
  !> @brief   Evaluate a space/time function
  !> @details Evaluate a scalar using space/time functions
  !
  !----------------------------------------------------------------------

  subroutine ker_space_time_function_scalar(ifunc,x,y,z,t,xvalu)
    integer(ip),  intent(in)  :: ifunc
    real(rp),     intent(in)  :: x
    real(rp),     intent(in)  :: y
    real(rp),     intent(in)  :: z
    real(rp),     intent(in)  :: t
    real(rp),     intent(out) :: xvalu
    integer(ip)               :: current_function

    variablesvalues(1) = x
    variablesvalues(2) = y
    variablesvalues(3) = z
    variablesvalues(4) = t

    current_function = function_position(ifunc) + 1
    xvalu = evalf(current_function,variablesvalues)
    if( EvalErrType > 0 ) write(*,*) 'Error was found while evaluating a space/time function: ',EvalErrMsg()

  end subroutine ker_space_time_function_scalar

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    25/02/2013
  !> @brief   Evaluate a space/time function
  !> @details Evaluate a scalar using space/time functions
  !
  !----------------------------------------------------------------------

  subroutine ker_space_time_function_scalar_1(ifunc,x,t,xvalu)
    integer(ip),  intent(in)  :: ifunc
    real(rp),     intent(in)  :: x(:)
    real(rp),     intent(in)  :: t
    real(rp),     intent(out) :: xvalu
    integer(ip)               :: current_function,ndime,ndim2,ndim3

    ndime = size(x)
    ndim2 = min(2_ip,ndime)
    ndim3 = min(3_ip,ndime)
    variablesvalues(1) = x(1)
    variablesvalues(2) = x(ndim2)
    variablesvalues(3) = x(ndim3)
    variablesvalues(4) = t

    current_function = function_position(ifunc) + 1
    xvalu = evalf(current_function,variablesvalues)
    if( EvalErrType > 0 ) write(*,*) 'Error was found while evaluating a space/time function: ',EvalErrMsg()

  end subroutine ker_space_time_function_scalar_1

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-01
  !> @brief   Parallelization
  !> @details Broadcast of data Parallelization
  !> 
  !-----------------------------------------------------------------------

  subroutine space_time_function_parall()

    integer(ip) :: ifunc,kdime,nexpr,parii,idime

    !call PAR_BROADCAST(number_space_time_function)                ! Space/Time functions
    !do ifunc = 1,size(space_time_function)
    !   call PAR_BROADCAST(space_time_function(ifunc) % ndime)
    !   call PAR_BROADCAST(space_time_function(ifunc) % nexpr)
    !   call PAR_BROADCAST(5_ip,space_time_function(ifunc) % name,'IN MY CODE')       
    !end do

    
    if( number_space_time_function > 0 ) then
       if( ISLAVE ) then
          do ifunc = 1,number_space_time_function
             igene = ifunc
             call ker_memory(7_ip)
          end do
       end if
       do parii = 1,2
          npari = 0
          nparr = 0
          nparc = 0
          do ifunc = 1,number_space_time_function
             nexpr = space_time_function(ifunc) % nexpr
             idime = space_time_function(ifunc) % ndime
             do kdime = 1,idime
                call cexcha(nexpr,space_time_function(ifunc) % expression(kdime))
             end do
          end do
          if( parii == 1 ) then
             call memory_alloca(mem_modul(1:2,modul),'PARIN','ker_parall',parin,npari,'DO_NOT_INITIALIZE')
             call memory_alloca(mem_modul(1:2,modul),'PARRE','ker_parall',parre,nparr,'DO_NOT_INITIALIZE')
             if( ISLAVE ) call par_broadc()
          end if
       end do
       if( IMASTER ) call par_broadc()
       call memory_deallo(mem_modul(1:2,modul),'PARRE','ker_parall',parre)
       call memory_deallo(mem_modul(1:2,modul),'PARIN','ker_parall',parin)
    end if

  end subroutine space_time_function_parall

end module mod_ker_space_time_function
