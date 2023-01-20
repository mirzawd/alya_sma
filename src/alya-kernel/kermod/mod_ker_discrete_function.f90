!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Kermod
!> @{
!> @file    mod_ker_discrete_functions.f90
!> @author  houzeaux
!> @date    2020-03-30
!> @brief   Apply discrete functions
!> @details Tools for discrete functions
!-----------------------------------------------------------------------

module mod_ker_discrete_function

  use def_kintyp,                  only : ip,rp
  use def_master
  use def_inpout
  use mod_memory,                  only : memory_alloca
  use mod_ecoute,                  only : ecoute
  use mod_exchange,                only : exchange_init
  use mod_exchange,                only : exchange_add
  use mod_exchange,                only : exchange_end
  implicit none
  private

  integer(ip), parameter  :: DISCRETE_FUNCTION_PIECEWISE_CONSTANT = 0
  integer(ip), parameter  :: DISCRETE_FUNCTION_LINEAR             = 1
  integer(ip), parameter  :: DISCRETE_FUNCTION_SMOOTH             = 2

  integer(ip)             :: num_discrete_functions

  type typ_discrete_function
     integer(ip)         :: kfl_type
     integer(ip)         :: kfl_periodic
     integer(ip)         :: dim
     integer(ip)         :: num_tim
     real(rp),   pointer :: tim(:)
     real(rp),   pointer :: val(:,:)
     character(5)        :: name
  end type typ_discrete_function
  type(typ_discrete_function), pointer :: fundi(:)

  interface ker_discrete_function
     module procedure ker_discrete_function_scalar,&
          &           ker_discrete_function_vector
  end interface ker_discrete_function

  public :: num_discrete_functions
  public :: ker_discrete_function
  public :: ker_discrete_function_parall
  public :: ker_discrete_function_initialization
  public :: ker_discrete_function_read
  public :: ker_discrete_function_number
  public :: ker_discrete_function_getdim
  
contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-05-29
  !> @brief   Read function
  !> @details Read discrete functions
  !> 
  !-----------------------------------------------------------------------

  subroutine ker_discrete_function_read()

    integer(ip) :: nauxi,ifunc,idime,istep
    real(rp)    :: rauxi(2)

    !----------------------------------------------------------
    !
    ! Functions
    !
    !----------------------------------------------------------
    !
    ! ADOC[d]> <ul> </ul>
    ! ADOC[d]> Functions:
    !
    ! ADOC[1]> DISCRETE_FUNCTIONS                                   $ Functions 
    ! ADOC[2]>   TOTAL_NUMBER: int1                                 $ Total number of functions
    ! ADOC[3]>   FUNCTION_NAME: char, NUMBER=int, DIMENSION=int     $ Function number
    ! ADOC[3]>     TIME_SHAPE: LINEAR , PIECEWISE, SMOOTH           $ Function type
    ! ADOC[3]>     SHAPE_DEFINITION
    ! ADOC[4]>       int3                                           $ Number of lines
    ! ADOC[4]>       rea1 rea2 rea3 rea4 ...                        $ rea1=time, rea2-4=DoFs
    ! ADOC[4]>       ...
    ! ADOC[3]>     END_SHAPE_DEFINITION
    ! ADOC[1]>   END_FUNCTION
    ! ADOC[1]> END_DISCRETE_FUNCTIONS
    ! ADOC[d]> <b>FUNCTIONS</b>:
    ! ADOC[d]> Space/time functions applied to those nodes/boundaries.
    ! ADOC[d]> <ul>
    ! ADOC[d]> <li>
    ! ADOC[d]>  LINEAR
    ! ADOC[d]>  Linear step (Ramp) is applied between discrete times.
    ! ADOC[d]> </li>
    ! ADOC[d]> <li>
    ! ADOC[d]>  SMOOTH
    ! ADOC[d]>  This function is such that first and second derivatives are zero between two discrete times.
    ! ADOC[d]>  This definition is intended to ramp up or down smoothly from one amplitude value to another.
    ! ADOC[d]> </li>
    ! ADOC[d]> </ul>

    if( exists('NUMBE') ) then
       num_discrete_functions = getint('FUNCT',1_ip,'#FUNCTION NUMBER')
    end if
    call ecoute('ker_discrete_function_read')

    ! Total number of functions
    if (words(1) == 'TOTAL') num_discrete_functions =  int(param(1))
    if( num_discrete_functions == 0 ) num_discrete_functions = 10
    call ker_dicsrete_function_allocate()

    if (num_discrete_functions == 0) then
       call runend('KER_DISCRETE_FUNCTION_READ: PROVIDE A TOTAL_NUMBER OF FUNCTIONS')
    else if (num_discrete_functions > 9_ip) then
       call runend('KER_DISCRETE_FUNCTION_READ: THE TOTAL_NUMBER OF FUNCTIONS MUST BE LOWER THAN 10')
    end if

    ! Initializations
    nauxi = 0
    ifunc = 0

    do while( words(1) /= 'ENDDI' )

       if( words(1) == 'FUNCT' ) then
          !
          ! Automatic function number
          !
          ifunc = ifunc + 1 
          if( exists('DIMEN') ) then
             idime = getint('DIMEN',1_ip,'#Dimension of the function')
          else
             idime = 1
          end if

          fundi(ifunc) % dim  = idime
          fundi(ifunc) % name = getcha('FUNCT','     ','#Space time function name')

          if( ifunc < 0 .or. ifunc > num_discrete_functions ) then
             call runend('KER_DISCRETE_FUNCTION_READ: WRONG FUNCION NUMBER IN FUNCTION '//trim(fundi(ifunc) % name))
          end if

          do while( words(1) /= 'ENDFU' ) 
             if ( words(1) == 'TIMES' ) then
                !
                ! Time shape
                !
                if( fundi(ifunc) % dim == 0 ) then             
                   call runend('KER_DISCRETE_FUNCTION_READ: DIMENSION OF THE DISCRETE FUNCTION '//trim(fundi(ifunc) % name)//' SHOULD BE GIVEN')
                end if
                if(      exists('SMOOT') ) then
                   fundi(ifunc) % kfl_type = DISCRETE_FUNCTION_SMOOTH 
                else if( exists('LINEA') ) then
                   fundi(ifunc) % kfl_type = DISCRETE_FUNCTION_LINEAR
                else if( exists('PIECE') ) then
                   fundi(ifunc) % kfl_type = DISCRETE_FUNCTION_PIECEWISE_CONSTANT 
                end if

                rauxi(:) = 0.0_rp

             else if ( words(1) == 'SHAPE' ) then   ! defining the shape by discrete points
                ! Reference time and value

                rauxi(1)= getrea('REFER', 0.0_rp, 'Reference value, sometimes useful')
                rauxi(2)= getrea('TSTAR', 0.0_rp, 'Startint time, sometimes useful')

                ! Save total number of tabular lines
                call ecoute('ker_discrete_function_read')
                fundi(ifunc) % num_tim = int(param(1)) ! number of data time

                ! Allocate the prescription time function vector for ifunc
                call ker_dicsrete_function_allocate(ifunc)

                ! Save discrete times and values
                call ecoute('ker_discrete_function_read')
                istep = 0
                do while(words(1)/='ENDSH')                   
                   istep = istep + 1
                   if(nnpar .ne. idime+1) call runend("KER_DISCRETE_FUNCTION_READ: FUNCTION "//trim(fundi(ifunc) % name)//" ROW "//trim(intost(istep))//&
                              " HAS "//trim(intost(nnpar-1))//" ELEMENTS (+1 MORE ELEMENT FOR THE TIME) INSTEAD OF DECLARED "//trim(intost(idime)))
                   if(istep > fundi(ifunc) % num_tim)  call runend("KER_DISCRETE_FUNCTION_READ: FUNCTION "//trim(fundi(ifunc) % name)//" HAS MORE ROWS ("&
                              //trim(intost(istep))//") THAN DECLARED ("//trim(intost(fundi(ifunc) % num_tim))//")")
                   fundi(ifunc) % tim(istep)         = param(1) - rauxi(2)         ! time
                   fundi(ifunc) % val(1:idime,istep) = param(2:idime+1) - rauxi(1) ! prescribed value
                   call ecoute('ker_discrete_function_read')
                end do
                if( istep .ne. fundi(ifunc) % num_tim) call runend("KER_DISCRETE_FUNCTION_READ: FUNCTION "//trim(fundi(ifunc) % name)//" HAS LESS ROWS ("&
                              //trim(intost(istep))//") THAN DECLARED ("//trim(intost(fundi(ifunc) % num_tim))//")")
             end if
             call ecoute('ker_discrete_function_read')

          end do

       end if

       call ecoute('ker_discrete_function_read')

    end do
    ! End functions

  end subroutine ker_discrete_function_read
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-05-29
  !> @brief   Initialization
  !> @details Module initialization
  !> 
  !-----------------------------------------------------------------------

  subroutine ker_discrete_function_initialization()

    num_discrete_functions = 0
    nullify(fundi)
    
  end subroutine ker_discrete_function_initialization
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-05-29
  !> @brief   Parallelization
  !> @details Parallelization
  !> 
  !-----------------------------------------------------------------------

  subroutine ker_dicsrete_function_allocate(ifunc)
 
    integer(ip), optional, intent(in) :: ifunc
    integer(ip)                       :: my_dim
    integer(ip)                       :: my_num_tim,ii

    if( num_discrete_functions > 0 ) then
       if( present(ifunc) ) then
          my_dim     = fundi(ifunc) % dim
          my_num_tim = fundi(ifunc) % num_tim
          call memory_alloca(mem_modul(1:2,modul),'FUNDI % DIM','ker_memall',fundi(ifunc) % tim,my_num_tim)
          call memory_alloca(mem_modul(1:2,modul),'FUNDI % VAL','ker_memall',fundi(ifunc) % val,my_dim,my_num_tim)
       else
          allocate(fundi(num_discrete_functions))
          do ii = 1,num_discrete_functions
             nullify(fundi(ii) % tim)
             nullify(fundi(ii) % val)
             fundi(ii) % dim          = 0
             fundi(ii) % num_tim      = 0
             fundi(ii) % kfl_type     = 0
             fundi(ii) % kfl_periodic = 0
             fundi(ii) % name         = ''
          end do
       end if
    end if
    
  end subroutine ker_dicsrete_function_allocate
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-05-29
  !> @brief   Parallelization
  !> @details Parallelization
  !> 
  !-----------------------------------------------------------------------
  
  subroutine ker_discrete_function_parall()

    integer(ip) :: ifunc

    call exchange_init()
    call exchange_add(num_discrete_functions) 
    call exchange_end()
    
    if( ISLAVE ) call ker_dicsrete_function_allocate()
    call exchange_init()
    do ifunc = 1,num_discrete_functions
       call exchange_add(fundi(ifunc) % kfl_type)
       call exchange_add(fundi(ifunc) % kfl_periodic)
       call exchange_add(fundi(ifunc) % dim)
       call exchange_add(fundi(ifunc) % num_tim)
    end do
    call exchange_end()
    if( ISLAVE ) then
       do ifunc = 1,num_discrete_functions
          call ker_dicsrete_function_allocate(ifunc)
       end do
    end if
    call exchange_init()
    do ifunc = 1,num_discrete_functions
       call exchange_add(fundi(ifunc) % val)
       call exchange_add(fundi(ifunc) % tim)
       call exchange_add(fundi(ifunc) % name)
    end do
    call exchange_end()

  end subroutine ker_discrete_function_parall
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-05-29
  !> @brief   Look for the time range
  !> @details Look for the current time range where
  !>          ctime \in [ tim(itime) : tim(itime+1) ]
  !> 
  !-----------------------------------------------------------------------

  subroutine ker_dicsrete_function_time_range(ctime,tim,kfl_periodic,itime,xfact) 

    real(rp),    intent(in)  :: ctime
    real(rp),    intent(in)  :: tim(:)
    integer(ip), intent(in)  :: kfl_periodic
    integer(ip), intent(out) :: itime
    real(rp),    intent(out) :: xfact

    integer(ip)              :: ntime,kk,istep
    real(rp)                 :: mintime,maxtime,cortime

    ntime    = size(tim)
    mintime  = tim(1)
    maxtime  = tim(ntime)
    !
    ! Periodic field, special treament
    !
    if( kfl_periodic == 1 ) then
       cortime = mod(ctime,maxtime)
    else
       cortime = ctime
    end if
    
    if(      cortime <= mintime ) then
       !
       ! Below minimum time
       !
       itime = 1
       xfact = 1.0_rp

    else if( cortime >= maxtime ) then
       !
       ! Above maximum time
       !
       itime = ntime-1 
       xfact = 0.0_rp       

    else
       !
       ! Right in the time range
       !
       kk = ntime
       steps: do istep = 1,ntime
          if( cortime <= tim(istep) ) then
             kk = istep - 1
             exit steps 
          end if
       end do steps
       kk = max(kk,1_ip)
       kk = min(kk,ntime-1)

       itime = kk              
       xfact = ( tim(kk+1) - cortime ) / ( tim(kk+1) - tim(kk) )

    end if

  end subroutine ker_dicsrete_function_time_range

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-05-29
  !> @brief   Discrete function
  !> @details Compute val from a discrete function
  !> 
  !-----------------------------------------------------------------------

  subroutine ker_discrete_function_scalar(ifunc,ctime,xx)
    integer(ip), intent(in)    :: ifunc
    real(rp),    intent(in)    :: ctime 
    real(rp),    intent(inout) :: xx
    real(rp)                   :: yy(1)

    yy(1) = xx
    call ker_discrete_function_go(ifunc,ctime,1_ip,yy)
    xx = yy(1)
    
  end subroutine ker_discrete_function_scalar

  subroutine ker_discrete_function_vector(ifunc,ctime,xx)
    integer(ip), intent(in)  :: ifunc
    real(rp),    intent(in)  :: ctime 
    real(rp),    intent(out) :: xx(:)
    integer(ip)              :: kdime

    kdime = size(xx)

    call ker_discrete_function_go(ifunc,ctime,kdime,xx)

  end subroutine ker_discrete_function_vector

  subroutine ker_discrete_function_go(ifunc,ctime,kdime,xx) 

    integer(ip), intent(in)  :: ifunc
    real(rp),    intent(in)  :: ctime
    integer(ip), intent(in)  :: kdime
    real(rp),    intent(out) :: xx(*)
    integer(ip)              :: itime    
    real(rp)                 :: xfact,t1,t2,m
    real(rp)                 :: x1(kdime),x2(kdime)

    call ker_dicsrete_function_time_range(ctime,fundi(ifunc) % tim,&
         fundi(ifunc) % kfl_periodic,itime,xfact)

    select case ( fundi(ifunc) % kfl_type )

    case ( DISCRETE_FUNCTION_PIECEWISE_CONSTANT )
       !
       ! Piecewise constant functions
       !
       if( xfact <= 0.0_rp ) then          
          xx(1:kdime) = fundi(ifunc) % val(1:kdime,itime+1)
       else
          xx(1:kdime) = fundi(ifunc) % val(1:kdime,itime)          
       end if
       
    case ( DISCRETE_FUNCTION_LINEAR )
       !
       ! Linear function
       !       
       xx(1:kdime) =  fundi(ifunc) % val(1:kdime,itime)   * xfact &
            &       + fundi(ifunc) % val(1:kdime,itime+1) * (1.0_rp-xfact)
       
    case ( DISCRETE_FUNCTION_SMOOTH )
       !
       ! Smooth function
       !
       t1           = fundi(ifunc) % tim(itime)
       t2           = fundi(ifunc) % tim(itime+1)
       m            = (ctime - t1)/(t2 - t1)
       x1(1:kdime)  = fundi(ifunc) % val(1:kdime,itime)
       x2(1:kdime)  = fundi(ifunc) % val(1:kdime,itime+1)
       xx(1:kdime)  = x1(1:kdime) + (x2(1:kdime) - x1(1:kdime))*(m**3.0_rp)*(10.0_rp - 15.0_rp*m + 6.0_rp*(m**2.0_rp))

    end select

  end subroutine ker_discrete_function_go

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    25/02/2013
  !> @brief   Get a function number
  !> @details Get a space and time function number given the function
  !>          name
  !
  !----------------------------------------------------------------------

  function ker_discrete_function_number(wfname)
    integer(ip)                :: ker_discrete_function_number
    character(*),  intent(in)  :: wfname
    integer(ip)                :: ifunc

    ker_discrete_function_number = 0
    do ifunc = 1,num_discrete_functions
       if( trim(wfname) == trim(fundi(ifunc) % name) ) then
          ker_discrete_function_number = ifunc
       end if
    end do
    if( ker_discrete_function_number == 0 ) &
         call runend('DISCRETE FUNCTION '//trim(wfname)//' DOES NOT EXIST')

  end function ker_discrete_function_number
  

  !----------------------------------------------------------------------
  !
  !> @author  Constantine Butakoff
  !> @date    26/10/2020
  !> @brief   Get a function's number of dimensions
  !> @details Get the dimension of a discrete  function given its number
  !
  !----------------------------------------------------------------------

  function ker_discrete_function_getdim(ifunc)
      integer(ip),  intent(in)   :: ifunc
      integer(ip)                :: ker_discrete_function_getdim

      ker_discrete_function_getdim = fundi(ifunc) % dim

  end function ker_discrete_function_getdim

end module mod_ker_discrete_function
!> @}
