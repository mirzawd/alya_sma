!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Output
!> @{
!> @file    mod_ansi_colors.f90
!> @author  houzeaux
!> @date    2020-05-03
!> @brief   Colors
!> @details Output with nice colors!
!-----------------------------------------------------------------------

module mod_ansi_colors

  use mod_strings, only : lower_case
  
  implicit none

  character(len=1), parameter :: c_esc           = achar(27)
  character(len=2), parameter :: c_start         = c_esc // '['
  character(len=1), parameter :: c_end           = 'm'
  character(len=*), parameter :: c_clear         = c_start // '0' // c_end
  
  character(len=*), parameter :: c_black         = '30'
  character(len=*), parameter :: c_red           = '31'
  character(len=*), parameter :: c_green         = '32'
  character(len=*), parameter :: c_yellow        = '33'
  character(len=*), parameter :: c_blue          = '34'
  character(len=*), parameter :: c_magenta       = '35'
  character(len=*), parameter :: c_cyan          = '36'
  character(len=*), parameter :: c_white         = '37'
  
  character(len=*), parameter :: c_grey          = '90'
  character(len=*), parameter :: c_light_red     = '91'
  character(len=*), parameter :: c_light_green   = '92'
  character(len=*), parameter :: c_light_yellow  = '93'
  character(len=*), parameter :: c_light_blue    = '94'
  character(len=*), parameter :: c_light_magenta = '95'
  character(len=*), parameter :: c_light_cyan    = '96'
  character(len=*), parameter :: c_light_white   = '97'

  private

  public :: ansi_colors
  public :: ansi_colors_name_to_code
  public :: ansi_colors_code_to_name

contains

  function ansi_colors_name_to_code(name_color) result(code_color)
    character(len=*), intent(in)  :: name_color
    character(len=:), allocatable :: code_color
    
    select case ( trim(lower_case(name_color)) )
       
    case( 'black'         ) ; code_color = c_black
    case( 'red'           ) ; code_color = c_red
    case( 'green'         ) ; code_color = c_green
    case( 'yellow'        ) ; code_color = c_yellow
    case( 'blue'          ) ; code_color = c_blue
    case( 'magenta'       ) ; code_color = c_magenta
    case( 'cyan'          ) ; code_color = c_cyan
    case( 'white'         ) ; code_color = c_white
       
    case( 'grey'          ) ; code_color = c_grey
    case( 'light red'     ) ; code_color = c_light_red
    case( 'light green'   ) ; code_color = c_light_green
    case( 'light yellow'  ) ; code_color = c_light_yellow
    case( 'light blue'    ) ; code_color = c_light_blue
    case( 'light magenta' ) ; code_color = c_light_magenta
    case( 'light cyan'    ) ; code_color = c_light_cyan
    case( 'light white'   ) ; code_color = c_light_white
       
    case( 'light_red'     ) ; code_color = c_light_red
    case( 'light_green'   ) ; code_color = c_light_green
    case( 'light_yellow'  ) ; code_color = c_light_yellow
    case( 'light_blue'    ) ; code_color = c_light_blue
    case( 'light_magenta' ) ; code_color = c_light_magenta
    case( 'light_cyan'    ) ; code_color = c_light_cyan
    case( 'light_white'   ) ; code_color = c_light_white
       
    case default      ;       code_color = ' ' ; return
       
    end select
    
  end function ansi_colors_name_to_code
  
  function ansi_colors_code_to_name(code_color) result(name_color)
    character(len=*), intent(in)  :: code_color
    character(len=:), allocatable :: name_color
    
    select case ( trim(code_color) ) 
       
    case( c_black         ) ; name_color = 'black'        
    case( c_red           ) ; name_color = 'red'          
    case( c_green         ) ; name_color = 'green'        
    case( c_yellow        ) ; name_color = 'yellow'       
    case( c_blue          ) ; name_color = 'blue'         
    case( c_magenta       ) ; name_color = 'magenta'      
    case( c_cyan          ) ; name_color = 'cyan'         
    case( c_white         ) ; name_color = 'white'        
                                                          
    case( c_grey          ) ; name_color = 'grey'         
    case( c_light_red     ) ; name_color = 'light red'    
    case( c_light_green   ) ; name_color = 'light green'  
    case( c_light_yellow  ) ; name_color = 'light yellow' 
    case( c_light_blue    ) ; name_color = 'light blue'   
    case( c_light_magenta ) ; name_color = 'light magenta'
    case( c_light_cyan    ) ; name_color = 'light cyan'   
    case( c_light_white   ) ; name_color = 'light white'  
                                                     
    case default      ;       name_color = ' ' ; return
       
    end select
    
  end function ansi_colors_code_to_name
  
  function ansi_colors(str, name_color, code_color) result(out)
    
    character(len=*), optional,   intent(in)  :: str
    character(len=*), optional,   intent(in)  :: name_color
    character(len=*), optional,   intent(in)  :: code_color
    character(len=:), allocatable             :: out
    character(2)                              :: code

    if( present(code_color) ) then

       code = code_color
       
    else if( present(name_color) ) then     

       code = ansi_colors_name_to_code(name_color) 
       
    end if

    if( trim(code) == '' ) then
       out = trim(str)
    else
       out = c_start // code // c_end // trim(str) // c_clear
    end if
    
  end function ansi_colors

end module mod_ansi_colors

!program ansi
!  use ansi_colors
!  implicit none
!  character(len=*), parameter :: endl = new_line('a')
!  print '(a)', &
!       color('Red',     c_red)     // endl // &
!       color('Green',   c_green)   // endl // &
!       color('Yellow',  c_yellow)  // endl // &
!       color('Blue',    c_blue)    // endl // &
!       color('Magenta', c_magenta) // endl // &
!       color('Cyan',    c_cyan)    // endl // &
!       color('White',   c_white)
!end program ansi
