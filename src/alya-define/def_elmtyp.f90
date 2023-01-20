!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Domain
!> @{
!> @file    def_elmtyp.f90
!> @author  houzeaux
!> @date    2020-05-08
!> @brief   Element type numbering
!> @details Element type naming and numbering
!>          Elements available in Alya
!>          The original numbering followed CGNS standrad:
!>          http://www.grc.nasa.gov/WWW/cgns/sids/conv.html
!>          As new elements have been introduced, CGNS standard is no 
!>          longer used.
!>          The numbering is:
!>          1D elements:  2 ->  9
!>          2D elements: 10 -> 29
!>          3D elements: 30 -> 50
!>         
!-----------------------------------------------------------------------

module def_elmtyp
  
  use def_kintyp_basic, only   :  ip
  !
  ! Starting according to dimension 0D, 1D, 2D and 3D
  !
  integer(ip), parameter :: element_num_ini(0:3) = (/ 1,2,10,30 /)
  integer(ip), parameter :: element_num_end(0:3) = (/ 1,9,29,99 /)
  integer(ip), parameter :: element_end          = element_num_end(3)
  integer(ip), parameter :: element_max          = 110
  !
  ! Different elements
  !
  integer(ip), parameter :: POINT =    1 ! 0D

  integer(ip), parameter :: BAR02 =    2 ! 1D
  integer(ip), parameter :: BAR03 =    3 ! 1D 
  integer(ip), parameter :: BAR04 =    4 ! 1D 

  integer(ip), parameter :: TRI03 =   10 ! 2D 
  integer(ip), parameter :: TRI06 =   11 ! 2D 
  integer(ip), parameter :: QUA04 =   12 ! 2D 
  integer(ip), parameter :: QUA08 =   13 ! 2D 
  integer(ip), parameter :: QUA09 =   14 ! 2D 
  integer(ip), parameter :: QUA16 =   15 ! 2D 
  integer(ip), parameter :: TRI10 =   16 ! 2D NEW TRI10

  integer(ip), parameter :: TET04 =   30 ! 3D 
  integer(ip), parameter :: TET10 =   31 ! 3D 
  integer(ip), parameter :: PYR05 =   32 ! 3D 
  integer(ip), parameter :: PYR14 =   33 ! 3D 
  integer(ip), parameter :: PEN06 =   34 ! 3D 
  integer(ip), parameter :: PEN15 =   35 ! 3D 
  integer(ip), parameter :: PEN18 =   36 ! 3D 
  integer(ip), parameter :: HEX08 =   37 ! 3D 
  !integer(ip), parameter :: HEX20 =   38 ! 3D 
  integer(ip), parameter :: HEX27 =   39 ! 3D 
  integer(ip), parameter :: HEX64 =   40 ! 3D 
  integer(ip), parameter :: TET20 =   41 ! 3D NEW TET20
  integer(ip), parameter :: SHELL =   51 ! 3D shell element
  integer(ip), parameter :: BAR3D =   52 ! 3D bar element
  integer(ip), parameter :: POI3D =   53 ! 3D point element

  integer(ip), parameter :: DDDNE = -100 ! DD: D/N element
  integer(ip), parameter :: DDESS = -101 ! DD 
  integer(ip), parameter :: DDNAT = -102 ! DD
  integer(ip), parameter :: DDROB = -103 ! DD
  integer(ip), parameter :: DDEXT = -104 ! DD
  !
  ! Element characteristics
  !
  integer(ip), parameter :: ELFEM =    0 ! Finite element
  integer(ip), parameter :: ELEXT =    1 ! Extension finite element
  integer(ip), parameter :: ELHOL =    2 ! Hole element
  integer(ip), parameter :: ELCNT =    3 ! Contact element
  integer(ip), parameter :: ELCUT =    4 ! Cut element
  integer(ip), parameter :: ELHAN =    5 ! Element with hanging node
  integer(ip), parameter :: ELCOH =    6 ! Cohesive element (Camacho & Ortiz)
  integer(ip), parameter :: ELINT =    7 ! Interface cohesive elements () 
  integer(ip), parameter :: ELVIR =    8 ! Virtual element
  !
  ! Boundary characteristics
  !
  integer(ip), parameter :: BOFEM =    0 ! Finite boundary
  integer(ip), parameter :: BOEXT =    1 ! Extension boundary
  integer(ip), parameter :: BOHOL =    2 ! Hole boundary
  integer(ip), parameter :: BOINT =    3 ! Internal boundary
  integer(ip), parameter :: BOFRI =    4 ! Fringe boundary
  !
  ! Node characteristics
  !
  integer(ip), parameter :: NOFEM =    0 ! Normal node
  integer(ip), parameter :: NOEXT =    1 ! Extension  node
  integer(ip), parameter :: NOHOL =    2 ! Hole node
  integer(ip), parameter :: NOFRI =    4 ! Fringe node
  integer(ip), parameter :: NODE_CONTACT_FLUID =    3 ! Contact node on the fluid side
  integer(ip), parameter :: NODE_CONTACT_SOLID =   -3 ! Contact node on the solid side

end module def_elmtyp

