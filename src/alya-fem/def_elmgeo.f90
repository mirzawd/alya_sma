!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!>
!> @defgroup Elemental_Geometric_Toolbox
!> ToolBox for elemental and general geometrical operations
!> @{
!> @file    mod_elmgeo.f90
!> @author  Guillaume Houzeaux
!> @brief   ToolBox for elements
!> @details Different functions useful in finite element implementations
!>
!>          To add an element, edit the following subroutines, mainly
!>          for compatibility reasons:
!>
!>          def_elmgeo.f90
!>          mod_elmgeo.f90 
!>          elmtyp.f90
!>          domfac.f90
!>          bouele.f90
!>
!------------------------------------------------------------------------

module def_elmgeo


#include "def_vector_size.inc"
  use def_kintyp_basic,      only : ip,rp,lg,i1p,i2p
#ifndef I_AM_NOT_ALYA
  use def_master,            only : kfl_paral
#endif
  use def_elmtyp,            only : POINT
  use def_elmtyp,            only : BAR02,BAR03,BAR04,TRI03,TRI06,QUA04,QUA08,QUA09,QUA16, TRI10
  use def_elmtyp,            only : TET04,TET10,PYR05,PYR14,PEN06,PEN15,PEN18,HEX08, TET20
  use def_elmtyp,            only : HEX27,HEX64,SHELL,BAR3D,POI3D
  use def_elmtyp,            only : element_num_ini 
  use def_elmtyp,            only : element_num_end 
  use def_elmtyp,            only : element_end 
  use def_elmtyp,            only : element_max 
  implicit none
  private

  integer(ip) :: i
  
  !----------------------------------------------------------------------
  !
  ! List of edges and faces
  !
  !----------------------------------------------------------------------
  
  integer(ip), parameter :: list_edges_BAR02(2,1)  = reshape ( (/ 1,2                                                        /), (/2,1 /) )
  integer(ip), parameter :: list_faces_BAR02(1,2)  = reshape ( (/ 1,       2                                                 /), (/1,2 /) )
  integer(ip), parameter :: type_faces_BAR02(2)    =           (/ POINT,   POINT                                             /)
  integer(ip), parameter :: node_faces_BAR02(2)    =           (/ 1,       1                                                 /)

  integer(ip), parameter :: list_edges_BAR03(2,2)  = reshape ( (/ 1,3,     3,2                                               /), (/2,2 /) )
  integer(ip), parameter :: list_faces_BAR03(1,2)  = reshape ( (/ 1,       2                                                 /), (/1,2 /) )
  integer(ip), parameter :: type_faces_BAR03(2)    =           (/ POINT,   POINT                                             /)
  integer(ip), parameter :: node_faces_BAR03(2)    =           (/ 1,       1                                                 /)

  integer(ip), parameter :: list_edges_BAR04(2,3)  = reshape ( (/ 1,3,     3,4,     4,2                                      /), (/2,3 /) )
  integer(ip), parameter :: list_faces_BAR04(1,2)  = reshape ( (/ 1,       2                                                 /), (/1,2 /) )
  integer(ip), parameter :: type_faces_BAR04(2)    =           (/ POINT,   POINT                                             /)
  integer(ip), parameter :: node_faces_BAR04(2)    =           (/ 1,       1                                                 /)

  integer(ip), parameter :: list_edges_TRI03(2,3)  = reshape ( (/ 1,2,     2,3,     3,1                                      /), (/2,3 /) )
  integer(ip), parameter :: type_edges_TRI03(3)    =           (/ BAR02,   BAR02,   BAR02                                    /)
  integer(ip), parameter :: node_edges_TRI03(3)    =           (/ 2,       2,       2                                        /)
  integer(ip), parameter :: list_faces_TRI03(2,3)  = reshape ( (/ 1,2,     2,3,     3,1                                      /), (/2,3 /) )
  integer(ip), parameter :: type_faces_TRI03(3)    =           (/ BAR02,   BAR02,   BAR02                                    /)
  integer(ip), parameter :: node_faces_TRI03(3)    =           (/ 2,       2,       2                                        /)

  integer(ip), parameter :: list_edges_TRI06(3,3)  = reshape ( (/ 1,2,4,   2,3,5,   3,1,6                                    /), (/3,3 /) )
  integer(ip), parameter :: type_edges_TRI06(3)    =           (/ BAR03,   BAR03,   BAR03                                    /)
  integer(ip), parameter :: node_edges_TRI06(3)    =           (/ 3,       3,       3                                        /)
  integer(ip), parameter :: list_faces_TRI06(3,3)  = reshape ( (/ 1,2,4,   2,3,5,   3,1,6                                    /), (/3,3 /) )
  integer(ip), parameter :: type_faces_TRI06(3)    =           (/ BAR03,   BAR03,   BAR03                                    /)
  integer(ip), parameter :: node_faces_TRI06(3)    =           (/ 3,       3,       3                                        /)

  integer(ip), parameter :: list_edges_TRI10(4,3)  = reshape ( (/ 1,2,4,5, 2,3,6,7, 3,1,8,9                                  /), (/4,3 /) )
  integer(ip), parameter :: type_edges_TRI10(3)    =           (/ BAR04,   BAR04,   BAR04                                    /)
  integer(ip), parameter :: node_edges_TRI10(3)    =           (/ 4,       4,       4                                        /)
  integer(ip), parameter :: list_faces_TRI10(4,3)  = reshape ( (/ 1,2,4,5, 2,3,6,7, 3,1,8,9                                  /), (/4,3 /) )
  integer(ip), parameter :: type_faces_TRI10(3)    =           (/ BAR04,   BAR04,   BAR04                                    /)
  integer(ip), parameter :: node_faces_TRI10(3)    =           (/ 4,       4,       4                                        /)

  integer(ip), parameter :: list_edges_QUA04(2,4)  = reshape ( (/ 1,2,     2,3,     3,4,     4,1                             /), (/2,4 /) )
  integer(ip), parameter :: type_edges_QUA04(4)    =           (/ BAR02,   BAR02,   BAR02,   BAR02                           /)
  integer(ip), parameter :: node_edges_QUA04(4)    =           (/ 2,       2 ,      2,       2                               /)
  integer(ip), parameter :: list_faces_QUA04(2,4)  = reshape ( (/ 1,2,     2,3,     3,4,     4,1                             /), (/2,4 /) )
  integer(ip), parameter :: type_faces_QUA04(4)    =           (/ BAR02,   BAR02,   BAR02,   BAR02                           /)
  integer(ip), parameter :: node_faces_QUA04(4)    =           (/ 2,       2 ,      2,       2                               /)

  integer(ip), parameter :: list_edges_QUA08(3,4)  = reshape ( (/ 1,2,5,   2,3,6,   3,4,7,   4,1,8                           /), (/3,4 /) )
  integer(ip), parameter :: type_edges_QUA08(4)    =           (/ BAR03,   BAR03,   BAR03,   BAR03                           /)
  integer(ip), parameter :: list_faces_QUA08(3,4)  = reshape ( (/ 1,2,5,   2,3,6,   3,4,7,   4,1,8                           /), (/3,4 /) )
  integer(ip), parameter :: type_faces_QUA08(4)    =           (/ BAR03,   BAR03,   BAR03,   BAR03                           /)
  integer(ip), parameter :: node_faces_QUA08(4)    =           (/ 3,       3 ,      3,       3                               /)

  integer(ip), parameter :: list_edges_QUA09(3,4)  = reshape ( (/ 1,2,5,   2,3,6,   3,4,7,   4,1,8                           /), (/3,4 /) )
  integer(ip), parameter :: type_edges_QUA09(4)    =           (/ BAR03,   BAR03,   BAR03,   BAR03                           /)
  integer(ip), parameter :: node_edges_QUA09(4)    =           (/ 3,       3 ,      3,       3                               /)
  integer(ip), parameter :: list_faces_QUA09(3,4)  = reshape ( (/ 1,2,5,   2,3,6,   3,4,7,   4,1,8                           /), (/3,4 /) )
  integer(ip), parameter :: type_faces_QUA09(4)    =           (/ BAR03,   BAR03,   BAR03,   BAR03                           /)
  integer(ip), parameter :: node_faces_QUA09(4)    =           (/ 3,       3 ,      3,       3                               /)

  integer(ip), parameter :: list_edges_QUA16(4,4)  = reshape ( (/ 1,2,5,6, 2,3,7,8, 3,4,9,10, 4,1,11,12                       /), (/4,4 /) )
  integer(ip), parameter :: type_edges_QUA16(4)    =           (/ BAR04,   BAR04,   BAR04,    BAR04                           /)
  integer(ip), parameter :: list_faces_QUA16(4,4)  = reshape ( (/ 1,2,5,6, 2,3,7,8, 3,4,9,10, 4,1,11,12                       /), (/4,4 /) )
  integer(ip), parameter :: type_faces_QUA16(4)    =           (/ BAR04,   BAR04,   BAR04,    BAR04                           /)
  integer(ip), parameter :: node_faces_QUA16(4)    =           (/ 4,       4 ,      4,        4                               /)

  integer(ip), parameter :: list_edges_TET04(2,6)  = reshape ( (/ 1,2, 1,3, 1,4, 2,3, 4,2, 3,4                               /), (/2,6 /) )
  integer(ip), parameter :: type_edges_TET04(6)    =           (/ BAR02, BAR02, BAR02, BAR02, BAR02, BAR02                   /)
  integer(ip), parameter :: list_faces_TET04(3,4)  = reshape ( (/ 1,3,2,   2,3,4,   1,2,4,   3,1,4                           /), (/3,4 /) )
  integer(ip), parameter :: type_faces_TET04(4)    =           (/ TRI03,   TRI03,   TRI03,   TRI03                           /)
  integer(ip), parameter :: node_faces_TET04(4)    =           (/ 3,       3,       3,       3                               /)

  integer(ip), parameter :: list_edges_TET10(3,6)  = reshape ( (/ 1,2,5, 2,3,6, 3,1,7, 1,4,8, 2,4,9, 3,4,10                  /), (/3,6 /) )
  integer(ip), parameter :: type_edges_TET10(6)    =           (/ (BAR03,i=1,6)                                              /)
  integer(ip), parameter :: list_faces_TET10(6,4)  = reshape ( (/ 1,2,4,5,9,8, 2,3,4,6,10,9, 3,1,4,7,8,10, 1,3,2,7,6,5       /), (/6,4 /) )
  integer(ip), parameter :: type_faces_TET10(4)    =           (/ TRI06,   TRI06,   TRI06,   TRI06                           /)
  integer(ip), parameter :: node_faces_TET10(4)    =           (/ 6,       6,       6,       6                               /)

  integer(ip), parameter :: list_edges_TET20(2,6)   = reshape ( (/ 1,2, 1,3, 1,4, 2,3, 4,2, 3,4                               /), (/2,6 /) )
  integer(ip), parameter :: list_faces_TET20(10,4)  = reshape &
     ( (/ 1,3,2,10,9,8,7,6,5,17, 1,2,4,5,6,16,15,11,12,18, 3,1,4,9,10,12,11,13,14,19, 2,3,4,7,8,14,13,15,16,20                /), (/10,4 /) )
  integer(ip), parameter :: type_faces_TET20(4)     =           (/ TRI10,   TRI10,   TRI10,   TRI10                           /)
  integer(ip), parameter :: node_faces_TET20(4)     =           (/ 10,       10,       10,       10                           /)

  integer(ip), parameter :: list_edges_PYR05(2,8)  = reshape ( (/ 1,2, 2,3, 3,4, 4,1, 1,5, 2,5, 3,5, 4,5                     /), (/2,8 /) )
  integer(ip), parameter :: type_edges_PYR05(8)    =           (/ (BAR02,i=1,8)                                              /)
  integer(ip), parameter :: list_faces_PYR05(4,5)  = reshape ( (/ 1,4,3,2, 1,2,5,0, 2,3,5,0, 3,4,5,0, 4,1,5,0                /), (/4,5 /) )
  integer(ip), parameter :: type_faces_PYR05(5)    =           (/ QUA04,   TRI03,   TRI03,   TRI03,   TRI03                  /)
  integer(ip), parameter :: node_faces_PYR05(5)    =           (/ 4,       3,       3,       3,       3                      /)

  integer(ip), parameter :: list_edges_PEN06(2,9)  = reshape ( (/ 1,2, 2,3, 3,1, 4,5, 5,6, 6,4, 1,4, 3,6, 2,5                /), (/2,9 /) )
  integer(ip), parameter :: type_edges_PEN06(9)    =           (/ (BAR02,i=1,9)                                              /)
  integer(ip), parameter :: list_faces_PEN06(4,5)  = reshape ( (/ 1,3,2,0, 4,5,6,0, 1,2,5,4, 2,3,6,5, 3,1,4,6                /), (/4,5 /) )
  integer(ip), parameter :: type_faces_PEN06(5)    =           (/ TRI03,   TRI03, QUA04,   QUA04,   QUA04                    /)
  integer(ip), parameter :: node_faces_PEN06(5)    =           (/ 3,       3,     4,       4,       4                        /)

  ! PEN15 not yet fully integrated!
  integer(ip), parameter :: list_edges_PEN15(2,9)  = reshape ( (/ 1,2, 2,3, 3,1, 4,5, 5,6, 6,4, 1,4, 3,6, 2,5                /), (/2,9 /) )
  integer(ip), parameter :: list_faces_PEN15(8,5)  = &
    reshape ( (/ 1,3,2,9,8,7,0,0, 4,5,6,13,14,15,0,0, 1,2,5,4,7,11,13,10, 2,3,6,5,8,12,14,11, 3,1,4,6,9,10,15,12             /), (/8,5 /) )
  integer(ip), parameter :: type_faces_PEN15(5)    =           (/ TRI06,   TRI06, QUA08,   QUA08,   QUA08                    /)
  integer(ip), parameter :: node_faces_PEN15(5)    =           (/ 6,       6,     8,       8,       8                        /)

  integer(ip), parameter :: list_edges_PEN18(2,9)  = reshape ( (/ 1,2, 2,3, 3,1, 4,5, 5,6, 6,4, 1,4, 3,6, 2,5                /), (/2,9 /) )
  integer(ip), parameter :: list_faces_PEN18(9,5)  = &
    reshape ( (/ 1,3,2,9,8,7,0,0,0, 4,5,6,13,14,15,0,0,0, 1,2,5,4,7,11,13,10,16, 2,3,6,5,8,12,14,11,17, 3,1,4,6,9,10,15,12,18/), (/9,5 /) )
  integer(ip), parameter :: type_faces_PEN18(5)    =           (/ TRI06,   TRI06, QUA09,   QUA09,   QUA09                    /)
  integer(ip), parameter :: node_faces_PEN18(5)    =           (/ 6,       6,     9,       9,       9                        /)

  integer(ip), parameter :: list_edges_HEX08(2,12) = reshape ( (/ 1,2, 1,4, 1,5, 2,3, 2,6, 3,4, 3,7, 4,8, 5,6, 5,8, 6,7, 7,8 /), (/2,12/) )
  integer(ip), parameter :: type_edges_HEX08(12)   =           (/ (BAR02,i=1,12)                                             /)
  integer(ip), parameter :: list_faces_HEX08(4,6)  = reshape ( (/ 1,4,3,2, 2,3,7,6, 5,6,7,8, 4,1,5,8, 1,2,6,5, 3,4,8,7       /), (/4,6 /) )
  integer(ip), parameter :: type_faces_HEX08(6)    =           (/ QUA04,   QUA04,   QUA04,   QUA04,   QUA04,   QUA04         /)
  integer(ip), parameter :: node_faces_HEX08(6)    =           (/ 4,       4,       4,       4,       4,       4             /)

  integer(ip), parameter :: list_edges_HEX27(3,12) = reshape ( (/  1, 2, 9,  2, 3,10,  3, 4,11,  4, 1,12,                    &
       &                                                           5, 6,17,  6, 7,18,  7, 8,19,  8, 5,20,                    &
       &                                                           1, 5,13,  2, 6,14,  3, 7,15,  4, 8,16                     /), (/3,12/) )
  integer(ip), parameter :: type_edges_HEX27(12)   =           (/ (BAR03,i=1,12)                                             /)
  integer(ip), parameter :: list_faces_HEX27(9,6)  = &
       reshape ( (/ 1,2,6,5,9,14,17,13,22, 2,3,7,6,10,15,18,14,23, 3,4,8,7,11,16,19,15,24,  4,1,5,8,12,13,20,16,25, 1,4,3,2,12,11,10,9,21, 5,6,7,8,17,18,19,20,26 /), (/9,6/) )
  ! abel says there are drawing of the elements in some subroutine but I have not been able to find them.
  ! I did them by hand looking at subroutine shape3
  ! for the fouth face more consistent with hex08 would have been to start by 4 but in any case both options have teh same orientation.
  integer(ip), parameter :: type_faces_HEX27(6)    =           (/ QUA09,   QUA09,   QUA09,   QUA09,   QUA09,   QUA09         /)
  integer(ip), parameter :: node_faces_HEX27(6)    =           (/ 9,       9,       9,       9,       9,       9             /)
  
  integer(ip), parameter :: list_edges_HEX64(2,12) = reshape ( (/ 1,2, 1,4, 1,5, 2,3, 2,6, 3,4, 3,7, 4,8, 5,6, 5,8, 6,7, 7,8 /), (/2,12/) )
  integer(ip), parameter :: list_faces_HEX64(16,6)  = &
       reshape ( (/ 1, 4, 3, 2, 16 ,15, 14, 13, 12, 11, 10,  9, 33, 36, 35, 34,  &
       &            2, 3, 7, 6, 11 ,12, 19, 23, 28, 27, 22, 18, 39, 40, 48, 47,  &
       &            5, 6, 7, 8, 25 ,26, 27, 28, 29, 30, 31, 32, 53, 54, 55, 56,  &
       &            1, 5, 8, 4, 17 ,21, 32, 31, 24, 20, 15, 16, 44, 52, 51, 43,  &
       &            1, 2, 6, 5,  9 ,10, 18, 22, 26, 25, 21, 17, 37, 38, 46, 45,  &
       &            3, 4, 8, 7, 13 ,14, 20, 24, 30, 29, 23, 19, 41, 42, 50, 49   /), (/16,6/) )
  integer(ip), parameter :: type_faces_HEX64(6)    =           (/ QUA16,   QUA16,   QUA16,   QUA16,   QUA16,   QUA16 /)
  integer(ip), parameter :: node_faces_HEX64(6)    =           (/    16,      16,      16,      16,      16,      16 /)

  !----------------------------------------------------------------------
  !
  ! Element linearization
  !
  !----------------------------------------------------------------------
  
  integer(ip), target                   :: BAR02_TO_BAR02(2,1) = reshape ( (/ &
       &                                                       1, 2           &
       &                                                      /),(/2,1/) )
  integer(ip), target                   :: BAR03_TO_BAR02(2,2) = reshape ( (/ &
       &                                                       1, 3,          &
       &                                                       3, 2           &
       &                                                      /),(/2,2/) )
  integer(ip), target                   :: BAR03_TO_BAR03(3,1) = reshape ( (/ &
       &                                                       1, 2, 3        &
       &                                                      /),(/3,1/) )
  integer(ip), target                   :: BAR04_TO_BAR02(2,3) = reshape ( (/ &
       &                                                       1, 3,          &
       &                                                       3, 4,          &
       &                                                       4, 2           &
       &                                                      /),(/2,3/) )
  integer(ip), target                   :: BAR04_TO_BAR04(4,1) = reshape ( (/ &
       &                                                       1, 2, 3, 4     &
       &                                                      /),(/4,1/) )
  integer(ip), target                   :: TRI03_TO_TRI03(3,1) = reshape ( (/ &
       &                                                       1, 2, 3        &
       &                                                      /),(/3,1/) )
  integer(ip), target                   :: TRI06_TO_TRI03(3,4) = reshape ( (/ &
       &                                                       1, 4, 6,       &
       &                                                       4, 2, 5,       & 
       &                                                       4, 5, 6,       & 
       &                                                       6, 5, 3        & 
       &                                                      /),(/3,4/) )
  integer(ip), target                   :: QUA04_TO_QUA04(4,1) = reshape ( (/ &
       &                                                       1, 2, 3, 4     &
       &                                                      /),(/4,1/) )
  integer(ip), target                   :: QUA09_TO_QUA04(4,4) = reshape ( (/ &
       &                                                       1, 5, 9, 8,    &
       &                                                       5, 2, 6, 9,    &
       &                                                       8, 9, 7, 4,    &
       &                                                       9, 6, 3, 7     &
       &                                                      /),(/4,4/) )
  integer(ip), target                   :: QUA16_TO_QUA04(4,9) = reshape ( (/ &
       &                                                       1, 5,13,12,    &
       &                                                       5, 6,14,13,    &
       &                                                       6, 2, 7,14,    &
       &                                                      12,13,16,11,    &  
       &                                                      13,14,15,16,    &
       &                                                      14, 7, 8,15,    &
       &                                                      11,16,10, 4,    &
       &                                                      16,15, 9,10,    &
       &                                                      15, 8, 3, 9     &
       &                                                      /),(/4,9/) )
  integer(ip), target                   :: HEX08_TO_HEX08(8,1) = reshape ( (/         &
       &                                                       1, 2, 3, 4, 5, 6, 7, 8 &
       &                                                      /),(/8,1/) )
  integer(ip), target                   :: HEX27_TO_HEX08(8,8) = reshape ( (/         &
       &                                                       1, 9,21,12,13,22,27,25,&
       &                                                       9, 2,10,21,22,14,23,27,&
       &                                                      12,21,11, 4,25,27,24,16,&
       &                                                      21,10, 3,11,27,23,15,24,&
       &                                                      13,22,27,25, 5,17,26,20,&
       &                                                      22,14,23,27,17, 6,18,26,&
       &                                                      25,27,24,16,20,26,19, 8,&
       &                                                      27,23,15,24,26,18, 7,19 &
       &                                                      /),(/8,8/) )
  integer(ip), target                   :: HEX08_TO_TET04(4,6) = reshape ( (/ &
       &                                                       1, 3, 4, 5,    &
       &                                                       3, 4, 5, 8,    &
       &                                                       1, 2, 3, 5,    &
       &                                                       2, 3, 5, 6,    &
       &                                                       3, 5, 6, 7,    &
       &                                                       3, 5, 8, 7     &
       &                                                       /),(/4,6/) ) 
  integer(ip), target                   :: PEN06_TO_PEN06(6,1) = reshape ( (/   &
       &                                                       1, 2, 3, 4, 5, 6 &
       &                                                      /),(/6,1/) )
  integer(ip), target                   :: PEN06_TO_TET04(4,3) = reshape ( (/ &
       &                                                       1, 2, 3, 4,    &
       &                                                       2, 3, 4, 5,    &
       &                                                       3, 4, 5, 6     &
       &                                                       /),(/4,3/) ) 
  integer(ip), target                   :: PYR05_TO_TET04(4,2) = reshape ( (/ &
       &                                                       1, 2, 3, 5,    &
       &                                                       1, 3, 4, 5     &
       &                                                       /),(/4,2/) )   
  integer(ip), target                   :: PYR05_TO_PYR05(5,1) = reshape ( (/ &
       &                                                       1, 2, 3, 4, 5  &
       &                                                       /),(/5,1/) )   
  integer(ip), target                   :: TET04_TO_TET04(4,1) = reshape ( (/ &
       &                                                       1, 2, 3, 4     &
       &                                                       /),(/4,1/) )   
  integer(ip), target                   :: TET10_TO_TET04(4,8) = reshape ( (/ &
       &                                                       1, 5, 7, 8,    &
       &                                                       5, 2, 6, 9,    &
       &                                                       6, 3, 7,10,    &
       &                                                       5, 9, 6, 8,    &
       &                                                       9,10, 6, 8,    &
       &                                                       7, 6,10, 8,    &
       &                                                       7, 5, 6, 8,    &
       &                                                      10, 8, 9, 4     &
       &                                                       /),(/4,8/) )   

  !----------------------------------------------------------------------
  !
  ! Element type
  !
  ! toplogy:
  !                     0D          1D         2D          3D
  !                  --------------------------------------------
  !              -2  Point         
  !              -1              Lines          
  !               0                     Quadrilateral   Hexahedra
  !               1                          Triangle  Tetrahedra
  !               2                                    Pentahedra
  !               3                                       Pyramid
  !
  !----------------------------------------------------------------------

  integer(ip), parameter  :: nelty = element_max
  type elem_typ
     integer(ip)          :: type                   !< Type
     character(5)         :: name                   !< Name of the element
     character(20)        :: nametopo               !< Name of element topology
     integer(ip)          :: dimensions             !< Dimension of the element
     integer(ip)          :: order                  !< Order of the interpolation
     integer(ip)          :: topology               !< Topologoy of the element
     integer(ip)          :: number_nodes           !< Number of nodes
     logical(lg)          :: hessian                !< If Hessain exists
     integer(ip)          :: max_face_nodes         !< Max number of node per face
     integer(ip)          :: number_faces           !< Number of faces
     integer(ip)          :: number_edges           !< Number of edges (excluding diagonal and repeated edges)
     real(rp)             :: natural_length         !< Natural length (length of isoparametric element)
     real(rp)             :: cog(3)                 !< Center of gravity natural coordinates
     integer(ip), pointer :: type_edges(:)          !< Type of edges
     integer(ip), pointer :: list_edges(:,:)        !< List of edges
     integer(ip), pointer :: type_faces(:)          !< Face types
     integer(ip), pointer :: node_faces(:)          !< Face number of nodes
     integer(ip), pointer :: list_faces(:,:)        !< List of faces (normal pointing inward)
     real(rp),    pointer :: coord(:,:)             !< Isoparametric coordinates
   contains
     procedure,   pass    :: init
     procedure,   pass    :: nodes_to_type          !< Nodes to type
  end type elem_typ 
  type(elem_typ) :: element_type(nelty)
  type(elem_typ) :: element_type_init = elem_typ(&
       0_ip,&
       'NULL',&
       'NULL',&
       0_ip,&
       0_ip,&
       0_ip,&
       0_ip,&
       .false.,&
       0_ip,&
       0_ip,&
       0_ip,&
       1.0_rp,&
       (/ 0.0_rp,0.0_rp,0.0_rp /),&
       null(),&
       null(),&
       null(),&
       null(),&
       null(),&
       null())
  
  public :: list_faces_BAR02,type_faces_BAR02,node_faces_BAR02,list_edges_BAR02
  public :: list_faces_BAR03,type_faces_BAR03,node_faces_BAR03,list_edges_BAR03
  public :: list_faces_BAR04,type_faces_BAR04,node_faces_BAR04,list_edges_BAR04
  public :: list_faces_TRI03,type_faces_TRI03,node_faces_TRI03,list_edges_TRI03,type_edges_TRI03
  public :: list_faces_TRI06,type_faces_TRI06,node_faces_TRI06,list_edges_TRI06,type_edges_TRI06
  public :: list_faces_TRI10,type_faces_TRI10,node_faces_TRI10,list_edges_TRI10,type_edges_TRI10
  public :: list_faces_QUA04,type_faces_QUA04,node_faces_QUA04,list_edges_QUA04,type_edges_QUA04
  public :: list_faces_QUA08,type_faces_QUA08,node_faces_QUA08,list_edges_QUA08,type_edges_QUA08
  public :: list_faces_QUA09,type_faces_QUA09,node_faces_QUA09,list_edges_QUA09,type_edges_QUA09
  public :: list_faces_QUA16,type_faces_QUA16,node_faces_QUA16,list_edges_QUA16,type_edges_QUA16
  public :: list_faces_TET04,type_faces_TET04,node_faces_TET04,list_edges_TET04,type_edges_TET04
  public :: list_faces_TET10,type_faces_TET10,node_faces_TET10,list_edges_TET10,type_edges_TET10
  public :: list_faces_TET20,type_faces_TET20,node_faces_TET20,list_edges_TET20
  public :: list_faces_PYR05,type_faces_PYR05,node_faces_PYR05,list_edges_PYR05,type_edges_PYR05
  public :: list_faces_PEN06,type_faces_PEN06,node_faces_PEN06,list_edges_PEN06,type_edges_PEN06
  public :: list_faces_PEN15,type_faces_PEN15,node_faces_PEN15,list_edges_PEN15
  public :: list_faces_PEN18,type_faces_PEN18,node_faces_PEN18,list_edges_PEN18
  public :: list_faces_HEX08,type_faces_HEX08,node_faces_HEX08,list_edges_HEX08,type_edges_HEX08
  public :: list_faces_HEX27,type_faces_HEX27,node_faces_HEX27,list_edges_HEX27,type_edges_HEX27
  public :: list_faces_HEX64,type_faces_HEX64,node_faces_HEX64,list_edges_HEX64
  public :: element_type
  public :: nelty
  public :: elem_typ
  public :: element_type_init
  public :: BAR02_TO_BAR02
  public :: BAR03_TO_BAR02
  public :: BAR03_TO_BAR03 
  public :: BAR04_TO_BAR02
  public :: BAR04_TO_BAR04 
  public :: TRI03_TO_TRI03
  public :: TRI06_TO_TRI03
  public :: QUA04_TO_QUA04
  public :: QUA09_TO_QUA04
  public :: QUA16_TO_QUA04
  public :: HEX08_TO_HEX08
  public :: HEX27_TO_HEX08
  public :: HEX08_TO_TET04
  public :: PEN06_TO_TET04
  public :: PEN06_TO_PEN06
  public :: PYR05_TO_TET04
  public :: PYR05_TO_PYR05
  public :: TET04_TO_TET04
  public :: TET10_TO_TET04

contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-03-19
  !> @brief   Initialization
  !> @details Initialization
  !> 
  !-----------------------------------------------------------------------

  subroutine init(self)

    class(elem_typ), intent(inout) :: self

    !self = element_type_init

  end subroutine init

  !-----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    12/01/2018
  !> @brief   Guess element type
  !> @details Find the element type as a function of number of nodes
  !>          and dimension
  !>
  !-----------------------------------------------------------------------

  subroutine nodes_to_type(self)

    class(elem_typ), intent(inout) :: self
    integer(ip)                    :: pdime
    integer(ip)                    :: pnode

    pnode = self % number_nodes
    pdime = self % dimensions
    
    select case ( pdime )
       
    case ( 3_ip )
       
       select case ( pnode )
       case (  3 ) ; self % type = SHELL
       case (  4 ) ; self % type = TET04
       case (  5 ) ; self % type = PYR05
       case (  6 ) ; self % type = PEN06
       case (  8 ) ; self % type = HEX08
       case ( 10 ) ; self % type = TET10
       case ( 14 ) ; self % type = PYR14
       case ( 15 ) ; self % type = PEN15
       case ( 18 ) ; self % type = PEN18
       case ( 20 ) ; self % type = TET20
       case ( 27 ) ; self % type = HEX27
       case ( 64 ) ; self % type = HEX64
       end select

    case ( 2_ip )
       
       select case ( pnode )
       case (  3 ) ; self % type = TRI03
       case (  4 ) ; self % type = QUA04
       case (  6 ) ; self % type = TRI06
       case (  8 ) ; self % type = QUA08
       case (  9 ) ; self % type = QUA09
       case ( 10 ) ; self % type = TRI10
       case ( 16 ) ; self % type = QUA16
       end select

    case ( 1_ip )
       
       select case ( pnode )
       case ( 2 ) ; self % type = BAR02
       case ( 3 ) ; self % type = BAR03
       case ( 4 ) ; self % type = BAR04
       end select
       
    case ( 0_ip )
       
       self % type = POINT
       
    end select  

  end subroutine nodes_to_type

end module def_elmgeo
!> @}
