!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine vocabu(itask,ndim1,ndim2)
  !------------------------------------------------------------------------
  !****f* kernel/vocabu
  ! NAME
  !    vocabu
  ! DESCRIPTION
  !    This routine sets Parall parameters according to the arrays
  !    to communicate (PARR1, PARR2, PARR3, PARI1, etc). ITASK is:
  !
  !    TYPE_KIND_DIM ..... TYPE = NELEM,NBOUN,NPOIN,NBOPO
  !                        KIND = INTE,REAL
  !                        DIM  = 1DIM,2DIM,3DIM,12DI
  !
  !    PARTY = 1,2,3,4 ... Vector dimensioned by 
  !    PARKI = 1,2,3 ..... Integer, real or character
  !    PARDI = 1,2,3 ..... Number of columns of a vector
  !    PARD1 = ........... Size of the first column
  !    PARD2 = ........... Size of the second column
  !          
  ! OUTPUT
  ! USED BY
  !    
  !***
  !------------------------------------------------------------------------
  use def_kintyp
  use def_master
  implicit none
  integer(ip), intent(in) :: itask,ndim1,ndim2

  select case(itask)

  case(-1_ip)
     !
     ! Translate kfl_paral to logical
     !
     if( kfl_paral == -1 ) then
        ISLAVE     = .false. 
        IMASTER    = .false.
        ISEQUEN    = .true.
        IPARALL    = .false.
        INOTMASTER = .true.
        INOTSLAVE  = .true.
        IEMPTY     = .false.
        INOTEMPTY  = .true.
     else if( kfl_paral == 0 ) then
        ISLAVE     = .false. 
        IMASTER    = .true.
        ISEQUEN    = .false.
        IPARALL    = .true.
        INOTMASTER = .false.
        INOTSLAVE  = .true.
        IEMPTY     = .true.
        INOTEMPTY  = .false.
     else if( kfl_paral >= 1 ) then
        ISLAVE     = .true. 
        IMASTER    = .false.
        ISEQUEN    = .false.
        IPARALL    = .true.
        INOTMASTER = .true.
        INOTSLAVE  = .false.
        IEMPTY     = .false.
        INOTEMPTY  = .true.
     end if

  case(NELEM_INTE_1DIM)
     !
     ! Element INT(NELEM)
     !
     party = 1
     parki = 1
     pardi = 1
     
  case(NELEM_INTE_2DIM)
     !
     ! Element INT(NDIM1,NELEM)
     !
     party = 1
     parki = 1
     pardi = 2
     pard1 = ndim1

  case(NBOUN_INTE_1DIM)
     !
     ! Boundary INT(NBOUN)
     !
     party =  2
     parki =  1
     pardi =  1

  case(NBOUN_INTE_2DIM)
     !
     ! Boundary INT(NDIM1,NBOUN)
     !
     party =  2
     parki =  1
     pardi =  2
     pard1 =  ndim1

  case(NBOUN_REAL_2DIM)
     !
     ! Boundary REAL(NDIM1,NBOUN)
     !
     party =  2
     parki =  2
     pardi =  2
     pard1 =  ndim1

  case(NBOUN_REAL_3DIM)
     !
     ! Boundary REAL(NDIM1,NDIM2,NBOUN)
     !
     party =  2
     parki =  2
     pardi =  3
     pard1 =  ndim1
     pard2 =  ndim2

  case(NPOIN_INTE_1DIM)
     !
     ! Node INT(NPOIN)
     !
     party =  3
     parki =  1
     pardi =  1
     
  case(NPOIN_INTE_2DIM)
     !
     ! Node INT(NDIM1,NPOIN)
     !
     party =  3
     parki =  1
     pardi =  2
     pard1 =  ndim1

  case(NPOIN_REAL_1DIM)
     !
     ! Node REAL(NPOIN)
     !
     party =  3 
     parki =  2
     pardi =  1 

  case(NPOIN_REAL_2DIM)
     !
     ! Node REAL(NDIM1,NPOIN)
     !
     party =  3
     parki =  2
     pardi =  2
     pard1 =  ndim1

  case(NPOIN_REAL_12DI)
     !
     ! Node REAL(NDIM1,NPOIN) => REAL(NDIM1*NPOIN)
     !
     party =  3 
     parki =  5  
     pardi =  1 
     pard1 =  ndim1
     
  case(NBOPO_INTE_1DIM)
     !
     ! Boundary node INT(NBOPO)
     !
     party =  4
     parki =  1
     pardi =  1
     
  case(NELEM_REAL_2DIM)
     !
     ! Element REAL(NDIM1,NPOIN)
     !
     party =  1
     parki =  2
     pardi =  2
     pard1 =  ndim1

  case(NELEM_REAL_3DIM)
     !
     ! Element REAL(NDIM1,NDIM2,NELEM)
     !
     party = 1
     parki = 2
     pardi = 3
     pard1 = ndim1
     pard2 = ndim2
     
  case(NPOIN_REAL_3DIM)
     !
     ! Node REAL(NDIM1,NDIM2,NPOIN)
     !
     party =  3
     parki =  2
     pardi =  3
     pard1 =  ndim1
     pard2 =  ndim2
     
  case(NBOPO_REAL_2DIM)
     !
     ! Boundary node REAL(NDIME,NBOPO)
     !
     party =  4
     parki =  2
     pardi =  2
     pard1 = ndim1

  case(NBOPO_REAL_3DIM)
     !
     ! Boundary node REAL(NDIME,NDIME,NBOPO)
     !
     party =  4
     parki =  2
     pardi =  3
     pard1 = ndim1
     pard2 = ndim2

  end select

end subroutine vocabu
