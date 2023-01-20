!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



module def_meteo
  use def_kintyp

  !----------------------------------------------------------------------
  !
  ! Meteo
  !
  !----------------------------------------------------------------------


  type meteo
     integer(ip)                    :: version, nx, ny, iproj
     real(rp)                       :: xfcst, xlvl, startlat, startlon, starti, startj, &
          deltalat, deltalon, dx, dy, xlonc, &
          truelat1, truelat2, earth_radius
     real(rp), pointer, dimension(:,:) :: slab
     logical                       :: is_wind_grid_rel
     character (len=8)             :: startloc
     character (len=9)             :: field
     character (len=24)            :: hdate
     character (len=25)            :: units
     character (len=32)            :: map_source
     character (len=46)            :: desc
  end type meteo

  ! Projection codes for proj_info structure:
  INTEGER, PUBLIC, PARAMETER  :: PROJ_LATLON = 0
  INTEGER, PUBLIC, PARAMETER  :: PROJ_LC = 1
  INTEGER, PUBLIC, PARAMETER  :: PROJ_PS = 2
  INTEGER, PUBLIC, PARAMETER  :: PROJ_PS_WGS84 = 102
  INTEGER, PUBLIC, PARAMETER  :: PROJ_MERC = 3
  INTEGER, PUBLIC, PARAMETER  :: PROJ_GAUSS = 4
  INTEGER, PUBLIC, PARAMETER  :: PROJ_CYL = 5
  INTEGER, PUBLIC, PARAMETER  :: PROJ_CASSINI = 6
  INTEGER, PUBLIC, PARAMETER  :: PROJ_ALBERS_NAD83 = 105
  INTEGER, PUBLIC, PARAMETER  :: PROJ_ROTLL = 203
  real(rp), parameter :: EARTH_RADIUS_M = 6370000.0_rp   ! same as MM5 system

  integer(8)                        ::     memor_meteo(2)
  !real(rp),pointer            :: uvel(:),vvel(:),temp(:),pres(:)




end module def_meteo
