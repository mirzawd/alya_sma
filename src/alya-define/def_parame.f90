!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



module def_parame

  !-----------------------------------------------------------------------
  !
  !     Parameters
  !
  !-----------------------------------------------------------------------  
  use def_kintyp, only : ip,rp
  !
  ! Frequently used mathematical constants (with precision to spare):
  !
  real(rp),     parameter :: pi      = 3.141592653589793_rp
                                      !3.141592653589793238462643383279502884197_rp
  real(rp),     parameter :: pio2    = 1.570796326794896_rp
                                      !1.57079632679489661923132169163975144209858_rp
  real(rp),     parameter :: twopi   = 6.283185307179586_rp
                                      !6.283185307179586476925286766559005768394_rp
  real(rp),     parameter :: in4pi   = 1.0_rp/(4.0_rp*pi)
  real(rp),     parameter :: sqrt2   = 1.414213562373095_rp
                                      !1.41421356237309504880168872420969807856967_rp
  real(rp),     parameter :: zero_rp = epsilon(0.0_rp)
  real(rp),     parameter :: kb      = 1.3806503e-23_rp
  real(rp),     parameter :: oneo3   = 1.0_rp/3.0_rp
  !
  ! Yes/No, Load/Unload and figures
  !
  integer(ip),  parameter :: yes=1_ip, no=0_ip
  integer(ip),  parameter :: mone=-1_ip, zero=0_ip, one=1_ip, two=2_ip, three=3_ip, four=4_ip
  integer(ip),  parameter :: five= 5_ip, six=6_ip, seven=7_ip, eight=8_ip, nine=9_ip, ten=10_ip
  integer(ip),  parameter :: eno=-1_ip
  integer(ip),  parameter :: load=100_ip, unload=200_ip
  !
  ! Null pointers
  !
  integer(ip),  target    :: nulip(1)
  real(rp),     target    :: nulir(1)
  !
  ! Special characters
  !
  character(1), parameter :: wback=achar(92)
  !
  ! Physical constants
  !
  real(rp),     parameter :: Stefan_Boltzmann = 5.670400e-08_rp ! W / (m^2 * K^4 )
  !
  ! Othres
  !
  real(rp),     parameter :: xmaxint4 = 2147483647.0_rp

end module def_parame
