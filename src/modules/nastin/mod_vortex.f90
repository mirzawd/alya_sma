!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-------------------------------------------------------------------------||---!
!-------------------------------------------------------------------| INIT |---!
!-------------------------------------------------------------------------||---!
module mod_vortex 
  use def_parame, only: ip, rp
  use def_domain, only: ndime, npoin
  use def_master, only: IMASTER, INOTMASTER, ISEQUEN, ID_NASTIN 
  use def_master, only: tempe, veloc, wmean, veloc  
  use def_nastin, only: kfl_regim_nsi
  use def_nastin, only: gamth_nsi   ! gamma = Cp/(Cp-R)
  use def_nastin, only: sphea_nsi   ! Specific heat Cp [J/K Kg]
  use def_nastin, only: prthe_nsi   ! Thermodynamic pressure = 101325.0_rp
  use def_nastin, only: bvess_nsi
  use def_nastin, only: nprev_nsi 
  use def_domain,   only: npoin, nboun, ndime, coord
  use mod_exchange, only :  exchange_add
  implicit none
  !-----------------------------------------------------------------------||---!
  !
  !-----------------------------------------------------------------------||---!
  real(rp)     :: parameters(6) = huge(1.0_rp) 
!  integer(ip)  :: initi
  integer(ip),  parameter :: VORTEX_ID   = 9   
  character(5), parameter :: VORTEX_NAME = 'VORTEX' 
  logical(ip),  parameter :: DEBUG       = .FALSE.
  !-----------------------------------------------------------------------||---!
  !
  !-----------------------------------------------------------------------||---!
  private
    public :: vortex_reabcs 
    public :: vortex_sendat  
    public :: vortex_parall  
    public :: vortex_iniunk  
    public :: VORTEX_NAME, VORTEX_ID  
!    public :: 
  !-----------------------------------------------------------------------||---!
  !
  !-----------------------------------------------------------------------||---!
  !=============================================================| contains |===!
contains
  !-----------------------------------------------------------------------||---! 
  subroutine vortex_reabcs( params, initi )
  implicit none
  real(rp),    intent(in   ) :: params(:)  
  integer(ip), intent(inout) :: initi  
  ! 
  ! vim nsi_reabcs.f90 +359 
  !         else if( words(1) == 'INITI' ) then
  !            ... 
  !            else if( exists('VORTE') ) then                                           ! 6: Poiseuille distribution 
  !               if (exists('BOUND')) then
  !  
  initi = VORTEX_ID  
  parameters(1:5) = params(3:)  !! Flow axis, Vin, alpha,   L,   n  
  !  
  if(DEBUG) print *, "[vortex_reabcs]", " axis,Vin,alpha,L,n:", parameters(:) 
  end subroutine 
  !-----------------------------------------------------------------------||---! 


  !-----------------------------------------------------------------------||---! 
  subroutine vortex_sendat()   
  implicit none
 !real(rp),    intent(in) :: poise_nsi(:) 
 !integer(ip), intent(in) :: kfl_initi_nsi   
  ! vim nsi_sendat.f90 +234
  !  
! call iexcha( initi )
  call rexcha( parameters(1) )
  call rexcha( parameters(2) )
  call rexcha( parameters(3) )
  call rexcha( parameters(4) )
  call rexcha( parameters(5) )
  ! 
  if(DEBUG) print *, "[vortex_sendat]", parameters(:) 
  !  
  end subroutine
  !-----------------------------------------------------------------------||---! 
  subroutine vortex_parall()   
  implicit none
 !real(rp),    intent(in) :: poise_nsi(:) 
 !integer(ip), intent(in) :: kfl_initi_nsi   
  ! vim nsi_sendat.f90 +234
  !  
  ! call iexcha( initi )
  call exchange_add(parameters)
  ! 
  if(DEBUG) print *, "[vortex_sendat]", parameters(:) 
  !  
  end subroutine
  !-----------------------------------------------------------------------||---! 


  !-----------------------------------------------------------------------||---! 
! vim nsi_iniunk.f90 +173 
  subroutine vortex_iniunk()
  implicit none
  integer(ip) :: axis 
  real(rp)    :: Vin, alpha,   L,   n    
  !
  ! L/2 = 20 * 7e-3  
  ! 1.0 + sin(2 * 3.14159 / (2 * 20 * 7e-3) * 1.0 * y) * 0.5 / ( 1 + exp( 1e6 * ( t - 2.0 * 1.773427e-05 ) ) )
  ! 
  ! axis, Vin, alpha,   L,   n  
  axis  = 1 !int(parameters(1)) 
  Vin   = parameters(2)  
  alpha = parameters(3)
  L     = parameters(4) 
  n     = parameters(5) 
  if(INOTMASTER) then 
      veloc(1:ndime,1:npoin,nprev_nsi) = 0.0 
      veloc(      1,1:npoin,nprev_nsi) = func( coord(2,1:npoin), L, n, alpha) * Vin 
  endif 
  !  
  if(DEBUG) print *, "[vortex_iniunk]"
  ! 
  end subroutine
  !-----------------------------------------------------------------------||---!
  !
  !-----------------------------------------------------------------------||---!


  !-----------------------------------------------------------------------||---!
  function func(x, D, n, alpha) result(f)
  implicit none
  real(rp), intent(in) :: x(npoin)
  real(rp), intent(in) :: D, n, alpha  ! input
  real(rp)             :: f(npoin)         ! output
  f = 1.0 + alpha * sin(2 * 3.14159 / D * n * x) 
  end function 
  !-----------------------------------------------------------------------||---!


  !=============================================================| contains |===!
end module 
!==============================================================================!
!==============================================================================!
