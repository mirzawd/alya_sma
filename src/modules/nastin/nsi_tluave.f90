!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_tluave(itask)

  !-----------------------------------------------------------------------
  !****f* Nastin/nsi_tluave
  ! NAME 
  !    nsi_tluave
  ! DESCRIPTION
  !    Calculate average velocity for two-layer model coupling
  !    U(n+1) = e*u(n+1) + (1-e)*U(n)
  !    where:
  !          e = dt/T is the weighting function
  !          T is the averaging period
  !          U(n) is the averaged velocity at timeG-step n
  !          u(n) is the instantaneous velocity at time-step n
  !
  !    ITASK = 1 ... Calculates time-averaged velocity at each time step
  !            2 ... Calculates initial value of time-averaged velocity
  ! USES
  ! USED BY
  !    nsi_averag
  !***
  !-----------------------------------------------------------------------

  use def_elmtyp
  use def_master
  use def_domain
  use def_nastin
  use def_kermod,         only : tlape_ker

  implicit none

  integer(ip),intent(in)    :: itask

  real(rp)                  :: avwei   ! e=dt/T, weighting function for averaging
  
  integer(ip)               :: ipoin


  if ( itask /= 2_ip ) then
     if ( abs(tlape_ker) > 1.0e-9_rp ) then
        avwei = dtime/tlape_ker    ! e=dt/T  where T is the time-period for averaging
     else
        avwei = 1.0_rp  ! If no averaging period is given in the input, the model uses the instantaneous velocity
     end if
  end if

  if ( itask == 1_ip ) then
  
     !-----------------------------------------------------------------------------
     !
     ! Calculate average velocity for the two-layer model
     !
     !-----------------------------------------------------------------------------
     
     do ipoin = 1, npoin
        tluav_nsi(1:ndime,ipoin) = avwei*veloc(1:ndime,ipoin,1)+(1.0_rp-avwei)*tluav_nsi(1:ndime,ipoin)
     end do

  else if ( itask == 2_ip ) then
     
     !-----------------------------------------------------------------------------
     !
     ! Calculate initial value of average velocity
     !
     !-----------------------------------------------------------------------------
 
     do ipoin = 1, npoin
        tluav_nsi(1:ndime,ipoin) = veloc(1:ndime,ipoin,1)
     end do

  end if

end subroutine nsi_tluave
