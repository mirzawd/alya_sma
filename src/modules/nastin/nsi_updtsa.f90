!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_updtsa(dtmin)
  !-----------------------------------------------------------------------
  !****f* Chemic/nsi_updtsa
  ! NAME 
  !    nsi_updtsa
  ! DESCRIPTION
  !    This routine computes next timestep based on the accuracy reached
  !    with previous timestep
  ! USED BY
  !    nsi_timste
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_nastin
  use mod_iofile, only : iofile_flush_unit
  use mod_communications, only : PAR_MAX
  implicit none
  real(rp),   intent(inout) :: dtmin
  integer(ip)               :: ipoin,idime
  real(rp)                  :: dtmax,dtht,dttr,dudt,d2udt2,facc
  real(rp)                  :: dt1,dt2,xfac1,xfac2,xfac3,xfac4,xfac5
  real(rp)                  :: dta,dtb,dudtm,d2udt2m,u1,u2,u3
  real(rp)                  :: dthtm,dttrm

  if( ittim <= 4 ) then

     dtmin = dtmin

  else

     ! Adaptive timestep
     ! Prediction of next timestep dtn+1 based on accuracy
     ! reached with previous timestep dtold.

     dta   = dtmin
     dtb   = 100.0_rp * dtmin
     dtmax = 1.0e6_rp
     dt1   = dtold(1)
     dt2   = dtold(2)
     xfac1 = dt1*dt2*(dt1+dt2)
     xfac2 = dt2*(dt2+2.0_rp*dt1)
     xfac3 = (dt1+dt2)*(dt1+dt2)
     xfac4 = dt1*dt1
     xfac5 = dt1+dt2

     dudtm   = 0.0_rp
     d2udt2m = 0.0_rp
     dthtm   = 0.0_rp
     dttrm   = 0.0_rp

     if( INOTMASTER ) then

        do ipoin = 1,npoin

           do idime = 1,ndime
              !
              ! Previous velocities
              !
              u1   =  veloc(idime,ipoin,3)
              u2   =  veloc(idime,ipoin,4)
              u3   =  veloc(idime,ipoin,5)
              !
              ! First order time derivative: du/dt
              !
              dudt   = ( xfac2*u1 - xfac3*u2 + xfac4*u3 ) / xfac1
              !
              ! Second order time derivative: d^2u/dt^2
              !
              d2udt2 = ( dt2*u1   - xfac5*u2 + dt1*u3 ) / (0.5_rp*xfac1)
              !
              ! Criterion of influence of higher order terms on u:
              ! DTHT <= 1/2 |d^u2/dt^2| * dt**2 = eps_TR * |u| 
              !
              if( abs(d2udt2) > 1.0e-6_rp .and. abs(u1) > 1.0e-6_rp ) then
                 dtht = sqrt(2.0_rp * epsht_nsi * abs(u1) / abs(d2udt2))
              else
                 dtht = strec_nsi * dtold(1)
              end if
              !
              ! Criterion of truncation error due to time discretization:
              ! DTTR <= 1/2 |d2^u/dt^2| * dt = eps_HT * |du/dt|
              !
              if( abs(d2udt2) > 1.0e-6_rp .and. abs(dudt) > 1.0e-6_rp ) then
                 dttr = 2.0_rp * epstr_nsi * abs(dudt) / abs(d2udt2)
              else
                 dttr = strec_nsi * dtold(1)
              end if
              !
              ! Maximum dt DTMAX over the nodes
              !
              dtmax   = min(dtmax,dtht,dttr)
              dudtm   = max(dudtm,dudt)
              d2udt2m = max(d2udt2m,d2udt2)
              dthtm   = max(dthtm,dtht)
              dttrm   = max(dttrm,dttr)
           end do
        end do
     end if
     !
     ! Look for maximum over subdomains
     !
     call PAR_MAX(dtmax)

     !   - Calculation of accuracy reached with precious timestep dtold:
     !      facc=max(dt_trunc,dt_HT)/dtold

     facc = dtmax / dtold(1)

     !   - Prediction of next time step dtn+1 to be used in the solver:
     !   dtn+1=dtold * min(strech,ln(1+(damp-1)*facc)/ln(damp))
     !  
     !   strech: strech factor to limit timestep increase such thate dtnew <= P2.dtold.
     !   strech_m to be defined by the user in INPUT. strech >1. (eg: strech=1.05 means 5%).
     !   damp: damping factor to avoid oscillations of the timestep.
     !   Given by user in INPUT. damp > 1.0.

     !if( facc > 0.5_rp ) then
     !dtmin = min( dtmin, dtold(1) * min(strec_nsi,log(1.0_rp + (dampi_nsi-1.0_rp)*facc) / log(dampi_nsi))/safet_nsi)
     dtmin = dtold(1) * min(strec_nsi,log(1.0_rp + (dampi_nsi-1.0_rp)*facc) / log(dampi_nsi))
     write(77,'(10(1x,e12.6))') facc,dtmin,dta,dtb,dt1,dt2,dthtm,dttrm
     !flush(77)
     write(78,'(10(1x,e12.6))') dtmin,dudtm,d2udt2m
     !flush(78)
     dtmin = max(dtmin,dta)
     dtmin = min(dtmin,dtb)
     !end if

     dtmin = dtmin / safet_nsi

  end if

end subroutine nsi_updtsa
