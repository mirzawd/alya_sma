!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_updtse(dtmin)
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
  implicit none
  real(rp),   intent(inout) :: dtmin
  integer(ip)               :: ipoin,idime,icrit
  real(rp)                  :: dtmax,dtht,dttr,dedt,d2edt2,facc
  real(rp)                  :: dt1,dt2,xfac1,xfac2,xfac3,xfac4,xfac5
  real(rp)                  :: dta,dtb,e1,e2,e3
  real(rp)                  :: dummr
  real(rp), save            :: energ_nsi(3)=0.0_rp

  icrit = 2

  if( INOTMASTER ) then

     if( icrit == 1 ) then
        energ_nsi(1) = 0.0_rp
        do ipoin = 1,npoin
           dummr = 0.0_rp
           do idime = 1,ndime
              dummr = dummr + veloc(idime,ipoin,3)*veloc(idime,ipoin,3)
           end do
           energ_nsi(1) = energ_nsi(1) + 0.5_rp*sqrt(dummr)
        end do
        energ_nsi(1) = energ_nsi(1) / real(npoin,rp)
     else
        !print*,lbsec(:)
        energ_nsi(1) = vbset(7,13) + vbset(7,14) + vbset(7,15) + vbset(7,16)
     end if

     if( ittim <= 4 ) then

        dtmin = dtmin 

     else

        e1    =  energ_nsi(1)
        e2    =  energ_nsi(2)
        e3    =  energ_nsi(3)

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
        !
        ! First order time derivative: du/dt
        !
        dedt   = ( xfac2*e1 - xfac3*e2 + xfac4*e3 ) / xfac1
        !
        ! Second order time derivative: d^2u/dt^2
        !
        d2edt2 = ( dt2*e1   - xfac5*e2 + dt1*e3 ) / (0.5_rp*xfac1)
        !
        ! Criterion of influence of higher order terms on u:
        ! DTHT <= 1/2 |d^e2/dt^2| * dt**2 = eps_TR * |u| 
        !
        if( abs(d2edt2) > 1.0e-6_rp .and. abs(e1) > 1.0e-6_rp ) then
           dtht = sqrt(2.0_rp * epsht_nsi * abs(e1) / abs(d2edt2))
        else
           !dtht = strec_nsi * dtold(1)
           dtht = 1.0e6_rp
        end if
        !
        ! Criterion of truncation error due to time discretization:
        ! DTTR <= 1/2 |d2^u/dt^2| * dt = eps_HT * |du/dt|
        !
        if( abs(d2edt2) > 1.0e-6_rp .and. abs(dedt) > 1.0e-6_rp ) then
           dttr = 2.0_rp * epstr_nsi * abs(dedt) / abs(d2edt2)
        else
           !dttr = strec_nsi * dtold(1)
           dttr = 1.0e6_rp
        end if
        !
        ! Maximum dt DTMAX over the nodes
        !
        dtmax = min(dtht,dttr)
        !
        ! Look for maximum over subdomains
        !
        !call pararr('MIN',0_ip,1_ip,dtmax)
        !
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
        dtmin = dtold(1) * min(strec_nsi,log(1.0_rp + (dampi_nsi-1.0_rp)*facc) / log(dampi_nsi))

        write(77,'(10(1x,e12.6))') facc,dt1,dttr,dtht,dedt,d2edt2
        flush(77)
        dtmin = max(dtmin,dta)
        dtmin = min(dtmin,dtb)
        !end if

        dtmin = dtmin / safet_nsi

     end if

  end if

  energ_nsi(3) = energ_nsi(2)
  energ_nsi(2) = energ_nsi(1)

end subroutine nsi_updtse
