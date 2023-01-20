!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_outdt2(itask,kdime,geunk)
  !-----------------------------------------------------------------------
  !****f* Chemic/nsi_outdt2
  ! NAME 
  !    nsi_outdt2
  ! DESCRIPTION
  !    This routine computes time derivatives:
  !    ITASK = 1 ... du/dt
  !          = 2 ... d^2u/dt^2
  ! USED BY
  !    nsi_outvar
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_nastin
  implicit none
  integer(ip),   intent(in) :: itask,kdime
  real(rp),      intent(in) :: geunk(kdime,npoin,*)
  integer(ip)               :: ipoin,idime
  real(rp)                  :: dt1,dt2,xfac1,xfac2,xfac3,xfac4,xfac5
  real(rp)                  :: u1,u2,u3,dudt,d2udt2,dtht,dttr

  if( IMASTER ) return

  dt1   = dtold(1)
  dt2   = dtold(2)
  xfac1 = dt1*dt2*(dt1+dt2)
  xfac2 = dt2*(dt2+2.0_rp*dt1)
  xfac3 = (dt1+dt2)*(dt1+dt2)
  xfac4 = dt1*dt1
  xfac5 = dt1+dt2
  if( xfac1 == 0.0_rp ) xfac1 = 1.0e12_rp

  if( kdime == ndime ) then

     do ipoin = 1,npoin

        do idime = 1,kdime
           !
           ! Previous velocities
           !
           u1     =  geunk(idime,ipoin,3)
           u2     =  geunk(idime,ipoin,4)
           u3     =  geunk(idime,ipoin,5)
           dudt   = abs(( xfac2*u1 - xfac3*u2 + xfac4*u3 ) / xfac1)
           d2udt2 = abs(( dt2*u1   - xfac5*u2 + dt1*u3 ) / (0.5_rp*xfac1))
           if( abs(d2udt2) > 1.0e-6_rp .and. abs(dudt) > 1.0e-6_rp ) then
              dttr = 2.0_rp * epstr_nsi * abs(dudt) / abs(d2udt2)
           else
              dttr = strec_nsi * dtold(1)
           end if
           if( abs(d2udt2) > 1.0e-6_rp .and. abs(u1) > 1.0e-6_rp ) then
              dtht = sqrt(2.0_rp * epsht_nsi * abs(u1) / abs(d2udt2)) 
           else
              dtht = strec_nsi * dtold(1)
           end if
           
           if( itask == 1 ) then
              !
              ! First order time derivative: du/dt
              !
              gevec(idime,ipoin) = dudt

           else if( itask == 2 ) then
              !
              ! Second order time derivative: d^2u/dt^2
              !
              gevec(idime,ipoin) = d2udt2

           else if( itask == 3 ) then
              !
              ! First criteria
              !
              gesca(ipoin) = dtht

           else if( itask == 4 ) then
              !
              ! Second criteria
              !
              gesca(ipoin) = dttr

           end if

        end do

     end do

  else

     do ipoin = 1,npoin
        !
        ! Previous velocities
        !
        u1   =  geunk(1,ipoin,3)
        u2   =  geunk(1,ipoin,4)
        u3   =  geunk(1,ipoin,5)

        if( itask == 1 ) then
           !
           ! First order time derivative: du/dt
           !
           gesca(ipoin) = ( xfac2*u1 - xfac3*u2 + xfac4*u3 ) / xfac1

        else if( itask == 2 ) then
           !
           ! Second order time derivative: d^2u/dt^2
           !
           gesca(ipoin) = ( dt2*u1   - xfac5*u2 + dt1*u3 ) / (0.5_rp*xfac1)
        end if

     end do

  end if

end subroutine nsi_outdt2
