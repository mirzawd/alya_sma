!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_usrbcs(iprob) 
  !------------------------------------------------------------------------
  !****f* Nastin/nsi_usrbcs
  ! NAME 
  !    nsi_usrbcs
  ! DESCRIPTION
  !    This routine sets user boundary contitions. There exists some
  !    predefined problems:
  !
  !    IPROB ... 0  Program yourself
  !          ... 1  Flow over cylinder (0,0) of radius 0.5
  !                 x=[-6,12], y=[-12,12]
  !          ... 2  Fully developed channel flow: [0,2]
  !          ... 3  NACA0015_6, chord=0.35
  ! USES
  ! USED BY
  !    nsi_reabcs
  !***
  !------------------------------------------------------------------------
  use def_parame
  use def_inpout
  use def_master
  use def_domain
  use def_nastin 
  use mod_memchk
  implicit none
  integer(ip), intent(in) :: iprob
  integer(ip)             :: ipoin,idime
  real(rp)                :: radio
!  real(rp)                :: ustar,nu,y,yplus

  select case(iprob)

  case(0_ip)
     !
     ! Program yourself
     !

  case(1_ip)
     !
     ! Flow over a cylinder
     !
     do ipoin = 1,npoin
        radio = 0.0_rp
        do idime = 1,ndime
           radio = radio + (coord(idime,ipoin)-0.0_rp)**2
        end do
        radio = sqrt(radio)
        if( radio <= 0.5_rp ) then                  ! Cylinder
           kfl_fixno_nsi(1,ipoin) = 1
           kfl_fixno_nsi(2,ipoin) = 1
           bvess_nsi(1,ipoin,1)   = 0.0_rp
           bvess_nsi(2,ipoin,1)   = 0.0_rp
        end if
        if(  coord(2,ipoin) >  11.9999_rp .or. &    ! Symmetry
             coord(2,ipoin) < -11.9999_rp ) then
           kfl_fixno_nsi(1,ipoin) = 0
           kfl_fixno_nsi(2,ipoin) = 1
           bvess_nsi(1,ipoin,1)   = 0.0_rp           
           bvess_nsi(2,ipoin,1)   = 0.0_rp           
        end if
        if( coord(1,ipoin) < -5.9999_rp ) then     ! Inflow
           kfl_fixno_nsi(1,ipoin) = 1
           kfl_fixno_nsi(2,ipoin) = 1
           bvess_nsi(1,ipoin,1)   = 1.0_rp           
           bvess_nsi(2,ipoin,1)   = 0.0_rp           
        end if        
     end do

  case(2_ip)
     !
     ! Fully developed channel flow: [0,2]
     !
     call runend('NSI_USRBCS: NOT CODED')
     !nu = visco_nsi(1,1)/densi_nsi(1,1)
     !do ipoin = 1,npoin
     !   if( coord(2,ipoin) < 1.0_rp ) then
     !      y = coord(2,ipoin)
     !   else
     !      y = 2.0_rp - coord(2,ipoin)
     !   end if
     !   ustar = 1.0_rp
     !   yplus = y*ustar/nu
     !   bvess_nsi(1,ipoin,1) = ustar*( log(1.0_rp+0.4_rp*yplus)/0.41_rp&
     !        +7.8_rp*(1.0_rp-exp(-yplus/11.0_rp)-yplus/11.0_rp*exp(-0.33_rp*yplus)))
     !end do

  case(3_ip)
     !
     ! NACA0015_6
     !
     do ipoin = 1,npoin
        if( lpoty(ipoin) /= 0 ) then

           if (   coord(1,ipoin) > -0.1_rp .and. coord(1,ipoin) < 0.5_rp .and. &
                & coord(2,ipoin) > -0.1_rp .and. coord(2,ipoin) < 0.1_rp ) then     
        
              kfl_fixno_nsi(1,ipoin) = 1
              kfl_fixno_nsi(2,ipoin) = 1
              bvess_nsi(1,ipoin,1)   = 0.0_rp           
              bvess_nsi(2,ipoin,1)   = 0.0_rp              

           else if(   coord(1,ipoin) > 0.296001_rp .and. &
                &   ( coord(2,ipoin) > 9.0_rp .or. coord(2,ipoin) < -9.0_rp ) ) then

              kfl_fixno_nsi(1,ipoin) = 0
              kfl_fixno_nsi(2,ipoin) = 1
              bvess_nsi(1,ipoin,1)   = 0.0_rp           
              bvess_nsi(2,ipoin,1)   = 0.0_rp

           else if( coord(1,ipoin) <= 0.296001_rp ) then

              kfl_fixno_nsi(1,ipoin) = 1
              kfl_fixno_nsi(2,ipoin) = 1
              bvess_nsi(1,ipoin,1)   = 1.0_rp           
              bvess_nsi(2,ipoin,1)   = 0.0_rp
              
           end if
        end if
     end do

  case(4_ip)
     !
     ! K1
     !
     do ipoin = 1,npoin
        if( lpoty(ipoin) /= 0 ) then

           if (   coord(1,ipoin) >  -3.69457_rp .and. coord(2,ipoin) < -3.7869_rp ) then     
        
              kfl_fixno_nsi(1,ipoin) = 0
              kfl_fixno_nsi(2,ipoin) = 1
              bvess_nsi(1,ipoin,1)   = 0.0_rp           
              bvess_nsi(2,ipoin,1)   = 0.0_rp              

           else if (   coord(1,ipoin) >  -3.69457_rp .and. coord(2,ipoin) > 3.2509_rp ) then     
        
              kfl_fixno_nsi(1,ipoin) = 0
              kfl_fixno_nsi(2,ipoin) = 1
              bvess_nsi(1,ipoin,1)   = 0.0_rp           
              bvess_nsi(2,ipoin,1)   = 0.0_rp              
              
           end if
        end if
     end do

  end select

end subroutine nsi_usrbcs
