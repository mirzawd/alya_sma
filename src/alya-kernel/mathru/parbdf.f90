!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine parbdf(order,pabdf)
!------------------------------------------------------------------------
!****f* mathru/parbdg
! NAME 
!    parbdf
! DESCRIPTION
!    This routine sets the coefficients for backward differentation
!    formula (BDF). The coefficients are:
!
!    k |      n   n-1   n-2    n-3   n-4   n-5  n-6
!    --+-------------------------------------------
!    1 |      1    -1
!    2 |    3/2    -2   1/2
!    3 |   11/6    -3   3/2   -1/3
!    4 |  25/12    -4     3   -4/3   1/4
!    5 | 137/60    -5     5  -10/3   5/4  -1/5
!    6 | 147/60    -6  15/2  -20/3  15/4  -6/5  1/6
!    --+-------------------------------------------
!
! USES
! USED BY
!    
!***
!------------------------------------------------------------------------
  use      def_kintyp
  use      def_master, only : dtime,dtold,kfl_timco
  implicit none
  integer(ip), intent(in)  :: order
  real(rp),    intent(out) :: pabdf(*)
  real(rp)                 :: ratio

  if( kfl_timco > 1 ) call runend('PARBDF: BDF MAKES NO SENSE WITH LOCAL TIME STEP')
  if( kfl_timco > 0 .and. order>2 ) call runend('PARBDF: BDF >2 ONLY WITH FIXED TIME STEP')

  if(order==1) then
     !
     ! 1/dt*( u^n - u^{n-1} )
     !
     pabdf(1)= 1.0_rp
     pabdf(2)=-1.0_rp

  else if(order==2) then
     ratio = dtime / dtold(1)
     if ( abs(ratio - 1.0_rp) < 1.0d-15) then 
        !
        ! 1/(2*dt)*( 3*u^n - 4*u^{n-1} + u^{n-2} )
        !
        pabdf(1)= 3.0_rp/2.0_rp
        pabdf(2)=-4.0_rp/2.0_rp
        pabdf(3)= 1.0_rp/2.0_rp
     else
        pabdf(1)= (1.0_rp + 2.0_rp*ratio) / (1.0_rp+ratio) 
        pabdf(2)=-(1.0_rp+ratio) 
        pabdf(3)= (ratio*ratio) / (1.0_rp+ratio)
     end if
     
  else if(order==3) then
     !
     ! 1/(6*dt)*( 11* u^n - 18*u^{n-1} + 9*u^{n-2} - 2*u^{n-3} )
     !
     pabdf(1)= 11.0_rp/6.0_rp
     pabdf(2)=-18.0_rp/6.0_rp
     pabdf(3)=  9.0_rp/6.0_rp
     pabdf(4)= -2.0_rp/6.0_rp
              
  else if(order==4) then
     !
     ! 1/(12*dt)*( 25*u^n - 48*u^{n-1} + 36*u^{n-2} - 16*u^{n-3} + 3*u^{n-4} )
     !
     pabdf(1)= 25.0_rp/12.0_rp
     pabdf(2)=-48.0_rp/12.0_rp
     pabdf(3)= 36.0_rp/12.0_rp
     pabdf(4)=-16.0_rp/12.0_rp
     pabdf(5)=  3.0_rp/12.0_rp

  else if(order==5) then
     !
     ! 1/(60*dt)*( 137*u^n - 300*u^{n-1} + 300*u^{n-2} - 200*u^{n-3} + 75*u^{n-4}
     !             -12*u^{n-5} )
     !
     pabdf(1)= 137.0_rp/60.0_rp
     pabdf(2)=-300.0_rp/60.0_rp
     pabdf(3)= 300.0_rp/60.0_rp
     pabdf(4)=-200.0_rp/60.0_rp
     pabdf(5)=  75.0_rp/60.0_rp
     pabdf(6)= -12.0_rp/60.0_rp

  else if(order==6) then
     !
     ! 1/(12*dt)*( 147*u^n - 360*u^{n-1} + 450*u^{n-2} - 400*u^{n-3} + 225*u^{n-4} 
     !             -72*u^{n-5} + 10*u^{n-6} ) 
     !
     pabdf(1)= 147.0_rp/60.0_rp
     pabdf(2)=-360.0_rp/60.0_rp
     pabdf(3)= 450.0_rp/60.0_rp
     pabdf(4)=-400.0_rp/60.0_rp
     pabdf(5)= 225.0_rp/60.0_rp
     pabdf(6)= -72.0_rp/60.0_rp
     pabdf(7)=  10.0_rp/60.0_rp
               
  end if

end subroutine parbdf
