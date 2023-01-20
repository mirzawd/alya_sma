!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tem_velfun(npoin,coord,vefun)
  !-----------------------------------------------------------------------
  !****f* Temper/tem_velfun
  ! NAME
  !   tem_velfun
  ! DESCRIPTION
  !   Compute velocity according to the function number 
  ! INPUT
  !   KFL_ADVEC_TEM ... Function
  !   COORD ........... Coordinates
  !   NPOIN ........... Number of nodes or Gauss points
  ! OUTPUT 
  !   VEFUN ........... Velocity
  ! USES
  ! USED BY
  !    tem_elmope
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp 
  use def_temper, only       :  kfl_advec_tem
  use def_master, only       :  veloc
  use def_domain, only       :  ndime
  implicit none
  integer(ip), intent(in)    :: npoin
  real(rp),    intent(in)    :: coord(ndime,npoin)
  real(rp),    intent(out)   :: vefun(ndime,npoin)
  integer(ip)                :: idime,ipoin
  real(rp)                   :: x,y,theta

  if(kfl_advec_tem==0) then
     do ipoin=1,npoin
        do idime=1,ndime
           vefun(idime,ipoin)= 0.0_rp
        end do
     end do
     
  else if(kfl_advec_tem==1) then
     do ipoin=1,npoin
        do idime=1,ndime
           vefun(idime,ipoin)= veloc(idime,ipoin,1)
        end do
     end do
     
  else if(kfl_advec_tem==2) then
     do ipoin=1,npoin
        x=coord(1,ipoin)
        y=coord(2,ipoin)
        vefun(1,ipoin)= 0.5_rp*(1.0_rp-x*x)*(1.0_rp+y)
        vefun(2,ipoin)=-0.5_rp*x*(4.0_rp-(1.0_rp+y)**2)
     end do

  else if(kfl_advec_tem==3) then
     do ipoin=1,npoin
        vefun(    1,ipoin) =  1.0_rp
        vefun(    2,ipoin) = -1.0_rp   
        vefun(ndime,ipoin) = -1.0_rp   
     end do

  else if(kfl_advec_tem==4) then
     do ipoin=1,npoin
        vefun(1,ipoin) = 1.0_rp
        if(ndime>1) vefun(2:ndime,ipoin) = 0.0_rp
     end do 

  else if(kfl_advec_tem==5) then
     theta = 1.1071487_rp  ! tan-1(2)
     do ipoin=1,npoin
        vefun(1,ipoin) =  cos(theta)
        vefun(2,ipoin) = -sin(theta)
     end do

  else if(kfl_advec_tem==6) then
     do ipoin=1,npoin
        vefun(1,ipoin) = 1.0_rp/sqrt(2.0_rp)
        vefun(2,ipoin) = 1.0_rp/sqrt(2.0_rp)
     end do

  else if(kfl_advec_tem==6) then
     do ipoin=1,npoin
        vefun(1,ipoin) = 1.0_rp/sqrt(2.0_rp)
        vefun(2,ipoin) = 1.0_rp/sqrt(2.0_rp)
     end do

  else if(kfl_advec_tem==7) then
     do ipoin=1,npoin
        vefun(1,ipoin)       = 1.0_rp
        vefun(2:ndime,ipoin) = 0.0_rp
     end do

  else if(kfl_advec_tem==8) then
     if(ndime == 2 ) then
        do ipoin=1,npoin
           vefun(1,ipoin) = 1.0_rp
           vefun(2,ipoin) = 1.0_rp
        end do
     else
        do ipoin=1,npoin
           vefun(1,ipoin) = 1.0_rp
           vefun(2,ipoin) = 0.0_rp
           vefun(3,ipoin) = 1.0_rp
        end do        
     end if

  else if(kfl_advec_tem==9) then
     do ipoin=1,npoin
        vefun(1,ipoin)     = 0.0_rp
        vefun(2,ipoin)     = 0.0_rp
        vefun(ndime,ipoin) = 1.0_rp
     end do

  else if(kfl_advec_tem==10) then

     do ipoin = 1,npoin
        vefun(1,ipoin) = -(coord(2,ipoin)-0.5_rp)
        vefun(2,ipoin) =  (coord(1,ipoin)-0.5_rp)
     end do

  else if(kfl_advec_tem==11) then

     do ipoin = 1,npoin
        vefun(1:ndime-1,ipoin) = 0.0_rp
        vefun(ndime,ipoin) = 1.0_rp
     end do

  end if

end subroutine tem_velfun
