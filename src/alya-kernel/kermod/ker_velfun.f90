!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Kermod
!> @{
!> @file    ker_velfun.f90
!> @author  Guillaume Houzeaux
!> @date    29/10/2014
!> @brief   Define advection velocity
!> @details Define the advection ADVEC:
!>
!>          ADVEC(1:NDIME,1:NPOIN,1) ... current advection
!>          ADVEC(1:NDIME,1:NPOIN,2) ... last coupling advection
!>          ADVEC(1:NDIME,1:NPOIN,3) ... last time step advection
!>
!>          According to KFL_VEFUN, it is computed as:
!>
!>          kfl_vefun = 0 ... ADVEC => VELOC. Nothing to do here.
!>                    < 0 ... ADVEC => FIELD. Constant in time and defined in ker_memall
!>                    > 0 ... ADVEC comes from a user defined function
!>
!>          Therefore, if kfl_vefun, there is nothing to do
!> @} 
!-----------------------------------------------------------------------
subroutine ker_velfun(itask)
  use def_parame
  use def_master
  use def_kermod
  use def_domain
  use mod_memory
  use mod_chktyp,                  only : check_type
  use mod_ker_space_time_function, only : ker_space_time_function

  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: ipoin,idime,ifunc,ifiel
  integer(ip)             :: kk,icomp
  real(rp)                :: AA,wx,wt,r,s,k,w,xx

  !----------------------------------------------------------------
  !
  ! ADVEC Computed only for user defined functions
  !
  !----------------------------------------------------------------

  if( INOTMASTER .and. kfl_vefun /= 0 ) then

     if( itask == ITASK_BEGSTE .or. itask == ITASK_INIUNK ) then
        
        if( kfl_vefun < 0 ) then
           !
           ! From Field
           !
           ifiel = -kfl_vefun
           if( kfl_field(4,ifiel) > 1 ) then
              kk = k_tran_fiel(ifiel)
              xx = x_tran_fiel(ifiel)

              do ipoin = 1,npoin
                 do idime = 1,ndime
                    advec(idime,ipoin,1) = xfiel(ifiel) % a(idime,ipoin,kk) * xx + &
                         xfiel(ifiel) % a(idime,ipoin,kk+1) * (1.0_rp-xx)
                 end do
              end do
           else
              do ipoin = 1,npoin
                 do idime = 1,ndime
                    advec(idime,ipoin,1) = xfiel(ifiel) % a(idime,ipoin,1) 
                 end do
              end do              
           end if
           
        else if( kfl_vefun > 1000 ) then
           !
           ! Space time function
           !
           do ipoin = 1,npoin
              ifunc = kfl_vefun  - 1000     
              call ker_space_time_function(&
                   ifunc,coord(1,ipoin),coord(2,ipoin),coord(ndime,ipoin),cutim,advec(1:ndime,ipoin,1))
           end do 

        else if( kfl_vefun > 0 ) then
           !
           ! Programmer functions
           !
           select case ( kfl_vefun )

           case ( 1_ip ) 

              wx = 2.0_rp
              wt = 2.0_rp
              AA = 0.0025_rp
              do ipoin = 1,npoin
                 !CALL RANDOM_NUMBER(dummr)
                 !advec(1,ipoin,1) = dummr*0.1_rp
                 !advec(1,ipoin,1) = 0.1_rp+(dummr-0.5_rp)*0.1_rp
                 advec(1,ipoin,1) = 0.1_rp
                 !advec(2,ipoin,1) = 0.03_rp*sin(2.0_rp*pi*cutim/4.0_rp)*sin(2.0_rp*pi*coord(1,ipoin)/(4.0_rp*0.1_rp))
                 advec(2,ipoin,1) = AA*sin(2.0_rp*pi*cutim/wt)*sin(2.0_rp*pi*coord(1,ipoin)/(wx*0.1_rp))
              end do

           case ( 2_ip )

              if( ndime == 3 ) then
                 do ipoin = 1,npoin
                    advec(1,ipoin,1) = -(coord(2,ipoin)-0.5_rp)
                    advec(2,ipoin,1) =  (coord(1,ipoin)-0.5_rp)
                    advec(3,ipoin,1) =  0.0_rp
                 end do
              else if( ndime == 2 ) then
                 do ipoin = 1,npoin
                    advec(1,ipoin,1) = -(coord(2,ipoin)-0.5_rp)
                    advec(2,ipoin,1) =  (coord(1,ipoin)-0.5_rp)
                 end do
              end if
              !do ipoin = 1,npoin
              !   advec(:,ipoin,3) = advec(:,ipoin,1)
              !end do
              
           case ( 3_ip ) 

              do ipoin = 1,npoin
                 advec(1:ndime,ipoin,1) = 0.0_rp
              end do

           case ( 4_ip ) 

              do ipoin = 1,npoin
                 advec(1,ipoin,1) = 0.1_rp
                 advec(2,ipoin,1) = 0.5_rp
              end do

           case ( 5_ip ) 

              do ipoin = 1,npoin
                 advec(    1,ipoin,1) = 1.0_rp
                 advec(    2,ipoin,1) = 0.0_rp
                 advec(ndime,ipoin,1) = 0.0_rp
              end do

           case ( 6_ip ) 

              do ipoin = 1,npoin
                 advec(1,ipoin,1) = 0.1_rp+0.1_rp*coord(1,ipoin)
                 if(coord(1,ipoin) > 0.1_rp) advec(1,ipoin,1) = 0.1_rp+0.1_rp*0.1_rp
              end do
              !-------- assign a helicoidal field-------------------   

           case ( 7_ip )

              r = 0.0_rp
              s = 5.0_rp
              k = 0.01_rp
              w = 0.01_rp
              do ipoin = 1,npoin
                 !--------- rotating vortex rope-------------
                 advec(1,ipoin,1) = -s * coord(2,ipoin) + r * s * sin((k * coord(3,ipoin)) + (w * cutim)) 
                 advec(2,ipoin,1) =  s * coord(1,ipoin) + r * s * cos((k * coord(3,ipoin)) + (w * cutim)) 
                 advec(3,ipoin,1) =  1.0_rp
                 !----------titling vortex ------------------
                 !advec(1,ipoin,1) = -coord(2,ipoin) + (cutim * coord(3,ipoin))
                 !advec(2,ipoin,1) =  coord(1,ipoin) - (cutim * coord(3,ipoin))
                 !advec(3,ipoin,1) =  coord(3,ipoin)
              end do

           case ( 8_ip ) 

              do ipoin = 1,npoin
                 advec(1:ndime,ipoin,1) = 0.0_rp
              end do

           case ( 9_ip ) 

              do ipoin = 1,npoin
                 advec(1,ipoin,1) = 0.01_rp * sin(1.0_rp*cutim)
              end do

           case ( 10_ip ) 
              !
              ! velocity field in a pipe (hadrien validation particle) 
              !
              !        
              !  U(r)=-R^2/4mu.(dP/dx).(1-r^2/R^2) => 2 * U * (1 - r^2/R^2)
              !  with U=>mean veloc, for a pipe U is Umax/2
              !       not the present case!!!  (pipe for R=3mm and U=0.5m/s (Re=200))
              !    present case : pipe for R=0.00225m and U=1.0m/s from
              !  (Numerical investigation of nano particles dispersion and deposition in fully developed laminar pipe flows)
              !

              !
              if(ndime==3)then
                 do ipoin = 1,npoin
                    advec(1,ipoin,1) = 2.0_rp*1.0_rp*(1-(coord(2,ipoin)*coord(2,ipoin)+coord(3,ipoin)*coord(3,ipoin))/(0.00225_rp*0.00225_rp))
                    advec(2,ipoin,1) = 0.0_rp
                    advec(3,ipoin,1) = 0.0_rp
                 end do
              else
                 do ipoin = 1,npoin
                    advec(1,ipoin,1) = -197530.86_rp*coord(2,ipoin)*coord(2,ipoin) + 888.88_rp * coord(2,ipoin) 
                    advec(2,ipoin,1) = 0.0_rp
                    advec(3,ipoin,1) = 0.0_rp
                 end do


              end if
           case ( 11_ip ) 

              do ipoin = 1,npoin
                 advec(1:ndime-1,ipoin,1) = 0.0_rp
                 advec(  ndime  ,ipoin,1) = 1.0_rp
              end do

           case ( 12_ip ) 

              do ipoin = 1,npoin
                 advec(1:ndime,ipoin,1) = real(ndime,rp) * coord(1:ndime,ipoin)
              end do

           case ( 13_ip ) 

              do ipoin = 1,npoin
                 r = sqrt( coord(1,ipoin)*coord(1,ipoin) + coord(2,ipoin)*coord(2,ipoin) )
                 ! Gamma = Re * 2* pi * nu, Gamma(Re=60,nu_air = 1.568e-5) = 5.911e-3
                 ! w_0 = Gamma / (2*pi*r) * ( 1 - exp( -r^2 / (4*nu*t) ) )
                 !
                 ! Pe = Gamma / (2*pi*nu) = 2
                 w = 0.0_rp
                 if (cutim > 0.0_rp) & 
                    w = 200.0_rp * 1.568e-5_rp / ( r + 0.0000001_rp) * ( 1.0_rp - exp( -r*r/ (4.0_rp * 1.568e-5_rp * cutim)) )
                    
                 advec(1,ipoin,1) =  -w * coord(2,ipoin) / r
                 advec(2,ipoin,1) =   w * coord(1,ipoin) / r
              end do

           case ( 14_ip ) 

              do ipoin = 1,npoin
                 r = sqrt( coord(1,ipoin)*coord(1,ipoin) + coord(2,ipoin)*coord(2,ipoin) )
                 ! Gamma = Re * 2* pi * nu, Gamma(Re=60,nu_air = 1.568e-5) = 5.911e-3
                 ! w_0 = Gamma / (2*pi*r) * ( 1 - exp( -r^2 / (4*nu*t) ) )
                 !
                 ! Pe = Gamma / (2*pi*nu) = 20
                 w = 0.0_rp
                 if (cutim > 0.0_rp) & 
                    w = 10.0_rp*200.0_rp * 1.568e-5_rp / ( r + 0.0000001_rp) * ( 1.0_rp - exp( -r*r/ (4.0_rp * 1.568e-5_rp * cutim)) )
                    
                 advec(1,ipoin,1) =  -w * coord(2,ipoin) / r
                 advec(2,ipoin,1) =   w * coord(1,ipoin) / r
              end do

           case ( 15_ip ) 

              do ipoin = 1,npoin
                 r = sqrt( coord(1,ipoin)*coord(1,ipoin) + coord(2,ipoin)*coord(2,ipoin) )
                 ! Gamma = Re * 2* pi * nu, Gamma(Re=60,nu_air = 1.568e-5) = 5.911e-3
                 ! w_0 = Gamma / (2*pi*r) * ( 1 - exp( -r^2 / (4*nu*t) ) )
                 !
                 ! Pe = Gamma / (2*pi*nu) = 200
                 w = 0.0_rp
                 if (cutim > 0.0_rp) & 
                    w = 100.0_rp*200.0_rp * 1.568e-5_rp / ( r + 0.0000001_rp) * ( 1.0_rp - exp( -r*r/ (4.0_rp * 1.568e-5_rp * cutim)) )
                    
                 advec(1,ipoin,1) =  -w * coord(2,ipoin) / r
                 advec(2,ipoin,1) =   w * coord(1,ipoin) / r
              end do

           case ( 16_ip ) 

              do ipoin = 1,npoin
                 r = sqrt( coord(1,ipoin)*coord(1,ipoin) + coord(2,ipoin)*coord(2,ipoin) )
                 ! Gamma = Re * 2* pi * nu, Gamma(Re=60,nu_air = 1.568e-5) = 5.911e-3
                 ! w_0 = Gamma / (2*pi*r) * ( 1 - exp( -r^2 / (4*nu*t) ) )
                 !
                 ! Pe = Gamma / (2*pi*nu) = 2000
                 w = 0.0_rp
                 if (cutim > 0.0_rp) & 
                    w = 1000.0_rp*200.0_rp * 1.568e-5_rp / ( r + 0.0000001_rp) * ( 1.0_rp - exp( -r*r/ (4.0_rp * 1.568e-5_rp * cutim)) )
                    
                 advec(1,ipoin,1) =  -w * coord(2,ipoin) / r
                 advec(2,ipoin,1) =   w * coord(1,ipoin) / r
              end do

           case ( 99_ip ) 
              !
              ! Do not do anything
              !
              continue

           end select

        end if
        !
        ! Assume constant initial advection
        !
        if( itask == ITASK_INIUNK ) then
           do icomp = 2,memory_size(advec,3_ip)
              do ipoin = 1,npoin 
                 do idime = 1,ndime
                    advec(idime,ipoin,icomp) = advec(idime,ipoin,1)
                 end do
              end do
           end do
        end if

     else if( itask == ITASK_ENDSTE ) then

        !----------------------------------------------------------------
        !
        ! Save previous advection
        ! KFL_VEFUN = 0, ADVEC point to VELOC which should not be modified
        !
        !----------------------------------------------------------------

        do ipoin = 1,npoin 
           do idime = 1,ndime
              advec(idime,ipoin,3) = advec(idime,ipoin,1)
           end do
        end do

     end if

  end if

end subroutine ker_velfun
