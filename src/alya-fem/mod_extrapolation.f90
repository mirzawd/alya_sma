!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Kinds_and_types
!> @{
!> @file    mod_extrapolation.f90
!> @author  houzeaux
!> @date    2020-04-04
!> @brief   Extrapolation arrays
!> @details Extrapolation shape function, derivative and Hessian
!-----------------------------------------------------------------------

module mod_extrapolation

  use def_kintyp_basic, only : ip,rp

  private

  public :: shaga0
  public :: shaga1
  public :: shaga2
  public :: shaga3

contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-02-13
  !> @brief   Extrapolation 0D
  !> @details Extrapolation 0D
  !> 
  !-----------------------------------------------------------------------

  pure subroutine shaga0(ngaus,shaga,ierro)

    integer(ip),           intent(in)  :: ngaus
    real(rp),              intent(out) :: shaga(ngaus)
    integer(ip), optional, intent(out) :: ierro

    if( present(ierro) ) ierro    = 0
    shaga(1) = 1.0_rp

  end subroutine shaga0

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-02-13
  !> @brief   Extrapolation 1D
  !> @details This routine evaluates shape functions associated to gauss points
  !>          for 1D : NGAUS = 1, 2, 3
  !> 
  !-----------------------------------------------------------------------

  pure subroutine shaga1(s,ngaus,shaga,ierro)

    real(rp),              intent(in)  :: s
    integer(ip),           intent(in)  :: ngaus
    real(rp),              intent(out) :: shaga(ngaus)
    integer(ip), optional, intent(out) :: ierro
    integer(ip)                        :: ierr
    
    ierr=0
    if(ngaus==1) then
       shaga(1)=1.0_rp
    else if(ngaus==2) then
       shaga(1)= 0.5_rp*sqrt(3.0_rp)*(1.0_rp/sqrt(3.0_rp)-s)
       shaga(2)= 0.5_rp*sqrt(3.0_rp)*(1.0_rp/sqrt(3.0_rp)+s)
    else if(ngaus==3) then
       shaga(1)= 5.0_rp/6.0_rp*(s-sqrt(0.6_rp))*s
       shaga(2)=-5.0_rp/3.0_rp*(s*s-0.60_rp)
       shaga(3)= 5.0_rp/6.0_rp*(s+sqrt(0.6_rp))*s
    else
       ierr=1
    end if

    if( present(ierro) ) ierro = ierr
    
  end subroutine shaga1

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-02-13
  !> @brief   Extrapolation 2D
  !> @details This routine evaluates shape functions associated to Gauss points
  !>          for 2D. The cases available so far are:
  !>
  !>          BRICKS    -->   NGAUS =   1   4   9  
  !>          SIMPLICES -->   NGAUS =   1   3   4  6  7  
  !> 
  !-----------------------------------------------------------------------

  pure subroutine shaga2(x,ltopo,ngaus,shaga,ierro)

    real(rp),              intent(in)  :: x(2)
    integer(ip),           intent(in)  :: ngaus
    integer(ip),           intent(in)  :: ltopo
    real(rp),              intent(out) :: shaga(ngaus)
    integer(ip), optional, intent(out) :: ierro
    integer(ip)                        :: ierr
    real(rp)                           :: s,t

    s = x(1)
    t = x(2)
    ierr = 0
    !
    ! Quadrilateral and hexahedral elements
    !
    if(ngaus==1) then
       shaga(1)=1.0_rp
    else if(ltopo==0 .and. ngaus==4) then
       shaga(1)= .75_rp*(s-1.0_rp/sqrt(3.0_rp))*(t-1.0_rp/sqrt(3.0_rp))
       shaga(2)=-.75_rp*(s-1.0_rp/sqrt(3.0_rp))*(t+1.0_rp/sqrt(3.0_rp))
       shaga(3)=-.75_rp*(s+1.0_rp/sqrt(3.0_rp))*(t-1.0_rp/sqrt(3.0_rp))
       shaga(4)= .75_rp*(s+1.0_rp/sqrt(3.0_rp))*(t+1.0_rp/sqrt(3.0_rp))
    else if(ltopo==0 .and. ngaus==9) then
       shaga(1)= 25.0_rp/36.0_rp*(s-sqrt(.6_rp))*(t-sqrt(.6_rp))*s*t
       shaga(2)=-25.0_rp/18.0_rp*(s-sqrt(.6_rp))*(t*t-.6_rp)*s
       shaga(3)= 25.0_rp/36.0_rp*(s-sqrt(.6_rp))*(t+sqrt(.6_rp))*s*t
       shaga(4)=-25.0_rp/18.0_rp*(s*s-.6_rp)*(t-sqrt(.6_rp))*t
       shaga(5)= 25.0_rp/9.0_rp*(s*s-.6_rp)*(t*t-.6_rp)
       shaga(6)=-25.0_rp/18.0_rp*(s*s-.6_rp)*(t+sqrt(.6_rp))*t
       shaga(7)= 25.0_rp/36.0_rp*(s+sqrt(.6_rp))*(t-sqrt(.6_rp))*s*t
       shaga(8)=-25.0_rp/18.0_rp*(s+sqrt(.6_rp))*(t*t-.6_rp)*s
       shaga(9)= 25.0_rp/36.0_rp*(s+sqrt(.6_rp))*(t+sqrt(.6_rp))*s*t
       !
       ! Triangular and tetrahedral elements
       !        
    else if (ltopo==1 .and. ngaus==3) then
       shaga(1)=2.0_rp*s-1.0_rp/3.0_rp
       shaga(2)=2.0_rp*t-1.0_rp/3.0_rp
       shaga(3)=2.0_rp*(1.0_rp-s-t)-1.0_rp/3.0_rp
    else if (ltopo==1 .and. ngaus==4) then
       shaga(1)=(-45.0_rp*(s+t)+225.0_rp*s*t+9.0_rp)/4.0_rp
       shaga(2)=(  5.0_rp*(s+t)- 75.0_rp*s*t+5.0_rp)/4.0_rp
       shaga(3)=(25.0_rp*s+15.0_rp*t-75.0_rp*s*t-5.0_rp)/4.0_rp
       shaga(4)=(15.0_rp*s+25.0_rp*t-75.0_rp*s*t-5.0_rp)/4.0_rp
    else if (ltopo==1 .and. ngaus==6) then
       shaga( 1)= 0.13855958741e+00_rp-0.19897337353e+01_rp*s&
            + 0.16597397329e+00_rp*t-0.16597397329e+00_rp*s*t&
            + 0.37248334214e+01_rp*s*s-0.16597397329e+00_rp*t*t
       shaga( 2)= 0.18736592735e+01_rp-0.54599331075e+01_rp*s&
            -0.54599331075e+01_rp*t+ 0.76156408161e+01_rp*s*t&
            + 0.37248334214e+01_rp*s*s+ 0.37248334214e+01_rp*t*t
       shaga( 3)= 0.13855958741e+00_rp+ 0.16597397329e+00_rp*s&
            -0.19897337353e+01_rp*t-0.16597397329e+00_rp*s*t&
            -0.16597397329e+00_rp*s*s+ 0.37248334214e+01_rp*t*t
       shaga( 4)=-0.63855958741e+00_rp+ 0.40859490520e+00_rp*s&
            + 0.79963036869e+01_rp*t-0.79963036869e+01_rp*s*t&
            + 0.35630540870e+00_rp*s*s-0.79963036869e+01_rp*t*t
       shaga( 5)= 0.12634072649e+00_rp-0.11212057226e+01_rp*s&
            -0.11212057226e+01_rp*t+ 0.87089145043e+01_rp*s*t&
            + 0.35630540870e+00_rp*s*s+ 0.35630540870e+00_rp*t*t
       shaga( 6)=-0.63855958741e+00_rp+ 0.79963036869e+01_rp*s&
            + 0.40859490520e+00_rp*t-0.79963036869e+01_rp*s*t&
            -0.79963036869e+01_rp*s*s+ 0.35630540870e+00_rp*t*t
    else if (ltopo==1 .and. ngaus==7 ) then
       shaga( 1)= 0.46660223518e-14_rp-0.22500000000e+01_rp*s&
            -0.22500000000e+01_rp*t+ 0.49500000000e+02_rp*s*t&
            + 0.22500000000e+01_rp*s*s+ 0.22500000000e+01_rp*t*t&
            -0.47250000000e+02_rp*(s*s*t+s*t*t)
       shaga( 2)=-0.76639777949e+00_rp+ 0.27797375097e+01_rp*s&
            + 0.87162291828e+01_rp*t-0.32857759237e+02_rp*s*t&
            -0.21106850116e+01_rp*s*s-0.87162291828e+01_rp*t*t&
            + 0.24141530054e+02_rp*(s*s*t+s*t*t)
       shaga( 3)=-0.76639777949e+00_rp+ 0.87162291828e+01_rp*s&
            + 0.27797375097e+01_rp*t-0.32857759237e+02_rp*s*t&
            -0.87162291828e+01_rp*s*s-0.21106850116e+01_rp*t*t&
            + 0.24141530054e+02_rp*(s*s*t+s*t*t)
       shaga( 4)=-0.97345281425e-01_rp+ 0.14416325135e+01_rp*s&
            + 0.14416325135e+01_rp*t-0.19646670894e+02_rp*s*t&
            -0.21106850116e+01_rp*s*s-0.21106850116e+01_rp*t*t&
            + 0.24141530054e+02_rp*(s*s*t+s*t*t)
       shaga( 5)= 0.26639777949e+00_rp-0.30297375097e+01_rp*s&
            -0.96622918276e+00_rp*t+ 0.93577592368e+01_rp*s*t&
            + 0.48606850116e+01_rp*s*s+ 0.96622918276e+00_rp*t*t&
            -0.83915300541e+01_rp*(s*s*t+s*t*t)
       shaga( 6)= 0.26639777949e+00_rp-0.96622918276e+00_rp*s&
            -0.30297375097e+01_rp*t+ 0.93577592368e+01_rp*s*t&
            + 0.96622918276e+00_rp*s*s+ 0.48606850116e+01_rp*t*t&
            -0.83915300541e+01_rp*(s*s*t+s*t*t)
       shaga( 7)= 0.20973452814e+01_rp-0.66916325135e+01_rp*s&
            -0.66916325135e+01_rp*t+ 0.17146670894e+02_rp*s*t&
            + 0.48606850116e+01_rp*s*s+ 0.48606850116e+01_rp*t*t&
            -0.83915300541e+01_rp*(s*s*t+s*t*t)
    else if (ltopo==1 .and. ngaus==13 ) then
       shaga(1:13)= 0.0_rp
    else
       ierr=1
    end if

    if( present(ierro) ) ierro = ierr
    
  end subroutine shaga2

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-02-13
  !> @brief   Extrapolation 3D
  !> @details This routine evaluates shape functions associated to gauss points
  !>          for 3D. The cases available so far are:
  !>
  !>          BRICKS    -->   NGAUS =   1   8  27 
  !>          SIMPLICES -->   NGAUS =   1   4
  !> 
  !-----------------------------------------------------------------------

  pure subroutine shaga3(x,ltopo,ngaus,shaga,ierro)

    real(rp),              intent(in)  :: x(3)
    integer(ip),           intent(in)  :: ngaus,ltopo
    real(rp),              intent(out) :: shaga(ngaus)
    integer(ip), optional, intent(out) :: ierro
    real(rp)                           :: s,t,z
    real(rp)                           :: r3,r6,a,b,c,d
    integer(ip)                        :: ierr

    s = x(1)
    t = x(2)
    z = x(3)
    ierr = 0
    !
    ! Non-available cases
    !
    if(ltopo==0) then
       if((ngaus/=1).and.(ngaus/=8).and.(ngaus/=27)) ierr = 1
    else if(ltopo==1) then
       if((ngaus/=1).and.(ngaus/=4)) ierr = 1
    end if
    !
    ! Available cases
    !
    if (ngaus==1) then

       shaga(1)=1.0_rp

    else if (ltopo==0 .and. ngaus==8) then           
       ! 
       ! 8 nodes Brick
       !
       r3=1.0_rp/sqrt(3.0_rp)
       shaga(1)=-1.125_rp*r3*(s-r3)*(t-r3)*(z-r3)
       shaga(2)= 1.125_rp*r3*(s-r3)*(t-r3)*(z+r3)
       shaga(3)= 1.125_rp*r3*(s-r3)*(t+r3)*(z-r3)
       shaga(4)=-1.125_rp*r3*(s-r3)*(t+r3)*(z+r3)
       shaga(5)= 1.125_rp*r3*(s+r3)*(t-r3)*(z-r3)
       shaga(6)=-1.125_rp*r3*(s+r3)*(t-r3)*(z+r3)
       shaga(7)=-1.125_rp*r3*(s+r3)*(t+r3)*(z-r3)
       shaga(8)= 1.125_rp*r3*(s+r3)*(t+r3)*(z+r3)
       !r3=sqrt(3.0_rp)
       !shaga(1)=-.046875_rp*r3*(s-r3)*(t-r3)*(z-r3)
       !shaga(2)= .046875_rp*r3*(s-r3)*(t-r3)*(z+r3)
       !shaga(3)= .046875_rp*r3*(s-r3)*(t+r3)*(z-r3)
       !shaga(4)=-.046875_rp*r3*(s-r3)*(t+r3)*(z+r3)
       !shaga(5)= .046875_rp*r3*(s+r3)*(t-r3)*(z-r3)
       !shaga(6)=-.046875_rp*r3*(s+r3)*(t-r3)*(z+r3)
       !shaga(7)=-.046875_rp*r3*(s+r3)*(t+r3)*(z-r3)
       !shaga(8)= .046875_rp*r3*(s+r3)*(t+r3)*(z+r3)

    else if (ltopo==0 .and. ngaus==27) then       
       !
       ! 27 node Brick
       !
       r6=sqrt(.6_rp)
       shaga( 1)= 125.0_rp/216.0_rp*(s-r6)*s*(t-r6)*t*(z-r6)*z
       shaga( 2)=-125.0_rp/108.0_rp*(s-r6)*s*(t-r6)*t*(z-r6)*(z+r6)
       shaga( 3)= 125.0_rp/216.0_rp*(s-r6)*s*(t-r6)*t*(z+r6)*z
       shaga( 4)= 125.0_rp/108.0_rp*(s-r6)*s*(t-r6)*(t+r6)*(z-r6)*z
       shaga( 5)=-125.0_rp/54.0_rp*(s-r6)*s*(t-r6)*(t+r6)*(z-r6)*(z+r6)
       shaga( 6)= 125.0_rp/108.0_rp*(s-r6)*s*(t-r6)*(t+r6)*(z+r6)*z
       shaga( 7)= 125.0_rp/216.0_rp*(s-r6)*s*(t+r6)*t*(z-r6)*z
       shaga( 8)=-125.0_rp/108.0_rp*(s-r6)*s*(t+r6)*t*(z-r6)*(z+r6)
       shaga( 9)= 125.0_rp/216.0_rp*(s-r6)*s*(t+r6)*t*(z+r6)*z
       shaga(10)=-125.0_rp/108.0_rp*(s-r6)*(s+r6)*(t-r6)*t*(z-r6)*z
       shaga(11)= 125.0_rp/54.0_rp*(s-r6)*(s+r6)*(t-r6)*t*(z-r6)*(z+r6)
       shaga(12)=-125.0_rp/108.0_rp*(s-r6)*(s+r6)*(t-r6)*t*(z+r6)*z
       shaga(13)=-125.0_rp/54.0_rp*(s-r6)*(s+r6)*(t-r6)*(t+r6)*(z-r6)*z
       shaga(14)= 125.0_rp/27.0_rp*(s-r6)*(s+r6)*(t-r6)*(t+r6)*(z-r6)*(z+r6)
       shaga(15)=-125.0_rp/54.0_rp*(s-r6)*(s+r6)*(t-r6)*(t+r6)*(z+r6)*z
       shaga(16)=-125.0_rp/108.0_rp*(s-r6)*(s+r6)*(t+r6)*t*(z-r6)*z
       shaga(17)= 125.0_rp/54.0_rp*(s-r6)*(s+r6)*(t+r6)*t*(z-r6)*(z+r6)
       shaga(18)=-125.0_rp/108.0_rp*(s-r6)*(s+r6)*(t+r6)*t*(z+r6)*z
       shaga(19)= 125.0_rp/216.0_rp*(s+r6)*s*(t-r6)*t*(z-r6)*z
       shaga(20)=-125.0_rp/108.0_rp*(s+r6)*s*(t-r6)*t*(z-r6)*(z+r6)
       shaga(21)= 125.0_rp/216.0_rp*(s+r6)*s*(t-r6)*t*(z+r6)*z
       shaga(22)= 125.0_rp/108.0_rp*(s+r6)*s*(t-r6)*(t+r6)*(z-r6)*z
       shaga(23)=-125.0_rp/54.0_rp*(s+r6)*s*(t-r6)*(t+r6)*(z-r6)*(z+r6)
       shaga(24)= 125.0_rp/108.0_rp*(s+r6)*s*(t-r6)*(t+r6)*(z+r6)*z
       shaga(25)= 125.0_rp/216.0_rp*(s+r6)*s*(t+r6)*t*(z-r6)*z
       shaga(26)=-125.0_rp/108.0_rp*(s+r6)*s*(t+r6)*t*(z-r6)*(z+r6)
       shaga(27)= 125.0_rp/216.0_rp*(s+r6)*s*(t+r6)*t*(z+r6)*z

    else if (ltopo==1 .and. ngaus==4) then 
       !
       ! 4 node tetrahedra
       !
       a= 1.927050983124842_rp
       b=-2.236067977499790_rp
       c=-2.236067977499789_rp  
       d=-2.236067977499789_rp
       shaga(1)=a+b*s+c*t+d*z
       a=-0.3090169943749473_rp
       b= 2.236067977499790_rp
       c= 0.0_rp
       d= 0.0_rp
       shaga(2)=a+b*s+c*t+d*z
       a=-0.3090169943749473_rp
       b= 0.0_rp
       c= 2.236067977499790_rp   
       d= 0.0_rp
       shaga(3)=a+b*s+c*t+d*z
       a=-0.3090169943749473_rp
       b= 0.0_rp
       c= 0.0_rp
       d= 2.236067977499790_rp
       shaga(4)=a+b*s+c*t+d*z

    else if (ltopo==2 .and. ngaus==6) then 
       !
       ! Prisms
       !
       if (s > 0.5_rp) then
          s= 2.0_rp - s
       else
          s= - s
       end if
       if (t > 0.5_rp) then
          t= 2.0_rp - t
       else
          t= - t
       end if
       if (z > 0.5_rp) then
          z= 2.0_rp - z
       else
          z= - z
       end if

       shaga(   1) = (1.0_rp-s-t)*(1.0_rp-z)
       shaga(   2) = s*(1.0_rp-z)
       shaga(   3) = t*(1.0_rp-z)
       shaga(   4) = (1.0_rp-s-t)*z
       shaga(   5) = s*z
       shaga(   6) = t*z

    else
       ierr=1
    end if
    
    if( present(ierro) ) ierro = ierr

  end subroutine shaga3

end module mod_extrapolation
!> @}

