!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Kinds_and_types
!> @{
!> @file    mod_lagranges.f90
!> @author  houzeaux
!> @date    2020-04-04
!> @brief   Iso-parametric arrays
!> @details Lagrange shape function, derivative and Hessian
!-----------------------------------------------------------------------

module mod_lagrange

  use def_kintyp_basic, only : ip,rp
  
  private

  public :: shape0
  public :: shape1
  public :: shape2
  public :: shape3

contains

  !-----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @brief   0D shape functions
  !> @details This routine evaluates shape functions and their first derivates 
  !>          for 0-d continuous with 1 nodes
  !
  !-----------------------------------------------------------------------

  pure subroutine shape0(nnode,shapf,ierro)
    
    integer(ip), intent(in)            :: nnode
    real(rp),    intent(out)           :: shapf(nnode)
    integer(ip), intent(out), optional :: ierro
    
    if( present(ierro) ) ierro    = 0 
    shapf(1) = 1.0_rp
    
  end subroutine shape0

  !-----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @brief   1D shape functions
  !> @details Evaluate shape functions and their first derivates
  !>          for 1-d continuos with 2 3 & 4 nodes
  !
  !-----------------------------------------------------------------------

  pure subroutine shape1(s,nnode,shapf,deriv,heslo,ierro)
    
    real(rp),    intent(in)            :: s
    integer(ip), intent(in)            :: nnode
    real(rp),    intent(out)           :: shapf(nnode)
    real(rp),    intent(out), optional :: deriv(1,nnode)
    real(rp),    intent(out), optional :: heslo(1,nnode)
    integer(ip), intent(out), optional :: ierro
    integer(ip)                        :: ierr

    ierr = 0
    if( present(heslo) ) then
       heslo(1,1:nnode) = 0.0_rp
    end if
    if( present(deriv) ) then
       deriv(1,1:nnode) = 0.0_rp
    end if

    if( nnode == 2 ) then

       shapf(1) =  0.5_rp*(1.0_rp-s)
       shapf(2) =  0.5_rp*(1.0_rp+s)
       if( present(deriv) ) then
          deriv(1,1) = -0.5_rp
          deriv(1,2) =  0.5_rp
       end if

    else if( nnode == 3 ) then

       shapf(1) =  0.5_rp*s*(s-1.0_rp)
       shapf(2) =  0.5_rp*s*(s+1.0_rp)
       shapf(3) = -(s+1.0_rp)*(s-1.0_rp)
       if( present(deriv) ) then
          deriv(1,1) =  s-0.5_rp
          deriv(1,2) =  s+0.5_rp
          deriv(1,3) = -2.0_rp*s
       end if
       if( present(heslo) ) then
          heslo(1,1) =   1.0_rp
          heslo(1,2) =   1.0_rp
          heslo(1,3) =  -2.0_rp
       end if

    else if( nnode == 4 ) then

       shapf(1) = -9.0_rp/16.0_rp*(s+1.0_rp/3.0_rp)*(s-1.0_rp/3.0_rp)*(s-1.0_rp)
       shapf(2) =  9.0_rp/16.0_rp*(s+1.0_rp)*(s+1.0_rp/3.0_rp)*(s-1.0_rp/3.0_rp)
       shapf(3) = 27.0_rp/16.0_rp*(s+1.0_rp)*(s-1.0_rp/3.0_rp)*(s-1.0_rp)
       shapf(4) = -27.0_rp/16.0_rp*(s+1.0_rp)*(s+1.0_rp/3.0_rp)*(s-1.0_rp)

       if( present(deriv) ) then
          deriv(1,1) = -9.0_rp/16.0_rp*((s-1.0_rp/3.0_rp)*(s-1.0_rp)+(s+1.0_rp/3.0_rp)&
               *(s-1.0_rp)+(s+1.0_rp/3.0_rp)*(s-1.0_rp/3.0_rp))
          deriv(1,2) =  9.0_rp/16.0_rp*((s+1.0_rp/3.0_rp)*(s-1.0_rp/3.0_rp)+(s+1.0_rp)&
               *(s-1.0_rp/3.0_rp)+(s+1.0_rp)*(s+1.0_rp/3.0_rp))
          deriv(1,3) = 27.0_rp/16.0_rp*((s-1.0_rp/3.0_rp)*(s-1.0_rp)+ (s+1.0_rp)&
               *(s-1.0_rp)+ (s+1.0_rp)*(s-1.0_rp/3.0_rp))
          deriv(1,4) = -27.0_rp/16.0_rp*((s+1.0_rp/3.0_rp)*(s-1.0_rp)+(s+1.0_rp)&
               *(s-1.0_rp)+(s+1.0_rp)*(s+1.0_rp/3.0_rp))
       end if
       if( present(heslo) ) then
          heslo(1,1) =  -9.0_rp/16.0_rp*(2.0_rp*(s-1.0_rp)+2.0_rp*(s+1.0_rp/3.0_rp)+2.0_rp*(s-1.0_rp/3.0_rp))
          heslo(1,2) =   9.0_rp/16.0_rp*(2.0_rp*(s-1.0_rp/3.0_rp)+2.0_rp*(s+1.0_rp/3.0_rp)+2.0_rp*(s+1.0_rp))
          heslo(1,3) =  27.0_rp/16.0_rp*(2.0_rp*(s-1.0_rp)+2.0_rp*(s-1.0_rp/3.0_rp)+2.0_rp*(s+1.0_rp))
          heslo(1,4) = -27.0_rp/16.0_rp*(2.0_rp*(s-1.0_rp)+2.0_rp*(s+1.0_rp/3.0_rp)+2.0_rp*(s+1.0_rp))
       end if
    else
       ierr = 1
    end if

    if( present(ierro) ) ierro = ierr
    
  end subroutine shape1

  !-----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @brief   2D shape functions
  !> @details Evaluate shape functions, derivates and Hessian
  !>
  !>          TRIANGLES       3   6  &  10  nodes
  !>          QUADRILATERALS  4   9  &  16  nodes
  !
  !-----------------------------------------------------------------------

  pure subroutine shape2(x,nnode,shapf,deriv,heslo,ierro)
   
    integer(ip), intent(in)            :: nnode
    real(rp),    intent(in)            :: x(2)
    real(rp),    intent(out)           :: shapf(nnode)
    real(rp),    intent(out), optional :: deriv(2,nnode)
    real(rp),    intent(out), optional :: heslo(3,nnode)
    integer(ip), intent(out), optional :: ierro
    integer(ip)                        :: ii,jj
    integer(ip)                        :: ierr
    real(rp)                           :: st,a1,a2,a3,ss,tt,s1,s2,s3,s4
    real(rp)                           :: t1,t2,t3,t4,s9,t9,c,a,s,t

    ierr = 0
    s = x(1)
    t = x(2)
    
    if( present(heslo) ) then
       do ii = 1,nnode
          do jj = 1,3
             heslo(jj,ii) = 0.0_rp
          end do
       end do
    end if
    if( present(deriv) ) then
       do ii = 1,nnode
          do jj = 1,2
             deriv(jj,ii) = 0.0_rp
          end do
       end do
    end if

    if( nnode == 3 ) then
       shapf(1)   =  1.0_rp-s-t
       shapf(2)   =  s
       shapf(3)   =  t
       if( present(deriv) ) then                              !  3
          deriv(1,1) = -1.0_rp                                !
          deriv(1,2) =  1.0_rp                                !
          deriv(1,3) =  0.0_rp                                !
          deriv(2,1) = -1.0_rp                                !  1       2
          deriv(2,2) =  0.0_rp
          deriv(2,3) =  1.0_rp
       end if

    else if( nnode == 4 ) then
       st         = s*t
       shapf(1)   = (1.0_rp-t-s+st)*0.25_rp                   !  4         3
       shapf(2)   = (1.0_rp-t+s-st)*0.25_rp                   !
       shapf(3)   = (1.0_rp+t+s+st)*0.25_rp                   !
       shapf(4)   = (1.0_rp+t-s-st)*0.25_rp                   !
       if( present(deriv) ) then                              !  1         2
          deriv(1,1) = (-1.0_rp+t)*0.25_rp
          deriv(1,2) = (+1.0_rp-t)*0.25_rp
          deriv(1,3) = (+1.0_rp+t)*0.25_rp
          deriv(1,4) = (-1.0_rp-t)*0.25_rp
          deriv(2,1) = (-1.0_rp+s)*0.25_rp
          deriv(2,2) = (-1.0_rp-s)*0.25_rp
          deriv(2,3) = (+1.0_rp+s)*0.25_rp
          deriv(2,4) = (+1.0_rp-s)*0.25_rp
       end if
       if( present(heslo) ) then
          heslo(3,1) =  0.25_rp
          heslo(3,2) = -0.25_rp
          heslo(3,3) =  0.25_rp
          heslo(3,4) = -0.25_rp
       end if

    else if( nnode == 6 ) then
       a1         = 1.0_rp-s-t
       a2         = s
       a3         = t
       shapf( 1)  = (2.0_rp*a1-1.0_rp)*a1                     !  3
       shapf( 2)  = (2.0_rp*a2-1.0_rp)*a2                     !
       shapf( 3)  = (2.0_rp*a3-1.0_rp)*a3                     !
       shapf( 4)  = 4.0_rp*a1*a2                              !  6      5
       shapf( 5)  = 4.0_rp*a2*a3                              !
       shapf( 6)  = 4.0_rp*a1*a3                              !
       if( present(deriv) ) then                              !  1     4     2
          deriv(1,1) = 1.0_rp-4.0_rp*a1
          deriv(1,2) = 4.0_rp*a2-1.0_rp
          deriv(1,3) = 0.0_rp
          deriv(1,4) = 4.0_rp*(a1-a2)
          deriv(1,5) = 4.0_rp*a3
          deriv(1,6) =-4.0_rp*a3
          deriv(2,1) = 1.0_rp-4.0_rp*a1
          deriv(2,2) = 0.0_rp
          deriv(2,3) = 4.0_rp*a3-1.0_rp
          deriv(2,4) =-4.0_rp*a2
          deriv(2,5) = 4.0_rp*a2
          deriv(2,6) = 4.0_rp*(a1-a3)
       end if
       if( present(heslo) ) then
          heslo(1,1) = 4.0_rp
          heslo(1,2) = 4.0_rp
          heslo(1,4) =-8.0_rp
          heslo(2,1) = 4.0_rp
          heslo(2,3) = 4.0_rp
          heslo(2,6) =-8.0_rp
          heslo(3,1) = 4.0_rp
          heslo(3,4) =-4.0_rp
          heslo(3,5) = 4.0_rp
          heslo(3,6) =-4.0_rp
       end if

    else if( nnode == 9 ) then
       ss         = s*s
       st         = s*t
       tt         = t*t
       s1         = s+1.0_rp
       t1         = t+1.0_rp
       s2         = s*2.0_rp
       t2         = t*2.0_rp
       s9         = s-1.0_rp
       t9         = t-1.0_rp                                   !  4      7      3
       shapf( 1)  = 0.25_rp*s9*st*t9                           !
       shapf( 2)  = 0.25_rp*s1*st*t9                           !
       shapf( 3)  = 0.25_rp*s1*st*t1                           !
       shapf( 4)  = 0.25_rp*s9*st*t1                           !  8      9      6
       shapf( 5)  = 0.5_rp*(1.0_rp-ss)*t*t9                    !
       shapf( 6)  = 0.5_rp*s*s1*(1.0_rp-tt)                    !
       shapf( 7)  = 0.5_rp*(1.0_rp-ss)*t*t1                    !
       shapf( 8)  = 0.5_rp*s*s9*(1.0_rp-tt)                    !  1      5      2
       shapf( 9)  = (1.0_rp-ss)*(1.0_rp-tt)
       if( present(deriv) ) then
          deriv(1,1) = 0.25_rp*t*t9*(-1.0_rp+s2)
          deriv(1,2) = 0.25_rp*(1.0_rp+s2)*t*t9
          deriv(1,3) = 0.25_rp*(1.0_rp+s2)*t*t1
          deriv(1,4) = 0.25_rp*(-1.0_rp+s2)*t*t1
          deriv(1,5) = -st*t9
          deriv(1,6) =  0.5_rp*(1.0_rp+s2)*(1.0_rp-tt)
          deriv(1,7) = -st*t1
          deriv(1,8) =  0.5_rp*(-1.0_rp+s2)*(1.0_rp-tt)
          deriv(1,9) = -s2*(1.0_rp-tt)
          deriv(2,1) =  0.25_rp*(-1.0_rp+t2)*s*s9
          deriv(2,2) =  0.25_rp*s*s1*(-1.0_rp+t2)
          deriv(2,3) =  0.25_rp*s*s1*(1.0_rp+t2)
          deriv(2,4) =  0.25_rp*s*s9*(1.0_rp+t2)
          deriv(2,5) =  0.5_rp*(1.0_rp-ss)*(-1.0_rp+t2)
          deriv(2,6) = -st*s1
          deriv(2,7) =  0.5_rp*(1.0_rp-ss)*(1.0_rp+t2)
          deriv(2,8) = -st*s9
          deriv(2,9) = -t2*(1.0_rp-ss)
       end if
       if( present(heslo) ) then
          heslo(1,1) = 0.5_rp*t*t9
          heslo(1,2) = 0.5_rp*t*t9
          heslo(1,3) = 0.5_rp*t*t1
          heslo(1,4) = 0.5_rp*t*t1
          heslo(1,5) =-t*t9
          heslo(1,6) = 1.0_rp-tt
          heslo(1,7) =-t*t1
          heslo(1,8) = 1.0_rp-tt
          heslo(1,9) =-2.0_rp*(1.0_rp-tt)
          heslo(2,1) = 0.5_rp*s*s9
          heslo(2,2) = 0.5_rp*s*s1
          heslo(2,3) = 0.5_rp*s*s1
          heslo(2,4) = 0.5_rp*s*s9
          heslo(2,5) = 1.0_rp-ss
          heslo(2,6) =-s*s1
          heslo(2,7) = 1.0_rp-ss
          heslo(2,8) =-s*s9
          heslo(2,9) =-2.0_rp*(1.0_rp-ss)
          heslo(3,1) = 0.25_rp*(-1.0_rp+t2)*(s9+s)
          heslo(3,2) = 0.25_rp*(-1.0_rp+t2)*(s1+s)
          heslo(3,3) = 0.25_rp*(1.0_rp+t2)*(s1+s)
          heslo(3,4) = 0.25_rp*(1.0_rp+t2)*(s9+s)
          heslo(3,5) =-s*(-1.0_rp+t2)
          heslo(3,6) =-t*s1-st
          heslo(3,7) =-s*(1.0_rp+t2)
          heslo(3,8) =-t*s9-st
          heslo(3,9) = s2*t2
       end if

    else if( nnode == 10 ) then
       c           = 9.0_rp/2.0_rp
       a1          = 1.0_rp-s-t
       a2          = 2.0_rp/3.0_rp-s-t
       a3          = 1.0_rp/3.0_rp-s-t
       shapf( 1)   = c*a1*a2*a3                                 !  3
       shapf( 2)   = c*(1.0_rp/3.0_rp-s)*(2.0_rp/3.0_rp-s)*s    !
       shapf( 3)   = c*(1.0_rp/3.0_rp-t)*(2.0_rp/3.0_rp-t)*t    !
       shapf( 4)   = 3.0_rp*c*a1*a2*s                           !  8    7
       shapf( 5)   =-3.0_rp*c*a1*(1.0_rp/3.0_rp-s)              !
       shapf( 6)   =-3.0_rp*c*(1.0_rp/3.0_rp-s)*s*t             !
       shapf( 7)   =-3.0_rp*c*s*(1.0_rp/3.0_rp-t)*t             !  9   10    6
       shapf( 8)   =-3.0_rp*c*a1*(1.0_rp/3.0_rp-t)*t            !
       shapf( 9)   = 3.0_rp*c*a1*a2*t                           !
       shapf(10)   = 6.0_rp*c*a1*s*t                            !  1    4    5    2
       if( present(deriv) ) then
          deriv(1, 1) =-c*(a1*a2+a1*a3+a2*a3)
          deriv(1, 2) =-c*((2.0_rp/3.0_rp-s)*s &
               &        + (1.0_rp/3.0_rp-s)*s-(1.0_rp/3.0_rp-s)*(2.0_rp/3.0_rp-s))
          deriv(1, 3) = 0.0_rp
          deriv(1, 4) = 3.0_rp*c*(a1*a2-a1*s-a2*s)
          deriv(1, 5) =-3.0_rp*c*(a1*(1.0_rp/3.0_rp-s)&
               &        -a1*s-(1.0_rp/3.0_rp-s)*s)
          deriv(1, 6) =-3.0_rp*c*((1.0_rp/3.0_rp-s)*t-s*t)
          deriv(1, 7) =-3.0_rp*c*((1.0_rp/3.0_rp-t)*t)
          deriv(1, 8) = 3.0_rp*c*((1.0_rp/3.0_rp-t)*t)
          deriv(1, 9) = 3.0_rp*c*(-a1*t-a2*t)
          deriv(1,10) = 6.0_rp*c*(a1*t-s*t)
          deriv(2, 1) =-c*(a1*a2+a1*a3+a2*a3)
          deriv(2, 2) = 0.0_rp
          deriv(2, 3) =-c*((2.0_rp/3.0_rp-t)*t&
               &        + (1.0_rp/3.0_rp-t)*t-(1.0_rp/3.0_rp-t)*(2.0_rp/3.0_rp-t))
          deriv(2, 4) = 3.0_rp*c*(-a1*s-a2*s)
          deriv(2, 5) =-3.0_rp*c*(-(1.0_rp/3.0_rp-s)*s)
          deriv(2, 6) =-3.0_rp*c*((1.0_rp/3.0_rp-s)*s)
          deriv(2, 7) =-3.0_rp*c*((1.0_rp/3.0_rp-t)*s-s*t)
          deriv(2, 8) =-3.0_rp*c*(-(1.0_rp/3.0_rp-t)*t&
               &        -a1*t+a1*(1.0_rp/3.0_rp-t))
          deriv(2, 9) = 3.0_rp*c*(-a1*t-a2*t+a1*a2)
          deriv(2,10) = 6.0_rp*c*(a1*s-s*t)
       end if
       if( present(heslo) ) then
          heslo(1, 1) = 2.0_rp*c*(a1+a2+a3)
          heslo(1, 2) =-c*(2.0_rp-6.0_rp*s)
          heslo(1, 3) = 0.0_rp
          heslo(1, 4) = 3.0_rp*c*(-2.0_rp*a1-2.0_rp*a2+2.0_rp*s)
          heslo(1, 5) =-3.0_rp*c*(-2.0_rp*a1-2.0_rp*(1.0_rp/3.0_rp-s)+2.0_rp*s)
          heslo(1, 6) = 3.0_rp*c*2.0_rp*t
          heslo(1, 7) = 0.0_rp
          heslo(1, 8) = 0.0_rp
          heslo(1, 9) = 3.0_rp*c*2.0_rp*t
          heslo(1,10) =-6.0_rp*c*2.0_rp*t
          heslo(2, 1) = c*(2.0_rp*a1+2.0_rp*a2+2.0_rp*a3)
          heslo(2, 2) = 0.0_rp
          heslo(2, 3) =-c*(2.0_rp-6.0_rp*t)
          heslo(2, 4) = 3.0_rp*c*2.0_rp*s
          heslo(2, 5) = 0.0_rp
          heslo(2, 6) = 0.0_rp
          heslo(2, 7) = 3.0_rp*c*2.0_rp*s
          heslo(2, 8) =-3.0_rp*c*(-2.0_rp*a1-2.0_rp*(1.0_rp/3.0_rp-t)+2.0_rp*t)
          heslo(2, 9) = 3.0_rp*c*(-2.0_rp*a1-2.0_rp*a2+2.0_rp*t)
          heslo(2,10) =-6.0_rp*c*2.0_rp*s
          heslo(3, 1) = 2.0_rp*c*(a1+a2+a3)
          heslo(3, 2) = 0.0_rp
          heslo(3, 3) = 0.0_rp
          heslo(3, 4) = 3.0_rp*c*(-a1-a2+2.0_rp*s)
          heslo(3, 5) =-3.0_rp*c*(-(1.0_rp/3.0_rp-s)+s)
          heslo(3, 6) =-3.0_rp*c*(1.0_rp/3.0_rp-2.0_rp*s)
          heslo(3, 7) =-3.0_rp*c*(1.0_rp/3.0_rp-2.0_rp*t)
          heslo(3, 8) = 3.0_rp*c*(1.0_rp/3.0_rp-2.0_rp*t)
          heslo(3, 9) = 3.0_rp*c*(-a1-a2+2.0_rp*t)
          heslo(3,10) = 6.0_rp*c*(a1-s-t)
       end if

    else if( nnode == 16 ) then
       a =81.0_rp/256.0_rp
       c =1.0_rp/3.0_rp
       s1=1.0_rp+s
       s2=c+s
       s3=c-s
       s4=1.0_rp-s
       t1=1.0_rp+t
       t2=c+t
       t3=c-t
       t4=1.0_rp-t
       shapf( 1) =   a*s2*s3*s4*t2*t3*t4                   ! 4    10    9    3
       shapf( 2) =   a*s1*s2*s3*t2*t3*t4                   !
       shapf( 3) =   a*s1*s2*s3*t1*t2*t3                   !
       shapf( 4) =   a*s2*s3*s4*t1*t2*t3                   ! 11   16   15    8
       shapf( 5) =-3.0_rp*a*s1*s3*s4*t2*t3*t4              !
       shapf( 6) =-3.0_rp*a*s1*s2*s4*t2*t3*t4              !
       shapf( 7) =-3.0_rp*a*s1*s2*s3*t1*t3*t4              ! 12   13   14    7
       shapf( 8) =-3.0_rp*a*s1*s2*s3*t1*t2*t4              !
       shapf( 9) =-3.0_rp*a*s1*s2*s4*t1*t2*t3              !
       shapf(10) =-3.0_rp*a*s1*s3*s4*t1*t2*t3              ! 1     5    6    2
       shapf(11) =-3.0_rp*a*s2*s3*s4*t1*t2*t4
       shapf(12) =-3.0_rp*a*s2*s3*s4*t1*t3*t4
       shapf(13) = 9.0_rp*a*s1*s3*s4*t1*t3*t4
       shapf(14) = 9.0_rp*a*s1*s2*s4*t1*t3*t4
       shapf(15) = 9.0_rp*a*s1*s2*s4*t1*t2*t4
       shapf(16) = 9.0_rp*a*s1*s3*s4*t1*t2*t4
       if( present(deriv) ) then
          deriv(1, 1)=  a *t2*t3*t4*(-s2*s3-s2*s4+s3*s4)
          deriv(1, 2)=  a *t2*t3*t4*(-s1*s2+s1*s3+s2*s3)
          deriv(1, 3)=  a *t1*t2*t3*(-s1*s2+s1*s3+s2*s3)
          deriv(1, 4)=  a *t1*t2*t3*(-s2*s3-s2*s4+s3*s4)
          deriv(1, 5)=-3.0_rp*a*t2*t3*t4*(-s1*s3-s1*s4+s3*s4)
          deriv(1, 6)=-3.0_rp*a*t2*t3*t4*(-s1*s2+s1*s4+s2*s4)
          deriv(1, 7)=-3.0_rp*a*t1*t3*t4*(-s1*s2+s1*s3+s2*s3)
          deriv(1, 8)=-3.0_rp*a*t1*t2*t4*(-s1*s2+s1*s3+s2*s3)
          deriv(1, 9)=-3.0_rp*a*t1*t2*t3*(-s1*s2+s1*s4+s2*s4)
          deriv(1,10)=-3.0_rp*a*t1*t2*t3*(-s1*s3-s1*s4+s3*s4)
          deriv(1,11)=-3.0_rp*a*t1*t2*t4*(-s2*s3-s2*s4+s3*s4)
          deriv(1,12)=-3.0_rp*a*t1*t3*t4*(-s2*s3-s2*s4+s3*s4)
          deriv(1,13)= 9.0_rp*a*t1*t3*t4*(-s1*s3-s1*s4+s3*s4)
          deriv(1,14)= 9.0_rp*a*t1*t3*t4*(-s1*s2+s1*s4+s2*s4)
          deriv(1,15)= 9.0_rp*a*t1*t2*t4*(-s1*s2+s1*s4+s2*s4)
          deriv(1,16)= 9.0_rp*a*t1*t2*t4*(-s1*s3-s1*s4+s3*s4)
          deriv(2, 1)=  a   *s2*s3*s4*(-t2*t3-t2*t4+t3*t4)
          deriv(2, 2)=  a   *s1*s2*s3*(-t2*t3-t2*t4+t3*t4)
          deriv(2, 3)=  a   *s1*s2*s3*(-t1*t2+t1*t3+t2*t3)
          deriv(2, 4)=  a   *s2*s3*s4*(-t1*t2+t1*t3+t2*t3)
          deriv(2, 5)= -3.0_rp*a *s1*s3*s4*(-t2*t3-t2*t4+t3*t4)
          deriv(2, 6)= -3.0_rp*a *s1*s2*s4*(-t2*t3-t2*t4+t3*t4)
          deriv(2, 7)= -3.0_rp*a *s1*s2*s3*(-t1*t3-t1*t4+t3*t4)
          deriv(2, 8)= -3.0_rp*a *s1*s2*s3*(-t1*t2+t1*t4+t2*t4)
          deriv(2, 9)= -3.0_rp*a *s1*s2*s4*(-t1*t2+t1*t3+t2*t3)
          deriv(2,10)= -3.0_rp*a *s1*s3*s4*(-t1*t2+t1*t3+t2*t3)
          deriv(2,11)= -3.0_rp*a *s2*s3*s4*(-t1*t2+t1*t4+t2*t4)
          deriv(2,12)= -3.0_rp*a *s2*s3*s4*(-t1*t3-t1*t4+t3*t4)
          deriv(2,13)=  9.0_rp*a *s1*s3*s4*(-t1*t3-t1*t4+t3*t4)
          deriv(2,14)=  9.0_rp*a *s1*s2*s4*(-t1*t3-t1*t4+t3*t4)
          deriv(2,15)=  9.0_rp*a *s1*s2*s4*(-t1*t2+t1*t4+t2*t4)
          deriv(2,16)=  9.0_rp*a *s1*s3*s4*(-t1*t2+t1*t4+t2*t4)
       end if
       if( present(heslo) ) then
          heslo(1, 1) =&
               a *t2*t3*t4*(2.0_rp*s2-2.0_rp*s3-2.0_rp*s4)
          heslo(1, 2) =&
               a *t2*t3*t4*(-2.0_rp*s1-2.0_rp*s2+2.0_rp*s3)
          heslo(1, 3) =&
               a *t1*t2*t3*(-2.0_rp*s1-2.0_rp*s2+2.0_rp*s3)
          heslo(1, 4) =&
               a *t1*t2*t3*(2.0_rp*s2-2.0_rp*s3-2.0_rp*s4)
          heslo(1, 5) =&
               -3.0_rp*a *t2*t3*t4*(2.0_rp*s1-2.0_rp*s3-2.0_rp*s4)
          heslo(1, 6) =&
               -3.0_rp*a *t2*t3*t4*(-2.0_rp*s1-2.0_rp*s2+2.0_rp*s4)
          heslo(1, 7) =&
               -3.0_rp*a *t1*t3*t4*(-2.0_rp*s1-2.0_rp*s2+2.0_rp*s3)
          heslo(1, 8) =&
               -3.0_rp*a *t1*t2*t4*(-2.0_rp*s1-2.0_rp*s2+2.0_rp*s3)
          heslo(1, 9) =&
               -3.0_rp*a *t1*t2*t3*(-2.0_rp*s1-2.0_rp*s2+2.0_rp*s4)
          heslo(1,10) =&
               -3.0_rp*a *t1*t2*t3*(2.0_rp*s1-2.0_rp*s3-2.0_rp*s4)
          heslo(1,11) =&
               -3.0_rp*a *t1*t2*t4*(2.0_rp*s2-2.0_rp*s3-2.0_rp*s4)
          heslo(1,12) =&
               -3.0_rp*a *t1*t3*t4*(2.0_rp*s2-2.0_rp*s3-2.0_rp*s4)
          heslo(1,13) =&
               9.0_rp*a *t1*t3*t4*(2.0_rp*s1-2.0_rp*s3-2.0_rp*s4)
          heslo(1,14) =&
               9.0_rp*a *t1*t3*t4*(-2.0_rp*s1-2.0_rp*s2+2.0_rp*s4)
          heslo(1,15) =&
               9.0_rp*a *t1*t2*t4*(-2.0_rp*s1-2.0_rp*s2+2.0_rp*s4)
          heslo(1,16) =&
               9.0_rp*a *t1*t2*t4*(2.0_rp*s1-2.0_rp*s3-2.0_rp*s4)
          heslo(2, 1) =&
               a *s2*s3*s4*(2.0_rp*t2-2.0_rp*t3-2.0_rp*t4)
          heslo(2, 2) =&
               a *s1*s2*s3*(2.0_rp*t2-2.0_rp*t3-2.0_rp*t4)
          heslo(2, 3) =&
               a *s1*s2*s3*(-2.0_rp*t1-2.0_rp*t2+2.0_rp*t3)
          heslo(2, 4) =&
               a *s2*s3*s4*(-2.0_rp*t1-2.0_rp*t2+2.0_rp*t3)
          heslo(2, 5) =&
               -3.0_rp*a *s1*s3*s4*(2.0_rp*t2-2.0_rp*t3-2.0_rp*t4)
          heslo(2, 6) =&
               -3.0_rp*a *s1*s2*s4*(2.0_rp*t2-2.0_rp*t3-2.0_rp*t4)
          heslo(2, 7) =&
               -3.0_rp*a *s1*s2*s3*(2.0_rp*t1-2.0_rp*t3-2.0_rp*t4)
          heslo(2, 8) =&
               -3.0_rp*a *s1*s2*s3*(-2.0_rp*t1-2.0_rp*t2+2.0_rp*t4)
          heslo(2, 9) =&
               -3.0_rp*a *s1*s2*s4*(-2.0_rp*t1-2.0_rp*t2+2.0_rp*t3)
          heslo(2,10) =&
               -3.0_rp*a *s1*s3*s4*(-2.0_rp*t1-2.0_rp*t2+2.0_rp*t3)
          heslo(2,11) =&
               -3.0_rp*a *s2*s3*s4*(-2.0_rp*t1-2.0_rp*t2+2.0_rp*t4)
          heslo(2,12) =&
               -3.0_rp*a *s2*s3*s4*(2.0_rp*t1-2.0_rp*t3-2.0_rp*t4)
          heslo(2,13) =&
               9.0_rp*a *s1*s3*s4*(2.0_rp*t1-2.0_rp*t3-2.0_rp*t4)
          heslo(2,14) =&
               9.0_rp*a *s1*s2*s4*(2.0_rp*t1-2.0_rp*t3-2.0_rp*t4)
          heslo(2,15) =&
               9.0_rp*a *s1*s2*s4*(-2.0_rp*t1-2.0_rp*t2+2.0_rp*t4)
          heslo(2,16) =&
               9.0_rp*a *s1*s3*s4*(-2.0_rp*t1-2.0_rp*t2+2.0_rp*t4)
          heslo(3, 1) =&
               a*(-s2*s3-s2*s4+s3*s4)*(-t2*t3-t2*t4+t3*t4)
          heslo(3, 2) =&
               a*(-s1*s2+s1*s3+s2*s3)*(-t2*t3-t2*t4+t3*t4)
          heslo(3, 3) =&
               a*(-s1*s2+s1*s3+s2*s3)*(-t1*t2+t1*t3+t2*t3)
          heslo(3, 4) =&
               a*(-s2*s3-s2*s4+s3*s4)*(-t1*t2+t1*t3+t2*t3)
          heslo(3, 5) =&
               -3.0_rp*a*(-s1*s3-s1*s4+s3*s4)*(-t2*t3-t2*t4+t3*t4)
          heslo(3, 6) =&
               -3.0_rp*a*(-s1*s2+s1*s4+s2*s4)*(-t2*t3-t2*t4+t3*t4)
          heslo(3, 7) =&
               -3.0_rp*a*(-s1*s2+s1*s3+s2*s3)*(-t1*t3-t1*t4+t3*t4)
          heslo(3, 8) =&
               -3.0_rp*a*(-s1*s2+s1*s3+s2*s3)*(-t1*t2+t1*t4+t2*t4)
          heslo(3, 9) =&
               -3.0_rp*a*(-s1*s2+s1*s4+s2*s4)*(-t1*t2+t1*t3+t2*t3)
          heslo(3,10) =&
               -3.0_rp*a*(-s1*s3-s1*s4+s3*s4)*(-t1*t2+t1*t3+t2*t3)
          heslo(3,11) =&
               -3.0_rp*a*(-s2*s3-s2*s4+s3*s4)*(-t1*t2+t1*t4+t2*t4)
          heslo(3,12) =&
               -3.0_rp*a*(-s2*s3-s2*s4+s3*s4)*(-t1*t3-t1*t4+t3*t4)
          heslo(3,13) =&
               9.0_rp*a*(-s1*s3-s1*s4+s3*s4)*(-t1*t3-t1*t4+t3*t4)
          heslo(3,14) =&
               9.0_rp*a*(-s1*s2+s1*s4+s2*s4)*(-t1*t3-t1*t4+t3*t4)
          heslo(3,15) =&
               9.0_rp*a*(-s1*s2+s1*s4+s2*s4)*(-t1*t2+t1*t4+t2*t4)
          heslo(3,16) =&
               9.0_rp*a*(-s1*s3-s1*s4+s3*s4)*(-t1*t2+t1*t4+t2*t4)
       end if
    else
       ierr = 1
    end if
    
    if( present(ierro) ) ierro = ierr

  end subroutine shape2

  !-----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @brief   3D shape functions
  !> @details Evaluate shape functions, derivative and Hessian
  !>          for 1-d continuos with 2 3 & 4 nodes
  !>
  !>          TETRAHEDRA:  4  10  &  20  nodes
  !>          HEXAHEDRA:   8  27  &  64  nodes
  !>          PRISM:       6  18         nodes
  !
  !
  !-----------------------------------------------------------------------

  pure subroutine shape3(x,nnode,shapf,deriv,heslo,ierro)
    
    integer(ip), intent(in)            :: nnode
    real(rp),    intent(in)            :: x(3)
    real(rp),    intent(out)           :: shapf(nnode)
    real(rp),    intent(out), optional :: deriv(3,nnode)
    real(rp),    intent(out), optional :: heslo(6,nnode)
    integer(ip), intent(out), optional :: ierro
    integer(ip)                        :: ii,jj, inode, idime, ind, ierr
    integer(ip)                        :: edg_perm(12,3), fac_perm(4,3)
    real(rp)                           :: L_vc(4), Laei(4,3), aux1, aux2
    real(rp)                           :: a1,a2,a3,a4,a,z1,z2,z3,z4,s1,s2,s3,s4
    real(rp)                           :: q1, q2, q3, q4, q5, q6, l13, l23, l33
    real(rp)                           :: dq1s, dq2s, dq3s, dq4s, dq5s, dq6s
    real(rp)                           :: dq1t, dq2t, dq3t, dq4t, dq5t, dq6t
    real(rp)                           :: dl1, dl2, dl3, s,t,z
    real(rp)                           :: t1,t2,t3,t4,sm,tm,zm,sq,tp,zp,s11,s21,s31,s41
    real(rp)                           :: t11,t21,t31,t41,z11,z21,z31,s12,s22,s32,s42
    real(rp)                           :: t12,t22,t32,t42,z41,z12,z22,z32,z42,sl,tl,zl
    real(rp)                           :: one8

    ierr = 0
    s = x(1)
    t = x(2)
    z = x(3)
    
    if( present(heslo) ) then
       do ii = 1,nnode
          do jj = 1,6
             heslo(jj,ii) = 0.0_rp
          end do
       end do
    end if
    if( present(deriv) ) then
       do ii = 1,nnode
          do jj = 1,3
             deriv(jj,ii) = 0.0_rp
          end do
       end do
    end if

    if( nnode == 4 ) then
       !
       ! Linear tetrahedron: s=[0:1], t=[0:1], z=[0:1]
       !
       shapf(   1) = 1.0_rp-s-t-z
       shapf(   2) = s
       shapf(   3) = t
       shapf(   4) = z
       if( present(deriv) ) then
          deriv(1, 1) =-1.0_rp
          deriv(2, 1) =-1.0_rp
          deriv(3, 1) =-1.0_rp
          deriv(1, 2) = 1.0_rp
          deriv(2, 3) = 1.0_rp
          deriv(3, 4) = 1.0_rp
       end if
       
    else if( nnode == 5 ) then
       !
       ! Linear Pyramid: s=[-1:1], t=[-1:1], z=[-1:1]
       !
       one8        =  0.125_rp
       shapf(   1) =  one8*(1.0_rp-s)*(1.0_rp-t)*(1.0_rp-z)
       shapf(   2) =  one8*(1.0_rp+s)*(1.0_rp-t)*(1.0_rp-z)
       shapf(   3) =  one8*(1.0_rp+s)*(1.0_rp+t)*(1.0_rp-z)
       shapf(   4) =  one8*(1.0_rp-s)*(1.0_rp+t)*(1.0_rp-z)
       shapf(   5) =  0.500_rp*(1.0_rp+z)
       if( present(deriv) ) then
          deriv(1, 1) = -one8*(1.0_rp-t)*(1.0_rp-z)
          deriv(2, 1) = -one8*(1.0_rp-s)*(1.0_rp-z)
          deriv(3, 1) = -one8*(1.0_rp-s)*(1.0_rp-t)
          deriv(1, 2) =  one8*(1.0_rp-t)*(1.0_rp-z)
          deriv(2, 2) = -one8*(1.0_rp+s)*(1.0_rp-z)
          deriv(3, 2) = -one8*(1.0_rp+s)*(1.0_rp-t)
          deriv(1, 3) =  one8*(1.0_rp+t)*(1.0_rp-z)
          deriv(2, 3) =  one8*(1.0_rp+s)*(1.0_rp-z)
          deriv(3, 3) = -one8*(1.0_rp+s)*(1.0_rp+t)
          deriv(1, 4) = -one8*(1.0_rp+t)*(1.0_rp-z)
          deriv(2, 4) =  one8*(1.0_rp-s)*(1.0_rp-z)
          deriv(3, 4) = -one8*(1.0_rp-s)*(1.0_rp+t)
          deriv(1, 5) =  0.0_rp
          deriv(2, 5) =  0.0_rp
          deriv(3, 5) =  0.5_rp
       end if
       if( present(heslo) ) then
          heslo(4, 1) =  one8*(1.0_rp-z)
          heslo(5, 1) =  one8*(1.0_rp-t)
          heslo(6, 1) =  one8*(1.0_rp-s)
          heslo(4, 2) = -one8*(1.0_rp-z)
          heslo(5, 2) = -one8*(1.0_rp-t)
          heslo(6, 2) =  one8*(1.0_rp+s)
          heslo(4, 3) =  one8*(1.0_rp-z)
          heslo(5, 3) = -one8*(1.0_rp+t)
          heslo(6, 3) = -one8*(1.0_rp+s)
          heslo(4, 4) = -one8*(1.0_rp-z)
          heslo(5, 4) =  one8*(1.0_rp+t)
          heslo(6, 4) = -one8*(1.0_rp-s)
       end if

    else if( nnode == 6 ) then
       !
       ! Linear Prism: s=[0:1], t=[0:1], z=[0:1]
       !
       shapf(   1) = (1.0_rp-s-t)*(1.0_rp-z)
       shapf(   2) = s*(1.0_rp-z)
       shapf(   3) = t*(1.0_rp-z)
       shapf(   4) = (1.0_rp-s-t)*z
       shapf(   5) = s*z
       shapf(   6) = t*z
       if( present(deriv) ) then
          deriv(1, 1) = z-1.0_rp
          deriv(2, 1) = z-1.0_rp
          deriv(3, 1) = s+t-1.0_rp
          deriv(1, 2) = 1.0_rp-z
          deriv(3, 2) = -s
          deriv(2, 3) = 1.0_rp-z
          deriv(3, 3) = -t
          deriv(1, 4) = -z
          deriv(2, 4) = -z
          deriv(3, 4) = 1.0_rp-s-t
          deriv(1, 5) = z
          deriv(3, 5) = s
          deriv(2, 6) = z
          deriv(3, 6) = t
       end if

  else if(nnode==18) then
     !
     ! Quadratic lagr. prism: s=[0:1], t=[0:1], z=[0:1]
     ! Product of TRI06(s,t) and O(2) Lagrange polynom.
     !
     a1 = 1.0_rp-s-t
     a2 = s
     a3 = t
     l13 = 2.0_rp*(z*z)-3.0_rp*z+1.0_rp
     l23 = 4.0_rp*(-(z*z)+z)
     l33 = 2.0_rp*(z*z)-z
     q1 = a1*(2.0_rp*a1-1.0_rp)
     q2 = a2*(2.0_rp*a2-1.0_rp)
     q3 = a3*(2.0_rp*a3-1.0_rp)
     q4 = 4.0_rp*a1*a2
     q5 = 4.0_rp*a2*a3
     q6 = 4.0_rp*a3*a1

     shapf(   1) = q1*l13
     shapf(   2) = q2*l13
     shapf(   3) = q3*l13
     shapf(   4) = q1*l33
     shapf(   5) = q2*l33
     shapf(   6) = q3*l33
     shapf(   7) = q4*l13
     shapf(   8) = q5*l13
     shapf(   9) = q6*l13
     shapf(  10) = q1*l23
     shapf(  11) = q2*l23
     shapf(  12) = q3*l23
     shapf(  13) = q4*l33
     shapf(  14) = q5*l33
     shapf(  15) = q6*l33
     shapf(  16) = q4*l23
     shapf(  17) = q5*l23
     shapf(  18) = q6*l23

     dq1s = -3.0_rp+4.0_rp*s+4.0_rp*t
     dq1t = -3.0_rp+4.0_rp*s+4.0_rp*t
     dq2s = 4.0_rp*s-1.0_rp
     dq2t = 0.0_rp
     dq3s = 0.0_rp
     dq3t = 4.0_rp*t-1.0_rp
     dq4s = 4.0_rp-8.0_rp*s-4.0_rp*t
     dq4t = -4.0_rp*s
     dq5s = 4.0_rp*t
     dq5t = 4.0_rp*s
     dq6s = -4.0_rp*t
     dq6t = 4.0_rp-4.0_rp*s-8.0_rp*t
     dl1 = 4.0_rp*z-3.0_rp
     dl2 = -8.0_rp*z+4.0_rp
     dl3 = 4.0_rp*z-1.0_rp

     deriv(1, 1) = dq1s*l13
     deriv(2, 1) = dq1t*l13
     deriv(3, 1) = q1*dl1
     deriv(1, 2) = dq2s*l13
     deriv(2, 2) = dq2t*l13
     deriv(3, 2) = q2*dl1
     deriv(1, 3) = dq3s*l13
     deriv(2, 3) = dq3t*l13
     deriv(3, 3) = q3*dl1
     deriv(1, 4) = dq1s*l33
     deriv(2, 4) = dq1t*l33
     deriv(3, 4) = q1*dl3
     deriv(1, 5) = dq2s*l33
     deriv(2, 5) = dq2t*l33
     deriv(3, 5) = q2*dl3
     deriv(1, 6) = dq3s*l33
     deriv(2, 6) = dq3t*l33
     deriv(3, 6) = q3*dl3
     deriv(1, 7) = dq4s*l13
     deriv(2, 7) = dq4t*l13
     deriv(3, 7) = q4*dl1
     deriv(1, 8) = dq5s*l13
     deriv(2, 8) = dq5t*l13
     deriv(3, 8) = q5*dl1
     deriv(1, 9) = dq6s*l13
     deriv(2, 9) = dq6t*l13
     deriv(3, 9) = q6*dl1
     deriv(1, 10) = dq1s*l23
     deriv(2, 10) = dq1t*l23
     deriv(3, 10) = q1*dl2
     deriv(1, 11) = dq2s*l23
     deriv(2, 11) = dq2t*l23
     deriv(3, 11) = q2*dl2
     deriv(1, 12) = dq3s*l23
     deriv(2, 12) = dq3t*l23
     deriv(3, 12) = q3*dl2
     deriv(1, 13) = dq4s*l33
     deriv(2, 13) = dq4t*l33
     deriv(3, 13) = q4*dl3
     deriv(1, 14) = dq5s*l33
     deriv(2, 14) = dq5t*l33
     deriv(3, 14) = q5*dl3
     deriv(1, 15) = dq6s*l33
     deriv(2, 15) = dq6t*l33
     deriv(3, 15) = q6*dl3
     deriv(1, 16) = dq4s*l23
     deriv(2, 16) = dq4t*l23
     deriv(3, 16) = q4*dl2
     deriv(1, 17) = dq5s*l23
     deriv(2, 17) = dq5t*l23
     deriv(3, 17) = q5*dl2
     deriv(1, 18) = dq6s*l23
     deriv(2, 18) = dq6t*l23
     deriv(3, 18) = q6*dl2

    else if( nnode == 10 ) then
       !
       ! Quadratic tetrahedron
       !
       a1= 1.0_rp-s-t-z
       a2=s
       a3=t
       a4=z
       shapf(   1) = (2.0_rp*a1-1.0_rp)*a1  ! [2 (1-s-t-z) - 1] (1-s-t-z)
       shapf(   2) = (2.0_rp*a2-1.0_rp)*a2  ! (2s-1) s
       shapf(   3) = (2.0_rp*a3-1.0_rp)*a3  ! (2t-1) t
       shapf(   4) = (2.0_rp*a4-1.0_rp)*a4  ! (2z-1) z
       shapf(   5) = 4.0_rp*a1*a2           ! 4 (1-s-t-z) s
       shapf(   6) = 4.0_rp*a2*a3           ! 4 s t
       shapf(   7) = 4.0_rp*a1*a3           ! 4 (1-s-t-z) t
       shapf(   8) = 4.0_rp*a1*a4           ! 4 (1-s-t-z) z
       shapf(   9) = 4.0_rp*a2*a4           ! 4 s z
       shapf(  10) = 4.0_rp*a3*a4           ! 4 t z
       if( present(deriv) ) then
          deriv(1, 1) = 1.0_rp-4.0_rp*a1
          deriv(2, 1) = 1.0_rp-4.0_rp*a1
          deriv(3, 1) = 1.0_rp-4.0_rp*a1
          deriv(1, 2) = 4.0_rp*a2-1.0_rp
          deriv(2, 3) = 4.0_rp*a3-1.0_rp
          deriv(3, 4) = 4.0_rp*a4-1.0_rp
          deriv(1, 5) = 4.0_rp*(a1-a2)
          deriv(2, 5) =-4.0_rp*a2
          deriv(3, 5) =-4.0_rp*a2
          deriv(1, 6) = 4.0_rp*a3
          deriv(2, 6) = 4.0_rp*a2
          deriv(1, 7) =-4.0_rp*a3
          deriv(2, 7) = 4.0_rp*(a1-a3)
          deriv(3, 7) =-4.0_rp*a3
          deriv(1, 8) =-4.0_rp*a4
          deriv(2, 8) =-4.0_rp*a4
          deriv(3, 8) = 4.0_rp*(a1-a4)
          deriv(1, 9) = 4.0_rp*a4
          deriv(3, 9) = 4.0_rp*a2
          deriv(2,10) = 4.0_rp*a4
          deriv(3,10) = 4.0_rp*a3
       end if
       if( present(heslo) ) then
          heslo(1:6,1) = 4.0_rp
          heslo(2, 3)  = 4.0_rp
          heslo(3, 4)  = 4.0_rp
          heslo(1, 2)  = 4.0_rp
          heslo(1, 5)  =-8.0_rp
          heslo(4, 5)  =-4.0_rp
          heslo(5, 5)  =-4.0_rp
          heslo(4, 6)  = 4.0_rp
          heslo(2, 7)  =-8.0_rp
          heslo(4, 7)  =-4.0_rp
          heslo(6, 7)  =-4.0_rp
          heslo(3, 8)  =-8.0_rp
          heslo(5, 8)  =-4.0_rp
          heslo(6, 8)  =-4.0_rp
          heslo(5, 9)  = 4.0_rp
          heslo(6,10)  = 4.0_rp
       end if

    else if(nnode==20) then
       !
       ! Cubic tetrahedron
       !        

       L_vc(1:4) = (/(1-s-t-z),s,t,z/)
       Laei(1,:) = (/-1.0_rp,-1.0_rp,-1.0_rp/)
       Laei(2,:) = (/ 1.0_rp, 0.0_rp, 0.0_rp/)
       Laei(3,:) = (/ 0.0_rp, 1.0_rp, 0.0_rp/)
       Laei(4,:) = (/ 0.0_rp, 0.0_rp, 1.0_rp/)

       edg_perm(1,:) = (/1,2,1/)
       edg_perm(2,:) = (/1,2,2/)
       edg_perm(3,:) = (/2,3,2/)
       edg_perm(4,:) = (/2,3,3/)
       edg_perm(5,:) = (/1,3,3/)
       edg_perm(6,:) = (/1,3,1/)
       edg_perm(7,:) = (/1,4,4/)
       edg_perm(8,:) = (/1,4,1/)
       edg_perm(9,:) = (/3,4,4/)
       edg_perm(10,:) = (/3,4,3/)
       edg_perm(11,:) = (/2,4,4/)
       edg_perm(12,:) = (/2,4,2/)

       fac_perm(1,:) = (/1,2,3/)
       fac_perm(2,:) = (/1,2,4/)
       fac_perm(3,:) = (/1,3,4/)
       fac_perm(4,:) = (/2,3,4/)

       ! Corner nodes

       do inode = 1,4
          shapf(inode) = 0.5_rp*(3.0_rp*L_vc(inode)-1.0_rp)*(3.0_rp*L_vc(inode)-2.0_rp)*L_vc(inode)
          do idime = 1,3
             deriv(idime,inode) = &
             0.5_rp*(27.0_rp*(L_vc(inode)**2)*Laei(inode,idime)-18*L_vc(inode)*Laei(inode,idime)+2.0_rp*Laei(inode,idime))
          end do
       end do

       ! Edge nodes

       do inode = 5,16
          ind = inode-4
          aux1 = L_vc(edg_perm(ind,1))*L_vc(edg_perm(ind,2))
          shapf(inode) = 4.5_rp*aux1*(3.0_rp*L_vc(edg_perm(ind,3))-1.0_rp)
          do idime = 1,3
             aux2 = Laei(edg_perm(ind,1),idime)*L_vc(edg_perm(ind,2))+Laei(edg_perm(ind,2),idime)*L_vc(edg_perm(ind,1))
             deriv(idime,inode) = 4.5_rp*(aux2*3.0_rp*L_vc(edg_perm(ind,3))+3.0_rp*aux1*Laei(edg_perm(ind,3),idime)-aux2)
          end do
       end do

       ! Face nodes

       do inode = 17,20
          ind = inode-16
          aux1 = L_vc(fac_perm(ind,1))*L_vc(fac_perm(ind,2))*L_vc(fac_perm(ind,3))
          shapf(inode) = 27.0_rp*aux1
          aux1 = L_vc(fac_perm(ind,1))*L_vc(fac_perm(ind,2))
          do idime = 1,3
             aux2 = Laei(fac_perm(ind,1),idime)*L_vc(fac_perm(ind,2))+Laei(fac_perm(ind,2),idime)*L_vc(fac_perm(ind,1))
             deriv(idime,inode) = 27.0_rp*(aux2*L_vc(fac_perm(ind,3))+aux1*Laei(fac_perm(ind,3),idime))
          end do
       end do

    else if( nnode == 8 ) then
       !
       ! Trilinear brick: s=[-1:1], t=[-1:1], z=[-1:1]
       !
       sm = 0.5_rp*(1.0_rp-s)
       tm = 0.5_rp*(1.0_rp-t)
       zm = 0.5_rp*(1.0_rp-z)
       sq = 0.5_rp*(1.0_rp+s)
       tp = 0.5_rp*(1.0_rp+t)
       zp = 0.5_rp*(1.0_rp+z)
       shapf(   1) = sm*tm*zm
       shapf(   2) = sq*tm*zm
       shapf(   3) = sq*tp*zm
       shapf(   4) = sm*tp*zm
       shapf(   5) = sm*tm*zp
       shapf(   6) = sq*tm*zp
       shapf(   7) = sq*tp*zp
       shapf(   8) = sm*tp*zp
       if( present(deriv) ) then
          deriv(1, 1) =-0.5_rp*tm*zm
          deriv(2, 1) =-0.5_rp*sm*zm
          deriv(3, 1) =-0.5_rp*sm*tm
          deriv(1, 2) = 0.5_rp*tm*zm
          deriv(2, 2) =-0.5_rp*sq*zm
          deriv(3, 2) =-0.5_rp*sq*tm
          deriv(1, 3) = 0.5_rp*tp*zm
          deriv(2, 3) = 0.5_rp*sq*zm
          deriv(3, 3) =-0.5_rp*sq*tp
          deriv(1, 4) =-0.5_rp*tp*zm
          deriv(2, 4) = 0.5_rp*sm*zm
          deriv(3, 4) =-0.5_rp*sm*tp
          deriv(1, 5) =-0.5_rp*tm*zp
          deriv(2, 5) =-0.5_rp*sm*zp
          deriv(3, 5) = 0.5_rp*sm*tm
          deriv(1, 6) = 0.5_rp*tm*zp
          deriv(2, 6) =-0.5_rp*sq*zp
          deriv(3, 6) = 0.5_rp*sq*tm
          deriv(1, 7) = 0.5_rp*tp*zp
          deriv(2, 7) = 0.5_rp*sq*zp
          deriv(3, 7) = 0.5_rp*sq*tp
          deriv(1, 8) =-0.5_rp*tp*zp
          deriv(2, 8) = 0.5_rp*sm*zp
          deriv(3, 8) = 0.5_rp*sm*tp
       end if
       if( present(heslo) ) then
          heslo(4, 1) = 0.25_rp*zm
          heslo(5, 1) = 0.25_rp*tm
          heslo(6, 1) = 0.25_rp*sm
          heslo(4, 2) =-0.25_rp*zm
          heslo(5, 2) =-0.25_rp*tm
          heslo(6, 2) = 0.25_rp*sq
          heslo(4, 3) = 0.25_rp*zm
          heslo(5, 3) =-0.25_rp*tp
          heslo(6, 3) =-0.25_rp*sq
          heslo(4, 4) =-0.25_rp*zm
          heslo(5, 4) = 0.25_rp*tp
          heslo(6, 4) =-0.25_rp*sm
          heslo(4, 5) = 0.25_rp*zp
          heslo(5, 5) =-0.25_rp*tm
          heslo(6, 5) =-0.25_rp*sm
          heslo(4, 6) =-0.25_rp*zp
          heslo(5, 6) = 0.25_rp*tm
          heslo(6, 6) =-0.25_rp*sq
          heslo(4, 7) = 0.25_rp*zp
          heslo(5, 7) = 0.25_rp*tp
          heslo(6, 7) = 0.25_rp*sq
          heslo(4, 8) =-0.25_rp*zp
          heslo(5, 8) =-0.25_rp*tp
          heslo(6, 8) = 0.25_rp*sm
       end if

    else if( nnode == 27 ) then
       !
       ! Triquadratic brick
       !
       sl=s*(s-1.0_rp)
       tl=t*(t-1.0_rp)
       zl=z*(z-1.0_rp)
       sq=s*(s+1.0_rp)
       tp=t*(t+1.0_rp)
       zp=z*(z+1.0_rp)
       s1= 2.0_rp*s-1.0_rp
       t1= 2.0_rp*t-1.0_rp
       z1= 2.0_rp*z-1.0_rp
       s2= 1.0_rp-s*s
       t2= 1.0_rp-t*t
       z2= 1.0_rp-z*z
       s3= 1.0_rp+2.0_rp*s
       t3= 1.0_rp+2.0_rp*t
       z3= 1.0_rp+2.0_rp*z
       s4=-2.0_rp*s
       t4=-2.0_rp*t
       z4=-2.0_rp*z
       shapf(   1) = 0.125_rp*sl*tl*zl
       shapf(   2) = 0.125_rp*sq*tl*zl
       shapf(   3) = 0.125_rp*sq*tp*zl
       shapf(   4) = 0.125_rp*sl*tp*zl
       shapf(   5) = 0.125_rp*sl*tl*zp
       shapf(   6) = 0.125_rp*sq*tl*zp
       shapf(   7) = 0.125_rp*sq*tp*zp
       shapf(   8) = 0.125_rp*sl*tp*zp
       shapf(   9) = 0.25_rp*s2*tl*zl
       shapf(  10) = 0.25_rp*sq*t2*zl
       shapf(  11) = 0.25_rp*s2*tp*zl
       shapf(  12) = 0.25_rp*sl*t2*zl
       shapf(  13) = 0.25_rp*sl*tl*z2
       shapf(  14) = 0.25_rp*sq*tl*z2
       shapf(  15) = 0.25_rp*sq*tp*z2
       shapf(  16) = 0.25_rp*sl*tp*z2
       shapf(  17) = 0.25_rp*s2*tl*zp
       shapf(  18) = 0.25_rp*sq*t2*zp
       shapf(  19) = 0.25_rp*s2*tp*zp
       shapf(  20) = 0.25_rp*sl*t2*zp
       shapf(  21) = 0.5_rp*s2*t2*zl
       shapf(  22) = 0.5_rp*s2*tl*z2
       shapf(  23) = 0.5_rp*sq*t2*z2
       shapf(  24) = 0.5_rp*s2*tp*z2
       shapf(  25) = 0.5_rp*sl*t2*z2
       shapf(  26) = 0.5_rp*s2*t2*zp
       shapf(  27) = s2*t2*z2
       if( present(deriv) ) then
          deriv(1, 1) = 0.125_rp*s1*tl*zl
          deriv(2, 1) = 0.125_rp*sl*t1*zl
          deriv(3, 1) = 0.125_rp*sl*tl*z1
          deriv(1, 2) = 0.125_rp*s3*tl*zl
          deriv(2, 2) = 0.125_rp*sq*t1*zl
          deriv(3, 2) = 0.125_rp*sq*tl*z1
          deriv(1, 3) = 0.125_rp*s3*tp*zl
          deriv(2, 3) = 0.125_rp*sq*t3*zl
          deriv(3, 3) = 0.125_rp*sq*tp*z1
          deriv(1, 4) = 0.125_rp*s1*tp*zl
          deriv(2, 4) = 0.125_rp*sl*t3*zl
          deriv(3, 4) = 0.125_rp*sl*tp*z1
          deriv(1, 5) = 0.125_rp*s1*tl*zp
          deriv(2, 5) = 0.125_rp*sl*t1*zp
          deriv(3, 5) = 0.125_rp*sl*tl*z3
          deriv(1, 6) = 0.125_rp*s3*tl*zp
          deriv(2, 6) = 0.125_rp*sq*t1*zp
          deriv(3, 6) = 0.125_rp*sq*tl*z3
          deriv(1, 7) = 0.125_rp*s3*tp*zp
          deriv(2, 7) = 0.125_rp*sq*t3*zp
          deriv(3, 7) = 0.125_rp*sq*tp*z3
          deriv(1, 8) = 0.125_rp*s1*tp*zp
          deriv(2, 8) = 0.125_rp*sl*t3*zp
          deriv(3, 8) = 0.125_rp*sl*tp*z3
          deriv(1, 9) = 0.25_rp*s4*tl*zl
          deriv(2, 9) = 0.25_rp*s2*t1*zl
          deriv(3, 9) = 0.25_rp*s2*tl*z1
          deriv(1,10) = 0.25_rp*s3*t2*zl
          deriv(2,10) = 0.25_rp*sq*t4*zl
          deriv(3,10) = 0.25_rp*sq*t2*z1
          deriv(1,11) = 0.25_rp*s4*tp*zl
          deriv(2,11) = 0.25_rp*s2*t3*zl
          deriv(3,11) = 0.25_rp*s2*tp*z1
          deriv(1,12) = 0.25_rp*s1*t2*zl
          deriv(2,12) = 0.25_rp*sl*t4*zl
          deriv(3,12) = 0.25_rp*sl*t2*z1
          deriv(1,13) = 0.25_rp*s1*tl*z2
          deriv(2,13) = 0.25_rp*sl*t1*z2
          deriv(3,13) = 0.25_rp*sl*tl*z4
          deriv(1,14) = 0.25_rp*s3*tl*z2
          deriv(2,14) = 0.25_rp*sq*t1*z2
          deriv(3,14) = 0.25_rp*sq*tl*z4
          deriv(1,15) = 0.25_rp*s3*tp*z2
          deriv(2,15) = 0.25_rp*sq*t3*z2
          deriv(3,15) = 0.25_rp*sq*tp*z4
          deriv(1,16) = 0.25_rp*s1*tp*z2
          deriv(2,16) = 0.25_rp*sl*t3*z2
          deriv(3,16) = 0.25_rp*sl*tp*z4
          deriv(1,17) = 0.25_rp*s4*tl*zp
          deriv(2,17) = 0.25_rp*s2*t1*zp
          deriv(3,17) = 0.25_rp*s2*tl*z3
          deriv(1,18) = 0.25_rp*s3*t2*zp
          deriv(2,18) = 0.25_rp*sq*t4*zp
          deriv(3,18) = 0.25_rp*sq*t2*z3
          deriv(1,19) = 0.25_rp*s4*tp*zp
          deriv(2,19) = 0.25_rp*s2*t3*zp
          deriv(3,19) = 0.25_rp*s2*tp*z3
          deriv(1,20) = 0.25_rp*s1*t2*zp
          deriv(2,20) = 0.25_rp*sl*t4*zp
          deriv(3,20) = 0.25_rp*sl*t2*z3
          deriv(1,21) = 0.5_rp*s4*t2*zl
          deriv(2,21) = 0.5_rp*s2*t4*zl
          deriv(3,21) = 0.5_rp*s2*t2*z1
          deriv(1,22) = 0.5_rp*s4*tl*z2
          deriv(2,22) = 0.5_rp*s2*t1*z2
          deriv(3,22) = 0.5_rp*s2*tl*z4
          deriv(1,23) = 0.5_rp*s3*t2*z2
          deriv(2,23) = 0.5_rp*sq*t4*z2
          deriv(3,23) = 0.5_rp*sq*t2*z4
          deriv(1,24) = 0.5_rp*s4*tp*z2
          deriv(2,24) = 0.5_rp*s2*t3*z2
          deriv(3,24) = 0.5_rp*s2*tp*z4
          deriv(1,25) = 0.5_rp*s1*t2*z2
          deriv(2,25) = 0.5_rp*sl*t4*z2
          deriv(3,25) = 0.5_rp*sl*t2*z4
          deriv(1,26) = 0.5_rp*s4*t2*zp
          deriv(2,26) = 0.5_rp*s2*t4*zp
          deriv(3,26) = 0.5_rp*s2*t2*z3
          deriv(1,27) = s4*t2*z2
          deriv(2,27) = s2*t4*z2
          deriv(3,27) = s2*t2*z4
       end if
       if( present(heslo) ) then
          heslo(1, 1) = 0.25_rp*tl*zl
          heslo(2, 1) = 0.25_rp*sl*zl
          heslo(3, 1) = 0.25_rp*sl*tl
          heslo(4, 1) = 0.125_rp*s1*t1*zl
          heslo(5, 1) = 0.125_rp*s1*tl*z1
          heslo(6, 1) = 0.125_rp*sl*t1*z1
          heslo(1, 2) = 0.25_rp*tl*zl
          heslo(2, 2) = 0.25_rp*sq*zl
          heslo(3, 2) = 0.25_rp*sq*tl
          heslo(4, 2) = 0.125_rp*s3*t1*zl
          heslo(5, 2) = 0.125_rp*s3*tl*z1
          heslo(6, 2) = 0.125_rp*sq*t1*z1
          heslo(1, 3) = 0.25_rp*tp*zl
          heslo(2, 3) = 0.25_rp*sq*zl
          heslo(3, 3) = 0.25_rp*sq*tp
          heslo(4, 3) = 0.125_rp*s3*t3*zl
          heslo(5, 3) = 0.125_rp*s3*tp*z1
          heslo(6, 3) = 0.125_rp*sq*t3*z1
          heslo(1, 4) = 0.25_rp*tp*zl
          heslo(2, 4) = 0.25_rp*sl*zl
          heslo(3, 4) = 0.25_rp*sl*tp
          heslo(4, 4) = 0.125_rp*s1*t3*zl
          heslo(5, 4) = 0.125_rp*s1*tp*z1
          heslo(6, 4) = 0.125_rp*sl*t3*z1
          heslo(1, 5) = 0.25_rp*tl*zp
          heslo(2, 5) = 0.25_rp*sl*zp
          heslo(3, 5) = 0.25_rp*sl*tl
          heslo(4, 5) = 0.125_rp*s1*t1*zp
          heslo(5, 5) = 0.125_rp*s1*tl*z3
          heslo(6, 5) = 0.125_rp*sl*t1*z3
          heslo(1, 6) = 0.25_rp*tl*zp
          heslo(2, 6) = 0.25_rp*sq*zp
          heslo(3, 6) = 0.25_rp*sq*tl
          heslo(4, 6) = 0.125_rp*s3*t1*zp
          heslo(5, 6) = 0.125_rp*s3*tl*z3
          heslo(6, 6) = 0.125_rp*sq*t1*z3
          heslo(1, 7) = 0.25_rp*tp*zp
          heslo(2, 7) = 0.25_rp*sq*zp
          heslo(3, 7) = 0.25_rp*sq*tp
          heslo(4, 7) = 0.125_rp*s3*t3*zp
          heslo(5, 7) = 0.125_rp*s3*tp*z3
          heslo(6, 7) = 0.125_rp*sq*t3*z3
          heslo(1, 8) = 0.25_rp*tp*zp
          heslo(2, 8) = 0.25_rp*sl*zp
          heslo(3, 8) = 0.25_rp*sl*tp
          heslo(4, 8) = 0.125_rp*s1*t3*zp
          heslo(5, 8) = 0.125_rp*s1*tp*z3
          heslo(6, 8) = 0.125_rp*sl*t3*z3
          heslo(1, 9) =-0.5_rp*tl*zl
          heslo(2, 9) = 0.5_rp*s2*zl
          heslo(3, 9) = 0.5_rp*s2*tl
          heslo(4, 9) = 0.25_rp*s4*t1*zl
          heslo(5, 9) = 0.25_rp*s4*tl*z1
          heslo(6, 9) = 0.25_rp*s2*t1*z1
          heslo(1,10) = 0.5_rp*t2*zl
          heslo(2,10) =-0.5_rp*sq*zl
          heslo(3,10) = 0.5_rp*sq*t2
          heslo(4,10) = 0.25_rp*s3*t4*zl
          heslo(5,10) = 0.25_rp*s3*t2*z1
          heslo(6,10) = 0.25_rp*sq*t4*z1
          heslo(1,11) =-0.5_rp*tp*zl
          heslo(2,11) = 0.5_rp*s2*zl
          heslo(3,11) = 0.5_rp*s2*tp
          heslo(4,11) = 0.25_rp*s4*t3*zl
          heslo(5,11) = 0.25_rp*s4*tp*z1
          heslo(6,11) = 0.25_rp*s2*t3*z1
          heslo(1,12) = 0.5_rp*t2*zl
          heslo(2,12) =-0.5_rp*sl*zl
          heslo(3,12) = 0.5_rp*sl*t2
          heslo(4,12) = 0.25_rp*s1*t4*zl
          heslo(5,12) = 0.25_rp*s1*t2*z1
          heslo(6,12) = 0.25_rp*sl*t4*z1
          heslo(1,13) = 0.5_rp*tl*z2
          heslo(2,13) = 0.5_rp*sl*z2
          heslo(3,13) =-0.5_rp*sl*tl
          heslo(4,13) = 0.25_rp*s1*t1*z2
          heslo(5,13) = 0.25_rp*s1*tl*z4
          heslo(6,13) = 0.25_rp*sl*t1*z4
          heslo(1,14) = 0.5_rp*tl*z2
          heslo(2,14) = 0.5_rp*sq*z2
          heslo(3,14) =-0.5_rp*sq*tl
          heslo(4,14) = 0.25_rp*s3*t1*z2
          heslo(5,14) = 0.25_rp*s3*tl*z4
          heslo(6,14) = 0.25_rp*sq*t1*z4
          heslo(1,15) = 0.5_rp*tp*z2
          heslo(2,15) = 0.5_rp*sq*z2
          heslo(3,15) =-0.5_rp*sq*tp
          heslo(4,15) = 0.25_rp*s3*t3*z2
          heslo(5,15) = 0.25_rp*s3*tp*z4
          heslo(6,15) = 0.25_rp*sq*t3*z4
          heslo(1,16) = 0.5_rp*tp*z2
          heslo(2,16) = 0.5_rp*sl*z2
          heslo(3,16) =-0.5_rp*sl*tp
          heslo(4,16) = 0.25_rp*s1*t3*z2
          heslo(5,16) = 0.25_rp*s1*tp*z4
          heslo(6,16) = 0.25_rp*sl*t3*z4
          heslo(1,17) =-0.5_rp*tl*zp
          heslo(2,17) = 0.5_rp*s2*zp
          heslo(3,17) = 0.5_rp*s2*tl
          heslo(4,17) = 0.25_rp*s4*t1*zp
          heslo(5,17) = 0.25_rp*s4*tl*z3
          heslo(6,17) = 0.25_rp*s2*t1*z3
          heslo(1,18) = 0.5_rp*t2*zp
          heslo(2,18) =-0.5_rp*sq*zp
          heslo(3,18) = 0.5_rp*sq*t2
          heslo(4,18) = 0.25_rp*s3*t4*zp
          heslo(5,18) = 0.25_rp*s3*t2*z3
          heslo(6,18) = 0.25_rp*sq*t4*z3
          heslo(1,19) =-0.5_rp*tp*zp
          heslo(2,19) = 0.5_rp*s2*zp
          heslo(3,19) = 0.5_rp*s2*tp
          heslo(4,19) = 0.25_rp*s4*t3*zp
          heslo(5,19) = 0.25_rp*s4*tp*z3
          heslo(6,19) = 0.25_rp*s2*t3*z3
          heslo(1,20) = 0.5_rp*t2*zp
          heslo(2,20) =-0.5_rp*sl*zp
          heslo(3,20) = 0.5_rp*sl*t2
          heslo(4,20) = 0.25_rp*s1*t4*zp
          heslo(5,20) = 0.25_rp*s1*t2*z3
          heslo(6,20) = 0.25_rp*sl*t4*z3
          heslo(1,21) =-t2*zl
          heslo(2,21) =-s2*zl
          heslo(3,21) = s2*t2
          heslo(4,21) = 0.5_rp*s4*t4*zl
          heslo(5,21) = 0.5_rp*s4*t2*z1
          heslo(6,21) = 0.5_rp*s2*t4*z1
          heslo(1,22) =-tl*z2
          heslo(2,22) = s2*z2
          heslo(3,22) =-s2*tl
          heslo(4,22) = 0.5_rp*s4*t1*z2
          heslo(5,22) = 0.5_rp*s4*tl*z4
          heslo(6,22) = 0.5_rp*s2*t1*z4
          heslo(1,23) = t2*z2
          heslo(2,23) =-sq*z2
          heslo(3,23) =-sq*t2
          heslo(4,23) = 0.5_rp*s3*t4*z2
          heslo(5,23) = 0.5_rp*s3*t2*z4
          heslo(6,23) = 0.5_rp*sq*t4*z4
          heslo(1,24) =-tp*z2
          heslo(2,24) = s2*z2
          heslo(3,24) =-s2*tp
          heslo(4,24) = 0.5_rp*s4*t3*z2
          heslo(5,24) = 0.5_rp*s4*tp*z4
          heslo(6,24) = 0.5_rp*s2*t3*z4
          heslo(1,25) = t2*z2
          heslo(2,25) =-sl*z2
          heslo(3,25) =-sl*t2
          heslo(4,25) = 0.5_rp*s1*t4*z2
          heslo(5,25) = 0.5_rp*s1*t2*z4
          heslo(6,25) = 0.5_rp*sl*t4*z4
          heslo(1,26) =-t2*zp
          heslo(2,26) =-s2*zp
          heslo(3,26) = s2*t2
          heslo(4,26) = 0.5_rp*s4*t4*zp
          heslo(5,26) = 0.5_rp*s4*t2*z3
          heslo(6,26) = 0.5_rp*s2*t4*z3
          heslo(1,27) =-2.0_rp*t2*z2
          heslo(2,27) =-2.0_rp*s2*z2
          heslo(3,27) =-2.0_rp*s2*t2
          heslo(4,27) = s4*t4*z2
          heslo(5,27) = s4*t2*z4
          heslo(6,27) = s2*t4*z4
       end if

    else if( nnode == 64 ) then
       !
       ! Tricubic brick
       !
       a  = 729.0_rp/4096.0_rp
       
       s1 = (1.0_rp/3.0_rp+s) * (1.0_rp/3.0_rp-s) * (1.0_rp-s)
       s2 = (1.0_rp+s)        * (1.0_rp/3.0_rp-s) * (1.0_rp-s)
       s3 = (1.0_rp+s)        * (1.0_rp/3.0_rp+s) * (1.0_rp-s)
       s4 = (1.0_rp+s)        * (1.0_rp/3.0_rp+s) * (1.0_rp/3.0_rp-s)
       
       t1 = (1.0_rp/3.0_rp+t) * (1.0_rp/3.0_rp-t) * (1.0_rp-t)
       t2 = (1.0_rp+t)        * (1.0_rp/3.0_rp-t) * (1.0_rp-t)
       t3 = (1.0_rp+t)        * (1.0_rp/3.0_rp+t) * (1.0_rp-t)
       t4 = (1.0_rp+t)        * (1.0_rp/3.0_rp+t) * (1.0_rp/3.0_rp-t)
       
       z1 = (1.0_rp/3.0_rp+z) * (1.0_rp/3.0_rp-z) * (1.0_rp-z)
       z2 = (1.0_rp+z)        * (1.0_rp/3.0_rp-z) * (1.0_rp-z)
       z3 = (1.0_rp+z)        * (1.0_rp/3.0_rp+z) * (1.0_rp-z)
       z4 = (1.0_rp+z)        * (1.0_rp/3.0_rp+z) * (1.0_rp/3.0_rp-z)
       
       s11 = (1.0_rp/3.0_rp-s)*(1.0_rp-s)        - (1.0_rp/3.0_rp+s)*(1.0_rp-s) - (1.0_rp/3.0_rp+s)*(1.0_rp/3.0_rp-s)
       s21 = (1.0_rp/3.0_rp-s)*(1.0_rp-s)        - (1.0_rp+s)*(1.0_rp-s)        - (1.0_rp+s)*(1.0_rp/3.0_rp-s)
       s31 = (1.0_rp/3.0_rp+s)*(1.0_rp-s)        + (1.0_rp+s)*(1.0_rp-s)        - (1.0_rp+s)*(1.0_rp/3.0_rp+s)
       s41 = (1.0_rp/3.0_rp+s)*(1.0_rp/3.0_rp-s) + (1.0_rp+s)*(1.0_rp/3.0_rp-s) - (1.0_rp+s)*(1.0_rp/3.0_rp+s)
       
       t11 = (1.0_rp/3.0_rp-t)*(1.0_rp-t)        - (1.0_rp/3.0_rp+t)*(1.0_rp-t) - (1.0_rp/3.0_rp+t)*(1.0_rp/3.0_rp-t)
       t21 = (1.0_rp/3.0_rp-t)*(1.0_rp-t)        - (1.0_rp+t)*(1.0_rp-t)        - (1.0_rp+t)*(1.0_rp/3.0_rp-t)
       t31 = (1.0_rp/3.0_rp+t)*(1.0_rp-t)        + (1.0_rp+t)*(1.0_rp-t)        - (1.0_rp+t)*(1.0_rp/3.0_rp+t)
       t41 = (1.0_rp/3.0_rp+t)*(1.0_rp/3.0_rp-t) + (1.0_rp+t)*(1.0_rp/3.0_rp-t) - (1.0_rp+t)*(1.0_rp/3.0_rp+t)
       
       z11 = (1.0_rp/3.0_rp-z)*(1.0_rp-z)        - (1.0_rp/3.0_rp+z)*(1.0_rp-z) - (1.0_rp/3.0_rp+z)*(1.0_rp/3.0_rp-z)
       z21 = (1.0_rp/3.0_rp-z)*(1.0_rp-z)        - (1.0_rp+z)*(1.0_rp-z)        - (1.0_rp+z)*(1.0_rp/3.0_rp-z)
       z31 = (1.0_rp/3.0_rp+z)*(1.0_rp-z)        + (1.0_rp+z)*(1.0_rp-z)        - (1.0_rp+z)*(1.0_rp/3.0_rp+z)
       z41 = (1.0_rp/3.0_rp+z)*(1.0_rp/3.0_rp-z) + (1.0_rp+z)*(1.0_rp/3.0_rp-z) - (1.0_rp+z)*(1.0_rp/3.0_rp+z)
       
       shapf(   1) =          -a * s1 * t1 * z1
       shapf(   2) =          -a * s4 * t1 * z1
       shapf(   3) =          -a * s4 * t4 * z1
       shapf(   4) =          -a * s1 * t4 * z1
       shapf(   5) =          -a * s1 * t1 * z4
       shapf(   6) =          -a * s4 * t1 * z4
       shapf(   7) =          -a * s4 * t4 * z4
       shapf(   8) =          -a * s1 * t4 * z4
       shapf(   9) =  3.0_rp * a * s2 * t1 * z1
       shapf(  10) =  3.0_rp * a * s3 * t1 * z1
       shapf(  11) =  3.0_rp * a * s4 * t2 * z1
       shapf(  12) =  3.0_rp * a * s4 * t3 * z1
       shapf(  13) =  3.0_rp * a * s3 * t4 * z1
       shapf(  14) =  3.0_rp * a * s2 * t4 * z1
       shapf(  15) =  3.0_rp * a * s1 * t3 * z1
       shapf(  16) =  3.0_rp * a * s1 * t2 * z1
       shapf(  17) =  3.0_rp * a * s1 * t1 * z2
       shapf(  18) =  3.0_rp * a * s4 * t1 * z2
       shapf(  19) =  3.0_rp * a * s4 * t4 * z2
       shapf(  20) =  3.0_rp * a * s1 * t4 * z2
       shapf(  21) =  3.0_rp * a * s1 * t1 * z3
       shapf(  22) =  3.0_rp * a * s4 * t1 * z3
       shapf(  23) =  3.0_rp * a * s4 * t4 * z3
       shapf(  24) =  3.0_rp * a * s1 * t4 * z3
       shapf(  25) =  3.0_rp * a * s2 * t1 * z4
       shapf(  26) =  3.0_rp * a * s3 * t1 * z4
       shapf(  27) =  3.0_rp * a * s4 * t2 * z4
       shapf(  28) =  3.0_rp * a * s4 * t3 * z4
       shapf(  29) =  3.0_rp * a * s3 * t4 * z4
       shapf(  30) =  3.0_rp * a * s2 * t4 * z4
       shapf(  31) =  3.0_rp * a * s1 * t3 * z4
       shapf(  32) =  3.0_rp * a * s1 * t2 * z4
       shapf(  33) = -9.0_rp * a * s2 * t2 * z1
       shapf(  34) = -9.0_rp * a * s3 * t2 * z1
       shapf(  35) = -9.0_rp * a * s3 * t3 * z1
       shapf(  36) = -9.0_rp * a * s2 * t3 * z1
       shapf(  37) = -9.0_rp * a * s2 * t1 * z2
       shapf(  38) = -9.0_rp * a * s3 * t1 * z2
       shapf(  39) = -9.0_rp * a * s4 * t2 * z2
       shapf(  40) = -9.0_rp * a * s4 * t3 * z2
       shapf(  41) = -9.0_rp * a * s3 * t4 * z2
       shapf(  42) = -9.0_rp * a * s2 * t4 * z2
       shapf(  43) = -9.0_rp * a * s1 * t3 * z2
       shapf(  44) = -9.0_rp * a * s1 * t2 * z2
       shapf(  45) = -9.0_rp * a * s2 * t1 * z3
       shapf(  46) = -9.0_rp * a * s3 * t1 * z3
       shapf(  47) = -9.0_rp * a * s4 * t2 * z3
       shapf(  48) = -9.0_rp * a * s4 * t3 * z3
       shapf(  49) = -9.0_rp * a * s3 * t4 * z3
       shapf(  50) = -9.0_rp * a * s2 * t4 * z3
       shapf(  51) = -9.0_rp * a * s1 * t3 * z3
       shapf(  52) = -9.0_rp * a * s1 * t2 * z3
       shapf(  53) = -9.0_rp * a * s2 * t2 * z4
       shapf(  54) = -9.0_rp * a * s3 * t2 * z4
       shapf(  55) = -9.0_rp * a * s3 * t3 * z4
       shapf(  56) = -9.0_rp * a * s2 * t3 * z4
       shapf(  57) = 27.0_rp * a * s2 * t2 * z2
       shapf(  58) = 27.0_rp * a * s3 * t2 * z2
       shapf(  59) = 27.0_rp * a * s3 * t3 * z2
       shapf(  60) = 27.0_rp * a * s2 * t3 * z2
       shapf(  61) = 27.0_rp * a * s2 * t2 * z3
       shapf(  62) = 27.0_rp * a * s3 * t2 * z3
       shapf(  63) = 27.0_rp * a * s3 * t3 * z3
       shapf(  64) = 27.0_rp * a * s2 * t3 * z3

       ! bad sign: 1,2,3,4,5,6,7,8
       if( present(deriv) ) then
          deriv(1, 1) =          -a * s11 * t1  * z1
          deriv(2, 1) =          -a * s1  * t11 * z1
          deriv(3, 1) =          -a * s1  * t1  * z11
          
          deriv(1, 2) =          -a * s41 * t1  * z1
          deriv(2, 2) =          -a * s4  * t11 * z1
          deriv(3, 2) =          -a * s4  * t1  * z11
          
          deriv(1, 3) =          -a * s41 * t4  * z1
          deriv(2, 3) =          -a * s4  * t41 * z1
          deriv(3, 3) =          -a * s4  * t4  * z11
          
          deriv(1, 4) =          -a * s11 * t4  * z1
          deriv(2, 4) =          -a * s1  * t41 * z1
          deriv(3, 4) =          -a * s1  * t4  * z11
          
          deriv(1, 5) =          -a * s11 * t1  * z4
          deriv(2, 5) =          -a * s1  * t11 * z4
          deriv(3, 5) =          -a * s1  * t1  * z41
          
          deriv(1, 6) =          -a * s41 * t1  * z4
          deriv(2, 6) =          -a * s4  * t11 * z4
          deriv(3, 6) =          -a * s4  * t1  * z41
          
          deriv(1, 7) =          -a * s41 * t4  * z4
          deriv(2, 7) =          -a * s4  * t41 * z4
          deriv(3, 7) =          -a * s4  * t4  * z41
          
          deriv(1, 8) =          -a * s11 * t4  * z4
          deriv(2, 8) =          -a * s1  * t41 * z4
          deriv(3, 8) =          -a * s1  * t4  * z41
          
          deriv(1, 9) =  3.0_rp * a * s21 * t1  * z1
          deriv(2, 9) =  3.0_rp * a * s2  * t11 * z1
          deriv(3, 9) =  3.0_rp * a * s2  * t1  * z11
          
          deriv(1,10) =  3.0_rp * a * s31 * t1  * z1
          deriv(2,10) =  3.0_rp * a * s3  * t11 * z1
          deriv(3,10) =  3.0_rp * a * s3  * t1  * z11
          
          deriv(1,11) =  3.0_rp * a * s41 * t2  * z1
          deriv(2,11) =  3.0_rp * a * s4  * t21 * z1
          deriv(3,11) =  3.0_rp * a * s4  * t2  * z11
          
          deriv(1,12) =  3.0_rp * a * s41 * t3  * z1
          deriv(2,12) =  3.0_rp * a * s4  * t31 * z1
          deriv(3,12) =  3.0_rp * a * s4  * t3  * z11
          
          deriv(1,13) =  3.0_rp * a * s31 * t4  * z1
          deriv(2,13) =  3.0_rp * a * s3  * t41 * z1
          deriv(3,13) =  3.0_rp * a * s3  * t4  * z11
          
          deriv(1,14) =  3.0_rp * a * s21 * t4  * z1
          deriv(2,14) =  3.0_rp * a * s2  * t41 * z1
          deriv(3,14) =  3.0_rp * a * s2  * t4  * z11
          
          deriv(1,15) =  3.0_rp * a * s11 * t3  * z1
          deriv(2,15) =  3.0_rp * a * s1  * t31 * z1
          deriv(3,15) =  3.0_rp * a * s1  * t3  * z11
          
          deriv(1,16) =  3.0_rp * a * s11 * t2  * z1
          deriv(2,16) =  3.0_rp * a * s1  * t21 * z1
          deriv(3,16) =  3.0_rp * a * s1  * t2  * z11
          
          deriv(1,17) =  3.0_rp * a * s11 * t1  * z2
          deriv(2,17) =  3.0_rp * a * s1  * t11 * z2
          deriv(3,17) =  3.0_rp * a * s1  * t1  * z21
          
          deriv(1,18) =  3.0_rp * a * s41 * t1  * z2
          deriv(2,18) =  3.0_rp * a * s4  * t11 * z2
          deriv(3,18) =  3.0_rp * a * s4  * t1  * z21
          
          deriv(1,19) =  3.0_rp * a * s41 * t4  * z2
          deriv(2,19) =  3.0_rp * a * s4  * t41 * z2
          deriv(3,19) =  3.0_rp * a * s4  * t4  * z21
          
          deriv(1,20) =  3.0_rp * a * s11 * t4  * z2
          deriv(2,20) =  3.0_rp * a * s1  * t41 * z2
          deriv(3,20) =  3.0_rp * a * s1  * t4  * z21
          
          deriv(1,21) =  3.0_rp * a * s11 * t1  * z3
          deriv(2,21) =  3.0_rp * a * s1  * t11 * z3
          deriv(3,21) =  3.0_rp * a * s1  * t1  * z31
          
          deriv(1,22) =  3.0_rp * a * s41 * t1  * z3
          deriv(2,22) =  3.0_rp * a * s4  * t11 * z3
          deriv(3,22) =  3.0_rp * a * s4  * t1  * z31
          
          deriv(1,23) =  3.0_rp * a * s41 * t4  * z3
          deriv(2,23) =  3.0_rp * a * s4  * t41 * z3
          deriv(3,23) =  3.0_rp * a * s4  * t4  * z31
          
          deriv(1,24) =  3.0_rp * a * s11 * t4  * z3
          deriv(2,24) =  3.0_rp * a * s1  * t41 * z3
          deriv(3,24) =  3.0_rp * a * s1  * t4  * z31
          
          deriv(1,25) =  3.0_rp * a * s21 * t1  * z4
          deriv(2,25) =  3.0_rp * a * s2  * t11 * z4
          deriv(3,25) =  3.0_rp * a * s2  * t1  * z41
          
          deriv(1,26) =  3.0_rp * a * s31 * t1  * z4
          deriv(2,26) =  3.0_rp * a * s3  * t11 * z4
          deriv(3,26) =  3.0_rp * a * s3  * t1  * z41
          
          deriv(1,27) =  3.0_rp * a * s41 * t2  * z4
          deriv(2,27) =  3.0_rp * a * s4  * t21 * z4
          deriv(3,27) =  3.0_rp * a * s4  * t2  * z41
          
          deriv(1,28) =  3.0_rp * a * s41 * t3  * z4
          deriv(2,28) =  3.0_rp * a * s4  * t31 * z4
          deriv(3,28) =  3.0_rp * a * s4  * t3  * z41
          
          deriv(1,29) =  3.0_rp * a * s31 * t4  * z4
          deriv(2,29) =  3.0_rp * a * s3  * t41 * z4
          deriv(3,29) =  3.0_rp * a * s3  * t4  * z41
          
          deriv(1,30) =  3.0_rp * a * s21 * t4  * z4
          deriv(2,30) =  3.0_rp * a * s2  * t41 * z4
          deriv(3,30) =  3.0_rp * a * s2  * t4  * z41
          
          deriv(1,31) =  3.0_rp * a * s11 * t3  * z4
          deriv(2,31) =  3.0_rp * a * s1  * t31 * z4
          deriv(3,31) =  3.0_rp * a * s1  * t3  * z41
          
          deriv(1,32) =  3.0_rp * a * s11 * t2  * z4
          deriv(2,32) =  3.0_rp * a * s1  * t21 * z4
          deriv(3,32) =  3.0_rp * a * s1  * t2  * z41
          
          deriv(1,33) = -9.0_rp * a * s21 * t2  * z1
          deriv(2,33) = -9.0_rp * a * s2  * t21 * z1
          deriv(3,33) = -9.0_rp * a * s2  * t2  * z11
          
          deriv(1,34) = -9.0_rp * a * s31 * t2  * z1
          deriv(2,34) = -9.0_rp * a * s3  * t21 * z1
          deriv(3,34) = -9.0_rp * a * s3  * t2  * z11
          
          deriv(1,35) = -9.0_rp * a * s31 * t3  * z1
          deriv(2,35) = -9.0_rp * a * s3  * t31 * z1
          deriv(3,35) = -9.0_rp * a * s3  * t3  * z11
          
          deriv(1,36) = -9.0_rp * a * s21 * t3  * z1
          deriv(2,36) = -9.0_rp * a * s2  * t31 * z1
          deriv(3,36) = -9.0_rp * a * s2  * t3  * z11
          
          deriv(1,37) = -9.0_rp * a * s21 * t1  * z2
          deriv(2,37) = -9.0_rp * a * s2  * t11 * z2
          deriv(3,37) = -9.0_rp * a * s2  * t1  * z21
          
          deriv(1,38) = -9.0_rp * a * s31 * t1  * z2
          deriv(2,38) = -9.0_rp * a * s3  * t11 * z2
          deriv(3,38) = -9.0_rp * a * s3  * t1  * z21
          
          deriv(1,39) = -9.0_rp * a * s41 * t2  * z2
          deriv(2,39) = -9.0_rp * a * s4  * t21 * z2
          deriv(3,39) = -9.0_rp * a * s4  * t2  * z21
          
          deriv(1,40) = -9.0_rp * a * s41 * t3  * z2
          deriv(2,40) = -9.0_rp * a * s4  * t31 * z2
          deriv(3,40) = -9.0_rp * a * s4  * t3  * z21
          
          deriv(1,41) = -9.0_rp * a * s31 * t4  * z2
          deriv(2,41) = -9.0_rp * a * s3  * t41 * z2
          deriv(3,41) = -9.0_rp * a * s3  * t4  * z21
          
          deriv(1,42) = -9.0_rp * a * s21 * t4  * z2
          deriv(2,42) = -9.0_rp * a * s2  * t41 * z2
          deriv(3,42) = -9.0_rp * a * s2  * t4  * z21
          
          deriv(1,43) = -9.0_rp * a * s11 * t3  * z2
          deriv(2,43) = -9.0_rp * a * s1  * t31 * z2
          deriv(3,43) = -9.0_rp * a * s1  * t3  * z21
          
          deriv(1,44) = -9.0_rp * a * s11 * t2  * z2
          deriv(2,44) = -9.0_rp * a * s1  * t21 * z2
          deriv(3,44) = -9.0_rp * a * s1  * t2  * z21
          
          deriv(1,45) = -9.0_rp * a * s21 * t1  * z3
          deriv(2,45) = -9.0_rp * a * s2  * t11 * z3
          deriv(3,45) = -9.0_rp * a * s2  * t1  * z31
          
          deriv(1,46) = -9.0_rp * a * s31 * t1  * z3
          deriv(2,46) = -9.0_rp * a * s3  * t11 * z3
          deriv(3,46) = -9.0_rp * a * s3  * t1  * z31
          
          deriv(1,47) = -9.0_rp * a * s41 * t2  * z3
          deriv(2,47) = -9.0_rp * a * s4  * t21 * z3
          deriv(3,47) = -9.0_rp * a * s4  * t2  * z31
          
          deriv(1,48) = -9.0_rp * a * s41 * t3  * z3
          deriv(2,48) = -9.0_rp * a * s4  * t31 * z3
          deriv(3,48) = -9.0_rp * a * s4  * t3  * z31
          
          deriv(1,49) = -9.0_rp * a * s31 * t4  * z3
          deriv(2,49) = -9.0_rp * a * s3  * t41 * z3
          deriv(3,49) = -9.0_rp * a * s3  * t4  * z31
          
          deriv(1,50) = -9.0_rp * a * s21 * t4  * z3
          deriv(2,50) = -9.0_rp * a * s2  * t41 * z3
          deriv(3,50) = -9.0_rp * a * s2  * t4  * z31
          
          deriv(1,51) = -9.0_rp * a * s11 * t3  * z3
          deriv(2,51) = -9.0_rp * a * s1  * t31 * z3
          deriv(3,51) = -9.0_rp * a * s1  * t3  * z31
          
          deriv(1,52) = -9.0_rp * a * s11 * t2  * z3
          deriv(2,52) = -9.0_rp * a * s1  * t21 * z3
          deriv(3,52) = -9.0_rp * a * s1  * t2  * z31
          
          deriv(1,53) = -9.0_rp * a * s21 * t2  * z4
          deriv(2,53) = -9.0_rp * a * s2  * t21 * z4
          deriv(3,53) = -9.0_rp * a * s2  * t2  * z41
          
          deriv(1,54) = -9.0_rp * a * s31 * t2  * z4
          deriv(2,54) = -9.0_rp * a * s3  * t21 * z4
          deriv(3,54) = -9.0_rp * a * s3  * t2  * z41
          
          deriv(1,55) = -9.0_rp * a * s31 * t3  * z4
          deriv(2,55) = -9.0_rp * a * s3  * t31 * z4
          deriv(3,55) = -9.0_rp * a * s3  * t3  * z41
          
          deriv(1,56) = -9.0_rp * a * s21 * t3  * z4
          deriv(2,56) = -9.0_rp * a * s2  * t31 * z4
          deriv(3,56) = -9.0_rp * a * s2  * t3  * z41
          
          deriv(1,57) = 27.0_rp * a * s21 * t2  * z2
          deriv(2,57) = 27.0_rp * a * s2  * t21 * z2
          deriv(3,57) = 27.0_rp * a * s2  * t2  * z21
          
          deriv(1,58) = 27.0_rp * a * s31 * t2  * z2
          deriv(2,58) = 27.0_rp * a * s3  * t21 * z2
          deriv(3,58) = 27.0_rp * a * s3  * t2  * z21
          
          deriv(1,59) = 27.0_rp * a * s31 * t3  * z2
          deriv(2,59) = 27.0_rp * a * s3  * t31 * z2
          deriv(3,59) = 27.0_rp * a * s3  * t3  * z21
          
          deriv(1,60) = 27.0_rp * a * s21 * t3  * z2
          deriv(2,60) = 27.0_rp * a * s2  * t31 * z2
          deriv(3,60) = 27.0_rp * a * s2  * t3  * z21
          
          deriv(1,61) = 27.0_rp * a * s21 * t2  * z3
          deriv(2,61) = 27.0_rp * a * s2  * t21 * z3
          deriv(3,61) = 27.0_rp * a * s2  * t2  * z31
          
          deriv(1,62) = 27.0_rp * a * s31 * t2  * z3
          deriv(2,62) = 27.0_rp * a * s3  * t21 * z3
          deriv(3,62) = 27.0_rp * a * s3  * t2  * z31
          
          deriv(1,63) = 27.0_rp * a * s31 * t3  * z3
          deriv(2,63) = 27.0_rp * a * s3  * t31 * z3
          deriv(3,63) = 27.0_rp * a * s3  * t3  * z31
          
          deriv(1,64) = 27.0_rp * a * s21 * t3  * z3
          deriv(2,64) = 27.0_rp * a * s2  * t31 * z3
          deriv(3,64) = 27.0_rp * a * s2  * t3  * z31
       end if

       if( present(heslo) ) then

          s12 = -2.0_rp * ( (1.0_rp-s)        + (1.0_rp/3.0_rp-s) - (1.0_rp/3.0_rp+s))
          s22 = -2.0_rp * ( (1.0_rp-s)        + (1.0_rp/3.0_rp-s) - (1.0_rp+s))
          s32 =  2.0_rp * ( (1.0_rp-s)        - (1.0_rp/3.0_rp+s) - (1.0_rp+s))
          s42 =  2.0_rp * ( (1.0_rp/3.0_rp-s) - (1.0_rp/3.0_rp+s) - (1.0_rp+s))
          
          t12 = -2.0_rp * ( (1.0_rp-t)        + (1.0_rp/3.0_rp-t) - (1.0_rp/3.0_rp+t))
          t22 = -2.0_rp * ( (1.0_rp-t)        + (1.0_rp/3.0_rp-t) - (1.0_rp+t))
          t32 =  2.0_rp * ( (1.0_rp-t)        - (1.0_rp/3.0_rp+t) - (1.0_rp+t))
          t42 =  2.0_rp * ( (1.0_rp/3.0_rp-t) - (1.0_rp/3.0_rp+t) - (1.0_rp+t))
           
          z12 = -2.0_rp * ( (1.0_rp-z)        + (1.0_rp/3.0_rp-z) - (1.0_rp/3.0_rp+z))
          z22 = -2.0_rp * ( (1.0_rp-z)        + (1.0_rp/3.0_rp-z) - (1.0_rp+z))
          z32 =  2.0_rp * ( (1.0_rp-z)        - (1.0_rp/3.0_rp+z) - (1.0_rp+z))
          z42 =  2.0_rp * ( (1.0_rp/3.0_rp-z) - (1.0_rp/3.0_rp+z) - (1.0_rp+z))

          heslo(1, 1) =  -a * s12 * t1  * z1
          heslo(2, 1) =  -a * s1  * t12 * z1
          heslo(3, 1) =  -a * s1  * t1  * z12
          heslo(4, 1) =  -a * s11 * t11 * z1
          heslo(5, 1) =  -a * s11 * t1  * z11
          heslo(6, 1) =  -a * s1  * t11 * z11

          heslo(1, 2) =  -a * s42 * t1  * z1
          heslo(2, 2) =  -a * s4  * t12 * z1
          heslo(3, 2) =  -a * s4  * t1  * z12
          heslo(4, 2) =  -a * s41 * t11 * z1
          heslo(5, 2) =  -a * s41 * t1  * z11
          heslo(6, 2) =  -a * s4  * t11 * z11

          heslo(1, 3) =  -a * s42 * t4  * z1
          heslo(2, 3) =  -a * s4  * t42 * z1
          heslo(3, 3) =  -a * s4  * t4  * z12
          heslo(4, 3) =  -a * s41 * t41 * z1
          heslo(5, 3) =  -a * s41 * t4  * z11
          heslo(6, 3) =  -a * s4  * t41 * z11

          heslo(1, 4) =  -a * s12 * t4  * z1
          heslo(2, 4) =  -a * s1  * t42 * z1
          heslo(3, 4) =  -a * s1  * t4  * z12
          heslo(4, 4) =  -a * s11 * t41 * z1
          heslo(5, 4) =  -a * s11 * t4  * z11
          heslo(6, 4) =  -a * s1  * t41 * z11

          heslo(1, 5) =  -a * s12 * t1  * z4
          heslo(2, 5) =  -a * s1  * t12 * z4
          heslo(3, 5) =  -a * s1  * t1  * z42
          heslo(4, 5) =  -a * s11 * t11 * z4
          heslo(5, 5) =  -a * s11 * t1  * z41
          heslo(6, 5) =  -a * s1  * t11 * z41

          heslo(1, 6) =  -a * s42 * t1  * z4
          heslo(2, 6) =  -a * s4  * t12 * z4
          heslo(3, 6) =  -a * s4  * t1  * z42
          heslo(4, 6) =  -a * s41 * t11 * z4
          heslo(5, 6) =  -a * s41 * t1  * z41
          heslo(6, 6) =  -a * s4  * t11 * z41

          heslo(1, 7) =  -a * s42 * t4  * z4
          heslo(2, 7) =  -a * s4  * t42 * z4
          heslo(3, 7) =  -a * s4  * t4  * z42
          heslo(4, 7) =  -a * s41 * t41 * z4
          heslo(5, 7) =  -a * s41 * t4  * z41
          heslo(6, 7) =  -a * s4  * t41 * z41

          heslo(1, 8) =  -a * s12 * t4  * z4
          heslo(2, 8) =  -a * s1  * t42 * z4
          heslo(3, 8) =  -a * s1  * t4  * z42
          heslo(4, 8) =  -a * s11 * t41 * z4
          heslo(5, 8) =  -a * s11 * t4  * z41
          heslo(6, 8) =  -a * s1  * t41 * z41

          heslo(1, 9) =  3.0_rp * a * s22 * t1  * z1
          heslo(2, 9) =  3.0_rp * a * s2  * t12 * z1
          heslo(3, 9) =  3.0_rp * a * s2  * t1  * z12
          heslo(4, 9) =  3.0_rp * a * s21 * t11 * z1
          heslo(5, 9) =  3.0_rp * a * s21 * t1  * z11
          heslo(6, 9) =  3.0_rp * a * s2  * t11 * z11

          heslo(1,10) =  3.0_rp * a * s32 * t1  * z1
          heslo(2,10) =  3.0_rp * a * s3  * t12 * z1
          heslo(3,10) =  3.0_rp * a * s3  * t1  * z12
          heslo(4,10) =  3.0_rp * a * s31 * t11 * z1
          heslo(5,10) =  3.0_rp * a * s31 * t1  * z11
          heslo(6,10) =  3.0_rp * a * s3  * t11 * z11

          heslo(1,11) =  3.0_rp * a * s42 * t2  * z1
          heslo(2,11) =  3.0_rp * a * s4  * t22 * z1
          heslo(3,11) =  3.0_rp * a * s4  * t2  * z12
          heslo(4,11) =  3.0_rp * a * s41 * t21 * z1
          heslo(5,11) =  3.0_rp * a * s41 * t2  * z11
          heslo(6,11) =  3.0_rp * a * s4  * t21 * z11

          heslo(1,12) =  3.0_rp * a * s42 * t3  * z1
          heslo(2,12) =  3.0_rp * a * s4  * t32 * z1
          heslo(3,12) =  3.0_rp * a * s4  * t3  * z12
          heslo(4,12) =  3.0_rp * a * s41 * t31 * z1
          heslo(5,12) =  3.0_rp * a * s41 * t3  * z11
          heslo(6,12) =  3.0_rp * a * s4  * t31 * z11

          heslo(1,13) =  3.0_rp * a * s32 * t4  * z1
          heslo(2,13) =  3.0_rp * a * s3  * t42 * z1
          heslo(3,13) =  3.0_rp * a * s3  * t4  * z12
          heslo(4,13) =  3.0_rp * a * s31 * t41 * z1
          heslo(5,13) =  3.0_rp * a * s31 * t4  * z11
          heslo(6,13) =  3.0_rp * a * s3  * t41 * z11

          heslo(1,14) =  3.0_rp * a * s22 * t4  * z1
          heslo(2,14) =  3.0_rp * a * s2  * t42 * z1
          heslo(3,14) =  3.0_rp * a * s2  * t4  * z12
          heslo(4,14) =  3.0_rp * a * s21 * t41 * z1
          heslo(5,14) =  3.0_rp * a * s21 * t4  * z11
          heslo(6,14) =  3.0_rp * a * s2  * t41 * z11

          heslo(1,15) =  3.0_rp * a * s12 * t3  * z1
          heslo(2,15) =  3.0_rp * a * s1  * t32 * z1
          heslo(3,15) =  3.0_rp * a * s1  * t3  * z12
          heslo(4,15) =  3.0_rp * a * s11 * t31 * z1
          heslo(5,15) =  3.0_rp * a * s11 * t3  * z11
          heslo(6,15) =  3.0_rp * a * s1  * t31 * z11

          heslo(1,16) =  3.0_rp * a * s12 * t2  * z1
          heslo(2,16) =  3.0_rp * a * s1  * t22 * z1
          heslo(3,16) =  3.0_rp * a * s1  * t2  * z12
          heslo(4,16) =  3.0_rp * a * s11 * t21 * z1
          heslo(5,16) =  3.0_rp * a * s11 * t2  * z11
          heslo(6,16) =  3.0_rp * a * s1  * t21 * z11

          heslo(1,17) =  3.0_rp * a * s12 * t1  * z2
          heslo(2,17) =  3.0_rp * a * s1  * t12 * z2
          heslo(3,17) =  3.0_rp * a * s1  * t1  * z22
          heslo(4,17) =  3.0_rp * a * s11 * t11 * z2
          heslo(5,17) =  3.0_rp * a * s11 * t1  * z21
          heslo(6,17) =  3.0_rp * a * s1  * t11 * z21

          heslo(1,18) =  3.0_rp * a * s42 * t1  * z2
          heslo(2,18) =  3.0_rp * a * s4  * t12 * z2
          heslo(3,18) =  3.0_rp * a * s4  * t1  * z22
          heslo(4,18) =  3.0_rp * a * s41 * t11 * z2
          heslo(5,18) =  3.0_rp * a * s41 * t1  * z21
          heslo(6,18) =  3.0_rp * a * s4  * t11 * z21

          heslo(1,19) =  3.0_rp * a * s42 * t4  * z2
          heslo(2,19) =  3.0_rp * a * s4  * t42 * z2
          heslo(3,19) =  3.0_rp * a * s4  * t4  * z22
          heslo(4,19) =  3.0_rp * a * s41 * t41 * z2
          heslo(5,19) =  3.0_rp * a * s41 * t4  * z21
          heslo(6,19) =  3.0_rp * a * s4  * t41 * z21

          heslo(1,20) =  3.0_rp * a * s12 * t4  * z2
          heslo(2,20) =  3.0_rp * a * s1  * t42 * z2
          heslo(3,20) =  3.0_rp * a * s1  * t4  * z22
          heslo(4,20) =  3.0_rp * a * s11 * t41 * z2
          heslo(5,20) =  3.0_rp * a * s11 * t4  * z21
          heslo(6,20) =  3.0_rp * a * s1  * t41 * z21

          heslo(1,21) =  3.0_rp * a * s12 * t1  * z3
          heslo(2,21) =  3.0_rp * a * s1  * t12 * z3
          heslo(3,21) =  3.0_rp * a * s1  * t1  * z32
          heslo(4,21) =  3.0_rp * a * s11 * t11 * z3
          heslo(5,21) =  3.0_rp * a * s11 * t1  * z31
          heslo(6,21) =  3.0_rp * a * s1  * t11 * z31

          heslo(1,22) =  3.0_rp * a * s42 * t1  * z3
          heslo(2,22) =  3.0_rp * a * s4  * t12 * z3
          heslo(3,22) =  3.0_rp * a * s4  * t1  * z32
          heslo(4,22) =  3.0_rp * a * s41 * t11 * z3
          heslo(5,22) =  3.0_rp * a * s41 * t1  * z31
          heslo(6,22) =  3.0_rp * a * s4  * t11 * z31

          heslo(1,23) =  3.0_rp * a * s42 * t4  * z3
          heslo(2,23) =  3.0_rp * a * s4  * t42 * z3
          heslo(3,23) =  3.0_rp * a * s4  * t4  * z32
          heslo(4,23) =  3.0_rp * a * s41 * t41 * z3
          heslo(5,23) =  3.0_rp * a * s41 * t4  * z31
          heslo(6,23) =  3.0_rp * a * s4  * t41 * z31

          heslo(1,24) =  3.0_rp * a * s12 * t4  * z3
          heslo(2,24) =  3.0_rp * a * s1  * t42 * z3
          heslo(3,24) =  3.0_rp * a * s1  * t4  * z32
          heslo(4,24) =  3.0_rp * a * s11 * t41 * z3
          heslo(5,24) =  3.0_rp * a * s11 * t4  * z31
          heslo(6,24) =  3.0_rp * a * s1  * t41 * z31

          heslo(1,25) =  3.0_rp * a * s22 * t1  * z4
          heslo(2,25) =  3.0_rp * a * s2  * t12 * z4
          heslo(3,25) =  3.0_rp * a * s2  * t1  * z42
          heslo(4,25) =  3.0_rp * a * s21 * t11 * z4
          heslo(5,25) =  3.0_rp * a * s21 * t1  * z41
          heslo(6,25) =  3.0_rp * a * s2  * t11 * z41

          heslo(1,26) =  3.0_rp * a * s32 * t1  * z4
          heslo(2,26) =  3.0_rp * a * s3  * t12 * z4
          heslo(3,26) =  3.0_rp * a * s3  * t1  * z42
          heslo(4,26) =  3.0_rp * a * s31 * t11 * z4
          heslo(5,26) =  3.0_rp * a * s31 * t1  * z41
          heslo(6,26) =  3.0_rp * a * s3  * t11 * z41

          heslo(1,27) =  3.0_rp * a * s42 * t2  * z4
          heslo(2,27) =  3.0_rp * a * s4  * t22 * z4
          heslo(3,27) =  3.0_rp * a * s4  * t2  * z42
          heslo(4,27) =  3.0_rp * a * s41 * t21 * z4
          heslo(5,27) =  3.0_rp * a * s41 * t2  * z41
          heslo(6,27) =  3.0_rp * a * s4  * t21 * z41

          heslo(1,28) =  3.0_rp * a * s42 * t3  * z4
          heslo(2,28) =  3.0_rp * a * s4  * t32 * z4
          heslo(3,28) =  3.0_rp * a * s4  * t3  * z42
          heslo(4,28) =  3.0_rp * a * s41 * t31 * z4
          heslo(5,28) =  3.0_rp * a * s41 * t3  * z41
          heslo(6,28) =  3.0_rp * a * s4  * t31 * z41

          heslo(1,29) =  3.0_rp * a * s32 * t4  * z4
          heslo(2,29) =  3.0_rp * a * s3  * t42 * z4
          heslo(3,29) =  3.0_rp * a * s3  * t4  * z42
          heslo(4,29) =  3.0_rp * a * s31 * t41 * z4
          heslo(5,29) =  3.0_rp * a * s31 * t4  * z41
          heslo(6,29) =  3.0_rp * a * s3  * t41 * z41

          heslo(1,30) =  3.0_rp * a * s22 * t4  * z4
          heslo(2,30) =  3.0_rp * a * s2  * t42 * z4
          heslo(3,30) =  3.0_rp * a * s2  * t4  * z42
          heslo(4,30) =  3.0_rp * a * s21 * t41 * z4
          heslo(5,30) =  3.0_rp * a * s21 * t4  * z41
          heslo(6,30) =  3.0_rp * a * s2  * t41 * z41

          heslo(1,31) =  3.0_rp * a * s12 * t3  * z4
          heslo(2,31) =  3.0_rp * a * s1  * t32 * z4
          heslo(3,31) =  3.0_rp * a * s1  * t3  * z42
          heslo(4,31) =  3.0_rp * a * s11 * t31 * z4
          heslo(5,31) =  3.0_rp * a * s11 * t3  * z41
          heslo(6,31) =  3.0_rp * a * s1  * t31 * z41

          heslo(1,32) =  3.0_rp * a * s12 * t2  * z4
          heslo(2,32) =  3.0_rp * a * s1  * t22 * z4
          heslo(3,32) =  3.0_rp * a * s1  * t2  * z42
          heslo(4,32) =  3.0_rp * a * s11 * t21 * z4
          heslo(5,32) =  3.0_rp * a * s11 * t2  * z41
          heslo(6,32) =  3.0_rp * a * s1  * t21 * z41

          heslo(1,33) = -9.0_rp * a * s22 * t2  * z1
          heslo(2,33) = -9.0_rp * a * s2  * t22 * z1
          heslo(3,33) = -9.0_rp * a * s2  * t2  * z12
          heslo(4,33) = -9.0_rp * a * s21 * t21 * z1
          heslo(5,33) = -9.0_rp * a * s21 * t2  * z11
          heslo(6,33) = -9.0_rp * a * s2  * t21 * z11

          heslo(1,34) = -9.0_rp * a * s32 * t2  * z1
          heslo(2,34) = -9.0_rp * a * s3  * t22 * z1
          heslo(3,34) = -9.0_rp * a * s3  * t2  * z12
          heslo(4,34) = -9.0_rp * a * s31 * t21 * z1
          heslo(5,34) = -9.0_rp * a * s31 * t2  * z11
          heslo(6,34) = -9.0_rp * a * s3  * t21 * z11

          heslo(1,35) = -9.0_rp * a * s32 * t3  * z1
          heslo(2,35) = -9.0_rp * a * s3  * t32 * z1
          heslo(3,35) = -9.0_rp * a * s3  * t3  * z12
          heslo(4,35) = -9.0_rp * a * s31 * t31 * z1
          heslo(5,35) = -9.0_rp * a * s31 * t3  * z11
          heslo(6,35) = -9.0_rp * a * s3  * t31 * z11

          heslo(1,36) = -9.0_rp * a * s22 * t3  * z1
          heslo(2,36) = -9.0_rp * a * s2  * t32 * z1
          heslo(3,36) = -9.0_rp * a * s2  * t3  * z12
          heslo(4,36) = -9.0_rp * a * s21 * t31 * z1
          heslo(5,36) = -9.0_rp * a * s21 * t3  * z11
          heslo(6,36) = -9.0_rp * a * s2  * t31 * z11

          heslo(1,37) = -9.0_rp * a * s22 * t1  * z2
          heslo(2,37) = -9.0_rp * a * s2  * t12 * z2
          heslo(3,37) = -9.0_rp * a * s2  * t1  * z22
          heslo(4,37) = -9.0_rp * a * s21 * t11 * z2
          heslo(5,37) = -9.0_rp * a * s21 * t1  * z21
          heslo(6,37) = -9.0_rp * a * s2  * t11 * z21

          heslo(1,38) = -9.0_rp * a * s32 * t1  * z2
          heslo(2,38) = -9.0_rp * a * s3  * t12 * z2
          heslo(3,38) = -9.0_rp * a * s3  * t1  * z22
          heslo(4,38) = -9.0_rp * a * s31 * t11 * z2
          heslo(5,38) = -9.0_rp * a * s31 * t1  * z21
          heslo(6,38) = -9.0_rp * a * s3  * t11 * z21

          heslo(1,39) = -9.0_rp * a * s42 * t2  * z2
          heslo(2,39) = -9.0_rp * a * s4  * t22 * z2
          heslo(3,39) = -9.0_rp * a * s4  * t2  * z22
          heslo(4,39) = -9.0_rp * a * s41 * t21 * z2
          heslo(5,39) = -9.0_rp * a * s41 * t2  * z21
          heslo(6,39) = -9.0_rp * a * s4  * t21 * z21

          heslo(1,40) = -9.0_rp * a * s42 * t3  * z2
          heslo(2,40) = -9.0_rp * a * s4  * t32 * z2
          heslo(3,40) = -9.0_rp * a * s4  * t3  * z22
          heslo(4,40) = -9.0_rp * a * s41 * t31 * z2
          heslo(5,40) = -9.0_rp * a * s41 * t3  * z21
          heslo(6,40) = -9.0_rp * a * s4  * t31 * z21

          heslo(1,41) = -9.0_rp * a * s32 * t4  * z2
          heslo(2,41) = -9.0_rp * a * s3  * t42 * z2
          heslo(3,41) = -9.0_rp * a * s3  * t4  * z22
          heslo(4,41) = -9.0_rp * a * s31 * t41 * z2
          heslo(5,41) = -9.0_rp * a * s31 * t4  * z21
          heslo(6,41) = -9.0_rp * a * s3  * t41 * z21

          heslo(1,42) = -9.0_rp * a * s22 * t4  * z2
          heslo(2,42) = -9.0_rp * a * s2  * t42 * z2
          heslo(3,42) = -9.0_rp * a * s2  * t4  * z22
          heslo(4,42) = -9.0_rp * a * s21 * t41 * z2
          heslo(5,42) = -9.0_rp * a * s21 * t4  * z21
          heslo(6,42) = -9.0_rp * a * s2  * t41 * z21

          heslo(1,43) = -9.0_rp * a * s12 * t3  * z2
          heslo(2,43) = -9.0_rp * a * s1  * t32 * z2
          heslo(3,43) = -9.0_rp * a * s1  * t3  * z22
          heslo(4,43) = -9.0_rp * a * s11 * t31 * z2
          heslo(5,43) = -9.0_rp * a * s11 * t3  * z21
          heslo(6,43) = -9.0_rp * a * s1  * t31 * z21

          heslo(1,44) = -9.0_rp * a * s12 * t2  * z2
          heslo(2,44) = -9.0_rp * a * s1  * t22 * z2
          heslo(3,44) = -9.0_rp * a * s1  * t2  * z22
          heslo(4,44) = -9.0_rp * a * s11 * t21 * z2
          heslo(5,44) = -9.0_rp * a * s11 * t2  * z21
          heslo(6,44) = -9.0_rp * a * s1  * t21 * z21

          heslo(1,45) = -9.0_rp * a * s22 * t1  * z3
          heslo(2,45) = -9.0_rp * a * s2  * t12 * z3
          heslo(3,45) = -9.0_rp * a * s2  * t1  * z32
          heslo(4,45) = -9.0_rp * a * s21 * t11 * z3
          heslo(5,45) = -9.0_rp * a * s21 * t1  * z31
          heslo(6,45) = -9.0_rp * a * s2  * t11 * z31

          heslo(1,46) = -9.0_rp * a * s32 * t1  * z3
          heslo(2,46) = -9.0_rp * a * s3  * t12 * z3
          heslo(3,46) = -9.0_rp * a * s3  * t1  * z32
          heslo(4,46) = -9.0_rp * a * s31 * t11 * z3
          heslo(5,46) = -9.0_rp * a * s31 * t1  * z31
          heslo(6,46) = -9.0_rp * a * s3  * t11 * z31

          heslo(1,47) = -9.0_rp * a * s42 * t2  * z3
          heslo(2,47) = -9.0_rp * a * s4  * t22 * z3
          heslo(3,47) = -9.0_rp * a * s4  * t2  * z32
          heslo(4,47) = -9.0_rp * a * s41 * t21 * z3
          heslo(5,47) = -9.0_rp * a * s41 * t2  * z31
          heslo(6,47) = -9.0_rp * a * s4  * t21 * z31

          heslo(1,48) = -9.0_rp * a * s42 * t3  * z3
          heslo(2,48) = -9.0_rp * a * s4  * t32 * z3
          heslo(3,48) = -9.0_rp * a * s4  * t3  * z32
          heslo(4,48) = -9.0_rp * a * s41 * t31 * z3
          heslo(5,48) = -9.0_rp * a * s41 * t3  * z31
          heslo(6,48) = -9.0_rp * a * s4  * t31 * z31

          heslo(1,49) = -9.0_rp * a * s32 * t4  * z3
          heslo(2,49) = -9.0_rp * a * s3  * t42 * z3
          heslo(3,49) = -9.0_rp * a * s3  * t4  * z32
          heslo(4,49) = -9.0_rp * a * s31 * t41 * z3
          heslo(5,49) = -9.0_rp * a * s31 * t4  * z31
          heslo(6,49) = -9.0_rp * a * s3  * t41 * z31

          heslo(1,50) = -9.0_rp * a * s22 * t4  * z3
          heslo(2,50) = -9.0_rp * a * s2  * t42 * z3
          heslo(3,50) = -9.0_rp * a * s2  * t4  * z32
          heslo(4,50) = -9.0_rp * a * s21 * t41 * z3
          heslo(5,50) = -9.0_rp * a * s21 * t4  * z31
          heslo(6,50) = -9.0_rp * a * s2  * t41 * z31

          heslo(1,51) = -9.0_rp * a * s12 * t3  * z3
          heslo(2,51) = -9.0_rp * a * s1  * t32 * z3
          heslo(3,51) = -9.0_rp * a * s1  * t3  * z32
          heslo(4,51) = -9.0_rp * a * s11 * t31 * z3
          heslo(5,51) = -9.0_rp * a * s11 * t3  * z31
          heslo(6,51) = -9.0_rp * a * s1  * t31 * z31

          heslo(1,52) = -9.0_rp * a * s12 * t2  * z3
          heslo(2,52) = -9.0_rp * a * s1  * t22 * z3
          heslo(3,52) = -9.0_rp * a * s1  * t2  * z32
          heslo(4,52) = -9.0_rp * a * s11 * t21 * z3
          heslo(5,52) = -9.0_rp * a * s11 * t2  * z31
          heslo(6,52) = -9.0_rp * a * s1  * t21 * z31

          heslo(1,53) = -9.0_rp * a * s22 * t2  * z4
          heslo(2,53) = -9.0_rp * a * s2  * t22 * z4
          heslo(3,53) = -9.0_rp * a * s2  * t2  * z42
          heslo(4,53) = -9.0_rp * a * s21 * t21 * z4
          heslo(5,53) = -9.0_rp * a * s21 * t2  * z41
          heslo(6,53) = -9.0_rp * a * s2  * t21 * z41

          heslo(1,54) = -9.0_rp * a * s32 * t2  * z4
          heslo(2,54) = -9.0_rp * a * s3  * t22 * z4
          heslo(3,54) = -9.0_rp * a * s3  * t2  * z42
          heslo(4,54) = -9.0_rp * a * s31 * t21 * z4
          heslo(5,54) = -9.0_rp * a * s31 * t2  * z41
          heslo(6,54) = -9.0_rp * a * s3  * t21 * z41

          heslo(1,55) = -9.0_rp * a * s32 * t3  * z4
          heslo(2,55) = -9.0_rp * a * s3  * t32 * z4
          heslo(3,55) = -9.0_rp * a * s3  * t3  * z42
          heslo(4,55) = -9.0_rp * a * s31 * t31 * z4
          heslo(5,55) = -9.0_rp * a * s31 * t3  * z41
          heslo(6,55) = -9.0_rp * a * s3  * t31 * z41

          heslo(1,56) = -9.0_rp * a * s22 * t3  * z4
          heslo(2,56) = -9.0_rp * a * s2  * t32 * z4
          heslo(3,56) = -9.0_rp * a * s2  * t3  * z42
          heslo(4,56) = -9.0_rp * a * s21 * t31 * z4
          heslo(5,56) = -9.0_rp * a * s21 * t3  * z41
          heslo(6,56) = -9.0_rp * a * s2  * t31 * z41

          heslo(1,57) = 27.0_rp * a * s22 * t2  * z2
          heslo(2,57) = 27.0_rp * a * s2  * t22 * z2
          heslo(3,57) = 27.0_rp * a * s2  * t2  * z22
          heslo(4,57) = 27.0_rp * a * s21 * t21 * z2
          heslo(5,57) = 27.0_rp * a * s21 * t2  * z21
          heslo(6,57) = 27.0_rp * a * s2  * t21 * z21

          heslo(1,58) = 27.0_rp * a * s32 * t2  * z2
          heslo(2,58) = 27.0_rp * a * s3  * t22 * z2
          heslo(3,58) = 27.0_rp * a * s3  * t2  * z22
          heslo(4,58) = 27.0_rp * a * s31 * t21 * z2
          heslo(5,58) = 27.0_rp * a * s31 * t2  * z21
          heslo(6,58) = 27.0_rp * a * s3  * t21 * z21

          heslo(1,59) = 27.0_rp * a * s32 * t3  * z2
          heslo(2,59) = 27.0_rp * a * s3  * t32 * z2
          heslo(3,59) = 27.0_rp * a * s3  * t3  * z22
          heslo(4,59) = 27.0_rp * a * s31 * t31 * z2
          heslo(5,59) = 27.0_rp * a * s31 * t3  * z21
          heslo(6,59) = 27.0_rp * a * s3  * t31 * z21

          heslo(1,60) = 27.0_rp * a * s22 * t3  * z2
          heslo(2,60) = 27.0_rp * a * s2  * t32 * z2
          heslo(3,60) = 27.0_rp * a * s2  * t3  * z22
          heslo(4,60) = 27.0_rp * a * s21 * t31 * z2
          heslo(5,60) = 27.0_rp * a * s21 * t3  * z21
          heslo(6,60) = 27.0_rp * a * s2  * t31 * z21

          heslo(1,61) = 27.0_rp * a * s22 * t2  * z3
          heslo(2,61) = 27.0_rp * a * s2  * t22 * z3
          heslo(3,61) = 27.0_rp * a * s2  * t2  * z32
          heslo(4,61) = 27.0_rp * a * s21 * t21 * z3
          heslo(5,61) = 27.0_rp * a * s21 * t2  * z31
          heslo(6,61) = 27.0_rp * a * s2  * t21 * z31

          heslo(1,62) = 27.0_rp * a * s32 * t2  * z3
          heslo(2,62) = 27.0_rp * a * s3  * t22 * z3
          heslo(3,62) = 27.0_rp * a * s3  * t2  * z32
          heslo(4,62) = 27.0_rp * a * s31 * t21 * z3
          heslo(5,62) = 27.0_rp * a * s31 * t2  * z31
          heslo(6,62) = 27.0_rp * a * s3  * t21 * z31

          heslo(1,63) = 27.0_rp * a * s32 * t3  * z3
          heslo(2,63) = 27.0_rp * a * s3  * t32 * z3
          heslo(3,63) = 27.0_rp * a * s3  * t3  * z32
          heslo(4,63) = 27.0_rp * a * s31 * t31 * z3
          heslo(5,63) = 27.0_rp * a * s31 * t3  * z31
          heslo(6,63) = 27.0_rp * a * s3  * t31 * z31

          heslo(1,64) = 27.0_rp * a * s22 * t3  * z3
          heslo(2,64) = 27.0_rp * a * s2  * t32 * z3
          heslo(3,64) = 27.0_rp * a * s2  * t3  * z32
          heslo(4,64) = 27.0_rp * a * s21 * t31 * z3
          heslo(5,64) = 27.0_rp * a * s21 * t3  * z31
          heslo(6,64) = 27.0_rp * a * s2  * t31 * z31
       end if

    else if( nnode == 3 ) then
       !
       ! SHELL
       !
       shapf(1) = 1.0_rp-s-t
       shapf(2) = s
       shapf(3) = t                                        
       if( present(deriv) ) then                             !  3
          deriv(1,1) =-1.0_rp                                !
          deriv(1,2) = 1.0_rp                                !
          deriv(1,3) = 0.0_rp                                !
          deriv(2,1) =-1.0_rp                                !  1       2
          deriv(2,2) = 0.0_rp
          deriv(2,3) = 1.0_rp
       end if
       
    else

       ierr = 1
       
    end if

    if( present(ierro) ) ierro = ierr
    
  end subroutine shape3
  
end module mod_lagrange
!> @}
