!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine pts_rungk45(                                       &
     ndime,kfl_drafo,kfl_extfo,grafo,buofo,gravi,veloc_fluid, &
     visco_fluid,densi_fluid,v,x,densi_parti,                 &
     diame_parti,spher_parti,dista,t,dt,vf,xf                 )
  !
  ! Returns final (position, velocity) tuple after time dt has passed
  !
  use def_kintyp, only : ip,rp
  implicit none
  integer(ip), intent(in)  :: ndime
  integer(ip), intent(in)  :: kfl_drafo
  integer(ip), intent(in)  :: kfl_extfo
  real(rp),    intent(in)  :: grafo
  real(rp),    intent(in)  :: buofo
  real(rp),    intent(in)  :: gravi(*)
  real(rp),    intent(in)  :: veloc_fluid(*)
  real(rp),    intent(in)  :: visco_fluid
  real(rp),    intent(in)  :: densi_fluid
  real(rp),    intent(in)  :: v(*)             !< Initial velocity
  real(rp),    intent(in)  :: x(*)             !< Initial position
  real(rp),    intent(in)  :: densi_parti
  real(rp),    intent(in)  :: diame_parti
  real(rp),    intent(in)  :: spher_parti
  real(rp),    intent(in)  :: dista               !< Distance
  real(rp),    intent(inout) :: t                !< Current time
  real(rp),    intent(inout) :: dt               !< Time step
  real(rp),    intent(out) :: vf(*)            !< Final velocity
  real(rp),    intent(out) :: xf(*)            !< Final position

  integer(ip)              :: idime
  real(rp)                 :: t1,x1(3),v1(3),k1(3)
  real(rp)                 :: t2,x2(3),v2(3),k2(3)
  real(rp)                 :: t3,x3(3),v3(3),k3(3)
  real(rp)                 :: t4,x4(3),v4(3),k4(3)
  real(rp)                 :: t5,x5(3),v5(3),k5(3)
  real(rp)                 :: t6,x6(3),v6(3),k6(3)
  real(rp)                 :: xf4(3),vf4(3)
  real(rp)                 :: xf5(3),vf5(3)
  real(rp)                 :: err,strec,dampi,facc,xdeno,toler

  ! def rkf45(t, dt, y, f, tolerance=1e-5):
  ! t = Current time.
  ! dt = Timestep.
  ! y = Initial value.
  ! f = Derivative function y' = f(t, y).
  ! tolerance = Error tolerance.
  ! Returns a (tf, dtf, yf) tuple after time dt has passed, where:
  ! tf = Final time.
  ! yf = Final value.
  ! dtf = Error-corrected timestep for next step."""
  ! Calculate the slopes at various points.
  ! Values taken from http://en.wikipedia.org/wiki/R....
  !

  toler = 1.0e-2_rp
  err   = 1.0_rp

  do while( err > toler )

     t1 = t
     t2 = t + (1.0_rp/4.0_rp)   * dt
     t3 = t + (3.0_rp/8.0_rp)   * dt
     t4 = t + (12.0_rp/13.0_rp) * dt
     t5 = t + dt
     t6 = t + (1.0_rp/2.0_rp)   * dt

     do idime = 1,ndime
        v1(idime) = v(idime)
        x1(idime) = x(idime)
        call pts_accele(                                               &
             ndime,kfl_drafo,kfl_extfo,grafo,buofo,gravi,veloc_fluid,  &
             visco_fluid,densi_fluid,v1,x1,densi_parti,                &
             diame_parti,spher_parti,dista,t,k1)
     end do

     do idime = 1,ndime
        v2(idime) = v(idime) + k1(idime) * (1.0_rp/4.0_rp)  * dt
        x2(idime) = x(idime) + v2(idime) * (1.0_rp/4.0_rp)  * dt
        call pts_accele(                                               &
             ndime,kfl_drafo,kfl_extfo,grafo,buofo,gravi,veloc_fluid,  &
             visco_fluid,densi_fluid,v2,x2,densi_parti,                &
             diame_parti,spher_parti,dista,t,k2)
     end do

     do idime = 1,ndime
        v3(idime) = v(idime) + k1(idime) * (3.0_rp/32.0_rp) * dt + k2(idime) * (9.0_rp/32.0_rp) * dt
        x3(idime) = x(idime) + v1(idime) * (3.0_rp/32.0_rp) * dt + v2(idime) * (9.0_rp/32.0_rp) * dt
        call pts_accele(                                               &
             ndime,kfl_drafo,kfl_extfo,grafo,buofo,gravi,veloc_fluid,  &
             visco_fluid,densi_fluid,v3,x3,densi_parti,                &
             diame_parti,spher_parti,dista,t,k3)
     end do

     do idime = 1,ndime
        v4(idime) = v(idime) + k1(idime) * (1932.0_rp/2197.0_rp) * dt + k2(idime) * (-7200.0_rp/2197.0_rp) * dt + k3(idime) * (7296.0_rp/2197.0_rp) * dt
        x4(idime) = x(idime) + v1(idime) * (1932.0_rp/2197.0_rp) * dt + v2(idime) * (-7200.0_rp/2197.0_rp) * dt + v3(idime) * (7296.0_rp/2197.0_rp) * dt
        call pts_accele(                                               &
             ndime,kfl_drafo,kfl_extfo,grafo,buofo,gravi,veloc_fluid,  &
             visco_fluid,densi_fluid,v4,x4,densi_parti,                &
             diame_parti,spher_parti,dista,t,k4)
     end do

     do idime = 1,ndime
        v5(idime) = v(idime) + k1(idime) * (439.0_rp/216.0_rp) * dt + k2(idime) * (-8.0_rp) * dt + k3(idime) * (3680.0_rp/513.0_rp) * dt + k4(idime) * (-845.0_rp/4104.0_rp) * dt
        x5(idime) = x(idime) + v1(idime) * (439.0_rp/216.0_rp) * dt + v2(idime) * (-8.0_rp) * dt + v3(idime) * (3680.0_rp/513.0_rp) * dt + v4(idime) * (-845.0_rp/4104.0_rp) * dt
        call pts_accele(                                               &
             ndime,kfl_drafo,kfl_extfo,grafo,buofo,gravi,veloc_fluid,  &
             visco_fluid,densi_fluid,v5,x5,densi_parti,                &
             diame_parti,spher_parti,dista,t,k5)
     end do

     do idime = 1,ndime
        v6(idime) = v(idime) + k1(idime) *(-8.0_rp/27.0_rp)*dt + k2(idime) * (2.0_rp) * dt + k3(idime) * (-3544.0_rp/2565.0_rp) * dt + k4(idime) * (1859.0_rp/4104.0_rp) * dt + k5(idime) * (-11.0_rp/40.0_rp) * dt
        x6(idime) = x(idime) + v1(idime) *(-8.0_rp/27.0_rp)*dt + v2(idime) * (2.0_rp) * dt + v3(idime) * (-3544.0_rp/2565.0_rp) * dt + v4(idime) * (1859.0_rp/4104.0_rp) * dt + v5(idime) * (-11.0_rp/40.0_rp) * dt
        call pts_accele(                                               &
             ndime,kfl_drafo,kfl_extfo,grafo,buofo,gravi,veloc_fluid,  &
             visco_fluid,densi_fluid,v6,x6,densi_parti,                &
             diame_parti,spher_parti,dista,t,k6)
     end do

     err   = 0.0_rp
     xdeno = 0.0_rp
     do idime = 1,ndime
        !
        ! 4th order approximation of the final value..0/
        !
        vf4(idime) = v(idime) + k1(idime)*(25.0_rp/216.0_rp)*dt + k3(idime)*(1408.0_rp/2565.0_rp)*dt + k4(idime)*(2197.0_rp/4104.0_rp)*dt + k5(idime)*(-1.0_rp/5.0_rp)*dt
        xf4(idime) = x(idime) + v1(idime)*(25.0_rp/216.0_rp)*dt + v3(idime)*(1408.0_rp/2565.0_rp)*dt + v4(idime)*(2197.0_rp/4104.0_rp)*dt + v5(idime)*(-1.0_rp/5.0_rp)*dt
        !
        ! 5th order approximation of the final value.
        !
        vf5(idime) = v(idime) + k1(idime)*(16.0_rp/135.0_rp)*dt + k3(idime)*(6656.0_rp/12825.0_rp)*dt + k4(idime)*(28561.0_rp/56430.0_rp)*dt + k5(idime)*(-9.0_rp/50.0_rp)*dt + k6(idime)*(2.0_rp/55.0_rp)*dt
        xf5(idime) = x(idime) + v1(idime)*(16.0_rp/135.0_rp)*dt + v3(idime)*(6656.0_rp/12825.0_rp)*dt + v4(idime)*(28561.0_rp/56430.0_rp)*dt + v5(idime)*(-9.0_rp/50.0_rp)*dt + v6(idime)*(2.0_rp/55.0_rp)*dt     
        !
        ! Timestep scaling factor. From http://math.fullerton.edu/math....
        !
        err   = err   + abs(vf5(idime)-vf4(idime))**2
        xdeno = xdeno + abs(vf5(idime)+vf4(idime))**2
        !
        !
        !
        vf(idime) = vf5(idime)
        xf(idime) = xf5(idime)
     end do

     if( xdeno > 0.0_rp ) then
        err = 2.0_rp * err / xdeno
     else
        err = 0.0_rp
     end if

     if( err > 0.0_rp ) then
        facc = toler / err
     else
        facc = 1.0e5_rp
     end if
     facc  = min(1.0e5_rp,facc)
     strec = 1.5_rp
     dampi = 1.2_rp
     if( err > toler ) t = t - dt
     dt    = dt * min(strec,log(1.0_rp + (dampi-1.0_rp)*facc) / log(dampi))
     dt    = max(dt,1.0e-12_rp)
     if( err > toler ) t = t + dt 
     !dt    = min(dt,dtime)
     write(400,*) err,dt,min(strec,log(1.0_rp + (dampi-1.0_rp)*facc) / log(dampi))

     ! s = 1 if err==0 else (tolerance*dt/(2*err))**(1/4.0_rp)

  end do

end subroutine pts_rungk45
