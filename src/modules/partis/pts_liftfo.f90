!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



function funReZeng(L,Re,d) result (x)
   use def_kintyp, only     :  ip,rp
   implicit none
   real(rp),  intent(in)  :: L
   real(rp),  intent(in)  :: Re
   real(rp),  intent(in)  :: d
   real(rp)               :: x
   real(rp)               :: f0, f1, Clt0, L2


   f0 = 1.0_rp + 0.329_rp*Re + 0.00485_rp*Re*Re
   f1 = -0.9_rp * tanh( 0.022_rp*Re )

   L2 = (L*Re)/d 

   if( L2 < 10_rp ) then

      Clt0 = ( 9.0_rp/8.0_rp + 5.78e-6 * L2 ) * exp( -0.292_rp*L2 )

   else

      Clt0 = 8.94_rp*L2**(-2.09_rp)

   end if


   x = f0 * Clt0 * L**f1

end function funReZeng

subroutine liftfo(imode,ur,wf,mu,rho,d,dista,Cls,Re,Res)

  !-----------------------------------------------------------------------
  !****f* liftfo/liftfo
  ! NAME
  !    liftfo
  ! DESCRIPTION
  ! 
  !  Lift coefficent for a particle immersed  into a viscous fluid
  !
  !
  ! USED BY
  !    pts_solite
  !***
  !-----------------------------------------------------------------------
 

   use def_kintyp, only     :  ip,rp
   implicit none
   integer(ip), intent(in)  :: imode
   real(rp),    intent(in)  :: ur          !< Relative velocity |u_fluid-u_particle|
   real(rp),    intent(in)  :: wf          !< Vorticity of velocity 
   real(rp),    intent(in)  :: mu         !< Fluid viscosity
   real(rp),    intent(in)  :: rho        !< Fluid density
   real(rp),    intent(in)  :: d          !< Particle diameter
   real(rp),    intent(in)  :: dista      !< Distance to the wall     
   real(rp),    intent(out) :: Cls        !< Lift coefficient
   real(rp),    intent(out) :: Re         !< Reynolds number
   real(rp),    intent(out) :: Res        !< Shear Reynolds number
   real(rp)                 :: alpha, beta, lambda, delta, Clsw, funcRe
   real(rp)                 :: zeror, g
   real(rp)                 :: funReZeng
   !
   ! Particle Reynolds number
   !
   zeror    = epsilon(1.0_rp)

   Re            =  rho * ur * d / mu + zeror
   Res           =  rho * d*d * wf / mu

   select case ( imode )
   case ( 1_ip )
      !
      ! Saffman Mei
      !
      if (Re == 0.0_rp) then
         beta   =  0.0_rp
      else
         beta   =  0.5_rp * Res/Re
      end if

      if(Re > 40.0_rp) then
         funcRe = 0.0524_rp * sqrt(beta*Re)
      else
         funcRe = (1.0_rp-0.3314_rp * sqrt(beta)) * exp(-Re/10.0_rp) + 0.3314_rp * sqrt(beta)
      end if

      Cls = 4.1126_rp * funcRe /sqrt(Res)

   case ( 2_ip )
      !
      ! Zeng, L., Najjar, F., Balachandar, S., Fischer, P., 2009. Forces on a finite-sized particle 
      ! located close to a wall in a linear shear flow. Phys. Fluids 21, 033302.
      !
      ! Is applicable in circumstances when a stationary particle is positioned in a wall-bounded linear shear flow for 1 < Res < 200 and
      ! even when the particle touches the wall. 
      !

      delta  = dista/d-0.5_rp      
      Clsw   = 3.663_rp / (Res*Res + 0.1173_rp)**0.22_rp 
      alpha  = - exp ( -0.3_rp + 0.025_rp*Res) 
      beta   = 0.8_rp + 0.01_rp * Res
      lambda = (1 - exp(-delta))**(2.5_rp)
      Cls    = Clsw * exp( -0.5_rp*delta  * (Res/250_rp)**(4.0_rp/3.0_rp) ) * (exp( alpha * delta**beta) - lambda) 

   case ( 3_ip )
      !
      ! Zeng, L., Najjar, F., Balachandar, S., Fischer, P., 2009. Forces on a finite-sized particle 
      ! located close to a wall in a linear shear flow. Phys. Fluids 21, 033302.
      !
      ! Is applicable with translation parallel to the nearby plane wall induced lift force, for 0 < Re < 100 and 0 < L* < 300. 
      !
      ! L* = LRe/d
      !
      ! where L is s the distance from the center of the particle to the nearby wall.
      !
      delta  = dista/d-0.5_rp      

      Clsw   = 0.313_rp + 0.812_rp*exp( -0.125_rp*Re**0.77_rp ) 
      g      = 3.0_rp * exp( -0.17_rp*Re**0.7_rp )
      Cls    = funReZeng(dista,Re,d) + (Clsw - funReZeng(0.5_rp,Re,d))*exp( -11_rp*(delta/g)**1.2_rp )

   case default
      !
      ! Others
      !
      call runend('LIFTFO: NON-EXISTING MODEL')

   end select


end subroutine liftfo
