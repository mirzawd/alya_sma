!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



module mod_tauadr
  implicit none
  private
  public :: tauadr
   
contains 
subroutine tauadr(taust,staco,a1,k1,s1,h1,h2,tau,dtinv)
  !------------------------------------------------------------------------
  !****f* mathru/tauadr
  ! NAME 
  !    tauadr
  ! DESCRIPTION
  !    Compute the stabilization parameter tau for the
  !    advection-diffusion-reaction equation
  !    -k*Lapl(u) + a.grad(u) + s*u = Q 
  !    accordining to the possible following strategies:
  !    TAUST = 1 ... Codina
  !          = 2 ... Average of exact 1D equation with constant residual 
  !          = 3 ... Shakib
  !          = 4 ... Directional (not implemented)
  !          = 5 ... Codina with time step ( 1/tau = 1/tau_codina + rho/dt) 
  !          = 6 ... Time step  (tau = dt/rho)
  ! USES
  ! USED BY
  !    *_elmsgs
  !***
  !------------------------------------------------------------------------
  use def_kintyp, only     :  ip,rp
  implicit none
  integer(ip), intent(in)           :: taust
  real(rp),    intent(in)           :: k1,a1,s1,h1,h2,staco(3)
  real(rp),    intent(in),optional  :: dtinv !should be affected by density
  real(rp),    intent(out)          :: tau
  real(rp)                 :: Dah,Dainv,Peh,Kh,Ah,alpha,PehDah
  real(rp)                 :: freq1,freq2,freq3,k,a,s

  if(h1==0.0_rp.and.h2==0.0_rp) then
     tau=0.0_rp
     return
  end if

  a=staco(2)*a1 ! Advection
  k=staco(1)*k1 ! Diffusion 
  s=staco(3)*s1 ! Reaction

  select case(taust) ! tau strategy

  case(0)
     !
     ! No stabilization
     !
     tau = 0.0_rp

  case(1)
     !
     ! Codina
     !
     freq1 = 4.0_rp*k/(h2*h2)
     freq2 = 2.0_rp*a/h1
     freq3 = abs(s)
     tau   = freq1+freq2+freq3
     if(tau/=0.0_rp) tau=1.0_rp/tau

  case(2)
     !
     ! Average of exact 1D equation with constant residual 
     !
     if(a/=0.0_rp.and.k/=0.0_rp.and.s==0.0_rp) then          ! AD
        Peh = a*h1/(2.0_rp*k)
        if(Peh<1.0e-3_rp) then
           alpha = Peh/3.0_rp
           tau   = h1*h1/(12.0_rp*k)
        else if(Peh>1.0e3_rp) then
           alpha = 1.0_rp
           tau   = h1/(2.0_rp*a)
        else
           alpha = 1.0_rp/tanh(Peh)-1.0_rp/(Peh)
           tau   = h1/(2.0_rp*a)*alpha
        end if

     else if(a==0.0_rp.and.k/=0.0_rp.and.s==0.0_rp) then     ! D
        tau  = h2*h2/(12.0_rp*k)

     else if(a/=0.0_rp.and.k==0.0_rp.and.s/=0.0_rp) then     ! AR

        Dah  = s*h1/a
        if(Dah>1e3_rp) then
           tau = 1.0_rp/s
        else if(Dah<-1.0e3_rp) then
           tau = (2.0_rp*exp(-Dah)/Dah)/Dah
        else if(abs(Dah)<1.0e-3_rp) then
           tau = h1/(2.0_rp*a)
        else
           Dainv = 1.0_rp/Dah
           tau   = h1/(2.0_rp*a)*(2.0_rp*Dainv*&
                &  (1.0_rp+Dainv*(exp(-Dah)-1.0_rp)))    
        end if

     else if(a/=0.0_rp.and.k==0.0_rp.and.s==0.0_rp) then     ! A
        tau   = h1/(2.0_rp*a)

     else if(a==0.0_rp.and.k==0.0_rp.and.s/=0.0_rp) then     ! R
        tau   = 1.0_rp/s

     else if(a==0.0_rp.and.k/=0.0_rp.and.s/=0.0_rp) then     ! DR
        Ah    = sqrt(s/k)*h2
        tau   = 2.0_rp*h2*(1.0_rp-cosh(Ah*h2))&
             &  /(Ah*s*sinh(Ah*h2))
        Kh    = sqrt(2.0_rp*s*h2*h2/(2.0_rp*k))
        tau   = 1.0_rp/s*(1.0_rp+2.0_rp*(1.0_rp-cosh(Kh))&
             &  /(Kh*sinh(Kh)))

     else if(a/=0.0_rp.and.k/=0.0_rp.and.s/=0.0_rp) then     ! ADR
 
        Peh    = a*h1/(2.0_rp*k)
        Dah    = s*h1/a
        PehDah = s*h1*h1/(2.0_rp*k)
        
        if(Dah<1.0e-6_rp) then
           alpha = 1.0_rp/tanh(Peh)-1.0_rp/(Peh)
           tau   = h1/(2.0_rp*a)*alpha
        else if(Dah<1.0e-6_rp*Peh) then
           freq2 = ( exp(-Dah)+exp(-2.0_rp*Peh-Dah) ) / (1.0_rp-exp(-2.0_rp*(Peh+Dah))) 
           freq3 = (Peh+Dah)/PehDah
           tau   = 1.0_rp/s*( 1.0_rp + (Peh+Dah)/PehDah*(freq2 - 1.0_rp/tanh(Peh+Dah)) )
        else
           Ah    = sqrt(Peh*Peh+2.0_rp*PehDah)
           freq2 = ( exp(Peh-Ah)+exp(-Peh-Ah) ) / (1.0_rp-exp(-2.0_rp*Ah)) 
           freq3 = Ah/PehDah
           tau   = 1.0_rp/s*( 1.0_rp + Ah/PehDah*(freq2 - 1.0_rp/tanh(Ah)) )
        end if

        if(s>0.0_rp) tau=max(0.0_rp,tau)

     else 
        tau   = 0.0_rp

     end if

  case(3)
     !
     ! Shakib
     !
     tau = 9.0_rp*(4.0_rp*k/(h2*h2))**2+(2.0_rp*a/h1)**2+s*s
     if(tau>0.0_rp) tau=1.0_rp/sqrt(tau)

  case(4)
     !
     ! Directional
     !
  case(5)
     !
     ! Include Dt in Tau - the rest idem Codina 
     ! 
     if( present(dtinv) ) then
        freq1 = 4.0_rp*k/(h2*h2)
        freq2 = 2.0_rp*a/h1
        freq3 = abs(s)
        tau   = freq1+freq2+freq3+dtinv
        if(tau/=0.0_rp) tau=1.0_rp/tau
     else
        call runend('TAUADR:.not.present(dtinv)') 
     end if

  case(6)
     !
     ! Tau =DT 
     ! 
     if( present(dtinv) ) then
        tau=1.0_rp/dtinv
     else
        call runend('TAUADR:.not.present(dtinv)') 
     end if
  end select

end subroutine tauadr
end module mod_tauadr
