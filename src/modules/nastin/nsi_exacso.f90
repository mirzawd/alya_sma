!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Nastin
!> @{
!> @file    nsi_exacso.f90
!> @author  Guillaume Houzeaux
!> @date    18/09/2012
!> @brief   Manufactured solutions
!> @details Manage manufactured solutions. They enable to check the correctness
!>          of the coding and to study the mesh convergence. If du/dt+L(u)=f is the
!>          Navier Stokes operator, this subroutine implements
!>          f = due/dt+L(ue), where ue is the manufactured (exact) solution
!>          programmed by the user.
!>          According to ITASK the subroutine performs the following:
!>          \verbatim
!>          ITASK=1 ... In order to calculate f.e.m. errors, compute:
!>                      EXVEL = u
!>                      EXVEG = grad(u)
!>                      EXPRE = p
!>                      EXPGR = grad(p) 
!>          ITASK=2 ... Compute force vector GPRHS=f to be applied as a load
!>          \endverbatim
!> @} 
!-----------------------------------------------------------------------
subroutine nsi_exacso(&
     itask,gpcod,gpden,gpvis,gppor,gpgvi,exvel,exveg,&
     expre,exprg,baloc,gprhs,gprhc,gprh2)
  use def_parame
  use def_master, only       :  cutim, zeror
  use def_domain, only       :  ndime
  use def_kermod, only       :  delta_dom,rough_dom
  use def_nastin, only       :  kfl_exacs_nsi,kfl_timei_nsi,kfl_advec_nsi,&
       &                        kfl_visco_nsi,corio_nsi,fvela_nsi,expar_nsi,&
       &                        fvins_nsi,kfl_convection_type_nsi,&
       &                        NSI_CONVECTION_SKEW,NSI_CONVECTION_EMAC,&
       &                        NSI_CONVECTION_CONSERVATIVE
  implicit none
  integer(ip), intent(in)    :: itask                  !< Task: element(1)/boundary(3) assembly, exact solution (2)
  real(rp),    intent(in)    :: gpden                  !< Gauss point density
  real(rp),    intent(in)    :: gpvis                  !< Gauss point viscosity
  real(rp),    intent(in)    :: gppor                  !< Gauss point porosity
  real(rp),    intent(in)    :: gpgvi(ndime)           !< Gauss point viscosity gradient
  real(rp),    intent(in)    :: gpcod(ndime)           !< Gauss point coordinate
  real(rp),    intent(out)   :: exvel(ndime)           !< Exact velocity
  real(rp),    intent(out)   :: exveg(ndime,ndime)     !< Exact velocity gradient
  real(rp),    intent(out)   :: expre                  !< Exact pressure
  real(rp),    intent(out)   :: exprg(ndime)           !< Exact pressure gradient
  real(rp),    intent(in)    :: baloc(ndime)           !< Boundary Gauss point normal
  real(rp),    intent(inout) :: gprhs(ndime)           !< Gauss point momentum RHS
  real(rp),    intent(inout) :: gprhc                  !< Gauss point continuity RHS 
  real(rp),    intent(inout) :: gprh2                  !< Gauss point continuity RHS using low mach
  real(rp)                   :: x,y,z,r,freq,diss,ft,dfdt
!  real(rp)                   :: n,omega,ri,ro,theta,vt
  real(rp)                   :: f,d1f,d2f,d3f,g,d1g,d2g,d3g,sx,sy,cx,cy
  real(rp)                   :: u,v,w,dudx,dudy,dudz,dvdx,dvdy,dvdz
  real(rp)                   :: dwdx,dwdy,dwdz,p,dpdx,dpdy,dpdz
  real(rp)                   :: d2udx,d2udy,d2udz,d2vdx,d2vdy,d2vdz
  real(rp)                   :: d2wdx,d2wdy,d2wdz,dudt,dvdt,dwdt
  real(rp)                   :: d2udxdy,d2vdxdy,d2wdxdy,divu,expt
  real(rp)                   :: d2udxdz,d2vdxdz,d2wdxdz,kap
  real(rp)                   :: d2udydz,d2vdydz,d2wdydz,ustar
  real(rp)                   :: mu,rho,sig,dmudx,dmudy,dmudz,wx,wy,wz
  !
  ! Initializations
  !
  x       = gpcod(1)
  y       = gpcod(2)
  z       = 0.0_rp
  rho     = gpden
  mu      = gpvis
  sig     = gppor
  dmudx   = gpgvi(1)
  dmudy   = gpgvi(2)
  dmudz   = 0.0_rp
  wx      = fvela_nsi(1)
  wy      = fvela_nsi(2)
  wz      = fvela_nsi(3)
  u       = 0.0_rp
  v       = 0.0_rp
  w       = 0.0_rp
  dudt    = 0.0_rp
  dvdt    = 0.0_rp
  dwdt    = 0.0_rp
  dudx    = 0.0_rp
  dudy    = 0.0_rp
  dudz    = 0.0_rp
  dvdx    = 0.0_rp
  dvdy    = 0.0_rp
  dvdz    = 0.0_rp
  dwdx    = 0.0_rp
  dwdy    = 0.0_rp
  dwdz    = 0.0_rp
  d2udx   = 0.0_rp
  d2udy   = 0.0_rp
  d2udz   = 0.0_rp
  d2vdx   = 0.0_rp
  d2vdy   = 0.0_rp
  d2vdz   = 0.0_rp
  d2wdx   = 0.0_rp
  d2wdy   = 0.0_rp
  d2wdz   = 0.0_rp
  p       = 0.0_rp
  dpdx    = 0.0_rp 
  dpdy    = 0.0_rp
  dpdz    = 0.0_rp
  d2udxdy = 0.0_rp 
  d2vdxdy = 0.0_rp
  d2wdxdy = 0.0_rp
  d2udxdz = 0.0_rp 
  d2vdxdz = 0.0_rp
  d2wdxdz = 0.0_rp
  d2udydz = 0.0_rp
  d2vdydz = 0.0_rp
  d2wdydz = 0.0_rp

  if(ndime==3 ) then
     z     = gpcod(ndime)
     dmudz = gpgvi(ndime)
  end if
  !
  ! Obtain unknowns and derivatives according to the exact solution.
  !
  if( kfl_exacs_nsi == 1 ) then
     !
     ! ux = a
     ! uy = b
     ! p  = 0
     ! 
     u = expar_nsi(1)
     v = expar_nsi(2)

  else if( kfl_exacs_nsi == 2 ) then
     !
     ! ux =  x
     ! uy = -y
     ! p  =  1
     !
     if(kfl_timei_nsi/=0 ) then
        diss = -expar_nsi(1)
        !       ft   = exp(diss*cutim)
        !       dfdt = diss*exp(diss*cutim)
        ft   = 1.0_rp+diss*cutim  ! First order series
        dfdt = diss
     else
        ft   = 1.0_rp
        dfdt = 0.0_rp
     end if
     p    =  1.0_rp
     u    =  x*ft
     v    = -y*ft
     dudt =  x*dfdt
     dvdt = -y*dfdt
     dudx =  ft
     dvdy = -ft

  else if( kfl_exacs_nsi == 3 ) then
     !
     ! ux = a*y
     ! uy = b*x
     ! p  = 4x+3y
     !
     u    = expar_nsi(1) * y
     v    = expar_nsi(2) * x
     dudy = expar_nsi(1)
     dvdx = expar_nsi(2)
     p    = 4.0_rp *x + 3.0_rp * y
     dpdx = 4.0_rp
     dpdy = 3.0_rp

  else if( kfl_exacs_nsi == 4 ) then
     !
     ! ux = y**n
     ! uy = 0
     ! p  = y
     !     
     u     =   y**expar_nsi(1)
     v     =   0.0_rp
     dudx  =   0.0_rp
     dudy  =   expar_nsi(1)*(y**(expar_nsi(1)-1))
     dvdx  =   0.0_rp
     dvdy  =   0.0_rp
     d2udy =   expar_nsi(1)*(expar_nsi(1)-1)*(y**(expar_nsi(1)-2))
     p     =   y
     dpdy  =   1.0_rp

  else if( kfl_exacs_nsi == 5 ) then
     !
     ! u =  2*y*(r-0.5)
     ! v = -2*x*(r-0.5)
     ! p = r-1
     !
     r     =   sqrt(x*x+y*y)
     u     =   2.0_rp*y*(r-0.50_rp)
     v     =  -2.0_rp*x*(r-0.50_rp)
     dudx  =   2.0_rp*x*y/r
     dudy  =   2.0_rp*(r-0.50_rp) + 2.0_rp*y*y/r
     dvdx  =  -2.0_rp*(r-0.50_rp) - 2.0_rp*x*x/r
     dvdy  =  -2.0_rp*x*y/r
     d2udy =   6.0_rp*y/r - 2.0_rp*y*y*y/(r*r*r)
     d2udx =   2.0_rp*y/r - 2.0_rp*x*x*y/(r*r*r)
     d2vdx =  -6.0_rp*x/r + 2.0_rp*x*x*x/(r*r*r)
     d2vdy =  -2.0_rp*x/r + 2.0_rp*x*y*y/(r*r*r)
     p     =  r-1.0_rp
     dpdx  =  x/r
     dpdy  =  y/r

  else if( kfl_exacs_nsi == 6 ) then
     !
     ! fx=a x^2 (1-x)^2 exp(b*x) cos(c*x)
     ! gy=d y^2 (1-y)^2 exp(e*y) cos(f*y)
     ! ft=cos(g*t)*exp(h*t)
     !
     ! u =  fx  * dgy * ft
     ! v = -dfx *  gy * ft
     ! p = i*x^2+j*y^2       ?
     !
     call spafun(expar_nsi(1),x,f,d1f,d2f,d3f) 
     call spafun(expar_nsi(4),y,g,d1g,d2g,d3g)
     freq = expar_nsi(7)*pi
     diss = -expar_nsi(8)
     ft   =       cos(freq*cutim)*exp(diss*cutim)
     dfdt = -freq*sin(freq*cutim)*exp(diss*cutim)          &
          &         +  diss*cos(freq*cutim)*exp(diss*cutim)
     u =  f  *d1g*ft
     v = -d1f*g  *ft
     dudx =  d1f*d1g*ft
     dudy =  f  *d2g*ft
     dvdx = -d2f*g  *ft
     dvdy = -d1f*d1g*ft
     d2udx = d2f*d1g*ft
     d2udy = f  *d3g*ft
     d2vdx =-d3f*g  *ft
     d2vdy =-d1f*d2g*ft
     !     p     = 100.0_rp*x*x
     !     dpdx  = 200.0_rp*x
     !     dpdy  = 0.0_rp
     if(kfl_timei_nsi/=0 ) then
        dudt =  f  *d1g*dfdt
        dvdt = -d1f*g  *dfdt
     end if

  else if( kfl_exacs_nsi == 7 ) then
     !
     !  Taylor vortex
     !
     if(kfl_timei_nsi/=0 ) then
        diss = expar_nsi(1)
        ft   =  exp(-2.0_rp*pi*pi*diss*cutim)
        dfdt = -2.0_rp*pi*pi*diss*ft
     else
        ft   =  1.0_rp
        dfdt =  0.0_rp       
     end if
     cx = cos(pi*x)
     cy = cos(pi*y)
     sx = sin(pi*x)
     sy = sin(pi*y)
     u = -cx*sy*ft
     v =  sx*cy*ft
     dudx =   pi*sx*sy*ft 
     dudy =  -pi*cx*cy*ft
     dvdx =   pi*cx*cy*ft
     dvdy =  -pi*sx*sy*ft
     d2udx =  pi*pi*cx*sy*ft
     d2udy =  pi*pi*cx*sy*ft
     d2vdx = -pi*pi*sx*cy*ft
     d2vdy = -pi*pi*sx*cy*ft
     p    = -0.25d0*(cos(2*pi*x) + cos(2*pi*y))*ft*ft
     dpdx =  0.5d0*pi*sin(2*pi*x)*ft*ft
     dpdy =  0.5d0*pi*sin(2*pi*y)*ft*ft
     if(kfl_timei_nsi/=0 ) then
        dudt = -cx*sy*dfdt
        dvdt =  sx*cy*dfdt
     end if

  else if( kfl_exacs_nsi == 8 ) then
     !
     !  Bochev
     !  P. Bochev, M. Gunzburger, and R. Lehoucq. On stabilized finite element methods for the Stokes problem in 
     !  the small time-step limit. International Journal for Numerical Methods in Fluids, 53:573–597, 2007. 
     !  Available from: http://dx.doi.org/10.1002/fld.1295. 
     !  sx(x) = sin(pi*x-0.7)
     !  sy(y) = sin(pi*y+0.2)
     !  cx(x) = cos(pi*x-0.7)
     !  cy(y) = cos(pi*y+0.2)
     !  splot sx(x)*sy(y) t 'u' , cx(x)*cy(y) t 'v', sin(x)*cos(y) t 'p'
     !
     sx    = sin(pi*x-0.7_rp)
     sy    = sin(pi*y+0.2_rp)

     cx    = cos(pi*x-0.7_rp)
     cy    = cos(pi*y+0.2_rp)

     u     =  sx*sy
     dudx  =  pi*cx*sy 
     d2udx = -pi*pi*sx*sy
     dudy  =  pi*sx*cy
     d2udy = -pi*pi*sx*sy

     v     =  cx*cy
     dvdx  = -pi*sx*cy
     d2vdx = -pi*pi*cx*cy
     dvdy  = -pi*cx*sy
     d2vdy = -pi*pi*cx*cy

     p    =  sin(x)*cos(y) 
     dpdx =  cos(x)*cos(y)
     dpdy = -sin(x)*sin(y)

  else if( kfl_exacs_nsi == -9 ) then
     !
     ! ux =  y*x^2
     ! uy = -x*y^2
     ! p  = x^2
     !
     if(kfl_timei_nsi/=0 ) then
        diss = -expar_nsi(1)
        !       ft   = exp(diss*cutim)
        !       dfdt = diss*exp(diss*cutim)
        ft   = 1.0_rp+diss*cutim  ! First order series
        dfdt = diss
     else
        ft   = 1.0_rp
        dfdt = 0.0_rp
     end if

     u    =  y*x*x*ft
     v    = -x*y*y*ft

     dudt =  y*x*x*dfdt
     dvdt = -x*y*y*dfdt

     dudx =  2.0_rp*y*x*ft
     dudy =  x*x*ft
     dvdx = -y*y*ft
     dvdy = -2.0_rp*x*y*ft

     d2udx =  2.0_rp*y*ft
     d2vdy = -2.0_rp*x*ft

     p     = x*x ! + 1.0_rp
     dpdx  = 2.0_rp*x
     dpdy  = 0.0_rp

  else if( kfl_exacs_nsi == 9 ) then
     !
     ! ux =  (1+x+y^2)
     ! uy = -(1+y+x^2)
     ! p  = 1+x+x^2
     !
     if(kfl_timei_nsi/=0 ) then
        diss = -expar_nsi(1)
        !       ft   = exp(diss*cutim)
        !       dfdt = diss*exp(diss*cutim)
        ft   = 1.0_rp+diss*cutim  ! First order series
        dfdt = diss
     else
        ft   = 1.0_rp
        dfdt = 0.0_rp
     end if

     u    =  (1.0_rp+x+y*y)*ft
     v    = -(1.0_rp+y+x*x)*ft

     dudt =  (1.0_rp+x+y*y)*dfdt
     dvdt = -(1.0_rp+y+x*x)*dfdt

     dudx =  ft
     dudy =  2.0_rp*y*ft
     dvdx = -2.0_rp*x*ft
     dvdy = -ft

     d2udx =  0.0_rp
     d2udy =  2.0_rp*ft
     d2vdx = -2.0_rp*ft
     d2vdy =  0.0_rp

     p     = 1.0_rp + x + x*x
     dpdx  = 1.0_rp + 2.0_rp*x
     dpdy  = 0.0_rp

  else if( kfl_exacs_nsi == 10 ) then
     !
     ! ux = y^2
     ! uy = x^2
     ! p  = 0
     !
     u     = y*y
     v     = x*x
     dudy  = 2.0_rp*y
     dvdx  = 2.0_rp*x
     d2udy = 2.0_rp
     d2vdx = 2.0_rp

  else if( kfl_exacs_nsi == 11 ) then
     !
     ! u =  u*/kap*ln[ (y+y0+delta)/y0 ]
     ! v =  0
     ! p =  0
     !
     if( itask == 1 ) then 
        ustar   =   0.23_rp
        kap     =   0.41_rp
        u       =   ustar / kap * log( (y+delta_dom+rough_dom)/rough_dom )
        dudy    =   ustar / ( kap * (y+delta_dom+rough_dom) )
        d2udy   =  -ustar / ( kap * (y+delta_dom+rough_dom)**2 )
     end if

  else if( kfl_exacs_nsi == 12 ) then
     !
     ! ux = 2 xy A
     ! uy = - A y^2
     ! p  = 7 x + B y
     !
     u     = x*y*2.0_rp*expar_nsi(1) 
     v     = - y**2 *expar_nsi(1) 
     dudx  = 2.0_rp * y*expar_nsi(1) 
     dudy  = 2.0_rp * x*expar_nsi(1) 
     dvdy  =  -2.0_rp*y*expar_nsi(1) 
     dvdx  = 0.0_rp
     d2udxdy = 2.0_rp*expar_nsi(1) 
     d2vdy = -2.0_rp*expar_nsi(1) 
     p     = 7.0_rp*x+expar_nsi(2) *y
     dpdx  = 7.0_rp
     dpdy  = expar_nsi(2) 

  else if( kfl_exacs_nsi == 13 ) then
     !
     ! ux =  2x +  y + 3z
     ! uy = - x - 3y - 4z
     ! uz =  6x + 7y - z
     ! p  = -5x + 6y + 2z => p(1,1,1)=3, p(0,0,0)=0
     !
     u    =   2.0_rp*x + y        + 3.0_rp*z
     v    =  -x        - 3.0_rp*y - 4.0_rp*z
     dudx =  2.0_rp
     dudy =  1.0_rp
     dvdx = -1.0_rp
     dvdy = -3.0_rp
     p    = -5.0_rp*x + 6.0_rp*y + 2.0_rp*z
     dpdx = -5.0_rp
     dpdy =  6.0_rp

     if( ndime == 3 ) then
        w    =   6.0_rp*x + 7.0_rp*y - z 
        dudz =  3.0_rp
        dvdz = -4.0_rp
        dwdx =  6.0_rp
        dwdy =  7.0_rp
        dwdz = -1.0_rp
        dpdz =  2.0_rp
     end if

  else if( kfl_exacs_nsi == 14 ) then
     !
     ! 
     !
     if(ndime==2)then
        u     =  x*x+y*y
        dudx  =  2.0_rp*x 
        d2udx =  2.0_rp
        dudy  =  2.0_rp*y
        d2udy =  2.0_rp

        v     =  x*x+y*y
        dvdx  =  2.0_rp*x 
        d2vdx =  2.0_rp
        dvdy  =  2.0_rp*y
        d2vdy =  2.0_rp

        p     =  x*x+y*y
        dpdx  =  2.0_rp*x 
        dpdy  =  2.0_rp*y
     else
        u     =  x*x + y*y + z*z
        dudx  =  2.0_rp*x 
        d2udx =  2.0_rp
        dudy  =  2.0_rp*y
        d2udy =  2.0_rp
        dudz  =  2.0_rp*z
        d2udz =  2.0_rp

        v     =  x*x + y*y + z*z
        dvdx  =  2.0_rp*x 
        d2vdx =  2.0_rp
        dvdy  =  2.0_rp*y
        d2vdy =  2.0_rp
        dvdz  =  2.0_rp*z
        d2vdz =  2.0_rp

        w     =  x*x + y*y + z*z
        dwdx  =  2.0_rp*x 
        d2wdx =  2.0_rp
        dwdy  =  2.0_rp*y
        d2wdy =  2.0_rp
        dwdz  =  2.0_rp*z
        d2wdz =  2.0_rp

        p     =  x*x + y*y + z*z
        dpdx  =  2.0_rp*x 
        dpdy  =  2.0_rp*y
        dpdz  =  2.0_rp*z


     end if

  else if( kfl_exacs_nsi == 15 ) then
     !
     ! u=2x+3y
     ! v=4x+5y
     ! p=6x+7y with p(1,1)=13
     !
     if( ndime == 2 ) then

        u     =  2.0_rp*x+3.0_rp*y
        dudx  =  2.0_rp
        dudy  =  3.0_rp

        v     =  4.0_rp*x+5.0_rp*y
        dvdx  =  4.0_rp
        dvdy  =  5.0_rp

        p     =  6.0_rp*x+7.0_rp*y
        dpdx  =  6.0_rp
        dpdy  =  7.0_rp

     else

        u     =  2.0_rp*x+3.0_rp*y+4.0_rp*z
        dudx  =  2.0_rp
        dudy  =  3.0_rp
        dudz  =  4.0_rp

        v     =  5.0_rp*x+6.0_rp*y+7.0_rp*z
        dvdx  =  5.0_rp
        dvdy  =  6.0_rp
        dvdz  =  7.0_rp

        w     =  8.0_rp*x+9.0_rp*y+10.0_rp*z
        dwdx  =  8.0_rp
        dwdy  =  9.0_rp
        dwdz  =  10.0_rp


        p     =  11.0_rp*x+12.0_rp*y+13.0_rp*z
        dpdx  =  11.0_rp
        dpdy  =  12.0_rp
        dpdz  =  13.0_rp

     end if

  else if (kfl_exacs_nsi==16)then
     !
     ! T=sin(pi*x)*sin(pi*y)*exp(x*y) in [0,1]x[0,1]
     !
     if( ndime == 2 ) then
        u      = sin(pi*x)*sin(pi*y)*exp(x*y)
        dudx   = sin(pi*y)*exp(x*y)*(pi*cos(pi*x)+y*sin(pi*x))
        dudy   = sin(pi*x)*exp(x*y)*(pi*cos(pi*y)+x*sin(pi*y))
        d2udx = sin(pi*y)*exp(x*y)*(2.0_rp*y*pi*cos(pi*x)&
             &   -pi*pi*sin(pi*x)+y*y*sin(pi*x))
        d2udy = sin(pi*x)*exp(x*y)*(2.0_rp*x*pi*cos(pi*y)&
             &   -pi*pi*sin(pi*y)+x*x*sin(pi*y))
        v      = sin(pi*x)*sin(pi*y)*exp(x*y)
        dvdx   = sin(pi*y)*exp(x*y)*(pi*cos(pi*x)+y*sin(pi*x))
        dvdy   = sin(pi*x)*exp(x*y)*(pi*cos(pi*y)+x*sin(pi*y))
        d2vdx = sin(pi*y)*exp(x*y)*(2.0_rp*y*pi*cos(pi*x)&
             &   -pi*pi*sin(pi*x)+y*y*sin(pi*x))
        d2vdy = sin(pi*x)*exp(x*y)*(2.0_rp*x*pi*cos(pi*y)&
             &   -pi*pi*sin(pi*y)+x*x*sin(pi*y))

        p      = sin(pi*x)*sin(pi*y)*exp(x*y)
        dpdx   = sin(pi*y)*exp(x*y)*(pi*cos(pi*x)+y*sin(pi*x))
        dpdy   = sin(pi*x)*exp(x*y)*(pi*cos(pi*y)+x*sin(pi*y))

     else
        !
        ! T=sin(pi*x)*sin(pi*y)*sin(pi*z)*exp(x*y*z) in [0,1]x[0,1]x[0,1]
        !
        u      = sin(pi*x)*sin(pi*y)*sin(pi*z)*exp(x*y*z)
        dudx   = sin(pi*y)*sin(pi*z)*exp(x*y*z)*(pi*cos(pi*x)+y*z*sin(pi*x))
        dudy   = sin(pi*x)*sin(pi*z)*exp(x*y*z)*(pi*cos(pi*y)+x*z*sin(pi*y))
        dudz   = sin(pi*x)*sin(pi*y)*exp(x*y*z)*(pi*cos(pi*z)+x*y*sin(pi*z))
        d2udx =  sin(pi*y)*sin(pi*z)*exp(x*y*z)*(2.0_rp*y*z*pi*cos(pi*x)+&
             &   (-pi*pi+y*y*z*z)*sin(pi*x))
        d2udy =  sin(pi*x)*sin(pi*z)*exp(x*y*z)*(2.0_rp*x*z*pi*cos(pi*y)+&
             &   (-pi*pi+x*x*z*z)*sin(pi*y))
        d2udz =  sin(pi*x)*sin(pi*y)*exp(x*y*z)*(2.0_rp*x*y*pi*cos(pi*z)+&
             &   (-pi*pi+x*x*y*y)*sin(pi*z)) 
        d2udxdy = (pi*cos(pi*x)+y*z*sin(pi*x))*sin(pi*z)*exp(x*y*z)*(pi*cos(pi*y)+x*z*sin(pi*y))+&
             &  z*sin(pi*x)*sin(pi*y)*sin(pi*z)*exp(x*y*z)
        d2udxdz = (pi*cos(pi*x)+y*z*sin(pi*x))*sin(pi*y)*exp(x*y*z)*(pi*cos(pi*z)+x*y*sin(pi*z))+&
             &  y*sin(pi*x)*sin(pi*y)*sin(pi*z)*exp(x*y*z)
        d2udydz = x*sin(pi*x)*sin(pi*y)*exp(x*y*z)*(y*(pi*cos(pi*y)+x*z*sin(pi*y))+sin(pi*y))

        v      = sin(pi*x)*sin(pi*y)*sin(pi*z)*exp(x*y*z)
        dvdx   = sin(pi*y)*sin(pi*z)*exp(x*y*z)*(pi*cos(pi*x)+y*z*sin(pi*x))
        dvdy   = sin(pi*x)*sin(pi*z)*exp(x*y*z)*(pi*cos(pi*y)+x*z*sin(pi*y))
        dvdz   = sin(pi*x)*sin(pi*y)*exp(x*y*z)*(pi*cos(pi*z)+x*y*sin(pi*z))
        d2vdx  = sin(pi*y)*sin(pi*z)*exp(x*y*z)*(2.0_rp*y*z*pi*cos(pi*x)+&
             &   (-pi*pi+y*y*z*z)*sin(pi*x))
        d2vdy = sin(pi*x)*sin(pi*z)*exp(x*y*z)*(2.0_rp*x*z*pi*cos(pi*y)+&
             &   (-pi*pi+x*x*z*z)*sin(pi*y))
        d2vdz = sin(pi*x)*sin(pi*y)*exp(x*y*z)*(2.0_rp*x*y*pi*cos(pi*z)+&
             &   (-pi*pi+x*x*y*y)*sin(pi*z)) 
        d2vdxdy = (pi*cos(pi*x)+y*z*sin(pi*x))*sin(pi*z)*exp(x*y*z)*(pi*cos(pi*y)+x*z*sin(pi*y))+&
             &  z*sin(pi*x)*sin(pi*y)*sin(pi*z)*exp(x*y*z)
        d2vdxdz = (pi*cos(pi*x)+y*z*sin(pi*x))*sin(pi*y)*exp(x*y*z)*(pi*cos(pi*z)+x*y*sin(pi*z))+&
             &  y*sin(pi*x)*sin(pi*y)*sin(pi*z)*exp(x*y*z)
        d2vdydz = x*sin(pi*x)*sin(pi*y)*exp(x*y*z)*(y*(pi*cos(pi*y)+x*z*sin(pi*y))+sin(pi*y))

        w      = sin(pi*x)*sin(pi*y)*sin(pi*z)*exp(x*y*z)
        dwdx   = sin(pi*y)*sin(pi*z)*exp(x*y*z)*(pi*cos(pi*x)+y*z*sin(pi*x))
        dwdy   = sin(pi*x)*sin(pi*z)*exp(x*y*z)*(pi*cos(pi*y)+x*z*sin(pi*y))
        dwdz   = sin(pi*x)*sin(pi*y)*exp(x*y*z)*(pi*cos(pi*z)+x*y*sin(pi*z))
        d2wdx  = sin(pi*y)*sin(pi*z)*exp(x*y*z)*(2.0_rp*y*z*pi*cos(pi*x)+&
             &   (-pi*pi+y*y*z*z)*sin(pi*x))
        d2wdy = sin(pi*x)*sin(pi*z)*exp(x*y*z)*(2.0_rp*x*z*pi*cos(pi*y)+&
             &   (-pi*pi+x*x*z*z)*sin(pi*y))
        d2wdz = sin(pi*x)*sin(pi*y)*exp(x*y*z)*(2.0_rp*x*y*pi*cos(pi*z)+&
             &   (-pi*pi+x*x*y*y)*sin(pi*z)) 
        d2wdxdy = (pi*cos(pi*x)+y*z*sin(pi*x))*sin(pi*z)*exp(x*y*z)*(pi*cos(pi*y)+x*z*sin(pi*y))+&
             &  z*sin(pi*x)*sin(pi*y)*sin(pi*z)*exp(x*y*z)
        d2wdxdz = (pi*cos(pi*x)+y*z*sin(pi*x))*sin(pi*y)*exp(x*y*z)*(pi*cos(pi*z)+x*y*sin(pi*z))+&
             &  y*sin(pi*x)*sin(pi*y)*sin(pi*z)*exp(x*y*z)
        d2wdydz = x*sin(pi*x)*sin(pi*y)*exp(x*y*z)*(y*(pi*cos(pi*y)+x*z*sin(pi*y))+sin(pi*y))

        p      = sin(pi*x)*sin(pi*y)*sin(pi*z)*exp(x*y*z)
        dpdx   = sin(pi*y)*sin(pi*z)*exp(x*y*z)*(pi*cos(pi*x)+y*z*sin(pi*x))
        dpdy   = sin(pi*x)*sin(pi*z)*exp(x*y*z)*(pi*cos(pi*y)+x*z*sin(pi*y))
        dpdz   = sin(pi*x)*sin(pi*y)*exp(x*y*z)*(pi*cos(pi*z)+x*y*sin(pi*z))
     end if

  else if (kfl_exacs_nsi==17)then
     !
     ! T=sin(pi*x)*sin(pi*y) in [0,1]x[0,1]
     !
     if( ndime == 2 ) then
        u      = sin(pi*x)*sin(pi*y)
        dudx   = sin(pi*y)*(pi*cos(pi*x)) 
        dudy   = sin(pi*x)*(pi*cos(pi*y))
        d2udx  = sin(pi*y)*(-pi*pi*sin(pi*x))
        d2udy  = sin(pi*x)*(-pi*pi*sin(pi*y))
        d2udxdy = pi*pi*cos(pi*y)*cos(pi*x)
        d2udxdz = pi*pi*cos(pi*x)*sin(pi*y)

        v      = sin(pi*x)*sin(pi*y)
        dvdx   = sin(pi*y)*(pi*cos(pi*x))
        dvdy   = sin(pi*x)*(pi*cos(pi*y))
        d2vdx  = sin(pi*y)*(-pi*pi*sin(pi*x))
        d2vdy  = sin(pi*x)*(-pi*pi*sin(pi*y))
        d2vdxdy = pi*pi*cos(pi*y)*cos(pi*x)
        d2vdxdz = pi*pi*cos(pi*x)*sin(pi*y)

        p      = sin(pi*x)*sin(pi*y)
        dpdx   = sin(pi*y)*(pi*cos(pi*x))
        dpdy   = sin(pi*x)*(pi*cos(pi*y))
     else
        !
        ! T=sin(pi*x)*sin(pi*y)*sin(pi*z) in [0,1]x[0,1]x[0,1]
        !
        u      = sin(pi*x)*sin(pi*y)*sin(pi*z)
        dudx   = pi*sin(pi*y)*sin(pi*z)*cos(pi*x)
        dudy   = pi*sin(pi*x)*sin(pi*z)*cos(pi*y)
        dudz   = pi*sin(pi*x)*sin(pi*y)*cos(pi*z)
        d2udx =  -pi*pi*sin(pi*x)*sin(pi*y)*sin(pi*z)
        d2udy =  -pi*pi*sin(pi*x)*sin(pi*y)*sin(pi*z)
        d2udz =  -pi*pi*sin(pi*x)*sin(pi*y)*sin(pi*z)
        d2udxdy = pi*pi*cos(pi*y)*cos(pi*x)*sin(pi*z)
        d2udxdz = pi*pi*cos(pi*z)*cos(pi*x)*sin(pi*y)
        d2udydz = pi*pi*cos(pi*z)*cos(pi*y)*sin(pi*x)

        v      = sin(pi*x)*sin(pi*y)*sin(pi*z)
        dvdx   = pi*sin(pi*y)*sin(pi*z)*cos(pi*x)
        dvdy   = pi*sin(pi*x)*sin(pi*z)*cos(pi*y)
        dvdz   = pi*sin(pi*x)*sin(pi*y)*cos(pi*z)
        d2vdxdy =pi*pi*cos(pi*y)*cos(pi*x)*sin(pi*z) 
        d2vdxdz =pi*pi*cos(pi*z)*cos(pi*x)*sin(pi*y) 
        d2vdydz =pi*pi*cos(pi*z)*cos(pi*y)*sin(pi*x) 
        d2vdx =  -pi*pi*sin(pi*x)*sin(pi*y)*sin(pi*z)
        d2vdy =  -pi*pi*sin(pi*x)*sin(pi*y)*sin(pi*z)
        d2vdz =  -pi*pi*sin(pi*x)*sin(pi*y)*sin(pi*z)

        w      = sin(pi*x)*sin(pi*y)*sin(pi*z)
        dwdx   = pi*sin(pi*y)*sin(pi*z)*cos(pi*x)
        dwdy   = pi*sin(pi*x)*sin(pi*z)*cos(pi*y)
        dwdz   = pi*sin(pi*x)*sin(pi*y)*cos(pi*z)
        d2wdxdy = pi*pi*cos(pi*y)*cos(pi*x)*sin(pi*z)
        d2wdxdz = pi*pi*cos(pi*z)*cos(pi*x)*sin(pi*y)
        d2wdydz = pi*pi*cos(pi*z)*cos(pi*y)*sin(pi*x)
        d2wdx =  -pi*pi*sin(pi*x)*sin(pi*y)*sin(pi*z)
        d2wdy =  -pi*pi*sin(pi*x)*sin(pi*y)*sin(pi*z)
        d2wdz =  -pi*pi*sin(pi*x)*sin(pi*y)*sin(pi*z)


        p      = sin(pi*x)*sin(pi*y)*sin(pi*z)
        dpdx   = pi*sin(pi*y)*sin(pi*z)*cos(pi*x)
        dpdy   = pi*sin(pi*x)*sin(pi*z)*cos(pi*y)
        dpdz   = pi*sin(pi*x)*sin(pi*y)*cos(pi*z)

     end if

  else if( kfl_exacs_nsi == 18 ) then
     !
     ! ux = y**n
     ! uy = 0
     ! p  = y
     !     
     if( ndime == 2 ) then
        u       = y*y
        v       = x*x
        dudy    = 2.0_rp*y
        dvdx    = 2.0_rp*x
        d2udy   = 2.0_rp
        d2vdx   = 2.0_rp
        p       = x*x+y*y
        dpdx    = 2.0_rp*x
        dpdy    = 2.0_rp*y
     else
        u     = y*y + z*z
        dudy  = 2.0_rp*y
        d2udy = 2.0_rp
        dudz  = 2.0_rp*z
        d2udz = 2.0_rp

        v     = x*x + z*z
        dvdx  = 2.0_rp*x
        d2vdx = 2.0_rp
        dvdz  = 2.0_rp*z
        d2vdz = 2.0_rp

        w     = x*x + y*y
        dwdx  = 2.0_rp*x
        d2wdx = 2.0_rp
        dwdy  = 2.0_rp*y
        d2wdy = 2.0_rp

        p    = x*x + y*y + z*z
        dpdx = 2.0_rp*x
        dpdy = 2.0_rp*y
        dpdz = 2.0_rp*z
     end if

  else if( kfl_exacs_nsi == 19 ) then
     !
     !  Bochev
     !  P. Bochev, M. Gunzburger, and R. Lehoucq. On stabilized finite element methods for the Stokes problem in 
     !  the small time-step limit. International Journal for Numerical Methods in Fluids, 53:573–597, 2007. 
     !  Available from: http://dx.doi.org/10.1002/fld.1295. 
     !  sx(x) = sin(pi*x-0.7)
     !  sy(y) = sin(pi*y+0.2)
     !  cx(x) = cos(pi*x-0.7)
     !  cy(y) = cos(pi*y+0.2)
     !  splot sx(x)*sy(y) t 'u' , cx(x)*cy(y) t 'v', sin(x)*cos(y) t 'p'
     !
     sx    = sin(pi*x-0.7_rp)
     sy    = sin(pi*y+0.2_rp)

     cx    = cos(pi*x-0.7_rp)
     cy    = cos(pi*y+0.2_rp)

     u     =  sy
     dudy  =  pi*cy
     d2udy = -pi*pi*sy

     v     =  cx
     dvdx  = -pi*sx
     d2vdx = -pi*pi*cx

     p    =  sin(x)*cos(y) 
     dpdx =  cos(x)*cos(y)
     dpdy = -sin(x)*sin(y)

  else if( kfl_exacs_nsi == 20 ) then
     !
     ! Green Taylor vortex
     !
     u       =       -cos(pi*x)*sin(pi*y)
     dudx    =     pi*sin(pi*x)*sin(pi*y)
     dudy    =    -pi*cos(pi*x)*cos(pi*y)
     d2udx   =  pi*pi*cos(pi*x)*sin(pi*y)
     d2udy   =  pi*pi*cos(pi*x)*sin(pi*y)
     d2udxdy =  pi*pi*sin(pi*x)*cos(pi*y)

     v       =         sin(pi*x)*cos(pi*y)
     dvdx    =      pi*cos(pi*x)*cos(pi*y)
     dvdy    =     -pi*sin(pi*x)*sin(pi*y)
     d2vdx   =  -pi*pi*sin(pi*x)*cos(pi*y)
     d2vdy   =  -pi*pi*sin(pi*x)*cos(pi*y)
     d2vdxdy =  -pi*pi*cos(pi*x)*sin(pi*y)

     f       =  sin(pi*0.0_rp)*sin(pi*0.0_rp) + sin(pi*0.5_rp)*sin(pi*0.5_rp)
     p       =  (sin(pi*x)*sin(pi*x) + sin(pi*y)*sin(pi*y))/f
     dpdx    =  (2.0_rp*pi*cos(pi*x)*sin(pi*x))/f
     dpdy    =  (2.0_rp*pi*cos(pi*y)*sin(pi*y))/f

  else if( kfl_exacs_nsi == 22 ) then

     if(ndime==2)then
        u       =              -cos(2.0_rp*pi*x)*sin(2.0_rp*pi*y)
        dudx    =     2.0_rp*pi*sin(2.0_rp*pi*x)*sin(2.0_rp*pi*y)
        dudy    =    -2.0_rp*pi*cos(2.0_rp*pi*x)*cos(2.0_rp*pi*y)
        d2udx   =  4.0_rp*pi*pi*cos(2.0_rp*pi*x)*sin(2.0_rp*pi*y)
        d2udy   =  4.0_rp*pi*pi*cos(2.0_rp*pi*x)*sin(2.0_rp*pi*y)
        d2udxdy =  4.0_rp*pi*pi*sin(2.0_rp*pi*x)*cos(2.0_rp*pi*y)

        v       =                sin(2.0_rp*pi*x)*cos(2.0_rp*pi*y)
        dvdx    =      2.0_rp*pi*cos(2.0_rp*pi*x)*cos(2.0_rp*pi*y)
        dvdy    =     -2.0_rp*pi*sin(2.0_rp*pi*x)*sin(2.0_rp*pi*y)
        d2vdx   =  -4.0_rp*pi*pi*sin(2.0_rp*pi*x)*cos(2.0_rp*pi*y)
        d2vdy   =  -4.0_rp*pi*pi*sin(2.0_rp*pi*x)*cos(2.0_rp*pi*y)
        d2vdxdy =  -4.0_rp*pi*pi*cos(2.0_rp*pi*x)*sin(2.0_rp*pi*y)

        p       =  2.0_rp*x + 3.0_rp*y 
        dpdx    =  2.0_rp
        dpdy    =  3.0_rp

     else

        u       =              -2.0_rp*cos(2.0_rp*pi*x)*sin(2.0_rp*pi*y)*sin(2.0_rp*pi*z)
        dudx    =     2.0_rp*2.0_rp*pi*sin(2.0_rp*pi*x)*sin(2.0_rp*pi*y)*sin(2.0_rp*pi*z)
        dudy    =    -2.0_rp*2.0_rp*pi*cos(2.0_rp*pi*x)*cos(2.0_rp*pi*y)*sin(2.0_rp*pi*z)
        dudz    =    -2.0_rp*2.0_rp*pi*cos(2.0_rp*pi*x)*sin(2.0_rp*pi*y)*cos(2.0_rp*pi*z)
        d2udx   =   8.0_rp*pi*pi*cos(2.0_rp*pi*x)*sin(2.0_rp*pi*y)*sin(2.0_rp*pi*z)
        d2udy   =   8.0_rp*pi*pi*cos(2.0_rp*pi*x)*sin(2.0_rp*pi*y)*sin(2.0_rp*pi*z)
        d2udz   =   8.0_rp*pi*pi*cos(2.0_rp*pi*x)*sin(2.0_rp*pi*y)*sin(2.0_rp*pi*z)
        d2udxdy =   8.0_rp*pi*pi*sin(2.0_rp*pi*x)*cos(2.0_rp*pi*y)*sin(2.0_rp*pi*z)
        d2udxdz =   8.0_rp*pi*pi*sin(2.0_rp*pi*x)*sin(2.0_rp*pi*y)*cos(2.0_rp*pi*z)
        d2udydz =  -8.0_rp*pi*pi*cos(2.0_rp*pi*x)*cos(2.0_rp*pi*y)*cos(2.0_rp*pi*z)

        v       =                sin(2.0_rp*pi*x)*cos(2.0_rp*pi*y)*sin(2.0_rp*pi*z)
        dvdx    =      2.0_rp*pi*cos(2.0_rp*pi*x)*cos(2.0_rp*pi*y)*sin(2.0_rp*pi*z)
        dvdy    =     -2.0_rp*pi*sin(2.0_rp*pi*x)*sin(2.0_rp*pi*y)*sin(2.0_rp*pi*z)
        dvdz    =      2.0_rp*pi*sin(2.0_rp*pi*x)*cos(2.0_rp*pi*y)*cos(2.0_rp*pi*z)
        d2vdx   =  -4.0_rp*pi*pi*sin(2.0_rp*pi*x)*cos(2.0_rp*pi*y)*sin(2.0_rp*pi*z)
        d2vdy   =  -4.0_rp*pi*pi*sin(2.0_rp*pi*x)*cos(2.0_rp*pi*y)*sin(2.0_rp*pi*z)
        d2vdz   =  -4.0_rp*pi*pi*sin(2.0_rp*pi*x)*cos(2.0_rp*pi*y)*sin(2.0_rp*pi*z)
        d2vdxdy =  -4.0_rp*pi*pi*cos(2.0_rp*pi*x)*sin(2.0_rp*pi*y)*sin(2.0_rp*pi*z)
        d2vdxdz =   4.0_rp*pi*pi*cos(2.0_rp*pi*x)*cos(2.0_rp*pi*y)*cos(2.0_rp*pi*z)
        d2vdydz =  -4.0_rp*pi*pi*sin(2.0_rp*pi*x)*sin(2.0_rp*pi*y)*cos(2.0_rp*pi*z)

        w       =                sin(2.0_rp*pi*x)*sin(2.0_rp*pi*y)*cos(2.0_rp*pi*z)
        dwdx    =      2.0_rp*pi*cos(2.0_rp*pi*x)*sin(2.0_rp*pi*y)*cos(2.0_rp*pi*z)
        dwdy    =      2.0_rp*pi*sin(2.0_rp*pi*x)*cos(2.0_rp*pi*y)*cos(2.0_rp*pi*z)
        dwdz    =     -2.0_rp*pi*sin(2.0_rp*pi*x)*sin(2.0_rp*pi*y)*sin(2.0_rp*pi*z)
        d2wdx   =  -4.0_rp*pi*pi*sin(2.0_rp*pi*x)*sin(2.0_rp*pi*y)*cos(2.0_rp*pi*z)
        d2wdy   =  -4.0_rp*pi*pi*sin(2.0_rp*pi*x)*sin(2.0_rp*pi*y)*cos(2.0_rp*pi*z)
        d2wdz   =  -4.0_rp*pi*pi*sin(2.0_rp*pi*x)*sin(2.0_rp*pi*y)*cos(2.0_rp*pi*z)
        d2wdxdy =   4.0_rp*pi*pi*cos(2.0_rp*pi*x)*cos(2.0_rp*pi*y)*cos(2.0_rp*pi*z)
        d2wdxdz =  -4.0_rp*pi*pi*cos(2.0_rp*pi*x)*sin(2.0_rp*pi*y)*sin(2.0_rp*pi*z)
        d2wdydz =  -4.0_rp*pi*pi*sin(2.0_rp*pi*x)*cos(2.0_rp*pi*y)*sin(2.0_rp*pi*z)

        p       = 2.0_rp*x+3.0_rp*y+4.0_rp*z
        dpdx    = 2.0_rp
        dpdy    = 3.0_rp  
        dpdz    = 4.0_rp

     end if

  else if( kfl_exacs_nsi == 23 ) then

     u       =               -cos(2.0_rp*pi*x)*sin(2.0_rp*pi*y)
     dudx    =      2.0_rp*pi*sin(2.0_rp*pi*x)*sin(2.0_rp*pi*y)
     dudy    =     -2.0_rp*pi*cos(2.0_rp*pi*x)*cos(2.0_rp*pi*y)
     d2udx   =   4.0_rp*pi*pi*cos(2.0_rp*pi*x)*sin(2.0_rp*pi*y)
     d2udy   =   4.0_rp*pi*pi*cos(2.0_rp*pi*x)*sin(2.0_rp*pi*y)
     d2udxdy =   4.0_rp*pi*pi*sin(2.0_rp*pi*x)*cos(2.0_rp*pi*y)

     v       =                sin(2.0_rp*pi*x)*cos(2.0_rp*pi*y)
     dvdx    =      2.0_rp*pi*cos(2.0_rp*pi*x)*cos(2.0_rp*pi*y)
     dvdy    =     -2.0_rp*pi*sin(2.0_rp*pi*x)*sin(2.0_rp*pi*y)
     d2vdx   =  -4.0_rp*pi*pi*sin(2.0_rp*pi*x)*cos(2.0_rp*pi*y)
     d2vdy   =  -4.0_rp*pi*pi*sin(2.0_rp*pi*x)*cos(2.0_rp*pi*y)
     d2vdxdy =  -4.0_rp*pi*pi*cos(2.0_rp*pi*x)*sin(2.0_rp*pi*y)

     p       = 2.0_rp*x+3.0_rp*y
     dpdx    = 2.0_rp
     dpdy    = 3.0_rp  

  else if( kfl_exacs_nsi == 24 ) then

     u       =  y*y
     dudy    =  2.0_rp*y 
     d2udy   =  2.0_rp
     p       =  1.0_rp
     dpdx    =  0.0_rp

  else if( kfl_exacs_nsi == 25 ) then

     u    =  (1.0_rp+x+y*y+z*z)
     v    = -(1.0_rp+y+x*x+z*z)
     w    =  (x*x + y*y)

     dudx =  1.0_rp
     dudy =  2.0_rp*y
     dudz =  2.0_rp*z

     dvdx = -2.0_rp*x
     dvdy = -1.0_rp
     dvdz = -2.0_rp*z

     dwdx =  2.0_rp*x
     dwdy =  2.0_rp*y
     dwdz =  0.0_rp

     d2udx =  0.0_rp
     d2udy =  2.0_rp
     d2udz =  2.0_rp

     d2vdx = -2.0_rp
     d2vdy =  0.0_rp
     d2vdz = -2.0_rp

     d2wdx =  2.0_rp
     d2wdy =  2.0_rp
     d2wdz =  0.0_rp

     p     =  x*x +y*y
     dpdx  =  2.0_rp*x
     dpdy  =  2.0_rp*y

  else if( kfl_exacs_nsi == 26 ) then
     !
     ! ux =  x
     ! uy = -y
     ! p  =  1
     !
     u     =  x
     v     = -y
     dudx  =  1.0_rp
     dvdy  = -1.0_rp
     p     =  1.0_rp

  else if( kfl_exacs_nsi == 27 ) then
     !
     ! ux =  2x^2 +  y^2 + 3z^2
     ! uy = - x^2 - 3y^2 - 4z^2
     ! uz =  6x^2 + 7y^2 - z^2
     ! p  = -5x^2 + 6y^2 + 2z^2 => p(1,1,1)=3
     !
     u    =   2.0_rp*x*x + y*y        + 3.0_rp*z*z
     v    =  -x*x        - 3.0_rp*y*y - 4.0_rp*z*z
     dudx =  4.0_rp*x
     dudy =  2.0_rp*y
     dvdx = -2.0_rp*x
     dvdy = -6.0_rp*y
     p    = -5.0_rp*x*x + 6.0_rp*y*y + 2.0_rp*z*z
     dpdx = -10.0_rp*x
     dpdy =  12.0_rp*y
     
     d2udx =  4.0_rp
     d2udy =  2.0_rp
     d2vdx = -2.0_rp
     d2vdy = -6.0_rp

     if( ndime == 3 ) then
        w     =  6.0_rp*x*x + 7.0_rp*y*y - z*z 
        dudz  =  6.0_rp*z
        dvdz  = -8.0_rp*z
        dwdx  = 12.0_rp*x
        dwdy  = 14.0_rp*y
        dwdz  = -2.0_rp*z
        dpdz  =  4.0_rp*z
        d2wdx = 12.0_rp
        d2wdy = 14.0_rp
        d2wdz = -2.0_rp
        d2udz =  6.0_rp
        d2vdz = -8.0_rp
     end if

  else if( kfl_exacs_nsi == 28 ) then
     !
     ! ux =  x * cos(f*t)
     ! uy = -y * cos(f*t)
     ! p  =  1 * cos(f*t)
     !
     f     =  5.0_rp
     u     =  x * cos(f*cutim) * cos(2.0_rp*f*cutim)
     v     = -y * cos(f*cutim) * cos(2.0_rp*f*cutim)
     dudx  =      cos(f*cutim) * cos(2.0_rp*f*cutim)
     dvdy  =     -cos(f*cutim) * cos(2.0_rp*f*cutim)
     dudt  =  - f * x * ( sin(f*cutim) * cos(2.0_rp*f*cutim) + cos(f*cutim) * sin(2.0_rp*f*cutim) )
     dvdt  =    f * y * ( sin(f*cutim) * cos(2.0_rp*f*cutim) + cos(f*cutim) * sin(2.0_rp*f*cutim) )
     p     =  0.0_rp 

  else if( kfl_exacs_nsi == 29 ) then

     f       =   5.0_rp

     u       =               -cos(2.0_rp*pi*x)*sin(2.0_rp*pi*y) * cos(f*cutim) 
     dudx    =      2.0_rp*pi*sin(2.0_rp*pi*x)*sin(2.0_rp*pi*y) * cos(f*cutim) 
     dudy    =     -2.0_rp*pi*cos(2.0_rp*pi*x)*cos(2.0_rp*pi*y) * cos(f*cutim) 
     d2udx   =   4.0_rp*pi*pi*cos(2.0_rp*pi*x)*sin(2.0_rp*pi*y) * cos(f*cutim) 
     d2udy   =   4.0_rp*pi*pi*cos(2.0_rp*pi*x)*sin(2.0_rp*pi*y) * cos(f*cutim) 
     d2udxdy =   4.0_rp*pi*pi*sin(2.0_rp*pi*x)*cos(2.0_rp*pi*y) * cos(f*cutim) 

     v       =                sin(2.0_rp*pi*x)*cos(2.0_rp*pi*y) * cos(f*cutim) 
     dvdx    =      2.0_rp*pi*cos(2.0_rp*pi*x)*cos(2.0_rp*pi*y) * cos(f*cutim) 
     dvdy    =     -2.0_rp*pi*sin(2.0_rp*pi*x)*sin(2.0_rp*pi*y) * cos(f*cutim) 
     d2vdx   =  -4.0_rp*pi*pi*sin(2.0_rp*pi*x)*cos(2.0_rp*pi*y) * cos(f*cutim) 
     d2vdy   =  -4.0_rp*pi*pi*sin(2.0_rp*pi*x)*cos(2.0_rp*pi*y) * cos(f*cutim) 
     d2vdxdy =  -4.0_rp*pi*pi*cos(2.0_rp*pi*x)*sin(2.0_rp*pi*y) * cos(f*cutim) 

     dudt    =   f * cos(2.0_rp*pi*x)*sin(2.0_rp*pi*y) * sin(f*cutim) 
     dvdt    =  -f * sin(2.0_rp*pi*x)*cos(2.0_rp*pi*y) * sin(f*cutim) 

     p       =   2.0_rp*x+3.0_rp*y
     dpdx    =   2.0_rp
     dpdy    =   3.0_rp  

  else if( kfl_exacs_nsi == 30 ) then

     expt    =  exp(max(-1.0e2_rp,-2.0_rp*cutim))

     u       = -cos(x)*sin(y)*expt
     dudx    =  sin(x)*sin(y)*expt
     dudy    = -cos(x)*cos(y)*expt
     d2udx   =  cos(x)*sin(y)*expt
     d2udy   =  cos(x)*sin(y)*expt
     d2udxdy =  sin(x)*cos(y)*expt

     v       =  sin(x)*cos(y)*expt
     dvdx    =  cos(x)*cos(y)*expt
     dvdy    = -sin(x)*sin(y)*expt
     d2vdx   = -sin(x)*cos(y)*expt
     d2vdy   = -sin(x)*cos(y)*expt 
     d2vdxdy = -cos(x)*sin(y)*expt

     dudt    =  2.0_rp*cos(x)*sin(y)*expt
     dvdt    = -2.0_rp*sin(x)*cos(y)*expt

     expt    =  exp(max(-1.0e2_rp,-4.0_rp*cutim))
     p       =  0.25_rp*(sin(2.0_rp*x)+sin(2.0_rp*y))*expt
     dpdx    =  0.25_rp*2.0_rp*(cos(2.0_rp*x))*expt
     dpdx    =  0.25_rp*2.0_rp*(cos(2.0_rp*y))*expt

  else if( kfl_exacs_nsi == 31 ) then
     !
     ! Taylor-Couette flow
     ! Omega= angular velocity
     ! n=     ratio inner/outer radius
     ! ri=    inner radius
     ! ro=    outer radius
     ! 
     !
!!$     ri      = 0.3_rp
!!$     ro      = 1.0_rp
!!$     omega   = 2.0_rp*pi
!!$     n       = ri/ro
!!$     r       = sqrt(x**2+y**2+z**2)
!!$     expt    = 1.0_rp/(1.0_rp-n**2)*(-n**2*r+ri**2/r)
!!$     vt      = omega*expt
!!$     if(      x == 0.0_rp ) then
!!$        if( y >= 0.0_rp ) then
!!$           u = -vt
!!$           y =  0.0_rp
!!$        else
!!$           u = vt
!!$           y = 0.0_rp
!!$        end if
!!$     else        
!!$        theta   = atan(abs(y/x))
!!$        if( x >= 0.0_rp .and. y >= 0.0_rp ) then
!!$           u    = -vt*cos(theta)
!!$           y    =  vt*sin(theta)
!!$        else if( x <= 0.0_rp .and. y >= 0.0_rp ) then
!!$           u    = -vt*cos(theta)
!!$           y    = -vt*sin(theta)
!!$        else if( x <= 0.0_rp .and. y <= 0.0_rp ) then
!!$           u    =  vt*cos(theta)
!!$           y    = -vt*sin(theta)
!!$        else if( x <= 0.0_rp .and. y <= 0.0_rp ) then
!!$           u    =  vt*cos(theta)
!!$           y    =  vt*sin(theta)
!!$        end if
!!$     end if
  else if( kfl_exacs_nsi == 32 ) then
     !
     ! u = 1.25*(0.4^2 - r^2)  
     ! v = 0.0
     ! w = 0.0 
     ! p = -10.31839621875*x+4.127358487500000
     !
     u     =   -1.25*(0.16-(x*x+y*y+z*z))
     !print*,'u',u,x,y,z
     dudx  =   -2.5*x
     dudy  =   -2.5*y
     dudz  =   -2.5*z
     d2udx =   -2.5
     d2udy =   -2.5
     d2udz =   -2.5
     p     =  0.175*x+0.7 
     dpdx  =  0.175
!     p     =  -10.318396264151*x+4.12735849056604  
     !     dpdx  =  -10.318396264151

  else if( kfl_exacs_nsi == 33 ) then
     u     = -7.8125*y*y+1.25
     dudy  = -15.625*y
     d2udy = -15.625
     p     = -2.0*x
     dpdx  = -2.0
  end if

  if( itask == 1 ) then 
     !
     ! Exact velocity, pressure and gradients
     !
     exvel(1)   = u
     exvel(2)   = v
     exveg(1,1) = dudx
     exveg(1,2) = dvdx
     exveg(2,1) = dudy
     exveg(2,2) = dvdy
     expre      = p
     exprg(1)   = dpdx
     exprg(2)   = dpdy
     if(ndime==3 ) then
        exvel(3)   = w
        exveg(1,3) = dwdx
        exveg(2,3) = dwdy
        exveg(3,3) = dwdz
        exveg(3,1) = dudz
        exveg(3,2) = dvdz
        exprg(3)   = dpdz
     end if

  else if( itask == 2 ) then
     !
     ! Force vector
     !
     gprhs(1) = gprhs(1)              + dpdx  +  sig * u
     gprhs(2) = gprhs(2)              + dpdy  +  sig * v
     if(ndime==3) gprhs(3) = gprhs(3) + dpdz  +  sig * w

     if( kfl_timei_nsi == 1 ) then
        !
        ! rho*du/dt
        !
        gprhs(1) = gprhs(1)              + rho * dudt
        gprhs(2) = gprhs(2)              + rho * dvdt
        if(ndime==3) gprhs(3) = gprhs(3) + rho * dwdt
     end if

     if( kfl_advec_nsi == 1 ) then
        !
        ! rho*(u.grad)u
        !
        divu     = dudx + dvdy + dwdz
        gprhs(1) =              gprhs(1) + rho * ( u * dudx  +  v * dudy  +  w * dudz )
        gprhs(2) =              gprhs(2) + rho * ( u * dvdx  +  v * dvdy  +  w * dvdz )
        if(ndime==3) gprhs(3) = gprhs(3) + rho * ( u * dwdx  +  v * dwdy  +  w * dwdz )
        
        if(     kfl_convection_type_nsi == NSI_CONVECTION_SKEW ) then
           !
           ! u.grad(u) + 1/2 (div u) u
           !
           gprhs(1) =              gprhs(1) + 0.5_rp * rho * u * divu
           gprhs(2) =              gprhs(2) + 0.5_rp * rho * v * divu
           if(ndime==3) gprhs(3) = gprhs(3) + 0.5_rp * rho * w * divu
           
        else if( kfl_convection_type_nsi == NSI_CONVECTION_CONSERVATIVE ) then
           !
           ! u.grad(u) + (div u) u
           !
           gprhs(1) =              gprhs(1) + rho * u * divu
           gprhs(2) =              gprhs(2) + rho * v * divu
           if(ndime==3) gprhs(3) = gprhs(3) + rho * w * divu
           
        else if( kfl_convection_type_nsi == NSI_CONVECTION_EMAC ) then
           !
           ! u.grad(u) + u.grad(u)^t + (div u) u - 0.5grad(u·u)
           !
           gprhs(1) =              gprhs(1) + rho * ( u * dudx + v * dvdx + w * dwdx ) + rho * u * divu - rho*(dudx*u+dvdx*v+dwdx*w)
           gprhs(2) =              gprhs(2) + rho * ( u * dudy + v * dvdy + w * dwdy ) + rho * v * divu - rho*(dudy*u+dvdy*v+dwdy*w)
           if(ndime==3) gprhs(3) = gprhs(3) + rho * ( u * dudz + v * dvdz + w * dwdz ) + rho * w * divu - rho*(dudz*u+dvdz*v+dwdz*w)
        end if
        
     end if
     !
     ! Viscous term: - div [ 2 mu eps(u) ] 
     !
     if( kfl_visco_nsi == 1  ) then
        !
        ! Laplacian form: - mu*Lapl(u) - dmu/dxj*dui/dxj (Laplacian form)
        !
        gprhs(1) = gprhs(1)              - mu * (d2udx + d2udy + d2udz)  -  dmudx * dudx - dmudy * dudy - dmudz * dudz 
        gprhs(2) = gprhs(2)              - mu * (d2vdx + d2vdy + d2vdz)  -  dmudx * dvdx - dmudy * dvdy - dmudz * dvdz 
        if(ndime==3) gprhs(3) = gprhs(3) - mu * (d2wdx + d2wdy + d2wdz)  -  dmudx * dwdx - dmudy * dwdy - dmudz * dwdz 

        if( fvins_nsi > 0.9_rp  ) then
           !
           ! Divergence form:  Laplacian form - mu*d^2uj/dxidxj - dmu/dxj duj/dxi
           !
           gprhs(1) = gprhs(1)              - mu * ( d2udx   + d2vdxdy + d2wdxdz )  -  dmudx * dudx - dmudy * dvdx - dmudz * dwdx 
           gprhs(2) = gprhs(2)              - mu * ( d2udxdy + d2vdy   + d2wdydz )  -  dmudx * dudy - dmudy * dvdy - dmudz * dwdy  
           if(ndime==3) gprhs(3) = gprhs(3) - mu * ( d2udxdz + d2vdydz + d2wdz   )  -  dmudx * dudz - dmudy * dvdz - dmudz * dwdz  

        end if

        if( fvins_nsi > 1.9_rp  ) then
           !
           ! Complete form: Laplacian + divergence + 2/3 dmu/dxi duj/dxj + 2/3 mu d^2uj/dxidxj
           !
           gprhs(1) = gprhs(1)              + 2.0_rp/3.0_rp * ( dmudx * (dudx+dvdy+dwdz)  +  mu * ( d2udx   + d2vdxdy + d2wdxdz ) )
           gprhs(2) = gprhs(2)              + 2.0_rp/3.0_rp * ( dmudy * (dudx+dvdy+dwdz)  +  mu * ( d2udxdy + d2vdy   + d2wdydz ) )
           if(ndime==3) gprhs(3) = gprhs(3) + 2.0_rp/3.0_rp * ( dmudz * (dudx+dvdy+dwdz)  +  mu * ( d2udxdz + d2vdydz + d2wdz   ) )
        end if

     end if

     if( corio_nsi > zeror ) then
        !
        ! 2*rho*(w x r)
        ! 
        gprhs(1) = gprhs(1)              + 2.0_rp * rho * ( -wz * v + wy * w ) 
        gprhs(2) = gprhs(2)              + 2.0_rp * rho * (  wz * u - wx * w ) 
        if(ndime==3) gprhs(3) = gprhs(3) + 2.0_rp * rho * ( -wy * u + wx * v ) 
     end if

     gprhc = gprhc + (dudx + dvdy + dwdz)
     gprh2 = gprhc

  else if( itask == 3 ) then
     !
     ! Traction
     !
     if( ndime == 2 ) then
        gprhs(1) = gprhs(1) - p * baloc(1) + mu * ( dudx * baloc(1) + dudy * baloc(2) )
        gprhs(2) = gprhs(2) - p * baloc(2) + mu * ( dvdx * baloc(1) + dvdy * baloc(2) )
     else
        gprhs(1) = gprhs(1) - p * baloc(1) + mu * ( dudx * baloc(1) + dudy * baloc(2) + dudz * baloc(3) )
        gprhs(2) = gprhs(2) - p * baloc(2) + mu * ( dvdx * baloc(1) + dvdy * baloc(2) + dvdz * baloc(3) )
        gprhs(3) = gprhs(3) - p * baloc(3) + mu * ( dwdx * baloc(1) + dwdy * baloc(2) + dwdz * baloc(3) )
     end if
     if( fvins_nsi > 0.9_rp ) then
        if( ndime == 2 ) then
           gprhs(1) = gprhs(1) + mu * ( dudx * baloc(1) + dvdx * baloc(2) )
           gprhs(2) = gprhs(2) + mu * ( dudy * baloc(1) + dvdy * baloc(2) )
        else
           gprhs(1) = gprhs(1) + mu * ( dudx * baloc(1) + dvdx * baloc(2) + dwdx * baloc(3) )
           gprhs(2) = gprhs(2) + mu * ( dudy * baloc(1) + dvdy * baloc(2) + dwdy * baloc(3) )
           gprhs(3) = gprhs(3) + mu * ( dudz * baloc(1) + dvdz * baloc(2) + dwdz * baloc(3) )
        end if
     end if
     if( fvins_nsi > 1.9_rp  ) then
        divu = dudx + dvdy + dwdz
        if( ndime == 2 ) then
           gprhs(1) = gprhs(1) - 2.0_rp/3.0_rp * mu * divu * baloc(1)
           gprhs(2) = gprhs(2) - 2.0_rp/3.0_rp * mu * divu * baloc(2)
        else
           gprhs(1) = gprhs(1) - 2.0_rp/3.0_rp * mu * divu * baloc(1)
           gprhs(2) = gprhs(2) - 2.0_rp/3.0_rp * mu * divu * baloc(2)
           gprhs(3) = gprhs(3) - 2.0_rp/3.0_rp * mu * divu * baloc(3)
        end if
     end if

  end if

end subroutine nsi_exacso

!
!  Subroutine to compute the derivatives of the spatial functions.
!   
subroutine spafun(param,x,f,d1f,d2f,d3f)
  use def_kintyp
  implicit none
  real(rp) param(3),x,f,d1f,d2f,d3f
  real(rp) a,b,c,      &
       &   f1,d1f1,d2f1,d3f1, &
       &   f2,d1f2,d2f2,d3f2, &
       &   f3,d1f3,d2f3,d3f3

  a = param(1)
  b = param(2)
  c = param(3)

  f1   = a*(      x*x -  2.0_rp*x*x*x +        x*x*x*x)
  d1f1 = a*(2.0_rp*x   -  6.0_rp*x*x   +  4.0_rp*x*x*x)
  d2f1 = a*(2.0_rp     - 12.0_rp*x     + 12.0_rp*x*x)
  d3f1 = a*(          - 12.0_rp       + 24.0_rp*x)
  f2   = exp(b*x)
  d1f2 = b*f2
  d2f2 = b*d1f2
  d3f2 = b*d2f2
  f3   =        cos(c*x)
  d1f3 =     -c*sin(c*x)
  d2f3 =    c*c*cos(c*x)
  d3f3 = -c*c*c*sin(c*x)
  f    =    f1*f2*f3                ! a*x*x*(1-x)*(1-x)*exp(b*x)*cos(c*x)
  d1f  =    d1f1*f2  *f3     &
       &   +       f1  *d1f2*f3     &
       &   +       f1  *f2  *d1f3
  d2f  =    d2f1*f2  *f3     &
       &   +       f1  *d2f2*f3     &
       &   +       f1  *f2  *d2f3   &
       &   + 2.0_rp*d1f1*d1f2*f3    &
       &   + 2.0_rp*d1f1*f2  *d1f3  &
       &   + 2.0_rp*f1  *d1f2*d1f3
  d3f  =    d3f1*f2  *f3     &
       &   +       f1  *d3f2*f3     &
       &   +       f1  *f2  *d3f3   &
       &   + 3.0_rp*d2f1*d1f2*f3    &
       &   + 3.0_rp*d1f1*d2f2*f3    &
       &   + 3.0_rp*d2f1*f2  *d1f3  &
       &   + 3.0_rp*d1f1*f2  *d2f3  &
       &   + 3.0_rp*f1  *d2f2*d1f3  &
       &   + 3.0_rp*f1  *d1f2*d2f3  &
       &   + 6.0_rp*d1f1*d1f2*d1f3

end subroutine spafun

