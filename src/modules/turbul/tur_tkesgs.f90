!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup TurbulModel
!> @ingroup    Turbul
!> @{
!> @file    tur_tkesgs.f90
!> @author  Matias Avila
!> @brief   tke for nu_t sgs
!> @details Compute coefficient of the TKE equation       \n
!!                                                                                \n
!!    - \f$ \rho Dk/Dt = P - rho* \varepsilon1 
!!      + div[ (\mu+\mu_t/\sigma_k) \nabla k  ] \f$                               \n
!!                                                                                \n 
!!
!> @} 
!-----------------------------------------------------------------------
subroutine tur_tkesgs(&
     gptur,gpden,gpvis,gpmut,gpgra,&
     gppr2, grnor_tur,hleng, gpwal, gprea,gpdif,gprhs,&
     gpcan ,sreac &
     )

  use def_kintyp, only     :  ip,rp 
  use def_domain, only     :  ndime
  use def_turbul, only     :  kfl_cotem_tur , prtur_tur
  
  implicit none
  real(rp),    intent(in)  :: gptur(3)                             !< Turb. varibales at Gauss points
  real(rp),    intent(inout)  :: gpgra                                !< Gauss point Gravity production GPGRA=beta*mut/Prt*g.grad(T)
  real(rp),    intent(in)  :: gpden                                !< Gauss point density
  real(rp),    intent(in)  :: gpvis                                !< Gauss point viscosity
  real(rp),    intent(inout)  :: gpmut(2)                             !< Gauss point mut
  real(rp),    intent(in)  :: gppr2                                !< Gauss point Shear production GPPR2 =      (dui/dxj+duj/dxi)*dui/dxj
  real(rp),    intent(out) :: gpdif                                !< Gauss point diffusion: k
  real(rp),    intent(out) :: gprea                                !< Gauss point reaction: r
  real(rp),    intent(in) ::  gpcan                                !< Gauss point canopy cd*LAD*|u|
  real(rp),    intent(out) :: sreac                                !< Gauss point canopy cd*LAD*|u|
  real(rp),    intent(inout) :: gprhs                              !< Gauss point RHS: f
  real(rp),    intent(in)    :: grnor_tur                          ! gravity norm
  real(rp),    intent(in)    :: hleng(ndime)                       ! gravity norm
  real(rp),    intent(in)    :: gpwal                              ! wall distance 
 
  real(rp)                 :: gpkin, gpkio
  real(rp)                 :: lmixi, lmixio, chale
  real(rp)                 :: gppro, sigka, react, prtur                             !< Gauss point Shear production GPPRO = mu_t*(dui/dxj+duj/dxi)*dui/dxj

  !
  ! Definitions
  !
  sigka = 1.0_rp
  gpkin = gptur(1) !max(gptur(1),0.0_rp) ! current tke
  gpkio = gptur(2) ! max(gptur(2),0.0_rp) ! old tke

! characteristic length
  chale = (hleng(1)*hleng(2)*hleng(3))**(0.3333333_rp)
  ! mixing length
  lmixi  = chale
  lmixio = lmixi
  gpgra  = gpgra*prtur_tur ! recover grad(temp) dot gravi /teref
  if (kfl_cotem_tur==1_ip.and.gpgra.lt.-1.0e-8_rp.and.gpwal.gt.70.0_rp) then  !stable stratification
     lmixi = max(  min(0.76_rp*sqrt(-gpkin/gpgra), chale), chale*0.01_rp)
     lmixio= max(  min(0.76_rp*sqrt(-gpkio/gpgra), chale), chale*0.01_rp)
  end if
  ! turbulent viscosity 
  gpmut(1) = max(0.1_rp*gpden*sqrt(gpkin)*lmixi , 0.0_rp) ! current
  gpmut(2) = max(0.1_rp*gpden*sqrt(gpkio)*lmixio, 0.0_rp) ! frozen

  ! thermal and mechanical production
  gppro    = gpmut(2)* gppr2  ! mechanical production of turbulence
  if (kfl_cotem_tur==1_ip) then ! thermal production
     ! turbulent prandtl number  (1/3 < Prt <= 1)
     prtur = 1.0_rp/(1.0_rp + 2.0_rp*lmixio/chale )
     gpgra    = gpmut(2)* gpgra/prtur  ! thermal production of turbulence
  end if

  ! reactive term, + rho*e
  react = max(gpden*(0.19_rp+ 0.51_rp*lmixio/chale)/lmixio*sqrt(gpkin), 0.0_rp)

  gprea = 1.5_rp*react + 8.0_rp/3.0_rp*gpcan

  sreac = gprea 
 
  gprhs = gprhs + max(gppro  +  gpgra, 0.0_rp) + 0.5_rp*react*gpkin    
  
  gpdif = gpvis + gpmut(2)/sigka
     
  
end subroutine tur_tkesgs


   
