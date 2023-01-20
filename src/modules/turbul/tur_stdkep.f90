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
!> @file    tur_stdkep.f90
!> @author  Guillaume Houzeaux
!> @brief   Standard k-eps model
!> @details Compute coefficient of the equation of standard k-eps model           \n
!!                                                                                \n
!!    - \f$ \rho Dk/Dt = P - rho* \varepsilon1 
!!      + div[ (\mu+\mu_t/\sigma_k) \nabla k  ] \f$                               \n
!!                                                                                \n 
!!    - \f$ \rho D\varepsilon/Dt = C_{\varepsilon1} \varepsilon/k P
!!       - \rho C_{\varepsilon2} \varepsilon^2 / k 
!!       + div[ (\mu+\mu_t/\sigma_\varepsilon) \nabla \varepsilon  ]\f$           \n
!!                                                                                \n
!!    - \f$ \mu_t = \rho C_{\mu} k^2/e  \f$                                       \n
!!                                                                                \n
!!    For Non-Neutral atmospheric boundary layer:                                 \n
!!    - \f$ C_{\varepsilon 1} <= [ C_{\varepsilon 1} 
!!      + (C_{\varepsilon 2}-C_{\varepsilon 1}) l_m/l_{max} ]  \f$                \n
!!    - \f$ l_m=C_\mu^{3/4}k^{3/2}/\varepsilon   \f$                           \n
!!    - if  (\f$ z/L < 0  \f$)  then
!!       \f$ \phi_m = (1 - 16  z/L )^{-1/4} \f$
!!      else
!!       \f$ \phi_m = 1 + 4.7  z/L \f$
!!      endif
!!
!> @} 
!-----------------------------------------------------------------------
subroutine tur_stdkep(&
     pnode,gptur,gpden,gpvis,gpmut,gpgra,&
     gpcar,gppr2,gprea,gpdif,gprhs,&
     kfl_kxmod_tur, eta, c_mu, gpcan ,sreac, &
     pmate, grunk, conve, gplmax)

  use def_kintyp, only     :  ip,rp 
  use def_domain, only     :  ndime
  use def_turbul, only     :  param_tur,nturb_tur,iunkn_tur, kfl_cotem_tur, ldiss_material_tur, &
       kfl_lmaxi_tur
  use def_turbul, only     :  inv_l_max, dtinv_tur
  use def_kermod, only     :  cmu_st
  use mod_ker_regularization, only : regul_k, regul_e, dregularization, d2regularization, kfl_regularization, kfl_second
  
  implicit none
  integer(ip), intent(in)  :: pnode                                !< Number of element nodes
  integer(ip), intent(in)  :: kfl_kxmod_tur                            !< Modification of ke model
  real(rp),    intent(in)  :: gptur(nturb_tur,3)                   !< Turb. varibales at Gauss points
  real(rp),    intent(inout)  :: gpgra                                !< Gauss point Gravity production GPGRA=beta*mut/Prt*g.grad(T)
  real(rp),    intent(in)  :: gpden                                !< Gauss point density
  real(rp),    intent(in)  :: gpvis                                !< Gauss point viscosity
  real(rp),    intent(inout)  :: gpmut(2)                             !< Gauss point mut
  real(rp),    intent(in)  :: gppr2                                !< Gauss point Shear production GPPR2 =      (dui/dxj+duj/dxi)*dui/dxj
  real(rp),    intent(in)  :: gpcar(ndime,pnode)                   !< Gauss point Shape function Cartesian derivatives
  real(rp),    intent(out) :: gpdif                                !< Gauss point diffusion: k
  real(rp),    intent(out) :: gprea                                !< Gauss point reaction: r
  real(rp),    intent(in) ::  gpcan                                !< Gauss point canopy cd*LAD*|u|
  real(rp),    intent(out) :: sreac                                !< Gauss point canopy cd*LAD*|u|
  real(rp),    intent(inout) :: gprhs                              !< Gauss point RHS: f
  real(rp),    intent(in)  :: eta, c_mu                            !< Eta parameter for RNG model
  real(rp),    intent(in)  :: grunk(ndime, 2)
  real(rp),    intent(inout) :: conve(ndime)                       !< conective term
  integer(ip), intent(in)  :: pmate
  real(rp),    intent(in)  :: gplmax 
  real(rp)                 :: gpkin,gpeps,Ce1,Ce2,Cmu, DCe1, Ce3, C1p, Ce4
  real(rp)                 :: lm, gpepo, gpkio, alpha, onerp, react,reac2, cmu0, reac3, canopy
  real(rp)                 :: reguk, regue, drdk, drde,drdko, drdeo, reguko, regueo, difco
  real(rp)                 :: d2rdk, d2rde, rfunc, dfdk, dfde, df2dk, df2de, df3de
  real(rp)                 :: facto,inv_lmax, epsam, keyam
  real(rp)                 :: gppro                               !< Gauss point Shear production GPPRO = mu_t*(dui/dxj+duj/dxi)*dui/dxj
 

  ! atmospheric values
  if (kfl_regularization.or.kfl_cotem_tur.ge.1_ip) then
     epsam = 7.208e-8_rp
     keyam = 1.0e-4_rp
  else
     epsam = 0.0_rp
     keyam = 1.0_rp
  end if
  !
  ! Definitions
  !
  gpkin = gptur(1,1)
  gpkio = gptur(1,2) ! old tke 
  gpeps = gptur(2,1)
  gpepo = gptur(2,2) ! old eps
  Ce1   = param_tur(3)
  Ce2   = param_tur(4)
  Ce3   = param_tur(5)
  Cmu   = param_tur(6) !constant cmu, from data files
  cmu0  = cmu ! for realizable
  DCe1 = 0.0_rp
  onerp = 1.0_rp
  if (kfl_kxmod_tur==2.or.kfl_kxmod_tur==3 ) cmu = c_mu ! gets modified cmu
  Ce4 = 0.37_rp
  ! inverse of maximum mixing length
  if (kfl_lmaxi_tur==0) then ! uniform max mixing length
     inv_lmax =  inv_l_max
  else                       ! Mellor yamada, non uniform lmax
     if (gplmax.gt.1.0_rp) then
        inv_lmax =  1.0_rp/gplmax
     else
        inv_lmax =  1.0_rp/1.0_rp ! minimum set to one meter
     end if
     
  end if
  if (kfl_regularization) then
     reguk  = regul_k(gpkin)
     reguko   = regul_k(gpkio)
     regue    = regul_e(gpeps)
     regueo   = regul_e(gpepo)
     gpmut(1) = cmu*gpden*reguk*reguk/regue    !+ gpvis
     gpmut(2) = cmu*gpden*reguko*reguko/regueo !+ gpvis
  else ! without regularization
     ! turbulent viscosity 
     gpmut(1) = max(cmu*gpkin*gpkin/gpeps*gpden, gpvis) !last iteration
     gpmut(2) = max(cmu*gpkio*gpkio/gpepo*gpden, gpvis) !frozen along inner iters
  end if
  gppro    = gpmut(2) * gppr2
  gpgra    = gpmut(2)*  gpgra
  !
  ! Maximum length scale limitation
  ! Ce1 <= [ Ce1 + (Ce2-Ce1) lm/lmax ]      \n
  !
  difco = 1.0_rp ! diffusion coefficient
    
     if( iunkn_tur == 1 ) then  ! TKE equation
        if (.not.kfl_regularization ) then
           gprea = max(2.0_rp*gpden*gpden*Cmu*gpkin/gpmut(2), 0.0_rp)
           sreac = max(2.0_rp*gpden*gpden*Cmu*gpkin/gpmut(2), 0.0_rp)
           if (gprea.lt.epsilon(1.0_rp)) onerp =0.0_rp
           gprhs = gprhs + max(gppro  +  gpgra,0.0_rp) + onerp*1.0_rp*gpden*gpeps + gpden*epsam    
           
        else if (.not.kfl_second) then ! regularization
           drdk  = dregularization (gpkin)
           d2rdk = d2regularization(gpkin)
           rfunc = reguk * reguk /drdk !           
           dfdk = (2.0_rp - reguk*d2rdk/(drdk*drdk))*reguk 
           !       r2fun = -1/drdk 
           df2dk = d2rdk/(drdk*drdk)
           react = gpden*gpden*Cmu/gpmut(2)* dfdk 
           reac2 = ( max(gppro  +  gpgra,0.0_rp) +  gpden*(epsam))   *df2dk
           gprea = react + reac2
           sreac = react + reac2
           gprhs = (  max(gppro  +  gpgra,0.0_rp)+  gpden*(epsam))/dregularization(gpkin)  &
          !      -  gpden*regue/dregularization(gpkio) &
                +  reac2*gpkin &
                +  react*gpkin  - gpden*gpden*Cmu/gpmut(2)*rfunc 
                
           !           conve(1:ndime) = conve (1:ndime) - 2.0_rp* gpmut(2)*delta*delta/(param_tur(1)*sqrtk*sqrtk*reguko )*grunk(1:ndime,2)/gpden
           conve(1:ndime) = conve (1:ndime) - gpmut(2)*d2regularization(gpkio)/(param_tur(1)*dregularization(gpkio))*grunk(1:ndime,2)/gpden
           
        else if (kfl_second) then    ! regularization second form
           drdk  = dregularization (gpkin)
           drdko  = dregularization (gpkio)
           conve(1:ndime) = conve (1:ndime) * drdko
           rfunc  = reguk * reguk 
           dfdk   = 2.0_rp*reguk*drdk
           react = gpden*gpden*Cmu/gpmut(2)* dfdk 
           gprea = react 
           sreac = react
           gprhs = max(gppro  +  gpgra,0.0_rp) - gpden*gpden*Cmu/gpmut(2)*rfunc &
                + react*gpkin &
                + gpden* dtinv_tur*(drdk*gpkin  + regul_k(gptur(1,3)) -reguk)
!                + gpden* dtinv_tur*(drdk*gptur(1,3) )
           
          difco = drdko
        end if
     else ! EPSILON EQUATION 
        if (.not.kfl_regularization ) then
           if (kfl_kxmod_tur==2) then  !Realizable model           
              Ce1 = max(0.43_rp, eta/(eta+5.0_rp))
              facto = sqrt(cmu_st/0.09_rp)
              Ce1 = facto*max(0.43_rp, eta/(eta+5.0_rp/facto))
              

              ! mixing length
              lm = (Cmu*gpkin*gpkin)**(0.75_rp)/gpepo    
              if (lm.lt.0.0_rp) lm=1.0_rp/inv_lmax
              
              Ce1 = Ce1 + (Ce2*sqrt(Cmu_st) - Ce1)*lm*inv_lmax         
              

              react = Ce2 * gpden * gpeps / (gpkin + sqrt(gpvis*gpepo))
              gprea = 2.0_rp*react
              sreac = 2.0_rp*react
              gprhs = gprhs    &
                   +  react * gpeps +  gpden*gpepo* sqrt(gppr2)*Ce1
              
           else   !normal model
              ! mixing length
              !           lm = min ((Cmu*gpkin*gpkin)**(0.75_rp)/gpepo, 1000.0/inv_l_max)
              lm = (Cmu*gpkin*gpkin)**(0.75_rp)/gpepo
              if (lm.lt.0.0_rp) lm=1.0_rp/inv_lmax
              ! l_max stores the inverse of the maximum (l_max=0 when no mixing model)
              DCe1  =  (Ce2-Ce1) * lm*inv_lmax    !  * min( lm/l_max, 1.0_rp)
              if (kfl_kxmod_tur==1) & ! RNG model
                   DCe1  = - eta*(1.0_rp-eta/4.377_rp)/(1.0_rp+0.012_rp*eta*eta*eta)      
              
              C1p  = Ce1 + DCe1  ! C1 modifidied by Apsley and Castro
              
              react = max(Ce2 * gpden * gpeps / gpkin, 0.0_rp) ! reactive term (linearized)
              gprea = 2.0_rp*react + 12.0_rp*(Ce1-Ce2)*sqrt(Cmu)*gpcan 
              sreac = 2.0_rp*react + 12.0_rp*(Ce1-Ce2)*sqrt(Cmu)*gpcan ! reaction in perturbation test function
              
              if (kfl_cotem_tur==2) then ! Thermal coupling using SOGACHEV's
                 if (gpgra.lt.0.0_rp) then  ! stable
                    alpha = 1.0_rp - lm*inv_lmax   
                    gpgra =  max(0.0d0,gpgra + gppro) - gppro
                 else             
                    alpha = 1.0_rp - ( 1.0_rp +(Ce2 -1.0_rp)/(Ce2-Ce1) )*lm*inv_lmax              
                 end if
                 Ce3 = 1.0_rp + (Ce1-Ce2)*alpha
              end if
              gprhs = gprhs +   max(C1p*Cmu*gpden*gpkin* gppr2 + Ce3*gpepo/gpkin*gpgra, 0.0_rp)&
                   + react*gpeps + gpden *Ce2*epsam*epsam/keyam
           end if
        else if(kfl_kxmod_tur.ne.2) then !  regularization but not realizable!!!!!
           
           ! mixing length
           lm = (Cmu*reguk*reguk)**(0.75_rp)/regueo
           
           ! l_max stores the inverse of the maximum (l_max=0 when no mixing model)
           DCe1  =  (Ce2-Ce1) * lm*inv_lmax    !  * min( lm/l_max, 1.0_rp)
           
           C1p  = Ce1 + DCe1  ! C1 modifidied by Apsley and Castro
           if (.not.kfl_second) then ! 
              drde  = dregularization(gpeps)
              d2rde = d2regularization(gpeps)
              !   r2fun = -1/drde
              df2de = d2rde/(drde*drde)
                    
              rfunc = regue * regue /drde !           
              dfde = (2.0_rp - regue*df2de)*regue
              
              df3de = 1.0_rp - regue*df2de ! linearization for canopy
              if (kfl_cotem_tur==2) then ! Thermal coupling using SOGACHEV's
                 if (gpgra.lt.0.0_rp) then  ! stable
                    alpha = 1.0_rp - lm*inv_lmax   
                    gpgra =  max(0.0d0,gpgra + gppro) - gppro
                 else             
                    alpha = 1.0_rp - ( 1.0_rp +(Ce2 -1.0_rp)/(Ce2-Ce1) )*lm*inv_lmax              
                 end if
                 Ce3 = 1.0_rp + (Ce1-Ce2)*alpha
              end if
              
              
              react = Ce2 * gpden / reguk* dfde 
              reac2 = ( max(C1p*Cmu*gpden*reguk* gppr2  + Ce3*regueo/reguk*gpgra, 0.0_rp) + Ce2*gpden *epsam*epsam/keyam)*df2de
              canopy= 12.0_rp*(Ce1-Ce2)*sqrt(Cmu)*gpcan
              reac3 = canopy*df3de
              ! reactive term (linearized)
              gprea = react + reac2 + reac3 
              sreac = react + reac2 + reac3
              gprhs = ( max(C1p*Cmu*gpden*reguk* gppr2  + Ce3*regueo/reguk*gpgra, 0.0_rp) + Ce2*gpden*epsam*epsam/keyam)/drde &
                   - Ce2 * gpden / reguk*rfunc &
                   - canopy*regue/drde &
                   + (reac2+react+reac3)*gpeps 
                  
             
              
              conve(1:ndime) = conve (1:ndime) &
                   - gpmut(2)*d2regularization(gpepo)/(param_tur(2)*dregularization(gpepo))*grunk(1:ndime,2)/gpden
           else if (kfl_second) then
              drde  = dregularization(gpeps)
              drdeo = dregularization(gpepo)
              rfunc = regue * regue         
              dfde =  2.0_rp* regue*drde
              react = Ce2 * gpden / reguk* dfde 
              canopy = 12.0_rp*(Ce1-Ce2)*sqrt(Cmu)*gpcan
              reac3 =  canopy*drde
              gprea = react + reac3
              sreac = react + reac3
              if (kfl_cotem_tur==2) then ! Thermal coupling using SOGACHEV's
                 if (gpgra.lt.0.0_rp) then  ! stable
                    alpha = 1.0_rp - lm*inv_lmax   
                    gpgra =  max(0.0d0,gpgra + gppro) - gppro
                 else             
                    alpha = 1.0_rp - ( 1.0_rp +(Ce2 -1.0_rp)/(Ce2-Ce1) )*lm*inv_lmax              
                 end if
                 Ce3 = 1.0_rp + (Ce1-Ce2)*alpha
              end if
              
              gprhs = max(C1p*Cmu*gpden*reguk* gppr2  + Ce3*regueo/reguk*gpgra, 0.0_rp) &
                   - Ce2 * gpden / reguk*rfunc &
                   +  (react + reac3) *gpeps - canopy*regue & 
                   + gpden*dtinv_tur*(drde*gpeps + regul_e(gptur(2,3)) - regue) &
                   + gpden*Ce2*epsam*epsam/keyam
!                 
 
              conve(1:ndime) = conve (1:ndime)* drdeo 
              difco = drdeo
           end if
        end if
           
        ! Disc dissipation model (dissipation equation)
        if (ldiss_material_tur(pmate)==1) & 
             gprhs = gprhs + max(Ce4*gppro*gppro/(gpden*gpkio), 0.0_rp)
        
     end if
     
  !
  ! GPDIF: diffusion and its gradient 
  !
     !  call tur_elmdif(pnode,gpvis,gpmut(iunkn_tur)*difco,gpcar,eledd, gpdif, gpgrd)
     gpdif = gpvis + gpmut(iunkn_tur)/param_tur(iunkn_tur)*difco
     
end subroutine tur_stdkep


   
