!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



module mod_tem_turbul

contains
subroutine tem_turbul(&
     pnode,pgaus,igaui,igauf,gpcon,gpsph, &
     gpdif,gpgrd,densi,gptur,gpgrt, gptke, hleng)   
  !-----------------------------------------------------------------------
  !****f* Temper/tem_turbul
  ! NAME 
  !    tem_turbul
  ! DESCRIPTION
  !  Coupling with TURBUL
  !    Compute effective diffusion coefficient rho*(D+D_t)
  ! USES
  ! USED BY
  !    tem_elmope_new
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only      : ip,rp
  use def_domain, only      : ndime
  use def_kermod, only      : turmu_ker, kfl_prtur_abl_ker
!  use def_kermod, only      : gravi, grnor
  use def_temper, only      : prtur_tem,kfl_grdif_tem, &
                              kfl_regim_tem, kfl_rhs_scal_tem
  use def_master, only      : div_enthalpy_transport,kfl_coupl,ID_TEMPER,ID_CHEMIC, &
                              kfl_htran

  implicit none
  integer(ip), intent(in)   :: pnode,pgaus,igaui,igauf
  real(rp),    intent(in)   :: gpcon(pgaus)
  real(rp),    intent(in)   :: gpsph(pgaus)
  real(rp),    intent(in)   :: densi(pgaus)
  real(rp),    intent(inout):: gptur(pgaus)
  real(rp),    intent(inout):: gpgrd(ndime,pgaus)
  real(rp),    intent(out)  :: gpdif(pgaus)
  real(rp),    intent(in), optional :: gpgrt(ndime, pgaus)
  real(rp),    intent(in), optional :: gptke(pgaus)
  real(rp),    intent(in), optional :: hleng(ndime)
  integer(ip)               :: igaus
  real(rp)                  :: prtur, length, facto, teref
!  real(rp)                  :: grtez
  !
  ! ENTHALPY EQUATION
  ! 
  if (kfl_regim_tem == 4) then
    
    !
    ! Laminar contribution (rho*D)
    !
    if (kfl_coupl(ID_TEMPER,ID_CHEMIC) >= 1 .and. associated(div_enthalpy_transport) .and. kfl_htran == 1) then
       gpdif(igaui:igauf) = 0.0_rp   ! In this case div_enthalpy_transport is computed in Chemic and no additional diffusion term is required for the energy eq.
    else
       gpdif(igaui:igauf) = gpcon(igaui:igauf) / gpsph(igaui:igauf)
    end if
    ! 
    ! Turbulent contribution (rho*D_t) for RANS & LES
    !
    if( turmu_ker % kfl_exist /= 0_ip ) then
      !
      ! Compute mu_t Compute mu_t -- now all are multiplied by densi outside
      !
      gptur(igaui:igauf) = gptur(igaui:igauf) * densi(igaui:igauf)
      !
      ! Effective diffusion coefficient: rho*(D+D_t)
      !
      gpdif(igaui:igauf) = gpdif(igaui:igauf) + gptur(igaui:igauf) / prtur_tem
 
    end if

  !
  ! TEMPERATURE EQUATION 
  !
  else
    !
    ! Laminar contribution (k)
    !
    gpdif(igaui:igauf) = gpcon(igaui:igauf)
    ! 
    ! Turbulent contribution (k_t) for RANS & LES
    !
    if( turmu_ker % kfl_exist /= 0_ip ) then
      !
      ! Compute mu_t -- now all are multiplied by densi outside
      !
      gptur(igaui:igauf) = gptur(igaui:igauf) * densi(igaui:igauf)
      !
      ! Effective diffusion coefficient: (k + k_t)
      !
      if (kfl_prtur_abl_ker==1_ip.and.present(gpgrt)) then ! turbulent prandtl depending on stable stratification 
         teref = 298.0  ! used for thermal expansion coefficient alpha = 1/Teref
         do igaus = igaui,igauf

            ! facto = mixlength/hleng
            if (gptke(igaus).lt.1.0e-10_rp) then
               facto = 1.0_rp
            else
               length =  (hleng(1)*hleng(2)*hleng(3))**(0.3333333_rp)          
               facto = 10.0_rp*gptur(igaus)/(densi(igaus)*sqrt(gptke(igaus))*length)
            end if
!!$            grtez = dot_product( gpgrt(1:ndime, igaus),gravi(1:ndime))
!!$            grtez=  max( -1.0_rp*grtez, 1.0e-9_rp)*grnor
!!$            length =  (hleng(1)*hleng(2)*hleng(3))**(0.3333333_rp)
!!$            facto =   min(1.0_rp,0.76_rp*sqrt(teref*gptke(igaus)/grtez)/length)
            prtur = 1.0_rp/(1.0_rp + 2.0_rp*facto ) ! between 1/3 (neutral) and 1 (stable stratification)
            gpdif(igaus) = gpdif(igaus) + gpsph(igaus)*gptur(igaus)/prtur
         end do
      else ! constant prandtl number
         gpdif(igaui:igauf) = gpdif(igaui:igauf) + gpsph(igaui:igauf)*gptur(igaui:igauf)/prtur_tem
      end if
    end if


  end if

  if ( kfl_rhs_scal_tem > 0 ) gpdif(igaui:igauf) = gpdif(igaui:igauf) / densi(igaui:igauf)  

  if(kfl_grdif_tem/=0) then
     !
     ! GPGRD=grad(k+kt)
     !
     if( turmu_ker % kfl_exist /= 0_ip ) call runend('TEM_TURBUL: GRADIENT OF TURBULENT VISCOSITY NOT CODED')

  endif

end subroutine tem_turbul
end module mod_tem_turbul
