!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine sld_elmcla(&
     itask,pgaus,pmate,gpvol,gpcau,gpgi0,gpgdi,gpene,gpstr,gppio,gpdet,gprat,&
     gpigd,gptmo,gpfib,ielem,elcod,pnode,lnods,gpsha,gpdds,&
     gpigd_eps,gpgdi_eps,gpdet_eps,gppio_eps,&
     gpmof)

!-----------------------------------------------------------------------
!****f* Solidz/sld_elmcla
! NAME
!    sld_elmcla
! DESCRIPTION
!    Call built-in or user materials interface
! INPUT
!    PGAUS ... Number of Gauss points
!    PMATE ... Material number
!    GPCAU ... Updated right Cauchy-Green tensor ............... C(n+1)
!    GPGDI ... Updated deformation gradient tensor ............. F(n+1)
!    GPDET ... Updated jacobian ................................ J = det(F(n+1))
!    GPRAT ... Rate of Deformation tensor ...................... Fdot = grad(phidot)
!    GPIGD ... Inverse of updated deformation gradient tensor .. F^{-1}(n+1)
!    IELEM ... Number of the element ........................... I_e
! INPUT/OUTPUT
!    GPENE ... Stored energy function .......................... W(n)/W(n+1)
!    GPSTR ... 2nd Piola-Kirchhoff stress tensor at t(n)/t(n+1)  S(n)/S(n+1)
!    GPPIO ... 1st Piola-Kirchhoff stress tensor at t(n)/t(n+1)  P(n)/P(n+1)
!    GPPIO_EPS ... Perturbed 1st Piola-Kirchhoff stress tensor at t(n)/t(n+1)  P(n)/P(n+1)
!    GPLEP ... Log strain in the {f s n} system
!
! FUTURE INPUT
!    TEMP0 ... Initial temperature ............................. temp(n)
!    TEMP1 ... Updated temperature ............................. temp(n+1)
!    GPGI0 ... Previous deformation gradient tensor ............ F(n)
!    FLAGT ... Flag for tangent moduli calculations ............ 1 if yes; 0 if no
!    FLAGS ... Flag for strain gradient calculations ........... 1 if yes; 0 if no
! FUTURE OUTPUT
!    GPTMO ... Tangent moduli at t(n+1) ........................ dP/dF(n+1)
!    GGSGO ... Deformation gradient gradient tensor ............ dF/dX(n+1)
!
! USES
!    sld_builtin_materials ... Interface subroutine for built-in materials
!    sld_user_materials ...... Interface subroutine for user materials
! USED BY
!    sld_elmope
!***
!-----------------------------------------------------------------------

  use def_kintyp,  only      : ip,rp
  use def_domain,  only      : ndime,mnode
  use mod_sld_atm, only      : kfl_therm_sld
  use mod_sld_atm, only      : ntmp_atm
  use mod_sld_atm, only      : sld_atm_get_temperature_at_GP
  use def_solidz

  implicit none
  integer(ip), intent(in)    :: itask,pgaus,pmate,ielem,pnode,lnods(pnode)
  integer(ip)                :: flagt,flags,idime,jdime,igaus
  real(rp)                   :: gptmp(ntmp_atm,pgaus)
  real(rp),    intent(in)    :: gpcau(ndime,ndime,pgaus)
  real(rp),    intent(in)    :: gpsha(pnode,pgaus)
  real(rp),    intent(in)    :: gpgi0(ndime,ndime,pgaus)
  real(rp),    intent(in)    :: gpgdi(ndime,ndime,pgaus)
  real(rp),    intent(in)    :: gpdet(pgaus),gpmof(pgaus)
  real(rp),    intent(in)    :: gprat(ndime,ndime,pgaus)
  real(rp),    intent(in)    :: gpvol(pgaus)
  real(rp),    intent(in)    :: gpigd(ndime,ndime,pgaus)
  real(rp),    intent(inout) :: gpene(pgaus)
  real(rp),    intent(inout) :: gpstr(ndime,ndime,pgaus,2)
  real(rp),    intent(inout) :: gpfib(ndime,pgaus)
  real(rp),    intent(inout) :: gppio(ndime,ndime,pgaus)
  real(rp),    intent(inout) :: gppio_eps(ndime,ndime,pgaus)
  real(rp),    intent(out)   :: gptmo(ndime,ndime,ndime,ndime,pgaus)
  real(rp)                   :: ggsgo(ndime,ndime)
  real(rp),    intent(in)    :: elcod(ndime,mnode)
  real(rp),    intent(out)   :: gpdds(ndime,ndime,ndime,ndime,pgaus)

  real(rp),    intent(in)        :: gpgdi_eps(ndime,ndime,pgaus) ! Displacement Gradient F for the perturbed state
  real(rp),    intent(in)        :: gpigd_eps(ndime,ndime,pgaus) ! Displacement Gradient Inverse F^{-1} for the perturbed state
  real(rp),    intent(in)        :: gpdet_eps(pgaus)

  !
  ! For now implicit methods and strain gradient calculations are not
  ! implemented
  !
  if( itask == 2 .or. kfl_exacs_sld /= 0 ) then
     flagt = 1_ip
  else
     flagt = 0_ip
  end if
  flags = 0_ip
  !
  ! Thermomechanical coupling
  !
  gptmp(:,:) = 0.0_rp
  if( kfl_therm_sld ) then
     call sld_atm_get_temperature_at_GP(pnode,pgaus,gpsha,lnods,gptmp)
  end if
  !
  ! Built-in materials
  !
  call sld_builtin_materials(pgaus,pmate,gpvol,gptmp,&
       gpgi0,gpgdi,gpigd,gpcau,gpdet,gprat,gppio,gpstr,gpene,&
       flagt,gptmo,flags,ggsgo,ielem,gpfib,elcod,pnode,lnods,gpsha,gpdds,&
       gpigd_eps,gpgdi_eps,gpdet_eps,gppio_eps,&
       gpmof)  
  !
  ! Stock F1 in gpgdi_sld and detF in dedef_sld
  !
  do igaus = 1,pgaus
     do idime = 1,ndime
        do jdime = 1,ndime
           gpgdi_sld(ielem) % a(idime,jdime,igaus) = gpgdi(idime,jdime,igaus)
        end do
     end do
     dedef_sld(ielem) % a(igaus) = gpdet(igaus)
  end do
  !
  ! Stock gpene in enden_sld
  !
  do igaus = 1, pgaus
     enden_sld(ielem) % a(igaus) = gpene(igaus)
  end do

end subroutine sld_elmcla
