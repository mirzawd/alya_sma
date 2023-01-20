!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine chm_elmprc_flamLet(&
  iclas,pgaus,gpden,gpmas,gphco,gpsph,gptur,gpdis,gpprd, &
  gpdif,gprhs)

  !-----------------------------------------------------------------------
  !****f* chemic/chm_elmprc
  ! NAME
  !    chm_elmprc
  ! DESCRIPTION
  !    Compute terms for each species ADR equation
  ! USES
  ! USED BY
  !    chm_element_operations
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only      :  ip,rp
  use def_chemic, only      :  nclas_chm
  use def_chemic, only      :  diffu_chm
  use def_kermod, only      :  turmu_ker

  implicit none

  integer(ip),  intent(in)  :: iclas
  integer(ip),  intent(in)  :: pgaus
  real(rp),     intent(in)  :: gpden(pgaus)
  real(rp),     intent(in)  :: gpmas(pgaus,nclas_chm)                ! Species source term (reaction rate)
  real(rp),     intent(in)  :: gphco(pgaus)                          ! heat conductivity
  real(rp),     intent(in)  :: gpsph(pgaus)                          ! specific heat capacity
  real(rp),     intent(in)  :: gptur(pgaus)                          ! turbulent viscosity
  real(rp),     intent(in)  :: gpdis(pgaus,nclas_chm)                ! Dissipation rate term in the Flamelet model
  real(rp),     intent(in)  :: gpprd(pgaus,nclas_chm)                ! Production term in the Flamelet model

  real(rp),     intent(out) :: gpdif(pgaus,nclas_chm)
  real(rp),     intent(out) :: gprhs(pgaus)

  integer(ip)               :: igaus

  do igaus = 1,pgaus
     !
     ! RHS terms: Source + Production term + Transport reaction rate fluctuations + Dissipation term
     !
     gprhs(igaus) = gprhs(igaus) + gpmas(igaus,iclas) + gpprd(igaus,iclas) - gpdis(igaus,iclas)

     !
     ! Diffusion coefficient
     !
     gpdif(igaus,iclas) = gphco(igaus) / gpsph(igaus)

     !
     ! Adding turbulent part
     !
     if(turmu_ker % kfl_exist /= 0_ip) then
        gpdif(igaus,iclas) = gpdif(igaus,iclas) + gptur(igaus) * gpden(igaus) / diffu_chm(1,1)
     end if

 enddo

end subroutine chm_elmprc_flamLet
