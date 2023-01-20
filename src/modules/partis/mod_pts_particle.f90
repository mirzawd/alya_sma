!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Partis
!> @{
!> @file    mod_pts_particle.f90
!> @author  houzeaux
!> @date    2019-03-18
!> @brief   Particles
!> @details Quantities related to particles
!-----------------------------------------------------------------------

module mod_pts_particle

  use def_partis
  use mod_physics
  implicit none

  private

  public :: pts_particle_diameter
  public :: pts_particle_mass
  
contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-03-18
  !> @brief   Particle diameter
  !> @details Compute particle diameter according to model used
  !> 
  !-----------------------------------------------------------------------

  real(rp) function pts_particle_diameter(itype,ilagr,rho)

    integer(ip),       intent(in) :: itype
    integer(ip),       intent(in) :: ilagr
    real(rp),optional, intent(in) :: rho

    if( parttyp(itype) % kfl_therm == 0 ) then
       pts_particle_diameter = parttyp(itype) % diame
    else
        if (present(rho)) then
            pts_particle_diameter = max(0.0_rp,physics_sphere_diameter(lagrtyp(ilagr) % mass_k, rho))
        else
            call physics_set_liquid_temperature( parttyp(itype) % liq , lagrtyp(ilagr) % tempe_k)
            pts_particle_diameter = max(0.0_rp,physics_sphere_diameter(lagrtyp(ilagr) % mass_k, parttyp(itype) % liq % rho))
        endif
    end if

  end function pts_particle_diameter
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-03-18
  !> @brief   Particle diameter
  !> @details Compute particle mass according to model used
  !> 
  !-----------------------------------------------------------------------

  real(rp) function pts_particle_mass(itype,ilagr)

    integer(ip), intent(in) :: itype
    integer(ip), intent(in) :: ilagr

    if( parttyp(itype) % kfl_therm == 0 ) then
       pts_particle_mass = physics_sphere_mass(parttyp(itype) % diame,parttyp(itype) % denpa)
    else
       pts_particle_mass = lagrtyp(ilagr) % mass_k 
    end if

  end function pts_particle_mass
  
end module mod_pts_particle
!> @}
