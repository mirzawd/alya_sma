!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Magnet
!> @{
!> @file    mag_updbcs.f90
!> @date    29/01/2018
!> @author  Guillaume Houzeaux
!> @brief   Update boundary conditions
!> @details This routine updates the magnetic field boundary conditions:
!>          1. Before a time step begins
!>          2. Before a global iteration begins
!>          3. Before an inner iteration begins
!>          4. Force boundary conditions in UNKNO
!> @} 
!-----------------------------------------------------------------------
subroutine mag_updbcs(itask)

  use def_magnet
  use def_master
  use def_domain
  use def_kermod, only: ndivi

  implicit none

  integer(ip), intent(in) :: itask

  integer(ip) :: iedge

  select case( itask )

  case( ITASK_BEGSTE )
    !
    ! Dirichlet boundary conditions and self field
    !
    if (.not. kfl_self_mag .or. (kfl_self_mag .and. kfl_prev_mag)) then
      !
      ! Self-Field calculation: Biot-Savart
      !
      if (kfl_self_mag) call mag_selfie()
      !
      ! External Applied Field
      !
      call mag_dirich()
      !
      ! Initialize Dirichlet boundary conditions
      !
      call mag_inidir()
      !
      ! Update Dirichlet boundary conditions
      !   
      call mag_upddir()
      ! 
      !**DIRBC**
      do iedge = 1_ip, meshe(ndivi) % nedge
        if (kfl_fixno_mag(1_ip, iedge) == 1_ip) then
          bvess_mag(1_ip, iedge, 1_ip) = He_mag(iedge)
        end if
      end do
      !
      ! Impose Dirichlet value on He_mag
      !
      do iedge = 1_ip, meshe(ndivi) % nedge
        if (kfl_fixno_mag(1_ip, iedge) == 1_ip) then
          He_mag(iedge) = bvess_mag(1_ip, iedge, 1_ip) 
        end if
      end do
      !**DIRBC**
      ! 
    end if
    
  case( ITASK_DOITER )
    !---------------------------------------------------------
    ! Self field calculation enabled
    !---------------------------------------------------------
    if (kfl_self_mag .and. .not. kfl_prev_mag) then
      !
      ! Self-Field calculation: Biot-Savart
      !
      call mag_selfie()
      !
      ! External Applied Field
      !
      call mag_dirich()
      !
      ! Initialize Dirichlet boundary conditions
      !
      call mag_inidir()
      !
      ! Update Dirichlet boundary conditions
      !
      call mag_upddir()
      !
      !**DIRBC**
      !
      do iedge = 1_ip, meshe(ndivi) % nedge
        if (kfl_fixno_mag(1_ip, iedge) == 1_ip) then
          bvess_mag(1_ip, iedge, 1_ip) = He_mag(iedge)
        end if
      end do
      !
      ! Impose Dirichlet value on He_mag
      !
      do iedge = 1_ip, meshe(ndivi) % nedge
        if (kfl_fixno_mag(1_ip, iedge) == 1_ip) then
          He_mag(iedge) = bvess_mag(1_ip, iedge, 1_ip)
        end if
      end do
      !
      !**DIRBC**
      !
    end if

  end select

end subroutine mag_updbcs
