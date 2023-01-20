!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Magnet
!> @{
!> @file    mag_inibcs.f90
!> @author  Guillaume Houzeaux
!> @brief   Impose boundary conditions
!> @details Impose boundary conditions
!> @} 
!-----------------------------------------------------------------------
subroutine mag_inibcs(itask)

  use def_parame
  use def_inpout
  use def_master
  use def_kermod
  use def_domain
  use def_magnet 
  use mod_boundary_conditions, only : boundary_conditions_impose_edge_codes 
  !use mod_boundary_conditions, only : boundary_conditions_impose_boundary_codes 

  implicit none

  integer(ip), intent(in) :: &
    itask

  integer(ip) :: &
    iedge,  &
    jedge

  logical(lg) :: &
    direnabled = .false.

  select case(itask)

    case (1_ip)

      !-------------------------------------------------------------
      !
      ! Allocate memory for the vectors needed to define the BC's 
      !
      !-------------------------------------------------------------

      call mag_membcs() 

      if (INOTEMPTY) then
        call boundary_conditions_impose_edge_codes(tecod_mag(1), kfl_fixno_mag, bvess_mag)
        !call boundary_conditions_impose_boundary_codes(tbcod_mag(1), kfl_fixbo_mag, bvnat_mag)
        !
        nedgdir_mag = 0_ip
        !
        do iedge = 1, meshe(ndivi) % nedge
          !
          kfl_fixno_mag(1_ip:ndime, iedge) = kfl_fixno_mag(1_ip, iedge)
          !
          if (kfl_fixno_mag(1_ip, iedge) > 1_ip) then
            !
            ! Wrong code
            !
            call runend("mag_inibcs: wrong code in mag.dat")
            !
          else if (kfl_fixno_mag(1_ip, iedge) == 1_ip) then
            !
            ! Dirichlet boundary conditions
            !
            direnabled = .true.
            nedgdir_mag = nedgdir_mag + 1_ip
            !
          end if
          !
        end do
        !
!        if (.not. direnabled) call runend("mag_inibcs: missing Dirichlet boundary conditions")
        !
      end if

    case(2_ip)

      if (INOTEMPTY) then
        !
        jedge = 0_ip
        !
        do iedge = 1, meshe(ndivi) % nedge
          if (kfl_fixno_mag(1_ip, iedge) == 1) then
            edgdir_mag(iedge) = .true.
            jedge = jedge + 1_ip
            diredg_mag(jedge) = iedge
          else
            edgdir_mag(iedge) = .false.  
          end if
        end do
        !
      end if

  end select
  
end subroutine mag_inibcs
