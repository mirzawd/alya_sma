!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup Chemic post scalar dissipation distribution
!> @{
!> @file    chm_post_scalar_dist.f90
!> @author  Guillermo Oyarzun
!> @brief   Distributes the post scalar dissipation rate using the packed data
!> @details Data reorganization to include vectorization and GPUs
!>
!> @}
!------------------------------------------------------------------------

subroutine chm_post_scalar_dist(ivari)

  use def_kintyp,                 only : ip,rp
  use mod_parall,                 only : num_subd_par
  use mod_parall,                 only : num_pack_par
  use def_domain,                 only : lnnod
  use def_domain,                 only : lgaus
  use mod_parall,                 only : list_elements_par
  use mod_parall,                 only : typ_list_elements_par
  use def_chemic,                 only : xYr_chm,xZr_chm,xYs_chm,xZs_chm
  use mod_solver,                 only : solver_lumped_mass_system


  implicit none
  integer(ip),  intent(in) :: ivari


  integer(ip)                          :: isubd,ipack,ielem
  integer(ip)                          :: pnode,pgaus
  integer(ip)                          :: num_subd
  integer(ip)                          :: VECTOR_DIM
  integer(ip),                 pointer :: num_pack(:)
  type(typ_list_elements_par), pointer :: list_elements(:)

  external                             :: chm_post_scalar_dissipation_rate_fast

  !
  num_subd      =  num_subd_par
  num_pack      => num_pack_par
  list_elements => list_elements_par

  !
  ! Initialization
  !
  select case (ivari)

    case (23_ip)
      xYr_chm = 0.0_rp

    case (24_ip)
      xZr_chm = 0.0_rp

    case (25_ip)
      xYs_chm = 0.0_rp

    case (26_ip)
      xZs_chm = 0.0_rp

  end select



  !-------------------------------------------------------------------
  !
  ! Packed version of gp_reatab: race condition
  !
  !-------------------------------------------------------------------

  do isubd = 1,num_subd
    do ipack = 1,num_pack(isubd)

        ielem = list_elements(isubd) % packs(ipack) % l(1)  ! Select first element
        pnode = lnnod(ielem)                                ! Number of nodes
        pgaus = lgaus(ielem)                                ! Number of Gauss points

        VECTOR_DIM = size(list_elements(isubd) % packs(ipack) % l,kind=ip)

        call chm_post_scalar_dissipation_rate_fast(VECTOR_DIM,pnode,pgaus,list_elements(isubd) % packs(ipack) % l,ivari)

     end do
  end do

  !
  ! Lumped mass matrix
  !
  select case (ivari)

    case (23_ip)
       call solver_lumped_mass_system(1_ip,xYr_chm)

    case (24_ip)
       call solver_lumped_mass_system(1_ip,xZr_chm)

    case (25_ip)
       call solver_lumped_mass_system(1_ip,xYs_chm)

    case (26_ip)
       call solver_lumped_mass_system(1_ip,xZs_chm)

  end select




end subroutine chm_post_scalar_dist

