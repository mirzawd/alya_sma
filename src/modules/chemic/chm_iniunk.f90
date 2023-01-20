!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine chm_iniunk()
  !-----------------------------------------------------------------------
  !****f* chemic/chm_iniunk
  ! NAME
  !    chm_iniunk
  ! DESCRIPTION
  !    This routine sets up the initial condition for the concentrations.
  ! USED BY
  !    chm_begste
  !***
  !-----------------------------------------------------------------------
  use def_master,                         only : ITASK_INIUNK, kfl_rstar, conce, cutim
  use def_domain,                         only : npoin
  use def_kintyp,                         only : ip
  use def_chemic,                         only : kfl_droplet_id_chm, kfl_model_chm, kfl_solve_cond_CMC_chm, kfl_soot_chm,&
                                                 nclas_chm, bvess_chm, avtim_chm
  use mod_chm_droplets,                   only : chm_droplet_id, chm_droplet_output
  use mod_chm_sectional_soot_model_fast,  only : chm_validation_test_ssm_fast
  use mod_chm_sectional_soot_model_fast,  only : soot_volume_discretization_ssm_fast
  use mod_chm_sectional_soot_model_fast,  only : initialize_constants_ssm_fast
!@
  use mod_chm_sectional_soot_model_fast,  only : chm_validation_test_ssm
  use mod_chm_sectional_soot_model_fast,  only : soot_volume_discretization_ssm
  use mod_chm_sectional_soot_model_fast,  only : initialize_constants_ssm
!@
  use def_kermod,                         only : kfl_soot_vect
  implicit none
  integer(ip) :: iclas,ipoin

  external    :: chm_fields
  external    :: chm_updbcs
  external    :: chm_updunk
  external    :: chm_update_model

  if( kfl_rstar == 0 ) then

     !
     ! Load initial conditions
     !
     if (.not. (kfl_model_chm == 4 .and. kfl_solve_cond_CMC_chm == 1)) then

        call chm_fields()  ! Initialization by fields activated

       do ipoin = 1,npoin
           do iclas = 1,nclas_chm
              conce(ipoin,iclas,1) = bvess_chm(iclas,ipoin)
           end do
        end do
        !
        ! Apply Dirichlet boundary conditions on (:,1)
        !
        call chm_updbcs(ITASK_INIUNK)
        !
        ! Droplet identification and output
        !
        if ( kfl_droplet_id_chm /= 0 ) then
           call chm_droplet_id
           call chm_droplet_output
        end if

        call chm_updunk(ITASK_INIUNK)

     end if

  end if

  !
  ! Update chemistry model
  !
  call chm_update_model()

  !
  ! Initialize sectional soot model
  !
  if (kfl_soot_chm /= 0) then
     if (kfl_soot_chm > 0) then
        if(kfl_soot_vect /=0) then
            !
            ! Determine location of gas phase species
            !
            call initialize_constants_ssm_fast()

            !
            ! Determine sections
            !
            call soot_volume_discretization_ssm_fast()
        else

            !
            ! Determine location of gas phase species
            !
            call initialize_constants_ssm()

            !
            ! Determine sections
            !
            call soot_volume_discretization_ssm()

        end if
     end if
  end if

  !
  ! Soot calculation validation test
  !
  if (kfl_soot_chm < 0) then
    if(kfl_soot_vect /= 0) then
         call chm_validation_test_ssm_fast()
    else
         call chm_validation_test_ssm()
    end if
  end if

  !
  ! Update the time with the current one for averages
  !
  avtim_chm = cutim

end subroutine chm_iniunk
