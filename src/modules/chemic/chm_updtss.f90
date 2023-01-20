!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine chm_updtss()
  !-----------------------------------------------------------------------
  !****f* Chemic/chm_updtss
  ! NAME
  !    chm_updtss
  ! DESCRIPTION
  !    This routine computes the time step size
  ! USED BY
  !    chm_timste
  !***
  !-----------------------------------------------------------------------
  use def_master,                only : dtinv, kfl_timco
  use def_kintyp,                only : ip, rp
  use def_chemic,                only : ADR_chm, dtcri_chm, dtinv_chm, kfl_dt_calc_CMC_chm, kfl_model_chm, kfl_solve_cond_CMC_chm,&
                                        kfl_start_CMC_chm, nclas_chm, safet_chm
  use mod_chm_finiteRate,        only : chm_updtcc_finiteRate
  use mod_chm_operations_CMC,    only : chm_updtcc_CMC
  use def_kermod,                only : kfl_chemic_vect
  implicit none
  integer(ip) :: iclas
  real(rp)    :: dtmin

  external    :: chm_updtcc_flamLet_dist
  external    :: chm_updtcc_flamLet

  if( ADR_chm(1) % kfl_time_integration /= 0 ) then

     dtmin = huge(1.0_rp)

     !
     ! Time step based on critical time step
     !
     if (kfl_model_chm == 1 .or. kfl_model_chm == 2_ip ) then
        if(kfl_chemic_vect == 1_ip) then
            call chm_updtcc_flamLet_dist(dtmin)
        else
            call chm_updtcc_flamLet(dtmin)
        end if

     else if (kfl_model_chm == 3) then
        call chm_updtcc_finiteRate(dtmin)

     else if (kfl_model_chm == 4) then
        if (kfl_solve_cond_CMC_chm == 1_ip) then
           if (kfl_dt_calc_CMC_chm == 1_ip)  call chm_updtcc_CMC(dtmin)
        else
           if (kfl_start_CMC_chm == 1) then
              dtmin = 1.0e-6_rp  !! Dummy dtmin just for the first time step
           else
              if(kfl_chemic_vect == 1_ip) then
                 call chm_updtcc_flamLet_dist(dtmin)
              else
                 call chm_updtcc_flamLet(dtmin)
              end if
           end if
        end if
     endif

     if (kfl_model_chm == 4 .and. kfl_solve_cond_CMC_chm == 1_ip .and. kfl_dt_calc_CMC_chm == 0_ip) then
        dtinv_chm = dtinv
        dtcri_chm = 1.0_rp / dtinv_chm
     else
        dtcri_chm = dtmin

        if( dtcri_chm /= 0.0_rp ) dtinv_chm = 1.0_rp/(dtcri_chm*safet_chm)
        if( kfl_timco == 1 )     then
           dtinv = max(dtinv,dtinv_chm)
        endif
     end if

     do iclas = 1,nclas_chm
        ADR_chm(iclas) % dtinv = dtinv_chm
     end do

  end if

end subroutine chm_updtss
