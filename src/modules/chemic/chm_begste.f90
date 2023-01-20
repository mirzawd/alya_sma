!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine chm_begste()
  !-----------------------------------------------------------------------
  !****f* partis/chm_begste
  ! NAME
  !    chm_begste
  ! DESCRIPTION
  !    This routine prepares a new time step
  ! USES
  !    chm_updunk
  ! USED BY
  !    partis
  !***
  !-----------------------------------------------------------------------
  use def_master,                 only : ITASK_BEGSTE, momod, modul
  use def_kintyp,                 only : ip
  use def_chemic,                 only : kfl_lookg_chm, kfl_model_chm, kfl_premix_chm, kfl_spray_chm, kfl_ufpv_chm,&
                                         kfl_solve_cond_CMC_chm, kfl_tab_fw_chm, kfl_multimod_chm
  use mod_chm_finiteRate,         only : chm_getProp_finiteRate
  use def_kermod,                 only : kfl_chemic_vect

  implicit none

  external :: chm_updunk
  external :: chm_post_scalar_dissipation_rate
  external :: chm_post_scalar_dist
  external :: chm_gp_reatab
  external :: chm_gp_multi_reatab
  external :: chm_reatab
  external :: chm_gp_ANN_eval

  if (kfl_model_chm /= 4) then
     if(  momod(modul) % kfl_stead /=1) then
        !
        ! Initial guess: c(n,0,*) <-- c(n-1,*,*).
        !
        call chm_updunk(ITASK_BEGSTE)
     end if
  end if

  if ((kfl_model_chm == 1  .or. kfl_model_chm == 2_ip) .and. kfl_tab_fw_chm>=0_ip) then

     !
     ! Read flamelet table for gas phase
     !
     if (kfl_spray_chm == 0 .or. ( kfl_spray_chm /= 0 .and. kfl_premix_chm == 0)) then
        if (kfl_ufpv_chm > 0) then
          if(kfl_chemic_vect /= 1_ip) then
             call chm_post_scalar_dissipation_rate(24_ip)
          else
             call chm_post_scalar_dist(24_ip)
          end if
        endif

        if (kfl_lookg_chm > 0) then
            if (kfl_multimod_chm > 0) then
               call chm_gp_multi_reatab(ITASK_BEGSTE)    
            else
               call chm_gp_reatab(ITASK_BEGSTE)
            endif  
            call chm_gp_ANN_eval(ITASK_BEGSTE)
        else
            call chm_reatab()
        endif
     end if

  elseif (kfl_model_chm == 3) then
       !
       ! Calculate transport properties
       !
       call chm_getProp_finiteRate()

  elseif (kfl_model_chm == 4 .and. kfl_solve_cond_CMC_chm == 0) then
       !
       ! Compute Xres for unconditional mixing variables from CMC model
       !
      if(kfl_chemic_vect /= 1_ip) then
         call chm_post_scalar_dissipation_rate(24_ip)
      else
         call chm_post_scalar_dist(24_ip)
      end if
  endif

end subroutine chm_begste

