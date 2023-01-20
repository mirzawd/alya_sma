!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine chm_endste()
  !-----------------------------------------------------------------------
  !****f* Chemic/chm_endste
  ! NAME
  !    chm_endste
  ! DESCRIPTION
  !    This routine ends a time step of the transport equation.
  ! USES
  !    chm_cvgunk
  !    chm_updunk
  !    chm_output
  !    chm_restar
  ! USED BY
  !    Chemic
  !***
  !-----------------------------------------------------------------------
  use def_master,                 only : INOTMASTER, ITASK_ENDSTE, momod, ittim, kfl_gotim, modul, mem_modul
  use def_kintyp,                 only : ip, rp, r1p
  use def_chemic,                 only : d32_gp_chm, sigma_gp_chm, sigma0_gp_chm, droplet_postprocess_frequency_chm, d32_chm,&
                                         kfl_droplet_id_chm, kfl_model_chm, kfl_spray_chm, kfl_timei_chm, Sigm0_chm, Sigma_chm
  use def_domain,                 only : ngaus, ltype, nelem
  use mod_memory,                 only : memory_alloca, memory_deallo
  use mod_chm_droplets,           only : chm_droplet_id, chm_droplet_output
  use def_kermod,                 only : kfl_chemic_vect

  implicit none
  integer(ip)             :: pelty,pgaus,ielem
  type(r1p),pointer       :: aux_r1p(:)

  external                :: chm_cvgunk
  external                :: chm_updunk
  external                :: smooth
  external                :: chm_averag
  external                :: chm_averag_fast

  !
  ! Compute convergence residual of the time evolution (that is,
  ! || c(n,*,*) - c(n-1,*,*)|| / ||c(n,*,*)||) and update unknowns
  ! c(n-1,*,*) <-- c(n,*,*)
  !
  if(momod(modul) % kfl_stead==0.and.kfl_timei_chm==1) then
     call chm_cvgunk(ITASK_ENDSTE)
     call chm_updunk(ITASK_ENDSTE)
     !
     ! Droplet identification and output
     !
     if ( kfl_droplet_id_chm /= 0 ) then
        if ( mod(ittim,droplet_postprocess_frequency_chm) == 0 ) then
           call chm_droplet_id
           call chm_droplet_output
        end if
     end if
  end if
  !
  ! Other updates
  !
!!DMM  call chm_upwmea(ITASK_ENDSTE)                          ! wmean(ipoin,1) ==> wmean(ipoin,3)
  !
  ! Read flamelet table for gas phase
  !
  if (kfl_model_chm == 1 .or. kfl_model_chm == 2) then
     if ( kfl_spray_chm > 0 ) then
        !
        ! For post-processing
        !
        if (INOTMASTER) then
           d32_chm   = 0.0_rp
           Sigma_chm = 0.0_rp
           Sigm0_chm = 0.0_rp

           nullify(aux_r1p)
           call memory_alloca(mem_modul(1:2,modul),'AUX_R1P','chm_endste',aux_r1p,nelem)
           do ielem = 1,nelem
              pelty = ltype(ielem)
              pgaus = ngaus(pelty)
              call memory_alloca(mem_modul(1:2,modul),'AUX_R1P % A','chm_endste',aux_r1p(ielem)%a,pgaus)
           end do

           do ielem = 1,nelem
              aux_r1p(ielem) % a = d32_gp_chm(ielem) % a(:,1,1)
           end do
           call smooth ( aux_r1p, d32_chm)

           do ielem = 1,nelem
              aux_r1p(ielem) % a = sigma_gp_chm(ielem) % a(:,1,1)
           end do
           call smooth ( aux_r1p, Sigma_chm)

           do ielem = 1,nelem
              aux_r1p(ielem) % a = sigma0_gp_chm(ielem) % a(:,1,1)
           end do
           call smooth ( aux_r1p, Sigm0_chm)

           call memory_deallo(mem_modul(1:2,modul),'AUX_R1P','chm_endste',aux_r1p)
        end if
     end if

     !
     ! Print performance of tab lookup
     !
     !call print_tab_interp_times()
  end if
  !
  ! Compute averaged variables
  !
  if(kfl_chemic_vect /= 1) then
      call chm_averag()
  else
      call chm_averag_fast()
  end if

  !
  ! If not steady, go on
  !
  if(momod(modul) % kfl_stead==0) then
     if(kfl_timei_chm==1) kfl_gotim = 1
  end if

end subroutine chm_endste
