!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine chm_reaous()
  !------------------------------------------------------------------------
  !****f* partis/chm_reaous
  ! NAME
  !    chm_reaphy
  ! DESCRIPTION
  !    This routine reads the output strategy
  ! USES
  ! USED BY
  !    chm_turnon
  !***
  !------------------------------------------------------------------------
!  use def_parame
  use def_inpout,             only : words, getint, param
  use def_master,             only : kfl_paral, kfl_timco
  use def_chemic,             only : ipara_chm, write_uncond_spec_CMC_chm, write_cond_spec_CMC_chm, write_iZ_spec_CMC_chm,&
                                     droplet_postprocess_frequency_chm, kfl_dt_calc_CMC_chm, kfl_model_chm, kfl_post_gp_CMC_chm,&
                                     kfl_solve_cond_CMC_chm, nclas_chm, nspec_cond_write_CMC_chm, nspec_uncond_write_CMC_chm,&
                                     nZ_CMC_chm, nZ_write_CMC_chm, rpara_chm, kfl_av_species_CMC_chm, kfl_av_enthalp_CMC_chm
  use def_kintyp,             only : ip, rp
  use mod_ecoute,             only : ecoute
  use mod_chm_operations_CMC, only : chm_initial_actions_reaous_CMC
  use mod_output_postprocess, only : output_postprocess_read
  implicit none
  integer(ip) :: ii

  external    :: runend
  external    :: chm_memous

  if(kfl_paral<=0) then
     !
     ! Initializations
     !
     ipara_chm     = 0       ! Integer parameters
     rpara_chm     = 0.0_rp  ! Real parameters
     droplet_postprocess_frequency_chm = 1_ip ! Droplet postprocess frequency
     kfl_post_gp_CMC_chm = 0_ip   ! 0: postprocessing on nodes, 1: postprocessing on Gauss points
     kfl_dt_calc_CMC_chm = 1_ip   ! 0: omit the calculation of dt if dt is prescribed, 1: otherwise
     kfl_av_species_CMC_chm = 0_ip  ! 0: do not post-process averages for Yk in CFD for CMC model; 1: otherwise
     kfl_av_enthalp_CMC_chm = 0_ip  ! 0: do not post-process averages for enthalpy in CFD for CMC model; 1: otherwise
     !
     ! Reach the section
     !
     call ecoute('chm_reaous')
     do while(words(1)/='OUTPU')
        call ecoute('chm_reaous')
     end do
     !
     ! Begin to read data.
     !
     do while(words(1)/='ENDOU')
        call ecoute('chm_reaous')

        call output_postprocess_read()

        if( words(1) == 'SPECY' ) then
           ipara_chm(1)=getint('SPECY',1_ip,'#Specy number')
        else if ( words(1) == 'FREQD' ) then
           droplet_postprocess_frequency_chm = getint('FREQD',1_ip,'#Droplet postprocess frequency')
           if ( droplet_postprocess_frequency_chm == 0 ) then
              call runend('CHEMIC REAOUS: Droplet postrprocess frequency should not be zero, look at your chm.dat file')
           end if

        else if ( words(1) == 'CMCMO' ) then
            !
            ! Data for CMC model
            !
            call ecoute('chm_reaous')

            CMC_model: do while(words(1)/='ENDCM')

               if ( words(1) == 'UNCON' ) then
                  if ( words(2) == 'TOTAL' ) then
                      nspec_uncond_write_CMC_chm = getint('TOTAL',0_ip,'#Number of species for unconditional mass fractions to be&
                          & post-processed')
                      if (nspec_uncond_write_CMC_chm >= nclas_chm)  call runend('CHEMIC REAOUS: number of species to be&
                          & post-processed is equal or exceeds the number of species in the mechanism')
                  else
                      call runend('CHEMIC REAOUS: number of species for unconditional mass fractions to be post-processed not&
                          & given')
                  end if
                  call chm_memous(1_ip)

                  if ( words(3) == 'POSIT' ) then
                     do ii = 1, nspec_uncond_write_CMC_chm
                        if (param(ii+2)<real(1,rp) .or. param(ii+2)>real(nclas_chm,rp)) then
                           call runend('CHEMIC REAOUS: number of species to be post-processed not valid')
                        else
                           write_uncond_spec_CMC_chm(ii) = int(param(ii+2_ip), kind=ip)
                        end if
                     end do
                  else
                     call runend('CHEMIC REAOUS: position of species for unconditional mass fractions to be post-processed not&
                        & given')
                  end if


               else if ( words(1) == 'CONDI' ) then
                  if ( words(2) == 'SPECI' ) then
                     if ( words(3) == 'TOTAL' ) then
                         nspec_cond_write_CMC_chm = getint('TOTAL',0_ip,'#Number of species for conditional mass fractions to be&
                            & post-processed')
                         if (nspec_cond_write_CMC_chm >= nclas_chm)  call runend('CHEMIC REAOUS: number of species to be&
                            & post-processed is equal or exceeds the number of species in the mechanism')
                     else
                         call runend('CHEMIC REAOUS: number of species for conditional mass fractions to be post-processed not&
                            & given')
                     end if
                     call chm_memous(2_ip)
                     if ( words(4) == 'POSIT' ) then
                        do ii = 1, nspec_cond_write_CMC_chm
                           if (param(ii+3)<real(1,rp) .or. param(ii+3)>real(nclas_chm,rp)) then
                              call runend('CHEMIC REAOUS: number of species to be post-processed not valid')
                           else
                              write_cond_spec_CMC_chm(ii) = int(param(ii+3_ip), kind=ip)
                           end if
                        end do
                     else
                        call runend('CHEMIC REAOUS: position of species for conditional mass fractions to be post-processed not&
                            & given')
                     end if

                  else if ( words(2) == 'MIXTU' ) then
                     if ( words(3) == 'TOTAL' ) then
                        nZ_write_CMC_chm = getint('TOTAL',0_ip,'#Number of mixture fractions to be written')
                        if (nZ_write_CMC_chm >= nZ_CMC_chm)  call runend('CHEMIC REAOUS: number of mixture fractions to be&
                            & post-processed is equal or exceeds the total number of mixture fractions')
                     else
                         call runend('CHEMIC REAOUS: number of mixture fractions for conditional mass fractions to be&
                            & post-processed not given')
                     end if
                     call chm_memous(3_ip)
                     if ( words(4) == 'POSIT' ) then
                        do ii = 1, nZ_write_CMC_chm
                           if (param(ii+3)<real(1,rp) .or. param(ii+3)>real(nZ_CMC_chm,rp)) then
                              call runend('CHEMIC REAOUS: position of mixture fractions to be post-processed not valid')
                           else
                              write_iZ_spec_CMC_chm(ii) = int(param(ii+3_ip), kind=ip)
                           end if
                        end do
                     else
                        call runend('CHEMIC REAOUS: position of mixture fractions for conditional mass fractions to be&
                            & post-processed not given')
                     end if
                  end if

               else if( words(1) == 'POSTP' ) then  ! Where to compute postprocessing fields
                  if (words(2) == 'GAUSS')   kfl_post_gp_CMC_chm = 1_ip

               else if( words(1) == 'TIMES' ) then  ! Omit the calculation of the time step if prescribed
                  if (words(2) == 'OMIT ' .and. kfl_timco == 0_ip)   kfl_dt_calc_CMC_chm = 0_ip
               end if

               call ecoute('chm_reaous')

            end do CMC_model
        end if

     end do

     if (kfl_model_chm == 4 .and. kfl_solve_cond_CMC_chm == 1) then
        if (.not. associated(write_uncond_spec_CMC_chm)) then
           nspec_uncond_write_CMC_chm = nclas_chm
           call chm_memous(1_ip)
           write_uncond_spec_CMC_chm = (/( ii, ii=1_ip, nspec_uncond_write_CMC_chm )/)
        end if

        if (.not. associated(write_cond_spec_CMC_chm)) then
           nspec_cond_write_CMC_chm = nclas_chm
           call chm_memous(2_ip)
           write_cond_spec_CMC_chm = (/( ii, ii=1_ip, nspec_cond_write_CMC_chm )/)
        end if

        if (.not. associated(write_iZ_spec_CMC_chm)) then
           nZ_write_CMC_chm = nZ_CMC_chm
           call chm_memous(3_ip)
           write_iZ_spec_CMC_chm = (/( ii, ii=1_ip, nZ_write_CMC_chm )/)
        end if

        call chm_initial_actions_reaous_CMC
     end if

  end if

end subroutine chm_reaous
