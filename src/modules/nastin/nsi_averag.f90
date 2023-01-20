!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_averag()
  !------------------------------------------------------------------------
  !****f* Nastin/nsi_averag
  ! NAME 
  !    nsi_averag
  ! DESCRIPTION
  !    Average velocity and pressure
  ! USES
  ! USED BY
  !    nsi_endste 
  !***
  !------------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_nastin
  use mod_gradie
  use mod_projec,             only : projec_elements_to_nodes
  use def_kermod,             only : kfl_wlaav_ker,nsteps_ensemble,kfl_twola_ker,kfl_noslw_ker,kfl_waexl_ker
  use mod_ker_proper,         only : ker_proper
  use mod_output_postprocess, only : output_postprocess_check_variable_postprocess
  use mod_solver,             only : solver_lumped_mass_system 
  use mod_memory,             only : memory_alloca, memory_deallo
  use mod_arrays,             only : arrays_number
  implicit none
  integer(ip)       :: idime,ipoin,dummi,iauxi
  real(rp)          :: auxii, minmod, maxmod, avemod
  real(rp), pointer :: velgr(:,:)
  real(rp), pointer :: prope_tmp(:,:)
  real(rp), pointer :: densi_tmp(:)

  if ( (kfl_wlaav_ker == 1_ip) .and. ( (kfl_waexl_ker == 0_ip) .or. (kfl_noslw_ker == 1_ip) ) .and. (ittim > 0_ip) ) then
     call nsi_wallav(2_ip)
  end if



  if( cutim > avtim_nsi .and. INOTEMPTY ) then

     if ( kfl_wlaav_ker == 1_ip ) then  
        call nsi_wallav(1_ip)      ! actualizes averaged velocity at gauss points (velav_ker)
        if (( (kfl_waexl_ker == 0_ip).or.(kfl_noslw_ker == 1_ip)) .and. (ittim .gt. 0_ip) ) &
             call nsi_wallav(2_ip) ! actualizes  averaged velocity at points (avupo_nsi)
     end if
     !AVVEL
     if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVVEL') ) then
        !if(kfl_regim_nsi == 3 .and. (kfl_coupl(ID_TEMPER,ID_CHEMIC) == 0) ) then
        !   do ipoin=1,npoin
        !      do idime=1,ndime
        !         avvel_nsi(idime,ipoin)=avvel_nsi(idime,ipoin)&
        !            +veloc(idime,ipoin,1)*dtime*(prthe(1)/gasco)*(1.0_rp/tempe(ipoin,1))
        !      end do
        !   end do
        !else 
        do ipoin=1,npoin
           do idime=1,ndime
              avvel_nsi(idime,ipoin)=avvel_nsi(idime,ipoin)&
                   +veloc(idime,ipoin,1)*dtime
           end do
        end do
        !end if
     end if

     ! AVPRE
     if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVPRE') ) then
        do ipoin = 1,npoin
           avpre_nsi(ipoin)=avpre_nsi(ipoin)&
                +press(ipoin,1)*dtime
        end do
     end if

     ! AVTAN
     if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVTAN') ) then
        call memgen(zero,ndime,npoin)
        call nsi_memphy(4_ip)  ! allocate prope_nsi(npoin, 2) 
        call nsi_denvis()      ! loads density and viscosity in prope_nsi(npoin,2) 
        call nsi_outtan(1_ip)      ! tangent traction in gevec, for no slip wall law, variational forces are obtained
        call nsi_memphy(-4_ip)
        do ipoin=1,npoin
           do idime=1,ndime
              avtan_nsi(idime,ipoin) = avtan_nsi(idime,ipoin)&
                   + gevec(idime,ipoin)*dtime
           end do
        end do
        call memgen(2_ip,ndime,npoin)
     end if

     ! V**2
     if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVVE2') ) then
        !if(kfl_regim_nsi == 3 .and. (kfl_coupl(ID_TEMPER,ID_CHEMIC) == 0) ) then
        !   do ipoin=1,npoin
        !      do idime=1,ndime
        !         auxii = veloc(idime,ipoin,1)                 
        !         avve2_nsi(idime,ipoin)=avve2_nsi(idime,ipoin)&
        !            +auxii*auxii*dtime*(prthe(1)/gasco)*(1.0_rp/tempe(ipoin,1))
        !      end do
        !   end do
        !else 
        do ipoin=1,npoin
           do idime=1,ndime
              auxii = veloc(idime,ipoin,1)                 
              avve2_nsi(idime,ipoin)=avve2_nsi(idime,ipoin)&
                   +auxii*auxii*dtime
           end do
        end do
        !end if

     end if

     !Vx*Vy
     if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVVXY') ) then
        !if(kfl_regim_nsi == 3 .and. (kfl_coupl(ID_TEMPER,ID_CHEMIC) == 0) ) then
        !   do ipoin=1,npoin
        !      avvxy_nsi(1,ipoin)=avvxy_nsi(1,ipoin)&
        !         +veloc(1,ipoin,1)*&
        !         veloc(2,ipoin,1)*dtime*(prthe(1)/gasco)*(1.0_rp/tempe(ipoin,1))
        !      if (ndime==3) then
        !         avvxy_nsi(2,ipoin)=avvxy_nsi(2,ipoin)&
        !            +veloc(2,ipoin,1)*&
        !            veloc(3,ipoin,1)*dtime*(prthe(1)/gasco)*(1.0_rp/tempe(ipoin,1))
        !         avvxy_nsi(3,ipoin)=avvxy_nsi(3,ipoin)&
        !            +veloc(3,ipoin,1)*&
        !            veloc(1,ipoin,1)*dtime*(prthe(1)/gasco)*(1.0_rp/tempe(ipoin,1))
        !      end if
        !   end do
        !else
        do ipoin=1,npoin
           avvxy_nsi(1,ipoin)=avvxy_nsi(1,ipoin)&
                +veloc(1,ipoin,1)*&
                veloc(2,ipoin,1)*dtime
           if (ndime==3) then
              avvxy_nsi(2,ipoin)=avvxy_nsi(2,ipoin)&
                   +veloc(2,ipoin,1)*&
                   veloc(3,ipoin,1)*dtime
              avvxy_nsi(3,ipoin)=avvxy_nsi(3,ipoin)&
                   +veloc(3,ipoin,1)*&
                   veloc(1,ipoin,1)*dtime
           end if
        end do
        !end if
     end if

     !P**2
     if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVPR2') ) then
        do ipoin = 1,npoin
           avpr2_nsi(ipoin)=avpr2_nsi(ipoin)&
                +press(ipoin,1)*press(ipoin,1)*dtime
        end do
     end if

     ! AVMUT ! is the averaged nu_t
     if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVMUT') ) then
        call memgen(zero,npoin,zero)
        call ker_proper('TURBU','NPOIN',dummi,dummi,gesca)

        !           call memgen(zero,npoin,zero)
        !           call projec_elements_to_nodes(turmu_nsi,gesca)
        do ipoin=1,npoin
           avmut_nsi(ipoin) = avmut_nsi(ipoin)&
                + gesca(ipoin)*dtime
        end do
        call memgen(2_ip,npoin,zero)
     end if

     ! AVSTX: nu_t * grad(u)
     if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVSTX') ) then
        call memory_alloca(mem_modul(1:2,modul),'VELGR','nsi_averag',velgr,ndime,npoin)
        call gradie(veloc(1,1:npoin,1),velgr)
        call memgen(zero,npoin,zero)
        call ker_proper('TURBU','NPOIN',dummi,dummi,gesca)
        !           call projec_elements_to_nodes(turmu_nsi,gesca)
        do ipoin=1,npoin
           do idime=1,ndime
              avstx_nsi(idime,ipoin) = avstx_nsi(idime,ipoin)&
                   + gesca(ipoin) * velgr(idime,ipoin) * dtime
           end do
        end do
        call memgen(2_ip,npoin,zero)
        call memory_deallo(mem_modul(1:2,modul),'VELGR','nsi_averag',velgr)
     end if

     ! AVSTY: nu_t * grad(v)
     if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVSTY') ) then
        call memory_alloca(mem_modul(1:2,modul),'VELGR','nsi_averag',velgr,ndime,npoin)
        call gradie(veloc(2,1:npoin,1),velgr)
        call memgen(zero,npoin,zero)
        call ker_proper('TURBU','NPOIN',dummi,dummi,gesca)
        !           call projec_elements_to_nodes(turmu_nsi,gesca)
        do ipoin=1,npoin
           do idime=1,ndime
              avsty_nsi(idime,ipoin) = avsty_nsi(idime,ipoin)&
                   + gesca(ipoin) * velgr(idime,ipoin) * dtime
           end do
        end do
        call memgen(2_ip,npoin,zero)
        call memory_deallo(mem_modul(1:2,modul),'VELGR','nsi_averag',velgr)
     end if

     ! AVSTZ: nu_t * grad(w)
     if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVSTZ') ) then
        call memory_alloca(mem_modul(1:2,modul),'VELGR','nsi_averag',velgr,ndime,npoin)
        call gradie(veloc(3,1:npoin,1),velgr)
        call memgen(zero,npoin,zero)
        call ker_proper('TURBU','NPOIN',dummi,dummi,gesca)
        !          call projec_elements_to_nodes(turmu_nsi,gesca)
        do ipoin=1,npoin
           do idime=1,ndime
              avstz_nsi(idime,ipoin) = avstz_nsi(idime,ipoin)&
                   + gesca(ipoin) * velgr(idime,ipoin) * dtime
           end do
        end do
        call memgen(2_ip,npoin,zero)
        call memory_deallo(mem_modul(1:2,modul),'VELGR','nsi_averag',velgr)
     end if

     ! AVMOS: average momentum source term from spray
     if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVMOS')  ) then

        nullify ( prope_tmp )
        call memory_alloca(mem_modul(1:2,modul),'PROPE_TMP','nsi_averag',prope_tmp,ndime,npoin)
        if (associated(momentum_sink)) then
           do ipoin=1,npoin
              do idime = 1,ndime
                 prope_tmp(idime,ipoin) = momentum_sink(idime,ipoin)
              enddo
           enddo
           call solver_lumped_mass_system(ndime,prope_tmp,EXCHANGE=.false.)
        else
           do ipoin=1,npoin
              do idime = 1,ndime
                 prope_tmp(idime,ipoin) = 0.0_rp 
              enddo
           end do
        endif

        do ipoin=1,npoin
           do idime = 1,ndime
              avmos_nsi(idime,ipoin) = avmos_nsi(idime,ipoin) + prope_tmp(idime,ipoin) * dtime
           enddo
        end do
        call memory_deallo(mem_modul(1:2,modul),'PROPE_TMP','nsi_averag',prope_tmp)

     end if

     ! AVMFL: average mass flux
     if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVMFL')  ) then

        nullify ( densi_tmp )
        call memory_alloca(mem_modul(1:2,modul),'DENSI_TMP','nsi_averag',densi_tmp,npoin)
        call ker_proper('DENSI','NPOIN',dummi,dummi,densi_tmp)
        do ipoin=1,npoin
           do idime = 1,ndime
              av_mass_flux_nsi(idime,ipoin) = av_mass_flux_nsi(idime,ipoin) + densi_tmp(ipoin) * veloc(idime,ipoin,1) * dtime
           enddo
        end do
        call memory_deallo(mem_modul(1:2,modul),'DENSI_TMP','nsi_averag',densi_tmp)

     end if

     ! AVRUU: average momentum flux, diagonal terms
     if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVRUU')  ) then

        nullify ( densi_tmp )
        call memory_alloca(mem_modul(1:2,modul),'DENSI_TMP','nsi_averag',densi_tmp,npoin)
        call ker_proper('DENSI','NPOIN',dummi,dummi,densi_tmp)
        do ipoin=1,npoin
           do idime = 1,ndime
              av_mom_flux_diag_nsi(idime,ipoin) = av_mom_flux_diag_nsi(idime,ipoin) + densi_tmp(ipoin) * veloc(idime,ipoin,1) * veloc(idime,ipoin,1)* dtime
           enddo
        end do
        call memory_deallo(mem_modul(1:2,modul),'DENSI_TMP','nsi_averag',densi_tmp)

     end if

     ! AVRUV: average momentum flux, off-diagonal terms
     if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVRUV')  ) then

        nullify ( densi_tmp )
        call memory_alloca(mem_modul(1:2,modul),'DENSI_TMP','nsi_averag',densi_tmp,npoin)
        call ker_proper('DENSI','NPOIN',dummi,dummi,densi_tmp)
        do ipoin=1,npoin
           av_mom_flux_off_nsi(1,ipoin) = av_mom_flux_off_nsi(1,ipoin) + densi_tmp(ipoin) * veloc(1,ipoin,1) * veloc(2,ipoin,1)* dtime
           if (ndime==3) then
              av_mom_flux_off_nsi(2,ipoin) = av_mom_flux_off_nsi(2,ipoin) + densi_tmp(ipoin) * veloc(2,ipoin,1) * veloc(3,ipoin,1)* dtime
              av_mom_flux_off_nsi(3,ipoin) = av_mom_flux_off_nsi(3,ipoin) + densi_tmp(ipoin) * veloc(3,ipoin,1) * veloc(1,ipoin,1)* dtime
           endif
        end do
        call memory_deallo(mem_modul(1:2,modul),'DENSI_TMP','nsi_averag',densi_tmp)

     end if


     !VMAXP
     if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='VMAXP') .or. &
         output_postprocess_check_variable_postprocess(VARIABLE_NAME='PINDE') ) then
           do ipoin = 1,npoin
             do idime=1,ndime
                 if ( abs(veloc(idime,ipoin,1)) .gt. abs(vmaxp_nsi(idime,ipoin)) ) then
                     vmaxp_nsi(idime,ipoin)=veloc(idime,ipoin,1)
                 endif
             end do
           end do
     end if
     !VMINP
     if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='VMINP') .or. &
         output_postprocess_check_variable_postprocess(VARIABLE_NAME='PINDE') ) then
           do ipoin = 1,npoin
             do idime=1,ndime
                 if ( abs(veloc(idime,ipoin,1)) .lt. abs(vminp_nsi(idime,ipoin)) ) then
                     vminp_nsi(idime,ipoin)=veloc(idime,ipoin,1)
                 endif
             end do
           end do
       end if
       !VAVEP
       if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='VAVEP') .or. &
           output_postprocess_check_variable_postprocess(VARIABLE_NAME='PINDE') ) then
             do ipoin = 1,npoin
                 do idime=1,ndime
                     vavep_nsi(idime,ipoin)=vavep_nsi(idime,ipoin) + veloc(idime,ipoin,1)*dtime
                 end do
               end do
         end if
         !PINDE
           if(output_postprocess_check_variable_postprocess(VARIABLE_NAME='PINDE') ) then
               do ipoin = 1,npoin
                 maxmod=sqrt(sum(vmaxp_nsi(1:ndime,ipoin)*vmaxp_nsi(1:ndime,ipoin)))
                 minmod=sqrt(sum(vminp_nsi(1:ndime,ipoin)*vminp_nsi(1:ndime,ipoin)))
                 avemod=sqrt(sum(vavep_nsi(1:ndime,ipoin)*vavep_nsi(1:ndime,ipoin)))
 
                 if(avemod .ne. 0.0_rp) then
                    pinde_nsi(1,ipoin)=(maxmod-minmod)/avemod
                 else
                    pinde_nsi(1,ipoin)=0.0_rp
                 endif
               end do
         end if


     if (nsteps_ensemble > 0 ) then

        !ENVEL
        if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='ENVEL') ) then
           iauxi = mod(ittim, postp(1) % npp_stepi(arrays_number('ENVEL'),0) )
           if ( iauxi > 0 .and. iauxi <= nsteps_ensemble ) then
              do ipoin=1,npoin
                 do idime=1,ndime
                    envel_nsi(idime,ipoin)=envel_nsi(idime,ipoin)&
                         +veloc(idime,ipoin,1)*dtime
                 end do
              end do
           end if
        end if

        ! ENPRE
        if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='ENPRE') ) then
           iauxi = mod(ittim, postp(1) % npp_stepi(arrays_number('ENPRE'),0) )
           if ( iauxi > 0 .and. iauxi <= nsteps_ensemble ) then
              do ipoin = 1,npoin
                 enpre_nsi(ipoin)=enpre_nsi(ipoin)&
                      +press(ipoin,1)*dtime
              end do
           end if
        end if

        ! V**2
        if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='ENVE2') ) then
           iauxi = mod(ittim, postp(1) % npp_stepi(arrays_number('ENVE2'),0) )
           if ( iauxi > 0 .and. iauxi <= nsteps_ensemble ) then
              do ipoin=1,npoin
                 do idime=1,ndime
                    auxii = veloc(idime,ipoin,1)                 
                    enve2_nsi(idime,ipoin)=enve2_nsi(idime,ipoin)&
                         +auxii*auxii*dtime
                 end do
              end do
           end if
        end if

        !Vx*Vy
        if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='ENVXY') ) then
           iauxi = mod(ittim, postp(1) % npp_stepi(arrays_number('ENVXY'),0) )
           if ( iauxi > 0 .and. iauxi <= nsteps_ensemble ) then
              do ipoin=1,npoin
                 envxy_nsi(1,ipoin)=envxy_nsi(1,ipoin)&
                      +veloc(1,ipoin,1)*&
                      veloc(2,ipoin,1)*dtime
                 if (ndime==3) then
                    envxy_nsi(2,ipoin)=envxy_nsi(2,ipoin)&
                         +veloc(2,ipoin,1)*&
                         veloc(3,ipoin,1)*dtime
                    envxy_nsi(3,ipoin)=envxy_nsi(3,ipoin)&
                         +veloc(3,ipoin,1)*&
                         veloc(1,ipoin,1)*dtime
                 end if
              end do
           end if
        end if

        !P**2
        if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='ENPR2') ) then
           iauxi = mod(ittim, postp(1) % npp_stepi(arrays_number('ENPR2'),0) )
           if ( iauxi > 0 .and. iauxi <= nsteps_ensemble ) then
              do ipoin = 1,npoin
                 enpr2_nsi(ipoin)=enpr2_nsi(ipoin)&
                      +press(ipoin,1)*press(ipoin,1)*dtime
              end do
           end if
        end if

        ! ENTAN
        if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='ENTAN') ) then
           iauxi = mod(ittim, postp(1) % npp_stepi(arrays_number('ENTAN'),0) )
           if ( iauxi > 0 .and. iauxi <= nsteps_ensemble ) then
              call memgen(zero,ndime,npoin)
              call nsi_memphy(4_ip)
              call nsi_denvis()
              call nsi_outtan(1_ip)  ! for no slip wall law, variational forces are obtained - comes in gevec
              call nsi_memphy(-4_ip)
              do ipoin=1,npoin
                 do idime=1,ndime
                    entan_nsi(idime,ipoin) = entan_nsi(idime,ipoin)&
                         + gevec(idime,ipoin)*dtime
                 end do
              end do
              call memgen(2_ip,ndime,npoin)
           end if
        end if

        ! ENMUT
        if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='ENMUT') ) then
           iauxi = mod(ittim, postp(1) % npp_stepi(arrays_number('ENMUT'),0) )
           if ( iauxi > 0 .and. iauxi <= nsteps_ensemble ) then
              call memgen(zero,npoin,zero)
              call projec_elements_to_nodes(turmu_nsi,gesca)
              do ipoin=1,npoin
                 enmut_nsi(ipoin) = enmut_nsi(ipoin)&
                      + gesca(ipoin)*dtime
              end do
              call memgen(2_ip,npoin,zero)
           end if
        end if

        ! ENSTX: mu_t * grad(u)
        if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='ENSTX') ) then
           iauxi = mod(ittim, postp(1) % npp_stepi(arrays_number('ENSTX'),0) )
           if ( iauxi > 0 .and. iauxi <= nsteps_ensemble ) then
              call memory_alloca(mem_modul(1:2,modul),'VELGR','nsi_averag',velgr,ndime,npoin)
              call gradie(veloc(1,1:npoin,1),velgr)
              call memgen(zero,npoin,zero)
              call projec_elements_to_nodes(turmu_nsi,gesca)
              do ipoin=1,npoin
                 do idime=1,ndime
                    enstx_nsi(idime,ipoin) = enstx_nsi(idime,ipoin)&
                         + gesca(ipoin) * velgr(idime,ipoin) * dtime
                 end do
              end do
              call memgen(2_ip,npoin,zero)
              call memory_deallo(mem_modul(1:2,modul),'VELGR','nsi_averag',velgr)
           end if
        end if

        ! ENSTY: mu_t * grad(v)
        if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='ENSTY') ) then
           iauxi = mod(ittim, postp(1) % npp_stepi(arrays_number('ENSTY'),0) )
           if ( iauxi > 0 .and. iauxi <= nsteps_ensemble ) then
              call memory_alloca(mem_modul(1:2,modul),'VELGR','nsi_averag',velgr,ndime,npoin)
              call gradie(veloc(2,1:npoin,1),velgr)
              call memgen(zero,npoin,zero)
              call projec_elements_to_nodes(turmu_nsi,gesca)
              do ipoin=1,npoin
                 do idime=1,ndime
                    ensty_nsi(idime,ipoin) = ensty_nsi(idime,ipoin)&
                         + gesca(ipoin) * velgr(idime,ipoin) * dtime
                 end do
              end do
              call memgen(2_ip,npoin,zero)
              call memory_deallo(mem_modul(1:2,modul),'VELGR','nsi_averag',velgr)
           end if
        end if

        ! ENSTZ: mu_t * grad(w)
        if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='ENSTZ') ) then
           iauxi = mod(ittim, postp(1) % npp_stepi(arrays_number('ENSTZ'),0) )
           if ( iauxi > 0 .and. iauxi <= nsteps_ensemble ) then
              call memory_alloca(mem_modul(1:2,modul),'VELGR','nsi_averag',velgr,ndime,npoin)
              call gradie(veloc(3,1:npoin,1),velgr)
              call memgen(zero,npoin,zero)
              call projec_elements_to_nodes(turmu_nsi,gesca)
              do ipoin=1,npoin
                 do idime=1,ndime
                    enstz_nsi(idime,ipoin) = enstz_nsi(idime,ipoin)&
                         + gesca(ipoin) * velgr(idime,ipoin) * dtime
                 end do
              end do
              call memgen(2_ip,npoin,zero)
              call memory_deallo(mem_modul(1:2,modul),'VELGR','nsi_averag',velgr)
           end if
        end if



     end if ! ensemble

     if ( kfl_twola_ker == 1_ip ) then
        call nsi_tluave(1_ip)
     end if

     if ( kfl_twola_ker == 2_ip ) then
        call nsi_coutan()
     end if

  end if

end subroutine nsi_averag
