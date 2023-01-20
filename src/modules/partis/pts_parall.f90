!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine pts_parall(order)
  !-----------------------------------------------------------------------
  !****f* Parall/pts_parall
  ! NAME
  !    pts_parall
  ! DESCRIPTION
  !    This routine exchange data 
  ! USES
  ! USED BY
  !    Reapro
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_kermod
  use def_domain
  use def_inpout
  use def_partis
  use mod_opebcs
  use mod_memory,         only : memory_alloca
  use mod_memory,         only : memory_deallo
  use mod_communications, only : PAR_BROADCAST
  use mod_communications, only : PAR_EXCHANGE
  use mod_pts_injection,  only : pts_injection_parallelization
  use mod_interp_tab,     only : tab_par_exchange
  use mod_interp_tab,     only : fw_par_exchange 
  use mod_interp_tab,     only : fw_allocate
  use mod_interp_tab,     only : tab_init_fw
  use mod_result_io,      only : result_io_exchange_add
  use mod_exchange,       only : exchange_init, exchange_end
  use mod_output_postprocess, only : output_postprocess_parall_old
  implicit none
  integer(ip), intent(in) :: order
  integer(ip)             :: ii,jj, kk

  if( ISEQUEN ) return


  select case ( order )

  case( 1_ip )

     !-------------------------------------------------------------------
     !
     ! Broadcast data read in *.pts.dat file
     !
     !-------------------------------------------------------------------
     call exchange_init()
     call result_io_exchange_add()
     call exchange_end()


     do parii = 1,2 
        npari = 0 
        nparr = 0
        nparc = 0
        nparl = 0
        !
        ! Physical problem
        !
        call PAR_EXCHANGE(mlagr,                parin,npari,parii)
        call PAR_EXCHANGE(dtmin_pts,            parre,nparr,parii)
        call PAR_EXCHANGE(dimin_pts,            parre,nparr,parii)
        call PAR_EXCHANGE(prthe_pts,            parre,nparr,parii)
        call PAR_EXCHANGE(mean_free_path_pts,   parre,nparr,parii)           
        call PAR_EXCHANGE(kfl_momentum_pts,     parin,npari,parii)           
        call PAR_EXCHANGE(kfl_thermo_pts,       parin,npari,parii)           
        call PAR_EXCHANGE(kfl_thermo_timsch_pts,parin,npari,parii) 
        call PAR_EXCHANGE(nclas_pts,            parin,npari,parii)           
        call PAR_EXCHANGE(kfl_momentum_sink_pts,parin,npari,parii)           
        call PAR_EXCHANGE(kfl_heat_sink_pts,    parin,npari,parii)           
        call PAR_EXCHANGE(kfl_mass_sink_pts,    parin,npari,parii)           
        !
        ! Type description
        !
        do ii = 1,mtyla
           call PAR_EXCHANGE(parttyp(ii) % kfl_exist,      parin,npari,parii)
           call PAR_EXCHANGE(parttyp(ii) % kfl_modla,      parin,npari,parii)
           call PAR_EXCHANGE(parttyp(ii) % kfl_therm,      parin,npari,parii)
           call PAR_EXCHANGE(parttyp(ii) % kfl_heattr_corr,parin,npari,parii)
           call PAR_EXCHANGE(parttyp(ii) % kfl_mass_pot,   parin,npari,parii)
           call PAR_EXCHANGE(parttyp(ii) % kfl_grafo,      parin,npari,parii)
           call PAR_EXCHANGE(parttyp(ii) % kfl_buofo,      parin,npari,parii)
           call PAR_EXCHANGE(parttyp(ii) % kfl_drafo,      parin,npari,parii)
           call PAR_EXCHANGE(parttyp(ii) % kfl_extfo,      parin,npari,parii)
           call PAR_EXCHANGE(parttyp(ii) % kfl_brown,      parin,npari,parii)
           call PAR_EXCHANGE(parttyp(ii) % kfl_turbu,      parin,npari,parii)
           call PAR_EXCHANGE(parttyp(ii) % kfl_saffm,      parin,npari,parii)
           call PAR_EXCHANGE(parttyp(ii) % kfl_schem,      parin,npari,parii)
           call PAR_EXCHANGE(parttyp(ii) % kfl_tstep,      parin,npari,parii)
           call PAR_EXCHANGE(parttyp(ii) % kfl_dmini,      parin,npari,parii)
           
           call PAR_EXCHANGE(parttyp(ii) % denpa,      parre,nparr,parii)           
           call PAR_EXCHANGE(parttyp(ii) % spher,      parre,nparr,parii)           
           call PAR_EXCHANGE(parttyp(ii) % diame,      parre,nparr,parii)           
           call PAR_EXCHANGE(parttyp(ii) % calor,      parre,nparr,parii)           
           call PAR_EXCHANGE(parttyp(ii) % emisi,      parre,nparr,parii)           
           call PAR_EXCHANGE(parttyp(ii) % scatt,      parre,nparr,parii)           
           call PAR_EXCHANGE(parttyp(ii) % diffu,      parre,nparr,parii)
           call PAR_EXCHANGE(parttyp(ii) % param_dmini,parre,nparr,parii) 
           call PAR_EXCHANGE(parttyp(ii) % n_drop,     parre,nparr,parii) 
           call PAR_EXCHANGE(parttyp(ii) % L_vapor,    parre,nparr,parii) 
           call PAR_EXCHANGE(parttyp(ii) % Cp,         parre,nparr,parii) 
           call PAR_EXCHANGE(parttyp(ii) % w,          parre,nparr,parii) 
           call PAR_EXCHANGE(parttyp(ii) % T_boiling,  parre,nparr,parii) 
           call PAR_EXCHANGE(parttyp(ii) % T_crit,     parre,nparr,parii) 
           call PAR_EXCHANGE(parttyp(ii) % P_boiling,  parre,nparr,parii) 
           call PAR_EXCHANGE(parttyp(ii) % weight_seen,parre,nparr,parii) 
           call PAR_EXCHANGE(parttyp(ii) % dtime,      parre,nparr,parii)     
           call PAR_EXCHANGE(parttyp(ii) % safet,      parre,nparr,parii)     
           call PAR_EXCHANGE(parttyp(ii) % chale,      parre,nparr,parii)     
           call PAR_EXCHANGE(parttyp(ii) % tursc,      parre,nparr,parii)      
           
           call PAR_EXCHANGE(parttyp(ii) % kfl_tab_fw, parin,npari,parii)
           call PAR_EXCHANGE(parttyp(ii) % kfl_spr_fw, parin,npari,parii)

           call PAR_EXCHANGE(parttyp(ii) % liq % name, parch,nparc,parii)
           do jj = 1,6
              do kk = 1,2
                 call PAR_EXCHANGE(parttyp(ii) % cpcoef_v_chm(jj,kk),parre,nparr,parii)
              end do
           end do

           call PAR_EXCHANGE(mlapr,parttyp(ii) % prope,parre,nparr,parii) 
           call PAR_EXCHANGE(mlapr,parttyp(ii) % prova,parin,npari,parii)           
        end do
        !
        ! Numerical problem
        !
        call PAR_EXCHANGE(kfl_adapt_pts,parin,npari,parii)
        call PAR_EXCHANGE(kfl_usbin_pts,parin,npari,parii)
        call PAR_EXCHANGE(kfl_order_pts,parin,npari,parii)
        call PAR_EXCHANGE(kfl_walld_pts,parin,npari,parii)
        call PAR_EXCHANGE(kfl_vesgs_pts,parin,npari,parii)
        
        call PAR_EXCHANGE(gamma_pts,parre,nparr,parii)
        call PAR_EXCHANGE(beta_pts ,parre,nparr,parii)
        call PAR_EXCHANGE(chale_pts,parre,nparr,parii)
        call PAR_EXCHANGE(safet_pts,parre,nparr,parii)
        call PAR_EXCHANGE(dtime_pts,parre,nparr,parii)
        solve_sol => solve(1:)
        call soldef(1_ip)
        !
        ! Output and postprocess
        !                  
        call output_postprocess_parall_old()

        call PAR_EXCHANGE(kfl_max_out_exist_state_pts ,parin,npari,parii)
        call PAR_EXCHANGE(kfl_min_out_exist_state_pts ,parin,npari,parii)


        call PAR_EXCHANGE(kfl_posla_pts        ,parin,npari,parii)
        call PAR_EXCHANGE(kfl_oudep_pts        ,parin,npari,parii)
        call PAR_EXCHANGE(kfl_depos_surface_pts,parin,npari,parii)
        call PAR_EXCHANGE(kfl_oufre_pts        ,parin,npari,parii)
        call PAR_EXCHANGE(kfl_dbfre_pts        ,parin,npari,parii)
        call PAR_EXCHANGE(kfl_exacs_pts        ,parin,npari,parii)
        call PAR_EXCHANGE(kfl_rstar_pts        ,parin,npari,parii)

#ifdef DBPARTICLES           
        call PAR_EXCHANGE(dbSettings % kfl_db_deep      ,parin,npari,parii)
        call PAR_EXCHANGE(dbSettings % maxInserts       ,parin,npari,parii)
        call PAR_EXCHANGE(dbSettings % port             ,parin,npari,parii)
        call PAR_EXCHANGE(dbSettings % numParticlesBatch,parin,npari,parii)
        call PAR_EXCHANGE(dbSettings % kfl_db_url_conn  ,parch,nparc,parii)           
#endif
        call PAR_EXCHANGE(mvarp_pts,postprocess_var_pts,parlo,nparl,parii)
        call PAR_EXCHANGE(mvarp_pts,deposition_var_pts ,parlo,nparl,parii)
        !
        ! Boundary conditions
        !
        call PAR_EXCHANGE(kfl_imax_injector_pts,parin,npari,parii)
        call PAR_EXCHANGE(kfl_boundary_injection,parin,npari,parii) 
        call PAR_EXCHANGE(pts_minj,codbo_pts,parin,npari,parii) 
        call PAR_EXCHANGE(tinla_pts,parre,nparr,parii)
        call PAR_EXCHANGE(tfila_pts,parre,nparr,parii)
        call PAR_EXCHANGE(tpela_pts,parre,nparr,parii)

        if( parii == 1 ) then
           call memory_alloca(mem_modul(1:2,modul),'PARIN','pts_parall',parin,npari)
           call memory_alloca(mem_modul(1:2,modul),'PARRE','pts_parall',parre,nparr)
           call memory_alloca(mem_modul(1:2,modul),'PARLO','pts_parall',parlo,nparl)
           if( ISLAVE  ) call PAR_BROADCAST(parin,      'IN MY CODE')
           if( ISLAVE  ) call PAR_BROADCAST(parre,      'IN MY CODE')
           if( ISLAVE  ) call PAR_BROADCAST(parlo,      'IN MY CODE')
           if( ISLAVE  ) call PAR_BROADCAST(nparc,parch,'IN MY CODE')
        else
           if( IMASTER ) call PAR_BROADCAST(parin,      'IN MY CODE')
           if( IMASTER ) call PAR_BROADCAST(parre,      'IN MY CODE')
           if( IMASTER ) call PAR_BROADCAST(parlo,      'IN MY CODE')
           if( IMASTER ) call PAR_BROADCAST(nparc,parch,'IN MY CODE')
        end if

     end do

     call memory_deallo(mem_modul(1:2,modul),'PARIN','pts_parall',parin)
     call memory_deallo(mem_modul(1:2,modul),'PARRE','pts_parall',parre)
     call memory_deallo(mem_modul(1:2,modul),'PARLO','pts_parall',parlo)

  end select

  npari = 0
  nparr = 0
  nparc = 0
  nparl = 0
  
  !
  ! Flamelet tables
  !
  do ii = 1,mtyla
     if (parttyp(ii) % kfl_therm == 2) then 
        if (parttyp(ii) % kfl_tab_fw > 0) &
           parttyp(ii) % table_fw => lookup_fw(parttyp(ii) % kfl_tab_fw)
        if (parttyp(ii) % kfl_spr_fw > 0) &
           parttyp(ii) % spr_tab_fw => lookup_fw(parttyp(ii) % kfl_spr_fw)
     endif
  enddo


  !
  ! Injector parallelization 
  !
  call pts_injection_parallelization()
  !
  ! Boundary conditions
  !
  call spbbcs(tbcod_pts)


end subroutine pts_parall

