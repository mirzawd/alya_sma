!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine par_sendat(order)
  !-----------------------------------------------------------------------
  !****f* Parall/par_sendat
  ! NAME
  !    par_sendat
  ! DESCRIPTION
  !    This routine exchange data
  ! USES
  ! USED BY
  !    Reapro
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_parall
  use def_parall,                 only : kfl_openacc_policy
  use def_parall,                 only : kfl_openacc_ranks_per_dev
  use def_parall,                 only : kfl_openacc_contiguous_ranks_per_dev
  use def_parall,                 only : kfl_openacc_streams_per_dev
  use def_parall,                 only : kfl_openacc_multithreading
  use def_inpout
  use def_solver
  use mod_memory
  use mod_parall
  use def_kermod
  use def_mpio
  use mod_messages
  use mod_iofile
  use mod_optimum_partition,      only : optimum_partition_parall
#ifdef ALYA_FTI
  use mod_alya2fti,               only : alya2fti_iexcha
#endif
  use mod_output_postprocess,     only : output_postprocess_parall_old
  implicit none
  integer(ip), intent(in) :: order
  integer(ip)             :: ii,ir,ic,ji,ki
  integer(ip)             :: ipass,kfl_ptask_old

  external                :: iexcha, rexcha, cexcha, lexcha

  select case (order)

  case(1_ip)
     !
     ! First communication always performed with MPI and not files
     !
     if( IMASTER .and. kfl_ptask == 0 .and. nproc_par == 1 ) return ! not need to broadcast
     kfl_ptask_old = kfl_ptask
     kfl_ptask     = 1
     call vocabu(-1_ip,0_ip,0_ip)
     !
     ! Exchange data read in Reapro (always using MPI)
     !
     if( IPARALL ) then

        strre = 'Reapro'
        strin = 'Reapro'
        strch = 'Reapro'
        do parii = 1,2
           npari = 0
           nparr = 0
           nparc = 0
           !
           ! Partition (just to check errors)
           !
           call iexcha(npart_par)
           !
           ! Read in rrudat
           !
           call iexcha(current_code)
           call iexcha(kfl_naked)
           call iexcha(kfl_rstar)
           call iexcha(kfl_ptask_old)
           call iexcha(npart_empty_par)
           call iexcha(par_omp_granularity)
           call iexcha(par_omp_coloring_alg)
           call iexcha(par_omp_nelem_chunk)
           call iexcha(par_omp_nboun_chunk)
           call iexcha(par_omp_npoin_chunk)
           call iexcha(par_omp_partition_alg)

           call iexcha(par_one_sided)
           call iexcha(par_topo_num_nodes)
           call iexcha(par_topo_num_cores_per_node)

           call iexcha(kfl_matri_par)

           call iexcha(kfl_partition_par)
           call iexcha(kfl_parseq_par)
           call iexcha(kfl_interface_parti_par)
           call iexcha(kfl_weigh_par)
           call iexcha(kfl_repart_par)
           call iexcha(kfl_repart_module_par)
           call iexcha(kfl_repart_windo_par)
           call iexcha(kfl_repart_criterion_par)
           call iexcha(kfl_repart_method_par)
           call iexcha(kfl_repart_freq_par)
           call iexcha(kfl_repart_post_par)
           call iexcha(kfl_parti_par)
           call iexcha(kfl_global_numbering_par)
           call iexcha(kfl_order_exchange_par)
           call iexcha(kfl_connectivity_par)

           call iexcha(kfl_openacc_policy)
           call iexcha(kfl_openacc_ranks_per_dev)
           call iexcha(kfl_openacc_contiguous_ranks_per_dev)
           call iexcha(kfl_openacc_streams_per_dev)
           call lexcha(kfl_openacc_multithreading)

           do ji = 1,size(weights_elements_par,KIND=ip)
              call rexcha(weights_elements_par(ji))
           end do
           do ji = 1,size(weights_materials_par,KIND=ip)
              call rexcha(weights_materials_par(ji))
           end do
           call rexcha(repart_threshold_lb)           
           call rexcha(repart_toler)           
           call rexcha(repart_min_weight)           
           call rexcha(repart_max_weight)           
           call rexcha(repart_wfact)           
           do ji = 1,3
              call iexcha(boxes_coarse_par(ji))
           end do
           do ji = 1,3
              call iexcha(boxes_fine_par(ji))
           end do
           do ji = 1,3
              call rexcha(vect_partition_par(ji))
           end do

           call iexcha(sfc_check)
           call iexcha(sfc_criteria)
           call iexcha(sfc_dim_bin_core)

           nparc=nparc+132
           if(parii==2.and.IMASTER) parch(1:66)   = title(1:66)
           if(parii==2.and.ISLAVE)  title(1:66)   = parch(1:66)
           if(parii==2.and.IMASTER) parch(67:132) = namda(1:66)
           if(parii==2.and.ISLAVE)  namda(1:66)   = parch(67:132)
           call cexcha(int(len(method_redistribution_par),ip),method_redistribution_par)

           call iexcha(kfl_output_partition)
           call iexcha(kfl_output_node_comm_arrays)
           call iexcha(kfl_output_edge_comm_arrays)
           !
           ! Read in readat and 'mod'_reapro of modules 'mod'
           !
           do ji=1,mblok
              call iexcha(micou(ji))
           end do
           call iexcha(nblok)
           call iexcha(kfl_timco)
           call iexcha(kfl_timei)
           call iexcha(kfl_timef)
           call iexcha(kfl_dtfun)
           call iexcha(mitim)
           call iexcha(mitsm)
           call iexcha(mitrf)
           call iexcha(kfl_postp_par)
           call iexcha(kfl_wwork)
           call iexcha(kfl_lumped)
           call rexcha(timei)
           call rexcha(timef)
           call rexcha(dtime)
           do ji=1,mmodu
              call iexcha(kfl_modul(ji))
           end do
           do ji=1,mmodu
              call iexcha(kfl_delay(ji))
           end do
           do ji=1,mmodu
              call iexcha(kfl_conve(ji))
           end do
           do ji=0,mmodu
              call iexcha(kfl_solve(ji))
           end do
           do ji=1,mmodu
              call iexcha(ndela(ji))
           end do
           do ji=1,mmodu
              do ki=1,mblok
                 call iexcha(lmord(ji,ki))
              end do
           end do
           do ji=0,mmodu
              call iexcha(lzone(ji))
           end do
           call output_postprocess_parall_old()
#ifdef ALYA_FTI
           call alya2fti_iexcha()
#endif
           !
           ! Automatic partitioning
           !
           call optimum_partition_parall()
           !
           ! Allocate memory for the first pass
           !
           if(parii==1) then
              call memory_alloca(par_memor,'PARIN','par_sendat',parin,npari)
              call memory_alloca(par_memor,'PARRE','par_sendat',parre,nparr)
              if(ISLAVE) call par_broadc()
           end if
        end do
        if(IMASTER) call par_broadc()

        call memory_deallo(par_memor,'PARIN','par_sendat',parin)
        call memory_deallo(par_memor,'PARRE','par_sendat',parre)

        kfl_ptask = kfl_ptask_old
        call vocabu(-1_ip,0_ip,0_ip)
     
     end if
     
  case(2_ip)

     if( IPARALL ) then

        kfl_ptask_old = kfl_ptask
        if( PART_AND_WRITE() .and. nproc_par > 1 ) kfl_ptask = 1

        strre='readim_reastr_reageo'
        strin='readim_reastr_reageo'
        strch='readim_reastr_reageo'
        do parii=1,2
           ipass = parii
           npari=0
           nparr=0
           nparc=0
           !
           ! Read in reageo
           !
           continue
           !
           ! Read in reaset
           !
           call iexcha(neset)
           call iexcha(nbset)
           call iexcha(nnset)           ! Re-computed in mod_reaset
           call iexcha(neset_origi)
           call iexcha(nbset_origi)
           call iexcha(nnset_origi)     ! Re-computed in mod_reaset
           !
           ! Read in reabcs
           !
           call iexcha(kfl_icodn)
           call iexcha(kfl_icodb)
           !
           ! Read in reafie
           !
           continue

           if( ipass == 1 ) then
              call memory_alloca(par_memor,'PARIN','par_sendat',parin,npari)
              call memory_alloca(par_memor,'PARRE','par_sendat',parre,nparr)
              if( ISLAVE ) call par_receiv()
           end if
        end do

        if( IMASTER ) call par_sendin()

        call memory_deallo(par_memor,'PARIN','par_sendat',parin)
        call memory_deallo(par_memor,'PARRE','par_sendat',parre)

        kfl_ptask = kfl_ptask_old
        
     end if

  case(5)

     if( kfl_paral >= 0 ) then

        strre='par_partit'
        strin='par_partit'
        strch='par_partit'
        do ipass=1,2
           ii=0
           ir=0
           ic=0
           !
           ! Calculated in partit
           !
           ii=ii+1
           if(ipass==2.and.IMASTER) parin(ii)   = ginde_par(3,kfl_desti_par)
           if(ipass==2.and.ISLAVE)  lni         = parin(ii)
           ii=ii+1
           if(ipass==2.and.IMASTER) parin(ii)   = ginde_par(4,kfl_desti_par)
           if(ipass==2.and.ISLAVE)  lnb         = parin(ii)
           ii=ii+1
           if(ipass==2.and.IMASTER) parin(ii)   = lneig_par(kfl_desti_par)
           if(ipass==2.and.ISLAVE)  nneig       = parin(ii)
           ii=ii+1
           if(ipass==2.and.IMASTER) parin(ii)   = slfbo_par(kfl_desti_par)
           if(ipass==2.and.ISLAVE)  slfbo       = parin(ii)
           if(ipass==1) then
              npari = ii
              nparr = ir
              nparc = ic
              call memory_alloca(par_memor,'PARIN','par_sendat',parin,npari)
              call memory_alloca(par_memor,'PARRE','par_sendat',parre,nparr)
              if(ISLAVE) call par_receiv()
           end if
        end do
        if(IMASTER) call par_sendin()

        call memory_deallo(par_memor,'PARIN','par_sendat',parin)
        call memory_deallo(par_memor,'PARRE','par_sendat',parre)

     end if

  case(7_ip)
     !
     ! Boundary conditions and solvers
     !
     if( IPARALL ) then

        strre='cderda_reaset'
        strin='cderda_reaset'
        strch='cderda_reaset'
        do parii=1,2
           npari=0
           nparr=0
           nparc=0
           
           call iexcha(kfl_schur)
           call iexcha(kfl_aiipr)
           
           if(parii==1) then
              call memory_alloca(par_memor,'PARIN','par_sendat',parin,npari)
              call memory_alloca(par_memor,'PARRE','par_sendat',parre,nparr)
              if( ISLAVE .or. kfl_ptask == 2 ) call par_broadc()
           end if
        end do

        if( IMASTER .and. kfl_ptask /= 2 ) call par_broadc()

        call memory_deallo(par_memor,'PARIN','par_sendat',parin)
        call memory_deallo(par_memor,'PARRE','par_sendat',parre)

     end if

  case(8_ip)

     if( IPARALL ) then

        strre='readim_reastr_reageo'
        strin='readim_reastr_reageo'
        strch='readim_reastr_reageo'
        do parii=1,2
           ipass = parii
           npari=0
           nparr=0
           nparc=0
           !
           ! Read in readim
           !
           npari=npari+1
           if(ipass==2.and.IMASTER) parin(npari) = npoin_par(kfl_desti_par)
           if(ipass==2.and.ISLAVE)  npoin        = parin(npari)
           npari=npari+1
           if(ipass==2.and.IMASTER) parin(npari) = nelem_par(kfl_desti_par)
           if(ipass==2.and.ISLAVE)  nelem        = parin(npari)
           npari=npari+1
           if(ipass==2.and.IMASTER) parin(npari) = nboun_par(kfl_desti_par)
           if(ipass==2.and.ISLAVE)  nboun        = parin(npari)

           if( ipass == 1 ) then
              call memory_alloca(par_memor,'PARIN','par_sendat',parin,npari)
              call memory_alloca(par_memor,'PARRE','par_sendat',parre,nparr)
              if( ISLAVE ) call par_receiv()
           end if
        end do

        if( IMASTER ) call par_sendin()

        call memory_deallo(par_memor,'PARIN','par_sendat',parin)
        call memory_deallo(par_memor,'PARRE','par_sendat',parre)

     end if

  end select

  npari=0
  nparr=0
  nparc=0

end subroutine par_sendat

