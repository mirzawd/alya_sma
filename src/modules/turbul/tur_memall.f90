!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tur_memall()
  !-----------------------------------------------------------------------
  !****f* Turbul/tur_memall
  ! NAME 
  !    tur_memall
  ! DESCRIPTION
  !    This routine allocates memory for the arrays needed to solve the
  !    turbulence equations      
  ! USES
  ! USED BY
  !    tur_turnon
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_inpout
  use def_master
  use def_domain
  use def_solver
  use def_turbul
  use mod_memory
  use mod_tur_arrays, only : tur_arrays
  use def_kermod, only : kfl_adj_prob,kfl_ndvars_opt
  use mod_ADR, only : FULL_OSS
  use mod_ADR, only : A_OSS  
  use mod_ADR, only : AR_OSS 
  use mod_ADR, only : BUBBLE
  use mod_ADR, only : ADR_initialize_type
  use mod_ADR, only : ADR_check_and_compute_data
  use mod_ADR, only : ADR_allocate_projections_bubble_sgs 
  use mod_output_postprocess, only : output_postprocess_check_variable_postprocess
  use mod_arrays,             only : arrays_number
  implicit none
  integer(ip) :: ielem,pelty,pnode

  call tur_arrays('ALLOCATE')
  !
  ! DUNKN_TUR: Delta unknown for Aitken relaxation strategy
  !
  if( kfl_relax_tur == 2 ) then
     call memory_alloca(mem_modul(1:2,modul),'DUNKN_TUR','tur_memall',dunkn_tur,nturb_tur*npoin)        
  end if
  !
  ! UNOLD_TUR: Old solution
  ! 
  if(    output_postprocess_check_variable_postprocess(arrays_number('RESI1')) .or.&
       & output_postprocess_check_variable_postprocess(arrays_number('RESI2')) .or.&
       & output_postprocess_check_variable_postprocess(arrays_number('RESIT')) ) then
     call memory_alloca(mem_modul(1:2,modul),'UNOLD_TUR','tur_memall',unold_tur,nturb_tur+1_ip,npoin)
  end if
  !
  ! UNPGR_TUR: Residual projections for Shock capturing methos
  !
  if (kfl_shock_tur/=0)  call memory_alloca(mem_modul(1:2,modul),'UNPGR_TUR','tur_memall',unpgr_tur,nturb_tur,ndime, npoin)

  !
  ! DETUR_TUR, VITUR_TUR: projected variable density and viscosity
  !
  if( kfl_colev_tur >= 1 ) then
     call memory_alloca(mem_modul(1:2,modul),'DETUR_TUR','tur_memall',detur_tur,npoin)
     call memory_alloca(mem_modul(1:2,modul),'VITUR_TUR','tur_memall',vitur_tur,npoin)
  end if
  !
  ! only for DDES, postprocessing value
  !
  if( kfl_ddesm_tur >= 1 ) then
     call memory_alloca(mem_modul(1:2,modul),'FDDES_TUR','tur_memall',fddes_tur,npoin)
     call memory_alloca(mem_modul(1:2,modul),'GDDES_TUR','tur_memall',gddes_tur,npoin)
  end if
  !
  ! only for SST, postprocessing values
  !
  if( TUR_SST_K_OMEGA) then
     call memory_alloca(mem_modul(1:2,modul),'SSTF1_TUR','tur_memall',sstf1_tur,npoin)
     call memory_alloca(mem_modul(1:2,modul),'SSTF2_TUR','tur_memall',sstf2_tur,npoin)
  end if
  if( TUR_SST_K_OMEGA .and. kfl_sasim_tur == 1) then
     call memory_alloca(mem_modul(1:2,modul),'SASSO_TUR','tur_memall',sasso_tur,npoin)
  end if
  if (TUR_FAMILY_K_EPS) &
       call memory_alloca(mem_modul(1:2,modul),'TUR_MAX_MIXLEN','tur_memall',tur_max_mixlen,npoin)
  !
  ! PRODU_TUR: production term
  !
  if( kfl_produ_tur == 1 ) then
     call memory_alloca(mem_modul(1:2,modul),'PRODU_TUR','tur_memall',produ_tur,npoin)
  end if
  ! provisional variable (MATIAS)
  call memory_alloca(mem_modul(1:2,modul),'TURVI_TUR','tur_memall',turvi_tur,2_ip, mgaus, nelem)
  !
  ! adjoint materials
  !
  if(kfl_adj_prob == 1 ) then     
     call memory_alloca(mem_modul(1:2,modul),'UNTUR_FORW','tur_memall',untur_forw,nturb_tur,npoin,ncomp_tur)
     call memory_alloca(mem_modul(1:2,modul),'RESDIFF_TUR','tur_memall',resdiff_tur,kfl_ndvars_opt,npoin)

     if (nturb_tur == 1) then
        call memory_alloca(mem_modul(1:2,modul),'Rhsadjtur_tur','tur_memall',Rhsadjtur_tur,nelem)
        do ielem=1,nelem
           pelty = ltype(ielem)
           pnode = nnode(pelty)
           call memory_alloca(mem_modul(1:2,modul),'Rhsadjtur_tur','tur_memall',Rhsadjtur_tur(ielem)%a,nturb_tur,pnode)
        end do
     endif

     call memory_alloca(mem_modul(1:2,modul),'RhsadjNas_tur','tur_memall',RhsadjNas_tur,nelem)

     do ielem=1,nelem
        pelty = ltype(ielem)
        pnode = nnode(pelty)
        call memory_alloca(mem_modul(1:2,modul),'RhsadjNas_tur','tur_memall',RhsadjNas_tur(ielem)%a,ndime,pnode)
     end do

  endif !kfl_adj_prob
 
end subroutine tur_memall
 
