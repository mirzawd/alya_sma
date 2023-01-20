!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Nastin
!> @{
!> @file    nsi_memall.f90
!> @date    06/06/1966
!> @author  Guillaume Houzeaux
!> @brief   Allocate memory
!> @details Allocate memory. Here are allocated variables which should
!>          redistributed. For example, main variables like the velocity
!>          but also postprcess variables like averaging.
!> @}
!-----------------------------------------------------------------------
subroutine nsi_memall()
  use def_parame
  use def_inpout
  use def_master
  use def_domain
  use def_solver
  use def_nastin
  use mod_nsi_arrays
  use mod_memory
  use def_kermod,              only : kfl_adj_prob
  use def_kermod,              only : kfl_ndvars_opt
  use def_kermod,              only : kfl_dvar_type
  use mod_communications,      only : PAR_MIN
  use mod_communications,      only : PAR_MAX
  use mod_output_postprocess,  only : output_postprocess_check_variable_postprocess
  use mod_arrays,              only : arrays
  implicit none
  integer(ip) :: ielem,pelty,pnode
  integer(ip) :: max_code,iboun,ncomp_loc

  !----------------------------------------------------------------------
  !
  ! Arrays
  !
  !----------------------------------------------------------------------

  call nsi_arrays('ALLOCATE')
  
  if (NSI_FRACTIONAL_STEP.and.kfl_stabi_nsi == NSI_ASGS) then
     ncomp_loc = ncomp_nsi+1
  else
     ncomp_loc = ncomp_nsi
  end if

  if( INOTMASTER ) then
     !
     ! Optimization materials for adjoint solution
     !
     if( kfl_adj_prob == 1 ) then
        if( kfl_dvar_type == 5 ) kfl_ndvars_opt = ndime*npoin
        call memory_alloca(mem_modul(1:2,modul),'VELOC_FORW','nsi_memall',veloc_forw,ndime,npoin,ncomp_nsi)
        call memory_alloca(mem_modul(1:2,modul),'PRESS_FORW','nsi_memall',press_forw,npoin,ncomp_nsi)
        ! call memory_alloca(mem_modul(1:2,modul),'RESDIFF_NSI','nsi_memall',resdiff_nsi,kfl_ndvars_opt, nzrhs)
        call memory_alloca(mem_modul(1:2,modul),'DCOST_DX_NSI','nsi_memall',dcost_dx_nsi,ndime*npoin)
        ! Coupling with temper
        if(kfl_coupl(ID_TEMPER,ID_NASTIN) == 1 ) then
           call memory_alloca(mem_modul(1:2,modul),'RHSADJTEM_NSI','nsi_memall',RhsadjTem_nsi,nelem)
           do ielem=1,nelem
              pelty = ltype(ielem)
              pnode = nnode(pelty)
              call memory_alloca(mem_modul(1:2,modul),'RHSADJTEM_NSI % A','nsi_memall',RhsadjTem_nsi(ielem)%a,pnode)
           end do
        end if
        ! Coupling with turbul
        nturb = 2_ip
        if( kfl_coupl(ID_TURBUL,ID_NASTIN) == 1 ) then
           call memory_alloca(mem_modul(1:2,modul),'RHSADJTUR_NSI','nsi_memall',RhsadjTur_nsi,nelem)
           do ielem=1,nelem
              pelty = ltype(ielem)
              pnode = nnode(pelty)
              call memory_alloca(mem_modul(1:2,modul),'RHSADJTUR_NSI % A','nsi_memall',RhsadjTur_nsi(ielem)%a,nturb,pnode)
           end do
        end if
     end if !adjoint
     !
     ! Coriolis stabilisation
     !
     call memory_alloca(mem_modul(1:2,modul),'PRDIVCOR_NSI','nsi_memall',prdivcor_nsi,npoin)

  else
     !
     ! Master: allocate minimum memory
     !
     call memory_alloca(mem_modul(1:2,modul),'VELOC','nsi_memall',veloc,ndime,1_ip,ncomp_loc)
     call memory_alloca(mem_modul(1:2,modul),'PRESS','nsi_memall',press,1_ip,ncomp_nsi)
     if( kfl_stabi_nsi > 0 ) then
        call memory_alloca(mem_modul(1:2,modul),'VEPRO_NSI','nsi_memall',vepro_nsi,1_ip,1_ip)
        call memory_alloca(mem_modul(1:2,modul),'PRPRO_NSI','nsi_memall',prpro_nsi,1_ip)
        if( kfl_stabi_nsi == 2 ) then
           call memory_alloca(mem_modul(1:2,modul),'GRPRO_NSI','nsi_memall',grpro_nsi,1_ip,1_ip)
        end if
     end if
     if(kfl_relax_nsi==2) then
        call memory_alloca(mem_modul(1:2,modul),'DUNKN_NSI','nsi_memall',dunkn_nsi,1_ip,1_ip)
     end if
     if(kfl_relap_nsi==2) then
        call memory_alloca(mem_modul(1:2,modul),'DUNKP_NSI','nsi_memall',dunkp_nsi,1_ip)
     end if
     if(kfl_regim_nsi==1.or.kfl_regim_nsi==2) then
        call memory_alloca(mem_modul(1:2,modul),'DENSI','nsi_memall',densi,1_ip,1_ip)
     end if
     !
     ! DRHODT_NSI: Time derivative of density  
     !
     if(kfl_regim_nsi == 3) then
        call memory_alloca(mem_modul(1:2,modul),'DRHODT_NSI','nsi_memall',drhodt_nsi,1_ip)
     end if
     !
     ! Optimization materials for adjoint solution
     !
     if( kfl_adj_prob == 1 ) then
        if( kfl_dvar_type == 5 ) kfl_ndvars_opt = ndime*npoin
        call memory_alloca(mem_modul(1:2,modul),'DCOST_DX_NSI','nsi_memall',dcost_dx_nsi,1)
     end if !adjoint
     !
     ! Coriolis Stab
     !
     call memory_alloca(mem_modul(1:2,modul),'PRDIVCOR_NSI','nsi_memall',prdivcor_nsi,1)

  end if
  !
  ! UNK2N_NSI, Nastin second variable: PRESS or DENSI
  !
  if( kfl_regim_nsi == 2 ) then
     unk2n_nsi => densi
  else
     unk2n_nsi => press
  end if

  !----------------------------------------------------------------------
  !
  ! For condition of 20 type, need to know how much outflows we have
  !
  !----------------------------------------------------------------------

  if( kfl_exist_fib20_nsi /= 0 ) then
     max_code = -huge(0_ip)
     if (INOTMASTER) then
        do iboun = 1,nboun
           if( kfl_fixbo_nsi(iboun) == 20 ) max_code = max(max_code,int(bvnat_nsi(4,iboun,1),ip))
        end do
     end if
     call PAR_MAX(max_code)
     if( max_code < 1 ) call runend('NSI_MEMALL: WRONG DEFINITION OF STABLE OUTFLOW CONDITION')
     call memory_alloca(mem_modul(1:2,modul),'OUTFLOW_MASS','nsi_memall',outflow_mass,max_code)
  end if

end subroutine nsi_memall
