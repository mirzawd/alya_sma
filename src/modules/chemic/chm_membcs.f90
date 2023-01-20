!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine chm_membcs(itask)
  !-----------------------------------------------------------------------
  !****f* Temper/chm_membcs
  ! NAME
  !    chm_membcs
  ! DESCRIPTION
  !    Allocate memory for the physical problem
  ! OUTPUT
  ! USES
  ! USED BY
  !    chm_reaphy
  !***
  !-----------------------------------------------------------------------
  use def_master,         only : kfl_rstar, mem_modul, modul
  use def_domain,         only : nboun, npoin
  use def_chemic,         only : bvess_chm, bvess_CMC_chm, bvess_ufield_CMC_chm, kfl_bc_alpha_CMC_chm, kfl_fixbo_chm,&
                                 kfl_fixno_chm, kfl_funno_chm, kfl_funtn_chm, kfl_initi_chm, kfl_model_chm, kfl_solve_cond_CMC_chm,&
                                 kfl_solve_enth_CMC_chm, kfl_weigh_in_eq_CMC_chm, nclas_chm, nvar_CMC_chm, xinit_chm, nZ_CMC_chm
  use def_kintyp,         only : ip, rp
  use def_parame,         only : two, zero
  use mod_memchk,         only : memchk
  use mod_memory,         only : memory_alloca
  implicit none
  integer(ip), intent(in) :: itask
  integer(4)              :: istat
  integer(ip)             :: nequs_loc, aux

  external                :: memerr

  !
  ! Number of equations solved
  !
  nequs_loc = nclas_chm
  if (kfl_model_chm == 4 .and. kfl_solve_cond_CMC_chm == 1)  nequs_loc = nvar_CMC_chm


  select case(itask)

  case(1_ip)
     !
     ! Fixity and boundary values
     !
     call memory_alloca(mem_modul(1:2,modul),'KFL_FIXNO_CHM','chm_membcs',kfl_fixno_chm,nequs_loc,npoin)
     call memory_alloca(mem_modul(1:2,modul),'KFL_FIXBO_CHM','chm_membcs',kfl_fixbo_chm,nequs_loc,nboun)

     if (kfl_model_chm == 4 .and. kfl_solve_cond_CMC_chm == 1 .and. kfl_rstar == 0 .and. kfl_weigh_in_eq_CMC_chm == 1) then
        call memory_alloca(mem_modul(1:2,modul),'BVESS_CHM','chm_membcs',bvess_chm,nequs_loc+1_ip,npoin)
     else
        call memory_alloca(mem_modul(1:2,modul),'BVESS_CHM','chm_membcs',bvess_chm,nequs_loc,npoin)
     endif


     if (kfl_model_chm == 4 .and. kfl_solve_cond_CMC_chm == 1) then
        allocate(bvess_CMC_chm(nZ_CMC_chm,npoin,nvar_CMC_chm),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'BVESS_CMC_CHM','chm_membcs',bvess_CMC_chm)

        if (kfl_bc_alpha_CMC_chm == 1_ip) then
           if (kfl_solve_enth_CMC_chm == 0_ip) then
              aux = 1_ip
           else
              aux = 2_ip
           end if
        else
           aux = nvar_CMC_chm
        end if
        allocate(bvess_ufield_CMC_chm(npoin,aux),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'BVESS_UNCOND_FIELD_CMC_CHM','chm_membcs',bvess_ufield_CMC_chm)
     end if

     if (npoin > 0) &
        kfl_fixno_chm(:,:)           = 0_ip
     if (nboun > 0) &
        kfl_fixbo_chm(:,:)           = 0_ip
     if (npoin > 0) &
        bvess_chm(:,:)               = 0.0_rp
     if (kfl_model_chm == 4 .and. kfl_solve_cond_CMC_chm == 1) then
        bvess_CMC_chm        = 0.0_rp
        bvess_ufield_CMC_chm = 0.0_rp
     end if

  case(2_ip)
     !
     ! Deallocate memory
     !
     call memchk(two,istat,mem_modul(1:2,modul),'KFL_FIXNO_CHM','chm_membcs',kfl_fixno_chm)
     deallocate(kfl_fixno_chm,stat=istat)
     if(istat/=0) call memerr(two,'KFL_FIXNO_CHM','nsi_membcs',0_ip)

     call memchk(two,istat,mem_modul(1:2,modul),'KFL_FIXBO_CHM','chm_membcs',kfl_fixbo_chm)
     deallocate(kfl_fixbo_chm,stat=istat)
     if(istat/=0) call memerr(two,'KFL_FIXBO_CHM','nsi_membcs',0_ip)

     call memchk(two,istat,mem_modul(1:2,modul),'BVESS_CHM','chm_membcs',bvess_chm)
     deallocate(bvess_chm,stat=istat)
     if(istat/=0) call memerr(two,'BVESS_CHM','nsi_membcs',0_ip)

     if (kfl_model_chm == 4 .and. kfl_solve_cond_CMC_chm == 1) then
        call memchk(two,istat,mem_modul(1:2,modul),'BVESS_CMC_CHM','chm_membcs',bvess_CMC_chm)
        deallocate(bvess_CMC_chm,stat=istat)
        if(istat/=0) call memerr(two,'BVESS_CMC_CHM','nsi_membcs',0_ip)

        call memchk(two,istat,mem_modul(1:2,modul),'BVESS_UNCOND_FIELD_CMC_CHM','chm_membcs',bvess_ufield_CMC_chm)
        deallocate(bvess_ufield_CMC_chm,stat=istat)
        if(istat/=0) call memerr(two,'BVESS_UNCOND_FIELD_CMC_CHM','nsi_membcs',0_ip)
     end if

  case(3_ip)
     !
     ! Allocate parameters
     !
     call memory_alloca(mem_modul(1:2,modul),'KFL_INITI_CHM','chm_membcs',kfl_initi_chm,nclas_chm)
     call memory_alloca(mem_modul(1:2,modul),'XINIT_CHM','chm_membcs',xinit_chm,nclas_chm,2_ip)

  case(-3_ip)
     !
     ! Deallocate parameters
     !
     call memchk(two,istat,mem_modul(1:2,modul),'XINIT_CHM','chm_membcs',xinit_chm)
     deallocate(xinit_chm,stat=istat)
     if(istat/=0) call memerr(two,'XINIT_CHM','nsi_membcs',0_ip)

     call memchk(two,istat,mem_modul(1:2,modul),'KFL_INITI_CHM','chm_membcs',kfl_initi_chm)
     deallocate(kfl_initi_chm,stat=istat)
     if(istat/=0) call memerr(two,'KFL_INITI_CHM','nsi_membcs',0_ip)

  case(4_ip)

     !
     ! Non-constant b.c.'s : Functions
     !
     call memory_alloca(mem_modul(1:2,modul),'KFL_FUNNO_CHM','chm_membcs',kfl_funno_chm,npoin,nequs_loc)
     call memory_alloca(mem_modul(1:2,modul),'KFL_FUNTN_CHM','chm_membcs',kfl_funtn_chm,npoin,nequs_loc)

  end select

end subroutine chm_membcs
