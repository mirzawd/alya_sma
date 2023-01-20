!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine chm_memous(itask)
  !-----------------------------------------------------------------------
  !****f* Chemic/chm_memous
  ! NAME
  !    chm_memous
  ! DESCRIPTION
  !    Allocate memory for the auxiliary variables in the output and post-
  !    processing part
  ! OUTPUT
  ! USES
  ! USED BY
  !    chm_reaous
  !***
  !-----------------------------------------------------------------------
  use def_master,       only : mem_modul, modul
  use def_chemic,       only : nspec_cond_write_CMC_chm, nspec_uncond_write_CMC_chm, nZ_write_CMC_chm, write_cond_spec_CMC_chm,&
                               write_iZ_spec_CMC_chm, write_uncond_spec_CMC_chm
  use def_kintyp,       only : ip
  use mod_memory,       only : memory_alloca
  implicit none
  integer(ip), intent(in) :: itask

  select case(itask)

  case(1_ip)
     !
     ! Allocate vector for <Y> species positions for post-processing for CMC model
     !

     nullify(write_uncond_spec_CMC_chm)
     call memory_alloca(mem_modul(1:2,modul),'WRITE_UNCOND_SPECIES_POS','chm_memous',write_uncond_spec_CMC_chm,&
         nspec_uncond_write_CMC_chm)
     write_uncond_spec_CMC_chm = 0_ip


  case(2_ip)
     !
     ! Allocate vector for <Y|Z> species for post-processing for CMC model
     !
     nullify(write_cond_spec_CMC_chm)
     call memory_alloca(mem_modul(1:2,modul),'WRITE_COND_SPECIES_POS','chm_memous',write_cond_spec_CMC_chm,nspec_cond_write_CMC_chm)
     write_cond_spec_CMC_chm = 0_ip

  case(3_ip)
     !
     ! Allocate vector for mixture fractions positions for post-processing for CMC model
     !
     nullify(write_iZ_spec_CMC_chm)
     call memory_alloca(mem_modul(1:2,modul),'WRITE_MIXT_FRAC_POS','chm_memous',write_iZ_spec_CMC_chm,nZ_write_CMC_chm)
     write_iZ_spec_CMC_chm = 0_ip

  end select

end subroutine chm_memous

