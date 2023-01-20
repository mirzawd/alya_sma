!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine chm_memnut(itask)
  !-----------------------------------------------------------------------
  !****f* Chemic/chm_memnut
  ! NAME
  !    chm_memous
  ! DESCRIPTION
  !    Allocate memory for the auxiliary variables in the numerical
  !    treatment part
  ! OUTPUT
  ! USES
  ! USED BY
  !    chm_reanut
  !***
  !-----------------------------------------------------------------------
  use def_master,       only : mem_modul, modul
  use def_chemic,       only : chem_int_iZ_CMC_chm, nspec_transf_CMC_chm, nZ_chm_int_CMC_chm, transf_spec_CMC_chm
  use def_kintyp,       only : ip
  use mod_memory,       only : memory_alloca
  implicit none
  integer(ip), intent(in) :: itask

  select case(itask)


  case(1_ip)
     !
     ! Allocate vector for mixture fractions positions for chemical integration for CMC model
     !
     nullify(chem_int_iZ_CMC_chm)
     call memory_alloca(mem_modul(1:2,modul),'CHEM_INT_MIXT_FRAC_POS','chm_memous',chem_int_iZ_CMC_chm,nZ_chm_int_CMC_chm)
     chem_int_iZ_CMC_chm = 0_ip

  case(2_ip)
     !
     ! Allocate vector for species to be transferred
     !
     nullify(transf_spec_CMC_chm)
     call memory_alloca(mem_modul(1:2,modul),'CHEM_SPEC_TRANSF','chm_memous',transf_spec_CMC_chm,nspec_transf_CMC_chm)
     transf_spec_CMC_chm = 0_ip
  end select

end subroutine chm_memnut

