!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tur_membcs(itask)
  !-----------------------------------------------------------------------
  !****f* Turbul/tur_membcs
  ! NAME 
  !    tur_membcs
  ! DESCRIPTION
  !    This routine reads TURBUL boundary conditions
  ! USES
  ! USED BY
  !    tur_turnon
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_turbul
  use mod_memchk
  use mod_memory
  implicit none 
  integer(ip), intent(in) :: itask
!  integer(4)              :: istat

  select case(itask)

  case(1_ip)
     !
     ! Allocate memory
     !
     call memory_alloca(mem_modul(1:2,modul),'KFL_FIXNO_TUR','tur_membcs',kfl_fixno_tur,1_ip,npoin,nturb_tur)
     call memory_alloca(mem_modul(1:2,modul),'BVESS_TUR'    ,'tur_membcs',bvess_tur,1_ip,npoin,nturb_tur)
     call memory_alloca(mem_modul(1:2,modul),'KFL_FUNNO_TUR','tur_membcs',kfl_funno_tur,npoin,nturb_tur)
     call memory_alloca(mem_modul(1:2,modul),'KFL_FUNTN_TUR','tur_membcs',kfl_funtn_tur,npoin,nturb_tur)
    
     !
     ! Conditions on boundaries
     !    
     call memory_alloca(mem_modul(1:2,modul),'KFL_FIXBO_TUR','tur_membcs',kfl_fixbo_tur,nboun)
     call memory_alloca(mem_modul(1:2,modul),'BVNAT_TUR','tur_membcs',bvnat_tur,npnat_tur,nboun,nturb_tur)

!  only case 1 is beeing called the erst are no longer beeing used for the moment I comment them

!  case(2_ip)
!     !
!     ! Deallocate memory
!     !
!     call memchk(two,istat,mem_modul(1:2,modul),'KFL_FIXNO_TUR','tur_membcs',kfl_fixno_tur)
!     deallocate(kfl_fixno_tur,stat=istat)
!     if(istat/=0) call memerr(two,'KFL_FIXNO_TUR','nsi_membcs',0_ip)
!     call memchk(two,istat,mem_modul(1:2,modul),'BVESS_TUR','tem_membcs',bvess_tur)
!     deallocate(bvess_tur,stat=istat)
!     if(istat/=0) call memerr(two,'BVESS_TUR','tur_membcs',0_ip)
!
!  case(3_ip)
!     !
!     ! WALLD_TUR: Wall distance
!     !
!     allocate(walld_tur(npoin),stat=istat)
!     call memchk(zero,istat,mem_modul(1:2,modul),'WALLD_TUR','tur_membcs',walld_tur)        call memory_alloca(mem_modul(1:2,modul),'VEOLD','nsi_memall',veold_nsi,ndime,npoin)

!
!  case(4_ip)
!     !
!     ! LWNEI_TUR: Nearest neighbor
!     !
!     allocate(lwnei_tur(npoin),stat=istat)
!     call memchk(zero,istat,mem_modul(1:2,modul),'LWNEI_TUR','tur_addarr',lwnei_tur)
!
  case(5_ip)
     !
     !
     !

  end select

end subroutine tur_membcs
