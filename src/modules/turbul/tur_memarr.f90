!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tur_memarr(itask)
  !-----------------------------------------------------------------------
  !****f* Turbul/tur_memarr
  ! NAME 
  !    tur_memarr
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
  implicit none
  integer(ip), intent(in) :: itask
  integer(4)              :: istat

  select case(itask)

  case(1_ip)
!     ! no onger beeing used - can erase this lines
!     !
!     ! WALLD_TUR: Wall distance
!     !
!     allocate(walld_tur(npoin),stat=istat)
!     call memchk(zero,istat,mem_modul(1:2,modul),'WALLD_TUR','tur_memarr',walld_tur)

  case(2_ip)
     !
     ! LWNEI_TUR: Nearest neighbor
     !
     allocate(lwnei_tur(npoin),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'LWNEI_TUR','tur_addarr',lwnei_tur)

  case(3_ip)
     !
     ! Empty
     !

  case(5_ip)
     !
     ! USTAR_TUR: Velocity gradient
     !
     if( kfl_ustar_tur/=0 .and. INOTMASTER ) then
        allocate(ustar_tur(nbopo),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'USTAR_TUR','tur_memarr',ustar_tur)
     end if

  case(6_ip)
     !
     ! LWNEI_TUR: deallocate
     !
     call memchk(two,istat,mem_modul(1:2,modul),'lwnei_tur','tur_memarr',lwnei_tur)
     deallocate(lwnei_tur,stat=istat)
     if(istat/=0) call memerr(two,'LWNEI_TUR','tur_memarr',0_ip)    
     
  case(7_ip)
!     this is no longer beeing used - can erase
!     !
!     ! WALLD_TUR: deallocate
!     !
!     call memchk(two,istat,mem_modul(1:2,modul),'lwnei_tur','tur_memarr',walld_tur)
!     deallocate(walld_tur,stat=istat)
!     if(istat/=0) call memerr(two,'WALLD_TUR','tur_memarr',0_ip)    

  case(8_ip)
     !
     ! KFL_GRK12_TUR, KFL_GRONO_TUR, KFL_GREPS_TUR
     !
     if( INOTMASTER ) then
        if(kfl_grk12_tur==1) then
           allocate(grk12_tur(nbopo),stat=istat)
           call memchk(zero,istat,mem_modul(1:2,modul),'GRK12_TUR','tur_memarr',grk12_tur)
        end if
        if(kfl_grono_tur==1) then
           allocate(grono_tur(nbopo),stat=istat)
           call memchk(zero,istat,mem_modul(1:2,modul),'GRONO_TUR','tur_memarr',grono_tur)
        end if
        if(kfl_greps_tur==1) then
           allocate(greps_tur(nbopo),stat=istat)
           call memchk(zero,istat,mem_modul(1:2,modul),'GREPS_TUR','tur_memarr',greps_tur)
        end if
     end if

  case(9_ip)
     !
     ! GRVE2_TUR: Velocity 2nd order gradients
     !
     if( kfl_grve2_tur==1 .and. INOTMASTER ) then
        allocate(grve2_tur(npoin),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'GRVE2_TUR','tur_memarr',grve2_tur)     
     end if

  case(10_ip)
     !
     ! GRSQK_TUR: grad(sqrt(k))
     !
     if( kfl_grsqk_tur==1 .and. INOTMASTER ) then
        allocate(grsqk_tur(npoin),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'GRSQK_TUR','tur_memarr',grsqk_tur)     
     end if

  case(11_ip)
     !
     ! GRPHI_TUR: grad(phi)
     !
     if( kfl_grphi_tur==1 .and. INOTMASTER ) then 
        allocate(grphi_tur(ndime,npoin),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'GRPHI_TUR','tur_memarr',grphi_tur)
     end if
     
  case(12_ip)
     !
     ! VORTI_TUR: W
     !
     if( kfl_vorti_tur == 1 .and. INOTMASTER ) then 
        allocate(vorti_tur(npoin),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'VORTI_TUR','tur_memarr',vorti_tur)
     end if
     
  end select

end subroutine tur_memarr
