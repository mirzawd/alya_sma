!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup NastinInput
!> @{
!> @file    nsi_membcs.f90
!> @author  Guillaume Houzeaux
!> @brief   Allocate memory for boudnary conditions
!> @details Allocate memory for boudnary conditions
!> @} 
!-----------------------------------------------------------------------
subroutine nsi_membcs(itask)
  use def_kintyp, only : ip
  use def_master, only : mem_modul
  use def_master, only : modul
  use def_master, only : INOTMASTER
  use def_domain, only : ndime
  use def_domain, only : npoin
  use def_domain, only : nboun
  use def_nastin, only : kfl_conbc_nsi
  use def_nastin, only : kfl_divcorrec_nsi
  use def_nastin, only : kfl_fixno_nsi
  use def_nastin, only : kfl_fixbo_nsi
  use def_nastin, only : kfl_fixrs_nsi
  use def_nastin, only : kfl_fixpr_nsi
  use def_nastin, only : kfl_fixpp_nsi
  use def_nastin, only : kfl_funno_nsi
  use def_nastin, only : kfl_funbo_nsi
  use def_nastin, only : kfl_funtn_nsi
  use def_nastin, only : kfl_funtb_nsi
  use def_nastin, only : kfl_wlawf_nsi
  use def_nastin, only : bvess_nsi
  use def_nastin, only : bvnat_nsi
  use def_nastin, only : bpess_nsi
  use def_nastin, only : kfl_fixno_div_nsi
  use mod_memory, only : memory_alloca
  use mod_memory, only : memory_alloca_min
  implicit none
  integer(ip), intent(in) :: itask

  select case ( itask )

  case( 1_ip )

     if( INOTMASTER ) then
        !
        ! Mandatory arrays
        !
        call memory_alloca( mem_modul(1:2,modul),'KFL_FIXNO_NSI','nsi_membcs',kfl_fixno_nsi,ndime,npoin      )   ! Velocity components
        call memory_alloca( mem_modul(1:2,modul),'KFL_FIXBO_NSI','nsi_membcs',kfl_fixbo_nsi,nboun            )   ! Momentum boundary fixity
        call memory_alloca( mem_modul(1:2,modul),'KFL_FIXRS_NSI','nsi_membcs',kfl_fixrs_nsi,npoin            )   ! Local axes
        call memory_alloca( mem_modul(1:2,modul),'KFL_FIXPR_NSI','nsi_membcs',kfl_fixpr_nsi,1_ip,npoin       )   ! Pressure Schur complement fixity
        call memory_alloca( mem_modul(1:2,modul),'KFL_FIXPP_NSI','nsi_membcs',kfl_fixpp_nsi,1_ip,npoin       )   ! Pressure fixity
        call memory_alloca( mem_modul(1:2,modul),'KFL_WLAWF_NSI','nsi_membcs',kfl_wlawf_nsi,npoin            )   ! Wall law flag

        if( kfl_conbc_nsi == 1 ) then 
           call memory_alloca( mem_modul(1:2,modul),'BVESS_NSI'     ,'nsi_membcs',bvess_nsi,ndime,npoin,1_ip )   ! Velocity values
           call memory_alloca( mem_modul(1:2,modul),'BVNAT_NSI'     ,'nsi_membcs',bvnat_nsi,5_ip ,nboun,1_ip )   ! Natural condition values
           call memory_alloca( mem_modul(1:2,modul),'BPESS_NSI'     ,'nsi_membcs',bpess_nsi,1_ip ,npoin,1_ip )   ! Pressure values
        else
           call memory_alloca( mem_modul(1:2,modul),'BVESS_NSI'     ,'nsi_membcs',bvess_nsi,ndime,npoin,2_ip )
           call memory_alloca( mem_modul(1:2,modul),'BVNAT_NSI'     ,'nsi_membcs',bvnat_nsi,5_ip, nboun,2_ip )
           call memory_alloca( mem_modul(1:2,modul),'BPESS_NSI'     ,'nsi_membcs',bpess_nsi,1_ip, npoin,2_ip )  
           call memory_alloca( mem_modul(1:2,modul),'KFL_FUNNO_NSI' ,'nsi_membcs',kfl_funno_nsi , npoin      )
           call memory_alloca( mem_modul(1:2,modul),'KFL_FUNBO_NSI' ,'nsi_membcs',kfl_funbo_nsi , nboun      )
           call memory_alloca( mem_modul(1:2,modul),'KFL_FUNTN_NSI' ,'nsi_membcs',kfl_funtn_nsi , npoin      )
           call memory_alloca( mem_modul(1:2,modul),'KFL_FUNTB_NSI' ,'nsi_membcs',kfl_funtb_nsi , nboun      )
        end if
        !
        ! Divergence free correction
        !
        if( kfl_divcorrec_nsi /= 0 ) then
           call memory_alloca(mem_modul(1:2,modul),'KFL_FIXNO_DIV_NSI','nsi_memall',kfl_fixno_div_nsi,1_ip,npoin)
        end if

     end if

     call memory_alloca_min( mem_modul(1:2,modul),'KFL_FIXNO_NSI','nsi_membcs',kfl_fixno_nsi)  
     call memory_alloca_min( mem_modul(1:2,modul),'KFL_FIXBO_NSI','nsi_membcs',kfl_fixbo_nsi)  
     call memory_alloca_min( mem_modul(1:2,modul),'KFL_FIXRS_NSI','nsi_membcs',kfl_fixrs_nsi)  
     call memory_alloca_min( mem_modul(1:2,modul),'KFL_FIXPR_NSI','nsi_membcs',kfl_fixpr_nsi)  
     call memory_alloca_min( mem_modul(1:2,modul),'KFL_FIXPP_NSI','nsi_membcs',kfl_fixpp_nsi)  
     call memory_alloca_min( mem_modul(1:2,modul),'KFL_WLAWF_NSI','nsi_membcs',kfl_wlawf_nsi)
     
     call memory_alloca_min( mem_modul(1:2,modul),'BVESS_NSI'    ,'nsi_membcs',bvess_nsi)  
     call memory_alloca_min( mem_modul(1:2,modul),'BVNAT_NSI'    ,'nsi_membcs',bvnat_nsi)  
     call memory_alloca_min( mem_modul(1:2,modul),'BPESS_NSI'    ,'nsi_membcs',bpess_nsi)  
     
  end select

end subroutine nsi_membcs

