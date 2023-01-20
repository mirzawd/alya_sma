!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Neutro
!> @{
!> @file    neu_membcs.f90
!> @author  Guillaume Houzeaux
!> @brief   Allocate memory for boundary conditions
!> @details Allocate memory for boundary conditions
!> @} 
!-----------------------------------------------------------------------
subroutine neu_membcs(itask)
  use def_kintyp, only : ip
  use def_master, only : mem_modul
  use def_master, only : modul
  use def_master, only : INOTMASTER
!   use def_domain, only : ndime
  use def_domain, only : npoin
  use def_domain, only : nboun
  use def_neutro, only : kfl_fixno_neu
  use def_neutro, only : kfl_fixbo_neu, kfl_funbo_neu, kfl_funtb_neu
  use def_neutro, only : bvess_neu
  use def_neutro, only : bvnat_neu
  use def_neutro, only : num_energies_neu
  use def_neutro, only : num_directions_neu
  use mod_memory, only : memory_alloca
  use mod_memory, only : memory_alloca_min
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: ienergy,idirection

  select case ( itask )

  case( 1_ip )

     if( INOTMASTER ) then
        !
        ! Codes
        !
        call memory_alloca( mem_modul(1:2,modul),'KFL_FIXNO_NEU','neu_membcs',kfl_fixno_neu,num_energies_neu,num_directions_neu)   
        call memory_alloca( mem_modul(1:2,modul),'KFL_FIXBO_NEU','neu_membcs',kfl_fixbo_neu,num_energies_neu,num_directions_neu) 
        call memory_alloca( mem_modul(1:2,modul),'KFL_FUNBO_NEU','neu_membcs',kfl_funbo_neu,num_energies_neu,num_directions_neu) 
        call memory_alloca( mem_modul(1:2,modul),'KFL_FUNTB_NEU','neu_membcs',kfl_funtb_neu,num_energies_neu,num_directions_neu) 
        do idirection = 1,num_directions_neu
           do ienergy = 1,num_energies_neu
              call memory_alloca( mem_modul(1:2,modul),'KFL_FIXNO_NEU % l','neu_membcs',&
                                  kfl_fixno_neu(ienergy,idirection) % l,1_ip,npoin )  
              call memory_alloca( mem_modul(1:2,modul),'KFL_FIXBO_NEU % l','neu_membcs',kfl_fixbo_neu(ienergy,idirection) % l,nboun )
              call memory_alloca( mem_modul(1:2,modul),'KFL_FUNBO_NEU % l','neu_membcs',kfl_funbo_neu(ienergy,idirection) % l,nboun )
              call memory_alloca( mem_modul(1:2,modul),'KFL_FUNTB_NEU % l','neu_membcs',kfl_funtb_neu(ienergy,idirection) % l,nboun )
           end do
        end do
        !
        ! Values
        !
        call memory_alloca( mem_modul(1:2,modul),'BVESS_NEU','neu_membcs',bvess_neu,num_energies_neu,num_directions_neu)  
        call memory_alloca( mem_modul(1:2,modul),'BVNAT_NEU','neu_membcs',bvnat_neu,num_energies_neu,num_directions_neu)
        do idirection = 1,num_directions_neu
           do ienergy = 1,num_energies_neu
              call memory_alloca(mem_modul(1:2,modul),'BVESS_NEU % a','neu_membcs',bvess_neu(ienergy,idirection) % a,1_ip,npoin,1_ip)
              call memory_alloca(mem_modul(1:2,modul),'BVNAT_NEU % a','neu_membcs',bvnat_neu(ienergy,idirection) % a,1_ip,nboun,1_ip)
           end do
        end do
        
     end if
     
  end select

end subroutine neu_membcs

