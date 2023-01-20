!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup NeutroTurnon
!> @{
!> @file    neu_memall.f90
!> @date    31/03/2016
!> @author  Guillaume Houzeaux
!> @brief   Allocate memory 
!> @details Allocate memory 
!> @} 
!-----------------------------------------------------------------------
subroutine neu_memall(itask)
! subroutine neu_memall()
  use def_kintyp, only : ip
  use def_master, only : neutr,heat_sink
  use def_master, only : solve
  use def_master, only : INOTMASTER, ISLAVE
  use def_master, only : mem_modul,modul
  use def_domain, only : npoin!,ndime
  use def_solver, only : solve_sol
  use def_neutro, only : num_energies_neu
  use def_neutro, only : num_directions_neu
  use def_neutro, only : num_materials_neu
  use def_neutro, only : num_legendre_lee
  use def_neutro, only : num_sources_neu
  use def_neutro, only : ncomp_neu
  use def_neutro, only : direc_neu,phi_neu,tita_neu
  use def_neutro, only : weigd_neu,ener_weigd_neu
  use def_neutro, only : scattering_neu
  use def_neutro, only : absor_neu
  use def_neutro, only : GRUPO_ENERGIAS
  use def_neutro, only : SCATT_NEU
  use def_neutro, only : FISS_NEU
  use def_neutro, only : EFECTIVOS_NEU_OUT
  use def_neutro, only : EFECTIVOS_NEU_in
  use def_neutro, only : source_bound_neu
  use def_neutro, only : aniso_neu, At_weight,Isotope, funsource, absor_neu_cte, scatt_neu_cte,Densidad_
  use def_neutro, only : kerma_poin_neu, kerma_neu
  use def_neutro, only : efectivos_neu_sparse, n_efectivos_neu_sparse, n_efectivos_neu_sparse_acum, max_val_acum_neu
  use def_neutro, only : resid_energy_group_neu

  use mod_memory, only : memory_alloca
  use mod_memory, only : memory_alloca_min
  implicit none
  integer(ip), intent(in) :: itask
  external :: soldef

select case (itask)
  case(1_ip)
    !----------------------------------------------------------------------
    !
    ! Solver
    !
    !----------------------------------------------------------------------

    solve_sol => solve(1:1)
    call soldef(4_ip)

    !----------------------------------------------------------------------
    !
    ! Arrays
    !
    !----------------------------------------------------------------------
    !
    ! Directions
    !
    call memory_alloca(mem_modul(1:2,modul),'DIREC_NEU'     ,'neu_memall',direc_neu     ,3_ip,num_directions_neu)
    call memory_alloca(mem_modul(1:2,modul),'SCATTERING_NEU','neu_memall',scattering_neu,num_directions_neu,num_directions_neu)
    call memory_alloca(mem_modul(1:2,modul),'ENER_WEIGD_NEU'     ,'neu_memall',ener_weigd_neu  ,num_energies_neu)
    call memory_alloca(mem_modul(1:2,modul),'WEIGD_NEU'     ,'neu_memall',weigd_neu     ,num_directions_neu)

    call memory_alloca(mem_modul(1:2,modul),'PHI_NEU'     ,'neu_memall',phi_neu     ,num_directions_neu)
    call memory_alloca(mem_modul(1:2,modul),'TITA_NEU'     ,'neu_memall',tita_neu     ,num_directions_neu)

    ! call memory_alloca(mem_modul(1:2,modul),'KERMA_POIN_NEU','neu_memall',kerma_poin_neu,npoin,num_energies_neu)
    
    if( INOTMASTER ) then
      !
      ! RADIA
      !
      call memory_alloca(mem_modul(1:2,modul),'RADIA','neu_memall',neutr,num_energies_neu,num_directions_neu,npoin,ncomp_neu)
      call memory_alloca(mem_modul(1:2,modul),'KERMA_POIN_NEU','neu_memall',kerma_poin_neu,npoin,num_energies_neu)
      call memory_alloca(mem_modul(1:2,modul),'HEAT_SINK','neu_memall',heat_sink,npoin)

    if(ISLAVE) then
      call memory_alloca(mem_modul(1:2,modul),'EFECTIVOS_NEU_IN' ,'neu_memall',efectivos_neu_in  ,num_energies_neu)
      call memory_alloca(mem_modul(1:2,modul),'efectivos_neu_out','neu_memall',efectivos_neu_out  ,num_energies_neu)
      call memory_alloca(mem_modul(1:2,modul),'absor_neu'      ,'neu_memall',absor_neu  ,num_materials_neu,num_energies_neu)
      call memory_alloca(mem_modul(1:2,modul),'grupo_energias' ,'neu_memall',grupo_energias  ,num_materials_neu,num_energies_neu)
      call memory_alloca(mem_modul(1:2,modul),'scatt_neu','neu_memall',scatt_neu,num_materials_neu,&
                          num_energies_neu,num_energies_neu,num_legendre_lee+1)
      call memory_alloca(mem_modul(1:2,modul),'fiss_neu','neu_memall',fiss_neu,num_materials_neu,num_energies_neu,num_energies_neu)
      call memory_alloca(mem_modul(1:2,modul),'source_bound_neu','neu_memall',source_bound_neu, num_sources_neu, num_energies_neu)

      call memory_alloca(mem_modul(1:2,modul),'Isotope'      ,'neu_memall',Isotope  ,num_materials_neu)
      call memory_alloca(mem_modul(1:2,modul),'At_weight'      ,'neu_memall',At_weight  ,num_materials_neu)
      call memory_alloca(mem_modul(1:2,modul),'Densidad_'      ,'neu_memall',Densidad_ ,num_materials_neu)
      call memory_alloca(mem_modul(1:2,modul),'aniso_neu'      ,'neu_memall',aniso_neu  ,num_materials_neu)
      call memory_alloca(mem_modul(1:2,modul),'absor_neu_cte'      ,'neu_memall',absor_neu_cte  ,num_materials_neu)
      call memory_alloca(mem_modul(1:2,modul),'scatt_neu_cte'      ,'neu_memall',scatt_neu_cte  ,num_materials_neu)
      call memory_alloca(mem_modul(1:2,modul),'funsource'      ,'neu_memall',funsource  ,num_materials_neu)

      call memory_alloca(mem_modul(1:2,modul),'KERMA_NEU','neu_memall',kerma_neu,num_materials_neu,num_energies_neu)

      call memory_alloca(mem_modul(1:2,modul),'n_efectivos_neu_sparse','neu_memall',&
                          n_efectivos_neu_sparse,num_materials_neu,num_energies_neu)
      call memory_alloca(mem_modul(1:2,modul),'n_efectivos_neu_sparse_acum','neu_memall',&
                          n_efectivos_neu_sparse_acum,num_materials_neu,num_energies_neu+1)

      call memory_alloca(mem_modul(1:2,modul),'resid_energy_group_neu','neu_memall',resid_energy_group_neu,num_energies_neu)

    endif

    else
      !
      ! RADIA
      !
      !call memory_alloca_min(neutr)
  ! call memory_alloca(mem_modul(1:2,modul),'EFECTIVOS_NEU_IN' ,'neu_memall',efectivos_neu_in  ,num_energies_neu)
  ! call memory_alloca(mem_modul(1:2,modul),'efectivos_neu_out','neu_memall',efectivos_neu_out  ,num_energies_neu)
  ! call memory_alloca(mem_modul(1:2,modul),'absor_neu'      ,'neu_memall',absor_neu  ,num_materials_neu,num_energies_neu)
  ! call memory_alloca(mem_modul(1:2,modul),'grupo_energias' ,'neu_memall',grupo_energias  ,num_materials_neu,num_energies_neu)
  ! call memory_alloca(mem_modul(1:2,modul),'scatt_neu','neu_memall',scatt_neu,num_materials_neu,&
  !                    num_energies_neu,num_energies_neu,num_legendre_lee+1)
  ! call memory_alloca(mem_modul(1:2,modul),'fiss_neu','neu_memall',fiss_neu,num_materials_neu,num_energies_neu,num_energies_neu)
  ! call memory_alloca(mem_modul(1:2,modul),'source_bound_neu','neu_memall',source_bound_neu,num_energies_neu)

  !call memory_alloca(mem_modul(1:2,modul),'Isotope'      ,'neu_memall',Isotope  ,num_materials_neu)
  !call memory_alloca(mem_modul(1:2,modul),'At_weight'      ,'neu_memall',At_weight  ,num_materials_neu)
  !call memory_alloca(mem_modul(1:2,modul),'aniso_neu'      ,'neu_memall',aniso_neu  ,num_materials_neu)
  !call memory_alloca(mem_modul(1:2,modul),'absor_neu_cte'      ,'neu_memall',absor_neu_cte  ,num_materials_neu)
  !call memory_alloca(mem_modul(1:2,modul),'scatt_neu_cte'      ,'neu_memall',scatt_neu_cte  ,num_materials_neu)
  !call memory_alloca(mem_modul(1:2,modul),'funsource'      ,'neu_memall',funsource  ,num_materials_neu)


    end if
  case (2_ip)

    if (ISLAVE) then
      call memory_alloca(mem_modul(1:2,modul),'efectivos_neu_sparse','neu_memall',&
                        efectivos_neu_sparse,num_materials_neu,max_val_acum_neu)
    end if

  end select 

  
end subroutine neu_memall
