!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine neu_parall(order)
  !-----------------------------------------------------------------------
  !****f* Parall/neu_parall
  ! NAME
  !    neu_parall
  ! DESCRIPTION
  !    This routine exchange data 
  ! USES
  ! USED BY
  !    Reapro
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_neutro
  use def_domain
  use def_inpout
  use def_solver
  use mod_memory
  use mod_opebcs
  use mod_ADR,                only : ADR_parallel_data
  use mod_exchange,           only : exchange_init
  use mod_exchange,           only : exchange_add
  use mod_exchange,           only : exchange_end
  use mod_solver,             only : solver_parall
  use mod_output_postprocess, only : output_postprocess_parall

  implicit none
  integer(ip), intent(in) :: order
  external :: iexcha, rexcha, par_broadc, soldef

  if( ISEQUEN ) return

  select case ( order )

  case ( 1_ip )

     !-------------------------------------------------------------------
     !
     ! Broadcast physical data read in *.neu.dat file
     !
     !-------------------------------------------------------------------

     call exchange_init()
     !
     ! Physical problem
     !           
     call exchange_add(num_energies_neu) 
     call exchange_add(num_directions_neu)
     call exchange_add(num_materials_neu)
     call exchange_add(num_sources_neu)
     call exchange_add(num_legendre_neu)
     call exchange_add(num_legendre_lee)

     call exchange_add(kfl_icosa_neu)
     call exchange_add(kfl_snord_neu)

     call exchange_add(albedo_neu)
     !
     ! ADR data
     !
     call ADR_parallel_data(ADR_neu)

     call exchange_end()

  case ( 2_ip ) 

     !-------------------------------------------------------------------
     !
     ! Other data
     !
     !-------------------------------------------------------------------

     call boundary_conditions_exchange(tncod_neu)
     call boundary_conditions_exchange(tbcod_neu)

     call exchange_init()
     !
     ! Numerical treatment
     !
     call exchange_add(miinn_neu) 
     call exchange_add(kfl_smobo_neu)
     call exchange_add(cotol_neu) 
     call exchange_add(relax_neu) 
     call exchange_add(nitsche_neu)
     !
     ! Solver and postprocess
     !
     call solver_parall()
     call output_postprocess_parall()
     call ADR_parallel_data(ADR_neu)
     
     call exchange_end()
 

  case ( 3_ip )

     ! debo enviar
     call exchange_init()
     !
     ! Physical problem
     !

     call exchange_add(aniso_neu)
     call exchange_add(At_weight)
     call exchange_add(Densidad_)
     call exchange_add(Isotope)
     call exchange_add(scatt_neu_cte)
     call exchange_add(absor_neu_cte)
     call exchange_add(funsource)

     call exchange_add(ener_weigd_neu)
     call exchange_add(efectivos_neu_in)
     call exchange_add(efectivos_neu_out)

     call exchange_add(grupo_energias)
     call exchange_add(absor_neu)
     call exchange_add(source_bound_neu)
     call exchange_add(fiss_neu)

     call exchange_add(scatt_neu)
     call exchange_add(kerma_neu)
     call exchange_add(n_efectivos_neu_sparse)
     call exchange_add(n_efectivos_neu_sparse_acum)
     call exchange_add(efectivos_neu_sparse)
     call exchange_add(resid_energy_group_neu)
     
     call exchange_end()

  case ( 4_ip )

     !-------------------------------------------------------------------
     !
     ! Broadcast variables needed for allocating
     !
     !-------------------------------------------------------------------

     call exchange_init()
     call exchange_add(max_val_acum_neu) 
     call exchange_end()

  end select

end subroutine neu_parall

