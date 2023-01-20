!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup Neutro 
!> @{
!> @file    neu_inivar.f90
!> @author  Guillaume Houzeaux
!> @brief   Initialize some variables
!> @details Initialize some variables\n
!>          ITASK=0 ... When starting the run (from Turnon)\n
!>          ITASK=1 ... After reading data\n
!> @} 
!------------------------------------------------------------------------
subroutine neu_inivar(itask)
  use def_parame
  use def_master
  use def_neutro
  use def_domain
  use def_solver
  use mod_ADR,   only : ADR_initialize_type
  implicit none
  integer(ip), intent(in) :: itask
  external soldef

  select case(itask)

  case(0_ip)   

      ! EG (17/08/2022): HAY QUE HACER ESTO O NO ES NECESARIO?????
   ! call arrays_register( 11_ip,(/'HEASK','SCALA','NPOIN','PRIMA'/),heat_sink    ,ENTITY_POSITION=1_ip)

     !
     ! Postprocess Variable 
     ! name:   '.....' 
     ! type:   'SCALA'/'VECTO' 
     ! entity: 'NPOIN'/'NELEM'
     !
     postp(1) % wopos (1:3, 1) = (/ 'CURRE' , 'VECTO' , 'NPOIN' /)
     postp(1) % wopos (1:3, 2) = (/ 'FLUX ' , 'SCALA' , 'NPOIN' /)
     postp(1) % wopos (1:3, 3) = (/ 'RADIA' , 'SCALA' , 'NPOIN' /)
     postp(1) % wopos (1:3, 4) = (/ 'DIREC' , 'SCALA' , 'NPOIN' /)
     postp(1) % wopos (1:3, 5) = (/ 'MATER' , 'SCALA' , 'NELEM' /)
     postp(1) % wopos (1:3, 6) = (/ 'HEAT ' , 'SCALA' , 'NPOIN' /)
     postp(1) % wopos (1:3, 7) = (/ 'SOURC' , 'SCALA' , 'NPOIN' /)
     postp(1) % wopos (1:3, 8) = (/ 'VACUU' , 'SCALA' , 'NPOIN' /)
     postp(1) % wopos (1:3, 9) = (/ 'REFLE' , 'SCALA' , 'NPOIN' /)
     postp(1) % wopos (1:3, 10) = (/ 'HEATG' , 'SCALA' , 'NPOIN' /)
     !
     ! Element set variables 
     !
     postp(1) % woese ( 1)     = 'VARI1'   ! Variable 1
     postp(1) % woese ( 2)     = 'VARI2'   ! Variable 2
     postp(1) % woese ( 3)     = 'VARI3'   ! Variable 3
     postp(1) % woese ( 4)     = 'VARI4'   ! Variable 4
     postp(1) % woese ( 5)     = 'VARI5'   ! Variable 5
     postp(1) % woese ( 6)     = 'VARI6'   ! Variable 6
     postp(1) % woese ( 7)     = 'VARI7'   ! Variable 7
     postp(1) % woese ( 8)     = 'VARI8'   ! Variable 8
     postp(1) % woese ( 9)     = 'VARI9'   ! Variable 9
     postp(1) % woese ( 10)     = 'VARI10'   ! Variable 10
     !
     ! Boundary set variables
     !
     postp(1) % wobse ( 1)     = 'VARI1'   ! Variable 1
     postp(1) % wobse ( 2)     = 'VARI2'   ! Variable 2
     postp(1) % wobse ( 3)     = 'VARI3'   ! Variable 3
     postp(1) % wobse ( 4)     = 'VARI4'   ! Variable 4
     postp(1) % wobse ( 5)     = 'VARI5'   ! Variable 5
     postp(1) % wobse ( 6)     = 'VARI6'   ! Variable 6
     postp(1) % wobse ( 7)     = 'VARI7'   ! Variable 7
     postp(1) % wobse ( 8)     = 'VARI8'   ! Variable 8
     postp(1) % wobse ( 9)     = 'VARI9'   ! Variable 9
     postp(1) % wobse ( 10)     = 'VARI10'   ! Variable 10
     !
     ! Node set variables
     !
     postp(1) % wonse ( 1)     = 'VARI1'   ! Variable 1
     postp(1) % wonse ( 2)     = 'VARI2'   ! Variable 2
     postp(1) % wonse ( 3)     = 'VARI3'   ! Variable 3
     postp(1) % wonse ( 4)     = 'VARI4'   ! Variable 4
     postp(1) % wonse ( 5)     = 'VARI5'   ! Variable 5
     postp(1) % wonse ( 6)     = 'VARI6'   ! Variable 6
     postp(1) % wonse ( 7)     = 'VARI7'   ! Variable 7
     postp(1) % wonse ( 8)     = 'VARI8'   ! Variable 8
     postp(1) % wonse ( 9)     = 'VARI9'   ! Variable 9
     postp(1) % wonse ( 10)     = 'VARI10'   ! Variable 10
     !
     ! Witness variables 
     !
     postp(1) % wowit ( 1)     = 'VARI1'   ! Variable 1
     postp(1) % wowit ( 2)     = 'VARI2'   ! Variable 2
     postp(1) % wowit ( 3)     = 'VARI3'   ! Variable 3
     postp(1) % wowit ( 4)     = 'VARI4'   ! Variable 4
     postp(1) % wowit ( 5)     = 'VARI5'   ! Variable 5
     postp(1) % wowit ( 6)     = 'VARI6'   ! Variable 6
     postp(1) % wowit ( 7)     = 'VARI7'   ! Variable 7
     postp(1) % wowit ( 8)     = 'VARI8'   ! Variable 8
     postp(1) % wowit ( 9)     = 'VARI9'   ! Variable 9
     postp(1) % wowit ( 10)     = 'VARI10'   ! Variable 10
     !
     ! Solvers
     !     
     call soldef(-1_ip)                           ! Allocate memory for NUM_SOLVERS solvers
     solve(1:1) % kfl_solve = 1                   ! Output flag
     solve(1) % ndofn       = 1                   ! dof
     solve(1) % wprob       = 'NEUTRON_RADIATION' ! Solver name
     !
     ! Nullify pointers
     !
     nullify(kfl_fixno_neu)
     nullify(kfl_fixbo_neu)
     nullify(kfl_funbo_neu)
     nullify(kfl_funtb_neu)
     nullify(bvess_neu)
     nullify(bvnat_neu)
     nullify(tncod_neu)    
     nullify(tbcod_neu)
     nullify(aniso_neu)

     nullify(efectivos_neu_in)
     nullify(efectivos_neu_out)
     nullify(source_bound_neu)
     nullify(absor_neu)
     nullify(scatt_neu)
     nullify(fiss_neu)
     nullify(funsource)
     nullify(At_weight)
     nullify(Densidad_)
     nullify(Isotope)
     nullify(aniso_neu)
     nullify(absor_neu_cte)
     nullify(scatt_neu_cte)
     nullify(grupo_energias)
     nullify(phi_neu)
     nullify(tita_neu) 

     nullify(num_isotopes_neu)
     nullify(At_weight_isotope)
     nullify(mass_percentage_isotope_neu)
     nullify(atom_percentage_isotope_neu)
     nullify(kerma_neu)
    !  nullify(kerma_isotope_neu)
     nullify(kerma_poin_neu)

     !
     ! ADR type
     !
     call ADR_initialize_type(ADR_neu)

  case(1_ip)   
     !
     ! Derived parameters
     !
     if( ADR_neu % kfl_time_integration /= 0 ) then
        ncomp_neu = 3
     else
        ncomp_neu = 2
     end if
     nprev_neu = min(3_ip,ncomp_neu)  ! Last time step or global iteration

     kfl_stead_neu=0

  case(2_ip)   
     !
     ! Dimensions
     !
     nunkn_neu = num_energies_neu * num_directions_neu
     
  end select

end subroutine neu_inivar
