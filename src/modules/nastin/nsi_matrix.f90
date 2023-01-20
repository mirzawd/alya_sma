!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup nastin
!> @{
!> @file    nsi_matrix.f90
!> @author  houzeaux
!> @date    2020-03-21
!> @brief   Assemble system
!> @details Assemble NS system according to the method used:
!>          Schur:        Auu, Aup, Apu, App, bu, bp, Q
!>          FS classical: bu, Aup, Apu, bu, bp, Q    
!>          FS fast:      bu, bp. Grad and Div computed apart   
!> @} 
!-----------------------------------------------------------------------

subroutine nsi_matrix()

  use def_solver
  use def_master
  use def_nastin
  use def_domain
  use mod_elmgeo
  use mod_ker_timeline,                only : ker_timeline
  use def_kermod,                      only : kfl_waexl_ker,velel_ker
  use def_kermod,                      only : kfl_noslw_ker
  use mod_interpolation,               only : COU_GET_INTERPOLATE_POINTS_VALUES
  use mod_solver,                      only : solver_preprocess
  use mod_solver,                      only : solver_initialize_matrix_and_rhs
  use mod_solver,                      only : solver_lumped_mass_system
  use mod_solver,                      only : solver_exchange
  use mod_communications,              only : PAR_MAX
  use mod_communications,              only : PAR_BARRIER
  use mod_nsi_dirichlet_global_system, only : nsi_dirichlet_global_system
  use mod_nsi_assembly_global_system,  only : nsi_assembly_algebraic_split_oss
  use mod_nsi_assembly_global_system,  only : nsi_assembly_scaling_pressure_equation
  use mod_timings,                     only : timings_assembly,timings_unity_test
  use mod_alya2dlb,                    only : alya2dlb_DLB_Enable
  use mod_alya2dlb,                    only : alya2dlb_DLB_Disable
  use mod_memory,                      only : memory_initia
  use mod_nsi_algebraic_forces,        only : nsi_algebraic_forces
  use mod_ker_updpro,                  only : ker_updpro
  use mod_nsi_low_Mach,                only : nsi_low_Mach_drho_dt
  use mod_nsi_spare_mesh,              only : nsi_spare_mesh_dirichlet
  use mod_wall_exchange,               only : ker_waexlo_getval
  use mod_wall_exchange,               only : ker_waexlo_getder
  use mod_nsi_elmope_all

  implicit none

  integer(ip) :: ierr, ipoin
  real(rp)    :: time1,time2,time3,time4
  real(rp)    :: timea,timeb

  call cputim(timea)

  call ker_timeline('INI_ASSEMBLY')

  !-----------------------------------------------------------------
  !
  ! Initializations
  !
  !-----------------------------------------------------------------

  time1     =  0.0_rp
  time2     =  0.0_rp
  time3     =  0.0_rp
  tamin_nsi =  huge(1.0_rp)
  tamax_nsi = -huge(1.0_rp)
  call memory_initia(bupor_nsi)
  porfo_nsi = 0.0_rp 

  !-----------------------------------------------------------------
  !
  ! Properties
  !
  !-----------------------------------------------------------------

  if ( kfl_noslw_ker /= 0_ip ) &
       call ker_updpro(ITASK_BEGINN)

  !-----------------------------------------------------------------
  !
  ! Matrix, preconditioner, RHS  and projections initializations
  !
  !-----------------------------------------------------------------

  if( INOTEMPTY ) then
     !
     ! Matrices, RHS and preconditioner
     !
     if( NSI_SEMI_IMPLICIT ) then
        call solver_initialize_matrix_and_rhs(solve(NSI_SOLVER_VISCOUS_TERM),visco_nsi)
     end if
     if( NSI_MONOLITHIC .or. kfl_grad_div_nsi /= 0 ) then
        if( kfl_tisch_nsi==3 .or. kfl_tisch_nsi==4) then
           call memory_initia(rhsid,LBOUN1=1_ip,UBOUN1=(ndime+1)*npoin)   ! Should be ndime+1  or ndime -- rethink 
        else
           call solver_initialize_matrix_and_rhs(solve,amatr,rhsid)
        end if
     else
        call solver_initialize_matrix_and_rhs(solve,amatr,rhsid,lapla_nsi)
     end if
     if( solve(ivari_nsi) % kfl_preco >= 3 ) pmatr(1:solve(ivari_nsi) % nzpre) = 0.0_rp

     if(  kfl_corre_nsi == 3 ) cmama_nsi = 0.0_rp  ! Consistent mass matrix
     if(  kfl_predi_nsi == 7 .or. &
          kfl_predi_nsi == 8 .or. &
          kfl_predi_nsi == 9 ) then
        call memory_initia(dt_rho_nsi)             ! rho / dt 
        call memory_initia(tau_nsi)                ! rho / tau
     end if
  end if
  if( associated(mass_rho_nsi) ) then
     do ipoin = 1,npoin
        mass_rho_nsi(ipoin,1) = 0.0_rp
     end do
  end if

  if( associated(rhsid_gravb) ) then
     do ipoin = 1,npoin * ndime
        rhsid_gravb(ipoin) = 0.0_rp
     end do
  end if
  !
  ! Accumulate CPU time not taken into account in elmope
  !
  call cputim(timeb)
  cputi_assembly_nsi(9) = cputi_assembly_nsi(9) + timeb - timea

  !-----------------------------------------------------------------
  !
  ! Update boundary conditions
  !
  !-----------------------------------------------------------------

  call nsi_updbcs(ITASK_BEGINN) 

  !-----------------------------------------------------------------
  !
  ! Element assembly: Matrix, preconditioner, RHS and projections
  !
  !-----------------------------------------------------------------

  if( INOTMASTER ) ierr = alya2dlb_DLB_Enable()

  call cputim(time1)

  if( kfl_vector_nsi == 0 ) then     
     call nsi_elmope_omp(1_ip) ! Classical version    
  else     
     call nsi_elmope_all(1_ip) ! Vectorized version     
  end if

  call cputim(time2)

  call cputim(timea)
#ifdef ALYA_DLB
  call PAR_BARRIER()
  if( INOTMASTER ) ierr = alya2dlb_DLB_Disable()
#endif
  
  !-----------------------------------------------------------------
  !
  ! Wall with exchange location
  !
  !-----------------------------------------------------------------

  if( kfl_waexl_ker == 1_ip ) then
     call times(2) % ini() 
     call ker_waexlo_getval(veloc,velel_ker)
     call times(2) % add()
  end if

  !-----------------------------------------------------------------
  !
  ! Boundary assembly
  !
  !-----------------------------------------------------------------

  call cputim(time3)
  call nsi_bouope_all(0_ip)
  call cputim(time4)

  !-----------------------------------------------------------------
  !
  ! Immersed Dirichlet boundary conditions 
  !
  !-----------------------------------------------------------------

  call nsi_spare_mesh_dirichlet()
  
  !-----------------------------------------------------------------
  !
  ! Internal force
  !
  !-----------------------------------------------------------------

  if( .not. NSI_MONOLITHIC ) then
     if(.not. NSI_SEMI_IMPLICIT) then
        call nsi_algebraic_forces(amatr(poauu_nsi:),amatr(poaup_nsi:),rhsid)
     end if
  end if

  !-----------------------------------------------------------------
  !
  ! Coupling with Immbou: compute force 
  !
  !-----------------------------------------------------------------

  call nsi_coupli(ITASK_MATRIX)

  !-----------------------------------------------------------------
  !
  ! Solver projection of dt/rho and tau
  !
  !-----------------------------------------------------------------

  if( kfl_predi_nsi == 7 .or. kfl_predi_nsi == 8 .or. kfl_predi_nsi == 9 ) then
     call solver_lumped_mass_system(1_ip,dt_rho_nsi)
     call solver_lumped_mass_system(1_ip,tau_nsi)
  end if
  ! if( associated(mass_rho_nsi) ) call solver_exchange(1_ip,mass_rho_nsi)

  !
  ! Compute drho/dt if necessary
  !
  if( associated(mass_rho_nsi) ) then
     call solver_exchange(1_ip,mass_rho_nsi) 
     if (associated(drhodt_nsi)) then
        if (ittim /= 1 .or. kfl_rstar /= 0 ) then
           call nsi_low_Mach_drho_dt(0_ip)
        else
           do ipoin = 1,npoin
              drhodt_nsi(ipoin) = 0.0_rp
           end do
        endif
     endif
  end if


  !-----------------------------------------------------------------
  !
  ! Solver pre-process and boundary conditions
  !
  !-----------------------------------------------------------------
  
  if( .not. NSI_FRACTIONAL_STEP ) then
     call solver_preprocess(momod(modul) % solve,amatr,rhsid,unkno,lapla_nsi)
  end if
  call nsi_dirichlet_global_system()

  !-----------------------------------------------------------------
  !
  ! Algebraic split OSS method for pressure
  ! ---------------------------------------
  !
  ! Pressure and velocity are not stabilized at all previously and
  ! pressure stabilization is added at the algebraic level:
  !
  ! App = Tau L
  ! bp  = bp - Apu Tau M^{-1} Aup p
  !
  ! so that the pressure stabilization matrix is
  !
  ! Tau L + Apu Tau M^{-1} Aup
  !
  !-----------------------------------------------------------------
  
  if( kfl_stabi_nsi == NSI_ALGEBRAIC_SPLIT_OSS .and. INOTEMPTY ) &
       call nsi_assembly_algebraic_split_oss(&
       amatr(poaup_nsi:),amatr(poapu_nsi:),&
       amatr(poapp_nsi:),lapla_nsi,&
       rhsid(ndbgs_nsi+1:),unkno(ndbgs_nsi+1:))

  !-----------------------------------------------------------------
  !
  ! Penalize pressure equation
  !
  !-----------------------------------------------------------------

  if( INOTMASTER .and. kfl_prepe_nsi == 1 ) &
       call nsi_penpre(&
       amatr(poapp_nsi),lapla_nsi,rhsid(ndbgs_nsi+1))

  !-----------------------------------------------------------------
  !
  ! Scaling of pressure equation
  ! ----------------------------
  !
  ! The objective is to force pressure Schur complement
  ! not to change => Q \approx Scal L, where L is the Laplacian
  ! and Scal is computed using the projections of tau and dt / rho.
  !
  !-----------------------------------------------------------------

  if( NSI_SCHUR_COMPLEMENT .and. INOTMASTER ) &
       call nsi_assembly_scaling_pressure_equation(&
       amatr(poapu_nsi:),amatr(poapp_nsi:),&
       rhsid(ndbgs_nsi+1:))
  !
  ! CPU times
  !
  cpu_ass_sol_nsi(1) = time2 - time1 
  cpu_ass_sol_nsi(4) = time4 - time3

  call timings_assembly(cpu_ass_sol_nsi(1),cpu_ass_sol_nsi(4),TYPE_OF_ASSEMBLY='ELEMENT, BOUNDARY')

  call ker_timeline('END_ASSEMBLY')
  !
  ! Accumulate CPU time not taken into account in elmope
  !
  call cputim(timeb)
  cputi_assembly_nsi(9) = cputi_assembly_nsi(9) + timeb - timea

end subroutine nsi_matrix
