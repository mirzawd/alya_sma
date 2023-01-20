!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Exmedi
!> @{
!> @file    exm_solite.f90
!> @author  Mariano Vazquez
!> @date    26/12/2016
!> @brief   Update
!> @details Update
!> @} 
!-----------------------------------------------------------------------
subroutine exm_solite
  use def_parame
  use def_master
  use def_domain
  use def_exmedi
  use mod_solver,  only : solver_solve
  use mod_parall,  only : PAR_COMM_MY_CODE_ARRAY
  use mod_memory,  only : memory_size
  use mod_exm_elmoperations
  use mod_exm_ionicurrents

  implicit none
  integer(ip) :: ipoin,ii
  real(rp)    :: cpu_refe1,cpu_refe2,cpu_solve,amauxi

  integer(ip), parameter :: EXPLICIT_SCHEME=1, IMPLICIT_SCHEME=2, CRANK=2


  !
  ! Update inner iteration counter and write headings in the solver file.
  !
  itinn(modul) = itinn(modul) + 1

  !
  ! Initializations
  !

  cpu_solve = 0.0_rp
  call cputim(cpu_refe1)

  if( INOTEMPTY ) ticel_exm = 0.0_rp
  if( INOTEMPTY ) jicel_exm = 0.0_rp
  !
  ! Compute per-node ionic currents 
  !

  call exm_nodalionicurrents
  !
  ! Perform per-element operations
  !
  call exm_elmoperations

  !
  ! Add nodal contributions
  !

   if( INOTEMPTY ) then
      call rhsmod(ndime,eflux_exm)
      call rhsmod(ndime,bipol_exm)
      do ipoin=1, npoin
         eflux_exm(1:ndime,ipoin)  = eflux_exm( 1:ndime,ipoin) / vmass(ipoin)
         bipol_exm(1:ndime,ipoin)  = bipol_exm( 1:ndime,ipoin) / vmass(ipoin)
      end do
   end if

   if (kfl_timet_exm == IMPLICIT_SCHEME) then     
     !
     ! Retrieve amatr K and compute rhsid using it as -Ku 
     !

     do ii = 1,npoin
        rhsid(ii) = 0.0_rp
     end do
     do ii = 1,memory_size(amatr_auxi_exm)
        amatr(ii) = amatr_auxi_exm(ii)
     end do
     
     if (kfl_tisch_exm == CRANK) then
        do ipoin = 1,npoin
           unkno(ipoin) = 0.5_rp *(elmag(ipoin,ITER_K) + elmag(ipoin,TIME_N) + ticel_exm(ipoin))
        end do
     else
        do ipoin = 1,npoin
           unkno(ipoin) = elmag(ipoin,ITER_K) + ticel_exm(ipoin)
        end do
     end if

     !
     ! Unkno=Phi+ticel is here an auxiliar value used to compute amatr*(Phi+ticel) and substract it to rhsid 
     !
     call exm_amavec(amatr, unkno, rhsid, -1.0_rp)
     
     do ipoin=1,npoin
        !
        ! Time matrix added to K (multiplicity required when running in parallel)
        !
        amauxi = 1.0_rp / dtime  
        if (ipoin > npoi1) amauxi = amauxi / real(PAR_COMM_MY_CODE_ARRAY(1) % bound_multiplicity(ipoin-npoi1),rp)
        amatr(idima_exm(ipoin)) = amatr(idima_exm(ipoin)) + vmass(ipoin) * amauxi

     end do
   end if

  !
  ! Solve the incremental delta form (DPhi will be stored in unkno)
  !

  unkno     = 0.0_rp              ! clean unkno
  if (kfl_nodif_exm == 0) then
     !
     ! diffusion computed, solve system
     !
     if (kfl_timet_exm == EXPLICIT_SCHEME) then
        do ipoin=1,npoin  
           vdiag_exm(ipoin) = dtime
        end do
        solve(1)%xdiag =  1.0_rp
        call solver_solve(momod(modul) % solve, amatr, rhsid, unkno, vdiag_exm)     
     else if (kfl_timet_exm == IMPLICIT_SCHEME) then
        call solver_solve(momod(modul) % solve, amatr, rhsid, unkno)     
     end if


  end if

  !
  ! Keep the number of iterations of the solver for this subiteration
  !
  last_iters_exm = solve_sol(1)%iters

  !
  ! Re-assign unkno to the unknown. At this point unkno results from solving the system and it is unkno = DPhi
  !
  do ipoin=1,npoin                       

!  ESTO DE ABAJO ESTA MAL, LA ACTUALIZACION ES LA MISMA PARA CN O BDF!!
!     if (kfl_tisch_exm == CRANK) then
!        unkno(ipoin) = unkno(ipoin) + 0.5_rp *(elmag(ipoin,ITER_K) + elmag(ipoin,TIME_N)) + ticel_exm(ipoin)
!     else
!        unkno(ipoin) =  unkno(ipoin) + elmag(ipoin,ITER_K) + ticel_exm(ipoin)
!     end if

     unkno(ipoin) =  unkno(ipoin) + elmag(ipoin,ITER_K) + ticel_exm(ipoin)


     if (kfl_appva_exm == 1) then
        ! when the stimuli is in voltage, it go straight to elmag
        if (abs(appfi_exm(ipoin)) > 0.0_rp) then
           unkno(ipoin)= appfi_exm(ipoin)
        end if
     end if
  end do



  call cputim(cpu_refe2)
  cpu_exmed(1) = cpu_exmed(1) + cpu_refe2 - cpu_refe1 
  cpu_exmed(2) = cpu_exmed(2) + cpu_solve

!!!  write(993,200) unkno(1:npoin)
!200 format(10(1x,f14.8))

end subroutine exm_solite
