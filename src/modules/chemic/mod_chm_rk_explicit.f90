!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



module mod_chm_rk_explicit

  use def_master,          only : ID_CHEMIC, ID_NASTIN, INOTEMPTY, INOTMASTER, ITASK_BEGINN, ITASK_ENDINN, ITASK_INNITE,&
                                  rhsid, unkno, soltyp, mem_modul, modul, postp, dtinv, conce, advec, ittim, itinn, kfl_coupl, solve
  use def_domain,          only : npoin, vmass
  use def_kintyp,          only : ip, rp
  use def_chemic,          only : comax_chm, comin_chm, iclaf_chm, iclai_chm, iclas_chm, imixf_rk, ittot_chm, kfl_avg_cond_CMC_chm,&
                                  kfl_entropy_chm, kfl_goit2_chm, kfl_goite_chm, kfl_mesh_interp_CMC_chm, kfl_model_chm,&
                                  kfl_solve_cond_CMC_chm, kfl_soot_chm, kfl_split_chm, kfl_start_CMC_chm, kfl_tiacc_chm,&
                                  kfl_trans_phs_spc_CMC_chm, nZ_CMC_chm, rtpts_chm, src_Yk_CMC_chm, dt_chm, dt_rho_chm,&
                                  kfl_spray_chm, kfl_entropy_chm, kfl_model_chm, kfl_solve_cond_CMC_chm, nclas_chm, nvar_CMC_chm,&
                                  kfl_fixno_chm, bvess_chm, veloc_CMC_chm
  use mod_solver,          only : solver_lumped_mass_system
  use mod_solver,          only : solver_explicit
  use mod_memory,          only : memory_alloca
  use mod_chm_finiteRate,  only : chm_IntegrateSource_finiteRate
  use mod_chm_entropy,     only : chm_entropy_solution, chm_entropy_memory, chm_entropy_initialization
  use mod_arrays,          only : arrays_number
  
  implicit none

  real(rp),   pointer :: bt(:)         ! Energy RHS
  real(rp),   pointer :: tt(:)         ! Unknown
  real(rp),   pointer :: tt0(:,:)      ! Unknown
  real(rp),   pointer :: rt(:,:,:)     ! Unknown residuals
  real(rp),   pointer :: aux_rt(:,:)   ! Unknown residuals
  real(rp),   pointer :: aux_tt(:,:)   ! Unknown residuals
  real(rp),   pointer :: Mass(:)       ! lumped mass
  real(rp),   pointer :: Tim(:,:,:)    ! projection dt / (rho)
  real(rp)            :: dt(3)

  real (rp)    :: a(5)
  real (rp)    :: b(5,5)
  real (rp)    :: sigma(5)
  integer (ip) :: nequs_loc

  type(soltyp), pointer  :: solve1(:)
  type(soltyp), pointer  :: solve2(:)
  real(rp),     pointer  :: x(:,:)

  private

  public :: chm_rk_explicit_solution
  public :: chm_rk_explicit_initialization
  public :: chm_rk_explicit_memory

contains

  subroutine chm_rk_explicit_initialization()

    nullify(rt)
    nullify(Tim)
    nullify(tt0)

    if (kfl_entropy_chm > 0 ) &
       call chm_entropy_initialization()

  end subroutine chm_rk_explicit_initialization



  subroutine chm_rk_explicit_memory()

    !
    ! Count number of equations solved
    !
    nequs_loc = nclas_chm
    if (kfl_model_chm == 4 .and. kfl_solve_cond_CMC_chm == 1)  nequs_loc = nvar_CMC_chm

    bt       => rhsid
    tt       => unkno
    Mass     => vmass

    if( INOTEMPTY ) then
       call memory_alloca(mem_modul(1:2,modul),'RT' ,'chm_rk_explicit_allocate',rt,nequs_loc,npoin,5_ip)
       call memory_alloca(mem_modul(1:2,modul),'AUX_RT','chm_rk_explicit_allocate',aux_rt,nequs_loc,npoin)
       call memory_alloca(mem_modul(1:2,modul),'AUX_TT','chm_rk_explicit_allocate',aux_tt,nequs_loc,npoin)
       call memory_alloca(mem_modul(1:2,modul),'TIM','chm_rk_explicit_allocate',Tim,nequs_loc,npoin,4_ip)
       call memory_alloca(mem_modul(1:2,modul),'TT0','chm_rk_explicit_allocate',tt0,nequs_loc,npoin)
    else
       call memory_alloca(mem_modul(1:2,modul),'RT' ,'chm_rk_explicit_allocate',rt,1_ip,1_ip,5_ip)
       call memory_alloca(mem_modul(1:2,modul),'AUX_RT','chm_rk_explicit_allocate',aux_rt,1_ip,1_ip)
       call memory_alloca(mem_modul(1:2,modul),'TIM','chm_rk_explicit_allocate',Tim,1_ip,1_ip,4_ip)
       call memory_alloca(mem_modul(1:2,modul),'TT0','chm_rk_explicit_allocate',tt0,1_ip,1_ip)
    end if

    solve1 => solve(1:)
    solve2 => solve(2:)
    x      => conce(:,:,1)

    dt    = 1.0_rp
    a     = 0.0_rp
    b     = 0.0_rp
    sigma = 0.0_rp

    if(kfl_tiacc_chm == 1) then !RK1
       a(4) = 1.0_rp
       sigma(4) = 1.0_rp
    else if(kfl_tiacc_chm == 2) then !Heun’s method or RK2 just for validation
       a(4) = 1.0_rp

       b(4,3) = 1.0_rp

       sigma(3) = 0.5_rp
       sigma(4) = 0.5_rp
    else if(kfl_tiacc_chm == 3) then ! standard 3 order RK
       a(3) = 1.0_rp
       a(4) = 0.5_rp

       b(3,2) = 1.0_rp
       b(4,2) = 1.0_rp/4.0_rp
       b(4,3) = 1.0_rp/4.0_rp

       sigma(2) = 1.0_rp/6.0_rp
       sigma(3) = 1.0_rp/6.0_rp
       sigma(4) = 2.0_rp/3.0_rp
    else if(kfl_tiacc_chm == 4) then ! standard 4 order RK
       a(2) = 0.5_rp
       a(3) = 0.5_rp
       a(4) = 1.0_rp

       b(2,1) = 0.5_rp
       b(3,2) = 0.5_rp
       b(4,3) = 1.0_rp

       sigma(1) = 1.0_rp/6.0_rp
       sigma(2) = 1.0_rp/3.0_rp
       sigma(3) = 1.0_rp/3.0_rp
       sigma(4) = 1.0_rp/6.0_rp
    end if

    if (kfl_entropy_chm > 0_ip ) &
       call chm_entropy_memory()

  end subroutine chm_rk_explicit_memory


  subroutine chm_multi_step_fs_eval(iclai,iclaf,istep)
    integer(ip) , intent(in) :: iclai
    integer(ip) , intent(in) :: iclaf
    integer(ip) , intent(in) :: istep
    integer(ip)              :: ipoin,iclass,kpoin,iclaf_gas,iclai_liq

    external                 :: chm_matrix
    external                 :: chm_partis
    external                 :: rhsmod

    dt_rho_chm = 0.0_rp

    if ( kfl_spray_chm /= 0_ip ) &
       dt_chm = 0.0_rp

    call chm_matrix()

    call solver_lumped_mass_system(1_ip,dt_rho_chm)

    if ( kfl_spray_chm /= 0_ip ) &
         call solver_lumped_mass_system(1_ip,dt_chm)

    !
    ! Coupling with partis
    !
    call chm_partis()

    !
    ! Solve explicit system
    !
    if ( kfl_spray_chm /= 0_ip ) then

       iclaf_gas = iclaf     - 2
       iclai_liq = iclaf_gas + 1

       do ipoin = 1,npoin
          do iclass = iclai,iclaf_gas
             Tim(iclass,ipoin,1) = dt_rho_chm(ipoin)/dt(1)
          end do
          do iclass = iclai_liq,iclaf
             Tim(iclass,ipoin,1) = dt_chm(ipoin)/dt(1)
          end do
       end do

    else
       do ipoin = 1,npoin
          do iclass = iclai,iclaf
             Tim(iclass,ipoin,1) = dt_rho_chm(ipoin)/dt(1)
          end do
       end do
    end if

    kpoin = 0
    do ipoin = 1,npoin
       do iclass = iclai,iclaf
          kpoin             = kpoin + 1
          aux_rt(iclass,ipoin) = bt(kpoin)
       end do
    end do

    call rhsmod(nequs_loc,aux_rt)

    if( INOTMASTER ) then
       do ipoin = 1,npoin
          do iclass = iclai,iclaf
             rt(iclass,ipoin,istep) = aux_rt(iclass,ipoin)
          end do
       end do
    endif

    if(istep == 4_ip) then
       call chm_multi_step_fs_solution_sij(iclai,iclaf,4_ip,sigma)
    else
       call chm_multi_step_fs_solution_sij(iclai,iclaf,istep,b(istep+1,:))
    endif

  end subroutine chm_multi_step_fs_eval


 subroutine chm_multi_step_fs_solution_sij(iclai,iclaf,istep,weight1)
   integer(ip) , intent(in) :: iclai
   integer(ip) , intent(in) :: iclaf
   integer(ip) , intent(in) :: istep
   real(rp) ,    intent(in) :: weight1(5)
   integer(ip)       :: ipoin,jstep,iclass,kpoin
   real(rp) :: aux2

   if( INOTMASTER ) then

      kpoin = 0
      do ipoin = 1,npoin
         do iclass = iclai,iclaf
            aux2 = 0.0_rp
            kpoin     = kpoin + 1
            do jstep=1,istep
               aux2   = aux2 + rt(iclass,ipoin,jstep)*weight1(jstep)
            end do
            !tt(kpoin) = tt0(iclass,ipoin) + ( Tim(iclass,ipoin,1)*aux2*dt(1) ) / Mass(ipoin)
            aux_tt(iclass,ipoin) = Tim(iclass,ipoin,1)*aux2*dt(1)
         end do
      end do

      !
      ! Solve system
      !
      call solver_explicit(solve1,aux_tt,EXCHANGE=.false.,solve_consistent=solve2)

      kpoin = 0
      do ipoin = 1,npoin
         do iclass = iclai,iclaf
            kpoin = kpoin+1
            tt(kpoin) = tt0(iclass,ipoin) + aux_tt(iclass,ipoin)
         end do
      end do

      !
      ! Dirichet bbcc
      !
      kpoin = 0
      do ipoin = 1,npoin
         do iclass = iclai,iclaf
            kpoin = kpoin + 1
            if( kfl_fixno_chm(iclass,ipoin) > 0 ) &
                 tt(kpoin) = bvess_chm(iclass,ipoin)
         end do
      end do

   end if

  end subroutine chm_multi_step_fs_solution_sij


 subroutine chm_rk_explicit_solution()

   use mod_chm_operations_CMC,     only : chm_integrate_chem_source_CMC, &
                                          chm_compute_interp_fields_mesh_CMC, &
                                          chm_local2global_CMC, chm_global2local_CMC, &
                                          chm_limit_Yk_mixt_fraction_slice_CMC, &
                                          chm_soot_transfer_source_terms_CMC


   integer(ip) ::       ipoin, kpoin, iter_rk
   integer(ip), save :: iter = 0
   real(rp)          :: dt_split

   external          :: chm_updbcs
   external          :: chm_updunk
   external          :: chm_endite

   iter = iter + 1

   !
   ! Update inner iteration counter
   !
   itinn(modul) = itinn(modul) + 1
   ittot_chm    = ittot_chm + 1
   rtpts_chm     =  0.0_rp
   comin_chm     =  1.0e9_rp
   comax_chm     = -1.0e9_rp
   kfl_goit2_chm =  0_ip

   !
   ! Allocate memory
   !
   bt       => rhsid
   tt       => unkno
   Mass     => vmass

   !
   ! Assemble equations:
   !
   dt(3) = dt(2)
   dt(2) = dt(1)
   dt(1) = 1.0_rp/dtinv


   iclai_chm = 1
   iclaf_chm = nequs_loc

   if( INOTMASTER) then
      !
      ! Integrate source terms
      !
      if(kfl_split_chm > 1_ip) then
         dt_split = dt(1) / dble(kfl_split_chm)
         if (kfl_model_chm == 3) then
             call chm_IntegrateSource_finiteRate(dt_split)
         else if (kfl_model_chm == 4 .and. kfl_solve_cond_CMC_chm == 1) then
             call chm_integrate_chem_source_CMC(dt_split)
         end if
      end if
   end if

   !
   ! Runge-Kutta loop
   !
   if (kfl_model_chm == 4 .and. kfl_solve_cond_CMC_chm == 1) then
      !
      ! Conditional CMC equations
      !

      if (INOTEMPTY)  src_Yk_CMC_chm = 0.0_rp

      ! Update boundary conditions
      if (kfl_start_CMC_chm == 0)  call chm_updbcs(1)

      if (kfl_trans_phs_spc_CMC_chm == 1_ip) then

         ! Assign velocity in case Nastin is run together with CMC
         if(kfl_coupl(ID_NASTIN,ID_CHEMIC) /= 0) then
            do ipoin=1,npoin
               veloc_CMC_chm(:,ipoin) = advec(:,ipoin,1)
            end do
         end if

         mixt_fr: do imixf_rk = 2,nZ_CMC_chm-1_ip

            if( INOTMASTER) then
               ! Define unkno(ipoin) from Yk_CMC_chm and enthalp_CMC_chm for CMC model
               call chm_global2local_CMC(imixf_rk)
               call chm_updbcs(2)
               call chm_updunk(ITASK_BEGINN)
               if (kfl_mesh_interp_CMC_chm == 1_ip .and. kfl_avg_cond_CMC_chm == 1_ip) then
                  call chm_compute_interp_fields_mesh_CMC(imixf_rk)
               end if

               kpoin = 0
               do ipoin = 1,npoin
                  do iclas_chm = iclai_chm, iclaf_chm
                     kpoin                = kpoin + 1
                     tt0(iclas_chm,ipoin) = tt(kpoin)
                  end do
               end do
            end if

            !
            ! Entropy prediction:
            !
            if(kfl_entropy_chm == 1_ip) call chm_entropy_solution()
            !!!!!!!!!!!!!!! VER SI AQUÍ HAY QUE CAMBIAR ALGO PARA EL CMC


            do iter_rk = 5_ip-kfl_tiacc_chm, 3_ip
               !
               ! Runge Kutta stage
               !
               call chm_multi_step_fs_eval(iclai_chm,iclaf_chm,iter_rk)

               !
               ! Clipping + update conce(:,:,1)
               !
               call chm_endite(ITASK_INNITE)
            end do
            !
            ! Last Runge Kutta stage
            !
            call chm_multi_step_fs_eval(iclai_chm,iclaf_chm,4_ip)


            call chm_endite(ITASK_ENDINN)
            ! Note: if enthalpy is not transported it does not change with time -> enthalpy
            ! matrix already filled during the initialization

            call chm_limit_Yk_mixt_fraction_slice_CMC
            call chm_local2global_CMC(imixf_rk)
            if (postp(1) % npp_stepi(arrays_number('MASUN'),0) /= 0) then
               ! arrays_number('MASUN') is for unconditional mass fractions
               if( kfl_soot_chm > 0 .and. mod(ittim, postp(1) % npp_stepi(arrays_number('MASUN'),0) ) == 0_ip ) &
                  call chm_soot_transfer_source_terms_CMC(imixf_rk)
            end if

         end do mixt_fr

      else
         kfl_goite_chm = 0_ip
      end if

   else

      !
      ! Normal equations
      !
      if( INOTMASTER) then
         !
         ! Define unkno(ipoin) from conce(:,iclas,1) for all iclas
         !
         call chm_updunk(ITASK_BEGINN)

         kpoin = 0
         do ipoin = 1,npoin
            do iclas_chm = iclai_chm, iclaf_chm
               kpoin                = kpoin + 1
               tt0(iclas_chm,ipoin) = tt(kpoin)
            end do
         end do
      end if

      !
      ! Entropy prediction:
      !
      if(kfl_entropy_chm == 1_ip) call chm_entropy_solution()


      if(     kfl_tiacc_chm == 2) then

         !
         ! Runge Kutta stage
         !
         call chm_multi_step_fs_eval(iclai_chm,iclaf_chm,3_ip)

         !
         ! Clipping + update conce(:,:,1)
         !
         call chm_endite(ITASK_INNITE)

      else if(kfl_tiacc_chm == 3) then

         if (kfl_split_chm == 2_ip .and. kfl_model_chm == 3_ip ) then
            dt_split = 0.5_rp * dt(1)
            call chm_IntegrateSource_finiteRate(dt_split)
         end if

         !
         ! Runge Kutta stage
         !
         call chm_multi_step_fs_eval(iclai_chm,iclaf_chm,2_ip)

         !
         ! Clipping + update conce(:,:,1)
         !
         call chm_endite(ITASK_INNITE)

         if (kfl_split_chm == 3_ip .and. kfl_model_chm == 3_ip ) then
            dt_split = 1.0_rp * dt(1)
            call chm_IntegrateSource_finiteRate(dt_split)
         end if

         !
         ! Runge Kutta stage
         !
         call chm_multi_step_fs_eval(iclai_chm,iclaf_chm,3_ip)

         !
         ! Clipping + update conce(:,:,1)
         !

         call chm_endite(ITASK_INNITE)

         if (kfl_split_chm == 3_ip .and. kfl_model_chm == 3_ip ) then
            dt_split = 0.5_rp * dt(1)
            call chm_IntegrateSource_finiteRate(dt_split)
         end if

      else ! order 4

         !
         ! Runge Kutta stage
         !
         call chm_multi_step_fs_eval(iclai_chm,iclaf_chm,1_ip)

         !
         ! Clipping + update conce(:,:,1)
         !
         call chm_endite(ITASK_INNITE)

         !
         ! Runge Kutta stage
         !
         call chm_multi_step_fs_eval(iclai_chm,iclaf_chm,2_ip)

         !
         ! Clipping + update conce(:,:,1)
         !
         call chm_endite(ITASK_INNITE)

         !
         ! Runge Kutta stage
         !
         call chm_multi_step_fs_eval(iclai_chm,iclaf_chm,3_ip)

         !
         ! Clipping + update conce(:,:,1)
         !
         call chm_endite(ITASK_INNITE)

      end if

      !
      ! Last Runge Kutta stage
      !
      call chm_multi_step_fs_eval(iclai_chm,iclaf_chm,4_ip)

      call chm_endite(ITASK_ENDINN)

   end if

 end subroutine chm_rk_explicit_solution

end module mod_chm_rk_explicit

