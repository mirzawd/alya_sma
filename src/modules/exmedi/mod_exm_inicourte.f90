!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup Exmedi
!> @{
!> @file    mod_exmedi_inicourte.f90
!> @author  Eva Casoni 
!> @date    10/12/2020
!> @brief   Initial condition setup for Courtemanche atrial  model
!> @details 
!> @} 
!!-----------------------------------------------------------------------
module mod_exm_inicourte

  use def_master
  use def_exmedi
  use mod_messages,      only: messages_live
  use mod_eccoupling,    only: kfl_exmsld_ecc, EXMSLD_EMEC_LAND, EXMSLD_EMEC_LAND_BIDIR
  use mod_parall,        only: PAR_MY_CODE_RANK
  use mod_memory,        only: memory_alloca, memory_deallo
  use mod_exm_cellmodel
  use mod_exm_courtemanche_model, only: exm_courtemanche_model

  implicit none

  private

  integer(ip), parameter :: &
     EXM_OCECOURTE_PURE = 0_ip, &   ! Run Courtemanche model
     EXM_OCECOURTE_LAND = 1_ip      ! Run Courtemanche coupled to Land

  real(rp) :: &
     elmlo_exm(3)         ! New local voltage

  real(rp), dimension(:,:), pointer :: &
     vaulo_exm => null(),&            ! Auxiliary variables for single cell solution
     vcolo_exm => null(),&            ! Concentrations for single cell solution
     viclo_exm => null()             ! Currents for single cell solution
     

  public :: exm_onecourte
 ! private :: exm_init_voltages, localvoltages2globalvoltage, exm_courte_allocate_locals, exm_courte_deallocate_locals


contains


subroutine exm_courte_allocate_locals()
   use mod_memory, only : memory_alloca
   implicit none

   call memory_alloca(mem_modul(1:2,modul),'vaulo_exm', 'mod_exm_onecourte', vaulo_exm,  nauxi_exm, 3_ip) 
   call memory_alloca(mem_modul(1:2,modul),'vcolo_exm', 'mod_exm_onecourte', vcolo_exm,  nconc_exm, 3_ip) 
   call memory_alloca(mem_modul(1:2,modul),'viclo_exm', 'mod_exm_onecourte', viclo_exm,  27_ip, 3_ip) 
   vaulo_exm = 0.0_rp
   vcolo_exm = 0.0_rp
   viclo_exm = 0.0_rp

end subroutine exm_courte_allocate_locals


subroutine exm_courte_deallocate_locals()
   use mod_memory, only : memory_deallo
   implicit none

   call memory_deallo(mem_modul(1:2,modul),'vaulo_exm', 'mod_exm_onecourte', vaulo_exm) 
   call memory_deallo(mem_modul(1:2,modul),'vcolo_exm', 'mod_exm_onecourte', vcolo_exm) 
   call memory_deallo(mem_modul(1:2,modul),'viclo_exm', 'mod_exm_onecourte', viclo_exm) 
end subroutine exm_courte_deallocate_locals



  !------------------------------------------------------------------------
  !> @addtogroup Exmedi
  !> @{
  !> @file    exm_onecourte.f90
  !> @date    11/12/2020
  !> @author  Eva Casoni
  !> @brief   Sets initial conditions for Courtemanche model
  !> @details 
  !> @}
  !------------------------------------------------------------------------
  ! icelltype -- epi, endo, mid (1-3)
  ! Returns:
  ! success_status == 0_ip -- sucess, >0 -- error, see EXM_ONEOHR_ERROR* for the meaning
   subroutine exm_onecourte(mat, icelltype, outputs)
      use      mod_eccoupling

      implicit none
      integer(ip), intent(in)                   :: mat, icelltype
      TYPE(CELL_MODEL_OUTPUTS), intent(inout)       :: outputs
      integer(ip)                               :: success_status
     ! real(rp), intent(inout)                   :: courte_stats(3_ip)  !saves Courtemanche stats: number of beats, tolerance
      !real(rp), dimension(1:7), intent(out)     :: land_variables

      !integer(ip)                               :: success_status_tmp
      logical                                   :: flag_land, initialized

      call exm_courte_allocate_locals()

      success_status = 0_ip

      flag_land = .FALSE.
      initialized = .FALSE. !flag to check if any of the IF was executed

      if (kfl_exmsld_ecc) then
         if (kfl_eccty(mat) == EXMSLD_EMEC_LAND .or. kfl_eccty(mat) == EXMSLD_EMEC_LAND_BIDIR) flag_land = .TRUE.
      end if

      call exm_init_voltages(mat, icelltype)

  
      if (kfl_exmsld_ecc .and. flag_land) then
 
        if (kfl_user_specified_celltypes_exm(mat, icelltype) .eq. 1_ip) then
            call exm_ocecourte(mat, icelltype, EXM_OCECOURTE_LAND, outputs)  
            success_status = success_status + outputs % success
         end if

      else
        if (kfl_user_specified_celltypes_exm(mat, icelltype) .eq. 1_ip) then
            call exm_ocecourte(mat, icelltype, EXM_OCECOURTE_PURE, outputs)  
            success_status = success_status + outputs % success
        end if

      end if

      call localvoltages2globalvoltages(mat, icelltype)  !EVA: revisar si hace falta

      initialized = .TRUE.

      if (success_status > 0_ip) then
         success_status = EXM_CELLMODEL_NOTCONVERGED
      end if

      if (.NOT. initialized) then
         success_status = EXM_CELLMODEL_NOTINITIALIZED
      end if 

      call exm_courte_deallocate_locals()

   end subroutine exm_onecourte


  
  subroutine localvoltages2globalvoltages(mat, ituss_exm)
     implicit none

     integer(ip), intent(in) :: mat, ituss_exm

     vauxi_exm_initial(:, ituss_exm, mat) = vaulo_exm(:, 1)
     vconc_initial(:, ituss_exm, mat) = vcolo_exm(:, 1)
     vminimate_exm(ituss_exm, mat) = elmlo_exm(2)

  end subroutine

  subroutine exm_init_voltages(mat, celltype)
     implicit none

     integer(ip), intent(in) :: mat, celltype

     vminimate_exm(celltype,mat) = -81.18_rp  !membrane [mV]

     ! Gate variables (dimensionless)
     vauxi_exm_initial(1,celltype,mat) = 0.0_rp        !Variable u 
     vauxi_exm_initial(2,celltype,mat) = 1.0_rp              !Variables v 
     vauxi_exm_initial(3,celltype,mat) = 0.9992_rp           !Variable w  
     vauxi_exm_initial(4,celltype,mat) = 1.367e-4_rp         !Variable d  
     vauxi_exm_initial(5,celltype,mat) = 7.755e-1_rp         !Variable fCa 
     vauxi_exm_initial(6,celltype,mat) = 9.996e-1_rp         !Variable f 
     vauxi_exm_initial(7,celltype,mat) = 9.649e-1_rp         !Variable h 
     vauxi_exm_initial(8,celltype,mat) = 9.775e-1_rp         !Variable j 
     vauxi_exm_initial(9,celltype,mat) = 2.908e-3_rp         !Variable m
     vauxi_exm_initial(10,celltype,mat) = 3.296e-5_rp        !Variable xr 
     vauxi_exm_initial(11,celltype,mat) = 1.869e-2_rp        !Variable xs
     vauxi_exm_initial(12,celltype,mat) = 3.043e-2_rp        !Variable oa
     vauxi_exm_initial(13,celltype,mat) = 9.992e-1_rp        !Variable oi
     vauxi_exm_initial(14,celltype,mat) = 4.966e-3_rp        !Variable ua
     vauxi_exm_initial(15,celltype,mat) = 9.986e-1_rp        !Variable ui
      
     ! Intracellular concentrations
     vconc_initial(1,celltype,mat) = 1.013e-4_rp         !Variable Cai [mM]
     vconc_initial(2,celltype,mat) = 1.488_rp            !Variable Carel [mM]
     vconc_initial(3,celltype,mat) = 1.488_rp            !Variable Caup [mM]
     vconc_initial(4,celltype,mat) = 1.39e2_rp           !Variable Ki [mM]
     vconc_initial(5,celltype,mat) = 1.117e1_rp          !Variable Nai [mM]

     ! Membrane currents
     viclo_exm(1:26,1:2) = 0.0_rp

     ! Initialise single cell variables
     elmlo_exm(1:2) = vminimate_exm(celltype,mat)


     ! Initialise single cell variables
     viclo_exm(:,1:2) = 0.0_rp
     vcolo_exm(:,1) = vconc_initial(:, celltype, mat)
     vcolo_exm(:,2) = vconc_initial(:, celltype, mat)
     elmlo_exm(1:2) = vminimate_exm(celltype,mat)
     vaulo_exm(:,1) = vauxi_exm_initial(:, celltype, mat)
     vaulo_exm(:,2) = vauxi_exm_initial(:, celltype, mat)

   end subroutine exm_init_voltages

   !-------------------------------------------------------------------
   ! Single cell run fo rinitial condition setup for Courtemanche
   ! This function initialize the single cell model until it achieves
   ! a steady state or a stable situation. Then the variables at this point
   ! are the initial conditions for the run (made in doiter)
   ! mat - material id
   ! ituss_exm: in this case is always 1
   ! -----------------------------------------------------------------

   subroutine exm_ocecourte(mat, ituss_exm, model_type, outputs)
      use mod_exm_cellmodel_convergence
      use mod_exm_drugs, only  : exm_drugs_ondrugs  !, exm_drugs_get_drugdmate
      use def_exmedi,    only  : kfl_isac_exm
      use mod_eccoupling, only : kfl_force_ca50_calculation

      ! definiton of variables 
      implicit none
      integer(ip), intent(in)                :: mat, ituss_exm, model_type
      TYPE(CELL_MODEL_OUTPUTS), intent(inout)    :: outputs

      integer(ip)   :: nbeat, j, itim, nsamples_per_beat
      real(rp)      :: qnet, cai
!      real(rp)      :: Inet
      real(rp)      :: dtimon
      real(rp)      :: aux1
      real(rp)      :: stim, batim, atim, dur, durcycle, numStim
!      integer(ip)   :: p
      real(rp)      :: Y(21)

      !Land variables 
      real(rp)    :: k_TRPN, n_TRPN, Ca50_ref, k_u, n_tm, TRPN_50, k_uw, k_ws, r_w, r_s, y_s, y_w, phi
      real(rp)    :: A_eff, beta_0, beta_1, T_ref, k_wu, k_su, A_s, c_s, c_w, k_u_2, lambda, lambda_rate
      real(rp)    :: S, W, CaTRPN, B, zeta_s, zeta_w, Ca50, Ca_max, Ca_min, y_su, y_wu, ydot(6)

      !number of voltages actually saved in the array, needed to verify that there was no error
      integer(ip) :: ss_array_current_id
      real(rp)    :: ss_epsilon ! threshold squared to declare steady state when comparing two beats
      real(rp)    :: ss_min_voltage_threshold
   
      integer(ip)                       :: post_frequency, post_file_handle
      logical                           :: save_variables, flag_land, flag_3D, flag_isac
      real(rp)                          :: statvar(7,2)

      type(ConvergenceChecker)          :: convergence

      !to save calcium and voltage and other variables
      save_variables = ( kfl_save_convergence_cellmodel == 1_ip )
      post_frequency = 500_ip
   
      post_file_handle = -1_ip
      Ca50 = 0.0_rp
      
      !sanity check
      SELECT CASE (model_type)
      CASE (EXM_OCECOURTE_PURE)
         ! nothing
      CASE (EXM_OCECOURTE_LAND)
         ! nothing
      CASE DEFAULT
         call runend("mod_exm_inicourte: Unknown cell model: "//trim(intost(model_type)))
      END SELECT
   
   
      !write file header
      !DEBUG-EVA
      save_variables = .true.
      if (save_variables) then
   
         post_file_handle = 1000000 + 100*mat + ituss_exm
   
         if (IPARALL) post_file_handle = post_file_handle + PAR_MY_CODE_RANK*10000

         SELECT CASE (model_type)
         CASE (EXM_OCECOURTE_PURE)
            WRITE (post_file_handle,*) "Global_Time Beat_Time Voltage Calcium Sodium"
         CASE (EXM_OCECOURTE_LAND)
            WRITE (post_file_handle,*) "Global_Time Beat_Time Calcium Voltage CaTRPN Qnet"
         CASE DEFAULT
            call runend("mod_exm_oceohr: Unknown cell model: "//trim(intost(model_type)))
         END SELECT
   
     end if

      if (model_type == EXM_OCECOURTE_LAND) then
         flag_land = .TRUE.
      else
         flag_land = .FALSE.
      end if
 
      flag_3D = .FALSE. 

      if (kfl_isac_exm(mat) == 1_ip) then
         flag_isac = .TRUE.
      else
         flag_isac = .FALSE.
      end if
 

      ! Varibles for beats and time integration of edos
      dtimon = 0.01_rp ! (ms)
      dur = 2.0_rp  ! Duration of stimulus (ms)  -it should be read from initial conditions
      nbeat = 0
      itim = 0
      batim = 0
      j = 0
      atim = 0.0_rp
      stim = 0.0_rp

      !---------------------------------------------------
      ! Land constants: REVIEW FOR COURTEMANCHE
      !-------------------------------------------------
      k_TRPN = 0.1_rp
      n_TRPN = 2.0_rp
      Ca50_ref = 0.86_rp!0.805_rp !solidz uses different units, Ca in solidz = Ca in exmedi * 1000
      k_u = 1.0_rp
      n_tm = 5.0_rp
      TRPN_50 = 0.35_rp
      k_uw = 3.0_rp*0.182_rp
      k_ws = 3.0_rp*0.012_rp
      r_w = 0.5_rp
      r_s = 0.25_rp
      y_s = 0.0085_rp
      y_w = 0.615_rp
      phi = 2.23_rp
      A_eff = 25.0_rp
      beta_0 = 2.3_rp
      beta_1 = -2.4_rp
      T_ref = 120.0_rp
      k_wu = k_uw*((1.0_rp/r_w) - 1.0_rp) - k_ws
      k_su = k_ws*((1.0_rp/r_s) - 1.0_rp)*r_w
      A_s = (A_eff*r_s)/((1.0_rp - r_s)*r_w + r_s)
      c_s = phi*k_ws*(1.0_rp - r_s)*r_w/r_s
      c_w = phi*k_uw*((1.0_rp - r_s)*(1.0_rp - r_w))/((1.0_rp - r_s)*r_w)
      k_u_2 = (k_u*(TRPN_50**n_tm))/(1.0_rp - r_s - (1.0_rp - r_s)*r_w)

      !Declare initial variables
      lambda = 1.0_rp
      lambda_rate = 0.0_rp
      Ca50 = Ca50_ref + beta_1*(lambda - 1.0_rp)

      ! Declare initial values of state variables
      S = 0.0_rp
      W = 0.0_rp
      CaTRPN = 0.000001_rp
      B = 1.0_rp
      zeta_s = 0.0_rp
      zeta_w = 0.0_rp

      !---------------------------------------------------
      ! End of Land constants
      !-------------------------------------------------
    

 
      !----------------------------------------------
      !
      ! begin: Initialize variables for steady state detection
      !
      !----------------------------------------------
      !allocate the arrays to store two full beats
      nsamples_per_beat = nint(moneclmate_exm(2, mat)/dtimon)
   
      if (steady_tol_cellm(mat) >= 0) then
         ss_epsilon = steady_tol_cellm(mat)
      else
         select case (kfl_steadystate_variable(mat))
         case (EXM_CELL_STEADY_VOLTAGE)
            ss_epsilon = 1.0_rp    !this is maximum allowed  mean voltage difference
         case (EXM_CELL_STEADY_CALCIUM)
            ss_epsilon = 1.0e-7_rp 
         case default
            call runend("EXM_OCEOHR: Unknown name of the variable to determine the cell convergence.")
         end select
      end if
      call convergence % Init( ss_epsilon, nsamples_per_beat )
      call convergence % SetMaxSteps( 3_ip )
   
      ss_min_voltage_threshold = -75.0_rp !the min volytage within the beat has to be lower than this, to consider that we are in steady state
      !----------------------------------------------
      !
      ! end: Initialize variables for steady state detection
      !
      !----------------------------------------------
   
      outputs % success = EXM_CELLMODEL_NOTCONVERGED !return error by default
   
      ! moneclmate_exm(1, mat) - number of beats
      ! moneclmate_exm(2, mat) - cycle length
      
      numStim = moneclmate_exm(1,mat)  !number of beats
      durcycle = moneclmate_exm(2,mat) !length in seconds   

      dtimon = 0.01_rp


      do while (batim < (moneclmate_exm(1, mat)*moneclmate_exm(2, mat)))

         Ca_max = 0.0_rp !reset Ca max for the next beat
         Ca_min = 1.0e+30_rp !reset
  
         aux1 = 0.0_rp
         atim = dtimon*aux1
         dur = 2.0_rp  !duration time (ms)
   
        !for steady state verification
         ss_array_current_id = 1_ip
         qnet = 0.0_rp
  !!       Inet = 0.0_rp
         do while (atim < moneclmate_exm(2, mat))
      
            aux1 = aux1 + 1.0_rp
            atim = dtimon*aux1

            if (atim < dur) then
                stim = -2.0_rp*754.0_rp  !Value of I-stim current [pA/pF], magnitud del estimulo
            else 
                stim = 0.0_rp
            end if 
            viclo_exm(18,1) = stim

            call  exm_courtemanche_model(ituss_exm,elmlo_exm(1),vcolo_exm,vaulo_exm,viclo_exm,dtimon,qnet,flag_land,flag_3D,flag_isac,statvar)

            Y(1:9) = vaulo_exm(1:9,2)
            Y(10:14) = vcolo_exm(1:5,2)
            Y(15) = elmlo_exm(2)
            Y(16:21) = vaulo_exm(10:15,2)


           !------------------------------------------------------
           ! 
           !   Ca50 stuff
           !
           if ( Ca_max < vcolo_exm(1, 1) * 1000.0_rp) then 
              Ca_max = vcolo_exm(1, 1) * 1000.0_rp !change of units
           elseif ( Ca_min > vcolo_exm(1, 1) * 1000.0_rp) then
              Ca_min = vcolo_exm(1, 1) * 1000.0_rp !change of units
           endif
           ! 
           !   Ca50 stuff
           !
           !------------------------------------------------------


            !
            ! Beginning of Land
            !
            if (model_type == EXM_OCECOURTE_LAND) then 
               !------------------------------------------------------
               ! 
               !   Beginning of Land
               !
               !----------------------------------------------------
   
               ! State variable dependent parameters (regularised parameters)
               if ((zeta_s + 1.0_rp) < 0.0_rp) then
                  y_su = y_s*(-zeta_s - 1.0_rp)
               elseif ((zeta_s + 1.0_rp) > 1.0_rp) then
                  y_su = y_s*zeta_s
               else
                  y_su = 0.0_rp
               end if
               y_wu = y_w*ABS(zeta_w)

               ! Define the RHS of the ODE system of Land
               ydot(1) = k_ws*W - k_su*S - y_su*S
               ydot(2) = k_uw*(1.0_rp - B - S - W) - k_wu*W - k_ws*W - y_wu*W
               ydot(3) = k_TRPN*(((vcolo_exm(1, 2)*1000.0_rp/Ca50)**n_TRPN)*(1.0_rp - CaTRPN) - CaTRPN)
               ydot(4) = k_u_2*MIN((CaTRPN**(-n_tm/2.0_rp)), 100.0_rp)*(1.0_rp - B - S - W) - k_u*(CaTRPN**(n_tm/2.0_rp))*B
               ydot(5) = A_s*lambda_rate - c_s*zeta_s
               ydot(6) = A_s*lambda_rate - c_w*zeta_w
   
               ! Update variables of Land
               S = S + dtimon*ydot(1)
               W = W + dtimon*ydot(2)
               CaTRPN = CaTRPN + dtimon*ydot(3)
               B = B + dtimon*ydot(4)
               zeta_s = zeta_s + dtimon*ydot(5)
               zeta_w = zeta_w + dtimon*ydot(6)

            end if
        !----------------------------------------------------------------------------------
        ! Outputs
        !----------------------------------------------------------------------------------
         cai = vcolo_exm(1,2)

        itim = itim + 1

        if (save_variables) then
           if (j >= post_frequency) then  !it was 10 for saving
   
              SELECT CASE (model_type)
                  CASE (EXM_OCECOURTE_PURE)
                     write (post_file_handle, "(6(e16.8E3, ' '))") itim*dtimon, atim, elmlo_exm(2), cai, viclo_exm(1,1)
                  CASE (EXM_OCECOURTE_LAND)
                     write (post_file_handle, "(6(e16.8E3, ' '))") itim*dtimon, atim, cai, elmlo_exm(2), CaTRPN, qnet
                  CASE DEFAULT
                     call runend("mod_exm_inicourte: Unknown cell model: "//trim(intost(model_type)))
                  END SELECT
          
                  j = 0
            end if
            j = j + 1
        end if !save_variables

        elmlo_exm(1) = elmlo_exm(2)
        vcolo_exm(:,1) = vcolo_exm(:,2)
        vaulo_exm(:,1) = vaulo_exm(:,2)

        select case (kfl_steadystate_variable(mat))
        case (EXM_CELL_STEADY_VOLTAGE)
           convergence % array (ss_array_current_id) = elmlo_exm(2)
        case (EXM_CELL_STEADY_CALCIUM)
           convergence % array (ss_array_current_id) = vcolo_exm(1, 2)
        case default
           call runend("EXM_OCEOHR: Unknown name of the variable to determine the cell convergence.")
        end select

        ss_array_current_id = ss_array_current_id + 1_ip

     end do !atim EVA: aqui veuleve a faltar TODO!!!!



     !batim = atim + durcycle*(p-1)
     batim = itim*dtimon


     !Update Ca50
     if( kfl_exmsld_ecc ) then ! TODO: the condition (.not. has_land) should disappear eventually
        if ( (.not. flag_land) .or. ( kfl_force_ca50_calculation(mat) ) ) then
          Ca50 = (Ca_min+Ca_max) * 0.5_rp 
        end if
     else
         Ca50 = (Ca_min+Ca_max) * 0.5_rp 
     end if



         !--------------------------------------------
         !
         !  begin: check if we reached steady state. If so, break
         !
         !----------------------------------------
         if (ss_array_current_id - 1 .NE. nsamples_per_beat) then
            call runend("EXM_OCECOURTE: Something went wrong in the loop of one beat. The number of timesteps executed is incorrect.")
         end if
      
         if ( convergence % converged() ) then
            outputs % success = EXM_CELLMODEL_CONVERGED !success
            exit
         end if
   
         !otherwise copy the beat and go to the beginning
         call convergence % EndIteration()         
         !--------------------------------------------
         !
         !  end: check if we reached steady state
         !
         !----------------------------------------

   end do !end of do while of batim < numberofbeats

      !------------------------------------------------
      !
      ! Save stats
      !
      !----------------------------------------------
      outputs % nbeats = itim/nsamples_per_beat
      outputs % toler = convergence % GetTolerance()
      outputs % rmse  = convergence % GetMetricValue()

      !-----------------------------------------------
      !
      !   Extra outputs
      !
      outputs % S      = S      ! land only
      outputs % W      = W      ! land only
      outputs % CaTRPN = CaTRPN * 1000.0_rp ! land only change from ms to s
      outputs % B      = B      ! land only
      outputs % zeta_s = zeta_s ! land only
      outputs % zeta_w = zeta_w ! land only
      outputs % Lambda = 1.0_rp ! not set here, introduced in TOR  
      outputs % Ca50   = Ca50 / 1000.0_rp   ! general. Denormalize, so that it has the same units as Ca0, and can be renormalized again in mod_exm_sld_eccoupling.f90
      outputs % dt     = dtimon / 1000.0_rp! ms -> s 


   end subroutine exm_ocecourte

end module mod_exm_inicourte
