!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Exmedi
!> @{
!> @file    mod_exm_onetor.f90
!> @author  mixed
!> @date    2019-11-16
!> @brief   mod_exm_onetor
!> @details mod_exm_onetor
!> @}
!-----------------------------------------------------------------------
module mod_exm_onetor

   use def_master
   use def_exmedi
   use mod_eccoupling
   use mod_messages, only: messages_live
   use mod_parall, only: PAR_MY_CODE_RANK
   use mod_memory, only: memory_alloca
   use mod_memory, only: memory_deallo
   use mod_exm_torord_model, only: exm_torord_model
   use mod_exm_cellmodel

   implicit none
   private

   integer(ip), parameter  :: &
      EXM_OCETOR_PURE = 0_ip, &  !run Pure ToR-ORd model
      EXM_OCETOR_LAND = 1_ip     !run ToR-ORd coupled to Land

   real(rp) :: &
      elmlo_exm(3)               ! New local voltage


   real(rp), dimension(:,:), pointer :: &
      vaulo_exm => null(),&            ! Auxiliary variables for single cell solution
      vcolo_exm => null(),&            ! Concentrations for single cell solution
      viclo_exm => null()              ! Currents for single cell solution


   public :: exm_onetor
   
contains



subroutine exm_tor_allocate_locals()
   use mod_memory, only : memory_alloca
   implicit none

   call memory_alloca(mem_modul(1:2,modul),'vaulo_exm', 'mod_exm_ohararudy', vaulo_exm,  nauxi_exm, 3_ip) 
   call memory_alloca(mem_modul(1:2,modul),'vcolo_exm', 'mod_exm_ohararudy', vcolo_exm,  nconc_exm, 3_ip) 
   call memory_alloca(mem_modul(1:2,modul),'viclo_exm', 'mod_exm_ohararudy', viclo_exm,  27_ip, 3_ip) 
   vaulo_exm = 0.0_rp
   vcolo_exm = 0.0_rp
   viclo_exm = 0.0_rp
end subroutine exm_tor_allocate_locals


subroutine exm_tor_deallocate_locals()
   use mod_memory, only : memory_deallo
   implicit none

   call memory_deallo(mem_modul(1:2,modul),'vaulo_exm', 'mod_exm_ohararudy', vaulo_exm) 
   call memory_deallo(mem_modul(1:2,modul),'vcolo_exm', 'mod_exm_ohararudy', vcolo_exm) 
   call memory_deallo(mem_modul(1:2,modul),'viclo_exm', 'mod_exm_ohararudy', viclo_exm) 
end subroutine exm_tor_deallocate_locals




subroutine exm_onetor(mat, icelltype, land_variables, success_status, torord_stats)

   implicit none
   integer(ip), intent(in)  :: mat, icelltype
   integer(ip), intent(out) :: success_status
   integer(ip)              :: success_status_tmp
   logical                  :: flag_land, initialized
   real(rp), dimension(1:7), intent(out) :: land_variables

   real(rp), intent(inout)  :: torord_stats(3_ip) ! saves ToR-ORd states: number of beats, tolerance


   call exm_tor_allocate_locals()

   success_status = 0_ip

   flag_land = .FALSE.
   initialized = .FALSE. !flag to check if any of the IF was executed


   if (kfl_exmsld_ecc) then
      if (kfl_eccty(mat) == 3_ip .or. kfl_eccty(mat) == 4_ip) flag_land = .TRUE.
   end if

   ! No stored initial states at the moment, steady-state needs to be calculated for every simulation.
   call exm_init_voltages(mat, icelltype)
   if (((coupling('SOLIDZ', 'EXMEDI') >= 1_ip) .or. (coupling('EXMEDI', 'SOLIDZ') >= 1_ip)) .and. flag_land) then

      if (kfl_user_specified_celltypes_exm(mat, icelltype) .eq. 1_ip) then
         call exm_ocetla(mat, icelltype, land_variables, success_status_tmp, torord_stats)
         success_status = success_status + success_status_tmp
      end if
   else
      if (kfl_user_specified_celltypes_exm(mat, icelltype) .eq. 1_ip) then
         call exm_ocetor(mat, icelltype, success_status_tmp, torord_stats)
         success_status = success_status + success_status_tmp
      end if
   end if

   call localvoltages2globalvoltages(mat, icelltype)
   initialized = .TRUE.

   if (success_status > 0_ip) then
      success_status = EXM_CELLMODEL_NOTCONVERGED
   end if

   if (.NOT. initialized) then
      success_status = EXM_CELLMODEL_NOTINITIALIZED
   end if

   call exm_tor_deallocate_locals()
end subroutine exm_onetor

subroutine localvoltages2globalvoltages(mat, ituss_exm)
   implicit none

   integer(ip), intent(in) :: mat, ituss_exm

   if (ituss_exm == EXM_CELLTYPE_EPI) then ! epi
      vauxi_exm_initial(:, ituss_exm, mat) = vaulo_exm(:, 1)
      vconc_initial(:, ituss_exm, mat) = vcolo_exm(:, 1)
      vminimate_exm(ituss_exm, mat) = elmlo_exm(2)
   else if (ituss_exm == EXM_CELLTYPE_ENDO) then ! endo
      vauxi_exm_initial(:, ituss_exm, mat) = vaulo_exm(:, 1)
      vconc_initial(:, ituss_exm, mat) = vcolo_exm(:, 1)
      vminimate_exm(ituss_exm, mat) = elmlo_exm(2)
   else if (ituss_exm == EXM_CELLTYPE_MID) then ! mid
      vauxi_exm_initial(:, ituss_exm, mat) = vaulo_exm(:, 1)
      vconc_initial(:, ituss_exm, mat) = vcolo_exm(:, 1)
      vminimate_exm(ituss_exm, mat) = elmlo_exm(2)
   end if
end subroutine

subroutine exm_init_voltages(mat, celltype)
   implicit none

   integer(ip), intent(in) :: mat, celltype

   vminimate_exm(celltype,mat) = -88.76_rp

   if (celltype == EXM_CELLTYPE_EPI) then
      ! Intracellular concentrations
      vconc_initial(1, celltype,mat) = 12.8363_rp   ! [Nai] = mM
      vconc_initial(2, celltype,mat) = 12.8366_rp   ! [Nass] = mM
      vconc_initial(3, celltype,mat) = 142.6951_rp  ! [Ki] = mM
      vconc_initial(4, celltype,mat) = 142.6951_rp  ! [Kss] = mM
      vconc_initial(5, celltype,mat) = 6.6309e-5_rp ! [Cai] = mM
      vconc_initial(6, celltype,mat) = 5.7672e-5_rp ! [Cass] = mM
      vconc_initial(7, celltype,mat) = 1.8119_rp    ! [Cansr] = mM
      vconc_initial(8, celltype,mat) = 1.8102_rp    ! [Cajsr] = mM
      vconc_initial(9, celltype,mat) = 24.0_rp      ! [Cli] = mM
      vconc_initial(10,celltype,mat) = 0.0129_rp   ! [CaMKt] = mM
      vconc_initial(11,celltype,mat) = 0.15_rp     ! [KmCaMK] = mM
      vconc_initial(12,celltype,mat) = 0.0015_rp   ! [KmCaM] = mM
      vconc_initial(13,celltype,mat) = 0.0_rp      ! [Jrel_np] = mM/ms
      vconc_initial(14,celltype,mat) = 0.0_rp      ! [Jrel_p] = mM/ms
   else if (celltype == EXM_CELLTYPE_MID) then
      vconc_initial(1, celltype,mat) = 15.0038_rp   ! [Nai] = mM
      vconc_initial(2, celltype,mat) = 15.0043_rp   ! [Nass] = mM
      vconc_initial(3, celltype,mat) = 143.0403_rp  ! [Ki] = mM
      vconc_initial(4, celltype,mat) = 143.0402_rp  ! [Kss] = mM
      vconc_initial(5, celltype,mat) = 8.166e-05_rp ! [Cai] = mM
      vconc_initial(6, celltype,mat) = 6.5781e-05_rp ! [Cass] = mM
      vconc_initial(7, celltype,mat) = 1.9557_rp    ! [Cansr] = mM
      vconc_initial(8, celltype,mat) = 1.9593_rp    ! [Cajsr] = mM
      vconc_initial(9, celltype,mat) = 24.0_rp      ! [Cli] = mM
      vconc_initial(10,celltype,mat) = 0.0192_rp   ! [CaMKt] = mM
      vconc_initial(11,celltype,mat) = 0.15_rp     ! [KmCaMK] = mM
      vconc_initial(12,celltype,mat) = 0.0015_rp   ! [KmCaM] = mM
      vconc_initial(13,celltype,mat) = 0.0_rp      ! [Jrel_np] = mM/ms
      vconc_initial(14,celltype,mat) = 0.0_rp      ! [Jrel_p] = mM/ms
   else if (celltype == EXM_CELLTYPE_ENDO) then 
      vconc_initial(1, celltype,mat) = 12.1025_rp   ! [Nai] = mM
      vconc_initial(2, celltype,mat) = 12.1029_rp   ! [Nass] = mM
      vconc_initial(3, celltype,mat) = 142.3002_rp  ! [Ki] = mM
      vconc_initial(4, celltype,mat) = 142.3002_rp  ! [Kss] = mM
      vconc_initial(5, celltype,mat) = 7.453481e-05_rp ! [Cai] = mM
      vconc_initial(6, celltype,mat) = 7.0305e-5_rp ! [Cass] = mM
      vconc_initial(7, celltype,mat) = 1.5211_rp    ! [Cansr] = mM
      vconc_initial(8, celltype,mat) = 1.5214_rp    ! [Cajsr] = mM
      vconc_initial(9, celltype,mat) = 24.0_rp      ! [Cli] = mM
      vconc_initial(10,celltype,mat) = 0.0111_rp   ! [CaMKt] = mM
      vconc_initial(11,celltype,mat) = 0.15_rp     ! [KmCaMK] = mM
      vconc_initial(12,celltype,mat) = 0.0015_rp   ! [KmCaM] = mM
      vconc_initial(13,celltype,mat) = 0.0_rp      ! [Jrel_np] = mM/ms
      vconc_initial(14,celltype,mat) = 0.0_rp      ! [Jrel_p] = mM/ms
   end if

   ! Gating state variables
   if (celltype == EXM_CELLTYPE_EPI) then 
      vauxi_exm_initial(1, celltype,mat) = 0.00074303_rp   ! m
      vauxi_exm_initial(2, celltype,mat) = 0.8360_rp      ! h
      vauxi_exm_initial(3, celltype,mat) = 0.8359_rp      ! j
      vauxi_exm_initial(4, celltype,mat) = 0.6828_rp      ! hp
      vauxi_exm_initial(5, celltype,mat) = 0.8357_rp      ! jp
      vauxi_exm_initial(6, celltype,mat) = 0.00015166_rp      ! mL
      vauxi_exm_initial(7, celltype,mat) = 0.5401_rp      ! hL
      vauxi_exm_initial(8, celltype,mat) = 0.3034_rp      ! hLp
      vauxi_exm_initial(9, celltype,mat) = 0.00092716_rp      ! a
      vauxi_exm_initial(10,celltype,mat) = 0.9996_rp        ! i_F
      vauxi_exm_initial(11,celltype,mat) = 0.9996_rp        ! iS
      vauxi_exm_initial(12,celltype,mat) = 0.0004724_rp     ! ap
      vauxi_exm_initial(13,celltype,mat) = 0.9996_rp        ! iFp
      vauxi_exm_initial(14,celltype,mat) = 0.9996_rp        ! iSp
      vauxi_exm_initial(15,celltype,mat) = 0.0_rp                ! d
      vauxi_exm_initial(16,celltype,mat) = 1.0_rp                   ! ff
      vauxi_exm_initial(17,celltype,mat) = 0.9485_rp        ! fs
      vauxi_exm_initial(18,celltype,mat) = 1.0_rp                   ! fcaf
      vauxi_exm_initial(19,celltype,mat) = 0.9999_rp        ! fcas
      vauxi_exm_initial(20,celltype,mat) = 1.0_rp                   ! jca
      vauxi_exm_initial(21,celltype,mat) = 1.0_rp                   ! ffp
      vauxi_exm_initial(22,celltype,mat) = 1.0_rp                   ! fcafp
      vauxi_exm_initial(23,celltype,mat) = 0.00030853_rp    ! nca_ss
      vauxi_exm_initial(24,celltype,mat) = 0.00053006_rp    ! nca_i
      vauxi_exm_initial(25,celltype,mat) = 0.00067941_rp   ! C1
      vauxi_exm_initial(26,celltype,mat) = 0.00082869_rp    ! C2
      vauxi_exm_initial(27,celltype,mat) = 0.9982_rp                   ! C3
      vauxi_exm_initial(28,celltype,mat) = 0.00027561_rp     ! O_IKr
      vauxi_exm_initial(29,celltype,mat) = 9.5416e-06_rp    ! I_IKr
      vauxi_exm_initial(30,celltype,mat) = 0.2309_rp        ! xs1
      vauxi_exm_initial(31,celltype,mat) = 0.00016975_rp     ! xs2
   else if (celltype == EXM_CELLTYPE_MID) then
      vauxi_exm_initial(1, celltype,mat) = 0.00073818_rp   ! m
      vauxi_exm_initial(2, celltype,mat) = 0.8365_rp      ! h
      vauxi_exm_initial(3, celltype,mat) = 0.8363_rp      ! j
      vauxi_exm_initial(4, celltype,mat) = 0.6838_rp      ! hp
      vauxi_exm_initial(5, celltype,mat) = 0.8358_rp      ! jp
      vauxi_exm_initial(6, celltype,mat) = 0.00015079_rp      ! mL
      vauxi_exm_initial(7, celltype,mat) = 0.5327_rp      ! hL
      vauxi_exm_initial(8, celltype,mat) = 0.2834_rp      ! hLp
      vauxi_exm_initial(9, celltype,mat) = 0.00092527_rp      ! a
      vauxi_exm_initial(10,celltype,mat) = 0.9996_rp        ! i_F
      vauxi_exm_initial(11,celltype,mat) = 0.5671_rp        ! iS
      vauxi_exm_initial(12,celltype,mat) = 0.00047143_rp     ! ap
      vauxi_exm_initial(13,celltype,mat) = 0.9996_rp        ! iFp
      vauxi_exm_initial(14,celltype,mat) = 0.6261_rp        ! iSp
      vauxi_exm_initial(15,celltype,mat) = 0.0_rp                ! d
      vauxi_exm_initial(16,celltype,mat) = 1.0_rp                   ! ff
      vauxi_exm_initial(17,celltype,mat) = 0.92_rp        ! fs
      vauxi_exm_initial(18,celltype,mat) = 1.0_rp                   ! fcaf
      vauxi_exm_initial(19,celltype,mat) = 0.9998_rp        ! fcas
      vauxi_exm_initial(20,celltype,mat) = 1.0_rp                   ! jca
      vauxi_exm_initial(21,celltype,mat) = 1.0_rp                   ! ffp
      vauxi_exm_initial(22,celltype,mat) = 1.0_rp                   ! fcafp
      vauxi_exm_initial(23,celltype,mat) = 0.00051399_rp    ! nca_ss
      vauxi_exm_initial(24,celltype,mat) = 0.0012_rp    ! nca_i
      vauxi_exm_initial(25,celltype,mat) = 0.00069560_rp   ! C1
      vauxi_exm_initial(26,celltype,mat) = 0.00082672_rp    ! C2
      vauxi_exm_initial(27,celltype,mat) = 0.9979_rp                   ! C3
      vauxi_exm_initial(28,celltype,mat) = 0.00054206_rp     ! O_IKr
      vauxi_exm_initial(29,celltype,mat) = 1.8784e-05_rp    ! I_IKr
      vauxi_exm_initial(30,celltype,mat) = 0.2653_rp        ! xs1
      vauxi_exm_initial(31,celltype,mat) = 0.00016921_rp     ! xs2
   else if(celltype == EXM_CELLTYPE_ENDO) then 
      vauxi_exm_initial(1, celltype,mat) = 8.0572e-4_rp   ! m
      vauxi_exm_initial(2, celltype,mat) = 0.8286_rp      ! h
      vauxi_exm_initial(3, celltype,mat) = 0.8284_rp      ! j
      vauxi_exm_initial(4, celltype,mat) = 0.6707_rp      ! hp
      vauxi_exm_initial(5, celltype,mat) = 0.8281_rp      ! jp
      vauxi_exm_initial(6, celltype,mat) = 1.629e-4_rp      ! mL
      vauxi_exm_initial(7, celltype,mat) = 0.5255_rp      ! hL
      vauxi_exm_initial(8, celltype,mat) = 0.2872_rp      ! hLp
      vauxi_exm_initial(9, celltype,mat) = 9.5098e-4_rp      ! a
      vauxi_exm_initial(10,celltype,mat) = 0.9996_rp        ! i_F
      vauxi_exm_initial(11,celltype,mat) = 0.5936_rp        ! iS
      vauxi_exm_initial(12,celltype,mat) = 4.8454e-4_rp     ! ap
      vauxi_exm_initial(13,celltype,mat) = 0.9996_rp        ! iFp
      vauxi_exm_initial(14,celltype,mat) = 0.6538_rp        ! iSp
      vauxi_exm_initial(15,celltype,mat) = 8.1084e-9_rp                ! d
      vauxi_exm_initial(16,celltype,mat) = 1.0_rp                   ! ff
      vauxi_exm_initial(17,celltype,mat) = 0.939_rp        ! fs
      vauxi_exm_initial(18,celltype,mat) = 1.0_rp                   ! fcaf
      vauxi_exm_initial(19,celltype,mat) = 0.9999_rp        ! fcas
      vauxi_exm_initial(20,celltype,mat) = 1.0_rp                   ! jca
      vauxi_exm_initial(21,celltype,mat) = 1.0_rp                   ! ffp
      vauxi_exm_initial(22,celltype,mat) = 1.0_rp                   ! fcafp
      vauxi_exm_initial(23,celltype,mat) = 6.6462e-4_rp    ! nca_ss
      vauxi_exm_initial(24,celltype,mat) = 0.0012_rp    ! nca_i
      vauxi_exm_initial(25,celltype,mat) = 7.0344e-4_rp   ! C1
      vauxi_exm_initial(26,celltype,mat) = 8.5109e-4_rp   ! C2
      vauxi_exm_initial(27,celltype,mat) = 0.9981_rp                   ! C3
      vauxi_exm_initial(28,celltype,mat) = 3.7585e-4_rp     ! O_IKr
      vauxi_exm_initial(29,celltype,mat) = 1.3289e-5_rp    ! I_IKr
      vauxi_exm_initial(30,celltype,mat) = 0.248_rp        ! xs1
      vauxi_exm_initial(31,celltype,mat) = 1.7707e-4_rp     ! xs2
   end if

   ! Initialise single cell variables
   viclo_exm(:,1:2) = 0.0_rp
   vcolo_exm(:,1) = vconc_initial(:, celltype, mat)
   vcolo_exm(:,2) = vconc_initial(:, celltype, mat)
   elmlo_exm(1:2) = vminimate_exm(celltype,mat)
   vaulo_exm(:,1) = vauxi_exm_initial(:, celltype, mat)
   vaulo_exm(:,2) = vauxi_exm_initial(:, celltype, mat)
end subroutine exm_init_voltages




subroutine exm_ocetor(mat, ituss_exm, success_status, torord_stats)
   implicit none
   integer(ip), intent(in)     :: mat, ituss_exm
   integer(ip), intent(out)    :: success_status
   real(rp), intent(inout)    :: torord_stats(3_ip)  !saves ToR-ORd stats: number of beats, tolerance
   real(rp), dimension(1:2)   :: dummy

   call exm_ocetor_general(mat, ituss_exm, success_status, torord_stats, EXM_OCETOR_PURE, dummy)
end subroutine exm_ocetor

subroutine exm_ocetla(mat, ituss_exm, land_variables, success_status, torord_stats)
   implicit none
   integer(ip), intent(in)     :: mat, ituss_exm
   integer(ip), intent(out)    :: success_status
   real(rp), intent(inout)    :: torord_stats(3_ip)  !saves ToR-ORd stats: number of beats, tolerance
   real(rp), dimension(1:7), intent(inout)  :: land_variables

   call exm_ocetor_general(mat, ituss_exm, success_status, torord_stats, EXM_OCETOR_LAND, land_variables)
end subroutine exm_ocetla

!For EXM_OCETOR_PURE extra_outputs - not used
subroutine exm_ocetor_general(mat, ituss_exm, success_status, torord_stats, model_type, extra_outputs)
   use mod_exm_drugs, only : exm_drugs_ondrugs, exm_drugs_get_drugdmate

   ! definition of variables
   implicit none
   integer(ip), intent(in)     :: mat, ituss_exm, model_type
   integer(ip), intent(out)    :: success_status
   real(rp), dimension(:), intent(inout)  :: extra_outputs
   real(rp), intent(inout)    :: torord_stats(3_ip)  !saves ToR-ORd stats: number of beats, tolerance

   integer(ip) :: j, i, itim, nsamples_per_beat

   !arrays to store voltages of the last 2 beats
   !voltages_last will always point ot the last elemnt of the array,
   !when reached the end of the array, the voltages_last restarts from the beginning
   real(rp), dimension(:), pointer :: ss_voltages_beat_current, ss_voltages_beat_previous

   !number of voltages actually saved in the array, needed to verify that there was no error
   integer(ip) :: ss_voltage_current_id
   real(rp) :: ss_voltages_epsilon, ss_calcium_epsilon ! threshold squared to declare steady state when comparing two beats
   real(rp) :: ss_rmse, ss_min_voltage_threshold, ss_min_voltage
   integer(ip) :: ss_rmse_counter, ss_rmse_counter_max !the rmse has to maintain low value for ss_rmse_counter_max number of consecutive beats

   integer(ip)                       :: post_frequency, post_file_handle, dump_file_handle
   logical                           :: save_variables, flag_land, flag_3D, flag_drug,flag_border
   real(rp)  :: aux1,dur,dtimon,h,lambda,r_s,T_ref,Tactf
   real(rp)  :: statvar(7,2),atim,batim,drugd(24_ip),gkr_scaling

   character(len=50)         :: number_storage

   !to save calcium and voltage and other variables
   save_variables = ( kfl_save_convergence_cellmodel == 1_ip )
   post_frequency = 200_ip

   nullify (ss_voltages_beat_current)
   nullify (ss_voltages_beat_previous)

   !sanity check
   SELECT CASE (model_type)
   CASE (EXM_OCETOR_PURE)
      ! nothing
   CASE (EXM_OCETOR_LAND)
      if ( size(extra_outputs,1,kind=ip)<7 ) then
         call runend("mod_exm_ocetor: Supplied array for Land variables needs to have at least 7 elements. The array has length = "//trim(intost(size(extra_outputs,1,kind=ip))))
      end if
   CASE DEFAULT
      call runend("mod_exm_ocetor: Unknown cell model: "//trim(intost(model_type)))
   END SELECT

   !write file header
   if (save_variables) then

      post_file_handle = 1000000 + 100*mat + ituss_exm

      if (IPARALL) post_file_handle = post_file_handle + PAR_MY_CODE_RANK*10000

      OPEN (unit=post_file_handle, file='torord_ode.m'//trim(intost(mat))//"c"//trim(intost(ituss_exm))//".csv", form="FORMATTED", status="REPLACE")


      SELECT CASE (model_type)
      CASE (EXM_OCETOR_PURE)
         WRITE (post_file_handle,*) "Global_Time Beat_Time Calcium Voltage"
      CASE (EXM_OCETOR_LAND)
         WRITE (post_file_handle,*) "Global_Time Beat_Time Calcium Voltage T_active"
      CASE DEFAULT
         call runend("mod_exm_ocetor: Unknown cell model: "//trim(intost(model_type)))
      END SELECT

   end if

   !if (kfl_borde_exm(mat) == 1) then
  !    flag_border = .true.
  !    gkr_scaling = border_gkrsc_exm(mat)
  ! else
   flag_border = .false.
   gkr_scaling = 1.0_rp
  ! end if

   if ( exm_drugs_ondrugs(mat) ) then
      flag_drug = .true.
      call exm_drugs_get_drugdmate(mat, drugd, size(drugd, 1_ip, kind=ip)/2_ip, 2_ip)
   else
      flag_drug = .false.
      drugd = 1.0_rp
   end if

   if (model_type == EXM_OCETOR_LAND) then
      flag_land = .true.
   else
      flag_land = .false.
   end if

   flag_3D = .false.

   dtimon = 0.02_rp ! [ms]


   !----------------------------------------------
   !
   ! begin: Initialize variables for steady state detection
   !
   !----------------------------------------------
   !allocate the arrays to store two full beats
   nsamples_per_beat = nint(moneclmate_exm(2, mat)/dtimon)

   call memory_alloca(mem_modul(1:2, modul), 'ss_voltages_beat_current', 'exm_ocetor', ss_voltages_beat_current, nsamples_per_beat)
   call memory_alloca(mem_modul(1:2, modul), 'ss_voltages_beat_previous', 'exm_ocetor', ss_voltages_beat_previous, nsamples_per_beat)

   ss_voltages_beat_current = 0.0_rp
   ss_voltages_beat_previous = 0.0_rp

   ss_voltages_epsilon = 1.0_rp !this is maximum allowed  mean voltage difference (squared)
   ss_calcium_epsilon = 1.0e-7_rp

   if (steady_tol_cellm(mat) >= 0) then
      ss_voltages_epsilon = steady_tol_cellm(mat)
      ss_calcium_epsilon  = steady_tol_cellm(mat)
   end if

   ss_rmse_counter = 0_ip
   ss_rmse_counter_max = 3_ip
   ss_min_voltage_threshold = -75.0_rp !the min volytage within the beat has to be lower than this, to consider that we are in steady state
   !----------------------------------------------
   !
   ! end: Initialize variables for steady state detection
   !
   !----------------------------------------------

   success_status = 1_ip !return error by default
   itim = 0_ip
   batim = 0.0_rp
   j = 0_ip
    ! moneclmate_exm(1, mat) - number of beats
   ! moneclmate_exm(2, mat) - cycle length
   do while (batim < (moneclmate_exm(1, mat)*moneclmate_exm(2, mat)))
      aux1 = 0.0_rp
      atim = dtimon*aux1

      !for steady state verification
      ss_voltage_current_id = 1_ip
      do while (atim < moneclmate_exm(2, mat))
         !if (atim < moneclmate_exm(2) )then
         aux1 = aux1 + 1.0_rp
         atim = dtimon*aux1

         dur = 1.0_rp ! duration time (ms)
         if (atim < dur) then
            ! External stimulus current I_stim [pA/pF]
            viclo_exm(20, 1) = -53.0_rp !-80.0_rp
            !viclo_exm(23,3) = -80.0_rp     ! Value of I_stim current [pA/pF]
         else
            viclo_exm(20, 1) = 0.0_rp
            !viclo_exm(23,3) = 0.0_rp
         end if

         if (flag_land) then
            call exm_torord_model(ituss_exm,elmlo_exm(1),vcolo_exm,vaulo_exm, &
             & viclo_exm,dtimon,1.0_rp,elmlo_exm(2),flag_land,flag_3D, &
             & flag_border,flag_drug,gkr_scaling,drugd,statvar)
         else
            call exm_torord_model(ituss_exm,elmlo_exm(1),vcolo_exm,vaulo_exm, &
             & viclo_exm,dtimon,1.0_rp,elmlo_exm(2),flag_land,flag_3D,&
             & flag_border,flag_drug,gkr_scaling,drugd)
         end if

         if (flag_land) then
            statvar(:,1) = statvar(:,2)
            lambda = 1.0_rp
            h = max(0.0_rp,1.0_rp+2.3_rp*(lambda+min(0.87_rp,lambda)-1.87_rp))
            T_ref = 120.0_rp
            r_s = 0.25_rp
            Tactf = h*(T_ref/r_s)*(statvar(1,2)*(statvar(5,2)+1.0_rp)+statvar(2,2)*statvar(6,2))
         end if

         itim = itim + 1

         if (save_variables) then
            if (j >= post_frequency) then  !it was 10 for saving


               SELECT CASE (model_type)
               CASE (EXM_OCETOR_PURE)
                  write (post_file_handle, "(4(e16.8E3, ' '))") itim*dtimon, atim, vcolo_exm(5, 2), elmlo_exm(2)
               CASE (EXM_OCETOR_LAND)
                  write (post_file_handle, "(5(e16.8E3, ' '))") itim*dtimon, atim, vcolo_exm(5, 2), elmlo_exm(2),Tactf
               CASE DEFAULT
                  call runend("mod_exm_ocetor: Unknown cell model: "//trim(intost(model_type)))
               END SELECT

               j = 0
            end if
            j = j + 1
         end if 
         elmlo_exm(1) = elmlo_exm(2)
         vcolo_exm(:,1) = vcolo_exm(:,2)
         vaulo_exm(:,1) = vaulo_exm(:,2)
         !save_variables

         select case (kfl_steadystate_variable(mat))
         case (EXM_CELL_STEADY_VOLTAGE)
            ss_voltages_beat_current(ss_voltage_current_id) = elmlo_exm(2)
         case (EXM_CELL_STEADY_CALCIUM)
            ss_voltages_beat_current(ss_voltage_current_id) = vcolo_exm(5, 2)
         case default
            call runend("EXM_OCETOR: Unknown name of the variable to determine the cell convergence.")
         end select

         ss_voltage_current_id = ss_voltage_current_id + 1_ip

      end do

      batim = itim*dtimon

      !--------------------------------------------
      !
      !  begin: check if we reached steady state. If so, break
      !
      !----------------------------------------
      if (ss_voltage_current_id - 1 .NE. nsamples_per_beat) then
         call runend("EXM_OCETOR: Something went wrong in the loop of one beat. The number of timesteps executed is incorrect.")
      end if

      select case (kfl_steadystate_variable(mat))
      case (EXM_CELL_STEADY_VOLTAGE)
         ss_min_voltage = minval(ss_voltages_beat_current)
         !print *, 'Beat ',itim/nsamples_per_beat,': min voltage: ',ss_min_voltage
         if (ss_min_voltage .LE. ss_min_voltage_threshold) then

            ss_rmse = sqrt(sum((ss_voltages_beat_current - ss_voltages_beat_previous)**2/nsamples_per_beat))
            !print *, 'Beat ',itim/nsamples_per_beat,': rmse between beats: ',ss_rmse
            if (ss_rmse < ss_voltages_epsilon) then
               ss_rmse_counter = ss_rmse_counter + 1
            else
               ss_rmse_counter = 0_ip
            end if
         end if
      case (EXM_CELL_STEADY_CALCIUM)
         ss_rmse = sqrt(sum((ss_voltages_beat_current - ss_voltages_beat_previous)**2/nsamples_per_beat))
         !print *, 'Beat ',itim/nsamples_per_beat,': rmse between beats: ',ss_rmse
         if (ss_rmse < ss_calcium_epsilon) then
            ss_rmse_counter = ss_rmse_counter + 1
         else
            ss_rmse_counter = 0_ip
         end if
      case default
         call runend("EXM_OCETOR: Unknown name of the variable to determine the cell convergence.")
      end select

      if (ss_rmse_counter .GE. ss_rmse_counter_max) then
         success_status = 0_ip !success
         exit
      end if

      !otherwise copy the beat and go to the beginning
      ss_voltages_beat_previous = ss_voltages_beat_current
      !--------------------------------------------
      !
      !  end: check if we reached steady state
      !
      !----------------------------------------

   end do

   !------------------------------------------------
   !
   ! Save stats
   !
   !----------------------------------------------
   torord_stats(1_ip) = itim/nsamples_per_beat
   select case (kfl_steadystate_variable(mat))
   case (EXM_CELL_STEADY_VOLTAGE)
      torord_stats(2_ip) = ss_voltages_epsilon
   case (EXM_CELL_STEADY_CALCIUM)
      torord_stats(2_ip) = ss_calcium_epsilon
   case default
      call runend("EXM_OCETOR: Unknown name of the variable to determine the cell convergence.")
   end select
   torord_stats(3_ip) = ss_rmse

   if (save_variables) close (post_file_handle)

   call memory_deallo(mem_modul(1:2, modul), 'ss_voltages_beat_current', 'exm_ocetor', ss_voltages_beat_current)
   call memory_deallo(mem_modul(1:2, modul), 'ss_voltages_beat_previous', 'exm_ocetor', ss_voltages_beat_previous)


   if (model_type == EXM_OCETOR_LAND) then
      !-----------------------------------------------
      !
      !   Land specific, start
      !
      !-----------------------------------------------------
      extra_outputs(:) = statvar(:,2)
      !-----------------------------------------------
      !
      !   Land specific, end
      !
      !-----------------------------------------------------
   end if

   !--------------------------------------------------------
   !
   ! This fragment dumps variables for predefining the conditions
   !
   !--------------------------------------------------------
   if ( kfl_save_init_cellmodel==1_ip ) then
      dump_file_handle = 2000000 + 100*mat + ituss_exm

      if (IPARALL) dump_file_handle = dump_file_handle + PAR_MY_CODE_RANK*10000

      OPEN (unit=dump_file_handle, file='torord_dump.m'//trim(intost(mat))//"c"//trim(intost(ituss_exm))//".csv", form="FORMATTED", status="REPLACE")

      WRITE(number_storage,*) elmlo_exm(2)
      WRITE(dump_file_handle,*) "elmlo_exm(1:2) = "//trim(number_storage)//"_rp"

      WRITE(dump_file_handle,*) "viclo_exm(1:26, 1) = 0.0_rp"
      do i=1,size(vcolo_exm, 1_ip, KIND=ip)
         WRITE(number_storage,*) vcolo_exm(i, 1)
         WRITE(dump_file_handle,*) "vcolo_exm("//trim(intost(i))//", 1:3) = "//trim(number_storage)//"_rp"
      end do
      do i=1,size(vaulo_exm, 1_ip, KIND=ip)
         WRITE(number_storage,*) vaulo_exm(i, 1)
         WRITE(dump_file_handle,*) "vaulo_exm("//trim(intost(i))//", 1:2) = "//trim(number_storage)//"_rp"
      end do
      CLOSE(dump_file_handle)
   end if



end subroutine exm_ocetor_general


end module mod_exm_onetor
