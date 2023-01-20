!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    mod_sld_rbo.f90
!> @author  Gerard Guillamet (gerard.guillamet@bsc.es)
!> @date    January, 2019
!>          - Subroutine written
!> @brief   Rigid body module
!>
!> @details Rigid body module
!>
!>  New flags for PDN-contact
!>             MATIAS PDN contact:
!>                UNILA: kfl_conta_sld=1, kfl_rigid_sld = 0
!>                BILAT: kfl_conta_sld=2, kfl_rigid_sld = 0
!>             NEW PDN contact:
!>                RBODE: kfl_conta_sld=3, kfl_rigid_sld = 1 o 0
!>
!>             NEW PDN contact for SHERLOC
!>                RBODE: kfl_conta_sld=4
!>              
!> @}
!------------------------------------------------------------------------

module mod_sld_rbo

  use def_kintyp,           only : ip, rp
  use def_master,           only : ittim
  use def_master,           only : INOTMASTER, INOTSLAVE, ITER_K, TIME_N
  use def_domain,           only : ndime, npoin
  use def_solidz,           only : kfl_rigid_sld, ncomp_sld
  use def_solidz,           only : nderi_sld
  use mod_array_operations, only : array_operations_residual_norm

  implicit none
  !
  ! Linear motion
  !
  real(rp)            :: rbdil_sld(3,3)
  real(rp)            :: rbvel_sld(3,3)
  real(rp)            :: rbacl_sld(3,3)
  real(rp)            :: rbbve_sld(3,3)
  real(rp)            :: uruk4_sld(6)
  
  public   :: sld_rbo_solrbo
  public   :: sld_rbo_updunk
  public   :: sld_rbo_updbcs 
  public   :: sld_rbo_doiter
  public   :: sld_rbo_inipro
  public   :: sld_rbo_iniunk 
  public   :: sld_rbo_solmem 
  public   :: sld_rbo_cvgunk 
  public   :: sld_rbo_restar

contains

  !-----------------------------------------------------------------------
  !> @author  Gerard Guillamet
  !> @date    March, 2019
  !> @brief   Solve RBO equations
  !> @details Solve RBO equations
  !-----------------------------------------------------------------------

  subroutine sld_rbo_solrbo()

    use def_master,               only : ITER_K, ITER_AUX, TIME_N
    use def_master,               only : displ
    use def_master,               only : dtime
    use def_domain,               only : ndime, npoin
    use mod_communications,       only : PAR_MAX
    use def_solidz,               only : kfl_rbfix_sld, kfl_timei_sld
    use def_solidz,               only : veloc_sld, accel_sld
    use def_solidz,               only : nderi_sld
    use def_solidz,               only : SLD_STATIC_PROBLEM, SLD_DYNAMIC_PROBLEM
    use mod_sld_solution_methods, only : sld_rk4, sld_rk4_derivs

    implicit none

    integer(ip)          :: ipoin, idime, icomp
    real(rp)             :: dU(nderi_sld)

    !-------------------------------------------------------------------
    !
    ! Solve RBO equations
    !
    !-------------------------------------------------------------------

    if( INOTMASTER ) then
       if( kfl_timei_sld == SLD_DYNAMIC_PROBLEM ) then
          !
          ! Get displacement and velocity
          !
          call sld_rk4(1_ip, dtime,  nderi_sld, uruk4_sld(:))
          rbdil_sld(1:ndime,ITER_K) = uruk4_sld(1:ndime)
          rbvel_sld(1:ndime,ITER_K) = uruk4_sld(ndime+1:ndime*2)
          !
          ! Get derivatives: velocity and acceleration
          !
          dU = 0.0_rp
          call sld_rk4_derivs(1_ip, nderi_sld, uruk4_sld(:), dU(:))
          rbacl_sld(1:ndime,ITER_K) = dU(ndime+1:ndime*2)
       end if
    end if
    !
    ! Parall service
    !
    ! Motion
    do icomp = 1,ncomp_sld
       call PAR_MAX(3_ip,rbdil_sld(:,icomp),'IN MY CODE')
       call PAR_MAX(3_ip,rbvel_sld(:,icomp),'IN MY CODE')
       call PAR_MAX(3_ip,rbacl_sld(:,icomp),'IN MY CODE')   
    end do
    call PAR_MAX(6_ip,uruk4_sld(:),'IN MY CODE')
    
    !-------------------------------------------------------------------
    !
    ! Move mesh according to the boundary conditions
    !
    !-------------------------------------------------------------------

    if( INOTMASTER ) then

       if(      kfl_timei_sld == SLD_STATIC_PROBLEM ) then

          do ipoin = 1,npoin
             do idime = 1,ndime
                if( kfl_rbfix_sld(idime) == 1 ) then
                   displ(idime,ipoin,ITER_K) = rbbve_sld(idime,ITER_K)
                end if
                veloc_sld(idime,ipoin,ITER_K) = 0.0_rp
                accel_sld(idime,ipoin,ITER_K) = 0.0_rp
             end do
          end do

       else if( kfl_timei_sld == SLD_DYNAMIC_PROBLEM ) then

          do ipoin = 1,npoin
             do idime = 1,ndime
                if( kfl_rbfix_sld(idime) == 1 ) then
                   displ(    idime,ipoin,ITER_K) = 0.0_rp
                   veloc_sld(idime,ipoin,ITER_K) = 0.0_rp
                   accel_sld(idime,ipoin,ITER_K) = 0.0_rp
                else
                   displ(    idime,ipoin,ITER_K) = rbdil_sld(idime,ITER_K)
                   veloc_sld(idime,ipoin,ITER_K) = rbvel_sld(idime,ITER_K)
                   accel_sld(idime,ipoin,ITER_K) = rbacl_sld(idime,ITER_K)
                end if
             end do
          end do

       end if

    end if

  end subroutine sld_rbo_solrbo

  !-----------------------------------------------------------------------
  !> @author  Gerard Guillamet
  !> @date    March, 2019
  !> @brief   Updates on rigid body unknowns
  !> @details
  !>          \verbatim
  !>          TIME STEP AND ITERATION LABELS (Defined in def_master)
  !>            ITER_K       = 1 ..... Current iteration at time step N
  !>            ITER_AUX     = 2 ..... Used for coupling iterations
  !>            TIME_N       = 3 ..... Time step (Converged)
  !>
  !>          PRE-PROCESS
  !>            sld_iniunk (itask=6) .......... if Restart
  !>
  !>          TIME LOOP
  !>            do time
  !>              sld_begste (itask=1) ........ (:,ITER_AUX) <= (:,TIME_N)
  !>              do outer
  !>                sld_begite (itask=2) ...... (:,ITER_K)   <= (:,ITER_AUX)
  !>                do inner
  !>                  sld_endite (itask=3) .... (:,ITER_K)   <= UNKNO
  !>                end do
  !>                sld_endite (itask=4p) ..... (:,ITER_AUX) <= (:,ITER_K)
  !>              end do
  !>              sld_endste (itask=5) ........ (:,TIME_N)   <= (:,ITER_K)
  !>            end do
  !>          \endverbatim
  !> @}
  !-----------------------------------------------------------------------

  subroutine sld_rbo_updunk(itask)

    use def_master, only : ITER_K, ITER_AUX, TIME_N
    use def_master, only : ITASK_INIUNK
    use def_master, only : ITASK_BEGSTE, ITASK_ENDSTE
    use def_master, only : ITASK_BEGITE, ITASK_ENDINN, ITASK_ENDITE
    use def_master, only : displ
    use def_domain, only : npoin
    use def_solidz, only : veloc_sld, accel_sld
    
    implicit none

    integer(ip), intent(in) :: itask !< What to do
    integer(ip)             :: icomp, ipoin

    select case (itask)

    case( ITASK_INIUNK )

       !-------------------------------------------------------------------
       !
       ! (:,*) <= (:,1) Initial solution
       !
       !-------------------------------------------------------------------
       !
       ! Rigid body (CoG)
       !
       ! Motion
       do icomp = 1,ncomp_sld
          rbdil_sld(1:ndime,icomp) = rbdil_sld(1:ndime,ITER_K)
          rbvel_sld(1:ndime,icomp) = rbvel_sld(1:ndime,ITER_K)
          rbacl_sld(1:ndime,icomp) = rbacl_sld(1:ndime,ITER_K)
       end do
       !
       ! Mesh
       !
       if( INOTMASTER ) then
          do icomp = 1,size(displ,DIM=3,KIND=ip)
             do ipoin = 1,npoin
                displ(:,ipoin,icomp)     = displ(:,ipoin,ITER_K)
                veloc_sld(:,ipoin,icomp) = veloc_sld(:,ipoin,ITER_K)
                accel_sld(:,ipoin,icomp) = accel_sld(:,ipoin,ITER_K)
             end do
          end do
       end if

    case( ITASK_BEGSTE )

       !-------------------------------------------------------------------
       !
       ! Initial guess for outer iterations (begin time step)
       !
       !-------------------------------------------------------------------
       !
       ! Rigid body (CoG)
       !
       ! Motion       
       rbdil_sld(1:ndime,ITER_AUX) = rbdil_sld(1:ndime,ITER_K)
       rbvel_sld(1:ndime,ITER_AUX) = rbvel_sld(1:ndime,ITER_K)
       rbacl_sld(1:ndime,ITER_AUX) = rbacl_sld(1:ndime,ITER_K)
       !
       ! Mesh
       !
       if( INOTMASTER ) then
          do ipoin = 1,npoin
             displ(    1:ndime,ipoin,ITER_AUX) = displ(    1:ndime,ipoin,ITER_K)
             veloc_sld(1:ndime,ipoin,ITER_AUX) = veloc_sld(1:ndime,ipoin,ITER_K)
             accel_sld(1:ndime,ipoin,ITER_AUX) = accel_sld(1:ndime,ipoin,ITER_K)
          end do
       end if
          
    case ( ITASK_BEGITE )

       !------------------------------------------------------------------
       !
       !  Initial guess for inner iterations (begin iterations)
       !
       !------------------------------------------------------------------
       !
       ! Rigid body (CoG)
       !
       ! Motion
       rbdil_sld(1:ndime,ITER_K) = rbdil_sld(1:ndime,ITER_AUX)
       rbvel_sld(1:ndime,ITER_K) = rbvel_sld(1:ndime,ITER_AUX)
       rbacl_sld(1:ndime,ITER_K) = rbacl_sld(1:ndime,ITER_AUX)
       !
       ! Mesh
       !
       if( INOTMASTER ) then
          do ipoin = 1,npoin
             displ(    1:ndime,ipoin,ITER_K) = displ(    1:ndime,ipoin,ITER_AUX)
             veloc_sld(1:ndime,ipoin,ITER_K) = veloc_sld(1:ndime,ipoin,ITER_AUX)
             accel_sld(1:ndime,ipoin,ITER_K) = accel_sld(1:ndime,ipoin,ITER_AUX)
          end do
       end if

    case( ITASK_ENDITE )

       !------------------------------------------------------------------
       !
       !  Updates (end iterations converged)
       !
       !------------------------------------------------------------------
       !
       ! Rigid body (CoG)
       !
       ! Motion
       rbdil_sld(1:ndime,ITER_AUX) = rbdil_sld(1:ndime,ITER_K)
       rbvel_sld(1:ndime,ITER_AUX) = rbvel_sld(1:ndime,ITER_K)
       rbacl_sld(1:ndime,ITER_AUX) = rbacl_sld(1:ndime,ITER_K)
       !
       ! Mesh
       !
       if( INOTMASTER ) then
          do ipoin = 1,npoin
             displ(    1:ndime,ipoin,ITER_AUX) = displ(    1:ndime,ipoin,ITER_K)
             veloc_sld(1:ndime,ipoin,ITER_AUX) = veloc_sld(1:ndime,ipoin,ITER_K)
             accel_sld(1:ndime,ipoin,ITER_AUX) = accel_sld(1:ndime,ipoin,ITER_K)
          end do
       end if

    case( ITASK_ENDSTE )

       !-------------------------------------------------------------------
       !
       ! End of time step
       !
       !-------------------------------------------------------------------
       !
       ! Rigid body (CoG)
       !
       ! Motion
       rbdil_sld(1:ndime,TIME_N) = rbdil_sld(1:ndime,ITER_K)
       rbvel_sld(1:ndime,TIME_N) = rbvel_sld(1:ndime,ITER_K)
       rbacl_sld(1:ndime,TIME_N) = rbacl_sld(1:ndime,ITER_K)
       !
       ! Mesh
       !
       if( INOTMASTER ) then
          do ipoin = 1,npoin
             displ(    1:ndime,ipoin,TIME_N) = displ(    1:ndime,ipoin,ITER_K)
             veloc_sld(1:ndime,ipoin,TIME_N) = veloc_sld(1:ndime,ipoin,ITER_K)
             accel_sld(1:ndime,ipoin,TIME_N) = accel_sld(1:ndime,ipoin,ITER_K)
          end do
       end if

    end select

  end subroutine sld_rbo_updunk

  !-----------------------------------------------------------------------
  !> @author  Gerard Guillamet
  !> @date    March, 2019
  !> @brief   Initialize unknowns at the begining of the run
  !> @details Initialize unknowns at the begining of the run
  !-----------------------------------------------------------------------

  subroutine sld_rbo_iniunk()

    use def_master, only : ITER_K
    use def_domain, only : npoin
    use def_solidz, only : veloc_sld
    use def_solidz, only : kfl_invel_sld, invel_sld

    implicit none

    integer(ip)     :: ipoin

    !----------------------------------------------------------------------
    !
    ! Initialization
    !
    !----------------------------------------------------------------------

    if( kfl_invel_sld /= 0_ip ) then
       ! CoG
       uruk4_sld(ndime+1:ndime*2) = invel_sld(1:ndime)
       rbvel_sld(1:ndime,ITER_K)  = invel_sld(1:ndime)
       ! Mesh
       if( INOTMASTER ) then
          do ipoin = 1,npoin
             veloc_sld(1:ndime,ipoin,ITER_K) = invel_sld(1:ndime)
          end do
       end if
    end if

  end subroutine sld_rbo_iniunk

  !-----------------------------------------------------------------------
  !> @author  Gerard Guillamet
  !> @date    March, 2019
  !> @brief   Calculation of RB variables at the begining of the run
  !> @details Calculation of RB variables at the begining of the run
  !-----------------------------------------------------------------------

  subroutine sld_rbo_inipro(itask)

    use def_solidz, only : rbfor_sld

    implicit none

    integer(ip), intent(in) :: itask

    select case (itask)

    case( 0_ip )

       !-------------------------------------------------------------------
       !
       ! Initialization
       !
       !-------------------------------------------------------------------
       !
       ! Properties
       !
       !
       ! Linear motion
       !
       rbdil_sld(:,:) = 0.0_rp
       rbvel_sld(:,:) = 0.0_rp
       rbacl_sld(:,:) = 0.0_rp
       !
       ! External force
       !
       rbfor_sld(:)   = 0.0_rp
       !
       ! Prescribed displacement
       !
       rbbve_sld(:,:) = 0.0_rp

    end select

  end subroutine sld_rbo_inipro

  !-----------------------------------------------------------------------
  !> @author  Gerard Guillamet
  !> @date    April, 2019
  !> @brief   General memory for the RBO variables
  !> @details General memory for the RBO variables
  !-----------------------------------------------------------------------

  subroutine sld_rbo_solmem()

    use def_master, only : mem_modul, modul
    use def_domain, only : nelem
    use mod_memory, only : memory_alloca
    use def_solidz, only : celen_sld

    implicit none

    if( INOTMASTER ) then
       !
       ! Element characteristic length
       !
       call memory_alloca(mem_modul(1:2,modul),'CELEN_SLD','sld_rbo_memall',celen_sld,nelem)

    else
       !
       ! Element characteristic length
       !
       call memory_alloca(mem_modul(1:2,modul),'CELEN_SLD','sld_rbo_memall',celen_sld,1_ip)

    end if

  end subroutine sld_rbo_solmem

  !-----------------------------------------------------------------------
  !> @author  Gerard Guillamet
  !> @date    April, 2019
  !> @brief   Iterations loop for RBO
  !> @details
  !-----------------------------------------------------------------------

  subroutine sld_rbo_doiter()

    use def_kintyp,     only : ip
    use def_master,     only : ITASK_BEGITE, ITASK_ENDITE, ITASK_ENDINN
    use def_master,     only : modul
    use def_master,     only : itcou, itinn
    use mod_messages,   only : livinf
    use def_solidz,     only : kfl_stead_sld, kfl_goite_sld

    implicit none

    if( kfl_stead_sld == 0_ip ) then
       !
       ! Initializations
       !
       kfl_goite_sld = 1
       itinn(modul)  = 0
       if ( itcou == 1 ) call sld_tistep()
       call livinf(15_ip,' ',modul)
       !
       ! Updates
       !
       call sld_rbo_updbcs(ITASK_BEGITE)
       call sld_rbo_updunk(ITASK_BEGITE)  ! Update unknowns: (:,ITER_K)   <= (:,ITER_AUX)
       !
       ! Iterate
       !
       do while( kfl_goite_sld == 1_ip ) ! if not converged
          call sld_solite()
          kfl_goite_sld = 0_ip
          !
          call livinf(-3_ip,' ',1_ip)
          call livinf(56_ip,' ',modul)
          !
          ! Convergence and Timings
          !
          call sld_rbo_cvgunk(0_ip)          
       end do
       !
       ! End iteration
       !
       call livinf(16_ip,' ',itinn(modul))
       call sld_rbo_updunk(ITASK_ENDITE)  ! Update unknowns: (:,ITER_AUX) <= (:,ITER_K)

    end if

  end subroutine sld_rbo_doiter

  !------------------------------------------------------------------------
  !> @author  Gerard Guillamet
  !> @date    June, 2019
  !> @brief   Convergence checks for Solidz RBO
  !> @details
  !------------------------------------------------------------------------

  subroutine sld_rbo_cvgunk(itask)

    use def_kintyp,         only : ip, rp
    use def_master,         only : zeror, isnain, INOTSLAVE, ITER_K, TIME_N
    use def_master,         only : ITASK_ENDSTE
    use def_master,         only : displ, modul, kfl_rstar
    use def_master,         only : cpu_initi, momod
    use def_master,         only : itinn, itcou, ittim, cutim, dtime
    use def_domain,         only : ndime, npoin
    use mod_array_operations, only : array_operations_residual_norm
    use mod_outfor,         only : outfor
    use mod_iofile,         only : iofile_flush_unit
    use def_solidz,         only : kfl_timei_sld, kfl_stead_sld
    use def_solidz,         only : rbfor_sld, rbmas_sld
    use def_solidz,         only : sstol_sld
    use def_solidz,         only : SLD_STATIC_PROBLEM

    implicit none

    integer(ip), intent(in) :: itask !< What to do
    integer(ip), save       :: kpass = 0
    real(rp),    save       :: cpuit_sld=0.0_rp
    real(rp),    save       :: resti=1.0_rp
    real(rp)                :: time1

    select case ( itask )

    case ( 0_ip )

       !-------------------------------------------------------------------
       !
       ! Timings and output
       !
       !-------------------------------------------------------------------
       !
       ! Write convergence
       !
       if( INOTSLAVE ) then
          call cputim(time1)
          if( kpass == 0 .and. kfl_rstar /= 2 ) then
             write(momod(modul) % lun_conve,100)
          end if
          if( kpass == 1 ) then
             time1 = time1 - cpuit_sld
          else
             time1 = time1 - cpu_initi
          end if
          !
          ! Write results
          !
          write(momod(modul) % lun_conve,101) &
               !   1     2            3     4     5
               ittim,itcou,itinn(modul),cutim,dtime,&
               !                 6                   7                   8
               rbdil_sld(1,ITER_K),rbdil_sld(2,ITER_K),rbdil_sld(3,ITER_K),&
               !                 9                  10                  11
               rbvel_sld(1,ITER_K),rbvel_sld(2,ITER_K),rbvel_sld(3,ITER_K),&
               !                12                  13                  14
               rbacl_sld(1,ITER_K),rbacl_sld(2,ITER_K),rbacl_sld(3,ITER_K),&
               !         15           16           17
               rbfor_sld(1),rbfor_sld(2),rbfor_sld(3),&
               !                          18                            19                            20
               rbacl_sld(1,ITER_K)*rbmas_sld,rbacl_sld(2,ITER_K)*rbmas_sld,rbacl_sld(3,ITER_K)*rbmas_sld,&
               !                                    21                                      22
               0.5_rp*rbmas_sld*rbvel_sld(1,ITER_K)**2,0.5_rp*rbmas_sld*rbvel_sld(2,ITER_K)**2,&
               ! 23
               0.5_rp*rbmas_sld*rbvel_sld(3,ITER_K)**2,&
               !  24
               time1

          call cputim(cpuit_sld)
          call iofile_flush_unit(momod(modul) % lun_conve)

       end if

       kpass = 1_ip

    case ( ITASK_ENDSTE )

       !-------------------------------------------------------------------
       !
       ! Check residual of the time iterations
       ! Only for dynamic (transient) problems reaching a steady state.
       ! RESTI: Time residual measured only using displacement
       !
       ! ||u^{n+1}-u^n|/||u^{n+1}||
       !-------------------------------------------------------------------

       if( kfl_timei_sld /= SLD_STATIC_PROBLEM ) then
          resti = array_operations_residual_norm(2_ip,ndime,ndime,displ,displ,0_ip,2_ip*npoin*ndime,ndime)
       end if
       !
       ! Some check to stop the code if things are going bad
       !
       if( isnain(resti)     ) kfl_stead_sld = 1_ip    ! NaN or +/- Inf
       if( resti > 1.0e10_rp ) kfl_stead_sld = 1_ip    ! Huge
       !
       ! Steady-state has been reached
       if( resti <= sstol_sld ) then
          kfl_stead_sld = 1_ip
          call outfor(28_ip,momod(modul) % lun_outpu,' ')
       end if

    end select

    !----------------------------------------------------------------------
    !
    ! Formats
    !
    !----------------------------------------------------------------------

100 format('# --| ALYA RBO Convergence'  ,/,&
         & '# --| Columns displayed:' ,/,&
         & '# --|                                                                          ',/,&
         & '# --|  1. Time step             2. Global Iteration      3. Inner Iteration    ',/,&
         & '# --|  4. Current time          5. Time increment                              ',/,&
         & '# --|  6. Linear displ. UX      7. Linear displ. UY      8. Linear displ. UZ   ',/,&
         & '# --|  9. Linear veloc. VX     10. Linear veloc. VY     11. Linear veloc. VZ   ',/,&
         & '# --| 12. Linear accel. AX     13. Linear accel. AY     14. Linear accel. AZ   ',/,&
         & '# --| 15. Linear force. FX     16. Linear force. FY     17. Linear force. FZ   ',/,&
         & '# --| 18. Force from accel. FX 19. Force from accel. FY 20. Force from acc. FZ ',/,&
         & '# --| 21. Kinetic energy EX    22. Kinetic energy KY    23. Kinetic energy KZ  ',/,&
         & '# --| 24. Elapsed CPU time                                                     ',/,&
         & '# --|                                                                          ')

101 format(4x,i9,2x,i9,2x,i9,50(2x,e16.8e3))

  end subroutine sld_rbo_cvgunk

  !------------------------------------------------------------------------
  !> @author  Gerard Guillamet
  !> @date    August, 2020
  !> @brief   Update boundary conditions for the rigid body
  !> @details
  !------------------------------------------------------------------------

  subroutine sld_rbo_updbcs(itask)

    use def_kintyp,                  only : ip, rp
    use def_master,                  only : INOTEMPTY, zeror
    use def_master,                  only : FUNCTION_TIME
    use def_master,                  only : FUNCTION_MODULE
    use def_master,                  only : FUNCTION_SPACE_TIME
    use def_master,                  only : ITER_K, ITER_AUX, TIME_N
    use def_master,                  only : ITASK_TURNON, ITASK_BEGSTE
    use def_master,                  only : ITASK_BEGITE, ITASK_INNITE
    use def_domain,                  only : ndime, npoin
    use mod_ker_functions,           only : ker_functions
    use mod_ker_space_time_function, only : ker_space_time_function
    use mod_communications,          only : PAR_MAX
    use mod_communications,          only : PAR_MIN
    use def_solidz,                  only : bvess_sld
    use def_solidz,                  only : kfl_conbc_sld, kfl_funno_sld
    use def_solidz,                  only : kfl_funtn_sld

    implicit none

    integer(ip), intent(in) :: itask      !< What to do
    integer(ip)             :: ipoin,ifunc,itype,idime
    real(rp)                :: dinew(ndime)
    real(rp)                :: dummin, dummax

    select case (itask)

    case( ITASK_TURNON )

    case( ITASK_BEGSTE )

       !-------------------------------------------------------------------
       !
       ! Before a time step
       !
       !-------------------------------------------------------------------

       if( INOTEMPTY .and. kfl_conbc_sld == 0 ) then
          !
          ! Update boundary conditions (Transient)
          !
          do ipoin = 1,npoin
             !
             ! Initialize
             !
             dinew(:) = 0.0_rp
             ifunc    = kfl_funno_sld(ipoin)
             itype    = kfl_funtn_sld(ipoin)

             if(      itype == FUNCTION_MODULE ) then
                !
                ! Solidz Functions
                !
                call sld_funbou(2_ip, ipoin, dinew(:))

             else if( itype == FUNCTION_SPACE_TIME ) then
                !
                ! Kernel Space/Time Functions
                !
                call ker_functions(ipoin,ifunc,itype,bvess_sld(:,ipoin,ITER_AUX),dinew)

             end if
             !
             ! Update Dirichlet
             !
             bvess_sld(1:ndime,ipoin,ITER_K) = dinew(1:ndime)
          end do
          !
          ! Save value
          !
          do idime = 1,ndime
             if( abs(minval(bvess_sld(idime,:,ITER_K))) > zeror ) then
                rbbve_sld(idime,ITER_K) = minval(bvess_sld(idime,:,ITER_K))
             else
                rbbve_sld(idime,ITER_K) = maxval(bvess_sld(idime,:,ITER_K))
             end if
          end do

       end if
       !
       ! Parall service
       !
       do idime = 1,3
          dummin = rbbve_sld(idime,ITER_K)
          dummax = rbbve_sld(idime,ITER_K)
          call PAR_MIN(dummin,'IN MY CODE')
          call PAR_MAX(dummax,'IN MY CODE')
          if( dummin == 0.0_rp ) rbbve_sld(idime,ITER_K) = dummax
          if( dummax == 0.0_rp ) rbbve_sld(idime,ITER_K) = dummin
       end do

    case( ITASK_BEGITE )

       !-------------------------------------------------------------------
       !
       ! Before a global iteration
       !
       !-------------------------------------------------------------------

    case( ITASK_INNITE )

       !-------------------------------------------------------------------
       !
       ! Before an inner iteration
       !
       !-------------------------------------------------------------------

    end select

  end subroutine sld_rbo_updbcs

  !------------------------------------------------------------------------
  !> @author  Gerard Guillamet
  !> @date    September, 2022
  !> @brief   Restart variables for rigid body
  !> @details Restart variables for rigid body
  !------------------------------------------------------------------------

  subroutine sld_rbo_restar()

    use mod_restart, only : restart_add
    use def_solidz,  only : rbmas_sld
    use def_solidz,  only : rbfor_sld
    
    implicit none
   
    call restart_add(rbmas_sld,'rbmas_sld')
    call restart_add(rbfor_sld,'rbfor_sld')
    call restart_add(rbdil_sld,'rbdil_sld')
    call restart_add(rbvel_sld,'rbvel_sld')
    call restart_add(rbacl_sld,'rbacl_sld')
    call restart_add(uruk4_sld,'uruk4_sld')

  end subroutine sld_rbo_restar
  
end module mod_sld_rbo
