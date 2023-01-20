!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    sld_cvgunk.f90
!> @author  Gerard Guillamet
!> @date    February, 2017
!>           - Subroutine written
!> @date    February, 2018
!>           - Add energies
!> @brief   Convergence checks for Solidz module
!> @details
!>          Convergence criteria implemented:
!>           - Residual error criterion (Belytschko p354)
!>           - Displ. increment error criterion (Belytschko p354)
!>           - Energy convergence criterion (Belytschko p354)
!>
!>          \verbatim
!>          ITASK = 0 ... Writting and timings output
!>                  1 ... Convergence inner iterations
!>                  2 ... Convergence outer iterations
!>                  3 ... Check residual time iterations
!>          \endverbatim
!>
!>          References:\n
!>          T. Belytschko, W. K. Liu, B. Moran, K. I. Elkhodary
!>          Nonlinear Finite elements for Continua and Structures.
!>
!> @todo    To do list::\n
!>           - Revision rhsid with and without fixity for force and energy
!>             criterions.
!>           - Include natural forces (from fsi problems) in the energy
!>             criterion.
!>
!> @}
!------------------------------------------------------------------------

subroutine sld_cvgunk(itask)

  use def_kintyp,            only : ip, rp
  use def_master,            only : zeror, isnain, INOTSLAVE, ITER_K, TIME_N
  use def_master,            only : ITASK_ENDINN, ITASK_ENDITE, ITASK_ENDSTE
  use def_master,            only : unkno, displ, rhsid, modul, kfl_rstar
  use def_master,            only : cpu_initi, momod, retost
  use def_master,            only : itinn, itcou, ittim, cutim, dtime
  use def_domain,            only : ndime, npoin
  use mod_outfor,            only : outfor
  use mod_iofile,            only : iofile_flush_unit
  use mod_messages,          only : messages_live
  use mod_array_operations,  only : array_operations_residual_norm
  use mod_array_operations,  only : array_operations_norm2
  use def_solidz,            only : last_iters_sld
  use def_solidz,            only : kfl_resid_sld, kfl_goite_sld
  use def_solidz,            only : kfl_timei_sld, kfl_stead_sld
  use def_solidz,            only : fint2_sld, fext2_sld, fine2_sld, fnatu_sld
  use def_solidz,            only : allie_sld, allwk_sld, allke_sld, etota_sld
  use def_solidz,            only : resid_sld, eener_sld
  use def_solidz,            only : miinn_sld, dtcri_sld, safet_sld
  use def_solidz,            only : cotol_sld, sstol_sld, kfl_damag_sld, kfl_cohes_sld
  use def_solidz,            only : SLD_STATIC_PROBLEM
  use mod_sld_cardiac_cycle, only : kfl_cardiac_cycle, sld_cardiac_cycle_write_cvg

  implicit none

  external                :: sld_minmax
  
  integer(ip), intent(in) :: itask !< What to do
  integer(ip), save       :: kpass = 0
  real(rp),    save       :: cpuit_sld=0.0_rp
  real(rp),    save       :: risld,rfsld,resld,resti=1.0_rp
  real(rp)                :: time1,dtdtc
  real(rp)                :: fmax2, rhsn2, emax
  real(rp)                :: cpu_ass_ave
  real(rp)                :: cpu_ass_max
  real(rp)                :: cpu_sol_ave
  real(rp)                :: cpu_sol_max
  real(rp)                :: displ_max

  select case ( itask )

  case ( 0_ip )

     !-------------------------------------------------------------------
     !
     ! Timings and output
     !
     !-------------------------------------------------------------------

     !
     ! Compute some minmax values
     !
     call sld_minmax(cpu_ass_ave,cpu_ass_max,cpu_sol_ave,cpu_sol_max,displ_max)
     !
     ! Write convergence
     !
     if( INOTSLAVE ) then
        call cputim(time1)
        !
        ! Write heading
        !
        if( kpass == 0 .and. kfl_rstar /= 2 ) then
           write(momod(modul) % lun_conve,100)
           if( kfl_cardiac_cycle ) call sld_cardiac_cycle_write_cvg(1_ip)

        end if

        if( kpass == 1 ) then
           time1 = time1 - cpuit_sld
        else
           time1 = time1 - cpu_initi
        end if
        
        if( dtcri_sld /= 0.0_rp ) then
           dtdtc = dtime / dtcri_sld
        else
           dtdtc = 0.0_rp
        end if
        !
        ! Write results
        !
        write(momod(modul) % lun_conve,101) &
             !   1     2            3              4     5         6
             ittim,itcou,itinn(modul),last_iters_sld,cutim,dtcri_sld,&
             !   7     8     9        10                11                12
             risld,rfsld,resld,displ_max,allie_sld(ITER_K),allwk_sld(ITER_K),&
             !              13                14
             allke_sld(ITER_K),etota_sld(ITER_K),&
             !  15          16          17          18          19    20         21  
             time1,cpu_ass_ave,cpu_ass_max,cpu_sol_ave,cpu_sol_max,dtdtc, safet_sld

        call cputim(cpuit_sld)
        call iofile_flush_unit(momod(modul) % lun_conve)
        !
        ! Write results cardiac cycle
        !
        if( kfl_cardiac_cycle ) call sld_cardiac_cycle_write_cvg(2_ip)
        
     end if

     kpass = 1_ip

  case ( ITASK_ENDINN )

     !-------------------------------------------------------------------
     !
     ! Check convergence of the inner iterations (N-R solidz iterations)
     !
     ! RISLD: Displacement residual
     ! RFSLD: Force residual
     ! RESLD: Energy residual
     ! 
     !-------------------------------------------------------------------
     !
     ! Displacement increment error (Belytschko p332)
     ! || d(n,i,j) - d(n,i,j-1)|| / ||d(n,i,j)|| (L2 norm)
     !
     risld = array_operations_residual_norm(2_ip,ndime,ndime,unkno,displ,0_ip,0_ip,ndime,1.0_rp)
     !
     ! Force residual error (Belytschko p332)
     ! ||rhs|| / max( ||fext||, ||fint||, ||M*a|| ) (L2 norm)
     !
     ! Force residual error (Calculated by solver)
     ! ||rhs|| <= cotol_sld 
     !
     rhsn2 = array_operations_norm2(rhsid,DOFS=ndime)
     fmax2 = max(fint2_sld, fext2_sld, fine2_sld, fnatu_sld)
     if( kfl_damag_sld == 1_ip .or. kfl_cohes_sld == 1_ip .or. &
         risld < zeror .or. fmax2 < zeror ) then
        rfsld = rhsn2 ! Algebraic value
     else
        rfsld = rhsn2/fmax2
     end if
     call messages_live(' |RES|= '//trim(retost(rfsld,REAL_FORMAT='(3E10.2)')), ADVANCE='no')
     !
     ! Energy convergence criterion (Belytschko p354)
     ! |du*rhs | / max( Wint, Wext, Wkin )
     !
     emax = max(abs(allie_sld(ITER_K)), abs(allwk_sld(ITER_K)), allke_sld(ITER_K))
     resld = eener_sld
     if( emax > zeror ) then
        resld = eener_sld/emax
     end if
     !resld = eener_sld/(emax+zeror)
     !
     ! KFL_GOITE_SLD: Check convergence
     !
     if( isnain(risld) ) then
        call runend("SOLIDZ: Nan or +/- Inf in displacement residual")
     end if
     if( isnain(rfsld) ) then
        call runend("SOLIDZ: Nan or +/- Inf in force residual")
     end if
     if( isnain(resld) ) then
        call runend("SOLIDZ: Nan or +/- Inf in energy residual")
     end if
     if( risld > 1.0e10_rp ) then
        call runend("SOLIDZ: Huge displacement residual")
     end if
     if( rfsld > 1.0e10_rp ) then
        call runend("SOLIDZ: Huge force residual")
     end if
     if( resld > 1.0e10_rp ) then
        !call runend("SOLIDZ: Huge energy residual")
     end if
     !
     ! Types of convergence criteria used to terminate iterations:
     !
     if (kfl_resid_sld == 1_ip) then
        !
        ! FORCE: Only force residual (Belytschko)
        !
        if ((rfsld <= cotol_sld) .or. &                       ! Inner tolerance achieved
             itinn(modul) >= miinn_sld) kfl_goite_sld = 0_ip  ! Number of iterations

     else if (kfl_resid_sld == 2_ip) then
        !
        ! ENERGY: Only energy residual (Belytschko)
        !
        if (resld <= cotol_sld .or. &                         ! Inner tolerance achieved
             itinn(modul) >= miinn_sld) kfl_goite_sld = 0_ip  ! Number of iterations

     else if (kfl_resid_sld == 3_ip) then
        !
        ! TOTAL: Displacement and force residuals (Belytschko)
        !
        if ((risld <= cotol_sld .and. rfsld <= cotol_sld &
             ) .or. &                                         ! Inner tolerance achieved
             itinn(modul) >= miinn_sld) kfl_goite_sld = 0_ip  ! Number of iterations
     else
        !
        ! DISPL: Only displacement residual (Belytschko)
        !
        if (risld <= cotol_sld .or. &                         ! Inner tolerance achieved
             itinn(modul) >= miinn_sld) kfl_goite_sld = 0_ip  ! Number of iterations

     end if

  case ( ITASK_ENDITE )

     !-------------------------------------------------------------------
     !
     ! Check convergence of the outer iterations
     ! External/global iterations in a coupling problem
     ! RESID_SLD: only displacement residual
     !
     ! ||u^{n+1,i+1}-u^{n+1,i-1}||/||u^{n+1,i+1}||
     !-------------------------------------------------------------------

     resid_sld = array_operations_residual_norm(2_ip,ndime,ndime,displ,displ,0_ip,npoin*ndime,ndime)

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

100 format('# --| ALYA Convergence',/,&
       &   '# --| Columns displayed:' ,/,&
       &   '# --|                                                                      ',/,&
       &   '# --|  1. Time step           2. Global Iteration    3. Inner Iteration    ',/,&
       &   '# --|  4. Solver iterations   5. Current time        6. Stable time incr.  ',/,&
       &   '# --|  7. RSI Displacement    8. RSI Force           9. RSI Energy         ',/,&
       &   '# --| 10. Max. Displacement  11. Internal energy    12. External work      ',/,&
       &   '# --| 13. Kinetic energy     14. Total energy       15. Elapsed CPU time   ',/,&
       &   '# --| 16. Ass. ave cpu time  17. Ass. max cpu time  18. Sol. ave cpu time  ',/,&
       &   '# --| 19. Sol. max cpu time  20. dt/dtc             21. Safety factor      ',/,&
       &   '# --|                                                                      ')

101 format(4x,i9,2x,i9,2x,i9,2x,i9,50(2x,e16.8e3))

end subroutine sld_cvgunk
