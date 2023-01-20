!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine chm_cvgunk(itask)
  !-----------------------------------------------------------------------
  !****f* Chemic/chm_cvgunk
  ! NAME
  !    chm_cvgunk
  ! DESCRIPTION
  !    Convergence criterion
  ! USES
  ! USED BY
  !    chm_doiter
  !***
  !-----------------------------------------------------------------------
  use def_parame,              only : one, zero
  use def_master,              only : ITASK_ENDINN, ITASK_ENDITE, ITASK_ENDSTE, momod, tempe_gp, conce, cpu_initi, cutim, IMASTER,&
                                      INOTMASTER, INOTSLAVE, IPARALL, itcou, ittim, kfl_rstar, modul, tempe, unkno, itinn
  use def_domain,              only : nelem, npoin, ltype, ngaus
  use def_kintyp,              only : ip, rp
  use def_chemic,              only : comax_chm, comin_chm, cotol_chm, hrr_int_chm, iclaf_chm, iclai_chm, kfl_goit2_chm,&
                                      kfl_goite_chm, kfl_model_chm, kfl_normc_chm, kfl_solve_cond_CMC_chm,&
                                      kfl_transfer_condField_CMC_chm, nclas_chm, nZ_CMC_chm, resid_chm, rtpts_chm, sstol_chm,&
                                      ripts_chm, Text_uncond_CMC_chm, Text_cond_CMC_chm
  use mod_outfor,              only : outfor
  use mod_array_operations,    only : array_operations_residual_norm
  use mod_array_operations,    only : array_operations_min_max
  use mod_communications,      only : PAR_MAX
  use mod_chm_operations_CMC,  only : chm_values_cvgunk_CMC
  use mod_communications,      only : PAR_MAX

  implicit none
  integer(ip), intent(in) :: itask
  integer(ip), save       :: ipass=0
  integer(ip)             :: iclas,kpoin,iZ
  integer(ip)             :: ielem, igaus, pelty, pgaus
  real(rp),    save       :: cpuit_chm=0.0_rp
  real(rp)                :: time1,comin,comax,temin,temax,tminmax(2)
  real(rp)                :: rtpts,dummr

  external                :: memerr
  external                :: minmax
  external                :: cputim

  select case(itask)

  case(ITASK_ENDINN)

     !
     ! Check convergence of the inner iterations:
     ! || c(n,i,j) - c(n,i,j-1)|| / ||c(n,i,j)||
     !
     do iclas = iclai_chm,iclaf_chm

        if( kfl_normc_chm /= 3 ) then
           ripts_chm(iclas) = array_operations_residual_norm(kfl_normc_chm,nclas_chm,1_ip,unkno,conce,iclas-1,npoin*(iclas-1),1_ip)
           !!!!!!!!!!!!!!!!!!!!!! CHECK IF THIS IS CORRECT FOR CMC MODEL
           !!!if (kfl_model_chm == 4 .and. kfl_solve_enth_CMC_chm /= 0) then
           !!!   if (iclas > nclas_chm) then
           !!!      call residu(&
           !!!           kfl_normc_chm,nvar_CMC_chm,one,unkno(1),therm(1,1),&
           !!!           iclas,one,one,1.0_rp,ripts_chm(iclas))
           !!!   else
           !!!      call residu(&
           !!!           kfl_normc_chm,nvar_CMC_chm,one,unkno(1),conce(1,iclas,1),&
           !!!           iclas,one,one,1.0_rp,ripts_chm(iclas))
           !!!   end if

           !!!else
           !!!   call residu(&
           !!!        kfl_normc_chm,nclas_chm,one,unkno(1),conce(1,iclas,1),&
           !!!        iclas,one,one,1.0_rp,ripts_chm(iclas))
           !!!end if
        end if

        if( ripts_chm(iclas) < cotol_chm ) then
           kfl_goit2_chm = kfl_goit2_chm + 1
        endif

        rtpts_chm = rtpts_chm + ripts_chm(iclas)

     end do

     !
     ! Stop or go on
     !
     if( kfl_goit2_chm == nclas_chm ) kfl_goite_chm = 0

     !!DMM if( kfl_spray_chm == 2 ) miinn_chm = 0

     !! DMM
     !if ( itinn(modul) >= miinn_chm ) kfl_goite_chm = 0
     kfl_goite_chm = 0
     !
     ! Compute min and max
     !
     kpoin    =  1
     do iclas = iclai_chm,iclaf_chm
        if( IMASTER ) then
           call minmax(one,npoin,zero,dummr,comin,comax)
        else
           call minmax(one,npoin,zero,unkno(kpoin),comin,comax)
        end if
        comin_chm = min(comin_chm,comin)
        comax_chm = max(comax_chm,comax)
        if( INOTMASTER ) kpoin = kpoin + npoin
     end do

     !
     ! Compute min and max temperature
     !
     temin = 0.0_rp
     temax = 0.0_rp
     if (associated(tempe)) then
         !
         ! Get limits of tempe on npoin
         !
         call array_operations_min_max(1_ip,npoin,0_ip,tempe,temin,temax)

     elseif (associated(tempe_gp)) then
         !
         ! Get limits of Gauss point temperatures
         !
         temin= huge(1.0_rp)
         temax=-huge(1.0_rp)
         do ielem = 1,nelem
            pelty = ltype(ielem)
            pgaus = ngaus(pelty)
            do igaus = 1,pgaus
                if (tempe_gp(ielem)%a(igaus,1,1) > temax) temax = tempe_gp(ielem)%a(igaus,1,1)
                if (tempe_gp(ielem)%a(igaus,1,1) < temin) temin = tempe_gp(ielem)%a(igaus,1,1)
            enddo
         enddo

         if( IPARALL ) then
            tminmax(1) = -temin
            tminmax(2) =  temax
            call PAR_MAX(2_ip,tminmax)
            temin = -tminmax(1)
            temax =  tminmax(2)
         end if

     endif
     !
     ! Write convergence
     !
     if( INOTSLAVE ) then
        ! Combustion saves species convergence by default
        if ( kfl_model_chm==1 .or. (kfl_model_chm==4 .and. kfl_solve_cond_CMC_chm == 0)) then
           !=============================================!
           ! FLAMELET MODEL AND CMC FOR MIXING VARIABLES !
           !=============================================!
           call cputim(time1)
           if( ipass == 0 ) then
              time1 = time1-cpu_initi
           else
              time1 = time1-cpuit_chm
           end if
           if (kfl_model_chm==4 .and. kfl_solve_cond_CMC_chm == 0 .and. &
                  kfl_transfer_condField_CMC_chm > 0) then
              if( ipass == 0 .and. kfl_rstar /= 2 ) write(momod(modul)%lun_conve,116)
              write(momod(modul)%lun_conve,111) &
                   ittim,itcou,itinn(modul),cutim,rtpts_chm,comin_chm,&
                   comax_chm,time1,hrr_int_chm,(ripts_chm(iclas),iclas = 1,nclas_chm)
           else
              if( ipass == 0 .and. kfl_rstar /= 2 ) write(momod(modul)%lun_conve,110)
              write(momod(modul)%lun_conve,111) &
                   ittim,itcou,itinn(modul),cutim,rtpts_chm,comin_chm,&
                   comax_chm,time1,(ripts_chm(iclas),iclas = 1,nclas_chm)
           end if
           call cputim(cpuit_chm)
           flush(momod(modul)%lun_conve)
           ipass = ipass + 1

        elseif ( kfl_model_chm==2 .or. kfl_model_chm==3 ) then
           !
           ! MIXED EQUATIONS MODEL OR FINITE RATE CHEMISTRY MODEL
           !
           call cputim(time1)
           if( ipass == 0 ) then
              time1 = time1-cpu_initi
           else
              time1 = time1-cpuit_chm
           end if
           if( ipass == 0 .and. kfl_rstar /= 2 ) write(momod(modul)%lun_conve,112)
           write(momod(modul)%lun_conve,113) &
                ittim,itcou,itinn(modul),cutim,rtpts_chm,comin_chm,&
                comax_chm,time1,hrr_int_chm,temin,temax,(ripts_chm(iclas),iclas = 1,nclas_chm)
           call cputim(cpuit_chm)
           flush(momod(modul)%lun_conve)
           ipass = ipass + 1

        endif

     end if

  case(ITASK_ENDITE)
     !
     ! Check convergence of the outer iterations:
     ! || c(n,i,*) - c(n,i-1,*)|| / ||c(n,i,*)||
     !
     resid_chm = 0.0_rp
     do iclas = 1,nclas_chm
        ripts_chm(iclas) = array_operations_residual_norm(kfl_normc_chm,1_ip,1_ip,conce,conce,0_ip,1_ip*npoin*nclas_chm,1_ip)
        resid_chm = resid_chm + ripts_chm(iclas)
     end do

  case(ITASK_ENDSTE)
     !
     ! Check residual of the time iterations:
     ! || c(n,*,*) - c(n-1,*,*)|| / ||c(n,*,*)||
     !
     rtpts = 0.0_rp
     do iclas = 1,nclas_chm
        ripts_chm(iclas) = array_operations_residual_norm(kfl_normc_chm,1_ip,1_ip,conce,conce,0_ip,2_ip*npoin*nclas_chm,1_ip)
        rtpts = rtpts + ripts_chm(iclas)
     end do

     rtpts = rtpts / real(nclas_chm,rp)
     if( rtpts <= sstol_chm ) then
        momod(modul) % kfl_stead = 1
        call outfor(28_ip,momod(modul)%lun_outpu,' ')
     end if

     if ( kfl_model_chm == 4_ip .and. kfl_solve_cond_CMC_chm == 1_ip ) then
        !
        !==================================!
        ! CONDITIONAL MOMENT CLOSURE MODEL !
        !==================================!
        !
        call chm_values_cvgunk_CMC

        if (INOTSLAVE) then
           call cputim(time1)
           if( ipass == 0 ) then
              time1 = time1-cpu_initi
           else
              time1 = time1-cpuit_chm
           end if
           if( ipass == 0 .and. kfl_rstar /= 2 ) write(momod(modul)%lun_conve,114)

           write(momod(modul)%lun_conve,115) &
                ittim, cutim, time1, Text_uncond_CMC_chm(1), Text_uncond_CMC_chm(2), &
                hrr_int_chm, (Text_cond_CMC_chm(iZ),iZ = 1,2*nZ_CMC_chm)  !!, &
                !! DESCOMENTAR
                !!(sum_Yk_ext_cond_CMC_chm(iZ),iZ = 1,2*nZ_CMC_chm)
           call cputim(cpuit_chm)
           flush(momod(modul)%lun_conve)
           ipass = ipass + 1
        end if

     endif


  end select
  !
  ! Formats
  !

110 format('# --| ALYA Convergence, Combustion model '  ,/,&
       &   '# --| Columns displayed:' ,/,&
       &   '# --| 1. Time step         2. Global Iteration   3. Inner Iteration   ' ,/,&
       &   '# --| 4. Current time      5. Concentration      6. Min. value        ' ,/,&
       &   '# --| 7. Max. value        8. Elapsed CPU Time   9. Residual Spec 1,2,... ' ,/,&
       &   '# ','          1','          2','          3',&
       &        '             4','             5','             6','             7',&
       &        '             8','             9...')
111 format(4x,i9,2x,i9,2x,i9,1000(2x,e13.6))

112 format('# --| ALYA Convergence, Combustion model '  ,/,&
       &   '# --| Columns displayed:' ,/,&
       &   '# --| 1. Time step         2. Global Iteration   3. Inner Iteration   ' ,/,&
       &   '# --| 4. Current time      5. Concentration      6. Min. value        ' ,/,&
       &   '# --| 7. Max. value        8. Elapsed CPU Time   9. Heat release (W)  ' ,/,&
       &   '# --|10. Min. temperature 11. Max. temperature  12. Residual Spec 1,2,... ' ,/,&
       &   '# ','          1','          2','          3',&
       &        '             4','             5','             6','             7',&
       &        '             8','             9','            10','            11',&
       &        '            12')
113 format(4x,i9,2x,i9,2x,i9,1000(2x,e13.6))

114 format('# --| ALYA Convergence, Combustion model '  ,/,&
       &   '# --| Columns displayed:' ,/,&
       &   '# --| 1. Time step         2. Current time       3. Elapsed CPU Time  ' ,/,&
       &   '# --| 4. Min. T            5. Max. T             6. Heat release (W)  ' ,/,&
       &   '# --| 7. Min/Max T|Z... '                                               ,/,&
       &   '# ','        Min T1','        Max T1','        Min T2',&
       &        '        Max T2','        Min T3','        Max T3',&
       &        '        Min T4','        Max T4','           ...')  !,&
!!       &   '# ','        Min Sum_Y1','    Max Sum_Y1','    Min Sum_Y2',&
!!       &        '        Max Sum_Y2','    Min Sum_Y3','    Max Sum_Y3',&
!!       &        '        Min Sum_Y4','    Max Sum_Y4','       ...')
115 format(4x,i9,1000(2x,e13.6))

116 format('# --| ALYA Convergence, Combustion model '  ,/,&
       &   '# --| Columns displayed:' ,/,&
       &   '# --| 1. Time step         2. Global Iteration   3. Inner Iteration   ' ,/,&
       &   '# --| 4. Current time      5. Concentration      6. Min. value        ' ,/,&
       &   '# --| 7. Max. value        8. Elapsed CPU Time   9. Heat release (W)  ' ,/,&
       &   '# --| 10. Residual Spec 1,2,... ' ,/,&
       &   '# ','          1','          2','          3',&
       &        '             4','             5','             6','             7',&
       &        '             8','             9...')

end subroutine chm_cvgunk

