!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine chm_updunk(itask)
  !-----------------------------------------------------------------------
  !****f* Chemic/chm_updunk
  ! NAME
  !    chm_updunk
  ! DESCRIPTION
  !    This routine performs several types of updates
  ! USED BY
  !    chm_begste (itask=1)
  !    chm_begite (itask=2)
  !    chm_endite (itask=3, inner loop)
  !    chm_endite (itask=4, outer loop)
  !    chm_endste (itask=5)
  !    chm_restar (itask=6)
  !***
  !-----------------------------------------------------------------------
  use def_master, only : INOTEMPTY, ITASK_INIUNK, ITASK_BEGITE, ITASK_BEGINN, ITASK_BEGSTE, ITASK_ENDINN, ITASK_ENDITE,&
                         ITASK_ENDSTE, ITASK_INIUNK, ITASK_INNITE, conce, unkno, therm
  use def_kintyp, only : ip
  use def_domain, only : npoin
  use def_chemic, only : iclaf_chm, iclai_chm, kfl_model_chm, kfl_solve_cond_CMC_chm, kfl_solve_enth_CMC_chm, nclas_chm, ncomp_chm,&
                         ADR_chm
  use mod_ADR,    only : ADR_end_time_step
  use mod_ADR,    only : ADR_end_inner_iteration
  use mod_ADR,    only : ADR_begin_time_step

  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: iclas,kpoin,ipoin,icomp,nprev_chm

  nprev_chm = min(3_ip,ncomp_chm)

  if( INOTEMPTY ) then

     select case (itask)

     case ( ITASK_INIUNK )
        !
        ! (:,3) <= (:,1): Initial solution
        !
        do icomp = 2,ncomp_chm
           do iclas = iclai_chm,iclaf_chm
              do ipoin = 1,npoin
                 conce(ipoin,iclas,icomp) = conce(ipoin,iclas,1)
              end do
           end do
        end do

     case( ITASK_BEGSTE )
        !
        ! (:,2) <= (:,1): before a time step
        !
        do iclas=1,nclas_chm
           call ADR_begin_time_step(ADR_chm(iclas),conce(:,iclas,:))
        end do

     case( ITASK_BEGITE, ITASK_BEGINN )
        !
        ! UNKNO <= CONCE(:,CLASS,1)
        !

        if (kfl_model_chm == 4 .and. kfl_solve_cond_CMC_chm == 1) then

           if (kfl_solve_enth_CMC_chm /= 0) then
              ! Assign unkno
              kpoin = 0
              do ipoin = 1,npoin
                 do iclas = 1,nclas_chm
                    kpoin = kpoin+1
                    unkno(kpoin) = conce(ipoin,iclas,1)
                 end do
                 kpoin = kpoin+1
                 unkno(kpoin) = therm(ipoin,1)
              end do

           else
              ! Assign unkno
              kpoin = 0
              do ipoin = 1,npoin
                 do iclas = 1,nclas_chm
                    kpoin = kpoin+1
                    unkno(kpoin) = conce(ipoin,iclas,1)
                 end do
              end do

           end if


        else
           kpoin = 0
           do ipoin = 1,npoin
              do iclas = iclai_chm,iclaf_chm
                 kpoin = kpoin+1
                 unkno(kpoin) = conce(ipoin,iclas,1)
              end do
           end do
        end if


     case( ITASK_INNITE, ITASK_ENDINN)
        !
        ! (:,1) <= UNKNO
        !
        if (kfl_model_chm == 4 .and. kfl_solve_cond_CMC_chm == 1) then
           if (kfl_solve_enth_CMC_chm /= 0) then
              ! Assign conce and therm
              kpoin = 0
              do ipoin = 1,npoin
                 do iclas = 1,nclas_chm
                    kpoin = kpoin+1
                    conce(ipoin,iclas,1) = unkno(kpoin)
                 end do
                 kpoin = kpoin+1
                 therm(ipoin,1) = unkno(kpoin)
              end do
           else
              ! Assign conce and therm
              kpoin = 0
              do ipoin = 1,npoin
                 do iclas = 1,nclas_chm
                    kpoin = kpoin+1
                    conce(ipoin,iclas,1) = unkno(kpoin)
                 end do
              end do
           end if

        else
           kpoin = 0
           do ipoin = 1,npoin
              do iclas = iclai_chm,iclaf_chm
                 kpoin = kpoin+1
                 conce(ipoin,iclas,1) = unkno(kpoin)
              end do
           end do
        end if


     case( ITASK_ENDITE )
        !
        ! (:,2) <= (:,1)
        !
        do iclas = 1,nclas_chm
           call ADR_end_inner_iteration(ADR_chm(iclas),conce(:,iclas,:))
        end do

     case( ITASK_ENDSTE )
        !
        ! (:,3) <= (:,1): End of time step
        ! (:,4) <= (:,3)
        ! (:,5) <= (:,4)
        ! ...
        !
        do iclas = 1,nclas_chm
           call ADR_end_time_step(ADR_chm(iclas),conce(:,iclas,:))
        end do

     end select

  end if

end subroutine chm_updunk

