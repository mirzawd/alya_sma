!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    sld_updunk.f90
!> @author  Guillaume Houzeaux
!> @date    August, 2006
!>          - Subroutine creation
!> @author  Gerard Guillamet and Mariano Vazquez
!> @date    July, 2018
!>          - Subroutine refactoring
!>          - OpenMP parallelization
!> @brief   This routine performs several types of updates for solid's
!>          unknowns.
!> @details
!>          \verbatim
!>          TIME STEP AND ITERATION LABELS (Defined in def_master)
!>            ITER_K       = 1 ..... Current iteration at time step N
!>            ITER_K_STATE = 1 ..... Current iteration at time step N for state variables
!>            ITER_AUX     = 2 ..... Used for coupling iterations
!>            TIME_N       = 3 ..... Time step (Converged)
!>            TIME_N_STATE = 2 ..... Time step (Converged) for state variables
!>
!>          PRE-PROCESS
!>            sld_iniunk (itask=6) .......... if Restart
!>
!>          TIME LOOP
!>            do time
!>              sld_begste (itask=1) ........ (:,ITER_AUX)     <= (:,TIME_N)
!>              do outer
!>                sld_begite (itask=2) ...... (:,ITER_K)       <= (:,ITER_AUX)
!>                do inner
!>                  sld_endite (itask=3) .... (:,ITER_K)       <= UNKNO
!>                end do
!>                sld_endite (itask=4p) ..... (:,ITER_AUX)     <= (:,ITER_K)
!>              end do
!>              sld_endste (itask=5) ........ (:,TIME_N)       <= (:,ITER_K)
!>              sld_endste (itask=7) ........ (:,TIME_N_STATE) <= (:,ITER_K_STATE)
!>            end do
!>          \endverbatim
!> @}
!-----------------------------------------------------------------------

subroutine sld_updunk(itask)

  use def_kintyp,             only : ip, rp
  use def_elmtyp,             only : ELINT
  use def_master,             only : INOTMASTER, displ, unkno
  use def_master,             only : ITER_AUX, ITER_K, TIME_N
  use def_master,             only : ITER_K_STATE, TIME_N_STATE
  use def_master,             only : ITASK_BEGSTE, ITASK_ENDSTE
  use def_master,             only : ITASK_BEGITE, ITASK_ENDITE
  use def_master,             only : ITASK_INNITE, ITASK_ENDINN
  use def_master,             only : ITASK_INIUNK
  use def_domain,             only : npoin, nelem, ndime
  use def_domain,             only : ltype, lelch, ngaus, nnode
  use def_solidz,             only : ndofn_sld, nsvar_sld
  use def_solidz,             only : veloc_sld, accel_sld
  use def_solidz,             only : kfl_sdvar_sld, svegm_sld
  use def_solidz,             only : SLD_DYNAMIC_PROBLEM

  implicit none

  integer(ip), intent(in)       :: itask  !> Update variables at selected case
  integer(ip)                   :: ipoin,idofn,idime,ielem,icomp
  integer(ip)                   :: pelty,pgaus,pnode

  if( INOTMASTER ) then

     select case (itask)

     case( ITASK_INIUNK )

        !-------------------------------------------------------------------
        !
        ! (:,*) <= (:,1) Initial solution
        !
        !-------------------------------------------------------------------
        do icomp = 1,size(displ,DIM=3,KIND=ip)
           do ipoin = 1,npoin
              displ(    1:ndime,ipoin,icomp) = displ(    1:ndime,ipoin,ITER_K)
              veloc_sld(1:ndime,ipoin,icomp) = veloc_sld(1:ndime,ipoin,ITER_K)
              accel_sld(1:ndime,ipoin,icomp) = accel_sld(1:ndime,ipoin,ITER_K)
           end do
        end do

     case( ITASK_BEGSTE )

        !------------------------------------------------------------------
        !
        ! (:,2) <= (:,1) Initial guess for outer iterations (begin time step)
        !
        !------------------------------------------------------------------
        do ipoin = 1,npoin
           displ(    1:ndime,ipoin,ITER_AUX) = displ(    1:ndime,ipoin,ITER_K)
           veloc_sld(1:ndime,ipoin,ITER_AUX) = veloc_sld(1:ndime,ipoin,ITER_K)
           accel_sld(1:ndime,ipoin,ITER_AUX) = accel_sld(1:ndime,ipoin,ITER_K)
        end do

     case( ITASK_BEGITE )

        !------------------------------------------------------------------
        !
        ! UNKNO <= (:,1) Initial guess for inner iterations (begin iterations)
        !
        !------------------------------------------------------------------
        do ipoin = 1,npoin
           do idime = 1,ndofn_sld
              idofn = (ipoin - 1)*ndofn_sld + idime
              unkno(idofn) = displ(idime,ipoin,ITER_AUX)
           end do
        end do

     case( ITASK_ENDINN )

        !------------------------------------------------------------------
        !
        ! (:,1) <= UNKNO Update
        !
        !------------------------------------------------------------------
        do ipoin = 1,npoin
           do idime = 1,ndofn_sld
              idofn = (ipoin - 1)*ndofn_sld + idime
              displ(idime,ipoin,ITER_K) = unkno(idofn)
           end do
        end do

     case( ITASK_ENDITE )

        !------------------------------------------------------------------
        !
        ! (:,2) <= (:,1) End of inner iterations (converged)
        !
        !------------------------------------------------------------------
        do ipoin = 1,npoin
           displ(    1:ndime,ipoin,ITER_AUX) = displ(    1:ndime,ipoin,ITER_K)
           veloc_sld(1:ndime,ipoin,ITER_AUX) = veloc_sld(1:ndime,ipoin,ITER_K)
           accel_sld(1:ndime,ipoin,ITER_AUX) = accel_sld(1:ndime,ipoin,ITER_K)
        end do

     case( ITASK_ENDSTE )

        !------------------------------------------------------------------
        !
        ! (:,3) <= (:,1) End of time step (converged)
        !
        !------------------------------------------------------------------
        do ipoin = 1,npoin
           displ(    1:ndime,ipoin,TIME_N) = displ(    1:ndime,ipoin,ITER_K)
           veloc_sld(1:ndime,ipoin,TIME_N) = veloc_sld(1:ndime,ipoin,ITER_K)
           accel_sld(1:ndime,ipoin,TIME_N) = accel_sld(1:ndime,ipoin,ITER_K)
        end do
        !
        ! Update State Dependent Variables (SDVs)
        !
        if( kfl_sdvar_sld == 1_ip ) then

           !--------------------------------------------------
           !$OMP PARALLEL DO SCHEDULE (STATIC)               &
           !$OMP DEFAULT  ( NONE )                           &
           !$OMP PRIVATE  ( ielem, pelty, pgaus, pnode )     &
           !$OMP SHARED   ( ltype, lelch, ngaus, nnode,      &
           !$OMP            nelem, nsvar_sld, svegm_sld )
           !--------------------------------------------------
           do ielem = 1,nelem
              pelty = ltype(ielem)
              pgaus = ngaus(pelty)
              pnode = nnode(pelty)
              if( lelch(ielem) == ELINT ) then
                 pgaus = pnode/2
              end if
              svegm_sld(ielem)%a(1:nsvar_sld,1:pgaus,TIME_N_STATE) = svegm_sld(ielem)%a(1:nsvar_sld,1:pgaus,ITER_K_STATE)
           end do
           !$OMP END PARALLEL DO
           !--------------------------------------------------

        end if

     end select

  end if

end subroutine sld_updunk

