!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!subroutine exm_rcpupd
!
!!-----------------------------------------------------------------------
!!
!! This routine updates the recovery potential by an explicit update
!!
!!
!!       OLD SUBROUTINE. Should not be used anymore. Replaced by
!!       exm_odeite. Will desapear in future versions
!!
!!
!!-----------------------------------------------------------------------
!  use      def_parame
!  use      def_master
!  use      def_domain
!
!  use      def_exmedi
!
!  implicit none
!  integer(ip) :: NCOMP_INTRA_LAST, ipoin
!  real(rp)    :: xhea1, xhea2, xphic
!  real(rp)    :: xmeps, xchic, xmgam, xfact, epsch, tauwn, tauwp, tauvp, tauvn, phivn, tauv1, tauv2
!
!
!        !
!        ! Compute the new value
!        NCOMP_INTRA_LAST = 2
!
!        nodeloop: do ipoin=1,npoin
!
!            if (kfl_cellmod(nodemat(ipoin) ) == CELL_NOMOD_EXMEDI) then
!            !        
!            ! No Ionic Curent model (NOMOD)
!                cycle nodeloop
!            
!            !
!            ! 
!
!            !elseif (kfl_cellmod(nodemat(ipoin) )  == 12) then
!            !
!            !FK ionic curren model
!            !
!                !tauwp = xmopa_exm( 17,nodemat(ipoin))
!                !tauwn = xmopa_exm( 18,nodemat(ipoin))
!                !xphic = xmopa_exm( 3, nodemat(ipoin))
!                !tauvp = xmopa_exm( 4, nodemat(ipoin))
!                !phivn = xmopa_exm( 5, nodemat(ipoin))
!                !tauv1 = xmopa_exm( 7, nodemat(ipoin))
!                !tauv2 = xmopa_exm( 8, nodemat(ipoin))
!
!
!                !xhea1 = heavis(elmag(ipoin,NCOMP_INTRA_LAST),xphic)
!                !xhea2 = heavis(xphic,elmag(ipoin,NCOMP_INTRA_LAST))
!                !tauvn = heavis( elmag(ipoin,NCOMP_INTRA_LAST), phivn) * tauv1&
!                !+ heavis(phivn, elmag(ipoin,NCOMP_INTRA_LAST))  * tauv2
!            !   xfacw = 1._rp / dtime + xhea2 / tauwn * 0.5 + xhea1 / tauwp * 0.5
!            !   xfacv = 1._rp / dtime + xhea2 / tauvn * 0.5 + xhea1 / tauvp * 0.5
!
!
!            !   v in fenton's model
!            !    vgate_exm(1,ipoin,2) =  vgate_exm(1,ipoin,1) &                                 !
!            !        * ( 1.0_rp  - dtime/ xmccm_exm* xhea2 / tauvn &                            ! Commented because vgate does not exist anymore
!            !        - dtime/xmccm_exm * xhea1 / tauvp ) + dtime/xmccm_exm * xhea2 / tauvn      !
!
!           !        ( vgate_exm(2,ipoin,1
!           !    vgate_exm(2,ipoin,2) = &
!           !         ( vgate_exm(2,ipoin,1) * ( 1.0_rp / dtime - xhea2 / tauvn &
!           !        - xhea1 / tauvp * 0.5_rp ) + xhea2 / tauvn ) / xfacv
!           
!           ! w in fenton's model
!           !    vgate_exm(2,ipoin,2) =  vgate_exm(2,ipoin,1) &                                  !
!           !        * ( 1.0_rp  - dtime/xmccm_exm * xhea2 / tauwn &                             ! Commented because vgate does not exist anymore
!           !        - dtime/xmccm_exm * xhea1 / tauwp ) + dtime/xmccm_exm * xhea2 / tauwn       !
!           !  vgate_exm(1,ipoin,2) = &
!           !       ( vgate_exm(1,ipoin,1) * ( 1.0_rp / dtime - xhea2 / tauwn &
!           !       - xhea1 / tauwp * 0.5_rp ) + xhea2 / tauwn ) / xfacw
!           
!
!           endif
!
!
!        enddo nodeloop
!
!        
!           
!  ! Update  U(r,i,2) <-- unkno       
!  ! elmag(3,1:npoin,2) = unkno(1:npoin)
!
!end subroutine exm_rcpupd
