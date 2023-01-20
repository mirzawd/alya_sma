!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine exm_nodset(insec,inset)
   !-----------------------------------------------------------------------
   !****f* Exmedi/exm_nodset
   ! NAME
   !    exm_nodset
   ! DESCRIPTION
   !    This routine computes variables on a node set
   ! USES
   ! USED BY
   !-----------------------------------------------------------------------
   use def_kintyp_basic,   only : ip
   use def_master,         only : postp, vnset
   use def_master,         only : elmag,vconc
   use mod_eccoupling,     only : kfl_cellmod,EXMSLD_CELL_TORORD
   use def_domain,         only : nodemat
   implicit none

   integer(ip), intent(in)  :: insec !< Node number
   integer(ip), intent(in)  :: inset !< Set code
   integer(ip)              :: ipoin,n,kmodel_imate

   if( insec .eq. 0 ) return

   !----------------------------------------------------------------------
   ! Initialisations
   !----------------------------------------------------------------------
   ipoin = insec
   !----------------------------------------------------------------------
   !
   ! Node sets
   !
   !----------------------------------------------------------------------
   !
   ! INTRA
   !

   if( postp(1) % npp_setsn(1) /= 0 ) then
           vnset(1,inset) =  elmag(ipoin,1)
   endif
   !
   ! CALCI
   !
   if( postp(1) % npp_setsn(2) /= 0 ) then
           n = nodemat(ipoin)
           kmodel_imate = kfl_cellmod(n)
           if (kmodel_imate == EXMSLD_CELL_TORORD) then
              vnset(2,inset) =  vconc(5,ipoin,1)
           else
              vnset(2,inset) =  vconc(1,ipoin,1)
           end if
   endif
   !
   ! SODI
   !
   if( postp(1) % npp_setsn(3) /= 0 ) then
           vnset(3,inset) =  vconc(3,ipoin,1)
   endif
   !
   ! POTAS
   !
   if( postp(1) % npp_setsn(4) /= 0 ) then
           vnset(4,inset) =  vconc(4,ipoin,1)
   endif


end subroutine exm_nodset
