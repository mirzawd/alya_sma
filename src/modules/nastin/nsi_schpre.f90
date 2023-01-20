!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_schpre()
  !----------------------------------------------------------------------
  !****f* Nastin/nsi_schpre
  ! NAME 
  !    nsi_schpre
  ! DESCRIPTION
  !    This routine imposes boundary conditions on the preconditioner
  ! USES
  ! USED BY
  !    nsi_solbgs
  !***
  !----------------------------------------------------------------------
  use def_master
  use def_domain
  use def_nastin
  use def_solver
  use mod_memchk
  implicit none
  integer(ip) :: icomp,ipoin

  if( IMASTER ) return

  !----------------------------------------------------------------------
  !
  ! Compute preconditioner
  !
  !----------------------------------------------------------------------

  if( kfl_predi_nsi==5 ) then
     !
     ! Lumped Mass 
     !
     call runend('NSI_SCHPRE: NOT CODED')
     if( solve(2)%kfl_symme==1 ) then
        do ipoin=1,npoin
           icomp=r_sym(ipoin+1)-1 
           !lapla_nsi(icomp)=vmass(ipoin)/visco_nsi(1,1)
        end do
     else
        do ipoin=1,npoin
           icomp=r_dom(ipoin)
           do while(c_dom(icomp)/=ipoin)
              icomp=icomp+1
           end do
           !lapla_nsi(icomp)=vmass(ipoin)/visco_nsi(1,1)
        end do
     end if
  end if

end subroutine nsi_schpre

