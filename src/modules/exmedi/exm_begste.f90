!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine exm_begste
!-----------------------------------------------------------------------
!****f* Exmedi/exm_begste
! NAME 
!    exm_begste
! DESCRIPTION
!    This routine prepares for a new time step 
! USES
!    exm_iniunk
!    exm_updtss
!    exm_updbcs
!    exm_updunk
! USED BY
!    Exmedi
!***
!-----------------------------------------------------------------------

  use def_parame
  use def_master
  use def_domain

  use def_exmedi


  implicit none
  integer(ip) :: ipoin

  !
  ! Initial guess fo the unknowns: U(n,0,*) <-- U(n-1,*,*).
  !
  call exm_updunk(ITASK_BEGSTE)     ! u(,ITER_AUX) <-- u(,TIME_N) 

  !reset the isochrones
  if (kfl_reset_fisoc) then
    if( INOTMASTER ) then
        do ipoin = 1,npoin
          fisoc(:,ipoin) = -1.0_rp
        end do
    end if
  
    kfl_reset_fisoc = .FALSE.
  end if


end subroutine exm_begste
