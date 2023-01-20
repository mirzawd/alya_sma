!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine lev_begste
!-----------------------------------------------------------------------
!****f* Wavequ/lev_begste
! NAME 
!    lev_begste
! DESCRIPTION
!    This routine prepares for a new time step of the level set
!    convection equation      
! USES
!    lev_updtss
!    lev_updbcs
!    lev_updunk
!    lev_radvuf
! USED BY
!    Temper
!***
!-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_levels
  implicit none
  integer(ip), save :: ipass=0

  if(ipass==0) then
     !
     ! First time we pass here: do not use ittim as this can be 
     ! a restart run 
     !
     ipass=1
     call lev_inivar(2_ip)
  end if
 
  if(kfl_stead_lev/=1) then     
     !
     ! Initial guess fo the temperature: u(n,0,*) <-- u(n-1,*,*).
     !
     call lev_updunk(1_ip)
  end if


end subroutine lev_begste

