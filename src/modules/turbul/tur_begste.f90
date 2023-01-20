!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tur_begste()
  !-----------------------------------------------------------------------
  !****f* Turbul/tur_begste
  ! NAME 
  !    tur_begste
  ! DESCRIPTION
  !    This routine prepares for a new time step of the turbulence
  !    equations
  ! USES
  !    tur_iniunk
  !    tur_updtss
  !    tur_updbcs
  !    tur_updunk
  !    tur_radvuf
  ! USED BY
  !    Temper
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_turbul
  implicit none


  if(kfl_stead_tur/=1_ip) then
     !
     ! Initial guess fo the turbulence variables: f(n,0,*) <-- f(n-1,*,*).
     !
     call tur_updunk(ITASK_BEGSTE)
     !
     ! Update boundary conditions
     !
     call tur_updbcs(TUR_BEFORE_TIME_STEP)

  end if
  !
  ! Calculate maximum mixing lenth for Apsey and Castro limitation (k eps)
  !
  if (kfl_lmaxi_tur==1_ip.and.(mod(ittim,10_ip).eq.0_ip.or.ittim==1_ip)) & ! .and.TUR_FAMILY_K_EPS) &
       call tur_maxlen()
end subroutine tur_begste
