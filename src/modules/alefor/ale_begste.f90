!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine ale_begste
  !-----------------------------------------------------------------------
  !****f* Alefor/ale_begste
  ! NAME 
  !    ale_begste
  ! DESCRIPTION
  !    This routine prepares for a new time step of the ALE formulation
  !    equation      
  ! USES
  !    ale_updunk
  ! USED BY
  !    Alefor
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_alefor
  implicit none
  integer(ip) :: idime

  if ( kfl_rigid_ale == 1 ) then
     !
     ! Obtain xrota_ale and xline_ale from nstro_ale and nstli_ale
     !
     xrota_ale = 0.0_rp
     xline_ale = 0.0_rp
     do idime =1,3_ip   
        if ( ittim >= nstro_ale(idime) )  xrota_ale(idime) = 1.0_rp
        if ( ittim >= nstli_ale(idime) )  xline_ale(idime) = 1.0_rp
     end do
     
  end if  
  !
  ! Initial guess for rigid body unknowns: a(n,0,*) <-- a(n-1,*,*).
  !
  call ale_updunk(ITASK_BEGSTE)
  !
  ! Use temporal predictor for zonal coupling (if any)
  !
  call ale_coupre()

!  if(ittim==6) then
!     if(kfl_paral==2) then
!        print*,'a=',dispm(:,1,1)    ! List of sets the rigid body is formed by
!        print*,'b=',dispm(:,1,2)    ! List of sets the rigid body is formed by
!        print*,'c=',dispm(:,1,3)    ! List of sets the rigid body is formed by
!     end if
!     call runend('O.K.!')
!  end if
  
end subroutine ale_begste

