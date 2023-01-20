!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup Gusano
!> @{
!> @file    gus_iniunk.f90
!> @author  Guillaume Houzeaux
!> @brief   Initial value of velocity and pressurfe
!> @details Set up the initial condition for velocity and pressure.
!> @} 
!------------------------------------------------------------------------
subroutine gus_iniunk()
  
  use def_parame 
  use def_master
  use def_elmtyp
  use def_kermod
  use def_domain
  use def_gusano 
  implicit none

  if( kfl_rstar == 0 ) then  

     !-------------------------------------------------------------------
     !
     ! Apply boundary conditions
     !
     !-------------------------------------------------------------------
     
     call gus_updbcs(ITASK_INIUNK)

     !-------------------------------------------------------------------
     !
     ! (:,3) <= (:,1)
     !
     !-------------------------------------------------------------------
     
     call gus_updunk( ITASK_INIUNK )
     
  end if
 
end subroutine gus_iniunk

