!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup Neutro 
!> @{
!> @file    neu_turnon.f90
!> @author  Guillaume Houzeaux
!> @brief   Turn on module
!> @details Read data and allocate memory
!> @} 
!------------------------------------------------------------------------
subroutine neu_turnon()

  use def_kintyp, only : ip
  implicit none
  external :: neu_inivar, neu_reaphy, neu_parall, neu_reanut, &
              neu_reaous, neu_reabcs, neu_inibcs, neu_memall, &
              neu_directions, neu_readXS, neu_inbcdir

  !
  ! Initial variables
  !
  call neu_inivar(0_ip)
  !
  ! Read the physical problem
  !
  call neu_reaphy()  
  !
  ! Service: Parall
  !
  call neu_parall(1_ip) ! * verificar envio de los datos nuevos leidos*
  !
  ! Initial variables
  !
  call neu_inivar(2_ip)
  !
  ! Read the numerical treatment
  !
  call neu_reanut()
  !
  ! Read the output strategy
  !
  call neu_reaous()
  !
  ! Read the boundary conditions (just codes)
  !
  call neu_reabcs()  
  !
  ! Service: Parall
  !
  call neu_parall(2_ip)

  !
  ! Modify boundary conditions. Create Kfl_fixbp_ne structure. Just the codes
  !
  call neu_inibcs() 
  !
  ! Initial variables
  !
  call neu_inivar(1_ip)
  !
  ! Allocate memory
  !
  call neu_memall(1_ip)  ! 
  ! call neu_memall()  ! 
  !
  ! Compute directions ans reflexive direction for boundaries
  !
  call neu_directions()
  !
  ! Compute scattering
  !
  ! call neu_scattering()  ! it need to be reprogramed  Untill now ::  scattering_neu==1
  ! E. Goldberg (13/09/2021): commented, does nothing (sets an unused variable to 1). Scattering is read from input
  !
  ! read total and scatering cross sections 
  !
  call neu_readXS()  !  
  call neu_parall(4_ip)
  call neu_memall(2_ip)  ! 
 
  call neu_inbcdir() 
!
  call neu_parall(3_ip) ! envio matrices y vectores leidos de archivos!!! * OJO FALTA* 
  !

  ! Warnings and errors
  !
  !  call neu_outerr(1_ip)

end subroutine neu_turnon




