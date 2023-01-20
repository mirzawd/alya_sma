!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup Alefor
!> @{
!> @file    ale_reanut.f90
!> @author  Guillaume Houzeaux
!> @date    February, 2017
!> @brief   Read numerical treatment
!>
!> @details Read numerical treatment
!>
!> @}
!------------------------------------------------------------------------
subroutine ale_reanut

 !.md<module>alefor
 !.md<input>case.ale.dat
 !.md<pos>1
 !.md<sec>

  use def_parame
  use def_inpout
  use def_master
  use def_solver
  use def_alefor
  use def_domain
  use mod_ecoute,             only : ecoute
  use mod_read_domain_arrays, only : read_domain_arrays_types
  implicit none

  if( INOTSLAVE ) then
     ! 
     !  Initializations (defaults)
     !
     kfl_smoot_ale     =  0                                      ! No mesh smoothing
     kfl_timef_ale     =  0                                      ! ALEFOR running coupled, with no timefunction
     nsmoo_ale         =  1                                      ! Number of smoothing steps (loading steps)
     kfl_defor_ale     =  0                                      ! No mesh deformation
     ndefo_ale         =  1                                      ! Number of deformation steps (loading steps)
     kfl_smobo_ale     =  0                                      ! Boundary smoothing
     kfl_fixsm_ale     =  1                                      ! Fix boundary nodes by default
     nsmob_ale         =  0                                      ! Number of boundary smoothing iterations
     kfl_crist_ale     =  0_ip                                   ! Cristobal 1 , new 0  ! RIGID BODY
     kfl_foexo_ale     =  2_ip                                   ! Force ( & Torque) extrapolation order  ! RIGID BODY
     kfl_disor_ale     =  2_ip                                   ! Integration order for the RB displacements  ! RIGID BODY
     kfl_nforc_ale     =  2_ip                                   ! Use the average of the previous 2 forces to calculate accel ! RB
     ansmo_ale         = -1.0_rp                                 ! Sharp edge detection
     resmo_ale         =  1.0_rp                                 ! Relaxation factor
     defor_param_ale   =  1.0_rp                                 ! Deformation parameters
     mnodb_ad          =  2
     !
     ! Reach the section
     !
     !.md<0># Numerical Treatment
     !.md<code>
     !.mdNUMERICAL_TREATMENT
     solve_sol => solve   
     call ecoute('ale_reanut')
     do while( words(1) /= 'NUMER' )
        call ecoute('ale_reanut')
     end do
     !
     ! Begin to read data
     !
     do while( words(1) /= 'ENDNU' )
        call ecoute('ale_reanut')

        if( words(1) == 'SMOOT' ) then
           !
           ! Mesh smoothing
           !
           if( exists('GAUSS') ) then
              kfl_smoot_ale = 2
           elseif( exists('CONST') ) then
              kfl_smoot_ale = 3
           else
              kfl_smoot_ale = 1
           end if
           if( exists('STEPS') ) nsmoo_ale = getint('STEPS',1_ip,  '#Number of smoothing steps') 
           if( exists('RELAX') ) resmo_ale = getrea('RELAX',1.0_rp,'#Smoothing relaxation factor') 

        else if( words(1) == 'DEFOR' ) then
           !
           ! Mesh deformation
           !
           if( exists('UNIFO') ) then
              kfl_defor_ale = 1
           else if( exists('SMALL') ) then
              kfl_defor_ale = 2
           else if( exists('ONLYB') ) then
              kfl_defor_ale = 4
           else if( exists('ISOTR') ) then
              kfl_defor_ale = 5
           else if( exists('BOUND') ) then
              kfl_defor_ale = 6
              if( exists('ALLEL') ) kfl_defor_ale = 7
              if( exists('TIMEF') ) kfl_timef_ale = 1
           else if( exists('DISTA') .or. exists('WALLD') ) then
              kfl_defor_ale = 8
           else
              kfl_defor_ale = 1
           end if
           if( exists('STEPS') ) ndefo_ale       = getint('STEPS',1_ip,'#Number of deformation steps') 
           if( exists('PARAM') ) defor_param_ale = getrea('PARAM',1.0_rp,'#Parameter for deformation algorithm') 
           
        else if( words(1) == 'BOUND' ) then
           !
           ! Boundary smoothing
           !
           kfl_smobo_ale = 1
           kfl_fixsm_ale = 0
           if( exists('STEPS') ) nsmob_ale = getint('STEPS',1_ip,'#Number of boundary smoothing steps') 
           if( exists('SHARP') ) ansmo_ale = getrea('SHARP',45.0_rp,'#Sharp edge angle') 

        else if( words(1) == 'FIXBO' ) then
           !
           ! Fix boundaries
           !
           if( words(2) == 'YES  ' .or. words(2) == 'ON   ' ) then
              kfl_fixsm_ale = 1
           else if( words(2) == 'NO   ' .or. words(2) == 'OFF  ' ) then
              kfl_fixsm_ale = 0
           end if

        else if( words(1) == 'ALGEB' ) then
           call reasol(1_ip)

        else if( words(1) == 'PRECO' ) then 
           call reasol(2_ip)

        else if( words(1) == 'RBTIM' ) then
           !
           ! Rigid Body time integration
           !
           if( words(2) == 'CRIST') then
              kfl_crist_ale = 1_ip
           else if( words(2) == 'SECON') then
              kfl_crist_ale = 0_ip
              kfl_foexo_ale = 2_ip                  ! Force ( & Torque) extrapolation order
              kfl_disor_ale = 2_ip                  ! Integration order for the RB displacements
           else if( words(2) == 'FIRST') then
              kfl_crist_ale = 0_ip
              kfl_foexo_ale = 1_ip                  ! Force ( & Torque) extrapolation order
              kfl_disor_ale = 1_ip                  ! Integration order for the RB displacements
           else if( words(2) == 'VSDF ') then
              kfl_crist_ale = 0_ip
              kfl_foexo_ale = 2_ip                  ! Force ( & Torque) extrapolation order
              kfl_disor_ale = 1_ip                  ! Integration order for the RB displacements
           else if( words(2) == 'VFDS ') then
              kfl_crist_ale = 0_ip
              kfl_foexo_ale = 1_ip                  ! Force ( & Torque) extrapolation order
              kfl_disor_ale = 2_ip                  ! Integration order for the RB displacements             
           else   ! default 2 2
              kfl_crist_ale = 0_ip
              kfl_foexo_ale = 2_ip                  ! Force ( & Torque) extrapolation order
              kfl_disor_ale = 2_ip                  ! Integration order for the RB displacements
           end if

        else if(words(1)=='FORCE') then
           if( words(2) == 'ORDIN' ) then
              kfl_nforc_ale = 1
           elseif( words(2) == 'AVERA' ) then
              kfl_nforc_ale = 2
           end if

        end if

     !.mdEND_NUMERICAL_TREATMENT
     !.md</code>

     end do
  end if

end subroutine ale_reanut
