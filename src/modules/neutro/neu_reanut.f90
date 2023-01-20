!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup NeutroInput
!> @ingroup    Neutro
!> @{
!> @file    neu_reanut.f90
!> @author  Guillaume Houzeaux
!> @brief   Read numerical data
!> @details Read numerical data
!> @} 
!-----------------------------------------------------------------------

subroutine neu_reanut()

  use def_parame
  use def_inpout
  use def_master
  use def_neutro
  use def_domain
  use def_solver
  use mod_ADR,   only : ADR_read_data
  use mod_ecoute, only :  ecoute
  implicit none
  external :: reasol

  if( INOTSLAVE ) then
     !
     !  Initializations (defaults)
     !
     miinn_neu     = 10         ! Maximum inner iterations
     cotol_neu     = 1.0e-6_rp  ! Tolerance inner iterations
     relax_neu     = 1.0_rp     ! Relaxation parameter
     nitsche_neu   = 1000.0_rp  ! Nitsche coefficient
     kfl_smobo_neu = 0          ! Boundary condition is not smoothed
     !
     ! Reach the section
     !
     call ecoute('neu_reanut')
     do while(words(1)/='NUMER')
        call ecoute('neu_reanut')
     end do
    
     do while( words(1) /= 'ENDNU' )
        if( words(1) == 'ALGEB' ) then
           !
           ! Algebraic solver
           !
           solve_sol => solve(1:)
           call reasol(1_ip)
           
        else if( words(1) == 'MAXIM' ) then
           !
           ! Maximum inner iterations
           !
           miinn_neu = int(param(1))

        else if( words(1) == 'CONVE' ) then
           !
           ! Tolerance inner iterations
           !
           cotol_neu = param(1)

        else if( words(1) == 'RELAX' ) then
           !
           ! Relaxation parameter
           !
           relax_neu = param(1)

        else if( words(1) == 'NITSC' ) then
           !
           ! Nitsche coefficient
           !
           nitsche_neu = param(1)
           if( exists('SMOOT') ) kfl_smobo_neu = 1  !nitsche b.c will be smoothed

        end if
        !
        ! Read ADR data
        !
        call ADR_read_data(ADR_neu) 

        call ecoute('neu_reanut')

     end do

  end if

end subroutine neu_reanut
 
