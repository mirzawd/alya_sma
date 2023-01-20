!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tem_reaous()

  !-----------------------------------------------------------------------
  !
  ! This routine reads the output strategy for the temperature
  ! equation.
  !
  !-----------------------------------------------------------------------
  use def_parame
  use def_inpout
  use def_master
  use def_temper
  use def_domain
  use mod_ecoute, only :  ecoute
  use mod_output_postprocess, only : output_postprocess_read
  implicit none

  if( INOTSLAVE ) then
     !
     ! Initializations
     !
     kfl_splot_tem          = 0           ! Surface plot
     kfl_psmat_tem          = 0           ! Postscript file of the matrix
     kfl_exacs_tem          = 0           ! Exact solution
     kfl_viewf_tem          = 0           ! View factors
     npp_bound_tem          = 0           ! Postprocess boundary conditions 
     avtim_tem              = 0.0_rp      ! Averaging initial time 
     expar_tem              = 0.0_rp      ! Exact solution parameters
   
     !
     ! Reach the section
     !
     call ecoute('tem_reaous')
     do while(words(1)/='OUTPU')
        call ecoute('tem_reaous')
     end do
     !
     ! Begin to read data.
     !
     do while(words(1)/='ENDOU')
        call ecoute('tem_reaous')

        call output_postprocess_read()

        if(words(1)=='OUTPU') then
           !
           ! Output
           !
           if(words(2)=='ERROR') then
              !
              ! Exact solution
              !
              kfl_exacs_tem=getint('SOLUT',1_ip,'#Exact solution')
              expar_tem=param(4:3+nexap_tem)

           else if(words(2)=='MATRI') then
              !
              ! Matrix profile
              !
              if( exists('GID  ') ) then
                 kfl_psmat_tem = -getint('ITERA',1_ip,'#Iteration to postprocess matrix') 
              else
                 kfl_psmat_tem =  getint('ITERA',1_ip,'#Iteration to postprocess matrix') 
              end if

           else if(words(2)=='SOLVE') then
              !
              ! Solver convergence
              !
              solve(1)%kfl_cvgso=1

           end if

        else if(words(1)=='SURFA') then
           !
           ! Draw 3D plots
           !
           kfl_splot_tem = 1

        else if(words(1)=='AVERA') then
           !
           ! Averaging starting time
           !
           avtim_tem = getrea('AVERA',0.0_rp,'#Averaging starting time') ! for turbulence
        end if

     end do

  end if

end subroutine tem_reaous
    
