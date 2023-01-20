!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!




!----------------------------------------------------------------------------------------
!> @addtogroup Alya
!> @{
!> @file    Alya.f90
!> @author  houzeaux
!> @date    2019-09-20
!> @brief   Alya main
!> @details A L Y A \n
!>          COMPUTATIONAL MECHANICS AND DESIGN  \n
!>          \n
!>          Contact and general info: \n
!>           \n
!>          guillaume.houzeaux@bsc.es   \n
!>          mariano.vazquez@bcs.es  \n
!>          \verbatim
!>                                                           
!>                @@        &@       @@,      @@@    ,@%           
!>              /@@@@       &@        ,@@   *@@     @@@@@          
!>             @@(  @@.     &@          &@@@@      @@   @@         
!>            @@@@@@@@@%    &@            @@      @@@@@@@@@,       
!>           @@       ,@@   &@            @@    ,@@       &@@      
!>         *@@          @@  &@@@@@@@@@,   @@   &@@         .@@    
!>
!>          \endverbatim
!> @} 
!----------------------------------------------------------------------------------------

program Alya
  
  use def_kintyp,         only : ip
  use def_master,         only : ITASK_INITIA
  use def_master,         only : ITASK_ENDTIM
  use def_master,         only : ITASK_ENDRUN
  use def_master,         only : ITASK_READ_RESTART
  use def_master,         only : ITASK_WRITE_RESTART
  use def_master,         only : kfl_gotim
  use def_master,         only : kfl_stop
  use def_master,         only : kfl_goblk
  use def_master,         only : kfl_gocou
  use def_coupli,         only : kfl_gozon
  use def_kermod,         only : kfl_reset
  use mod_AMR,            only : AMR
  implicit none

  call Initia()                                               ! Initialization of the run
  call Readom()                                               ! Domain reading
  call Partit()                                               ! Domain partitioning
  call Reaker()                                               ! Read Kermod
  call Domtra()                                               ! Domain transformation
  call Domain()                                               ! Domain construction
  call Turnon()                                               ! Read modules
  call Solmem()                                               ! Solver and output memory
  call Begrun()                                               ! Initial computations

  call Restar(ITASK_READ_RESTART)                             ! Read restart
  call Iniunk()                                               ! Initial solution
  call Output(ITASK_INITIA)                                   ! Initial output
 
  time: do while ( kfl_gotim == 1_ip .and. kfl_stop == 0_ip )

      call Timste()                                           ! Compute time step

      reset: do
        call Begste()                                         ! Begin time step

         block: do while ( kfl_goblk == 1_ip .and. kfl_stop == 0_ip )

           zone_coupling: do while ( kfl_gozon == 1_ip .and. kfl_stop == 0_ip )

              call Begzon()                                   ! Start multicode block coupling

              coupling_modules: do while ( kfl_gocou == 1_ip .and. kfl_stop == 0_ip .and. kfl_reset /= 1_ip )
                 call Doiter()                                ! Inner loop of modules
                 call Concou()                                ! Check coupling convergence
              end do coupling_modules

              call Endzon()

           end do zone_coupling

           call Conblk()                                      ! End multicode block coupling

        end do block
        if( kfl_reset /= 1_ip .or. kfl_stop /= 0_ip ) exit reset

     enddo reset

     call Endste()                                            ! End of time step

     call Output(ITASK_ENDTIM)
     call EndTim(ITASK_ENDTIM)

     call Restar(ITASK_WRITE_RESTART)                         ! Write restart
     call AMR()                                               ! Adaptive mesh refinement
     call Repart()                                            ! Repartitioning
     
  end do time
  
  call Output(ITASK_ENDRUN)

  call Turnof()                                               ! Finish Alya
  
end program Alya
