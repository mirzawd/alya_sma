!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tur_reaous()
  !-----------------------------------------------------------------------
  !****f* Turbul/tur_reaous
  ! NAME 
  !    tur_reaous
  ! DESCRIPTION
  !    This routine reads the output strategy for the turbulence
  !    equation.
  ! USES
  ! USED BY
  !    tur_turnon
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_inpout
  use def_master
  use def_turbul
  use def_domain
  use mod_memchk
  use mod_ecoute, only :  ecoute
  use mod_output_postprocess, only : output_postprocess_read
  implicit none
 
  if( INOTSLAVE ) then
     !
     ! Initializations
     !
     kfl_exacs_tur = 0           ! Manufactured solution
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
              ! Manufactured solution
              !
              kfl_exacs_tur=getint('SOLUT',1_ip,'#Exact solution')
           end if
        end if

     end do

  end if


end subroutine tur_reaous

    
