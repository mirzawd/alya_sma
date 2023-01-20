!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine read_module_options()
  !------------------------------------------------------------------------
  !****f* kernel/read_module_options
  ! NAME 
  !    read_module_options
  ! DESCRIPTION
  !    This routine reads the generic module data
  ! USES
  ! USED BY
  !    ***_reapro
  !***
  !------------------------------------------------------------------------
  use def_master
  use def_inpout
  use mod_ecoute, only :  ecoute
  implicit none
  character(10) :: messa

  if( INOTSLAVE ) then
     ! 
     ! Reach the section
     !      
     messa=exmod(modul)//'_REAPRO'
     rewind(lisda)
     do while(words(1)/='RUNDA')
        call ecoute(messa)
     end do
     do while(words(1)/='ENDRU')
        call ecoute(messa)
     end do
     do while(words(1)/='PROBL')
        call ecoute(messa)
     end do
     !
     ! Read data
     !      
     do while(words(1)/='ENDPR')
        call ecoute(messa)
        if(words(1)==namod(modul)(1:5)) then   
           if(exists('ON   ').or.exists('ACTIV')) then
              kfl_modul(modul) = 1
              do while(words(1)/='END'//namod(modul)(1:2))
                 call ecoute(messa)
                 if(words(1)=='DELAY') then
                    kfl_delay(modul) = 1
                    ndela(modul)     = int(param(2),ip)
                 else if(words(1)=='CONVE') then
                    if(exists('YES  ').or. exists('ON   ')) kfl_conve(modul) = 1
                    if(exists('NO   ').or. exists('OFF  ')) kfl_conve(modul) = 0
                 else if(words(1)=='SOLVE') then
                    if( words(2) == 'ATEAC' ) kfl_solve(modul) = AT_EACH_TIME_STEP
                    if( words(2) == 'ATBEG' ) kfl_solve(modul) = AT_BEGINNING
                 end if
              end do
           end if
        end if
     end do
  end if

end subroutine read_module_options
