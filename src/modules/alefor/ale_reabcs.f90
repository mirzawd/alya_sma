!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine ale_reabcs()
  !------------------------------------------------------------------------
  !****f* Alefor/ale_reabcs
  ! NAME 
  !    ale_reabcs
  ! DESCRIPTION
  !    This routine reads the boundary conditions 
  ! USES
  !    ale_bcntoe
  !    ecoute
  !    memchk
  !    runend
  !    ale_bounod
  !    ale_autbcs
  ! USED BY
  !    ale_turnon
  !***
  !------------------------------------------------------------------------
  use def_parame
  use def_inpout
  use def_master
  use def_domain
  use def_alefor 
  use mod_memchk
  use mod_opebcs
  use mod_ecoute, only :  ecoute
  implicit none
  integer(ip) :: dummi
  !
  ! Allocate memory
  !
  if( kfl_icodn > 0 ) then
     call opnbcs(1_ip,1_ip,dummi,dummi,tncod_ale) ! Memory for structure
     call opnbcs(2_ip,1_ip,ndime, 0_ip,tncod_ale) ! Memory for variable
 end if
  if( kfl_icodb > 0 ) then
     call opbbcs(0_ip,1_ip,1_ip,tbcod_ale)      
  end if
  if( kfl_geome > 0 ) then
     call opnbcs(0_ip,1_ip,ndime, 0_ip,tgcod_ale)
  end if

  if( INOTSLAVE ) then
     !
     ! Reach the nodal-wise section
     !
     call ecoute('ale_reabcs')
     do while(words(1)/='BOUND')
        call ecoute('ale_reabcs')
     end do
     !
     ! Loop over nodes and or boundaries
     !
     call ecoute('ale_reabcs')
     do while(words(1)/='ENDBO')

        if( words(1) == 'CODES' .and. exists('NODES') ) then

           !-------------------------------------------------------------
           !
           ! User-defined codes on nodes
           !
           !-------------------------------------------------------------

           if(exists('GEOME')) then
              !
              ! Geometrical node code
              !              
              tgcod => tgcod_ale(1:)
              call reacod(4_ip)

           else
              !
              ! Velocity: node codes
              !
              tncod => tncod_ale(1:)
              call reacod(1_ip)

           end if

        else if( words(1) == 'CODES' .and. exists('BOUND') ) then

           !-------------------------------------------------------------
           !
           ! User-defined codes on boundaries
           !          
           !-------------------------------------------------------------

           tbcod => tbcod_ale(1:)
           call reacod(2_ip)

        else if( words(1) == 'FUNCT' ) then

           !-------------------------------------------------------------
           !
           ! Function definitions
           !
           !-------------------------------------------------------------

           call runend('ALRE_REABCS: OBSOLETE OPTION')

        end if
        call ecoute('ale_reabcs')
     end do

  end if

end subroutine ale_reabcs
