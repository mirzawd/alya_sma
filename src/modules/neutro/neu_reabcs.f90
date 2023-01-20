!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup NeutroInput
!> @{
!> @file    nsi_reabcs.f90
!> @author  Guillaume Houzeaux
!> @brief   Read boundary conditions 
!> @details Read boundary conditions, initial conditions and parameters
!> @} 
!-----------------------------------------------------------------------
subroutine neu_reabcs()
  use def_parame
  use def_inpout
  use def_master
  use def_domain
  use def_neutro
  use mod_opebcs, only : opebcs_initialization_variable
  use mod_opebcs, only : opebcs_initialization_structure
  use mod_opebcs, only : cpybcs_boundaries
  use mod_ecoute, only : ecoute
  implicit none
  integer(ip) :: ienergy,idirection,iunkn,iunk1
  external :: reacod
  !
  ! Allocate memory
  !
  if( kfl_icodn > 0 ) then
     call opebcs_initialization_structure(nunkn_neu,tncod_neu)
     call opebcs_initialization_variable(1_ip,tncod_neu)
  end if
  if( kfl_icodb > 0 ) then
     call opebcs_initialization_structure(nunkn_neu,tbcod_neu)
     call opebcs_initialization_variable(1_ip,tbcod_neu)
   end if
 
   if( INOTSLAVE ) then
     !
     ! Initialization global variables
     !
     !
     ! Reach the boundary condition section
     !
     call ecoute('neu_reabcs')
     do while( words(1) /= 'BOUND' )
        call ecoute('neu_reabcs')
     end do
     call ecoute('neu_reabcs')

     do while( words(1) /= 'ENDBO' )

        if(      words(1) == 'CODES' .and. exists('NODES') ) then

           !-------------------------------------------------------------
           !
           ! User-defined codes on nodes
           !
           !-------------------------------------------------------------

           ienergy    = getint('ENERG',1_ip,'#Node code for current variable ENERGY')
           idirection = getint('DIREC',1_ip,'#Node code for current variable DIRECTION')
           
           if( ienergy    < 1 .or. ienergy    > num_energies_neu   ) call runend('NEU_REABCS: WRONG ENERGY NUMBER')
           if( idirection < 1 .or. idirection > num_directions_neu ) call runend('NEU_REABCS: WRONG DIRECTION NUMBER')

           iunkn =  (idirection-1)*num_energies_neu + ienergy
           tncod => tncod_neu(iunkn:)
           call reacod(READ_NODE_CODES)

        else if( words(1) == 'CODES' .and. exists('BOUND') ) then

           !-------------------------------------------------------------
           !
           ! User-defined codes on boundaries
           !          
           !-------------------------------------------------------------

           if( exists('ALL  ' ) ) then
              ienergy    =  1
              idirection =  1
              iunk1      =  (idirection-1)*num_energies_neu + ienergy
              tbcod      => tbcod_neu(iunk1:)
              call reacod(READ_BOUNDARY_CODES)
              do ienergy = 1,num_energies_neu
                 do idirection = 1,num_directions_neu
                    if( ienergy /= 1 .or. idirection /= 1 ) then
                       iunkn =  (idirection-1)*num_energies_neu + ienergy
                       call cpybcs_boundaries(1_ip,iunkn,tbcod_neu)
                       !tbcod_neu(iunkn:) => tbcod_neu(iunk1:)
                    end if
                 end do
              end do
           else 
               ienergy    = getint('ENERG',1_ip,'#Node code for current variable ENERGY')
               idirection = getint('DIREC',1_ip,'#Node code for current variable DIRECTION')
               if( ienergy    < 1 .or. ienergy    > num_energies_neu   ) call runend('NEU_REABCS: WRONG ENERGY NUMBER')
               if( idirection < 1 .or. idirection > num_directions_neu ) call runend('NEU_REABCS: WRONG DIRECTION NUMBER')
               iunkn =  (idirection-1)*num_energies_neu + ienergy
               tbcod => tbcod_neu(iunkn:)
               call reacod(READ_BOUNDARY_CODES)
           end if

        end if

        call ecoute('neu_reabcs')

     end do

  end if


end subroutine neu_reabcs
