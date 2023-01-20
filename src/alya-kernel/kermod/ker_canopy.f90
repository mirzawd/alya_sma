!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine ker_canopy()
  !-----------------------------------------------------------------------
  !****f* kermod/ker_canopy
  ! NAME
  !    ker_canopy
  ! DESCRIPTION
  !   Calculate canopy height & height over terrain
  ! OUTPUT
  !    canopy height & height over terrain
  ! USED BY
  !    ker_iniunk
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_kermod
  use mod_memchk
  use def_domain
  use mod_communications_global, only : PAR_SUM
  implicit none
  integer(ip) :: ierro

  if( kfl_canhe /= -1 ) then  !exists canopy

     ierro = 0

     if( INOTMASTER ) then

        if( kfl_canhe > 0 ) then   ! for the constant case it is done inside mod_ker_proper
           !
           ! Given by a field
           !
           if( nfiel < kfl_canhe ) then
              ierro =  1
           else
              canhe => xfiel(kfl_canhe) % a(1,:,1)
           end if

        end if

     end if
     !
     ! Check errors
     !
     call PAR_SUM(ierro)
     if( ierro /= 0 ) then
        call runend('FIELD FOR CANOPY HEIGHT HAS NOT BEEN DEFINED')
     end if
 

  end if


  if( kfl_heiov /= -1 ) then ! heiov not defined (error)

     ierro = 0

     if( INOTMASTER ) then

        if( kfl_heiov == 0 ) then ! flag when walld
           !
           ! height over ground is approximated by the wall distance - was what we were using initially - poor with slope
           !
           heiov => walld

        else if( kfl_heiov > 0 ) then
           !
           ! Given by a field
           !
           if( nfiel < kfl_heiov ) then
              ierro =  1
           else
              heiov => xfiel(kfl_heiov) % a(1,:,1)
           end if

        end if

     end if
     !
     ! Check errors
     !
     call PAR_SUM(ierro)
     if( ierro /= 0 ) then
        call runend('FIELD FOR HEIGHT OVER GROUND HAS NOT BEEN DEFINED')
     end if
 

  end if

  if( kfl_canla /= 0 ) then  ! leaf area density
     
     ierro = 0
     
     if( INOTMASTER ) then
        
           !
           ! Given by a field
           !
           if( nfiel < kfl_canla ) then
              ierro =  1
           else
              canla => xfiel(kfl_canla) % a(1,:,1)
           end if
           
        
     end if
     !
     ! Check errors
     !
     call PAR_SUM(ierro)
     if( ierro /= 0 ) then
        call runend('FIELD FOR CANOPY LEAF AREA DENSITY HAS NOT BEEN DEFINED')
     end if
 

  end if


end subroutine ker_canopy

