!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tur_doiter
!-----------------------------------------------------------------------
!****f* Turbul/tur_doiter
! NAME 
!    tur_doiter
! DESCRIPTION
!    This routine controls the internal loop of the temperature equation.
! USES
!    tur_begite
!    tur_solite
!    tur_endite
! USED BY
!    Turbul
!***
!-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_solver
  use def_turbul
  use mod_timings, only : timings_ini
  use mod_timings, only : timings_end
  implicit none
  
  
  if( kfl_stead_tur == 0 ) then
     call timings_ini()
     call tur_begite()
     call timings_end(ITASK_BEGITE)

     do while( kfl_goite_tur == 1 )
        if( kfl_algor_tur == 1 ) then
           do iunkn_tur = 1,nturb_tur
              do itera_tur =1,niter_tur
                 call tur_solite()
                 call tur_endite(ITASK_ENDINN)
              end do              
              ! acualizes untur(:,2) with relaxation factor
              call tur_updunk(ITASK_DOITER)
           end do
        else
           call tur_solite()
           call tur_endite(ITASK_ENDINN)
        end if
     end do
     call timings_ini()
     call tur_endite(ITASK_ENDITE)
     call timings_end(ITASK_ENDITE)
  end if

end subroutine tur_doiter
