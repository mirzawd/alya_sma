!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Gusano
!> @{
!> @file    gus_updunk.f90
!> @author  houzeaux
!> @date    2018-12-28
!> @brief   This routine performs several types of updates for the
!>          thermal variable (T or h)
!> @details Solution updates
!> @}
!-----------------------------------------------------------------------

subroutine gus_updunk(itask)

  use def_parame
  use def_master
  use def_domain
  use def_gusano
  use mod_gus_projections, only : gus_projections_updunk
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: ipoin
!  integer(ip)             :: itime
  integer(ip)             :: nprev_gus,ncomp_gus
  
  ncomp_gus = 3
  nprev_gus = 3 !min(3_ip,ncomp_gus)
  
  select case ( itask )

  case ( ITASK_INIUNK )
     !
     ! (:,all) <= (:,1): Initial solution
     !
     do ipoin = 1,npoin
        press(ipoin,2:ncomp_gus) = press(ipoin,1)
        flowr(ipoin,2:ncomp_gus) = flowr(ipoin,1)
     end do
     
  case ( ITASK_BEGSTE )
     !
     ! (:,2) <= (:,1): Initial guess for outer iterations
     !
     do ipoin = 1,npoin
        press(ipoin,2) = press(ipoin,1)
        flowr(ipoin,2) = flowr(ipoin,1)
     end do

  case ( ITASK_BEGITE )
     !
     ! UNKNO <= (:,1)
     !
     if( solve(1) % kfl_block == 0 ) then
        do ipoin = 1,npoin
           unkno((ipoin-1)*2+1) = flowr(ipoin,1)
           unkno((ipoin-1)*2+2) = press(ipoin,1)
        end do
     else
        do ipoin = 1,npoin
           unkno(ipoin)       = flowr(ipoin,1)
           unkno(ipoin+npoin) = press(ipoin,1)
        end do        
     end if

  case ( ITASK_ENDINN )
     !
     ! (:,1) <= UNKNO
     !
     if( solve(1) % kfl_block == 0 ) then
        do ipoin = 1,npoin
           flowr(ipoin,1)     = unkno((ipoin-1)*2+1)      
           press(ipoin,1)     = unkno((ipoin-1)*2+2) 
           vel1d(ipoin,1)     = flowr(ipoin,1) / areas(ipoin,1)
        end do
     else
        do ipoin = 1,npoin
           flowr(ipoin,1)     = unkno(ipoin)      
           press(ipoin,1)     = unkno(ipoin+npoin) 
           vel1d(ipoin,1)     = flowr(ipoin,1) / areas(ipoin,1)
        end do        
     end if
     call gus_projections_updunk(ITASK_ENDINN)
     
  case ( ITASK_ENDITE )
     !
     ! (:,2) <= (:,1): End of inner iteration
     !        
     do ipoin = 1,npoin
        flowr(ipoin,2) = flowr(ipoin,1)
        press(ipoin,2) = press(ipoin,1)
     end do
     
  case ( ITASK_ENDSTE )
     !
     ! (:,3) <= (:,1): End of time step
     ! (:,4) <= (:,3)
     ! (:,5) <= (:,4)
     ! ...
     !        
     !if( kfl_tisch_gus == 2 ) then
     !   do itime = 2+kfl_tiaor_gus,4,-1
     !      do ipoin=1,npoin
     !         tempe(ipoin,itime) = tempe(ipoin,itime-1)
     !      end do
     !   end do
     !end if
     do ipoin = 1,npoin
        flowr(ipoin,nprev_gus) = flowr(ipoin,1)
        press(ipoin,nprev_gus) = press(ipoin,1)
     end do
             
  end select


end subroutine gus_updunk
