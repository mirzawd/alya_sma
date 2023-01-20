!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine setgts(order)
  !-----------------------------------------------------------------------
  !****f* master/setgts
  ! NAME
  !    setgts
  ! DESCRIPTION
  !    This routine computes the time step
  !    ITASK     = 1 ... Initialization (Turnon)
  !              = 2 ... On each time step (Timste)
  !
  !    KFL_TIMCO = 0 ... Prescribed time step
  !              = 1 ... From critical time step
  !              = 2 ... local/per module time step
  ! USES
  ! USED BY
  !    Turnon, Timste
  !***
  !-----------------------------------------------------------------------
  
  use def_master
  use mod_outfor,   only : outfor
  use def_kermod,   only : reset_factor
  use def_kermod,   only : kfl_reset
  use def_kermod,   only : kfl_reset_steps
  use def_kermod,   only : ittim_reset
  use def_kermod,   only : itti2_reset
  use def_kermod,   only : cutim_reset
  use def_kermod,   only : dtime_reset
  use def_kermod,   only : dtinv_reset
  use mod_messages, only : livinf
  
  implicit none
  
  integer(ip), intent(in) :: order
  integer(ip)             :: ii
  
  select case (order)

  case ( ITASK_INITIA )
     
     !-------------------------------------------------------------------
     !
     ! Initialize time variables
     !
     !-------------------------------------------------------------------
     !
     ! Compute DTINV and the time for the first time step, as well as the
     ! maximum number of time steps allowed
     !
     dtinv = 0.0_rp
     cutim = timei
     !
     ! Only execute if time is prescribed
     !
     if( kfl_timco == 0 ) call timfun()

     if( dtime >  0.0_rp ) dtinv = 1.0_rp / dtime

     dtold(1)     = dtime
     dtinv_old(1) = 1.0_rp / max(dtold(1),zeror)
     oltim        = 0.0_rp

  case ( ITASK_TIMSTE )

     !-------------------------------------------------------------------
     !
     ! Actualize time variables
     !
     !-------------------------------------------------------------------
     !
     ! Use minimum of critical time steps
     !
     dtold(1) = dtime
     !
     ! Only execute if time is prescribed
     !
     if( kfl_timco == 0 ) call timfun()
     !
     ! From critical time step
     !
     if( kfl_timco == 1 ) then
        if( dtinv /= 0.0_rp ) dtime = 1.0_rp / dtinv
     end if

     where( dtold == 0.0_rp ) dtold = dtime
     
     dtinv_old  = 1.0_rp / max(dtold,zeror)
     oltim      = cutim
     cutim      = cutim + dtime
     ioutp(1)   = ittim
     routp(1)   = dtime
     routp(2)   = cutim
     call outfor(15_ip,lun_outpu,' ')

  case ( ITASK_ENDSTE )

     !-------------------------------------------------------------------
     !
     ! Ending a time step: save old values
     !
     !-------------------------------------------------------------------

     do ii = size(dtold),2,-1
        dtold(ii) = dtold(ii-1)
     end do

  case ( ITASK_SAVE_RESET )
     
     !-------------------------------------------------------------------
     !
     ! Save for reset
     !
     !-------------------------------------------------------------------

     ittim_reset = ittim
     itti2_reset = itti2
     cutim_reset = cutim
     dtime_reset = dtime
     dtinv_reset = dtinv

  case ( ITASK_RECOVER_RESET )
     
     !-------------------------------------------------------------------
     !
     ! Reset and go back in time
     !
     !-------------------------------------------------------------------

     if( kfl_reset == 1 ) then
        !
        ! Recover old values
        !
        ittim = ittim_reset 
        itti2 = itti2_reset 
        cutim = cutim_reset
        dtime = dtime_reset
        dtinv = dtinv_reset
        !
        ! Modify time steps
        !
        dtinv = dtinv / reset_factor
        dtime = dtime * reset_factor
        !
        ! Recover old time step sizes
        !
        do ii = 1,size(dtold)-kfl_reset_steps
           dtold(ii) = dtold(ii+kfl_reset_steps-1)
        end do
        !
        ! Inverse time
        !
        dtinv_old(1:10) = 1.0_rp / max(dtold(1:10),zeror)
        oltim           = cutim
       
        call livinf(201_ip, ' ',1_ip)
     end if

  end select

end subroutine setgts
