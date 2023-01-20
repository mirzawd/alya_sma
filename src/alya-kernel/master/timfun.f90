!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine timfun()
  !------------------------------------------------------------------------
  !****f* master/timfun
  ! NAME
  !    timfun
  ! DESCRIPTION
  !    This routine computes the time step
  ! OUTPUT
  ! USES
  ! USED BY
  !    setgts
  !***
  !------------------------------------------------------------------------
  use def_master
  use mod_parall
  use mod_communications, only : PAR_BROADCAST
  implicit none
  integer(ip) :: idata
  real(rp)    :: tol

  if(kfl_timef.gt.0) then
     !
     ! User defined time function
     !
     ! THIS IS REALLY OLD AND SHOULD BE DEPRECATED
     !

     !
     !USER DEFINED TIME FUNCTION
     !
     select case(kfl_timef)

     case(1)
        !
        ! Marek
        !
        if(cutim<535.0_rp) then
           dtime=(397.91_rp*cutim+2435.9_rp)/10000.0_rp
        else
           dtime=(0.1513_rp*(cutim*cutim)-381.18_rp*cutim+367561.0_rp)/10000.0_rp
        end if

     case(2)
        !
        ! Marek
        !
        if(cutim<375.84_rp) then
           dtime=(397.91_rp*cutim+2435.9_rp)/10000.0_rp
        else
           dtime=15.0_rp
        end if

     case(3)
        !
        ! Temporal test for catamaran, actually it would me much flexible to give a table in the .dat file
        !
        if(cutim<(80.0_rp*0.005_rp)) then
           dtime=0.005_rp
        else
           dtime=0.02_rp
        end if

     case(4)
        !
        ! Increasing time step for Planet Boundary Layer Diurnal Cycle 
        ! From 1s to 5s time step using approx. 800 time steps
        if (ittim<=100) then
           dtime=0.2_rp
        else if (ittim<=3000) then
           dtime=0.2_rp+(0.002_rp*real((ittim-100_ip),rp))
        else 
           dtime=6.0_rp
        end if
        
     case(5)
        !
        ! Increasing time step for Planet Boundary Layer Diurnal Cycle
        ! From 1s to 7.5s time step using apporx. 3000 time steps
        ! At 500 steps -> 1.5 s  / at 1000 steps -> 4 s 
        dtime= 3.45_rp+4.0_rp*TANH(((cutim)/900.0_rp)-1.0_rp)
        if (dtime>6.0_rp) dtime = 6.0_rp
     end select

     if(dtime>0.0_rp) dtinv = 1.0_rp/dtime

  else if( kfl_timef < 0 ) then
     !
     ! Discrete time evolution
     !
     if( INOTSLAVE ) then

        if( kfl_dtfun /= 1 ) then
           call runend('kernel/master/timefun: DISCRETE TIME FUNCTION WAS REQUESTED BUT THERE IS NO DISCRETE TIME FUNCTION DEFINED')
        end if
        !
        ! Read data table
        !
        tol = 1.0E-15_rp ! Double precision
        do idata = 1, nfundt
           if ( cutim < dtfun(1,idata) - tol ) then
              dtime = dtfun(2,idata)              
              exit
           end if
        end do
        if ( dtime /= 0.0_rp ) dtinv = 1.0_rp / dtime

     endif
     call PAR_BROADCAST(dtime,'IN MY CODE')
     call PAR_BROADCAST(dtinv,'IN MY CODE')

  end if

end subroutine timfun
