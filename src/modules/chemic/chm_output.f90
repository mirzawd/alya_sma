!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine chm_output()
  !------------------------------------------------------------------------
  !****f* Temper/chm_output
  ! NAME
  !    chm_output
  ! DESCRIPTION
  !    End of a time step
  !
  ! ITASK = -1 ... Initial solution
  ! ITASK =  0 ... When timemarching is true. There is output or
  !                post-process of results if required.
  ! ITASK =  1 ... When timemarching is false. Output and/or
  !                post-process of results is forced if they have not
  !                been written previously.
  ! USES
  ! USED BY
  !    chm_iniunk
  !    chm_endite
  !    chm_endste
  !***
  !------------------------------------------------------------------------
  use def_master,             only : ITASK_ENDTIM, ITASK_INITIA, postp, cutim, ittyp, ittim
  use def_kintyp,             only : ip, nvarp
  use def_chemic,             only : avtim_chm
  use mod_output_postprocess, only : output_postprocess_variables
  use def_kermod,             only : kfl_chemic_vect
  implicit none
  integer(ip) :: ivarp

  external    :: chm_outvar
  external    :: chm_outset
  external    :: chm_outwit
  external    :: chm_post_scalar_dissipation_rate
  external    :: chm_post_scalar_dist

  !
  ! Calculate scalar dissipation rates for output on volume or witness points
  !
  if(kfl_chemic_vect /= 1_ip) then
      if( postp(1) % npp_witne(8+3) == 1 ) then
         call chm_post_scalar_dissipation_rate(23_ip) ! XYR
      endif
      if( postp(1) % npp_witne(8+1) == 1 ) then
         call chm_post_scalar_dissipation_rate(24_ip) ! XZR
      endif
      if( postp(1) % npp_witne(8+4) == 1 ) then
         call chm_post_scalar_dissipation_rate(25_ip) ! XYS
      endif
      if( postp(1) % npp_witne(8+2) == 1 ) then
         call chm_post_scalar_dissipation_rate(26_ip) ! XZS
      endif
  else
      if( postp(1) % npp_witne(8+3) == 1 ) then
         call chm_post_scalar_dist(23_ip) ! XYR
      endif
      if( postp(1) % npp_witne(8+1) == 1 ) then
         call chm_post_scalar_dist(24_ip) ! XZR
      endif
      if( postp(1) % npp_witne(8+4) == 1 ) then
         call chm_post_scalar_dist(25_ip) ! XYS
      endif
      if( postp(1) % npp_witne(8+2) == 1 ) then
         call chm_post_scalar_dist(26_ip) ! XZS
      endif
  end if


  !
  ! Initial solution, end of a time step and and of run
  !
  call output_postprocess_variables(chm_outvar)
  !
  ! Update reference time for time-averaging
  !
  do ivarp = 1,nvarp
     if( postp(1) % npp_stepi(ivarp,0) /= 0 ) then
        if( mod(ittim, postp(1) % npp_stepi(ivarp,0) ) == 0 ) then
           avtim_chm = cutim
        endif
     endif
  end do


  if( ittyp == ITASK_INITIA .or. ittyp == ITASK_ENDTIM ) then
     !
     ! Calculations on sets
     !
     call chm_outset()
     !
     ! Calculations on witness points
     !
     call chm_outwit()

  end if

end subroutine chm_output
