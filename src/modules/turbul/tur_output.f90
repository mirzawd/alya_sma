!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tur_output()
  !-----------------------------------------------------------------------
  !****f* Turbul/tur_output
  ! NAME
  !   tur_output
  ! DESCRIPTION
  !   End of a temperature time step: 
  !   itask = 0  When timemarching is true. There is output or post-process
  !              of results if required.
  !   itask = 1  When timemarching is false. Output and/or post-process of
  !              results is forced if they have not been written previously.
  ! USES
  !    tur_outvar
  ! USED BY
  !    tur_endste
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_turbul
  use mod_iofile
  use mod_output_postprocess, only : output_postprocess_variables
  use mod_arrays,             only : arrays_number
  implicit none
  external    :: tur_outvar
!  integer(ip) :: ivari
!  real(rp)    :: time1,time2
  !
  ! Initial solution, end of a time step and and of run
  !
  ! call cputim(time1)
  call output_postprocess_variables(tur_outvar)  ! makes call tur_outvar(ivari,imesh)

! call cputim(time2)
! print*,'TUR_OUTPUT:',time2-time1

  if( postp(1) % npp_stepi(arrays_number('AVTVI'),0) /= 0 ) then
     if( mod(ittim, postp(1) % npp_stepi(arrays_number('AVTVI'),0) ) == 0 ) then              ! AVTVI_TUR
        avtim_tur = cutim  ! Update reference time for time-averaging
     endif
  endif

  if( ittyp == ITASK_INITIA .or. ittyp == ITASK_ENDTIM ) then
     !
     ! Postprocess on witness points
     !
     call tur_outwit()

  end if

end subroutine tur_output

