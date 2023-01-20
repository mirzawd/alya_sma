!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tur_averag()
  !------------------------------------------------------------------------
  !****f* Turbul/tur_averag
  ! NAME 
  !    tur_averag
  ! DESCRIPTION
  !    Average key, omega, turvi
  ! USES
  ! USED BY
  !    tur_endste 
  !***
  !------------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_turbul
  use mod_output_postprocess, only : output_postprocess_check_variable_postprocess
  use mod_arrays,             only : arrays_number
  implicit none
  integer(ip) :: ipoin
  real(rp)    :: xfact

  if( INOTMASTER ) then

     if( cutim > avtim_tur ) then
        xfact     = 0.5_rp*dtime 

        ! AVKEY: average key

        if( output_postprocess_check_variable_postprocess(arrays_number('AVKEY')) ) then
          do ipoin = 1,npoin
              avkey_tur(ipoin)=avkey_tur(ipoin)&
                   +(untur(1,ipoin,1)+untur(1,ipoin,3)) * xfact
           end do
        end if

        ! AVOME: average omega

        if( output_postprocess_check_variable_postprocess(arrays_number('AVOME')) ) then
          do ipoin = 1,npoin
              avome_tur(ipoin)=avome_tur(ipoin)&
                   +(untur(2,ipoin,1)+untur(2,ipoin,3)) * xfact
           end do
        end if

        ! AVTVI: average turbulent viscosity

        if( output_postprocess_check_variable_postprocess(arrays_number('AVTVI')) ) then
          do ipoin = 1,npoin
              avtvi_tur(ipoin)=avtvi_tur(ipoin)&
                   +(turmu(ipoin)+olded_tur(ipoin)) * xfact

           end do
        end if

     end if

  end if

end subroutine tur_averag

