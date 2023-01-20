!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    sld_wetno.f90
!> @date    18/05/2016
!> @author  Alfonso Santiago
!> @brief   Postprocess wet nodes
!> @details Postprocess wet nodes
!> @}
!------------------------------------------------------------------------
subroutine sld_wetno(touched)

  use def_kintyp,    only : ip,rp
  use def_master,    only : TIME_N, ITER_AUX, ITER_K, inotmaster, current_code
  use def_domain,    only : npoin
  use def_coupli,    only :  coupling_type,mcoup

  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  implicit none
  integer(ip)                       :: icoup

  logical(ip)                       :: code_j
  integer(ip)                       :: n_wets
  real(rp), intent(inout)           :: touched(npoin)
  integer(ip), pointer              :: wets(:) => null()
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!

  touched = 0.0_rp
  n_wets  = -1
  !
  do icoup=1_ip,  mcoup
    code_j = current_code == coupling_type(icoup) % code_target
    if(code_j) then
        n_wets =  coupling_type(icoup)%wet%npoin_wet
        wets   => coupling_type(icoup)%wet%lpoin_wet(:)
        if(inotmaster .and. associated(wets)) touched( wets(1:n_wets) ) = real(icoup,rp)
    endif
  enddo


end subroutine sld_wetno

