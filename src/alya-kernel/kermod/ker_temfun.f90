!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Kermod
!> @{
!> @file    ker_temfun.f90
!> @author  Guillaume Houzeaux
!> @date    18/03/2019
!> @brief   Define advection temperature
!> @details Define the advection TEMPE
!> @} 
!-----------------------------------------------------------------------
subroutine ker_temfun(itask)
  use def_parame
  use def_master
  use def_kermod
  use def_domain
  use mod_memory
  use mod_chktyp,                  only : check_type
  use mod_ker_functions,           only : ker_functions
  use mod_ker_space_time_function, only : ker_space_time_function

  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: ipoin,ifunc

  !----------------------------------------------------------------
  !
  ! ADVEC Computed only for user defined functions
  !
  !----------------------------------------------------------------

  if( INOTMASTER .and. kfl_tefun > 0 ) then

     if( itask == ITASK_BEGSTE .or. itask == ITASK_INIUNK ) then        

        if( kfl_tefun > 1000 .and. kfl_tefun < 6000 ) then
           !
           ! Space and time functions
           !
           do ipoin = 1,npoin
              ifunc = kfl_tefun - 1000
              call ker_space_time_function(&
                   ifunc,coord(1,ipoin),coord(2,ipoin),coord(ndime,ipoin),cutim,tempe(ipoin,1))
              therm(ipoin,1) = tempe(ipoin,1)
           end do
           
        else if( kfl_tefun > 6000 ) then
           !
           ! Discrete functions
           !
           do ipoin = 1,npoin
              ifunc = kfl_tefun - 6000
              call ker_functions(ipoin,ifunc,6_ip,1.0_rp,therm(ipoin,1),WHEREIN='ON NODES')
           end do

        else
           ! 
           ! User values
           !
           select case ( kfl_tefun )

           case ( 1_ip ) 
              do ipoin = 1,npoin
                 tempe(ipoin,1) = 273.15_rp
                 therm(ipoin,1) = tempe(ipoin,1)
              end do

           case ( 3_ip ) 
              
              do ipoin = 1,npoin
                 tempe(ipoin,1) = 0.0_rp
                 therm(ipoin,1) = tempe(ipoin,1)
              end do
              
           case ( 99_ip ) 
              !
              ! Do not do anything
              !
              continue

           case default
              
              call runend('KER_TEFUN: NOT CODED')

           end select
           
        end if
        !
        ! Assume constant initial advection
        !
        if( itask == ITASK_INIUNK ) then
           do ipoin = 1,npoin 
              tempe(ipoin,2) = tempe(ipoin,1)
              tempe(ipoin,3) = tempe(ipoin,1)
              therm(ipoin,2) = therm(ipoin,1)
              therm(ipoin,3) = therm(ipoin,1)
           end do
        else if( itask == ITASK_BEGSTE ) then
           do ipoin = 1,npoin 
              tempe(ipoin,2) = tempe(ipoin,1)
              therm(ipoin,2) = therm(ipoin,1)
           end do
        end if

     else if( itask == ITASK_ENDSTE ) then

        !----------------------------------------------------------------
        !
        ! Save previous advection
        ! KFL_TEFUN = 0, TEMPE point to TEMPE which should not be modified
        !
        !----------------------------------------------------------------

        if( kfl_tefun /= 0 ) then
           do ipoin = 1,npoin 
              tempe(ipoin,3) = tempe(ipoin,1)
              therm(ipoin,3) = therm(ipoin,1)
           end do
        end if

     end if

  end if

end subroutine ker_temfun
