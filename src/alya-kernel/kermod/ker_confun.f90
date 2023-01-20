!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Kermod
!> @{
!> @file    ker_confun.f90
!> @author  Guillaume Houzeaux
!> @date    18/03/2019
!> @brief   Define concentration
!> @details Define concentration
!> @} 
!-----------------------------------------------------------------------
subroutine ker_confun(itask)
  use def_parame
  use def_master
  use def_kermod
  use def_domain
  use mod_memory
  use mod_ker_space_time_function, only : ker_space_time_function

  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: ipoin,iclas,ifunc,ifiel

  !----------------------------------------------------------------
  !
  ! ADVEC Computed only for user defined functions 
  !
  !----------------------------------------------------------------

  if( INOTMASTER .and. kfl_cofun /= 0 ) then

     if( itask == ITASK_BEGSTE .or. itask == ITASK_INIUNK ) then        

        if( kfl_cofun < 0 ) then
           !
           ! From Field
           !
           ifiel = -kfl_cofun
           do ipoin = 1,npoin
              do iclas = 1,kfl_field(1,ifiel)
                 conce(ipoin,iclas,1) = xfiel(ifiel) % a(iclas,ipoin,1) 
              end do
           end do              
           
        elseif( kfl_cofun > 1000 ) then
           
           do ipoin = 1,npoin
              ifunc = kfl_cofun  - 1000     
              call ker_space_time_function(&
                   ifunc,coord(1,ipoin),coord(2,ipoin),coord(ndime,ipoin),cutim,conce(ipoin,:,1))
           end do 

        else

           select case ( kfl_cofun )

           case ( 3_ip ) 
              
              do ipoin = 1,npoin
                 conce(ipoin,:,1) = 0.0_rp
              end do

           case ( 666_ip)
              !
              ! Relative humidity
              !
              do ipoin = 1,npoin
                 conce(ipoin,:,1) = conce_relhu
              end do
 
           case ( 99_ip ) 
              !
              ! Do not do anything
              !
              continue

           case default
              
              call runend('KER_COFUN: NOT CODED')

           end select
           
        end if
        !
        ! Assume constant initial concentration
        !
        if( itask == ITASK_INIUNK ) then
           do ipoin = 1,npoin 
              conce(ipoin,:,2) = conce(ipoin,:,1)
              conce(ipoin,:,3) = conce(ipoin,:,1)
           end do
        else if( itask == ITASK_BEGSTE ) then
           do ipoin = 1,npoin 
              conce(ipoin,:,2) = conce(ipoin,:,1)
           end do
        end if

     else if( itask == ITASK_ENDSTE ) then

        !----------------------------------------------------------------
        !
        ! Save previous concentration
        ! KFL_COFUN = 0, CONCE point to CONCE which should not be modified
        !
        !----------------------------------------------------------------

        if( kfl_cofun /= 0 ) then
           do ipoin = 1,npoin 
              conce(ipoin,:,3) = conce(ipoin,:,1)
           end do
        end if

     end if

  end if

end subroutine ker_confun
