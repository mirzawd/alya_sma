!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Temper
!> @{
!> @file    tem_coupli.f90
!> @author  Guillaume Houzeaux
!> @brief   Coupling of tmper with other modules
!> @details Coupling of tmper with other modules\n
!>          - Coupling with CHEMIC: water vapor model\n
!>            qL = w/Area * Dh \n
!>            w/Area = rho_wv k_wv * gradc_wv.n \n
!> @} 
!-----------------------------------------------------------------------
subroutine tem_coupli(itask)
  
  use def_master
  use def_domain
  use def_elmtyp
  use def_temper
  use mod_gradie
  
  implicit none
  
  integer(ip), intent(in) :: itask

  select case ( itask )

  case ( ITASK_CONCOU )


  case ( ITASK_BEGSTE )


  case ( ITASK_INIUNK )


  case ( ITASK_BEGITE )

     if( kfl_coupl(ID_TEMPER,ID_CHEMIC) ==2 ) then
        !
        ! Water vapor concentration gradients
        !
        if( INOTMASTER ) then
           call gradie(conce(1:npoin,1,1),gradc_tem)
        end if

     end if

  case(ITASK_DOITER) 


  end select

end subroutine tem_coupli
