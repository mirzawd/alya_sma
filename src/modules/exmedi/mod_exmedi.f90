!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Exmedi
!> @{
!> @file    mod_exmedi.f90
!> @author  Alfonso Santiago (BSC/ELEM)
!> @date    2020-04-20
!> @brief   Exmedi main module
!> @details 
!-----------------------------------------------------------------------

module mod_exmedi

   use def_kintyp, only :  ip, rp, lg
   use def_master
   use def_exmedi
   implicit none
   private

   public :: exmedi_main
   
   contains


      !-----------------------------------------------------------------------
      !> 
      !> @author  Alfonso Santiago (BSC/ELEM)
      !> @date    2020-04-20
      !> @brief   Main subroutine
      !> @details Main subroutine
      !> 
      !-----------------------------------------------------------------------
   
      subroutine exmedi_main(order)
         implicit none
         integer(ip), intent(in) :: order
         call Exmedi(order)
      end subroutine exmedi_main

end module mod_exmedi
