!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    mod_tem_memory.f90
!> @author  Adrià Quintanas-Corominas
!> @date    2020-04-20
!> @brief   Solidz main module
!> @details Manage memory of solidz, etc.
!-----------------------------------------------------------------------

module mod_solidz

   use def_kintyp, only :  ip, rp, lg
   use def_master
   use def_solidz
   implicit none
   private

   public :: solidz_main
   
   contains


      !-----------------------------------------------------------------------
      !> 
      !> @author  Adrià Quintanas-Corominas
      !> @date    2020-04-20
      !> @brief   IMain subroutine
      !> @details Main subroutine
      !> 
      !-----------------------------------------------------------------------
   
      subroutine solidz_main(order)
         implicit none
         integer(ip), intent(in) :: order
         call Solidz(order)
      end subroutine solidz_main

end module mod_solidz
