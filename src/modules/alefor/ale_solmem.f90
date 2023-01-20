!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Alefor
!> @{
!> @file    ale_memall.f90
!> @author  Guillaume Houzeaux
!> @date    16/11/1966
!> @brief   Arrays allocation subroutine
!> @details Arrays allocation subroutine
!> @} 
!-----------------------------------------------------------------------
subroutine ale_solmem
  
  use      def_parame
  use      def_inpout
  use      def_master
  use      def_domain
  use      def_solver
  use      def_alefor
  use      mod_memory
  use      mod_ale_arrays
  
  implicit none
  
  !
  ! Solver memory
  !
  solve_sol => solve(1:1)
  call soldef(4_ip)  
  solve(1) % bvess     => bvess_ale(:,:,1)
  solve(1) % kfl_fixno => kfl_fixno_ale
  solve(1) % kfl_fixrs => kfl_fixrs_ale

end subroutine ale_solmem
      
