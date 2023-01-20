!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Chemic
!> @{
!> @file    chm_massma.f90
!> @date    23/10/2013
!> @author  Guillaume Houzeaux
!> @brief   Compute weighted mass matrix
!> @details Assemble the preconditioner diagonal coefficients
!>          Used for Richardson solvers:
!>          x^{i+1} = x^i + P.M^-1.r^i
!> @}
!-----------------------------------------------------------------------

subroutine chm_massma(pmatr)
  use def_kintyp,     only  :  ip,rp
  use def_solver,     only  :  solve_sol,SOL_LOCAL_DIAGONAL
  use def_master,     only  :  gesca,INOTMASTER
  use def_domain,     only  :  npoin
  use def_chemic,     only  :  dtinv_chm
  use mod_ker_proper, only  :  ker_proper
  implicit none
  real(rp),   intent(inout) :: pmatr(*)     !< Diagonal preconditioner P
  integer(ip)               :: ipoin,dummi

  external                  :: memgen

  if( solve_sol(1) % kfl_preco == SOL_LOCAL_DIAGONAL .and. INOTMASTER ) then

     call memgen(0_ip,npoin,0_ip)
     call ker_proper('DENSI','NPOIN',dummi,dummi,gesca)
     do ipoin = 1,npoin
        pmatr(ipoin) = 1.0_rp/(dtinv_chm*gesca(ipoin))
     end do
     call memgen(2_ip,npoin,0_ip)

  end if

end subroutine chm_massma
