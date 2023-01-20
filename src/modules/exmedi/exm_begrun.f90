!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Exmper
!> @{
!> @file    exm_begrun.f90
!> @date    14/06/2019
!> @author  Guillaume Houzeaux
!> @brief   Beginning the run... 
!> @details Beginning the run... we can compute matrices!
!> @}
!-----------------------------------------------------------------------

subroutine exm_begrun()

  use def_master
  use def_exmedi
  use def_solver
  use def_domain
  use mod_memory
  implicit none
  !
  ! Additional arrays
  !
  call exm_addarr()
  !
  ! auxiliary system matrix
  !
  call memory_alloca(mem_modul(1:2,modul),'AMATR_AUXI_EXM','exm_memall',amatr_auxi_exm,memory_size(amatr))     
  !
  ! Compute diagonal indices
  !
  if( INOTEMPTY ) then
     if( solve_sol(1)%kfl_symme == 0 ) then
        call diagoi(npoin,1_ip,solve_sol(1) % kfl_symme,r_dom,c_dom,idima_exm)
     else
        call diagoi(npoin,1_ip,solve_sol(1) % kfl_symme,r_sym,c_sym,idima_exm)
     end if
  end if

end subroutine exm_begrun
 


