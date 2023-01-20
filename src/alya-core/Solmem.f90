!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Solmem
!> @{
!> @file    Solmem.f90
!> @author  houzeaux
!> @date    2019-06-17
!> @brief   Solver memory
!> @details Compute solver requirements
!> @} 
!-----------------------------------------------------------------------

subroutine Solmem()

  use def_kintyp,             only : ip
  use def_master,             only : iblok
  use def_master,             only : nblok
  use def_master,             only : modul
  use def_master,             only : mmodu
  use def_master,             only : momod
  use def_master,             only : kfl_modul
  use def_master,             only : ITASK_SOLMEM
  use def_domain,             only : nzmat,nzrhs
  use mod_output_postprocess, only : output_postprocess_allocate_sets_and_witness
  use mod_moduls,             only : moduls
  use mod_solver,             only : solver_matrix_RHS_sizes
  
  implicit none
  integer(ip) :: ivari,imodu

  call Kermod(ITASK_SOLMEM)
  do iblok = 1_ip,nblok
     call moduls(ITASK_SOLMEM)
  end do
  modul = 0_ip
  !
  ! Postprocess NOT GOOD, CHECK ORDER OF THINGS... put open file in turnon after moddef...
  !
  do imodu = 1_ip,mmodu
     if( kfl_modul(imodu) /= 0_ip ) then
        call output_postprocess_allocate_sets_and_witness(imodu,momod(imodu) % postp(1))
     end if
  end do
  !
  ! Matrix and RHS sizes
  !
  do imodu = 1_ip,mmodu
     if( kfl_modul(imodu) /= 0_ip ) then
        if( associated(momod(imodu) % solve) ) then
           do ivari = 1_ip,size(momod(imodu) % solve,KIND=ip)
              call solver_matrix_RHS_sizes(momod(imodu) % solve(ivari),nzmat,nzrhs)
           end do
        end if
     end if
  end do
  !
  ! Allocate memory for all unknowns of the problem and coefficients
  !
  call memunk(1_ip)

end subroutine Solmem
