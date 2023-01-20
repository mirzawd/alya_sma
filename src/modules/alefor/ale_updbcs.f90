!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Alefor
!> @{
!> @file    ale_updbcs.f90
!> @author  houzeaux
!> @date    2020-04-09
!> @brief   Update boundary conditions
!> @details This routine sets up the boundary condition for the mesh 
!>          displacement.
!>          The boundary conditions is the increment in displacement
!>          Delta d = d^{n+1}-d^n. Fixity code 1 should not be processed 
!>          as the boundary condition could be provided by a coupling
!> @} 
!-----------------------------------------------------------------------

subroutine ale_updbcs()

  use def_master
  use def_domain
  use def_alefor
  use def_kermod
  use mod_ker_functions, only : ker_functions
  use mod_local_basis,   only : local_basis_global_to_local
  use mod_local_basis,   only : local_basis_local_to_global
  implicit none

  integer(ip) :: ipoin,ifunc,itype,idime
  real(rp)    :: dispm_new(ndime),dispm_old(ndime)
  real(rp)    :: dispm_delta(3)
  !
  ! Non-constant boundary conditions
  !

  do ipoin = 1,npoin
     
     ifunc       = kfl_funno_ale(ipoin)
     itype       = kfl_funtn_ale(ipoin)

     if( itype /= 0 ) then    
        !
        ! Displacement function
        !
        call ker_functions(ipoin,ifunc,itype,bvess_ale(:,ipoin,2),dispm_new)              ! d^{n+1}
        !
        ! Old displacement
        !
        dispm_old(1:ndime)   = coord(1:ndime,ipoin) - coord_ori(1:ndime,ipoin)            ! d^n
        if( associated(kfl_fixrs_ale) ) &
             call local_basis_global_to_local(kfl_fixrs_ale(ipoin),dispm_old,NODE=ipoin)  ! Rotate d^n
        dispm_delta(1:ndime) = dispm_new(1:ndime)   - dispm_old(1:ndime)                  ! Delta d
        do idime = 1,ndime
           if( kfl_fixno_ale(idime,ipoin) == 1 ) then
              !
              ! Total displacement is prescribed
              !
              bvess_ale(idime,ipoin,1) = dispm_delta(idime)
           else if( kfl_fixno_ale(idime,ipoin) == 4 ) then
              !
              ! Coordinates are prescribed
              !
              call ker_functions(ipoin,ifunc,itype,bvess_ale(:,ipoin,2),dispm_new,COORDINATES=coord_ori(:,ipoin)) ! d^{n+1}
              bvess_ale(idime,ipoin,1) = dispm_new(idime) - coord(idime,ipoin)
           end if
        end do 
     else
        !
        ! Constant velocity is prescribed
        !
        do idime = 1,ndime
           if( kfl_fixno_ale(idime,ipoin) == 2 ) then
              dispm_delta(idime)       = bvess_ale(idime,ipoin,2) * cutim
              bvess_ale(idime,ipoin,1) = dispm_delta(idime)
           end if
        end do
     end if

  end do
  !
  ! Impose displacement to prescribed value
  !
  call local_basis_global_to_local(kfl_fixrs_ale,dispm)
  do ipoin = 1,npoin
     do idime = 1,ndime
        if( kfl_fixno_ale(idime,ipoin) > 0 ) then
           dispm(idime,ipoin,1) = bvess_ale(idime,ipoin,1)
        end if
     end do
  end do
  call local_basis_local_to_global(kfl_fixrs_ale,dispm)

end subroutine ale_updbcs
