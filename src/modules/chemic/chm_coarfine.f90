!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



  !----------------------------------------------------------------------
  !> @addtogroup Nastin
  !> @{
  !> @file    chm_coarfin.f90
  !> @author  J.C. Cajas
  !> @date    09/04/2014
  !> @brief   Interpolate concentrations from coarse mesh to fine mesh
  !> @details Interpolate concentrations from coarse mesh to fine mesh using
  !> @        the coupling structures and functions
  !> @}
  !----------------------------------------------------------------------

subroutine chm_coarfine(itask)
  use def_domain,        only :  npoin
  use def_kintyp,        only :  ip,rp
  use def_master,        only :  current_code
  use def_master,        only :  conce
  use def_master,        only :  nspec
  use def_master,        only :  INOTMASTER
  use mod_couplings,     only :  COU_INTERPOLATE_NODAL_VALUES
  implicit none
  integer(ip), intent(in)     :: itask
  real(rp),    pointer        :: xvalu(:)
  real(rp),    pointer        :: current_conce(:)
  integer(ip)                 :: i

  nullify(xvalu)
  nullify(current_conce)

  if ( INOTMASTER ) then
     allocate(xvalu(npoin))
     allocate(current_conce(npoin))
  else
     allocate(xvalu(1_ip))
     allocate(current_conce(1_ip))
  end if

  if ( itask == 1 .and. current_code == 2 ) then !

     do i = 1_ip, nspec

        current_conce(:) = conce(:,i,1)

        ! Arguments for the coupling function:
        ! COU_INTERPOLATE_NODAL_VALUES(coupling label, number of dimensions of the variable to interpolate,
        ! array to store the results, variable to interpolate )
        ! The setup of the communicators, coordinates and everything else is done in COU_INITIALIZE_COUPLING
        call COU_INTERPOLATE_NODAL_VALUES(3_ip,1_ip,xvalu,current_conce)
        if( associated(xvalu) ) then

           conce(:,i,1) = xvalu
           conce(:,i,2) = conce(:,i,1)
           conce(:,i,3) = conce(:,i,1)

        end if

     end do

  else if ( itask == 2 .and. current_code == 1 ) then

     do i = 1_ip, nspec

        current_conce(:) = conce(:,i,1)

        call COU_INTERPOLATE_NODAL_VALUES(3_ip,1_ip,xvalu,current_conce)

     end do
  end if

  if( associated(xvalu) ) deallocate( xvalu )
  if( associated(current_conce) ) deallocate( current_conce )

end subroutine chm_coarfine
!>
