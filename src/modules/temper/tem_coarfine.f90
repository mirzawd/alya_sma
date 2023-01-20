!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



  !----------------------------------------------------------------------
  !> @addtogroup Nastin
  !> @{
  !> @file    tem_coarfin.f90
  !> @author  J.C. Cajas
  !> @date    09/04/2014
  !> @brief   Interpolate temperature from coarse mesh to fine mesh
  !> @details Interpolate temperature from coarse mesh to fine mesh using 
  !> @        the coupling structures and functions
  !> @} 
  !----------------------------------------------------------------------

subroutine tem_coarfine(itask)
  use def_domain,        only :  npoin
  use def_kintyp,        only :  ip,rp
  use def_master,        only :  current_code
  use def_master,        only :  tempe
  use def_master,        only :  INOTMASTER
  use mod_couplings,     only :  COU_INTERPOLATE_NODAL_VALUES
  implicit none
  integer(ip), intent(in)     :: itask
  real(rp),    pointer        :: xvalu(:)

  nullify(xvalu) 

  if ( INOTMASTER ) then
     allocate(xvalu(npoin))
  else
     allocate(xvalu(1_ip))
  end if

  if ( itask == 1 .and. current_code == 2 ) then ! 
     
     ! Arguments for the coupling function: 
     ! COU_INTERPOLATE_NODAL_VALUES(coupling label, number of dimensions of the variable to interpolate, 
     ! array to store the results, variable to interpolate )
     ! The setup of the communicators, coordinates and everything else is done in COU_INITIALIZE_COUPLING
     call COU_INTERPOLATE_NODAL_VALUES(2_ip,1_ip,xvalu,tempe)
     if ( associated(xvalu) ) then

        tempe(:,1) = xvalu
        tempe(:,2) = tempe(:,1) 
        tempe(:,3) = tempe(:,1)

     end if
  
  else if ( itask == 2 .and. current_code == 1 ) then

     call COU_INTERPOLATE_NODAL_VALUES(2_ip,1_ip,xvalu,tempe)
     
  end if

  if( associated(xvalu) ) deallocate( xvalu )

end subroutine tem_coarfine
!>
