!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



  !----------------------------------------------------------------------
  !> @addtogroup Nastin
  !> @{
  !> @file    nsi_coarfin.f90
  !> @author  J.C. Cajas
  !> @date    09/04/2014
  !> @brief   Interpolate velocity from coarse mesh to fine mesh
  !> @details Interpolate velocity from coarse mesh to fine mesh using 
  !> @        the coupling structures and functions
  !> @} 
  !----------------------------------------------------------------------

subroutine nsi_coarfine(itask)
  use def_domain,         only :  ndime
  use def_domain,         only :  npoin
  use def_kintyp,         only :  ip,rp
  use def_master,         only :  current_code
  use def_master,         only :  veloc
  use def_master,         only :  press
  use def_master,         only :  INOTMASTER
  use mod_couplings,      only :  COU_INTERPOLATE_NODAL_VALUES
  use mod_communications, only :  PAR_BARRIER
  implicit none
  integer(ip), intent(in)     :: itask
  real(rp),    pointer        :: xvalu(:,:,:)
  real(rp),    pointer        :: pvalu(:,:)

  nullify(xvalu) 
  nullify(pvalu)

  if ( INOTMASTER ) then
     allocate(xvalu(ndime,npoin,1))
     allocate(pvalu(npoin,1))
  else
     allocate(xvalu(1,1,1))
     allocate(pvalu(1,1))
  end if

  if ( itask == 1 .and. current_code == 2 ) then ! 

     ! Arguments for the coupling function: 
     ! COU_INTERPOLATE_NODAL_VALUES(coupling label, number of dimensions of the variable to interpolate, 
     ! array to store the results, variable to interpolate )
     ! The setup of the communicators, coordinates and everything else is done in COU_INITIALIZE_COUPLING

     ! velocity
     call COU_INTERPOLATE_NODAL_VALUES(1_ip,ndime,xvalu,veloc)  

     if( associated(xvalu) ) then

        veloc(:,:,1) = xvalu(:,:,1)
        veloc(:,:,2) = veloc(:,:,1) 
        veloc(:,:,3) = veloc(:,:,1) 

     end if

     ! pressure
     call COU_INTERPOLATE_NODAL_VALUES(1_ip,1_ip,pvalu,press)

     if( associated(pvalu) ) then

        press(:,1) = pvalu(:,1)
        press(:,2) = press(:,1) 
        press(:,3) = press(:,1) 

     end if

     call PAR_BARRIER("IN THE WORLD")

  else if ( itask == 2 .and. current_code == 1 ) then

     ! velocity
     call COU_INTERPOLATE_NODAL_VALUES(1_ip,ndime,xvalu,veloc)

     ! pressure
     call COU_INTERPOLATE_NODAL_VALUES(1_ip,1_ip,pvalu,press)    

     call PAR_BARRIER("IN THE WORLD")

  end if

  if( associated(xvalu) ) deallocate( xvalu )
  if( associated(pvalu) ) deallocate( pvalu )

end subroutine nsi_coarfine
!> 
