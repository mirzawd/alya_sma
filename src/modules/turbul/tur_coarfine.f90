!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



  !----------------------------------------------------------------------
  !> @addtogroup Nastin
  !> @{
  !> @file    tur_coarfin.f90
  !> @author  J.C. Cajas
  !> @date    09/04/2014
  !> @brief   Interpolate concentrations from coarse mesh to fine mesh
  !> @details Interpolate concentrations from coarse mesh to fine mesh using 
  !> @        the coupling structures and functions
  !> @} 
  !----------------------------------------------------------------------

subroutine tur_coarfine(itask)
  use def_domain,        only :  npoin
  use def_kintyp,        only :  ip,rp
  use def_master,        only :  current_code
  use def_master,        only :  untur
  use def_turbul,        only :  nturb_tur
  use def_master,        only :  INOTMASTER
  use mod_couplings,     only :  COU_INTERPOLATE_NODAL_VALUES
  implicit none
  integer(ip), intent(in)     :: itask
  real(rp),    pointer        :: xvalu(:)
  real(rp),    pointer        :: current_untur(:)
  integer(ip)                 :: i

  nullify(xvalu)
  nullify(current_untur)

  if ( INOTMASTER ) then
     allocate(xvalu(npoin))
     allocate(current_untur(npoin))
  else
     allocate(xvalu(1_ip))
     allocate(current_untur(1_ip))
  end if

  if ( itask == 1 .and. current_code == 2 ) then ! 
     
     do i = 1_ip, nturb_tur   

        current_untur(:) = untur(i,:,1)

        ! Arguments for the coupling function: 
        ! COU_INTERPOLATE_NODAL_VALUES(coupling label, number of dimensions of the variable to interpolate, 
        ! array to store the results, variable to interpolate )
        ! The setup of the communicators, coordinates and everything else is done in COU_INITIALIZE_COUPLING
        call COU_INTERPOLATE_NODAL_VALUES(2_ip,1_ip,xvalu,current_untur)  
        if( associated(xvalu) ) then

           untur(i,:,1) = xvalu
           untur(i,:,2) = untur(i,:,1) 
           untur(i,:,3) = untur(i,:,1)           
           
        end if

     end do

  else if ( itask == 2 .and. current_code == 1 ) then

     do i = 1_ip, nturb_tur

        current_untur(:) = untur(i,:,1)
        
        call COU_INTERPOLATE_NODAL_VALUES(2_ip,1_ip,xvalu,current_untur)
        
     end do
  end if

  if( associated(xvalu) ) deallocate( xvalu )
  if( associated(current_untur) ) deallocate( current_untur )

end subroutine tur_coarfine
!>
