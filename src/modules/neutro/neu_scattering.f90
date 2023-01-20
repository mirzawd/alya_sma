!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Neutro
!> @{
!> @file    neu_scatterings.f90
!> @date    31/03/2016
!> @author  Guillaume Houzeaux
!> @brief   Scattering
!> @details Here we compute the scattering coefficients for each pair of directions. 
!>          This should be changed by a function that would read the parameters from a table.
!> @} 
!-----------------------------------------------------------------------

subroutine neu_scattering()

  use def_kintyp, only : ip,rp
  use def_master, only : INOTMASTER
  use def_domain, only : ndime
  use def_neutro, only : scattering_neu
  use def_neutro, only : num_directions_neu
  use def_neutro, only : direc_neu
  !use def_neutro, only : aniso_neu
  implicit none
  integer(ip) :: idire,jdire
  real(rp)    :: dotpr

  if( INOTMASTER ) then
     do idire = 1,num_directions_neu   
        do jdire = 1,num_directions_neu
           dotpr = dot_product(direc_neu(1:ndime,idire),direc_neu(1:ndime,jdire))
           scattering_neu(idire,jdire) = 1.0_rp  ! + aniso_neu * dotpr
        end do
     end do
  end if

end subroutine neu_scattering
