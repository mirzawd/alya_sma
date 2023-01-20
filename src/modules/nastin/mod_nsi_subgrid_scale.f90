!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Nastin
!> @{
!> @file    mod_nsi_subgrid_scale.f90
!> @author  Guillaume Houzeaux
!> @date    23/07/2015
!> @brief   Subrgid scale
!> @details Subrgid scale operation
!> @} 
!-----------------------------------------------------------------------
module mod_nsi_subgrid_scale

  use def_kintyp, only : ip,rp,r3p
  use def_nastin, only : kfl_sgsco_nsi
  use def_nastin, only : kfl_sgsti_nsi
  implicit none  
  private

  interface nsi_subgrid_scale_residual_and_update
     module procedure nsi_subgrid_scale_residual_and_update_scalar,&
          &           nsi_subgrid_scale_residual_and_update_vector
  end interface nsi_subgrid_scale_residual_and_update

  public :: nsi_subgrid_scale_gather
  public :: nsi_subgrid_scale_residual_and_update

contains 

  !-----------------------------------------------------------------------
  !
  !> @brief   Gather subgrid scale
  !> @details Gather the subgrid scale at the element level
  !> @author  Guillaume Houzeaux
  !
  !-----------------------------------------------------------------------

  subroutine nsi_subgrid_scale_gather(ndime,pgaus,ielem,vesgs,gpsgs)
    integer(ip), intent(in)            :: ndime
    integer(ip), intent(in)            :: pgaus
    integer(ip), intent(in)            :: ielem
    type(r3p),   intent(in),   pointer :: vesgs(:)
    real(rp),    intent(out)           :: gpsgs(ndime,pgaus,2)
    
    if( kfl_sgsco_nsi == 1 .or. kfl_sgsti_nsi == 1 ) then

       gpsgs(1:ndime,1:pgaus,1) = vesgs(ielem) % a(1:ndime,1:pgaus,1)

       if( kfl_sgsti_nsi == 1 ) then
          gpsgs(1:ndime,1:pgaus,2) = vesgs(ielem) % a(1:ndime,1:pgaus,2)
       else
          gpsgs(1:ndime,1:pgaus,2) = 0.0_rp
       end if

    else

       gpsgs = 0.0_rp

    end if

  end subroutine nsi_subgrid_scale_gather

  !-----------------------------------------------------------------------
  !
  !> @brief   Gather subgrid scale
  !> @details Compute subgrid scale residual and update its value
  !> @author  Guillaume Houzeaux
  !
  !-----------------------------------------------------------------------

  subroutine nsi_subgrid_scale_residual_and_update_scalar(ndime,pgaus,ielem,gpsgs,vesgs,resgs)
    integer(ip), intent(in)              :: ndime
    integer(ip), intent(in)              :: pgaus
    integer(ip), intent(in)              :: ielem
    real(rp),    intent(in)              :: gpsgs(ndime,pgaus,*)
    type(r3p),   intent(inout),  pointer :: vesgs(:)
    real(rp),    intent(out)             :: resgs(2)
    integer(ip)                          :: idime,igaus

    if( kfl_sgsco_nsi == 1 .or. kfl_sgsti_nsi == 1 ) then

       do igaus = 1,pgaus
          do idime = 1,ndime 
             resgs(1) = resgs(1) + ( gpsgs(idime,igaus,1) - vesgs(ielem) % a(idime,igaus,1) ) ** 2
             resgs(2) = resgs(2) +   gpsgs(idime,igaus,1) * gpsgs(idime,igaus,1)
             vesgs(ielem) % a(idime,igaus,1) = gpsgs(idime,igaus,1)
          end do
       end do

    end if

  end subroutine nsi_subgrid_scale_residual_and_update_scalar

  subroutine nsi_subgrid_scale_residual_and_update_vector(nsize,ndime,pgaus,list_elements,gpsgs,vesgs,resgs)
    integer(ip), intent(in)              :: nsize
    integer(ip), intent(in)              :: ndime
    integer(ip), intent(in)              :: pgaus
    integer(ip), intent(in)              :: list_elements(nsize)
    real(rp),    intent(in)              :: gpsgs(nsize,ndime,pgaus,*)
    type(r3p),   intent(inout),  pointer :: vesgs(:)
    real(rp),    intent(out)             :: resgs(2)
    integer(ip)                          :: idime,igaus,isize,ielem

    if( kfl_sgsco_nsi == 1 .or. kfl_sgsti_nsi == 1 ) then

       do isize = 1,nsize
          ielem = list_elements(isize)
          if( ielem > 0 ) then
             do igaus = 1,pgaus
                do idime = 1,ndime 
                   resgs(1) = resgs(1) + ( gpsgs(isize,idime,igaus,1) - vesgs(ielem) % a(idime,igaus,1) ) ** 2
                   resgs(2) = resgs(2) +   gpsgs(isize,idime,igaus,1) * gpsgs(isize,idime,igaus,1)
                   vesgs(ielem) % a(idime,igaus,1) = gpsgs(isize,idime,igaus,1)
                end do
             end do
          end if
       end do

    end if

  end subroutine nsi_subgrid_scale_residual_and_update_vector

end module mod_nsi_subgrid_scale
 
