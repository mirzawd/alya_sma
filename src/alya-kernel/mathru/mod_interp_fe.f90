!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Domain
!> @{
!> @file    mod_interpolation.f90
!> @author  houzeaux
!> @date    2020-09-04
!> @brief   Interpolation
!> @details Interpolate variables
!-----------------------------------------------------------------------

module mod_interp_fe

  use def_kintyp_basic, only : ip,rp 
  use def_kintyp_basic, only : typ_interp
  use def_domain,       only : memor_dom
  use mod_memory,       only : memory_alloca
  use mod_memory,       only : memory_deallo
  use mod_elmgeo,       only : element_type
  use def_domain,       only : lnods
  use def_domain,       only : ltype
  implicit none
  private
  
  interface interp_fe
     module procedure interp_fe_RP1,&
          &           interp_fe_RP2
  end interface interp_fe
  
  interface interp_fe_deallocate
     module procedure interp_fe_deallocate_RP1,&
          &           interp_fe_deallocate_RP2
  end interface interp_fe_deallocate
  
  public :: interp_fe
  public :: interp_fe_deallocate

contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-09-04
  !> @brief   Interpolation
  !> @details Interpolation of variables 
  !> 
  !-----------------------------------------------------------------------

  subroutine interp_fe_RP1(xx_in,inte,xx_out,nn)

    real(rp),                   pointer, intent(in)    :: xx_in(:)
    real(rp),                   pointer, intent(inout) :: xx_out(:)
    type(typ_interp),                    intent(in)    :: inte
    integer(ip),                         intent(in)    :: nn
    integer(ip)                                        :: ii,ielem,inode,pelty,ipoin

    if( .not. associated(xx_out) ) &
         call memory_alloca(memor_dom,'XX_OUT','outvar',xx_out,nn)
    do ii = 1,nn
       ielem = inte % lelem(ii)
       pelty = abs(ltype(ielem))
       do inode = 1,element_type(pelty) % number_nodes
          ipoin = lnods(inode,ielem)
          xx_out(ii) = xx_out(ii) + inte % shapf(inode,ii) * xx_in(ipoin)
       end do
    end do

  end subroutine interp_fe_RP1
  
  subroutine interp_fe_RP2(xx_in,inte,xx_out,nn)

    real(rp),                  pointer, intent(in)    :: xx_in(:,:)
    real(rp),                  pointer, intent(inout) :: xx_out(:,:)
    type(typ_interp),                   intent(in)    :: inte
    integer(ip),                        intent(in)    :: nn
    integer(ip)                                       :: ii,ielem,inode,pelty,kdime,ipoin

    kdime = size(xx_in,1_ip)
    if( .not. associated(xx_out) ) &
         call memory_alloca(memor_dom,'XX_OUT','outvar',xx_out,kdime,nn)
    do ii = 1,nn
       ielem = inte % lelem(ii)
       pelty = abs(ltype(ielem))
       do inode = 1,element_type(pelty) % number_nodes
          ipoin = lnods(inode,ielem)
          xx_out(1:kdime,ii) = xx_out(1:kdime,ii) + inte % shapf(inode,ii) * xx_in(1:kdime,ipoin)
       end do
    end do

  end subroutine interp_fe_RP2
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-09-04
  !> @brief   Deallocate 
  !> @details Deallocate interp_fe variables 
  !> 
  !-----------------------------------------------------------------------

  subroutine interp_fe_deallocate_RP1(xx_out)

    real(rp), pointer, intent(inout) :: xx_out(:)

    call memory_deallo(memor_dom,'XX_OUT','outvar',xx_out)

  end subroutine interp_fe_deallocate_RP1
  
  subroutine interp_fe_deallocate_RP2(xx_out)

    real(rp), pointer, intent(inout) :: xx_out(:,:)

    call memory_deallo(memor_dom,'XX_OUT','outvar',xx_out)
 
  end subroutine interp_fe_deallocate_RP2
  
end module mod_interp_fe
!> @}
