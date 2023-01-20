!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    mod_sld_arrays.f90
!> @author  houzeaux
!> @date    2019-11-16
!> @brief   Solidz arrays
!> @details Solidz arrays
!> @}
!-----------------------------------------------------------------------

module mod_sld_arrays

  use def_kintyp,              only : ip
  use def_master,              only : modul, mem_modul
  use def_master,              only : displ
  use def_domain,              only : ndime, nelem, npoin
  use def_domain,              only : ltype, ngaus
  use mod_arrays,              only : arrays, arrays_number
  use mod_memory,              only : memory_alloca

  private

  public :: sld_arrays

contains

  !-----------------------------------------------------------------------
  !>
  !> @author  houzeaux
  !> @date    2019-11-16
  !> @brief   Solidz arrays
  !> @details Do what you have to do with nastin arrays
  !>
  !-----------------------------------------------------------------------

  subroutine sld_arrays(wtask)

    use def_solidz, only : kfl_rigid_sld
    use def_solidz, only : kfl_sdvar_sld, nsvar_sld, ncomp_sld
    use def_solidz, only : veloc_sld, accel_sld, svegm_sld
    use def_solidz, only : kfl_conta_sld, fcont_sld

    implicit none

    character(len=*), intent(in) :: wtask
    integer(ip)                  :: ielem
    integer(ip)                  :: pelty,pgaus
    !
    ! Displacement, velocity and acceleration
    !
    call arrays(arrays_number('DISPL'),trim(wtask),    displ,ndime,npoin,ncomp_sld)
    call arrays(arrays_number('VELOC'),trim(wtask),veloc_sld,ndime,npoin,ncomp_sld)
    call arrays(arrays_number('ACCEL'),trim(wtask),accel_sld,ndime,npoin,ncomp_sld)
    !
    ! State dependent variables (stored at Gauss Points for each Element)
    !
    if( kfl_sdvar_sld == 1 .and. kfl_rigid_sld == 0 ) then
       call arrays(arrays_number('SVEGM'),trim(wtask),svegm_sld,nelem)
       if( trim(wtask) == 'ALLOCATE' ) then
          do ielem = 1,nelem
             pelty = abs(ltype(ielem))
             pgaus = ngaus(pelty)
             call memory_alloca(mem_modul(1:2,modul),'SVEGM_SLD % A','sld_arrays',svegm_sld(ielem)%a,nsvar_sld,pgaus,2_ip)
          end do
       end if
    end if
    !
    ! PDN Contact 
    !
    if( kfl_conta_sld /= 0_ip .and. kfl_rigid_sld == 0 ) then
      call arrays(arrays_number('FCONT'),trim(wtask),fcont_sld,ndime,npoin)
    end if
    
  end subroutine sld_arrays

end module mod_sld_arrays
