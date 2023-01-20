!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Chemic
!> @{
!> @file    mod_pts_arrays.f90
!> @author  houzeaux
!> @date    2019-11-16
!> @brief   Nastin arrays
!> @details Nastin arrays
!-----------------------------------------------------------------------

module mod_pts_arrays

  use def_parame
  use def_domain
  use def_master
  use def_kermod
  use def_partis
  use def_inpout
  use mod_memory
  use mod_output_postprocess, only : output_postprocess_check_variable_postprocess
  use mod_arrays,             only : arrays
  use mod_arrays,             only : arrays_number
  implicit none

  private

  public :: pts_arrays

contains

  subroutine pts_arrays(wtask)

    character(len=*), intent(in) :: wtask    
    !
    ! Two-way coupling: momentum, hear or vapor mass creatiomn rate is accumulated
    !
    if( kfl_momentum_sink_pts /= 0 ) then
       if( wtask == 'RESIZE' ) then
          call memory_resize(mem_modul(1:2,modul),'MOMSK','pts_inibcs',momentum_sink,ndime,npoin_2)
       else
          call arrays(arrays_number('MOMSK'),wtask,momentum_sink,ndime,npoin)
       end if
    end if
    if( kfl_heat_sink_pts /= 0 ) then
       if( wtask == 'RESIZE' ) then
          call memory_resize(mem_modul(1:2,modul),'HEASK','pts_inibcs',heat_sink,npoin_2)
       else
          call arrays(arrays_number('HEASK'),wtask,heat_sink,npoin)          
       end if
    end if
    if( kfl_mass_sink_pts /= 0 ) then
       if( wtask == 'RESIZE' ) then
          call memory_resize(mem_modul(1:2,modul),'MASSK','pts_inibcs',mass_sink,npoin_2)
       else        
          call arrays(arrays_number('MASSK'),wtask,mass_sink,npoin)
       end if
    end if
    !
    ! Residence time
    !
    if(  output_postprocess_check_variable_postprocess(VARIABLE_NAME='RESID') .or. &
         output_postprocess_check_variable_postprocess(VARIABLE_NAME='NRESI') ) then
       kfl_resid_pts = 1
       if( wtask == 'RESIZE' ) then
          call memory_resize(mem_modul(1:2,modul),'RESID','pts_inibcs',resid_pts,ntyla_pts,nelem_2)
       else
          call arrays(arrays_number('RESID'),wtask,resid_pts,ntyla_pts,nelem)          
       end if
    end if
    !
    ! Deposition
    !
    if(     output_postprocess_check_variable_postprocess(VARIABLE_NAME='DEPOE') .or. &
         &  output_postprocess_check_variable_postprocess(VARIABLE_NAME='DEPOB') ) then
       kfl_depos_pts = 1
       if( ntyla_pts > 0 ) then
          if( wtask == 'RESIZE' ) then
             call memory_resize(mem_modul(1:2,modul),'DEPOE','pts_inibcs',depoe_pts,ntyla_pts,nelem_2)
             call memory_resize(mem_modul(1:2,modul),'DEPOB','pts_inibcs',depob_pts,ntyla_pts,nboun_2)
          else
             call arrays(arrays_number('DEPOE'),wtask,depoe_pts,ntyla_pts,nelem)
             call arrays(arrays_number('DEPOB'),wtask,depob_pts,ntyla_pts,nboun)             
          end if
       end if
    end if
    !
    ! Deposition surface
    !
    if( kfl_depos_pts == 0 .and. kfl_depos_surface_pts == 1 .and. kfl_oudep_pts /= 0 ) then
       if( ntyla_pts > 0 ) then
          if( wtask == 'RESIZE' ) then
             call memory_resize(mem_modul(1:2,modul),'DEPOB','pts_inibcs',depob_pts,ntyla_pts,nboun_2)
          else
             call arrays(arrays_number('DEPOB'),wtask,depob_pts,ntyla_pts,nboun)
          end if
       end if
    end if
    !
    ! Mie scattering equivalent
    !
    if(     output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVMIE') ) then
       if( ntyla_pts > 0 ) then
          if( wtask == 'RESIZE' ) then
             call memory_resize(mem_modul(1:2,modul),'AVMIE','pts_inibcs',avg_mie_pts,npoin_2)
          else
             call arrays(arrays_number('AVMIE'),wtask,avg_mie_pts,npoin)
          end if
       end if
    end if
    
  end subroutine pts_arrays

end module mod_pts_arrays
!> @}
