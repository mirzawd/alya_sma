!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup Nastin
!> @{
!> @file    nsi_gpu_destroy.f90
!> @author  Guillermo Oyarzun
!> @brief   Destroy data in GPU
!> @details Deallocates GPU data for the repartitioning
!> @}
!------------------------------------------------------------------------
subroutine nsi_gpu_destroy()

  use def_kintyp
  use def_master
  use def_domain
  use def_nastin

  use mod_memory_basic

  integer(ip)::  ndim1_old
  integer(ip)::  ndim2_old
  integer(ip)::  ndim3_old



  ndim1_old = memory_size(veloc,1_ip)
  ndim2_old = memory_size(veloc,2_ip)
  ndim3_old = memory_size(veloc,3_ip)



#ifdef OPENACCHHH
  if (kfl_paral /= 0 ) then
    
     if( kfl_savda == 0 ) then
        !$acc exit data delete (coord,ltype,lnods,lnodb,gravi_nsi,  &
        !$acc                   dt_rho_nsi,mass_rho_nsi, rhsid,     &
        !$acc                   veloc(1:ndim1_old, 1:ndim2_old, 1:ndim3_old))
     else
        !$acc exit data delete (coord,ltype,lnods,lnodb,gravi_nsi,   &
        !$acc                   elmda_gpvol, elmda_gpcar ,           &
        !$acc                   dt_rho_nsi, mass_rho_nsi, rhsid,     &
        !$acc                   veloc(1:ndim1_old, 1:ndim2_old, 1:ndim3_old))
     end if

     if( associated(velom)) then
        !$acc exit data delete(veloc(1:ndim1_old, 1:ndim2_old, 1:ndim3_old))
     end if
  end if  
#endif

end subroutine nsi_gpu_destroy

