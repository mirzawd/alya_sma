!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    sld_begrun.f90
!> @date    14/06/2019
!> @author  Guillaume Houzeaux
!> @brief   Beginning the run...
!> @details Beginning the run... we can compute matrices!
!> @}
!-----------------------------------------------------------------------

subroutine sld_begrun()

  use def_kintyp,    only : ip, rp
  use def_elmtyp,    only : SHELL, BAR3D, ELINT
  use def_master,    only : INOTMASTER
  use def_master,    only : INOTEMPTY
  use def_master,    only : fiber
  use def_master,    only : kfl_timco, dtinv
  use def_domain,    only : ltype, lmate, lelch
  use def_domain,    only : xfiel, npoin, ndime, nelem
#if defined(_OPENMP)
  use mod_parall,    only : par_omp_nelem_chunk
#endif
  use def_solidz,    only : kfl_local_sld, kfl_fiber_sld
  use def_solidz,    only : kfl_vofor_sld, vofor_sld
  use def_solidz,    only : kfl_conta_stent
  use def_solidz,    only : fibts_sld, fibtn_sld
  use def_solidz,    only : kfl_vecto_sld
  use mod_sld_csys,  only : sld_csys_build_jacrot
  use mod_sld_csys,  only : sld_csys_assign_material_axes
  use mod_sld_stent, only : sld_set_boundaries
  use def_solidz,    only : kfl_rigid_sld
  use def_solidz,    only : kfl_savdt_sld
  use def_solidz,    only : safet_sld, safex_sld, safma_sld
  use def_solidz,    only : dtcri_sld, dtinv_sld
  use def_solidz,    only : kfl_timei_sld, SLD_STATIC_PROBLEM
  
  implicit none

  external               :: vecpro
  external               :: sld_updtss
  external               :: sld_velsnd
  external               :: sld_inimat
  
  integer(ip)            :: ipoin, ielem, pmate, pelty
  real(rp)               :: dummr(3),vauxi(3),vmodu
  
  !-------------------------------------------------------------------
  !
  ! Time step options and speed of sound for materials
  !
  !-------------------------------------------------------------------
  !
  ! Speed of sound for materials (at this moment only for vectorized version)
  !
  if( kfl_vecto_sld .and. kfl_timei_sld /= SLD_STATIC_PROBLEM ) then

     if( INOTMASTER ) then
        
        !---------------------------------------------------
        !$OMP PARALLEL  DO                                 &
        !$OMP SCHEDULE  ( DYNAMIC , par_omp_nelem_chunk )  &
        !$OMP DEFAULT   (NONE)                             &
        !$OMP PRIVATE   (ielem,pelty,pmate)                &
        !$OMP SHARED    (ltype,lmate,nelem,lelch,          &
        !$OMP            par_omp_nelem_chunk)        
        !---------------------------------------------------
        !
        ! Loop over elements
        !
        elements: do ielem = 1,nelem

           pelty = ltype(ielem)

           if( pelty > 0      .and. &
               pelty /= SHELL .and. &
               pelty /= BAR3D .and. &
               lelch(ielem) /= ELINT ) then
              pmate = lmate(ielem)
              !
              ! Speed of sound according to the material
              !
              call sld_velsnd(pmate,ielem)
           end if

        end do elements
        !$OMP END PARALLEL DO
        !---------------------------------------------------
     end if
     
  end if
  !
  ! Save time step at the begining of the run
  !
  if( kfl_savdt_sld == 1_ip ) then
     if(      kfl_rigid_sld == 0_ip ) then        
        safet_sld = min(safet_sld*safex_sld, safma_sld)
        call sld_updtss()
     else if( kfl_rigid_sld == 1_ip ) then
        dtcri_sld = 1.0_rp
        dtinv_sld = 1.0_rp
        if( kfl_timco == 1 ) dtinv = max(dtinv,dtinv_sld)
     end if
  end if
  
  !-------------------------------------------------------------------
  !
  ! Local axes
  !
  !-------------------------------------------------------------------
  !
  ! Build rotation matrix for local axes prescription
  !
  if( INOTMASTER ) then
     if( kfl_local_sld /= 0_ip ) call sld_csys_build_jacrot()
  end if

  !-------------------------------------------------------------------
  !
  ! Material axes
  !
  !-------------------------------------------------------------------
  !
  ! Material axes at element level (fibers)
  !
  if( INOTEMPTY ) then
     if( kfl_fiber_sld > 3 ) call sld_csys_assign_material_axes()
  end if
  !
  ! Material axes at node level (fibers)
  !
  if (kfl_fiber_sld == 3 .and. INOTMASTER) then
     dummr(:) = 0.0_rp
     do ipoin = 1,npoin

        vauxi   = 0.0_rp
        vauxi(1)= 1.0_rp
        call vecpro(fiber(1,ipoin),vauxi,fibts_sld(1,ipoin),ndime)
        dummr(1:ndime) = fibts_sld(1:ndime,ipoin)
        vmodu= sqrt(dummr(1)*dummr(1) + dummr(2)*dummr(2) + dummr(3)*dummr(3))
        if (vmodu < 1.0e-10_rp) then
           vauxi   = 0.0_rp
           vauxi(2)= 1.0_rp
           call vecpro(fiber(1,ipoin),vauxi,fibts_sld(1,ipoin),ndime)
           dummr(1:ndime) = fibts_sld(1:ndime,ipoin)
           vmodu= sqrt(dummr(1)*dummr(1) + dummr(2)*dummr(2) + dummr(3)*dummr(3))
        end if

        fibts_sld(1:ndime,ipoin) = fibts_sld(1:ndime,ipoin)/vmodu

        call vecpro(fiber(1,ipoin),fibts_sld(1,ipoin),fibtn_sld(1,ipoin),ndime)
        dummr(1:ndime) = fibtn_sld(1:ndime,ipoin)
        vmodu= sqrt(dummr(1)*dummr(1) + dummr(2)*dummr(2) + dummr(3)*dummr(3))

        fibtn_sld(1:ndime,ipoin) = fibtn_sld(1:ndime,ipoin)/vmodu

     end do

  end if

  !-------------------------------------------------------------------
  !
  ! Material model initializations
  !
  !-------------------------------------------------------------------
  !
  ! Pre-calculation of some material properties (at this moment only for vecto)
  !
  if( kfl_vecto_sld ) then
     call sld_inimat()
  end if
  
  !-------------------------------------------------------------------
  !
  ! Concentrated Loads
  !
  !-------------------------------------------------------------------
  !
  ! Concentrated loads using fields
  !
  if( kfl_vofor_sld > 0 ) then
     vofor_sld => xfiel(kfl_vofor_sld) % a(:,:,1)
  end if
  
  !-------------------------------------------------------------------
  !
  ! Stent
  !
  !-------------------------------------------------------------------
  ! 
  ! Set boundaries for stent case 
  !
  if (INOTMASTER) then
     if (kfl_conta_stent /= 0_ip) then
        call sld_set_boundaries()
        call sld_csys_build_jacrot()
     end if
  end if
  
end subroutine sld_begrun
