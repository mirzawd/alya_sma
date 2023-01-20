!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Partis
!> @{
!> @file    pts_outwig.f90
!> @author  houzeaux
!> @date    2020-01-29
!> @brief   Witness output
!> @details Output of witness geometries
!> @} 
!-----------------------------------------------------------------------

subroutine pts_outwig()

  use def_parame
  use def_master
  use def_kermod
  use def_domain
  use def_partis
  use mod_ker_proper,     only : ker_proper
  use mod_communications, only : PAR_SUM
  use mod_witness,        only : witness_in_geometry
  use mod_witness,        only : WITNESS_RING 
  use mod_witness,        only : witness_ring_transform
  use mod_witness,        only : witness_volume_geometry
  use mod_physics,        only : physics_set_liquid_temperature
  use mod_pts_particle,   only : pts_particle_diameter
  implicit none
  integer(ip) :: iwitg,ilagr,ilagr_local,ivarg
  real(rp)    :: denpa,diame,VolCV
  real(rp)    :: n_drop, vel_loc(ndime)
  real(rp)    :: tempnumer(kfl_max_prop_gwit_pts) 
  real(rp)    :: tempdenomSauter 

  !
  ! Loop through geometric witness points
  !
  do iwitg = 1,nwitg
     !
     ! Initialize local sums 
     !
     tempnumer       = 0.0_rp
     tempdenomSauter = 0.0_rp

     !
     ! Control volume
     ! 
     if( postp(1) % npp_witng(14) > 0 ) then
        VolCV   = witness_volume_geometry( gewit(iwitg) % kfl_geometry,gewit(iwitg) % param)
     endif

     !
     ! Loop through particles
     !
     do ilagr_local = 1,nlagr_local_pts
        !
        ! Get global particle index
        !
        ilagr = permu_nlagr_pts(ilagr_local)

        !
        ! Check if particle is within control colume of geometric witness
        !
        if( witness_in_geometry( lagrtyp(ilagr) % coord,      &
           &                     gewit(iwitg) % kfl_geometry, &
           &                     gewit(iwitg) % param         )) then
           
           !==================!
           ! Process particle !
           !==================!
           !
           ! Properties common to multiple methods
           !
           n_drop = parttyp(lagrtyp(ilagr) % itype) % n_drop

           !
           ! Diameter
           !
           if(  postp(1) % npp_witng(2) > 0 .or. &
              & postp(1) % npp_witng(3) > 0 .or. &
              & postp(1) % npp_witng(8) > 0 ) then
              if( parttyp(lagrtyp(ilagr) % itype) % kfl_therm == 0 ) then
                 denpa  = parttyp(lagrtyp(ilagr) % itype) % denpa
              else
                 call physics_set_liquid_temperature( parttyp(lagrtyp(ilagr) % itype) % liq , lagrtyp(ilagr) % tempe_k)
                 denpa  = parttyp(lagrtyp(ilagr) % itype) % liq % rho 
              endif
              diame     = pts_particle_diameter(lagrtyp(ilagr) % itype,ilagr,denpa)
           endif

           !
           ! Velocity
           !
           if(  postp(1) % npp_witng(5)  > 0 .or. &
              & postp(1) % npp_witng(10) > 0 .or. &
              & postp(1) % npp_witng(14) > 0 ) then
              vel_loc   = lagrtyp(ilagr) % veloc
              !
              ! Transform velocity to local coordinate system of ring
              !
              if (gewit(iwitg) % kfl_geometry == WITNESS_RING) call witness_ring_transform(lagrtyp(ilagr) % coord,gewit(iwitg) % param,vec=vel_loc)
           endif

           !
           ! 1: NUMBE 
           !    Number of particles is always processed, becuase others need it too
           !                      \int_t n(t) dt
           ! Number of particles: ------------
           !                        \int_t dt
           !
           !
           tempnumer(1) = tempnumer(1) + n_drop

           !
           ! 2: DIA10 
           !                   \sum_n d
           ! Average diameter:----------
           !                   \sum_n 1
           !
           if( postp(1) % npp_witng(2) > 0 ) then
              tempnumer(2) = tempnumer(2) + n_drop * diame
           endif

           !
           ! 3: DIA32 
           !                       \sum_n d^3
           ! Sauter mean diameter:------------
           !                       \sum_n d^2
           !
           if( postp(1) % npp_witng(3) > 0 ) then
              tempnumer(3)    = tempnumer(3)    + n_drop * diame**3
              tempdenomSauter = tempdenomSauter + n_drop * diame**2
           endif 

           !
           ! 4: TEMPE 
           !                      \sum_n T
           ! Average temperature:----------
           !                      \sum_n 1
           !
           if( postp(1) % npp_witng(4) > 0 ) then
              tempnumer(4) = tempnumer(4) + n_drop * lagrtyp(ilagr) % tempe_k
           endif

           !
           ! 5-7: VELOC
           !                   \sum_n U
           ! Average velocity:----------
           !                   \sum_n 1
           !
           if( postp(1) % npp_witng(5) > 0 ) then
              tempnumer(5:4+ndime) = tempnumer(5:4+ndime) + n_drop * vel_loc
           endif

           !
           ! 8: DIA20
           !                          \sum_n d^2
           ! Average diameter square:------------
           !                           \sum_n 1
           !
           if( postp(1) % npp_witng(8) > 0 ) then
              tempnumer(8) = tempnumer(8) + n_drop * diame**2
           endif

           !
           ! 9: TEMP2
           !                             \sum_n T^2
           ! Average temperature square:------------
           !                              \sum_n 1
           !
           if( postp(1) % npp_witng(9) > 0 ) then
              tempnumer(9) = tempnumer(9) + n_drop * lagrtyp(ilagr) % tempe_k**2
           endif

           !
           ! 10-12: VELO2
           !                          \sum_n U^2
           ! Average velocity square:------------
           !                           \sum_n 1
           !
           if( postp(1) % npp_witng(10) > 0 ) then
              tempnumer(10:9+ndime) = tempnumer(10:9+ndime) + n_drop * vel_loc**2
           endif 

           !
           ! 13: MASS 
           !               \sum_n m
           ! Average mass:-----------
           !               \sum_n 1
           !
           if( postp(1) % npp_witng(13) > 0 ) then
              tempnumer(13) = tempnumer(13) + n_drop * lagrtyp(ilagr) % mass_k
           endif

           !
           ! 14-16: MFLUX
           !                           \sum_n m U             \sum_n m U  
           ! Average massu flux:-------------------------- = ------------
           !                     Area*thickness  \sum_n 1     V \sum_n 1
           !
           if( postp(1) % npp_witng(14) > 0 ) then
              tempnumer(14:13+ndime) = tempnumer(14:13+ndime) + n_drop * vel_loc * lagrtyp(ilagr) % mass_k / VolCV
           endif
        endif
     enddo

     !
     ! Sum results of one geometric witness point over subdomains
     !
     call PAR_SUM(kfl_max_prop_gwit_pts, tempnumer)
     call PAR_SUM(tempdenomSauter)


     !
     ! Give final results back
     ! 
     if (IMASTER) then
         do ivarg = 1,kfl_max_prop_gwit_pts
            if( postp(1) % npp_witng(ivarg) > 0 ) then 
                witng(ivarg,iwitg)                     = witng(ivarg,iwitg) + tempnumer(ivarg)
                if ( ivarg == 3 ) then
                   postp(1) % witng_deldenom(ivarg,iwitg) = tempdenomSauter
                elseif ( ivarg == 1 ) then
                   postp(1) % witng_deldenom(ivarg,iwitg) = 0.0_rp
                else
                   postp(1) % witng_deldenom(ivarg,iwitg) = tempnumer(1)
                endif 
            endif
         enddo
     endif
  enddo

end subroutine pts_outwig

