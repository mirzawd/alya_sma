!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup Partis
!! @{
!> @name    Partis memory
!! @file    pts_outvar.f90
!> @author  Guillaume Houzeaux
!> @date    28/06/2012
!! @brief   Postprocess particle results on nodes
!! @details Postprocess particle results:
!!          - Number of particles per type
!!          - Desposited particles per type
!> @} 
!------------------------------------------------------------------------

subroutine pts_outvar(ivari,imesh)
  use def_parame
  use def_master
  use def_domain
  use def_partis
  use mod_memory
  use mod_gradie
  use mod_postpr
  use mod_postpr_tools
  use def_kermod,      only : ndivi
  use mod_projec,      only : projec_elements_to_nodes
  use mod_projec,      only : projec_boundaries_to_nodes
  use mod_solver,      only : solver_lumped_mass_system
  use mod_arrays,      only : arrays
  use mod_arrays,      only : arrays_name
  use mod_outvar,      only : outvar
  implicit none
  integer(ip), intent(in) :: ivari
  integer(ip), intent(in) :: imesh   !< Mesh to postprocess on
  integer(ip)             :: ipoin,itype,iboun,ii
  integer(ip)             :: ilagr,ielem
  real(rp),    pointer    :: gelem(:,:)
  real(rp),    pointer    :: depob_tmp(:)
  character(5)            :: wopos_loc(5)
  
  if( ivari == 0 ) return

  select case ( arrays_name(ivari) )  

  case( 'PARTI' )
     !
     ! Number of particles/node
     !
     nullify(gelem)
     call memgen(0_ip,ntyla_pts,max(1_ip,npoin))
     call memory_alloca(mem_modul(1:2,modul),'GELEM','pts_outvar',gelem,ntyla_pts,max(1_ip,nelem))
     if( INOTEMPTY ) then
        !
        ! Count particles inside each element
        !
        do ilagr = 1,mlagr 
           if( lagrtyp(ilagr) % kfl_exist == -1 ) then 
              itype = lagrtyp(ilagr) % itype
              ielem = lagrtyp(ilagr) % ielem
              if( ielem > 0 .and. ielem <= nelem ) then
                 gelem(itype,ielem) = gelem(itype,ielem) + 1.0_rp
              end if
           end if 
        end do
        call projec_elements_to_nodes(gelem,gevec) 
     end if 
     call memory_deallo(mem_modul(1:2,modul),'GELEM','pts_outvar',gelem)
     call arrays(ivari,'POSTPROCESS',gevec,FORCE=.true.)
     call memgen(2_ip,ntyla_pts,1_ip)     
     return
     
  case( 'DEPOS' )
     !
     ! Deposition map
     !
     if( kfl_depos_pts == 0 ) return
     call memgen(0_ip,ntyla_pts,max(npoin,1_ip))
     if( INOTEMPTY ) then
        call projec_boundaries_to_nodes(depob_pts,meshe(ndivi),elmar,gevec)
     end if
     call arrays(ivari,'POSTPROCESS',gevec,FORCE=.true.)
     return
     
  case( 'SLIPW' )
     !
     ! Slip wall distance
     !     
     gesca => walld_slip_pts

  case( 'FRICT' )
     !
     ! Friction
     !     
     gesca => friction_pts

  case( 'RESID' )
     !
     ! Residence
     !
     if( kfl_resid_pts == 0 ) return
     postp(1) % wopos(2,5) = 'MULTI'
     gevec => resid_pts
     call outvar(&
       ivari,ittim,cutim,postp(1) % wopos(:,ivari))
     postp(1) % wopos(2,5) = 'SCALA'
     return

  case( 'DEPOB' )
     !
     ! Deposition map on boundaries
     !
     call memgen(0_ip,nboun,0_ip)
     ii = len(postp(1) % wopos,1)
     wopos_loc(1:ii) = postp(1) % wopos(1:ii,ivari)
     wopos_loc(1)(2:2) = 'B'
     do itype = 1,ntyla_pts
        call postpr_components(itype,wopos_loc(1))        
        do iboun = 1,nboun
           gesca(iboun) = depob_pts(itype,iboun)
        end do
        call postpr(gesca,wopos_loc,ittim,cutim)
     end do
     call memgen(2_ip,nboun,0_ip)
     return
     
  case( 'DEPOT' )
     !
     ! Total deposition map
     !
     if( kfl_depos_pts == 0 ) return
     
     if( INOTMASTER ) then
        nullify(depob_tmp)
        call memory_alloca(mem_modul(1:2,modul),'DEPOB_TMP','pts_outvar',depob_tmp,max(1_ip,nboun))
        do iboun = 1,nboun
           do itype = 1,ntyla_pts
              if( parttyp(itype) % kfl_exist /= 0 ) then
                 depob_tmp(iboun) = depob_tmp(iboun) + depob_pts(itype,iboun)
              end if
           end do
        end do        
        call memgen(0_ip,npoin,0_ip)
        call projec_boundaries_to_nodes(depob_tmp,meshe(ndivi),elmar,gesca)
        call memory_deallo(mem_modul(1:2,modul),'DEPOB_TMP','pts_outvar',depob_tmp)
     end if
     
  case( 'BOUNC' )
     !
     ! Bouncing wall distance
     !     
     gesca => walld_bouncing_pts

  case( 'TEMPE' )
     !
     ! TEMPE: Fluid temperature
     !
     gesca => therm(:,1)
  
 case( 'MASSK' )
     !
     ! MASSK: Mass sink
     !
     call memgen(zero,max(1_ip,npoin),zero)
     if (associated(mass_sink)) then
        do ipoin=1,npoin
           gesca(ipoin) = mass_sink(ipoin)
        enddo
        call solver_lumped_mass_system(1_ip,gesca,EXCHANGE=.false.)
     else
        do ipoin=1,npoin
           gesca(ipoin) = 0.0_rp 
        end do
     endif

 case( 'HEASK' )
     !
     ! HEASK: Heat sink
     !
     call memgen(zero,max(1_ip,npoin),zero)
     if (associated(heat_sink)) then
         do ipoin=1,npoin
            gesca(ipoin) = heat_sink(ipoin)
        enddo
         call solver_lumped_mass_system(1_ip,gesca,EXCHANGE=.false.)
     else
        do ipoin=1,npoin
            gesca(ipoin) = 0.0_rp 
        end do
     endif

 case( 'MOMSK' )
     !
     ! MOMSK: Momentum sink
     !
     call memgen(zero,ndime,max(1_ip,npoin))
     if (associated(momentum_sink)) then
         do ipoin=1,npoin
            gevec(1:ndime,ipoin) =  momentum_sink(1:ndime,ipoin)
         enddo
         call solver_lumped_mass_system(ndime,gevec,EXCHANGE=.false.)
     else
        do ipoin=1,npoin
            gevec(1:ndime,ipoin) = 0.0_rp 
        end do
     endif

  case( 'DEPOE' )
     !
     ! Deposition on elements
     !
     call memgen(0_ip,nelem,0_ip)
     ii = len(postp(1) % wopos,1)
     wopos_loc(1:ii) = postp(1) % wopos(1:ii,ivari)
     do itype = 1,ntyla_pts
        wopos_loc(1) = postp(1) % wopos(1,ivari)
        call postpr_components(itype,wopos_loc(1))        
        do ielem = 1,nelem
           gesca(ielem) = depoe_pts(itype,ielem)
        end do
        call postpr(gesca,wopos_loc,ittim,cutim)
     end do
     call memgen(2_ip,nelem,0_ip)
     return

  case( 'NRESI' )
     !
     ! Residence time on nodes
     !
     if( kfl_resid_pts == 0 ) return
     call memgen(0_ip,ntyla_pts,npoin)
     call projec_elements_to_nodes(resid_pts,gevec)
     call arrays(ivari,'POSTPROCESS',gevec,FORCE=.true.)
     call memgen(2_ip,ntyla_pts,npoin)
     return

  case( 'AVMIE' )
     !
     ! Average Mie scattering equivalent
     !
     call arrays(ivari,'POSTPROCESS',avg_mie_pts)
     !call memgen(0_ip,max(npoin,1_ip),0_ip)
     !call projec_elements_to_nodes(avg_mie_pts,gesca)
     !call arrays(ivari,'POSTPROCESS',gesca,FORCE=.true.)
     !call memgen(2_ip,max(npoin,1_ip),0_ip)
     !do ielem = 1,nelem
     !   avg_mie_pts(ielem) = 0.0_rp
     !enddo
     do ipoin = 1,npoin
        avg_mie_pts(ipoin) = 0.0_rp
     enddo
     return

  end select

  call outvar(&
       ivari,&
       ittim,cutim,postp(1) % wopos(:,ivari),MESH_ID=imesh)

end subroutine pts_outvar
