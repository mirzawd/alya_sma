!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Coupling
!> @{
!> @file    cou_parall.f90
!> @author  Guillaume Houzeaux
!> @date    03/03/2014
!> @brief   Parallelization
!> @details Parallelization
!> @} 
!-----------------------------------------------------------------------

subroutine cou_parall()
  use def_kintyp
  use def_parame
  use def_master
  use def_kermod
  use def_inpout
  use def_coupli
  use mod_communications
  use mod_parall

  implicit none
  integer(ip)          :: icoup
  integer(ip)          :: code_target,zone_target,subdomain_target
  integer(ip)          :: code_source,zone_source,subdomain_source
  integer(ip)          :: ii,jj
  integer(ip), pointer :: all_zone_source(:)

  !----------------------------------------------------------------------
  !
  ! Broadcast general data
  !
  !----------------------------------------------------------------------

  do parii = 1,2 
     npari = 0 ; nparr = 0 ; nparc = 0
     !
     ! Data
     !
     call PAR_EXCHANGE(mcoup,                        parin,npari,parii)
     call PAR_EXCHANGE(coudt,                        parin,npari,parii)
     call PAR_EXCHANGE(kfl_lost_wet_point_cou,       parin,npari,parii)
     call PAR_EXCHANGE(kfl_absolute_cou,             parin,npari,parii)
     call PAR_EXCHANGE(toler_absolute_cou,           parre,nparr,parii)
     call PAR_EXCHANGE(toler_relative_cou,           parre,nparr,parii)
     call PAR_EXCHANGE(kfl_timco_cou,                parin,npari,parii)
     !
     ! Driver
     !
     do ii = 1,max_block_cou
        do jj = 1,max_coupl_cou
           call PAR_EXCHANGE(coupling_driver_couplings(jj,ii),parin,npari,parii)        
        end do
        call PAR_EXCHANGE(coupling_driver_number_couplings(ii),parin,npari,parii)
        call PAR_EXCHANGE(coupling_driver_max_iteration(ii),   parin,npari,parii)
        call PAR_EXCHANGE(coupling_driver_iteration(ii),       parin,npari,parii)
        call PAR_EXCHANGE(coupling_driver_tolerance(ii),       parre,nparr,parii)
     end do
     !
     ! Allocate memory for the first pass
     !
     if( parii == 1 ) then
        allocate(parin(npari))
        allocate(parre(nparr))
        if( ISLAVE  ) call PAR_BROADCAST(parin,'IN MY CODE')
        if( ISLAVE  ) call PAR_BROADCAST(parre,'IN MY CODE')
     else
        if( IMASTER ) call PAR_BROADCAST(parin,'IN MY CODE')
        if( IMASTER ) call PAR_BROADCAST(parre,'IN MY CODE')
     end if
  end do

  deallocate(parin)
  deallocate(parre)

  if( mcoup > 0 ) then

     !----------------------------------------------------------------------
     !
     ! Exchange source zones. This is because if i am not the source code, 
     ! I cannot know the zone it corresponds to
     !
     !----------------------------------------------------------------------

     nullify(all_zone_source)
     allocate( all_zone_source(mcoup) )
     all_zone_source = 0
     if( IMASTER ) then
        do icoup = 1,mcoup
           if( coupling_type(icoup) % code_source == current_code .and. coupling_type(icoup) % module_source /= 0 ) then 
              all_zone_source(icoup) = lzone(coupling_type(icoup) % module_source) 
           end if
        end do
     end if
     call PAR_MAX_ALL(all_zone_source)
     if( IMASTER ) then
        do icoup = 1,mcoup
           coupling_type(icoup) % zone_source = all_zone_source(icoup) 
        end do
     end if
     deallocate( all_zone_source )
     !
     ! Allocate memory for coupling 
     !
     if( ISLAVE ) call cou_memory(1_ip)

     !----------------------------------------------------------------------
     !
     ! Broadcast coupling definition
     !
     !----------------------------------------------------------------------

     do parii = 1,2 
        npari = 0 ; nparr = 0 ; nparc = 0
        !
        ! Data
        !
        do icoup = 1,mcoup
           call PAR_EXCHANGE(coupling_type(icoup) % number,                  parin,npari,parii)  
           call PAR_EXCHANGE(coupling_type(icoup) % itype,                   parin,npari,parii)  
           call PAR_EXCHANGE(coupling_type(icoup) % code_target,             parin,npari,parii)  
           call PAR_EXCHANGE(coupling_type(icoup) % zone_target,             parin,npari,parii) 
           call PAR_EXCHANGE(coupling_type(icoup) % module_target,           parin,npari,parii) 
           call PAR_EXCHANGE(coupling_type(icoup) % subdomain_target,        parin,npari,parii) 
           call PAR_EXCHANGE(coupling_type(icoup) % code_source,             parin,npari,parii)  
           call PAR_EXCHANGE(coupling_type(icoup) % zone_source,             parin,npari,parii)  
           call PAR_EXCHANGE(coupling_type(icoup) % module_source,           parin,npari,parii)
           call PAR_EXCHANGE(coupling_type(icoup) % subdomain_source,        parin,npari,parii) 
           call PAR_EXCHANGE(coupling_type(icoup) % where_type,              parin,npari,parii) 
           call PAR_EXCHANGE(coupling_type(icoup) % where_number,            parin,npari,parii)
           call PAR_EXCHANGE(coupling_type(icoup) % where_number_source,     parin,npari,parii)
           call PAR_EXCHANGE(coupling_type(icoup) % kfl_par_transm,          parin,npari,parii)
           call PAR_EXCHANGE(coupling_type(icoup) % kfl_check_exha,          parin,npari,parii)
           call PAR_EXCHANGE(coupling_type(icoup) % kfl_toda_costa,          parin,npari,parii)
           call PAR_EXCHANGE(coupling_type(icoup) % what,                    parin,npari,parii)
           call PAR_EXCHANGE(coupling_type(icoup) % scheme,                  parin,npari,parii)
           call PAR_EXCHANGE(coupling_type(icoup) % relax,                   parre,nparr,parii)
           call PAR_EXCHANGE(coupling_type(icoup) % scaling_iqnls,           parre,nparr,parii)
           call PAR_EXCHANGE(coupling_type(icoup) % efilter_iqnls,           parre,nparr,parii)
           call PAR_EXCHANGE(coupling_type(icoup) % ranku_iqnls,             parin,npari,parii)
           call PAR_EXCHANGE(coupling_type(icoup) % history_iqnls,           parin,npari,parii)
           call PAR_EXCHANGE(coupling_type(icoup) % itera,                   parin,npari,parii)
           call PAR_EXCHANGE(coupling_type(icoup) % task_compute_and_send,   parin,npari,parii)
           call PAR_EXCHANGE(coupling_type(icoup) % when_compute_and_send,   parin,npari,parii)
           call PAR_EXCHANGE(coupling_type(icoup) % frequ_send,              parin,npari,parii) 
           call PAR_EXCHANGE(coupling_type(icoup) % task_recv_and_assemble,  parin,npari,parii) 
           call PAR_EXCHANGE(coupling_type(icoup) % when_recv_and_assemble,  parin,npari,parii)
           call PAR_EXCHANGE(coupling_type(icoup) % frequ_recv,              parin,npari,parii)
           call PAR_EXCHANGE(coupling_type(icoup) % variable,                parch,nparc,parii)
           call PAR_EXCHANGE(coupling_type(icoup) % conservation,            parin,npari,parii)
           call PAR_EXCHANGE(coupling_type(icoup) % kfl_symmetry,            parin,npari,parii)
           call PAR_EXCHANGE(coupling_type(icoup) % overlap,                 parin,npari,parii)
           call PAR_EXCHANGE(coupling_type(icoup) % ngaus,                   parin,npari,parii)
           call PAR_EXCHANGE(coupling_type(icoup) % kfl_multi_source,        parin,npari,parii)
           call PAR_EXCHANGE(coupling_type(icoup) % kfl_lost_wet_points,     parin,npari,parii)
           call PAR_EXCHANGE(coupling_type(icoup) % when_update_wet_geometry,parin,npari,parii)
           call PAR_EXCHANGE(coupling_type(icoup) % when_update_coupling,    parin,npari,parii)
           call PAR_EXCHANGE(coupling_type(icoup) % temporal_predictor,      parin,npari,parii)
           call PAR_EXCHANGE(coupling_type(icoup) % temporal_predictor_order,parin,npari,parii)
        end do
        !
        ! Allocate memory for the first pass
        !
        if( parii == 1 ) then
           allocate(parin(npari))
           allocate(parre(nparr))
           if( ISLAVE  ) call PAR_BROADCAST(parin,      'IN MY CODE')
           if( ISLAVE  ) call PAR_BROADCAST(parre,      'IN MY CODE')
           if( ISLAVE  ) call PAR_BROADCAST(nparc,parch,'IN MY CODE')
        else
           if( IMASTER ) call PAR_BROADCAST(parin,      'IN MY CODE')
           if( IMASTER ) call PAR_BROADCAST(parre,      'IN MY CODE')
           if( IMASTER ) call PAR_BROADCAST(nparc,parch,'IN MY CODE')
        end if
     end do

     deallocate(parin)
     deallocate(parre)

     !----------------------------------------------------------------------
     ! 
     ! Define colors of the coupling
     !
     !----------------------------------------------------------------------

     do icoup = 1,mcoup
        code_target                         = coupling_type(icoup) % code_target
        zone_target                         = coupling_type(icoup) % zone_target
        subdomain_target                    = coupling_type(icoup) % subdomain_target
        code_source                         = coupling_type(icoup) % code_source
        zone_source                         = coupling_type(icoup) % zone_source
        subdomain_source                    = coupling_type(icoup) % subdomain_source
        coupling_type(icoup) % color_target = par_code_zone_subd_to_color(code_target,zone_target,subdomain_target)
        coupling_type(icoup) % color_source = par_code_zone_subd_to_color(code_source,zone_source,subdomain_source)
        color_target = coupling_type(icoup) % color_target
        color_source = coupling_type(icoup) % color_source
     end do
     
  end if

  npari = 0 ; nparr = 0 ; nparc = 0

end subroutine cou_parall
