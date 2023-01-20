!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Coupling
!> @{
!> @file    cou_initialize_coupling.f90
!> @author  Guillaume Houzeaux
!> @date    03/03/2014
!> @brief   Read coupling data
!> @details Initialize coupling data
!>
!>          ISTAT(in)  =   1 ... First iteration
!>                     =   2 ... Force interpolation in all couplings
!>          ISTAT(out) =   0 ... No null rows in transmat
!>                         2 ... Null rows in the transpose couplings
!> @} 
!-----------------------------------------------------------------------

subroutine cou_initialize_coupling(istat)

  use def_kintyp,                   only : ip,rp,lg
  use def_master,                   only : INOTMASTER
  use def_master,                   only : intost,momod
  use mod_couplings,                only : couplings_select_boundaries
  use def_domain,                   only : ndime,npoin
  use def_domain,                   only : coord
  use def_domain,                   only : nelem
  use def_domain,                   only : lesub
  use def_coupli,                   only : mcoup
  use def_coupli,                   only : kfl_graph_cou
  use def_coupli,                   only : ON_MIRROR
  use def_coupli,                   only : PROJECTION
  use def_coupli,                   only : STRESS_PROJECTION
  use def_coupli,                   only : TRANSPOSE_MIRROR
  use def_coupli,                   only : coupling_type
  use def_coupli,                   only : memor_cou
  use def_coupli,                   only : BETWEEN_SUBDOMAINS
  use def_coupli,                   only : BOUNDARY_INTERPOLATION
  use def_coupli,                   only : nboun_cou
  use def_coupli,                   only : lnodb_cou
  use def_coupli,                   only : ltypb_cou
  use def_coupli,                   only : lnnob_cou
  use def_coupli,                   only : lelbo_cou
  use def_coupli,                   only : SCALAR
  use def_coupli,                   only : STOP_ALYA     
  use def_coupli,                   only : STOP_ALYA_WITH_WARNINGS 
  use def_coupli,                   only : kfl_lost_wet_point_cou
  use def_coupli,                   only : ON_FLOATING_POINTS
  use def_coupli,                   only : ON_IMMERSED_MESH
  use def_coupli,                   only : kfl_immer 
  use def_coupli,                   only : kfl_dimbou
  use mod_parall,                   only : color_target
  use mod_parall,                   only : color_source
  use mod_communications,           only : PAR_BARRIER
  use mod_communications,           only : PAR_MAX
  use mod_interpolation,            only : COU_SOURCE_INTERFACE_MASS_MATRIX
  use mod_interpolation,            only : COU_TARGET_INTERFACE_MASS_MATRIX
  use mod_couplings,                only : COU_INIT_INTERPOLATE_POINTS_VALUES
  use mod_couplings,                only : I_AM_IN_COUPLING
  use mod_kdtree,                   only : kdtree_construct
  use mod_memory,                   only : memory_alloca
  use mod_memory,                   only : memory_deallo
  use mod_messages,                 only : livinf
  use mod_couplings_communications, only : COU_GENERATE_LOCAL_TRANSMISSION_MATRICES
  use mod_couplings_communications, only : COU_GENERATE_TRANSPOSED_LOCAL_TRANSMISSION_MATRICES
  use mod_couplings_communications, only : COU_PARALLELIZE_TRANSMISSION_MATRICES
  use mod_domain,                   only : domain_memory_deallocate
  use mod_couplings,                only : couplings_elements_graph
  use mod_couplings,                only : couplings_are_exhaustive
  use mod_couplings,                only : cou_wet_points_from_mirror
  use mod_immersed,                 only : cou_check_if_fringe_nodes_are_necessary
  use mod_run_config,               only : run_config
  use def_master,                   only : times

  implicit none
  integer(ip),           intent(inout) :: istat
  integer(ip)                          :: itype, jtype
  integer(ip)                          :: icoup
  integer(ip)                          :: ipoin,kpoin,jcoup,zone_source
  integer(ip)                          :: subdomain_source,ierro,code_target,thereis_float_point
  integer(ip)                          :: pnodb,kboun,iboun,ielem,jwhere
  logical(lg),    pointer              :: lesou(:)
  integer(ip),    pointer              :: lbsou(:)
  logical(lg),    pointer              :: boundary_mask(:)
  real(rp)                             :: time1,time2
  logical(lg)                          :: all_interp   !force all couplings to be interpolation?
  logical(lg)                          :: kfl_advec = .false.
  logical(lg)                          :: kfl_remsk = .false.
  character(100), PARAMETER            :: vacal = "cou_initialize_coupling"
  logical(lg)                          :: check_exha


  nullify(lesou)
  nullify(lbsou)
  nullify(boundary_mask)

  all_interp = .false.
  if(istat == 2_ip ) all_interp = .true.
  istat = 0_ip

  if( mcoup > 0 ) then

#ifdef ALYA_EXTRAE
     call extrae_eventandcounters(900,98_8) 
#endif

     call momod(0) % times(10) % ini()
     
     if(.not. all_interp) then
        call livinf(0_ip,'COUPLI: INITIALIZE COUPLING',0_ip)
     else
        call livinf(0_ip,'COUPLI: INITIALIZE COUPLING SWITCHING TO ALL INTERPOLATION',0_ip)
     endif

     !-------------------------------------------------------------------
     !     
     ! Loop over couplings
     !
     !-------------------------------------------------------------------
     do icoup = 1,mcoup
                
        if( I_AM_IN_COUPLING(icoup) .and. coupling_type(icoup) % what /= SCALAR ) then

           code_target      = coupling_type(icoup) % code_target
           zone_source      = coupling_type(icoup) % zone_source
           subdomain_source = coupling_type(icoup) % subdomain_source
           color_target     = coupling_type(icoup) % color_target
           color_source     = coupling_type(icoup) % color_source
           jcoup            = coupling_type(icoup) % mirror_coupling
           itype            = coupling_type(icoup) % itype
           if( jcoup /= 0 ) then
              jtype         = coupling_type(jcoup) % itype
              jwhere        = coupling_type(jcoup) % where_type
           else
              jtype         = 0
              jwhere        = 0
           end if
           !
           ! Check: wet surface is valid
           !
           if(    all_interp                      .and. &  
                & itype /= BOUNDARY_INTERPOLATION .and. &
                & itype /= TRANSPOSE_MIRROR               )then
              call runend("Cou_initialize_coupling: case not implemented")    
           endif

           !
           ! Compute KDTree of the sources if necessary
           !
           if( run_config%timing ) call PAR_BARRIER('IN CURRENT COUPLING')
           call cputim(time1)           
           !
           ! Check if fringe wetnodes are necessary
           !
           call cou_check_if_fringe_nodes_are_necessary(icoup)
           !
           ! Get wet nodes from mirror coupling
           !              
           if( coupling_type(icoup) % where_type == ON_MIRROR ) then
              call cou_wet_points_from_mirror(coupling_type(icoup))
           end if
           !
           ! Identify if any immersed coupling exists
           !
           if( coupling_type(icoup) % where_type == ON_IMMERSED_MESH) kfl_immer = .true.
           !
           ! Identify if variables are "advection" or "residual momentum sink" to identify if
           ! I'm using the deformable body immersed boundary method 
           !
           if( coupling_type(icoup) % variable == 'ADVEC' ) kfl_advec = .true.
           if( coupling_type(icoup) % variable == 'REMSK' ) kfl_remsk = .true.

           if( INOTMASTER ) then

              if(  itype == STRESS_PROJECTION .or. itype == PROJECTION        ) then
                 !
                 ! Get mirror kdtree
                 !
                 if( jcoup /= 0 ) then
                    call kdtree_construct(&
                       nboun_cou,npoin,lnodb_cou,ltypb_cou,coord,coupling_type(icoup) % geome % kdtree,&
                       coupling_type(jcoup) % wet % lboun_wet)
                 end if

              else if( itype == BOUNDARY_INTERPOLATION .or. all_interp ) then
                 !
                 ! Get mirror kdtree
                 !
                 if( jcoup /= 0 .and. jtype == TRANSPOSE_MIRROR .and. jwhere /= ON_MIRROR ) then
                    call kdtree_construct(&
                       nboun_cou,npoin,lnodb_cou,ltypb_cou,coord,coupling_type(icoup) % geome % kdtree,&
                       coupling_type(jcoup) % wet % lboun_wet)
                 else
                    !
                    ! Consider only boundaries involved as a source
                    !
                    call memory_alloca(memor_cou,'LESOU',vacal,lesou,nelem)
                    if( coupling_type(icoup) % kind == BETWEEN_SUBDOMAINS ) then
                       do ielem = 1,nelem
                          if( lesub(ielem) == coupling_type(icoup) % subdomain_source ) lesou(ielem) = .true.
                       end do
                    else         
                       do ielem = 1,nelem
                          lesou(ielem) = .true.
                       end do
                    end if
                    kboun = 0
                    do iboun = 1,nboun_cou
                       pnodb = lnnob_cou(iboun)
                       ielem = lelbo_cou(iboun)
                       if( lesou(ielem) ) kboun = kboun + 1
                    end do
                    call memory_alloca(memor_cou,'LBSOU',vacal,lbsou,kboun)
                    kboun = 0
                    do iboun = 1,nboun_cou
                       pnodb = lnnob_cou(iboun)
                       ielem = lelbo_cou(iboun)
                       if( lesou(ielem) ) then
                          kboun = kboun + 1
                          lbsou(kboun) = iboun
                       end if
                    end do
                    call kdtree_construct(nboun_cou,npoin,lnodb_cou,ltypb_cou,coord,coupling_type(icoup) % geome % kdtree,lbsou)
                    call memory_deallo(memor_cou,'LBSOU',vacal,lbsou)
                    call memory_deallo(memor_cou,'LESOU',vacal,lesou)
                 endif
              end if
           end if

           call cputim(time2) 
           coupling_type(icoup) % cputim(2) = coupling_type(icoup) % cputim(2) + time2 - time1
           !
           ! Get the CPUs in charge of my wet points
           !
           call livinf(0_ip,'COUPLI: DEFINE INTERPOLATION OF WET NODES FOR COUPLING '//trim(intost(icoup)),0_ip)
           ierro = 0
                      
           if( .not. coupling_type(icoup) % itype == TRANSPOSE_MIRROR .or. all_interp ) then 
              !
              ! When transpose fails
              !
              if( coupling_type(icoup) % where_type_source /= 0 ) then
                 call couplings_select_boundaries(coupling_type(icoup),boundary_mask,&
                      coupling_type(icoup) % where_type_source,&
                      coupling_type(icoup) % where_number_source)
                 call COU_INIT_INTERPOLATE_POINTS_VALUES(coupling_type(icoup) % wet % coord_wet,color_target,&
                      color_source,coupling_type(icoup) , CANDIDATE_SOURCE_BOUNDARIES=boundary_mask,all_interp_=all_interp)
                 call memory_deallo(memor_cou,'LBSOU',vacal,boundary_mask)
                                  
              else
                 !
                 ! Normal coupling
                 !
                 call COU_INIT_INTERPOLATE_POINTS_VALUES(coupling_type(icoup) % wet % coord_wet,color_target,&
                      color_source,coupling_type(icoup) ,all_interp_=all_interp)
                 
              end if
              !
              ! Check if all points have been found
              !
              if( coupling_type(icoup) % kfl_lost_wet_points == 0 ) then
                 do kpoin = 1,coupling_type(icoup) % wet % number_wet_points
                    ipoin = coupling_type(icoup) % geome % status(kpoin)
                    if( ipoin == 0 ) then
                       if( INOTMASTER .and. kfl_lost_wet_point_cou == STOP_ALYA_WITH_WARNINGS ) &
                            write(*,*) icoup,coupling_type(icoup) % wet % coord_wet(1:ndime,kpoin)
                       ierro = ierro + 1
                    end if
                 end do
              end if
           end if
           if( kfl_lost_wet_point_cou == STOP_ALYA_WITH_WARNINGS .or. kfl_lost_wet_point_cou == STOP_ALYA ) then
              call PAR_MAX(ierro,'IN CURRENT COUPLING')
              if( ierro > 0 ) call runend('COU_INITIALIZE_COUPLING: SOME WET NODES ARE LOST')
           end if
        end if

        ! call PAR_SCHEDULING(coupling_type(icoup))
     end do
     !
     ! Compute mass matrix if required
     !
     if( INOTMASTER ) then
        do icoup = 1,mcoup
           if( coupling_type(icoup) % itype == STRESS_PROJECTION ) then
              jcoup = coupling_type(icoup) % mirror_coupling 
              if( jcoup /= 0 ) then
                 call COU_SOURCE_INTERFACE_MASS_MATRIX(coupling_type(icoup),coupling_type(jcoup))
              end if
           else if( coupling_type(icoup) % itype == PROJECTION ) then
              call COU_TARGET_INTERFACE_MASS_MATRIX(coupling_type(icoup))
           end if
        end do
     end if

     !-------------------------------------------------------------------
     !     
     ! Generater Transmission Matrices
     !
     !-------------------------------------------------------------------
     ! 1. Generate local transmission matrices
     do icoup = 1,mcoup
        if( coupling_type(icoup) % what /= SCALAR ) then
           if( I_AM_IN_COUPLING(icoup) .and. (.not. coupling_type(icoup) % itype == TRANSPOSE_MIRROR &
                .or. all_interp) ) then
              call COU_GENERATE_LOCAL_TRANSMISSION_MATRICES(coupling_type(icoup),all_interp_=all_interp)
           end if
        end if
     end do
     ! 2. Generate trans. matrices for transpose couplings
     do icoup = 1,mcoup
        if( coupling_type(icoup) % what /= SCALAR ) then
           if( I_AM_IN_COUPLING(icoup) .and. ( coupling_type(icoup) % itype == TRANSPOSE_MIRROR &
                .and. .not. all_interp) )  then
              call COU_GENERATE_TRANSPOSED_LOCAL_TRANSMISSION_MATRICES(coupling_type(icoup))
           end if
        end if
     end do
     !
     ! 3. Expand matrices for parallel executions
     !
     do icoup = 1,mcoup
        if( coupling_type(icoup) % what /= SCALAR ) then
           if( I_AM_IN_COUPLING(icoup) .and.  coupling_type(icoup) % kfl_par_transm == 1 ) then
              call COU_PARALLELIZE_TRANSMISSION_MATRICES(coupling_type(icoup))
           end if
        end if
     end do
     !
     ! 4. Chech exhaustivity of couplings (currently for all of them or none)
     !
     check_exha = .false.
     do icoup = 1,mcoup
        if(coupling_type(icoup) % kfl_check_exha == 1) check_exha = .true. 
     end do
     thereis_float_point = 0_ip
     do icoup = 1,mcoup
        if ( coupling_type(icoup) % where_type == ON_FLOATING_POINTS )  thereis_float_point = 1_ip
     end do
     if ( thereis_float_point == 0_ip ) then
        if(check_exha .and. (.not.couplings_are_exhaustive())) istat = 2
     end if
     
#ifdef ALYA_EXTRAE
     call extrae_eventandcounters(900,0_8) 
#endif
     
     call momod(0) % times(10) % add()

  endif

  !
  ! If ADVEC and REMSK couplings are active => Activate deformable immmersed boundary flag
  !
  if (kfl_advec .and. kfl_remsk) kfl_dimbou = .true.

  if(kfl_graph_cou == 1_ip .and. mcoup>0) then

     call couplings_elements_graph()
     call domain_memory_deallocate('PELEL')
     call domain_memory_deallocate('LELEL')

  endif

end subroutine cou_initialize_coupling

