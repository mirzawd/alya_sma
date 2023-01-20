!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Coupling
!> @{
!> @file    cou_define_wet_geometry.f90
!> @author  Guillaume Houzeaux
!> @date    03/03/2014
!> @brief   Read coupling data
!> @details Initialize coupling data
!>
!>          \verbatim
!>
!>          --------
!>          Wet mesh
!>          --------
!>
!>               wet boundaries
!>          <----------------------->
!>          o-----o-----o-----o-----o <= wet nodes (o) 
!>
!>          COUPLING % GEOME % NBOUN_WET ........... Number of wet boundaries
!>          COUPLING % GEOME % NPOIN_WET ........... Number of wet nodes
!>          COUPLING % GEOME % LBOUN_WET(:) ........ List of wet boundaries
!>          COUPLING % GEOME % LPOIN_WET(:) ........ List of wet nodes
!> 
!>          ----------
!>          Wet points
!>          ----------
!>
!>          They can be the wet nodes or the gauss points of the 
!>          wet boundaries:
!>
!>          COUPLING % GEOME % NUMBER_WET_POINTS ... Number of wet points
!>          COUPLING % GEOME % COORD_WET(:,:) ...... Coordinates of wet points
!>
!>          wet points are wet nodes:     x-----x-----x-----x-----x
!>          wet points are gauss points:  o-x-x-o-x-x-o-x-x-o-x-x-o
!>
!>          \endverbatim
!> @} 
!-----------------------------------------------------------------------

subroutine cou_define_wet_geometry()
  use def_kintyp,         only :  ip,rp,lg 
  use def_elmtyp,         only :  BOFRI,NOFRI
  use def_master,         only :  INOTMASTER 
  use def_master,         only :  current_code
  use def_master,         only :  kfl_async,gisca
  use def_master,         only :  NBOUN_TYPE
  use def_master,         only :  mmodu
  use def_master,         only :  id_nastin, id_solidz
  use def_domain,         only :  lnnob,ltype,lnods
  use def_domain,         only :  nboun,lbset,nelem
  use def_domain,         only :  xfiel,kfl_codbo
  use def_domain,         only :  npoin,coord
  use def_domain,         only :  lnoch
  use def_domain,         only :  ndime,kfl_field
  use def_domain,         only :  lelbo,lesub
  use def_domain,         only :  coord_cloud_gp
  use def_domain,         only :  gp_total_cloud_gp
  use def_coupli,         only :  coupling_type
  use def_coupli,         only :  mcoup
  use def_coupli,         only :  memor_cou
  use def_coupli,         only :  ON_SET
  use def_coupli,         only :  ON_FIELD
  use def_coupli,         only :  ON_CODE
  use def_coupli,         only :  ON_WHOLE_MESH
  use def_coupli,         only :  ON_IMMERSED_MESH
  use def_coupli,         only :  ON_FLOATING_POINTS
  use def_coupli,         only :  PROJECTION
  use def_coupli,         only :  STRESS_PROJECTION
  use def_coupli,         only :  BETWEEN_SUBDOMAINS
  use def_coupli,         only :  BETWEEN_ZONES
  use def_coupli,         only :  RESIDUAL
  use def_coupli,         only :  ON_CHIMERA_MESH
  use def_coupli,         only :  ON_EMBEDDED_MESH
  use def_coupli,         only :  ON_MIRROR
  use def_coupli,         only :  FLOATING_TARGET_ENTITY
  use def_coupli,         only :  nboun_cou
  use def_coupli,         only :  lnodb_cou
  use def_coupli,         only :  ltypb_cou
  use def_coupli,         only :  lboch_cou
  use def_coupli,         only :  lnnob_cou
  use mod_parall,         only :  par_memor
  use mod_parall,         only :  I_AM_IN_COLOR
  use mod_parall,         only :  color_target 
  use mod_memory,         only :  memory_alloca 
  use mod_memory,         only :  memory_deallo
  use mod_memory,         only :  memory_copy
  use mod_communications, only :  PAR_BARRIER
  use mod_communications, only :  PAR_INTERFACE_NODE_EXCHANGE      
  use mod_couplings,      only :  COU_LIST_SOURCE_NODES
  use mod_couplings,      only :  MIRROR_COUPLING
  use mod_interpolation,  only :  COU_PROJECTION_TYPE
  use mod_holcut,         only :  cou_holcut
  use mod_holcut,         only :  cou_holcut_multicode
  use mod_messages,       only :  livinf
  use mod_cloud_gp,       only :  generate_cloud_gp
  use mod_elmgeo,         only :  elmgeo_shapf_deriv_heslo
  use def_elmgeo,         only :  element_type

  implicit none
  integer(ip)                  :: icoup,ipoin,iboun,inodb,kpoin,ielem
  integer(ip)                  :: set_target,zone_target,code_target,bcode_target
  integer(ip)                  :: nboun_wet,pblty,igaub,ipass,inode
  integer(ip)                  :: field_target,number_wet_nodes,iboun_wet
  integer(ip)                  :: number_wet_points,ngaub,pnodb
  integer(ip)                  :: number_overlap_nodes
  integer(ip), pointer         :: kweight_wet(:)
  integer(ip), pointer         :: list_wet_nodes(:)
  integer(ip), pointer         :: list_overlap_nodes(:)
  logical(lg)                  :: I_must_compute_something
  real(rp)                     :: time1,time5,time6

  nullify(kweight_wet)
  nullify(list_wet_nodes)
  nullify(list_overlap_nodes)
  call cputim(time1)

  if( mcoup > 0 ) then

     call livinf(0_ip,'COUPLI: DEFINE WET GEOMETRY',0_ip)

     if( kfl_async == 1 ) then
        do icoup = 1,mcoup
           if( coupling_type(icoup) % kind == BETWEEN_SUBDOMAINS ) then  
              call runend('COUPLI: ASYNCHRONOUS COMMUNICATION WITH COUPLING IS NOT YET POSSIBLE')
           end if
        end do
     end if
     if( INOTMASTER ) then 
        call memory_alloca(memor_cou,'LIST_WET_NODES'    ,'cou_initialize_coupling',list_wet_nodes,    npoin)
        call memory_alloca(memor_cou,'LIST_OVERLAP_NODES','cou_initialize_coupling',list_overlap_nodes,npoin)
     end if
     
     !-------------------------------------------------------------------
     !     
     ! Loop over couplings
     !
     !-------------------------------------------------------------------

     do ipass = 1,2
        !
        ! Hole cut
        !
        if( ipass == 2 ) then
           call cou_holcut()
!           call cou_holcut_multicode() ! ### david debug ###
        end if

        do icoup = 1,mcoup
           
           code_target  = coupling_type(icoup) % code_target
           color_target = coupling_type(icoup) % color_target
           call cputim(time5)
           !
           ! I should do something only if I am in target
           !
           I_must_compute_something = & 
                &  (   I_AM_IN_COLOR(color_target) .and. code_target == current_code .and. INOTMASTER ) .and. &
                &  ( ( ipass == 1 .and. coupling_type(icoup) % where_type /= ON_CHIMERA_MESH ) .or. &
                &    ( ipass == 2 .and. coupling_type(icoup) % where_type == ON_CHIMERA_MESH ) )


           if( coupling_type(icoup) % where_type == ON_FLOATING_POINTS .and. ipass == 1 .and. &
               I_must_compute_something ) then
              !-----------------------------------------------------------
              !
              ! Floating points
              !
              !-----------------------------------------------------------

              call generate_cloud_gp
              number_wet_points = gp_total_cloud_gp
              if (number_wet_points > 0_ip) then
                 call memory_alloca(memor_cou,'COUPLING % WET % COORD_WET','cou_initialize_coupling',&
                         coupling_type(icoup) % wet % coord_wet,ndime,number_wet_points)
                 coupling_type(icoup) % wet % coord_wet      = coord_cloud_gp
              end if
              coupling_type(icoup) % wet % number_wet_points = number_wet_points
              coupling_type(icoup) % wet % npoin_wet         = number_wet_points
              coupling_type(icoup) % target_entity           = FLOATING_TARGET_ENTITY

           else if( I_must_compute_something ) then

              !-----------------------------------------------------------
              !
              ! My code is the target. Compute:
              ! - COUPLING_TYPE(ICOUP) % WET % NBOUN_WET
              ! - COUPLING_TYPE(ICOUP) % WET % LBOUN_WET(:,:)
              ! - LIST_WET_NODES(1:NPOIN)
              !
              !-----------------------------------------------------------

              if( coupling_type(icoup) % where_type == ON_SET ) then
                 !
                 ! Coupling on a set
                 !
                 set_target  = coupling_type(icoup) % where_number
                 zone_target = coupling_type(icoup) % zone_target
                 nboun_wet = 0
                 !do kboun = 1,nbouz(lzone(zone_target))
                 !   iboun = lbouz(lzone(zone_target)) % l(kboun)
                 do iboun = 1,nboun
                    if( lbset(iboun) == set_target ) then
                       nboun_wet = nboun_wet + 1
                    end if
                 end do
                 coupling_type(icoup) % wet % nboun_wet = nboun_wet
                 call memory_alloca(memor_cou,'COUPLING % WET % LBOUN_WET','cou_initialize_coupling',coupling_type(icoup) % wet % lboun_wet,nboun_wet)

                 nboun_wet = 0
                 !do kboun = 1,nbouz(lzone(zone_target))
                 !   iboun = lbouz(lzone(zone_target)) % l(kboun)
                 do iboun = 1,nboun
                    if( lbset(iboun) == set_target ) then
                       nboun_wet = nboun_wet + 1
                       coupling_type(icoup) % wet % lboun_wet(nboun_wet) = iboun
                       do inodb = 1,lnnob(iboun)
                          ipoin = lnodb_cou(inodb,iboun)
                          list_wet_nodes(ipoin) = 1
                       end do
                    end if
                 end do

              else if( coupling_type(icoup) % where_type == ON_MIRROR ) then
                 !
                 ! On mirror: at this stage, we do not know which are the source points of the mirror coupling
                 !
                 continue
                 
              else if( coupling_type(icoup) % where_type == ON_FIELD ) then
                 !
                 ! Coupling on a field
                 !
                 field_target = coupling_type(icoup) % where_number 
                 nboun_wet    = 0
                 do iboun = 1,nboun
                    if( int(xfiel(field_target) % a(1,iboun,1),ip) /= 0 ) then
                       nboun_wet = nboun_wet + 1
                    end if
                 end do
                 coupling_type(icoup) % wet % nboun_wet = nboun_wet
                 call memory_alloca(memor_cou,'COUPLING % WET % LBOUN_WET','cou_initialize_coupling',coupling_type(icoup) % wet % lboun_wet,nboun_wet)
                 if( kfl_field(2,field_target) == NBOUN_TYPE ) then
                    nboun_wet = 0
                    do iboun = 1,nboun
                       if( int(xfiel(field_target) % a(1,iboun,1),ip) /= 0 ) then
                          nboun_wet = nboun_wet + 1
                          coupling_type(icoup) % wet % lboun_wet(nboun_wet) = iboun
                          do inodb = 1,lnnob(iboun)
                             ipoin = lnodb_cou(inodb,iboun)
                             list_wet_nodes(ipoin) = 1
                          end do
                       end if
                    end do
                 end if

              else if( coupling_type(icoup) % where_type == ON_CODE ) then
                 !
                 ! Coupling on a boundary code
                 !
                 bcode_target = coupling_type(icoup) % where_number 
                 nboun_wet    = 0
                 do iboun = 1,nboun
                    if( kfl_codbo(iboun) == bcode_target ) then
                       nboun_wet = nboun_wet + 1
                    end if
                 end do
                 coupling_type(icoup) % wet % nboun_wet = nboun_wet
                 call memory_alloca(memor_cou,'COUPLING % WET % LBOUN_WET','cou_initialize_coupling',coupling_type(icoup) % wet % lboun_wet,nboun_wet)

                 nboun_wet = 0
                 do iboun = 1,nboun
                    if( kfl_codbo(iboun) == bcode_target ) then
                       nboun_wet = nboun_wet + 1
                       coupling_type(icoup) % wet % lboun_wet(nboun_wet) = iboun
                       do inodb = 1,lnnob(iboun)
                          ipoin = lnodb_cou(inodb,iboun)
                          list_wet_nodes(ipoin) = 1
                       end do
                    end if
                 end do

              else if( coupling_type(icoup) % where_type == ON_WHOLE_MESH ) then
                 !
                 ! Coupling on the whole mesh
                 !
                 do ipoin = 1,npoin
                    list_wet_nodes(ipoin) = 1
                 end do

              else if( coupling_type(icoup) % where_type == ON_CHIMERA_MESH ) then
                 !
                 ! Chimera
                 !
                 do ipoin = 1,npoin
                    if( lnoch(ipoin) == NOFRI ) list_wet_nodes(ipoin) = 1
                 end do
                 nboun_wet = 0
                 do iboun = 1,nboun_cou
                    if( lboch_cou(iboun) == BOFRI ) nboun_wet = nboun_wet + 1
                 end do
                 coupling_type(icoup) % wet % nboun_wet = nboun_wet
                 call memory_alloca(memor_cou,'COUPLING % WET % LBOUN_WET','cou_initialize_coupling',coupling_type(icoup) % wet % lboun_wet,nboun_wet)

                 nboun_wet = 0
                 do iboun = 1,nboun_cou
                    if( lboch_cou(iboun) == BOFRI ) then
                       nboun_wet = nboun_wet + 1
                       coupling_type(icoup) % wet % lboun_wet(nboun_wet) = iboun
                       do inodb = 1,lnnob_cou(iboun)
                          ipoin = lnodb_cou(inodb,iboun)
                          list_wet_nodes(ipoin) = 1
                       end do
                    end if
                 end do

              else if( coupling_type(icoup) % where_type == ON_IMMERSED_MESH .or. &
                     & coupling_type(icoup) % where_type == ON_EMBEDDED_MESH ) then

                 ! Target boundary code of immersed/embedded mesh
                 bcode_target = coupling_type(icoup) % where_number

                 ! Set all wet nodes to 1
                 do ipoin = 1,npoin
                    list_wet_nodes(ipoin) = 1
                 end do

                 ! Set overlapping nodes to 1 for elements containing boundary code
                 if ( bcode_target > 0 ) then
                    nboun_wet = 0
                    do iboun = 1,nboun
                       if( kfl_codbo(iboun) == bcode_target ) then
                          nboun_wet = nboun_wet + 1
                       end if
                    end do
                    coupling_type(icoup) % wet % nboun_wet = nboun_wet
                    call memory_alloca(memor_cou,'COUPLING % WET % LBOUN_WET','cou_initialize_coupling',coupling_type(icoup) % wet % lboun_wet,nboun_wet)
                    nboun_wet = 0
                    do iboun = 1,nboun
                       if( kfl_codbo(iboun) == bcode_target ) then
                          nboun_wet = nboun_wet + 1
                          coupling_type(icoup) % wet % lboun_wet(nboun_wet) = iboun
                          do inodb = 1,lnnob(iboun)
                             ipoin = lnodb_cou(inodb,iboun)
                             list_overlap_nodes(ipoin) = 1
                          end do
                       end if
                    end do
                 end if


              end if
              !
              ! LIST_WET_NODES: Exchange list of wet nodes
              !
              call PAR_INTERFACE_NODE_EXCHANGE(list_wet_nodes,    'MAX','IN CURRENT TARGET COLOR')
              call PAR_INTERFACE_NODE_EXCHANGE(list_overlap_nodes,'MAX','IN CURRENT TARGET COLOR')
              !
              ! NUMBER_WET_NODES: compute number of wet nodes
              ! NUMBER_OVERLAP_NODES: compute number of wet nodes
              !
              number_wet_nodes     = 0
              number_overlap_nodes = 0
              do ipoin = 1,npoin
                 ! Compute number of wet nodes
                 if( list_wet_nodes(ipoin) /= 0 ) then
                    number_wet_nodes = number_wet_nodes + 1
                    list_wet_nodes(ipoin) = number_wet_nodes
                 end if
                 ! Compute number of overlap nodes
                 if( list_overlap_nodes(ipoin) /= 0 ) then
                    number_overlap_nodes = number_overlap_nodes + 1
                    list_overlap_nodes(ipoin) = number_overlap_nodes
                 end if
              end do
              !
              ! Allocate memory
              !
              call memory_alloca(memor_cou,'KWEIGHT_WET'                    ,'cou_initialize_coupling',kweight_wet,npoin)
              call memory_alloca(memor_cou,'COUPLING % WET % LPOIN_WET'     ,'cou_initialize_coupling',coupling_type(icoup) % wet % lpoin_wet     ,number_wet_nodes)
              call memory_alloca(memor_cou,'COUPLING % WET % WEIGHT_WET'    ,'cou_initialize_coupling',coupling_type(icoup) % wet % weight_wet    ,number_wet_nodes) 
              call memory_alloca(memor_cou,'COUPLING % WET % WEIGHT_WET_IMP','cou_initialize_coupling',coupling_type(icoup) % wet % weight_wet_imp,number_wet_nodes) 
              call memory_alloca(memor_cou,'COUPLING % WET % LOLAP_WET'     ,'cou_initialize_coupling',coupling_type(icoup) % wet % lolap_wet     ,number_wet_nodes)
              !
              ! LPOIN_WET(1:NUMBER_WET_NODES): Permutation array 
              ! LOLAP_WET(1:NUMBER_WET_NODES): Overlapping nodes array wrt LPOIN_WET
              ! KWEIGHT_WET: weights
              !
              number_wet_nodes = 0
              do ipoin = 1,npoin
                 if( list_wet_nodes(ipoin) /= 0 ) then
                    number_wet_nodes   = list_wet_nodes(ipoin)
                    kweight_wet(ipoin) = 1
                    coupling_type(icoup) % wet % lpoin_wet(number_wet_nodes) = ipoin
                    coupling_type(icoup) % wet % lolap_wet(number_wet_nodes) = list_overlap_nodes(ipoin)
                 end if
              end do
              coupling_type(icoup) % wet % npoin_wet = number_wet_nodes
              coupling_type(icoup) % wet % nolap_wet = number_overlap_nodes


              !----------------------------------------------------------
              !
              ! Get wet points (can be nodes, gauss points or knots)
              !
              ! COORD_WET(:,:) ... Define wet points coordinates
              ! NPOIN_WET ........ Number of wet points             
              !
              !----------------------------------------------------------

              if(    coupling_type(icoup) % itype == STRESS_PROJECTION .or. &
                   & coupling_type(icoup) % itype == PROJECTION        ) then
                 !
                 ! Wet points are Gauss points
                 !
                 call COU_PROJECTION_TYPE(coupling_type(icoup))
                 number_wet_points = 0
                 do iboun_wet = 1,coupling_type(icoup) % wet % nboun_wet
                    number_wet_points = number_wet_points &
                         + coupling_type(icoup) % wet % proje_target(iboun_wet) % pgaub
                 end do
                 coupling_type(icoup) % wet % number_wet_points = number_wet_points          
                 call memory_alloca(memor_cou,'COUPLING % WET % COORD_WET','cou_initialize_coupling',&
                      coupling_type(icoup) % wet % coord_wet,ndime,number_wet_points)
                 call memory_alloca(memor_cou,'COUPLING % WET % LSUBD'    ,'cou_initialize_coupling',&
                      coupling_type(icoup) % wet % lsubd,number_wet_points) 
                 number_wet_points = 0
                 do iboun_wet = 1,coupling_type(icoup) % wet % nboun_wet
                    iboun = coupling_type(icoup) % wet % lboun_wet(iboun_wet)
                    ielem = lelbo(iboun)
                    ngaub = coupling_type(icoup) % wet % proje_target(iboun_wet) % pgaub
                    pblty = abs(ltypb_cou(iboun))
                    pnodb = lnnob_cou(iboun)
                    do igaub = 1,ngaub
                       number_wet_points = number_wet_points + 1
                       coupling_type(icoup) % wet % coord_wet(1:ndime,number_wet_points) = 0.0_rp
                       coupling_type(icoup) % wet % lsubd(number_wet_points) = lesub(ielem)
                       do inodb = 1,pnodb
                          ipoin = lnodb_cou(inodb,iboun)
                          coupling_type(icoup) % wet % coord_wet(1:ndime,number_wet_points) = &
                               &   coupling_type(icoup) % wet % coord_wet(1:ndime,number_wet_points) &
                               & + coupling_type(icoup) % wet % proje_target(iboun_wet) % shapb(inodb,igaub) &
                               & * coord(1:ndime,ipoin)
                       end do
                    end do
                 end do
              else
                 !
                 ! Wet points are nodes
                 !
                 number_wet_points = number_wet_nodes
                 coupling_type(icoup) % wet % number_wet_points = number_wet_points
                 call memory_alloca(memor_cou,'COUPLING % WET % COORD_WET','cou_initialize_coupling',&
                      coupling_type(icoup) % wet % coord_wet,ndime,number_wet_points)   
                 call memory_alloca(memor_cou,'COUPLING % WET % LSUBD'    ,'cou_initialize_coupling',&
                      coupling_type(icoup) % wet % lsubd,number_wet_points) 
                 number_wet_points = 0
                 call memgen(1_ip,npoin,0_ip)
                 do ielem = 1,nelem
                    do inode = 1,element_type(abs(ltype(ielem))) % number_nodes
                       ipoin = lnods(inode,ielem)
                       gisca(ipoin) = lesub(ielem)
                    end do                    
                 end do
                 do ipoin = 1,npoin
                    if( list_wet_nodes(ipoin) /= 0 ) then
                       number_wet_points = list_wet_nodes(ipoin)
                       coupling_type(icoup) % wet % lsubd(number_wet_points) = gisca(ipoin)
                       coupling_type(icoup) % wet % coord_wet(1:ndime,number_wet_points) &
                            = coord(1:ndime,ipoin)
                    end if
                 end do
                 call memgen(3_ip,npoin,0_ip)
              end if
              !
              ! Copy initial conditions
              !
              call memory_copy(memor_cou,'COUPLING % WET % COORD_WET_INI','cou_initialize_coupling', coupling_type(icoup) % wet % coord_wet,coupling_type(icoup) % wet % coord_wet_ini,'DO_NOT_DEALLOCATE')
              !
              ! Define weights of wet nodes             
              !
              call PAR_INTERFACE_NODE_EXCHANGE(kweight_wet,'SUM','IN CURRENT TARGET COLOR')
              if( coupling_type(icoup) % kind == BETWEEN_SUBDOMAINS ) then
                 !
                 ! Subdomain coupling: weight is one, because coupling is carried out AFTER exchange
                 !
                 do kpoin = 1,coupling_type(icoup) % wet % npoin_wet
                    coupling_type(icoup) % wet % weight_wet(kpoin) = 1.0_rp
                 end do
              else
                 !
                 ! Weight depends on the number of neighbors, because coupling is carried out BEFORE exchange
                 !   
                 do kpoin = 1,coupling_type(icoup) % wet % npoin_wet
                    ipoin = coupling_type(icoup) % wet % lpoin_wet(kpoin)
                    coupling_type(icoup) % wet % weight_wet(kpoin) = 1.0_rp / real(kweight_wet(ipoin),rp)
                 end do
              end if
              !
              ! Weight used when transmision matrices are not parallelized
              !
              do kpoin = 1,coupling_type(icoup) % wet % npoin_wet
                 ipoin = coupling_type(icoup) % wet % lpoin_wet(kpoin)
                 coupling_type(icoup) % wet % weight_wet_imp(kpoin) = 1.0_rp / real(kweight_wet(ipoin),rp)
              end do
              call memory_deallo(par_memor,'KWEIGHT_WET','cou_initialize_coupling',kweight_wet)
              !
              ! LIST_WET_NODES: reinitialize
              !
              do ipoin = 1,npoin
                 list_wet_nodes(ipoin) = 0
              end do
           end if

           call cputim(time6)
           coupling_type(icoup) % cputim(1) = coupling_type(icoup) % cputim(1) + time6 - time5
        end do
     end do

     if ( INOTMASTER ) then
        call memory_deallo(memor_cou,'LIST_WET_NODES',    'cou_initialize_coupling',list_wet_nodes)
        call memory_deallo(memor_cou,'LIST_OVERLAP_NODES','cou_initialize_coupling',list_overlap_nodes)
     end if

  end if

end subroutine cou_define_wet_geometry

subroutine cou_trapezoidal_rule(ndime,ngaus,posgp,weigp) 
  use def_kintyp, only : ip,rp
  implicit none
  integer(ip), intent(in)            :: ndime
  integer(ip), intent(in)            :: ngaus
  real(rp),    intent(out)           :: posgp(ndime,ngaus)
  real(rp),    intent(out)           :: weigp(ngaus)
  !integer(ip), intent(out), optional :: ierro
  real(rp)                           :: delta,xpanels
  integer(ip)                        :: igaus,npanels
  
  !if( present(ierro) ) ierro = 0
  
  if( ndime == 1 ) then
     npanels = ngaus - 1
     xpanels = 1.0_rp / real(npanels,rp)
     delta   = 2.0_rp / real(ngaus-1_ip,rp)
     do igaus = 1,ngaus
        posgp(1,igaus) = -1.0_rp + delta * real(igaus-1_ip,rp)
        weigp(igaus)   =  2.0_rp * xpanels
     end do
     weigp(1)     = xpanels
     weigp(ngaus) = xpanels
  else
     call runend('NOT CODED')
  end if

end subroutine cou_trapezoidal_rule
