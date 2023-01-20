!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Partis
!> @{
!> @file    pts_moving_mesh.f90
!> @author  houzeaux
!> @date    2020-07-06
!> @brief   Moving mesh
!> @details Recompute host elements if necessary
!> @} 
!-----------------------------------------------------------------------

subroutine pts_moving_mesh()

  use def_master
  use def_domain
  use def_partis
  use mod_elmgeo
  use def_kermod
  use mod_memory
  use mod_messages,                      only : messages_live
  use def_kintyp_comm,                   only : comm_data_par
  use mod_pts_host_element,              only : pts_host_element
  use mod_par_bin_structure,             only : par_bin_structure_points_in_partitions
  use mod_communications_global,         only : PAR_MAX
  use mod_communications_global,         only : PAR_SUM
  use mod_communications_global,         only : PAR_ALLTOALL
  use mod_communications_point_to_point, only : PAR_SEND_RECEIVE
  use mod_pts_parallelization,           only : pts_parallelization_migration
  use mod_strings,                       only : integer_to_string
  use mod_parall,                        only : typ_bin_structure
  use mod_parall,                        only : PAR_COMM_MY_CODE 
  use mod_par_bin_structure,             only : par_bin_structure
  use mod_par_bin_structure,             only : par_bin_structure_partition_bounding_box
  use mod_par_bin_structure,             only : par_bin_structure_deallocate
  use mod_par_bin_structure,             only : par_bin_structure_initialization
  implicit none

  integer(ip)              :: ilagr,ifoun,pelty
  integer(ip)              :: pnode,ielem,nlost
  integer(ip)              :: ipoin,inode,ilost
  integer(ip)              :: nfoun,ii,nn,cpart
  integer(ip)              :: pp,jj,iimin(1),ineig
  real(rp)                 :: elcod(ndime,mnode)
  real(rp)                 :: shapf(mnode)
  real(rp)                 :: deriv(ndime,mnode)
  real(rp)                 :: coloc(3)
  character(20), parameter :: vacal = 'pts_moving_mesh'
  type(comm_data_par)      :: comm
  
  type(typ_bin_structure)  :: bin_structure
  real(rp)                 :: comin(3),comax(3),delta(3)
  
  integer(ip),   pointer   :: list_lost(:)
  real(rp),      pointer   :: xx_lost(:,:)
  real(rp),      pointer   :: host_shapf(:,:)          
  integer(ip),   pointer   :: host_element(:)
  integer(ip),   pointer   :: npoin_in_part(:)
  integer(ip),   pointer   :: npoin_recv(:)
  type(i1p),     pointer   :: lpoin_in_part(:)
  type(i1p),     pointer   :: decision_send(:)
  type(i1p),     pointer   :: decision_recv(:)
  type(r2p),     pointer   :: xx_send(:)
  type(r2p),     pointer   :: xx_recv(:) 
  real(rp),      pointer   :: xx_shapf(:,:)       ! => used to point
  integer(ip),   pointer   :: xx_element(:)       ! => used to point
  integer(ip),   pointer   :: list_owner(:,:)
  integer(ip),   pointer   :: nlagr_migrating(:)  ! Number of particles to send
  integer(ip),   pointer   :: nlagr_recv(:)       ! Number of particles to send

  nullify(list_lost)
  nullify(xx_lost)
  nullify(host_shapf)
  nullify(host_element) 
  nullify(npoin_in_part)
  nullify(npoin_recv)
  nullify(lpoin_in_part)
  nullify(decision_send)
  nullify(decision_recv)
  nullify(xx_send)
  nullify(xx_recv)
  nullify(xx_shapf)
  nullify(xx_element)
  nullify(list_owner)
  nullify(nlagr_migrating)
  nullify(nlagr_recv)

  nlost               = 0
  nfoun               = 0
  
  if( if_moving_mesh_pts ) then
     !
     ! List of existing particles, permutation array permu_nlagr
     !
     call pts_compute_permutation()
     !
     ! Check bin element searh 
     !
     if (kfl_usbin_pts == 1 ) call runend('PARTIS: BECAUSE OF MOVING_MESH BIN ELEMENT SEARCH STRATEGY IS NOT AVAILABLE')
     !
     ! Check if this is not a multicode simulation
     !
     if( kfl_modul(ID_PARTIS) == 0 .and. kfl_difun == 0 ) &
          call runend('PTS_MOVING_MESH: PARTIS AND ALEFOR SHOULD BE COMPUTED IN THE SAME ALYA CODE... BECAUSE OF BIN STRUCTURE')
     !
     ! List of lost particles: LIST_LOST
     !
     call memory_alloca(mem_modul(1:2,modul),'LIST_LOST',vacal,list_lost,mlagr)
     nlost = 0
     do ilagr = 1,mlagr
        if( lagrtyp(ilagr) % kfl_exist == -1 ) then
           ielem = lagrtyp(ilagr) % ielem
           pelty = ltype(ielem)
           pnode = element_type(pelty) % number_nodes
           do inode = 1,pnode
              ipoin = lnods(inode,ielem)
              elcod(1:ndime,inode) = coord(1:ndime,ipoin)
           end do
           call elmgeo_natural_coordinates(             &
                ndime,pelty,pnode,elcod,shapf,deriv,    &
                lagrtyp(ilagr) % coord,coloc,ifoun,relse(1))
           if( ifoun == 0 ) then
              nlost            = nlost + 1
              list_lost(nlost) = ilagr
           end if
        end if
     end do
     call memory_alloca(mem_modul(1:2,modul),'XX_LOST',vacal,xx_lost,ndime,nlost)
     do ilost = 1,nlost
        ilagr = list_lost(ilost)
        xx_lost(1:ndime,ilost) = lagrtyp(ilagr) % coord(1:ndime)
     end do
     ii = nlost
     call PAR_MAX(ii)     
     !
     ! Check in my mesh first
     !
     if( ii > 0 ) then
        call memory_alloca(mem_modul(1:2,modul),'HOST_ELEMENT'  ,vacal,host_element  ,nlost)
        call memory_alloca(mem_modul(1:2,modul),'HOST_SHAPF'    ,vacal,host_shapf    ,mnode,nlost)  
        call pts_host_element(&
             nlost,nfoun,xx_lost,&
             host_shapf,host_element)
        ii = 0
        do ilost = 1,nlost
           if( host_element(ilost) <= 0 ) then 
              ii = ii + 1
              xx_lost(1:ndime,ii)    = xx_lost(1:ndime,ilost)
              list_lost(ii)          = list_lost(ilost)
           else
              ilagr                  = list_lost(ilost)
              lagrtyp(ilagr) % ielem = host_element(ilost)
           end if
        end do
        nlost = ii
        call memory_deallo(mem_modul(1:2,modul),'HOST_ELEMENT'  ,vacal,host_element)
        call memory_deallo(mem_modul(1:2,modul),'HOST_SHAPF'    ,vacal,host_shapf)  
        !
        ! Still have lost particles!
        !
        call PAR_MAX(ii)

        if( ii > 0 ) then 

           call par_bin_structure_initialization        (bin_structure)
           call par_bin_structure_partition_bounding_box(comin,comax,delta)
           call par_bin_structure                       (bin_structure,comin,comax,COMM=PAR_COMM_MY_CODE,VERBOSE=.false.)       
           call par_bin_structure_points_in_partitions  (xx_lost,npoin_in_part,lpoin_in_part,NUMBER_POINTS=nlost,BIN_STRUCTURE=bin_structure)
           call par_bin_structure_deallocate            (bin_structure)

           call memory_alloca(mem_modul(1:2,modul),'XX_SEND'      ,vacal,xx_send      ,npart)
           call memory_alloca(mem_modul(1:2,modul),'NPOIN_RECV'   ,vacal,npoin_recv   ,npart+1,lboun=0_ip)
           call memory_alloca(mem_modul(1:2,modul),'DECISION_RECV',vacal,decision_recv,npart)
           do cpart = 1,npart
              nn = memory_size(lpoin_in_part(cpart) % l)
              call memory_alloca(mem_modul(1:2,modul),'XX_SEND % A'      ,vacal,xx_send(cpart) % a,ndime,nn)
              call memory_alloca(mem_modul(1:2,modul),'DECISION_RECV % L',vacal,decision_recv(cpart) % l,nn)
              do ii = 1,nn
                 ilost = lpoin_in_part(cpart) % l(ii)
                 xx_send(cpart) % a(1:ndime,ii) = xx_lost(1:ndime,ilost)
              end do
           end do
           !
           ! Send particles to selected neighbors
           !
           call PAR_ALLTOALL(1_ip,1_ip,npoin_in_part,npoin_recv)
           !
           ! Allocate for receiving particles
           !
           call memory_alloca(mem_modul(1:2,modul),'XX_RECV'      ,vacal,xx_recv      ,npart)
           call memory_alloca(mem_modul(1:2,modul),'DECISION_SEND',vacal,decision_send,npart)
           do cpart = 1,npart
              call memory_alloca(mem_modul(1:2,modul),'XX_RECV       % A',vacal,xx_recv(cpart) % a      ,ndime,npoin_recv(cpart))
              call memory_alloca(mem_modul(1:2,modul),'DECISION_SEND % L',vacal,decision_send(cpart) % l,npoin_recv(cpart))
              call PAR_SEND_RECEIVE(xx_send(cpart) % a,xx_recv(cpart) % a,'IN MY CODE',cpart,'BLOCKING')
           end do
           !
           ! Check if I have the NN particles I received
           !
           nn = sum(npoin_recv)
           if( nn > 0 ) then
              call memory_alloca(mem_modul(1:2,modul),'HOST_ELEMENT'  ,vacal,host_element  ,nn)
              call memory_alloca(mem_modul(1:2,modul),'HOST_SHAPF'    ,vacal,host_shapf    ,mnode,nn)  
              ii = 1
              jj = 0
              do cpart = 1,npart
                 if( npoin_recv(cpart) > 0 ) then
                    xx_shapf   => host_shapf(1:,ii:)
                    xx_element => host_element(ii:)
                    call pts_host_element(&
                         npoin_recv(cpart),nfoun,xx_recv(cpart) % a,&
                         xx_shapf,xx_element)
                    ii = ii + npoin_recv(cpart)
                 end if
                 do pp = 1,npoin_recv(cpart)
                    jj = jj + 1
                    decision_send(cpart) % l(pp) = host_element(jj)
                 end do
              end do
           end if
           do cpart = 1,npart
              call PAR_SEND_RECEIVE(decision_send(cpart) % l,decision_recv(cpart) % l,'IN MY CODE',cpart,'BLOCKING')
           end do
           !
           ! Choose the first one
           !           
           call memory_alloca(mem_modul(1:2,modul),'LIST_OWNER'     ,vacal,list_owner     ,2_ip,nlost)
           call memory_alloca(mem_modul(1:2,modul),'NLAGR_MIGRATING',vacal,nlagr_migrating,npart+1,lboun=0_ip)
           call memory_alloca(mem_modul(1:2,modul),'NLAGR_RECV'     ,vacal,nlagr_recv     ,npart+1,lboun=0_ip)
           jj = 0
           do cpart = 1,npart
              do ii = 1,memory_size(decision_recv(cpart) % l)
                 pp    = decision_recv(cpart) % l(ii)
                 ilost = lpoin_in_part(cpart) % l(ii)
                 ilagr = list_lost(ilost)
                 if( pp /= 0 ) then
                    if( list_owner(1,ilost) == 0 ) then
                       nlagr_migrating(cpart)     = nlagr_migrating(cpart) + 1
                       list_owner(1,ilost)        = cpart
                       list_owner(2,ilost)        = ii
                       lagrtyp(ilagr) % kfl_exist = cpart
                    end if
                 end if
              end do
           end do
           !
           ! Check lost particles
           !
           nlagr_going_out_pts = 0
           kfl_lost_moving_pts = 0
           do ilost = 1,nlost
              if( list_owner(1,ilost) == 0 ) then
                 nlagr_going_out_pts        = nlagr_going_out_pts + 1
                 ilagr                      = list_lost(ilost)
                 lagrtyp(ilagr) % kfl_exist = PTS_PARTICLE_MOVING_MESH 
                 !print*,'kk-------corrd-lost',lagrtyp(ilagr) % coord(1:ndime)
                 kfl_lost_moving_pts        = 1
              end if
           end do
           ii = nlagr_going_out_pts
           call PAR_SUM(ii)
           if( ii > 0 ) then 
              call messages_live(&
                   'PARTICLES GOING OUT OF COMPUTATIONAL DOMAIN DUE TO MESH MOTION= '//integer_to_string(ii),&
                   'MODULE')
           end if
           !
           ! Tell subdomain he has to add me has a neighbor
           !
           call PAR_ALLTOALL(1_ip,1_ip,nlagr_migrating,nlagr_recv)
           !
           ! Migrate particles
           !
           allocate(comm % neights(npart))
           comm % nneig = 0
           do cpart = 1,npart
              if( nlagr_migrating(cpart) /= 0 ) then
                 comm % nneig                  = comm % nneig + 1
                 comm % neights(comm % nneig)  = cpart
                 nlagr_migrating(comm % nneig) = nlagr_migrating(cpart)
              else if( nlagr_recv(cpart) /= 0 ) then
                 comm % nneig                  = comm % nneig + 1
                 comm % neights(comm % nneig)  = cpart                 
              end if
           end do
           do ilost = 1,nlost
              if( list_owner(1,ilost) /= 0 ) then
                 cpart                      = list_owner(1,ilost)
                 ilagr                      = list_lost(ilost)
                 iimin                      = minloc(comm % neights,comm % neights==cpart)
                 ineig                      = iimin(1)
                 lagrtyp(ilagr) % kfl_exist = ineig
                 lagrtyp(ilagr) % ielem     = decision_recv(cpart) % l(list_owner(2,ilost))
              end if
           end do

           call pts_parallelization_migration(comm,nlagr_migrating(1:))  

           deallocate(comm % neights)

           nlost = nlost - nlagr_going_out_pts
           call PAR_SUM(nlost)
           if( nlost > 0 ) & 
                call messages_live('MIGRATING '//integer_to_string(nlost)//' PARTICLES DUE TO MESH MOTION','MODULE')
        end if
     end if
     !
     ! Deallocate everything
     !
     call memory_deallo(mem_modul(1:2,modul),'LIST_LOST'      ,vacal,list_lost)
     call memory_deallo(mem_modul(1:2,modul),'XX_LOST'        ,vacal,xx_lost)
     call memory_deallo(mem_modul(1:2,modul),'XX_RECV'        ,vacal,xx_recv)
     call memory_deallo(mem_modul(1:2,modul),'HOST_SHAPF'     ,vacal,host_shapf)  
     call memory_deallo(mem_modul(1:2,modul),'HOST_ELEMENT'   ,vacal,host_element)
     call memory_deallo(mem_modul(1:2,modul),'NPOIN_IN_PART'  ,vacal,npoin_in_part)
     call memory_deallo(mem_modul(1:2,modul),'NPOIN_RECV'     ,vacal,npoin_recv)
     call memory_deallo(mem_modul(1:2,modul),'LPOIN_IN_PART'  ,vacal,lpoin_in_part)
     call memory_deallo(mem_modul(1:2,modul),'DECISION_SEND'  ,vacal,decision_send)
     call memory_deallo(mem_modul(1:2,modul),'DECISION_RECV'  ,vacal,decision_recv)
     call memory_deallo(mem_modul(1:2,modul),'XX_SEND'        ,vacal,xx_send)
     call memory_deallo(mem_modul(1:2,modul),'XX_RECV'        ,vacal,xx_recv)
     call memory_deallo(mem_modul(1:2,modul),'LIST_OWNER'     ,vacal,list_owner)
     call memory_deallo(mem_modul(1:2,modul),'NLAGR_MIGRATING',vacal,nlagr_migrating)
     call memory_deallo(mem_modul(1:2,modul),'NLAGR_RECV'     ,vacal,nlagr_recv)

  end if

end subroutine pts_moving_mesh
