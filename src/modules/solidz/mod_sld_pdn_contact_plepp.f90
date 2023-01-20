!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    mod_sld_pdn_contact_plepp.f90
!> @author  gguillamet
!> @date    2022-05-25
!> @brief   PDN contact algorithm for rigid-to-deformable 
!> @details PLE++ version for Explicit Analysis  
!>
!> @{
!-----------------------------------------------------------------------

module mod_sld_pdn_contact_plepp

  use def_kintyp_basic,      only : ip, rp, lg
  use def_kintyp,            only : soltyp
  use def_master,            only : INOTMASTER,ISEQUEN
  use def_master,            only : momod, intost
  use def_master,            only : mem_modul, modul
  use def_domain,            only : ndime, npoin
  use mod_memory,            only : memory_alloca
  use mod_memory,            only : memory_deallo
  use mod_communications,    only : PAR_SUM, PAR_MAX, PAR_MIN
  use mod_communications,    only : PAR_BROADCAST
  use mod_messages,          only : messages_live
#if defined COMMDOM && COMMDOM == 2
  use mod_plepp_pdn_contact, only : CNT_CPLNG
  use mod_plepp_pdn_contact, only : PLEPP_CPLNG
  use mod_plepp_pdn_contact, only : CNT_SENDRECV
#endif
  use def_solidz,            only : kfl_fixno_sld, kfl_fixrs_sld
  use def_solidz,            only : kfl_goite_sld
  use def_solidz,            only : frxid_sld, fcont_sld
  use def_solidz,            only : rbfor_sld
  use def_solidz,            only : jacrot_du_dq_sld, jacrot_dq_du_sld
  use def_solidz,            only : bvess_sld

  implicit none
  
  !----------------------------------------------------------------------
  !
  ! Parameters
  !
  !----------------------------------------------------------------------

  type(soltyp), pointer   :: solve(:)
  real(rp)                :: force_dbo(3)
#if defined COMMDOM && COMMDOM == 2
  integer(ip),  private   :: MPI_RANK_MASTER_DEF
  integer(ip),  private   :: MPI_RANK_MASTER_RBO
#endif
  integer(ip),  parameter :: SLD_PDN_RIGID_TO_DEFOR = 4_ip
  integer(ip),  parameter :: lun_pdnco_sld          = 118_ip

  logical(lg)                         :: &
       kfl_relea_sld,                    & ! PDN contact release node flag
       kfl_pdncf_sld                       ! PDN contact file convergence flag

  integer(ip)                         :: &
       kfl_pdnco_sld,                    & ! PDN Contact type flag
       bcont_sld,                        & ! Identificator for contact boundary
       itrel_sld,                        & ! Sub-iterations for release nodes
       ntoco_sld,                        & ! Total contact nodes counter (artificial + contact)
       nacco_sld,                        & ! Active contact nodes counter (detected)
       narco_sld                           ! Artificial contact nodes counter

  real(rp)                            :: &
       tolpe_sld                           ! Penetration tolerance

  integer(ip), pointer                :: &
       kfl_rnode_sld(:)                    ! List of marked nodes to be released

#if defined COMMDOM && COMMDOM == 2
  
  public  :: sld_pdn_contact_cvgunk
  public  :: sld_pdn_contact_reabcs
  public  :: sld_pdn_contact_parall
  public  :: sld_pdn_contact_begste
  public  :: sld_pdn_contact_plugin
  public  :: sld_pdn_contact_releas

  private :: sld_pdn_contact_proje_send
  private :: sld_pdn_contact_proje_recv
  private :: sld_pdn_contact_force_send
  private :: sld_pdn_contact_projection
  private :: sld_pdn_contact_master_rank

#endif
  
contains
  
#if defined COMMDOM && COMMDOM == 2
  
  !-----------------------------------------------------------------------
  !>
  !> @author  Miguel Zavala and Gerard Guillamet
  !> @date    March, 2019
  !> @brief   Main subroutine for PDN contact coupling (Generalized)
  !> @details Main subroutine for PDN contact coupling (Generalized)
  !>
  !-----------------------------------------------------------------------

  subroutine sld_pdn_contact_plugin

    implicit none

    integer(ip)          :: nsend            !< Localized nodes sent
    integer(ip)          :: nrecv            !< Localized nodes recieved
    integer(ip)          :: ndata            !< Dimensions data
    integer(ip)          :: kpoin
    real(rp),    pointer :: data_send(:,:) 
    real(rp),    pointer :: data_recv(:,:)
    logical(lg), save    :: first_time = .true.

    nullify(data_send)
    nullify(data_recv)

    if( ISEQUEN ) call runend('SLD_PDN_CONTACT_PLUGIN: SEND FORCE IS ONLY PROGRAMMED FOR PARALL')
    
    if( any(CNT_SENDRECV) ) then
       !
       ! Master ranks
       !
       if( first_time ) then
          call sld_pdn_contact_master_rank()
          first_time = .false.
       end if
       !
       ! Data from PLEPP
       !
       nsend = PLEPP_CPLNG%n_ij ! Nodos localizados del deformable respeto al rigido
       nrecv = PLEPP_CPLNG%n_ji ! Nodos localizados del rigido respeto al deformable
       !
       ! Dimensions
       !
       if( ndime == 2_ip ) then
          ndata = 5_ip
       else if( ndime == 3_ip ) then
          ndata = 10_ip
       end if
       !
       ! Allocate memory
       !
       allocate(data_send(ndata,nsend))
       allocate(data_recv(ndata,nrecv))
       !
       ! Initialize
       !
       do kpoin = 1,nsend
          data_send(1:ndata,kpoin) = 0.0_rp
       end do
       do kpoin = 1,nrecv
          data_recv(1:ndata,kpoin) = 0.0_rp
       end do

       code_i: if(CNT_CPLNG%current_code == CNT_CPLNG%code_i) then

          !------------------------------------------------------------------
          !
          ! Rigid body (code_i==1)
          !
          !------------------------------------------------------| U--> |---!

          if( CNT_SENDRECV(7) ) then
             !
             ! Projections: Dirichlet
             !
             call sld_pdn_contact_proje_send(ndata,nsend,data_send(:,:))
             !
             ! Exchange: SEND
             !
             call commdom_locator_exchange_double_stride(data_send(:,:),data_recv(:,:),ndata)
          end if
          !
          ! Bilateral
          !---------------------------------------------------| dUdn<-- |---!
          if( CNT_SENDRECV(12) ) then
             !
             ! Compatibility checks
             !
             solve(1:) => momod(modul) % solve(1:)
             if(solve(1)%kfl_react==1) call runend("solve(1)%kfl_react==1")
             !
             ! Exchange: RECIEVE
             !
             call PAR_BROADCAST(ndime,force_dbo,wherein='IN THE UNIVERSE',root_rank=MPI_RANK_MASTER_DEF)
             rbfor_sld(:) = -force_dbo(:)
             call messages_live('<- RECV '//trim('FORCE')//' AS A TARGET FOR COUPLING '//trim(intost(1_ip)))
          end if

       end if code_i

       code_j: if(CNT_CPLNG%current_code == CNT_CPLNG%code_j) then

          !------------------------------------------------------------------
          !
          ! Deformable body (code_j==2)
          !
          !------------------------------------------------------| U<-- |---!

          if( CNT_SENDRECV(7) ) then
             !
             ! Exchange: RECIEVE
             !
             call commdom_locator_exchange_double_stride(data_send(:,:),data_recv(:,:),ndata)
             !
             ! Projections: Apply Dirichlet
             !
             call sld_pdn_contact_proje_recv(ndata,nrecv,data_recv(:,:))

          end if
          !
          ! Bilateral
          !---------------------------------------------------| dUdn--> |---!
          if( CNT_SENDRECV(12) ) then
             !
             ! Compatibility checks
             !
             solve(1_ip:) => momod(modul) % solve(1_ip:)
             if(solve(1)%kfl_bvnat==1) call runend("solve(1)%kfl_bvnat==1")  ! si es Neumaan sale!!
             !
             ! Force: Neumann
             !
             call sld_pdn_contact_force_send()
             !
             ! Exchange: SEND
             !
             call PAR_BROADCAST(ndime,force_dbo,wherein='IN THE UNIVERSE',root_rank=MPI_RANK_MASTER_DEF)
             call messages_live('-> SEND '//trim('FORCE')//' AS A SOURCE FOR COUPLING '//trim(intost(1_ip)))
             
          end if

       end if code_j
       !
       ! Deallocate memory arrays
       !
       deallocate(data_send)
       deallocate(data_recv)

    end if

  end subroutine sld_pdn_contact_plugin

  !-----------------------------------------------------------------------
  !>
  !> @author  gguillamet
  !> @date    2021-07-21
  !> @brief   PDN-contact send projections
  !> @details PDN-contact send projections
  !>
  !-----------------------------------------------------------------------

  subroutine sld_pdn_contact_proje_send(ndata,nsend,data_send)

    implicit none

    integer(ip), intent(in)    :: ndata
    integer(ip), intent(in)    :: nsend
    real(rp),    intent(inout) :: data_send(ndata,nsend) !< Data array to send

    integer(ip)                :: kpoin
    real(rp), pointer          :: local_basis(:,:,:)
    real(rp), pointer          :: norma_dista(:)

    !------------------------------------------------------------------
    !
    ! Body i
    !
    !------------------------------------------------------------------

    nullify(local_basis)
    nullify(norma_dista)
    !
    ! Allocate memory for arrays
    !
    call memory_alloca(mem_modul(1:2,modul),'LOCAL_BASIS','sld_pdn_contact_proje_send',local_basis,ndime,ndime,nsend)
    call memory_alloca(mem_modul(1:2,modul),'NORMA_DISTA','sld_pdn_contact_proje_send',norma_dista,nsend)
    !
    ! Get projections: Nodes body j which are inside body i are projected to boundary i
    !
    call sld_pdn_contact_projection(nsend,local_basis,norma_dista)
    !
    ! Store data
    !
    if(      ndime == 2_ip ) then
       do kpoin = 1,nsend
          data_send(1:2,kpoin) = local_basis(1:2,1,kpoin)
          data_send(3:4,kpoin) = local_basis(1:2,2,kpoin)
          data_send(  5,kpoin) = norma_dista(kpoin)
       end do
    else if( ndime == 3_ip ) then
       do kpoin = 1,nsend
          data_send(1:3,kpoin) = local_basis(1:3,1,kpoin)
          data_send(4:6,kpoin) = local_basis(1:3,2,kpoin)
          data_send(7:9,kpoin) = local_basis(1:3,3,kpoin)
          data_send( 10,kpoin) = norma_dista(kpoin)
       end do
    end if
    !
    ! Deallocate memory for arrays
    !
    call memory_deallo(mem_modul(1:2,modul),'LOCAL_BASIS','sld_pdn_contact_proje_send',local_basis)
    call memory_deallo(mem_modul(1:2,modul),'NORMA_DISTA','sld_pdn_contact_proje_send',norma_dista)

  end subroutine sld_pdn_contact_proje_send

  !-----------------------------------------------------------------------
  !>
  !> @author  gguillamet
  !> @date    2021-07-21
  !> @brief   PDN-contact projection recieve
  !> @details PDN-contact projection recieve
  !>
  !-----------------------------------------------------------------------

  subroutine sld_pdn_contact_proje_recv(ndata,nrecv,data_recv)

    use def_master,       only : displ, TIME_N, ITER_K
    use def_master,       only : kfl_rstar
    use def_domain,       only : kfl_codno, lpoty

    implicit none

    integer(ip),intent(in)    :: ndata
    integer(ip),intent(in)    :: nrecv
    real(rp),   intent(inout) :: data_recv(ndata,nrecv)

    integer(ip)               :: kpoin, ipoin, idime, ibopo

    !------------------------------------------------------------------
    !
    ! Deformable body
    !
    !------------------------------------------------------------------

    if( INOTMASTER ) then
       !
       ! Initializations
       !
       !if( kfl_rstar == 0 ) then
          do ipoin = 1,npoin
             do idime = 1,ndime
                if( kfl_fixno_sld(idime,ipoin) == 3_ip ) kfl_fixno_sld(idime,ipoin) = 0_ip
             end do
          end do
       !end if
       !
       ! Loop nodes recieved from PLEPP
       !
       do kpoin = 1,nrecv

          ipoin = PLEPP_CPLNG%interior_list_j(kpoin)
          ibopo = lpoty(ipoin)

          if( ibopo > 0 .and. any(kfl_codno(:,ipoin) == bcont_sld) ) then   ! nodes that belong to the boundary
             !
             ! Fill rotation matrix
             !
             jacrot_du_dq_sld(1:ndime,1,ipoin) = data_recv(1:ndime,kpoin)   ! Normal vector
             jacrot_du_dq_sld(1:ndime,2,ipoin) = data_recv(ndime+1:ndime*2,kpoin) ! Tangent vector
             if( ndime == 3_ip ) jacrot_du_dq_sld(1:ndime,3,ipoin) = data_recv(ndime*2+1:ndime*3,kpoin)
             ! Transpose of R
             jacrot_dq_du_sld(:,:,ipoin) = transpose(jacrot_du_dq_sld(:,:,ipoin))
             !
             ! Contact nodes with fixity in the first dof
             !
             if( kfl_fixno_sld(1,ipoin) == 1_ip ) then
                do idime = 1,ndime-1
                   kfl_fixno_sld(idime+1,ipoin) = 1_ip
                end do
             end if
             !
             ! Mark contact node
             !
             kfl_fixno_sld(1,ipoin) = 3_ip
             ntoco_sld              = ntoco_sld + 1
             nacco_sld              = nacco_sld + 1
             !
             ! Prescribed Dirichlet
             !
             !displ(1:ndime,ipoin,TIME_N) =  displ(1:ndime,ipoin,TIME_N) + & !! For Implicit scheme
             !        data_recv(ndime*ndime+1,kpoin)*data_recv(1:ndime,kpoin) - abs(tolpe_sld)*data_recv(1:ndime,kpoin)
             bvess_sld(1:ndime,ipoin,ITER_K) = displ(1:ndime,ipoin,TIME_N) + &
                  data_recv(ndime*ndime+1,kpoin)*data_recv(1:ndime,kpoin) - abs(tolpe_sld)*data_recv(1:ndime,kpoin)
             !bvess_sld(1:ndime,ipoin,ITER_K) = displ(1:ndime,ipoin,TIME_N) + &
             !     (abs(data_recv(ndime*ndime+1,kpoin)) - abs(tolpe_sld))*data_recv(1:ndime,kpoin)
          end if

       end do
       !
       ! Allocate memory for release nodes array
       !
       if( .not. associated(kfl_rnode_sld)) then
          call memory_alloca(mem_modul(1:2,modul),'KFL_RNODE_SLD','sld_pdn_conta_projection_recv',kfl_rnode_sld,nrecv)
          do kpoin = 1,nrecv
             kfl_rnode_sld(kpoin) = 0_ip
          end do
          kfl_relea_sld = .false.
       end if

    end if

  end subroutine sld_pdn_contact_proje_recv

  !-----------------------------------------------------------------------
  !>
  !> @author  gguillamet
  !> @date    2021-07-21
  !> @brief   PDN-contact send force
  !> @details PDN-contact send force
  !>
  !-----------------------------------------------------------------------

  subroutine sld_pdn_contact_force_send()

    use def_domain,         only : lpoty
    use mod_sld_csys,       only : sld_csys_rotuni

    implicit none

    integer(ip)       :: ipoin,kpoin,itott,ibopo
    integer(ip)       :: nrecv

    !------------------------------------------------------------------
    !
    ! Deformable body
    !
    !------------------------------------------------------------------
    !
    ! Rotate to global
    !
    if( INOTMASTER ) then
       nrecv = PLEPP_CPLNG%n_ji
       do kpoin = 1,nrecv
          ipoin = PLEPP_CPLNG%interior_list_j(kpoin)
          ibopo = lpoty(ipoin)
          itott = (ipoin - 1) * ndime
          if( ibopo > 0 .and. kfl_fixno_sld(1,ipoin) == 3_ip ) then
             ! Local to global
             call sld_csys_rotuni(2_ip,ndime,ipoin,fcont_sld(:,ipoin)) ! Not exchanged yet!
          end if
       end do
    end if
    !
    ! Get total force (sum) to send to rigid
    !
    force_dbo = 0.0_rp
    if( INOTMASTER ) then
       force_dbo(1:ndime) = sum(fcont_sld(1:ndime,:), dim=2)
       !
       ! Parall exchange
       !
       call rhsmod(ndime,fcont_sld)
    end if
    !
    ! Parall service
    !
    call PAR_SUM(ndime,force_dbo(1:ndime),'IN MY CODE')     
    !
    ! Deallocate release node vector
    !
    if( associated(kfl_rnode_sld) ) then
       call memory_deallo(mem_modul(1:2,modul),'KFL_RNODE_SLD','sld_pdn_contact_force_send',kfl_rnode_sld)
    end if

  end subroutine sld_pdn_contact_force_send

  !-----------------------------------------------------------------------
  !>
  !> @author  gguillamet
  !> @date    2021-07-21
  !> @brief   PDN-contact projection
  !> @details PDN-contact projection
  !>
  !-----------------------------------------------------------------------

  subroutine sld_pdn_contact_projection(nsend,local_basis,norma_dista)

    use def_master,         only : displ,IPARALL,ITER_K
    use def_elmtyp,         only : BAR02,TRI03,QUA04
    use def_domain,         only : coord,lnodb,kfl_codbo
    use def_domain,         only : nboun,ndime,mnodb,nnode
    use def_domain,         only : ltypb
    use mod_parall,         only : PAR_COMM_MY_CODE
    use mod_communications, only : PAR_ALLGATHER
    use mod_communications, only : PAR_ALLGATHERV
    use mod_communications, only : PAR_SUM,PAR_MAX
    use mod_communications, only : PAR_COMM_RANK_AND_SIZE
    use mod_maths,          only : maths_normalize_vector
    use mod_maths_geometry, only : maths_point_segment_distance
    use mod_maths_geometry, only : maths_point_in_triangle
    use mod_maths_geometry, only : maths_point_face_distance
    use mod_maths_geometry, only : maths_normal_to_triangle
    use mod_maths,          only : maths_vectorial_product
    use mod_elmgeo,         only : elmgeo_projection_on_a_face
    use mod_elmgeo,         only : elmgeo_natural_coordinates_on_boundaries

    implicit none

    integer(ip),       intent(in)    :: nsend
    real(rp), pointer, intent(inout) :: local_basis(:,:,:)
    real(rp), pointer, intent(inout) :: norma_dista(:)

    integer(ip), pointer :: recvcount1(:)
    integer(ip), pointer :: recvcount2(:)
    integer(ip), pointer :: pblty_glo(:)
    integer(ip), pointer :: pblty_loc(:)
    integer(ip), pointer :: pnodb_glo(:)
    integer(ip), pointer :: pnodb_loc(:)
    real(rp),    pointer :: bocod_glo(:,:,:)
    real(rp),    pointer :: bocod_loc(:,:,:)
    real(rp),    pointer :: coord_j(:,:)

    integer(ip)          :: iboun,kboun
    integer(ip)          :: inodb,ifoun,idime,ifoun2
    integer(ip)          :: ipoin,kpoin
    integer(ip)          :: pnodb,pblty
    integer(ip)          :: nboun_loc_i,nboun_glo_i,ntria,mnodb_max
    integer(ip)          :: rankempi,sizempi
    integer(ip)          :: sizearray1,sizearray2

    real(rp)             :: bari1,bari2
    real(rp)             :: N(ndime),T1(ndime),T2(ndime),rho
    real(rp)             :: dista,proje(ndime)

    !--------------------------------------------------------------------
    !
    ! PLEPP stuff
    !
    !--------------------------------------------------------------------

    nullify(coord_j)

    if( INOTMASTER .and. nsend > 0 ) then
       !
       ! Coordinates of ALL nodes of the other instance which have penetrated my current instance
       !
       call memory_alloca(mem_modul(1:2,modul),'COORD_J','sld_pdn_contact_projection',coord_j,ndime,nsend)
       coord_j = 0.0_rp
       do kpoin = 1,nsend
          do idime = 1,ndime
             coord_j(idime,kpoin) = PLEPP_CPLNG%dist_coords_j(ndime*(kpoin-1)+idime)
          end do
       end do
    end if

    !--------------------------------------------------------------------
    !
    ! Total contact boundary list (rigid body)
    !
    !--------------------------------------------------------------------

    nullify(recvcount1)
    nullify(recvcount2)
    nullify(pblty_loc)
    nullify(pblty_glo)
    nullify(pnodb_loc)
    nullify(pnodb_glo)
    nullify(bocod_loc)
    nullify(bocod_glo)
    !
    ! Calculate total boundaries body i (rigid)
    !
    nboun_loc_i = 0_ip
    mnodb_max = 0_ip
    if( INOTMASTER ) then
       nboun_loc_i = count(kfl_codbo(1:nboun) == bcont_sld)
       mnodb_max = mnodb
    end if
    !
    ! Parall
    !
    nboun_glo_i = nboun_loc_i
    call PAR_SUM(nboun_glo_i,'IN MY CODE')
    call PAR_MAX(mnodb_max,  'IN MY CODE')
    !
    ! Get total boundaries and their nodal coordinates from body i
    !
    call memory_alloca(mem_modul(1:2,modul),'PBLTY_GLO','sld_pdn_contact_projection',pblty_glo,nboun_glo_i)
    call memory_alloca(mem_modul(1:2,modul),'PNODB_GLO','sld_pdn_contact_projection',pnodb_glo,nboun_glo_i)
    call memory_alloca(mem_modul(1:2,modul),'BOCOD_GLO','sld_pdn_contact_projection',bocod_glo,ndime,mnodb_max,nboun_glo_i)
    pblty_glo = 0_ip
    pnodb_glo = 0_ip
    bocod_glo = 0.0_rp
    if( INOTMASTER ) then
       call memory_alloca(mem_modul(1:2,modul),'PBLTY_LOC','sld_pdn_contact_projection',pblty_loc,nboun_loc_i)
       call memory_alloca(mem_modul(1:2,modul),'PNODB_LOC','sld_pdn_contact_projection',pnodb_loc,nboun_loc_i)
       call memory_alloca(mem_modul(1:2,modul),'BOCOD_LOC','sld_pdn_contact_projection',bocod_loc,ndime,mnodb_max,nboun_loc_i)
       kboun = 0
       do iboun = 1,nboun
          pblty = ltypb(iboun)
          pnodb = nnode(pblty)
          if( kfl_codbo(iboun) == bcont_sld ) then
             kboun = kboun + 1_ip
             pblty_loc(kboun) = pblty
             pnodb_loc(kboun) = pnodb
             do inodb = 1,pnodb
                ipoin = lnodb(inodb,iboun)
                bocod_loc(1:ndime,inodb,kboun) = coord(1:ndime,ipoin) + displ(1:ndime,ipoin,ITER_K)
             end do
          end if
       end do
    else
       call memory_alloca(mem_modul(1:2,modul),'PBLTY_LOC','sld_pdn_contact_projection',pblty_loc,0_ip)
       call memory_alloca(mem_modul(1:2,modul),'PNODB_LOC','sld_pdn_contact_projection',pnodb_loc,0_ip)
       call memory_alloca(mem_modul(1:2,modul),'BOCOD_LOC','sld_pdn_contact_projection',bocod_loc,ndime,mnodb_max,0_ip)
    end if
    !
    ! Parall
    !
    sizearray1 = nboun_loc_i
    sizearray2 = ndime*mnodb_max*nboun_loc_i
    call PAR_COMM_RANK_AND_SIZE(PAR_COMM_MY_CODE,rankempi,sizempi)
    call memory_alloca(mem_modul(1:2,modul),'RECVCOUNT1','sld_pdn_contact_projection',recvcount1,sizempi)
    call memory_alloca(mem_modul(1:2,modul),'RECVCOUNT2','sld_pdn_contact_projection',recvcount2,sizempi)
    recvcount1 = -666_ip
    recvcount2 = -666_ip
    if( IPARALL ) then
       call PAR_ALLGATHER(sizearray1,recvcount1,1_4,       'IN MY CODE')
       call PAR_ALLGATHER(sizearray2,recvcount2,1_4,       'IN MY CODE')
       call PAR_ALLGATHERV(pblty_loc,pblty_glo,recvcount1, 'IN MY CODE')
       call PAR_ALLGATHERV(pnodb_loc,pnodb_glo,recvcount1, 'IN MY CODE')
       call PAR_ALLGATHERV(bocod_loc,bocod_glo,recvcount2, 'IN MY CODE')
    else
       recvcount1 = sizearray1
       recvcount2 = sizearray2
       pblty_glo  = pblty_loc
       pnodb_glo  = pnodb_loc
       bocod_glo  = bocod_loc
    end if

    !--------------------------------------------------------------------
    !
    ! Projections
    !
    !--------------------------------------------------------------------

    if( INOTMASTER ) then

       if( ndime == 2 ) then

          do kpoin = 1,nsend
             norma_dista(kpoin)     = 0.0_rp
             local_basis(:,:,kpoin) = 0.0_rp
             iboun_loop2d: do kboun = 1,nboun_glo_i
                pblty = pblty_glo(kboun) 
                pnodb = pnodb_glo(kboun)
                ifoun = 0_ip
                if( pblty /= BAR02 ) then
                   call runend("SLD_PDN_CONTACT_PROJECTION: projections other than BAR02 not implemented")
                end if
                call maths_point_segment_distance(ndime,coord_j(:,kpoin),bocod_glo(:,1,kboun),bocod_glo(:,2,kboun),dista,proje,ifoun)
                if( ifoun == 1_ip ) then ! Projected point is inside the boundary element
                   !
                   ! Node projection on the element segment
                   !
                   ! Normal
                   N(1:ndime) = proje(1:ndime)-coord_j(1:ndime,kpoin)
                   call maths_normalize_vector(ndime,N(:),rho)
                   local_basis(1:ndime,1,kpoin) = N(1:ndime)
                   ! Tangent
                   local_basis(1,2,kpoin) = -N(2)
                   local_basis(2,2,kpoin) =  N(1)
                   ! Distance
                   norma_dista(kpoin) = dista
                   exit iboun_loop2d
                end if
             end do iboun_loop2d
          end do
          
       else if( ndime == 3 ) then

          do kpoin = 1,nsend
             iboun_loop3d: do kboun = 1,nboun_glo_i
                pblty = pblty_glo(kboun) 
                pnodb = pnodb_glo(kboun)  
                ifoun = 0_ip
                if( pblty == TRI03 .or. pblty == QUA04 ) then
                   call maths_point_in_triangle(pnodb,coord_j(:,kpoin),bocod_glo(:,:,kboun),ifoun,bari1,bari2,ntria)
                   call instri(pnodb,coord_j(:,kpoin),bocod_glo(:,:,kboun),ifoun2)
                   if( ifoun /= ifoun2) call runend('PROJECTION IS DIFERENT BETWEEN SUBRUS')
                else
                   call runend("SLD_PDN_CONTACT_PROJECTION: projections other than TRI03 or QUA04 not implemented")
                end if
                if( ifoun == 1_ip ) then ! Projected point is inside the boundary element
                   !
                   ! Node projection on the element face
                   !
                   call elmgeo_projection_on_a_face(ndime,pblty,bocod_glo(:,:,kboun),coord_j(:,kpoin),proje)
                   dista = sqrt(dot_product(proje(1:ndime)-coord_j(1:ndime,kpoin),proje(1:ndime)-coord_j(1:ndime,kpoin)))
                   ! Normal
                   N(1:ndime) = proje(1:ndime)-coord_j(1:ndime,kpoin)
                   call maths_normalize_vector(ndime,N(:),rho)
                   local_basis(1:ndime,1,kpoin) = N(1:ndime)
                   ! Tangents
                   T1(:) = bocod_glo(1:ndime,2,kboun) - bocod_glo(1:ndime,1,kboun)
                   call maths_normalize_vector(ndime,T1(:),rho)
                   local_basis(1:ndime,2,kpoin) = T1(1:ndime)
                   call maths_vectorial_product(N(:), T1(:), T2(:), ndime)
                   local_basis(1:ndime,3,kpoin) = T2(1:ndime)
                   ! Dista
                   norma_dista(kpoin) = dista
                   exit iboun_loop3d
                end if
             end do iboun_loop3d
          end do

       end if

    end if

    !--------------------------------------------------------------------
    !
    ! Deallocate memory
    !
    !--------------------------------------------------------------------

    call memory_deallo(mem_modul(1:2,modul),'RECVCOUNT1','sld_pdn_contact_projection',recvcount1)
    call memory_deallo(mem_modul(1:2,modul),'RECVCOUNT2','sld_pdn_contact_projection',recvcount2)
    call memory_deallo(mem_modul(1:2,modul),'PBLTY_GLO', 'sld_pdn_contact_projection',pblty_glo )
    call memory_deallo(mem_modul(1:2,modul),'PNODB_GLO', 'sld_pdn_contact_projection',pnodb_glo )
    call memory_deallo(mem_modul(1:2,modul),'BOCOD_GLO', 'sld_pdn_contact_projection',bocod_glo )
    if( INOTMASTER ) then
       call memory_deallo(mem_modul(1:2,modul),'COORD_J',  'sld_pdn_contact_projection',coord_j  )
       call memory_deallo(mem_modul(1:2,modul),'PBLTY_LOC','sld_pdn_contact_projection',pblty_loc)
       call memory_deallo(mem_modul(1:2,modul),'PNODB_LOC','sld_pdn_contact_projection',pnodb_loc)
       call memory_deallo(mem_modul(1:2,modul),'BOCOD_LOC','sld_pdn_contact_projection',bocod_loc)
    end if

  end subroutine sld_pdn_contact_projection

  !-----------------------------------------------------------------------
  !>
  !> @author  gguillamet
  !> @date    2021-07-21
  !> @brief   PDN-contact release nodes
  !> @details PDN-contact release nodes
  !>
  !-----------------------------------------------------------------------

  subroutine sld_pdn_contact_releas()

    use def_master,  only : displ
    use def_master,  only : ITER_K,TIME_N
    use def_solidz,  only : veloc_sld,accel_sld
    
    implicit none

    integer(ip) :: ipoin, kpoin, itott
    integer(ip) :: nrecv

    if( CNT_CPLNG%current_code == CNT_CPLNG%code_j ) then
       !
       ! Check existance of sticky nodes
       !
       if( kfl_goite_sld /= 1_ip ) then     ! if converged
          !
          ! Compute contact reaction and mark nodes in traction
          !
          if( INOTMASTER ) then
             nrecv = PLEPP_CPLNG%n_ji
             !
             ! Get contact force and mark artificial nodes
             !
             do kpoin = 1,nrecv
                ipoin = PLEPP_CPLNG%interior_list_j(kpoin)
                itott = (ipoin - 1) * ndime
                !
                ! Nodes in contact
                !
                if( kfl_fixno_sld(1,ipoin) == 3_ip ) then
                   if( frxid_sld(itott+1) > 0.0_rp .and. &
                       kfl_rnode_sld(kpoin) == 0_ip ) then ! not marked
                      kfl_rnode_sld(kpoin) = 1_ip
                      kfl_relea_sld        = .true.
                   end if
                end if
             end do
             if( associated(kfl_rnode_sld) ) narco_sld = sum(kfl_rnode_sld)
          end if
          !
          ! Release nodes (if any) and repeat iteration
          !
          if( INOTMASTER .and. kfl_relea_sld ) then
             nrecv = PLEPP_CPLNG%n_ji
             do kpoin = 1,nrecv
                ipoin = PLEPP_CPLNG%interior_list_j(kpoin)
                if( kfl_fixno_sld(1,ipoin) == 3_ip .and. &
                    kfl_rnode_sld(kpoin)   == 1_ip ) then
                   ! Set previous state and set node free
                   kfl_fixno_sld(1,ipoin) = 0_ip
                   kfl_rnode_sld(kpoin)   = 0_ip
                   bvess_sld(1:ndime,ipoin,ITER_K) = 0.0_rp
                   jacrot_du_dq_sld(:,:,ipoin)     = 0.0_rp
                   jacrot_dq_du_sld(:,:,ipoin)     = 0.0_rp
                   displ(1:ndime,ipoin,ITER_K)     = displ(1:ndime,ipoin,TIME_N)
                   veloc_sld(1:ndime,ipoin,ITER_K) = veloc_sld(1:ndime,ipoin,TIME_N)
                   accel_sld(1:ndime,ipoin,ITER_K) = accel_sld(1:ndime,ipoin,TIME_N) 
                   nacco_sld                       = nacco_sld - 1_ip
                end if
             end do
             ! Init flag
             kfl_relea_sld = .false.
             ! Repeat iteration (not converged)
             kfl_goite_sld = 1_ip
             itrel_sld = itrel_sld + 1_ip
          end if

       end if
       !
       ! Check if all subdomains are converged to solution
       !
       call PAR_SUM(kfl_goite_sld, 'IN MY CODE')
       if( kfl_goite_sld >= 1_ip ) kfl_goite_sld = 1_ip
       !
       ! Parall of some convergence parameters
       !
       call PAR_SUM(nacco_sld, 'IN MY CODE')
       call PAR_SUM(ntoco_sld, 'IN MY CODE')
       call PAR_SUM(narco_sld, 'IN MY CODE')
       call PAR_MAX(itrel_sld, 'IN MY CODE')

    end if

  end subroutine sld_pdn_contact_releas

  !-----------------------------------------------------------------------
  !>
  !> @author  gguillamet
  !> @date    2021-07-21
  !> @brief   PDN-contact read data
  !> @details PDN-contact read data
  !>
  !-----------------------------------------------------------------------

  subroutine sld_pdn_contact_reabcs()

    use def_inpout,   only : words, exists, getint, getrea, getcha
    use mod_ecoute,   only : ecoute
    use def_solidz,   only : kfl_conta_sld

    implicit none

    call ecoute('sld_reabcs')
    do while( words(1) /= 'ENDCO' )
       if( words(1) == 'RBODE' ) then
          kfl_pdnco_sld = SLD_PDN_RIGID_TO_DEFOR
          kfl_conta_sld = SLD_PDN_RIGID_TO_DEFOR ! Old one required for testing

       else if( words(1) == 'BOUND' ) then
          bcont_sld = getint('BOUND',1_ip,'#Contact boundary')

       else if( words(1) == 'PENET' ) then
          tolpe_sld = getrea('PENET',1.0e-4_rp,'#Penetration tolerance')

       else if( words(1) == 'CONVE' ) then
          if( words(2) =='ON   ') then
             kfl_pdncf_sld = .true.
          else if( words(2) == 'OFF  ' ) then
             kfl_pdncf_sld = .false.
          else
             kfl_pdncf_sld = .false.
          end if
       end if
       call ecoute('sld_reabcs')
    end do

  end subroutine sld_pdn_contact_reabcs

  !-----------------------------------------------------------------------
  !>
  !> @author  gguillamet
  !> @date    2021-07-21
  !> @brief   PDN-contact parall send data
  !> @details PDN-contact parall send dat
  !>
  !-----------------------------------------------------------------------

  subroutine sld_pdn_contact_parall()

    use mod_exchange, only : exchange_add
    use mod_exchange, only : exchange_end

    call exchange_add(kfl_pdnco_sld)
    call exchange_add(kfl_pdncf_sld)
    call exchange_add(bcont_sld)
    call exchange_add(tolpe_sld)
    call exchange_end()

  end subroutine sld_pdn_contact_parall

  !-----------------------------------------------------------------------
  !>
  !> @author  gguillamet
  !> @date    2021-07-21
  !> @brief   PDN-contact convergence
  !> @details PDN-contact convergence
  !>
  !-----------------------------------------------------------------------

  subroutine sld_pdn_contact_cvgunk()

    use def_master, only : INOTSLAVE
    use def_master, only : modul
    use def_master, only : solve_sol, kfl_rstar
    use def_master, only : itinn, itcou, ittim
    use mod_iofile, only : iofile_flush_unit

    implicit none

    integer(ip), save   :: kpass = 0

    if( INOTSLAVE ) then
       !
       ! Write heading
       !
       if( kpass == 0 .and. kfl_rstar /= 2 ) then
          write(lun_pdnco_sld,400)
       end if
       !
       ! Write convergence
       !
       write(lun_pdnco_sld,401) &
            !   1     2            3                  4
            ittim,itcou,itinn(modul),solve_sol(1)%iters, &
            !       5         6         7
            itrel_sld,ntoco_sld,narco_sld,nacco_sld

       call iofile_flush_unit(lun_pdnco_sld)
    end if

    kpass = 1_ip

400 format('# --| ALYA PDN-contact convergence'  ,/,&
         & '# --| Columns displayed:' ,/,&
         & '# --|                                                                            ',/,&
         & '# --|  1. Time step           2. Global Iteration          3. Inner Iteration    ',/,&
         & '# --|  4. Solver iterations   5. Iteration release nodes   6. Total contact nodes',/,&
         & '# --|  7. Art. contact nodes  8. Act. contact nodes                              ',/,&
         & '# --|                                                                            ')

401 format(4x,i9,2x,i9,2x,i9,2x,i9,2x,i9,2x,i9,2x,i9,2x,i9)

  end subroutine sld_pdn_contact_cvgunk
  
  !-----------------------------------------------------------------------
  !>
  !> @author  gguillamet
  !> @date    2021-07-21
  !> @brief   PDN-contact begste
  !> @details PDN-contact begste
  !>
  !-----------------------------------------------------------------------
  
  subroutine sld_pdn_contact_begste
    
    use def_master, only : iblok
    use def_coupli, only : coupling_driver_max_iteration
    use def_solidz, only : coupling_contact_its, kfl_conta_sld
    use def_solidz, only : SLD_PDN_BILATERAL
    
    implicit none
    !
    ! Initializations
    !
    nacco_sld = 0
    narco_sld = 0
    ntoco_sld = 0
    itrel_sld = 0
    
    if( kfl_conta_sld == SLD_PDN_BILATERAL ) then
       !
       ! Bilateral (Matias)
       !
       coupling_driver_max_iteration(iblok) = coupling_contact_its
    else
       !
       ! Unilateral (Matias) and Rigid-Deformable (Gerard)
       !
       coupling_driver_max_iteration(iblok) = 1_ip
    end if
    
  end subroutine sld_pdn_contact_begste
  
  !-----------------------------------------------------------------------
  !>
  !> @author  aquintan and gguillam
  !> @date    2021-07-21
  !> @brief   PDN-contact get rank master for deformable body
  !> @details PDN-contact get rank master for deformable body
  !>
  !-----------------------------------------------------------------------
  
  subroutine sld_pdn_contact_master_rank()
    
    use def_master,         only : IMASTER
    use mod_parall,         only : PAR_MY_UNIVERSE_RANK
    use mod_communications, only : PAR_SUM

    implicit none

    integer(ip), pointer         :: ranks(:)

    allocate(ranks(2))
    ranks(:) = 0_ip
    if(      CNT_CPLNG%current_code == CNT_CPLNG%code_i ) then
       if( IMASTER ) then
          ranks(1) = PAR_MY_UNIVERSE_RANK
       end if
    else if( CNT_CPLNG%current_code == CNT_CPLNG%code_j ) then
       if( IMASTER ) then
          ranks(2) = PAR_MY_UNIVERSE_RANK
       endif
    end if
    call PAR_SUM(ranks,'IN THE UNIVERSE')
    MPI_RANK_MASTER_RBO = ranks(1)
    MPI_RANK_MASTER_DEF = ranks(2)
    deallocate(ranks)

  end subroutine sld_pdn_contact_master_rank

#endif
  
end module mod_sld_pdn_contact_plepp

