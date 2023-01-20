!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



module mod_meshes_gather_submesh

  use def_kintyp_basic,   only : lg
  use def_domain
  use mod_parall,         only : PAR_WORLD_SIZE
  use def_master,         only : IMASTER
  use def_master,         only : INOTMASTER
  use def_master,         only : lninv_loc
  use def_master,         only : title
  use def_master,         only : INOTSLAVE
  use def_master,         only : NPOIN_TYPE
  use def_master,         only : NELEM_TYPE
  use def_master,         only : NBOUN_TYPE
  use mod_communications, only : PAR_GATHERV
  use mod_communications, only : PAR_GATHER
  use mod_communications, only : PAR_MAX
  use mod_output,         only : output_mesh_gid_format 
  use mod_output,         only : output_result_gid_format 
  use mod_output,         only : output_domain_alya_format
  use mod_iofile,         only : iofile_open_unit
  use mod_mesh_type,      only : mesh_type_initialize
  use mod_messages,       only : messages_live
  use mod_memory,         only : memory_alloca
  use mod_memory,         only : memory_deallo
  use mod_iofile,         only : iofile_open_unit
  use mod_iofile,         only : iofile_close_unit
  
  implicit none
  private

  public :: meshes_gather_submesh
  
contains

  subroutine meshes_gather_submesh(xx_in,diameter,NODE_NUMBER)

    real(rp),       intent(in)           :: xx_in(3)
    real(rp),       intent(in)           :: diameter
    integer(ip),    intent(in), optional :: NODE_NUMBER
    integer(ip)                          :: ipoin,idime,knode,ielem,kpoin,kboun
    integer(ip)                          :: inode,ineig,iboun,kelem,nneig_gat
    integer(ip)                          :: jpoin,ndim1,ndim2,ifiel,nfiel_sav
    integer(ip)                          :: nelem_sub,npoin_sub,nboun_sub
    integer(ip)                          :: kneig
    real(rp)                             :: dummr,xx(3)
    type(r3p)                            :: xfiel_sub(mfiel)

    integer(ip),    pointer              :: mark_node(:)
    logical(lg),    pointer              :: mark_elem(:)
    logical(lg),    pointer              :: mark_boun(:)
    integer(ip),    pointer              :: permu_node(:)
    integer(ip),    pointer              :: permu_elem(:)
    integer(ip),    pointer              :: permu_boun(:)
    integer(ip),    pointer              :: invpe_node(:)
    integer(ip),    pointer              :: invpe_elem(:)

    integer(ip),    pointer              :: nelem_gat(:)
    integer(ip),    pointer              :: npoin_gat(:)
    integer(ip),    pointer              :: nboun_gat(:)
    integer(4),     pointer              :: lsiz4_gat(:)

    integer(ip),    pointer              :: lnods_sub(:,:)
    integer(ip),    pointer              :: ltype_sub(:)
    integer(ip),    pointer              :: lesub_sub(:)
    integer(ip),    pointer              :: lmate_sub(:)
    integer(ip),    pointer              :: lnnod_sub(:)
    integer(ip),    pointer              :: leset_sub(:)

    integer(ip),    pointer              :: lnodb_sub(:,:)
    integer(ip),    pointer              :: ltypb_sub(:)
    integer(ip),    pointer              :: lnnob_sub(:)
    integer(ip),    pointer              :: lelbo_sub(:)
    integer(ip),    pointer              :: codbo_sub(:)
    integer(ip),    pointer              :: lbset_sub(:)

    real(rp),       pointer              :: coord_sub(:,:)
    integer(ip),    pointer              :: lninv_sub(:)

    real(rp),       pointer              :: xarray(:)

    integer(ip)                          :: nunit_msh
    integer(ip)                          :: nunit_res
    integer(ip)                          :: nunit_dat
    type(mesh_type)                      :: meshe_sub

    call mesh_type_initialize(meshe_sub)
    nfiel_sav = nfiel

    nullify(mark_node)
    nullify(mark_elem)
    nullify(mark_boun)
    nullify(permu_node)
    nullify(permu_elem)
    nullify(permu_boun)
    nullify(invpe_node)
    nullify(invpe_elem)

    nullify(nelem_gat)
    nullify(npoin_gat)
    nullify(nboun_gat)
    nullify(lsiz4_gat)

    nullify(lnods_sub)
    nullify(ltype_sub)
    nullify(lesub_sub)
    nullify(lmate_sub)
    nullify(lnnod_sub)
    nullify(leset_sub)

    nullify(lnodb_sub)
    nullify(ltypb_sub)
    nullify(lnnob_sub)
    nullify(lelbo_sub)
    nullify(codbo_sub)
    nullify(lbset_sub)

    nullify(coord_sub)
    nullify(lninv_sub)

    nullify(xarray)

    npoin_sub = 0
    nelem_sub = 0
    nboun_sub = 0
    allocate(mark_node(npoin))
    allocate(mark_elem(nelem))
    allocate(mark_boun(nboun))

    if( present(NODE_NUMBER) ) then
       ipoin_loop: do ipoin = 1,npoin_own
          if( lninv_loc(ipoin) == NODE_NUMBER ) then
             xx(1:ndime) = coord(1:ndime,ipoin)
             exit ipoin_loop
          end if
       end do ipoin_loop
       call PAR_MAX(ndime,xx)
    else
       xx(1:ndime) = xx_in(1:ndime)
    end if

    !--------------------------------------------------------------------
    !
    ! Create submeshes from marked nodes
    !
    !--------------------------------------------------------------------

    if( INOTMASTER ) then
       !
       ! Mark nodes
       !
       do ipoin = 1,npoin
          dummr = 0.0_rp
          do idime = 1,ndime
             dummr = dummr + (coord(idime,ipoin)-xx(idime))*(coord(idime,ipoin)-xx(idime))
          end do
          if( sqrt(dummr) <= diameter ) then
             mark_node(ipoin) = 1
             npoin_sub     = 1
          else
             mark_node(ipoin) = 0
          end if
       end do
       if( npoin_sub == 1 ) then
          !
          ! Mark elements
          !
          nelem_sub = 0
          do ielem = 1,nelem
             mark_elem(ielem) = .false.
             loop_inode: do inode = 1,lnnod(ielem)
                ipoin = lnods(inode,ielem)
                if( mark_node(ipoin) == 1 ) then
                   do knode = 1,lnnod(ielem)
                      kpoin = lnods(knode,ielem)
                      if( mark_node(kpoin) == 0 ) mark_node(kpoin) = 2
                   end do
                   nelem_sub = nelem_sub + 1
                   mark_elem(ielem) = .true.
                   exit loop_inode
                end if
             end do loop_inode
          end do
          !
          ! Mark boundaries
          !
          nboun_sub = 0
          do iboun = 1,nboun
             ielem = lelbo(iboun)
             mark_boun(iboun) = mark_elem(ielem)
             if( mark_boun(iboun) ) nboun_sub = nboun_sub + 1
          end do
          !
          ! Renumber nodes
          !
          npoin_sub = 0
          do ipoin = 1,npoin
             if( mark_node(ipoin) /= 0 ) then
                npoin_sub = npoin_sub + 1
                mark_node(ipoin) = npoin_sub
             end if
          end do
          !
          ! Permutation
          !
          allocate(permu_node(npoin_sub))
          allocate(permu_elem(nelem_sub))
          allocate(permu_boun(nelem_sub))
          allocate(invpe_node(npoin))
          allocate(invpe_elem(nelem))
          npoin_sub = 0
          do ipoin = 1,npoin
             if( mark_node(ipoin) /= 0 ) then
                npoin_sub             = npoin_sub + 1
                permu_node(npoin_sub) = ipoin
                invpe_node(ipoin)     = npoin_sub
             end if
          end do
          nelem_sub = 0
          do ielem = 1,nelem
             if( mark_elem(ielem) ) then
                nelem_sub             = nelem_sub + 1
                permu_elem(nelem_sub) = ielem
                invpe_elem(ielem)     = nelem_sub
             end if
          end do
          nboun_sub = 0
          do iboun = 1,nboun
             if( mark_boun(iboun) ) then
                nboun_sub             = nboun_sub + 1
                permu_boun(nboun_sub) = iboun
             end if
          end do
          !
          ! Create submesh
          !
          call memory_alloca(memor_dom,'LNODS','mod_meshes',lnods_sub,mnode,nelem_sub)
          call memory_alloca(memor_dom,'LTYPE','mod_meshes',ltype_sub,nelem_sub)
          call memory_alloca(memor_dom,'LESUB','mod_meshes',lesub_sub,nelem_sub)
          call memory_alloca(memor_dom,'LMATE','mod_meshes',lmate_sub,nelem_sub)
          call memory_alloca(memor_dom,'LNNOD','mod_meshes',lnnod_sub,nelem_sub)
          if( neset > 0 ) call memory_alloca(memor_dom,'LESET','mod_meshes',leset_sub,nelem_sub)

          call memory_alloca(memor_dom,'LNODS','mod_meshes',lnodb_sub,mnodb,nboun_sub)
          call memory_alloca(memor_dom,'LTYPE','mod_meshes',ltypb_sub,nboun_sub)
          call memory_alloca(memor_dom,'LNNOB','mod_meshes',lnnob_sub,nboun_sub)
          call memory_alloca(memor_dom,'LELBO','mod_meshes',lelbo_sub,nboun_sub)
          call memory_alloca(memor_dom,'CODBO','mod_meshes',codbo_sub,nboun_sub)
          if( nbset > 0 ) call memory_alloca(memor_dom,'LBSET','mod_meshes',lbset_sub,nboun_sub)

          call memory_alloca(memor_dom,'COORD','mod_meshes',coord_sub,ndime,npoin_sub)
          call memory_alloca(memor_dom,'LNINV','mod_meshes',lninv_sub,npoin_sub)

          do ifiel = 1,nfiel
             ndim1 = kfl_field(1,ifiel)
             ndim2 = kfl_field(4,ifiel)
             if(      kfl_field(2,ifiel) == NELEM_TYPE ) then
                call memory_alloca(memor_dom,'XFIEL','mod_meshes',xfiel_sub(ifiel) % a,ndim1,nelem_sub,ndim2)
             else if( kfl_field(2,ifiel) == NPOIN_TYPE ) then                                                                        
                call memory_alloca(memor_dom,'XFIEL','mod_meshes',xfiel_sub(ifiel) % a,ndim1,npoin_sub,ndim2)
             else if( kfl_field(2,ifiel) == NBOUN_TYPE ) then                                                                        
                call memory_alloca(memor_dom,'XFIEL','mod_meshes',xfiel_sub(ifiel) % a,ndim1,nboun_sub,ndim2)
             end if
          end do

          do ielem = 1,nelem_sub
             kelem              = permu_elem(ielem)
             lnods_sub(:,ielem) = invpe_node(lnods(:,kelem))
             ltype_sub(ielem)   = ltype(kelem)
             lesub_sub(ielem)   = lesub(kelem)
             lmate_sub(ielem)   = lmate(kelem)
             lnnod_sub(ielem)   = lnnod(kelem)
             if( neset > 0 ) leset_sub(ielem)   = leset(kelem)
          end do
          do iboun = 1,nboun_sub
             kboun              = permu_boun(iboun)
             lnodb_sub(:,iboun) = invpe_node(lnodb(:,kboun))
             ltypb_sub(iboun)   = ltypb(kboun)
             lnnob_sub(iboun)   = lnnob(kboun)
             lelbo_sub(iboun)   = invpe_elem(lelbo(kboun))
             codbo_sub(iboun)   = kfl_codbo(kboun)
             if( nbset > 0 ) lbset_sub(iboun)   = lbset(kboun)
          end do
          do ipoin = 1,npoin_sub
             kpoin              = permu_node(ipoin)
             coord_sub(:,ipoin) = coord(:,kpoin)
             lninv_sub(ipoin)   = lninv_loc(kpoin)
          end do
          do ifiel = 1,nfiel
             if(      kfl_field(2,ifiel) == NELEM_TYPE ) then
                do ielem = 1,nelem_sub
                   kelem = permu_elem(ielem)
                   xfiel_sub(ifiel) % a(:,ielem,:) = xfiel(ifiel) % a(:,kelem,:) 
                end do
             else if( kfl_field(2,ifiel) == NPOIN_TYPE ) then                                                                        
                do ipoin = 1,npoin_sub
                   kpoin = permu_node(ipoin)
                   xfiel_sub(ifiel) % a(:,ipoin,:) = xfiel(ifiel) % a(:,kpoin,:) 
                end do
             else if( kfl_field(2,ifiel) == NBOUN_TYPE ) then                                                                        
                do iboun = 1,nboun_sub
                   kboun = permu_boun(iboun)
                   xfiel_sub(ifiel) % a(:,iboun,:) = xfiel(ifiel) % a(:,kboun,:) 
                end do
             end if
          end do
          if( kfl_ngrou /= 0 ) then
             nfiel              = nfiel + 1
             ifiel              = nfiel
             kfl_field(1,ifiel) = 1
             kfl_field(2,ifiel) = NPOIN_TYPE
             kfl_field(4,ifiel) = 1
             ndim1              = kfl_field(1,ifiel)
             ndim2              = kfl_field(4,ifiel)
             call memory_alloca(memor_dom,'XFIEL','mod_meshes',xfiel_sub(ifiel) % a,ndim1,npoin_sub,ndim2)
             do ipoin = 1,npoin_sub
                kpoin = permu_node(ipoin)
                xfiel_sub(ifiel) % a(1,ipoin,1) = real(lgrou_dom(kpoin),rp) 
             end do
          end if
       end if

    end if

    !--------------------------------------------------------------------
    !
    ! Gather geometry
    !
    !--------------------------------------------------------------------

    if( INOTSLAVE ) then
       nneig_gat = PAR_WORLD_SIZE-1
       allocate( nelem_gat(0:PAR_WORLD_SIZE-1) )
       allocate( nboun_gat(0:PAR_WORLD_SIZE-1) )
       allocate( npoin_gat(0:PAR_WORLD_SIZE-1) )
       allocate( lsiz4_gat(0:PAR_WORLD_SIZE-1) )
    else
       nneig_gat = -1
    end if

    call PAR_GATHER(nelem_sub,nelem_gat,'IN THE WORLD') 
    call PAR_GATHER(nboun_sub,nboun_gat,'IN THE WORLD') 
    call PAR_GATHER(npoin_sub,npoin_gat,'IN THE WORLD')

    if( INOTSLAVE ) then
       nelem_sub = sum(nelem_gat)
       nboun_sub = sum(nboun_gat)
       npoin_sub = sum(npoin_gat)
       call memory_alloca(memor_dom,'LNODS','mod_meshes',lnods_sub,mnode,nelem_sub)
       call memory_alloca(memor_dom,'LTYPE','mod_meshes',ltype_sub,nelem_sub)
       call memory_alloca(memor_dom,'LESUB','mod_meshes',lesub_sub,nelem_sub)
       call memory_alloca(memor_dom,'LMATE','mod_meshes',lmate_sub,nelem_sub)
       call memory_alloca(memor_dom,'LNNOD','mod_meshes',lnnod_sub,nelem_sub)
       if( neset > 0 ) call memory_alloca(memor_dom,'LESET','mod_meshes',leset_sub,nelem_sub)

       call memory_alloca(memor_dom,'LNODB','mod_meshes',lnodb_sub,mnodb,nboun_sub)
       call memory_alloca(memor_dom,'LTYPB','mod_meshes',ltypb_sub,nboun_sub)
       call memory_alloca(memor_dom,'LNNOB','mod_meshes',lnnob_sub,nboun_sub)
       call memory_alloca(memor_dom,'LELBO','mod_meshes',lelbo_sub,nboun_sub)
       call memory_alloca(memor_dom,'CODBO','mod_meshes',codbo_sub,nboun_sub)
       if( nbset > 0 ) call memory_alloca(memor_dom,'LBSET','mod_meshes',lbset_sub,nboun_sub)

       call memory_alloca(memor_dom,'COORD','mod_meshes',coord_sub,ndime,npoin_sub)
       call memory_alloca(memor_dom,'LNINV','mod_meshes',lninv_sub,npoin_sub)

       do ifiel = 1,nfiel
          ndim1 = kfl_field(1,ifiel)
          ndim2 = kfl_field(4,ifiel)
          if(      kfl_field(2,ifiel) == NELEM_TYPE ) then
             call memory_alloca(memor_dom,'XFIEL','mod_meshes',xfiel_sub(ifiel) % a,ndim1,nelem_sub,ndim2)
          else if( kfl_field(2,ifiel) == NPOIN_TYPE ) then                                                                        
             call memory_alloca(memor_dom,'XFIEL','mod_meshes',xfiel_sub(ifiel) % a,ndim1,npoin_sub,ndim2)
          else if( kfl_field(2,ifiel) == NBOUN_TYPE ) then                                                                        
             call memory_alloca(memor_dom,'XFIEL','mod_meshes',xfiel_sub(ifiel) % a,ndim1,nboun_sub,ndim2)
          end if
       end do
       !
       ! Add field for groups
       !
       if( kfl_ngrou /= 0 ) then 
          nfiel = nfiel + 1
          ifiel = nfiel
          if( ifiel > mfiel ) call runend('MOD_MESHES: TOO MANY FIELDS!')
          kfl_field(1,ifiel) = 1
          kfl_field(2,ifiel) = NPOIN_TYPE
          kfl_field(4,ifiel) = 1
          ndim1              = kfl_field(1,ifiel)
          ndim2              = kfl_field(4,ifiel)
          call memory_alloca(memor_dom,'XFIEL','mod_meshes',xfiel_sub(ifiel) % a,ndim1,npoin_sub,ndim2)
       end if
    end if
    !
    ! ELement arrays
    !
    do ineig = 0,nneig_gat 
       lsiz4_gat(ineig) = int(mnode * nelem_gat(ineig),4)
    end do
    call PAR_GATHERV(lnods_sub,lnods_sub,lsiz4_gat,'IN THE WORLD')

    do ineig = 0,nneig_gat 
       lsiz4_gat(ineig) = int(nelem_gat(ineig),4)
    end do
    call PAR_GATHERV(ltype_sub,ltype_sub,lsiz4_gat,'IN THE WORLD')
    call PAR_GATHERV(lesub_sub,lesub_sub,lsiz4_gat,'IN THE WORLD')
    call PAR_GATHERV(lmate_sub,lmate_sub,lsiz4_gat,'IN THE WORLD')
    call PAR_GATHERV(lnnod_sub,lnnod_sub,lsiz4_gat,'IN THE WORLD')
    if( neset > 0 ) call PAR_GATHERV(leset_sub,leset_sub,lsiz4_gat,'IN THE WORLD')
    !
    ! Boundary arrays
    !
    do ineig = 0,nneig_gat 
       lsiz4_gat(ineig) = int(mnodb * nboun_gat(ineig),4)
    end do
    call PAR_GATHERV(lnodb_sub,lnodb_sub,lsiz4_gat,'IN THE WORLD')
    do ineig = 0,nneig_gat 
       lsiz4_gat(ineig) = int(nboun_gat(ineig),4)
    end do
    call PAR_GATHERV(ltypb_sub,ltypb_sub,lsiz4_gat,'IN THE WORLD')
    call PAR_GATHERV(lnnob_sub,lnnob_sub,lsiz4_gat,'IN THE WORLD')
    call PAR_GATHERV(lelbo_sub,lelbo_sub,lsiz4_gat,'IN THE WORLD')
    call PAR_GATHERV(codbo_sub,codbo_sub,lsiz4_gat,'IN THE WORLD')
    if( nbset > 0 ) call PAR_GATHERV(lbset_sub,lbset_sub,lsiz4_gat,'IN THE WORLD')
    !
    ! Node arrays
    !
    do ineig = 0,nneig_gat 
       lsiz4_gat(ineig) = int(ndime * npoin_gat(ineig),4)
    end do
    call PAR_GATHERV(coord_sub,coord_sub,lsiz4_gat,'IN THE WORLD')
    do ineig = 0,nneig_gat 
       lsiz4_gat(ineig) = int(npoin_gat(ineig),4)
    end do
    call PAR_GATHERV(lninv_sub,lninv_sub,lsiz4_gat,'IN THE WORLD')
    !
    ! Fields
    !
    do ifiel = 1,nfiel
       ndim1 = kfl_field(1,ifiel)
       ndim2 = kfl_field(4,ifiel)
       if(      kfl_field(2,ifiel) == NELEM_TYPE ) then
          do ineig = 0,nneig_gat 
             lsiz4_gat(ineig) = int(ndim1*ndim2*nelem_gat(ineig),4)
          end do
          call PAR_GATHERV(xfiel_sub(ifiel) % a,xfiel_sub(ifiel) % a,lsiz4_gat,'IN THE WORLD')
       else if( kfl_field(2,ifiel) == NPOIN_TYPE ) then                                                                        
          do ineig = 0,nneig_gat 
             lsiz4_gat(ineig) = int(ndim1*ndim2*npoin_gat(ineig),4)
          end do
          call PAR_GATHERV(xfiel_sub(ifiel) % a,xfiel_sub(ifiel) % a,lsiz4_gat,'IN THE WORLD')
       else if( kfl_field(2,ifiel) == NBOUN_TYPE ) then                                                                        
          do ineig = 0,nneig_gat 
             lsiz4_gat(ineig) = int(ndim1*ndim2*nboun_gat(ineig),4)
          end do
          call PAR_GATHERV(xfiel_sub(ifiel) % a,xfiel_sub(ifiel) % a,lsiz4_gat,'IN THE WORLD')
       end if
    end do
    !
    ! Renumber nodes
    !
    if( INOTSLAVE ) then
       ielem = 0
       iboun = 0
       kpoin = 0
       do ineig = 0,nneig_gat
          do kelem = 1,nelem_gat(ineig)
             ielem = ielem + 1
             lnods_sub(:,ielem) = lnods_sub(:,ielem) + kpoin
          end do
          do kboun = 1,nboun_gat(ineig)
             iboun = iboun + 1
             lnodb_sub(:,iboun) = lnodb_sub(:,iboun) + kpoin
          end do
          kpoin = kpoin + npoin_gat(ineig)
       end do
    end if
    !
    ! Renumber boundaries
    !
    if( INOTSLAVE ) then
       iboun = 0
       kelem = 0
       do ineig = 0,nneig_gat
          do kboun = 1,nboun_gat(ineig)
             iboun = iboun + 1
             lelbo_sub(iboun) = lelbo_sub(iboun) + kelem
          end do
          kelem = kelem + nelem_gat(ineig)
       end do
    end if
    !
    ! Eliminate duplicated nodes
    !
    if( INOTSLAVE ) then

       allocate(permu_node(npoin_sub))
       do ipoin = 1,npoin_sub
          permu_node(ipoin) = 0
       end do
       kpoin = 0
       do ipoin = 1,npoin_sub
          if( permu_node(ipoin) == 0 ) then
             kpoin = kpoin + 1
             permu_node(ipoin) = kpoin
             do jpoin = ipoin+1,npoin_sub
                if( lninv_sub(jpoin) == lninv_sub(ipoin) ) then
                   permu_node(jpoin) = -permu_node(ipoin)
                end if
             end do
          end if
       end do

       jpoin = 0
       do ipoin = 1,npoin_sub
          if( permu_node(ipoin) > 0 ) then
             jpoin = jpoin + 1
             coord_sub(:,jpoin) = coord_sub(:,ipoin)
          end if
       end do
       if( kfl_ngrou /= 0 ) then
          jpoin = 0
          do ipoin = 1,npoin_sub
             if( permu_node(ipoin) > 0 ) then
                jpoin = jpoin + 1
                xfiel_sub(nfiel) % a(1,jpoin,1) = xfiel_sub(nfiel) % a(1,ipoin,1)
             end if
          end do
       end if

       npoin_sub = kpoin
       do ielem = 1,nelem_sub
          lnods_sub(:,ielem) = abs(permu_node(lnods_sub(:,ielem)))
       end do
       do iboun = 1,nboun_sub
          lnodb_sub(:,iboun) = abs(permu_node(lnodb_sub(:,iboun)))
       end do

    end if

    !--------------------------------------------------------------------
    !
    ! Add a field for parallelization
    !
    !--------------------------------------------------------------------

    if( INOTSLAVE ) then
       nfiel = nfiel + 1
       ifiel = nfiel
       if( ifiel > mfiel ) call runend('MOD_MESHES: TOO MANY FIELDS!')
       kfl_field(1,ifiel) = 1
       kfl_field(2,ifiel) = NELEM_TYPE
       kfl_field(4,ifiel) = 1
       ndim1              = kfl_field(1,ifiel)
       ndim2              = kfl_field(4,ifiel)
       call memory_alloca(memor_dom,'XFIEL','mod_meshes',xfiel_sub(ifiel) % a,ndim1,nelem_sub,ndim2)
       ielem = 0
       kneig = 0
       do ineig = 0,nneig_gat
          if( nelem_gat(ineig) > 0 ) then
             kneig = kneig + 1
             do kelem = 1,nelem_gat(ineig)
                ielem = ielem + 1
                xfiel_sub(ifiel) % a(1,ielem,1) = real(kneig,rp)
             end do
          end if
       end do
    end if
    if( IMASTER ) then
       call messages_live('NUMBER OF PARTITIONS INVOLVED IN THIS SUBMESH= ',INT_NUMBER=kneig)       
    end if

    !--------------------------------------------------------------------
    !
    ! Output mesh
    !
    !--------------------------------------------------------------------

    nunit_msh = 90
    nunit_res = 91
    nunit_dat = 92
    call iofile_open_unit(nunit_msh,trim(title)//'-submesh.post.msh','POSTPROCESS SUBMESH')
    call iofile_open_unit(nunit_res,trim(title)//'-submesh.post.res','POSTPROCESS SUBMESH')
    call iofile_open_unit(nunit_dat,trim(title)//'-submesh.dom.dat', 'POSTPROCESS SUBMESH')

    if( INOTSLAVE ) then
       meshe_sub % ndime     =  ndime
       meshe_sub % mnode     =  mnode
       meshe_sub % mnodb     =  mnodb
       meshe_sub % nfiel     =  nfiel
       meshe_sub % kfl_field =  kfl_field
       meshe_sub % nelem     =  nelem_sub
       meshe_sub % nboun     =  nboun_sub
       meshe_sub % npoin     =  npoin_sub
       if( kfl_ngrou /= 0 ) meshe_sub % kfl_ngrou =  nfiel-1
       meshe_sub % lnods     => lnods_sub
       meshe_sub % ltype     => ltype_sub
       meshe_sub % lesub     => lesub_sub
       meshe_sub % lmate     => lmate_sub
       meshe_sub % lnnod     => lnnod_sub
       meshe_sub % leset     => leset_sub
       meshe_sub % lnodb     => lnodb_sub
       meshe_sub % ltypb     => ltypb_sub
       meshe_sub % lnnob     => lnnob_sub
       meshe_sub % lelbo     => lelbo_sub
       meshe_sub % kfl_codbo => codbo_sub
       meshe_sub % lbset     => lbset_sub
       meshe_sub % coord     => coord_sub
       do ifiel = 1,nfiel
          meshe_sub % xfiel(ifiel) % a => xfiel_sub(ifiel) % a
       end do
       call output_mesh_gid_format   (meshe_sub,'SUBMESH',nunit_msh)
       call output_domain_alya_format(meshe_sub,nunit_dat)

       allocate(xarray(nelem_sub))
       do ielem = 1,nelem_sub
          xarray(ielem) = meshe_sub % xfiel(nfiel)%a(1,ielem,1)
       end do
       call output_result_gid_format (nunit_res,nelem_sub,xarray,NAME='PARTITION',wherein='ON ELEMENTS')
       deallocate(xarray)
       if( meshe_sub % kfl_ngrou > 0 ) then
          allocate(xarray(npoin_sub))
          do ipoin = 1,npoin_sub
             xarray(ipoin) = meshe_sub % xfiel(nfiel-1)%a(1,ipoin,1)
          end do
          call output_result_gid_format (nunit_res,npoin_sub,xarray,NAME='GROUPS')
          deallocate(xarray)
       end if
    end if

    !--------------------------------------------------------------------
    !
    ! Deallocate 
    !
    !--------------------------------------------------------------------

    call memory_deallo(memor_dom,'LNODS','mod_meshes',lnods_sub)
    call memory_deallo(memor_dom,'LTYPE','mod_meshes',ltype_sub)
    call memory_deallo(memor_dom,'LESUB','mod_meshes',lesub_sub)
    call memory_deallo(memor_dom,'LMATE','mod_meshes',lmate_sub)
    call memory_deallo(memor_dom,'LNNOD','mod_meshes',lnnod_sub)
    call memory_deallo(memor_dom,'LESET','mod_meshes',leset_sub)

    call memory_deallo(memor_dom,'LNODB','mod_meshes',lnodb_sub)
    call memory_deallo(memor_dom,'LTYPB','mod_meshes',ltypb_sub)
    call memory_deallo(memor_dom,'LNNOB','mod_meshes',lnnob_sub)
    call memory_deallo(memor_dom,'LELBO','mod_meshes',lelbo_sub)
    call memory_deallo(memor_dom,'CODBO','mod_meshes',codbo_sub)
    call memory_deallo(memor_dom,'LBSET','mod_meshes',lbset_sub)

    call memory_deallo(memor_dom,'COORD','mod_meshes',coord_sub)
    call memory_deallo(memor_dom,'LNINV','mod_meshes',lninv_sub)
    do ifiel = 1,mfiel
       call memory_deallo(memor_dom,'XFIEL','mod_meshes',xfiel_sub(ifiel) % a)
    end do

    if( associated(mark_node ) ) deallocate(mark_node )
    if( associated(mark_elem ) ) deallocate(mark_elem )
    if( associated(mark_boun ) ) deallocate(mark_boun )
    if( associated(permu_node) ) deallocate(permu_node)
    if( associated(permu_elem) ) deallocate(permu_elem)
    if( associated(permu_boun) ) deallocate(permu_boun)
    if( associated(invpe_node) ) deallocate(invpe_node)
    if( associated(invpe_elem) ) deallocate(invpe_elem)

    call iofile_close_unit(nunit_msh)
    call iofile_close_unit(nunit_dat)
    call iofile_close_unit(nunit_res)
    nfiel = nfiel_sav 

    call runend('O.K.!')

  end subroutine meshes_gather_submesh
  
end module mod_meshes_gather_submesh
!> @}
