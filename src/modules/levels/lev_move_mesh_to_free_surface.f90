!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



 subroutine lev_move_mesh_to_free_surface()
!!$  use def_kintyp,     only : ip,rp
!!$  use def_master
!!$  use def_domain
!!$  use def_elmtyp
!!$  use mod_cutele
!!$  use def_kermod,         only : kfl_cutel 
!!$  use def_kermod,         only : ndivi
!!$  use mod_ker_proper,     only : ker_updpro
!!$  use mod_kdtree,         only : typ_kdtree
!!$  use mod_meshes,         only : meshes_surface_from_nodal_array
!!$  use mod_meshes,         only : meshes_surface_from_nodal_array_deallocate
!!$  use mod_communications, only : PAR_INTERFACE_NODE_EXCHANGE
!!$  use mod_communications, only : PAR_MAX
!!$  use mod_communications, only : PAR_ALLGATHER
!!$  use mod_communications, only : PAR_ALLGATHERV
!!$  use mod_parall,         only : PAR_CODE_SIZE
!!$  use mod_ker_deform,     only : deform_deform
!!$  use def_levels 
!!$  use mod_memory
!!$  use mod_postpr
!!$  implicit none
!!$  integer(ip)            :: ielem,inode,ipoin
!!$  integer(ip)            :: knode,pnode,jpoin
!!$  integer(ip)            :: idime,jnode,pelty
!!$  integer(ip)            :: iboun,ibopo,ipart
!!$  integer(ip)            :: inodb
!!$  real(rp)               :: rsign,norma(3)
!!$  real(rp)               :: fi,fj,t,xx(ndime,mnode)
!!$  real(rp)               :: pcoor(3),dummr,dimax
!!$  integer(ip), pointer   :: linte(:)
!!$  type(typ_kdtree)       :: kdtree_lev
!!$  integer(ip), parameter :: AIR_NODE=1
!!$  integer(ip), parameter :: WATER_NODE=2
!!$  integer(ip)            :: FRINGE_NODE
!!$
!!$  integer(ip)            :: istar
!!$  integer(ip)            :: mnodb_sur
!!$  integer(ip)            :: npoin_sur
!!$  integer(ip)            :: nboun_sur
!!$  integer(ip)            :: npoin_all
!!$  integer(ip)            :: nboun_all
!!$  integer(ip), pointer   :: npoin_nboun_sur(:)
!!$  integer(ip), pointer   :: npoin_nboun_glo(:)
!!$  integer(ip), pointer   :: npoin_glo(:)
!!$  integer(ip), pointer   :: nboun_glo(:)
!!$  integer(ip), pointer   :: lnodb_sur(:,:) , lnodb_glo(:,:)
!!$  integer(ip), pointer   :: ltypb_sur(:)   , ltypb_glo(:)
!!$  real(rp),    pointer   :: coord_sur(:,:) , coord_glo(:,:)
!!$  real(rp),    pointer   :: dista_sur(:)   , dista_glo(:)
!!$  real(rp),    pointer   :: proje_sur(:,:) , proje_glo(:,:)
!!$
!!$  real(rp),    pointer   :: dista_cpy(:)
!!$  real(rp)               :: dista,proje(3)
!!$
!!$  integer(ip)            :: ndefo,kfl_dmeth
!!$  integer(ip), pointer   :: kfl_fixno_sur(:,:) 
!!$  real(rp),    pointer   :: bvess_sur(:,:)
!!$
!!$  character(100)         :: wmess
!!$  nullify(linte)
!!$  nullify(dista_cpy)
!!$
!!$  nullify(npoin_nboun_sur)
!!$  nullify(dista_sur)
!!$  nullify(proje_sur)
!!$  nullify(ltypb_sur)
!!$  nullify(lnodb_sur)
!!$  nullify(ltypb_sur)
!!$  nullify(coord_sur)
!!$  
!!$  nullify(npoin_nboun_glo)
!!$  nullify(npoin_glo)
!!$  nullify(nboun_glo)
!!$  nullify(dista_glo)
!!$  nullify(proje_glo)
!!$  nullify(ltypb_glo)
!!$  nullify(lnodb_glo)
!!$  nullify(ltypb_glo)
!!$  nullify(coord_glo)
!!$  
!!$  nullify(kfl_fixno_sur)
!!$  nullify(bvess_sur)
!!$
!!$  npoin_sur = 0
!!$  nboun_sur = 0
!!$  mnodb_sur = ndime
!!$
!!$  !----------------------------------------------------------------------
!!$  !
!!$  ! Boundary condition
!!$  !
!!$  !----------------------------------------------------------------------
!!$
!!$  if( INOTMASTER ) then
!!$     call memory_alloca(mem_modul(1:2,modul),'KFL_FIXNO_SUR','lev_move_mesh_to_free_surface',kfl_fixno_sur,ndime,npoin)
!!$     call memory_alloca(mem_modul(1:2,modul),'BVESS_SUR'    ,'lev_move_mesh_to_free_surface',bvess_sur,ndime,npoin)
!!$     do ipoin = 1,npoin
!!$        ibopo = lpoty(ipoin)
!!$        if( ibopo > 0 ) then
!!$           kfl_fixno_sur(1:ndime,ipoin) = 1
!!$           bvess_sur(1:ndime,ipoin)     = 0.0_rp
!!$        end if
!!$     end do
!!$     !
!!$     ! Recover old coordinates
!!$     !
!!$     do ipoin = 1,npoin
!!$        do idime = 1,ndime
!!$           if( kfl_fixno_sur(idime,ipoin) == 2 ) then
!!$              kfl_fixno_sur(idime,ipoin) = 0
!!$              bvess_sur(idime,ipoin)     = 0.0_rp
!!$           end if
!!$           !coord(idime,ipoin) = coord_ori(idime,ipoin)
!!$        end do
!!$        dispm(1:ndime,ipoin,1) = 0.0_rp
!!$        velom(1:ndime,ipoin)   = 0.0_rp        
!!$     end do
!!$  end if
!!$
!!$  !----------------------------------------------------------------------
!!$  !
!!$  ! Detect interface nodes
!!$  !
!!$  !----------------------------------------------------------------------
!!$
!!$  if( INOTMASTER ) then
!!$
!!$     call memory_alloca(mem_modul(1:2,modul),'LINTE'    ,'lev_move_mesh_to_free_surface',linte,npoin)
!!$     call memory_alloca(mem_modul(1:2,modul),'DISTA_SUR','lev_move_mesh_to_free_surface',dista_sur,npoin)
!!$     call memory_alloca(mem_modul(1:2,modul),'DISTA_CPY','lev_move_mesh_to_free_surface',dista_cpy,npoin)
!!$     call memory_alloca(mem_modul(1:2,modul),'PROJE_SUR','lev_move_mesh_to_free_surface',proje_sur,ndime,npoin)
!!$
!!$     do ielem = 1,nelem
!!$        if( lelch(ielem) /= ELHOL ) then   
!!$           ipoin = lnods(1,ielem)
!!$           rsign = fleve(ipoin,1)
!!$           pnode = lnnod(ielem)
!!$           loop_inode: do inode = 2,pnode
!!$              ipoin = lnods(inode,ielem)
!!$              if( rsign*fleve(ipoin,1) < 0.0_rp ) then
!!$                 do jnode = 1,pnode
!!$                    jpoin = lnods(jnode,ielem)
!!$                    if( fleve(jpoin,1) < 0.0_rp ) then
!!$                       linte(jpoin) = AIR_NODE
!!$                    else
!!$                       linte(jpoin) = WATER_NODE
!!$                    end if
!!$                 end do
!!$                 exit loop_inode
!!$              end if
!!$           end do loop_inode
!!$        end if
!!$     end do
!!$  end if
!!$
!!$  !----------------------------------------------------------------------
!!$  !
!!$  ! Construct local level-set surface mesh using current level-set function
!!$  !
!!$  !----------------------------------------------------------------------
!!$
!!$  if( INOTMASTER ) then
!!$     call meshes_surface_from_nodal_array(&
!!$          fleve,meshe(ndivi),npoin_sur,nboun_sur,lnodb_sur,coord_sur,ltypb_sur)
!!$  end if
!!$
!!$  !----------------------------------------------------------------------
!!$  !
!!$  ! Gather global level-set surface mesh from local ones 
!!$  !
!!$  !----------------------------------------------------------------------
!!$  !
!!$  ! Gather level set mesh parameters
!!$  ! npoin_sur of kfl_paral = npoin_nboun_glo( 2*ipart+0 )
!!$  ! nboun_sur of kfl_paral = npoin_nboun_glo( 2*ipart+1 )
!!$  ! NB: master is included
!!$  !
!!$  call memory_alloca(mem_modul(1:2,modul),'NBOUN_SUR','lev_move_mesh_to_free_surface',npoin_nboun_sur,2_ip)   
!!$  call memory_alloca(mem_modul(1:2,modul),'NBOUN_ALL','lev_move_mesh_to_free_surface',npoin_nboun_glo,2_ip*PAR_CODE_SIZE,'INITIALIZE',0_ip)   
!!$  call memory_alloca(mem_modul(1:2,modul),'NPOIN_GLO','lev_move_mesh_to_free_surface',npoin_glo,PAR_CODE_SIZE,'INITIALIZE',0_ip)
!!$  call memory_alloca(mem_modul(1:2,modul),'NBOUN_GLO','lev_move_mesh_to_free_surface',nboun_glo,PAR_CODE_SIZE,'INITIALIZE',0_ip)
!!$  npoin_nboun_sur(1) = npoin_sur
!!$  npoin_nboun_sur(2) = nboun_sur
!!$  call PAR_ALLGATHER(npoin_nboun_sur,npoin_nboun_glo,2_ip,'IN MY CODE')
!!$  npoin_all = 0
!!$  nboun_all = 0  
!!$  do ipart = 0,PAR_CODE_SIZE-1
!!$     npoin_glo(ipart) = npoin_nboun_glo( 2*ipart+0 )
!!$     nboun_glo(ipart) = npoin_nboun_glo( 2*ipart+1 )
!!$     npoin_all        = npoin_all + npoin_nboun_glo( 2*ipart+0 )
!!$     nboun_all        = nboun_all + npoin_nboun_glo( 2*ipart+1 )
!!$  end do
!!$  !
!!$  ! Slave renumber their boundaries to obtain a global node numbering
!!$  !
!!$  istar = 0
!!$  do ipart = 1,kfl_paral-1
!!$     istar = istar + npoin_glo(ipart)
!!$  end do
!!$  do iboun = 1,nboun_sur
!!$     lnodb_sur(1:mnodb_sur,iboun) = lnodb_sur(1:mnodb_sur,iboun) + istar
!!$  end do
!!$  !
!!$  ! Gather level set mesh arrays
!!$  !
!!$  call memory_alloca(mem_modul(1:2,modul),'LNODB_GLO','lev_move_mesh_to_free_surface',lnodb_glo,mnodb_sur,nboun_all)   
!!$  call memory_alloca(mem_modul(1:2,modul),'LTYPB_GLO','lev_move_mesh_to_free_surface',ltypb_glo,nboun_all)   
!!$  call memory_alloca(mem_modul(1:2,modul),'COORD_GLO','lev_move_mesh_to_free_surface',coord_glo,ndime,npoin_all)
!!$  do ipart = 1,PAR_CODE_SIZE-1
!!$     npoin_glo(ipart) = npoin_glo(ipart) * ndime
!!$  end do
!!$  call PAR_ALLGATHERV(coord_sur,coord_glo,npoin_glo,'IN MY CODE')
!!$  call PAR_ALLGATHERV(ltypb_sur,ltypb_glo,nboun_glo,'IN MY CODE')
!!$  do ipart = 1,PAR_CODE_SIZE-1
!!$     nboun_glo(ipart) = nboun_glo(ipart) * mnodb_sur
!!$  end do
!!$  call PAR_ALLGATHERV(lnodb_sur,lnodb_glo,nboun_glo,'IN MY CODE')
!!$  
!!$  wmess = 'LEVELS: LEVET-SET SURFACE MESH: NPOIN, NBOUN= '//trim(intost(npoin_all))//', '//trim(intost(nboun_all))
!!$  call livinf(0_ip,trim(wmess),0_ip)
!!$
!!$  if(kfl_paral==0) then
!!$     call gid_mesh(99_4,ndime,npoin_all,nboun_all,mnodb_sur,ltypb_glo,coord_glo,lnodb_glo)
!!$  end if
!!$
!!$  !----------------------------------------------------------------------
!!$  !
!!$  ! Construct global KDTree
!!$  !
!!$  !----------------------------------------------------------------------
!!$
!!$  if( INOTMASTER ) then
!!$     call kdtree_construct(&
!!$          nboun_all,npoin_all,lnodb_glo,ltypb_glo,coord_glo,kdtree_lev)
!!$  end if
!!$
!!$  !----------------------------------------------------------------------
!!$  !
!!$  ! Compute projection of fringe nodes on level-set surface mesh
!!$  !
!!$  !----------------------------------------------------------------------
!!$
!!$  !if( kfl_paral==1) then
!!$  !   pcoor(1)= 0.15_rp
!!$  !   pcoor(2)= 1.9134
!!$  !   ipoin = 1
!!$  !   call kdtree_nearest_boundary(pcoor,kdtree_lev,iboun,dista_sur(ipoin),proje_sur(:,ipoin))
!!$  !   print*,'TESTEO=',iboun,dista_sur(ipoin),proje_sur(1:2,ipoin)
!!$  !end if
!!$
!!$  FRINGE_NODE = AIR_NODE
!!$
!!$  if( INOTMASTER ) then
!!$     do ipoin = 1,npoin
!!$        if( linte(ipoin) == FRINGE_NODE ) then
!!$           call kdtree_nearest_boundary(coord(:,ipoin),kdtree_lev,iboun,dista_sur(ipoin),proje_sur(:,ipoin))
!!$           do idime = 1,ndime
!!$              if( kfl_fixno_sur(idime,ipoin) == 0 ) then
!!$                 kfl_fixno_sur(idime,ipoin) = 2
!!$                 bvess_sur(idime,ipoin)     = proje_sur(idime,ipoin)-coord(idime,ipoin)
!!$                 dispm(idime,ipoin,1)       = bvess_sur(idime,ipoin)
!!$              end if
!!$           end do
!!$        end if
!!$     end do
!!$  end if
!!$
!!$  !----------------------------------------------------------------------
!!$  !
!!$  ! Deform mesh
!!$  !
!!$  !----------------------------------------------------------------------
!!$
!!$  solve_sol => solve(4:)
!!$  ndefo = 1
!!$  kfl_dmeth = 6
!!$  call deform_deform(&
!!$       ndefo,kfl_dmeth,ID_LEVELS,kfl_fixno_sur,bvess_sur,&
!!$       coord,amatr,unkno,rhsid,solve_sol)
!!$  solve_sol => solve(1:)
!!$  do ipoin = 1,npoin 
!!$     do idime = 1,ndime
!!$        dispm(idime,ipoin,1) = unkno( (ipoin-1)*ndime + idime ) 
!!$        coord(idime,ipoin)   = coord(idime,ipoin) + dispm(idime,ipoin,1)
!!$        velom(idime,ipoin)   = dtinv * dispm(idime,ipoin,1)
!!$     end do
!!$  end do
!!$  !
!!$  ! Maximum distance
!!$  !
!!$  !call PAR_MAX(dimax,'IN MY CODE')
!!$  !if( INOTSLAVE ) print*,'MAX DISTANCE=',dimax
!!$  !if(kfl_paral==0) then
!!$  !   !call gid_mesh(99_4,ndime,npoin,nboun_all,2_ip,coord_glo,lnodb_glo)
!!$  !   call gid_mesh(99_4,ndime,npoin,nelem,mnode,ltype,coord,lnods)
!!$  !end if
!!$  !call runend('O.K.!')
!!$  !call memgen(0_ip,ndime,npoin)
!!$  !do ipoin = 1,npoin
!!$  !   gevec(1:ndime,ipoin) = dispm(1:ndime,ipoin,1)
!!$  !end do
!!$  !call postpr_right_now('XXXXX','VECTO','NPOIN',gevec)
!!$  !call runend('O.K.!')
!!$
!!$  !----------------------------------------------------------------------
!!$  !
!!$  ! Redistanciation: Recompute geometrical level set
!!$  !
!!$  !----------------------------------------------------------------------
!!$
!!$  do ipoin = 1,npoin
!!$     call kdtree_nearest_boundary(coord(:,ipoin),kdtree_lev,iboun,dista_sur(ipoin),proje_sur(:,ipoin))
!!$     fleve(ipoin,1) = dista_sur(ipoin)
!!$     if( fleve(ipoin,1) >= 0.0_rp ) fleve(ipoin,1) = max(fleve(ipoin,1),1.0e-8_rp)
!!$  end do
!!$  !call memgen(0_ip,npoin,0_ip)
!!$  !do ipoin = 1,npoin
!!$  !   gesca(ipoin) = fleve(ipoin,1)
!!$  !end do
!!$  !call postpr_right_now('XXXXX','SCALA','NPOIN',gesca)
!!$  !call runend('O.K.!')
!!$
!!$  !
!!$  ! Deallocate
!!$  !
!!$  call kdtree_deallocate(kdtree_lev)
!!$  call meshes_surface_from_nodal_array_deallocate(lnodb_sur,coord_sur,ltypb_sur)
!!$  call memory_deallo(mem_modul(1:2,modul),'LINTE'        ,'lev_move_mesh_to_free_surface',linte)
!!$  call memory_deallo(mem_modul(1:2,modul),'DISTA_SUR'    ,'lev_move_mesh_to_free_surface',dista_sur)
!!$  call memory_deallo(mem_modul(1:2,modul),'DISTA_CPY'    ,'lev_move_mesh_to_free_surface',dista_cpy)
!!$  call memory_deallo(mem_modul(1:2,modul),'PROJE_SUR'    ,'lev_move_mesh_to_free_surface',proje_sur)
!!$  call memory_deallo(mem_modul(1:2,modul),'KFL_FIXNO_SUR','lev_move_mesh_to_free_surface',kfl_fixno_sur)
!!$  call memory_deallo(mem_modul(1:2,modul),'BVESS_SUR'    ,'lev_move_mesh_to_free_surface',bvess_sur)
!!$  
!!$  call memory_deallo(mem_modul(1:2,modul),'NBOUN_SUR'    ,'lev_move_mesh_to_free_surface',npoin_nboun_sur)   
!!$  call memory_deallo(mem_modul(1:2,modul),'NBOUN_ALL'    ,'lev_move_mesh_to_free_surface',npoin_nboun_glo)   
!!$  call memory_deallo(mem_modul(1:2,modul),'NPOIN_GLO'    ,'lev_move_mesh_to_free_surface',npoin_glo)
!!$  call memory_deallo(mem_modul(1:2,modul),'NBOUN_GLO'    ,'lev_move_mesh_to_free_surface',nboun_glo)
!!$  call memory_deallo(mem_modul(1:2,modul),'LNODB_GLO'    ,'lev_move_mesh_to_free_surface',lnodb_glo)   
!!$  call memory_deallo(mem_modul(1:2,modul),'LTYPB_GLO'    ,'lev_move_mesh_to_free_surface',ltypb_glo)   
!!$  call memory_deallo(mem_modul(1:2,modul),'COORD_GLO'    ,'lev_move_mesh_to_free_surface',coord_glo)
!!$
!!$  !----------------------------------------------------------------------
!!$  !
!!$  ! Recompute geometrical arrays
!!$  !
!!$  !----------------------------------------------------------------------
!!$
!!$  call domarr(2_ip)
!!$
!!$  !----------------------------------------------------------------------
!!$  !
!!$  ! Force to update properties
!!$  !
!!$  !----------------------------------------------------------------------
!!$
!!$  call ker_updpro() 

end subroutine lev_move_mesh_to_free_surface




subroutine gid_mesh(lunit,ndime,npoin,nelem,mnode,ltype,coord,lnods)
  use def_kintyp, only : ip,rp
  integer(4),  intent(in) :: lunit
  integer(ip), intent(in) :: ndime
  integer(ip), intent(in) :: npoin
  integer(ip), intent(in) :: nelem
  integer(ip), intent(in) :: mnode
  integer(ip), intent(in) :: ltype(nelem)
  integer(ip), intent(in) :: lnods(mnode,nelem)
  real(rp),    intent(in) :: coord(ndime,npoin)
  integer(ip)             :: ipoin,ielem,pnode
  
  open(unit=lunit,file='mesh_on_the_fly.post.msh',status='unknown')
  !write(lunit,'(a,1x,i1,1x,a)') 'MESH ON_THE_FLY dimension',ndime,'Elemtype Triangle Nnode 3'
  write(lunit,'(a,1x,i1,1x,a)') 'MESH ON_THE_FLY dimension',ndime,'Elemtype Linear Nnode 2'
  write(lunit,*) 'coordinates'
  do ipoin = 1,npoin
     write(lunit,*) ipoin,coord(1:ndime,ipoin)
  end do
  write(lunit,*) 'end coordinates'
  write(lunit,*) 'elements'
  do ielem = 1,nelem     
     pnode = ltype(ielem)
     write(lunit,'(10(1x,i7))') ielem,lnods(1:pnode,ielem)
  end do
  write(lunit,*) 'end elements'
  close(lunit)

end subroutine gid_mesh
