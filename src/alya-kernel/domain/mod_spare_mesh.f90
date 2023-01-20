!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Domain
!> @{
!> @file    mod_sparce_mesh.f90
!> @author  houzeaux
!> @date    2020-01-28
!> @brief   Spare meshes 
!> @details Spare meshes
!-----------------------------------------------------------------------

module mod_spare_mesh

  use def_kintyp_basic,                  only : ip,i1p,r2p,lg,rp
  use def_kintyp_spare_mesh,             only : typ_spare_mesh 
  use def_elmtyp,                        only : POI3D
  use def_master,                        only : INOTMASTER
  use def_master,                        only : IPARALL
  use def_kermod,                        only : ndivi
  use def_kermod,                        only : spare_meshes
  use def_kermod,                        only : num_spare_meshes
  use def_kermod,                        only : ielse,relse
  use def_domain,                        only : meshe
  use def_kintyp_mesh_basic,             only : mesh_type_basic
  use def_kintyp_mesh,                   only : mesh_type
  use mod_memory_basic,                  only : memory_alloca
  use mod_memory_basic,                  only : memory_deallo
  use mod_mesh_type_basic,               only : mesh_type_basic_broadcast
  use mod_elmgeo,                        only : element_type
  use mod_communications_point_to_point, only : PAR_SEND_RECEIVE
  use mod_elsest,                        only : elsest_host_element
  use mod_kdtree,                        only : kdtree_construct
  use mod_kdtree,                        only : typ_kdtree
  use mod_kdtree,                        only : kdtree_nearest_boundary
  use mod_kdtree,                        only : kdtree_deallocate
  use mod_kdtree,                        only : kdtree_initialize
  use def_mpi
  
  implicit none

  private
  public :: spare_mesh_setup
  public :: spare_mesh_parall
  public :: spare_mesh_alloca

contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-02-14
  !> @brief   Set up spare meshes info
  !> @details Set up spare meshes info
  !> 
  !-----------------------------------------------------------------------

  subroutine spare_mesh_alloca()

    integer(ip) :: ii

    if( .not. associated(spare_meshes) ) then
       allocate(spare_meshes(num_spare_meshes))
       do ii = 1,num_spare_meshes
          call spare_meshes(ii) % init()
       end do
    end if
    
  end subroutine spare_mesh_alloca
           
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-02-14
  !> @brief   Set up spare meshes info
  !> @details Set up spare meshes info
  !> 
  !-----------------------------------------------------------------------

  subroutine spare_mesh_setup()

    integer(ip) :: ii

    do ii = 1,num_spare_meshes
       call spare_mesh_gauss_points(spare_meshes(ii),spare_meshes(ii) % mesh,meshe(ndivi))
       call spare_mesh_distance    (spare_meshes(ii),spare_meshes(ii) % mesh,meshe(ndivi))
    end do

  end subroutine spare_mesh_setup

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-02-14
  !> @brief   Set up spare meshes info
  !> @details Set up spare meshes info
  !> 
  !-----------------------------------------------------------------------

  subroutine spare_mesh_parall()

    integer(ip) :: ii

    if( num_spare_meshes > 0 ) then
       if( INOTMASTER ) call spare_mesh_alloca()
       do ii = 1,num_spare_meshes
          call mesh_type_basic_broadcast(spare_meshes(ii) % mesh)
       end do
    end if

  end subroutine spare_mesh_parall

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-05-14
  !> @brief   Check who is in charge of Gauss-point
  !> @details Check which rank is in charge of Gauss-point in a 
  !>          background mesh mesh_bak
  !> 
  !-----------------------------------------------------------------------
  
  subroutine spare_mesh_gauss_points(spare,mesh,mesh_bak)

    type(typ_spare_mesh),           intent(inout) :: spare
    type(mesh_type_basic),          intent(inout) :: mesh
    type(mesh_type),                intent(in)    :: mesh_bak
    integer(ip)                                   :: igaus,tgaus,ineig
    integer(ip)                                   :: ielem,kgaus,pnode
    integer(ip)                                   :: pelty,pgaus,ipoin
    integer(ip)                                   :: inode,ii,jelem,jnode
    integer(ip)                                   :: ndime,dom_i,my_rank
    real(rp)                                      :: coglo(3),coloc(3)
    real(rp)                                      :: deriv(mesh_bak % ndime,mesh_bak % mnode)
    real(rp)                                      :: shapf(mesh_bak % mnode)
    real(rp)                                      :: dista
    integer(ip),           pointer                :: list_gauss(:)
    logical(lg),           pointer                :: list_frine(:)
    integer(ip),           pointer                :: lsend_fring(:)
    type(i1p),             pointer                :: lrecv_fring(:)
    integer(ip)                                   :: nsend_fring
    integer(ip),           pointer                :: nrecv_fring(:)
    
    if( mesh % nelem <= 0 .or. mesh_bak % nelem <= 0 ) return    

    nullify(list_gauss)
    nullify(list_frine)
    nullify(lsend_fring)
    nullify(lrecv_fring)
    nullify(nrecv_fring)

    tgaus = 0
    do ielem = 1,mesh % nelem
       pelty = mesh % ltype(ielem)
       pgaus = mesh % quad(pelty) % ngaus
       tgaus = tgaus + pgaus
    end do
    call memory_alloca(spare % memor,'LIST_GAUSS' ,'mesh_type_basic_gauss_points',list_gauss,tgaus)        
    call memory_alloca(spare % memor,'SPARE % ELEME' ,'mesh_type_basic_gauss_points',spare % eleme,mesh % nelem)        
    call memory_alloca(spare % memor,'SPARE % SHAPF' ,'mesh_type_basic_gauss_points',spare % shapf,mesh % nelem)
    !
    ! List of fringe elements
    !
    if( IPARALL ) then
       call memory_alloca(spare % memor,'LIST_FRING','mesh_type_basic_gauss_points',list_frine,mesh_bak % nelem)        
       do ielem = 1,mesh_bak % nelem
          pnode = mesh_bak % lnnod(ielem)
          if( any( mesh_bak % lnods(1:pnode,ielem) > mesh_bak % npoi1 ) ) list_frine(ielem) = .true.
       end do
       call memory_alloca(spare % memor,'LSEND_FRING','mesh_type_basic_gauss_points',lsend_fring,max(1_ip,tgaus))        
       call memory_alloca(spare % memor,'NRECV_FRING','mesh_type_basic_gauss_points',nrecv_fring,mesh_bak % comm % nneig)        
       call memory_alloca(spare % memor,'LRECV_FRING','mesh_type_basic_gauss_points',lrecv_fring,mesh_bak % comm % nneig)        
       nsend_fring = 0
       nrecv_fring = 0
    end if
    
    kgaus   = 0  
    ndime   = mesh_bak % ndime
    my_rank = int(mesh_bak % comm % RANK4,ip)
    
    do ielem = 1,mesh % nelem
       pelty    = mesh % ltype(ielem)
       pgaus    = mesh % quad(pelty) % ngaus

       pnode    = element_type(pelty) % number_nodes
       do igaus = 1,pgaus
          kgaus = kgaus + 1
          coglo = 0.0_rp
          do inode = 1,pnode
             ipoin = mesh % lnods(inode,ielem)
             coglo(1:mesh % ndime) = coglo(1:ndime) &
                  + mesh % iso(pelty) % shape(inode,igaus) &
                  * mesh % coord(1:ndime,ipoin)
          end do
          call elsest_host_element(&
               ielse,relse,1_ip,mesh_bak,coglo,jelem,&
               shapf,deriv,coloc,dista)
          if( jelem /= 0 ) then
             jnode = mesh_bak % lnnod(jelem)
             if( .not. associated(spare % shapf(ielem)%a) ) &
                  call memory_alloca(spare % memor,'SPARE % SHAPF % A',&
                  'mesh_type_basic_gauss_points',spare % shapf(ielem) % a,jnode,pgaus)
             if( .not. associated(spare % eleme(ielem)%l) ) &
                  call memory_alloca(spare % memor,'SPARE % ELEME % L',&
                  'mesh_type_basic_gauss_points',spare % eleme(ielem) % l,pgaus)
             list_gauss(kgaus)                       = my_rank
             spare % eleme(ielem) % l(igaus)         = jelem
             spare % shapf(ielem) % a(1:jnode,igaus) = shapf(1:jnode)
             if( IPARALL ) then
                if( list_frine(jelem) )then
                   nsend_fring  = nsend_fring + 1
                   lsend_fring(nsend_fring) = kgaus
                end if
             end if
          end if
       end do
    end do
    
    if( IPARALL ) then
       do ineig = 1,mesh_bak % comm % nneig
          dom_i = mesh_bak % comm % neights(ineig)
          call PAR_SEND_RECEIVE(nsend_fring,nrecv_fring(ineig),DOM_I=dom_i,&
               PAR_COMM_IN=mesh_bak % comm % PAR_COMM_WORLD)
          call memory_alloca(spare % memor,'LRECV_FRING','mesh_type_basic_gauss_points',lrecv_fring(ineig) % l,max(1_ip,nrecv_fring(ineig)))
       end do

       do ineig = 1,mesh_bak % comm % nneig
          dom_i = mesh_bak % comm % neights(ineig)
          call PAR_SEND_RECEIVE(nsend_fring,nrecv_fring(ineig),lsend_fring,lrecv_fring(ineig) % l,DOM_I=dom_i,&
               PAR_COMM_IN=mesh_bak % comm % PAR_COMM_WORLD)
       end do

       do ineig = 1,mesh_bak % comm % nneig
          do ii = 1,nrecv_fring(ineig)
             kgaus = lrecv_fring(ineig) % l(ii)
             if( list_gauss(kgaus) > 0 ) then
                if( my_rank < mesh_bak % comm % neights(ineig) ) then
                   list_gauss(kgaus) = 0
                end if
             end if
          end do
       end do

       kgaus = 0
       do ielem = 1,mesh % nelem
          pgaus = mesh % quad(mesh % ltype(ielem)) % ngaus
          if( associated(spare % eleme(ielem) % l) ) then
             do igaus = 1,pgaus
                kgaus = kgaus + 1
                if( list_gauss(kgaus) == 0 ) then
                   spare % eleme(ielem) % l(igaus) = 0
                end if
             end do
          else
             kgaus = kgaus + pgaus
          end if 
       end do

    end if
    
    call memory_deallo(spare % memor,'LIST_GAUSS' ,'mesh_type_basic_gauss_points',list_gauss)               
    call memory_deallo(spare % memor,'NRECV_FRING','mesh_type_basic_gauss_points',nrecv_fring)               
    call memory_deallo(spare % memor,'LRECV_FRING','mesh_type_basic_gauss_points',lrecv_fring)               
    call memory_deallo(spare % memor,'LSEND_FRING','mesh_type_basic_gauss_points',lsend_fring)               
    
  end subroutine spare_mesh_gauss_points

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-05-14
  !> @brief   Distance to mesh
  !> @details Distance to mesh
  !> 
  !-----------------------------------------------------------------------
  
  subroutine spare_mesh_distance(spare,mesh,mesh_bak)

    type(typ_spare_mesh),  intent(inout) :: spare
    type(mesh_type_basic), intent(inout) :: mesh
    type(mesh_type),       intent(in)    :: mesh_bak
    integer(ip)                          :: ipoin,iboun
    real(rp)                             :: proje(3)
    type(typ_kdtree)                     :: kdtree

    if( associated(mesh % ltype) ) then
       if( mesh % ltype(1) /= POI3D ) then
          call kdtree_initialize(kdtree)
          call memory_alloca(spare % memor,'SPARE % DISTA' ,'mesh_type_basic_gauss_points',spare % dista,mesh_bak % npoin)
          if( mesh % npoin > 0 ) call kdtree_construct(mesh % nelem,mesh % npoin,mesh % lnods,mesh % ltype,mesh % coord,kdtree)
          do ipoin = 1,mesh_bak % npoin
             call kdtree_nearest_boundary(mesh_bak % coord(:,ipoin),kdtree,iboun,spare % dista(ipoin),proje)
          end do
          call kdtree_deallocate(kdtree)
       end if
    end if
    
  end subroutine spare_mesh_distance
    
end module mod_spare_mesh
!> @}
