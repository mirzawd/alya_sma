!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Domain
!> @{
!> @file    finite_volume.f90
!> @date    19/04/2016
!> @author  Guillaume Houzeaux
!> @brief   Finite volume
!> @details Compute arrays required by FV moethods
!>          1) FV_CENTER_COORD(:,:) ....... Coord of centers
!>          2) FV_CENTER_DISTANCE(:) ...... Distance between centers
!>          3) FV_CENTER_VECTOR(:) ........ Vector between centers
!>          4) FV_FACE_AREA(:) ............ Surface of faces
!>          5) FV_FACE_NORMAL(:,:) ........ Normales de las caras (order is whatever)
!>          6) FV_CENTER_FACE(:,:,:) ...... Vector from centroid to face
!>          7) FV_CELL_VOLUME(:) .......... Volume of elements
!>          8) FV_FACE_BOUNDARY(:) ........ Boundary to face correspondance
!>          9) FV_FACE_GRAPH(:,:) ......... Position of the face in the graph
!>
!> @} 
!-----------------------------------------------------------------------

module mod_finite_volume

  use def_kintyp,   only : ip,rp,lg
  use def_master,   only : INOTMASTER,kfl_paral
  use def_master,   only : nfacg
  use def_master,   only : lfacg
  use def_kermod,   only : kfl_fv_data
  use def_domain,   only : mesh_type,memor_dom
  use def_domain,   only : vmass
  use def_domain,   only : ndime
  use def_domain,   only : mnode
  use def_domain,   only : mnodb
  use def_domain,   only : fv_center_coord
  use def_domain,   only : fv_cell_volume  
  use def_domain,   only : fv_face_area    
  use def_domain,   only : fv_face_normal
  use def_domain,   only : fv_face_orientation
  use def_domain,   only : fv_center_distance
  use def_domain,   only : fv_center_vector
  use def_domain,   only : fv_center_face
  use def_domain,   only : fv_face_boundary
  use def_domain,   only : fv_face_graph
  use def_domain,   only : fv_graph_diag
  use mod_elmgeo,   only : elmgeo_element_volume
  use mod_elmgeo,   only : elmgeo_face_area
  use mod_elmgeo,   only : element_type
  use mod_memory,   only : memory_alloca
  use mod_memory,   only : memory_deallo
  use mod_messages, only : livinf
  implicit none

  private

  public :: finite_volume_arrays
  public :: finite_volume_element_to_nodes

contains

  subroutine finite_volume_arrays(meshe)

    type(mesh_type), intent(inout) :: meshe                 !< Mesh type
    integer(ip)                    :: ipoin,pnode,inodf
    integer(ip)                    :: ifacg,pelty,pflty
    integer(ip)                    :: inode,iface,jelem
    integer(ip)                    :: ielem,iz,mface,kz
    integer(ip)                    :: idime,pnodf,iboun
    integer(ip)                    :: kelem,jiz
    real(rp)                       :: xcoor(3),phi
    real(rp)                       :: elcod(ndime,mnode)
    real(rp)                       :: bocod(ndime,mnodb)
    logical(lg),     pointer       :: kfl_marked(:)

    if( kfl_fv_data == 1 ) call livinf(0_ip,'COMPUTE FINITE VOLUME ARRAYS',0_ip)
    if( kfl_fv_data == 1 .and. INOTMASTER ) then
       !
       ! Maximum nuber of faces
       !
       mface = 0
       do ielem = 1,meshe % nelem                                      
          pelty = abs(meshe % ltype(ielem))
          mface = max(mface,element_type(pelty) % number_faces)
       end do
       !
       ! Allocate memory
       !
       nullify(kfl_marked)
       call memory_alloca( memor_dom , 'FV_CENTER_COORD'     , 'finite_volume' , fv_center_coord     , ndime,meshe % nelem_2       )
       call memory_alloca( memor_dom , 'FV_VOLUME'           , 'finite_volume' , fv_cell_volume      , meshe % nelem_2             )
       call memory_alloca( memor_dom , 'FV_FACE_AREA'        , 'finite_volume' , fv_face_area        , nfacg                       )
       call memory_alloca( memor_dom , 'FV_FACE_NORMAL'      , 'finite_volume' , fv_face_normal      , ndime,nfacg                 )
       call memory_alloca( memor_dom , 'FV_FACE_ORIENTATION' , 'finite_volume' , fv_face_orientation , meshe % nzelm_2             )
       call memory_alloca( memor_dom , 'FV_CENTER_DISTANCE'  , 'finite_volume' , fv_center_distance  , meshe % nzelm_2             )
       call memory_alloca( memor_dom , 'FV_CENTER_VECTOR'    , 'finite_volume' , fv_center_vector    , ndime,meshe % nzelm_2       )
       call memory_alloca( memor_dom , 'FV_CENTER_FACE'      , 'finite_volume' , fv_center_face      , ndime,mface,meshe % nelem_2 ) 
       call memory_alloca( memor_dom , 'FV_FACE_BOUNDARY'    , 'finite_volume' , fv_face_boundary    , meshe % nboun               )
       call memory_alloca( memor_dom , 'FV_FACE_GRAPH'       , 'finite_volume' , fv_face_graph       , 4_ip,nfacg                  )
       call memory_alloca( memor_dom , 'FV_GRAPH_DIAG'       , 'finite_volume' , fv_graph_diag       , meshe % nelem               )
       !
       ! FV_CELL_VOLUME, FV_CENTER_COORD
       !
       do ielem = 1,meshe % nelem_2
          pelty = abs(meshe % ltype(ielem))
          pnode = meshe % lnnod(ielem)
          elcod(1:ndime,1:pnode) = meshe % coord(1:ndime,meshe % lnods(1:pnode,ielem))
          call elmgeo_element_volume(ndime,pelty,elcod,fv_cell_volume(ielem),fv_center_coord(1:ndime,ielem))
       end do

       do ielem = 1,meshe % nelem
          do iz = meshe % r_elm_2(ielem),meshe % r_elm_2(ielem+1)-1
             jelem = meshe % c_elm_2(iz)
             fv_center_vector(1:ndime,iz) = fv_center_coord(1:ndime,jelem) - fv_center_coord(1:ndime,ielem)
             fv_center_distance(iz)       = sqrt(dot_product(fv_center_vector(1:ndime,iz),fv_center_vector(1:ndime,iz)))
          end do
       end do
       !
       ! FV_FACE_AREA, FV_FACE_NORMAL
       !
       do ifacg = 1,nfacg
          ielem = lfacg(1,ifacg)  
          iface = lfacg(3,ifacg) 
          pelty = abs(meshe % ltype(ielem))
          pflty = element_type(pelty) % type_faces(iface) 
          do inodf = 1,element_type(pelty) % node_faces(iface) 
             inode = element_type(pelty) % list_faces(inodf,iface)
             ipoin = meshe % lnods(inode,ielem)
             bocod(1:ndime,inodf) = meshe % coord(1:ndime,ipoin)
          end do
          call elmgeo_face_area(ndime,pflty,bocod,fv_face_area(ifacg),fv_face_normal(1:ndime,ifacg))
          fv_face_normal(1:ndime,ifacg) = fv_face_normal(1:ndime,ifacg) &
               / sqrt(dot_product(fv_face_normal(1:ndime,ifacg),fv_face_normal(1:ndime,ifacg)))
       end do
       !
       ! FV_CENTER_FACE(NDIME,MFACE,NELEM_2): vector from element centroide to face centroide
       !
       do ielem = 1,meshe % nelem_2
          pelty = abs(meshe % ltype(ielem))
          do iface = 1,element_type(pelty) % number_faces
             pnodf = element_type(pelty) % node_faces(iface)
             do idime = 1,ndime
                xcoor(idime) = sum(meshe % coord(idime,meshe % lnods(element_type(pelty) % list_faces(1:pnodf,iface),ielem)))
             end do
             xcoor = xcoor / real(pnodf,rp)
             fv_center_face(1:ndime,iface,ielem) = xcoor(1:ndime) - fv_center_coord(1:ndime,ielem)
          end do
       end do
       !
       ! FV_FACE_BOUNDARY
       !
       call memory_alloca(memor_dom,'FV_FACE_GRAPH','finite_volume',kfl_marked,nfacg)       
       do iboun = 1,meshe % nboun
          ielem = meshe % lelbo(iboun)
          ifacg_loop: do ifacg = 1,nfacg
             if( lfacg(1,ifacg) == ielem .and. lfacg(2,ifacg) == 0 .and. .not. kfl_marked(ifacg) ) then 
                kfl_marked(ifacg) = .true.
                fv_face_boundary(iboun) = ifacg
                exit ifacg_loop
             end if
          end do ifacg_loop
       end do
       call memory_deallo(memor_dom,'FV_FACE_GRAPH','finite_volume',kfl_marked)
       !
       ! FV_FACE_GRAPH
       ! 
       do ifacg = 1,nfacg
          ielem = lfacg(1,ifacg)
          jelem = lfacg(2,ifacg)
          if( jelem == 0 ) then
             kz = 1
          else
             kz = 2
          end if
          iz = meshe % r_elm_2(ielem)
          do while( iz <= meshe % r_elm_2(ielem+1)-1 )
             kelem = meshe % c_elm_2(iz)
             if(      kelem == ielem ) then
                fv_face_graph(1,ifacg) = iz
                if( kz == 1 ) iz = meshe % r_elm_2(ielem+1)
             else if( kelem == jelem ) then
                fv_face_graph(2,ifacg) = iz
                iz = meshe % r_elm_2(ielem+1)
             end if
             iz = iz + 1
          end do

          if( jelem > 0 ) then
             iz = meshe % r_elm_2(jelem)
             do while( iz <= meshe % r_elm_2(jelem+1)-1 )
                kelem = meshe % c_elm_2(iz)
                if(      kelem == jelem ) then
                   fv_face_graph(3,ifacg) = iz
                else if( kelem == ielem ) then
                   fv_face_graph(4,ifacg) = iz
                end if
                iz = iz + 1
             end do
          end if

       end do
       !
       ! FV_GRAPH_DIAG
       ! 
       do ielem = 1,meshe % nelem
          iz = meshe % r_elm_2(ielem)
          do while( iz <= meshe % r_elm_2(ielem+1)-1 )
             if( meshe % c_elm_2(iz) == ielem ) then
                fv_graph_diag(ielem) = iz
                iz = meshe % r_elm_2(ielem+1)
             end if
             iz = iz + 1
          end do
       end do
       !
       ! FV_FACE_ORIENTATION
       ! 
       do ifacg = 1,nfacg
          ielem = lfacg(1,ifacg) 
          jelem = lfacg(2,ifacg) 
          iface = lfacg(3,ifacg)
          
          if( jelem > 0 ) then
             iz = fv_face_graph(2,ifacg)
          else
             iz = fv_face_graph(1,ifacg) 
          end if

          phi = dot_product(fv_center_face(1:ndime,iface,ielem),fv_face_normal(1:ndime,ifacg))
          fv_face_orientation(iz) = sign(1.0_rp,phi)
          if( jelem > 0 ) then
             jiz = fv_face_graph(4,ifacg)
             fv_face_orientation(jiz) = -fv_face_orientation(iz) 
          end if

       end do

    end if

  end subroutine finite_volume_arrays

  subroutine finite_volume_element_to_nodes(meshe,ndofn,xin,xout)

    type(mesh_type), intent(in)  :: meshe         !< Mesh type
    integer(ip),     intent(in)  :: ndofn
    real(rp),        intent(in)  :: xin(ndofn,*)
    real(rp),        intent(out) :: xout(ndofn,*)
    integer(ip)                  :: ielem,inode
    integer(ip)                  :: ipoin,pnode
    real(rp)                     :: xfact
    !
    ! Loop over elements
    !
    xout(1:ndofn,1:meshe % npoin) = 0.0_rp
    elements: do ielem = 1,meshe % nelem
       if( meshe % ltype(ielem) > 0 ) then
          pnode = meshe % lnnod(ielem)
          xfact = fv_cell_volume(ielem) / real(pnode,rp)
          do inode = 1,pnode
             ipoin                = meshe % lnods(inode,ielem)
             xout(1:ndofn,ipoin)  = xout(1:ndofn,ipoin) + xfact * xin(1:ndofn,ielem)
          end do
       end if
    end do elements
    !
    ! Parallelization
    !
    call rhsmod(ndofn,xout) 
    !
    ! Solve system
    !     
    do ipoin = 1,meshe % npoin
       xout(1:ndofn,ipoin) = xout(1:ndofn,ipoin) / vmass(ipoin)
    end do

  end subroutine finite_volume_element_to_nodes

end module mod_finite_volume
