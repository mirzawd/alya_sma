!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Domain
!> @{
!> @file    renumber_elements.f90
!> @author  Guillaume Houzeaux
!> @date    02/11/2015
!> @brief   Renumber the elements of the mesh
!> @details Renumber the elements of the mesh
!>          NEW = PERMR(OLD)
!>          OLD = INVPR(NEW)
!>          and the following variables which depend on the element
!>          numbering:
!>          LTYPE, LELCH, LNNOD, LESUB, LMATE, LNODS, 
!>          LEINV_LOC, LELBO, LESET, XFIEL, LELPO, PELPO,
!>          OMPS_DOMAINS
!> @} 
!-----------------------------------------------------------------------

module mod_renumbering

  use def_kintyp,            only : ip,rp,lg
  use def_master,            only : IMASTER
  use def_master,            only : IPARALL
  use def_master,            only : ISEQUEN
  use def_master,            only : INOTMASTER,ISLAVE
  use def_master,            only : npart,lun_outpu
  use def_master,            only : kfl_paral,ioutp
  use def_master,            only : intost
  use def_master,            only : zeror
  use def_master,            only : npoi1
  use def_master,            only : npoi2
  use def_master,            only : npoi3
  use def_domain,            only : mesh_type
  use def_domain,            only : memor_dom
  use def_domain,            only : nsteps_fiel_ondemand
  use def_domain,            only : npoin
  use def_domain,            only : npoin_own  
  use def_domain,            only : npoin_halo 
  use mod_mpio_config,       only : mpio_config
  use mod_graphs,            only : graphs_elepoi
  use mod_graphs,            only : graphs_elepoi_deallocate
  use mod_memory,            only : memory_alloca
  use mod_memory,            only : memory_deallo
  use mod_memory,            only : memory_copy
  use mod_memory,            only : memory_size
  use mod_mesh_type,         only : mesh_type_update_last_mesh
  use mod_parall,            only : PAR_COMM_MY_CODE_ARRAY
  use mod_communications,    only : PAR_ALLGATHER
  use mod_communications,    only : PAR_INTERFACE_NODE_EXCHANGE
  use mod_communications,    only : PAR_INTERFACE_OWN_NODE_EXCHANGE
  use mod_communications,    only : PAR_MAX
  use mod_communications,    only : PAR_SUM
  use mod_optional_argument, only : optional_argument
  use def_AMR,               only : interp_AMR_nelem
  use def_AMR,               only : interp_AMR_npoin
  use mod_std

#ifdef __PGI
#define MEMPGI )
#else
#define MEMPGI ,memor=memor_dom)
#endif
  implicit none
  private

  interface renumbering_lexical_order
     module procedure renumbering_lexical_order_node,&
          &           renumbering_lexical_order_type
  end interface renumbering_lexical_order

  public :: renumbering_elements
  public :: renumbering_element_arrays
  public :: renumbering_node_arrays
  public :: renumbering_nodes
  public :: renumbering_boundary_arrays
  public :: renumbering_lexical_order
  public :: renumbering_lexical_order_type
  public :: renumbering_temporal_locality
  public :: renumbering_sfc
  public :: renumbering_sfc_recursive
  public :: renumbering_reverse_cuthill_mckee
  public :: renumbering_update
   
contains

  !-----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    27/09/2016
  !> @brief   Update
  !> @details Some required updates after renumbering the mesh
  !
  !-----------------------------------------------------------------------

  subroutine renumbering_update()

    if( ISLAVE ) then
       !
       ! Useful parameters
       !
       npoin_own  = npoi3
       npoin_halo = npoin_own
       !
       ! Update main communicator
       !
       PAR_COMM_MY_CODE_ARRAY(1) % npoi1 = npoi1
       PAR_COMM_MY_CODE_ARRAY(1) % npoi2 = npoi2
       PAR_COMM_MY_CODE_ARRAY(1) % npoi3 = npoi3
       PAR_COMM_MY_CODE_ARRAY(1) % npoin = npoin

    else if( ISEQUEN ) then

       npoi3      = npoin
       npoin_own  = npoin
       npoin_halo = npoin

    else if( IMASTER ) then

       npoi3      = -1
       npoin_own  = -1
       npoin_halo = -1

    end if
    !
    ! Point to mesh structure => MESHE(NDIVI)
    !
    call mesh_type_update_last_mesh()

  end subroutine renumbering_update

  !-----------------------------------------------------------------------
  !
  !> @brief   Permutation arrays uusing element numbering
  !> @author  Guillaume Houzeaux
  !> @date    27/09/2016
  !> @details Compute the permutaion PERMR array using some 
  !>          renumering techniques
  !
  !-----------------------------------------------------------------------

  subroutine renumbering_elements(itask,meshe,permr,list_elements,offset_opt,NAME)

    integer(ip),          intent(in)             :: itask             !< Renumbering strategy
    type(mesh_type),      intent(inout)          :: meshe             !< Mesh type
    integer(ip), pointer, intent(inout)          :: permr(:)          !< Permutation
    integer(ip), pointer, intent(in),   optional :: list_elements(:)  !< List elements (old numbering)
    integer(ip),          intent(in),   optional :: offset_opt        !< Starting new element (new numbering)
    character(LEN=*),     intent(in),   optional :: NAME
    integer(ip)                                  :: ielem,ipmin
    integer(ip)                                  :: ipoin,ipmax
    integer(ip)                                  :: inode,ienew
    integer(ip)                                  :: mepoi,ielpo
    integer(ip)                                  :: kelem
    integer(ip)                                  :: jelem,izdom
    integer(ip)                                  :: jpoin
!    integer(ip)                                  :: jelpo
    integer(ip)                                  :: offset,ifjpoin(1)
    integer(ip), pointer                         :: lpoim(:)  
    integer(ip), pointer                         :: ellis(:,:) 
    integer(ip), pointer                         :: lren2(:) 
    integer(ip), pointer                         :: lnods_tmp(:,:) 
    integer(ip), pointer                         :: lnods_cpy(:,:) 

    integer(ip), pointer                         :: pelpo(:) 
    integer(ip), pointer                         :: lelpo(:)  

    logical(lg), pointer                         :: consider(:)
    logical(lg), pointer                         :: consider_ipoin(:)
    character(LEN=:), allocatable                :: my_name

    my_name = optional_argument('PERMR',NAME)
    
    if( IMASTER ) return
    !
    ! Allocate memory
    !
    !if( .not. associated(permr) ) then
    call memory_alloca(memor_dom,my_name,'renumber_elements',permr,meshe % nelem)
    !end if

    nullify(lelpo)
    nullify(pelpo)
    nullify(consider)
    nullify(consider_ipoin)

    if( itask == 0 ) then

       !----------------------------------------------------------------------
       !
       ! No renumbering
       !
       !----------------------------------------------------------------------

       return

    else if( itask == 1 ) then

       !----------------------------------------------------------------------
       !
       ! Renumber elements in the order of appearance of nodes 
       !
       !----------------------------------------------------------------------

       call memory_alloca(memor_dom,'CONSIDER','renumber_elements',consider,meshe % nelem)
       if( present(list_elements) ) then
          if( present(offset_opt) ) then
             offset = offset_opt
          else
             call runend('OFFSET SHOULD BE PRESCRIBED')
          end if
          do ielem = 1,size(list_elements,KIND=ip)
             kelem = list_elements(ielem)
             consider(kelem) = .true.
          end do
       else
          offset = 0
          consider = .true.
       end if
 
       call graphs_elepoi(meshe % npoin,meshe % nelem,meshe % mnode,meshe % lnods,meshe % lnnod,mepoi,pelpo,lelpo MEMPGI

       kelem = offset
       do ipoin = 1,meshe % npoin
          do ielpo = pelpo(ipoin),pelpo(ipoin+1)-1
             ielem = lelpo(ielpo)
             if( permr(ielem) == 0 .and. consider(ielem) ) then
                kelem        = kelem + 1
                permr(ielem) = kelem
             end if
          end do
       end do

!!$       call memory_alloca(memor_dom,'CONSIDER_IPOIN','renumber_elements',consider_ipoin,meshe % npoin)
!!$       consider_ipoin = .true.       
!!$       kelem = offset
!!$       do izdom = 1,meshe % nzdom
!!$          ipoin = meshe % c_dom(izdom)
!!$          if( consider_ipoin(ipoin) ) then
!!$             do ielpo = pelpo(ipoin),pelpo(ipoin+1)-1
!!$                ielem = lelpo(ielpo)
!!$                if( permr(ielem) == 0 .and. consider(ielem) ) then
!!$                   kelem        = kelem + 1
!!$                   permr(ielem) = kelem
!!$                end if
!!$             end do
!!$          end if
!!$          consider_ipoin(ipoin) = .false.
!!$       end do

       call graphs_elepoi_deallocate(pelpo,lelpo,memor=memor_dom)
       call memory_deallo(memor_dom,'CONSIDER','renumber_elements',consider)

    else if( itask == 2 ) then

       !----------------------------------------------------------------------
       !
       ! Renumber elements in the order of appearance of edges 
       !
       !----------------------------------------------------------------------

       call memory_alloca(memor_dom,'CONSIDER','renumber_elements',consider,meshe % nelem)
       if( present(list_elements) ) then
          if( present(offset_opt) ) then
             offset = offset_opt
          else
             call runend('OFFSET SHOULD BE PRESCRIBED')
          end if
          do ielem = 1,size(list_elements,KIND=ip)
             consider(ielem) = .true.
          end do
       else
          offset = 0
          consider = .true.
       end if

       call graphs_elepoi(meshe % npoin,meshe % nelem,meshe % mnode,meshe % lnods,meshe % lnnod,mepoi,pelpo,lelpo MEMPGI
       !
       ! Select elements with edge IPOIN-JPOIN
       !
!!$       kelem = offset
!!$       do ipoin = 1,meshe % npoin
!!$          do izdom = meshe % r_dom(ipoin),meshe % r_dom(ipoin+1)-1
!!$             jpoin = meshe % c_dom(izdom)
!!$             do ielpo = pelpo(ipoin),pelpo(ipoin+1)-1
!!$                ielem = lelpo(ielpo)
!!$                if( permr(ielem) == 0 ) then
!!$                   jelpo = pelpo(jpoin)
!!$                   do while( jelpo <= pelpo(jpoin+1)-1 )
!!$                      jelem = lelpo(jelpo)
!!$                      if( ielem == jelem .and. consider(ielem) ) then
!!$                         jelpo        = pelpo(jpoin+1)
!!$                         kelem        = kelem + 1
!!$                         permr(ielem) = kelem
!!$                      end if
!!$                      jelpo = jelpo + 1
!!$                   end do
!!$                end if
!!$             end do
!!$          end do
!!$       end do

       kelem = offset
       do ipoin = 1,meshe % npoin
          do izdom = meshe % r_dom(ipoin),meshe % r_dom(ipoin+1)-1
             jpoin = meshe % c_dom(izdom)
             do ielpo = pelpo(ipoin),pelpo(ipoin+1)-1
                ielem = lelpo(ielpo)
                ifjpoin = count(meshe % lnods(1:meshe % lnnod(ielem),ielem)==jpoin,KIND=ip)
                if( ifjpoin(1) /= 0 .and. permr(ielem) == 0 .and. consider(ielem) ) then
                   kelem        = kelem + 1
                   permr(ielem) = kelem
                end if
             end do
          end do
       end do

       call graphs_elepoi_deallocate(pelpo,lelpo,memor=memor_dom)
       call memory_deallo(memor_dom,'CONSIDER','renumber_elements',consider)

    else if( itask == 3 ) then

       !----------------------------------------------------------------------
       !
       ! Old renumbering
       !
       !----------------------------------------------------------------------

       nullify(lpoim)   
       nullify(ellis)     
       nullify(lren2)  
       nullify(lnods_tmp)  
       nullify(lnods_cpy)  

       !----------------------------------------------------------------------
       !
       ! Compute permutation array: PERMR
       !
       !----------------------------------------------------------------------
       ! 
       ! This subroutine renumber the elements
       ! It assumes that the first point in the element is the smallest
       !
       call memory_alloca(memor_dom,'LPOIM'    ,'renumber_elements',lpoim,meshe % npoin+1_ip)
       call memory_alloca(memor_dom,'LNODS_TMP','renumber_elements',lnods_tmp,meshe % mnode,meshe % nelem)
       call memory_alloca(memor_dom,'ELLIS'    ,'renumber_elements',ellis,2_ip,meshe % nelem)
       call memory_alloca(memor_dom,'LREN2'    ,'renumber_elements',lren2,meshe % nelem)
       call memory_copy(  memor_dom,'LNODS_CPY','renumber_elements',meshe % lnods,lnods_cpy,'DO_NOT_DEALLOCATE')
       !
       ! Order LNODS in increasing order
       !
       do ielem = 1,meshe % nelem
          ipmin = huge(1_ip)
          ipmax = 0
          do inode = 1,meshe % lnnod(ielem)
             ipoin = lnods_cpy(inode,ielem)
             if( ipoin > ipmax ) ipmax = ipoin
             if( ipoin < ipmin ) ipmin = ipoin
          end do
          ellis(1,ielem) = ipmin
          ellis(2,ielem) = ipmax
       end do
       !
       ! First reorder with respect to the max
       !
       do ielem = 1,meshe % nelem
          ipoin = ellis(2,ielem) + 1
          lpoim(ipoin) = lpoim(ipoin) + 1
       end do

       lpoim(1) = 1
       do ipoin = 2,meshe % npoin+1
          lpoim(ipoin) = lpoim(ipoin) + lpoim(ipoin-1)
       end do

       do ielem = 1,meshe % nelem
          ipoin        = ellis(2,ielem)
          ienew        = lpoim(ipoin)
          permr(ielem) = ienew
          lpoim(ipoin) = lpoim(ipoin) + 1
       end do

       do ielem = 1,meshe % nelem
          ienew = permr(ielem)
          do inode = 1,meshe % mnode
             lnods_tmp(inode,ienew) = lnods_cpy(inode,ielem)
          end do
       end do

       do ielem = 1,meshe % nelem
          do inode = 1,meshe % mnode
             lnods_cpy(inode,ielem) = lnods_tmp(inode,ielem)
          end do
       end do
       !
       ! Then with respect to the min
       !
       do ipoin = 1,meshe % npoin+1
          lpoim(ipoin) = 0
       end do

       do ielem = 1,meshe % nelem
          ipoin = ellis(1,ielem) + 1
          lpoim(ipoin) = lpoim(ipoin) + 1
       end do

       lpoim(1) = 1
       do ipoin = 2,meshe % npoin+1
          lpoim(ipoin) = lpoim(ipoin) + lpoim(ipoin-1)
       end do

       do ielem = 1,meshe % nelem
          ipoin        = ellis(1,ielem)
          ienew        = lpoim(ipoin)
          lren2(ielem) = ienew
          lpoim(ipoin) = lpoim(ipoin) + 1
       end do

       do ielem = 1,meshe % nelem
          ienew = lren2(ielem)
          do inode = 1,meshe % mnode
             lnods_tmp(inode,ienew) = lnods_cpy(inode,ielem) 
          end do
       end do

       do ielem = 1,meshe % nelem 
          do inode = 1,meshe % mnode
             lnods_cpy(inode,ielem) = lnods_tmp(inode,ielem) 
          end do
          jelem = permr(ielem)        ! JELEM (new) = PERMR(IELEM)
          permr(ielem) = lren2(jelem)
       end do

       call memory_deallo(memor_dom,'LREN2'    ,'renumber_elements',lren2)
       call memory_deallo(memor_dom,'ELLIS'    ,'renumber_elements',ellis)
       call memory_deallo(memor_dom,'LNODS_TMP','renumber_elements',lnods_tmp)
       call memory_deallo(memor_dom,'LNODS_CPY','renumber_elements',lnods_cpy)
       call memory_deallo(memor_dom,'LPOIM'    ,'renumber_elements',lpoim)

    end if

    if( allocated(my_name) ) deallocate(my_name)
    
  end subroutine renumbering_elements

  !-----------------------------------------------------------------------
  !
  !> @brief   Permutation arrays using element numbering
  !> @author  Guillaume Houzeaux
  !> @date    27/09/2016
  !> @details Renu,ber the following element arrays:
  !>          LTYPE, LELCH, LNNOD, LESUB, LMATE, LNODS, 
  !>          LEINV_LOC, LELBO, LESET, XFIEL, LELPO, PELPO,
  !>          OMPS_DOMAINS
  !>
  !>          NEW = PERMR(OLD)
  !
  !-----------------------------------------------------------------------

  subroutine renumbering_element_arrays(permr,NAME)

    use def_master,         only : leinv_loc
    use def_master,         only : NELEM_TYPE
    use def_domain,         only : nelem,ltype
    use def_domain,         only : lelch,lnnod
    use def_domain,         only : lesub,lmate
    use def_domain,         only : lelbo,leset
    use def_domain,         only : nboun,nfiel
    use def_domain,         only : lemsh,lelev
    use def_domain,         only : lgaus
    use def_domain,         only : kfl_field
    use def_domain,         only : neset,lnods
    use def_domain,         only : npoin,lelpo
    use def_domain,         only : pelpo,xfiel
    use def_domain,         only : ompss_domains
    use mod_memory,         only : memory_copy
    use mod_memory,         only : memory_deallo
    use mod_memory,         only : memory_renumber
    use mod_redistribution, only : commd_nelem

    integer(ip), pointer, intent(inout) :: permr(:)             !< Inverse permutation
    character(LEN=*),     intent(in),   optional :: NAME

    integer(ip)                         :: iboun,kelem
    integer(ip)                         :: ipoin,ielem
    integer(ip)                         :: ielpo,isubd,ifiel,istep
    integer(ip)                         :: ielem_old,ielem_new
    integer(ip)                         :: nsteps,ii               
    integer(ip), pointer                :: ltype_tmp(:)         ! Element type
    integer(ip), pointer                :: lelch_tmp(:)         ! Element characteristic
    integer(ip), pointer                :: lnnod_tmp(:)         ! Element number of nodes
    integer(ip), pointer                :: lnods_tmp(:,:)       ! Element connectivity
    integer(ip), pointer                :: lesub_tmp(:)         ! Element subdomain
    integer(ip), pointer                :: lmate_tmp(:)         ! Element material
    integer(ip), pointer                :: leinv_tmp(:)         ! Element original numbering
    integer(ip), pointer                :: leset_tmp(:)         ! Element set
    integer(ip), pointer                :: lemsh_tmp(:)         ! Postprocess on original mesh
    integer(ip), pointer                :: lelev_tmp(:)         ! Postprocess on original mesh
    integer(ip), pointer                :: lgaus_tmp(:)         ! Postprocess on original mesh

    real(rp),    pointer                :: xfiel_tmp(:,:)
    character(LEN=:), allocatable       :: my_name

    if( IMASTER ) return
    
    my_name = optional_argument('PERMR',NAME)
    !
    ! Check permutation array
    !
    if( .not. associated(permr) ) call runend('RENUMBERING_ELEMENT_ARRAYS: PERMR NOT ASSOCIATED')
    if( minval(permr) < 1     )   call runend('RENUMBERING_ELEMENT_ARRAYS: PERMR NOT CORRECT 1')
    if( maxval(permr) > nelem )   call runend('RENUMBERING_ELEMENT_ARRAYS: PERMR NOT CORRECT 2')

    nullify(ltype_tmp)
    nullify(lelch_tmp)
    nullify(lnnod_tmp)
    nullify(lnods_tmp)
    nullify(lesub_tmp)
    nullify(lmate_tmp)
    nullify(leinv_tmp)
    nullify(leset_tmp)
    nullify(lemsh_tmp)
    nullify(lelev_tmp)
    nullify(lgaus_tmp)

    !----------------------------------------------------------------------
    !
    ! LTYPE, LNNOD, LNODS, LELCH, LESUB, LMATE and LEINV_LOC
    !
    !----------------------------------------------------------------------

    call memory_copy(memor_dom,'LTYPE_TMP','renumbering_element_arrays',ltype    ,ltype_tmp,'DO_NOT_DEALLOCATE',COPY_NAME='LTYPE_TMP')
    call memory_copy(memor_dom,'LELCH_TMP','renumbering_element_arrays',lelch    ,lelch_tmp,'DO_NOT_DEALLOCATE',COPY_NAME='LELCH_TMP')
    call memory_copy(memor_dom,'LNNOD_TMP','renumbering_element_arrays',lnnod    ,lnnod_tmp,'DO_NOT_DEALLOCATE',COPY_NAME='LNNOD_TMP')
    call memory_copy(memor_dom,'LNODS_TMP','renumbering_element_arrays',lnods    ,lnods_tmp,'DO_NOT_DEALLOCATE',COPY_NAME='LNODS_TMP')
    call memory_copy(memor_dom,'LESUB_TMP','renumbering_element_arrays',lesub    ,lesub_tmp,'DO_NOT_DEALLOCATE',COPY_NAME='LESUB_TMP')
    call memory_copy(memor_dom,'LMATE_TMP','renumbering_element_arrays',lmate    ,lmate_tmp,'DO_NOT_DEALLOCATE',COPY_NAME='LMATE_TMP')
    call memory_copy(memor_dom,'LEINV_TMP','renumbering_element_arrays',leinv_loc,leinv_tmp,'DO_NOT_DEALLOCATE',COPY_NAME='LEINV_TMP')

    do ielem_old = 1,nelem
       ielem_new            = permr(ielem_old)       ! IELEM_NEW: new numbering
       ltype(ielem_new)     = ltype_tmp(ielem_old)
       lelch(ielem_new)     = lelch_tmp(ielem_old)
       lnods(:,ielem_new)   = lnods_tmp(:,ielem_old)
       lesub(ielem_new)     = lesub_tmp(ielem_old)
       lmate(ielem_new)     = lmate_tmp(ielem_old)
       leinv_loc(ielem_new) = leinv_tmp(ielem_old)
    end do

    if( associated(lnnod) ) then
       do ielem_old = 1,nelem
          ielem_new            = permr(ielem_old)       ! IELEM_NEW: new numbering
          lnnod(ielem_new)     = lnnod_tmp(ielem_old)
       end do
    end if
    
    call memory_deallo(memor_dom,'LEINV_TMP','renumbering_element_arrays',leinv_tmp)
    call memory_deallo(memor_dom,'LMATE_TMP','renumbering_element_arrays',lmate_tmp)
    call memory_deallo(memor_dom,'LESUB_TMP','renumbering_element_arrays',lesub_tmp)
    call memory_deallo(memor_dom,'LNODS_TMP','renumbering_element_arrays',lnods_tmp)
    call memory_deallo(memor_dom,'LNNOD_TMP','renumbering_element_arrays',lnnod_tmp)
    call memory_deallo(memor_dom,'LELCH_TMP','renumbering_element_arrays',lelch_tmp)
    call memory_deallo(memor_dom,'LTYPE_TMP','renumbering_element_arrays',ltype_tmp)

    !----------------------------------------------------------------------
    !
    ! LGAUS
    !
    !----------------------------------------------------------------------

    if( associated(lgaus) ) then

       call memory_copy(memor_dom,'LGAUS_TMP','renumbering_element_arrays',lgaus,lgaus_tmp,'DO_NOT_DEALLOCATE')
       do ielem_old = 1,nelem
          ielem_new        = permr(ielem_old)
          lgaus(ielem_new) = lgaus_tmp(ielem_old)
       end do
       call memory_deallo(memor_dom,'LGAUS_TMP','renumbering_element_arrays',lgaus_tmp)

    end if

    !----------------------------------------------------------------------
    !
    ! LELBO
    !
    !----------------------------------------------------------------------

    if( associated(lelbo) ) then

       do iboun = 1,nboun
          ielem_old = lelbo(iboun)
          if( ielem_old /= 0 ) then
             lelbo(iboun) = permr(ielem_old)
          end if
       end do

    end if

    !--------------------------------------------------------------------
    !
    ! Mesh level: LELEV 
    !
    !--------------------------------------------------------------------

    if( associated(lelev) ) then

       call memory_copy(memor_dom,'LELEV_TMP','renumbering_element_arrays',lelev,lelev_tmp,'DO_NOT_DEALLOCATE')
       do ielem_old = 1,nelem
          ielem_new        = permr(ielem_old)
          lelev(ielem_new) = lelev_tmp(ielem_old)
       end do
       call memory_deallo(memor_dom,'LELEV_TMP','renumbering_element_arrays',lelev_tmp)

    end if

    !--------------------------------------------------------------------
    !
    ! Postprocess on original mesh: LEMSH 
    !
    !--------------------------------------------------------------------

    if( associated(lemsh) ) then

       do ielem_new = 1,memory_size(lemsh)
          ielem_old = lemsh(ielem_new)
          if( ielem_old /= 0 ) lemsh(ielem_new) = permr(ielem_old)
       end do

    end if

    !----------------------------------------------------------------------
    !
    ! Sets: LESET
    !
    !----------------------------------------------------------------------

    if( neset > 0 .and. associated(leset) ) then

       call memory_copy(memor_dom,'LESET_TMP','renumbering_element_arrays',leset,leset_tmp,'DO_NOT_DEALLOCATE')
       do ielem_old = 1,nelem
          ielem_new        = permr(ielem_old)
          leset(ielem_new) = leset_tmp(ielem_old)
       end do
       call memory_deallo(memor_dom,'LESET_TMP','renumbering_element_arrays',leset_tmp)

    end if

    !----------------------------------------------------------------------
    !
    ! Fields: XFIEL
    !
    !----------------------------------------------------------------------

    do ifiel = 1,nfiel
       if( kfl_field(2,ifiel) == NELEM_TYPE ) then
          if ( (kfl_field(6,ifiel) /= 1) .OR. (mpio_config%output%post_process%export_only) ) then 
             nsteps = kfl_field(4,ifiel)
          else
             nsteps = nsteps_fiel_ondemand
          end if
          do istep = 1,nsteps
             nullify(xfiel_tmp)
             call memory_alloca(memor_dom,'XFIEL_TMP','renumbering_element_arrays',xfiel_tmp,kfl_field(1,ifiel),nelem)
             do ielem = 1,nelem
                xfiel_tmp(:,ielem) = xfiel(ifiel) % a(:,ielem,istep)
             end do
             call memory_renumber(memor_dom,'XFIEL','renumbering_element_arrays',xfiel_tmp,permr)
             do ielem = 1,nelem
               xfiel(ifiel) % a(:,ielem,istep) = xfiel_tmp(:,ielem)
             end do
             call memory_deallo(memor_dom,'XFIEL_TMP','renumbering_element_arrays',xfiel_tmp) 
          end do
       end if
    end do

    !----------------------------------------------------------------------
    !
    ! Node to element graph: PELPO, LELPO
    !
    !----------------------------------------------------------------------

    if( associated(pelpo) .and. associated(lelpo) ) then
       do ipoin = 1,npoin
          do ielpo = pelpo(ipoin),pelpo(ipoin+1)-1
             ielem_old    = lelpo(ielpo)
             ielem_new    = permr(ielem_old)
             lelpo(ielpo) = ielem_new
          end do
       end do
    end if

    !----------------------------------------------------------------------
    !
    ! Partitions have been subdivided into sub-subdomains: OMPSS_DOMAINS
    !
    !----------------------------------------------------------------------

    if( associated(ompss_domains) ) then 
       do isubd = 1,size(ompss_domains,KIND=ip)
          do kelem = 1,size(ompss_domains(isubd) % elements,KIND=ip)
             ielem_old = ompss_domains(isubd) % elements(kelem)
             ielem_new = permr(ielem_old)
             ompss_domains(isubd) % elements(kelem) = ielem_new
          end do
       end do
    end if

    !----------------------------------------------------------------------
    !
    ! Repartitioning communicator
    !
    !----------------------------------------------------------------------

    if( associated(commd_nelem % lrecv_perm) ) then
       do ii = 1,commd_nelem % lrecv_dim         
          commd_nelem % lrecv_perm(ii) = permR(commd_nelem % lrecv_perm(ii))
       end do
    end if

    !----------------------------------------------------------------------
    !
    ! AMR communicator
    !
    !----------------------------------------------------------------------

    if( associated(interp_AMR_nelem) ) then       
       call interp_AMR_nelem % permutation(permr)
    end if
    
    !----------------------------------------------------------------------
    !
    ! Deallocate memory
    !
    !----------------------------------------------------------------------

    call memory_deallo(memor_dom,my_name,'renumbering_element_arrays',permr)
    if( allocated(my_name) ) deallocate(my_name)

  end subroutine renumbering_element_arrays

  !-----------------------------------------------------------------------
  !>
  !> @date    07/11/2013
  !> @author  Guillaume Houzeaux
  !> @brief   Renumber geometrical arrays
  !> @details Renumber geometrical arrays, read from files,
  !>          using a node permutation array PERMR
  !>          \verbatim
  !>            LNODS ........
  !>            LNODB ........
  !>            COORD ........
  !>            LNSEC ........
  !>            LNOCH ........
  !>            LMAST ........
  !>            LPERI ........
  !>            XFIEL ........
  !>            LGROU_DOM ....
  !>            KFL_CODNO ....
  !>            BOUND_PERM ...
  !>            LNINV_LOC ....
  !>            BVCOD ........
  !>          \endverbatim
  !>
  !>          The bounds can be specified if for example, arrays are to be
  !>          renumbered including the halo nodes.
  !>
  !>          Special care on arrays that are not defiend on halos.
  !>
  !-----------------------------------------------------------------------

  subroutine renumbering_node_arrays(permr,npoin_opt,nelem_opt,nboun_opt,PERMUTE_LMAST)

    use def_kintyp
    use def_parame
    use def_domain
    use def_master
    use def_kermod
    use mod_parall,         only : commd
    use mod_memory,         only : memory_alloca
    use mod_memory,         only : memory_deallo
    use mod_elmgeo,         only : element_type,elmgeo_element_type_initialization
    use mod_redistribution, only : commd_npoin
    implicit none

    integer(ip),      pointer,  intent(inout) :: permr(:)             !< Inverse permutation
    integer(ip),      optional, intent(in)    :: npoin_opt            !< Optional # nodes
    integer(ip),      optional, intent(in)    :: nelem_opt            !< Optional # elements
    integer(ip),      optional, intent(in)    :: nboun_opt            !< Optional # boundaries
    logical(lg),      optional, intent(in)    :: PERMUTE_LMAST        !< If periodicity should be permuted

    integer(ip)                               :: ipoin,ielem,inode,iboun,idime,jpoin
    integer(ip)                               :: icono,inset,inodb,pblty,ifiel
    integer(ip)                               :: nn,ii,kpoin,npoin_min,iperi,istep
    integer(ip)                               :: ipoin_new,jpoin_new,pelty
    integer(ip)                               :: npoin_end
    integer(ip)                               :: nelem_end
    integer(ip)                               :: nboun_end
    integer(ip)                               :: nsteps
    logical(lg)                               :: if_permute_lmast
    integer(ip),      pointer                 :: kfl_codno_tmp(:,:) 
    integer(ip),      pointer                 :: lnoch_tmp(:) 
    integer(ip),      pointer                 :: lllll_tmp(:) 
    integer(ip),      pointer                 :: lgrou_dom_tmp(:)   
    integer(ip),      pointer                 :: lninv_loc_tmp(:)   
    integer(ip),      pointer                 :: bound_perm(:) 
    integer(ip),      pointer                 :: bound_invp(:) 
    real(rp),         pointer                 :: coord_tmp(:,:)     
    type(r3p),        pointer                 :: xfiel_tmp(:)       
    type(r2p),        pointer                 :: bvcod_tmp(:)       
    type(i1pp),       pointer                 :: linno_tmp(:)       

    nullify(kfl_codno_tmp)
    nullify(lnoch_tmp)
    nullify(lllll_tmp)
    nullify(lgrou_dom_tmp)
    nullify(lninv_loc_tmp)
    nullify(bound_perm)
    nullify(bound_invp)
    nullify(coord_tmp)
    nullify(xfiel_tmp)
    nullify(bvcod_tmp)
    nullify(linno_tmp)

    if( present(npoin_opt) ) then
       npoin_end = npoin_opt
    else
       npoin_end = npoin
    end if
    if( present(nelem_opt) ) then
       nelem_end = nelem_opt
    else
       nelem_end = nelem
    end if
    if( present(nboun_opt) ) then
       nboun_end = nboun_opt
    else
       nboun_end = nboun
    end if

    if( present(PERMUTE_LMAST) ) then
       if_permute_lmast = PERMUTE_LMAST
    else
       if_permute_lmast = .true.       
    end if

    !----------------------------------------------------------------------
    !
    ! LNODS, LNODB, COORD: geometrical arrays
    !
    !----------------------------------------------------------------------

    if( associated(lnods) ) then

       if( associated(lnnod) ) then
          do ielem = 1,nelem_end
             do inode = 1,lnnod(ielem)
                lnods(inode,ielem) = permr(lnods(inode,ielem))
             end do
          end do
       else
          do ielem = 1,nelem_end
             pelty = abs(ltype(ielem))
             do inode = 1,element_type(pelty) % number_nodes
                lnods(inode,ielem) = permr(lnods(inode,ielem))
             end do
          end do
       end if

    end if

    if( associated(lnodb) ) then

       do iboun = 1,nboun_end
          pblty = ltypb(iboun)
          do inodb = 1,nnode(pblty)
             lnodb(inodb,iboun) = permr(lnodb(inodb,iboun))
          end do
       end do

    end if

    if( associated(coord) ) then

       call memory_alloca(memor_dom,'COORD_TMP','mod_renumbering',coord_tmp,ndime,npoin_end)
       do ipoin = 1,npoin_end
          coord_tmp(1:ndime,ipoin) = coord(1:ndime,ipoin)
       end do
       do ipoin = 1,npoin_end
          jpoin = permr(ipoin)
          coord(1:ndime,jpoin) = coord_tmp(1:ndime,ipoin)
       end do
       call memory_deallo(memor_dom,'COORD_TMP','mod_renumbering',coord_tmp)

    end if
    
    !--------------------------------------------------------------------
    !
    ! Sets: LNSET
    !
    !--------------------------------------------------------------------

    if( associated(lnset) ) then

       npoin_min = min(memory_size(lnset),npoin_end)
       call memory_alloca(memor_dom,'LLLLL_TMP','mod_renumbering',lllll_tmp,npoin_min)
       lllll_tmp(1:npoin_min) = lnset(1:npoin_min)
       do ipoin = 1,npoin_min
          jpoin = permr(ipoin)
          lnset(jpoin) = lllll_tmp(ipoin)
       end do
       call memory_deallo(memor_dom,'LLLLL_TMP','mod_renumbering',lllll_tmp)

    end if

    !--------------------------------------------------------------------
    !
    ! Sets: LNSEC
    !
    !--------------------------------------------------------------------

    if( associated(lnsec) ) then

       do inset = 1,nnset
          ipoin = lnsec(inset)
          lnsec(inset) = permr(ipoin)
       end do

    end if

    !--------------------------------------------------------------------
    !
    ! LNOCH
    !
    !--------------------------------------------------------------------

    if( associated(lnoch) ) then

       npoin_min = min(memory_size(lnoch),npoin_end)
       call memory_alloca(memor_dom,'LNOCH_TMP','mod_renumbering',lnoch_tmp,npoin_min)
       lnoch_tmp(1:npoin_min) = lnoch(1:npoin_min)
       do ipoin = 1,npoin_min
          jpoin = permr(ipoin)
          lnoch(jpoin) = lnoch_tmp(ipoin)
       end do
       call memory_deallo(memor_dom,'LNOCH_TMP','mod_renumbering',lnoch_tmp)

    end if

    !--------------------------------------------------------------------
    !
    ! LPOTY
    !
    !--------------------------------------------------------------------

    if( associated(lpoty) ) then

       npoin_min = min(memory_size(lpoty),npoin_end)
       call memory_alloca(memor_dom,'LLLLL_TMP','mod_renumbering',lllll_tmp,npoin_min)
       lllll_tmp(1:npoin_min) = lpoty(1:npoin_min)
       do ipoin = 1,npoin_min
          jpoin = permr(ipoin)
          lpoty(jpoin) = lllll_tmp(ipoin)
       end do
       call memory_deallo(memor_dom,'LLLLL_TMP','mod_renumbering',lllll_tmp)

    end if

    !--------------------------------------------------------------------
    !
    ! LNLEV
    !
    !--------------------------------------------------------------------

    if( associated(lnlev) ) then

       npoin_min = min(memory_size(lnlev),npoin_end)
       call memory_alloca(memor_dom,'LLLLL_TMP','mod_renumbering',lllll_tmp,npoin_min)
       lllll_tmp(1:npoin_min) = lnlev(1:npoin_min)
       do ipoin = 1,npoin_min
          jpoin = permr(ipoin)
          lnlev(jpoin) = lllll_tmp(ipoin)
       end do
       call memory_deallo(memor_dom,'LLLLL_TMP','mod_renumbering',lllll_tmp)

    end if

    !--------------------------------------------------------------------
    !
    ! LPMSH: postprocess on original mesh
    !
    !--------------------------------------------------------------------

    if( associated(lpmsh) ) then

       do ipoin = 1,memory_size(lpmsh)
          jpoin = lpmsh(ipoin)
          lpmsh(ipoin) = permr(jpoin)
       end do

    end if

    !--------------------------------------------------------------------
    !
    ! Fields
    !
    !--------------------------------------------------------------------

    if( nfiel > 0 ) then

       call memory_alloca(memor_dom,'XFIEL_TMP','mod_renumbering',xfiel_tmp,nfiel)
       do ifiel = 1,nfiel
          if( kfl_field(1,ifiel) > 0 ) then

             if ( (kfl_field(6,ifiel) /= 1) .OR. (mpio_config%output%post_process%export_only) ) then 
                nsteps = kfl_field(4,ifiel)
             else
                nsteps = nsteps_fiel_ondemand
             end if


             if( kfl_field(2,ifiel) == NPOIN_TYPE ) then
                npoin_min = min( size(xfiel(ifiel) % a,2,KIND=ip),npoin_end)                
                call memory_alloca(memor_dom,'XFIEL_TMP % A','mod_renumbering',xfiel_tmp(ifiel)%a,kfl_field(1,ifiel),npoin_min,kfl_field(4,ifiel))
                do istep = 1,nsteps
                   do ipoin = 1,npoin_min
                      do idime = 1,kfl_field(1,ifiel)
                         xfiel_tmp(ifiel) % a(idime,ipoin,istep) = xfiel(ifiel) % a(idime,ipoin,istep)
                      end do
                   end do
                end do !istep
                do istep = 1,nsteps
                   do ipoin = 1,npoin_min
                      jpoin = permr(ipoin)
                      do idime = 1,kfl_field(1,ifiel)
                         xfiel(ifiel) % a(idime,jpoin,istep) = xfiel_tmp(ifiel) % a(idime,ipoin,istep)
                      end do
                   end do
                end do !istep
             end if ! kfl_field


          end if
       end do
       call memory_deallo(memor_dom,'XFIEL_TMP','mod_renumbering',xfiel_tmp)

    end if

    !--------------------------------------------------------------------
    !
    ! Groups
    !
    !--------------------------------------------------------------------

    if( ngrou_dom > 0 .and. associated(lgrou_dom) ) then

       npoin_min = min(memory_size(lgrou_dom),npoin_end)

       call memory_alloca(memor_dom,'LGROU_DOM_TMP','mod_renumbering',lgrou_dom_tmp,npoin_min)       
       do ipoin = 1,npoin_min
          lgrou_dom_tmp(ipoin) = lgrou_dom(ipoin)
       end do
       do ipoin = 1,npoin_min
          jpoin = permr(ipoin)
          lgrou_dom(jpoin) = lgrou_dom_tmp(ipoin)
       end do
       call memory_deallo(memor_dom,'LGROU_DOM_TMP','mod_renumbering',lgrou_dom_tmp)

    end if

    !--------------------------------------------------------------------
    !
    ! Boundary conditions: KFL_CODNO
    !
    !--------------------------------------------------------------------

    if( kfl_icodn > 0 .and. associated(kfl_codno) ) then

       npoin_min = min(memory_size(kfl_codno,2_ip),npoin_end)
       call memory_alloca(memor_dom,'KFL_CODNO_TMP','mod_renumbering',kfl_codno_tmp,mcono,npoin_min)
       do ipoin = 1,npoin_min
          do icono = 1,mcono
             kfl_codno_tmp(icono,ipoin) = kfl_codno(icono,ipoin)
          end do
       end do
       do ipoin = 1,npoin_min
          jpoin = permr(ipoin)
          do icono = 1,mcono
             kfl_codno(icono,jpoin) = kfl_codno_tmp(icono,ipoin)
          end do
       end do

       call memory_deallo(memor_dom,'KFL_CODNO_TMP','mod_renumbering',kfl_codno_tmp)

    end if

    !--------------------------------------------------------------------
    !
    ! Periodicity: LPERI
    !
    !--------------------------------------------------------------------

    if( nperi > 0 .and. associated(lperi) ) then

       do iperi = 1,nperi
          ipoin          = lperi(1,iperi)
          jpoin          = lperi(2,iperi)
          lperi(1,iperi) = permr(ipoin)
          lperi(2,iperi) = permr(jpoin)
       end do

    end if

    ! PERDIOCITY TODAY
    if( new_periodicity == 1 ) then
       if( associated(lmast) ) then
          npoin_min = min(memory_size(lmast),npoin_end)
          call memory_alloca(memor_dom,'LLLLL_TMP','mod_renumbering',lllll_tmp,npoin_min)
          lllll_tmp(1:npoin_min) = lmast(1:npoin_min)
          lmast = 0
          do ipoin = 1,npoin
             jpoin = lllll_tmp(ipoin)
             if( abs(jpoin) > 0 ) then
                ipoin_new = permr(ipoin)
                lmast(ipoin_new) = jpoin
             end if
          end do
          call memory_deallo(memor_dom,'LLLLL_TMP','mod_renumbering',lllll_tmp)
       end if
    else if( new_periodicity == 0 ) then
       if( associated(lmast) ) then
          npoin_min = min(memory_size(lmast),npoin_end)
          call memory_alloca(memor_dom,'LLLLL_TMP','mod_renumbering',lllll_tmp,npoin_min)
          lllll_tmp(1:npoin_min) = lmast(1:npoin_min)
          lmast = 0
          do ipoin = 1,npoin
             jpoin = lllll_tmp(ipoin)
             if( jpoin > 0 ) then
                ipoin_new = permr(ipoin)
                if( if_permute_lmast ) then
                   jpoin_new = permr(jpoin)
                   lmast(ipoin_new) = jpoin_new
                else
                   lmast(ipoin_new) = jpoin
                end if
             end if
          end do
          call memory_deallo(memor_dom,'LLLLL_TMP','mod_renumbering',lllll_tmp)
       end if
    end if

    !--------------------------------------------------------------------
    !
    ! BOUND_PERM: Parallel communication arrays
    !
    !--------------------------------------------------------------------

    if( associated(commd) .and. ISLAVE ) then
       if( associated(commd % bound_perm) ) then

          call memory_alloca(memor_dom,'BOUND_PERM_1','mod_renumbering',bound_perm,commd % bound_dim)
          do ii = 1,commd % bound_dim         
             bound_perm(ii) = commd % bound_perm(ii)
          end do
          do ii = 1,commd % bound_dim
             commd % bound_perm(ii) = permR( bound_perm(ii) )
          end do
          call memory_deallo(memor_dom,'BOUND_PERM_1','mod_renumbering',bound_perm)

       end if
       if( associated(commd % bound_invp) ) then

          call memory_alloca(memor_dom,'BOUND_INVP_1','mod_renumbering',bound_invp,commd % bound_dim)
          do ii = 1,commd % bound_dim         
             bound_invp(ii) = commd % bound_invp(ii)
          end do
          do ii = 1,commd % bound_dim 
             commd % bound_invp(ii) = permR( bound_invp(ii) )
          end do
          call memory_deallo(memor_dom,'BOUND_INVP_1','mod_renumbering',bound_invp)

       end if
    end if

    !--------------------------------------------------------------------
    !
    ! Redistribution communication arrays
    !
    !--------------------------------------------------------------------

    if( associated(commd_npoin % lrecv_perm) ) then
       do ii = 1,commd_npoin % lrecv_dim         
          commd_npoin % lrecv_perm(ii) = permR(commd_npoin % lrecv_perm(ii))
       end do
    end if
    
    !--------------------------------------------------------------------
    !
    ! AMR
    !
    !--------------------------------------------------------------------
    
    if( associated(interp_AMR_npoin) ) then
       call interp_AMR_npoin % permutation(permr)
    end if

    !--------------------------------------------------------------------
    !
    ! LNINV_LOC: global numbering
    !
    !--------------------------------------------------------------------

    if( associated(lninv_loc) ) then

       npoin_min = min(memory_size(lninv_loc),npoin_end)
       call memory_alloca(memor_dom,'LNINV_LOC_TMP','mod_renumbering',lninv_loc_tmp,npoin_min)
       lninv_loc_tmp(1:npoin_min) = lninv_loc(1:npoin_min)
       lninv_loc = 0
       do ipoin = 1,npoin_min
          kpoin = permR(ipoin) 
          lninv_loc(kpoin) = lninv_loc_tmp(ipoin)
       end do
       call memory_deallo(memor_dom,'LNINV_LOC_TMP','mod_renumbering',lninv_loc_tmp)

    end if

    !--------------------------------------------------------------------
    !
    ! Halo communication
    !
    !--------------------------------------------------------------------

    if( associated(commd) .and. IPARALL ) then
       if( associated(commd % ghost_recv_node_perm) ) then

          do ii = 1,commd % ghost_recv_node_dim
             ipoin = commd % ghost_recv_node_perm(ii)
             if( ipoin <= size(permr) ) &
                  commd % ghost_recv_node_perm(ii) = permr(ipoin)
          end do

       end if

       if( associated(commd % ghost_send_node_perm) ) then

          do ii = 1,commd % ghost_send_node_dim
             ipoin = commd % ghost_send_node_perm(ii)
             if( ipoin <= size(permr) ) &
                  commd % ghost_send_node_perm(ii) = permr(ipoin)
          end do
       end if
    end if

    !--------------------------------------------------------------------
    !
    ! Interpolation
    !
    !--------------------------------------------------------------------

    if( 1 == 2 ) then

       allocate(linno_tmp(npoin_end))
       do ipoin = 1,npoin_end
          nn    = meshe(ndivi) % linno(ipoin) % n
          linno_tmp(ipoin) % n = nn
          allocate( linno_tmp(ipoin) % l(nn))
          do ii = 1,nn
             linno_tmp(ipoin) % l(ii) = meshe(ndivi) % linno(ipoin) % l(ii)
          end do
       end do
       do ipoin = 1,npoin_end
          jpoin = permr(ipoin)
          nn    = linno_tmp(ipoin) % n
          meshe(ndivi) % linno(jpoin) % n = nn
          allocate( linno_tmp(ipoin) % l(nn))
          do ii = 1,nn
             meshe(ndivi) % linno(jpoin) % l(ii) = linno_tmp(ipoin) % l(ii)
          end do
          deallocate(linno_tmp(ipoin) % l)
       end do
       deallocate(linno_tmp)

    end if

  end subroutine renumbering_node_arrays

  !-----------------------------------------------------------------------
  !>
  !> @date    07/11/2013
  !> @author  Guillaume Houzeaux
  !> @brief   Renumber nodes
  !> @details Renumber nodes
  !>
  !-----------------------------------------------------------------------

  subroutine renumbering_nodes(what,kpoin,ia,ja,permr)

    character(*),             intent(in)    :: what
    integer(ip)                             :: kpoin             !< Number of nodes
    integer(ip),    pointer,  intent(in)    :: ia(:)             !< Node graph
    integer(ip),    pointer,  intent(in )   :: ja(:)             !< Node graph
    integer(ip),    pointer,  intent(inout) :: permr(:)          !< Permutation

    call runend('GENERIC SUBROUTINE FOR NODE RENUMBERING UNDER CONSTRUCTION')

  end subroutine renumbering_nodes

  !-----------------------------------------------------------------------
  !>
  !> @date    07/11/2013
  !> @author  Guillaume Houzeaux
  !> @brief   Renumber nodes in lexical order
  !> @details Renumber nodes in global lexical order using mesh type
  !>
  !-----------------------------------------------------------------------

  subroutine renumbering_lexical_order_deallocate(permr,what)

    integer(ip),     pointer, intent(inout) :: permr(:)          !< Permutation
    character(*),    optional               :: what              !< Should we deallocate

    if( INOTMASTER ) call memory_deallo(memor_dom,'permr','renumbering_lexical_order',permr)

  end subroutine renumbering_lexical_order_deallocate

  !-----------------------------------------------------------------------
  !>
  !> @date    07/11/2013
  !> @author  Guillaume Houzeaux
  !> @brief   Renumber nodes in lexical order
  !> @details Renumber nodes in global lexical order using mesh type
  !>          \verbatim
  !>          o Own nodes
  !>          b Boundary nodes
  !>          h Halo nodes (not in my subdomain)
  !>
  !>          o----o----o----o----o----b----b----b----h----h----h
  !>          <------------------->                               NPOIN_OWN
  !>          <---------------------------------->                NPOIN
  !>          <-------------------------------------------------> NPOIN_HALO
  !>          \endverbatim
  !>
  !-----------------------------------------------------------------------

  subroutine renumbering_lexical_order_node(npoin,npoin_halo,npoin_own,permr,what,ndofn_opt)

    integer(ip),              intent(in)    :: npoin             !< Nodes
    integer(ip),              intent(in)    :: npoin_own         !< Nodes
    integer(ip),              intent(in)    :: npoin_halo        !< Nodes
    integer(ip),     pointer, intent(inout) :: permr(:)          !< Permutation
    character(*),    optional               :: what              !< Should we deallocate
    integer(ip),     optional               :: ndofn_opt         !< If dof > 1
    integer(ip),     pointer                :: npoin_own_tot(:)
    integer(ip)                             :: ipoin,ndofn
    integer(ip)                             :: npoin_offset
    logical(lg)                             :: with_halo
    !
    ! What to do
    !
    if( present(what) ) then
       if( trim(what) == 'DEALLOCATE' ) then
          if( INOTMASTER ) call memory_deallo(memor_dom,'permr','renumbering_lexical_order',permr)
          return
       end if
    end if
    with_halo = .true.
    if( present(what) ) then
       if( trim(what) == 'WITHOUT HALOS' ) then
          with_halo = .false.
       end if
    end if
    !
    ! NDOFN > 1
    !
    if( present(ndofn_opt) ) then
       ndofn = ndofn_opt
    else
       ndofn = 1
    end if
    !
    ! Permute
    !
    nullify(npoin_own_tot)

    if( with_halo ) then
       !
       ! Renumber nodes including halos (default)
       !
       if( .not. associated(permr) .and. INOTMASTER ) then
          call memory_alloca(memor_dom,'permr','renumbering_lexical_order',permr,npoin_halo*ndofn)
       end if
       call memory_alloca(memor_dom,'npoin_own_tot','renumbering_lexical_order',npoin_own_tot,npart+1,'INITIALIZE',0_ip)
       call PAR_ALLGATHER(npoin_own,npoin_own_tot)       
       if( INOTMASTER ) then
          npoin_offset = sum(npoin_own_tot(0:kfl_paral-1))*ndofn
          do ipoin = 1,npoin_own * ndofn          
             permr(ipoin) = ipoin + npoin_offset
          end do
          call PAR_INTERFACE_OWN_NODE_EXCHANGE(ndofn,permr)
       end if
       call memory_deallo(memor_dom,'npoin_own_tot','renumbering_lexical_order',npoin_own_tot)

    else
       !
       ! Renumber nodes without including halos
       !       
       if( .not. associated(permr) .and. INOTMASTER ) then
          call memory_alloca(memor_dom,'permr','renumbering_lexical_order',permr,npoin*ndofn)
       end if
       call memory_alloca(memor_dom,'npoin_own_tot','renumbering_lexical_order',npoin_own_tot,npart+1,'INITIALIZE',0_ip)
       call PAR_ALLGATHER(npoin_own,npoin_own_tot)
       if( INOTMASTER ) then
          npoin_offset = sum(npoin_own_tot(0:kfl_paral-1))*ndofn
          do ipoin = 1,npoin_own * ndofn
             permr(ipoin) = ipoin + npoin_offset
          end do
          call PAR_INTERFACE_NODE_EXCHANGE(ndofn,permr,'MAX','IN MY CODE')
       end if
       call memory_deallo(memor_dom,'npoin_own_tot','renumbering_lexical_order',npoin_own_tot)

    end if

  end subroutine renumbering_lexical_order_node

  !-----------------------------------------------------------------------
  !>
  !> @date    07/11/2013
  !> @author  Guillaume Houzeaux
  !> @brief   Renumber nodes in lexical order
  !> @details Renumber nodes in global lexical order using mesh type
  !>
  !-----------------------------------------------------------------------

  subroutine renumbering_lexical_order_type(meshe,permr,what,ndofn_opt)

    type(mesh_type),          intent(in)    :: meshe             !< Mesh type
    integer(ip),     pointer, intent(inout) :: permr(:)          !< Permutation
    character(*),    optional               :: what              !< Should we deallocate
    integer(ip),     optional               :: ndofn_opt         !< If dof > 1
    integer(ip),     pointer                :: npoin_own_tot(:)
    integer(ip)                             :: ipoin,ndofn
    integer(ip)                             :: npoin_offset
    logical(lg)                             :: with_halo

    call renumbering_lexical_order_node(meshe % npoin,meshe % npoin_halo,meshe % npoin_own,permr,what,ndofn_opt)
    return


    ! ESTO ES MUY FEO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!! todas las lineas de abaj son al PEDO!!!!!  -- las borro
    ! Pasar un optional que noe sta definido tampoco es lindo -- ric dice que se al banca pero no debería
    ! I thought it was incorrect but --- jimdempseyatthecove
    ! I believe that an incomming optional argument(s) can be passed as optional argumentsin a call to a nested subroutine
    ! as long as it is (they are) last in the argument list of the nested subroutine. 

    !
    ! What to do
    !
    if( present(what) ) then
       if( trim(what) == 'DEALLOCATE' ) then
          if( INOTMASTER ) call memory_deallo(memor_dom,'permr','renumbering_lexical_order',permr)
          return
       end if
    end if
    with_halo = .true.
    if( present(what) ) then
       if( trim(what) == 'WITHOUT HALOS' ) then
          with_halo = .false.
       end if
    end if
    !
    ! NDOFN > 1
    !
    if( present(ndofn_opt) ) then
       ndofn = ndofn_opt
    else
       ndofn = 1
    end if
    !
    ! Permute
    !
    nullify(npoin_own_tot)

    if( with_halo ) then

       if( .not. associated(permr) .and. INOTMASTER ) then
          call memory_alloca(memor_dom,'permr','renumbering_lexical_order',permr,meshe % npoin_halo*ndofn)
       end if
       call memory_alloca(memor_dom,'npoin_own_tot','renumbering_lexical_order',npoin_own_tot,npart+1,'INITIALIZE',0_ip)
       call PAR_ALLGATHER(meshe % npoin_own,npoin_own_tot)       
       if( INOTMASTER ) then
          npoin_offset = sum(npoin_own_tot(0:kfl_paral-1))*ndofn
          do ipoin = 1,meshe % npoin_own * ndofn          
             permr(ipoin) = ipoin + npoin_offset
          end do
          call PAR_INTERFACE_OWN_NODE_EXCHANGE(permr)
       end if
       call memory_deallo(memor_dom,'npoin_own_tot','renumbering_lexical_order',npoin_own_tot)

    else

       if( .not. associated(permr) .and. INOTMASTER ) then
          call memory_alloca(memor_dom,'permr','renumbering_lexical_order',permr,meshe % npoin*ndofn)
       end if
       call memory_alloca(memor_dom,'npoin_own_tot','renumbering_lexical_order',npoin_own_tot,npart+1,'INITIALIZE',0_ip)
       call PAR_ALLGATHER(meshe % npoin_own,npoin_own_tot)
       if( INOTMASTER ) then
          npoin_offset = sum(npoin_own_tot(0:kfl_paral-1))*ndofn
          do ipoin = 1,meshe % npoin_own * ndofn
             permr(ipoin) = ipoin + npoin_offset
          end do
          call PAR_INTERFACE_NODE_EXCHANGE(permr,'MAX','IN MY CODE')
       end if
       call memory_deallo(memor_dom,'npoin_own_tot','renumbering_lexical_order',npoin_own_tot)

    end if

  end subroutine renumbering_lexical_order_type


  !-----------------------------------------------------------------------
  !
  !> @brief   Estimate temporal locality
  !> @author  Guillaume Houzeaux
  !> @date    17/10/2017
  !> @details Estimate the temporal locality of matrix scatter obtained
  !>          during the element assbmly
  !
  !-----------------------------------------------------------------------

  subroutine renumbering_temporal_locality(meshe)

    use mod_graphs, only : graphs_find_edge

    type(mesh_type),      intent(inout) :: meshe             !< Mesh type
    integer(ip)                         :: ielem,inode,ipoin
    integer(ip)                         :: jnode,izdom
    integer(ip)                         :: jpoin,i1,i2
    integer(ip)                         :: max_nelem,ave_nelem,tot_nelem
    integer(ip)                         :: max_nedge,ave_nedge,tot_nedge
    integer(ip)                         :: max_npoin,ave_npoin,tot_npoin
    integer(ip)                         :: spatial_graph_max 
    integer(ip)                         :: spatial_graph_ave 
    integer(ip),          pointer       :: las_time(:)
    integer(ip),          pointer       :: max_time(:)
    integer(ip),          pointer       :: num_time(:)
    integer(ip),          pointer       :: las_time_edge(:)
    integer(ip),          pointer       :: max_time_edge(:)
    integer(ip),          pointer       :: num_time_edge(:)

    nullify(las_time)
    nullify(max_time)
    nullify(num_time)
    nullify(las_time_edge)
    nullify(max_time_edge)
    nullify(num_time_edge)

    if( INOTMASTER ) then
       !
       ! Spatial and temporal locality
       !
       call memory_alloca(memor_dom,'LAS_TIME','renumbering_temporal_locality',las_time,meshe % npoin)
       call memory_alloca(memor_dom,'MAX_TIME','renumbering_temporal_locality',max_time,meshe % npoin)
       call memory_alloca(memor_dom,'NUM_TIME','renumbering_temporal_locality',num_time,meshe % npoin)

       max_npoin = 0
       ave_npoin = 0
       tot_npoin = 0

       do ielem = 1,meshe % nelem
          do inode = 1,meshe % lnnod(ielem)
             ipoin = meshe % lnods(inode,ielem)
             if( num_time(ipoin) > 0 ) max_time(ipoin) = max(max_time(ipoin),ielem-las_time(ipoin))
             num_time(ipoin) = num_time(ipoin) + 1
             las_time(ipoin) = ielem
             do jnode = inode+1,meshe % lnnod(ielem)
                jpoin = meshe % lnods(jnode,ielem)
                if( ipoin > jpoin ) then
                   tot_npoin = tot_npoin + 1
                   max_npoin = max(max_npoin,ipoin-jpoin)
                   ave_npoin = ave_npoin + ipoin-jpoin
                end if
             end do
          end do
       end do

       max_nelem = maxval(max_time)
       ave_nelem = sum(max_time(1:meshe % npoin))
       tot_nelem = meshe % npoin

       call memory_deallo(memor_dom,'LAS_TIME','renumbering_temporal_locality',las_time)
       call memory_deallo(memor_dom,'MAX_TIME','renumbering_temporal_locality',max_time)
       call memory_deallo(memor_dom,'NUM_TIME','renumbering_temporal_locality',num_time)
       !
       ! Temporal locality of edges
       !
       call memory_alloca(memor_dom,'LAS_TIME_EDGE','renumbering_temporal_locality',las_time_edge,meshe % nzdom)
       call memory_alloca(memor_dom,'MAX_TIME_EDGE','renumbering_temporal_locality',max_time_edge,meshe % nzdom)
       call memory_alloca(memor_dom,'NUM_TIME_EDGE','renumbering_temporal_locality',num_time_edge,meshe % nzdom)

       do ielem = 1,meshe % nelem
          do inode = 1,meshe % lnnod(ielem)
             ipoin = meshe % lnods(inode,ielem)
             do jnode = 1,meshe % lnnod(ielem)
                jpoin = meshe % lnods(jnode,ielem)
                call graphs_find_edge(ipoin,jpoin,meshe % r_dom,meshe % c_dom,izdom)
                if( num_time_edge(izdom) > 0 ) &                   
                     max_time_edge(izdom) = max(max_time_edge(izdom),abs(ielem-las_time_edge(izdom)))
                num_time_edge(izdom) = num_time_edge(izdom) + 1
                las_time_edge(izdom) = ielem
             end do
          end do
       end do

       max_nedge = maxval(max_time_edge)
       ave_nedge = sum(max_time_edge(1:meshe % nzdom))
       tot_nedge = meshe % nzdom

       call memory_deallo(memor_dom,'LAS_TIME_EDGE','renumbering_temporal_locality',las_time_edge)
       call memory_deallo(memor_dom,'MAX_TIME_EDGE','renumbering_temporal_locality',max_time_edge)
       call memory_deallo(memor_dom,'NUM_TIME_EDGE','renumbering_temporal_locality',num_time_edge)
       !
       ! Spatial locality on graph (SpMV)
       !
       spatial_graph_max = 0
       spatial_graph_ave = 0 

       do ipoin = 1,meshe % npoin
          i1                = minval(meshe % c_dom(meshe % r_dom(ipoin):meshe % r_dom(ipoin+1)-1))
          i2                = maxval(meshe % c_dom(meshe % r_dom(ipoin):meshe % r_dom(ipoin+1)-1))
          spatial_graph_max = max(spatial_graph_max,i2-i1)
          spatial_graph_ave = spatial_graph_ave + i2-i1
       end do

    end if

    call PAR_MAX(max_npoin)
    call PAR_SUM(ave_npoin)
    call PAR_SUM(tot_npoin)

    call PAR_MAX(max_nelem)
    call PAR_SUM(ave_nelem)
    call PAR_SUM(tot_nelem)

    call PAR_MAX(spatial_graph_max)
    call PAR_SUM(spatial_graph_ave)

    call PAR_MAX(max_nedge)
    call PAR_SUM(ave_nedge)
    call PAR_SUM(tot_nedge)

    ioutp(1) = ave_npoin / tot_npoin
    ioutp(2) = max_npoin 

    ioutp(3) = ave_nelem / tot_nelem
    ioutp(4) = max_nelem

    ioutp(5) = spatial_graph_ave / tot_npoin
    ioutp(6) = spatial_graph_max

    ioutp(7) = ave_nedge / tot_nedge
    ioutp(8) = max_nedge
    !
    ! We cannot depend on outfor because of ciruclar dependency
    !
    !call outfor(86_ip,lun_outpu,' ')

  end subroutine renumbering_temporal_locality

  !-----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    24/01/2018
  !> @brief   SFC
  !> @details Renumbering using SFC
  !>
  !-----------------------------------------------------------------------

  subroutine renumbering_sfc(nsfc,coord,perm,NPOIN)

    use mod_maths,  only : maths_sfc_1d_to2d3d_tab
    use mod_maths,  only : maths_mapping_coord_to_3d
    use mod_graphs, only : graphs_number_to_linked_list

    integer(ip),          intent(in)    :: nsfc
    real(rp),    pointer, intent(in)    :: coord(:,:)
    integer(ip), pointer, intent(inout) :: perm(:)
    integer(ip),          optional      :: NPOIN
    integer(ip), pointer                :: i2dto1d(:,:,:)
    integer(ip), pointer                :: numd(:)
    integer(ip)                         :: ndime,kpoin,d,x,y,z
    integer(ip)                         :: ipoin,ii,jj,kk,ll,idime,nd
    integer(8)                          :: memor(2)
    integer(ip)                         :: boxes(3)
    real(rp)                            :: comin(3)
    real(rp)                            :: comax(3)

    nullify(i2dto1d)
    nullify(numd)

    ndime = memory_size(coord,1_ip)
    if( present(NPOIN) ) then
       kpoin = NPOIN
    else
       kpoin = memory_size(coord,2_ip)
    end if
    memor = 0_8
    nd    = nsfc**ndime

    if( .not. associated(perm) ) then
       call memory_alloca(memor,'direct_solver % mumps_par % rhs','direct_solver_factorization_mumps',perm,kpoin)
    end if

    if(      ndime == 2 ) then
       call memory_alloca(memor,'direct_solver % mumps_par % rhs','direct_solver_factorization_mumps',i2dto1d,nsfc,nsfc,1_ip)
    else if( ndime == 3 ) then
       call memory_alloca(memor,'direct_solver % mumps_par % rhs','direct_solver_factorization_mumps',i2dto1d,nsfc,nsfc,nsfc)
    end if

    call memory_alloca(memor,'direct_solver % mumps_par % rhs','direct_solver_factorization_mumps',numd,nd+1_ip)

    if( ndime == 2 ) then
       do d = 1,nd
          call maths_sfc_1d_to2d3d_tab(nsfc,d,x,y)
          i2dto1d(x,y,1) = d
       end do
    else
       do d = 1,nd
          call maths_sfc_1d_to2d3d_tab(nsfc,d,x,y,z)
          i2dto1d(x,y,z) = d
       end do
    end if

    boxes = nsfc

    if( ndime == 2 ) then
       do idime = 1,2
          comin(idime) = minval(coord(idime,:))
          comax(idime) = maxval(coord(idime,:))
       end do
       do ipoin = 1,kpoin
          call maths_mapping_coord_to_3d(ndime,boxes,comin,comax,coord(:,ipoin),ii,jj,kk)
          perm(ipoin) =  i2dto1d(ii,jj,1)
       end do
    else
       do idime = 1,3
          comin(idime) = minval(coord(idime,:))
          comax(idime) = maxval(coord(idime,:))
       end do
       do ipoin = 1,kpoin
          call maths_mapping_coord_to_3d(ndime,boxes,comin,comax,coord(:,ipoin),ii,jj,kk)
          perm(ipoin) =  i2dto1d(ii,jj,kk)
       end do
    end if

    do ipoin = 1,kpoin
       d       = perm(ipoin)
       numd(d) = numd(d) + 1   
    end do


    kk      = numd(1)
    numd(1) = 1 
    do ii = 2,nd+1
       ll       = numd(ii)
       numd(ii) = numd(ii-1) + kk
       kk       = ll
    end do

    do ipoin = 1,kpoin
       d           = perm(ipoin)
       perm(ipoin) = numd(d)
       numd(d)     = numd(d)+1
    end do

    !do ipoin = 1,kpoin
    !   pp = 0
    !   do jpoin = 1,kpoin
    !      if( perm(jpoin) == perm(ipoin) ) pp=pp+1
    !   end do
    !   if( pp /= 1 .or. perm(ipoin) <= 0 .or. perm(ipoin) > kpoin ) print*,'popo=',ipoin,pp
    !end do

    call memory_deallo(memor,'direct_solver % mumps_par % rhs','direct_solver_factorization_mumps',i2dto1d)
    call memory_deallo(memor,'direct_solver % mumps_par % rhs','direct_solver_factorization_mumps',numd)

  end subroutine renumbering_sfc

  subroutine renumbering_sfc_recursive(nsfc,coord,perm,NUMBER_OF_NODES)

    use mod_maths,  only : maths_mapping_coord_to_3d
    use mod_maths,  only : maths_sfc_1d_to2d3d_tab,maths_sfc_d2xy_tab
    use mod_graphs, only : graphs_number_to_linked_list

    type sfcbox
       integer(ip)               :: id          ! My global ID
       integer(ip)               :: level       ! Generation
       integer(ip)               :: npoinbox    ! Number of nodes
       integer(ip)               :: childid     ! Child ID 
       integer(2)                :: orientation ! Orientation  
       integer(ip)               :: whoiam      ! Father or have nodes
       integer(ip),  pointer     :: nodes(:)    ! nodes
       real(rp)                  :: minc(3)     ! Min coordinates
       real(rp)                  :: maxc(3)     ! Max coordinates
       type(sfcbox), pointer     :: parent      ! Pointer to parent
       type(sfcbox), pointer     :: children(:) ! Pointer to children
    end type sfcbox

    type sfctype
       integer(ip)               :: iallo
       type(sfcbox),pointer      :: tree_root
       integer(ip), pointer      :: kstat(:)  
       integer(ip)               :: divmax
       real(rp)                  :: comin(3)
       real(rp)                  :: comax(3)
       real(rp),    pointer      :: cputi(:)    
    end type sfctype

    type(sfcbox),   parameter    :: sfcbox_init = sfcbox(&
         0_ip,&                                             ! id         
         0_ip,&                                             ! level      
         0_ip,&                                             ! npoinbox   
         0_ip,&                                             ! childid    
         0_2,&                                              ! orientation    
         0_ip,&                                             ! whoiam     
         null(),&                                           ! nodes(:)   
         (/0.0_rp,0.0_rp,0.0_rp/),&                         ! minc(3)    
         (/0.0_rp,0.0_rp,0.0_rp/),&                         ! maxc(3)    
         null(),&                                           ! parent     
         null())                                            ! children(:)
    type(sfctype),   parameter   :: sfc_struc_init = sfctype(&
         0_ip,&                                             ! iallo      
         null(),&                                           ! tree_root              
         null(),&                                           ! kstat(:)             
         0_ip,&                                             ! divmax             
         0.0_rp,&                                           ! comin(3)        
         0.0_rp,&                                           ! comax(3)     
         null())                                            ! cputi(:)      

    integer(ip),          intent(in)    :: nsfc
    real(rp),    pointer, intent(in)    :: coord(:,:)
    integer(ip), pointer, intent(inout) :: perm(:)
    integer(ip),          optional      :: NUMBER_OF_NODES

    type(sfctype),pointer               :: sfc_struc(:) => null()
    type(sfcbox), pointer               :: tree_root
    type(sfcbox), pointer               :: current_o

    type(sfcbox), pointer               :: old_pointer
    type(sfcbox), pointer               :: tm1_pointer
    type(sfcbox), pointer               :: tm2_pointer

    integer(ip)                         :: counter
    integer(ip)                         :: divmax,limit,npoin
    integer(ip)                         :: ndime,kpoin,d,jpoin
    integer(ip)                         :: ipoin,ii,jj,kk,idime
    integer(ip)                         :: boxes(3)
    integer(2)                          :: orientation
    integer(ip), pointer                :: i2dto1d(:,:,:)

    real(rp)                            :: comin(3)
    real(rp)                            :: comax(3)

    logical(lg)                         :: conti
    !
    ! Nullify
    !
    nullify(sfc_struc)
    nullify(old_pointer)
    nullify(tm1_pointer)
    nullify(tm2_pointer)
    nullify(tree_root)
    nullify(i2dto1d)
    !
    ! Initialization
    !
    if( present(NUMBER_OF_NODES) ) then
       npoin = NUMBER_OF_NODES
    else
       npoin = memory_size(coord,2_ip)
    end if
    ndime  = memory_size(coord,1_ip)
    divmax = nsfc**ndime
    boxes  = nsfc
    limit  = 1
    comin  = 0.0_rp
    comax  = 0.0_rp
    !
    ! Mapping 3D to 1D
    !
    if(     ndime == 2 ) then
       call memory_alloca(memor_dom,'i2dto1d','renumbering_sfc_recursive',i2dto1d,nsfc,nsfc,1_ip)
    else if( ndime == 3 ) then
       call memory_alloca(memor_dom,'i2dto1d','renumbering_sfc_recursive',i2dto1d,nsfc,nsfc,nsfc)
    end if

    !--------------------------------------------------------------------
    !
    ! Init sfc_struc
    !
    !--------------------------------------------------------------------

    allocate(sfc_struc(1))
    sfc_struc(1)          = sfc_struc_init
    sfc_struc(1) % iallo  = 1 
    sfc_struc(1) % divmax = divmax
    allocate( sfc_struc(1) % cputi(10) )
    allocate( sfc_struc(1) % kstat(10) )

    !--------------------------------------------------------------------
    !
    ! Init tree_root values, fill it in with all nodes
    !
    !--------------------------------------------------------------------

    allocate( sfc_struc(1) % tree_root )
    sfc_struc(1) % tree_root =  sfcbox_init
    tree_root                => sfc_struc(1) % tree_root
    current_o                => tree_root

    call memory_alloca(memor_dom,'NODES','elsest_sfc_preprocess',current_o % nodes,npoin)

    current_o % id          =  0
    current_o % level       =  0
    current_o % whoiam      =  0
    current_o % childid     =  0
    current_o % orientation =  0_2
    current_o % npoinbox    =  npoin

    do ipoin = 1,npoin
       current_o % nodes(ipoin) = ipoin
    end do
    do idime = 1,ndime
       current_o % minc(idime) = minval(coord(idime,1:npoin))
       current_o % maxc(idime) = maxval(coord(idime,1:npoin))
    end do

    conti   = .true.
    counter = 0
    jpoin   = 0

    !--------------------------------------------------------------------
    !
    ! Fill in tree
    !
    !--------------------------------------------------------------------

    do while( conti )
       !
       ! If maximum number of points inside current box is exceeded, subdivide
       !     
       if( current_o % npoinbox > limit ) then
          allocate( current_o % children(divmax) )
          current_o % children(1:divmax) = sfcbox_init
          !
          ! Give birth to my DIVMAX children
          !
          do d = 1,divmax
             counter                               =  counter+1
             current_o % children(d) % id          =  counter
             current_o % children(d) % childid     =  d                      
             current_o % children(d) % level       =  current_o % level + 1 
             current_o % children(d) % whoiam      =  0
             current_o % children(d) % npoinbox    =  0
             current_o % children(d) % parent      => current_o
          end do
          comin = current_o % minc
          comax = current_o % maxc
          !
          ! Compute SFC coordinates according to father's orientation
          !
          if( ndime == 2 ) then
             do d = 1,divmax
                orientation = current_o % orientation
                call maths_sfc_1d_to2d3d_tab(nsfc,d,ii,jj,orientation)
                i2dto1d(ii,jj,1) = d
                current_o % children(d) % orientation = orientation
             end do
          else
             do d = 1,divmax
                orientation = current_o % orientation
                call maths_sfc_1d_to2d3d_tab(nsfc,d,ii,jj,kk,orientation)
                i2dto1d(ii,jj,kk) = d
                current_o % children(d) % orientation = orientation
             end do

          end if
          !
          ! Offer my nodes to my children
          !
          do ipoin = 1,current_o % npoinbox            
             kpoin = current_o % nodes(ipoin)
             call maths_mapping_coord_to_3d(ndime,boxes,comin,comax,coord(:,kpoin),ii,jj,kk)
             d = i2dto1d(ii,jj,kk)
             current_o % children(d) % npoinbox = current_o % children(d) % npoinbox + 1
          end do

          do d = 1,divmax
             call memory_alloca(memor_dom,'NODES','elsest_sfc_preprocess',current_o % children(d) % nodes,current_o % children(d) % npoinbox)
             current_o % children(d) % npoinbox = 0
          end do

          do ipoin = 1,current_o % npoinbox
             kpoin = current_o % nodes(ipoin)
             call maths_mapping_coord_to_3d(ndime,boxes,comin,comax,coord(:,kpoin),ii,jj,kk)
             d = i2dto1d(ii,jj,kk)
             current_o % children(d) % npoinbox = current_o % children(d) % npoinbox + 1
             current_o % children(d) % nodes(current_o % children(d) % npoinbox) = kpoin
          end do
          !
          ! Compute children bounding boxes
          !
          do d = 1,divmax
             comin(1:ndime) =  huge(1.0_rp)
             comax(1:ndime) = -huge(1.0_rp)
             if( current_o % children(d) % npoinbox > 0 ) then
                do ipoin = 1,current_o % children(d) % npoinbox
                   kpoin = current_o % children(d) % nodes(ipoin)
                   do idime = 1,ndime
                      comin(idime) = min( comin(idime) , coord(idime,kpoin) )
                      comax(idime) = max( comax(idime) , coord(idime,kpoin) )
                   end do
                end do
             end if
             current_o % children(d) % minc = comin - zeror
             current_o % children(d) % maxc = comax + zeror     
          end do

          call memory_deallo(memor_dom,'CURRENT % NODES','elsest_sfc_preprocess',current_o % nodes)


          current_o % whoiam   =  0
          current_o % npoinbox =  0
          current_o            => current_o % children(1)

       else if( current_o % id == 0 .and. current_o % npoinbox <= limit ) then
          !
          ! If the Padrino has too few elements
          !
          call memory_alloca(memor_dom,'TREE_ROOT % NODES','elsest_sfc_preprocess',current_o % nodes,npoin)
          do ipoin = 1,npoin
             current_o % nodes(ipoin) = ipoin
          end do
          conti = .false.
          current_o % npoinbox =  npoin
          current_o % whoiam   =  1
          current_o            => old_pointer

       else 
          !
          ! We have nodes
          !
          if( current_o % npoinbox /= 0 ) then
             current_o % whoiam  = 1
             do ipoin = 1,current_o % npoinbox
                kpoin = current_o % nodes(ipoin)
                jpoin = jpoin + 1
                perm(kpoin) = jpoin
             end do
          else
             current_o % whoiam  = 2
          end if
          !
          ! if limit of points inside box is not exceeded, assign elements
          !
          if( current_o % childid < divmax .and. current_o % id /= 0 ) then
             !
             ! Go to next children
             !
             tm1_pointer => current_o
             tm2_pointer => tm1_pointer%parent%children(tm1_pointer%childid+1)
             current_o   => tm2_pointer
             goto 10

          else if(current_o % childid == 0 ) then  
             !
             ! Padrino
             !
             goto 10

          else if(current_o % childid == divmax ) then
             !
             ! Last children
             !
             noparent: do while( current_o % id > 0 )
                if(current_o % parent%id == 0) then
                   conti=.false.
                   exit noparent
                else
                   if(current_o % parent%childid /=divmax) then 
                      tm1_pointer => current_o
                      tm2_pointer => tm1_pointer%parent%parent%children(tm1_pointer%parent%childid+1)
                      current_o   => tm2_pointer
                      exit
                   else 
                      current_o => current_o % parent
                   end if
                end if
             end do noparent

          else 
             !
             ! Wrong child ID
             !
             call runend('WRONG CHILD ID: '//trim(intost(current_o % childid)))    
          end if

       end if

10     continue
       old_pointer => current_o

    end do

    !--------------------------------------------------------------------
    !
    ! Deallocate
    !
    !--------------------------------------------------------------------

    call memory_deallo(memor_dom,'i2dto1d','renumbering_sfc_recursive',i2dto1d)

    divmax    =  sfc_struc(1) % divmax
    current_o => sfc_struc(1) % tree_root
    conti     =  .true.

    do while( conti )
       !
       ! First go to deepest level in first branch
       !
       do while( current_o % whoiam == 0 )
          current_o => current_o % children(1)
       end do
       !
       ! Deallocate list of elements
       !
       if( current_o % whoiam == 1 ) then
          call memory_deallo(memor_dom,'CURRENT % NODES','elsest_sfc_deallocate',current_o % nodes)
       end if

       if( current_o % childid < divmax .and. current_o % childid /= 0 ) then
          !
          ! I'm not the last child neither the Padrino
          !
          current_o => current_o % parent % children(current_o % childid+1)

       else if( current_o % childid == divmax ) then
          !
          ! I'm the last child
          !
          current_o => current_o % parent 
          deallocate( current_o % children)
          current_o % whoiam = 3

       else if( current_o % id == 0 ) then
          !
          ! I'm the Padrino: end of deallocation
          !
          deallocate( current_o)
          conti = .false.

       end if

    end do

    deallocate( sfc_struc(1) % cputi )
    deallocate( sfc_struc(1) % kstat )
    deallocate( sfc_struc )

  end subroutine renumbering_sfc_recursive

  
!> @author  Guillermo Oyarzun
!> @date    01/07/2020
!> Update of the renumbering using cuthill
  subroutine renumbering_reverse_cuthill_mckee(nv,ia,ja,permR)
    integer(ip), intent(in)             :: nv       ! #vertices
    integer(ip), pointer, intent(in)    :: ia(:)    ! ia, ja edjes grahp in csr format 
    integer(ip), pointer, intent(in)    :: ja(:)  
    integer(ip), pointer, intent(inout) :: permR(:) ! reverse permuation
    integer(ip), pointer                :: degree(:)
    integer(ip), pointer                :: neworder(:)
    integer(ip), pointer                :: treejA(:)
    integer(ip), pointer                :: mkdnode(:)
    integer(ip), pointer                :: conv(:)
    integer(ip)                         :: ntree, lvl
    integer(ip)                         :: i , j , k , l, nodeid
    integer(ip)                         :: mindeg, maxdeg, idmin
    integer(ip)                         :: suma, empty

    character(100), PARAMETER :: vacal = "numbering_reverse_cuthill_mckee"

    ntree = nv + 1
    
    nullify(degree)
    nullify(neworder)
    nullify(treejA)
    nullify(mkdnode)
    nullify(conv)


    call memory_alloca(memor_dom,'degree',vacal,degree,nv)
    call memory_alloca(memor_dom,'neworder',vacal,neworder,nv)
    call memory_alloca(memor_dom,'treejA',vacal,treejA,ntree)
    call memory_alloca(memor_dom,'mkdnode',vacal,mkdnode,nv)
    call memory_alloca(memor_dom,'conv',vacal,conv,nv)


    !Find the minimun and max degree
    mindeg = 99999999
    idmin = 1
    maxdeg = 0

    do i= 1, nv

        degree(i) = ia(i+1) - ia(i) - 1

        if(mindeg > degree(i)) then
            mindeg = degree(i)
            idmin = i
        end if

        if(maxdeg < degree(i)) then
            maxdeg = degree(i)
        end if
    end do


    !Array to mark with 1 the nodes that have been added to the tree
    mkdnode(:) = 0

    !Adding the root of the tree
    lvl = 1
    l=1
    neworder(l) = idmin
    conv(idmin) = l

    treejA(l)  = 1
    l=l+1
    treejA(l)  = 2

    mkdnode(idmin) = 1

    !Loop adding the dependencies until all nodes are in the tree
    do while ( l <= nv )

        empty = 1

        ! Seach for the connected with the previous level
        do k = treejA(lvl) , treejA(lvl+1) - 1

            nodeid = neworder(k)
            
            !Seach all the connected with nodeid except from itself  
            do j = ia(nodeid), ia(nodeid+1) - 1
                !If still not in the tree, then add it
                if(mkdnode(ja(j)) == 0) then
                    neworder(l) = ja(j)
                    conv(ja(j)) = l
                    mkdnode(ja(j)) = 1             
                    l = l + 1
                    empty = 0
                end if 
            end do
        end do

        !If new level dont provide new insertions in the tree
        if( empty == 1 .and. l < nv ) then

            !Find new min degree between the remaining nodes
            mindeg=1000000
            do i = 1, nv
                if(mindeg > degree(i) .and. mkdnode(i) < 1) then
                    mindeg = degree(i)
                    idmin  = i
                end if 
            end do

            !Add a new level from that node with min degree 
            neworder(l) = idmin
            conv(idmin) = l
            mkdnode(idmin) = 1             
            l = l + 1
            empty = 0
        end if

        lvl = lvl + 1
        treejA(lvl+1) = l 

    end do

    !Checksum if all nodes have been visited 
    suma = 0
    do i = 1, nv
        suma = suma + mkdnode(i)
    end do

    if( suma /= nv) then
         call runend("Incongruent data at subro: renumbering_reverse_cuthill_mckee")
    end if

    !Copying the permutation vector
    do i = 1, nv
        permR(i) = conv(i)
    end do

    call memory_deallo(memor_dom,'degree',vacal,degree)
    call memory_deallo(memor_dom,'mkdnode',vacal,mkdnode)
    call memory_deallo(memor_dom,'neworder',vacal,neworder)
    call memory_deallo(memor_dom,'treejA',vacal,treejA)
    call memory_deallo(memor_dom,'conv',vacal,conv)

  end subroutine renumbering_reverse_cuthill_mckee


  subroutine renumbering_reverse_cuthill_mckee_old(nv,ia,ja,permR)

    use mod_maths, only : maths_heap_sort

    integer(ip), intent(in)             :: nv       ! #vertices
    integer(ip), pointer, intent(in)    :: ia(:)    ! ia, ja edjes grahp in csr format 
    integer(ip), pointer, intent(in)    :: ja(:)  
    integer(ip), pointer, intent(inout) :: permR(:) ! reverse permuation

    integer(ip)                         :: ne       ! #edjes
    integer(ip)                         :: ii,jj, iperm
    integer(ip)                         :: ivert
    integer(ip)                         :: maxdeg, nadj  
    integer(ip), pointer                :: degree(:)
    integer(ip), pointer                :: isin(:)
    integer(ip), pointer                :: adj(:)
    integer(ip), pointer                :: adj_deg(:)

    character(100), PARAMETER :: vacal = "numbering_reverse_cuthill_mckee"

    !
    ! Evaluate degree of each vertex
    !
    ne = ia(nv+1)-1
    maxdeg = 0_ip
    if(ne /= memory_size(ja)) call runend("Incongruent data at subro: renumbering_reverse_cuthill_mckee")
    nullify(degree)
    call memory_alloca(memor_dom,'degree',vacal,degree,nv)
    degree = 0_ip
    do ii = 1,nv
       degree(ii) = degree(ii) + ia(ii+1)-ia(ii)
       if( degree(ii) > maxdeg ) maxdeg=degree(ii)
       do jj = ia(ii),ia(ii+1)-1 
          degree(ja(jj)) = degree(ja(jj)) + 1
          if( degree(ja(jj)) > maxdeg ) maxdeg=degree(ja(jj))
       enddo
    enddo
    !
    ! Find a seed,  permR(1) := one of the nodes with minimum degree
    !
    nullify(permR,isin)
    call memory_alloca(memor_dom,'permR',vacal,permR,nv)
    call memory_alloca(memor_dom,'isin',vacal,isin,nv)
    permR = 0_ip
    isin = 0_ip
    permR(1) = minloc(degree,1)
    isin(permR(1)) = 1_ip
    iperm = 1_ip
    ! 
    ! Fill permR
    !
    nullify(adj,adj_deg)
    call memory_alloca(memor_dom,'adj',vacal,adj,maxdeg)
    call memory_alloca(memor_dom,'adj_deg',vacal,adj_deg,maxdeg)

    do ivert = 1, nv
       adj = 0_ip
       adj_deg = 0_ip
       nadj = 0_ip
       !
       ! adj:=Adjacency set of permR(ivert) excluding elements already in permR (isin=1)
       !
       do jj=ia(permR(ivert)), ia(permR(ivert)+1)-1
          if(isin(ja(jj))==0_ip) then
             isin(ja(jj)) = 1_ip
             nadj = nadj + 1_ip
             adj(nadj) = ja(jj)
             adj_deg(nadj) = degree(ja(jj)) 
          endif
       enddo
       !
       ! Sort Ai
       !
       call maths_heap_sort(2_ip, nadj, adj_deg(1:nadj),ivo1=adj(1:nadj))
       !
       ! Append Ai to R
       !
       do ii = 1,nadj
          iperm = iperm + 1_ip
          permR(iperm) = adj(ii)
       enddo
       !
       ! If no more elements to continue find a new seed
       !
       if( iperm==ivert .and. ivert < nv ) then
          iperm = iperm + 1_ip
          permR( iperm) = minloc(degree,dim=1,mask = degree < 1)
       endif

    enddo
    !
    ! Deallocate memory
    !
    call memory_deallo(memor_dom,'degree',vacal,degree)
    call memory_deallo(memor_dom,'isin',vacal,isin)
    call memory_deallo(memor_dom,'adj',vacal,adj)
    call memory_deallo(memor_dom,'adj_deg',vacal,adj_deg)

  end subroutine renumbering_reverse_cuthill_mckee_old

  subroutine degree ( root, adj_num, adj_row, adj, mask, deg, iccsze, ls, &
       node_num )

    !*****************************************************************************80
    !
    !! DEGREE computes the degrees of the nodes in the connected component.
    !
    !  Discussion:
    !
    !    The connected component is specified by MASK and ROOT.
    !    Nodes for which MASK is zero are ignored.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license. 
    !
    !  Modified:
    !
    !    05 January 2003
    !
    !  Author:
    !
    !    Original FORTRAN77 version by Alan George, Joseph Liu.
    !    FORTRAN90 version by John Burkardt.
    !
    !  Reference:
    !
    !    Alan George, Joseph Liu,
    !    Computer Solution of Large Sparse Positive Definite Systems,
    !    Prentice Hall, 1981.
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) ROOT, the node that defines the connected 
    !    component.
    !
    !    Input, integer ( kind = 4 ) ADJ_NUM, the number of adjacency entries.
    !
    !    Input, integer ( kind = 4 ) ADJ_ROW(NODE_NUM+1).  Information about 
    !    row I is stored in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
    !
    !    Input, integer ( kind = 4 ) ADJ(ADJ_NUM), the adjacency structure.
    !    For each row, it contains the column indices of the nonzero entries.
    !
    !    Input, integer ( kind = 4 ) MASK(NODE_NUM), is nonzero for those nodes 
    !    which are to be considered.
    !
    !    Output, integer ( kind = 4 ) DEG(NODE_NUM), contains, for each  node in 
    !    the connected component, its degree.
    !
    !    Output, integer ( kind = 4 ) ICCSIZE, the number of nodes in the 
    !    connected component.
    !
    !    Output, integer ( kind = 4 ) LS(NODE_NUM), stores in entries 1 through 
    !    ICCSIZE the nodes in the connected component, starting with ROOT, and 
    !    proceeding by levels.
    !
    !    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
    !
    implicit none

    integer ( kind = ip ) adj_num
    integer ( kind = ip ) node_num

    integer ( kind = ip ) adj(adj_num)
    integer ( kind = ip ) adj_row(node_num+1)
    integer ( kind = ip ) deg(node_num)
    integer ( kind = ip ) i
    integer ( kind = ip ) iccsze
    integer ( kind = ip ) ideg
    integer ( kind = ip ) j
    integer ( kind = ip ) iz2
    integer ( kind = ip ) iz1
    integer ( kind = ip ) lbegin
    integer ( kind = ip ) ls(node_num)
    integer ( kind = ip ) lvlend
    integer ( kind = ip ) lvsize
    integer ( kind = ip ) mask(node_num)
    integer ( kind = ip ) nbr
    integer ( kind = ip ) node
    integer ( kind = ip ) root
    !
    !  The sign of ADJ_ROW(I) is used to indicate if node I has been considered.
    !
    ls(1) = root
    adj_row(root) = -adj_row(root)
    lvlend = 0
    iccsze = 1
    !
    !  LBEGIN is the pointer to the beginning of the current level, and
    !  LVLEND points to the end of this level.
    !
    do

       lbegin = lvlend + 1
       lvlend = iccsze
       !
       !  Find the degrees of nodes in the current level,
       !  and at the same time, generate the next level.
       !
       do i = lbegin, lvlend

          node = ls(i)
          iz1 = -adj_row(node)
          iz2 = abs ( adj_row(node+1) ) - 1
          ideg = 0

          do j = iz1, iz2

             nbr = adj(j)

             if ( mask(nbr) /= 0 ) then

                ideg = ideg + 1

                if ( 0 <= adj_row(nbr) ) then
                   adj_row(nbr) = -adj_row(nbr)
                   iccsze = iccsze + 1
                   ls(iccsze) = nbr
                end if

             end if

          end do

          deg(node) = ideg

       end do
       !
       !  Compute the current level width.
       !
       lvsize = iccsze - lvlend
       !
       !  If the current level width is nonzero, generate another level.
       !
       if ( lvsize == 0 ) then
          exit
       end if

    end do
    !
    !  Reset ADJ_ROW to its correct sign and return.
    !
    do i = 1, iccsze
       node = ls(i)
       adj_row(node) = -adj_row(node)
    end do

    return
  end subroutine degree

  !-----------------------------------------------------------------------
  !
  !> @brief   Permutation arrays using element numbering
  !> @author  Guillaume Houzeaux
  !> @date    27/09/2016
  !> @details Renu,ber the following element arrays:
  !>          LTYPB, LBOCH, LNNOB, LESUB, LNODB, 
  !>          LBINV_LOC...
  !>
  !>          NEW = PERMR(OLD)
  !
  !-----------------------------------------------------------------------

  subroutine renumbering_boundary_arrays(permr)

    use def_master,         only : lbinv_loc
    use def_master,         only : NBOUN_TYPE
    use def_domain,         only : nboun,ltypb
    use def_domain,         only : lboch,lnnob
    use def_domain,         only : lbset,mnodb
    use def_domain,         only : ltypb
    use def_domain,         only : nboun,nfiel
    use def_domain,         only : lbmsh,lblev
    use def_domain,         only : lelbo,kfl_codbo
    use def_domain,         only : kfl_field,nboun_2
    use def_domain,         only : nbset,lnodb
    use def_domain,         only : xfiel,lboel
    use def_domain,         only : ompss_boundaries
    use mod_memory,         only : memory_copy
    use mod_memory,         only : memory_deallo
    use mod_memory,         only : memory_renumber
    use mod_redistribution, only : commd_nboun

    integer(ip), pointer, intent(inout) :: permr(:)             !< Inverse permutation

    integer(ip)                         :: iboun,kk
    integer(ip)                         :: ifiel,istep
!    integer(ip)                         :: isubd
    integer(ip)                         :: iboun_old,iboun_new
    integer(ip)                         :: nsteps,ii,nboun_new              
    integer(ip), pointer                :: ltypb_tmp(:)         ! Element type
    integer(ip), pointer                :: lboch_tmp(:)         ! Element characteristic
    integer(ip), pointer                :: lnnob_tmp(:)         ! Element number of nodes
    integer(ip), pointer                :: lnodb_tmp(:,:)       ! Element connectivity
    integer(ip), pointer                :: lelbo_tmp(:)         ! Element connectivity
    integer(ip), pointer                :: kfl_codbo_tmp(:)     ! Boundary code
    integer(ip), pointer                :: lboel_tmp(:,:)       ! Element connectivity
    integer(ip), pointer                :: lbinv_tmp(:)         ! Element original numbering
    integer(ip), pointer                :: lbset_tmp(:)         ! Element set
    integer(ip), pointer                :: lbmsh_tmp(:)         ! Postprocess on original mesh
    integer(ip), pointer                :: lblev_tmp(:)         ! Postprocess on original mesh
    integer(ip), pointer                :: lrecv_tmp(:)         ! Receive permutation
    real(rp),    pointer                :: xfiel_tmp(:,:)       ! Fields

    if( IMASTER .or. nboun == 0 ) return
    !
    ! Check permutation array
    !
    if( .not. associated(permr) ) call runend('RENUMBERING_BOUNDARY_ARRAYS: PERMR NOT ASSOCIATED')
    if( maxval(permr) > nboun )   call runend('RENUMBERING_BOUNDARY_ARRAYS: PERMR NOT CORRECT 2')

    nullify(ltypb_tmp)
    nullify(lboch_tmp)
    nullify(lnnob_tmp)
    nullify(lnodb_tmp)
    nullify(lelbo_tmp)
    nullify(kfl_codbo_tmp)
    nullify(lboel_tmp)
    nullify(lbinv_tmp)
    nullify(lbset_tmp)
    nullify(lbmsh_tmp)
    nullify(lblev_tmp)
    nullify(lrecv_tmp)

    nboun_new = 0
    do iboun_old = 1,nboun
       if( permr(iboun_old) > 0 ) nboun_new = nboun_new + 1
    end do
    
    !----------------------------------------------------------------------
    !
    ! LTYPB, LNNOB, LNODB, LBOCH, and LEBNV_LOC
    !
    !----------------------------------------------------------------------

    call memory_copy  (memor_dom,'LTYPB_TMP','renumbering_element_arrays',ltypb    ,ltypb_tmp,COPY_NAME='LTYPB_TMP')
    call memory_copy  (memor_dom,'LBOCH_TMP','renumbering_element_arrays',lboch    ,lboch_tmp,COPY_NAME='LBOCH_TMP')
    call memory_copy  (memor_dom,'LNODB_TMP','renumbering_element_arrays',lnodb    ,lnodb_tmp,COPY_NAME='LNODB_TMP')
    call memory_copy  (memor_dom,'LELBO_TMP','renumbering_element_arrays',lelbo    ,lelbo_tmp,COPY_NAME='LELBO_TMP')
    call memory_copy  (memor_dom,'LBINV_TMP','renumbering_element_arrays',lbinv_loc,lbinv_tmp,COPY_NAME='LBINV_TMP')
    
    call memory_alloca(memor_dom,'LTYPB'    ,'renumbering_element_arrays',ltypb    ,nboun_new)
    call memory_alloca(memor_dom,'LBOCH'    ,'renumbering_element_arrays',lboch    ,nboun_new)
    call memory_alloca(memor_dom,'LNODB'    ,'renumbering_element_arrays',lnodb    ,mnodb,nboun_new)
    call memory_alloca(memor_dom,'LELBO'    ,'renumbering_element_arrays',lelbo    ,nboun_new)
    call memory_alloca(memor_dom,'LBINV'    ,'renumbering_element_arrays',lbinv_loc,nboun_new)


    do iboun_old = 1,nboun
       iboun_new = permr(iboun_old)       ! IBOUN_NEW: new numbering
       if( iboun_new > 0 ) then
          ltypb(iboun_new)     = ltypb_tmp(iboun_old)
          lboch(iboun_new)     = lboch_tmp(iboun_old)
          lnodb(:,iboun_new)   = lnodb_tmp(:,iboun_old)
          lelbo(iboun_new)     = lelbo_tmp(iboun_old)
          lbinv_loc(iboun_new) = lbinv_tmp(iboun_old)
       end if
    end do

    call memory_deallo(memor_dom,'LBINV_TMP','renumbering_element_arrays',lbinv_tmp)
    call memory_deallo(memor_dom,'LELBO_TMP','renumbering_element_arrays',lelbo_tmp)
    call memory_deallo(memor_dom,'LNODB_TMP','renumbering_element_arrays',lnodb_tmp)
    call memory_deallo(memor_dom,'LBOCH_TMP','renumbering_element_arrays',lboch_tmp)
    call memory_deallo(memor_dom,'LTYPB_TMP','renumbering_element_arrays',ltypb_tmp)

    !--------------------------------------------------------------------
    !
    ! LNNOB
    !
    !--------------------------------------------------------------------

    if( associated(lnnob) ) then

       call memory_copy  (memor_dom,'LNNOB_TMP','renumbering_element_arrays',lnnob,lnnob_tmp)
       call memory_alloca(memor_dom,'LNNOB'    ,'renumbering_element_arrays',lnnob,nboun_new)
       do iboun_old = 1,nboun
          iboun_new = permr(iboun_old)
          if( iboun_new > 0 ) lnnob(iboun_new) = lnnob_tmp(iboun_old)
       end do
       call memory_deallo(memor_dom,'LNNOB_TMP','renumbering_element_arrays',lnnob_tmp)

    end if

    !--------------------------------------------------------------------
    !
    ! Codes
    !
    !--------------------------------------------------------------------

    if( associated(kfl_codbo) ) then

       call memory_copy  (memor_dom,'KFL_CODBO_TMP','renumbering_element_arrays',kfl_codbo,kfl_codbo_tmp)
       call memory_alloca(memor_dom,'KFL_CODBO'    ,'renumbering_element_arrays',kfl_codbo,nboun_new)
       do iboun_old = 1,nboun
          iboun_new = permr(iboun_old)
          if( iboun_new > 0 ) kfl_codbo(iboun_new) = kfl_codbo_tmp(iboun_old)
       end do
       call memory_deallo(memor_dom,'KFL_CODBO_TMP','renumbering_element_arrays',kfl_codbo_tmp)

    end if

    !--------------------------------------------------------------------
    !
    ! Boundary to element connectivity
    !
    !--------------------------------------------------------------------

    if( associated(lboel) ) then

       call memory_copy  (memor_dom,'LBOEL_TMP','renumbering_element_arrays',lboel,lboel_tmp)
       call memory_alloca(memor_dom,'LBOEL'    ,'renumbering_element_arrays',lboel,mnodb,nboun_new)
       do iboun_old = 1,nboun
          iboun_new = permr(iboun_old)
          if( iboun_new > 0 ) lboel(:,iboun_new) = lboel_tmp(:,iboun_old)
       end do
       call memory_deallo(memor_dom,'LBOEL_TMP','renumbering_element_arrays',lboel_tmp)

    end if

    !--------------------------------------------------------------------
    !
    ! Mesh level: LBLEV 
    !
    !--------------------------------------------------------------------

    if( associated(lblev) ) then

       call memory_copy  (memor_dom,'LBLEV_TMP','renumbering_element_arrays',lblev,lblev_tmp)
       call memory_alloca(memor_dom,'LBLEV'    ,'renumbering_element_arrays',lblev,nboun_new)
       do iboun_old = 1,nboun
          iboun_new = permr(iboun_old)
          if( iboun_new > 0 ) lblev(iboun_new) = lblev_tmp(iboun_old)
       end do
       call memory_deallo(memor_dom,'LBLEV_TMP','renumbering_element_arrays',lblev_tmp)

    end if

    !--------------------------------------------------------------------
    !
    ! Postprocess on original mesh: LBMSH 
    !
    !--------------------------------------------------------------------

    if( associated(lbmsh) ) then

       do iboun_new = 1,memory_size(lbmsh)
          iboun_old = lbmsh(iboun_new)
          if( iboun_old /= 0 ) lbmsh(iboun_new) = permr(iboun_old)
       end do

    end if

    !----------------------------------------------------------------------
    !
    ! Sets: LBSET
    !
    !----------------------------------------------------------------------

    if( nbset > 0 .and. associated(lbset) ) then

       call memory_copy  (memor_dom,'LBSET_TMP','renumbering_element_arrays',lbset,lbset_tmp)
       call memory_alloca(memor_dom,'LBSET'    ,'renumbering_element_arrays',lbset,nboun_new)
       do iboun_old = 1,nboun
          iboun_new        = permr(iboun_old)
          if( iboun_new > 0 ) lbset(iboun_new) = lbset_tmp(iboun_old)
       end do
       call memory_deallo(memor_dom,'LBSET_TMP','renumbering_element_arrays',lbset_tmp)

    end if

    !----------------------------------------------------------------------
    !
    ! Fields: XFIEL
    !
    !----------------------------------------------------------------------

    do ifiel = 1,nfiel
       if( kfl_field(2,ifiel) == NBOUN_TYPE ) then
          if ( (kfl_field(6,ifiel) /= 1) .OR. (mpio_config%output%post_process%export_only) ) then 
             nsteps = kfl_field(4,ifiel)
          else
             nsteps = nsteps_fiel_ondemand
          end if
          call memory_deallo(memor_dom,'XFIEL','renumbering_element_arrays',xfiel(ifiel)%a)
          call memory_alloca(memor_dom,'XFIEL','renumbering_element_arrays',xfiel(ifiel)%a,kfl_field(1,ifiel),nboun_new,nsteps)
          do istep = 1,nsteps
             nullify(xfiel_tmp)
             call memory_alloca(memor_dom,'XFIEL_TMP','renumbering_element_arrays',xfiel_tmp,kfl_field(1,ifiel),nboun)
             xfiel_tmp(:,:) = xfiel(ifiel) % a(:,:,istep)
             do iboun_old = 1,nboun
                iboun_new        = permr(iboun_old)
                if( iboun_new > 0 ) xfiel(ifiel) % a(:,iboun_new,istep) = xfiel_tmp(:,iboun_old)
             end do
             call memory_deallo(memor_dom,'XFIEL_TMP','renumbering_element_arrays',xfiel_tmp) 
          end do
       end if
    end do

    !----------------------------------------------------------------------
    !
    ! Partitions have been subdivided into sub-subdomains: OMPSS_DOMAINS
    !
    !----------------------------------------------------------------------

    if( associated(ompss_boundaries) ) then
       call runend('RENUMBER_BOUNDARY_ARRAYS: OMPSS_BOUNDARIES NOT CODED')
       !do isubd = 1,size(ompss_boundaries,KIND=ip)
       !   do kboun = 1,size(ompss_boundaries(isubd) % elements,KIND=ip)
       !      iboun_old = ompss_boundaries(isubd) % elements(kboun)
       !      iboun_new = permr(iboun_old)
       !      ompss_boundaries(isubd) % elements(kboun) = iboun_new
       !   end do
       !end do
    end if

    !----------------------------------------------------------------------
    !
    ! Repartitioning communicator
    !
    !----------------------------------------------------------------------

    if( associated(commd_nboun % lrecv_perm) ) then
       call memory_copy  (memor_dom,'COMMD_NBOUN % LRECV_PERM','renumbering_element_arrays',commd_nboun % lrecv_perm,lrecv_tmp)
       call memory_alloca(memor_dom,'COMMD_NBOUN % LRECV_PERM','renumbering_element_arrays',commd_nboun % lrecv_perm,nboun_new)
       kk = 0
       do ii = 1,commd_nboun % lrecv_dim
          iboun     = lrecv_tmp(ii)
          iboun_new = permr(iboun)
          if( iboun_new > 0 ) then
             kk = kk + 1
             commd_nboun % lrecv_perm(kk) = iboun_new
          end if
       end do
       commd_nboun % lrecv_dim = kk
       call memory_deallo(memor_dom,'COMMD_NBOUN % LRECV_PERM','renumbering_element_arrays',lrecv_tmp)
    end if

    !----------------------------------------------------------------------
    !
    ! Deallocate memory
    !
    !----------------------------------------------------------------------

    call memory_deallo(memor_dom,'PERMR','renumbering_element_arrays',permr)

    !----------------------------------------------------------------------
    !
    ! NEW NUMBER OF BOUDARIES!
    !
    !----------------------------------------------------------------------

    nboun   = nboun_new
    nboun_2 = nboun
    
  end subroutine renumbering_boundary_arrays
  
end module mod_renumbering
!> @}

