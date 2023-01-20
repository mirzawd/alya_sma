!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Domain
!> @{
!> @file    mod_materials.f90
!> @author  houzeaux
!> @date    2018-11-13
!> @brief   Module for materials
!> @details Several tools to treat materials
!-----------------------------------------------------------------------

module mod_materials

  use def_kintyp
  use def_elmtyp
  use def_parame
  use def_master
  use def_domain
  use def_inpout
  use mod_memchk
  use mod_memory,         only : memory_alloca
  use mod_memory,         only : memory_deallo
  use mod_communications, only : PAR_MAX
  use mod_communications, only : PAR_INTERFACE_NODE_EXCHANGE
  use mod_domain,         only : domain_memory_allocate
  use mod_domain,         only : domain_memory_deallocate
 use mod_std
  implicit none
  private

  public :: materials_on_nodes
  public :: materials_from_boundaries
  public :: materials_destructor

contains
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-11-13
  !> @brief   Destructor
  !> @details Destroy materials
  !> 
  !-----------------------------------------------------------------------

  subroutine materials_destructor()

    call domain_memory_deallocate('MATERIALS')
    
  end subroutine materials_destructor
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-11-13
  !> @brief   ???
  !> @details Compute the list of node materials LMATN(IPOIN) % L
  !> 
  !-----------------------------------------------------------------------

  subroutine materials_on_nodes()

    integer(ip)          :: ipoin,ielem,inode,imate,kmate,kpoin,ierro(2),pnode,pmate,pelty
    integer(ip)          :: kfl_mater(100)
    integer(ip), pointer :: num_materials(:)
    integer(ip), pointer :: list_materials(:)

    nullify(num_materials)
    nullify(list_materials)
    !
    ! Check non-assigned materials
    !
    if( INOTMASTER ) then
       ierro = 0
       do ielem = 1,nelem
          if(         lmate(ielem) == 0 ) then
             ierro(1) = ierro(1) + 1
          else if( lmate(ielem) > nmate ) then
             ierro(2) = lmate(ielem) 
          end if
       end do
    end if
    call PAR_MAX(2_ip,ierro)
    if( ierro(1) /= 0 ) call runend('SOME ELEMENTS WERE NOT ASSIGNED A MATERIAL')
    if( ierro(2) /= 0 ) call runend('SOME ELEMENTS HAVE A MATERIAL NUMBER= '//trim(intost(ierro(2)))//&
         ' > NUMBER OF MATERIALS= '//trim(intost(nmate)))
    !
    ! Contact elements should be declared as solid elements
    ! This contact is for ALEFOR
    !
    do ielem = 1,nelem
       if( lelch(ielem) == ELCNT ) then
          lmate(ielem) = nmate
       end if
    end do
    !
    ! Count number of materials
    !
    do imate = 1,size(kfl_mater)
       kfl_mater(imate) = 0
    end do
    kmate = 0
    do ielem = 1,nelem
       imate = lmate(ielem)
       if( kfl_mater(imate) == 0 ) then
          kmate = kmate + 1
          kfl_mater(imate) = 1
       end if
    end do
    !if( nmate <= 1 ) nmate = kmate
    !
    ! Allocate memory for LMATN and NMATN
    !
    call domain_memory_allocate('MATERIALS')
    call memgen(1_ip,npoin,0_ip)

!!$     call memory_alloca(memor_dom,'NUM_MATERIALS','materi',num_materials,npoin)
!!$     call memory_alloca(memor_dom,'NUM_MATERIALS','materi',list_materials,nmate,npoin)
!!$     do imate = 1,nmate
!!$        do ielem = 1,nelem
!!$           if( lmate(ielem) == imate ) then
!!$              do inode = 1,lnnod(ielem)
!!$                 ipoin = lnods(inode,ielem)
!!$                 num_materials(ipoin) = num_materials(ipoin) + 1
!!$                 list_materials(num_materials(ipoin),ipoin) = imate
!!$              end do
!!$           end if
!!$        end do
!!$     end do
!!$     call PAR_INTERFACE_NODE_EXCHANGE(list_materials,'MERGE')
!!$     do ipoin = 1,npoin
!!$        igene = ipoin
!!$        igen2 = num_materials(ipoin)
!!$        call domain_memory_allocate('LMATN % L',NUMBER1=ipoin,NUMBER2=igen2)
!!$        lmatn(ipoin)%l(1:igen2) = list_materials(1:igen2,ipoin)
!!$     end do
!!$     return
    !
    ! Loop over elements to create a material-per-node vector: NODEMAT
    !
    do ielem = 1,nelem
       pelty = ltype(ielem)        
       if( pelty > 0 ) then
          pnode = nnode(pelty)           
          pmate = 1
          if( nmate > 1 ) then
             pmate = lmate(ielem)
          end if
          do inode= 1,pnode
             ipoin= lnods(inode,ielem)
             ! convention: material n-1 supersedes material n
             if ((nodemat(ipoin) == -1) .or. (pmate < nodemat(ipoin)) ) then
                nodemat(ipoin) = pmate                 
             end if
          end do
       end if
    end do

    do imate = 1,nmate
       !
       ! GISCA(IPOIN) = 1: nodes with material IMATE
       !
       do ielem = 1,nelem
          if( lmate(ielem) == imate ) then
             do inode = 1,lnnod(ielem)
                ipoin = lnods(inode,ielem)
                gisca(ipoin) = 1
             end do
          end if
       end do
       call PAR_INTERFACE_NODE_EXCHANGE(gisca,'SUM')
       !
       ! KPOIN = number of nodes with material IMATE
       !
       if( npoin > 0 ) then
          kpoin = count( gisca(1:npoin) > 0 ,KIND=ip)
       else
          kpoin = 0
       end if
       igene = imate
       igen2 = kpoin
       call domain_memory_allocate('LMATN % L',NUMBER1=imate,NUMBER2=kpoin)
       !
       ! LMATN(IMATE) % L(:) = List of node with material IMATE
       !
       kpoin = 0
       do ipoin = 1,npoin
          if( gisca(ipoin) > 0 ) then
             kpoin = kpoin + 1
             lmatn(imate) % l(kpoin) = ipoin
             gisca(ipoin) = 0
          end if
       end do
       nmatn(imate) = kpoin
    end do
    call memgen(3_ip,npoin,0_ip)

    call PAR_INTERFACE_NODE_EXCHANGE(nodemat,'MIN','IN MY CODE')

  end subroutine materials_on_nodes

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-11-13
  !> @brief   Generate materials from boundaries
  !> @details Mark elements recursively to assign material MATERIALS_IMATE
  !>          to MATERIALS_NLAYE layers of elements, starting from boundaries
  !>          with code MATERIALS_ICODE
  !>
  !-----------------------------------------------------------------------

  subroutine materials_from_boundaries()
    use mod_memory, only : memory_size

    integer(ip)          :: nchange,jnode,ipoin,imate
    integer(ip)          :: ielem,inode,ilaye,iboun,inodb
    logical(lg)          :: communicate
    integer(ip), pointer :: layer_cur(:)
    integer(ip), pointer :: layer_old(:)
    logical(lg), pointer :: check_element(:)

    if( INOTMASTER .and. associated(materials_icode) ) then

       do imate = 1,memory_size(materials_icode)

          if( materials_icode(imate) > 0 ) then
             nullify(layer_cur)
             nullify(layer_old)
             nullify(check_element)
             call memory_alloca(memor_dom,'LAYER_CUR','membcs',layer_cur,npoin)
             call memory_alloca(memor_dom,'LAYER_CUR','membcs',layer_old,npoin)
             call memory_alloca(memor_dom,'LAYER_CUR','membcs',check_element,nelem)
             !
             ! First layer: mark boundary node for boundaries with code MATERIALS_ICODE 
             !
             do iboun = 1,nboun
                if( kfl_codbo(iboun) == materials_icode(imate) ) then
                   do inodb = 1,nnode(ltypb(iboun))
                      layer_cur(lnodb(inodb,iboun)) = 1
                   end do
                   ielem = lelbo(iboun)
                   lmate(ielem) = materials_imate(imate)
                   do inode = 1,nnode(ltype(ielem))
                      layer_cur(lnods(inode,ielem)) = 1
                   end do
                end if
             end do

             call PAR_INTERFACE_NODE_EXCHANGE(layer_cur,'MAX')

             check_element = .true.
             !
             ! Go recursively through the layers
             !
             do ilaye = 1,materials_nlaye(imate)-1

                communicate   = .true.

                do while ( communicate ) 
                   layer_old(1:npoin) = layer_cur(1:npoin)
                   do ielem = 1,nelem
                      if( check_element(ielem) ) then
                         inode = 1
                         do while( inode <= nnode(ltype(ielem)) )
                            if( layer_cur(lnods(inode,ielem)) == ilaye ) then
                               do jnode = 1,nnode(ltype(ielem))
                                  ipoin = lnods(jnode,ielem)
                                  if( layer_cur(ipoin) == 0 ) layer_cur(ipoin) = ilaye + 1
                               end do
                               lmate(ielem) = materials_imate(imate)
                               inode = nnode(ltype(ielem))
                               check_element(ielem) = .false.
                            end if
                            inode = inode + 1
                         end do
                      end if
                   end do
                   call PAR_INTERFACE_NODE_EXCHANGE(layer_cur,'MAX')
                   nchange = 0
                   do ipoin = 1,npoin
                      if( layer_cur(ipoin) /= layer_old(ipoin) ) nchange = nchange + 1
                   end do
                   call PAR_MAX(nchange,'IN MY CODE WITHOUT MASTER',INCLUDE_ROOT=.true.)
                   if( nchange > 0 ) then
                      communicate = .true.
                   else
                      communicate = .false.
                   end if
                end do

             end do

             call memory_deallo(memor_dom,'LAYER_CUR','membcs',layer_cur)
             call memory_deallo(memor_dom,'LAYER_CUR','membcs',layer_old)
             call memory_deallo(memor_dom,'LAYER_CUR','membcs',check_element)

          end if

       end do

    end if

  end subroutine materials_from_boundaries

end module mod_materials
!> @}
