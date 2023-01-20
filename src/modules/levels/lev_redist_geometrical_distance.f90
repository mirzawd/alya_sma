!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Levels
!> @{
!> @file    lev_redist_generalized_distance.f90
!> @author  houzeaux
!> @date    2020-09-07
!> @brief   Redistantiation
!> @details Redistantiation based on the generlized distance
!> @} 
!-----------------------------------------------------------------------

subroutine lev_redist_geometrical_distance()

  use def_kintyp_mesh_basic, only : mesh_type_basic
  use def_elmtyp
  use def_master
  use def_domain 
  use def_kermod
  use def_levels
  use mod_elmgeo
  use mod_memory
  use mod_solver
  use mod_elsest,                only : elsest_host_element
  use mod_generalized_distance,  only : generalized_distance_assembly
  use mod_generalized_distance,  only : generalized_distance_solution
  use mod_generalized_distance,  only : generalized_distance_update
  use mod_messages,              only : messages_live
  use mod_communications_global, only : PAR_SUM
  use mod_communications_global, only : PAR_ALLGATHER
  use mod_communications_global, only : PAR_ALLGATHERV
  use mod_mesh_type_basic,       only : mesh_type_basic_parallel
  use mod_mesh_type_basic,       only : mesh_type_basic_output
  use mod_parall,                only : PAR_COMM_MY_CODE
  use mod_kdtree,                only : kdtree_construct
  use mod_kdtree,                only : typ_kdtree
  use mod_kdtree,                only : kdtree_nearest_boundary
  use mod_kdtree,                only : kdtree_deallocate
  use mod_kdtree,                only : kdtree_initialize
  use mod_maths,                 only : maths_cross_product
  implicit none

  type(mesh_type_basic) :: mesh
  integer(ip)           :: ipoin,ielem
  integer(ip)           :: inode
  integer(ip)           :: iboun
  integer(ip)           :: kpoin,kelem
  real(rp)              :: den
  real(rp)              :: proje(3),v1(3),v2(3)

  integer(ip)           :: mnode_sum
  integer(ip)           :: nelem_sum
  integer(ip)           :: npoin_sum
  real(rp),     pointer :: dista(:)
  integer(ip),  pointer :: nelem_gat(:)
  integer(ip),  pointer :: npoin_gat(:)
  integer(ip),  pointer :: lnods_gat(:,:)
  integer(ip),  pointer :: ltype_gat(:)
  real(rp),     pointer :: coord_gat(:,:)
  integer(4),   pointer :: lnods_count_gat4(:)
  integer(4),   pointer :: coord_count_gat4(:)
  logical(lg),  pointer :: lmask(:)
  type(typ_kdtree)      :: kdtree
  
  nullify(nelem_gat)
  nullify(npoin_gat)
  nullify(lnods_gat)
  nullify(ltype_gat)
  nullify(coord_gat)
  nullify(lnods_count_gat4)
  nullify(coord_count_gat4)
  nullify(lmask)
  !
  ! Create mesh
  !
  !do ipoin = 1,npoin
  !   fleve(ipoin,1)=cos(coord(1,ipoin))+coord(2,ipoin)-2.0_rp
  !end do
  dista => fleve(:,1)
  call mesh % init()
  call mesh % cut_level(meshe(ndivi),dista,mem_modul(1:2,modul))

  nelem_sum = mesh % nelem
  call PAR_SUM(nelem_sum)
  call mesh_type_basic_parallel(mesh,PAR_COMM_MY_CODE,kfl_paral)
  if( 1 == 2 ) then
     mesh % name = 'LEVEL'
     call mesh_type_basic_output(mesh)
  end if

  if( nelem_sum > 0 ) then
     
     call messages_live('REDISTANCIATION USING GEOMETRICAL DISTANCE')
     !
     ! Gather geometry
     !
     call memory_alloca(mem_modul(1:2,modul),'NELEM_GAT'       ,'lev_redist',nelem_gat       ,npart+1_ip       ,LBOUN=0_ip)
     call memory_alloca(mem_modul(1:2,modul),'NPOIN_GAT'       ,'lev_redist',npoin_gat       ,npart+1_ip       ,LBOUN=0_ip)
     call memory_alloca(mem_modul(1:2,modul),'COORD_COUNT_GAT4','lev_redist',coord_count_gat4,int(npart+1_ip,4),LBOUN=0_4)
     call memory_alloca(mem_modul(1:2,modul),'LNODS_COUNT_GAT4','lev_redist',lnods_count_gat4,int(npart+1_ip,4),LBOUN=0_4)
     call PAR_ALLGATHER(mesh % nelem,nelem_gat,1_4)
     call PAR_ALLGATHER(mesh % npoin,npoin_gat,1_4)
     
     mnode_sum        = ndime
     npoin_sum        = sum(npoin_gat)
     lnods_count_gat4 = int(nelem_gat,4) * int(mnode_sum,4)
     coord_count_gat4 = int(npoin_gat,4) * int(ndime,4)
     !
     ! LNODS
     !
     kpoin = sum(npoin_gat(0:kfl_paral-1))
     do ielem = 1,mesh % nelem
        do inode = 1,mnode_sum
           mesh % lnods(inode,ielem) = mesh % lnods(inode,ielem) + kpoin
        end do
     end do
     call memory_alloca(mem_modul(1:2,modul),'LNODS_GAT','lev_redist',lnods_gat,mnode_sum,nelem_sum)     
     call PAR_ALLGATHERV(mesh % lnods,lnods_gat,lnods_count_gat4)
     !
     ! COORD
     !
     call memory_alloca(mem_modul(1:2,modul),'COORD_GAT','lev_redist',coord_gat,ndime,npoin_sum)     
     call PAR_ALLGATHERV(mesh % coord,coord_gat,coord_count_gat4)
     !
     ! LTYPE
     !
     call memory_alloca(mem_modul(1:2,modul),'LTYPE_GAT','lev_redist',ltype_gat,nelem_sum)     
     if( ndime == 2 ) then
        ltype_gat = BAR02
     else
        ltype_gat = TRI03
     end if
     !
     ! Remove null elements
     !
     call memory_alloca(mem_modul(1:2,modul),'LMASK','lev_redist',lmask,nelem_sum)
     proje = 0.0_rp
     if( ndime == 2 ) then
        do ielem = 1,nelem_sum
           v1(1:ndime) = coord_gat(1:ndime,lnods_gat(2,ielem))-coord_gat(1:ndime,lnods_gat(1,ielem)) 
           den         = sqrt(dot_product(v1(1:2),v1(1:2)))
           if( den > zeror ) then
              lmask(ielem) = .true.
           else
              lmask(ielem) = .false.
           end if
        end do
     else
        kelem = 0
        do ielem = 1,nelem_sum
           v1(1:ndime) = coord_gat(1:ndime,lnods_gat(2,ielem))-coord_gat(1:ndime,lnods_gat(1,ielem)) 
           v2(1:ndime) = coord_gat(1:ndime,lnods_gat(3,ielem))-coord_gat(1:ndime,lnods_gat(1,ielem)) 
           proje       = maths_cross_product(v1,v2,ndime)
           den         = sqrt(dot_product(proje,proje))
           if( den > zeror ) then
              lmask(ielem) = .true.
           else
              lmask(ielem) = .false.
           end if
        end do
     end if
     !
     ! Find minimum distance using kdtree
     !
     call kdtree_initialize(kdtree)
     if(npoin>0) call kdtree_construct(nelem_sum,npoin_sum,lnods_gat,ltype_gat,coord_gat,kdtree)
     do ipoin = 1,npoin
        call kdtree_nearest_boundary(coord(:,ipoin),kdtree,iboun,dista_lev(ipoin),proje,LMASK=lmask)
        dista_lev(ipoin) = abs(dista_lev(ipoin))
     end do
     call kdtree_deallocate(kdtree)
     !
     ! Correct level set
     !
     do ipoin = 1,npoin
        if( kfl_fixno_lev(1,ipoin) < 1 ) then
           if( fleve(ipoin,1) >= 0.0_rp ) then
              fleve(ipoin,1) = dista_lev(ipoin)
           else
              fleve(ipoin,1) =-dista_lev(ipoin)
           endif
        end if
     end do
  end if
  !
  ! Deallocate 
  !
  call mesh % deallo()
  call memory_deallo(mem_modul(1:2,modul),'LMASK'           ,'lev_redist',lmask           )
  call memory_deallo(mem_modul(1:2,modul),'NELEM_GAT'       ,'lev_redist',nelem_gat       )
  call memory_deallo(mem_modul(1:2,modul),'NPOIN_GAT'       ,'lev_redist',npoin_gat       )
  call memory_deallo(mem_modul(1:2,modul),'COORD_COUNT_GAT4','lev_redist',coord_count_gat4)
  call memory_deallo(mem_modul(1:2,modul),'LNODS_COUNT_GAT4','lev_redist',lnods_count_gat4)
  call memory_deallo(mem_modul(1:2,modul),'LTYPE_GAT'       ,'lev_redist',ltype_gat       )     
  call memory_deallo(mem_modul(1:2,modul),'LNODS_GAT'       ,'lev_redist',lnods_gat       )     
  call memory_deallo(mem_modul(1:2,modul),'COORD_GAT'       ,'lev_redist',coord_gat       )     

end subroutine lev_redist_geometrical_distance
