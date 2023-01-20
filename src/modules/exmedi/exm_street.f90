!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup Exmedi
!> @{
!> @file    exm_street.f90
!> @date    12/03/2013
!> @author  Mariano Vazquez
!> @brief   Computes a rule-based fiber streeter model
!> @details Computes a rule-based fiber streeter model
!> @}
!
!
!   Changelog.
!       Code 1. The input and subroutine were modified to accept a generic boundary,
!               not a hardcoded one
!       Code 2. nstrb_exm added to know if it has to generate fiber field for one or
!               two cavities
!       Code 3. Variable stran_exm added to make the angles an user-input variable.
!
!------------------------------------------------------------------------
subroutine exm_street

  use def_kintyp
  use def_master
  use def_kermod
  use def_domain
  use def_exmedi
  use mod_memory
  use mod_kdtree
  use mod_gradie
  use mod_memory,         only : memory_alloca
  use mod_memory,         only : memory_deallo
  use mod_messages,       only : messages_live


  implicit none
  integer(ip)  :: ipoin,iboun,ielem,inode,inodb,pnode,pelty,pnodb,peltb,idime,npoib_perset,nboun_perset
  real(rp)     :: dummr(3),dista,proje(3),prome,vmodu,pival,veaux(3,2)

  real(rp)     :: fibers_angle, fang_endo, fang_epi

  integer(ip), dimension(:), pointer    :: ltypb_perset
  integer(ip), dimension(:,:), pointer  :: lnodb_perset
  integer(ip), dimension(:), pointer    :: npoin_tmp
  integer(ip), dimension(:), pointer    :: local2global
  integer(ip), dimension(:), pointer    :: npoin_tmp2

  real(rp), dimension(:,:), pointer     :: disau
  real(rp), dimension(:,:), pointer     :: disau_tmp
  real(rp), dimension(:,:,:), pointer   :: fibers_vector
  real(rp), dimension(:,:), pointer     :: coord_perset

  type(typ_kdtree)                      :: kdtree

  nullify(ltypb_perset)
  nullify(lnodb_perset)
  nullify(npoin_tmp)
  nullify(local2global)
  nullify(npoin_tmp2)
  nullify(disau)
  nullify(disau_tmp)
  nullify(fibers_vector)
  nullify(coord_perset)




  if (kfl_stree_exm == 0) return


  pival=3.141592653589793_rp


  open(unit=191,file='fibers_field.dat')
  open(unit=192,file='cell_types.dat')

  if (IPARALL)           call runend('EXM_STREET: STREETER IS ONLY FOR SEQUENTIAL PREPROCESS')
  if (nbset < 2)         call runend('EXM_STREET: AT LEAST TWO SETS ARE REQUIRED')
  if (nstrb_exm.le.1_ip) call runend('EXM_STREET: More than one set is needed') !Changelog_code 2 (in)

  call memory_alloca(mem_modul(1:2,modul),'DISAU'         ,'exm_street', disau, npoin, 3_ip)
  call memory_alloca(mem_modul(1:2,modul),'DISAU_TMP'     ,'exm_street', disau_tmp, npoin, 3_ip)
  call memory_alloca(mem_modul(1:2,modul),'FIBERS_VECTOR' ,'exm_street', fibers_vector, ndime, npoin, 3_ip)
  call memory_alloca(mem_modul(1:2,modul),'LTYPB_PERSET'  ,'exm_street', ltypb_perset, nboun)
  call memory_alloca(mem_modul(1:2,modul),'LNODB_PERSET'  ,'exm_street', lnodb_perset, mnodb, nboun)
  call memory_alloca(mem_modul(1:2,modul),'COORD_PERSET'  ,'exm_street', coord_perset, ndime, npoin)
  call memory_alloca(mem_modul(1:2,modul),'NPOIN_TMP'     ,'exm_street', npoin_tmp, npoin)
  call memory_alloca(mem_modul(1:2,modul),'LOCAL2GLOBAL'  ,'exm_street', local2global, npoin)
  call memory_alloca(mem_modul(1:2,modul),'NPOIN_TMP2'    ,'exm_street', npoin_tmp2, npoin)
  call memory_alloca(mem_modul(1:2,modul),'FIBE2_EXM'     ,'exm_street', fibe2_exm, ndime, npoin)

  !
  ! Rock 'n' roll nena...
  !--------------------------------------------------------------------------     
  ! Construct KD-Tree for set 1
  !
  npoin_tmp = 0
  do iboun = 1,nboun
     !if (lbset(iboun) == 1) then                                            ! Changelog_code 1 (out)
     if (lbset(iboun) == strbo_exm(1)) then                                  ! Changelog_code 1 (in)
        peltb = ltypb(iboun)
        pnodb = nnode(peltb)
        do inodb = 1,pnodb
           npoin_tmp(lnodb(inodb,iboun)) = 1
        end do
     end if
  end do

  local2global = 0
  npoib_perset = 0
  do ipoin = 1,npoin
     if (npoin_tmp(ipoin) == 1) then
        npoib_perset = npoib_perset + 1
        local2global(npoib_perset) = ipoin 
     end if
  end do

  npoin_tmp2 = 0
  do ipoin = 1,npoib_perset
     npoin_tmp2(local2global(ipoin)) = ipoin
  end do

  nboun_perset = 0
  do iboun =1,nboun
     !if (lbset(iboun) == 1) then                                            ! Changelog_code 1 (out)
     if (lbset(iboun) == strbo_exm(1)) then                                  ! Changelog_code 1 (in)
        nboun_perset = nboun_perset + 1
        ltypb_perset(nboun_perset) = ltypb(iboun)
        peltb = ltypb(iboun)
        pnodb = nnode(peltb)
        do inodb = 1,pnodb
           lnodb_perset(inodb,nboun_perset) = npoin_tmp2(lnodb(inodb,iboun))
        end do
     end if
  end do

  do idime = 1,ndime
     do ipoin = 1,npoib_perset
        coord_perset(idime,ipoin) = coord(idime,local2global(ipoin))
     end do
  end do

  call kdtree % init()
  !call kdtree_construct(NDIME=ndime,NBOUN=nboun_perset,NPOIN=npoib_perset,LNODB=lnodb_perset,LTYPB=ltypb_perset,COORD=coord_perset,kdtree,NDIME=ndime)
  call kdtree % construct(NDIME=ndime,NBOUN=nboun_perset,NPOIN=npoib_perset,LNODB=lnodb_perset,LTYPB=ltypb_perset,COORD=coord_perset)
  !
  ! Look for minimum distance to the surface for all boundary nodes
  ! Displacement is BVESS_SUPPO_KER on these nodes
  !

  do ipoin = 1,npoin

     if(npoin_tmp(ipoin) /= 1) then
        !
        ! Nodes which doesn't belong to set 1
        !
        call kdtree % find(coord(1:ndime,ipoin),iboun,dista,proje)
        disau_tmp(ipoin,1) = dista
     else
        !
        ! Boundary nodes
        !
        disau_tmp(ipoin,1) = 0.0_rp
     end if
  end do
  !
  ! Deallocate KD-Tree for set 1
  !
  call kdtree % deallo()
  !
  !--------------------------------------------------------------------------     
  !
  ! Construct KD-Tree for set 2
  !
  npoin_tmp = 0
  do iboun = 1,nboun
     !if (lbset(iboun) == 2) then                                            ! Changelog_code 1 (out)
     if (lbset(iboun) == strbo_exm(2)) then                                  ! Changelog_code 1 (in)
        peltb = ltypb(iboun)
        pnodb = nnode(peltb)
        do inodb = 1,pnodb
           npoin_tmp(lnodb(inodb,iboun)) = 1
        end do
     end if
  end do

  npoib_perset = 0
  do ipoin = 1,npoin
     if (npoin_tmp(ipoin) == 1) then
        npoib_perset = npoib_perset + 1
        local2global(npoib_perset) = ipoin 
     end if
  end do

  npoin_tmp2 = 0
  do ipoin = 1,npoib_perset
     npoin_tmp2(local2global(ipoin)) = ipoin
  end do

  nboun_perset = 0
  do iboun =1,nboun
     !if (lbset(iboun) == 2) then                                            ! Changelog_code 1 (out)
     if (lbset(iboun) == strbo_exm(2)) then                                  ! Changelog_code 1 (in)
        nboun_perset = nboun_perset + 1
        ltypb_perset(nboun_perset) = ltypb(iboun)
        peltb = ltypb(iboun)
        pnodb = nnode(peltb)
        do inodb = 1,pnodb
           lnodb_perset(inodb,nboun_perset) = npoin_tmp2(lnodb(inodb,iboun))
        end do
     end if
  end do

  do idime = 1,ndime
     do ipoin = 1,npoib_perset
        coord_perset(idime,ipoin) = coord(idime,local2global(ipoin))
     end do
  end do

  call kdtree % init()
  call kdtree % construct(NDIME=ndime,NBOUN=nboun_perset,NPOIN=npoib_perset,LNODB=lnodb_perset,LTYPB=ltypb_perset,COORD=coord_perset)
  !call kdtree_construct(nboun_perset,npoib_perset,lnodb_perset,ltypb_perset,coord_perset,kdtree,NDIME=ndime)
  !
  ! Look for minimum distance to the surface for all boundary nodes
  ! Displacement is BVESS_SUPPO_KER on these nodes
  !
  do ipoin = 1,npoin


     if(npoin_tmp(ipoin) /= 1) then
        !
        ! Nodes which doesn't belong to set 2 
        !
        call kdtree % find(coord(1:ndime,ipoin),iboun,dista,proje)
        disau_tmp(ipoin,2) = dista
     else
        !
        ! Boundary nodes
        !
        disau_tmp(ipoin,2) = 0.0_rp
     end if
  end do
  !
  ! Deallocate KD-Tree for set 2
  !
  call kdtree % deallo()

  if (nstrb_exm .eq. 3_ip) then ! 3 sets - bi-ventricular case           ! Changelog_code 2 (in)
     !if (nbset == 3) then ! 3 sets - bi-ventricular case                   ! Changelog_code 2 (out)
     !
     !
     !--------------------------------------------------------------------------     
     !
     ! Construct KD-Tree for set 3
     !
     npoin_tmp = 0
     do iboun = 1,nboun
        !if (lbset(iboun) == 3) then                                          ! Changelog_code 1 (out)
        if (lbset(iboun) == strbo_exm(3)) then                                ! Changelog_code 1 (in)
           peltb = ltypb(iboun)
           pnodb = nnode(peltb)
           do inodb = 1,pnodb
              npoin_tmp(lnodb(inodb,iboun)) = 1
           end do
        end if
     end do

     npoib_perset = 0
     do ipoin = 1,npoin
        if (npoin_tmp(ipoin) == 1) then
           npoib_perset = npoib_perset + 1
           local2global(npoib_perset) = ipoin 
        end if
     end do

     npoin_tmp2 = 0
     do ipoin = 1,npoib_perset
        npoin_tmp2(local2global(ipoin)) = ipoin
     end do

     nboun_perset = 0
     do iboun =1,nboun
        !if (lbset(iboun) == 3) then                                          ! Changelog_code 1 (out)
        if (lbset(iboun) == strbo_exm(3)) then                                ! Changelog_code 1 (in)
           nboun_perset = nboun_perset + 1
           ltypb_perset(nboun_perset) = ltypb(iboun)
           peltb = ltypb(iboun)
           pnodb = nnode(peltb)
           do inodb = 1,pnodb
              lnodb_perset(inodb,nboun_perset) = npoin_tmp2(lnodb(inodb,iboun))
           end do
        end if
     end do

     do idime = 1,ndime
        do ipoin = 1,npoib_perset
           coord_perset(idime,ipoin) = coord(idime,local2global(ipoin))
        end do
     end do

     call kdtree % init()
     call kdtree % construct(NDIME=ndime,NBOUN=nboun_perset,NPOIN=npoib_perset,LNODB=lnodb_perset,LTYPB=ltypb_perset,COORD=coord_perset)
     !call kdtree_construct(nboun_perset,npoib_perset,lnodb_perset,ltypb_perset,coord_perset,kdtree,NDIME=ndime)
     !
     ! Look for minimum distance to the surface for all boundary nodes
     ! Displacement is BVESS_SUPPO_KER on these nodes
     !
     do ipoin = 1,npoin


        if(npoin_tmp(ipoin) /= 1) then
           !
           ! Nodes which doesn't belong to set 3 
           !
           call kdtree % find(coord(1:ndime,ipoin),iboun,dista,proje)
           disau_tmp(ipoin,3) = dista
        else
           !
           ! Boundary nodes
           !
           disau_tmp(ipoin,3) = 0.0_rp
        end if
     end do
     !
     ! Deallocate KD-Tree for set 3
     !
     call kdtree % deallo()

  end if !end if for nbset == 3
  !
  !
  !--------------------------------------------------------------------------     
  !
  ! Compute e
  !
  if (nstrb_exm .eq. 2_ip) then ! 2 sets - one-ventricular case           ! Changelog_code 2 (in)
     !if (nbset == 2) then ! 2 sets - 1 ventricle case                       ! Changelog_code 2 (out)

     do ipoin = 1,npoin

        dista = abs(disau_tmp(ipoin,1)) / (abs(disau_tmp(ipoin,2)) + abs(disau_tmp(ipoin,1)))        
        disau(ipoin,1) = abs(dista)
        disau(ipoin,2) = 0.0_rp
        disau(ipoin,3) = 0.0_rp  

        ! cell types definition (per node)
        if(dista .le. 0.0_rp) then !endo  was lt 0.33333
           write(192,101) ipoin,1
        else if((dista .gt. 0.0_rp) .and. (dista .le. 0.60_rp)) then !mid  was ge 0.33333  and le 0.66666
           write(192,101) ipoin,2
        else if(dista .gt. 0.60_rp) then !epi  was gt 0.66666
           write(192,101) ipoin,3
        end if

     end do

  else if (nstrb_exm .eq. 3_ip) then ! 2 sets - bi-ventricular case           ! Changelog_code 2 (in)
     !else if (nbset == 3) then ! 3 sets - bi-ventricular case                   ! Changelog_code 2 (out)

     do ipoin = 1,npoin

        !node region definition (the region where each node belongs to):
        if ((abs(disau_tmp(ipoin,3)) >= abs(disau_tmp(ipoin,1))) .and. (abs(disau_tmp(ipoin,3)) >= abs(disau_tmp(ipoin,2)))) then !left ventricle
           disau(ipoin,1) = disau_tmp(ipoin,2) !dist_endo = dist set 2 (left vent.)
           disau(ipoin,2) = disau_tmp(ipoin,1) !dist_epi = dist set 1 (epi)
        else if ((abs(disau_tmp(ipoin,2)) >= abs(disau_tmp(ipoin,1))) .and. (abs(disau_tmp(ipoin,2)) >= abs(disau_tmp(ipoin,3)))) then !right ventricle
           disau(ipoin,1) = disau_tmp(ipoin,3) !dist_endo = dist set 3 (right vent.)
           disau(ipoin,2) = disau_tmp(ipoin,1) !dist_epi = dist set 1 (epi)
        else if ((abs(disau_tmp(ipoin,1)) >= abs(disau_tmp(ipoin,2))) .and. (abs(disau_tmp(ipoin,1)) >= abs(disau_tmp(ipoin,3)))) then !septum
           if (abs(disau_tmp(ipoin,3)) > abs(disau_tmp(ipoin,2))) then !left septum
              disau(ipoin,1) = disau_tmp(ipoin,2) !dist_endo = dist set 2 (left vent.)
              disau(ipoin,2) = disau_tmp(ipoin,3) !dist_epi = dist set 3 (right vent.)
           else !right septum
              disau(ipoin,1) = disau_tmp(ipoin,3) !dist_endo = dist set 3 (right vent.)
              disau(ipoin,2) = disau_tmp(ipoin,2) !dist_epi = dist set 2 (left vent.)
           end if
        else
           disau(ipoin,1) = 0.0_rp !dist_endo
           disau(ipoin,2) = 1.0_rp !dist_epi
        end if

        dista = abs(disau(ipoin,1)) / (abs(disau(ipoin,2)) + abs(disau(ipoin,1)))        
        disau(ipoin,1) = dista
        disau(ipoin,2) = 0.0_rp
        disau(ipoin,3) = 0.0_rp   

        ! cell types definition (per node)
        if(dista .le. 0.0_rp) then !endo was .lt. 0.33333
           write(192,101) ipoin,1
        else if((dista .gt. 0.0_rp) .and. (dista .le. 0.60_rp)) then !mid was .ge. 0.33333  .le. 0.66666
           write(192,101) ipoin,2
        else if(dista .gt. 0.60_rp) then !epi  was .gt. 0.66666
           write(192,101) ipoin,3
        end if

     end do

  end if ! end of calculation 

  !
  ! Compute mean e
  !
  do ielem = 1,nelem
     pelty = ltype(ielem)
     pnode = nnode(pelty)
     prome = 0.0_rp
     do inode = 1,pnode
        ipoin = lnods(inode,ielem)
        prome = prome + disau(ipoin,1)
     end do
     prome = prome / real(pnode,rp)
     do inode = 1,pnode
        ipoin = lnods(inode,ielem)
        disau(ipoin,2) = disau(ipoin,2) + prome
        disau(ipoin,3) = disau(ipoin,3) + 1.0_rp
     end do
  end do

  do ipoin= 1,npoin
     disau(ipoin,2) = disau(ipoin,2)/disau(ipoin,3)
  end do

  !
  !  Compute mean e gradient 
  !
  call grasca(disau(1:npoin,2),fibers_vector(1:ndime,1:npoin,1))

  !
  !  Compute the local basis (u,v,w)
  !
  do ipoin= 1,npoin

     dummr(1:ndime) = fibers_vector(1:ndime,ipoin,1)
     vmodu= sqrt(dummr(1)*dummr(1) + dummr(2)*dummr(2) + dummr(3)*dummr(3))
     fibers_vector(1,ipoin,1) = fibers_vector(1,ipoin,1)/vmodu
     fibers_vector(2,ipoin,1) = fibers_vector(2,ipoin,1)/vmodu
     fibers_vector(3,ipoin,1) = fibers_vector(3,ipoin,1)/vmodu
     !
     ! fibers_vector(:,ipoin,1) = grad(e) -- u
     ! fibers_vector(:,ipoin,1) CAN NOT BE COLINEAL with fiaxe_exm (the apex-to-base vector)
     !

     call vecpro(fibers_vector(1,ipoin,1),fiaxe_exm,fibers_vector(1,ipoin,2),ndime)
     dummr(1:ndime) = fibers_vector(1:ndime,ipoin,2)
     vmodu= sqrt(dummr(1)*dummr(1) + dummr(2)*dummr(2) + dummr(3)*dummr(3))
     fibers_vector(1,ipoin,2) = fibers_vector(1,ipoin,2)/vmodu
     fibers_vector(2,ipoin,2) = fibers_vector(2,ipoin,2)/vmodu
     fibers_vector(3,ipoin,2) = fibers_vector(3,ipoin,2)/vmodu
     !
     ! fibers_vector(:,ipoin,2) = v
     !

     call vecpro(fibers_vector(1,ipoin,1),fibers_vector(1,ipoin,2),fibers_vector(1,ipoin,3),ndime)
     dummr(1:ndime) = fibers_vector(1:ndime,ipoin,3)
     vmodu= sqrt(dummr(1)*dummr(1) + dummr(2)*dummr(2) + dummr(3)*dummr(3))
     fibers_vector(1,ipoin,3) = fibers_vector(1,ipoin,3)/vmodu
     fibers_vector(2,ipoin,3) = fibers_vector(2,ipoin,3)/vmodu
     fibers_vector(3,ipoin,3) = fibers_vector(3,ipoin,3)/vmodu
     !
     ! fibers_vector(:,ipoin,3) = w
     !


     !fibers_factor= stran_exm*pival/180.0_rp !pival/3.0_rp ! R value = pi/3 !Changelog_code 3 (in) ! Changelog_code 4 (out)

     fang_endo = stran_endo_exm*pival/180.0_rp
     fang_epi = stran_epi_exm*pival/180.0_rp

     !fibers_angle= fibers_factor*(1.0_rp - 2.0_rp * disau(ipoin,1))**kfl_stree_exm ! alpha = R*(1-2e)**kfl_stree_exm  !Changelog_code 4 (out)


     fibers_angle= ( fang_endo * ( 1.0_rp - disau(ipoin,1) ) + fang_epi * disau(ipoin,1) )**kfl_stree_exm ! variation from endo(0) to epi (1)  !Changelog_code 4 (in)



     veaux(1,1) = fibers_vector(1,ipoin,2) * cos(fibers_angle) + fibers_vector(1,ipoin,3) * sin(fibers_angle)
     veaux(2,1) = fibers_vector(2,ipoin,2) * cos(fibers_angle) + fibers_vector(2,ipoin,3) * sin(fibers_angle)
     veaux(3,1) = fibers_vector(3,ipoin,2) * cos(fibers_angle) + fibers_vector(3,ipoin,3) * sin(fibers_angle)

     veaux(1,2) = -fibers_vector(1,ipoin,2) * sin(fibers_angle) + fibers_vector(1,ipoin,3) * cos(fibers_angle)
     veaux(2,2) = -fibers_vector(2,ipoin,2) * sin(fibers_angle) + fibers_vector(2,ipoin,3) * cos(fibers_angle)
     veaux(3,2) = -fibers_vector(3,ipoin,2) * sin(fibers_angle) + fibers_vector(3,ipoin,3) * cos(fibers_angle)

     fibers_vector(1:ndime,ipoin,2) = veaux(1:ndime,1)
     fibers_vector(1:ndime,ipoin,3) = veaux(1:ndime,2)

     dummr(1:ndime) = fibers_vector(1:ndime,ipoin,3)
     vmodu= sqrt(dummr(1)*dummr(1) + dummr(2)*dummr(2) + dummr(3)*dummr(3))
     fibers_vector(1,ipoin,3) = fibers_vector(1,ipoin,3)/vmodu
     fibers_vector(2,ipoin,3) = fibers_vector(2,ipoin,3)/vmodu
     fibers_vector(3,ipoin,3) = fibers_vector(3,ipoin,3)/vmodu

     dummr(1:ndime) = fibers_vector(1:ndime,ipoin,2)
     vmodu= sqrt(dummr(1)*dummr(1) + dummr(2)*dummr(2) + dummr(3)*dummr(3))
     fibers_vector(1,ipoin,2) = fibers_vector(1,ipoin,2)/vmodu
     fibers_vector(2,ipoin,2) = fibers_vector(2,ipoin,2)/vmodu
     fibers_vector(3,ipoin,2) = fibers_vector(3,ipoin,2)/vmodu

     write(191,100) ipoin, fibers_vector(1:3,ipoin,2)  ! along fiber

  end do

100 format(i6,3(2x,e16.8e3))
101 format(i8,i2)

  close(unit=191)  
  close(unit=192)

  do idime = 1,3
     do ipoin = 1,npoin
        fibe2_exm(idime,ipoin) = fibers_vector(idime,ipoin,2) !along fiber, fiber postprocessing in Alya (FIBE2)
     end do
  end do


  call memory_deallo(mem_modul(1:2,modul),'DISAU'         ,'exm_street', disau)
  call memory_deallo(mem_modul(1:2,modul),'DISAU_TMP'     ,'exm_street', disau_tmp)
  call memory_deallo(mem_modul(1:2,modul),'FIBERS_VECTOR' ,'exm_street', fibers_vector)
  call memory_deallo(mem_modul(1:2,modul),'LTYPB_PERSET'  ,'exm_street', ltypb_perset)
  call memory_deallo(mem_modul(1:2,modul),'LNODB_PERSET'  ,'exm_street', lnodb_perset)
  call memory_deallo(mem_modul(1:2,modul),'COORD_PERSET'  ,'exm_street', coord_perset)
  call memory_deallo(mem_modul(1:2,modul),'NPOIN_TMP'     ,'exm_street', npoin_tmp)
  call memory_deallo(mem_modul(1:2,modul),'LOCAL2GLOBAL'  ,'exm_street', local2global)
  call memory_deallo(mem_modul(1:2,modul),'NPOIN_TMP2'    ,'exm_street', npoin_tmp2)
  call memory_deallo(mem_modul(1:2,modul),'FIBE2_EXM'     ,'exm_street', fibe2_exm)


  call messages_live('MODULE DATA','END SECTION')
  call runend('O.K.!')

end subroutine exm_street
