!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Kernel
!> @{
!> @file    mod_ker_noslwa.f90
!> @author  Herbert Owen
!> @brief   Preliminary operations to apply wall law increasing viscosity in the first element and using no slip 
!> @details - obtains: kount_nsw_ele_ker, normal_nsw_ker(ndime,kount_nsw_ele), kfl_nswel_ker(nelem), kfl_nswpo_ker(npoin) 
!>          - kount_nsw_ele_ker    : Number of elements where additional viscosity for no slip wall law must be added.
!>          - normal_nsw_ker       : normal of those elements. If more than 1 no slip wall law boundary reaches that element we use the average of all such boundaries. 
!>          - avta1_nsw_ker        : running time average of (mu+mut) du_t/dn (initialized to 0 here)
!>          - kfl_nswel_ker(nelem) : pointer - if /=0 indicates the number of number it occupies in the list of non slip wall law elements 
!>          - kfl_nswpo_ker(npoin) : Nodal flag used in nsi_wallav to calculate avupo_ker only where needed.
!> @} 
!-----------------------------------------------------------------------
module mod_ker_noslwa
  use def_kintyp,                   only : ip,rp
  implicit none
  private
  public :: ker_noslwa

  contains
subroutine ker_noslwa()
  
  use def_kintyp,                   only : ip,rp
  use def_master,                   only : INOTMASTER
  use def_kermod,                   only : kfl_noslw_ker, kfl_nswel_ker, kfl_nswpo_ker, kount_nsw_ele_ker, normal_nsw_ker
  use def_kermod,                   only : kfl_fixbo_nsw_ker,lnsw_exch
  use def_master,                   only : modul,mem_modul
  use def_domain,                   only : lnods,lmate,nelem,npoin,nboun,ltypb,nnode,lelbo,ltype,ngaus,coord,lnodb
  use def_domain,                   only : ndimb,elmar,ndime,mnodb,mnode
  use mod_memory,                   only : memory_alloca,memory_deallo
  use mod_communications,           only : PAR_INTERFACE_NODE_EXCHANGE
  use mod_bouder
  implicit none

  integer(ip)              :: ielem,inode,ipoin,idime
  integer(ip)              :: pnode,iboun,igaub,inodb
  integer(ip)              :: pelty,pblty,pnodb,pgaub,pmate
  real(rp)                 :: bocod(ndime,mnodb),elcod(ndime,mnode)
  real(rp)                 :: baloc(ndime,ndime),eucta,normal_aux(ndime)
  integer(ip), pointer     :: kount_boun(:)   ! number of boundaries that touch an element

  nullify(kfl_nswel_ker)
  nullify(kfl_nswpo_ker)
  nullify(normal_nsw_ker)
  nullify(kount_boun)

  if ( kfl_noslw_ker /= 0_ip ) then 

     if( INOTMASTER ) then     ! pensar que pasa con el master
        call memory_alloca(mem_modul(1:2,modul),'KFL_NSWEL_KER','ker_noslwa',kfl_nswel_ker,nelem)
        call memory_alloca(mem_modul(1:2,modul),'KFL_NSWPO_KER','ker_noslwa',kfl_nswpo_ker,npoin)
        !
        ! Loop over boundaries - preliminary just to obtain kount
        !
        kfl_nswel_ker = 0_ip
        lnsw_exch(1:nelem)%nbogp = 0_ip
        boun0: do iboun = 1,nboun

           if(  kfl_fixbo_nsw_ker(iboun) ==  1_ip ) then      ! Wall law with no slip
              !
              ! Element properties and dimensions
              !
              pblty = ltypb(iboun) 
              pnodb = nnode(pblty)
              ielem = lelbo(iboun)
              kfl_nswel_ker(ielem) = 1_ip
              lnsw_exch(ielem)%nbogp = lnsw_exch(ielem)%nbogp + ngaus(pblty)
           end if

        end do boun0
        
        kount_nsw_ele_ker = 0_ip
        elem0: do ielem = 1,nelem

           if(  kfl_nswel_ker(ielem) /= 0_ip ) then
              !
              ! Element properties and dimensions
              !
              kount_nsw_ele_ker = kount_nsw_ele_ker + 1_ip
              kfl_nswel_ker(ielem)  = ielem  ! lo reuso 
              pelty = ltype(ielem) 
              pnode = nnode(pelty)
              do inode = 1,pnode
                 ipoin = lnods(inode,ielem)
                 kfl_nswpo_ker(ipoin) = 1_ip 
              end do
              lnsw_exch(ielem)%fact = 1.0_rp / real(lnsw_exch(ielem)%nbogp,rp)
           end if

        end do elem0

        call PAR_INTERFACE_NODE_EXCHANGE(kfl_nswpo_ker,'MAX','IN MY CODE')
        

     end if

     !
     ! Loop over boundaries
     !
     if( INOTMASTER ) then
        call memory_alloca(mem_modul(1:2,modul),'NORMAL_NSW_KER','ker_noslwa',normal_nsw_ker,ndime,nelem)
        call memory_alloca(mem_modul(1:2,modul),'KOUNT_BOUN','ker_noslwa',    kount_boun,nelem)

        boundaries: do iboun =1, nboun
           
           if(  kfl_fixbo_nsw_ker(iboun) ==  1_ip ) then      ! Wall law with no slip
              !
              ! Element properties and dimensions
              !
              pblty = ltypb(iboun) 
              pnodb = nnode(pblty)
              ielem = lelbo(iboun)
              pelty = ltype(ielem)
              pnode = nnode(pelty)
              pgaub = ngaus(pblty) 
              pmate = lmate(ielem)

              if ( ( pmate /= -1 ) )  then 
                 !
                 ! Gather operations: ELCOD, BOCOD
                 !
                 do inode = 1,pnode
                    ipoin = lnods(inode,ielem)
                    do idime = 1,ndime
                       elcod(idime,inode) = coord(idime,ipoin)             
                    end do
                 end do

                 do inodb = 1,pnodb     ! obtain bocod for bouder
                    ipoin = lnodb(inodb,iboun)
                    do idime = 1,ndime
                       bocod(idime,inodb) = coord(idime,ipoin)
                    end do
                 end do

                 normal_aux = 0_rp
                 gauss_points: do igaub = 1,pgaub
                    !
                    ! Obtain normal (baloc(:,ndime) to the surface (following nsi_bouset)
                    !                 
                    call bouder(&
                         pnodb,ndime,ndimb,elmar(pblty)%deriv(:,:,igaub),&    ! Cartesian derivative
                         bocod,baloc,eucta)                                   ! and Jacobian
                    call chenor(pnode,baloc,bocod,elcod)                      ! Check normal

                    normal_aux = normal_aux + baloc(:,ndime)

                 end do gauss_points
                 !
                 ! Just the average over all gauss points - might not be the best for high order elements
                 !
                 kount_boun(ielem)       = kount_boun(ielem) + 1_ip
                 normal_nsw_ker(:,ielem) = normal_nsw_ker (:,ielem) + ( normal_aux / real(pgaub,rp))

                 ! beware check taht it is normalized!!!
              end if

           end if

        end do boundaries

        do ielem = 1,nelem
           if( kount_boun(ielem) > 0 ) &
                normal_nsw_ker(:,ielem) = normal_nsw_ker(:,ielem) / real(kount_boun(ielem),rp)
        end do
        call memory_deallo(mem_modul(1:2,modul),'KOUNT_BOUN','ker_noslwa',kount_boun)
        
     end if  ! not master

  end if 

end subroutine ker_noslwa

end module mod_ker_noslwa
