!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Domain
!> @{
!> @file    lbouel.f90
!> @author  houzeaux
!> @date    2018-07-21
!> @brief   Compute LELBO
!> @details This routine associates a boundary element with a volume element
!>          in array LELBO
!> @} 
!-----------------------------------------------------------------------

subroutine lbouel()

  use def_kintyp
  use def_master
  use def_domain
  use mod_memory
  use mod_maths,                 only : maths_heap_sort
  use mod_messages,              only : messages_live
  use mod_elmgeo,                only : elmgeo_number_nodes
  use mod_communications_global, only : PAR_MAX
  implicit none

  integer(ip)              :: ipoin,inodb,pnode,pnodb
  integer(ip)              :: jpoin,icoin,inode,ielem
  integer(ip)              :: kboib,mnoel,iboun,kelem
!  integer(ip)              :: dummi
  integer(ip), allocatable :: lboib(:)
  type(i1p),   pointer     :: lpoel(:)
  integer(ip), pointer     :: lnoel(:)
  real(rp)                 :: time1,time2
  
  nullify(lpoel)
  nullify(lnoel)
  !
  ! Check LELBO to prevent user from doing strange things
  !
  !dummi = memory_size(lelbo)
  !call PAR_MAX(dummi,INCLUDE_ROOT=.true.)
  !if( dummi > 0 ) then
  !   if( memory_size(lelbo)> 0 ) then
  !      dummi = maxval(lelbo)
  !   else
  !      dummi = 0_ip
  !   end if
  !   call PAR_MAX(dummi,INCLUDE_ROOT=.true.)
  !   if( dummi == 0 .and. kfl_bouel == 1 ) then
  !      call runend('LBOUEL: YOU HAVE NOT PROVIDED THE BOUNDARY/ELEMENT CONNECTIVITY FILE (LELBO)')
  !   end if
  !end if
  
  call cputim(time1)

  if( kfl_bouel == 0 .and. nboun > 0 ) then
     !
     ! Allocate memory
     !
     call messages_live('COMPUTE BOUNDARY/ELEMENT CONNECTIVITY')

     call memory_alloca(memor_dom,'LPOEL','lbouel',lpoel,npoin)
     call memory_alloca(memor_dom,'LNOEL','lbouel',lnoel,npoin)
     !
     ! Compute LNOEL
     !
     do ielem = 1,nelem
        do inode = 1,elmgeo_number_nodes(ltype(ielem),lnods(:,ielem))
           ipoin = lnods(inode,ielem)
           lnoel(ipoin) = lnoel(ipoin) + 1
        end do
     end do
     !
     ! Maxmimum number of elements connected to nodes
     !
     mnoel = maxval(lnoel)
     !
     ! Construct the list of elements connected to nodes LPOEL(:) % L
     !
     do ipoin = 1,npoin
        kelem = lnoel(ipoin)
        call memory_alloca(memor_dom,'LPOEL % L','lbouel',lpoel(ipoin) % l,kelem)
        lnoel(ipoin) = 0
     end do
     
     do ielem = 1,nelem
        pnode = elmgeo_number_nodes(ltype(ielem),lnods(:,ielem))
        do inode = 1,pnode
           ipoin = lnods(inode,ielem)
           lnoel(ipoin) = lnoel(ipoin) + 1
           lpoel(ipoin) % l(lnoel(ipoin)) = ielem
        end do
     end do

     allocate(lboib(mnodb*mnoel))
     lboib = 0
     
     !$OMP PARALLEL   DO                                          &
     !$OMP SCHEDULE ( STATIC )                                    &
     !$OMP DEFAULT  ( NONE )                                      &
     !$OMP SHARED   ( nboun,nnode,ltypb,lnodb,lnoel,lpoel,lelbo,  &
     !$OMP            lnods,ltype                              )  &
     !$OMP PRIVATE  ( iboun,kboib,inodb,ipoin,kelem,icoin,inode,  &
     !$OMP            jpoin,pnode,pnodb,ielem,lboib               )

     do iboun = 1,nboun
        kboib = 0
        pnodb = nnode(ltypb(iboun))
        do inodb = 1,pnodb
           ipoin = lnodb(inodb,iboun)
           do kelem = 1,lnoel(ipoin)
              kboib        = kboib + 1
              lboib(kboib) = lpoel(ipoin) % l(kelem)
           end do
        end do 
        !
        ! Remove replicates
        !
        call maths_heap_sort(2_ip,kboib,lboib,ELIMINATE_REPLICATES=.true.)
        !
        ! Loop over elements
        !
        kelem = 0
        do while( kelem < kboib )
           kelem = kelem + 1
           icoin = 0
           ielem = lboib(kelem)
           pnode = elmgeo_number_nodes(ltype(ielem),lnods(:,ielem))
           do inodb = 1,pnodb
              ipoin = lnodb(inodb,iboun)
              inode = 0
              do while( inode < pnode )
                 inode = inode + 1
                 jpoin = lnods(inode,ielem)
                 if( ipoin == jpoin ) then
                    icoin = icoin + 1
                    inode = pnode
                 end if
              end do
           end do
           !
           ! Element IELEM has been found
           !
           if( icoin == pnodb ) then
              lelbo(iboun) = ielem
              kelem = kboib
           end if
        end do
     end do
     
     !$OMP END PARALLEL DO
     !
     ! Check errors
     !
     do iboun = 1,nboun
        if(lelbo(iboun)==0) print*,iboun,lnodb(:,iboun)
     end do
     if( nboun > 0 ) then
        kelem = minval(lelbo)
        if( kelem == 0 ) call runend('LBOUEL: SOME BOUNDARY ELEMENTS HAVE NOT BEEN FOUND')
     end if
     !
     ! Deallocate memory
     !
     call memory_deallo(memor_dom,'LNOEL','lbouel',lnoel)
     call memory_deallo(memor_dom,'LPOEL','lbouel',lpoel)
     deallocate(lboib)

     kfl_bouel = 1

  end if
  call cputim(time2)
  cpu_start(CPU_ADDTIONAL_ARRAYS) = cpu_start(CPU_ADDTIONAL_ARRAYS) + time2 - time1

end subroutine lbouel
