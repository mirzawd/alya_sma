!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine bougra()
  !-----------------------------------------------------------------------
  !****f* Domain/bougra
  ! NAME
  !    bougra
  ! DESCRIPTION
  !    This routine creates the boundary graph.
  ! OUTPUT
  !    nzbou:   Number of nonzero coefficients of the graph
  !    p_r_dom: Pointer to the array of rows r_dom(npoin+1) (r_dom(ipoin) = 
  !             coefficient of the graph where row ipoin starts)
  !    p_c_dom: Pointer to the array of columns c_dom(nzbou) (c_dom (izdom)
  !             = column of the izdom coefficient of mesh graph)
  ! USED BY
  !    nsi_autint
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use mod_memchk
  use mod_memory
  use mod_domain, only : domain_memory_allocate
  implicit none
  integer(ip)          :: ncoef,inodb,ipoin,iboun,nlbop,icoef
  integer(ip)          :: mbpoi,izbou,mbopo,mtouc
  integer(ip), pointer :: nbpoi(:) 
  integer(ip), pointer :: lbopo(:) 
  integer(ip), pointer :: pbopo(:)   
  logical(lg), pointer :: touch(:) 
  !
  ! Nullify pointers
  !
  nullify(nbpoi)
  nullify(lbopo)
  nullify(pbopo)  
  nullify(touch)
  !
  ! Boundary graph has been computed
  !
  !if( kfl_crbou == 1 ) return
  !kfl_crbou = 1
  !
  ! Allocate memory for NBPOI and compute it
  !
  call memory_alloca(memor_dom,'NBPOI','bougra',nbpoi,npoin+1_ip)
  do iboun = 1,nboun
     do inodb = 1,nnode(ltypb(iboun))
        ipoin = lnodb(inodb,iboun)
        nbpoi(ipoin) = nbpoi(ipoin)+1
     end do
  end do
  !  
  ! Allocate memory for PBOPO and compute it
  !
  call memory_alloca(memor_dom,'PBOPO','bougra',pbopo,npoin+1_ip)
  pbopo(1) = 1
  do ipoin = 1,npoin
     pbopo(ipoin+1) = pbopo(ipoin) + nbpoi(ipoin)
  end do
  !
  ! Allocate memory for LBOPO and construct the list.
  !
  nlbop = pbopo(npoin+1)
  call memory_alloca(memor_dom,'lbopo','bougra',Lbopo,nlbop)
  mbopo=0
  do iboun = 1,nboun
     mbopo = mbopo + nnode(ltypb(iboun)) * nnode(ltypb(iboun))
     do inodb = 1,nnode(ltypb(iboun))
        ipoin = lnodb(inodb,iboun)
        lbopo(pbopo(ipoin)) = iboun
        pbopo(ipoin) = pbopo(ipoin) + 1
     end do
  end do
  !
  ! Recompute PELPO and maximum number of element neighbors MEPOI
  !
  pbopo(1) =  1
  mbpoi    = -1
  do ipoin = 1,npoin
     pbopo(ipoin+1) = pbopo(ipoin) + nbpoi(ipoin)
     mbpoi = max(mbpoi,nbpoi(ipoin))
  end do
  !
  ! Deallocate memory for temporary node/element connectivity
  !
  call memory_deallo(memor_dom,'NBPOI','bougra',nbpoi)
  !
  ! Total number of nonzero coefficients of the graph.
  !
  mtouc=0
  do ipoin = 1,npoin
     mtouc = max(mtouc,(pbopo(ipoin+1)-pbopo(ipoin))*mnodb)
  end do
  nzbou = 0
  call memory_alloca(memor_dom,'TOUCH','bougra',touch,max(1_ip,mtouc))
  do ipoin = 1,npoin
     nlbop = pbopo(ipoin+1) - pbopo(ipoin)
     if( nlbop > 0 ) then
        ncoef = nlbop * mnodb
        do icoef = 1,ncoef
           touch(icoef) = .false.
        end do
        call nzecob(&
             nlbop,ncoef,nzbou,lbopo(pbopo(ipoin)),touch)
     end if
  end do
  !
  ! Allocate memory for the arrays of indexes of the mesh graph. These are 
  ! organized as follows: R_BOU(IPOIN) = coefficient of the graph where IPOIN
  ! starts, C_BOU(IZDOM) = column of the IZDOM coefficient of matrix graph.
  !
  call domain_memory_allocate('R_BOU AND C_BOU')
  !
  ! Construct the array of indexes
  !
  izbou = 1
  do ipoin = 1,npoin
     nlbop = pbopo(ipoin+1) - pbopo(ipoin)
     ncoef = nlbop * mnodb
     do icoef = 1,ncoef
        touch(icoef) = .false.
     end do
     call arrinb(&
          nlbop,ncoef,lbopo(pbopo(ipoin)),touch,izbou,ipoin)
  end do
  r_bou(npoin+1) = nzbou + 1
  !
  ! Deallocate volatile memory
  !
  call memory_deallo(memor_dom,'LBOPO','bougra',lbopo)
  call memory_deallo(memor_dom,'PBOPO','bougra',pbopo)
  call memory_deallo(memor_dom,'TOUCH','bougra',touch)

end subroutine bougra
