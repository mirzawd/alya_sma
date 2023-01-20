!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine ale_bodfit()
  !-----------------------------------------------------------------------
  !****f* Alefor/ale_bodfit
  ! NAME
  !    ale_bodfit
  ! DESCRIPTION
  !    This routines computes the body fitted contribution ! borrowed from ibm_bodfit and adapted
  ! USED BY
  !    domain
  !***
  !-----------------------------------------------------------------------
  use def_master
  use def_domain
  use def_alefor
  use mod_memchk
  use mod_alefor
  implicit none
  integer(ip)          :: pblty,iimbo,kauxi
  integer(ip)          :: isets,kpoin,kboun,inodb,iboun,ipoin,jpoin,iperi

  if( IMASTER ) then
     !
     ! No boundaries nor nodes
     !
     do iimbo = 1,nrbod
        rbbou(iimbo) % nboib = 0
        rbbou(iimbo) % npoib = 0
     end do
     
  else

     call memgen(1_ip,npoin,0_ip)

     do iimbo = 1,nrbod

        igene                = iimbo   ! se usaba para que le llegue a ibm_memall ! ahoar que no hay iimbo no tiene mayor interes
        mnoib                = max(mnoib,mnodb)
        mgaib                = max(mgaib,mgaub)
        rbbou(iimbo) % nboib = 0
        rbbou(iimbo) % npoib = 0
        !
        ! Count number of boundaries and mark nodes
        !
        kboun = 0
        do iboun = 1,nboun
           kauxi = 0
           do isets = 1,rbbou(iimbo) % nrbse
              if (lbset(iboun) == rbbou(iimbo) % lrbse(isets) ) kauxi = 1 
           end do
           if( kauxi == 1 ) then
              lexib(ltypb(iboun)) = 2
              kboun = kboun + 1
              do inodb = 1,nnode(ltypb(iboun))
                 ipoin = lnodb(inodb,iboun)
                 gisca(ipoin) = 1
              end do
           end if
        end do
        rbbou(iimbo) % nboib = kboun
        !
        ! Count number of nodes
        !
        ! Problemas en el loop de boundaries, solucion temporal
        !
        do isets = 1,rbbou(iimbo) % nrbse
           do ipoin = 1,npoin
              if ( kfl_codno(1_ip,ipoin) == rbbou(iimbo) % lrbse(isets) .or. kfl_codno(2_ip,ipoin) == rbbou(iimbo) % lrbse(isets) &
                   .or. kfl_codno(3_ip,ipoin) == rbbou(iimbo) % lrbse(isets) ) then
                 gisca(ipoin) = 1_ip
              endif
           end do
        end do
        !
        ! Special check for periodic nodes
        !
        do iperi = 1,nperi
           ipoin = lperi(1,iperi)
           jpoin = lperi(2,iperi)
           if ( ipoin > 0 .and. jpoin > 0 .and. kfl_codno(1,ipoin) == rbbou(1_ip) % lrbse(1_ip))  then
              gisca(ipoin) = 1_ip
           end if
        end do 
        !
        ! Final mark of rigid body nodes
        !
        kpoin = 0
        do ipoin = 1,npoin
           if( gisca(ipoin) == 1 ) then
              kpoin = kpoin + 1
              gisca(ipoin) = kpoin
           end if
        end do
       
        rbbou(iimbo) % npoib = kpoin
        !
        ! Allocate memory
        !
        call alefor_memory_allocate('RBBOU TOPOLOGY',NUMBER1=iimbo)
        !
        ! Permutation with mesh nu,bering
        !
        do ipoin = 1,npoin
           if( gisca(ipoin) /= 0 ) then
              kpoin = gisca(ipoin)
              rbbou(iimbo) % lninv(kpoin) = ipoin
           end if
        end do
        !
        ! Boundaries
        !
        kboun = 0
        do iboun = 1,nboun
           kauxi = 0
           do isets = 1,rbbou(iimbo) % nrbse
              if (lbset(iboun) == rbbou(iimbo) % lrbse(isets) ) kauxi = 1 
           end do
           if( kauxi == 1 ) then
              kboun = kboun + 1
              rbbou(iimbo) % ltyib(kboun) = ltypb(iboun)
              rbbou(iimbo) % lbinv(kboun) = iboun
              do inodb = 1,nnode(ltypb(iboun))
                 rbbou(iimbo) % lnoib(inodb,kboun) = gisca(lnodb(inodb,iboun))
              end do
           end if
        end do
        ! 
        ! Put GISCA to 0
        !
        do ipoin = 1,npoin
           gisca(ipoin) = 0
        end do
        !
        ! Gauss points and integration rule
        !
        do pblty = ibsta_dom,ibsto_dom
           if( lexib(pblty) == 2 ) then
              ngaib(pblty) = ngaus(pblty)
              lruib(pblty) = lrule(pblty)
           end if
        end do
     end do
     !
     ! Shape functions
     !
     call ale_cshder(2_ip)
     do pblty = ibsta_dom,ibsto_dom
        if( lexib(pblty) == 2 ) lexib(pblty) = 1
     end do
     !
     ! GISCA: Deallocate
     !
     call memgen(3_ip,npoin,0_ip)

  end if

end subroutine ale_bodfit
