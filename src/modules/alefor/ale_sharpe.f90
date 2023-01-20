!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine ale_sharpe()
  !-----------------------------------------------------------------------
  !****f* alefor/ale_sharpe
  ! NAME 
  !    ale_sharpe
  ! DESCRIPTION
  !    This routine detect sharp edges
  ! USES
  !    ecoute
  !    memchk
  !    runend
  ! USED BY
  !    opebcs
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_alefor
  use mod_memchk
  use mod_outfor,        only : outfor
  use mod_domain,        only : domain_memory_deallocate
  use mod_communications, only : PAR_MAX
  implicit none
  integer(ip), pointer  :: ledbo(:,:)
  real(rp),    pointer  :: bouno(:,:)
  integer(ip)           :: inodb,ipoin,jnodb,jpoin,ista0,ista1
  integer(ip)           :: iboun,icolu,pnodb,ielem,kpoin,idime
  integer(ip)           :: ibopo,jboun,kboun,ierro,ipoi1,ipoi2
  integer(ip)           :: ibou1,nboun_tmp,npoin_tmp,dummi
  integer(4)            :: istat
  real(rp)              :: vnorm,cosim,extno,exwor(2,2)

  ierro = 0

  if( INOTMASTER ) then 
     !
     ! Save geometry dimensions: we are going to consider fringe geometry
     !
     nboun_tmp = nboun
     npoin_tmp = npoin
     !
     ! Boundary graph required C_BOU, R_BOU: must include fringe geometry
     !
     nboun = nboun_2
     npoin = npoin_2
     call domain_memory_deallocate('R_BOU AND C_BOU')
     call bougra()
     nboun = nboun_tmp
     npoin = npoin_tmp     
     !
     ! Allocate temporal memory
     !
     allocate(ledbo(2,nzbou),stat=istat)
     call memchk(zero,istat,memor_dom,'LEDBO','ale_sharpe',ledbo)
     allocate(bouno(ndime,max(1_ip,nboun_2)),stat=istat)
     call memchk(zero,istat,memor_dom,'BOUNO','ale_sharpe',bouno)

     !----------------------------------------------------------------------
     !
     ! Compute boundary edges IEDGE
     ! LEDBO(1,iedge) = boundary 1
     ! LEDBO(2,iedge) = boundary 2 / 0
     !
     !----------------------------------------------------------------------

     do iboun = 1,nboun_2
        pnodb = nnode(ltypb(iboun))
        ielem = lelbo(iboun)
        if( ielem >= 1 ) then  ! Volume element
           do inodb = 1,pnodb
              ipoin = lnodb(inodb,iboun)
              do jnodb = 1,pnodb
                 if( inodb /= jnodb ) then
                    jpoin = lnodb(jnodb,iboun)                
                    ista0 = r_bou(ipoin  )-1
                    ista1 = r_bou(ipoin+1)-1
                    icolu = ista0
                    do while( icolu /= ista1 )
                       icolu = icolu + 1
                       if( c_bou(icolu) == jpoin ) then
                          if( ledbo(1,icolu) == 0 ) then
                             ledbo(1,icolu) = iboun
                          else if( ledbo(2,icolu) == 0 ) then
                             ledbo(2,icolu) = iboun
                          else
                             ierro = ierro + 1
                             ipoi1 = ipoin
                             ipoi2 = jpoin
                             ibou1 = iboun
                             print*,'error in ale_sharpe'
                             print*,ipoin,jpoin,ledbo(1,icolu),ledbo(2,icolu),iboun
                             stop
                          end if
                          icolu = ista1
                       end if
                    end do
                 end if
              end do
           end do
        end if
     end do
  end if

  !----------------------------------------------------------------------
  !
  ! Check errors
  !
  !----------------------------------------------------------------------

  call PAR_MAX(ierro)
  if( ierro /= 0 ) then
     ioutp(1) = ierro
     ioutp(2) = ipoi1
     ioutp(3) = ipoi2
     ioutp(4) = ibou1
     call outfor(45_ip,lun_outpu,' ')
  end if

  if( INOTMASTER ) then

     !----------------------------------------------------------------------
     !
     ! Initialization
     !
     !----------------------------------------------------------------------

     !
     ! Compute the normals to the boundaries
     !
     if( nboun_2 > 0 ) then
        call bounor(nboun_2,lnodb,ltypb,lelbo,dummi,bouno)
     end if
     !
     ! Maximum angle and initialization
     !
     cosim = cos(ansmo_ale * pi / 180.0_rp)

     if( ndime == 2 ) then

        !-------------------------------------------------------------------    
        !
        ! 2D
        !
        !-------------------------------------------------------------------

        do ipoin = 1,npoin
           ibopo = lpoty(ipoin)
           kpoin = 0
           if( ibopo /= 0 ) then
              !
              ! Check if point is on a slip-wall or wall-law boundary
              !
              exwor = 0.0_rp

              do iboun = 1,nboun_2
                 pnodb = nnode(ltypb(iboun))
                 nodes_2d: do inodb = 1,pnodb
                    if( lnodb(inodb,iboun) == ipoin ) then
                       kpoin                = kpoin + 1
                       exwor(1:ndime,kpoin) = bouno(1:ndime,iboun)
                       exit nodes_2d
                    end if
                 end do nodes_2d
              end do
              !
              ! Check if point is on an obtuse angle or if it is a corner
              !
              if( kpoin == 2 ) then
                 vnorm = dot_product(exwor(1:ndime,1),exwor(1:ndime,2))
                 exwor(1:ndime,1) = exwor(1:ndime,1) + exwor(1:ndime,2)
                 call vecuni(ndime,exwor,extno)
                 if( vnorm < cosim ) then    
                    do idime = 1,ndime
                       kfl_fixno_ale(idime,ipoin) = 2
                    end do
                 end if
              end if
           end if
        end do

     else

        !-------------------------------------------------------------------
        !
        ! 3D
        !
        !-------------------------------------------------------------------

        do ipoin = 1,npoin

           if( lpoty(ipoin) /= 0 ) then

              loop_segment_ipoin_jpoin: do icolu = r_bou(ipoin),r_bou(ipoin+1)-1
                 !
                 ! Angle between jboun-iboun
                 !
                 jboun = ledbo(1,icolu)
                 kboun = ledbo(2,icolu)
                 if( jboun /= 0 .and. kboun /= 0 ) then
                    vnorm = dot_product(bouno(1:ndime,jboun),bouno(1:ndime,kboun))
                    if( vnorm < cosim ) then  
                       do idime = 1,ndime
                          kfl_fixno_ale(idime,ipoin) = 2
                       end do
                       exit loop_segment_ipoin_jpoin
                    end if
                 end if
                 
              end do loop_segment_ipoin_jpoin

           end if
        end do

     end if
     !
     ! Deallocate volatile memory
     !
     call domain_memory_deallocate('R_BOU AND C_BOU')

     call memchk(two,istat,memor_dom,'LEDBO','geonor',ledbo)
     deallocate(ledbo,stat=istat)
     if(istat/=0)  call memerr(two,'LEDBO','geonor',0_ip)

     call memchk(two,istat,memor_dom,'BOUNO','geonor',bouno)
     deallocate(bouno,stat=istat)
     if(istat/=0)  call memerr(two,'BOUNO','geonor',0_ip)

  end if

end subroutine ale_sharpe
