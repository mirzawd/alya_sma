!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine ale_smobou()
  !-----------------------------------------------------------------------
  !****f* domain/ale_smobou
  ! NAME
  !    domain
  ! DESCRIPTION
  !    This routines smoothes the boundary mesh.
  ! USED BY
  !    Turnon
  !***
  !-----------------------------------------------------------------------
  use def_master
  use def_parame
  use def_elmtyp
  use def_domain
  use def_alefor
  use mod_memchk
  use mod_bouder
  use mod_communications_global, only : PAR_SUM
  implicit none
  integer(ip)          :: iboun,pnodb,pblty,inodb,ielem
  integer(ip)          :: ii,ipoin,jpoin,idime,iters,izdom
  integer(4)           :: istat
  real(rp)             :: lamda,mu,alpha,asele,ri(3),ro(3)
  real(rp)             :: volui,voluf,gbcod(3),gbcod_tmp(3)
  real(rp)             :: eucta,gbsur
  real(rp)             :: baloc(ndime,ndime)
  real(rp)             :: bocod(ndime,mnodb)
  real(rp)             :: bocod_tmp(ndime,mnodb)
  real(rp)             :: elmat(mnodb,mnodb)
  real(rp)             :: elrhs(mnodb)
  real(rp),    pointer :: amatr_tmp(:)
  real(rp),    pointer :: coord_tmp(:,:)
  real(rp),    pointer :: rhsid_tmp(:,:)
  real(rp),    pointer :: unkno_tmp(:,:)
  real(rp),    pointer :: invdiag(:)

  !-------------------------------------------------------------------
  !
  ! Change fixity
  !
  !-------------------------------------------------------------------

  if( INOTMASTER ) then
     do ipoin = 1,npoin
        if( lpoty(ipoin) == 0 ) then
           do idime = 1,ndime
              kfl_fixno_ale(idime,ipoin) = -1
           end do
        else
           do idime = 1,ndime
              kfl_fixno_ale(idime,ipoin) =  0
           end do
        end if
     end do
     call memgen(1_ip,npoin,0_ip)
     do iboun = 1,nboun
        if( lboch(iboun) == BOEXT ) then
           do inodb = 1,nnode(abs(ltypb(iboun)))           
              ipoin = lnodb(inodb,iboun)
              gisca(ipoin) = 1
           end do
        end if
     end do
     call parari('SLX',NPOIN_TYPE,npoin,0_ip)
     do ipoin = 1,npoin
        if( gisca(ipoin) /= 0 ) then
           do idime = 1,ndime
              kfl_fixno_ale(idime,ipoin) = 0
           end do
        end if
     end do
     call memgen(3_ip,npoin,0_ip)
  end if

  !----------------------------------------------------------------------
  !
  ! Detect sharp edges
  !
  !----------------------------------------------------------------------

  if( ansmo_ale > 0.0_rp ) call ale_sharpe()

  !-------------------------------------------------------------------
  !
  ! Allocate memory
  !
  !-------------------------------------------------------------------

  if( INOTMASTER ) then

     allocate( coord_tmp(ndime,npoin) , stat = istat ) 
     call memchk(zero,istat,memor_dom,'COORD_TMP','ale_smobou',coord_tmp) 

     allocate( rhsid_tmp(ndime,npoin) , stat = istat )
     call memchk(zero,istat,memor_dom,'RHSID_TMP','ale_smobou',rhsid_tmp) 

     allocate( amatr_tmp(solve_sol(1)%nzmat) , stat = istat )
     call memchk(zero,istat,memor_dom,'AMATR_TMP','ale_smobou',amatr_tmp) 

     allocate( unkno_tmp(ndime,npoin)  , stat = istat )
     call memchk(zero,istat,memor_dom,'UNKNO_TMP','ale_smobou',unkno_tmp) 

     allocate( invdiag(npoin)          , stat = istat )
     call memchk(zero,istat,memor_dom,'INVDIAG','ale_smobou',invdiag) 

  end if

  !-------------------------------------------------------------------
  !
  ! Assemble matrix
  !
  !-------------------------------------------------------------------

  if( INOTMASTER ) then

     do ipoin = 1,npoin
        do idime = 1,ndime
           rhsid_tmp(idime,ipoin) = 0.0_rp
           coord_tmp(idime,ipoin) = coord(idime,ipoin)
        end do
     end do
     do izdom = 1,nzdom
        amatr_tmp(izdom) = 0.0_rp
     end do

     do iboun = 1,nboun
        pblty = ltypb(iboun)
        pnodb = nnode(pblty)
        ielem = lelbo(iboun)
        do inodb = 1,pnodb
           ipoin = lnodb(inodb,iboun)
           do idime = 1,ndime
              bocod(idime,inodb) = coord(idime,ipoin)
           end do
        end do
        asele = -1.0_rp
        call bouder(&
             pnodb,ndime,ndimb,elmar(pblty)%dercg,&    ! Cartesian derivative
             bocod,baloc,eucta)                        ! and Jacobian
        !
        ! Assemble equation
        !
        call ale_elmpno(&
             pnodb,pblty,lnodb(1,iboun),lelch(ielem),bocod,&
             asele,elmat,elrhs,rhsid_tmp) 
        call assmat(&
             1_ip,pnodb,pnodb,npoin,solve_sol(1) % kfl_algso,&
             -iboun,lnodb(1,iboun),elmat,amatr_tmp) 
     end do

     !-------------------------------------------------------------------
     !
     ! Prescribe boundary conditions
     !
     !-------------------------------------------------------------------

     do ipoin = 1,npoin
        if( lpoty(ipoin) /= 0 .and. kfl_fixno_ale(1,ipoin) == 1 ) then
           do izdom = r_dom(ipoin),r_dom(ipoin+1)-1
              jpoin = c_dom(izdom)
              if( ipoin == jpoin ) then
                 amatr_tmp(izdom) = 1.0_rp
              else
                 amatr_tmp(izdom) = 0.0_rp
              end if
           end do
        end if
     end do

     !-------------------------------------------------------------------
     !
     ! Diagonal and matrix scaling
     !
     !-------------------------------------------------------------------

     do ipoin = 1,npoin
        if( lpoty(ipoin) /= 0 ) then
           izdom = r_dom(ipoin)
           jpoin = c_dom(izdom)
           do while( jpoin /= ipoin )
              izdom = izdom + 1
              jpoin = c_dom(izdom)
           end do
           invdiag(ipoin) = amatr_tmp(izdom)
        else
           invdiag(ipoin) = 0.0_rp
        end if
     end do
     call pararr('SLX',NPOIN_TYPE,npoin,invdiag)
     call pararr('SLX',NPOIN_TYPE,ndime*npoin,rhsid_tmp)
     do ipoin = 1,npoin
        if( lpoty(ipoin) /= 0 ) invdiag(ipoin) = 1.0_rp / invdiag(ipoin)
     end do

     do ipoin = 1,npoin
        if( lpoty(ipoin) /= 0 ) then
           do izdom = r_dom(ipoin),r_dom(ipoin+1)-1
              jpoin = c_dom(izdom)
              if( ipoin == jpoin ) then
                 amatr_tmp(izdom) = invdiag(ipoin) * amatr_tmp(izdom) 
              else
                 amatr_tmp(izdom) = invdiag(ipoin) * amatr_tmp(izdom) 
              end if
           end do
        end if
     end do

     !-------------------------------------------------------------------
     !
     ! Laplacian smoothing using Gauss-Seidel
     !
     !-------------------------------------------------------------------

     lamda = 0.6307_rp
     mu    = 1.0_rp / ( 0.1_rp - 1.0_rp / lamda )

     do iters = 1,nsmob_ale
        do ii = 1,2
           if( ii == 1 ) then
              alpha = lamda
           else
              alpha = mu
           end if

           do ipoin = 1,npoi1
              if( kfl_fixno_ale(1,ipoin) == 0 .and. lpoty(ipoin) /= 0 ) then
                 do idime = 1,ndime
                    ri(idime) = coord_tmp(idime,ipoin) + invdiag(ipoin) * rhsid_tmp(idime,ipoin)
                    ro(idime) = coord_tmp(idime,ipoin)
                 end do
                 do izdom = r_dom(ipoin),r_dom(ipoin+1)-1
                    jpoin = c_dom(izdom)
                    if( lpoty(jpoin) /= 0 ) then
                       do idime = 1,ndime
                          ri(idime) = ri(idime) - amatr_tmp(izdom) * coord_tmp(idime,jpoin) 
                       end do
                    end if
                 end do

                 do idime = 1,ndime
                    coord_tmp(idime,ipoin) = (1.0_rp-alpha) * ro(idime) + alpha * ri(idime)
                 end do
              end if
           end do
           !
           ! Update boundary nodes
           !
           do ipoin = npoi1+1,npoin
              if( kfl_fixno_ale(1,ipoin) == 0 .and. lpoty(ipoin) /= 0 ) then
                 do idime = 1,ndime
                    unkno_tmp(idime,ipoin) = 0.0_rp
                 end do
                 do izdom = r_dom(ipoin),r_dom(ipoin+1)-1
                    jpoin = c_dom(izdom)
                    if( lpoty(jpoin) /= 0 ) then
                       do idime = 1,ndime
                          unkno_tmp(idime,ipoin) = unkno_tmp(idime,ipoin) - amatr_tmp(izdom) * coord_tmp(idime,jpoin) 
                       end do
                    end if
                 end do
              end if
           end do

           call pararr('SLX',NPOIN_TYPE,ndime*npoin,unkno_tmp)
           do ipoin = npoi1+1,npoin
              if( kfl_fixno_ale(1,ipoin) == 0 .and. lpoty(ipoin) /= 0 ) then
                 do idime = 1,ndime
                    unkno_tmp(idime,ipoin) = unkno_tmp(idime,ipoin) &
                         + coord_tmp(idime,ipoin) + invdiag(ipoin) * rhsid_tmp(idime,ipoin)
                    coord_tmp(idime,ipoin) = (1.0_rp-alpha) * coord_tmp(idime,ipoin) + alpha * unkno_tmp(idime,ipoin)
                 end do
              end if
           end do
        end do
     end do

     !-------------------------------------------------------------------
     !
     ! Compute displacement and change fixity for volume deformation
     !
     !-------------------------------------------------------------------
     do ipoin = 1,npoin
        if( lpoty(ipoin) == 0 ) then
           kfl_fixno_ale(1,ipoin) = 0
        else        
           kfl_fixno_ale(1,ipoin) = max(kfl_fixno_ale(1,ipoin),1_ip)
           do idime = 1,ndime
              kfl_fixno_ale(idime,ipoin) = kfl_fixno_ale(1,ipoin) 
              bvess_ale(idime,ipoin,1)   = coord_tmp(idime,ipoin) - coord(idime,ipoin)
           end do
        end if
     end do
     !-------------------------------------------------------------------
     !
     ! Compute initial and final volumes
     !
     !-------------------------------------------------------------------

     volui = 0.0_rp
     voluf = 0.0_rp
     do iboun = 1,nboun
        pblty = ltypb(iboun)
        pnodb = nnode(pblty)
        do inodb = 1,pnodb
           ipoin = lnodb(inodb,iboun)
           do idime = 1,ndime
              bocod(idime,inodb)     = coord(idime,ipoin)
              bocod_tmp(idime,inodb) = coord_tmp(idime,ipoin)
           end do
        end do
        call bouder(&
             pnodb,ndime,ndimb,elmar(pblty)%dercg,&    ! Cartesian derivative
             bocod,baloc,eucta)                        ! and Jacobian
        gbsur = elmar(pblty) % weicg * eucta 
        do idime = 1,ndime 
           gbcod(idime)     = 0.0_rp
           gbcod_tmp(idime) = 0.0_rp
           do inodb = 1,pnodb
              gbcod(idime)     = gbcod(idime)     + bocod(idime,inodb)     * elmar(pblty) % shacg(inodb)
              gbcod_tmp(idime) = gbcod_tmp(idime) + bocod_tmp(idime,inodb) * elmar(pblty) % shacg(inodb)
           end do
        end do
        do idime = 1,ndime
           volui = volui + gbsur * baloc(idime,ndime) * gbcod(idime)     
           voluf = voluf + gbsur * baloc(idime,ndime) * gbcod_tmp(idime) 
        end do

     end do

  end if

  !-------------------------------------------------------------------
  !
  ! Deallocate memory
  !
  !-------------------------------------------------------------------

  if( INOTMASTER ) then

     deallocate( amatr_tmp , stat = istat ) 
     if(istat/=0) call memerr(two,'AMATR_TMP','ale_smobou',0_ip)
     call memchk(two,istat,memor_dom,'amatr_tmp','ale_smobou',amatr_tmp) 

     deallocate( rhsid_tmp , stat = istat ) 
     if(istat/=0) call memerr(two,'RHSID_TMP','ale_smobou',0_ip)
     call memchk(two,istat,memor_dom,'rhsid_tmp','ale_smobou',rhsid_tmp) 

     deallocate( coord_tmp , stat = istat ) 
     if(istat/=0) call memerr(two,'COORD_TMP','ale_smobou',0_ip)
     call memchk(two,istat,memor_dom,'coord_tmp','ale_smobou',coord_tmp) 

     deallocate( invdiag , stat = istat ) 
     if(istat/=0) call memerr(two,'INVDIAG','ale_smobou',0_ip)
     call memchk(two,istat,memor_dom,'invdiag','ale_smobou',invdiag) 

     deallocate( unkno_tmp , stat = istat ) 
     if(istat/=0) call memerr(two,'UNKNO_TMP','ale_smobou',0_ip)
     call memchk(two,istat,memor_dom,'unkno_tmp','ale_smobou',unkno_tmp) 

  end if

  !-------------------------------------------------------------------
  !
  ! Output information
  !
  !-------------------------------------------------------------------

  call PAR_SUM(voluf) 
  call PAR_SUM(volui) 
  volui = volui / real(ndime,rp)
  voluf = voluf / real(ndime,rp)

  if( INOTSLAVE ) then
     write(momod(modul) % lun_outpu,'(a,1x,e12.6)') 'Initial volume= ',volui
     write(momod(modul) % lun_outpu,'(a,1x,e12.6)') 'Final   volume= ',voluf
  end if

end subroutine ale_smobou
