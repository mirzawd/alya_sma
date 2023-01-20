!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine par_mygather()
  !-----------------------------------------------------------------------
  !****f* Parall/par_mygather
  ! NAME
  !    par_mygather
  ! DESCRIPTION
  !    This subroutine exchange arrays between master and slaves
  ! USED BY
  !    Parall
  !***
  !-----------------------------------------------------------------------
  use def_master
  use def_parame
  use def_parall
  use def_domain
  use mod_memchk
  use mod_parall, only : par_memor
  implicit none
  integer(ip)          :: ielem,iboun,ipoin,kdime
  integer(ip)          :: jelem,jboun,jpoin,idofn
!  integer(ip)          :: ibopo,jbopo
  integer(ip)          :: ii,ipar1,ipar2,jdofn
  integer(4)           :: istat
  integer(ip), pointer :: loc_pari1(:),loc_pari2(:,:)
  real(rp),    pointer :: loc_parr1(:),loc_parr2(:,:),loc_parr3(:,:,:)
  !
  ! Initialize
  !
  npari = 0
  nparr = 0
  nparc = 0
  if( ISLAVE ) kfl_desti_par = 0

  if( party == 1 ) then

     !-------------------------------------------------------------------
     !
     ! Element: NELEM
     !
     !-------------------------------------------------------------------

     if( pardi == 1 .and. parki == 1 ) then
        !
        ! INTEGER(NELEM)
        !
        if( IMASTER .and. .not. READ_AND_RUN() ) then   
           allocate(loc_pari1(nelem),stat=istat)
           call memchk(zero,istat,par_memor,'LOC_PARI1','par_mygather',loc_pari1)
           do ielem = 1, nelem
              jelem = leinv_par(ielem)
              loc_pari1(ielem) = pari1(jelem)
           enddo
           ii = 1
           do kfl_desti_par = 1,npart_par
              kdime = nelem_par(kfl_desti_par)
              if( kdime > 0 ) call par_srinte(1_ip,kdime,loc_pari1(ii:))
              ii = ii + nelem_par(kfl_desti_par)
           enddo
           call memchk(two,istat,par_memor,'LOC_PARI1','par_mygather',loc_pari1)
           deallocate(loc_pari1,stat=istat)
           if(istat/=0) call memerr(two,'par_mygather','LOC_PARI1',0_ip)           

        else if( ISLAVE ) then 

           kdime = nelem
           if( kdime > 0 ) call par_srinte(2_ip,kdime,pari1)

        end if

     else if( pardi == 2 .and. parki == 1 ) then
        !
        ! INTEGER(KDIME,NELEM)
        !
        if( IMASTER .and. .not. READ_AND_RUN() ) then

           allocate(loc_pari2(pard1,nelem),stat=istat)
           call memchk(zero,istat,par_memor,'LOC_PARI2','par_mygather',loc_pari2)
           do ielem = 1, nelem
              jelem = leinv_par(ielem)
              do ipar1 = 1,pard1
                 loc_pari2(ipar1,ielem) = pari2(ipar1,jelem)
              end do
           end do
           ii = 1
           do kfl_desti_par = 1, npart_par
              kdime = nelem_par(kfl_desti_par)*pard1
              if( kdime > 0 ) call par_srinte(1_ip,kdime,loc_pari2(1,ii:))
              ii = ii + nelem_par(kfl_desti_par)
           end do
           call memchk(two,istat,par_memor,'LOC_PARI2','par_mygather',loc_pari2)
           deallocate(loc_pari2,stat=istat)
           if(istat/=0) call memerr(two,'LOC_PARI2','par_mygather',0_ip)

        else if( ISLAVE ) then

           kdime = pard1*nelem
           if( kdime > 0 ) call par_srinte(2_ip,kdime,pari2)

        end if        

      else if( pardi == 1 .and. parki == 2 ) then
        !
        ! REAL(NELEM)
        !
        if( IMASTER .and. .not. READ_AND_RUN() ) then   
           allocate(loc_parr1(nelem),stat=istat)
           call memchk(zero,istat,par_memor,'LOC_PARR1','par_mygather',loc_parr1)
           do ielem = 1, nelem
              jelem = leinv_par(ielem)
              loc_parr1(ielem) = parr1(jelem)
           enddo
           ii = 1
           do kfl_desti_par = 1,npart_par
              kdime = nelem_par(kfl_desti_par)
              if( kdime > 0 ) call par_srreal(1_ip,kdime,loc_parr1(ii))
              ii = ii + nelem_par(kfl_desti_par)
           enddo
           call memchk(two,istat,par_memor,'LOC_PARR1','par_mygather',loc_parr1)
           deallocate(loc_parr1,stat=istat)
           if(istat/=0) call memerr(two,'par_mygather','LOC_PARR1',0_ip)           

        else if( ISLAVE ) then 

           kdime = nelem
           if( kdime > 0 ) call par_srreal(2_ip,kdime,parr1)

        end if

    else if( pardi == 1 .and. parki == 6 ) then
        !
        ! REAL(PARD1*NELEM)
        !
        if( IMASTER .and. .not. READ_AND_RUN() ) then
           allocate(loc_parr1(pard1*nelem),stat=istat)
           call memchk(zero,istat,par_memor,'LOC_PARR1','par_mygather',loc_parr1)
           do ielem = 1, nelem
              jelem = leinv_par(ielem)
              idofn = ( ielem - 1 ) * pard1
              jdofn = ( jelem - 1 ) * pard1
              do ipar1 = 1,pard1
                 idofn = idofn + 1
                 jdofn = jdofn + 1
                 loc_parr1(idofn) = parr1(jdofn)
              end do 
           end do
           ii = 1
           do kfl_desti_par = 1,npart_par
              kdime = pard1 * nelem_par(kfl_desti_par)
              if( kdime > 0 ) call par_srreal(1_ip,kdime,loc_parr1(ii:))
              ii = ii + kdime
           end do
           call memchk(two,istat,par_memor,'LOC_PARR1','par_mygather',loc_parr1)
           deallocate(loc_parr1,stat=istat)
           if(istat/=0) call memerr(two,'LOC_PARR1','par_mygather',0_ip)

        else if( ISLAVE ) then

           kdime = pard1*nelem
           if( kdime > 0 ) call par_srreal(2_ip,kdime,parr1)

        end if

     else if( pardi == 3 .and. parki == 2 ) then
        !
        ! REAL(KDIME,NDIM2,NELEM)
        !
        if( IMASTER .and. .not. READ_AND_RUN() ) then

           allocate(loc_parr3(pard1,pard2,nelem),stat=istat)
           call memchk(zero,istat,par_memor,'LOC_PARR3','par_mygather',loc_parr3)
           do ielem = 1,nelem
              jelem = leinv_par(ielem)
              do ipar2 = 1,pard2
                 do ipar1 = 1,pard1
                    loc_parr3(ipar1,ipar2,ielem) = parr3(ipar1,ipar2,jelem)
                 end do
              end do
           end do
           ii = 1
           do kfl_desti_par = 1, npart_par
              kdime = nelem_par(kfl_desti_par)*pard1*pard2
              if( kdime > 0 ) call par_srreal(1_ip,kdime,loc_parr3(1,1,ii))
              ii = ii + nelem_par(kfl_desti_par)
           end do
           call memchk(two,istat,par_memor,'LOC_PARR3','par_mygather',loc_parr3)
           deallocate(loc_parr3,stat=istat)
           if(istat/=0) call memerr(two,'LOC_PARR3','par_mygather',0_ip)

        else if( ISLAVE ) then

           kdime = pard1*pard2*nelem
           if( kdime > 0 ) call par_srreal(2_ip,kdime,parr3)

        end if        

     else

        call runend('PAR_MYGATHER: THIS GATHER DOES NOT EXIST')

     end if

  else if( party == 2 ) then

     !-------------------------------------------------------------------
     !
     ! Boundary: NBOUN
     !
     !-------------------------------------------------------------------

     if( pardi == 1 .and. parki == 1 ) then
        !
        ! INTEGER(NBOUN)
        !
        if( IMASTER .and. .not. READ_AND_RUN() ) then     
  
           allocate(loc_pari1(nboun),stat=istat)
           call memchk(zero,istat,par_memor,'LOC_PARI1','par_mygather',loc_pari1)
           do iboun = 1, nboun
              jboun = lbper_par(iboun)
              loc_pari1(jboun) = pari1(iboun)
           end do
           ii = 1
           do kfl_desti_par = 1, npart_par
              kdime =  nboun_par(kfl_desti_par)
              if( kdime > 0 ) call par_srinte(1_ip,kdime,loc_pari1(ii))
              ii = ii + nboun_par(kfl_desti_par)
           end do
           call memchk(two,istat,par_memor,'LOC_PARI1','par_mygather',loc_pari1)
           deallocate(loc_pari1,stat=istat)
           if(istat/=0) call memerr(two,'par_mygather','LOC_PARI1',0_ip)     
      
        else if( ISLAVE ) then 

           kdime = nboun
           if( kdime > 0 ) call par_srinte(2_ip,kdime,pari1)

        end if

     else if( pardi == 1 .and. parki == 2 ) then
        !
        ! REAL(NBOUN)
        !
        if( IMASTER .and. .not. READ_AND_RUN() ) then     
  
           allocate(loc_parr1(nboun),stat=istat)
           call memchk(zero,istat,par_memor,'LOC_PARR1','par_mygather',loc_parr1)
           do iboun = 1, nboun
              jboun = lbper_par(iboun)
              loc_parr1(jboun) = parr1(iboun)
           end do
           ii = 1
           do kfl_desti_par = 1, npart_par
              kdime = nboun_par(kfl_desti_par)
              if( kdime > 0 ) call par_srreal(1_ip,kdime,loc_parr1(ii))
              ii = ii + nboun_par(kfl_desti_par)
           end do
           call memchk(two,istat,par_memor,'LOC_PARR1','par_mygather',loc_parr1)
           deallocate(loc_parr1,stat=istat)
           if(istat/=0) call memerr(two,'par_mygather','LOC_PARR1',0_ip)    
       
        else if( ISLAVE ) then 

           kdime = nboun
           if( kdime > 0 ) call par_srreal(2_ip,kdime,parr1)

        end if

     else if( pardi == 2 .and. parki == 2 ) then
        !
        ! REAL(:,NBOUN)
        !
        if( IMASTER .and. .not. READ_AND_RUN() ) then

           allocate(loc_parr2(pard1,nboun),stat=istat)
           call memchk(zero,istat,par_memor,'LOC_PARR2','par_mygather',loc_parr2)
           do iboun = 1, nboun
              jboun = lbper_par(iboun)
              do ipar1=1,pard1
                 loc_parr2(ipar1,jboun) = parr2(ipar1,iboun)
              end do
           end do

           ii = 1
           do kfl_desti_par = 1, npart_par
              kdime = nboun_par(kfl_desti_par)*pard1
              if( kdime > 0 ) call par_srreal(1_ip,kdime,loc_parr2(:,ii))
              ii = ii + nboun_par(kfl_desti_par)
           end do
           call memchk(two,istat,par_memor,'LOC_PARR2','par_mygather',loc_parr2)
           deallocate(loc_parr2,stat=istat)
           if(istat/=0) call memerr(two,'par_mygather','LOC_PARR2',0_ip)    
       
        else if( ISLAVE ) then

           kdime = pard1*nboun
           if( kdime > 0 ) call par_srreal(2_ip,kdime,parr2)

        end if

     else if( pardi == 2 .and. parki == 1 ) then
        !
        ! INTEGER(:,NBOUN)
        !
        if( IMASTER .and. .not. READ_AND_RUN() ) then
           allocate(loc_pari2(pard1,nboun),stat=istat)
           call memchk(zero,istat,par_memor,'LOC_PARI2','par_mygather',loc_pari2)
           do iboun = 1, nboun
              jboun = lbper_par(iboun)
              do ipar1=1,pard1
                 loc_pari2(ipar1,jboun) = pari2(ipar1,iboun)
              end do
           end do

           ii = 1
           do kfl_desti_par = 1, npart_par
              kdime = nboun_par(kfl_desti_par)*pard1
              if( kdime > 0 ) call par_srinte(1_ip,kdime,loc_pari2(:,ii))
              ii = ii + nboun_par(kfl_desti_par)
           end do
           call memchk(two,istat,par_memor,'LOC_PARI2','par_mygather',loc_pari2)
           deallocate(loc_pari2,stat=istat)
           if(istat/=0) call memerr(two,'par_mygather','LOC_PARI2',0_ip)           

        else if( ISLAVE ) then

           kdime = pard1*nboun
           if( kdime > 0 ) call par_srinte(2_ip,kdime,pari2)

        end if

     else if( pardi == 3 .and. parki == 2 ) then
        !
        ! REAL(:,:,NBOUN)
        !
        if( IMASTER .and. .not. READ_AND_RUN() ) then

           allocate(loc_parr3(pard1,pard2,nboun),stat=istat)
           call memchk(zero,istat,par_memor,'LOC_PARR3','par_mygather',loc_parr3)
           do iboun = 1, nboun
              jboun = lbper_par(iboun)
              do ipar2 = 1,pard2
                 do ipar1 = 1,pard1
                    loc_parr3(ipar1,ipar2,jboun) = parr3(ipar1,ipar2,iboun)
                 end do
              end do
           end do

           ii = 1
           do kfl_desti_par = 1, npart_par
              kdime = nboun_par(kfl_desti_par)*pard1*pard2
              if( kdime > 0 ) call par_srreal(1_ip,kdime,loc_parr3(1,1,ii))
              ii = ii + nboun_par(kfl_desti_par)
           end do
           call memchk(two,istat,par_memor,'LOC_PARR3','par_mygather',loc_parr3)
           deallocate(loc_parr3,stat=istat)
           if(istat/=0) call memerr(two,'par_mygather','LOC_PARR3',0_ip) 
          
        else if( ISLAVE ) then

           kdime = pard1*pard2*nboun
           if( kdime > 0 ) call par_srreal(2_ip,kdime,parr3)

        end if

     else if( pardi == 2 .and. parki == 1 ) then
        !
        ! INTEGER(:,NBOUN)
        !
        if( IMASTER .and. .not. READ_AND_RUN() ) then
           allocate(loc_pari2(pard1,nboun),stat=istat)
           call memchk(zero,istat,par_memor,'LOC_PARI2','par_mygather',loc_pari2)
           do iboun = 1, nboun
              jboun = lbper_par(iboun)
              do ipar1=1,pard1
                 loc_pari2(ipar1,jboun) = pari2(ipar1,iboun)
              end do
           end do

           ii = 1
           do kfl_desti_par = 1, npart_par
              kdime = nboun_par(kfl_desti_par)*pard1
              if( kdime > 0 ) call par_srinte(1_ip,kdime,loc_pari2(1,ii))
              ii = ii + nboun_par(kfl_desti_par)
           end do
           call memchk(two,istat,par_memor,'LOC_PARI2','par_mygather',loc_pari2)
           deallocate(loc_pari2,stat=istat)
           if(istat/=0) call memerr(two,'par_mygather','LOC_PARI2',0_ip)           

        else if( ISLAVE ) then

           kdime = pard1*nboun
           if( kdime > 0 ) call par_srinte(2_ip,kdime,pari2)

        end if

     else if( pardi == 1 .and. parki == 6 ) then
        !
        ! REAL(PARD1*NBOUN)
        !
        if( IMASTER .and. .not. READ_AND_RUN() ) then
           allocate(loc_parr1(pard1*nboun),stat=istat)
           call memchk(zero,istat,par_memor,'LOC_PARR1','par_mygather',loc_parr1)
           do iboun = 1, nboun
              jboun = lbper_par(iboun)
              idofn = (iboun-1) * pard1
              jdofn = (jboun-1) * pard1
              do ipar1 = 1,pard1
                 idofn = idofn + 1
                 jdofn = jdofn + 1
                 loc_parr1(jdofn) = parr1(idofn)
              end do
           end do

           ii = 1
           do kfl_desti_par = 1, npart_par
              kdime = nboun_par(kfl_desti_par)*pard1
              if( kdime > 0 ) call par_srreal(1_ip,kdime,loc_parr1(ii))
              ii = ii + kdime
           end do
           call memchk(two,istat,par_memor,'LOC_PARR1','par_mygather',loc_parr1)
           deallocate(loc_parr1,stat=istat)
           if(istat/=0) call memerr(two,'par_mygather','LOC_PARR1',0_ip)           

        else if( ISLAVE ) then

           kdime = pard1*nboun
           if( kdime > 0 ) call par_srreal(2_ip,kdime,parr1)

        end if

     else

        call runend('PAR_MYGATHER: THIS GATHER DOES NOT EXIST')

     end if

  else if( party == 3 ) then

     !-------------------------------------------------------------------
     !
     ! Node: NPOIN
     !
     !-------------------------------------------------------------------

     if( pardi == 1 .and. parki == 1 ) then
        !
        ! INT(NPOIN)
        !
        if( IMASTER .and. .not. READ_AND_RUN() ) then

           allocate(loc_pari1(npoin_total),stat=istat)
           call memchk(zero,istat,par_memor,'LOC_PARI1','par_mygather',loc_pari1)
           do ipoin = 1, npoin_total
              jpoin = lninv_loc(ipoin)
              loc_pari1(ipoin) = pari1(jpoin)
           end do
           ii = 1
           do kfl_desti_par = 1, npart_par
              kdime = npoin_par(kfl_desti_par)
              if( kdime > 0 ) call par_srinte(1_ip,kdime,loc_pari1(ii))
              ii = ii + npoin_par(kfl_desti_par)
           end do
           call memchk(two,istat,par_memor,'LOC_PARI1','par_mygather',loc_pari1)
           deallocate(loc_pari1,stat=istat)
           if(istat/=0) call memerr(two,'par_mygather','LOC_PARI1',0_ip)

        else if( ISLAVE ) then

           kdime = npoin
           if( kdime > 0 ) call par_srinte(2_ip,kdime,pari1)

        end if

     else if( pardi == 1 .and. parki == 5 ) then
        !
        ! INT(PARD1*NPOIN)
        !
        if( IMASTER .and. .not. READ_AND_RUN() ) then
           allocate(loc_pari1(pard1*npoin_total),stat=istat)
           call memchk(zero,istat,par_memor,'LOC_PARI1','par_mygather',loc_pari1)
           do ipoin = 1, npoin_total
              jpoin = lninv_loc(ipoin)
              idofn = (ipoin-1)*pard1
              jdofn = (jpoin-1)*pard1
              do ipar1 = 1,pard1
                 idofn = idofn + 1
                 jdofn = jdofn + 1
                 loc_pari1(idofn) = pari1(jdofn)
              end do 
           end do
           ii = 1
           do kfl_desti_par = 1, npart_par
              kdime = pard1 * npoin_par(kfl_desti_par)
              if( kdime > 0 ) call par_srinte(1_ip,kdime,loc_pari1(ii))
              ii = ii + kdime
           end do
           call memchk(two,istat,par_memor,'LOC_PARI1','par_mygather',loc_pari1)
           deallocate(loc_pari1,stat=istat)
           if(istat/=0) call memerr(two,'LOC_PARI1','par_mygather',0_ip)

        else if( ISLAVE ) then

           kdime = pard1*npoin
           if( kdime > 0 ) call par_srinte(2_ip,kdime,pari1)

        end if

     else if( pardi == 1 .and. parki == 6 ) then
        !
        ! REAL(PARD1*NPOIN)
        !
        if( IMASTER .and. .not. READ_AND_RUN() ) then
           allocate(loc_parr1(pard1*npoin_total),stat=istat)
           call memchk(zero,istat,par_memor,'LOC_PARR1','par_mygather',loc_parr1)
           do ipoin = 1, npoin_total
              jpoin = lninv_loc(ipoin)
              idofn = (ipoin-1)*pard1
              jdofn = (jpoin-1)*pard1
              do ipar1 = 1,pard1
                 idofn = idofn + 1
                 jdofn = jdofn + 1
                 loc_parr1(idofn) = parr1(jdofn)
              end do 
           end do
           ii = 1
           do kfl_desti_par = 1, npart_par
              kdime = pard1 * npoin_par(kfl_desti_par)
              if( kdime > 0 ) call par_srreal(1_ip,kdime,loc_parr1(ii))
              ii = ii + kdime
           end do
           call memchk(two,istat,par_memor,'LOC_PARR1','par_mygather',loc_parr1)
           deallocate(loc_parr1,stat=istat)
           if(istat/=0) call memerr(two,'LOC_PARR1','par_mygather',0_ip)

        else if( ISLAVE ) then

           kdime = pard1*npoin
           if( kdime > 0 ) call par_srreal(2_ip,kdime,parr1)

        end if

     else if( pardi == 2 .and. parki == 1 ) then
        !
        ! INT(:,NPOIN)
        !
        if( IMASTER .and. .not. READ_AND_RUN() ) then

           allocate(loc_pari2(pard1,npoin_total),stat=istat)
           call memchk(zero,istat,par_memor,'LOC_PARI2','par_mygather',loc_pari2)
           do ipoin = 1, npoin_total
              jpoin = lninv_loc(ipoin)
              do ipar1 = 1,pard1
                 loc_pari2(ipar1,ipoin) = pari2(ipar1,jpoin)
              end do
           end do
           ii = 1
           do kfl_desti_par = 1, npart_par
              kdime = pard1* npoin_par(kfl_desti_par)
              if( kdime > 0 ) call par_srinte(1_ip,kdime,loc_pari2(1,ii))
              ii = ii + npoin_par(kfl_desti_par)
           end do
           call memchk(two,istat,par_memor,'LOC_PARI2','par_mygather',loc_pari2)
           deallocate(loc_pari2,stat=istat)
           if(istat/=0) call memerr(two,'LOC_PARI2','par_mygather',0_ip)

        else if( ISLAVE ) then

           kdime = pard1*npoin
           if( kdime > 0 ) call par_srinte(2_ip,kdime,pari2)

        end if

     else if( pardi == 1 .and. parki == 2 ) then
        !
        ! REAL(NPOIN)
        !
        if( IMASTER .and. .not. READ_AND_RUN() ) then

           allocate(loc_parr1(npoin_total),stat=istat)
           call memchk(zero,istat,par_memor,'LOC_PARR1','par_mygather',loc_parr1)
           do ipoin = 1, npoin_total
              jpoin = lninv_loc(ipoin)
              loc_parr1(ipoin) = parr1(jpoin)
           end do
           ii = 1
           do kfl_desti_par = 1, npart_par
              nparr = npoin_par(kfl_desti_par)
              parre => loc_parr1(ii:)
              call par_sendin()
              ii = ii + npoin_par(kfl_desti_par)
           end do
           call memchk(two,istat,par_memor,'LOC_PARR1','par_mygather',loc_parr1)
           deallocate(loc_parr1,stat=istat)
           if(istat/=0) call memerr(two,'par_mygather','LOC_PARR1',0_ip)

        else if( ISLAVE ) then

           nparr =  npoin
           parre => parr1
           call par_receiv()

        end if

     else if( pardi == 2 .and. parki == 2 ) then
        !
        ! REAL(:,NPOIN)
        !
        if( IMASTER .and. .not. READ_AND_RUN() ) then

           allocate(loc_parr2(pard1,npoin_total),stat=istat)
           call memchk(zero,istat,par_memor,'LOC_PARR2','par_mygather',loc_parr2)
           do ipoin = 1, npoin_total
              jpoin = lninv_loc(ipoin)
              do ipar1 = 1,pard1
                 loc_parr2(ipar1,ipoin) = parr2(ipar1,jpoin)
              end do
           end do
           ii = 1
           do kfl_desti_par = 1, npart_par
              kdime = pard1*npoin_par(kfl_desti_par)
              if( kdime > 0 ) call par_srreal(1_ip,kdime,loc_parr2(1,ii))
              ii = ii + npoin_par(kfl_desti_par)
           end do
           call memchk(two,istat,par_memor,'LOC_PARR2','par_mygather',loc_parr2)
           deallocate(loc_parr2,stat=istat)
           if(istat/=0) call memerr(two,'LOC_PARR2','par_mygather',0_ip)

        else

           kdime = pard1*npoin
           if( kdime > 0 ) call par_srreal(2_ip,kdime,parr2)

        end if


     else if( pardi == 3 .and. parki == 2 ) then
        !
        ! REAL(:,:,NPOIN)
        !
        if( IMASTER .and. .not. READ_AND_RUN() ) then

           allocate(loc_parr3(pard1,pard2,npoin_total),stat=istat)
           call memchk(zero,istat,par_memor,'LOC_PARR3','par_mygather',loc_parr3)
           do ipoin = 1, npoin_total
              jpoin = lninv_loc(ipoin)
              do ipar2 = 1,pard2
                 do ipar1 = 1,pard1
                    loc_parr3(ipar1,ipar2,ipoin) = parr3(ipar1,ipar2,jpoin)
                 end do
              end do
           end do
           ii = 1
           do kfl_desti_par = 1, npart_par
              kdime = pard1*pard2*npoin_par(kfl_desti_par)
              if( kdime > 0 ) call par_srreal(1_ip,kdime,loc_parr3(1,1,ii))
              ii = ii + npoin_par(kfl_desti_par)
           end do
           call memchk(two,istat,par_memor,'LOC_PARR3','par_mygather',loc_parr3)
           deallocate(loc_parr3,stat=istat)
           if(istat/=0) call memerr(two,'LOC_PARR3','par_mygather',0_ip)

        else if( ISLAVE ) then

           kdime = pard1*pard2*npoin
           call par_srreal(2_ip,kdime,parr3)

        end if

     else

        call runend('PAR_MYGATHER: THIS GATHER DOES NOT EXIST')

     end if

  else if( party == 4 ) then

     !-------------------------------------------------------------------
     !
     ! Boundary node: NBOPO
     !
     !-------------------------------------------------------------------

     if( pardi == 1 .and. parki == 1 ) then
        !
        ! INTEGER(NBOPO)
        !
        call runend('PAR_MYGATHER: NO LONGER SUPPORTED')
!!$        if( IMASTER .and. .not. READ_AND_RUN() ) then
!!$           allocate(loc_pari1(gnbop_loc+1),stat=istat)
!!$           call memchk(zero,istat,par_memor,'LOC_PARI1','par_mygather',loc_pari1)
!!$           do ibopo = 1, gnbop_loc
!!$              jbopo = exnpe_loc(ibopo)
!!$              loc_pari1(ibopo) = pari1(jbopo)
!!$           end do
!!$           ii = 1
!!$           do kfl_desti_par = 1, npart_par
!!$              kdime =  nbopo_par(kfl_desti_par)
!!$              if( kdime > 0 ) call par_srinte(1_ip,kdime,loc_pari1(ii))
!!$              ii = ii + nbopo_par(kfl_desti_par)
!!$           end do
!!$           call memchk(two,istat,par_memor,'LOC_PARI1','par_mygather',loc_pari1)
!!$           deallocate(loc_pari1,stat=istat)
!!$           if(istat/=0) call memerr(two,'LOC_PARI1','par_mygather',0_ip)
!!$        else if( ISLAVE ) then
!!$           kdime =  nbopo
!!$           if(kdime > 0 ) call par_srinte(2_ip,kdime,pari1)
!!$        end if

     else if(pardi==3.and.parki==2) then
        !
        ! REAL(:,:,NBOPO)
        !
        call runend('PAR_MYGATHER: NO LONGER SUPPORTED')
!!$        if( IMASTER .and. .not. READ_AND_RUN() ) then 
!!$           allocate(loc_parr3(pard1,pard2,gnbop_loc+1),stat=istat)
!!$           call memchk(zero,istat,par_memor,'LOC_PARR3','par_mygather',loc_parr3)
!!$           do ibopo = 1, gnbop_loc
!!$              jbopo = exnpe_loc(ibopo)
!!$              do ipar2 = 1,pard2
!!$                 do ipar1 = 1,pard1
!!$                    loc_parr3(ipar1,ipar2,ibopo) = parr3(ipar1,ipar2,jbopo)
!!$                 end do
!!$              end do
!!$           end do
!!$           ii = 1
!!$           do kfl_desti_par = 1, npart_par
!!$              kdime = pard1*pard2*nbopo_par(kfl_desti_par)
!!$              if( kdime > 0 ) call par_srreal(1_ip,kdime,loc_parr3(1,1,ii))
!!$              ii = ii + nbopo_par(kfl_desti_par)            
!!$           end do
!!$           call memchk(two,istat,par_memor,'LOC_PARR3','par_mygather',loc_parr3)
!!$           deallocate(loc_parr3,stat=istat)
!!$           if(istat/=0) call memerr(two,'LOC_PARR3','par_mygather',0_ip)
!!$        else if( ISLAVE ) then
!!$           kdime = pard1*pard2*nbopo
!!$           if( kdime > 0 ) call par_srreal(2_ip,kdime,parr3)
!!$        end if

     else if( pardi == 2 .and. parki == 2 ) then
        !
        ! REAL(:,NBOPO)
        !
        call runend('PAR_MYGATHER: NO LONGER SUPPORTED')
!!$        if( IMASTER .and. .not. READ_AND_RUN() ) then
!!$           allocate(loc_parr2(pard1,gnbop_loc+1),stat=istat)
!!$           call memchk(zero,istat,par_memor,'LOC_PARR2','par_mygather',loc_parr2)
!!$           do ibopo = 1, gnbop_loc
!!$              jbopo = exnpe_loc(ibopo)
!!$              do ipar1 = 1,pard1
!!$                 loc_parr2(ipar1,ibopo) = parr2(ipar1,jbopo)
!!$              end do
!!$           end do
!!$           ii = 1
!!$           do kfl_desti_par = 1, npart_par
!!$              kdime =  pard1*nbopo_par(kfl_desti_par)
!!$              if( kdime > 0 ) call par_srreal(1_ip,kdime,loc_parr2(1,ii))
!!$              ii = ii + nbopo_par(kfl_desti_par)            
!!$           end do
!!$           call memchk(two,istat,par_memor,'LOC_PARR2','par_mygather',loc_parr2)
!!$           deallocate(loc_parr2,stat=istat)
!!$           if(istat/=0) call memerr(two,'LOC_PARR2','par_mygather',0_ip)
!!$        else if( ISLAVE ) then
!!$           kdime = pard1*nbopo
!!$           if( kdime > 0 ) call par_srreal(2_ip,kdime,parr2)
!!$        end if

     else if( pardi == 2 .and. parki == 1 ) then
        !
        ! INTEGER(:,NBOPO)
        !
        call runend('PAR_MYGATHER: NO LONGER SUPPORTED')
!!$        if( IMASTER .and. .not. READ_AND_RUN() ) then
!!$           allocate(loc_pari2(pard1,gnbop_loc+1),stat=istat)
!!$           call memchk(zero,istat,par_memor,'LOC_PARI2','par_mygather',loc_pari2)
!!$           do ibopo = 1, gnbop_loc
!!$              jbopo = exnpe_loc(ibopo)
!!$              do ipar1 = 1,pard1
!!$                 loc_pari2(ipar1,ibopo) = pari2(ipar1,jbopo)
!!$              end do
!!$           end do
!!$           ii = 1
!!$           do kfl_desti_par = 1, npart_par
!!$              kdime =  pard1*nbopo_par(kfl_desti_par)
!!$              if( kdime > 0 ) call par_srinte(1_ip,kdime,loc_pari2(1,ii))
!!$              ii = ii + nbopo_par(kfl_desti_par)            
!!$           end do
!!$           call memchk(two,istat,par_memor,'LOC_PARI2','par_mygather',loc_pari2)
!!$           deallocate(loc_pari2,stat=istat)
!!$           if(istat/=0) call memerr(two,'LOC_PARI2','par_mygather',0_ip)
!!$        else if( ISLAVE ) then
!!$           kdime = pard1*nbopo
!!$           if( kdime > 0 ) call par_srinte(2_ip,kdime,pari2)
!!$        end if

     else

        call runend('PAR_MYGATHER: THIS GATHER DOES NOT EXIST')

     end if

  end if

  npari = 0
  nparr = 0
  nparc = 0

end subroutine par_mygather
