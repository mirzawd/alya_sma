!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine par_sengeo(itask)
  !------------------------------------------------------------------------
  !****f* Parall/par_sengeo
  ! NAME
  !    par_sengeo
  ! DESCRIPTION
  !    Send 
  !    LNINV_LOC:  1->npoin to 1->npoin_total
  !    XLNIN_LOC:  Pointer to LNINV_LOC   
  ! OUTPUT
  ! USED BY
  !    Domain
  !***
  !------------------------------------------------------------------------
  use def_kintyp
  use def_parame
  use def_parall
  use def_domain
  use def_master
  use mod_memory
  use mod_elmgeo,                        only : elmgeo_number_nodes
  use mod_domain,                        only : domain_memory_deallocate     
  use mod_domain,                        only : domain_memory_allocate     
  use mod_parall,                        only : par_memor
  use mod_communications_point_to_point, only : PAR_SEND
  use mod_mpio_config,                   only : mpio_config

  implicit none  
  integer(ip), intent(in) :: itask
  integer(ip)             :: ielem,jelem,inode,jpoin,ipoin,inodb,jboun,istep
  integer(ip)             :: indice0,indice1,indice2,indice3,indice4,indice5
  integer(ip)             :: ii,jj,iboun_loc,ielem_loc,iposi,iboun,nnodb
  integer(ip)             :: domai,idime,kdime,ifiel,ipoin_local,ipoin_global
  integer(ip)             :: nsteps,what_to_do,pnode

  integer(ip), pointer    :: indice_dom(:)  
  integer(ip), pointer    :: lnods_loc(:)   
  integer(ip), pointer    :: lnodb_loc(:)   
  integer(ip), pointer    :: ltype_loc(:)   
  integer(ip), pointer    :: ltypb_loc(:)   
  integer(ip), pointer    :: lelbo_loc(:) 
  integer(ip), pointer    :: lelch_loc(:)   
  integer(ip), pointer    :: leinv_tmp(:)   
  integer(ip), pointer    :: lninv_tmp(:)
  integer(ip), pointer    :: lgrou_loc(:)     
  integer(ip), pointer    :: lnoch_loc(:)   
  integer(ip), pointer    :: lmast_loc(:)
  integer(ip), pointer    :: lbinv_tmp(:)   
  integer(ip), pointer    :: lboch_loc(:)   
  integer(ip), pointer    :: lesub_loc(:)   
  integer(ip), pointer    :: lmate_loc(:)   
  real(rp),    pointer    :: coord_loc(:)   
  logical(lg)             :: found
  logical(lg)             :: if_parallel_preprocess
  !
  ! Are we in a parallel preprocess?
  !
  if( kfl_ptask == 0 .and. nproc_par > 1 ) then
     if_parallel_preprocess = .true.
     what_to_do             = 1
  else
     if_parallel_preprocess = .false.
     what_to_do             = 2
  end if

  nullify(indice_dom)    
  nullify(lnods_loc)   
  nullify(lnodb_loc)   
  nullify(ltype_loc)   
  nullify(ltypb_loc)   
  nullify(lelbo_loc) 
  nullify(lelch_loc)   
  nullify(leinv_tmp)   
  nullify(lninv_tmp)   
  nullify(lgrou_loc)   
  nullify(lnoch_loc)   
  nullify(lmast_loc)   
  nullify(lbinv_tmp)   
  nullify(lboch_loc)   
  nullify(lesub_loc)   
  nullify(lmate_loc)   
  nullify(coord_loc)   

  select case( itask ) 

  case( 1_ip )

     if( IMASTER .and. ( .not. if_parallel_preprocess ) ) then 

        !----------------------------------------------------------------
        !
        ! Master
        !
        !----------------------------------------------------------------
        !
        ! Nodal arrays
        !
        call memory_alloca(par_memor,'LNOCH_LOC','par_sengeo',lnoch_loc,npoin_total)
        call memory_alloca(par_memor,'LMAST_LOC','par_sengeo',lmast_loc,npoin_total)
        call memory_alloca(par_memor,'LNINV_TMP','par_sengeo',lninv_tmp,npoin_total)
        call memory_alloca(par_memor,'COORD_LOC','par_sengeo',coord_loc,npoin_total*ndime)
        !
        ! LNINV, LNOCH, COORD
        !
        ii = 0
        jj = 0
        do ipoin = 1,npoin_total
           jj            = jj + 1
           jpoin         = lninv_loc(ipoin)
           lninv_tmp(jj) = jpoin
           lnoch_loc(jj) = lnoch(jpoin)
           do idime = 1,ndime
              ii = ii + 1
              coord_loc(ii) = coord(idime,jpoin)
           end do
        end do
        !
        ! LMAST: transform to local numbering
        !
        if( new_periodicity == 0 ) then
           jj = 0
           indice0 = 0
           do domai = 1,npart_par
              do ipoin = 1,npoin_par(domai)
                 jj           = jj + 1
                 jpoin        = lninv_loc(jj)
                 ipoin_global = lmast(jpoin)              
                 if( ipoin_global /= 0 ) then
                    found        =.true.
                    ipoin_local  = 0
                    do while( ipoin_local < npoin_par(domai) .and. found )
                       ipoin_local = ipoin_local + 1
                       if( lninv_loc(indice0+ipoin_local) == ipoin_global ) then
                          found         = .false.
                          lmast_loc(jj) = ipoin_local
                          ipoin_local   = npoin_par(domai)
                       end if
                    end do
                 else
                    lmast_loc(jj) = 0
                 end if
              end do
              indice0 = indice0 + npoin_par(domai)
           end do
        else
           jj = 0
           do ipoin = 1,npoin_total
              jj            = jj + 1
              jpoin         = lninv_loc(ipoin)
              lmast_loc(jj) = lmast(jpoin)
           end do
        end if
        !
        ! Elements arrays
        !
        call memory_alloca(par_memor,'LNODS_LOC','par_sengeo',lnods_loc,mnode*nelem)
        call memory_alloca(par_memor,'LTYPE_LOC','par_sengeo',ltype_loc,nelem)
        call memory_alloca(par_memor,'LELCH_LOC','par_sengeo',lelch_loc,nelem)
        call memory_alloca(par_memor,'LEINV_TMP','par_sengeo',leinv_tmp,nelem)
        call memory_alloca(par_memor,'LESUB_LOC','par_sengeo',lesub_loc,nelem)
        call memory_alloca(par_memor,'LMATE_LOC','par_sengeo',lmate_loc,nelem)

        do ielem = 1, nelem
           jelem            = leinv_par(ielem)
           domai            = lepar_par(jelem)
           ltype_loc(ielem) = ltype(jelem)
           lelch_loc(ielem) = lelch(jelem)
           leinv_tmp(ielem) = jelem
           lesub_loc(ielem) = lesub(jelem)
           lmate_loc(ielem) = lmate(jelem)
           !
           ! Special and finite elements
           !
           pnode = elmgeo_number_nodes(ltype(jelem),lnods(:,jelem))
           !
           ! LNODS
           !
           do inode = 1,pnode
              jpoin = lnper_par( lnods(inode,jelem) )
              iposi = ( ielem - 1 ) * mnode + inode
              if( jpoin <= gni ) then 
                 !
                 ! Internal node
                 !
                 lnods_loc(iposi) = jpoin - ginde_par(1,domai) + 1
              else  
                 !
                 ! Boundary node
                 !
                 jpoin = jpoin - gni
                 ii    = badj(jpoin)
                 do while( bdom(ii) /= domai )
                    ii = ii + 1
                    if( ii== badj(jpoin+1) ) &
                         call runend('PAR_SENGEO (LNODB): '&
                         //' GLOBAL NODE '//trim(intost(lnods(inode,ielem)))&
                         //' ('//trim(intost(jpoin+gni))//' IN I/B NUMBERING)'&
                         //' (LOCAL NODE '//trim(intost(inode))&
                         //' OF BOUNDARY '//trim(intost(ielem))//')'&
                         //' DOES NOT BELONG TO DOMAIN '&
                         //trim(intost(domai)))
                 end do
                 lnods_loc(iposi) = bpoin(ii)                 
              endif
           end do

        end do
        !
        ! Boundary arrays
        !
        call memory_alloca(par_memor,'LNODB_LOC','par_sengeo',lnodb_loc,nboun*mnodb)
        call memory_alloca(par_memor,'LTYPB_LOC','par_sengeo',ltypb_loc,nboun)
        call memory_alloca(par_memor,'LBOCH_LOC','par_sengeo',lboch_loc,nboun)
        call memory_alloca(par_memor,'LBINV_TMP','par_sengeo',lbinv_tmp,nboun)

        do iboun = 1,nboun
           jboun            = lbper_par(iboun)
           domai            = lbpar_par(iboun)
           ltypb_loc(jboun) = ltypb(iboun)
           lboch_loc(jboun) = lboch(iboun)
           lbinv_tmp(jboun) = iboun
           do inodb = 1,nnode(ltypb(iboun))
              jpoin = lnper_par( lnodb(inodb,iboun) )
              iposi = ( jboun - 1 ) * mnodb + inodb
              if( jpoin <= gni ) then
                 !
                 ! Internal node
                 !
                 lnodb_loc(iposi) = jpoin - ginde_par(1,domai) + 1
              else 
                 !
                 ! Boundary node
                 !
                 jpoin = jpoin - gni
                 ii    = badj(jpoin)
                 do while( bdom(ii) /= domai )
                    ii = ii + 1
                    if( ii == badj(jpoin+1) ) &
                         call runend('PAR_SENGEO (LNODB): '&
                         //' GLOBAL NODE '//trim(intost(lnodb(inodb,iboun)))&
                         //' ('//trim(intost(jpoin+gni))//' IN I/B NUMBERING)'&
                         //' (LOCAL NODE '//trim(intost(inodb))&
                         //' OF BOUNDARY '//trim(intost(iboun))//')'&
                         //' DOES NOT BELONG TO DOMAIN '&
                         //trim(intost(domai)))
                 end do
                 lnodb_loc(iposi) = bpoin(ii)
              endif
           end do
        end do
        !
        ! Send data
        !
        nparc   = 0
        indice0 = 1
        indice1 = 1
        indice2 = 1
        indice3 = 1
        indice4 = 1
        indice5 = 1
        do domai= 1, npart_par
           kfl_desti_par = domai
           !
           ! Send LNOCH_LOC
           !
           npari =  npoin_par(domai)
           parin => lnoch_loc(indice0:)
           strin =  'LNOCH_LOC'
           call par_sendin()
           !
           ! Send LMAST_LOC
           !
           npari =  npoin_par(domai)
           parin => lmast_loc(indice0:)
           strin =  'LMAST_LOC'
           call par_sendin()
           !
           ! Send LNINV_TMP
           !
           npari =  npoin_par(domai)
           parin => lninv_tmp(indice0:)
           strin =  'LNINV_TMP'
           call par_sendin()
           !
           ! Send COORD_LOC
           !
           nparr =  ndime*npoin_par(domai)
           parre => coord_loc(indice1:)
           strre =  'COORD_LOC'
           call par_sendin()
           !
           ! Send LNODS_LOC
           !
           npari =  mnode*nelem_par(domai)
           parin => lnods_loc(indice4:)
           strin =  'LNODS_LOC'
           call par_sendin()
           !
           ! Send LTYPE_LOC
           !
           npari =  nelem_par(domai)
           parin => ltype_loc(indice2:)
           strin =  'LTYPE_LOC'
           call par_sendin()
           !
           ! Send LELCH_LOC
           !
           npari =  nelem_par(domai)
           parin => lelch_loc(indice2:)
           strin =  'LELCH_LOC'
           call par_sendin()
           !
           ! Send LEINV_TMP
           !
           npari =  nelem_par(domai)
           parin => leinv_tmp(indice2:)
           strin =  'LEINV_TMP'
           call par_sendin()
           !
           ! Send LESUB_LOC
           !
           npari =  nelem_par(domai)
           parin => lesub_loc(indice2:)
           strin =  'LESUB_LOC'
           call par_sendin()
           !
           ! Send LMATE_LOC
           !
           npari =  nelem_par(domai)
           parin => lmate_loc(indice2:)
           strin =  'LMATE_LOC'
           call par_sendin()

           if( nboun_par(domai) > 0 ) then
              !
              ! Send LNODB_LOC
              !
              npari =  mnodb*nboun_par(domai)
              parin => lnodb_loc(indice5:)
              strin =  'LNODB_LOC'
              call par_sendin()
              !
              ! Send LTYPB_LOC
              !
              npari =  nboun_par(domai)
              parin => ltypb_loc(indice3:)
              strin =  'LTYPB_LOC'
              call par_sendin()
              !
              ! Send LBOCH_LOC
              !
              npari =  nboun_par(domai)
              parin => lboch_loc(indice3:)
              strin =  'LBOCH_LOC'
              call par_sendin()
              !
              ! Send LBINV_TMP
              !
              npari =  nboun_par(domai)
              parin => lbinv_tmp(indice3:)
              strin =  'LBINV_TMP'
              call par_sendin()
           end if

           indice0 = indice0 + npoin_par(domai)
           indice1 = indice1 + npoin_par(domai) * ndime
           indice2 = indice2 + nelem_par(domai)
           indice3 = indice3 + nboun_par(domai)
           indice4 = indice4 + nelem_par(domai) * mnode
           indice5 = indice5 + nboun_par(domai) * mnodb

        end do
        !
        ! Deallocate memory
        !
        call memory_deallo(par_memor,'LBINV_TMP','par_sengeo',lbinv_tmp)
        call memory_deallo(par_memor,'LBOCH_LOC','par_sengeo',lboch_loc)
        call memory_deallo(par_memor,'LTYPB_LOC','par_sengeo',ltypb_loc)
        call memory_deallo(par_memor,'LNODB_LOC','par_sengeo',lnodb_loc)

        call memory_deallo(par_memor,'LESUB_LOC','par_sengeo',lesub_loc)
        call memory_deallo(par_memor,'LMATE_LOC','par_sengeo',lmate_loc)
        call memory_deallo(par_memor,'LEINV_TMP','par_sengeo',leinv_tmp)
        call memory_deallo(par_memor,'LELCH_LOC','par_sengeo',lelch_loc)
        call memory_deallo(par_memor,'LTYPE_LOC','par_sengeo',ltype_loc)
        call memory_deallo(par_memor,'LNODS_LOC','par_sengeo',lnods_loc)

        call memory_deallo(par_memor,'COORD_LOC','par_sengeo',coord_loc)
        call memory_deallo(par_memor,'LNINV_TMP','par_sengeo',lninv_tmp)
        call memory_deallo(par_memor,'LMAST_LOC','par_sengeo',lmast_loc)
        call memory_deallo(par_memor,'LNOCH_LOC','par_sengeo',lnoch_loc)

        if( kfl_ngrou == -1 ) then
           !
           ! Groups
           !
           call memory_alloca(par_memor,'LGROU_LOC','par_sengeo', lgrou_loc,npoin_total)
           jj = 0
           do ipoin = 1, npoin_total
              jj    = jj + 1
              jpoin = lninv_loc(ipoin)
              lgrou_loc(jj) = lgrou_dom(jpoin)
           end do
           !
           ! Send data
           !
           nparc   = 0
           indice0 = 1
           do domai= 1, npart_par
              kfl_desti_par = domai
              !
              ! Send LGROU_LOC
              !
              npari =  npoin_par(domai)
              parin => lgrou_loc(indice0:)
              strin =  'LGROU_LOC'
              call par_sendin()           
              indice0 = indice0 + npoin_par(domai)              
           end do
           call domain_memory_deallocate('LGROU_DOM')
           call memory_deallo(par_memor,'LGROU_LOC','par_sengeo', lgrou_loc)

        end if

     else if( ISLAVE ) then

        !----------------------------------------------------------------
        !
        ! Slaves
        !
        !----------------------------------------------------------------

        if( .not. associated(lninv_loc) ) call memory_alloca(par_memor,'LNINV_LOC','par_sengeo',lninv_loc,npoin) 
        if( .not. associated(leinv_loc) ) call memory_alloca(par_memor,'LEINV_LOC','par_sengeo',leinv_loc,nelem) 
        if( .not. associated(lbinv_loc) ) call memory_alloca(par_memor,'LBINV_LOC','par_sengeo',lbinv_loc,max(1_ip,nboun))

        if( npoin > 0 ) then
           !
           ! Receive LNOCH_LOC
           !
           call par_srinte(what_to_do,npoin,lnoch)
           !
           ! Receive LMAST_LOC
           !
           call par_srinte(what_to_do,npoin,lmast)
           !
           ! Receive LNINV_LOC
           !
           call par_srinte(what_to_do,npoin,lninv_loc) 
           !
           ! Receive COORD
           !
           call par_srreal(what_to_do,ndime*npoin,coord)

        end if
        
        if( nelem > 0 ) then
           !
           ! Receive LNODS
           !
           call par_srinte(what_to_do,mnode*nelem,lnods)
           !
           ! Receive LTYPE
           !
           call par_srinte(what_to_do,nelem,ltype)
           !
           ! Receive LELCH
           !
           call par_srinte(what_to_do,nelem,lelch)
           !
           ! Receive LEINV_LOC
           !
           call par_srinte(what_to_do,nelem,leinv_loc) 
           !
           ! Receive LESUB_LOC
           !
           call par_srinte(what_to_do,nelem,lesub) 
           !
           ! Receive LMATE_LOC
           !
           call par_srinte(what_to_do,nelem,lmate) 

        end if

        if( nboun > 0 ) then
           !
           ! Receive LNODB
           !          
           call par_srinte(what_to_do,mnodb*nboun,lnodb)
           !
           ! Receive LTYPB
           !
           call par_srinte(what_to_do,nboun,ltypb)
           !
           ! Receive LBOCH
           !
           call par_srinte(what_to_do,nboun,lboch)
           !
           ! Receive LBINV_LOC
           !
           call par_srinte(what_to_do,nboun,lbinv_loc)
        end if
        !
        ! Groups
        !
        if( kfl_ngrou == -1 .and. npoin > 0 ) then
           call domain_memory_allocate('LGROU_DOM')
           call par_srinte(what_to_do,npoin,lgrou_dom)
        end if

     end if

  case( 2_ip )
     !
     ! Special case of LELBO
     !
     if( IMASTER .and. kfl_ptask /= 2 ) then

        call memory_alloca(par_memor,'LELBO_LOC','par_sengeo', lelbo_loc,nboun)
        call memory_alloca(par_memor,'INDICE_DOM','par_sengeo',indice_dom,npart_par+1_ip)
        indice_dom(1) = 1
        do domai= 1,npart_par
           indice_dom(domai+1)  = indice_dom(domai) + nelem_par(domai)
        end do

        do iboun = 1, nboun
           iboun_loc            = lbper_par(iboun)
           nnodb                = nnode(ltypb(iboun))
           domai                = lbpar_par(iboun)
           ielem                = lelbo(iboun)
           ielem_loc            = leper_par(ielem)-indice_dom(domai)+1
           lelbo_loc(iboun_loc) = ielem_loc
        end do
        nparc   = 0
        nparr   = 0
        indice3 = 1
        do domai= 1, npart_par
           kfl_desti_par = domai
           !
           ! Send LELBO_LOC
           !
           if( nboun_par(domai) > 0 ) then
              kdime =  nboun_par(domai)
              strin =  'LELBO_LOC'
              call par_srinte(1_ip,kdime,lelbo_loc(indice3))
           end if
           indice3 = indice3 + nboun_par(domai)
        end do
        !
        ! Deallocate memory
        !
        call memory_deallo(par_memor,'INDICE_DOM','par_sengeo',indice_dom)
        call memory_deallo(par_memor,'LELBO_LOC' ,'par_sengeo', lelbo_loc)

     else if(ISLAVE) then
        !
        ! Receive LELBO
        !
        if( nboun > 0 ) then
           call par_srinte(2_ip,nboun,lelbo)
        end if

     end if

  case( 3_ip )

     !----------------------------------------------------------------------
     !
     ! Fields
     !
     !----------------------------------------------------------------------

     if( nfiel > 0 ) then

        do ifiel = 1,nfiel
           if( kfl_field(1,ifiel) /= 0 ) then
              !
              ! This field has been defined
              !
              if( ISLAVE ) call domain_memory_allocate('XFIEL % A',NUMBER1=ifiel)

              if(      kfl_field(2,ifiel) == NPOIN_TYPE ) then
                 kfl_field(5,ifiel) = npoin
              else if( kfl_field(2,ifiel) == NELEM_TYPE ) then
                 kfl_field(5,ifiel) = nelem
              else
                 kfl_field(5,ifiel) = nboun
              end if


              if ( (kfl_field(6,ifiel) /= 1) .OR. (mpio_config%output%post_process%export_only) ) then 
                 nsteps = kfl_field(4,ifiel)
              else
                 nsteps = nsteps_fiel_ondemand
              end if

              if(      kfl_field(2,ifiel) == NPOIN_TYPE ) then
                 do istep = 1,nsteps
                    call par_pararr('GAT',NPOIN_TYPE,kfl_field(1,ifiel)*kfl_field(5,ifiel),xfiel(ifiel) % a(1,1,istep))
                 end do
              else if( kfl_field(2,ifiel) == NELEM_TYPE ) then
                 do istep = 1,nsteps
                    call par_pararr('GAT',NELEM_TYPE,kfl_field(1,ifiel)*kfl_field(5,ifiel),xfiel(ifiel) % a(1,1,istep))
                 end do
              else
                 do istep = 1,nsteps
                    call par_pararr('GAT',NBOUN_TYPE,kfl_field(1,ifiel)*kfl_field(5,ifiel),xfiel(ifiel) % a(1,1,istep))
                 end do
              end if

              if( IMASTER ) call domain_memory_deallocate('XFIEL % A',NUMBER1=ifiel)

           end if

        end do

     end if

  end select

end subroutine par_sengeo
