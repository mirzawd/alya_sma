!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine par_comset(itask)
  !------------------------------------------------------------------------
  !****f* Parall/par_comset
  ! NAME
  !    par_senset
  ! DESCRIPTION
  !    Send/receive sets
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
  use def_kermod
  use mod_memchk  
  use mod_parall,         only : PAR_COMM_MY_CODE,PAR_INTEGER,commd
  use mod_communications, only : PAR_SEND_RECEIVE
  use mod_communications, only : PAR_SEND_RECEIVE_RP
  use mod_parall,         only : par_memor
  use def_mpi
  implicit none  

  integer(ip), intent(in) :: itask
  integer(ip)             :: inset,domai,nvari,ivars,knset,ii
  integer(ip)             :: iwitn,ivawi,kwitn,icoun,ineig,inset_total
  integer(ip)             :: ini,bsize,isize,isiz1,isiz2  
  integer(ip), pointer    :: loc_spari1(:),loc_rpari1(:)
  integer                 :: istat,my_sendcount,nsize4,bsize4,dom_i
  integer(ip)             :: ifina
  integer(ip), pointer    :: my_recvbuf(:)
  integer(4),  pointer    :: my_displs(:),my_recvcounts(:)
  integer(ip)             :: my_sendbuf(4)
  integer(ip), pointer    :: mark_witness(:)
!  real(rp)                :: rnul2(1,1)
  
  nullify(loc_spari1)
  nullify(loc_rpari1)
  nullify(my_recvbuf)
  nullify(my_displs)
  nullify(my_recvcounts)
  nullify(mark_witness)
  
  select case( itask )

  case ( 1_ip )

     !-------------------------------------------------------------------
     !
     ! Node sets
     !
     !-------------------------------------------------------------------

     if( IMASTER ) then     
        !
        ! Master: receive set values 
        !
        nvari       = nparr/nnset
        isize       = nparr
        inset_total = 0
        call memgen(0_ip,isize,0_ip)
        
        do domai =  1,npart_par
           nparr =  nvari * nnset_par(domai)
           if( nparr > 0 ) then
              parre => gesca
              kfl_desti_par = domai
              call par_receiv()
              
              do inset = 1,nnset_par(domai)
                 inset_total = inset_total + 1 
                 knset = lnsec_par(inset_total,2)
                 do ivars = 1,nvari
                    vnset(ivars,knset) = parre((inset-1)*nvari+ivars)
                 end do
              end do
              
              nullify(parre)
           end if
        end do 
        call memgen(2_ip,isize,0_ip)

     else if( ISLAVE ) then
        !
        ! Slaves: end set values
        !
        kfl_desti_par = 0
        if( nparr > 0 ) call par_srreal(1_ip,nparr,vnset)

     end if

!!$    if( IMASTER ) then     
!!$       !
!!$       ! Master: receive set values 
!!$       !
!!$       nvari = nparr/nnset
!!$       do domai =  1,npart_par
!!$          if( nnset_par(domai) > 0 ) then
!!$             
!!$             nparr = nvari * nnset_par(domai)
!!$             call memgen(0_ip,nvari,nnset_par(domai))
!!$             
!!$             call PAR_SEND_RECEIVE_RP(0_ip,nparr,rnul2,gevec,'IN MY CODE',domai,'SYNCHRONOUS')
!!$             
!!$             do inset = 1,nnset
!!$                if( lnsec_par(inset,1) == domai ) then
!!$                   knset = lnsec_par(inset,2)
!!$                   print*,'recv=',domai,inset,knset 
!!$                   do ivars = 1,nvari
!!$                      vnset(ivars,inset) = gevec(ivars,knset)
!!$                   end do
!!$                end if
!!$             end do
!!$             call memgen(2_ip,nvari,nnset_par(domai))
!!$             
!!$          end if
!!$       end do
!!$       
!!$     else if( ISLAVE ) then
!!$        !
!!$        ! Slaves: end set values
!!$        !
!!$        kfl_desti_par = 0
!!$        if( nparr > 0 ) call PAR_SEND_RECEIVE_RP(nparr,0_ip,vnset,rnul2,'IN MY CODE',0_ip,'SYNCHRONOUS')
!!$
!!$     end if

  case ( 2_ip )

     !-------------------------------------------------------------------
     !
     ! Witness points
     !
     !-------------------------------------------------------------------

     !do kfl_desti_par = 1,npart_par
     !   nwitn_par(kfl_desti_par) = nwitn_par(kfl_desti_par) * postp(1) % nvawi
     !end do
     !call PAR_GATHERV(witne,witne,nwitn_par,'IN MY CODE')
     !do kfl_desti_par = 1,npart_par
     !   nwitn_par(kfl_desti_par) = nwitn_par(kfl_desti_par) / postp(1) % nvawi
     !end do
     !return

     if( IMASTER ) then      
        !
        ! Master: receive witness values 
        !        
        isiz1 = postp(1)%nvawi
        isiz2 = nwitn
        call memgen(0_ip,isiz1,isiz2)
        kwitn = 0
        do kfl_desti_par = 1,npart_par
           if( nwitn_par(kfl_desti_par) > 0 ) then
              nparr =  nwitn_par(kfl_desti_par) * postp(1) % nvawi
              call par_srreal(2_ip,nparr,gevec)
              icoun = 0
              do iwitn = 1,nwitn_par(kfl_desti_par) 
                 kwitn = kwitn + 1
                 do ivawi = 1,postp(1) % nvawi
                    icoun = icoun + 1
                    witne( ivawi , lnwit(kwitn) ) = gevec(ivawi,iwitn)
                 end do
              end do
           end if
        end do
        call memgen(2_ip,isiz1,isiz2)

     else if( ISLAVE ) then
        !
        ! Slaves: send witness values 
        !         
        kfl_desti_par = 0
        nparr         = postp(1) % nvawi * nwitn
        if( nparr > 0 ) call par_srreal(1_ip,nparr,postp(1) % witne)

     end if

  case ( 3_ip )

     !-------------------------------------------------------------------
     ! 
     ! Witness points: define who owns it
     !
     !-------------------------------------------------------------------

     if( ISLAVE ) then

        allocate(loc_spari1(nwitn),stat=istat)
        call memchk(zero,istat,par_memor,'LOC_SPARI1','par_slexch',loc_spari1)

        allocate(loc_rpari1(nwitn*nneig),stat=istat)
        call memchk(zero,istat,par_memor,'LOC_RPARI1','par_slexch',loc_rpari1)

        do iwitn = 1,nwitn
           loc_spari1(iwitn) = lewit(iwitn)
        end do

        do ineig = 1,nneig

           dom_i = commd%neights(ineig)
           ini   = (ineig-1) * nwitn + 1
           bsize = nwitn

#ifdef MPI_OFF
#else
           bsize4 = int(bsize,4)
           call MPI_Sendrecv( loc_spari1(1:), bsize4, &
                PAR_INTEGER,  dom_i, 0,               &
                loc_rpari1(ini:), bsize4,             &
                PAR_INTEGER, dom_i, 0,                &
                PAR_COMM_MY_CODE, status, istat )
#endif
        end do

        do iwitn = 1,nwitn
           do ineig = 1,nneig
              ii    = (ineig-1) * nwitn + iwitn
              dom_i = commd % neights(ineig)
              if( lewit(iwitn) > 0 .and. loc_rpari1(ii) > 0 ) then
                 if( kfl_paral > dom_i ) then
                    lewit(iwitn) = 0
                 end if
              end if
           end do
        end do

        call memchk(two,istat,par_memor,'LOC_RPARI1','par_slexch',loc_rpari1)
        deallocate(loc_rpari1,stat=istat)
        if(istat/=0) call memerr(two,'LOC_RPARI1','par_slexch',0_ip)

        call memchk(two,istat,par_memor,'LOC_SPARI1','par_slexch',loc_spari1)
        deallocate(loc_spari1,stat=istat)
        if(istat/=0) call memerr(two,'LOC_SPARI1','par_slexch',0_ip)

     end if

  case ( 4_ip )

     !----------------------------------------------------------------------
     !
     ! NWITN_PAR: number witness point for each subdomain
     !
     !----------------------------------------------------------------------

     my_sendcount = 1_4
     ifina     = npart_par+1
     nsize4    = int(ifina,4)

     allocate( my_displs(npart_par+1) )
     allocate( my_recvcounts(npart_par+1) )
     allocate( my_recvbuf(ifina) )
     
     do ii = 1,npart_par+1
        my_displs(ii)     = 1*(int(ii,4)-1)
        my_recvcounts(ii) = 1
     end do
     do ii = 1,ifina
        my_recvbuf(ii) = 0
     end do
     my_sendbuf(1) = nwitn
     !
     ! MPI all gather
     !
#ifdef MPI_OFF
#else
     CALL MPI_ALLGATHERV(                             &
          my_sendbuf(1:1), my_sendcount, PAR_INTEGER, &
          my_recvbuf(1:ifina), my_recvcounts ,        &
          my_displs(1:npart_par+1), PAR_INTEGER,      &
          PAR_COMM_MY_CODE, istat                     )
#endif

     if( IMASTER ) then
        call par_memset(2_ip)        
        do ii = 1,npart_par
           nwitn_par(ii) = my_recvbuf(ii+1)
        end do
     end if

     deallocate( my_recvcounts )
     deallocate( my_displs     ) 
     deallocate( my_recvbuf    )

     if( ISLAVE ) then
        kfl_desti_par = 0
        if( nwitn > 0 ) call par_parari('SND',0_ip,nwitn,lnwit)
     else
        kwitn = 1
        do kfl_desti_par = 1,npart
           if( nwitn_par(kfl_desti_par) > 0 ) then
              call par_parari('RCV',0_ip,nwitn_par(kfl_desti_par),lnwit(kwitn))
              kwitn = kwitn + nwitn_par(kfl_desti_par)
           end if
        end do
        mwitn = kwitn - 1
     end if

  end select

end subroutine par_comset
