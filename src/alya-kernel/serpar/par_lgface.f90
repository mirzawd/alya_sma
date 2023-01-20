!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine par_lgface(itask)
  !-----------------------------------------------------------------------
  !****f* Parall/par_lgface
  ! NAME
  !    par_lgface
  ! DESCRIPTION
  !    This subroutine 
  !    commd%bound_perm(jj):      permutation array
  !    loc_sparr1:                my local values
  !    loc_rparr1:                values given by neighbor ii
  !    commd%bound_dim:           size of communication array
  !    nneig:                     number of neighbors that share same group
  !    commd%neights(ii):         number of subdomain ii
  !    commd%bound_size(ii):      where my local arrays sart to exchange with ii
  !    commd%bound_size(ii+1)
  !        -commd%bound_size(ii): number of groups to exchange with ii
  ! USED BY
  !    Parall
  !***
  !-----------------------------------------------------------------------
  use def_master
  use def_parame
  use def_parall
  use def_domain
  use mod_memchk
  use mod_parall, only : PAR_INTEGER,commd
  use mod_parall, only : PAR_REAL
  use mod_parall, only : PAR_COMM_MY_CODE
  use mod_parall, only : par_memor
  use mod_memory
  use def_mpi
  implicit none

  integer(ip), intent(in)  :: itask
  integer(ip)              :: ii,jj,bsize,ji,poin,ini
  integer(ip)              :: ielem,jelem,ifacg,kk,isize,kfacg
  integer(ip)              :: ielem_global,jelem_global
  integer(4)               :: istat,bsize4,dom_i
  integer(ip), allocatable :: lorf1(:),lorf2(:)
  integer(ip), allocatable :: loc_spari1(:),   loc_rpari1(:)
  real(rp),    allocatable :: loc_sparr1(:),   loc_rparr1(:)
  complex(rp), allocatable :: loc_sparx1(:),   loc_rparx1(:)

  if( ISLAVE ) then

     select case ( itask )

     case ( 1_ip ) 

        !----------------------------------------------------------------
        !
        ! Construct face communication arrays
        !
        !----------------------------------------------------------------

        call memgen(1_ip,nneig,0_ip)

        commd % bface_dim = 0

        do ifacg = 1,nfacg
           ielem = lfacg(1,ifacg)
           jelem = lfacg(2,ifacg)
           ii    = 0
           if( ielem <= nelem .and. jelem > nelem ) then
              ii = leldo(1,jelem-nelem) 
           else if( ielem > nelem .and. jelem <= nelem ) then
              ii = leldo(1,ielem-nelem) 
           end if
           if( ii /= 0 ) then
              gisca(ii) = gisca(ii) + 1
              commd % bface_dim = commd % bface_dim + 1
           end if
        end do
        !
        ! Allocate permutation arrays
        !
        call memory_alloca(par_memor,'COMMD % BFACE_PERM','par_lgface',commd % bface_perm,commd % bface_dim)
        call memory_alloca(par_memor,'COMMD % BFACE_SIZE','par_lgface',commd % bface_size,nneig+1) 
        commd % bface_size(1) = 1
        do ii = 2,nneig+1
           commd % bface_size(ii) = commd % bface_size(ii-1) + gisca(ii-1)
        end do
        !
        ! LORF1, LORF2: Temporal arrays
        !
        allocate( lorf1(commd % bface_dim) , stat = istat )
        call memchk(zero,istat,par_memor,'LORF1','par_lgface',lorf1)
        allocate( lorf2(commd % bface_dim) , stat = istat )
        call memchk(zero,istat,par_memor,'LORF2','par_lgface',lorf2)
        !
        ! List of faces
        !
        kk = 0
        do ii = 1,nneig
           gisca(ii) = 0 
        end do

        do ifacg = 1,nfacg
           ielem = lfacg(1,ifacg)
           jelem = lfacg(2,ifacg)
           ii    = 0
           if( ielem <= nelem .and. jelem > nelem ) then
              ii = leldo(1,jelem-nelem) 
           else if( ielem > nelem .and. jelem <= nelem ) then
              ii = leldo(1,ielem-nelem) 
           end if
           if( ii /= 0 ) then
              gisca(ii)    = gisca(ii) + 1
              ini          = commd % bface_size(ii)
              kk           = ini + gisca(ii) - 1
              ielem_global = leinv_loc(ielem)
              jelem_global = leinv_loc(jelem)
              commd % bface_perm(kk) = ifacg
              if( jelem_global < ielem_global ) then
                 lorf1(kk)  = ielem_global
                 lorf2(kk)  = jelem_global
              else
                 lorf1(kk)  = jelem_global
                 lorf2(kk)  = ielem_global
              end if
           end if
        end do
        !
        ! Order the faces COMMD % BFACE_PERM using LORF1 and LORF2
        !
        do ii = 1,nneig
           ini   = commd % bface_size(ii)
           isize = commd % bface_size(ii+1) - commd % bface_size(ii)
           if( isize /= 0 ) call par_sorti2(isize,lorf1(ini),lorf2(ini),commd % bface_perm(ini))
        end do
        !
        ! Deallocate memory
        !
        call memchk(two,istat,par_memor,'LORF1','par_lgface',lorf1)
        deallocate( lorf1 , stat = istat )
        if(istat/=0) call memerr(two,'LORF1','par_lgface',0_ip)
        call memchk(two,istat,par_memor,'LORF2','par_lgface',lorf2)
        deallocate( lorf2 , stat = istat )
        if(istat/=0) call memerr(two,'LORF2','par_lgface',0_ip)

        call memgen(3_ip,nneig,0_ip)

        !do ii = 1, nneig
        !   do kk = commd % bface_size(ii),commd % bface_size(ii+1)-1
        !      commd % bface_perm(kk) = lorfa(kk,1)
        !      ifacg = lorfa(kk,1)
        !      ifacg = commd % bface_perm(kk)
        !      !write(100+kfl_paral,*) commd % bface_perm(kk),leinv_loc(lfacg(1,ifacg)),leinv_loc(lfacg(2,ifacg))
        !      !write(100+kfl_paral,*) lorf1(kk),lorf2(kk)
        !   end do
        !end do
        !do ii = 1,nneig
        !   print*,commd % neights(ii),gisca(ii)
        !end do

     case ( 2_ip ) 

        !----------------------------------------------------------------
        !
        ! Exchange face arrays
        !
        !----------------------------------------------------------------

        if( pardi == 1 .and. parki == 1 ) then

           !-------------------------------------------------------------
           !
           ! INT(NFACG)
           !
           !-------------------------------------------------------------

           allocate(loc_spari1(commd % bface_dim),stat=istat)
           call memchk(zero,istat,par_memor,'LOC_SPARI1','par_slexch',loc_spari1)

           allocate(loc_rpari1(commd % bface_dim),stat=istat)
           call memchk(zero,istat,par_memor,'LOC_RPARI1','par_slexch',loc_rpari1)

           do jj = 1,commd % bface_dim
              ifacg = commd % bface_perm(jj)
              loc_spari1(jj) = pari1(ifacg)
           enddo

           do ii= 1, nneig

              dom_i = commd % neights(ii)
              ini   = commd % bface_size(ii)
              bsize = commd % bface_size(ii+1) - ini

#ifdef MPI_OFF
#else
              if( bsize > 0 ) then
                 !
                 ! Check if I have to send something. dom_i could be a neighbor
                 ! with no common faces with me!
                 !
                 bsize4 = int(bsize,4)
                 call MPI_Sendrecv( loc_spari1(ini:), bsize4,&
                      PAR_INTEGER,  dom_i, 0_4,     &
                      loc_rpari1(ini:), bsize4,              &
                      PAR_INTEGER, dom_i, 0_4,      &
                      PAR_COMM_MY_CODE, status, istat )
              end if
#endif
           enddo

           do jj = 1,commd % bface_dim
              ifacg = commd % bface_perm(jj)
              pari1(ifacg) = pari1(ifacg) + loc_rpari1(jj)
           enddo

           call memchk(two,istat,par_memor,'LOC_RPARI1','par_slexch',loc_rpari1)
           deallocate(loc_rpari1,stat=istat)
           if(istat/=0) call memerr(two,'LOC_RPARI1','par_slexch',0_ip)

           call memchk(two,istat,par_memor,'LOC_SPARI1','par_slexch',loc_spari1)
           deallocate(loc_spari1,stat=istat)
           if(istat/=0) call memerr(two,'LOC_SPARI1','par_slexch',0_ip)

        else if( pardi >= 1 .and. parki == 6 ) then

           !-------------------------------------------------------------
           !
           ! INT(PARD1,NPOIN) => INT(PARD1*NPOIN)
           !
           !-------------------------------------------------------------

           allocate(loc_spari1(pard1*commd % bface_dim),stat=istat)
           call memchk(zero,istat,par_memor,'LOC_SPARI1','par_slexch',loc_spari1)

           allocate(loc_rpari1(pard1*commd % bface_dim),stat=istat)
           call memchk(zero,istat,par_memor,'LOC_RPARI1','par_slexch',loc_rpari1)

           do jj = 1,commd % bface_dim
              ifacg = commd % bface_perm(jj)
              kfacg = pard1*(ifacg-1)
              kk    = pard1*(jj-1)
              do ii = 1,pard1
                 kk    = kk + 1
                 kfacg = kfacg + 1
                 loc_spari1(kk) = pari1(kfacg)
              enddo
           enddo

           do ii = 1,nneig

              dom_i = commd%neights(ii)
              ini   = pard1*(commd % bface_size(ii)-1)   + 1
              bsize = pard1*(commd % bface_size(ii+1)-1) - ini + 1

#ifdef MPI_OFF
#else
              if( bsize > 0 ) then
                 !
                 ! Check if I have to send something. dom_i could be a neighbor
                 ! with no common faces with me!
                 !
                 bsize4 = int(bsize,4)
                 call MPI_Sendrecv( loc_spari1(ini:), bsize4,&
                      PAR_INTEGER,  dom_i, 0_4,     &
                      loc_rpari1(ini:), bsize4,              &
                      PAR_INTEGER, dom_i, 0_4,      &
                      PAR_COMM_MY_CODE, status, istat )
              end if
#endif
           end do

           do jj = 1,commd % bface_dim
              ifacg = commd % bface_perm(jj)
              kfacg = pard1*(ifacg-1)
              kk    =  pard1*(jj-1)
              do ii = 1,pard1
                 kk    = kk + 1
                 kfacg = kfacg + 1
                 pari1(kfacg) = pari1(kfacg) + loc_rpari1(kk)
              enddo
           enddo
           call memchk(two,istat,par_memor,'LOC_RPARI1','par_slexch',loc_rpari1)
           deallocate(loc_rpari1,stat=istat)
           if(istat/=0) call memerr(two,'LOC_RPARI1','par_slexch',0_ip)

           call memchk(two,istat,par_memor,'LOC_SPARI1','par_slexch',loc_spari1)
           deallocate(loc_spari1,stat=istat)
           if(istat/=0) call memerr(two,'LOC_SPARI1','par_slexch',0_ip)

        else if( pardi == 1 .and. parki == 2 ) then

           !-------------------------------------------------------------
           !
           ! REAL(NFACG)
           !
           !-------------------------------------------------------------

           allocate(loc_sparr1(commd % bface_dim),stat=istat)
           call memchk(zero,istat,par_memor,'LOC_SPARR1','par_slexch',loc_sparr1)

           allocate(loc_rparr1(commd % bface_dim),stat=istat)
           call memchk(zero,istat,par_memor,'LOC_RPARR1','par_slexch',loc_rparr1)

           do jj = 1,commd % bface_dim
              ifacg = commd % bface_perm(jj)
              loc_sparr1(jj) = parr1(ifacg)
           enddo

           do ii = 1,nneig

              dom_i = commd%neights(ii)
              ini   = commd % bface_size(ii)
              bsize = commd % bface_size(ii+1) - ini

#ifdef MPI_OFF
#else
              if( bsize > 0 ) then
                 !
                 ! Check if I have to send something. dom_i could be a neighbor
                 ! with no common faces with me!
                 !
                 bsize4=int(bsize,4)
                 call MPI_Sendrecv( loc_sparr1(ini:), bsize4,&
                      PAR_REAL,  dom_i, 0_4,     &
                      loc_rparr1(ini:), bsize4,              &
                      PAR_REAL, dom_i, 0_4,      &
                      PAR_COMM_MY_CODE, status, istat )
              end if

#endif
           enddo

           do jj = 1,commd % bface_dim
              ifacg = commd % bface_perm(jj)
              parr1(ifacg) = parr1(ifacg) + loc_rparr1(jj)
           enddo

           call memchk(two,istat,par_memor,'LOC_RPARR1','par_slexch',loc_rparr1)
           deallocate(loc_rparr1,stat=istat)
           if(istat/=0) call memerr(two,'LOC_RPARR1','par_slexch',0_ip)

           call memchk(two,istat,par_memor,'LOC_SPARR1','par_slexch',loc_sparr1)
           deallocate(loc_sparr1,stat=istat)
           if(istat/=0) call memerr(two,'LOC_SPARR1','par_slexch',0_ip)

        else if( pardi >= 1 .and. parki == 5 ) then

           !-------------------------------------------------------------
           !
           ! REAL(PARD1,NPOIN) => REAL(PARD1*NPOIN)
           !
           !-------------------------------------------------------------

           allocate(loc_sparr1(pard1*commd % bface_dim),stat=istat)
           call memchk(zero,istat,par_memor,'LOC_SPARR1','par_slexch',loc_sparr1)

           allocate(loc_rparr1(pard1*commd % bface_dim),stat=istat)
           call memchk(zero,istat,par_memor,'LOC_RPARR1','par_slexch',loc_rparr1)

           do jj = 1,commd % bface_dim
              ifacg = commd % bface_perm(jj)
              do ii = 1,pard1
                 loc_sparr1(pard1*(jj-1)+ii) = parr1(pard1*(ifacg-1)+ii)
              enddo
           enddo

           do ii= 1,nneig

              dom_i = commd%neights(ii)
              ini   = pard1*(commd % bface_size(ii)-1)   + 1
              bsize = pard1*(commd % bface_size(ii+1)-1) - ini + 1

#ifdef MPI_OFF
#else
              if( bsize > 0 ) then
                 !
                 ! Check if I have to send something. dom_i could be a neighbor
                 ! with no common faces with me!
                 !
                 bsize4 = int(bsize,4)
                 call MPI_Sendrecv( loc_sparr1(ini:), bsize4,&
                      PAR_REAL,  dom_i, 0_4,     &
                      loc_rparr1(ini:), bsize4,              &
                      PAR_REAL, dom_i, 0_4,      &
                      PAR_COMM_MY_CODE, status, istat )
              end if
#endif
           end do

           do jj = 1,commd % bface_dim
              ifacg = commd % bface_perm(jj)
              do ii = 1,pard1
                 parr1(pard1*(ifacg-1)+ii) = parr1(pard1*(ifacg-1)+ii) + loc_rparr1(pard1*(jj-1)+ii)
              enddo
           enddo
           call memchk(two,istat,par_memor,'LOC_RPARR1','par_slexch',loc_rparr1)
           deallocate(loc_rparr1,stat=istat)
           if(istat/=0) call memerr(two,'LOC_RPARR1','par_slexch',0_ip)

           call memchk(two,istat,par_memor,'LOC_SPARR1','par_slexch',loc_sparr1)
           deallocate(loc_sparr1,stat=istat)
           if(istat/=0) call memerr(two,'LOC_SPARR1','par_slexch',0_ip)

        else if( pardi == 2 .and. parki == 1 ) then

           call runend('PAR_LGFACE: OPTION NOT AVAILABLE')

        else if( pardi == 2 .and. parki == 2 ) then

           !-------------------------------------------------------------
           !
           ! REAL(PARD1,NPOIN)
           !
           !-------------------------------------------------------------

           allocate(loc_sparr1(pard1*commd % bface_dim),stat=istat)
           call memchk(zero,istat,par_memor,'LOC_SPARR1','par_slexch',loc_sparr1)

           allocate(loc_rparr1(pard1*commd % bface_dim),stat=istat)
           call memchk(zero,istat,par_memor,'LOC_RPARR1','par_slexch',loc_rparr1)

           do jj = 1,commd % bface_dim
              ifacg = commd % bface_perm(jj)
              do ii = 1,pard1
                 loc_sparr1(pard1*(jj-1)+ii) = parr2(ii,ifacg)
              enddo
           enddo

           do ii = 1,nneig

              dom_i = commd%neights(ii)
              ini   = pard1*(commd % bface_size(ii)-1)   + 1
              bsize = pard1*(commd % bface_size(ii+1)-1) - ini + 1

#ifdef MPI_OFF
#else
              if( bsize > 0 ) then
                 !
                 ! Check if I have to send something. dom_i could be a neighbor
                 ! with no common faces with me!
                 !
                 bsize4 = int(bsize,4)
                 call MPI_Sendrecv( loc_sparr1(ini:), bsize4,&
                      PAR_REAL,  dom_i, 0_4,     &
                      loc_rparr1(ini:), bsize4,              &
                      PAR_REAL, dom_i, 0_4,      &
                      PAR_COMM_MY_CODE, status, istat )
              end if
#endif
           enddo

           do jj = 1,commd % bface_dim
              ifacg = commd % bface_perm(jj)
              do ii = 1,pard1
                 parr2(ii,ifacg) = parr2(ii,ifacg) + loc_rparr1(pard1*(jj-1)+ii)
              enddo
           enddo
           call memchk(two,istat,par_memor,'LOC_RPARR1','par_slexch',loc_rparr1)
           deallocate(loc_rparr1,stat=istat)
           if(istat/=0) call memerr(two,'LOC_RPARR1','par_slexch',0_ip)

           call memchk(two,istat,par_memor,'LOC_SPARR1','par_slexch',loc_sparr1)
           deallocate(loc_sparr1,stat=istat)
           if(istat/=0) call memerr(two,'LOC_SPARR1','par_slexch',0_ip)

        else if ( pardi == 1 .and. parki == 4 ) then

           !-------------------------------------------------------------
           !
           ! COMPLEX(NPOIN)
           !
           !-------------------------------------------------------------

           allocate(loc_sparx1(commd % bface_dim),stat=istat)
           call memchk(zero,istat,par_memor,'LOC_SPARX1','par_slexch',loc_sparx1)

           allocate(loc_rparx1(commd % bface_dim),stat=istat)
           call memchk(zero,istat,par_memor,'LOC_RPARX1','par_slexch',loc_rparx1)

           do jj = 1,commd % bface_dim
              ifacg = commd % bface_perm(jj)
              loc_sparx1(jj) = parx1(ifacg)
           enddo

           do ii = 1,nneig

              dom_i = commd%neights(ii)
              ini   = commd % bface_size(ii)
              bsize = commd % bface_size(ii+1) - ini

#ifdef MPI_OFF
#else
              if( bsize > 0 ) then
                 !
                 ! Check if I have to send something. dom_i could be a neighbor
                 ! with no common faces with me!
                 !
                 bsize4 = int(bsize,4)
                 call MPI_Sendrecv( loc_sparx1(ini:), bsize4,&
                      MPI_DOUBLE_COMPLEX, dom_i, 0_4,        &
                      loc_rparx1(ini:), bsize4,              &
                      MPI_DOUBLE_COMPLEX, dom_i, 0_4,        &
                      PAR_COMM_MY_CODE, status, istat )
              end if
#endif
           enddo

           do jj = 1,commd % bface_dim
              ifacg = commd % bface_perm(jj)
              parx1(ifacg) = parx1(ifacg) + loc_rparx1(jj)
           enddo

           call memchk(two,istat,par_memor,'LOC_RPARX1','par_slexch',loc_rparx1)
           deallocate(loc_rparx1,stat=istat)
           if(istat/=0) call memerr(two,'LOC_RPARX1','par_slexch',0_ip)

           call memchk(two,istat,par_memor,'LOC_SPARX1','par_slexch',loc_sparx1)
           deallocate(loc_sparx1,stat=istat)
           if(istat/=0) call memerr(two,'LOC_SPARX1','par_slexch',0_ip)

        else if ( pardi == 1 .and. parki == 7 ) then

           !-------------------------------------------------------------
           !
           ! COMPLEX(PARD1,NPOIN) => COMPLEX(PARD1*NPOIN)
           !
           !-------------------------------------------------------------

           allocate(loc_sparx1(pard1*commd % bface_dim),stat=istat)
           call memchk(zero,istat,par_memor,'LOC_SPARX1','par_slexch',loc_sparx1)

           allocate(loc_rparx1(pard1*commd % bface_dim),stat=istat)
           call memchk(zero,istat,par_memor,'LOC_RPARX1','par_slexch',loc_rparx1)

           do jj = 1,commd % bface_dim

              ji    = pard1 * (jj - 1)
              ifacg = commd % bface_perm(jj)
              poin  = pard1 * (ifacg - 1)
              do ii = 1,pard1

                 loc_sparx1(ji+ii) = parx1(poin+ii)

              enddo

           enddo

           do ii = 1,nneig

              dom_i = commd%neights(ii)
              ini   = pard1 * (commd % bface_size(ii) - 1) + 1
              bsize = pard1 * (commd % bface_size(ii+1) - 1) - ini + 1

#ifdef MPI_OFF
#else
              if( bsize > 0 ) then
                 !
                 ! Check if I have to send something. dom_i could be a neighbor
                 ! with no common faces with me!
                 !
                 bsize4 = int(bsize,4)
                 call MPI_Sendrecv(loc_sparx1(ini:),bsize4,&
                      MPI_DOUBLE_COMPLEX,dom_i,0_4,        &
                      loc_rparx1(ini:),bsize4,             &
                      MPI_DOUBLE_COMPLEX,dom_i,0_4,        &
                      PAR_COMM_MY_CODE,status,istat)
              end if
#endif
           enddo

           do jj = 1,commd % bface_dim

              ji    = pard1 * (jj - 1)
              ifacg = commd % bface_perm(jj)
              poin  = pard1 * (ifacg - 1)
              do ii = 1,pard1

                 parx1(poin+ii) = parx1(poin+ii) + loc_rparx1(ji+ii)

              enddo

           enddo

           call memchk(two,istat,par_memor,'LOC_RPARX1','par_slexch',loc_rparx1)
           deallocate(loc_rparx1,stat=istat)
           if(istat/=0) call memerr(two,'LOC_RPARX1','par_slexch',0_ip)

           call memchk(two,istat,par_memor,'LOC_SPARX1','par_slexch',loc_sparx1)
           deallocate(loc_sparx1,stat=istat)
           if(istat/=0) call memerr(two,'LOC_SPARX1','par_slexch',0_ip)


        end if

     end select

  end if

end subroutine par_lgface
