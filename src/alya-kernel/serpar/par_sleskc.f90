!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine par_sleskc()
  !-----------------------------------------------------------------------
  !****f* Parall/par_sleskc
  ! NAME
  !    par_sleskc
  ! DESCRIPTION
  !    This subroutine exchange arrays between master and slaves
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
  use mod_parall, only : commd,PAR_COMM_MY_CODE
  use mod_parall, only : PAR_INTEGER
  use mod_parall, only : PAR_REAL
  use mod_parall, only : par_memor
  use def_mpi
  implicit none

  integer(ip)              :: ipoin,ii,jj,bsize,ini,ibopo,kk,mm
  integer                  :: istat,bsize4,dom_i
  integer(ip), allocatable :: loc_spari1(:),loc_rpari1(:)
  integer(ip), allocatable :: loc_sparii(:),loc_rparii(:)
  real(rp),    allocatable :: loc_sparr1(:),loc_rparr1(:)

  if( ISLAVE ) then

     !-------------------------------------------------------------------
     !
     ! Take from the neighbor's who has the mastr's node:
     ! SKCOS
     ! KFL_GEONO
     ! In fact, they can be different because of round of errors
     ! or ordering....
     !
     !-------------------------------------------------------------------
     !
     ! LOC_PARII = 1 when my neighbor owns the node
     !
     allocate(loc_sparii(commd%bound_dim),stat=istat)
     call memchk(zero,istat,par_memor,'LOC_SPARII','par_slexch',loc_sparii)

     allocate(loc_rparii(commd%bound_dim),stat=istat)
     call memchk(zero,istat,par_memor,'LOC_RPARII','par_slexch',loc_rparii)

     do jj= 1, commd%bound_dim
        ipoin = commd%bound_perm(jj)
        if( ipoin >= npoi2 .and. ipoin <= npoi3 ) then
           loc_sparii(jj) = 1                             ! I own this node
        end if
     end do

     do ii= 1, nneig
        dom_i = commd%neights(ii)

        ini   = commd%bound_size(ii)
        bsize = commd%bound_size(ii+1) - ini

#ifdef MPI_OFF
#else
        bsize4=int(bsize,4)
        call MPI_Sendrecv( loc_sparii(ini:), bsize4,&
             PAR_INTEGER,  dom_i, 0_4,     &
             loc_rparii(ini:), bsize4,              &
             PAR_INTEGER, dom_i, 0_4,      &
             PAR_COMM_MY_CODE, status, istat )
#endif
     end do
     !
     ! SKCOS
     !
     pard1 = ndime*ndime

     allocate(loc_sparr1(pard1*commd%bound_dim),stat=istat)
     call memchk(zero,istat,par_memor,'LOC_SPARR1','par_sleskc',loc_sparr1)

     allocate(loc_rparr1(pard1*commd%bound_dim),stat=istat)
     call memchk(zero,istat,par_memor,'LOC_RPARR1','par_sleskc',loc_rparr1)

     do jj = 1,commd%bound_dim
        ipoin = commd%bound_perm(jj)
        ibopo = lpoty(ipoin)
        if( ibopo > 0 ) then
           mm = pard1*(jj-1) 
           do ii = 1,ndime
              do kk = 1,ndime
                 mm = mm + 1
                 loc_sparr1(mm) = skcos(kk,ii,ibopo)
              end do
           end do
        end if
     end do

     do ii= 1, nneig
        dom_i = commd%neights(ii)

        ini   = pard1*(commd%bound_size(ii)-1)   + 1
        bsize = pard1*(commd%bound_size(ii+1)-1) - ini + 1

#ifdef MPI_OFF
#else
        bsize4 = int(bsize,4)
        call MPI_Sendrecv( loc_sparr1(ini:), bsize4,&
             PAR_REAL,  dom_i, 0,                  &
             loc_rparr1(ini:), bsize4,              &
             PAR_REAL, dom_i, 0,                   &
             PAR_COMM_MY_CODE, status, istat )
#endif
     end do

     do jj = 1,commd%bound_dim
        if( loc_rparii(jj) == 1 ) then     ! Node is owned by this neighbor
           ipoin = commd%bound_perm(jj)
              ibopo = lpoty(ipoin)
              if( ibopo > 0 ) then
                 mm = pard1*(jj-1)
                 do ii = 1,ndime
                    do kk = 1,ndime
                       mm = mm + 1
                       skcos(kk,ii,ibopo) = loc_rparr1(mm) 
                    end do
                 end do
              end if
        end if
     end do

     call memchk(two,istat,par_memor,'LOC_RPARR1','par_sleskc',loc_rparr1)
     deallocate(loc_rparr1,stat=istat)
     if(istat/=0) call memerr(two,'LOC_RPARR1','par_sleskc',0_ip)

     call memchk(two,istat,par_memor,'LOC_SPARR1','par_sleskc',loc_sparr1)
     deallocate(loc_sparr1,stat=istat)
     if(istat/=0) call memerr(two,'LOC_SPARR1','par_sleskc',0_ip)
     !
     ! KFL_GEONO
     !
     allocate(loc_spari1(commd%bound_dim),stat=istat)
     call memchk(zero,istat,par_memor,'LOC_SPARI1','par_slexch',loc_spari1)

     allocate(loc_rpari1(commd%bound_dim),stat=istat)
     call memchk(zero,istat,par_memor,'LOC_RPARI1','par_slexch',loc_rpari1)

     do jj= 1, commd%bound_dim
        ipoin = commd%bound_perm(jj)
        ibopo = lpoty(ipoin)
        if( ibopo > 0 ) then
           loc_spari1(jj) = kfl_geono(ibopo)
        end if
     enddo

     do ii= 1, nneig
        dom_i = commd%neights(ii)

        ini   = commd%bound_size(ii)
        bsize = commd%bound_size(ii+1) - ini

#ifdef MPI_OFF
#else
        bsize4=int(bsize,4)
        call MPI_Sendrecv( loc_spari1(ini:), bsize4,&
             PAR_INTEGER,  dom_i, 0,                &
             loc_rpari1(ini:), bsize4,              &
             PAR_INTEGER, dom_i, 0,                 &
             PAR_COMM_MY_CODE, status, istat )
#endif
     enddo

     do jj= 1, commd%bound_dim
        if( loc_rparii(jj) == 1 ) then
           ipoin = commd%bound_perm(jj)
           ibopo = lpoty(ipoin)
           if( ibopo > 0 ) then
              kfl_geono(ibopo) = loc_spari1(jj) 
           end if
        end if
     end do

     call memchk(two,istat,par_memor,'LOC_RPARI1','par_slexch',loc_rpari1)
     deallocate(loc_rpari1,stat=istat)
     if(istat/=0) call memerr(two,'LOC_RPARI1','par_slexch',0_ip)

     call memchk(two,istat,par_memor,'LOC_SPARI1','par_slexch',loc_spari1)
     deallocate(loc_spari1,stat=istat)
     if(istat/=0) call memerr(two,'LOC_SPARI1','par_slexch',0_ip)

     call memchk(two,istat,par_memor,'LOC_SPARII','par_slexch',loc_sparii)
     deallocate(loc_sparii,stat=istat)
     if(istat/=0) call memerr(two,'LOC_SPARII','par_slexch',0_ip)

     call memchk(two,istat,par_memor,'LOC_RPARII','par_slexch',loc_rparii)
     deallocate(loc_rparii,stat=istat)
     if(istat/=0) call memerr(two,'LOC_RPARII','par_slexch',0_ip)

  end if

end subroutine par_sleskc
