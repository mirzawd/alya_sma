!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine par_minxch
  !-----------------------------------------------------------------------
  !****f* Parall/par_minxch
  ! NAME
  !    par_minxch
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
  use mod_parall, only : commd
  use mod_parall, only : PAR_REAL
  use mod_parall, only : PAR_COMM_MY_CODE
  use mod_parall, only : par_memor
  use def_mpi
  implicit none

  integer(ip)           :: ipoin,ii,jj,bsize,ini
  integer(4)            :: istat,bsize4,bsizepard14, dom_i
  real(rp), allocatable :: loc_sparr1(:),   loc_rparr1(:)
  real(rp), allocatable :: loc_sparr2(:,:), loc_rparr2(:,:)


  if(kfl_paral>0) then

     if(party==1) then
        !
        ! Element
        !
     else if(party==2) then
        !
        ! Boundary
        !
     else if(party==3) then
        !
        ! Node
        !
        if(pardi==1.and.parki==1) then

        else if(pardi==1.and.parki==2) then
           !
           ! commd%bound_perm(jj):                        permutation array
           ! loc_sparr1:                                  my local values
           ! loc_rparr1:                                  values given by neighbor ii
           ! commd%bound_dim:                             size of communication array
           ! nneig:                                       number of neighbors that share same group
           ! commd%neights(ii):                           number of subdomain ii
           ! commd%bound_size(ii):                        where my local arrays sart to exchange with ii
           ! commd%bound_size(ii+1)-commd%bound_size(ii): number of groups to exchange with ii
           !
           allocate(loc_sparr1(commd%bound_dim),stat=istat)
           call memchk(zero,istat,par_memor,'LOC_SPARR1','par_minxch',loc_sparr1)

           allocate(loc_rparr1(commd%bound_dim),stat=istat)
           call memchk(zero,istat,par_memor,'LOC_RPARR1','par_minxch',loc_rparr1)

           do jj= 1, commd%bound_dim
              ipoin = commd%bound_perm(jj)
              loc_sparr1(jj) = parr1(ipoin)
           enddo

           do ii= 1, nneig
              dom_i = commd%neights(ii)

              ini   = commd%bound_size(ii)
              bsize = commd%bound_size(ii+1) - ini

#ifdef MPI_OFF
#else
              bsize4=int(bsize,4)
              call MPI_Sendrecv( loc_sparr1(ini:), bsize4,&
                   PAR_REAL,  dom_i, 0,       &
                   loc_rparr1(ini:), bsize4,              &
                   PAR_REAL, dom_i, 0,        &
                   PAR_COMM_MY_CODE, status, istat )
#endif
           enddo
           do jj= 1, commd%bound_dim
              ipoin = commd%bound_perm(jj)
              parr1(ipoin) = min(parr1(ipoin),loc_rparr1(jj))
           enddo

           call memchk(two,istat,par_memor,'LOC_RPARR1','par_minxch',loc_rparr1)
           deallocate(loc_rparr1,stat=istat)
           if(istat/=0) call memerr(two,'LOC_RPARR1','par_minxch',0_ip)

           call memchk(two,istat,par_memor,'LOC_SPARR1','par_minxch',loc_sparr1)
           deallocate(loc_sparr1,stat=istat)
           if(istat/=0) call memerr(two,'LOC_SPARR1','par_minxch',0_ip)

        else if(pardi==1.and.parki==5) then

           allocate(loc_sparr1(pard1*commd%bound_dim),stat=istat)
           call memchk(zero,istat,par_memor,'LOC_SPARR1','par_minxch',loc_sparr1)

           allocate(loc_rparr1(pard1*commd%bound_dim),stat=istat)
           call memchk(zero,istat,par_memor,'LOC_RPARR1','par_minxch',loc_rparr1)

           do jj= 1, commd%bound_dim
              ipoin = commd%bound_perm(jj)
              do ii= 1, pard1
                 loc_sparr1(pard1*(jj-1)+ii) = parr1(pard1*(ipoin-1)+ii)
              enddo
           enddo

           do ii= 1, nneig
              dom_i = commd%neights(ii)

              ini   = pard1*(commd%bound_size(ii)-1)   + 1
              bsize = pard1*(commd%bound_size(ii+1)-1) - ini + 1

#ifdef MPI_OFF
#else
              bsize4 = int(bsize,4)
              call MPI_Sendrecv( loc_sparr1(ini:), bsize4,&
                   PAR_REAL,  dom_i, 0,       &
                   loc_rparr1(ini:), bsize4,              &
                   PAR_REAL, dom_i, 0,        &
                   PAR_COMM_MY_CODE, status, istat )
#endif
           enddo

           do jj= 1, commd%bound_dim
              ipoin = commd%bound_perm(jj)
              do ii= 1, pard1
                 parr1(pard1*(ipoin-1)+ii) = min(parr1(pard1*(ipoin-1)+ii) ,loc_rparr1(pard1*(jj-1)+ii))
              enddo
           enddo
           call memchk(two,istat,par_memor,'LOC_RPARR1','par_minxch',loc_rparr1)
           deallocate(loc_rparr1,stat=istat)
           if(istat/=0) call memerr(two,'LOC_RPARR1','par_minxch',0_ip)

           call memchk(two,istat,par_memor,'LOC_SPARR1','par_minxch',loc_sparr1)
           deallocate(loc_sparr1,stat=istat)
           if(istat/=0) call memerr(two,'LOC_SPARR1','par_minxch',0_ip)

        else if(pardi==2.and.parki==1) then

        else if(pardi==2.and.parki==2) then
           allocate(loc_sparr2(pard1,commd%bound_dim),stat=istat)
           call memchk(zero,istat,par_memor,'LOC_SPARR2','par_minxch',loc_sparr2)

           allocate(loc_rparr2(pard1,commd%bound_dim),stat=istat)
           call memchk(zero,istat,par_memor,'LOC_RPARR2','par_minxch',loc_rparr2)

           do jj= 1, commd%bound_dim
              ipoin = commd%bound_perm(jj)
              loc_sparr2(1:pard1,jj) = parr2(1:pard1,ipoin)
           enddo

           do ii= 1, nneig
              dom_i = commd%neights(ii)

              ini   = commd%bound_size(ii)
              bsize = commd%bound_size(ii+1) - ini

#ifdef MPI_OFF
#else
              bsize4      = int(bsize,4)
              bsizepard14 = bsize4*int(pard1,4)
              call MPI_Sendrecv( loc_sparr2(1:,ini), bsizepard14,&
                   PAR_REAL,  dom_i, 0_4,            &
                   loc_rparr2(1:,ini), bsizepard14,              &
                   PAR_REAL, dom_i, 0_4,             &
                   PAR_COMM_MY_CODE, status, istat )
#endif
           enddo

           do jj= 1, commd%bound_dim
              ipoin = commd%bound_perm(jj)
              parr2(1:pard1,ipoin) = min(parr2(1:pard1,ipoin),loc_rparr2(1:pard1,jj))
           enddo

           call memchk(two,istat,par_memor,'LOC_RPARR2','par_minxch',loc_rparr2)
           deallocate(loc_rparr2,stat=istat)
           if(istat/=0) call memerr(two,'LOC_RPARR2','par_minxch',0_ip)

           call memchk(two,istat,par_memor,'LOC_SPARR2','par_minxch',loc_sparr2)
           deallocate(loc_sparr2,stat=istat)
           if(istat/=0) call memerr(two,'LOC_SPARR2','par_minxch',0_ip)

        else if(pardi==1.and.parki==6) then

           allocate(loc_sparr1(pard1*commd%bound_dim),stat=istat)
           call memchk(zero,istat,par_memor,'LOC_SPARR1','par_minxch',loc_sparr1)

           allocate(loc_rparr1(pard1*commd%bound_dim),stat=istat)
           call memchk(zero,istat,par_memor,'LOC_RPARR1','par_minxch',loc_rparr1)

           do jj= 1, commd%bound_dim
              ipoin = commd%bound_perm(jj)
              do ii= 1, pard1
                 loc_sparr1(pard1*(jj-1)+ii) = parr2(ii,ipoin)
              enddo
           enddo

           do ii= 1, nneig
              dom_i = commd%neights(ii)

              ini   = pard1*(commd%bound_size(ii)-1)   + 1
              bsize = pard1*(commd%bound_size(ii+1)-1) - ini + 1

#ifdef MPI_OFF
#else
              bsize4 = int(bsize,4)
              call MPI_Sendrecv( loc_sparr1(ini:), bsize4,&
                   PAR_REAL,  dom_i, 0_4,     &
                   loc_rparr1(ini:), bsize4,              &
                   PAR_REAL, dom_i, 0_4,      &
                   PAR_COMM_MY_CODE, status, istat )
#endif
           enddo

           do jj= 1, commd%bound_dim
              ipoin = commd%bound_perm(jj)
              do ii= 1, pard1
                 parr2(ii,ipoin) = min(parr2(ii,ipoin),loc_rparr1(pard1*(jj-1)+ii))
              enddo
           enddo
           call memchk(two,istat,par_memor,'LOC_RPARR1','par_minxch',loc_rparr1)
           deallocate(loc_rparr1,stat=istat)
           if(istat/=0) call memerr(two,'LOC_RPARR1','par_minxch',0_ip)

           call memchk(two,istat,par_memor,'LOC_SPARR1','par_minxch',loc_sparr1)
           deallocate(loc_sparr1,stat=istat)
           if(istat/=0) call memerr(two,'LOC_SPARR1','par_minxch',0_ip)
        end if
        !
        ! Boundary node: REAL ARRAY(NBOPO)
        !
     else if(party==4) then
     end if
  endif


end subroutine par_minxch
