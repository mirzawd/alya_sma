!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine par_slesch()
  !-----------------------------------------------------------------------
  !****f* Parall/par_slesch
  ! NAME
  !    par_slexch
  ! DESCRIPTION
  !    This subroutine exchange arrays between slaves
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
  use mod_parall, only : PAR_REAL
  use mod_parall, only : par_memor
  use def_mpi
  implicit none

  integer(ip)              :: ipoin,ii,jj,bsize,ini
  integer(4)               :: istat,bsize4,dom_i
  real(rp)                 :: time1,time2
  real(rp),    allocatable :: loc_sparr1(:),   loc_rparr1(:)

  call cputim(time1)

  if(kfl_paral>0) then

     !-------------------------------------------------------------
     !
     ! REAL(NPOIN)
     !
     !-------------------------------------------------------------

     allocate(loc_sparr1(commd%bound_dim),stat=istat)
     call memchk(zero,istat,par_memor,'LOC_SPARR1','par_slexch',loc_sparr1)

     allocate(loc_rparr1(commd%bound_dim),stat=istat)
     call memchk(zero,istat,par_memor,'LOC_RPARR1','par_slexch',loc_rparr1)

     do jj= 1, commd%bound_dim
        ipoin = commd%bound_perm(jj)
        loc_sparr1(jj) = parr1(ipoin-npoi1)
     enddo

     do ii= 1, nneig
        dom_i = commd%neights(ii)

        ini   = commd%bound_size(ii)
        bsize = commd%bound_size(ii+1) - ini

#ifdef MPI_OFF
#else
        bsize4=int(bsize,4)
        call MPI_Sendrecv( loc_sparr1(ini:), bsize4,&
             PAR_REAL,  dom_i, 0_4,     &
             loc_rparr1(ini:), bsize4,              &
             PAR_REAL, dom_i, 0_4,      &
             PAR_COMM_MY_CODE, status, istat )

#endif
     enddo

     do jj= 1, commd%bound_dim
        ipoin = commd%bound_perm(jj)
        parr1(ipoin-npoi1) = parr1(ipoin-npoi1) + loc_rparr1(jj)
     enddo

     call memchk(two,istat,par_memor,'LOC_RPARR1','par_slexch',loc_rparr1)
     deallocate(loc_rparr1,stat=istat)
     if(istat/=0) call memerr(two,'LOC_RPARR1','par_slexch',0_ip)

     call memchk(two,istat,par_memor,'LOC_SPARR1','par_slexch',loc_sparr1)
     deallocate(loc_sparr1,stat=istat)
     if(istat/=0) call memerr(two,'LOC_SPARR1','par_slexch',0_ip)

  endif

  call cputim(time2)
  cpu_paral(25)=cpu_paral(25)+time2-time1

end subroutine par_slesch
