!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine par_slesca
  !-----------------------------------------------------------------------
  !****f* Parall/par_slesca
  ! NAME
  !    par_slesca
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
  use mod_parall, only : commd,PAR_COMM_MY_CODE
  use mod_parall, only : PAR_INTEGER
  use mod_parall, only : par_memor
  use def_mpi
  implicit none

  integer(ip)              :: ii,bsize,bsire,jj,kk
  integer(4)               :: istat,bsize4,bsire4,dom_i
  real(rp)                 :: time1,time2
  integer(ip), allocatable :: loc_spari1(:),loc_rpari1(:)

  call cputim(time1)

  if( ISLAVE ) then

     allocate(loc_spari1(nneig),stat=istat)
     call memchk(zero,istat,par_memor,'LOC_SPARI1','par_slesca',loc_spari1)

     do ii= 1, nneig

        dom_i  = commd%neights(ii)
        bsize  = 1
        bsize4 = int(bsize,4)

#ifdef MPI_OFF
#else
        !
        ! Send parin to my neighbor DOM_I
        ! Receive parin from my neighbor
        !
        call MPI_Sendrecv( parin(1:1), bsize4 ,&
             PAR_INTEGER,  dom_i, 0_4,         &
             loc_spari1(ii:), bsize4,          &
             PAR_INTEGER, dom_i, 0_4,          &
             PAR_COMM_MY_CODE, status, istat )
#endif
     end do

     do ii= 1, nneig

        allocate(loc_rpari1(max(1_ip,loc_spari1(ii))),stat=istat)
        call memchk(zero,istat,par_memor,'LOC_SPARI1','par_slesca',loc_rpari1)

        dom_i  = commd%neights(ii)
        bsize  = npari
        bsize4 = int(bsize,4)
        bsire  = loc_spari1(ii)
        bsire4 = int(bsire,4)

        if( bsize > 0 .and. bsire > 0 ) then
#ifdef MPI_OFF
#else
           !
           ! Send parin to my neighbor DOM_I
           ! Receive parin from my neighbor
           !
           call MPI_Sendrecv( pari1(1:bsize), bsize4 , &
                PAR_INTEGER,  dom_i, 0_4,              &
                loc_rpari1(1:), bsire4,                &
                PAR_INTEGER, dom_i, 0_4,               &
                PAR_COMM_MY_CODE, status, istat )
#endif
           do kk = 1,bsize
              do jj = 1,bsire
                 if( abs(pari1(kk)) ==  abs(loc_rpari1(jj)) ) then
                    if( kfl_paral > dom_i ) pari1(kk) = -abs(pari1(kk))
                 end if
              end do
           end do
        end if

        call memchk(two,istat,par_memor,'LOC_RPARI1','par_slesca',loc_rpari1)
        deallocate(loc_rpari1,stat=istat)
        if(istat/=0) call memerr(two,'LOC_RPARI1','par_slesca',0_ip)

     end do

     call memchk(two,istat,par_memor,'LOC_SPARI1','par_slesca',loc_spari1)
     deallocate(loc_spari1,stat=istat)
     if(istat/=0) call memerr(two,'LOC_SPARI1','par_slesca',0_ip)


  end if

  call cputim(time2)
  cpu_paral(25)=cpu_paral(25)+time2-time1

end subroutine par_slesca
