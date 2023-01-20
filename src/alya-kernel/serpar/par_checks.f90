!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine par_checks()
  !-----------------------------------------------------------------------
  !****f* Parall/par_slexch
  ! NAME
  !    par_slexch
  ! DESCRIPTION
  !    This subroutine performs some checks
  ! USED BY
  !    Parall
  !***
  !-----------------------------------------------------------------------
  use def_master
  use def_parame
  use def_parall
  use def_domain
  use mod_memchk
  use mod_parall, only    :  PAR_COMM_MY_CODE,PAR_INTEGER,commd
  use mod_parall, only    :  PAR_REAL
  use mod_parall, only : par_memor
  use def_mpi
  implicit none

  integer(ip)              :: ipoin,ii,jj,bsize,ini,jdime,jdofn
  integer(4)               :: istat,bsize4,dom_i
  integer(ip), allocatable :: loc_spari1(:),loc_rpari1(:)
  real(rp),    allocatable :: loc_sparr1(:),loc_rparr1(:)

  if( ISLAVE ) then
     !
     ! Check if OWN boundary nodes are not repeated
     !
     allocate(loc_spari1(commd%bound_dim),stat=istat)
     call memchk(zero,istat,par_memor,'LOC_SPARI1','par_slexch',loc_spari1)

     allocate(loc_rpari1(commd%bound_dim),stat=istat)
     call memchk(zero,istat,par_memor,'LOC_RPARI1','par_slexch',loc_rpari1)

     call memgen(1_ip,npoin,0_ip)

     do jj= 1, commd%bound_dim
        ipoin = commd%bound_perm(jj)
        if( ipoin >= npoi2 .and. ipoin <= npoi3 ) then
           loc_spari1(jj) = 1
           gisca(ipoin)   = 1
        end if
     end do

#ifdef MPI_OFF
#else
     do ii = 1, nneig

        dom_i  = commd%neights(ii)
        ini    = commd%bound_size(ii)
        bsize  = commd%bound_size(ii+1) - ini
        bsize4 = int(bsize,4)

        call MPI_Sendrecv( &
             loc_spari1(ini:), bsize4,     &
             PAR_INTEGER,  dom_i, 0_4,     &
             loc_rpari1(ini:), bsize4,     &
             PAR_INTEGER,  dom_i, 0_4,     &
             PAR_COMM_MY_CODE, status, istat )
     end do
#endif

     do jj = 1, commd%bound_dim
        ipoin = commd%bound_perm(jj)
        gisca(ipoin) = gisca(ipoin) + loc_rpari1(jj)
     end do

     do ii = 1, nneig
        dom_i  = commd%neights(ii)
        do jj = commd%bound_size(ii),commd%bound_size(ii+1)-1
           ipoin = commd%bound_perm(jj)
           if( gisca(ipoin) > 1 ) then
              write(*,'(a,4(1x,i6))') 'OWN BOUNDARY NODE REPEATED=',kfl_paral,dom_i,ipoin,lninv_loc(ipoin)
           else if( gisca(ipoin) == 0 ) then
              write(*,'(a,4(1x,i6))') 'NO OWN BOUNDARY OWNER=     ',kfl_paral,dom_i,ipoin,lninv_loc(ipoin)
           end if
        end do
     end do

     call memgen(3_ip,npoin,0_ip)

     call memchk(two,istat,par_memor,'LOC_RPARI1','par_slexch',loc_rpari1)
     deallocate(loc_rpari1,stat=istat)
     if(istat/=0) call memerr(two,'LOC_RPARI1','par_slexch',0_ip)

     call memchk(two,istat,par_memor,'LOC_SPARI1','par_slexch',loc_spari1)
     deallocate(loc_spari1,stat=istat)
     if(istat/=0) call memerr(two,'LOC_SPARI1','par_slexch',0_ip)
     !
     ! Check if boundary node coordinates coincide
     !
     allocate(loc_sparr1(ndime*commd%bound_dim),stat=istat)
     call memchk(zero,istat,par_memor,'LOC_SPARR1','par_slexch',loc_sparr1)

     allocate(loc_rparr1(ndime*commd%bound_dim),stat=istat)
     call memchk(zero,istat,par_memor,'LOC_RPARR1','par_slexch',loc_rparr1)

     jdofn = 0
     do jj = 1,commd%bound_dim
        ipoin = commd%bound_perm(jj)
        do jdime = 1,ndime
           jdofn = jdofn + 1
           loc_sparr1(jdofn) = coord(jdime,ipoin)
        end do
     end do

#ifdef MPI_OFF
#else
     do ii = 1 , nneig

        dom_i  = commd%neights(ii)
        ini    = (commd%bound_size(ii)-1)*ndime+1
        bsize  = (commd%bound_size(ii+1)-1)*ndime+1 - ini
        bsize4 = int(bsize,4)

        call MPI_Sendrecv( &
             loc_sparr1(ini:), bsize4,              &
             PAR_REAL,  dom_i, 0_4,     &
             loc_rparr1(ini:), bsize4,              &
             PAR_REAL, dom_i, 0_4,      &
             PAR_COMM_MY_CODE, status, istat )
     end do
#endif

     jdofn = 0
     do jj = 1 , commd%bound_dim
        ipoin = commd%bound_perm(jj)
        ii = 0
        do jdime = 1,ndime
           jdofn = jdofn + 1
           if( abs( coord(jdime,ipoin) - loc_rparr1(jdofn)) / (coord(jdime,ipoin) + zeror) > 1.0e-15_rp ) ii = ii + 1
        end do
        if( ii /= 0 ) then
           print*,''
           print*,'PAR_CHECKS: BAD COORDINATES OF BOUNDARY NODE='
           print*,'KFL_PARAL    =',kfl_paral
           print*,'LNINV NODE 1 =',lninv_loc(ipoin)
           print*,'COORD NODE 1 =',coord(1:ndime,ipoin)
           print*,'COORD NODE 2 =',loc_rparr1( (jj-1)*ndime+1:(jj-1)*ndime+ndime)
           print*,'DIFF COORD   =',abs(coord(1:ndime,ipoin) - loc_rparr1( (jj-1)*ndime+1:(jj-1)*ndime+ndime))
        end if
     end do

     call memchk(two,istat,par_memor,'LOC_RPARR1','par_slexch',loc_rparr1)
     deallocate(loc_rparr1,stat=istat)
     if(istat/=0) call memerr(two,'LOC_RPARR1','par_slexch',0_ip)

     call memchk(two,istat,par_memor,'LOC_SPARR1','par_slexch',loc_sparr1)
     deallocate(loc_sparr1,stat=istat)
     if(istat/=0) call memerr(two,'LOC_SPARR1','par_slexch',0_ip)

  end if

end subroutine par_checks
