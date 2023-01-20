!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine par_cregrp(itask)
  !-----------------------------------------------------------------------
  !****f* Parall/par_cregrp
  ! NAME
  !    par_cregrp
  ! DESCRIPTION
  !    This routine computes the groups for the deflated CG.
  ! USED BY
  !    domain
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_solver
  use def_parall
  use mod_memchk
  use mod_postpr
  use mod_parall, only    :  PAR_COMM_MY_CODE,PAR_INTEGER
  use mod_parall, only    :  PAR_REAL
  use mod_parall, only    :  commd
  use mod_parall, only    : par_memor
  use def_mpi
  implicit none

  integer(ip), intent(in)  :: itask
  integer(ip)              :: ipoin,ii,jj,bsize,ini
  integer(ip)              :: igrou,pgrou,ngrou
  integer(ip), pointer     :: ngros(:)
  integer                  :: istat,bsize4,dom_i
  integer(ip), allocatable :: loc_spari1(:),loc_rpari1(:)
  real(rp),    allocatable :: loc_sparr1(:),loc_rparr1(:)

  select case (itask)

  case(1_ip)
     !
     ! Only one subdomain has the number
     !
     if( IMASTER ) return

     allocate(loc_spari1(commd%bound_dim),stat=istat)
     call memchk(zero,istat,par_memor,'LOC_SPARI1','par_cregrp',loc_spari1)

     allocate(loc_rpari1(commd%bound_dim),stat=istat)
     call memchk(zero,istat,par_memor,'LOC_RPARI1','par_cregrp',loc_rpari1)

     do jj= 1, commd%bound_dim
        ipoin = commd%bound_perm(jj)
        loc_spari1(jj) = pari1(ipoin)
     enddo

     do ii= 1, nneig
        dom_i = commd%neights(ii)

        ini   = commd%bound_size(ii)
        bsize = commd%bound_size(ii+1) - ini

#ifdef MPI_OFF
#else
        bsize4=int(bsize,4)
        call MPI_Sendrecv( loc_spari1(ini:), bsize4,&
             PAR_INTEGER,  dom_i, 0,     &
             loc_rpari1(ini:), bsize4,              &
             PAR_INTEGER, dom_i, 0,      &
             PAR_COMM_MY_CODE, status, istat )
#endif
     enddo

     do jj = 1, commd%bound_dim
        ipoin = commd%bound_perm(jj)
        pari1(ipoin) = min(pari1(ipoin),loc_rpari1(jj))
     enddo

     call memchk(two,istat,par_memor,'LOC_RPARI1','par_cregrp',loc_rpari1)
     deallocate(loc_rpari1,stat=istat)
     if(istat/=0) call memerr(two,'LOC_RPARI1','par_cregrp',0_ip)

     call memchk(two,istat,par_memor,'LOC_SPARI1','par_cregrp',loc_spari1)
     deallocate(loc_spari1,stat=istat)
     if(istat/=0) call memerr(two,'LOC_SPARI1','par_cregrp',0_ip)

  case(2_ip)
     !
     ! Determine who's propietary of the group
     !
     if( IMASTER ) return

     ngrou = solve_sol(1) % ngrou
     allocate(lessl(2,ngrou),stat=istat)
     call memchk(zero,istat,par_memor,'LESSL','par_cregrp',lessl)
 
     do ipoin = 1,npoin
        igrou = solve_sol(1) % lgrou(ipoin)
        if( igrou > 0 ) lessl(1,igrou) = 1
        if( igrou > 0 ) lessl(2,igrou) = 1
     end do

     do ii = 1,nneig
        dom_i = commd % neights(ii)
        jj    = commd % bound_size(ii)
        do while( jj <= commd % bound_size(ii+1) - 1 )
           ipoin = commd % bound_perm(jj)
           igrou = solve_sol(1) % lgrou(ipoin)
           if( igrou > 0 ) lessl(1,igrou) = 0
           if( igrou > 0 .and. kfl_paral > dom_i ) lessl(2,igrou) = 0
           jj    = jj + 1
        end do
     end do

  case(3_ip)
     !
     ! 
     !
     if( IMASTER ) return

     ngrou = solve_sol(1) % ngrou
     !allocate( lgros(nneig), stat = istat )
     !allocate( lgros(nneig), stat = istat )
     allocate( ngros(ngrou), stat = istat )

     allocate(loc_sparr1(ngrou),stat=istat)
     call memchk(zero,istat,par_memor,'LOC_SPARR1','par_cregrp',loc_sparr1)
     
     allocate(loc_rparr1(ngrou),stat=istat)
     call memchk(zero,istat,par_memor,'LOC_RPARR1','par_cregrp',loc_rparr1)

     do ii = 1,nneig

        do igrou = 1,ngrou
           ngros(igrou) = 0
        end do
        dom_i = commd%neights(ii)
        jj    = commd%bound_size(ii)
        do while( jj <= commd%bound_size(ii+1) - 1 )
           ipoin = commd%bound_perm(jj)
           igrou = solve_sol(1) % lgrou(ipoin)
           if( igrou > 0 ) then
              ngros(igrou) = ngros(igrou) + 1
           end if
           jj    = jj + 1
        end do

        pgrou = 0
        do igrou = 1,ngrou
           if( ngros(igrou) > 1 ) pgrou = pgrou + 1
        end do

        if( pgrou > 0 ) then
           continue
        else
           write(*,*) 'PAR_CERGRP: ERROOOOOOOOOOOOOOOOOORR'
        end if

        pgrou = 0
        do igrou = 1,ngrou
           if( ngros(igrou) > 1 ) then
              pgrou = pgrou + 1
              loc_sparr1(pgrou) = parre(igrou)
              loc_rparr1(pgrou) = 0.0_rp
              !lgros(ii) % l(pgrou) = igrou
           end if
        end do
       
        dom_i = commd%neights(ii)
        ini   = 1
        bsize = pgrou
        
#ifdef MPI_OFF
#else
        bsize4 = int(bsize,4)
        call MPI_Sendrecv( loc_sparr1(ini:), bsize4,&
             PAR_REAL,  dom_i, 0,     &
             loc_rparr1(ini:), bsize4,              &
             PAR_REAL, dom_i, 0,      &
             PAR_COMM_MY_CODE, status, istat )
#endif

     end do
        pgrou = 0
        do igrou = 1,ngrou
           if( ngros(igrou) > 1 ) then
              pgrou = pgrou + 1
              parre(igrou) = parre(igrou) + loc_rparr1(pgrou)
              !lgros(ii) % l(pgrou) = igrou
           end if
        end do

    
     !deallocate( lgros, stat = istat )
     deallocate( ngros, stat = istat )

     call memchk(two,istat,par_memor,'LOC_RPARR1','par_cregrp',loc_rparr1)
     deallocate(loc_rparr1,stat=istat)
     if(istat/=0) call memerr(two,'LOC_RPARR1','par_cregrp',0_ip)
     
     call memchk(two,istat,par_memor,'LOC_SPARR1','par_cregrp',loc_sparr1)
     deallocate(loc_sparr1,stat=istat)
     if(istat/=0) call memerr(two,'LOC_SPARR1','par_cregrp',0_ip)

  end select

end subroutine par_cregrp
