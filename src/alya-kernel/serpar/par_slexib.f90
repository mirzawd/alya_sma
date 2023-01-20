!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine par_slexib(itask)
  !-----------------------------------------------------------------------
  !****f* Parall/par_slexib
  ! NAME
  !    par_slexib
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
  use mod_parall, only : PAR_REAL
  use mod_parall, only : par_memor
  use def_mpi
  implicit none

  integer(ip), intent(in)  :: itask 
  integer(ip)              :: ipoin,ii,jj,bsize,ini,kk,ll
  integer(ip)              :: idime,idofi,idofj,nneno
  integer                  :: istat,bsize4,dom_i
  real(rp)                 :: time1,time2
  real(rp),    pointer     :: loc_sparr1(:),loc_rparr1(:)
  integer(ip), pointer     :: loc_spari1(:),loc_rpari1(:)

  call cputim(time1)

  if( ISLAVE ) then

     if( itask == 1 ) then

        !-------------------------------------------------------------
        !
        ! FIND WHICH SUBDOMAINS HAVE THE LOVERS
        !
        !-------------------------------------------------------------

        allocate(loc_sparr1(commd%bound_dim),stat=istat)
        call memchk(zero,istat,par_memor,'LOC_SPARR1','par_slexib',loc_sparr1)

        allocate(loc_rparr1(commd%bound_dim),stat=istat)
        call memchk(zero,istat,par_memor,'LOC_RPARR1','par_slexib',loc_rparr1)
        !
        ! Gather
        !
        do jj= 1, commd%bound_dim
           ipoin = commd%bound_perm(jj)
           loc_sparr1(jj) = parr1(ipoin)
        end do

        do ii= 1, nneig
           dom_i = commd%neights(ii)

           ini   = commd%bound_size(ii)
           bsize = commd%bound_size(ii+1) - ini

#ifdef MPI_OFF
#else
           bsize4=int(bsize,4)
           call MPI_Sendrecv( loc_sparr1(ini:), bsize4, &
                PAR_REAL,  dom_i, 0_4,      &
                loc_rparr1(ini:), bsize4,               &
                PAR_REAL, dom_i, 0_4,       &
                PAR_COMM_MY_CODE, status, istat )
#endif
        end do

        ini = 0 ! Trade off
        do ii= 1, nneig
           dom_i = commd%neights(ii)
           
           do kk = 1,commd%bound_size(ii+1)-commd%bound_size(ii)
              jj    = kk + ini

              ipoin = commd%bound_perm(jj)

              if( pari1(ipoin) /= 0 ) then
                 !
                 ! I am a travesty
                 !
                 if( loc_rparr1(jj) < parr1(ipoin) ) then
                    !
                    ! If found a travesty with a greater lover
                    !
                    pari1(ipoin) = -dom_i  
                    parr1(ipoin) = loc_rparr1(jj)

                 else if( loc_rparr1(jj) == parr1(ipoin) ) then
                    !
                    ! We have both a nice lover
                    !
                    if( pari1(ipoin) > 0 ) then
                       if( kfl_paral < dom_i ) pari1(ipoin) = -dom_i
                    else
                       pari1(ipoin) = min(-int(dom_i,ip),pari1(ipoin))
                    end if

                 end if

              end if

           end do
           ini = ini + commd%bound_size(ii+1)-commd%bound_size(ii)
        end do

        call memchk(two,istat,par_memor,'LOC_RPARR1','par_slexib',loc_rparr1)
        deallocate(loc_rparr1,stat=istat)
        if(istat/=0) call memerr(two,'LOC_RPARR1','par_slexib',0_ip)

        call memchk(two,istat,par_memor,'LOC_SPARR1','par_slexib',loc_sparr1)
        deallocate(loc_sparr1,stat=istat)
        if(istat/=0) call memerr(two,'LOC_SPARR1','par_slexib',0_ip)

     else if( itask == 2 ) then

        !-------------------------------------------------------------
        !
        ! FIND WHICH SUBDOMAINS HAVE THE LOVERS
        !
        !-------------------------------------------------------------

        allocate(loc_sparr1(pard1*commd%bound_dim),stat=istat)
        call memchk(zero,istat,par_memor,'LOC_SPARR1','par_slexib',loc_sparr1)

        allocate(loc_rparr1(pard1*commd%bound_dim),stat=istat)
        call memchk(zero,istat,par_memor,'LOC_RPARR1','par_slexib',loc_rparr1)

        do jj= 1, commd%bound_dim
           ipoin = commd%bound_perm(jj)
           idofi = (ipoin-1)*pard1
           idofj = (jj-1)*pard1
           do idime = 1, pard1
              idofi = idofi + 1
              idofj = idofj + 1
              loc_sparr1(idofj) = parr1(idofi)
           end do
        end do

        ini = 1
        do ii= 1, nneig

           dom_i = commd%neights(ii)
           bsize = (commd%bound_size(ii+1)-commd%bound_size(ii))*pard1

#ifdef MPI_OFF
#else
           bsize4=int(bsize,4)
           call MPI_Sendrecv( loc_sparr1(ini:), bsize4, &
                PAR_REAL,  dom_i, 0_4,      &
                loc_rparr1(ini:), bsize4,               &
                PAR_REAL, dom_i, 0_4,       &
                PAR_COMM_MY_CODE, status, istat )
#endif
           ini = ini + bsize
        end do

        ini = 0
        do ii = 1, nneig
           dom_i = commd%neights(ii)
           bsize = (commd%bound_size(ii+1)-commd%bound_size(ii))

           do jj = 1,bsize

              kk    = ini + jj
              idofj = (kk-1)*pard1
              ipoin = commd%bound_perm(kk)
              idofi = (ipoin-1)*pard1 

              if( pari1(ipoin) < 0 .and. pari1(ipoin) == -dom_i ) then
                 do idime = 1,pard1
                    idofj = idofj + 1
                    idofi = idofi + 1
                    parr1(idofi) = loc_rparr1(idofj)
                 end do
              end if

           end do

           ini = ini + bsize
        end do

        call memchk(two,istat,par_memor,'LOC_RPARR1','par_slexib',loc_rparr1)
        deallocate(loc_rparr1,stat=istat)
        if(istat/=0) call memerr(two,'LOC_RPARR1','par_slexib',0_ip)

        call memchk(two,istat,par_memor,'LOC_SPARR1','par_slexib',loc_sparr1)
        deallocate(loc_sparr1,stat=istat)
        if(istat/=0) call memerr(two,'LOC_SPARR1','par_slexib',0_ip)  

     else if( itask == 3 ) then

        !-------------------------------------------------------------
        !
        ! Identify fringe nodes
        !
        !-------------------------------------------------------------

        allocate(loc_spari1(commd%bound_dim),stat=istat)
        call memchk(zero,istat,par_memor,'LOC_SPARI1','par_slexib',loc_spari1)

        allocate(loc_rpari1(commd%bound_dim),stat=istat)
        call memchk(zero,istat,par_memor,'LOC_RPARR1','par_slexib',loc_rpari1)
        !
        ! Gather
        !
        do jj= 1, commd%bound_dim
           ipoin = commd%bound_perm(jj)
           loc_spari1(jj) = pari1(ipoin)
        end do

        do ii= 1, nneig
           dom_i = commd%neights(ii)

           ini   = commd%bound_size(ii)
           bsize = commd%bound_size(ii+1) - ini

#ifdef MPI_OFF
#else
           bsize4=int(bsize,4)
           call MPI_Sendrecv( loc_spari1(ini:), bsize4, &
                PAR_INTEGER,  dom_i, 0_4,               &
                loc_rpari1(ini:), bsize4,               &
                PAR_INTEGER, dom_i, 0_4,                &
                PAR_COMM_MY_CODE, status, istat )
#endif
        end do

        do jj= 1, commd%bound_dim
           ipoin = commd%bound_perm(jj)
           if( loc_rpari1(jj) < 0 ) pari1(ipoin) = loc_rpari1(jj)
        enddo

        call memchk(two,istat,par_memor,'LOC_RPARI1','par_slexib',loc_rpari1)
        deallocate(loc_rpari1,stat=istat)
        if(istat/=0) call memerr(two,'LOC_RPARI1','par_slexib',0_ip)

        call memchk(two,istat,par_memor,'LOC_SPARI1','par_slexib',loc_spari1)
        deallocate(loc_spari1,stat=istat)
        if(istat/=0) call memerr(two,'LOC_SPARI1','par_slexib',0_ip)

     else if( itask == 4 ) then

        nneno = 0_ip
        do ipoin = 1,npoin
           if (ipoin > npoi1) then
              nneno = nneno + pari1(ipoin)
           end if
        end do

        allocate(loc_spari1(nneno),stat=istat)
        call memchk(zero,istat,par_memor,'LOC_SPARI1','par_slexch',loc_spari1)
        
        allocate(loc_rpari1(nneno),stat=istat)
       call memchk(zero,istat,par_memor,'LOC_RPARI1','par_slexch',loc_rpari1)

       !
       ! Loop over neighbor subdomains
       !
       do ii = 1,nneig
          dom_i = commd%neights(ii)
          ini = 1
          bsize = 0
          kk    = 0
          !
          ! Loop over shared nodes with the current neighbor subdomain 
          !
          do jj = commd%bound_size(ii),commd%bound_size(ii+1) - 1
             ipoin = commd%bound_perm(jj)
             !
             ! Consider the fringe and travesti nodes
             !
             if (lntib(ipoin) < 0 .or. (lntib(ipoin) == 0 .and. lnti2(ipoin) > 0)) then                     
                do ll = 1,pari1(ipoin)
                   kk = kk + 1
                   !
                   ! List of fringe values sorted according the global numbering
                   !
                   loc_spari1(kk) = pari2(ll,ipoin)
                end do
                bsize = bsize + pari1(ipoin)
             end if
          end do
          !
          ! Parallel interchange
          !          
#ifdef MPI_OFF
#else
          bsize4 = int(bsize,4)
          call MPI_Sendrecv( loc_spari1(ini:), bsize4,&
               PAR_INTEGER,  dom_i, 0_4, &
               loc_rpari1(ini:), bsize4, &
               PAR_INTEGER, dom_i, 0_4, &
               PAR_COMM_MY_CODE, status, istat )
#endif
          ini = ini + bsize          
       end do


       do ii = 1,nneig
          kk    = 0
          do jj = commd%bound_size(ii),commd%bound_size(ii+1) - 1
             ipoin = commd%bound_perm(jj)
             !
             ! Consider the fringe and travesti nodes
             !
             if (lntib(ipoin) < 0 .or. (lntib(ipoin) == 0 .and. lnti2(ipoin) > 0)) then           
                do ll = 1,pari1(ipoin)
                   kk = kk + 1
                   !
                   ! List of fringe values sorted according the global numbering
                   !
                   if (loc_rpari1(kk) > -9e8 ) pari2(ll,ipoin) = loc_rpari1(kk)
                end do
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


  else if( itask == 5 ) then
     
     
  end if

  call cputim(time2)
  cpu_paral(25)=cpu_paral(25)+time2-time1

end subroutine par_slexib
