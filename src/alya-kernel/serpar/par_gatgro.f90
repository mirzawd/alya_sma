!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine par_gatgro()
  !-------------------------------------------------------------------------------
  !****f* parall/par_grogro
  ! NAME
  !    par_grogro
  ! DESCRIPTION
  !    Set up the data base COML(ICOML) for the allgatherv used in the 
  !    deflated CG. Simple Example:
  !    a. 2 subdomains
  !    b. Subdomain 1 has groups 1,2,3.
  !    c. Subdomain 2 has groups 1 and 3
  !    DISPL  = [ 0 0 3 ]
  !    LCOUN  = [ 0 3 2 ]
  !    LBIG   = [ 2 1 3 | 1 3 ]
  !    NSMALL = 3 or 2
  ! INPUT
  ! OUTPUT
  ! USED BY
  !***
  !-------------------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_solver
  use def_parall
  use mod_memchk
  use mod_parall, only : par_memor
  use mod_memory, only : memory_copy
  use mod_communications_global, only : PAR_SUM
  implicit none
  integer(ip)          :: igrou, jgrou, kk, ipoin, ipart
  integer(ip)          :: ngrou  
  integer(4)           :: istat
  integer(ip), pointer :: ngloc(:),domli(:,:),lgrou(:)

  ngrou =  solve_sol(1) % ngrou
  lgrou => solve_sol(1) % lgrou

  allocate(ngloc(npart_par),stat=istat)
  call memchk(zero,istat,par_memor,'ngloc','par_domgro',ngloc)
  allocate(domli(ngrou,npart_par),stat=istat)
  call memchk(zero,istat,par_memor,'domli','par_memgro',domli)

  call par_memory(11_ip) ! DISPL and LCOUN

  if( INOTMASTER ) then

     !-------------------------------------------------------------------
     !
     ! NDOMI: number of subdomain to which igrou belongs
     ! DOMLI: list of subdomains to which igrou belongs
     !
     !-------------------------------------------------------------------

     call memgen(1_ip,ngrou,0_ip)
     do ipoin = 1,npoin 
        igrou   = lgrou(ipoin)
        if( igrou > 0 ) then
           gisca(igrou) = 1 
        end if
     end do
     do igrou = 1,ngrou
        if( gisca(igrou) == 1 ) then
           ngloc(kfl_paral) = ngloc(kfl_paral) + 1
           domli(ngloc(kfl_paral),kfl_paral) = igrou
        end if
     end do
     call memgen(3_ip,ngrou,0_ip)

  end if

  call PAR_SUM(ngloc)
  call PAR_SUM(domli)
  !
  ! LCOUN(IPART): Size of array subdomain IPART has to exchange (0 is master)
  ! DISPL(IPART): Displacement of array in XBIG for subdomain IPART
  ! LBIG(IPART):  Permutation array for subdomain IPART
  !
  solve_sol(1) % displ(1) = 0
  solve_sol(1) % displ(2) = 0
  solve_sol(1) % lcoun(1) = 0
  solve_sol(1) % nbig     = 0
  
  do ipart = 2,npart_par
     solve_sol(1) % displ(ipart+1) = solve_sol(1) % displ(ipart) + ngloc(ipart-1)
  end do
  
  do ipart = 1,npart_par
     solve_sol(1) % lcoun(ipart+1) = ngloc(ipart)
     solve_sol(1) % nbig           = solve_sol(1) % nbig + ngloc(ipart)
  end do
  solve_sol(1) % nsmall = solve_sol(1) % lcoun(kfl_paral+1)

  call par_memory(12_ip) ! LBIG and XBIG

  kk = 0
  do ipart = 1,npart_par
     do igrou = 1,ngloc(ipart)
        jgrou = domli(igrou,ipart)
        kk    = kk + 1
        solve_sol(1) % lbig(kk) = jgrou
     end do
  end do
  solve_sol(1) % nsmall = 0

  call memchk(two,istat,par_memor,'ngloc','par_gatgro',ngloc)
  deallocate(ngloc,stat=istat)
  if(istat/=0) call memerr(two,'ngloc','par_gatgro',0_ip) 

  call memchk(two,istat,par_memor,'domli','par_gatgro',domli)
  deallocate(domli,stat=istat)
  if(istat/=0) call memerr(two,'domli','par_gatgro',0_ip) 
  !
  ! Master does not need LBIG
  !
  if( IMASTER ) call par_memory(-15_ip)

  !----------------------------------------------------------------------
  !
  ! DISPL and LCOUN must be Integer(4) for MPI_ALLGATHERV
  !
  !----------------------------------------------------------------------

#ifdef I8
  call par_memory( 13_ip)
  do kk = 1,npart_par + 1
     solve_sol(1) % disp4(kk) = int(solve_sol(1) % displ(kk),4_4)
     solve_sol(1) % lcou4(kk) = int(solve_sol(1) % lcoun(kk),4_4)
  end do
  call par_memory(-11_ip)
#else
  call memory_copy(par_memor,'SOLVE % DISP4','par_gatgro',solve_sol(1) % displ,solve_sol(1) % disp4,'DO_NOT_DEALLOCATE')
  call memory_copy(par_memor,'SOLVE % LCOU4','par_gatgro',solve_sol(1) % lcoun,solve_sol(1) % lcou4,'DO_NOT_DEALLOCATE')
#endif

  !icoml = icoml + 1

end subroutine par_gatgro
