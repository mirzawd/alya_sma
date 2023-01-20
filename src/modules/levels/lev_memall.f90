!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine lev_memall(order)
  !-----------------------------------------------------------------------
  !****f* levels/lev_memall
  ! NAME 
  !    lev_memall
  ! DESCRIPTION
  !    This routine allocates memory for the arrays needed to solve the
  !    wave equation equation
  ! USES
  ! USED BY
  !    lev_turnon
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_inpout
  use def_master
  use def_domain
  use def_solver
  use def_levels 
  use mod_memory, only : memory_alloca
  use mod_memchk
  implicit none
  integer(ip), intent(in) :: order 
  integer(4)              :: istat

  select case(order)
     !
     !
     !
  case(1_ip)
     !
     ! Number of unknown components
     !
     if(kfl_tisch_lev==1) then
        !
        ! Trapezoidal rule
        !
        ncomp_lev=3

     else if(kfl_tisch_lev==2) then
        !
        ! BDF scheme
        !
        ncomp_lev=2+kfl_tiaor_lev
     end if

     if(kfl_timet_lev==1) then
        !
        ! Explicit treatment
        !
        ncomp_lev=4
     end if

     if( INOTMASTER ) then
        !
        ! level set equation unknown: amplitude LEVELS
        ! 
        call memory_alloca(mem_modul(1:2,modul),'FLEVE'    ,'lev_memall',fleve,npoin,ncomp_lev)
        call memory_alloca(mem_modul(1:2,modul),'DISTA'    ,'lev_memall',dista_lev,npoin)
        call memory_alloca(mem_modul(1:2,modul),'FLEV0_LEV','lev_memall',flev0_lev,npoin)
        call memory_alloca(mem_modul(1:2,modul),'NORML_LEV','lev_memall',norml_lev,ndime,npoin)
        call memory_alloca(mem_modul(1:2,modul),'NORMV_LEV','lev_memall',normv_lev,npoin)
        call memory_alloca(mem_modul(1:2,modul),'ICUOT_LEV','lev_memall',icupt_lev,npoin)
        call memory_alloca(mem_modul(1:2,modul),'ELCRO_LEV','lev_memall',elcro_lev,nelem)
        call memory_alloca(mem_modul(1:2,modul),'FLSEX_LEV','lev_memall',flsex_lev,npoin)
        !
        ! Used in reinitialization
        !
        if(( tyred_lev == 3 ).or.( tyred_lev == 5 ).or.( tyred_lev == 6 )) then
           call memory_alloca(mem_modul(1:2,modul),'PSBEL_LEV','lev_memall',psbel_lev,nelwh+1)  
        end if

        !
        ! Allocate ! TESTEO
        ! 
        !if( .not. associated(dispm) ) call memory_alloca(mem_modul(1:2,modul),'DISPM','lev_memall',dispm,ndime,npoin,3_ip)
        !if( .not. associated(velom) ) call memory_alloca(mem_modul(1:2,modul),'VELOM','lev_memall',velom,ndime,npoin)
        !
        ! Solver
        !
        if( tyred_lev == 1 ) then
           solve_sol => solve(1:1)
           call soldef(4_ip)
        else
          solve_sol => solve(1:)
          call soldef(4_ip)
        end if
        !solve_sol => solve(4:4)  ! TESTEO
        !call soldef(4_ip)        ! TESTEO

     else
        allocate(fleve(1,3),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'LEVEL','lev_memall',fleve) !!!! PORQUE DECIA WAVAM ACA !?!?!?
        allocate(norml_lev(1,1),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),' GRADPHI ','lev_memall',norml_lev)
     end if

  case(2_ip)

     if( ISEQUEN ) then

        kfl_alloc_lev = 1

        allocate(lnodb_lev(ndime,nboun_lev),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'LNODB_LEV','lev_memall',lnodb_lev)

        allocate(coord_lev(ndime,npoin_lev),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'COORD_LEV','lev_memall',coord_lev)

        allocate(lebsu_lev(nboun_lev),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'LEBSU_LEV','lev_memall',lebsu_lev)

     else if( ISLAVE ) then

        allocate(lnodp_lev(ndime,nredi_lev(1)),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'LNODP_LEV','lev_memall',lnodp_lev)

        allocate(coorp_lev(ndime,nredi_lev(2)),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'COORP_LEV','lev_memall',coorp_lev)

        allocate(lebsp_lev(nredi_lev(1)),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'LEBSP_LEV','lev_memall',lebsp_lev)

     end if

  case(3_ip)
     
     kfl_alloc_lev = 1

     allocate(lnodb_lev(ndime,nboun_lev),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'LNODB_LEV','lev_memall',lnodb_lev)

     allocate(coord_lev(ndime,npoin_lev),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'COORD_LEV','lev_memall',coord_lev)

     allocate(lebsu_lev(nboun_lev),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'LEBSU_LEV','lev_memall',lebsu_lev)

  case(5_ip)

     if( IMASTER ) then
        allocate(nredm_lev(2,pard1),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'NREDM_LEV','lev_memall',nredm_lev)     
     end if

  case(7_ip) 
     !
     ! Deallocate what was allocated in case 5
     !
     if ( IMASTER ) then
        
        call memchk(two,istat,mem_modul(1:2,modul),'NREDM_LEV','lev_memall',nredm_lev)
        deallocate(nredm_lev,stat=istat)
        if(istat/=0) call memerr(two,'NREDM_LEV','lev_memall',0_ip)

     else if( ISLAVE ) then

        call memchk(two,istat,mem_modul(1:2,modul),'LNODP_LEV','lev_memall',lnodp_lev)
        deallocate(lnodp_lev,stat=istat)
        if(istat/=0) call memerr(two,'LNODP_LEV','lev_memall',0_ip)

        call memchk(two,istat,mem_modul(1:2,modul),'COORP_LEV','lev_memall',coorp_lev)
        deallocate(coorp_lev,stat=istat)
        if(istat/=0) call memerr(two,'COORP_LEV','lev_memall',0_ip)

        call memchk(two,istat,mem_modul(1:2,modul),'LEBSP_LEV','lev_memall',lebsp_lev)
        deallocate(lebsp_lev,stat=istat)
        if(istat/=0) call memerr(two,'LEBSP_LEV','lev_memall',0_ip)

     end if

  case(8_ip) 
     !
     ! Deallocate
     !
     if( ISLAVE ) then

        call memchk(two,istat,mem_modul(1:2,modul),'LNODB_LEV','lev_memall',lnodb_lev)
        deallocate(lnodb_lev,stat=istat)
        if(istat/=0) call memerr(two,'LNODB_LEV','lev_memall',0_ip)
        
        call memchk(two,istat,mem_modul(1:2,modul),'COORD_LEV','lev_memall',coord_lev)
        deallocate(coord_lev,stat=istat)
        if(istat/=0) call memerr(two,'COORD_LEV','lev_memall',0_ip)

        call memchk(two,istat,mem_modul(1:2,modul),'LEBSU_LEV','lev_memall',lebsu_lev)
        deallocate(lebsu_lev,stat=istat)
        if(istat/=0) call memerr(two,'LEBSU_LEV','lev_memall',0_ip)

     end if

  case(9_ip)
     !
     ! Allocate bvcod_lev(:) analogously to bvcod in kernel/domain/membcs.f90
     !   
     call memory_alloca(mem_modul(1:2,modul),'BVCOD_LEV','lev_memall',bvcod_lev,npoin)
     
  case(10_ip) 
     !
     ! Deallocate 
     !
     if( kfl_alloc_lev == 1 ) then

        kfl_alloc_lev = 0

        if( INOTSLAVE ) then

           call memchk(two,istat,mem_modul(1:2,modul),'LNODB_LEV','lev_memall',lnodb_lev)
           deallocate(lnodb_lev,stat=istat)
           if(istat/=0) call memerr(two,'LNODB_LEV','lev_memall',0_ip)
           
           call memchk(two,istat,mem_modul(1:2,modul),'COORD_LEV','lev_memall',coord_lev)
           deallocate(coord_lev,stat=istat)
           if(istat/=0) call memerr(two,'COORD_LEV','lev_memall',0_ip)

           call memchk(two,istat,mem_modul(1:2,modul),'LEBSU_LEV','lev_memall',lebsu_lev)
           deallocate(lebsu_lev,stat=istat)
           if(istat/=0) call memerr(two,'LEBSU_LEV','lev_memall',0_ip)

        end if

     end if

  end select

end subroutine lev_memall
      
