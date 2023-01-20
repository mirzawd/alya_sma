!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine par_cregro()
  !-----------------------------------------------------------------------
  !****f* Parall/par_cregro
  ! NAME
  !    par_cregro
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
  use mod_parall, only : par_memor
  implicit none
  integer(ip)          :: ipoin,jpoin,domai,indice1,ngrou,nzgro
  integer(4)           :: istat
  integer(ip), target  :: dummi(2) 
  integer(ip), pointer :: lgrou_loc(:)

  if( IPARALL ) then

     !-------------------------------------------------------------------
     !
     ! Broadcast dimensions
     !
     !-------------------------------------------------------------------

     call par_livinf(14_ip,' ',modul)
     !
     ! Broadcast NGROU, NZGRO
     !
     npari    =  2
     nparr    =  0
     dummi(1) =  solve_sol(1)%ngrou
     dummi(2) =  solve_sol(1)%nskyl
     parin    => dummi
     strin    =  'NGROU,NSKYL'
     call par_broadc()
     solve_sol(1)%ngrou = parin(1)
     solve_sol(1)%nskyl = parin(2)
     ngrou              = solve_sol(1)%ngrou
     nzgro              = solve_sol(1)%nskyl
     !
     ! Slaves allocate memory for ISKYL (master does it only for broacdast)
     !
     if( ISLAVE .or. ( IMASTER .and. READ_AND_RUN() ) ) then
        allocate(solve_sol(1)%iskyl(ngrou+1),stat=istat)
        call memchk(0_ip,istat,memit,'ISKYL','par_cregro',solve_sol(1)%iskyl)
     end if
     !
     ! Broadcast ISKYL(NGROU)
     !
     npari =  ngrou+1
     parin => solve_sol(1)%iskyl
     strin =  'ISKYL'
     call par_broadc()
     !
     ! Master deallocate ISKYL
     !
     if( IMASTER ) then
        call memchk( two, istat, par_memor, 'ISKYL','par_sengeo', solve_sol(1)%iskyl)
        deallocate( solve_sol(1)%iskyl, stat=istat )
        if(istat/=0) call memerr( two, 'ISKYL', 'par_sengeo',0_ip)
     end if

  end if

  if( IMASTER .and. .not. READ_AND_RUN() ) then

     !-------------------------------------------------------------------
     !
     ! Master: construct local-to-global group
     !
     !-------------------------------------------------------------------

     allocate( lgrou_loc(npoin_total),stat=istat)
     call memchk( zero, istat, par_memor, 'lgrou_loc', 'par_sengeo', lgrou_loc )
     do ipoin = 1, npoin_total
        jpoin            = lninv_loc(ipoin)
        lgrou_loc(ipoin) = solve_sol(1)%lgrou(jpoin)
     end do
     !
     ! Send LGROU_LOC
     !
     nparc   = 0
     indice1 = 1
     do domai= 1, npart_par
        kfl_desti_par = domai
        npari =  npoin_par(domai)
        parin => lgrou_loc(indice1:)
        strin =  'LGROU_LOC'
        call par_sendin()
        indice1 = indice1 + npoin_par(domai)
     end do
     !
     ! Deallocate memory
     !
     call memchk( two, istat, par_memor, 'LGROU_LOC','par_sengeo', lgrou_loc)
     deallocate( lgrou_loc, stat=istat ) 
     if(istat/=0) call memerr( two, 'LGROU_LOC', 'par_sengeo',0_ip)
     call memchk( two, istat, par_memor, 'LGROU','par_sengeo', solve_sol(1)%lgrou)
     deallocate( solve_sol(1)%lgrou, stat=istat )
     if(istat/=0) call memerr( two, 'LGROU', 'par_sengeo',0_ip)


  else if( ISLAVE ) then

     !-------------------------------------------------------------------
     !
     ! Slaves: allocate memory and receive groups (LGROU)
     !
     !-------------------------------------------------------------------

     allocate(solve_sol(1)%lgrou(npoin),stat=istat)
     call memchk(zero,istat,memit,'LGROU','cregro',solve_sol(1)%lgrou)
     !
     ! Receive LGROU
     !
     nparr = 0
     nparc = 0
     kfl_desti_par = 0
     npari =  npoin
     parin => solve_sol(1)%lgrou
     call par_receiv()

  end if

end subroutine par_cregro
