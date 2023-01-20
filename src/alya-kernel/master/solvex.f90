!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine solvex(rhsix,unknx,amatx,pmatx)
  !-----------------------------------------------------------------------
  !****f* master/solvex
  ! NAME 
  !    solvex
  ! DESCRIPTION
  !    This routine calls the solvexs
  !    For diagonal solve which uses vmass, amatx must NOT be modified
  !
  !    About residual RHSIX:
  !
  !    solve_sol(1)%kfl_recov = 0 ... Do not recover local residual
  !                           = 1 ... Recover local residual
  !                           = 2 ... Residual is already global
  !
  !    About solve tolerance SOLCO:
  !
  !    solve_sol(1)%kfl_adres = 0 ... Solvex tolerance is given by user
  !                           = 1 ... Solvex tolerance is adaptive
  !
  ! USES
  !    memchk
  !    mediso
  ! USED BY
  !***
  !-----------------------------------------------------------------------
  use def_kintyp
  use def_master, only       :  NPOIN_REAL_12DI,INOTSLAVE,NPOIN_TYPE
  use def_domain, only       :  npoin,c_dom,r_dom
  use def_solver, only       :  solve_sol,cpu_solve,resi1,resi2,&
       &                        resin,resfi,iters
  use mod_outfor, only       :  outfor
  use mod_memchk
  implicit none
  complex(rp), intent(inout) :: unknx(*)
  complex(rp), intent(in)    :: amatx(*)
  complex(rp), intent(in)    :: pmatx(*)
  complex(rp), intent(inout) :: rhsix(*)
  real(rp)                   :: time1,time2
  real(rp)                   :: cpre1,cpre2

  call cputim(time1)
  !
  ! Headers
  !
  if(solve_sol(1)%kfl_algso/=9) then
     if(solve_sol(1)%heade==0) then
        solve_sol(1)%heade=1
        call outfor(39_ip,0_ip,' ')
     end if
     call outfor(5_ip,0_ip,' ')
  end if
  !
  ! Modify RHS due to periodicity and Parall service: RHSIX
  !
  if( solve_sol(1)%kfl_recov /= 2 ) call pararx('SLX',NPOIN_TYPE,npoin*solve_sol(1)%ndofn,rhsix)
  !
  ! Initial algebraic residual: RESIN = ||b-Ax||/||b||
  !
  !!call algres(rhsix,unknx,amatx)
  !
  ! Adaptive residual tolerance: SOLCO = max ( alpha*r_0 , eps_min )
  !
  !!if( solve_sol(1)%kfl_adres == 1 ) then
  !!   solve_sol(1)%solco = max( solve_sol(1)%resin*solve_sol(1)%adres , solve_sol(1)%solmi )
  !!end if
  !
  ! Algebraic solver
  !

  call cputim(cpre1)

  if( solve_sol(1)%kfl_algso == 5 ) then

     call bcgplx(&
          npoin,solve_sol(1)%ndofn,&
          solve_sol(1)%miter,solve_sol(1)%kfl_preco,&
          solve_sol(1)%solco,amatx,c_dom,r_dom,&
          rhsix,unknx)

  elseif( solve_sol(1)%kfl_algso == 17 ) then

     call bcgptx(&
          npoin,solve_sol(1)%ndofn,&
          solve_sol(1)%miter,solve_sol(1)%kfl_preco,&
          solve_sol(1)%solco,amatx,c_dom,r_dom,&
          rhsix,unknx)


  else

     call runend('SOLVEX: COMPLEX SOLVER NOT CODED')

  end if

  call cputim(cpre2)
  call cputim(time2)
  solve_sol(1)%cputi = solve_sol(1)%cputi + (time2-time1)
  cpu_solve          = cpu_solve          + (time2-time1)
  solve_sol(1)%nsolv = solve_sol(1)%nsolv + 1

  if( INOTSLAVE ) then
     if( solve_sol(1)%kfl_cvgso == 1 ) then
        write(solve_sol(1)%lun_cvgso,100) iters,resfi
        flush(solve_sol(1)%lun_cvgso)
     end if
     if( solve_sol(1)%kfl_solve == 1 ) then
        if( resi2>0.0_rp ) then
           if( resi1/=0.0_rp ) then
              write(solve_sol(1)%lun_solve,110)&
                   iters,resin,resfi,log10(resi2/resi1),cpre2-cpre1
           else
              write(solve_sol(1)%lun_solve,110)&
                   iters,resin,resfi,0.0_rp,cpre2-cpre1    
           end if
        else
           write(solve_sol(1)%lun_solve,110)&
                iters,resin,resfi,1.0_rp,cpre2-cpre1    
        end if
        flush(solve_sol(1)%lun_solve)
     end if
  end if
 
100 format(i7,1x,e12.6)
110 format(i5,18(2x,e12.6))

end subroutine solvex
