!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine solite(rhsid,unkno,amatr,pmatr)
  !------------------------------------------------------------------------
  !****f* solite/solite
  ! NAME
  !    solite
  ! DESCRIPTION
  !    This routine drives the library to solve a linear system of equations
  !    by iterative methods.
  ! USES
  ! USED BY
  !***
  !------------------------------------------------------------------------
  use def_master,                        only : INOTSLAVE,IMASTER,ISLAVE,NPOIN_TYPE
  use def_master,                        only : IEMPTY
  use def_domain,                        only : r_sym,c_sym,npoin
  use def_solver
  use mod_iterative_solver
  use mod_alya2maphys,                   only : alya2maphys_solution
  use mod_alya2agmg,                     only : alya2agmg_solution
  use mod_iofile,                        only : iofile_flush_unit
  use mod_communications_point_to_point, only : PAR_INTERFACE_NODE_EXCHANGE
#ifdef INC_PSBLAS
  use mod_alya2psblas,                   only : alya2psblas_solution
#endif

  implicit none
  real(rp),    target  :: rhsid(*)
  real(rp),    target  :: unkno(*)
  real(rp),    target  :: amatr(*)
  real(rp)             :: pmatr(*)
  real(rp)             :: cpre1,cpre2,logres
  integer(ip), pointer :: ia(:),ja(:)
  integer(ip), target  :: pointer_tmp(2)

  call cputim(cpre1)
  
  !----------------------------------------------------------------------
  !
  ! Force solution continuity across subdomain boundaries
  !
  !----------------------------------------------------------------------
 
  if( solve_sol(1) % kfl_force == 1 .and. ISLAVE .and. npoin > 0 ) then
     call PAR_INTERFACE_NODE_EXCHANGE(solve_sol(1) % ndofn,unkno,'TAKE MIN','IN MY CODE')
  end if

  !----------------------------------------------------------------------
  !
  ! Check if full row format should used and select
  ! the appropriate graph
  !
  !----------------------------------------------------------------------

  if( solve_sol(1) % kfl_full_rows == 1 ) then
     ia => solve_sol(1) % ia_full
     ja => solve_sol(1) % ja_full
  else
     ia => solve_sol(1) % ia
     ja => solve_sol(1) % ja
  end if
  if( IEMPTY ) then
     if( .not. associated(ia) ) ia => pointer_tmp
     if( .not. associated(ja) ) ja => pointer_tmp
  end if

  !----------------------------------------------------------------------
  !
  ! Iterative solvers
  !
  !----------------------------------------------------------------------

  if( solve_sol(1) % kfl_schur == 1 ) then
     !
     ! Schur complement solver
     !
     if( IMASTER ) then
        call solshu(1_ip,1_ip,1_ip,rhsid,unkno,amatr)
     else
        call solshu(solve_sol(1) % nzmat,solve_sol(1) % ndofn,npoin,rhsid,unkno,amatr)
     end if

  else if( solve_sol(1) % kfl_algso == SOL_SOLVER_MINRES ) then
     !
     ! MINRES
     !
     call minres( &
          solve_sol(1) % nequa,     solve_sol(1) % ndofn, &
          solve_sol(1) % kfl_preco, solve_sol(1) % miter, &
          solve_sol(1) % solco, &
          amatr, pmatr, solve_sol(1) % kfl_cvgso, &
          solve_sol(1) % lun_cvgso, solve_sol(1) % kfl_solve, &
          solve_sol(1) % lun_solve, ja, ia, rhsid, unkno )

  else if( solve_sol(1) % kfl_algso == SOL_SOLVER_CG ) then
     !
     ! CG from PLS
     !
     call cgrpls( &
          solve_sol(1) % nequa,     solve_sol(1) % ndofn, &
          solve_sol(1) % kfl_preco, solve_sol(1) % miter, &
          solve_sol(1) % solco, amatr, pmatr, &
          solve_sol(1) % kfl_cvgso, solve_sol(1) % lun_cvgso,&
          solve_sol(1) % kfl_solve, solve_sol(1) % lun_solve,&
          ja, ia, rhsid, unkno )


  else if( solve_sol(1) % kfl_algso == SOL_SOLVER_DEFLATED_CG .or. &
       &   solve_sol(1) % kfl_algso == SOL_SOLVER_A_DEF2      ) then
     !
     ! Deflated CG
     !
     if( solve_sol(1) % kfl_symme == 1 ) then
        call deflcg( &
             solve_sol(1) % nequa,     solve_sol(1) % ndofn, &
             solve_sol(1) % kfl_preco, solve_sol(1) % miter, &
             solve_sol(1) % solco, &
             amatr, pmatr, solve_sol(1) % kfl_cvgso, &
             solve_sol(1) % lun_cvgso, solve_sol(1) % kfl_solve, &
             solve_sol(1) % lun_solve, c_sym, r_sym, rhsid, unkno )
     else
        call deflcg( &
             solve_sol(1) % nequa,     solve_sol(1) % ndofn, &
             solve_sol(1) % kfl_preco, solve_sol(1) % miter, &
             solve_sol(1) % solco, &
             amatr, pmatr, solve_sol(1) % kfl_cvgso, &
             solve_sol(1) % lun_cvgso, solve_sol(1) % kfl_solve, &
             solve_sol(1) % lun_solve, ja, ia, rhsid, unkno )
     end if
#ifdef NINJA
  else if( solve_sol(1) % kfl_algso == SOL_GDECG ) then
     !
     ! GPU Deflated CG
     !
     call GPUDEFLCG( &
          solve_sol(1) % nequa, solve_sol(1) % ndofn, &
          solve_sol(1) % ngrou, solve_sol(1) % miter, &
          amatr, ia, ja, rhsid, unkno)
#endif
  else if( solve_sol(1) % kfl_algso == SOL_SOLVER_BICGSTAB ) then
     !
     ! Bi-CGSTAB from PLS
     !
     call bcgpls( &
          solve_sol(1) % nequa,     solve_sol(1) % ndofn, &
          solve_sol(1) % kfl_preco, solve_sol(1) % miter, &
          solve_sol(1) % solco,     amatr, pmatr, &
          solve_sol(1) % kfl_cvgso, solve_sol(1) % lun_cvgso,&
          solve_sol(1) % kfl_solve, solve_sol(1) % lun_solve,&
          ja, ia, rhsid, unkno )

  else if( solve_sol(1) % kfl_algso == SOL_SOLVER_GMRES ) then
     !
     ! GMRES from PLS
     !
     call gmrpls( &
          solve_sol(1) % nequa,     solve_sol(1) % ndofn,&
          solve_sol(1) % kfl_preco, solve_sol(1) % miter,&
          solve_sol(1) % nkryd,     solve_sol(1) % kfl_ortho, &
          solve_sol(1) % kfl_cvgso, solve_sol(1) % lun_cvgso, &
          solve_sol(1) % kfl_solve, solve_sol(1) % lun_solve, &
          solve_sol(1) % solco, &
          amatr, pmatr, ja , ia , rhsid, unkno )
#ifdef NINJA
  else if( solve_sol(1) % kfl_algso == SOL_GGMRS ) then
     !
     ! GMRES GPU
     !
     call GPUGMRES( &
          solve_sol(1) % nequa,     solve_sol(1) % ndofn,&
          solve_sol(1) % miter,&
          solve_sol(1) % nkryd,&
          amatr, ia, ja,rhsid, unkno )
#endif

#ifdef AMGX
  else if( solve_sol(1) % kfl_algso == SOL_GAMGX ) then
     !
     ! AMGX GPU
     !
     if( ISEQUEN ) then

        call GPUAMGX( &
             solve_sol(1) % nequa,     solve_sol(1) % ndofn,&
             solve_sol(1) % miter,&
             amatr, ia, ja,rhsid, unkno )

     else

        !call GPUAMGXPAR( &
        !     solve_sol(1) % nequa,     solve_sol(1) % ndofn,&
        !     solve_sol(1) % miter,&
        !     amatr, ia, ja,rhsid, unkno )

        call alya2wsmp_solution(solve_sol(1),amatr,rhsid,unkno)

     end if


#endif

  else if( solve_sol(1) % kfl_algso == SOL_SOLVER_RICHARDSON .or. &
       &   solve_sol(1) % kfl_algso == SOL_SOLVER_MATRIX_RICHARDSON ) then
     !
     ! Richardson solver
     !
     call solric(solve_sol(1) % nequa,solve_sol(1) % ndofn,rhsid,unkno,amatr,pmatr)

  else if( solve_sol(1) % kfl_algso == 11 ) then
     !
     ! Richardson
     !
     call runend('SOLITE: OBSOLETE')
     if( solve_sol(1) % kfl_symme == 1 ) then
        call richar( &
             solve_sol(1) % nequa,     solve_sol(1) % ndofn, &
             solve_sol(1) % kfl_preco, solve_sol(1) % miter, &
             solve_sol(1) % solco,     amatr, pmatr, &
             solve_sol(1) % kfl_cvgso, solve_sol(1) % lun_cvgso,&
             solve_sol(1) % kfl_solve, solve_sol(1) % lun_solve,&
             c_sym, r_sym, rhsid, unkno )
     else
        call richar( &
             solve_sol(1) % nequa,     solve_sol(1) % ndofn, &
             solve_sol(1) % kfl_preco, solve_sol(1) % miter, &
             solve_sol(1) % solco,     amatr, pmatr, &
             solve_sol(1) % kfl_cvgso, solve_sol(1) % lun_cvgso,&
             solve_sol(1) % kfl_solve, solve_sol(1) % lun_solve,&
             ja, ia, rhsid, unkno )
     end if

  else if( solve_sol(1) % kfl_algso == SOL_SOLVER_DEFLATED_BICGSTAB ) then
     !
     ! Deflated Bi-CGSTAB
     !
     call defbcg( &
          solve_sol(1) % nequa,     solve_sol(1) % ndofn, &
          solve_sol(1) % kfl_preco, solve_sol(1) % miter, &
          solve_sol(1) % solco,     amatr, pmatr, &
          solve_sol(1) % kfl_cvgso, solve_sol(1) % lun_cvgso,&
          solve_sol(1) % kfl_solve, solve_sol(1) % lun_solve,&
          ja, ia, rhsid, unkno )

  else if( solve_sol(1) % kfl_algso == SOL_SOLVER_DEFLATED_GMRES ) then
     !
     ! Deflated GMRES
     !
     call defgmr( &
          solve_sol(1) % nequa,     solve_sol(1) % ndofn,&
          solve_sol(1) % kfl_preco, solve_sol(1) % miter, solve_sol(1) % nkryd,&
          solve_sol(1) % kfl_cvgso, solve_sol(1) % lun_cvgso,&
          solve_sol(1) % kfl_solve, solve_sol(1) % lun_solve,&
          solve_sol(1) % solco, amatr, pmatr, ja, ia, rhsid, unkno )

  else if( solve_sol(1) % kfl_algso == SOL_SOLVER_STEEPEST_DESCENT ) then
     !
     ! Steepest descent
     !
     if( solve_sol(1) % kfl_symme == 1 ) then
        call steepe( &
             solve_sol(1) % nequa,     solve_sol(1) % ndofn, &
             solve_sol(1) % kfl_preco, solve_sol(1) % miter, &
             solve_sol(1) % solco, amatr, pmatr, &
             solve_sol(1) % kfl_cvgso, solve_sol(1) % lun_cvgso,&
             solve_sol(1) % kfl_solve, solve_sol(1) % lun_solve,&
             c_sym, r_sym, rhsid, unkno )
     else
        call steepe( &
             solve_sol(1) % nequa,     solve_sol(1) % ndofn, &
             solve_sol(1) % kfl_preco, solve_sol(1) % miter, &
             solve_sol(1) % solco, amatr, pmatr, &
             solve_sol(1) % kfl_cvgso, solve_sol(1) % lun_cvgso,&
             solve_sol(1) % kfl_solve, solve_sol(1) % lun_solve,&
             ja, ia, rhsid, unkno )
     end if
#ifdef NINJA
  else if(solve_sol(1) % kfl_algso == SOL_GCG ) then
     call GPUCG( &
          solve_sol(1) % nequa,solve_sol(1) % ndofn, &
          solve_sol(1) % miter,&
          amatr,ia, ja, rhsid, unkno )
  else if(solve_sol(1) % kfl_algso == SOL_GCGNOPREC ) then
     call GPUCGNOPREC( &
          solve_sol(1) % nequa,solve_sol(1) % ndofn, &
          solve_sol(1) % miter,&
          amatr,ia, ja, rhsid, unkno )


  else if(solve_sol(1) % kfl_algso == SOL_GPCG ) then
     call GPUPIPECG( &
          solve_sol(1) % nequa,solve_sol(1) % ndofn, &
          solve_sol(1) % miter,&
          amatr,ia, ja, rhsid, unkno )

#endif
  else if( solve_sol(1) % kfl_algso == SOL_SOLVER_PIPELINED_CG         .or. &
       &   solve_sol(1) % kfl_algso == SOL_SOLVER_PIPELINED_DEFLATED_CG ) then
     !
     ! Pipelined CG
     !
     call pipelined_CG_rr( &
          solve_sol(1) % nequa,     solve_sol(1) % ndofn, &
          solve_sol(1) % kfl_preco, solve_sol(1) % miter, &
          solve_sol(1) % solco,     amatr, pmatr, &
          solve_sol(1) % kfl_cvgso, solve_sol(1) % lun_cvgso,&
          solve_sol(1) % kfl_solve, solve_sol(1) % lun_solve,&
          ja, ia, rhsid, unkno )

  else if( solve_sol(1) % kfl_algso == SOL_SOLVER_MAPHYS_UNSYMMETRIC .or. &
       &   solve_sol(1) % kfl_algso == SOL_SOLVER_MAPHYS_SYMMETRIC   .or. &
       &   solve_sol(1) % kfl_algso == SOL_SOLVER_MAPHYS_SPD         ) then
     !
     ! MAPHYS solver
     !
     call alya2maphys_solution(solve_sol(1),amatr,rhsid,unkno)

  else if( solve_sol(1) % kfl_algso == SOL_SOLVER_POLYNOMIAL ) then
     !
     ! Polynomial solver
     !
     call iterative_solver_polynomial(solve_sol,amatr,rhsid,unkno)
     
  else if( solve_sol(1) % kfl_algso == SOL_SOLVER_AGMG ) then
     !
     ! AGMG solver
     !
#ifdef PARAL_AGMG
     call alya2agmg_solution(solve_sol(1),amatr,rhsid,unkno)
#endif

  else if( solve_sol(1) % kfl_algso == SOL_SOLVER_PSBLAS ) then
     !
     ! PSBLAS solver
     !
#ifdef INC_PSBLAS
     call alya2psblas_solution(solve_sol(1),amatr,rhsid,unkno)
#endif

else if( solve_sol(1) % kfl_algso == SOL_SOLVER_PETSC ) then
   !
   ! PSETSc solver
   !
#ifdef PETSC
   block
   use mod_alya2petsc
   call alya2petsc_solution(solve_sol(1) % petsc, amatr, rhsid, unkno)
   call alya2petsc_interfaceNodeExchange(solve_sol(1) % petsc, rhsid, unkno)
   call alya2petsc_getOutInfo(solve_sol(1) % petsc, solve_sol(1) % iters)
   call alya2petsc_solverInfoExchange(solve_sol(1))
   endblock
#endif

  end if

  !----------------------------------------------------------------------
  !
  ! Force solution continuity across subdomain boundaries
  !
  !----------------------------------------------------------------------

  if( solve_sol(1) % kfl_force == 1 .and. ISLAVE ) then
     call PAR_INTERFACE_NODE_EXCHANGE(solve_sol(1) % ndofn,unkno,'TAKE MIN','IN MY CODE')
  end if

  !----------------------------------------------------------------------
  !
  ! Output solver information
  !
  !----------------------------------------------------------------------

  call cputim(cpre2)
  if( solve_sol(1) % kfl_cvgso == 1 .and. INOTSLAVE ) call iofile_flush_unit(solve_sol(1) % lun_cvgso)
  if( INOTSLAVE ) then
     if( solve_sol(1) % kfl_solve == 1 ) then
        if(solve_sol(1) % kfl_algso == SOL_SOLVER_AGMG   )  resi2 = -10.0_rp  ! I set it to negative so that it goes to the third print
        if(solve_sol(1) % kfl_algso == SOL_SOLVER_PSBLAS )  resi2 = -10.0_rp  ! missing check if this is corect - otherwise resi2 was NaN
        if(solve_sol(1) % kfl_algso == SOL_SOLVER_PETSC )   resi2 = -10.0_rp  ! missing check if this is corect - otherwise resi2 was NaN
        ! I need to clarify with guillaume why there are solve_sol(1) % resi2  & resi2   -- isn't it a repetition - also clarify resi1
       if( resi2 > 0.0_rp ) then
           if( resi1 /= 0.0_rp ) then
              logres = log10(resi2/resi1)
           else
              logres = 0.0_rp
           end if
        else
           logres = 1.0_rp
        end if
        write(solve_sol(1) % lun_solve,110)&
             !                  1                    2                    3
             solve_sol(1) % iters,solve_sol(1) % resin,solve_sol(1) % resfi,&
             !                  4                    5        6
             solve_sol(1) % resi2,solve_sol(1) % resf2,logres,&
             !         7                    8                    9
             cpre2-cpre1,solve_sol(1) % xorth,solve_sol(1) % bnorm,&
             !                              10                    11
             solve_sol(1) % kfl_update_precond, solve_sol(1) % kappa, &
             !                      12                         13
             solve_sol(1) % lambda_min, solve_sol(1) % lambda_max

        call iofile_flush_unit(solve_sol(1) % lun_solve)
     end if
  end if

110 format(i7,8(2x,e16.8E3),2x,i12,13(2x,e16.8E3))

end subroutine solite
