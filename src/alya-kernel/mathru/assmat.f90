!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine assmat(&
     ndofn,pnode,pevat,nunkn,kfl_algso,&
     ielem,lnods,elmat,amatr)
  !-----------------------------------------------------------------------
  !****f* mathru/assmat
  ! NAME 
  !    assmat
  ! DESCRIPTION
  !    Assembly an elemental matrix ELMAT in global matrix AMATR
  ! USES
  ! USED BY
  !    ***_elmope
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp
  use def_solver, only       :  solve_sol
  use def_solver, only       :  nzdom_aii,nzdom_aib,nzdom_abi
  implicit none
  integer(ip), intent(in)    :: ndofn,pnode,pevat,nunkn
  integer(ip), intent(in)    :: kfl_algso
  integer(ip), intent(in)    :: ielem
  integer(ip), intent(in)    :: lnods(pnode)
  real(rp),    intent(in)    :: elmat(ndofn,ndofn,pnode,pnode)
  real(rp),    intent(inout) :: amatr(*)
  integer(ip)                :: poaii,poaib,poabi,poabb

  solve_sol(1) % kfl_assem = 1 ! Equation is assembled

  if( kfl_algso == 0 ) then
     !
     ! Direct LDU solver
     ! 
     call runend('ASSMAT: OBSOLETE OPTION')

  else if( kfl_algso /= -333333 ) then
     !
     ! Sparse matrix based solvers
     !
     if( solve_sol(1)%kfl_symme == 1 ) then
        !
        ! Symmetric assembly
        !
        call csrase(elmat,amatr,ndofn,pnode,pevat,&
             lnods,3_ip)
     else
        !
        ! General case
        !
        if( solve_sol(1) % kfl_schur == 1 ) then
           poaii = 1
           poaib = poaii + nzdom_aii
           poabi = poaib + nzdom_aib
           poabb = poabi + nzdom_abi
           call csrshu(elmat,Amatr(poaii),Amatr(poaib),Amatr(poabi),Amatr(poabb),&
                ndofn,pnode,pnode,lnods,2_ip)
        else
           !
           ! No element to CSR straucture
           !
           call csrase(elmat,amatr,ndofn,pnode,pevat,lnods,2_ip)
        end if
     end if
  end if

end subroutine assmat
