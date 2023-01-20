!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine reaeig(itask)
  !-----------------------------------------------------------------------
  !****f* outrut/reaeig
  ! NAME 
  !    reaeig
  ! DESCRIPTION
  !    This routine reads the eigen solver data
  ! INPUT
  ! OUTPUT
  ! USES
  ! USED BY
  !    ***_reanut
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only    :  ip,rp
  use def_solver, only    :  eigen_sol 
  use def_master, only    :  kfl_symgr
  use def_inpout
  implicit none
  integer(ip), intent(in) :: itask

  if( itask == 1 ) then
     !
     ! Solver type
     !
     if(exists('DIREC')) eigen_sol(1)%kfl_algso = 0
     if(exists('EIGE1')) eigen_sol(1)%kfl_algso = 1
     if(exists('EIGE2')) eigen_sol(1)%kfl_algso = 2
     if(exists('EIGE3')) eigen_sol(1)%kfl_algso = 3
     !
     ! Solver name
     !
     if(exists('DIREC')) eigen_sol(1)%wsolv = 'DIRECT'
     if(exists('EIGE1')) eigen_sol(1)%wsolv = 'ITERATIVE 1'
     if(exists('EIGE2')) eigen_sol(1)%wsolv = 'ITERATIVE 2'
     if(exists('EIGE3')) eigen_sol(1)%wsolv = 'ITERATIVE 3'
     !
     ! Solver parameters
     !
     if(exists('ITERA')) eigen_sol(1)%miter = getint('ITERA', 1_ip,   '#Max iterations')
     if(exists('TOLER')) eigen_sol(1)%solco = getrea('TOLER', 1e-4_rp,'#Tolerance')
     if(exists('NUMBE')) eigen_sol(1)%neiva = getint('NUMBE', 10_ip,  '#Number of eigenvalues')
     if(exists('SHIFT')) eigen_sol(1)%shift = getrea('SHIFT',-1.0_rp,   '#amount of shift')
     !
     ! Mass matrix
     !
     if(exists('MASS ')) then
        if(exists('DIAGO').or.exists('LUMPE')) then
           eigen_sol(1)%kfl_massm = 0
        else
           eigen_sol(1)%kfl_massm = 1
           kfl_symgr = 1
        end if
     end if

  end if

end subroutine reaeig
