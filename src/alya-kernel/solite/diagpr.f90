!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine diagpr(npoin,nbvar,an,ja,ia,wa1,wa2)

  !--------------------------------------------------------------------------------------------
  ! Sources/kernel/solite/diagpr.f90
  ! NAME
  !    diagpr
  ! DESCRIPTION
  !    This routine calculates simple complex Jacobi (diagonal) preconditioning: 
  !   
  !                             [M]^-1 = 1/diag([A])
  !                            
  ! INPUT ARGUMENTS
  !    NPOIN ..... Number of nodes in the interior of a mesh
  !    NBVAR ..... Number of unknowns in each node
  !    AN ........ Sparse matrix of the original system in BCSR (Blocked Compressed Sparse Row) format
  !    JA ........ Compressed Sparse format: index vector for column numbers
  !    IA ........ Compressed Sparse format: index vector for beginning of a row block
  ! OUTPUT ARGUMENTS 
  !    WA1 ....... Matrix of diagonal preconditioner, [M]
  !    WA2 ....... Inverse matrix of diagonal preconditioner, [M]^-1 
  !--------------------------------------------------------------------------------------------

  use def_kintyp, only       :  ip,rp
  use def_master, only       :  INOTMASTER
  use def_solver, only       :  solve_sol
  use mod_memchk                               

 !Declaration statements
  implicit none

 !Dummy arguments
  integer(ip), intent(in)     :: npoin,nbvar
  complex(rp), intent(in)     :: an(*)
  integer(ip), intent(in)     :: ja(*),ia(*)
  complex(rp), intent(inout)  :: wa1(*),wa2(*)

 !Local variables
  integer(ip)              :: ii,total
!  complex(rp)              :: norm

  if (INOTMASTER) then
    total = npoin * nbvar
  	!wa1 = M = diag(A)
  	call diagox(npoin,nbvar,solve_sol(1)%kfl_symme,ia,ja,an,wa1)
    !wa2 = M^(-1) = 1/diag(A)
	  do ii = 1,total
   	  if (wa1(ii) == (0.0_rp,0.0_rp)) then
        write (*,*) 'Muy mal!!!!!!!!!!!!!!! Zero element in diagonal!'
        wa1(ii) = (1.0_rp,0.0_rp)
        wa2(ii) = (1.0_rp,0.0_rp)      
     	else
     		wa2(ii) = (1.0_rp,0.0_rp)/wa1(ii)
     	endif
  	enddo
  endif

!  write (*,*) 'Subroutine DIAGPR: ', total
!  write (*,*) ' '
!  call no2plx(npoin,nbvar,wa1,norm)
!  norm = cmplx(0.0,0.0)
!  do ii = 1,total
!    norm = norm + wa1(ii)
!  end do
!  write (*,*) norm
!  call no2plx(npoin,nbvar,wa2,norm)
!  norm = cmplx(0.0,0.0)
!  do ii = 1,total
!    norm = norm + wa2(ii)
!  end do
!  write (*,*) norm

end subroutine diagpr

