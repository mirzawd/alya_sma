!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine inisol()
  !-----------------------------------------------------------------------
  !****f* master/inisol
  ! NAME
  !    inisol
  ! DESCRIPTION
  !    This subroutine initializes the solver arrays
  !    AMATR ... Matrix
  !    RHSID ... Right-hand side
  !    PMATR ... Preconditioner matrix
  ! USES
  ! USED BY
  !    Alya
  !***
  !-----------------------------------------------------------------------
  use def_master 
  use def_solver
  implicit none  
  integer(ip)  :: izmat,izrhs,izpre

  if( INOTMASTER ) then 

     solve_sol(1) % kfl_assem = 0 ! Matrix is not assembled

     if( solve_sol(1) % kfl_cmplx == 0 ) then

        !----------------------------------------------------------------
        !
        ! Algebraic REAL solver
        !
        !----------------------------------------------------------------

        do izmat = 1,solve_sol(1) % nzmat 
           amatr(izmat) = 0.0_rp
        end do
        do izrhs = 1,solve_sol(1) % nzrhs 
           rhsid(izrhs) = 0.0_rp
        end do

        if( solve_sol(1) % kfl_preco == 3 ) then
           do izpre = 1,solve_sol(1) % nzpre 
              pmatr(izpre) = 0.0_rp
           end do
        end if

     else

        !----------------------------------------------------------------
        !
        ! Algebraic COMPLEX solver
        !
        !----------------------------------------------------------------

           do izmat = 1,solve_sol(1) % nzmat * solve_sol(1) % nrhss 
              amatx(izmat) = CMPLX( 0.0_rp, 0.0_rp, kind=rp )
           end do
           do izrhs = 1,solve_sol(1) % nzrhs * solve_sol(1) % nrhss
              rhsix(izrhs) = CMPLX( 0.0_rp, 0.0_rp, kind=rp ) 
           end do

           if( solve_sol(1) % kfl_preco == 3 ) then
              do izpre = 1,solve_sol(1) % nzpre 
                 pmatx(izpre) = CMPLX( 0.0_rp, 0.0_rp, kind=rp ) 
              end do
           end if

     end if

  end if

end subroutine inisol
 
