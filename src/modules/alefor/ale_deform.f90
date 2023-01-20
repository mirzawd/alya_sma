!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Alefor
!> @{
!> @file    ale_deform.f90
!> @author  houzeaux
!> @date    2020-04-21
!> @brief   Deform a mesh
!> @details This routines ale_deformes the mesh. Idea: Put a node at the middle
!>          of its neighbors. A Gauss Seidel in 1D would give:
!>
!>          X^{k+1}(i) = 1/2 [  X^k(i-1) + X^k(i+1) ]
!>          X^{k+1}(i) = X^k(i) + 1/2 [ X^k(i-1) -2 X^k(i) + X^k(i+1) ]
!>          X^{k+1}(i) = X^k(i) + 1/2 h^2 [ X^k(i-1) -2 X^k(i) + X^k(i+1) ]/h^2
!>          X^{k+1}(i) = X^k(i) + h^2 [ 0 - Lapl(X) ]
!>       
!>          Let X = x0+d. x0 is initial solution. d is displacement.
!>          Solve div[ a * grad(X) ] = 0 with X=x0 on boundary
!>        
!>          In 1D:
!>      
!>          i-1   i   i+1
!>           o----o----o
!>              h    h
!>          +-
!>          |  a* grad(X).grad(v) dx =   a*h* [X(i)-X(i-1)]/h (1/h) 
!>         -+                          + a*h* [X(i+1)-X(i)]/h (-1/h) 
!>                                   = -alpha*h* Lapla(X)
!>         => a = h
!> @} 
!-----------------------------------------------------------------------

subroutine ale_deform()
  use def_master
  use def_domain
  use def_alefor
  use mod_ker_deform, only : deform_deform
  implicit none
  
  real(rp),  pointer :: coord_loc(:,:)
  real(rp),  pointer :: bvess_loc(:,:)

  if( INOTEMPTY ) then
     coord_loc => coord_ale(:,:,1)
     bvess_loc => bvess_ale(:,:,1)
  else
     nullify(coord_loc)
     nullify(bvess_loc)
  end if
  
  call deform_deform(&  
       ndefo_ale,kfl_defor_ale,defor_param_ale,ID_ALEFOR,kfl_fixno_ale,bvess_loc,&
       coord_loc,amatr,unkno,rhsid,solve_sol)
  
end subroutine ale_deform
