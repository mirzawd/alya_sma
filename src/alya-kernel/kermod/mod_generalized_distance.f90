!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Kermod
!> @{
!> @file    mod_ker_generalized_distance.f90
!> @author  houzeaux
!> @date    2020-09-07
!> @brief   Generalized distance
!> @details Toolbox for computing generlized distance
!>          Compute the generalized distance to the wall via a 
!>          Poisson equation:
!>          1. Solve Lapl(f)=-1, with f=0 on wall
!>          2. d=sqrt[ grad(f)^2 +2*f ] - sqrt[grad(f)^2]
!>          See the following references:
!>          P.G. Tucker, Differential equation-based wall distance computation for
!>            DES and RANS, J. Comp. Phys. 190 (2003) 229-248.
!>          P.G. Tucker, Int. J. Numer. Fluids 33 (2000) 869.
!>          P.G. Tucker, Appl. Math. Model. 22 (1998) 293.
!>          
!-----------------------------------------------------------------------

module mod_generalized_distance

  use def_kintyp_basic,   only : ip,rp
  use def_kintyp_solvers, only : soltyp
  use def_master,         only : modul
  use def_master,         only : mem_modul
  use def_kermod,         only : ndivi
  use def_domain,         only : ndime,npoin
  use def_domain,         only : meshe
  use def_domain,         only : elmar
  use mod_gradie,         only : gradie
  use mod_solver,         only : solver_solve
  use mod_solver,         only : solver_initialize_matrix_and_rhs
  use mod_ADR,            only : ADR_assemble_laplacian
  use mod_memory,         only : memory_alloca
  use mod_memory,         only : memory_deallo
  implicit none
  private

  public :: generalized_distance_update
  public :: generalized_distance_assembly
  public :: generalized_distance_solution
  
contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-09-07
  !> @brief   Update
  !> @details Compute generalized distance from the solution of Lu=-1
  !> 
  !-----------------------------------------------------------------------

  subroutine generalized_distance_assembly(solve,A,b)

    type(soltyp), pointer, intent(inout) :: solve(:)
    real(rp),     pointer, intent(inout) :: A(:)
    real(rp),     pointer, intent(inout) :: b(:)
    integer(ip)                          :: ipoin
    real(rp),     pointer                :: rhs(:)

    nullify(rhs)
    call memory_alloca(mem_modul(1:2,modul),'RHS','lev_memall',rhs,npoin)
    do ipoin = 1,npoin
       rhs(ipoin) = -1.0_rp
    end do
    
    call solver_initialize_matrix_and_rhs(solve,A,b)
    call ADR_assemble_laplacian(meshe(ndivi),elmar,A,rhs,b)

    call memory_deallo(mem_modul(1:2,modul),'RHS','lev_memall',rhs)
     
  end subroutine generalized_distance_assembly
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-09-07
  !> @brief   Update
  !> @details Compute generalized distance from the solution of Lu=-1
  !> 
  !-----------------------------------------------------------------------

  subroutine generalized_distance_solution(solve,A,b,x)

    type(soltyp), pointer,             intent(inout) :: solve(:)
    real(rp),     pointer, contiguous, intent(inout) :: A(:)
    real(rp),     pointer, contiguous, intent(inout) :: b(:)
    real(rp),     pointer,             intent(inout) :: x(:)
    
    call solver_solve(solve,A,b,x)
    
  end subroutine generalized_distance_solution
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-09-07
  !> @brief   Update
  !> @details Compute generalized distance from the solution of Lf=-1
  !>          d=sqrt[ grad(f)^2 +2*f ] - sqrt[grad(f)^2]
  !> 
  !-----------------------------------------------------------------------

  subroutine generalized_distance_update(unkno,dista)

    real(rp),   pointer, intent(in)    :: unkno(:)
    real(rp),   pointer, intent(inout) :: dista(:)
    integer(ip)                        :: ipoin
    real(rp)                           :: fact1,fact2
    real(rp),   pointer                :: gradu(:,:)

    nullify(gradu)
    call memory_alloca(mem_modul(1:2,modul),'GRADU','lev_memall',gradu,ndime,npoin)
    call gradie(unkno,gradu)
    do ipoin = 1,npoin
       fact1        = dot_product(gradu(:,ipoin),gradu(:,ipoin))
       fact2        = fact1 + 2.0_rp*abs(unkno(ipoin))
       dista(ipoin) = sqrt(fact2) - sqrt(fact1) 
    end do
    call memory_deallo(mem_modul(1:2,modul),'GRADU','lev_memall',gradu)
    
  end subroutine generalized_distance_update
  
end module mod_generalized_distance
!> @}
