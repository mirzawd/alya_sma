!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



recursive subroutine tem_elmadr()
  !------------------------------------------------------------------------
  !****f* Temper/tem_elmadr
  ! NAME 
  !    tem_elmadr
  ! DESCRIPTION
  !    ORDER=1:
  !      Temperature equation, elemental operations:
  !      1. Compute elemental matrix and RHS 
  !      2. Impose Dirichlet boundary conditions
  !      3. Assemble them
  !    ORDER=2:
  !      Update the subgrid scale
  ! USES
  ! USED BY
  !    tem_matrix
  !------------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_temper
  implicit none

end subroutine tem_elmadr
