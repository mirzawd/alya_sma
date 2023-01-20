!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine exm_endste
  !-----------------------------------------------------------------------
  !****f* Exmedi/exm_endste
  ! NAME 
  !    exm_endste
  ! DESCRIPTION
  !    This routine ends a time step of the incompressible NS equations.
  ! USES
  !    exm_output
  ! USED BY
  !    exm_endste (itask=1)
  !    exm_turnof (itask=2)
  !***
  !-----------------------------------------------------------------------
  use      def_parame
  use      def_master
  use      def_domain

  use      def_exmedi
  use      mod_exm_activation, only : exm_stim_synchronize

  implicit none


  call exm_updunk(ITASK_ENDSTE)         
  call exm_isochr( 1_ip)

  kfl_gotim = 1
  
  call exm_stim_synchronize()

end subroutine exm_endste
