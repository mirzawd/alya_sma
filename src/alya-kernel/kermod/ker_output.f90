!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine ker_output()
  !------------------------------------------------------------------------
  !****f* Kermod/ker_output
  ! NAME 
  !    ker_output
  ! DESCRIPTION
  !    Output module variables
  ! USES
  !    ker_outvar
  ! USED BY
  !    ker_endste (itask=1)
  !    ker_turnof (itask=2)
  !***
  !------------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use mod_witness
  use mod_output_postprocess
  implicit none
  external    :: ker_outvar
  !
  ! Create and output witness meshes
  !
  if( ittyp == ITASK_INITIA ) then
     call witness_mesh_create()   
  end if
  !
  ! Postprocess variables
  !
  call output_postprocess_variables(ker_outvar)
    
end subroutine ker_output
