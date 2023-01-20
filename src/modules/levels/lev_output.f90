!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine lev_output()
  !-----------------------------------------------------------------------
  !****f* Levels/lev_output
  ! NAME 
  !    lev_output
  ! DESCRIPTION
  !    End of a time step 
  !    ITASK=0 ... When timemarching is true. There is output or post-process
  !                of results if required.
  !    ITASK=1 ... When timemarching is false. Output and/or post-process of
  !                results is forced if they have not been written previously.
  ! USED BY
  !    lev_iniunk
  !    lev_endste
  !***
  !-----------------------------------------------------------------------
  use  def_parame
  use  def_master
  use  def_domain
  use  def_levels
  use  mod_output_postprocess, only : output_postprocess_variables
  implicit none
  external    :: lev_outvar
  !
  ! Postprocess on mesh
  !
  call output_postprocess_variables(lev_outvar)
  !
  ! Output interface mesh
  !
  call lev_output_interface_mesh()
  !
  ! Deallocate redistantiation
  !
  call lev_memall(10_ip)

  if( ittyp == ITASK_INITIA .or. ittyp == ITASK_ENDTIM ) then
     !
     ! Calculations on sets
     !
     call lev_outset()

  end if

end subroutine lev_output
