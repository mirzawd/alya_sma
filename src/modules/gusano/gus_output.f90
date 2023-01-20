!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Gusano
!> @{
!> @file    gus_output.f90
!> @author  houzeaux
!> @date    2020-10-20
!> @brief   Output
!> @details Gusano output
!> @} 
!-----------------------------------------------------------------------

subroutine gus_output()

  use def_parame
  use def_master
  use def_domain
  use def_gusano
  use mod_output_postprocess

  implicit none
  external          :: gus_outvar
  !
  ! Postprocess variables
  !
  call output_postprocess_variables(gus_outvar)
  !
  ! Calculations on sets
  !
  if( ittyp == ITASK_INITIA .or. ittyp == ITASK_ENDTIM ) then
     call gus_outset()
  end if

end subroutine gus_output
