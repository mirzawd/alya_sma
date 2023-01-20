!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Kermod
!> @{
!> @file    ker_outerr.f90
!> @author  Guillaume Houzeaux
!> @date    19/02/2016
!> @brief   Check errors
!> @details Check errors
!> @} 
!-----------------------------------------------------------------------

subroutine ker_outerr()

#include "def_vector_size.inc"
  use def_master
  use def_kermod
  use def_domain
  use def_solver
  use mod_communications, only : PAR_SUM
  use mod_outfor,         only : outfor
  use mod_mpio_config,    only : mpio_config
  implicit none
  integer(ip) :: ierro,iwarn
  integer(ip) :: ifw, id

  ierro = 0
  iwarn = 0
  !
  ! Direct solver
  ! 
  if( kfl_direct_solver == SOL_DIRECT_SOLVER_PASTIX ) then
#ifndef PASTIX 
     ierro = ierro + 1
     call outfor(1_ip,momod(modul)%lun_outpu,'COMPILE ALYA WITH -DPASTIX AND LINK IT WITH PASTIX')     
#endif
  end if
  if( kfl_direct_solver == SOL_DIRECT_SOLVER_MUMPS ) then
#ifndef MUMPS
     ierro = ierro + 1 
     call outfor(1_ip,momod(modul)%lun_outpu,'COMPILE ALYA WITH -DMUMPS AND LINK IT WITH MUMPS')     
#endif
  end if
  if( kfl_direct_solver == SOL_DIRECT_SOLVER_WSMP ) then
#ifndef WSMP
     ierro = ierro + 1
     print*,'ierro',ierro
     call outfor(1_ip,momod(modul)%lun_outpu,'COMPILE ALYA WITH -DWSMP AND LINK IT WITH WSMP')
#endif
  end if
  if( kfl_direct_solver == SOL_DIRECT_SOLVER_PWSMP) then
#ifndef PWSMP
     ierro = ierro + 1
     call outfor(1_ip,momod(modul)%lun_outpu,'COMPILE ALYA WITH -DPWSMP AND LINK IT WITH PWSMP')
#endif
  end if

  !
  ! Periodicity does not work with SFC 
  !
  if( nperi/=0_ip .and. kfl_ngrou /= -1_ip .and. kfl_ngrou /= 0_ip .and. new_periodicity == 0 ) then
     ierro = ierro + 1
     call outfor(1_ip,momod(modul)%lun_outpu,&
          'FOR PERIODIC CASES YOU NEED TO USE GROUPS = int, SEQUENTIAL_FRONTAL IN THE DOM.DAT')
  end if
  !
  ! Exchange location needs elsest 
  !
  if( kfl_waexl_ker/=0_ip .and. kfl_elses == 0_ip  ) then
     ierro = ierro + 1
     call outfor(1_ip,momod(modul)%lun_outpu,&
          'EXCHANGE LOCATION NEEDS ELSEST')
  end if
  !
  ! MPIO incompatibilities 
  !
  if( kfl_posdi == 0 .and. mpio_config%output%merge .and. ndivi > 0 .and.&
       ( mpio_config%output%post_process%enabled .or. mpio_config%output%restart%enabled ) ) then
     call runend('KER_OUTERR: CANNOT POSTPROCESS OR RESTART ON ORIGINAL MESH WHEN USING MPIO WITH MERGE OPTION')
  end if

  !----------------------------------!
  ! ERRORS RELATED LOOKUP FRAMEWORKS !
  !----------------------------------!
  !
  ! Undefined scaling type 
  !
  do ifw = 1,max_lookup_fw
     !
     ! Check if this framwork exists
     !
     if (lookup_fw(ifw) % kfl_tab_main > 0_ip) then
        !
        ! Check every dimension
        ! 
        do id = 1,lookup_fw(ifw) % main_table % ndim
           !
           ! Check if un-initialized
           !
           if (lookup_fw(ifw) % kfl_scale(id) == -2_ip) then
              ierro = ierro + 1
              call outfor(1_ip,momod(modul)%lun_outpu, &
                 'IN LOOKUP FRAMEWORK: '&
                 //trim(intost(ifw))//&
                 ' THE SCALING LAW OF DIMENSION '&
                 //trim(intost(id))//&
                 '('//lookup_fw(ifw) % main_table % coords(id) % name//') IS UNINITIALIZED')
           endif
        enddo
     endif
  enddo










  call PAR_SUM(ierro,'IN MY CODE',INCLUDE_ROOT=.true.)
  if( ierro /= 0 ) call runend('AN ERROR HAS BEEN FOUND IN KERMOD')

end subroutine ker_outerr
