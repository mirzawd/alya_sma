!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine Endste()
  !-----------------------------------------------------------------------
  !****f* master/Endste
  ! NAME
  !    Endste
  ! DESCRIPTION
  !    This routine closes a time step.
  ! USES
  !    Nastin
  !    Temper
  !    Codire
  !    Alefor
  ! USED BY
  !    Alya
  !***
  !-----------------------------------------------------------------------

  use def_kintyp,              only : ip,rp
  use mod_bourgogne_pinotnoir, only : bourgogne
  use def_master,              only : iblok
  use def_master,              only : nblok
  use def_master,              only : ITASK_ENDSTE
  use def_master,              only : kfl_gotim
  use def_master,              only : cutim
  use def_master,              only : timef
  use def_master,              only : ittim
  use def_master,              only : mitim
  use def_master,              only : momod
  use def_master,              only : mmodu
  use def_master,              only : kfl_modul
  use mod_messages,            only : livinf
  use mod_moduls,              only : moduls
  use mod_optimum_partition,   only : optimum_partition_check
  use mod_state,               only : state_output
  
  implicit none
  integer(ip) :: imodu

  call bourgogne(1_ip)
  !
  ! Initializations
  !
  kfl_gotim = 0
  !
  ! End a time step for each module
  !
  do iblok = 1,nblok
    call moduls(ITASK_ENDSTE)
  end do
  call Kermod(ITASK_ENDSTE)
  !
  ! Save old time steps
  !
  call setgts(ITASK_ENDSTE)
  !
  ! Postprocess ppm
  !
  call posppm()
  !
  ! Live information
  !
  call livinf(9_ip,' ',0_ip)
  !
  ! Check if the time evolution has to be stopped or not
  !  
  if( kfl_gotim /= 0 ) then !Some module is still running, then wake up all modules
     do imodu = 1,mmodu-1
        if( kfl_modul(imodu) /= 0 ) then
           momod(imodu) % kfl_stead = 0
        end if
     end do
  end if
!!!!!!!kfl_gotim = 1!!!OJO Cambiado para que alefor funcione solo y haga iters de t

  if( cutim >= timef-epsilon(1.0_rp) )  kfl_gotim = 0

  if( timef < 0.0_rp )                  kfl_gotim = 1
  if( ittim >= mitim )                  kfl_gotim = 0
  !
  ! Write restart files
  !
  call restar(2_ip) ! General data
  !
  ! Output memory evolution
  !
  call output_memory_evolution()

  !
  ! Check if partition is optimum
  !
  call optimum_partition_check()
  !
  ! State output
  !
  call state_output('RUNNING')
  
end subroutine Endste
