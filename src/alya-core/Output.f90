!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine Output(itask)
  !-----------------------------------------------------------------------
  !****f* master/Output
  ! NAME
  !    Output
  ! DESCRIPTION
  !    This routine output and postprocess the solution
  ! USES
  !    moduls
  ! USED BY
  !    Alya
  !***
  !-----------------------------------------------------------------------
  use def_kintyp,             only : ip
  use def_master,             only : iblok
  use def_master,             only : nblok
  use def_master,             only : ittyp
  use def_master,             only : kfl_modul
  use def_master,             only : modul
  use def_master,             only : mmodu
  use def_master,             only : ITASK_OUTPUT
  use def_master,             only : postp
  use def_master,             only : momod
  use def_master,             only : ITASK_INITIA
  use def_master,             only : ITASK_ENDTIM 
  use def_master,             only : INOTSLAVE
  use mod_output_postprocess, only : output_postprocess_allocate_sets_and_witness
  use mod_output_postprocess, only : output_postprocess_header_sets_and_witness
  use mod_output_postprocess, only : output_postprocess_output_witness
  use mod_output_postprocess, only : output_postprocess_witness_parall
  use mod_moduls,             only : moduls
  use mod_witness,            only : witness_value_ini
  use mod_witness,            only : witness_value_end
  use mod_witness,            only : witness_value_increment_denom
#ifdef CATA
  use def_domain, only : coord, npoin, lnods, nelem, ltype
  use tcp
#endif

  implicit none

  integer(ip), intent(in) :: itask
  !
  ! Where we are
  ! 1 = initial solution
  ! 2 = end of a time step
  ! 3 = end of the run
  !
  ittyp = itask
  !
  ! Set and witness point headers
  !
  if( ittyp == ITASK_INITIA .and. INOTSLAVE ) then
     do modul = 1,mmodu
        if( kfl_modul(modul) /= 0 ) then
           postp => momod(modul) % postp
           call output_postprocess_header_sets_and_witness(momod(modul) % postp)
        end if
     end do
     modul = 0
  end if
  !
  ! Set and witness point headers
  !
  do modul = 1,mmodu
     if( kfl_modul(modul) /= 0 ) then
        call witness_value_ini(momod(modul) % postp)
     end if
  end do
  modul = 0
  !
  ! Output and postprocess of modules
  !
  call Kermod(-ITASK_OUTPUT)
  do iblok = 1,nblok
     call moduls(ITASK_OUTPUT)
  end do
  !
  ! Increment denominators
  !
  do modul = 1,mmodu
     if( kfl_modul(modul) /= 0 ) then
        call witness_value_increment_denom(momod(modul) % postp)
     end if
  end do
  modul = 0
  !
  ! Set and witness point output
  !
  if( ittyp ==  ITASK_INITIA .or. ittyp == ITASK_ENDTIM ) then
     do modul = 1,mmodu
        if( kfl_modul(modul) /= 0 ) then
           call output_postprocess_witness_parall(momod(modul) % postp) 
           call output_postprocess_output_witness(momod(modul) % postp)
        end if
     end do
     modul = 0
  end if
  !
  ! Averaging for witness
  !
  do modul = 1,mmodu
     if( kfl_modul(modul) /= 0 ) then
        call witness_value_end(momod(modul) % postp)
     end if
  end do
  modul = 0
  !
  !CATALYST
  !
#ifdef CATA
  call testcoprocessor(ittim,cutim,npoin,coord,nelem,lnods,ltype,kfl_paral,veloc,press)
#endif

end subroutine Output
