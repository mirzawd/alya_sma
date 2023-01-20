!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



module mod_bourgogne_pinotnoir

contains
  !DEC$ ATTRIBUTES FORCEINLINE :: bourgogne
  subroutine bourgogne(itask)

    use def_kintyp, only : ip
#ifdef BOURGOGNE
    use def_kintyp, only : rp
    use def_master, only : INOTSLAVE,cpu_initi
#endif
    use def_parall
    implicit none
    integer(ip), intent(in) :: itask
    !CHARACTER(50)           :: s = "Alya-License.lic"
#ifdef BOURGOGNE
    CHARACTER(500)          :: licenseName
    INTEGER(4)              :: hourCounter = 50
    REAL(rp)                :: caca
    LOGICAL                 :: license_exist

    if( INOTSLAVE ) then
      !CALL get_command_argument(1, licenseName
      !licenseName = s

      CALL get_command_argument(2, licenseName)

      INQUIRE(FILE=licenseName, EXIST=license_exist)
      if(license_exist) then


        if( itask == 1 ) then
          call check(licenseName)

        else
          !Actualize License
          call cputim(caca)
          caca = ((caca - cpu_initi))
          call actualize(licenseName, caca)     
        end if
      else
        STOP "License FILE do not exist"
      end if
    end if

#endif

  end subroutine bourgogne

end module mod_bourgogne_pinotnoir
