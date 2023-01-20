!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



module mod_module_interface

  use def_kintyp,          only : ip
  implicit none
  public

contains

#ifdef CMAKE

  subroutine nastin(itask)
#ifdef NASTIN_MODULE
    NASTIN_MODULE mod_nastin
#endif
    integer(ip), intent(in) :: itask
#ifdef NASTIN_MODULE
    call nastin_main(itask)
#else
    call runend("NASTIN IS NOT AVAILABLE! PLEASE COMPILE THE NASTIN LIBRARY")
#endif
  end subroutine

  subroutine temper(itask)
#ifdef TEMPER_MODULE
    TEMPER_MODULE mod_temper
#endif
    integer(ip), intent(in) :: itask
#ifdef TEMPER_MODULE
    call temper_main(itask)
#else
    call runend("TEMPER IS NOT AVAILABLE! PLEASE COMPILE THE TEMPER LIBRARY")
#endif
  end subroutine
  
  subroutine partis(itask)
#ifdef PARTIS_MODULE
    PARTIS_MODULE mod_partis
#endif
    integer(ip), intent(in) :: itask
#ifdef PARTIS_MODULE
    call partis_main(itask)
#else
    call runend("PARTIS IS NOT AVAILABLE! PLEASE COMPILE THE PARTIS LIBRARY")
#endif
  end subroutine

  subroutine alefor(itask)
#ifdef ALEFOR_MODULE
    ALEFOR_MODULE mod_alefor
#endif
    integer(ip), intent(in) :: itask
#ifdef ALEFOR_MODULE
    call alefor_main(itask)
#else
    call runend("ALEFOR IS NOT AVAILABLE! PLEASE COMPILE THE ALEFOR LIBRARY")
#endif
  end subroutine

  subroutine solidz(itask)
#ifdef SOLIDZ_MODULE
    SOLIDZ_MODULE mod_solidz
#endif
    integer(ip), intent(in) :: itask
#ifdef SOLIDZ_MODULE
    call solidz_main(itask)
#else
    call runend("SOLIDZ IS NOT AVAILABLE! PLEASE COMPILE THE SOLIDZ LIBRARY")
#endif
  end subroutine
  
  subroutine gusano(itask)
#ifdef GUSANO_MODULE
    GUSANO_MODULE mod_gusano
#endif
    integer(ip), intent(in) :: itask
#ifdef GUSANO_MODULE
    call gusano_main(itask)
#else
    call runend("GUSANO IS NOT AVAILABLE! PLEASE COMPILE THE GUSANO LIBRARY")
#endif
  end subroutine

  subroutine exmedi(itask)
#ifdef EXMEDI_MODULE
    EXMEDI_MODULE mod_exmedi
#endif
    integer(ip), intent(in) :: itask
#ifdef EXMEDI_MODULE
    call exmedi_main(itask)
#else
    call runend("EXMEDI IS NOT AVAILABLE! PLEASE COMPILE THE EXMEDI LIBRARY")
#endif
  end subroutine

  subroutine chemic(itask)
#ifdef CHEMIC_MODULE
    CHEMIC_MODULE mod_chemic
#endif
    integer(ip), intent(in) :: itask
#ifdef CHEMIC_MODULE
    call chemic_main(itask)
#else
    call runend("CHEMIC IS NOT AVAILABLE! PLEASE COMPILE THE CHEMIC LIBRARY")
#endif
  end subroutine

  subroutine turbul(itask)
#ifdef TURBUL_MODULE
    TURBUL_MODULE mod_turbul
#endif
    integer(ip), intent(in) :: itask
#ifdef TURBUL_MODULE
    call turbul_main(itask)
#else
    call runend("TURBUL IS NOT AVAILABLE! PLEASE COMPILE THE TURBUL LIBRARY")
#endif
  end subroutine

  subroutine neutro(itask)
#ifdef NEUTRO_MODULE
    NEUTRO_MODULE mod_neutro
#endif
    integer(ip), intent(in) :: itask
#ifdef NEUTRO_MODULE
    call neutro_main(itask)
#else
    call runend("NEUTRO IS NOT AVAILABLE! PLEASE COMPILE THE NEUTRO LIBRARY")
#endif
  end subroutine

  subroutine levels(itask)
#ifdef LEVELS_MODULE
    LEVELS_MODULE mod_levels
#endif
    integer(ip), intent(in) :: itask
#ifdef LEVELS_MODULE
    call levels_main(itask)
#else
    call runend("LEVELS IS NOT AVAILABLE! PLEASE COMPILE THE LEVELS LIBRARY")
#endif
  end subroutine

  subroutine magnet(itask)
#ifdef MAGNET_MODULE
    MAGNET_MODULE mod_magnet
#endif
    integer(ip), intent(in) :: itask
#ifdef MAGNET_MODULE
    call magnet_main(itask)
#else
    call runend("MAGNET IS NOT AVAILABLE! PLEASE COMPILE THE MAGNET LIBRARY")
#endif
  end subroutine

#endif

end module mod_module_interface
