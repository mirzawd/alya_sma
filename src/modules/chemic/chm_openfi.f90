!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine chm_openfi(itask)
  !------------------------------------------------------------------------
  !****f* partis/chm_openfi
  ! NAME
  !    chm_openfi
  ! DESCRIPTION
  !    This subroutine gets ALL the file names and open them to be used by
  !    the module in two possible ways:
  !
  !    1. Recalling them from the environment, when Alya is launched
  !    encapsulated in a shell script, or
  !
  !    2. Composing the names out of the problem name which is given as argument
  !    when the binary file Alya is launched "naked".
  ! USES
  ! USED BY
  !    chm_turnon
  !------------------------------------------------------------------------
  use def_chemic,     only : fil_droplet_chm, fil_resou_chm, fil_time2_chm, fil_times_chm, lun_droplet_chm, lun_times_chm
  use def_parame,     only : zero, two
  use def_master,     only : INOTSLAVE, kfl_naked, kfl_rstar, modul, namda, exmod, namod
  use def_kintyp,     only : ip
  use mod_iofile,     only : iofile

  implicit none
  integer(ip), intent(in) :: itask
  character(150)          :: fil_spcvg_chm
  character(7)            :: statu
  character(11)           :: forma
  character(6)            :: posit

  if( INOTSLAVE ) then

     if( kfl_rstar == 2 ) then
        statu = 'old'
        forma = 'formatted'
        posit = 'append'
     else
        statu = 'unknown'
        forma = 'formatted'
        posit = 'asis'
     end if

     select case (itask)

     case (2_ip)
        !
        ! Open files needed occasionally
        !
        if( kfl_naked == 0 ) then
           call GET_ENVIRONMENT_VARIABLE('FOR1911',fil_times_chm)
           call GET_ENVIRONMENT_VARIABLE('FOR1912',fil_time2_chm)
           call GET_ENVIRONMENT_VARIABLE('FOR1931',fil_resou_chm)
           call GET_ENVIRONMENT_VARIABLE('FOR1932',fil_spcvg_chm)
           call GET_ENVIRONMENT_VARIABLE('FOR1941',fil_droplet_chm)
        else
           fil_times_chm   = adjustl(trim(namda))//'-time-step.'        //exmod(modul)//'.res'
           fil_time2_chm   = adjustl(trim(namda))//'-time-step-target.' //exmod(modul)//'.res'
           fil_resou_chm   = adjustl(trim(namda))//'.'                  //exmod(modul)//'.src'
           fil_spcvg_chm   = adjustl(trim(namda))//'-species.'          //exmod(modul)//'.cvg'
           fil_droplet_chm = adjustl(trim(namda))//'-droplet.'          //exmod(modul)//'.res'
        end if

        !
        ! Time step strategy
        !
        call iofile(zero,lun_times_chm,fil_times_chm,namod(modul)//' TIME STEP INFORMATION')
        call iofile(zero,lun_droplet_chm,fil_droplet_chm,namod(modul)//' DROPLETS INFORMATION')

     case (4_ip)

        call iofile(two,lun_times_chm,fil_times_chm,namod(modul)//' TIME STEP')
        call iofile(two,lun_droplet_chm,fil_droplet_chm,namod(modul)//' DROPLETS INFORMATION')

     end select

  end if

end subroutine chm_openfi

