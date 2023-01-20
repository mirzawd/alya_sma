!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine lev_openfi(itask)

  !-----------------------------------------------------------------------
  !    
  ! This subroutine gets ALL the file names and open them to be used by 
  ! the module in two possible ways:
  ! 
  ! 1. Recalling them from the environment, when Alya is launched
  ! encapsulated in a shell script, or
  ! 
  ! 2. Composing the names out of the problem name which is given as argument
  ! when the binary file Alya is launched "naked". 
  !
  !-----------------------------------------------------------------------
  use      def_levels
  use      def_parame
  use      def_master
  use      mod_iofile
  implicit none
  integer(ip), intent(in) :: itask 
  character(150)          :: fil_capte_lev,fil_volum_lev,fil_reave_lev
  character(150)          :: fil_inter_lev, fil_gauge_lev, fil_cored_lev
  character(150)          :: outputname
  character(7)            :: statu
  character(11)           :: forma
  character(6)            :: posit

  if( INOTSLAVE ) then

     if(kfl_rstar==2) then
        statu='old'
        forma='formatted'
        posit='append'
     else
        statu='unknown'
        forma='formatted'
        posit='asis'
     end if

     select case (itask)

     case (2)
        !
        ! Open files needed occasionally
        !
        if (kfl_naked==0) then
           call GET_ENVIRONMENT_VARIABLE('FOR1411',fil_reave_lev)
           call GET_ENVIRONMENT_VARIABLE('FOR1409',fil_volum_lev)
           call GET_ENVIRONMENT_VARIABLE('FOR1410',fil_capte_lev)
           if(npp_gauge_lev==1) then
              call GET_ENVIRONMENT_VARIABLE('FOR1413',fil_gauge_lev)
           end if
           if(inred_lev==1) then
              call GET_ENVIRONMENT_VARIABLE('FOR1413',fil_cored_lev)
           end if
        else
           fil_reave_lev = adjustl(trim(namda))//'-velocity.' //exmod(modul)//'.res'
           if(npp_gauge_lev==1) then
              fil_gauge_lev = adjustl(trim(namda))//'-gauges.' //exmod(modul)//'.res'
           end if
           if(tyred_lev>0) then
              fil_cored_lev = adjustl(trim(namda))//'-redist.' //exmod(modul)//'.cvg'
           end if
           fil_volum_lev = adjustl(trim(namda))//'-volume.'//exmod(modul)//'.res'
           fil_capte_lev = adjustl(trim(namda))//'-gauge.' //exmod(modul)//'.res'
        end if
        ! Read velocity initialisation file
        if(kfl_reave_lev==1) &
             call iofile(zero,lun_reave_lev,fil_reave_lev,namod(modul)//' READ INITIAL VELOCITY ',statu,forma,posit)

        ! Postprocess Interface Gauges
        if(npp_gauge_lev==1) then
           call iofile(zero,lun_gauge_lev,fil_gauge_lev,namod(modul)//' POSTPROCESS INTERFACE GAUGES ',statu,forma,posit)
        end if

        ! Postprocess convergence of redistanciation
        if(tyred_lev>0) then
           call iofile(zero,lun_cored_lev,fil_cored_lev,namod(modul)//' POSTPROCESS REDISTANCIATION CONVERGENCE ',statu,forma,posit)
        end if

        ! Volume
        call iofile(zero,lun_volum_lev,fil_volum_lev,namod(modul)//' LIQUID VOLUME',statu,forma,posit)

        ! Capte
        call iofile(zero,lun_capte_lev,fil_capte_lev,namod(modul)//' INTERPHASE GAUGE',statu,forma,posit)

     case(4)
        !
        ! Close result file
        !
        close(lun_volum_lev)
        close(lun_capte_lev)
        if(npp_gauge_lev==1) then
           close(lun_gauge_lev)
        endif
        if(tyred_lev>0) then
           close(lun_cored_lev)
        end if

     case(5)

        statu='unknown'
        forma='formatted'
        posit='asis'
        Write(outputname,"('-interf.',i4.4,'.lev.res')") ittim
        call GET_ENVIRONMENT_VARIABLE('FOR1412',fil_inter_lev)
        fil_inter_lev = adjustl(trim(namda))//outputname 
        call iofile(zero,lun_inter_lev,fil_inter_lev,namod(modul)//' POSTPROCESS INTERFACE ',statu,forma,posit)

     case (8_ip)

        close(lun_inter_lev)

     end select

  end if

end subroutine lev_openfi

