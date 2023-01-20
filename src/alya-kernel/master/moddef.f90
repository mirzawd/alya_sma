!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine moddef(itask)
  !-----------------------------------------------------------------------
  !****f* master/moddef
  ! NAME
  !    moddef
  ! DESCRIPTION
  !    This subroutine define units, file names, open and close files
  !
  !    ------------------------------------------------------------
  !    Units #      Unit name     File name     Description
  !    ------------------------------------------------------------
  !
  !    MOMOD(:) %
  !
  !    1 .......... lun_pdata ... fil_pdata ... Data
  !    2 .......... lun_outpu ... fil_outpu ... Output
  !    3 .......... lun_conve ... fil_conve ... Convergence
  !    5 .......... lun_rstar ... fil_rstar ... Restart
  !    11.......... lun_time  ... fil_time  ... Time values
  !
  !    POSTP(1) %
  !
  !    6 .......... lun_setse ... fil_setse ... Element set
  !    7 .......... lun_setsb ... fil_setsb ... Boundary set
  !    8 .......... lun_setsn ... fil_setsn ... Node set
  !    9 .......... lun_setsi ... fil_setsi ... IB set
  !    26 ......... lun_witne ... fil_witne ... Witness point
  !    27 ......... lun_witng ... fil_witng ... Witness geometry
  !
  !    SOLVE(:) %
  !
  !    40 -> 49 ... lun_exsol ... fil_exsol ... External solver output
  !    50 -> 69 ... lun_solve ... fil_solve ... Solver output
  !    70 -> 89 ... lun_cvgso ... fil_cvgso ... Solver convergence
  !
  !    ------------------------------------------------------------
  !
  ! USES
  ! USED BY
  !    Alya
  !***
  !-----------------------------------------------------------------------
  
  use def_parame
  use def_domain
  use def_master
  use def_kermod
  use def_solver
  use def_inpout
  use mod_memchk
  use mod_memory
  use mod_iofile
  use mod_parall,         only : PAR_MY_CODE_RANK_WM
  use mod_communications, only : PAR_AVERAGE
  use mod_communications, only : PAR_MAX
  use mod_outfor,         only : outfor
  use mod_messages,       only : messages_live
  use mod_moduls_conf,    only : moduls_set_current_module
  use mod_strings,        only : integer_to_string
  use mod_run_config,     only : run_config

  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: ivari,jtask,imod0,imod1,imod2,nsolv,ierror
  character(5)            :: wmodu,wfort(mmodu)
  character(20)           :: wtmp
  character(7)            :: statu
  character(11)           :: forma
  character(6)            :: posit
  character(150)          :: fil_rstar_new
  character(150)          :: fil_rstib_new
  character(150)          :: fil_exsol          ! Extrenal solver file
#ifdef _CRAYFTN
  integer(ip),parameter   :: lun_base = 200_ip
#else
  integer(ip),parameter   :: lun_base = 100_ip
#endif

  jtask = abs(itask)
  if( itask > 0 ) then
     !
     ! Treat all modules
     !
     imod0 = 0
     imod1 = 1
     imod2 = mmodu
  else
     !
     ! Treat only module MODUL
     !
     imod0 = modul
     imod1 = modul
     imod2 = modul
  end if

  if ( ( jtask == 1_ip .or. jtask == 2 ) ) then
     !
     ! Define unit opening option if this is a restart run
     !
     if ( kfl_rstar == 2 ) then
        statu = 'old'
        forma = 'formatted'
        posit = 'append'
     else
        statu = 'unknown'
        forma = 'formatted'
        posit = 'asis'
     end if
     do modul = imod1,imod2
        wtmp         = intost(modul)
        wmodu        = wtmp(1:5)
        wfort(modul) = 'FOR'//trim(wmodu)
     end do
  end if

  if ( jtask == 1_ip .and. INOTSLAVE ) then

     !-------------------------------------------------------------------
     !
     ! Open files that are always needed and get file names
     !
     !-------------------------------------------------------------------

     do modul = imod1,imod2

        if( kfl_modul(modul) /= 0 ) then

           postp => momod(modul) % postp
           solve => momod(modul) % solve
           eigeg => momod(modul) % eigen
           !
           ! Memory check
           !
           if( .not. associated(postp) ) call runend('MODDEF: POSTP NOT ASSOCIATED FOR MODULE '//integer_to_string(modul))
           !
           ! Define all units
           !
           momod(modul) % lun_pdata = modul * lun_base +  1     ! Data
           momod(modul) % lun_outpu = modul * lun_base +  2     ! Output
           momod(modul) % lun_conve = modul * lun_base +  3     ! Convergence
           momod(modul) % lun_rstpo = modul * lun_base +  4     ! Restart variables
           momod(modul) % lun_rstar = modul * lun_base +  5     ! Restart
           momod(modul) % lun_timin = modul * lun_base + 10     ! Timings
           postp(1)     % lun_setse = modul * lun_base +  6     ! Element set
           postp(1)     % lun_setsb = modul * lun_base +  7     ! Boundary set
           postp(1)     % lun_setsn = modul * lun_base +  8     ! Node set
           postp(1)     % lun_setsi = modul * lun_base +  9     ! IB set
           postp(1)     % lun_witne = modul * lun_base + 26     ! Witness point
           postp(1)     % lun_witng = modul * lun_base + 27     ! Witness geometry

           if ( kfl_naked == 0 ) then
              !
              !  Encapsulated, then get names from the environment
              !
              call GET_ENVIRONMENT_VARIABLE( trim(wfort(modul)) // '01', momod(modul) % fil_pdata )
              call GET_ENVIRONMENT_VARIABLE( trim(wfort(modul)) // '02', momod(modul) % fil_outpu )
              call GET_ENVIRONMENT_VARIABLE( trim(wfort(modul)) // '03', momod(modul) % fil_conve )
              call GET_ENVIRONMENT_VARIABLE( trim(wfort(modul)) // '05', momod(modul) % fil_rstar )
              call GET_ENVIRONMENT_VARIABLE( trim(wfort(modul)) // '10', momod(modul) % fil_timin )
              !call GET_ENVIRONMENT_VARIABLE( trim(wfort(modul)) // '11', momod(modul) % fil_time  )

           else if ( kfl_naked == 1 ) then
              !
              !  Naked, then compose the names
              !
              momod(modul) % fil_pdata = adjustl(trim(namda))//'.'        //exmod(modul)//'.dat'
              momod(modul) % fil_outpu = adjustl(trim(namda))//'.'        //exmod(modul)//'.log'
              momod(modul) % fil_conve = adjustl(trim(namda))//'.'        //exmod(modul)//'.cvg'
              momod(modul) % fil_rstar = adjustl(trim(namda))//'.'        //exmod(modul)//'.rst'
              momod(modul) % fil_timin = adjustl(trim(namda))//'-timings.'//exmod(modul)//'.res'
              !momod(modul) % fil_time  = adjustl(trim(namda))//'.'//exmod(modul)//'.tim'

           end if

           if( modul == mmodu ) then
              kfl_reawr = 1
              call iofile(4_ip,momod(modul) % lun_pdata,momod(modul) % fil_pdata,namod(modul)//' DATA',       'old')
              if( kfl_reawr == -1 ) momod(modul) % lun_pdata = 0
              kfl_reawr = 0
           end if
           if( momod(modul) % lun_pdata > 0 ) then
              call iofile(zero,momod(modul) % lun_pdata,momod(modul) % fil_pdata,namod(modul)//' DATA',       'old')
           end if
           !
           ! Output (log) and convergence (cvg) and time (tim) files
           !
           if( kfl_rstar == 2 ) then
              kfl_reawr = 1
              call iofile(4_ip,momod(modul) % lun_outpu,momod(modul) % fil_outpu,namod(modul)//' OUTPUT' ,    statu,forma,posit)
              if( kfl_reawr == 1 ) then
                 call iofile(zero,momod(modul) % lun_outpu,momod(modul) % fil_outpu,namod(modul)//' OUTPUT', 'old','formatted','append')
              else
                 call iofile(zero,momod(modul) % lun_outpu,momod(modul) % fil_outpu,namod(modul)//' OUTPUT')
              end if
              kfl_reawr = 1
              call iofile(4_ip,momod(modul) % lun_conve,momod(modul) % fil_conve,namod(modul)//' CONVERGENCE',     statu,forma,posit)
              if( kfl_reawr == 1 ) then
                 call iofile(zero,momod(modul) % lun_conve,momod(modul) % fil_conve,namod(modul)//' CONVERGENCE','old','formatted','append')
              else
                 call iofile(zero,momod(modul) % lun_conve,momod(modul) % fil_conve,namod(modul)//' CONVERGENCE')
              end if
              kfl_reawr = 1
              call iofile(4_ip,momod(modul) % lun_timin,momod(modul) % fil_timin,namod(modul)//' TIMINGS',    statu,forma,posit)
              if( kfl_reawr == 1 ) then
                 call iofile(zero,momod(modul) % lun_timin,momod(modul) % fil_timin,namod(modul)//' TIMINGS','old','formatted','append')
              else
                 call iofile(zero,momod(modul) % lun_timin,momod(modul) % fil_timin,namod(modul)//' TIMINGS')
              end if
              kfl_reawr = 0
           else
              call iofile(zero,momod(modul) % lun_outpu,momod(modul) % fil_outpu,namod(modul)//' OUTPUT',     statu,forma,posit)
              call iofile(zero,momod(modul) % lun_conve,momod(modul) % fil_conve,namod(modul)//' CONVERGENCE',statu,forma,posit)
              call iofile(zero,momod(modul) % lun_timin,momod(modul) % fil_timin,namod(modul)//' TIMINGS'    ,statu,forma,posit)
           end if
           !
           ! Write header
           !
           call outfor(12_ip,momod(modul) % lun_outpu)
           call outfor(98_ip,momod(modul) % lun_timin)

        end if
     end do

     call moduls_set_current_module(0_ip)

  else if ( jtask == 2 ) then

     !-------------------------------------------------------------------
     !
     ! Close data file and open files needed occasionally
     !
     !-------------------------------------------------------------------

     if( INOTSLAVE ) then

        do modul = imod1,imod2
           if ( kfl_modul(modul) /= 0 ) then
              if( momod(modul) % lun_pdata /= 0 ) &
                   call iofile(two,momod(modul) % lun_pdata,momod(modul) % fil_pdata,namod(modul)//' DATA')
           end if
        end do

     end if

     do modul = imod1,imod2

        if ( kfl_modul(modul) /= 0 ) then

           postp => momod(modul) % postp
           solve => momod(modul) % solve
           eigeg => momod(modul) % eigen

           if( associated(momod(modul) % solve) ) then
              nsolv = size(solve,kind=ip)
           else
              nsolv = 0
           end if

           if( associated(momod(modul) % solve) ) then
              do ivari = 1,nsolv
                 if ( solve(ivari) % kfl_algso /= -999 ) then
                    if ( solve(ivari) % kfl_solve /= 0 ) then
                       solve(ivari) % lun_solve = modul*lun_base + 49 + ivari ! 50 -> 69
                    end if
                    if ( solve(ivari) % kfl_cvgso /= 0 ) then
                       solve(ivari) % lun_cvgso = modul*lun_base + 69 + ivari ! 70 -> 89
                    end if
                    if(      solve(ivari) % kfl_algso == SOL_SOLVER_MAPHYS_UNSYMMETRIC .or. &
                         &   solve(ivari) % kfl_algso == SOL_SOLVER_MAPHYS_SYMMETRIC   .or. &
                         &   solve(ivari) % kfl_algso == SOL_SOLVER_MAPHYS_SPD         ) then
                       solve(ivari) % lun_exsol = modul*lun_base + 40 + ivari ! 40 -> 49
                    end if

                 end if
              end do
              if ( size(solve,kind=ip) > 20 ) call runend('MODDEF: TOO MANY UNITS FOR SOLVER')
           end if

           if ( kfl_naked == 0 ) then
              !
              !  encapsulated, then get names from the environment
              !
              call GET_ENVIRONMENT_VARIABLE( trim(wfort(modul))//'06', postp(1) % fil_setse )
              call GET_ENVIRONMENT_VARIABLE( trim(wfort(modul))//'07', postp(1) % fil_setsb )
              call GET_ENVIRONMENT_VARIABLE( trim(wfort(modul))//'08', postp(1) % fil_setsn )
              call GET_ENVIRONMENT_VARIABLE( trim(wfort(modul))//'09', postp(1) % fil_setsi )
              call GET_ENVIRONMENT_VARIABLE( trim(wfort(modul))//'26', postp(1) % fil_witne )
              call GET_ENVIRONMENT_VARIABLE( trim(wfort(modul))//'27', postp(1) % fil_witng )

              do ivari = 1,nsolv
                 if ( solve(ivari) % kfl_algso /= -999 ) then
                    if ( solve(ivari) % kfl_solve /= 0 ) then
                       wtmp  = intost(solve(ivari)%lun_solve)
                       wmodu = wtmp(1:5)
                       call GET_ENVIRONMENT_VARIABLE( trim(wfort(modul))//trim(wmodu), solve(ivari) % fil_solve )
                    end if
                    if ( solve(ivari) % kfl_cvgso /= 0 ) then
                       wtmp  = intost(solve(ivari)%lun_cvgso)
                       wmodu = wtmp(1:5)
                       call GET_ENVIRONMENT_VARIABLE( trim(wfort(modul))//trim(wmodu), solve(ivari) % fil_cvgso )
                    end if
                    if(      solve(ivari) % kfl_algso == SOL_SOLVER_MAPHYS_UNSYMMETRIC .or. &
                         &   solve(ivari) % kfl_algso == SOL_SOLVER_MAPHYS_SYMMETRIC   .or. &
                         &   solve(ivari) % kfl_algso == SOL_SOLVER_MAPHYS_SPD         ) then
                       wtmp  = intost(solve(ivari)%lun_exsol)
                       wmodu = wtmp(1:5)
                       call GET_ENVIRONMENT_VARIABLE( trim(wfort(modul))//trim(wmodu), fil_exsol )
                    end if
                 end if
              end do

           else

              postp(1) % fil_setse = adjustl(trim(namda)) // '-element.'  // exmod(modul) // '.set'
              postp(1) % fil_setsb = adjustl(trim(namda)) // '-boundary.' // exmod(modul) // '.set'
              postp(1) % fil_setsn = adjustl(trim(namda)) // '-node.'     // exmod(modul) // '.set'
              postp(1) % fil_setsi = adjustl(trim(namda)) // '-IB.'       // exmod(modul) // '.set'
              postp(1) % fil_witne = adjustl(trim(namda)) // '.'          // exmod(modul) // '.wit'
              postp(1) % fil_witng = adjustl(trim(namda)) // '-geometry.' // exmod(modul) // '.wit'

              do ivari = 1,nsolv
                 if ( solve(ivari)%kfl_algso /= -999 ) then
                    if ( solve(ivari)%kfl_solve /= 0 ) then
                       solve(ivari) % fil_solve = adjustl(trim(namda))//'-'//trim(solve(ivari) % wprob)//'.'//exmod(modul)//'.sol'
                       if( INOTSLAVE ) &
                            call iofile(&
                            zero,solve(ivari) % lun_solve,solve(ivari) % fil_solve,namod(modul)//'  SOLVER', statu,forma,posit)
                    end if
                    if ( solve(ivari) % kfl_cvgso /= 0 ) then
                       solve(ivari) % fil_cvgso = adjustl(trim(namda))//'-'//trim(solve(ivari) % wprob)//'.'//exmod(modul)//'.cso'
                       if( INOTSLAVE ) &
                            call iofile(&
                            zero,solve(ivari) % lun_cvgso,solve(ivari) % fil_cvgso,namod(modul)//'  SOLVER', statu,forma,posit)
                    end if

                    if(      solve(ivari) % kfl_algso == SOL_SOLVER_MAPHYS_UNSYMMETRIC .or. &
                         &   solve(ivari) % kfl_algso == SOL_SOLVER_MAPHYS_SYMMETRIC   .or. &
                         &   solve(ivari) % kfl_algso == SOL_SOLVER_MAPHYS_SPD         ) then
                       fil_exsol = adjustl(trim(namda))//'-'//trim(solve(ivari) % wprob)//'-maphys.'//exmod(modul)//'.org'
                       if( PAR_MY_CODE_RANK_WM == 0 ) then
                          call iofile(&
                               zero,solve(ivari) % lun_exsol,trim(fil_exsol),namod(modul)//'  MAPHYS SOLVER', statu,forma,posit)
                       end if
                    end if
                 end if
              end do

           end if

           if( INOTSLAVE ) then
              !
              ! Element set file
              !
              if ( maxval(postp(1) % npp_setse) > 0 ) &
                   call iofile(zero,postp(1) % lun_setse,postp(1) % fil_setse,namod(modul)//' ELEMENT SETS ',statu,forma,posit)
              !
              ! Boundary set file
              !
              if ( maxval(postp(1) % npp_setsb) > 0 ) &
                   call iofile(zero,postp(1) % lun_setsb,postp(1) % fil_setsb,namod(modul)//' BOUNDARY SETS',statu,forma,posit)
              !
              ! Node set file
              !
              if ( maxval(postp(1) % npp_setsn) > 0) &
                   call iofile(zero,postp(1) % lun_setsn,postp(1) % fil_setsn,namod(modul)//' NODE SETS',    statu,forma,posit)
              !
              ! IB set file
              !
              if ( maxval(postp(1) % npp_setsb) > 0 .and. nboib > 0 ) &
                   call iofile(zero,postp(1) % lun_setsi,postp(1) % fil_setsi,namod(modul)//' IB SETS',statu,forma,posit)
              !
              ! Witness point
              !
              if ( maxval(postp(1) % npp_witne) > 0 .and. nwitn > 0 ) &
                   call iofile(zero,postp(1) % lun_witne,postp(1) % fil_witne,namod(modul)//' WITNESS POINT',statu,forma,posit)
              !
              ! Witness geometry
              !
              if ( maxval(postp(1) % npp_witng) > 0 .and. nwitg > 0 ) &
                   call iofile(zero,postp(1) % lun_witng,postp(1) % fil_witng,namod(modul)//' WITNESS GEOMETRY',statu,forma,posit)
           end if
        end if
     end do

     modul     =  0
     postp     => momod(modul) % postp
     solve     => momod(modul) % solve
     eigeg     => momod(modul) % eigen
     solve_sol => momod(modul) % solve
     veset     => postp(1)     % veset
     vbset     => postp(1)     % vbset
     vnset     => postp(1)     % vnset
     witne     => postp(1)     % witne
     witng     => postp(1)     % witng
     tncod     => momod(modul) % tncod
     tgcod     => momod(modul) % tgcod
     tbcod     => momod(modul) % tbcod

  else if ( jtask == 4_ip ) then

     !-------------------------------------------------------------------
     !
     ! Close result files
     !
     !-------------------------------------------------------------------

     if( INOTSLAVE ) then
        do modul = imod1,imod2
           if ( kfl_modul(modul) /= 0 ) then

              postp => momod(modul) % postp
              solve => momod(modul) % solve
              eigeg => momod(modul) % eigen

              if( associated(momod(modul) % solve) ) then
                 !
                 ! Close result files
                 !
                 call iofile(two,momod(modul) % lun_outpu,' ',namod(modul)//' OUTPUT')
                 call iofile(two,momod(modul) % lun_conve,' ',namod(modul)//' CONVERGENCE')
                 !if (modul==10_ip) then
                 !   call iofile(two,momod(modul) % lun_time,' ',namod(modul)//' TIME VALUES')
                 !end if
                 do ivari = 1,size(solve,kind=ip)
                    if ( solve(ivari) % kfl_algso /= -999 ) then
                       if ( solve(ivari) % kfl_solve /= 0 ) then
                          call iofile(two,solve(ivari) % lun_solve,solve(ivari) % fil_solve,namod(modul)//' SOLVER')
                       end if
                       if ( solve(ivari) % kfl_cvgso /= 0 ) then
                          call iofile(two,solve(ivari) % lun_cvgso,solve(ivari) % fil_cvgso,namod(modul)//' SOLVER CONVERGENCE')
                       end if
                       if(      solve(ivari) % kfl_algso == SOL_SOLVER_MAPHYS_UNSYMMETRIC .or. &
                            &   solve(ivari) % kfl_algso == SOL_SOLVER_MAPHYS_SYMMETRIC   .or. &
                            &   solve(ivari) % kfl_algso == SOL_SOLVER_MAPHYS_SPD         ) then
                          call iofile(two,solve(ivari) % lun_exsol,fil_exsol,namod(modul)//' MAPHYS SOLVER')
                       end if
                    end if
                 end do
              end if

              if ( maxval(postp(1) % npp_setse) > 0 ) &
                   call iofile(two,postp(1) % lun_setse,postp(1) % fil_setse,namod(modul)//' ELEMENT SETS')
              if ( maxval(postp(1) % npp_setsb) > 0 ) &
                   call iofile(two,postp(1) % lun_setsb,postp(1) % fil_setsb,namod(modul)//' BOUNDARY SETS')
              if ( maxval(postp(1) % npp_setsn) > 0) &
                   call iofile(two,postp(1) % lun_setsn,postp(1) % fil_setsn,namod(modul)//' NODE SETS')
              if ( maxval(postp(1) % npp_setsb) > 0 .and. nboib > 0 ) &
                   call iofile(two,postp(1) % lun_setsi,postp(1) % fil_setsi,namod(modul)//' IB SETS')
              if ( maxval(postp(1) % npp_witne) > 0 .and. nwitn > 0 ) &
                   call iofile(two,postp(1) % lun_witne,postp(1) % fil_witne,namod(modul)//' WITNESS POINT')
              if ( maxval(postp(1) % npp_witng) > 0 .and. nwitg > 0 ) &
                   call iofile(two,postp(1) % lun_witng,postp(1) % fil_witng,namod(modul)//' WITNESS GEOMETRY')

           end if
        end do
     end if

  else if ( jtask == 6_ip .and. INOTSLAVE ) then

     !-------------------------------------------------------------------
     !
     ! Close restart file
     !
     !-------------------------------------------------------------------

     if( kfl_modul(modul) /= 0 ) then
        if( modul == -1 ) then
           call iofile(two,lun_rstib,fil_rstib,'IB RESTART')
        else if( modul == 0 ) then
           call iofile(two,lun_rstar,fil_rstar,'RESTART')
        else
           call iofile(two,momod(modul) % lun_rstar,momod(modul) % fil_rstar,namod(modul)//' RESTART')
        end if
     end if

  else if ( jtask == 7_ip .and. INOTSLAVE ) then

     !-------------------------------------------------------------------
     !
     ! Open restart file for reading
     !
     !-------------------------------------------------------------------

     if( kfl_modul(modul) /= 0 ) then
        if( modul == -1 ) then
           call iofile_open_unit(lun_rstib,fil_rstib,'IB RESTART','old','unformatted')
        else if( modul == 0 ) then
           call iofile_open_unit(lun_rstar,fil_rstar,'RESTART', 'old','unformatted')
        else
           call iofile_open_unit(momod(modul) % lun_rstar,momod(modul) % fil_rstar,namod(modul)//' RESTART','old','unformatted',IOSTAT=ierror)
           if( ierror /= 0 ) call messages_live('COULD NOT OPEN RESTART FILE '//trim(momod(modul) % fil_rstar)//' OF '&
                //namod(modul)//'. SOME RESTART INFO MAY MISSING...','WARNING')
        end if
     end if

  else if ( jtask == 8_ip .and. INOTSLAVE ) then

     !-------------------------------------------------------------------
     !
     ! Open restart file for writing
     !
     !-------------------------------------------------------------------

     if( kfl_modul(modul) /= 0 ) then
        if( modul == -1 ) then
           fil_rstib_new = fil_rstib
           if( run_config%restart%append_timestep ) call appnum(ittim,fil_rstib_new)
           call iofile(zero,lun_rstib,fil_rstib_new,'IB RESTART','unknown','unformatted')
        else if( modul == 0 ) then
           fil_rstar_new = fil_rstar
           if( run_config%restart%append_timestep ) call appnum(ittim,fil_rstar_new)
           call iofile(zero,lun_rstar,fil_rstar_new,'RESTART', 'unknown','unformatted')
        else
           fil_rstar_new = momod(modul) % fil_rstar
           if( run_config%restart%append_timestep ) call appnum(ittim,fil_rstar_new)
           call iofile(zero,momod(modul) % lun_rstar,fil_rstar_new,namod(modul)//' RESTART','unknown','unformatted')
        end if
     end if

  else if ( jtask == 9_ip ) then

     !-------------------------------------------------------------------
     !
     ! Pointer to current module structures
     !
     !-------------------------------------------------------------------

     call moduls_set_current_module()

  end if

end subroutine moddef
