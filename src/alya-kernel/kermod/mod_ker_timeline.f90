!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!
!> @addtogroup Timeline_Toolbox
!> @{
!> @name    Timeline toolbox... figure out what Alya is doing graphically
!> @file    mod_ker_timeline.f90
!> @author  Guillaume Houzeaux
!> @brief   ToolBox for timeline
!> @details ToolBox for timeline
!
!-----------------------------------------------------------------------

module mod_ker_timeline
  
  use def_kintyp,         only : ip,rp
  use mod_communications, only : PAR_MIN
  use mod_communications, only : PAR_MAX
  use mod_communications, only : PAR_BARRIER
  use mod_iofile,         only : iofile
  use mod_iofile,         only : iofile_flush_unit
  use def_master
  use def_kermod
  use mod_run_config,     only : run_config
  implicit none
  
  real(rp) :: timeline_ini

  private 

  character(13) :: task_name(20) = (/ &
       'REAPRO', &
       'TURNON', &
       'INIUNK', &
       'TIMSTE', &
       'BEGSTE', &
       'DOITER', &
       'CONCOU', &
       'CONBLK', &
       '------', &
       'ENDSTE', &
       'FILTER', &
       'OUTPUT', &
       'TURNOF', &
       'BEGITE', &
       'ENDITE', &
       'MATRIX', &
       'DOOPTI', &
       'ENDOPT', &
       'BEGZON', &
       'ENDZON' /)     

  interface ker_timeline
     module procedure ker_timeline_module, &
          &           ker_timeline_kernel
  end interface ker_timeline
  
  public :: ker_timeline
  public :: ker_timeline_synchronization
  public :: ker_timeline_open_file

  public :: task_name
  
contains
  
  !-----------------------------------------------------------------------
  !
  !> @date    09/01/2017
  !> @author  Guillaume Houzeaux
  !> @brief   Synchronization
  !> @details Synchronization of all Alya codes
  !
  !-----------------------------------------------------------------------

  subroutine ker_timeline_synchronization()

    if( run_config%timeline ) then
       call PAR_BARRIER('IN THE WORLD')
       call cputim(timeline_ini)
    end if

  end subroutine ker_timeline_synchronization

  !-----------------------------------------------------------------------
  !
  !> @date    16/09/2016
  !> @author  Guillaume Houzeaux
  !> @brief   Write in file
  !> @details Output events from kernel in timeline file
  !
  !-----------------------------------------------------------------------

  subroutine ker_timeline_kernel(kmodu,event,intnumber,intcode)

    integer(ip),  intent(in)           :: kmodu
    character(*), intent(in)           :: event
    integer(ip),  intent(in), optional :: intnumber
    integer(ip),  intent(in), optional :: intcode
    
    call ker_timeline_all(kmodu,event,intnumber,intcode)

  end subroutine ker_timeline_kernel

  !-----------------------------------------------------------------------
  !
  !> @date    16/09/2016
  !> @author  Guillaume Houzeaux
  !> @brief   Write in file
  !> @details Output events from a module in timeline file
  !
  !-----------------------------------------------------------------------

  subroutine ker_timeline_module(event,intnumber)

    character(*), intent(in)           :: event
    integer(ip),  intent(in), optional :: intnumber
    
    call ker_timeline_all(modul,event,intnumber)

  end subroutine ker_timeline_module

  !-----------------------------------------------------------------------
  !>
  !> @date    16/09/2016
  !> @author  Guillaume Houzeaux
  !> @brief   Write in file
  !> @details Output events in timeline file
  !>          Universal codes to plot directly on gnuplot:
  !>           
  !> 
  !-----------------------------------------------------------------------

  subroutine ker_timeline_all(kmodu,event,intnumber,intcode)

    integer(ip),  intent(in)           :: kmodu
    character(*), intent(in)           :: event 
    integer(ip),  intent(in), optional :: intnumber
    integer(ip),  intent(in), optional :: intcode
    integer(ip),  save                 :: ipass=0
    real(rp)                           :: timeline
    integer(ip)                        :: inumb,icode
    integer(ip),  parameter            :: length_message=33
    character(length_message)          :: message

    if( run_config%timeline ) then

       call cputim(timeline)
       timeline = timeline-timeline_ini

       if(      trim(event(1:3)) == 'INI' ) then
          call PAR_MIN(timeline)
       else if( trim(event(1:3)) == 'END' ) then
          call PAR_MAX(timeline)
       end if

       if( INOTSLAVE ) then
          !
          ! Header
          !          
          if( ipass == 0 ) then
             write(lun_timeline,1) 
             ipass = 1
          end if
          !
          ! Input number
          !
          if( present(intnumber) ) then
             inumb = intnumber
          else
             inumb = -999
          end if
          !
          ! Automatic universal code
          !
          if(      trim(event) == 'INI_ASSEMBLY' ) then          ! icode = 1001,2001...
             icode = 10 * modul + 1
          else if( trim(event) == 'END_ASSEMBLY' ) then          ! icode = 1001,2001...
             icode = 10 * modul + 1
          else if( event(1:10) == 'INI_SOLVER'   ) then          ! icode = 1002,2002...             
             icode = 10 * modul + 2
             if( modul == ID_NASTIN ) then
                if( event(12:12) == 'C' ) icode = 10 * modul + 3
             end if
          else if( event(1:10) == 'END_SOLVER'   ) then          ! icode = 1002,2002...             
             icode = 10 * modul + 2
             if( modul == ID_NASTIN ) then
                if( event(12:12) == 'C' ) icode = 10 * modul + 3
             end if
          else if( trim(event) == 'INI_OUTPUT' ) then            ! icode = 1004,2004...
             icode = 10 * modul + 4
          else if( trim(event) == 'END_OUTPUT' ) then            ! icode = 1004,2004...
             icode = 10 * modul + 4
          else if( trim(event) == 'INI_READ_MESH' ) then         ! icode = 1
             icode = 1
          else if( trim(event) == 'END_READ_MESH' ) then
             icode = 1
          else if( trim(event) == 'INI_PARTITION_MESH' ) then    ! icode = 2
             icode = 2
          else if( trim(event) == 'END_PARTITION_MESH' ) then
             icode = 2
          else if( trim(event) == 'INI_COUPLING' ) then
             if( present(intcode) ) then
                icode = 10 * intcode + 1
             else
                call runend('WRONG CALL TO TIMELINE')
             end if
         else if( trim(event) == 'END_COUPLING' ) then
             if( present(intcode) ) then
                icode = 10 * intcode + 1
             else
                call runend('WRONG CALL TO TIMELINE')
             end if
          else
             icode = -999
          end if

          message(1:min(len(event,kind=ip),length_message)) = event(1:min(len(event,kind=ip),length_message))
          message(len(event)+1:length_message) = ' '
          write(lun_timeline,2) timeline,kmodu,adjustl(message),icode,inumb
          call iofile_flush_unit(lun_timeline)

       end if

    end if

1   format(&
         & '# Time        # module # event                           # code # numbers')

2   format(&
         & (e13.6),1x,(i8),1x,(a33),1x,(i6),1x,(i9))

  end subroutine ker_timeline_all

  !-----------------------------------------------------------------------
  !>
  !> @date    16/09/2016
  !> @author  Guillaume Houzeaux
  !> @brief   Open timeline file
  !> @details Open timeline file
  !> 
  !-----------------------------------------------------------------------

  subroutine ker_timeline_open_file()

    character(150) :: fil_timel

    if( run_config%timeline ) then
       
       if (INOTSLAVE) then

           if ( kfl_naked == 0 ) then
              call GET_ENVIRONMENT_VARIABLE('FOR042',fil_timel) 
           else if ( kfl_naked == 1 ) then
              fil_timel = adjustl(trim(namda))//'-timeline.res'     ! Unit 42
           end if
           
           if( kfl_rstar == 2 ) then
              call iofile(0_ip,lun_timeline,fil_timel,'TIMELINE','old','formatted','append')
           else
              call iofile(0_ip,lun_timeline,fil_timel,'TIMELINE')
           end if

       end if

    end if

  end subroutine ker_timeline_open_file

end module mod_ker_timeline
!> @}
