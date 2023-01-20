!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Output
!> @{
!> @file    mod_status.f90
!> @author  houzeaux
!> @date    2020-09-11
!> @brief   Alya state
!> @details Output Alya state in file
!-----------------------------------------------------------------------

module mod_state

  use def_kintyp_basic, only : ip,rp
  use mod_strings,      only : string_to_integer
  use mod_strings,      only : integer_to_string
  use def_master,       only : namda
  use def_master,       only : ittim
  use def_master,       only : mitim
  use def_master,       only : cutim
  use def_master,       only : timef
  use def_master,       only : dtime
  use def_master,       only : cpu_times
  use def_master,       only : cpu_initi
  use def_master,       only : ittim_ini
  use def_master,       only : ittim_rst
  use def_master,       only : lun_state
  use def_master,       only : INOTSLAVE
  use mod_io_config,    only : io_config
  use mod_iofile,       only : iofile_open_unit
  use mod_iofile,       only : iofile_close_unit

  implicit none
  private

  integer                  :: cmd_s
  character(200)           :: cmd_m
  character(200)           :: job_name
  character(200)           :: job_clustername
  integer(ip)              :: job_id
  integer(ip)              :: num_times
  integer(ip),   parameter :: ave_times=3
  
  public :: state_output

contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-09-11
  !> @brief   JOBID
  !> @details Get the JOBID
  !> 
  !-----------------------------------------------------------------------

  function jobinfo_id() result(id) 

    integer(ip)   :: id
    character(20) :: wid
    integer(ip)   :: stat

    id = 0
    call get_environment_variable("SLURM_JOB_ID",wid)
    id = string_to_integer(wid,stat) 
    if( stat /= 0 ) id = 0

  end function jobinfo_id

  subroutine jobinfo_name(name) 

    character(*) :: name

    !call execute_command_line("squeue -o %j > foo",WAIT=.true.,CMDSTAT=cmd_s,CMDMSG=cmd_m)
    
    call execute_command_line("bash -c ''squeue -o %j > foo 2> foo.err''",WAIT=.true.,CMDSTAT=cmd_s,CMDMSG=cmd_m)

    if( cmd_s == 0 ) then
       call execute_command_line("squeue -o %j > foo",CMDSTAT=cmd_s)
       open(unit=99,file='foo',status='unknown')  
       read(99,'(a)',end=999,err=999) name
       read(99,'(a)',end=999,err=999) name
       close(99)
       name = trim(name)
       return
    else
       call execute_command_line("rm -rf foo; rm -rf foo.err",CMDSTAT=cmd_s)       
    end if

999 name = 'null'

  end subroutine jobinfo_name

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-09-11
  !> @brief   Cluster name
  !> @details Name of the cluster where the job is being executed
  !> 
  !-----------------------------------------------------------------------

  subroutine jobinfo_clustername(name)

    character(*) :: name

    !call execute_command_line("squeue -O cluster > foo",WAIT=.true.,CMDSTAT=cmd_s,CMDMSG=cmd_m)
    call execute_command_line("bash -c ''squeue -O cluster > foo 2> foo.err''",WAIT=.true.,CMDSTAT=cmd_s,CMDMSG=cmd_m)

    if( cmd_s == 0 ) then
       open(unit=99,file='foo',status='unknown')  
       read(99,'(a)',end=999,err=999) name
       read(99,'(a)',end=999,err=999) name
       close(99)
       name = trim(name)
       return
    else
       call execute_command_line("rm -rf foo; rm -rf foo.err",CMDSTAT=cmd_s)       
    end if

999 name = ' '

  end subroutine jobinfo_clustername

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-09-11
  !> @brief   Cluster name
  !> @details Name of the cluster where the job is being executed
  !> 
  !-----------------------------------------------------------------------

  subroutine jobinfo_timeleft(name)

    character(*) :: name

    !call execute_command_line("squeue -o %L > foo",WAIT=.true.,CMDSTAT=cmd_s,CMDMSG=cmd_m)
    call execute_command_line("bash -c ''squeue -o %L > foo 2> foo.err''",WAIT=.true.,CMDSTAT=cmd_s,CMDMSG=cmd_m)

    if( cmd_s == 0 ) then
       open(unit=99,file='foo',status='unknown')  
       read(99,'(a)',end=999,err=999) name
       read(99,'(a)',end=999,err=999) name
       close(99)
       name = trim(name)
       return
    end if

999 name = ' '

  end subroutine jobinfo_timeleft

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-09-11
  !> @brief   Cluster name
  !> @details Name of the cluster where the job is being executed
  !> 
  !-----------------------------------------------------------------------

  subroutine jobinfo_state(name)

    character(*) :: name

    !call execute_command_line("squeue -o %T > foo",WAIT=.true.,CMDSTAT=cmd_s,CMDMSG=cmd_m)
    call execute_command_line("bash -c ''squeue -o %T > foo 2> foo.err''",WAIT=.true.,CMDSTAT=cmd_s,CMDMSG=cmd_m)

    if( cmd_s == 0 ) then
       open(unit=99,file='foo',status='unknown')  
       read(99,'(a)',end=999,err=999) name
       read(99,'(a)',end=999,err=999) name
       close(99)
       name = trim(name)
       return
    end if

999 name = ' '

  end subroutine jobinfo_state

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-09-11
  !> @brief   State
  !> @details Output the state of Alya
  !> 
  !-----------------------------------------------------------------------

  subroutine state_output(MY_STATE)

    character(*),     optional, intent(in) :: MY_STATE
    integer                                :: values(8)
    character(2)                           :: wh
    character(2)                           :: wm
    character(2)                           :: ws
    real(rp)                               :: cpu_refer
    real(rp)                               :: cpu_estim
    character(200)                         :: job_state
    character(200)                         :: job_date
    character(200)                         :: job_time
    integer(ip)                            :: num_steps_left
    integer(ip)                            :: num_steps_done
    integer(ip),      save                 :: ipass=0

    if( INOTSLAVE .and. io_config%state ) then
       !
       ! JOB_ID: Constant info
       !
       if( ipass == 0 ) then
          ipass = 1
          call jobinfo_name       (job_name       )
          call jobinfo_clustername(job_clustername)
          job_id = jobinfo_id()
          num_times = 0          
       end if
       !
       ! JOB_DATE= Job date
       !
       call DATE_AND_TIME(VALUES=values)
       job_date = &
            &      integer_to_string(int(values(3),ip))//'/'// &
            &      integer_to_string(int(values(2),ip))//'/'// &
            &      integer_to_string(int(values(1),ip))
       !
       ! JOB_TIME= Job time
       !
       if( values(5) < 10 ) then
          wh = '0'//integer_to_string(int(values(5),ip))
       else
          wh = integer_to_string(int(values(5),ip))
       end if
       if( values(6) < 10 ) then
          wm = '0'//integer_to_string(int(values(6),ip))
       else
          wm = integer_to_string(int(values(6),ip))
       end if
       if( values(7) < 10 ) then
          ws = '0'//integer_to_string(int(values(7),ip))
       else
          ws = integer_to_string(int(values(7),ip))
       end if
       job_time = &
            wh//':'//wm//':'//ws
       !
       ! CPU_REFER = CPU per time step
       !
       call cputim(cpu_refer)
       if( ittim-ittim_ini /= 0 ) then
          cpu_refer = (cpu_refer-cpu_times) / real(ittim-ittim_ini,rp)
       else
          cpu_refer = 0.0_rp
       end if
       !
       ! NUM_STEPS_DONE: Number of time steps done
       ! NUM_STEPS_LEFT: Number of remaining time steps and ETA
       !
       num_steps_done = ittim-ittim_ini
       num_steps_left = min( int((timef-cutim)/dtime,ip) , mitim-ittim )
       cpu_estim      = cpu_refer * real(num_steps_left,rp)
       !
       ! Current state
       !
       if( present(MY_STATE) ) then
          job_state = trim(MY_STATE)
       else
          job_state = 'RUNNING'
       end if

       !call jobinfo_state   (job_state      )
       !call jobinfo_timeleft(job_timeleft   )

       call iofile_open_unit(lun_state,'state.log','STATE FILE')

       write(lun_state,'(a)') &
            &     'name, job name, date, time, job id, job state, time steps done, '// &
            &     ' time steps remaining, estimated time left'
       write(lun_state,1)         &
            trim(namda),          & ! 1
            trim(job_name),       & ! 2
            trim(job_date),       & ! 3
            trim(job_time),       & ! 4
            job_id,               & ! 5
            trim(job_state),      & ! 6
            num_steps_done,       & ! 7
            num_steps_left,       & ! 8
            cpu_estim               ! 9       

       call iofile_close_unit(lun_state,'state.log','STATE FILE')

    end if

1   format(a, ', ', a, ', ', a, ', ', a, ', ' , i6, ', ', a, ', '  , i6, ', ', i6, ',', e15.8e3)

  end subroutine state_output

end module mod_state
