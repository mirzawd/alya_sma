!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @defgroup IO_Toolbox
!> Toolbox for IO, like open and close files
!> @{
!> @file    mod_iofile.f90
!> @author  houzeaux
!> @date    2018-04-11
!> @brief   Do things with file
!> @details Do some operations with units and files
!>
!-----------------------------------------------------------------------

module mod_iofile

  use def_kintyp_basic, only : ip,rp,lg
  use def_master,       only : file_opened
  use def_master,       only : kfl_reawr
  use def_master,       only : lun_outpu
  use mod_iofile_basic, only : iofile_available_unit   
  use mod_iofile_basic, only : iofile_flush_unit       
  use mod_iofile_basic, only : iofile_opened           
  use mod_iofile_basic, only : iofile_append_tag       
  use mod_iofile_basic, only : iofile_number_of_lines  
  use mod_iofile_basic, only : iofile_file_exists      
  use mod_iofile_basic, only : iofile_restart_run      
  use mod_iofile_basic, only : iofile_normal_run       
  use mod_iofile_basic, only : iofile_create_directory 
  use mod_iofile_basic, only : iofile_delete_file      
  implicit none
  private
  
  character(LEN=256) :: my_iomsg
  
  public :: iofile                  ! Do plenty of things with lots of arguments!
  public :: iofile_open_unit        ! Open a unit
  public :: iofile_close_unit       ! Close a unit
  public :: iofile_available_unit   ! Look for an available unit
  public :: iofile_flush_unit       ! Flush a unit
  public :: iofile_opened           ! Inquire if unit is opened
  public :: iofile_append_tag       ! Append a tag to a file (e.g. MPI rank)
  public :: iofile_number_of_lines  ! Computes the number of lines in a file until blank line is found
  public :: iofile_file_exists      ! If a file exists
  public :: iofile_restart_run      ! Define options for a restart run
  public :: iofile_normal_run       ! Define options for a normal run
  public :: iofile_create_directory ! Create a directory
  public :: iofile_delete_file      ! Delete a file if it exists

contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-04-16
  !> @brief   Open unit
  !> @details Open a file
  !> 
  !-----------------------------------------------------------------------

  subroutine iofile_open_unit(nunit,files,messa,stato,formo,posio,conver,crash_if_cannot_open,IOSTAT)

    integer(ip),      intent(in)              :: nunit 
    character(len=*), intent(in)              :: files
    character(len=*), intent(in),    optional :: messa
    character(len=*), intent(in),    optional :: stato
    character(len=*), intent(in),    optional :: formo
    character(len=*), intent(in),    optional :: posio
    character(len=*), intent(in),    optional :: conver
    logical(lg),      intent(in),    optional :: crash_if_cannot_open
    integer(ip),      intent(inout), optional :: IOSTAT
    logical(lg)                               :: crash

    if( present(IOSTAT) ) then
       crash = .false.
    else
       crash = .true.
    end if
    if( present(crash_if_cannot_open) ) crash = crash_if_cannot_open

    if( present(IOSTAT) ) then
       if( crash ) then
          call iofile(0_ip,nunit,files,messa,stato,formo,posio,conver,IOSTAT=IOSTAT)
       else
          call iofile(7_ip,nunit,files,messa,stato,formo,posio,conver,IOSTAT=IOSTAT)
       end if
    else
       if( crash ) then
          call iofile(0_ip,nunit,files,messa,stato,formo,posio,conver)
       else
          call iofile(7_ip,nunit,files,messa,stato,formo,posio,conver)
       end if
    end if

  end subroutine iofile_open_unit

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-04-11
  !> @brief   Close a file unit
  !> @details Close a file unit
  !> 
  !-----------------------------------------------------------------------

  subroutine iofile_close_unit(nunit,files,messa,IOSTAT,STATUS)

    integer(ip),      intent(in)              :: nunit 
    character*(*),    intent(in),    optional :: files
    character*(*),    intent(in),    optional :: messa
    integer(ip),      intent(inout), optional :: IOSTAT
    character(len=*), intent(in),    optional :: STATUS
    integer(4)                                :: iostat4
    integer(4)                                :: nunit4
    integer(4)                                :: lun_outpu4
    
    nunit4 = int(nunit,4)
    
    if(present(STATUS)) then
       close(UNIT=nunit4,status=trim(STATUS),IOSTAT=iostat4,ERR=100)
    else
       close(UNIT=nunit4,IOSTAT=iostat4,ERR=100)
    end if

    if( iostat4 /= 0 ) then
       if( present(messa) ) then
          if( trim(messa) == 'RESTART' ) then
             lun_outpu4 = int(lun_outpu,4)
             if( present(files) ) then
                write(lun_outpu4,101) 'COULD NOT OPEN THEN FOLLOWING RESTART FILE: '//trim(files)
             else
                write(lun_outpu4,101) 'COULD NOT OPEN THEN FOLLOWING RESTART FILE (UNKNOWN NAME)'
             end if
          else
             if( present(files) ) then
                call runend('ERROR WHEN CLOSING THE '//trim(messa)//' FILE: '//adjustl(trim(files)))
             else
                call runend('ERROR WHEN CLOSING FILE')
             end if
          end if
       else
          call runend('ERROR WHEN CLOSING FILE')          
       end if
    end if

    return
    
100 continue
    if( present(IOSTAT) ) then
       IOSTAT = int(iostat4,ip)
    else 
       !call runend('IOFILE_FLUSH_UNIT: COULD NOT CLOSE UNIT '//integer_to_string(nunit)//' WITH ERROR= '//integer_to_string(int(iostat4,ip)))
    end if
    
101 format(5x,'WARNING: ',a)
    
  end subroutine iofile_close_unit
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-04-11
  !> @brief   Main routine of this module
  !> @details Do all the operations with units and files
  !> 
  !-----------------------------------------------------------------------

  subroutine iofile(itask,nunit,files,messa,stato,formo,posio,conver,lgexi,IOSTAT)
    use def_master, only : intost
    implicit none
    integer(ip),      intent(in)            :: itask  !< Task
    integer(ip),      intent(in),  optional :: nunit  !< Unit
    character(len=*), intent(in),  optional :: files  !< File name
    character(len=*), intent(in),  optional :: messa  !< Message for warnings and errors
    character(len=*),              optional :: stato  !< Status
    character(len=*),              optional :: formo  !< Form
    character(len=*),              optional :: posio  !< Position
    character(len=*),              optional :: conver !< Convert
    logical(lg),                   optional :: lgexi
    integer(ip),      intent(out), optional :: IOSTAT
    integer(4)                              :: ioerr4
    integer(4)                              :: nunit4
    character(7)                            :: statu
    character(6)                            :: posit
    character(11)                           :: forma
    integer(4)                              :: lun_outpu4

    if( present(nunit) ) nunit4=int(nunit,4)

    select case (itask)

    case ( 0_ip )
       !
       ! Open unit
       !
       if( present(files) ) then
          if(present(stato)) then
             statu=stato
          else
             statu='unknown'
          end if
          if(present(formo)) then
             forma=formo
          else
             forma='formatted'
          end if
          if(present(posio)) then
             posit=posio
          else
             posit='asis'
          end if

#ifdef BIG_ENDIAN
          !
          ! Forces all to big endian
          !
          open(nunit4,file=adjustl(trim(files)),status=statu,form=forma,iostat=ioerr4,iomsg=my_iomsg,position=posit,convert='BIG_ENDIAN')
#else
          open(nunit4,file=adjustl(trim(files)),status=statu,form=forma,iostat=ioerr4,iomsg=my_iomsg,position=posit)
#endif   
          !
          ! Error when opening the file
          !
          if(ioerr4/=0) then
             file_opened = .false.
             if( present(messa) ) then
                call runend('ERROR WHEN OPENING THE '//trim(messa)//&
                     ' FILE: ' //adjustl(trim(files))//&
                     ', UNIT ' //trim(intost(nunit4))//&
                     ', ERROR '//trim(intost(ioerr4))//&
                     ', ERROR MESSAGE: '//trim(my_iomsg))
             else
                call runend('ERROR WHEN OPENING FILE: '//adjustl(trim(files)//&
                     ', UNIT ' //trim(intost(nunit4)))//&
                     ', ERROR '//trim(intost(ioerr4))//&
                     ', ERROR MESSAGE: '//trim(my_iomsg))
             end if
          else
             file_opened = .true.
          end if
       else
          call runend('IOFILE: FILE NAME NOT GIVEN!')
       end if

    case ( 2_ip )
       !
       ! Close unit
       !
       if(present(stato)) then
          statu=stato
          close(nunit4,status=statu,iostat=ioerr4)
       else
          close(nunit4,iostat=ioerr4)
       end if

       if(ioerr4/=0) then
          if( present(messa) .and. present(files) ) then
             if(trim(messa)=='RESTART') then
                lun_outpu4 = int(lun_outpu,4)
                write(lun_outpu4,101) 'COULD NOT OPEN THEN FOLLOWING RESTART FILE: '//trim(files)
             else
                call runend('ERROR WHEN CLOSING THE '//trim(messa)//' FILE: '//adjustl(trim(files)))
             end if
          else
             call runend('ERROR WHEN CLOSING A FILE')
          end if
       end if

    case ( 3_ip )
       !
       ! Delete file
       !
       close(nunit4,status='delete',iostat=ioerr4)
       if(ioerr4/=0) then
          if( present(messa) .and. present(files) ) then
             call runend('ERROR WHEN CLOSING THE '//trim(messa)//' FILE: '//adjustl(trim(files)))
          else
             call runend('ERROR WHEN CLOSING A FILE')
          end if
       end if

    case ( 4_ip )
       !
       ! Check if file exist
       !
       if(present(formo)) then
          forma=formo
       else
          forma='formatted'
       end if
       if(present(posio)) then
          posit=posio
       else
          posit='asis'
       end if
       if( present(files) ) & 
            open(nunit4,file=adjustl(trim(files)),status='old',form=forma,iostat=ioerr4,position=posit)

       if( ioerr4 /= 0 ) then
          file_opened = .false.
          kfl_reawr   = -abs(kfl_reawr)
       else
          file_opened = .true.
          close(nunit4)
       end if

    case ( 7_ip )
       !
       ! Open unit but do not crash if cannot
       !
       if(present(stato)) then
          statu=stato
       else
          statu='unknown'
       end if
       if(present(formo)) then
          forma=formo
       else
          forma='formatted'
       end if
       if(present(posio)) then
          posit=posio
       else
          posit='asis'
       end if

       if( present(files) ) then
#ifdef BIG_ENDIAN
          ! forces all to big endian
          open(nunit4,file=adjustl(trim(files)),status=statu,form=forma,iostat=ioerr4,position=posit,convert='BIG_ENDIAN')
#else
          if( present(conver) ) then
             open(nunit4,file=adjustl(trim(files)),status=statu,form=forma,iostat=ioerr4,position=posit)!,convert=conver)
          else
             open(nunit4,file=adjustl(trim(files)),status=statu,form=forma,iostat=ioerr4,position=posit)
          end if
#endif   
          if(ioerr4/=0) then
             file_opened = .false.
          else
             file_opened = .true.
          end if
       end if

    end select
    !
    ! Format
    !
    if( present(IOSTAT) ) IOSTAT = int(ioerr4,ip)

101 format(5x,'WARNING: ',a)

  end subroutine iofile
  
end module mod_iofile
!> @}
