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
!> @file    mod_iofile_basic.f90
!> @author  houzeaux
!> @date    2018-04-11
!> @brief   Do things with file
!> @details Do some operations with units and files
!>
!-----------------------------------------------------------------------

module mod_iofile_basic
  use, intrinsic :: iso_c_binding, only: c_int, c_char, c_int16_t, c_null_char
  use def_kintyp_basic, only : ip,rp,lg
  use mod_strings,      only : integer_to_string

  implicit none
  private
  
  interface iofile_flush_unit
     module procedure iofile_flush_unit_4,&
          &           iofile_flush_unit_8
  end interface iofile_flush_unit


  
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


  interface
      !  should be unsigned int ... not available in Fortran
      !  OK until highest bit gets set.
      function FortSleep (seconds)  bind ( C, name="sleep" )
            import
            integer (c_int) :: FortSleep
            integer (c_int), VALUE :: seconds
      end function FortSleep


      function FortMkdir(path,mode) bind(c,name="mkdir")
         import
         integer(c_int) :: FortMkdir
         character(kind=c_char,len=1) :: path(*)
         integer(c_int16_t), value :: mode
      end function FortMkdir

      !function FortChmod(path,mode) bind(c,name="chmod")
      !   import
      !   integer(c_int) :: FortChmod
      !   character(kind=c_char,len=1) :: path(*)
      !   integer(c_int16_t), value :: mode
      !end function FortChmod


  end interface

contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-04-16
  !> @brief   Delete a file
  !> @details Delete a file if it exists
  !> 
  !-----------------------------------------------------------------------

  subroutine iofile_delete_file(wfile,nunit)

    character(len=*),           intent(in) :: wfile
    integer(ip),      optional, intent(in) :: nunit
    integer(4)                             :: nunit4
    
    if( iofile_file_exists(wfile) ) then
       if( present(nunit) ) then
          nunit4 = int(nunit,4)
       else
          nunit4 = int(iofile_available_unit(10000_ip),4)
       end if
       open(UNIT=nunit4,file=trim(wfile),status='old')
       close(UNIT=nunit4,status='DELETE')
    end if
    
  end subroutine iofile_delete_file
 
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-04-16
  !> @brief   Inquire if file exists
  !> @details Inquire if file exists
  !> 
  !-----------------------------------------------------------------------

  logical(lg) function iofile_file_exists(wfile)

    character(len=*), intent(in) :: wfile

    inquire(FILE=trim(wfile),EXIST=iofile_file_exists)

  end function iofile_file_exists

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-04-16
  !> @brief   Inquire if unit is opened
  !> @details Inquire if unit is opened
  !> 
  !-----------------------------------------------------------------------

  logical(lg) function iofile_opened(lunit)

    integer(ip), intent(in) :: lunit
    integer(4)              :: unit4

    unit4=int(lunit,4)
    inquire(unit=unit4,opened=iofile_opened)

  end function iofile_opened

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-04-16
  !> @brief   Flush a unit
  !> @details Flush a unit
  !> 
  !-----------------------------------------------------------------------

  subroutine iofile_flush_unit_ip(nunit)

    integer(ip), intent(in) :: nunit
    integer(4)              :: iostat4,nunit4
    integer(ip)             :: iostat

    nunit4 = int(nunit,4)
    flush(UNIT=nunit4,IOSTAT=iostat4,ERR=100)
    return

    iostat = int(iostat4,ip)
100 call runend('IOFILE_FLUSH_UNIT: COULD NOT FLUSH UNIT '//integer_to_string(nunit)//' WITH ERROR= '//integer_to_string(iostat))
    
  end subroutine iofile_flush_unit_ip
  
  subroutine iofile_flush_unit_4(nunit)

    integer(4),  intent(in) :: nunit
    integer(ip)             :: nunit_ip

    nunit_ip = int(nunit,ip)
    call iofile_flush_unit_ip(nunit_ip)

  end subroutine iofile_flush_unit_4
 
  subroutine iofile_flush_unit_8(nunit)

    integer(8),  intent(in) :: nunit
    integer(ip)             :: nunit_ip

    nunit_ip = int(nunit,ip)
    call iofile_flush_unit_ip(nunit_ip)

  end subroutine iofile_flush_unit_8
 
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-04-16
  !> @brief   Available unit
  !> @details Look for an available unit
  !> 
  !-----------------------------------------------------------------------

  integer(ip) function iofile_available_unit(START_AT)

    integer(ip), optional, intent(in) :: START_AT
    integer(4)                        :: ioerr,unit4,munit4
    integer(4)                        :: unit4_start
    logical(lg)                       :: opened

    munit4 = 10000_4

    if( present(START_AT) ) then
       unit4_start = int(START_AT,4)
    else
       unit4_start = 90_4
    end if
    
    do unit4 = unit4_start,munit4
       inquire(unit=unit4,opened=opened,iostat=ioerr)
       if( ioerr /= 0 )  cycle
       if( .not. opened ) exit
    end do

    if( unit4 > munit4 ) then
       iofile_available_unit = 0_ip
    else
       iofile_available_unit = int(unit4,ip)
    end if

  end function iofile_available_unit

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-04-16
  !> @brief   Append a tag to a file
  !> @details Append a tag to a file (e.g. MPI rank)
  !> 
  !-----------------------------------------------------------------------

  subroutine iofile_append_tag(filen,tag,fileo)

    character(len=*), intent(inout)          :: filen
    integer(ip),      intent(in)             :: tag
    character(len=*), intent(out),  optional :: fileo
    integer(ip)                              :: fdot,flen
    
    fdot = index(filen,'.')
    flen = len_trim(filen)
    if( fdot == 0 ) fdot = flen
    if( present(fileo) ) then
       fileo = trim(filen(1:fdot-1)//'-'//trim(integer_to_string(tag))//filen(fdot:flen))
    else
       filen = trim(filen(1:fdot-1)//'-'//trim(integer_to_string(tag))//filen(fdot:flen))
    end if
    
  end subroutine iofile_append_tag

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-04-16
  !> @brief   Number of lines
  !> @details Gives the number of lines of a file
  !> 
  !-----------------------------------------------------------------------

  function iofile_number_of_lines(nunit) result(nlines)

    integer(ip), intent(in)  :: nunit
    integer(4)               :: nunit4
    integer(4)               :: ierro
    integer(ip)              :: nlines
    character(10)            :: str
    
    nunit4 = int(nunit,4)
    nlines = 0
    rewind(nunit)
    do
       read(nunit,FMT='(A)',iostat=ierro,end=10) str
       if( len_trim(str) == 0 ) goto 10
       if( ierro/= 0 ) exit
       nlines = nlines + 1
    end do
    
10  rewind(nunit)

  end function iofile_number_of_lines
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-04-16
  !> @brief   File option for a restart run
  !> @details File option for a restart run
  !> 
  !-----------------------------------------------------------------------

  subroutine iofile_restart_run(stato,formo,posio,FORMATTED,UNFORMATTED)

    character(len=*), intent(out)          :: stato
    character(len=*), intent(out)          :: formo
    character(len=*), intent(out)          :: posio
    logical(lg),      intent(in), optional :: FORMATTED
    logical(lg),      intent(in), optional :: UNFORMATTED
    logical(lg)                            :: if_formatted

    if_formatted = .true.
    if( present(FORMATTED) )   if_formatted = FORMATTED
    if( present(UNFORMATTED) ) if_formatted = .not. UNFORMATTED

    stato = 'old'
    posio = 'append'

    if( if_formatted ) then
       formo = 'formatted'
    else
       formo = 'unformatted'
    end if
   
  end subroutine iofile_restart_run
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-04-16
  !> @brief   File option for a normal run
  !> @details File option for a normal run
  !> 
  !-----------------------------------------------------------------------

  subroutine iofile_normal_run(stato,formo,posio,FORMATTED,UNFORMATTED)

    character(len=*), intent(out)          :: stato
    character(len=*), intent(out)          :: formo
    character(len=*), intent(out)          :: posio
    logical(lg),      intent(in), optional :: FORMATTED
    logical(lg),      intent(in), optional :: UNFORMATTED
    logical(lg)                            :: if_formatted

    if_formatted = .true.
    if( present(FORMATTED) )   if_formatted = FORMATTED
    if( present(UNFORMATTED) ) if_formatted = .not. UNFORMATTED

    stato = 'unknown'
    posio = 'asis'
    
    if( if_formatted ) then
       formo = 'formatted'
    else
       formo = 'unformatted'
    end if

  end subroutine iofile_normal_run


  !-----------------------------------------------------------------------
  !> 
  !> @author  cbutakoff
  !> @date    2021-11-29
  !> @brief   find the location of the next "/" in the path
  !> @details find the location of the next "/" in the path
  !> 
  !-----------------------------------------------------------------------
  integer(4) function find_next_dir_delimiter( str2, pos )
      implicit none
      character(len=*)  , intent(in)  :: str2            
      integer(4)        , intent(in)  :: pos
      integer(4)                      :: i
      character(len=1)                :: delimiter

         
      delimiter = '/'
#ifdef _WIN32
      delimiter = '\\'
#endif

      find_next_dir_delimiter = -1
      i = pos-1
      do while ( i<len(str2) )
         i = i+1
         if (str2(i:i)==delimiter) then
               find_next_dir_delimiter = i
               i = len(str2) + 1 
         end if
      end do

   end function find_next_dir_delimiter


   !-----------------------------------------------------------------------
   !> 
   !> @author  cbutakoff
   !> @date    2021-11-29
   !> @brief   make directory tree akin to mkdir -p BOU1/test/x
   !> @details make directory tree akin to mkdir -p BOU1/test/x
   !> 
   !-----------------------------------------------------------------------
   function iofile_maketree(path)
      implicit none
      integer(4)                   :: iofile_maketree
      character(len=*), intent(in) :: path
      integer(4)                   :: pos, pos_next
      integer (c_int)              :: mkdir_status
      integer (c_int16_t)          :: folder_mode

      folder_mode = int(509, c_int16_t) !chmod 775 (decimal 509)
               
      iofile_maketree = 1
      pos = 1
      do while ( pos > 0 )
         pos_next = find_next_dir_delimiter(path,pos+1)    
         if (pos_next>0) then
            mkdir_status = FortMkdir( path(1:pos_next-1) // c_null_char, mode=folder_mode)
            iofile_maketree = iofile_maketree * mkdir_status
         end if

         pos = pos_next
      end do

      mkdir_status = FortMkdir( trim(path) // c_null_char, mode=folder_mode) 
      iofile_maketree = iofile_maketree * mkdir_status

   end function iofile_maketree




  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-04-16
  !> @brief   Create a directory
  !> @details Create a directory
  !> 
  !-----------------------------------------------------------------------

  subroutine iofile_create_directory(dirname)
    use def_master,                  only : ip, intost
    use, intrinsic :: iso_c_binding, only : c_int, c_null_char

    character(len=*), intent(in) :: dirname
    !integer(4)                   :: MY_EXITSTAT
    !integer(4)                   :: MY_CMDSTAT
    !character(len=256)           :: MY_MESSAGE
    !integer(4)                   :: attempt_count
    integer (c_int)              :: status   

    !MY_EXITSTAT = 1
    !attempt_count        = 0

    !call execute_command_line ('mkdir -p '//trim(dirname),wait=.true.,EXITSTAT=MY_EXITSTAT,CMDSTAT=MY_CMDSTAT, CMDMSG=MY_MESSAGE)
    
    !do while( MY_EXITSTAT /= 0 .and. attempt_count < 10 ) 
    !   attempt_count = attempt_count + 1
       !with EXTRAE, mkdir seems to be executed asynchronously and EXITSTAT is always==1, therefore need this additional test
    !   status = FortSleep(1)
    !   if ( direxist(trim(dirname)) ) MY_EXITSTAT = 0
    !end do

    !if( MY_EXITSTAT /= 0 ) then !if system call failed, try POSIX 
    !  status = iofile_maketree( trim(dirname) )
    !   if ( direxist(trim(dirname)) ) MY_EXITSTAT = 0
    !end if

    status = iofile_maketree( trim(dirname) )
    if ( .not. direxist(trim(dirname)) ) call runend('Failed to create folder '//trim(dirname))


    !if ( MY_EXITSTAT /= 0 ) then
    !   call runend('Failed to create folder '//trim(dirname)//&
    !               ', EXITSTAT='//trim(intost(int(MY_EXITSTAT,ip)))//&
    !               ', CMDSTAT='//trim(intost(int(MY_CMDSTAT,ip)))//&
    !               ', CMDMSG='//trim(MY_MESSAGE) )
    !end if

    contains

         logical(lg) function direxist(dir)
            implicit none
            character(len=*), intent(in) :: dir
            integer(ip)                  :: unitno

            ! Test whether the directory exists
            
            unitno = iofile_available_unit()
            open(newunit=unitno,file='./'//trim(dir)//'/deleteme.txt',status='unknown',err=1234)
            !do not delete the file, if several alyas are running, and one alya deletes the file before the other - everything dies
            !TODO: make a better direxist() or call something from C/C++
            close(unitno) 
            direxist = .true.
            return

            ! If doesn't exist, end gracefully
1234        direxist = .false.
            return

         end function

  end subroutine iofile_create_directory
  
end module mod_iofile_basic
!> @}
