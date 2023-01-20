!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Kernel
!> @{
!> @file    mod_ecoute.f90
!> @author  houzeaux
!> @date    2022-12-19
!> @brief   Sparse a file
!> @details Sparse a file and interpret words and parameters
!-----------------------------------------------------------------------

module mod_ecoute

  use def_kintyp,  only : ip,rp
  use mod_strings, only : upper_case
  use mod_strings, only : integer_to_string
  use def_inpout

  implicit none

  private
  
  abstract interface
     subroutine message_generic(messa)
       character(LEN=*), intent(in) :: messa
     end subroutine message_generic
  end interface
  
  procedure(message_generic), pointer :: message

  public :: ecoute
  public :: ecoute_reach_section
  public :: ecoute_initialization
  public :: ecoute_set_read_unit
  public :: ecoute_set_write_unit
  
contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-01-24
  !> @brief   Set the unit to read from
  !> @details Set the unit to read from
  !> 
  !-----------------------------------------------------------------------

  subroutine ecoute_initialization(messa)
     
    procedure(message_generic) :: messa

    message => messa

  end subroutine ecoute_initialization
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-01-24
  !> @brief   Set the unit to read from
  !> @details Set the unit to read from
  !> 
  !-----------------------------------------------------------------------

  subroutine ecoute_set_read_unit(current_unit)

    integer(ip), intent(in) :: current_unit

    call set_read_unit(current_unit)
    
  end subroutine ecoute_set_read_unit

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-01-24
  !> @brief   Set the unit to write
  !> @details Set the unit to write
  !> 
  !-----------------------------------------------------------------------

  subroutine ecoute_set_write_unit(current_unit)

    integer(ip), intent(in) :: current_unit

    call set_write_unit(current_unit)
    
  end subroutine ecoute_set_write_unit

  !-----------------------------------------------------------------------
  !
  !> @author  houzeaux
  !> @date    2022-01-24
  !> @brief   Interpret a line
  !> @details Reads a string and interprets it as words and parameters.
  !>
  !>          - Maximum number of words and parameters = maxwp.
  !>          - Only the first five characters of each word are decoded.
  !>          - The underline characters '_' are discarted.
  !>          - Lower case letters are converted to upper case.
  !>          - Each word or parameter must be separated by ' ', '=', ':' or ','
  !>          - A "comment" begins with '$', '!', '/'  or '/' in any place.
  !>          - A line that has a comment beginning with '/' or '/'
  !>            continuates in the next line.
  !>          - A line beginning with title is not decoded. title directive.
  !>          - A line beginning with include is an include directive.
  !>          - A line beginning with echo turns 'on' or 'off' the echo.
  !
  !-----------------------------------------------------------------------

  subroutine ecoute(subna,STOP_END_OF_FILE, DO_NOT_READ_INCLUDE)

    character(*), intent(in)           :: subna
    logical(lg),  intent(in), optional :: STOP_END_OF_FILE
    logical(lg),  intent(in), optional :: DO_NOT_READ_INCLUDE
    real(rp)                           :: digit
    integer(ip)                        :: first,firsp,i,last,lastp,ptrwo,npptr,nwptr
    integer(ip)                        :: leng,flag,resum
    logical(lg)                        :: newline=.false.
    logical(lg), save                  :: echo=.false.     ! default echo off. to change it use: echo on
    logical(lg)                        :: stop_end_of_file_opt
    logical(lg)                        :: do_not_read_include_opt
    !
    ! Options
    !
    stop_end_of_file_opt = .true.
    if( present(STOP_END_OF_FILE) ) then
       stop_end_of_file_opt = STOP_END_OF_FILE
    end if
    do_not_read_include_opt = .false.
    if( present(DO_NOT_READ_INCLUDE) ) then
       do_not_read_include_opt = DO_NOT_READ_INCLUDE
    end if
    !
    ! Set unit
    !
    if(lispa==0) then
       call set_unit(lisda,RESET=.true.)               !  initial data file
    end if
    !
    ! Begin
    !
    ccard = ' '
    nnwor = 0                                           ! initialize.
    nnpar = 0
    nwptr = 0
    npptr = 0
    resum = 0
    do i = 1,maxwp
       words(i) = ' '
       param(i) = 0.0_rp
    end do
    !
    ! Binary ecoute reading
    !
100 continue
    do while(((nnwor==0).and.(nnpar==0)&              ! don't return without answer
         .or.newline).or.resum==1)                    ! continue reading if / or \

       if (resum==0) then
          newline=.false.                             ! initialize.
          last=0
          lastp=0
       end if
       firsp=1
       resum=0
       read(nunit,10,end=101,err=1) ccard             ! read a card
       if( subna(1:4) == 'DONT' ) return
       !     leng=lnbln1(ccard)                             ! calculate the length.
       leng=len_trim(ccard)                             ! calculate the length.

       decode_card: do while(last<leng)               ! decode all the card.
          first=last+1
          loop_first: do while(            &
               ccard(first:first)=='_'.or. &          ! jump null character (_)
               ccard(first:first)==' '.or. &          ! jump separators ( =:,)
               ccard(first:first)=='='.or. &
               ccard(first:first)==':'.or. &
               ccard(first:first)==','.or. &
               iachar(ccard(first:first)) == 9 )         ! Tab character ASCII code is 9
             first=first+1
             if(first>leng) exit loop_first
          end do loop_first
          if(last==0) firsp=first                     ! save first to print card
          last=first
          loop_last: do while(             &
               ccard(last:last)/=' ' .and. &          ! look for separator ( =:,).
               ccard(last:last)/='=' .and. &
               ccard(last:last)/=':' .and. &
               ccard(last:last)/=',' .and. &
               ccard(last:last)/='$' .and. &          ! look for coment ($!/\).
               ccard(last:last)/='!' .and. &
               ccard(last:last)/='\' .and. &
               ccard(last:last)/='\' .and. &
               iachar(ccard(last:last)) /= 9 )        ! Tab character ASCII code is 9     )
             last=last+1
             if(last>leng) exit loop_last
          end do loop_last

          if(last<=251            .and.(&
               ccard(last:last)=='$'.or.&
               ccard(last:last)=='!')) leng=last-1
          if(last<=251            .and.(&
               ccard(last:last)=='\'.or.&
               ccard(last:last)=='\')) then
             leng=last-1                              ! deal with continuation
             newline=.true.                           ! set new line flag.
          end if
          last=last-1
          if(last>=first) then                        ! is it a word or a parameter
             lastp=last                               ! save last to print cardlogic
             call decod1(last-first+1_ip,&
                  ccard(first:last),flag,digit)
             wname= ccard(first:last)                 ! keep the useful ccard part in wname
             if(flag==0) then                         ! it is a parameter.
                nnpar=nnpar+1                         ! # of parameters
                npptr=npptr+1                         ! integer :: to next parameter
                if(npptr>maxwp) go to 4               ! error.
                if(nwptr>npptr) npptr=nwptr
                param(npptr)=digit
             else                                     ! it is a word.
                nnwor=nnwor+1                         ! # of words
                nwptr=nwptr+1                         ! integer :: to next word
                if(nwptr>maxwp) go to 5               ! error.
                if(npptr>=nwptr) nwptr=npptr+1
                ptrwo=1
                do while ((first<=last).and.(ptrwo<=5))
                   words(nwptr)(ptrwo:ptrwo)=ccard(first:first)
                   ptrwo=ptrwo+1
                   first=first+1
                   do while (ccard(first:first)=='_') ! jump null character
                      first=first+1
                   end do
                end do
                words(nwptr) = upper_case(words(nwptr)) ! convert to upper case.
             end if
          end if ! (last>=first)
          if((nnwor==1).and.(nnpar==0).and.&          ! deal with title or include
               ((words(1)=='TITLE').or.((words(1)=='INCLU').and.&
               (subna/='NOREAD')))) then
             if(echo.and.(subna/='noecho'))&
                  write(lisre,20) 'ecoute',ccard(firsp:leng)
             last=last+2
             do while(ccard(last:last)==' ')          ! remove blank spaces
                last=last+1
             end do
             if(leng<last) go to 6                    ! error
             ccard=ccard(last:leng)                   ! remove words(1) from ccard
             leng=leng-last+1
             if(words(1)=='TITLE') then               ! deal with titles
                firsp=1
                lastp=leng
             else if (.not.do_not_read_include_opt) then     ! deal with include directive
                !if(nunit==lisin) go to 3             ! error
                last=1                                ! remove tail comments
                do while((last<=leng).and.&           ! look for end (last=leng)
                     ccard(last:last)/=' ')           ! look for separator ( )
                   last=last+1
                end do
                ccard=ccard(1:last-1)                 ! remove tail comments
                if(nunit==lisin) then
                   call set_unit(lisi1)               !  initial data file
                else
                   call set_unit(lisin)               !  initial data file
                end if
                call opincl()
                lastp=0                               ! to ignore the echo
                nnwor=0                               ! forget all
                nwptr=0
                words(1)=' '
             end if
             last=leng                                ! to break the do while
          else if(subna=='NOREAD') then
             newline=.false.
             nnwor=1
          end if
       end do decode_card

       if((words(1)=='ECHO') .and.&
            (nnpar==0).and.(nnwor==2)) then           ! deal with echo
          if(words(2)=='OFF') then
             echo=.false.
             if(subna/='noecho') write(lisre,20) 'ECOUTE','ECHO OFF'
          else
             echo=.true.
             if(subna/='noecho') write(lisre,20) 'ECOUTE','ECHO ON'
          endif
          nnwor=0                                     ! forget all
          nwptr=0
          do i=1,maxwp
             words(i)=' '
          end do
       else                                           ! print card
          if((echo).and.(firsp<=lastp).and.(subna/='noecho')) then
             if(newline) then
                lastp=lastp+2
                ccard(lastp-1:lastp)=' /'
             end if
             write(lisre,20) subna,ccard(firsp:lastp)
          end if
       end if

    end do ! while ((nnwor==0).and.(nnpar==0).or.newline)

    nwopa=max(npptr,nwptr)
    return
    !
    ! End of include
    !
101 continue
    if(nunit/=lisin.and.nunit/=lisi1) then
       if( endst == 1 .and. stop_end_of_file_opt ) then
          goto 2 ! error
       else
          words(1)='EOF  '
          return
       end if
    end if
    close(unit=nunit)
    if(nunit==lisi1) then
       call set_unit(lisin)  
    else
       call set_unit(lisda)  
    end if
    if(echo.and.(subna/='noecho'))&
         write(lisre,20) 'ECOUTE','END OF INCLUDE FILE'
    resum=1
    go to 100 ! for resume the error return to the same place.
    !
    ! Errors:
    !
1   call message('ECOUTE: ERROR DETECTED WHEN READING')
2   call message('ECOUTE: END OF FILE DETECTED IN SUBROUTINE '//trim(subna))
3   call message('ECOUTE: ERROR: INCLUDE FROM INCLUDE')
4   call message('ECOUTE: TOO MANY PARAM IN COMMAND  ')
5   call message('ECOUTE: TOO MANY WORDS IN COMMAND  ')
6   call message('ECOUTE: BLANK IS ILEGAL HERE       ')
7   call message('ECOUTE: COULD NOT OPEN FILE: '//adjustl(trim(ccard)))
    !
    ! Format
    !
10  format(a250)
20  format(1x,a6,' <-- ',a)

  end subroutine ecoute

  subroutine opincl()
    !
    ! Include file cannot be opened: look the directory of the data file
    !
    integer(ip)              :: istat
    integer(4)               :: istat4,nunit4
    character(20)            :: wstat

    nunit4 = int(nunit,4)
    open(unit=nunit4,file=adjustl(trim(ccard)),err=7,form='FORMATTED',status='OLD',IOSTAT=istat4)
    return

7   call message('ecoute: failed to open '//adjustl(trim(ccard)))

    istat = int(istat4,ip)
    wstat = integer_to_string(istat)
    call message('ECOUTE: COULD NOT OPEN FILE: '//trim(ccard)//'. ERROR '//trim(wstat))

    istat = int(istat4,ip)
    wstat = integer_to_string(istat)

  end subroutine opincl

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-10-20
  !> @brief   Reach a particular section
  !> @details Reach a specific section in the input file
  !> 
  !-----------------------------------------------------------------------

  subroutine ecoute_reach_section(fword,vacal)

    character(5),               intent(in) :: fword
    character(len=*), optional, intent(in) :: vacal !< Calling name

    words(1) = 'XXXXX'

    if( present(vacal) ) then
       do while(words(1)/=fword)
          call ecoute(trim(vacal)) 
       end do
    else
       do while(words(1)/=fword)
          call ecoute('mod_ecoute')
       end do
    end if
    
  end subroutine ecoute_reach_section
   
end module mod_ecoute
!> @}
