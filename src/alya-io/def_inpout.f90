!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



module def_inpout

  !-----------------------------------------------------------------------
  !
  !    This common block contains the parameters needed for some input-
  !    output operations
  !
  !-----------------------------------------------------------------------
  use def_kintyp_basic, only : ip,rp,lg
  
  implicit none
  !
  ! Listen files
  !
  integer(ip), protected :: nunit
  integer(ip), protected :: lisda
  integer(ip), protected :: lisre
  integer(ip), protected :: lispa = 0
  integer(ip), parameter :: lisin = 10
  integer(ip), parameter :: lisi1 = 9
  !
  ! Listen parameters:
  !
  integer(ip), parameter :: maxwp=60
  !
  ! Listen:
  !
  character(251)         :: ccard,wname             ! it must be 250+1
  integer(ip)            :: nwopa,nnwor,nnpar,endst
  real(rp)               :: param(maxwp)
  character(5)           :: words(maxwp)

contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-01-24
  !> @brief   Set the unit to read from
  !> @details Set the unit to read from
  !> 
  !-----------------------------------------------------------------------

  subroutine set_unit(current_unit,RESET)

    integer(ip),           intent(in) :: current_unit
    logical(lg), optional, intent(in) :: RESET
    
    nunit = current_unit
    if( present(RESET) ) then
       if( RESET ) lispa = 1
    end if
    
  end subroutine set_unit

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-01-24
  !> @brief   Set the unit to read from
  !> @details Set the unit to read from
  !> 
  !-----------------------------------------------------------------------

  subroutine set_read_unit(current_unit)

    integer(ip), intent(in) :: current_unit

    lispa = 0
    lisda = current_unit
    if( lisda == lisin ) call runend('MOD_ECOUTE: THIS UNIT IS RESERVED BY DEF_INPOUT')

  end subroutine set_read_unit

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-01-24
  !> @brief   Set the unit to write to 
  !> @details Set the unit to write to
  !> 
  !-----------------------------------------------------------------------

  subroutine set_write_unit(current_unit)

    integer(ip), intent(in) :: current_unit

    lisre = current_unit
    
  end subroutine set_write_unit

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-05-04
  !> @brief   Reda a real value
  !> @details This routine gets the integer value associated with fword
  !>          - texts(1:1).eq.'!' => compulsory parameter
  !>          - texts(1:1).eq.'*' => compulsory parameter
  !> 
  !-----------------------------------------------------------------------

  function getint(fword,defau,textvari)

    character(5)  :: fword
    integer(ip)   :: defau
    character(len = *)                :: textvari
    character(len = len(textvari))    :: texts
    integer(ip)   :: valui
    integer(ip)   :: iword
    logical(lg)   :: found
    character(1)  :: markc
    integer(ip)   :: getint

    valui = defau
    found = .false.
    iword = 0
    texts = textvari
    markc = texts(1:1)

    do while((iword<nwopa).and.(.not.found))
       iword = iword+1
       if(words(iword)==fword) then
          found = .true.
          valui = int(param(iword),ip)
       end if
    end do

    if(((markc=='!').or.(markc=='*')).and.(.not.found)) then
       write(lisre,1) fword,texts(2:35)
       call runend('GETINT: COMPULSORY PARAM. NOT FOUND ('//trim(fword)//', '//trim(textvari)//')')
    end if

    getint = valui

1   format(/,4x,'*** ERROR: ',a5,' = ',a34,/,&
         & 15x,'IS A COMPULSORY PARAMETER. SPECIFY IT !')

  end function getint

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-05-04
  !> @brief   Reda a real value
  !> @details This routine gets the real value associated with fword
  !>          - texts(1:1).eq.'!' => compulsory parameter
  !>          - texts(1:1).eq.'*' => compulsory parameter
  !> 
  !-----------------------------------------------------------------------

  function getrea(fword,defau,textvari)
    
    character(5)  :: fword
    real(rp)  :: defau
    character(len = *)                :: textvari
    character(len = len(textvari))    :: texts
    integer(ip)   :: iword
    logical(lg)   :: found
    real(rp)      :: valur
    character(1)  :: markc
    real(rp)      :: getrea

    valur = defau
    found = .false.
    iword = 0
    texts = textvari
    markc = texts(1:1)

    do while((iword<nwopa).and.(.not.found))
       iword = iword+1
       if(words(iword)==fword) then
          found = .true.
          valur = param(iword)
       end if
    end do

    if(((markc=='!').or.(markc=='*')).and.(.not.found)) then
       write(lisre,1) fword,texts(2:35)
       call runend('GETREA: COMPULSORY PARAM. NOT FOUND')
    end if

    getrea = valur

1   format(/,4x,'*** ERROR: ',a5,' = ',a34,/,&
         & 15x,'IS A COMPULSORY PARAMETER. SPECIFY IT !')

  end function getrea

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-05-04
  !> @brief   Get a character
  !> @details This routine gets the character valuc associated with fword
  !>          - texts(1:1).eq.'!' => compulsory parameter
  !>          - texts(1:1).eq.'*' => compulsory parameter
  !> 
  !-----------------------------------------------------------------------

  function getcha(fword,defau,textvari)
    
    character(5)  :: fword
    character(5)  :: defau
    character(len = *)                :: textvari
    character(len = len(textvari))    :: texts
    character(5)  :: valuc
    integer(ip)   :: iword
    logical(lg)   :: found
    character(1)  :: markc
    character(5)  :: getcha

    valuc=defau
    found=.false.
    iword=0
    texts=textvari
    markc=texts(1:1)

    do while((iword<nwopa).and.(.not.found))
       iword=iword+1
       if(words(iword)==fword) then
          found=.true.
          valuc=words(iword+1)
       end if
    end do

    if(((markc=='!').or.(markc=='*')).and.(.not.found)) then
       write(lisre,1) fword,texts(2:35)
       call runend('GETCHA: COMPULSORY PARAM. NOT FOUND')
    end if

    getcha=valuc

1   format(/,4x,'*** ERROR: ',a5,' = ',a34,/,&
         & 15x,'IS A COMPULSORY PARAMETER. SPECIFY IT !')

  end function getcha

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-05-04
  !> @brief   Get a long character
  !> @details Get a long character
  !> 
  !-----------------------------------------------------------------------

  function getcha_long(fword,defau,textvari,lenstr) result(getcha)
    
    character(5),               intent(in) :: fword
    character(5),               intent(in) :: defau
    character(len=*),           intent(in) :: textvari
    integer(ip),      optional, intent(in) :: lenstr
    integer(ip)                            :: iword
    character(len=:), allocatable          :: getcha
    
    iword = 1
    do while( trim(ccard(iword:iword+4)) /= fword )
       iword = iword + 1
    end do

    iword = iword + 1
    do while( ccard(iword:iword) /= ' ' )
       iword = iword + 1
    end do
    do while( ccard(iword:iword) == ' ' .or. ccard(iword:iword) == ':' .or. ccard(iword:iword) == '=' )
       iword = iword + 1
    end do

    if( present(lenstr) ) then
       getcha = trim(ccard(iword:iword+lenstr-1))
    else
       getcha = trim(ccard(iword:))
    end if
    
  end function getcha_long

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-05-04
  !> @brief   Check if a word exists
  !> @details Check if a word exists
  !> 
  !-----------------------------------------------------------------------

  function exists(fword)
    
    character(5) :: fword
    integer(ip)  :: iword
    logical(lg)  :: exists

    exists=.false.
    iword=0
    do while((iword<nwopa).and.(.not.exists))
       iword=iword+1
       if(words(iword)==fword) then
          exists=.true.
       end if
    end do

  end function exists

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-05-04
  !> @brief   If an option is on
  !> @details This routine chescks if an option is on or off
  !> 
  !-----------------------------------------------------------------------
  
  function option(fword)
    
    character(5) :: fword
    integer(ip)  :: iword
    logical(lg)  :: option

    option = .false.
    iword  = 0
    do while( iword < nwopa )
       iword = iword+1
       if( words(iword) == fword ) then
          if( words(iword+1) == 'ON   ' .or.  words(iword+1) == 'YES  ' ) then
             option = .true.
             return
          else if( words(iword+1) == 'OFF  ' .or.  words(iword+1) == 'NO   ' ) then
             option = .false.
             return
          else
             option = .false.
             call runend('BAD OPTION FOR FIELD '//fword//': SHOULD BE ON OR OFF')
          end if
       end if
    end do
    if( iword == nwopa ) call runend('COULD NOT FIND FIELD '//fword)

  end function option

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-05-04
  !> @brief   If an option is on
  !> @details This routine chescks if an option is not off
  !> 
  !-----------------------------------------------------------------------
  
  function option_not_off(fword)
    
    character(5) :: fword
    integer(ip)  :: iword
    logical(lg)  :: option_not_off

    option_not_off = .true.
    iword  = 0
    do while( iword < nwopa )
       iword = iword+1
       if( words(iword) == fword ) then
          if( words(iword+1) == 'OFF  ' .or.  words(iword+1) == 'NO   ' ) then
             option_not_off = .false.
             return
          end if
       end if
    end do

  end function option_not_off

  subroutine decod1(nstri,strin,lflag,digit)
    !-----------------------------------------------------------------------
    !
    !     This routine decodified a string
    !
    !-----------------------------------------------------------------------
    implicit none
    integer(ip),intent(in)     :: nstri
    character(*),intent(in)    :: strin                            ! (nstri+1)
    integer(ip),intent(out)    :: lflag
    real(rp),intent(out)       :: digit
    character(251) :: stri1
    integer(ip)    :: istri,decim
    !
    ! Special treatmenet if it starts with E,D,e,D,.
    !
    lflag=0
    decim = ichar(strin(1:1))                       ! decimal value
    if( decim < 48 .or. decim > 57 ) then           ! It is not a num.
       if( &
            decim ==  68 .or. decim ==  69 .or. &   ! starts by D or E
            decim == 100 .or. decim == 101 .or. &   ! starts by d or e
            decim ==  46 ) then                     ! starts by .
          lflag = 1
       end if
    end if

    if(lflag/=1) then ! if it does not start by  E,D,e,D,.

       lflag = 0                                          ! It is a parameter
       istri = 1
       do while( istri <= nstri )
          decim = ichar(strin(istri:istri))               ! decimal value
          if( decim < 48 .or. decim > 57 ) then           ! It is not a num.
             if( &
                  decim /=  43 .and. decim /=  45 .and. & ! discard + -
                  decim /=  68 .and. decim /=  69 .and. & ! discard D E
                  decim /= 100 .and. decim /= 101 .and. & ! discard d e
                  decim /=  46 ) then                     ! discard .
                lflag = 1
                istri = nstri
             end if
          end if
          istri = istri+1
       end do

       if( lflag == 0 ) then                              ! It's a number
          do istri = 1,nstri
             stri1(istri:istri) = strin(istri:istri)
          end do
          stri1(nstri+1:nstri+1) = ' '
          read(stri1(1:nstri),*, err=1983) digit                    ! DIGIT <- 'STRIN'
       end if

    end if

    return

1983 call runend('DEF_INPOUT.DECOD1: Error converting string to number:' // stri1(1:nstri) //'; larger string:'//stri1)
!1 call runend('DECOD1: ERROR WHILE READING '//adjustl(trim(strin)))

  end subroutine decod1

  function getpos(fword)
    !-----------------------------------------------------------------------
    !
    !     This routine gets the integer value associated with fword
    !
    !      - texts(1:1).eq.'!' => compulsory parameter
    !      - texts(1:1).eq.'*' => compulsory parameter
    !
    !-----------------------------------------------------------------------
    implicit         none
    character(5)  :: fword
    integer(ip)   :: iword
    logical(lg)   :: found
    integer(ip)   :: getpos

    found=.false.
    iword=0

    do while((iword<nwopa).and.(.not.found))
       iword=iword+1
       if(words(iword)==fword) then
          found=.true.
       end if
    end do

    if( found ) then
       getpos=iword
    else
       call runend('GETPOS: CANNOT FIND WORD FILTER')
    end if

  end function getpos

  subroutine getbigcha(fresu,fword,defau,textvari)
    !-----------------------------------------------------------------------
    !
    !     This routine gets the character value associated with fword
    !
    !      - texts(1:1).eq.'!' => compulsory parameter
    !      - texts(1:1).eq.'*' => compulsory parameter
    !
    !-----------------------------------------------------------------------
    implicit         none
    character(*)  :: fresu
    character(5)  :: fword
    character(5)  :: defau
    character(len = *)                :: textvari
    integer(ip)   :: iword
    !character(*)  :: getbigcha

    iword = 1
    do while( trim(ccard(iword:iword+4)) /= fword )
       iword = iword + 1
    end do

    iword = iword + 1
    do while( ccard(iword:iword) /= ' ' )
       iword = iword + 1
    end do
    do while( ccard(iword:iword) == ' ' .or. ccard(iword:iword) == ':' .or. ccard(iword:iword) == '=' )
       iword = iword + 1
    end do
    fresu = trim(ccard(iword:))

  end subroutine getbigcha


end module def_inpout

