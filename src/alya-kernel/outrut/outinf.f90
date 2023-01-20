!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine outinf
  !-----------------------------------------------------------------------
  !****f* Master/outinf
  ! NAME
  !    outinf
  ! DESCRIPTION
  !    This routine writes some information about the run
  ! OUTPUT
  ! USED BY
  !    Domain
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_elmtyp
  use def_master
  use def_domain
  use mod_memchk
  use mod_communications, only : PAR_SUM
  use mod_communications, only : PAR_AVERAGE
  use mod_outfor,         only : outfor
  use mod_iofile,         only : iofile_flush_unit
  use def_elmgeo,         only : element_type
#if !defined __PGI && defined ISOFORTRANENV
  use iso_fortran_env,    only : compiler_version
  use iso_fortran_env,    only : compiler_options
#endif
  use mod_std
  implicit none
  integer(ip)                      :: ielty,imodu,ibloc,iorde,nchar
  integer(ip)                      :: nelem_exte,nelem_hole,nelem_tota,ielem
  integer(ip)                      :: nboun_exte,nboun_hole,nboun_tota,iboun
  integer(ip)                      :: npoin_hole,npoin_tota,ipoin,imate
  character(100)                   :: celem,ceint,cblem,cbint,cmodu,cbloc,ctimc
  character(7)                     :: crule
  character(20)                    :: clnut
  character(len=:),    allocatable :: mychar(:)
  character(7)                     :: c_version, c_options
  
  if( kfl_rstar /= 2 .and. .not. READ_AND_RUN() ) then

     if( INOTMASTER ) then
        !
        ! Count elements
        !
        nelem_exte = 0
        nelem_hole = 0
        nelem_tota = nelem
        do ielem = 1,nelem
           if( lelch(ielem) == ELEXT ) nelem_exte = nelem_exte + 1
           if( lelch(ielem) == ELHOL ) nelem_hole = nelem_hole + 1
        end do
        nboun_exte = 0
        nboun_hole = 0
        nboun_tota = nboun
        do iboun = 1,nboun
           if( lboch(iboun) == BOEXT ) nboun_exte = nboun_exte + 1
           if( lboch(iboun) == BOHOL ) nboun_hole = nboun_hole + 1
        end do
        npoin_hole = 0
        if( ISEQUEN ) then
           npoin_tota = npoin
        else
           npoin_tota = npoi3
        end if
        do ipoin = 1,npoi3
           if( lnoch(ipoin) == NOHOL ) npoin_hole = npoin_hole + 1
        end do
     end if

     call PAR_SUM(nelem_tota)
     call PAR_SUM(nelem_exte)
     call PAR_SUM(nelem_hole)
     call PAR_SUM(nboun_tota)
     call PAR_SUM(nboun_exte)
     call PAR_SUM(nboun_hole)
     call PAR_SUM(npoin_tota)
     call PAR_SUM(npoin_hole)
     call PAR_AVERAGE(bandw_dom)
     call PAR_AVERAGE(profi_dom)
     
     !-------------------------------------------------------------------
     !
     ! GENERAL DATA
     !
     !-------------------------------------------------------------------
 
     if( INOTSLAVE ) then
        !
        ! System info
        !
        nchar = 100
        allocate( character(len=nchar) :: mychar(7) )
        call get_environment_variable("HOST",         coutp(1))
        call get_environment_variable("USER",         coutp(2))
        call get_environment_variable("BSC_MACHINE",  coutp(3))
        call get_environment_variable("OSTYPE",       coutp(4))
        call get_environment_variable("MACHTYPE",     coutp(5))
        call get_environment_variable("HOSTTYPE",     coutp(6))
        call get_environment_variable("LOADEDMODULES",coutp(7))
        if( len_trim(coutp(3)) == 0 ) coutp(3) = 'No'
        if( len_trim(coutp(7)) == 0 ) coutp(3) = 'None'
                
        mychar(2) = 'USER:             '//trim(coutp(2))
        mychar(1) = 'HOST:             '//trim(coutp(1))
        mychar(3) = 'BSC_MACHINE:      '//trim(coutp(3))
        mychar(4) = 'OSTYPE:           '//trim(coutp(4))
        mychar(5) = 'MACHTYPE:         '//trim(coutp(5))
        mychar(6) = 'HOSTTYPE:         '//trim(coutp(6))
        mychar(7) = 'LOADMODULES:      '//trim(coutp(7))
        
        call outfor(40_ip, lun_syste, CHA_LIST=mychar)
        deallocate(mychar)
        !
        ! Compiler info
        !

#if !defined __PGI && defined ISOFORTRANENV
        c_version = compiler_version()
        c_options = compiler_options()
#else
        c_version = "UNKNOWN"
        c_options = "UNkNWON"
#endif

        nchar = max(len_trim('COMPILER VERSION: '//c_version),&
             &      len_trim('COMPILER OPTIONS: '//c_options),100)
        allocate( character(len=nchar) :: mychar(2) )
        
        
        mychar( 1) = 'COMPILER VERSION: '//c_version
        mychar( 2) = 'COMPILER OPTIONS: '//c_options
        
        call outfor(41_ip, lun_syste, CHA_LIST=mychar)
        deallocate(mychar)
        !
        ! Git revision
        !
        call infrev()        
        !
        ! General data
        !
        cmodu = ''
        cbloc = ''
        do imodu = 1,mmodu
           if( kfl_modul(imodu) /= 0 ) then
              cmodu = trim(cmodu)//trim(namod(imodu))//','
           end if
        end do
        cmodu = trim(cmodu(1:len(trim(cmodu))-1))

        do ibloc = 1,nblok
           cbloc = trim(cbloc)//trim(intost(ibloc))//': '
           do iorde = 1,mmodu
              imodu = lmord(iorde,ibloc)
              if( imodu /= 0 ) then
                 if( kfl_modul(imodu) /= 0 ) then
                    cbloc = trim(cbloc) // trim(namod(imodu)) // ','
                 end if
              end if
           end do
           cbloc = trim(cbloc(1:len(trim(cbloc))))
        end do
        if( len(trim(cbloc))-1 > 0 ) cbloc = trim(cbloc(1:len(trim(cbloc))-1))

        if(kfl_timco==0) then
           ctimc='PRESCRIBED'
        else if(kfl_timco==1) then
           ctimc='MINIMUM OF CRITICAL TIME STEPS'
        else
           ctimc='LOCAL TIME STEP'
        end if
     end if

     !-------------------------------------------------------------------
     !
     ! DOMAIN DATA
     !
     !-------------------------------------------------------------------

     call PAR_SUM(nelty,lnuty) 
     
     if( INOTSLAVE ) then
        celem=''
        ceint=''
        do ielty = iesta_dom,iesto_dom
           if( lexis(ielty) == 1 ) then
              clnut = intost(lnuty(ielty))
              celem = trim(celem)//trim(clnut)//' '//trim(element_type(ielty) % name)//','
              if( lquad(ielty) == 0 ) then
                 crule = '(open)'
              else
                 crule = '(close)'
              end if
              ceint = trim(ceint)//trim(intost(ngaus(ielty)))//trim(crule)//','
           end if
        end do
        celem = trim(celem(1:len(trim(celem))-1))
        ceint = trim(ceint(1:len(trim(ceint))-1))
        cblem = ''
        cbint = ''
        do ielty = ibsta_dom,ibsto_dom
           if( lexis(ielty) == 1 ) then
              clnut = intost(lnuty(ielty))
              cblem = trim(cblem)//trim(clnut)//' '//trim(element_type(ielty) % name)//','
              if( lquad(ielty) == 0 ) then
                 crule = '(open)'
              else
                 crule = '(close)'
              end if
              cbint = trim(cbint)//trim(intost(ngaus(ielty)))&
                   //trim(crule)//','
           end if
        end do
        cblem = trim(cblem(1:len(trim(cblem))-1))
        cbint = trim(cbint(1:len(trim(cbint))-1))
        !
        ! Write information
        !
        coutp(1) = adjustl(trim(cmodu))
        coutp(2) = adjustl(trim(cbloc))
        coutp(3) = adjustl(trim(ctimc))
        if( kfl_timco == -1 ) then
           routp(5) = 0.0_rp
           routp(6) = 0.0_rp
           routp(7) = 0.0_rp
        else
           routp(5) = dtime
           routp(6) = timei
           routp(7) = timef
        end if
        ioutp(1)  = ndime
        ioutp(2)  = npoin_tota
        ioutp(11) = npoin_hole
        routp(1)  = vodom
        routp(2)  = voave
        routp(3)  = vomin
        ioutp(3)  = elmin
        routp(4)  = vomax
        ioutp(4)  = elmax
        ioutp(5)  = nelem_tota
        ioutp(6)  = nelem_exte
        ioutp(7)  = nelem_hole
        coutp(4)  = adjustl(trim(celem))
        coutp(5)  = adjustl(trim(ceint))
        ioutp(8)  = nboun_tota
        ioutp(9)  = nboun_exte
        ioutp(10) = nboun_hole
        coutp(6)  = adjustl(trim(cblem))
        coutp(7)  = adjustl(trim(cbint))
        call outfor(17_ip,lun_outpu,' ')

        !-------------------------------------------------------------------
        !
        ! MATERIAL DATA
        !
        !-------------------------------------------------------------------

        call outfor(106_ip,lun_outpu,' ')
        do imate = 1,nmate
           ioutp(1) = imate
           routp(1) = vomat(imate)
           call outfor(107_ip,lun_outpu,' ')
        end do
           
     end if

  end if
  
  if( INOTSLAVE ) call iofile_flush_unit(lun_outpu)

end subroutine outinf
