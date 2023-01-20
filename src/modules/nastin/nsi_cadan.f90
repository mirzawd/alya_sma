!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_cadan(itask)
  !-----------------------------------------------------------------------
  !****f* Nastin/nsi_cadan
  ! NAME 
  !    nsi_cadan !Coupling with ADAN
  ! DESCRIPTION
  !     This subroutine executes several task for coupling with ADAN
  !
  ! USES
  !   -
  ! USED BY
  !    nsi_openfi
  !    nsi_endstei
  !    nsi_begste
  !    nsi_turnof
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_elmtyp
  use def_master
  use def_kermod
  use def_domain
  use def_nastin
  use mod_ker_proper
  use mod_communications, only : PAR_BROADCAST
  use mod_output_postprocess, only : output_postprocess_boundary_sets_parall
  implicit none

  integer(ip), intent(in)                :: itask
  real(rp)                               :: aflux, auxpres
  real(rp), save                         :: q0, p0 ! Flux and pressure in the conection boundary
  integer(ip)                            :: kbset, auxntran, nn
  !  integer(ip)                            :: CheckConv, CheckConvm1
  integer(ip),save                       :: sz, oldsz, ntran, oldnt
  integer(ip), save                      :: old_If_couit ! Last value of IFace iteration
  character(len=18)                      :: stringb
  logical, save                          :: temp_on !Existence of temporary file
  logical, parameter                     :: DEBUG_MODE=.False. !Long file debug mode



  if (kfl_cadan_nsi==1) then


     select case (itask)

     case (1_ip)

        !-------------------------------------------------------------------!
        !
        ! Open
        !
        !-------------------------------------------------------------------!

        if(INOTSLAVE .and. DEBUG_MODE) then
           open(unit=314, file="../../ALYA2COUPLING.send",status='unknown', &
                access='sequential', position='append', action='write') 
        endif

     case (2_ip)

        !-------------------------------------------------------------------!
        !
        ! Read
        !
        !-------------------------------------------------------------------!

        if(INOTSLAVE) then

           if (ittim==1 .and. itcou==1) then
              inquire(file='../../C2AL.temp',exist=temp_on)
              if(temp_on) then 
                 write(6,*) 'Temporary files for distributed execution detected.'
                 call ReadTempFile(.True.)
              elseif(DEBUG_MODE) then
                 call ReadLongFile(.True.)
              else
                 call RUNEND('NSI_CADAN: NOT TEMP FILE NOR DEBUG_MODE DETECTED')
              endif
           else
              if(temp_on) then
                 call ReadTempFile(.False.)
              elseif(DEBUG_MODE) then 
                 call ReadLongFile(.False.)
              else
                 call RUNEND('NSI_CADAN: NOT TEMP FILE NOR DEBUG_MODE DETECTED')
              endif
           endif

        endif

        call PAR_BROADCAST(press_cadan_nsi,'IN MY CODE')

     case (3_ip)

        !-------------------------------------------------------------------!
        !
        ! Write
        !    
        !-------------------------------------------------------------------!

        if(ittim==1 .and. itcou==1) then
           if(INOTMASTER) then
              do kbset = 1,nbset
                if(lbsec(kbset).eq.kfl_aiobo_nsi) call nsi_bouset(lbsec(kbset),kbset)
                ! call nsi_bouset(lbsec(kbset),kbset)
              enddo
           endif
           call output_postprocess_boundary_sets_parall()

           if(INOTSLAVE) then
              if(kfl_aiobo_nsi==0) call RUNEND('NSI_CADAN: A CONNECTION BOUNDARY WITH ADAN MUST BE SPECIFIED')

              if(DEBUG_MODE) then
                    allocate (Q_cadan_nsi(nbset))
                    allocate (P_cadan_nsi(nbset))
              endif

              nn    =  postp(1) % nvabs + 1 !Not used anymore

              do kbset= 1, nbset
                if(lbsec(kbset).eq.kfl_aiobo_nsi) then
                    q0=vbset(16,kbset) !*vbset(nn,kbset)
                    p0=vbset(1,kbset)
                endif
                if(DEBUG_MODE) then
                    Q_cadan_nsi(kbset)=vbset(16,kbset) !*vbset(nn,kbset) 
                    P_cadan_nsi(kbset)=vbset(1,kbset)                  
                endif
              enddo
           endif      
        endif


        if(INOTSLAVE) then
           if(.not.DEBUG_MODE) then
                call WriteTempFile(ittim==1 .and. itcou==1)
           else
                call WriteLongFile(ittim==1 .and. itcou==1)
           endif
           write(6,*) 'Info enviada a ADAN! |---->'
        endif

        if( INOTMASTER ) then
           do kbset = 1,nbset 
              if(lbsec(kbset).eq.kfl_aiobo_nsi) call nsi_bouset(lbsec(kbset),kbset)
              !call nsi_bouset(lbsec(kbset),kbset)
           end do
        end if
        call output_postprocess_boundary_sets_parall()

        if(INOTSLAVE) then
            nn    =  postp(1) % nvabs + 1
            do kbset= 1, nbset
                if(lbsec(kbset).eq.kfl_aiobo_nsi) then
                    q0=vbset(16,kbset) !*vbset(nn,kbset)
                    p0=vbset(1,kbset)
                endif
                if(DEBUG_MODE) then
                    Q_cadan_nsi(kbset)=vbset(16,kbset) !*vbset(nn,kbset) 
                    P_cadan_nsi(kbset)=vbset(1,kbset) 
                endif                 
            enddo
        endif
    
     case (4_ip)
        !-------------------------------------------------------------------!
        !
        ! To be defined
        !
        !-------------------------------------------------------------------!

     case (7_ip)
        !-------------------------------------------------------------------!
        !
        ! Close
        !
        !-------------------------------------------------------------------!
        if (INOTSLAVE) then
           close(314)
           close(413)
           deallocate(Q_cadan_nsi)
           deallocate(P_cadan_nsi)
        endif


     end select

  endif

  !------------------------------------------------------------------------!
  !
  !                             CONTAINS
  !
  !------------------------------------------------------------------------!
contains 

  !------------------------|
  !write the temporary file|
  !------------------------|
  subroutine WriteTempFile(IterOneFlag)
    implicit none
    logical, intent(in)              :: IterOneFlag
    character(len=*), parameter      :: TempOutputFormat="(i8, ' ', i8, ' ', ES15.8,' ',ES15.8,' ', ES15.8/)"

    open(unit=320, file="../../AL2C.temp",status='replace', &
         access='sequential', position='asis', action='write') 


    if(IterOneFlag) then
       do kbset=1, nbset
          write(320,*) '#', kbset, lbsec(kbset)
       enddo
       write(320,TempOutputFormat,advance="no") ittim, itcou, cutim, q0*0.6_rp, p0*0.6_rp
    else
       write(320,TempOutputFormat,advance="no") ittim, itcou, cutim, q0, p0
    end if

    close(320)

  end subroutine WriteTempFile

  !-------------------|
  !write the long file|
  !-------------------|
  subroutine WriteLongFile(IterOneFlag)
    implicit none
    logical                          :: IterOneFlag
    character(len=*), parameter      :: LongOutputFormat="(ES15.8,' ', ES15.8,' ', ES15.8,' ', ES15.8,' ', ES15.8,' ', ES15.8,' ', ES15.8/)"

    if(IterOneFlag) then
       call WriteHeader(314_ip)
       write(314,"(ES15.8, 6ES15.8/)",advance="no") cutim, Q_cadan_nsi(:)*0.6_rp, P_cadan_nsi(:)*0.6_rp
    else
       write(314,"(ES15.8, 6ES15.8/)",advance="no") cutim, Q_cadan_nsi(:), P_cadan_nsi(:)
    end if

  end subroutine WriteLongFile

  !-----------------|
  !write the header |
  !-----------------|
  subroutine WriteHeader(FileCode)
    implicit none
    integer(ip)     :: FileCode
    write(FileCode,'(a76)') "# C1 timestep | C2-4 meanvelocity (W-IN-OUT) | C5-7 mean pressure (W-IN-OUT)"
    do kbset=1, nbset
       write(FileCode,*) '#', kbset, lbsec(kbset)
    enddo
  end subroutine WriteHeader

  !----------------------------------------|
  !Obtains the file size value accurrately |
  !----------------------------------------|
  function AccurrateSize() result(sze)
    implicit none
    integer(ip)         :: sze

    close(413)
    open (413, file='../../COUPLING2ALYA.send', status='unknown', &
         access='sequential', action='read' ) 

    if(sz.gt.10) then
       read(413,'(a18)') stringb
       auxntran=0
       do while (auxntran.lt.oldnt)
          read(413,'(1i8, ES15.8, ES15.8)') auxntran, auxpres, aflux
       enddo
    endif

    inquire(unit=413, size=sze)

  endfunction AccurrateSize
  !--------------------|
  !Reads the long file |
  !--------------------|
  subroutine ReadLongFile(itone)
    implicit none
    logical, intent(in)     :: itone


    open(unit=413, file="../../COUPLING2ALYA.send",status='old',  &
         access='sequential', action='read')
    inquire(unit=413, size=sz)

    if(itone) then
       oldsz=sz
       oldnt=0
       ntran=0
    endif

    write(6,*) "Esperando info de ADAN..."
    !        write(6,*) 'DEBUGDEBUGDEBUG ittim: ', ittim, '|itcou:', itcou,'|itinn: ', itinn(modul)
    !        write(6,*) 'DEBUGDEBUGDEBUG sz: ', sz, '| oldsz: ', oldsz

    if (.not.itone) then
       waiting: do while ( (sz.lt.(10+oldsz)) )
          sz=AccurrateSize()
          write(6,*) '             Durmiendo...'
          !call sleep(1)
          write(6,*) '             Up!'
       enddo waiting
    else
       read(413,'(a18)') stringb
    endif

    read(413,'(1i8, ES15.8, ES15.8)') ntran, press_cadan_nsi, aflux

    oldnt=ntran
    oldsz=sz

    !write(6,'(1i8, ES15.8,ES15.8)') ntran, press_cadan_nsi, aflux
    write(6,*) 'Info conseguida! <---|'

  end subroutine ReadLongFile

  !-------------------------|
  !Reads the temporary file |
  !-------------------------|
  subroutine ReadTempFile(itone)
    implicit none
    logical, intent(in)             :: itone !Iteration One Flag
    logical                         :: cont  ! Continuation flag
    logical                         :: ex    ! File existence flag
    integer(ip)                     :: If_timit ! IFace time Iteration
    integer(ip)                     :: If_couit ! IFace coupling Iteration
    integer(ip)                     :: sleepcount ! Time code goes to sleep
    integer(ip), parameter          :: maxsleepcount=1000_ip
    integer(ip)                     :: ierr, ierr2
    real(rp)                        :: p0_read, q0_read   ! Iface recieved pressure and flux              
    character(len=*), parameter     :: TempInputFormat="(i8, ' ', i8, ' ',ES15.8,' ', ES15.8/)"

    cont=.False.
    ex=.False.
    sleepcount=0_ip

    if(itone) old_If_couit=0

    write(6,*) 'Waiting for IFace temporary file...'

!
!If the file is not existent yet, wait for it
!
        do while(.not.ex)
            inquire(file='../../C2AL.temp', exist=ex)
                if(ex) then
                    open (323, file='../../C2AL.temp', status='old', &
                          access='sequential', action='read')
                    !call sleep(1)   ! Este PAUSE esta aca para asegurarnos que el archivo se escriba.
                else
                    if(sleepcount == 0) write(6,*) 'Still waiting...'
                    sleepcount=sleepcount+1
                    !call sleep(1)
                    !!! if(old_Al_couit.eq.0_ip) iteroneflag=.True. !!! ESTE ITERONEFLAG ME PARECE QUE NO LO NECESITO
                endif
        enddo

        if (sleepcount .gt. 0) write(6,*) '|-------------- Up after ', sleepcount, ' naps.'
        if(sleepcount.eq.maxsleepcount) stop 'TEMPORARY INPUT FILE NOT FOUND. PROGRAM FINISHED.'
!
! End of existence of file waiting loop (the file exists)
!

       read(323,TempInputFormat, iostat=ierr) If_timit, If_couit, p0_read, q0_read

!!!!
!!!!  ERROR DETECTION CODE
!!!!
        if(ierr.eq.0_ip) then
            close(323,status='delete')
        elseif(ierr.lt.0_ip) then
            write(6,*) 'WARNING: EOF detected! Trying again...'
            !call sleep(1)
            close(323)
            open (323, file='../../C2AL.temp', status='old', &
                       access='sequential', action='read')

            read(323,TempInputFormat, iostat=ierr2) If_timit, If_couit, p0_read, q0_read

            if (ierr2.eq.0_ip) then
                write(6,*) 'Woof! Made it this time!'
                close(323,status='delete')
            else
                write(6,*) 'PROGRAM STOPED: The file has a problem, wont be deleted'
                stop 'PROGRAM STOPPED DUE TO ERROR IN READ IN TEMP FILE'
            endif

        elseif(ierr.gt.0_ip) then
            write (6,*) 'WARNING: shit was readed from temp file!'
            !call sleep(1)
            close(323)
            open (323, file='../../C2AL.temp', status='old', &
                       access='sequential', action='read')

            read(323,TempInputFormat, iostat=ierr2) If_timit, If_couit, p0_read, q0_read

            if (ierr2.eq.0_ip) then
                write(6,*) 'Woof! Made it this time!'
                close(323,status='delete')
            else
                write(6,*) 'PROGRAM STOPED: The file has a problem, wont be deleted'
                stop 'PROGRAM STOPPED DUE TO ERROR IN READ IN TEMP FILE'
            endif

        endif
!!!!
!!!! END  ERROR DETECTION CODE
!!!!


    if (itone.and.((If_couit .ne. 1) .or. (If_timit .ne. 1))) then
        call RUNEND('temp file not of the first iteration.')
    endif

    press_cadan_nsi=p0_read
    old_If_couit=If_couit

  end subroutine ReadTempFile
end subroutine nsi_cadan
