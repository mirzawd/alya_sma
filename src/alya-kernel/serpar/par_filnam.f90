!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine par_filnam(itask,ipart,wfili,wfilo)
  !-----------------------------------------------------------------------
  !****f* Parall/par_filnam
  ! NAME
  !    par_filnam
  ! DESCRIPTION
  !    Compose restart file name according to file hierarchy. For 
  !    process "n", the file names are:
  !
  !    KFL_FILEH_PAR = 0 ... No hierarchy:    NAME.par.rst"n"
  !    KFL_FILEH_PAR = 1 ... Hierarchy is on: PAR*****/NAME.par.rst"n"
  !                  
  !    where ***** is the range of the process using step of 100:
  !    PAR00000: master restart file
  !    PAR00100:   1 to 100 process restart files
  !    PAR00200: 101 to 200 process restart files, etc.
  !
  !    There also can be a prefix:
  !
  !    KFL_PREFI_PAR = 0 ... No prefix, file name as is.
  !    KFL_PREFI_PAR = 1 ... Files are nbames WPREF_PAR/PAR00000...
  !
  ! USED BY
  !    par_openfi
  !*** 
  !-----------------------------------------------------------------------
  use def_kintyp
  use def_parame
  use def_master
  use def_domain
  use def_parall
  implicit none
  integer(ip),    intent(in)    :: ipart,itask
  character(150), intent(in)    :: wfili
  character(150), intent(inout) :: wfilo
  character(6)                  :: wdire
  character(20)                 :: winte
  integer(ip)                   :: ii
  !
  ! Add process number to tail
  !  
  winte=intost(ipart)
  wfilo=trim(wfili)//trim(winte)  
  !
  ! File hierarchy: add directory as prefix
  !
  if(kfl_fileh_par/=0) then
     if(ipart==0) then
        wdire='000000'
     else
        ii=100*(int(real(ipart-1)/100.0)+1)
        if(ii<1000) then
           write(wdire,'(a,i3)') '000',ii
        else if(ii<10000) then
           write(wdire,'(a,i4)') '00', ii
        else if(ii<100000) then
           write(wdire,'(a,i5)') '0', ii
        else 
           write(wdire,'(i6)')        ii
        end if
     end if
     if(kfl_fileh_par==1) then
        !
        ! One-level hierarchy
        !
        wfilo='PAR'//wdire//'/'//trim(wfilo)
     else
        !
        ! Two-level hierarchy
        !
        if(itask==1) then
           wfilo='PAR'//wdire//'/data/'//trim(wfilo)
        else
           wfilo='PAR'//wdire//'/results/'//trim(wfilo)
        end if
     end if
  end if

end subroutine par_filnam
