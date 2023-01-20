!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine par_receiv()
!------------------------------------------------------------------------
!****f* Parall/par_receiv
! NAME
!    par_receiv
! DESCRIPTION
!    This routine Receive all buffers from process 'kfl_desti'
! OUTPUT
!   
! USED BY
!    Parall
!***
!------------------------------------------------------------------------
  use def_parall
  use def_master
  use mod_iofile
  use mod_communications_point_to_point, only : PAR_RECEIVE
  use mod_memory
  use def_mpi
  implicit none

  integer(4)     :: iunit4
  integer(ip)    :: ipari,iparr,iunit
  real(rp)       :: time1,time2
  character(300) :: messa
  character(20)  :: cnume

  call cputim(time1)

  if( PART_AND_RUN() ) then

     !-------------------------------------------------------------------
     !
     ! MPI communication
     !
     !-------------------------------------------------------------------

#ifdef MPI_OFF
#else
     if( npari > 0 ) then
        call PAR_RECEIVE(npari,parin,DOM_I=kfl_desti_par)
        npari = 0
     end if
     
     if( nparr > 0 ) then
        call PAR_RECEIVE(nparr,parre,DOM_I=kfl_desti_par)
        nparr = 0
     end if

     if( nparc > 0 ) then
        call PAR_RECEIVE(nparc,parch,DOM_I=kfl_desti_par)
        nparc = 0
     end if
#endif

  else 

     !-------------------------------------------------------------------
     !
     ! Read from file
     !
     !-------------------------------------------------------------------

     iunit  = lun_aonlp_par !+ kfl_paral
     iunit4 = int(iunit,4)

     if(kfl_ascii_par==0) then
        read(iunit4,err=1) npari,nparr,nparc
        if( npari > 0 ) read(iunit4,err=1,end=1)  ( parin(ipari),   ipari=1,npari )
        if( nparr > 0 ) read(iunit4,err=1,end=1)  ( parre(iparr),   iparr=1,nparr )
        if( nparc > 0 ) read(iunit4,err=1,end=1)    parch(1:nparc)
     else
        read(iunit4,*,err=1) npari,nparr,nparc
        read(iunit4,*,err=1) strin,strre,strch
        if( npari > 0 ) read(iunit4,*,err=1,end=1) ( parin(ipari),  ipari=1,npari )
        if( nparr > 0 ) read(iunit4,*,err=1,end=1) ( parre(iparr),  iparr=1,nparr )
        if( nparc > 0 ) read(iunit4,*,err=1,end=1)   parch(1:nparc)      
     end if
     npari=0
     nparr=0
     nparc=0

  end if

  call cputim(time2)
  cpu_paral(22)=cpu_paral(22)+time2-time1
  return

1 cnume=intost(kfl_paral)
  messa='PARALL: ERROR WHILE SLAVE '//trim(cnume)&
       //' IS READING RESTART FILE. CHECK FILE FORMAT.'&
       //'TRYING TO READ: '//trim(strin)&
       //', '//trim(strre)//', '//trim(strch)
  call runend(messa)

end subroutine par_receiv
