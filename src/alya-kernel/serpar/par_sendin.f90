!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine par_sendin()
  !------------------------------------------------------------------------
  !****f* Parall/par_send
  ! NAME
  !    par_send
  ! DESCRIPTION
  !    This routine Send all buffers to process 'kfl_desti'
  ! OUTPUT
  !   
  ! USED BY
  !    Parall
  !***
  !------------------------------------------------------------------------
  use def_parall
  use def_master
  use mod_memchk
  use mod_iofile
  use mod_par_virfil
  use mod_communications_point_to_point, only : PAR_SEND
  use def_mpi
  use mod_memory
  implicit none
  
  integer(4)             :: iunit4
  integer(ip)            :: iunit,ivari,nparr2,iparx
  real(rp)               :: time1,time2  
  character(150)         :: cfile
  real(rp),      pointer :: parrx(:)

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
        call PAR_SEND(npari,parin,DOM_I=kfl_desti_par)
        npari=0
     end if

     if( nparr > 0 ) then
        call PAR_SEND(nparr,parre,DOM_I=kfl_desti_par)
        nparr=0
     end if

     if( nparc > 0 ) then
        call PAR_SEND(nparc,parch,DOM_I=kfl_desti_par)
        nparc=0
     end if

     if( nparx > 0 ) then
        nparr2  = 2*nparx
        allocate( parrx(nparr2) )
        ivari = 0
        do iparx = 1,nparx
           ivari = ivari + 1
           parrx(ivari) = real  ( parcx(iparx) )
           ivari = ivari + 1
           parrx(ivari) = aimag ( parcx(iparx) )
        end do
        call PAR_SEND(nparr2,parrx,DOM_I=kfl_desti_par)
        nparx=0
        deallocate( parrx )
     end if

#endif

  else if( ( IMASTER .or. nproc_par > 1 ) .and. PART_AND_WRITE() ) then

     !-------------------------------------------------------------------
     !
     ! File communication
     !
     !-------------------------------------------------------------------
     
     iunit  = lun_aonlp_par + kfl_desti_par
     iunit4 = int(iunit,4) 

     if( kfl_virfi_par == 1 ) then
        !
        ! Virtual file
        !
        call par_wribuf(kfl_desti_par)

     else 
        !
        ! Binary/ASCII format
        !
        if( kfl_filio_par == 1 ) then
           call par_filnam(1_ip,kfl_desti_par,fil_rstar_par,cfile)
           if( kfl_ascii_par == 0 ) then
              call iofile(zero,iunit,trim(cfile),'PARALL RESTART','old','unformatted','append')
           else
              call iofile(zero,iunit,trim(cfile),'PARALL RESTART','old','formatted',  'append')      
           end if
        end if

        call par_wrifil(iunit)

        if( kfl_filio_par == 1 )  close(iunit4)

     end if

     npari = 0
     nparr = 0
     nparc = 0
     nparx = 0

  end if

  call cputim(time2)
  cpu_paral(21)=cpu_paral(21)+time2-time1

end subroutine par_sendin
