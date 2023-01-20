!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine par_broadc()
  !------------------------------------------------------------------------
  !****f* Parall/par_broadc
  ! NAME
  !    par_broadc
  ! DESCRIPTION
  !    This routine Send/Receive all buffers
  ! OUTPUT
  !   
  ! USED BY
  !    Parall
  !***
  !------------------------------------------------------------------------
  use def_parall
  use def_parame
  use def_master
  use mod_iofile
  use mod_par_virfil
  use mod_parall, only : PAR_COMM_MY_CODE,PAR_INTEGER
  use mod_parall, only : PAR_REAL
  use mod_memory
  use def_mpi
  implicit none

  integer(4)             :: istat,npari4
  integer(4)             :: nparr4,nparc4,nparx4,nparl4
  integer(ip)            :: iunit
  integer(4)             :: iunit4
  integer(ip)            :: npart_loc
  real(rp)               :: time1,time2
  character(300)         :: messa
  character(150)         :: cfile
  character(20)          :: cnume

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
        npari4 = int(npari,4)
        call MPI_Bcast( parin(1:npari), npari4, PAR_INTEGER, 0_4, PAR_COMM_MY_CODE, istat )
        npari = 0
     end if

     if( nparr > 0 ) then
        nparr4 = int(nparr,4)
        call MPI_Bcast( parre(1:nparr), nparr4, PAR_REAL, 0_4, PAR_COMM_MY_CODE, istat )
        nparr = 0
     end if

     if( nparc > 0 ) then
        nparc4 = int(nparc,4)
        call MPI_Bcast( parch, nparc4, MPI_CHARACTER, 0_4, PAR_COMM_MY_CODE, istat )
        nparc = 0
     end if

     if( nparx > 0 ) then
        nparx4 = int(nparx,4)
        call MPI_Bcast( parcx(1:nparx), nparx4, MPI_DOUBLE_COMPLEX, 0_4, PAR_COMM_MY_CODE, istat )
        nparx = 0
     end if

     if( nparl > 0 ) then
        nparl4 = int(nparl,4)
        call MPI_Bcast( parlo(1:nparl), nparl4, MPI_LOGICAL, 0_4, PAR_COMM_MY_CODE, istat )
        nparl = 0
     end if

#endif

  else if( ( IMASTER .or. nproc_par > 1 ) .and. PART_AND_WRITE() ) then

     !-------------------------------------------------------------------
     !
     ! Sequential preprocess: Master parts and writes restart files
     ! Parallel preprocess: Everybody does it!
     !
     !-------------------------------------------------------------------

     if( nproc_par > 1 ) then
        npart_loc = 0
     else
        npart_loc = npart_par
     end if
     
     if( kfl_virfi_par == 1 ) then
        !
        ! Virtual files
        !
        do kfl_desti_par = 0,npart_loc
           call par_wribuf(kfl_desti_par)
        end do
        
     else
        !
        ! ASCII/Binary format
        !
        
        do kfl_desti_par = 0,npart_loc
           iunit  = lun_aonlp_par + kfl_desti_par
           iunit4 = int(iunit,4)
           if( kfl_filio_par == 1 ) then
              call par_filnam(1_ip,kfl_desti_par,fil_rstar_par,cfile)
              if (kfl_ascii_par == 0 ) then
                 call iofile(zero,iunit,trim(cfile),'PARALL RESTART','old','unformatted','append')
              else
                 call iofile(zero,iunit,trim(cfile),'PARALL RESTART','old','formatted','append')
              end if
           end if
           call par_wrifil(iunit)
           if( kfl_filio_par == 1 )  close(iunit4)
        end do

     end if
     npari=0
     nparr=0
     nparc=0

  else if( READ_AND_RUN() ) then

     !-------------------------------------------------------------------
     !
     ! Read restart files
     !
     !-------------------------------------------------------------------

     iunit  = 1000 + kfl_paral
     iunit4 = int(iunit,4)
     if(kfl_ascii_par == 0 ) then
        read(iunit4,err=2) npari,nparr,nparc
        if( npari > 0 ) read(iunit4,err=2,end=2) parin(1:npari)
        if( nparr > 0 ) read(iunit4,err=2,end=2) parre(1:nparr)
        if( nparc > 0 ) read(iunit4,err=2,end=2) parch(1:nparc)
        if( nparl > 0 ) read(iunit4,err=2,end=2) parlo(1:nparl)
     else
        read(iunit4,*,err=2) npari,nparr,nparc
        read(iunit4,*,err=2) strin,strre,strch
        if( npari > 0 ) read(iunit4,*,err=2,end=2) parin(1:npari)
        if( nparr > 0 ) read(iunit4,*,err=2,end=2) parre(1:nparr)
        if( nparc > 0 ) read(iunit4,*,err=2,end=2) parch(1:nparc)
        if( nparl > 0 ) read(iunit4,*,err=2,end=2) parlo(1:nparl)
     end if
     npari = 0
     nparr = 0
     nparc = 0
     nparl = 0

  end if

  call cputim(time2)
  cpu_paral(23)=cpu_paral(23)+time2-time1
  return

1 cnume=intost(kfl_desti_par)
  messa='PARALL: ERROR WHILE MASTER IS WRITING SLAVE '//trim(cnume)&
       //' RESTART FILE. CHECK FILE FORMAT'
  call runend(trim(messa))
2 cnume=intost(kfl_paral)   
  messa='PARALL: ERROR WHILE SLAVE '//trim(cnume)&
       //' IS READING RESTART FILE. CHECK FILE FORMAT.'&
       //'TRYING TO READ: '//trim(strin)&
       //', '//trim(strre)//', '//trim(strch)
  call runend(trim(messa))

end subroutine par_broadc
