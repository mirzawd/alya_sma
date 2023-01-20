!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine par_operat(itask)
  !------------------------------------------------------------------------
  !****f* Parall/par_operat
  ! NAME
  !    par_operat
  ! DESCRIPTION
  !    This routine operates between slave arrays and broadcast the result
  !    to all slaves and master
  !    ITASK=1 ... Minimum 
  !    ITASK=2 ... Maximum
  !    ITASK=3 ... Sum
  ! OUTPUT
  !    NPARI ..... Integer array
  !    NPARR ..... Real array
  ! USED BY
  !    Parall
  !***
  !------------------------------------------------------------------------
  use def_parall
  use def_master
  use mod_parall, only : PAR_INTEGER
  use mod_parall, only : PAR_REAL
  use mod_parall, only : PAR_COMM_MY_CODE
  use def_mpi
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: ii 
  integer(4)              :: istat=0,npari4,nparr4,nparx4
  real(rp)                :: time1,time2
  integer(ip), pointer    :: iwa(:)
  real(rp),    pointer    :: rwa(:)
  complex(rp), pointer    :: cwa(:)
  !
  ! Operations on real arrays
  !

  if(nparr>0) then 

#ifdef MPI_OFF
#else
     call cputim(time1)     
     allocate(rwa(nparr))
     if(ISLAVE) then
        do ii=1,nparr
           rwa(ii)=parre(ii)
        end do
     end if
     nparr4=int(nparr,4)
     if(itask==1) then
        !
        ! Minimum
        !
        if(IMASTER) then
           do ii=1,nparr
              rwa(ii)=huge(1.0_rp)
           end do
        end if
        call MPI_AllReduce(rwa,parre,nparr4,PAR_REAL,MPI_MIN,PAR_COMM_MY_CODE,istat) 
     else if(itask==2) then
        !
        ! Maximum
        !
        if(IMASTER) then
           do ii=1,nparr
              rwa(ii)=-huge(1.0_rp)
           end do
        end if
        call MPI_AllReduce(rwa,parre,nparr4,PAR_REAL,MPI_MAX,PAR_COMM_MY_CODE,istat)        
     else if(itask==3) then
        !
        ! Sum
        !
        if(IMASTER) then
           do ii=1,nparr
              rwa(ii)=0.0_rp
           end do
        end if
        call MPI_AllReduce(rwa,parre,nparr4,PAR_REAL,MPI_SUM,PAR_COMM_MY_CODE,istat)        
     end if
     deallocate(rwa)
     call cputim(time2)
     cpu_paral(22)=cpu_paral(22)+time2-time1
     nparr=0
#endif

  end if
  !
  ! Operations on integer arrays
  !
  if(npari>0) then 

#ifdef MPI_OFF
#else
     call cputim(time1)     
     allocate(iwa(npari))
     if(ISLAVE) then
        do ii=1,npari
           iwa(ii)=parin(ii)
        end do
     end if
     npari4=int(npari,4)
     if(itask==1) then
        !
        ! Minimum
        !
        if(IMASTER) then
           do ii=1,npari
              iwa(ii)=huge(1_ip)
           end do
        end if
        call MPI_AllReduce(iwa,parin,npari4,PAR_INTEGER,MPI_MIN,PAR_COMM_MY_CODE,istat) 

     else if(itask==2) then
        !
        ! Maximum
        !
        if(IMASTER) then
           do ii=1,npari
              iwa(ii)=-huge(1_ip)
           end do
        end if
        call MPI_AllReduce(iwa,parin,npari4,PAR_INTEGER,MPI_MAX,PAR_COMM_MY_CODE,istat)  

     else if(itask==3) then
        !
        ! Sum
        !
        if(IMASTER) then
           do ii=1,npari
              iwa(ii)=0
           end do
        end if
        call MPI_AllReduce(iwa,parin,npari4,PAR_INTEGER,MPI_SUM,PAR_COMM_MY_CODE,istat)        
     end if
     deallocate(iwa)
     call cputim(time2)
     cpu_paral(24)=cpu_paral(24)+time2-time1
     npari=0
#endif

  end if

  !
  ! Operations on complex arrays
  !
  if(nparx>0) then 

#ifdef MPI_OFF
#else
     call cputim(time1)     
     allocate(cwa(nparx))
     if(ISLAVE) then
        do ii=1,nparx
           cwa(ii)=parcx(ii)
        end do
     end if
     nparx4=int(nparx,4)
     if(itask==1) then
        !
        ! Minimum
        !
        if(IMASTER) then
           do ii=1,nparx
              cwa(ii)=huge(1.0_rp)
           end do
        end if
        call MPI_AllReduce(cwa,parcx,nparx4,MPI_DOUBLE_COMPLEX,MPI_MIN,PAR_COMM_MY_CODE,istat) 
     else if(itask==2) then
        !
        ! Maximum
        !
        if(IMASTER) then
           do ii=1,nparx
              cwa(ii)=-huge(1.0_rp)
           end do
        end if
        call MPI_AllReduce(cwa,parcx,nparx4,MPI_DOUBLE_COMPLEX,MPI_MAX,PAR_COMM_MY_CODE,istat)        
     else if(itask==3) then
        !
        ! Sum
        !
        if(IMASTER) then
           do ii=1,nparx
              cwa(ii)=(0.0_rp,0.0_rp)
           end do
        end if
        call MPI_AllReduce(cwa,parcx,nparx4,MPI_DOUBLE_COMPLEX,MPI_SUM,PAR_COMM_MY_CODE,istat)        
     end if
     deallocate(cwa)
     call cputim(time2)
     cpu_paral(22)=cpu_paral(22)+time2-time1
     nparx=0
#endif

  end if
  
  if(istat/=0) call runend('PARALL: FUNCTION MPI_ALLREDUCE HAS FAILED')

end subroutine par_operat
