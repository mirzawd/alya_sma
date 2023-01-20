!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine par_slaves(itask)
  !-----------------------------------------------------------------------
  !****f* Parall/par_slaves
  ! NAME
  !    par_slaves
  ! DESCRIPTION
  !    This subroutine exchange arrays between slaves
  ! USED BY
  !    Parall
  !***
  !-----------------------------------------------------------------------
  use def_master
  use def_parame
  use def_parall
  use def_domain
  use mod_memchk
  use mod_parall, only : commd,PAR_COMM_MY_CODE
  use mod_parall, only : PAR_INTEGER
  use mod_parall, only : PAR_REAL
  use mod_parall, only : par_memor
  use def_mpi
  implicit none
  
  integer(ip), intent(in)  :: itask
  integer(4)               :: istat,npasi4,npasr4,npari4,nparr4
  integer(4)               :: dom_i
  real(rp)                 :: time1,time2

  call cputim(time1)

  call runend('PAR_SLAVES: OUT OF DATE SINCE 28/10/2014 FOR TASK '//intost(itask))

  if( ISLAVE ) then

     if( itask == 1_ip .or. itask == 4_ip ) then

        call runend('PAR_SLAVES 1-4: OUT OF DATE SINCE 28/10/2014')
        !-------------------------------------------------------------
        !
        ! INTEGER:                REAL:
        !
        ! Send:    PARIS(NPASI)   PARRS(NPASR)
        ! Receive: PARIN(NPARI)   PARRE(NPARR)
        ! 
        !-------------------------------------------------------------

#ifdef MPI_OFF
#else
        if( itask == 1_ip ) then
           dom_i  = int(kfl_desti_par,4)
        else
           dom_i  = int(commd % neights(kfl_desti_par),4)
        end if
        npasi4 = int(npasi,4)
        npari4 = int(npari,4)
        npasr4 = int(npasr,4)
        nparr4 = int(nparr,4)
        !
        ! Integer
        !
        if( npari /= 0 .and. npasi == 0 ) then
           call MPI_Recv(                          &
                parin(1:npari), npari4,            &
                PAR_INTEGER, dom_i, 0,             &
                PAR_COMM_MY_CODE, status, istat   )

        else if( npari == 0 .and. npasi /= 0 ) then
           call MPI_Send(                          &
                paris(1:npasi), npasi4,            &
                PAR_INTEGER, dom_i, 0,             &
                PAR_COMM_MY_CODE, istat           )

        else if( npari /= 0 .and. npasi /= 0 ) then
           call MPI_Sendrecv(                      &
                paris(1:npasi), npasi4,            &
                PAR_INTEGER,  dom_i, 0,            &
                parin(1:npari), npari4,            &
                PAR_INTEGER, dom_i, 0,             &
                PAR_COMM_MY_CODE, status, istat   )
        end if
        !
        ! Real
        !
        if( nparr /= 0 .and. npasr == 0 ) then
           call MPI_Recv(                          &
                parre(1:nparr), nparr4,            &
                PAR_REAL, dom_i, 0,    &
                PAR_COMM_MY_CODE, status, istat   )

        else if( nparr == 0 .and. npasr /= 0 ) then
           call MPI_Send(                          &
                parrs(1:npasr), npasr4,            &
                PAR_REAL, dom_i, 0,    &
                PAR_COMM_MY_CODE, istat           )

        else if( nparr /= 0 .and. npasr /= 0 ) then
           call MPI_Sendrecv(                      &
                parrs(1:), npasr4,                 &
                PAR_REAL,  dom_i, 0_4, &
                parre(1:), nparr4,                 &
                PAR_REAL, dom_i, 0,    &
                PAR_COMM_MY_CODE, status, istat   )
        end if

        npasi = 0_ip
        npari = 0_ip
        npasr = 0_ip
        nparr = 0_ip

#endif

     end if

  end if

  call cputim(time2)
  cpu_paral(25) = cpu_paral(25) + time2 - time1

end subroutine par_slaves
