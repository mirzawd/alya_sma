!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine par_algath()
  !-------------------------------------------------------------------------------
  !****f* parall/par_algath
  ! NAME
  !    par_algath
  ! DESCRIPTION
  !    All gather
  ! INPUT
  ! OUTPUT
  ! USED BY
  !***
  !-------------------------------------------------------------------------------
  use def_master
  use def_domain
  use def_solver
  use def_parall  
  use mod_parall, only : PAR_INTEGER,PAR_COMM_MY_CODE
  use mod_memchk
  use def_mpi
  implicit none

  integer(4)           :: istat,nsize4
  integer(ip)          :: ii,kk,ifina,nsmall
  integer(4)           :: my_sendcount
  integer(ip), pointer :: my_recvbuf(:)
  integer(4),  pointer :: my_displs(:),my_recvcounts(:)
  integer(ip)          :: my_sendbuf(4)

  nsmall    = 4_ip
  my_sendcount = 4_4
  ifina     = 4*(npart_par+1)
  nsize4    = int(ifina,4)

  allocate( my_displs(npart_par+1) , stat = istat )
  allocate( my_recvcounts(npart_par+1) , stat = istat )
  allocate( my_recvbuf(ifina)       , stat = istat )

  do ii = 1,npart_par+1
     my_displs(ii)     = 4_4*(int(ii,4)-1_4)
     my_recvcounts(ii) = 4
     npoin_tot(ii)  = 0
     npoia_tot(ii)  = 0
     nelem_tot(ii)  = 0
     nboun_tot(ii)  = 0
  end do
  do ii = 1,ifina
     my_recvbuf(ii) = 0
  end do

  !----------------------------------------------------------------------
  !
  ! NPOIN: Number of interior+own nodes
  !
  !----------------------------------------------------------------------

  if( IMASTER ) then
     my_sendbuf(1) = 0
     my_sendbuf(2) = 0
     my_sendbuf(3) = 0
     my_sendbuf(4) = 0
  else
     my_sendbuf(1) = npoi1+(npoi3-npoi2)+1
     my_sendbuf(2) = npoin
     my_sendbuf(3) = nelem
     my_sendbuf(4) = nboun
  end if
  !
  ! MPI all gather
  !
#ifdef MPI_OFF
#else
  CALL MPI_ALLGATHERV(&
       my_sendbuf(1:4), my_sendcount, PAR_INTEGER,&
       my_recvbuf(1:ifina), my_recvcounts ,&
       my_displs(1:npart_par+1), PAR_INTEGER,&
       PAR_COMM_MY_CODE, istat)
#endif

  kk = 0
  do ii = 1,npart_par+1
     kk            = kk + 1
     npoin_tot(ii) = my_recvbuf(kk)
     kk            = kk + 1
     npoia_tot(ii) = my_recvbuf(kk)
     kk            = kk + 1
     nelem_tot(ii) = my_recvbuf(kk)
     kk            = kk + 1
     nboun_tot(ii) = my_recvbuf(kk)
  end do

  deallocate( my_recvcounts , stat = istat )
  deallocate( my_displs     , stat = istat )
  deallocate( my_recvbuf    , stat = istat )

end subroutine par_algath
