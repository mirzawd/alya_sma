!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine par_lagran(itask)
  !-----------------------------------------------------------------------
  !****f* Parall/par_lagran
  ! NAME
  !    par_lagran
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
  use mod_parall, only : PAR_INTEGER,commd
  use mod_parall, only : PAR_COMM_MY_CODE
  use mod_parall, only : par_memor
  use mod_memory, only : memory_size
  use mod_memory, only : memory_alloca
  use mod_memory, only : memory_deallo
  use def_mpi
#include "def_mpi.inc"
  implicit none

  integer(ip), intent(in)  :: itask
  integer(ip)              :: ii,jj
  integer(4)               :: bsize4,dom_i
  integer(4)               :: istat
  integer(ip), pointer     :: loc_rpari1(:)
  
  if( IPARALL ) then

     nullify(loc_rpari1)
     
     select case ( itask )

     case ( 3_ip )
        !
        ! Check repeated particles
        !
        if( ISLAVE ) then
           call memory_alloca(par_memor,'LOC_RPARI1','par_lagran',loc_rpari1,npari)
           pard1 = 0

           do ii = 1, nneig

              dom_i = int(commd % neights(ii),4)            

#ifdef MPI_OFF
#else
              bsize4 = int(npari,4)
              call MPI_Sendrecv(                 &
                   pari1(1:),  bsize4,           &
                   PAR_INTEGER,  dom_i, 0_4,     &
                   loc_rpari1(1:), bsize4,       &
                   PAR_INTEGER, dom_i, 0_4,      &
                   PAR_COMM_MY_CODE, status, istat )
#endif

              do jj = 1,npari
                 if( pari1(jj) > 0 .and. loc_rpari1(jj) > 0 ) then
                    if( kfl_paral > dom_i ) then
                       pard1 = 1
                       pari1(jj) = -abs(pari1(jj))
                    end if
                 end if
              end do

           end do

           call memory_deallo(par_memor,'LOC_RPARI1','par_lagran',loc_rpari1)

        end if

     case default

        call runend('PAR_ALLGAT: DOES NOT EXIST')
        
      end select

  end if

end subroutine par_lagran
