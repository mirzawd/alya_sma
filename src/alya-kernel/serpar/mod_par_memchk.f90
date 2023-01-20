!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



module mod_par_memchk

  use def_parame
  use def_parall
  use def_kintyp_comm, only : tAdj_par
  use mod_memory,      only : lbytm
  use mod_memory,      only : memory_output_info
  !------------------------------------------------------------------------
  !****f* Parall/mod_par_memchk
  ! NAME
  !    mod_par_memchk
  ! DESCRIPTION
  !    Check some memory operations:
  !    - allocation;
  !    - deallocation; 
  !    - reallocation.
  ! OUTPUT
  ! USES
  !    memerr
  !    memctr
  ! USED BY
  !***
  !------------------------------------------------------------------------

  interface par_memchk
     module procedure mem_comm_data_par, &
          &           mem_tAdj_par
  end interface

contains

  subroutine mem_comm_data_par(itask,istat,memor,vanam,vacal,varia)
    !
    ! Type(comm_data_par)(:)
    !
    implicit none
    character(*), intent(in)    :: vanam,vacal
    integer(4),   intent(in)    :: istat
    integer(ip),  intent(in)    :: itask
    integer(8),   intent(inout) :: memor(2)
    type(comm_data_par)         :: varia(:)
    integer(ip)                 :: isize,nsize
    if(itask==0) then
       if(istat==0) then
          nsize=size(varia)
          lbytm=int(nsize*ip,KIND=8)
          do isize=1,nsize
             varia(isize)%bound_dim=0
             nullify(varia(isize)%neights) 
             nullify(varia(isize)%bound_size)
             nullify(varia(isize)%bound_perm)
          end do
       else
          call memerr(itask,vanam,vacal,istat)
       end if
    else
       lbytm=-size(varia,KIND=8)*int(ip,KIND=8)
    end if
    call memory_output_info(memor,vanam,vacal,'type')
  end subroutine mem_comm_data_par

  subroutine mem_tAdj_par(itask,istat,memor,vanam,vacal,varia)
    !
    ! Type(tAdj_par)(:)
    !
    implicit none
    character(*), intent(in)    :: vanam,vacal
    integer(4),   intent(in)    :: istat
    integer(ip),  intent(in)    :: itask
    integer(8),   intent(inout) :: memor(2)
    type(tAdj_par)              :: varia(:)
    integer(ip)                 :: isize,nsize
    if(itask==0) then
       if(istat==0) then
          nsize=size(varia)
          lbytm=int(nsize*ip,KIND=8)
          do isize=1,nsize
             varia(isize)%node1=0
             varia(isize)%node2=0
          end do
       else
          call memerr(itask,vanam,vacal,istat)
       end if
    else
       lbytm=-size(varia,KIND=8)*int(ip,KIND=8)
    end if
    call memory_output_info(memor,vanam,vacal,'type')
  end subroutine mem_tAdj_par

end module mod_par_memchk


