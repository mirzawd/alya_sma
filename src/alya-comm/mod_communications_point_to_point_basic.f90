!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @defgroup Communication_Toolbox
!> Toolbox for MPI communication, bridge to MPI
!> @{
!> @name    Parallelization toolbox
!> @file    mod_communications_point_to_point_basic.f90
!> @author  Guillaume Houzeaux
!> @date    28/06/2012
!> @brief   ToolBox for parallel communications
!> @details ToolBox for parallel communications
!------------------------------------------------------------------------

module mod_communications_point_to_point_basic

  use def_kintyp_basic,   only : ip,rp,lg
  use def_kintyp_comm,    only : comm_data_par_basic
  use def_communications
  use def_mpi
#include "def_mpi.inc"
  implicit none

  private

  integer(ip),            allocatable :: tmp_isend(:)
  integer(ip),            allocatable :: tmp_irecv(:)
  real(rp),               allocatable :: tmp_rsend(:)
  real(rp),               allocatable :: tmp_rrecv(:)  
  MY_MPI_REQUEST,         allocatable :: ireq4(:)

  interface PAR_INTERFACE_NODE_EXCHANGE
     module procedure PAR_INTERFACE_NODE_EXCHANGE_BASIC_IP_1,&
          &           PAR_INTERFACE_NODE_EXCHANGE_BASIC_RP_1
  end interface PAR_INTERFACE_NODE_EXCHANGE

  public :: PAR_INTERFACE_NODE_EXCHANGE
  
contains
  
  subroutine PAR_INTERFACE_NODE_EXCHANGE_BASIC_IP_1(xx,what,COMM,wsynch)
    
    integer(ip),                pointer,  intent(inout) :: xx(:)
    character(*),                         intent(in)    :: what
    type(comm_data_par_basic),            intent(in)    :: COMM
    character(*),               optional, intent(in)    :: wsynch
    integer(ip)                                         :: ndofn

    if( comm % RANK4 >= 0 .and. associated(xx) ) then
       ndofn = 1
       call PAR_INTERFACE_NODE_EXCHANGE_BASIC_IP(ndofn,xx,what,comm,wsynch)
    end if
    
  end subroutine PAR_INTERFACE_NODE_EXCHANGE_BASIC_IP_1  

  subroutine PAR_INTERFACE_NODE_EXCHANGE_BASIC_RP_1(xx,what,COMM,wsynch)
    
    real(rp),                   pointer,  intent(inout) :: xx(:)
    character(*),                         intent(in)    :: what
    type(comm_data_par_basic),            intent(in)    :: COMM
    character(*),               optional, intent(in)    :: wsynch
    integer(ip)                                         :: ndofn

    if( comm % RANK4 >= 0 .and. associated(xx) ) then
       ndofn = 1
       call PAR_INTERFACE_NODE_EXCHANGE_BASIC_RP(ndofn,xx,what,comm,wsynch)
    end if
    
  end subroutine PAR_INTERFACE_NODE_EXCHANGE_BASIC_RP_1

  subroutine PAR_INTERFACE_NODE_EXCHANGE_BASIC_IP(ndofn,xx,what,commu,wsynch)

    integer(ip),                    intent(in)    :: ndofn
    integer(ip),                    intent(inout) :: xx(ndofn,*)
    character(*),                   intent(in)    :: what
    type(comm_data_par_basic),      intent(in)    :: commu
    character(*),         optional, intent(in)    :: wsynch
    integer(ip)                                   :: ii,nsize,jj,dom_i
    integer(ip)                                   :: ipoin,ini,kk,idofn
    integer(4)                                    :: istat4,nsize4,count4
    integer(4)                                    :: dom_i4,my_rank4
    logical(lg)                                   :: asynch
    integer(ip),                       save       :: ipass = 0
    MY_MPI_STATUS,                     pointer    :: status4(:)
    integer(ip)                                   :: bound_dim
    integer(ip),                       pointer    :: bound_perm(:)
    integer(ip),                       pointer    :: bound_invp(:)
    integer(ip),                       pointer    :: bound_size(:)
    MY_MPI_COMM                                   :: PAR_COMM_TO_USE

    PAR_COMM_TO_USE = commu % PAR_COMM_WORLD
    my_rank4        = commu % RANK4
 
#ifndef MPI_OFF

    if( my_rank4 >= 0 .and. associated(commu % bound_perm) ) then
       !
       ! Passes
       !
       ipass = ipass + 1
       !
       ! On nodes or edges
       ! 
       bound_dim  =  size(commu % bound_perm)
       bound_perm => commu % bound_perm
       bound_invp => commu % bound_perm
       bound_size => commu % bound_size
       !
       ! Synchronous or asynchronous
       !
       if( present(wsynch) ) then
          if( trim(wsynch) == 'SYNCHRONOUS' .or. trim(wsynch) == 'BLOCKING' ) then
             asynch = .false.
          else if( trim(wsynch) == 'ASYNCHRONOUS' .or. trim(wsynch) == 'NON BLOCKING' ) then
             asynch = .true.
          else
             call runend('PAR_NODE_ASSMEMBLY: UNKNOWN COMMUNICATION TYPE')
          end if
       else
          asynch = .false.
       end if

       if( ipass == 1 ) then
          !
          ! Allocate memory
          !
          if( asynch ) allocate(ireq4(commu % nneig*2))
          allocate(tmp_isend(bound_dim * ndofn))
          allocate(tmp_irecv(bound_dim * ndofn))

          kk = 0
          do jj = 1,bound_dim
             ipoin = bound_perm(jj)
             do idofn = 1,ndofn
                kk = kk + 1
                tmp_isend(kk) = xx(idofn,ipoin)
                tmp_irecv(kk) = 0
             end do
          end do
          !
          ! Send    temp_send
          ! Receive temp_recv
          !
          istat4 = 0_4
          kk = 0
          do ii = 1,commu % nneig

             dom_i  = commu % neights(ii)
             dom_i4 = int(dom_i,4)

             ini   = ndofn * ( bound_size(ii)   - 1 ) + 1
             nsize = ndofn * ( bound_size(ii+1) - 1 ) + 1 - ini
             nsize4 = int(nsize,4)

#if MPI_VERSION==3

#else
             if( asynch ) then
                kk = kk + 1
                call MPI_Isend(&
                     tmp_isend(ini:ini+nsize-1), nsize4, &
                     PAR_INTEGER,  dom_i4, 0_4,          &
                     PAR_COMM_TO_USE, ireq4(kk), istat4 )
                kk = kk + 1
                call MPI_Irecv(&
                     tmp_irecv(ini:ini+nsize-1), nsize4, &
                     PAR_INTEGER,  dom_i4, 0_4,          &
                     PAR_COMM_TO_USE, ireq4(kk), istat4 )
             else
                call MPI_Sendrecv(                       &
                     tmp_isend(ini:), nsize4,            &
                     PAR_INTEGER, dom_i4, 0_4,           &
                     tmp_irecv(ini:), nsize4,            &
                     PAR_INTEGER, dom_i4, 0_4,           &
                     PAR_COMM_TO_USE, status, istat4    )
             end if
#endif
             
             if( istat4 /= 0_4 ) call runend('PAR_INTERFACE_NODE_EXCHANGE_IP: MPI ERROR')

          end do

       end if
       !
       ! sum,max,min on temp_recv
       !
       if( asynch .and. ipass == 2 ) then
          count4 = 2*int(commu % nneig,4)
          allocate( status4(MPI_STATUS_SIZE*2*commu % nneig) )
          CALL MPI_WAITALL(count4,ireq4,status4,istat4)
          deallocate( status4 )
          deallocate(ireq4)
       end if

       if( ( asynch .and. ipass == 2 ) .or. ( .not. asynch .and. ipass == 1 ) ) then

          if( trim(what) == 'SUM' .or. trim(what) == 'ASSEMBLY' .or. trim(what) == 'DISTRIBUTE' ) then
             !
             ! SUM
             !
             kk = 0
             do jj = 1,bound_dim
                ipoin = bound_invp(jj)
                do idofn = 1,ndofn
                   kk = kk + 1
                   xx(idofn,ipoin) = xx(idofn,ipoin) + tmp_irecv(kk)
                end do
             end do

          else if( trim(what) == 'MAX' .or. trim(what) == 'MAXIMUM' ) then
             !
             ! MAX
             !
             kk = 0
             do jj = 1,bound_dim
                ipoin = bound_invp(jj)
                do idofn = 1,ndofn
                   kk = kk + 1
                   xx(idofn,ipoin) = max(xx(idofn,ipoin),tmp_irecv(kk))
                end do
             end do

          else if( trim(what) == 'MIN' .or. trim(what) == 'MINIMUM' ) then
             !
             ! MIN
             !
             kk = 0
             do jj = 1,bound_dim
                ipoin = bound_invp(jj)
                do idofn = 1,ndofn
                   kk = kk + 1
                   xx(idofn,ipoin) = min(xx(idofn,ipoin),tmp_irecv(kk))
                end do
             end do

          else if( trim(what) == 'TAKE MIN' ) then
             !
             ! Take value of minimum neighbor rank
             !
             kk = 0
             do ii = 1,commu % nneig
                dom_i = commu % neights(ii)
                if( dom_i < my_rank4 ) then
                   ini   = ndofn * ( bound_size(ii)   - 1 ) + 1
                   nsize = ndofn * ( bound_size(ii+1) - 1 ) + 1 - ini
                   do jj = bound_size(ii),bound_size(ii+1)-1
                      ipoin = bound_perm(jj)
                      do idofn = 1,ndofn
                         kk = kk + 1
                         xx(idofn,ipoin) = tmp_irecv(kk)
                      end do
                   end do
                else
                   kk = kk + (bound_size(ii+1)-bound_size(ii))*ndofn
                end if
             end do 

          else

             call runend('UNKNOWN ORDER: '//trim(what))

          end if

          ipass = 0
          deallocate(tmp_irecv)
          deallocate(tmp_isend)

       end if

    end if
#endif

  end subroutine PAR_INTERFACE_NODE_EXCHANGE_BASIC_IP

  subroutine PAR_INTERFACE_NODE_EXCHANGE_BASIC_RP(ndofn,xx,what,commu,wsynch)

    integer(ip),                    intent(in)    :: ndofn
    real(rp),                       intent(inout) :: xx(ndofn,*)
    character(*),                   intent(in)    :: what
    type(comm_data_par_basic),      intent(in)    :: commu
    character(*),         optional, intent(in)    :: wsynch
    integer(ip)                                   :: ii,nsize,jj,dom_i
    integer(ip)                                   :: ipoin,ini,kk,idofn
    integer(4)                                    :: istat4,nsize4,count4
    integer(4)                                    :: dom_i4,my_rank4
    logical(lg)                                   :: asynch
    integer(ip),                       save       :: ipass = 0
    MY_MPI_STATUS,                     pointer    :: status4(:)
    integer(ip)                                   :: bound_dim
    integer(ip),                       pointer    :: bound_perm(:)
    integer(ip),                       pointer    :: bound_invp(:)
    integer(ip),                       pointer    :: bound_size(:)
    MY_MPI_COMM                                   :: PAR_COMM_TO_USE

    PAR_COMM_TO_USE = commu % PAR_COMM_WORLD
    my_rank4        = commu % RANK4

#ifndef MPI_OFF

    if( my_rank4 >= 0 .and. associated(commu % bound_perm) ) then
       !
       ! Passes
       !
       ipass = ipass + 1
       !
       ! On nodes or edges
       ! 
       bound_dim  =  size(commu % bound_perm)
       bound_perm => commu % bound_perm
       bound_invp => commu % bound_perm
       bound_size => commu % bound_size
       !
       ! Synchronous or asynchronous
       !
       if( present(wsynch) ) then
          if( trim(wsynch) == 'SYNCHRONOUS' .or. trim(wsynch) == 'BLOCKING' ) then
             asynch = .false.
          else if( trim(wsynch) == 'ASYNCHRONOUS' .or. trim(wsynch) == 'NON BLOCKING' ) then
             asynch = .true.
          else
             call runend('PAR_NODE_ASSMEMBLY: UNKNOWN COMMUNICATION TYPE')
          end if
       else
          asynch = .false.
       end if

       if( ipass == 1 ) then
          !
          ! Allocate memory
          !
          if( asynch ) allocate(ireq4(commu % nneig*2))
          allocate(tmp_rsend(bound_dim * ndofn))
          allocate(tmp_rrecv(bound_dim * ndofn))

          kk = 0
          do jj = 1,bound_dim
             ipoin = bound_perm(jj)
             do idofn = 1,ndofn
                kk = kk + 1
                tmp_rsend(kk) = xx(idofn,ipoin)
                tmp_rrecv(kk) = 0
             end do
          end do
          !
          ! Send    temp_send
          ! Receive temp_recv
          !
          istat4 = 0_4
          kk = 0
          do ii = 1,commu % nneig

             dom_i  = commu % neights(ii)
             dom_i4 = int(dom_i,4)

             ini   = ndofn * ( bound_size(ii)   - 1 ) + 1
             nsize = ndofn * ( bound_size(ii+1) - 1 ) + 1 - ini
             nsize4 = int(nsize,4)
             
             if( asynch ) then
                kk = kk + 1
                call MPI_Isend(&
                     tmp_rsend(ini:ini+nsize-1), nsize4, &
                     PAR_REAL, dom_i4, 0_4,  &
                     PAR_COMM_TO_USE, ireq4(kk), istat4  )
                kk = kk + 1
                call MPI_Irecv(&
                     tmp_rrecv(ini:ini+nsize-1), nsize4, &
                     PAR_REAL,  dom_i4, 0_4, &
                     PAR_COMM_TO_USE, ireq4(kk), istat4  )
             else
                call MPI_Sendrecv(                       &
                     tmp_rsend(ini:), nsize4,            &
                     PAR_REAL, dom_i4, 0_4,  &
                     tmp_rrecv(ini:), nsize4,            &
                     PAR_REAL, dom_i4, 0_4,  &
                     PAR_COMM_TO_USE, status, istat4     )
             end if
             if( istat4 /= 0_4 ) call runend('PAR_INTERFACE_NODE_EXCHANGE_IP: MPI ERROR')

          end do

       end if
       !
       ! sum,max,min on temp_recv
       !
       if( asynch .and. ipass == 2 ) then
          count4 = 2*int(commu % nneig,4)
          allocate( status4(MPI_STATUS_SIZE*2*commu % nneig) )
          CALL MPI_WAITALL(count4,ireq4,status4,istat4)
          !CALL MPI_WAITALL(ireq4,status4,istat4)
          deallocate( status4 )
          deallocate(ireq4)
       end if

       if( ( asynch .and. ipass == 2 ) .or. ( .not. asynch .and. ipass == 1 ) ) then

          if( trim(what) == 'SUM' .or. trim(what) == 'ASSEMBLY' .or. trim(what) == 'DISTRIBUTE' ) then
             !
             ! SUM
             !
             kk = 0
             do jj = 1,bound_dim
                ipoin = bound_invp(jj)
                do idofn = 1,ndofn
                   kk = kk + 1
                   xx(idofn,ipoin) = xx(idofn,ipoin) + tmp_rrecv(kk)
                end do
             end do

          else if( trim(what) == 'MAX' .or. trim(what) == 'MAXIMUM' ) then
             !
             ! MAX
             !
             kk = 0
             do jj = 1,bound_dim
                ipoin = bound_invp(jj)
                do idofn = 1,ndofn
                   kk = kk + 1
                   xx(idofn,ipoin) = max(xx(idofn,ipoin),tmp_rrecv(kk))
                end do
             end do

          else if( trim(what) == 'MIN' .or. trim(what) == 'MINIMUM' ) then
             !
             ! MIN
             !
             kk = 0
             do jj = 1,bound_dim
                ipoin = bound_invp(jj)
                do idofn = 1,ndofn
                   kk = kk + 1
                   xx(idofn,ipoin) = min(xx(idofn,ipoin),tmp_rrecv(kk))
                end do
             end do
             
          else if( trim(what) == 'TAKE MIN' ) then
             !
             ! Take value of minimum neighbor rank
             !
             kk = 0
             do ii = 1,commu % nneig
                dom_i = commu % neights(ii)
                if( dom_i < my_rank4 ) then
                   ini   = ndofn * ( bound_size(ii)   - 1 ) + 1
                   nsize = ndofn * ( bound_size(ii+1) - 1 ) + 1 - ini
                   do jj = bound_size(ii),bound_size(ii+1)-1
                      ipoin = bound_perm(jj)
                      do idofn = 1,ndofn
                         kk = kk + 1
                         xx(idofn,ipoin) = tmp_rrecv(kk)
                      end do
                   end do
                else
                   kk = kk + (bound_size(ii+1)-bound_size(ii))*ndofn
                end if
             end do
             
          else

             call runend('UNKNOWN ORDER: '//trim(what))

          end if

          ipass = 0
          deallocate(tmp_rrecv)
          deallocate(tmp_rsend)

       end if

    end if
#endif

  end subroutine PAR_INTERFACE_NODE_EXCHANGE_BASIC_RP

end module mod_communications_point_to_point_basic
!> @}

