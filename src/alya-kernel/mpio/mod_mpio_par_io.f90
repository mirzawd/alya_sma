!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup mpio
!> @{
!> @file    mod_mpio_par_io.f90
!> @author  Damien Dosimont
!> @date    27/09/2017
!> @brief   Parallel I/O interface
!> @details This module offers a high level interface to I/O parallel operations
!>          It is responsible for managing MPI-IO file format metadata and file operations
!> @}
!-----------------------------------------------------------------------

module mod_mpio_par_io

    use def_master,             only : IIOSLAVE, IMASTER, INOTMASTER, IIOSLAVE
    use def_master,             only : namda
    use def_master,             only : intost, retost
    use def_master,             only : optional_argument
    use def_kintyp,             only : ip, rp, lg
    use mod_parall,             only : PAR_COMM_WORLD
    use mod_parall,             only : PAR_CODE_SIZE, PAR_MY_CODE_RANK, PAR_COMM_MPIO, PAR_COMM_MPIO_RANK_WM, PAR_COMM_MPIO_WM, PAR_COMM_MPIO_WM_SIZE
    use mod_communications,     only : PAR_COMM_SPLIT, PAR_GATHERV, PAR_GATHER, PAR_SEND_RECEIVE
    use mod_communications,     only : PAR_BROADCAST, PAR_BARRIER
    use mod_memory,             only : memory_alloca, memory_deallo
    use mod_mpio_seq_io
    use mod_mpio_par_mpiwrapper
    use def_mpio
    use mod_mpio_par_log
    use mod_mpio_par_async_io
    use mod_mpio_config,        only : mpio_config
    use def_mpi
#include "def_mpi.inc"
    
    implicit none

    integer(ip)                             :: i, j

    private

    character(150)                          ::  wherein_code="IN MY CODE"
    character(150)                          ::  wherein_p="IN MPIO WITH MASTER"

#ifndef MPI_OFF
!    integer                                 ::  SMALL_FILE_SIZE = 104857600 ! 100 MB
#endif


    interface PAR_FILE_WRITE_ALL
        module procedure                        PAR_FILE_WRITE_ALL_INT_V,&
                                                PAR_FILE_WRITE_ALL_INT_M,&
                                                PAR_FILE_WRITE_ALL_REAL_V,&
                                                PAR_FILE_WRITE_ALL_REAL_M
    end interface

    interface PAR_FILE_READ_ALL
        module procedure                        PAR_FILE_READ_ALL_INT_V,&
                                                PAR_FILE_READ_ALL_INT_M,&
                                                PAR_FILE_READ_ALL_REAL_V,&
                                                PAR_FILE_READ_ALL_REAL_M
    end interface

    public                                  ::  PAR_FILE_WRITE_ALL, PAR_FILE_READ_ALL, PAR_FILE_READ_HEADER, PAR_FILE_OPTION_ENABLED, par_broadcast_header
    public                                  ::  PAR_FILE_WRITE_HEADER
    
    contains

    subroutine par_broadcast_header(header)
        type(mpio_header), intent(inout)               :: header
#ifndef MPI_OFF
        call PAR_BROADCAST(header%magic_number, wherein_code)
        call PAR_BROADCAST(string_count, header%format, wherein_code)
        call PAR_BROADCAST(string_count, header%version, wherein_code)
        call PAR_BROADCAST(string_count, header%object, wherein_code)
        call PAR_BROADCAST(string_count, header%dimension, wherein_code)
        call PAR_BROADCAST(string_count, header%resultson, wherein_code)
        call PAR_BROADCAST(string_count, header%type, wherein_code)
        call PAR_BROADCAST(string_count, header%size, wherein_code)
        call PAR_BROADCAST(string_count, header%par, wherein_code)
        call PAR_BROADCAST(string_count, header%filter, wherein_code)
        call PAR_BROADCAST(string_count, header%sorting, wherein_code)
        call PAR_BROADCAST(string_count, header%id, wherein_code)
        call PAR_BROADCAST(string_count, header%align_chars, wherein_code)
        call PAR_BROADCAST(header%columns, wherein_code)
        call PAR_BROADCAST(header%lines, wherein_code)
        call PAR_BROADCAST(header%ittim, wherein_code)
        call PAR_BROADCAST(header%nsubd, wherein_code)
        call PAR_BROADCAST(header%divi, wherein_code)
        call PAR_BROADCAST(header%tag1, wherein_code)
        call PAR_BROADCAST(header%tag2, wherein_code)
        call PAR_BROADCAST(header%time, wherein_code)
        call PAR_BROADCAST(string_count, header%align_chars, wherein_code)
        do i=1, option_size
            call PAR_BROADCAST(string_count, header%options%opt(i), wherein_code)
        end do
        call PAR_BROADCAST(header%item_size, wherein_code)
        call PAR_BROADCAST(header%file_size, wherein_code)
#endif
    end subroutine

    subroutine par_check_header(header)
        type(mpio_header), intent(inout)               :: header
#ifndef MPI_OFF
        if (IMASTER) then
            call check_header(header)
        end if
#endif
    end subroutine

    subroutine par_fill_header(header, divi, tag1, tag2, ittim, cutim, opt)
        type(mpio_header), intent(inout)               :: header
        integer(ip),           intent(in)              :: divi
        integer(ip),           intent(in)              :: tag1
        integer(ip),           intent(in)              :: tag2
        integer(ip),            intent(in), optional    :: ittim
        real(rp),               intent(in), optional    :: cutim
        type(mpio_header_options), intent(in), optional:: opt
#ifndef MPI_OFF
        call FILL_HEADER(header, divi, tag1, tag2, ittim, cutim, opt)
#endif
    end subroutine

    subroutine PAR_FILE_READ_HEADER(filename, header, wherein)
        character(*)                                        ::  filename
        type(mpio_header), intent(inout)                    ::  header
        MY_MPI_FILE                                         ::  fh
        character(150)   , intent(in)                       ::  wherein
#ifndef MPI_OFF
        call PAR_FILE_OPEN_READ(fh, filename, wherein)
        if (IMASTER) then
            call PAR_FILE_READ(fh, header%magic_number)
            call PAR_FILE_READ(fh, header%format)
            call PAR_FILE_READ(fh, header%version)
            call PAR_FILE_READ(fh, header%object)
            call PAR_FILE_READ(fh, header%dimension)
            call PAR_FILE_READ(fh, header%resultson)
            call PAR_FILE_READ(fh, header%type)
            call PAR_FILE_READ(fh, header%size)
            call PAR_FILE_READ(fh, header%par)
            call PAR_FILE_READ(fh, header%filter)
            call PAR_FILE_READ(fh, header%sorting)
            call PAR_FILE_READ(fh, header%id)
            call PAR_FILE_READ(fh, header%align_chars)
            call PAR_FILE_READ(fh, header%columns)
            call PAR_FILE_READ(fh, header%lines)
            call PAR_FILE_READ(fh, header%ittim)
            call PAR_FILE_READ(fh, header%nsubd)
            call PAR_FILE_READ(fh, header%divi)
            call PAR_FILE_READ(fh, header%tag1)
            call PAR_FILE_READ(fh, header%tag2)
            call PAR_FILE_READ(fh, header%time)
            call PAR_FILE_READ(fh, header%align_chars)
            call PAR_FILE_READ(fh, header%options%opt, option_size)
            header%item_size=0_ip
            if (header%size(1:7)=='8BYTE00') then
               header%item_size=8
            elseif (header%size(1:7)=='4BYTE00') then
               header%item_size=4_ip
            end if
            header%file_size=int(header_size,8)+(int(header%lines,8)*int(header%columns,8)*int(header%item_size,8))
        end if
        call par_broadcast_header(header)
        call par_check_header(header)
        call PAR_FILE_CLOSE(fh)
        call PAR_BARRIER(wherein)
#endif
    end subroutine

    function PAR_FILE_OPTION_ENABLED(header, opt) result(enabled)
        type(mpio_header), intent(in)   ::  header
        character(7), intent(in)        ::  opt
        logical                         ::  enabled
        enabled=file_option_enabled(header, opt)
    end function

    subroutine PAR_FILE_WRITE_HEADER(fh, header,rank4)
        MY_MPI_FILE   ,       intent(in) :: fh
        type(mpio_header),    intent(in) :: header
        integer(4), optional, intent(in) :: rank4
#ifndef MPI_OFF

        if( IMASTER .or. optional_argument(1_4,rank4) == 0_4 ) then
            call PAR_FILE_WRITE(fh, header % magic_number)
            call PAR_FILE_WRITE(fh, header % format)
            call PAR_FILE_WRITE(fh, header % version)
            call PAR_FILE_WRITE(fh, header % object)
            call PAR_FILE_WRITE(fh, header % dimension)
            call PAR_FILE_WRITE(fh, header % resultson)
            call PAR_FILE_WRITE(fh, header % type)
            call PAR_FILE_WRITE(fh, header % size)
            call PAR_FILE_WRITE(fh, header % par)
            call PAR_FILE_WRITE(fh, header % filter)
            call PAR_FILE_WRITE(fh, header % sorting)
            call PAR_FILE_WRITE(fh, header % id)
            call PAR_FILE_WRITE(fh, header % align_chars)
            call PAR_FILE_WRITE(fh, header % columns)
            call PAR_FILE_WRITE(fh, header % lines)
            call PAR_FILE_WRITE(fh, header % ittim)
            call PAR_FILE_WRITE(fh, header % nsubd)
            call PAR_FILE_WRITE(fh, header % divi)
            call PAR_FILE_WRITE(fh, header % tag1)
            call PAR_FILE_WRITE(fh, header % tag2)
            call PAR_FILE_WRITE(fh, header % time)
            call PAR_FILE_WRITE(fh, header % align_chars)
            call PAR_FILE_WRITE(fh, header % options%opt, option_size)
        end if
        if( .not. present(rank4) ) call PAR_BARRIER(wherein_code)
#endif
    end subroutine

    subroutine PAR_FILE_WRITE_ALL_INT_V(buf, filename, object, resultson, dim, divi, tag1, tag2, ittim, cutim, opt, nsubd)
      integer(ip), pointer,  intent(in)              :: buf(:)
      character(*),          intent(in)              :: filename
      character(8),          intent(in)              :: object
      character(8),          intent(in)              :: resultson
      integer(ip),           intent(in)              :: dim
      integer(ip),           intent(in)              :: divi
      integer(ip),           intent(in)              :: tag1
      integer(ip),           intent(in)              :: tag2
      integer(ip),           intent(in), optional    :: ittim
      real(rp),              intent(in), optional    :: cutim
      type(mpio_header_options), intent(in), optional:: opt
      integer(ip),           intent(in), optional    :: nsubd
      integer(ip), pointer                           :: vec(:)
      integer(ip)                                    :: myoffset
      MY_MPI_FILE                                    :: fh
      integer(8)                                     :: offset                              ! Current file offset
      type(mpio_header)                              :: header
      character(100), PARAMETER :: vacal = "PAR_FILE_WRITE_INT_V"
      nullify(vec)
#ifndef MPI_OFF
      header%columns=1_header_ip
      if (ip==4) then
         header%size='4BYTE00'//char(0)
      else
         header%size='8BYTE00'//char(0)
      end if
      header%type='INTEG00'//char(0)
      call memory_alloca(mpio_memor,'vec',vacal,vec, PAR_CODE_SIZE)
      call PAR_GATHER(dim, vec, wherein_code)
      call PAR_BROADCAST(vec, wherein_code)
      vec(1)=0_ip
      header%lines=int(sum(vec), header_ip)
      myoffset=0_ip
      if (INOTMASTER) then
         myoffset=sum(vec(2:PAR_MY_CODE_RANK))
      end if
      header%object=object
      header%resultson=resultson
      if (present(nsubd)) then
         header%nsubd=int(nsubd, header_ip)
      else
         header%nsubd=1_header_ip 
      end if

      call par_fill_header(header, divi, tag1, tag2, ittim, cutim, opt)
      call PAR_INFO_CREATE()
      call par_start_timer(wherein_code, header%object(1:5), 'WRITE', 'MPIAL', PAR_CODE_SIZE-1, PAR_CODE_SIZE, header%file_size)
      call PAR_FILE_OPEN_WRITE(fh, filename)
      call PAR_INFO_FREE()
      call PAR_FILE_SET_SIZE(fh, header%file_size)
      call PAR_BARRIER(wherein_code)
      call PAR_FILE_WRITE_HEADER(fh, header)
      call PAR_FILE_SET_VIEW(fh, buf, int(header_size,8))
      offset = int(myoffset,8)
      call PAR_FILE_SEEK_SET(fh, offset)
      if (.not. mpio_config%output%asynchronous) then
         call PAR_FILE_WRITE_ALL(fh, buf, dim)
         call PAR_FILE_CLOSE(fh)
      else
         call PAR_ASYNCHRONOUS_WRITE_POP_ALL(mpio_config%output%flush, ittim)
         call PAR_ASYNCHRONOUS_WRITE_PUSH(fh, buf, dim)
      end if
      call par_end_timer(wherein_code)
      call memory_deallo(mpio_memor,'vec',vacal,vec)
      call PAR_BARRIER(wherein_code)
#endif
    end subroutine PAR_FILE_WRITE_ALL_INT_V

    subroutine PAR_FILE_WRITE_ALL_INT_M(buf, filename, object, resultson, dim1, dim2, divi, tag1, tag2, ittim, cutim, opt, nsubd)
        integer(ip), pointer,   intent(in)             :: buf(:,:)
        character(*),          intent(in)              :: filename
        character(8),          intent(in)              :: object
        character(8),          intent(in)              :: resultson
        integer(ip),           intent(in)              :: dim1
        integer(ip),           intent(in)              :: dim2
        integer(ip),           intent(in)              :: divi
        integer(ip),           intent(in)              :: tag1
        integer(ip),           intent(in)              :: tag2
        integer(ip),           intent(in), optional    :: ittim
        real(rp),              intent(in), optional    :: cutim
        type(mpio_header_options), intent(in), optional:: opt
        integer(ip),           intent(in), optional    :: nsubd
        integer(ip), pointer                           :: vec(:)
        integer(ip)                                    :: myoffset
        MY_MPI_FILE                                    :: fh
        integer(8)                                     :: offset                              ! Current file offset
        type(mpio_header)                              :: header
        character(100), PARAMETER :: vacal = "PAR_FILE_WRITE_INT_M"
        nullify(vec)
#ifndef MPI_OFF
        ! specific to this subroutine
        header%columns=int(dim1, header_ip)
        if (ip==4) then
            header%size='4BYTE00'//char(0)
        else
            header%size='8BYTE00'//char(0)
        end if
        header%type='INTEG00'//char(0)
        call memory_alloca(mpio_memor,'vec',vacal,vec, PAR_CODE_SIZE)
        call PAR_GATHER(dim2, vec, wherein_code)
        call PAR_BROADCAST(vec, wherein_code)
        vec(1)=0_ip
        header%lines=int(sum(vec), header_ip)
        myoffset=0_ip
        if (INOTMASTER) then
            myoffset=sum(vec(2:PAR_MY_CODE_RANK))
        end if
        header%object=object
        header%resultson=resultson
        if (present(nsubd)) then
            header%nsubd=int(nsubd, header_ip)
        else
            header%nsubd=1_header_ip
        end if
        call par_fill_header(header, divi, tag1, tag2, ittim, cutim, opt)
        call PAR_INFO_CREATE()
        call par_start_timer(wherein_code, header%object(1:5), 'WRITE', 'MPIAL', PAR_CODE_SIZE-1, PAR_CODE_SIZE, header%file_size)
        call PAR_FILE_OPEN_WRITE(fh, filename)
        call PAR_INFO_FREE()
        call PAR_FILE_SET_SIZE(fh, header%file_size)
        call PAR_BARRIER(wherein_code)
        call PAR_FILE_WRITE_HEADER(fh, header)
        call PAR_FILE_SET_VIEW(fh, buf, int(header_size,8))
        offset = int(myoffset,8) * int(header%columns,8)
        call PAR_FILE_SEEK_SET(fh, offset)
        if (.not. mpio_config%output%asynchronous) then
            call PAR_FILE_WRITE_ALL(fh, buf, int(header%columns,ip)*dim2)
            call PAR_FILE_CLOSE(fh)
        else
            call PAR_ASYNCHRONOUS_WRITE_POP_ALL(mpio_config%output%flush, ittim)
            call PAR_ASYNCHRONOUS_WRITE_PUSH(fh, buf, int(header%columns, ip), dim2)
        end if
        call par_end_timer(wherein_code)
        call memory_deallo(mpio_memor,'vec',vacal,vec)
        call PAR_BARRIER(wherein_code)
#endif
    end subroutine

    subroutine PAR_FILE_WRITE_ALL_REAL_V(buf, filename, object, resultson, dim, divi, tag1, tag2, ittim, cutim, opt, nsubd)
        real(rp), pointer,     intent(in)              :: buf(:)
        character(*),          intent(in)              :: filename
        character(8),          intent(in)              :: object
        character(8),          intent(in)              :: resultson
        integer(ip),           intent(in)              :: dim
        integer(ip),           intent(in)              :: divi
        integer(ip),           intent(in)              :: tag1
        integer(ip),           intent(in)              :: tag2
        integer(ip),           intent(in), optional    :: ittim
        real(rp),              intent(in), optional    :: cutim
        type(mpio_header_options), intent(in), optional:: opt
        integer(ip),           intent(in), optional    :: nsubd
        integer(ip), pointer                           :: vec(:)
        integer(ip)                                    :: myoffset
        MY_MPI_FILE                                    :: fh
        integer(8)                                     :: offset                              ! Current file offset
        type(mpio_header)                              :: header
        character(100), PARAMETER :: vacal = "PAR_FILE_WRITE_REAL_V"
        nullify(vec)
#ifndef MPI_OFF
        ! specific to this subroutine
        header%columns=1_ip
        header%size='8BYTE00'//char(0)
        header%type='REAL000'//char(0)
        call memory_alloca(mpio_memor,'vec',vacal,vec, PAR_CODE_SIZE)
        call PAR_GATHER(dim, vec, wherein_code)
        call PAR_BROADCAST(vec, wherein_code)
        vec(1)=0_ip
        header%lines=int(sum(vec), header_ip)
        myoffset=0_ip
        if (INOTMASTER) then
            myoffset=sum(vec(2:PAR_MY_CODE_RANK))
        end if
        header%object=object
        header%resultson=resultson
        if (present(nsubd)) then
            header%nsubd=int(nsubd, header_ip)
        else
            header%nsubd=1_header_ip
        end if
        call par_fill_header(header, divi, tag1, tag2, ittim, cutim, opt)
        call PAR_INFO_CREATE()
        call par_start_timer(wherein_code, header%object(1:5), 'WRITE', 'MPIAL', PAR_CODE_SIZE-1, PAR_CODE_SIZE, header%file_size)
        call PAR_FILE_OPEN_WRITE(fh, filename)
        call PAR_INFO_FREE()
        call PAR_FILE_SET_SIZE(fh, header%file_size)
        call PAR_BARRIER(wherein_code)
        call PAR_FILE_WRITE_HEADER(fh, header)
        call PAR_FILE_SET_VIEW(fh, buf, int(header_size,8))
        offset = int(myoffset,8)
        call PAR_FILE_SEEK_SET(fh, offset)
        if (.not. mpio_config%output%asynchronous) then
            call PAR_FILE_WRITE_ALL(fh, buf, dim)
            call PAR_FILE_CLOSE(fh)
        else
            call PAR_ASYNCHRONOUS_WRITE_POP_ALL(mpio_config%output%flush, ittim)
            call PAR_ASYNCHRONOUS_WRITE_PUSH(fh, buf, dim)
        end if
        call par_end_timer(wherein_code)
        call memory_deallo(mpio_memor,'vec',vacal,vec)
        call PAR_BARRIER(wherein_code)
#endif
    end subroutine

    subroutine PAR_FILE_WRITE_ALL_REAL_M(buf, filename, object, resultson, dim1, dim2, divi, tag1, tag2, ittim, cutim, opt, nsubd)
        real(rp), pointer,      intent(in)             :: buf(:,:)
        character(*),          intent(in)              :: filename
        character(8),          intent(in)              :: object
        character(8),          intent(in)              :: resultson
        integer(ip),           intent(in)              :: dim1
        integer(ip),           intent(in)              :: dim2
        integer(ip),           intent(in)              :: divi
        integer(ip),           intent(in)              :: tag1
        integer(ip),           intent(in)              :: tag2
        integer(ip),           intent(in), optional    :: ittim
        real(rp),              intent(in), optional    :: cutim
        type(mpio_header_options), intent(in), optional:: opt
        integer(ip),           intent(in), optional    :: nsubd
        integer(ip), pointer                           :: vec(:)
        integer(ip)                                    :: myoffset
        MY_MPI_FILE                                    :: fh
        integer(8)                                     :: offset                              ! Current file offset
        type(mpio_header)                              :: header
        character(100), PARAMETER :: vacal = "PAR_FILE_WRITE_REAL_M"
        nullify(vec)
#ifndef MPI_OFF
        ! specific to this subroutine
        header%columns = int(dim1, header_ip)
        header%size    = '8BYTE00'//char(0)
        header%type    = 'REAL000'//char(0)
        call memory_alloca(mpio_memor,'vec',vacal,vec, PAR_CODE_SIZE)
        call PAR_GATHER(dim2, vec, wherein_code)
        call PAR_BROADCAST(vec, wherein_code)
        vec(1)=0_ip
        header%lines=int(sum(vec), header_ip)
        myoffset=0_ip
        if (INOTMASTER) then
            myoffset=sum(vec(2:PAR_MY_CODE_RANK))
        end if
        header%object=object
        header%resultson=resultson
        if (present(nsubd)) then
            header%nsubd=int(nsubd, header_ip)
        else
            header%nsubd=1_header_ip
        end if
        call par_fill_header(header, divi, tag1, tag2, ittim, cutim, opt)
        call PAR_INFO_CREATE()
        call par_start_timer(wherein_code, header%object(1:5), 'WRITE', 'MPIAL', PAR_CODE_SIZE-1, PAR_CODE_SIZE, header%file_size)
        call PAR_FILE_OPEN_WRITE(fh, filename)
        call PAR_INFO_FREE()
        call PAR_FILE_SET_SIZE(fh, header%file_size)
        call PAR_BARRIER(wherein_code)
        call PAR_FILE_WRITE_HEADER(fh, header)
        call PAR_FILE_SET_VIEW(fh, buf, int(header_size,8))
        offset = int(myoffset,8) * int(header%columns,8)
        call PAR_FILE_SEEK_SET(fh, offset)
        if (.not. mpio_config%output%asynchronous) then
            call PAR_FILE_WRITE_ALL(fh, buf, int(header%columns,ip)*dim2)
            call PAR_FILE_CLOSE(fh)
        else
            call PAR_ASYNCHRONOUS_WRITE_POP_ALL(mpio_config%output%flush, ittim)
            call PAR_ASYNCHRONOUS_WRITE_PUSH(fh, buf, int(header%columns, ip), dim2)
        end if
        call par_end_timer(wherein_code)
        call memory_deallo(mpio_memor,'vec',vacal,vec)
        call PAR_BARRIER(wherein_code)
#endif
    end subroutine

    subroutine PAR_FILE_READ_ALL_INT_V(buf, filename, dim, wherein, header, force_subd)
        integer(ip), pointer,  intent(inout)           :: buf(:)
        character(*),          intent(in)              :: filename
        integer(ip),           intent(in)              :: dim
        integer(ip), pointer                           :: vec(:)
        integer(ip)                                    :: myoffset
        MY_MPI_FILE                                    :: fh
        integer(8)                                     :: offset                              ! Current file offset
        integer(ip)                                    :: lines
        logical                                        :: IREADER
        integer(ip)                                    :: SIZE
        integer(ip)                                    :: RANK
        type(mpio_header),     intent(inout), optional :: header
        type(mpio_header)                              :: h
        character(150),        intent(in)              :: wherein
        logical,               intent(in), optional    :: force_subd
        logical                                        :: force_subd_loc=.false.
        logical                                        :: convert4to8=.false.
        integer(4), pointer                            :: buf4(:)
        character(100), PARAMETER :: vacal = "PAR_FILE_READ_HEADER_INT_V"
        nullify(vec)
        nullify(buf4)
#ifndef MPI_OFF
        if (wherein==wherein_code) then
            IREADER=INOTMASTER
            SIZE=PAR_CODE_SIZE-1
            RANK=PAR_MY_CODE_RANK-1
        else if (wherein==wherein_p) then
            IREADER=IIOSLAVE
            SIZE=PAR_COMM_MPIO_WM_SIZE
            RANK=PAR_COMM_MPIO_RANK_WM
        else
            call runend("Communicator not compatible with MPI-IO service")
        endif
        if( present(header)) then
           h=header
        else 
           call PAR_FILE_READ_HEADER(filename, h, wherein)
        end if
        if (present(force_subd)) then
            force_subd_loc=force_subd
        end if
        if ((.not.force_subd_loc).and.h%nsubd/=1_header_ip .and. h%nsubd/=int(SIZE, header_ip)) then
            call runend ("Trying to read a subdomain file with a different process number")
        end if
        call memory_alloca(mpio_memor,'vec',vacal,vec, SIZE+1)
        call PAR_GATHER(dim, vec, wherein)
        call PAR_BROADCAST(vec, wherein)
        vec(1)=0_ip
        lines=sum(vec)
        if (IREADER .or. IMASTER) then
            if (h%lines/=int(lines, header_ip)) then
                call runend ("File line number ("//intost(int(h%lines, ip))//") and dimension 1 ("//intost(lines)//") are inconsistent")
            end if
            if (h%columns/=1_ip) then
                call runend ("Trying to fill a vector from a multi column file")
            end if
            if (h%type(1:7)/='INTEG00') then
                call runend ("Data contained in the file are not integers")
            end if
            if (h%item_size/=ip) then
                if (h%item_size==4.and.ip==8) then
                    convert4to8=.true.
                else
                    call runend ("File contains integer "//intost(int(h%item_size, ip))//" but Alya is compiled with integer "//intost(int(ip, ip)))
                end if
            end if
            myoffset=0_ip
            if (IREADER) then
                myoffset=sum(vec(2:RANK+1))
            end if
        end if
        call PAR_INFO_CREATE()
        call par_start_timer(wherein, h%object(1:5), 'READ ', 'MPIAL', SIZE, PAR_CODE_SIZE, h%file_size)
        call PAR_FILE_OPEN_READ(fh, filename, wherein)
        call PAR_INFO_FREE()
        call PAR_BARRIER(wherein)
        if (convert4to8) then
            call memory_alloca(mpio_memor,'buf4',vacal,buf4,int(dim,4))
            call PAR_FILE_SET_VIEW(fh, buf4, int(header_size,header_ip))
        else
            call PAR_FILE_SET_VIEW(fh, buf, int(header_size,header_ip))
        end if
        if (IREADER .or. IMASTER) then
            offset = int(myoffset,8)
            call PAR_FILE_SEEK_SET(fh, offset)
            if (convert4to8) then
                call PAR_FILE_READ_ALL(fh, buf4, dim)
                do i=1,dim
                    buf(i)=int(buf4(i),ip)
                end do
                call memory_deallo(mpio_memor,'buf4',vacal,buf4)
            else
                call PAR_FILE_READ_ALL(fh, buf, dim)
            end if
        end if
        call PAR_FILE_CLOSE(fh)
        call par_end_timer(wherein)
        call memory_deallo(mpio_memor,'vec',vacal,vec)
        call PAR_BARRIER(wherein)
#endif
    end subroutine

    subroutine PAR_FILE_READ_ALL_INT_M(buf, filename, dim1, dim2, wherein, header, force_subd)
        integer(ip), pointer,   intent(inout)          :: buf(:,:)
        integer(ip), pointer                           :: temp(:,:)
        character(*),          intent(in)              :: filename
        integer(ip),           intent(in)              :: dim1
        integer(ip),           intent(in)              :: dim2
        integer(ip), pointer                           :: vec(:)
        integer(ip)                                    :: myoffset
        MY_MPI_FILE                                    :: fh
        integer(8)                                     :: offset                              ! Current file offset
        integer(ip)                                    :: lines
        logical                                        :: IREADER
        integer(ip)                                    :: SIZE
        integer(ip)                                    :: RANK
        type(mpio_header),     intent(inout), optional :: header
        type(mpio_header)                              :: h
        character(150),        intent(in)              :: wherein
        logical,               intent(in), optional    :: force_subd
        logical                                        :: force_subd_loc=.false.
        logical                                        :: convert4to8=.false.
        integer(4), pointer                            :: buf4(:,:)
        character(100), PARAMETER :: vacal = "PAR_FILE_READ_HEADER_INT4_M"
        nullify(vec)
        nullify(buf4)
        nullify(temp)
#ifndef MPI_OFF
        if (wherein==wherein_code) then
            IREADER=INOTMASTER
            SIZE=PAR_CODE_SIZE-1
            RANK=PAR_MY_CODE_RANK-1
        else if (wherein==wherein_p) then
            IREADER=IIOSLAVE
            SIZE=PAR_COMM_MPIO_WM_SIZE
            RANK=PAR_COMM_MPIO_RANK_WM
        else
            call runend("Communicator not compatible with MPI-IO service")
        endif
        if (.not.present(header)) then
            call PAR_FILE_READ_HEADER(filename, h, wherein)
        else
            h=header
        end if
        if (present(force_subd)) then
            force_subd_loc=force_subd
        end if
        if ((.not.force_subd_loc).and.h%nsubd/=1_header_ip.and. h%nsubd/=int(SIZE, header_ip)) then
            call runend ("Trying to read a subdomain file with a different process number")
        end if
        call memory_alloca(mpio_memor,'vec',vacal,vec, SIZE+1)
        call PAR_GATHER(dim2, vec, wherein)
        call PAR_BROADCAST(vec, wherein)
        vec(1)=0_ip
        lines=sum(vec)
        if (IREADER .or. IMASTER) then
            if (h%lines/=int(lines, header_ip)) then
                call runend ("File line number ("//intost(int(h%lines, ip))//") and dimension 2 ("//intost(lines)//") are inconsistent")
            end if
            if (h%type(1:7)/='INTEG00') then
                call runend ("Data contained in the file are not integers")
            end if
            if (h%item_size/=ip) then
                if (h%item_size==4.and.ip==8) then
                    convert4to8=.true.
                else
                    call runend ("File contains integer "//intost(int(h%item_size, ip))//" but Alya is compiled with integer "//intost(int(ip, ip)))
                end if
            end if
            myoffset=0_ip
            if (IREADER) then
                myoffset=sum(vec(2:RANK+1))
            end if
        end if
        call PAR_INFO_CREATE()
        call par_start_timer(wherein, h%object(1:5), 'READ ', 'MPIAL', SIZE, PAR_CODE_SIZE, h%file_size)
        call PAR_FILE_OPEN_READ(fh, filename, wherein)
        call PAR_INFO_FREE()
        call PAR_BARRIER(wherein)
        if (convert4to8) then
            call memory_alloca(mpio_memor,'buf4',vacal,buf4, int(h%columns,4), int(dim2, 4))
            call PAR_FILE_SET_VIEW(fh, buf4, int(header_size,header_ip))
        else
            call PAR_FILE_SET_VIEW(fh, buf, int(header_size,header_ip))
        end if
        if (IREADER .or. IMASTER) then
            offset = int(myoffset,8) * int(h%columns,8)
            call PAR_FILE_SEEK_SET(fh, offset)
            if (h%columns==int(dim1, header_ip)) then
                if (convert4to8) then
                    call PAR_FILE_READ_ALL(fh, buf4, int(h%columns,ip)*dim2)
                    do j=1,int(h%columns,4)
                        do i=1,dim2
                            buf(j,i)=int(buf4(j,i),ip)
                        end do
                    end do
                    call memory_deallo(mpio_memor,'buf4',vacal,buf4)
                else
                    call PAR_FILE_READ_ALL(fh, buf, int(h%columns,ip)*dim2)
                end if
            else if (h%columns>int(dim1, header_ip)) then
                if (convert4to8) then
                    call PAR_FILE_READ_ALL(fh, buf4, int(h%columns,ip)*dim2)
                    do j=1,dim1
                        do i=1,dim2
                            buf(j,i)=int(buf4(j,i),ip)
                        end do
                    end do
                    call memory_deallo(mpio_memor,'buf4',vacal,buf4)
                else
                    allocate(temp(h%columns, dim2))
                    call PAR_FILE_READ_ALL(fh, temp, int(h%columns,ip)*dim2)
                    if (IREADER .and. dim1*dim2 > 0) then
                        buf(1:dim1, 1:dim2)=temp(1:dim1, 1:dim2)
                    end if
                    deallocate(temp)
                end if
            else
                call runend ("File column number and dimension 1 are inconsistent")
            end if
        end if
        call PAR_FILE_CLOSE(fh)
        call par_end_timer(wherein)
        call memory_deallo(mpio_memor,'vec',vacal,vec)
        call PAR_BARRIER(wherein)
#endif
    end subroutine

    subroutine PAR_FILE_READ_ALL_REAL_V(buf, filename, dim, wherein, header, force_subd)
        real(rp), pointer,     intent(inout)           :: buf(:)
        character(*),          intent(in)              :: filename
        integer(ip),           intent(in)              :: dim
        integer(ip), pointer                           :: vec(:)
        integer(ip)                                    :: myoffset
        MY_MPI_FILE                                    :: fh
        integer(8)                                     :: offset                              ! Current file offset
        integer(ip)                                    :: lines
        logical                                        :: IREADER
        integer(ip)                                    :: SIZE
        integer(ip)                                    :: RANK
        type(mpio_header),     intent(inout), optional :: header
        type(mpio_header)                              :: h
        character(150),        intent(in)              :: wherein
        logical,               intent(in), optional    :: force_subd
        logical                                        :: force_subd_loc=.false.
        character(100), PARAMETER :: vacal = "PAR_FILE_READ_HEADER_REAL_V"
        nullify(vec)
#ifndef MPI_OFF
        if (wherein==wherein_code) then
            IREADER=INOTMASTER
            SIZE=PAR_CODE_SIZE-1
            RANK=PAR_MY_CODE_RANK-1
        else if (wherein==wherein_p) then
            IREADER=IIOSLAVE
            SIZE=PAR_COMM_MPIO_WM_SIZE
            RANK=PAR_COMM_MPIO_RANK_WM
        else
            call runend("Communicator not compatible with MPI-IO service")
        endif
        if (.not.present(header)) then
            call PAR_FILE_READ_HEADER(filename, h, wherein)
        else
            h=header
        end if
        if (present(force_subd)) then
            force_subd_loc=force_subd
        end if
        if ((.not.force_subd_loc).and.h%nsubd/=1_header_ip.and. h%nsubd/=int(SIZE, header_ip)) then
            call runend ("Trying to read a subdomain file with a different process number")
        end if
        call memory_alloca(mpio_memor,'vec',vacal,vec, SIZE+1)
        call PAR_GATHER(dim, vec, wherein)
        call PAR_BROADCAST(vec, wherein)
        vec(1)=0
        lines=sum(vec)
        if (IREADER .or. IMASTER) then
            if (h%lines/=int(lines, header_ip)) then
                call runend ("File line number ("//intost(int(h%lines, ip))//") and dimension 1 ("//intost(lines)//") are inconsistent")
            end if
            if (h%columns/=1_header_ip) then
                call runend ("Trying to fill a vector from a multi column file")
            end if
            if (h%type(1:7)/='REAL000') then
                call runend ("Data contained in the file are not reals")
            end if
            if (h%item_size/=rp) then
                call runend ("Data contained in the file have the wrong size")
            end if
            myoffset=0_ip
            if (IREADER) then
                myoffset=sum(vec(2:RANK+1))
            end if
        end if
        call PAR_INFO_CREATE()
        call par_start_timer(wherein, h%object(1:5), 'READ ', 'MPIAL', SIZE, PAR_CODE_SIZE, h%file_size)
        call PAR_FILE_OPEN_READ(fh, filename, wherein)
        call PAR_INFO_FREE()
        call PAR_BARRIER(wherein)
        call PAR_FILE_SET_VIEW(fh, buf, int(header_size,8))
        if (IREADER .or. IMASTER) then
            offset = int(myoffset,8)
            call PAR_FILE_SEEK_SET(fh, offset)
            call PAR_FILE_READ_ALL(fh, buf, dim)
        end if
        call PAR_FILE_CLOSE(fh)
        call par_end_timer(wherein)
        call memory_deallo(mpio_memor,'vec',vacal,vec)
        call PAR_BARRIER(wherein)
#endif
    end subroutine

    subroutine PAR_FILE_READ_ALL_REAL_M(buf, filename, dim1, dim2, wherein, header, force_subd)
        real(rp), pointer,   intent(inout)             :: buf(:,:)
        real(rp), pointer                              :: temp(:,:)
        character(*),              intent(in)          :: filename
        integer(ip),               intent(in)          :: dim1
        integer(ip),               intent(in)          :: dim2
        integer(ip), pointer                           :: vec(:)
        integer(ip)                                    :: myoffset
        MY_MPI_FILE                                    :: fh
        integer(8)                                     :: offset                              ! Current file offset
        integer(ip)                                    :: lines
        logical                                        :: IREADER
        integer(ip)                                    :: SIZE
        integer(ip)                                    :: RANK
        type(mpio_header),     intent(inout), optional :: header
        type(mpio_header)                              :: h
        character(150),        intent(in)              :: wherein
        logical,               intent(in), optional    :: force_subd
        logical                                        :: force_subd_loc=.false.
        character(100), PARAMETER :: vacal = "PAR_FILE_READ_HEADER_REAL_M"
        nullify(vec)
        nullify(temp)
#ifndef MPI_OFF
        if (wherein==wherein_code) then
            IREADER=INOTMASTER
            SIZE=PAR_CODE_SIZE-1
            RANK=PAR_MY_CODE_RANK-1
        else if (wherein==wherein_p) then
            IREADER=IIOSLAVE
            SIZE=PAR_COMM_MPIO_WM_SIZE
            RANK=PAR_COMM_MPIO_RANK_WM
        else
            call runend("Communicator not compatible with MPI-IO service")
        endif
        if (.not.present(header)) then
            call PAR_FILE_READ_HEADER(filename, h, wherein)
        else
            h=header
        end if
        if (present(force_subd)) then
            force_subd_loc=force_subd
        end if
        if ((.not.force_subd_loc).and.(h%nsubd/=1_header_ip.and. h%nsubd/=int(SIZE, header_ip))) then
            call runend ("Trying to read a subdomain file with a different process number")
        end if
        call memory_alloca(mpio_memor,'vec',vacal,vec, SIZE+1)
        call PAR_GATHER(dim2, vec, wherein)
        call PAR_BROADCAST(vec, wherein)
        vec(1)=0
        lines=sum(vec)
        if (IREADER .or. IMASTER) then
            if (h%lines/=int(lines, header_ip)) then
                call runend ("File line number ("//intost(int(h%lines, ip))//") and dimension 2 ("//intost(lines)//") are inconsistent")
            end if
            if (h%type(1:7)/='REAL000') then
                call runend ("Data contained in the file are not reals")
            end if
            if (h%item_size/=rp) then
                call runend ("Data contained in the file have the wrong size")
            end if
            myoffset=0
            if (IREADER) then
                myoffset=sum(vec(2:RANK+1))
            end if
        end if
        call PAR_INFO_CREATE()
        call par_start_timer(wherein, h%object(1:5), 'READ ', 'MPIAL', SIZE, PAR_CODE_SIZE, h%file_size)
        call PAR_FILE_OPEN_READ(fh, filename, wherein)
        call PAR_INFO_FREE()
        call PAR_BARRIER(wherein)
        call PAR_FILE_SET_VIEW(fh, buf, int(header_size,8))
        if (IREADER .or. IMASTER) then
            offset = int(myoffset,8) * int(h%columns,8)
            call PAR_FILE_SEEK_SET(fh, offset)
            if (h%columns==int(dim1, header_ip)) then
                call PAR_FILE_READ_ALL(fh, buf, int(h%columns,ip)*dim2)
            else if (h%columns>int(dim1, header_ip)) then
                nullify(temp)
                allocate(temp(h%columns, dim2))
                call PAR_FILE_READ_ALL(fh, temp, int(h%columns,ip)*dim2)
                if (IREADER .and. dim1*dim2 > 0 ) then
                    buf(1:dim1, 1:dim2)=temp(1:dim1, 1:dim2)
                end if
                deallocate(temp)
            else
                call runend ("File column number and dimension 1 are inconsistent")
            end if
        end if
        call PAR_FILE_CLOSE(fh)
        call par_end_timer(wherein)
        call memory_deallo(mpio_memor,'vec',vacal,vec)
        call PAR_BARRIER(wherein)
#endif
    end subroutine

end module
