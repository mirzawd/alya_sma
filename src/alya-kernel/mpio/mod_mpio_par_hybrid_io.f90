!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup mpio
!> @{
!> @file    mod_mpio_par_hybrid_io.f90
!> @author  Damien Dosimont
!> @date    27/09/2017
!> @brief   Parallel to sequential I/O interface
!> @details This module offers a high level interface to I/O sequential operations
!>          It is responsible for managing MPI-IO file format metadata and file operations
!>          as well as inter-process communications to gather/scatter the file data
!> @}
!-----------------------------------------------------------------------

module mod_mpio_par_hybrid_io

    use def_master,             only : IIOSLAVE, IMASTER, INOTMASTER, IIOSLAVE, ISEQUEN
    use def_master,             only : namda
    use def_master,             only : intost, retost
    use def_kintyp,             only : ip, rp, lg
    use mod_memory,             only : memory_alloca, memory_deallo
    use mod_mpio_seq_iowrapper
    use mod_mpio_seq_io
    use mod_mpio_par_io,        only : par_broadcast_header
    use mod_parall,             only : PAR_COMM_WORLD
    use mod_parall,             only : PAR_CODE_SIZE, PAR_MY_CODE_RANK, PAR_COMM_MPIO, PAR_COMM_MPIO_RANK_WM, PAR_COMM_MPIO_WM, PAR_COMM_MPIO_WM_SIZE
    use mod_communications,     only : PAR_COMM_SPLIT, PAR_GATHERV, PAR_GATHER, PAR_SUM
    use mod_communications,     only : PAR_BROADCAST, PAR_BARRIER
    use mod_communications
    use def_mpio
    use mod_mpio_seq_log

    implicit none

    integer(ip)                             :: i, j

    private

    character(150)                          ::  wherein_code="IN MY CODE"


    interface PARSEQ_FILE_WRITE_ALL
        module procedure                        PARSEQ_FILE_WRITE_ALL_INT_V,&
                                                PARSEQ_FILE_WRITE_ALL_INT_M,&
                                                PARSEQ_FILE_WRITE_ALL_REAL_V,&
                                                PARSEQ_FILE_WRITE_ALL_REAL_M
    end interface

    interface PARSEQ_FILE_READ_ALL
        module procedure                        PARSEQ_FILE_READ_ALL_INT_V,&
                                                PARSEQ_FILE_READ_ALL_INT_M,&
                                                PARSEQ_FILE_READ_ALL_REAL_V,&
                                                PARSEQ_FILE_READ_ALL_REAL_M
    end interface

    public                                  ::  PARSEQ_FILE_READ_ALL, PARSEQ_FILE_WRITE_ALL


    contains

    subroutine PARSEQ_FILE_WRITE_ALL_INT_V(buf, filename, object, resultson, dim, divi, tag1, tag2, ittim, cutim, opt, nsubd)
        integer(ip), pointer,  intent(in)             :: buf(:)
        character(*),          intent(in)              :: filename
        character(8),          intent(in)              :: object
        character(8),          intent(in)              :: resultson
        integer(ip),           intent(in)              :: dim
        integer(ip),           intent(in)              :: divi
        integer(ip),           intent(in)              :: tag1
        integer(ip),           intent(in)              :: tag2
        integer(ip),            intent(in), optional   :: ittim
        real(rp),               intent(in), optional   :: cutim
        type(mpio_header_options), intent(in), optional:: opt
        integer(ip),            intent(in), optional   :: nsubd
        integer(ip)                                    :: fh
        type(mpio_header)                              :: header
        integer(ip)                                    :: dimtmp
        integer(ip), pointer                           :: bufrcv(:)
        integer(ip), pointer                           :: bufsend(:)
        integer(4), pointer                            :: vec(:)
        character(100), PARAMETER :: vacal = "PARSEQ_FILE_WRITE_ALL_INT_V"
        dimtmp=dim
        if (IMASTER) then
            dimtmp=0
        end if
        nullify(vec)
        nullify(bufrcv)
        call memory_alloca(mpio_memor,'vec',vacal,vec,PAR_CODE_SIZE)
        call PAR_GATHER(int(dimtmp,4), vec, wherein_code)
        call PAR_BROADCAST(vec, wherein_code)
        call PAR_SUM(dimtmp, wherein_code)
        if (IMASTER) then
            header%columns=1_header_ip
            if (ip==4) then
                header%size='4BYTE00'//char(0)
            else
                header%size='8BYTE00'//char(0)
            end if
            header%type='INTEG00'//char(0)
            header%lines=int(dimtmp, header_ip)
            header%object=object
            header%resultson=resultson
            if (present(nsubd)) then
                header%nsubd=int(nsubd,header_ip)
            else
                header%nsubd=int(1,header_ip)
            end if
            call fill_header(header, divi, tag1, tag2, ittim, cutim, opt)
            call memory_alloca(mpio_memor,'bufrcv',vacal,bufrcv,dimtmp)
        end if
        call par_broadcast_header(header)
        call start_timer(header%object(1:5), 'WRITE', 'MPIAL', 1_ip, 1_ip, header%file_size)
        if (associated(buf)) then
            call PAR_GATHERV(buf, bufrcv, int(dim,4), vec, wherein_code)
        else
            nullify(bufsend)
            call memory_alloca(mpio_memor,'bufsend',vacal,bufsend,1_ip)
            call PAR_GATHERV(bufsend, bufrcv, 0_4, vec, wherein_code)
            call memory_deallo(mpio_memor,'bufsend',vacal,bufsend)
        end if
        if (IMASTER) then
            call SEQ_FILE_OPEN_WRITE(fh, filename)
            call SEQ_FILE_WRITE_HEADER(fh, header)
            call SEQ_FILE_WRITE(fh, bufrcv, dimtmp)
            call SEQ_FILE_CLOSE(fh)
            call memory_deallo(mpio_memor,'bufrcv',vacal,bufrcv)
        end if
        call memory_deallo(mpio_memor,'vec',vacal,vec)
        call end_timer()
    end subroutine

    subroutine PARSEQ_FILE_WRITE_ALL_INT_M(buf, filename, object, resultson, dim1, dim2, divi, tag1, tag2, ittim, cutim, opt, nsubd)
        integer(ip), pointer,   intent(in)             :: buf(:,:)
        character(*),          intent(in)              :: filename
        character(8),          intent(in)              :: object
        character(8),          intent(in)              :: resultson
        integer(ip),           intent(in)              :: dim1, dim2
        integer(ip),           intent(in)              :: divi
        integer(ip),           intent(in)              :: tag1
        integer(ip),           intent(in)              :: tag2
        integer(ip),           intent(in), optional    :: ittim
        real(rp),              intent(in), optional    :: cutim
        type(mpio_header_options), intent(in), optional:: opt
        integer(ip),           intent(in), optional    :: nsubd
        integer(ip)                                    :: fh
        type(mpio_header)                              :: header
        integer(ip)                                    :: dimtmp
        integer(ip), pointer                           :: bufrcv(:,:)
        integer(ip), pointer                           :: bufsend(:,:)
        integer(4), pointer                            :: vec(:)
        character(100), PARAMETER :: vacal = "PARSEQ_FILE_WRITE_ALL_INT_M"
        dimtmp=dim2
        if (IMASTER) then
            dimtmp=0
        end if
        nullify(vec)
        nullify(bufrcv)
        call memory_alloca(mpio_memor,'vec',vacal,vec,PAR_CODE_SIZE)
        call PAR_GATHER(int(dimtmp,4), vec, wherein_code)
        call PAR_BROADCAST(vec, wherein_code)
        vec=vec*int(dim1, 4)
        call PAR_SUM(dimtmp, wherein_code)
        if (IMASTER) then
            header%columns=int(dim1, header_ip)
            if (ip==4) then
                header%size='4BYTE00'//char(0)
            else
                header%size='8BYTE00'//char(0)
            end if
            header%type='INTEG00'//char(0)
            header%lines=int(dimtmp, header_ip)
            header%object=object
            header%resultson=resultson
            if (present(nsubd)) then
                header%nsubd=int(nsubd, header_ip)
            else
                header%nsubd=1_header_ip
            end if
            call fill_header(header, divi, tag1, tag2, ittim, cutim, opt)
            call memory_alloca(mpio_memor,'bufrcv',vacal,bufrcv,dim1,dimtmp)
        end if
        call par_broadcast_header(header)
        call start_timer(header%object(1:5), 'WRITE', 'MPIAL', 1_ip, 1_ip, header%file_size)
        if (associated(buf)) then
            call PAR_GATHERV(buf, bufrcv, int(dim2*dim1,4), vec, wherein_code)
        else
            nullify(bufsend)
            call memory_alloca(mpio_memor,'bufsend',vacal,bufsend,1_ip,1_ip)
            call PAR_GATHERV(bufsend, bufrcv, 0_4, vec, wherein_code)
            call memory_deallo(mpio_memor,'bufsend',vacal,bufsend)
        end if
        if (IMASTER) then
            call SEQ_FILE_OPEN_WRITE(fh, filename)
            call SEQ_FILE_WRITE_HEADER(fh, header)
            call SEQ_FILE_WRITE(fh, bufrcv, int(header%columns, ip),dimtmp)
            call SEQ_FILE_CLOSE(fh)
            call memory_deallo(mpio_memor,'bufrcv',vacal,bufrcv)
        end if
        call memory_deallo(mpio_memor,'vec',vacal,vec)
        call end_timer()
    end subroutine

    subroutine PARSEQ_FILE_WRITE_ALL_REAL_V(buf, filename, object, resultson, dim, divi, tag1, tag2, ittim, cutim, opt, nsubd)

        real(rp), pointer,   intent(in)                 :: buf(:)
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
        integer(ip)                                    :: fh
        type(mpio_header)                              :: header
        integer(ip)                                    :: dimtmp
        real(rp), pointer                               :: bufsend(:)
        real(rp), pointer                               :: bufrcv(:)
        integer(4), pointer                            :: vec(:)
        character(100), PARAMETER :: vacal = "PARSEQ_FILE_WRITE_ALL_REAL_V"

        dimtmp=dim
        if (IMASTER) then
            dimtmp=0
        end if
        nullify(vec)
        nullify(bufrcv)
        call memory_alloca(mpio_memor,'vec',vacal,vec,PAR_CODE_SIZE)
        call PAR_GATHER(int(dimtmp,4), vec, wherein_code)
        call PAR_BROADCAST(vec, wherein_code)
        call PAR_SUM(dimtmp, wherein_code)
        if (IMASTER) then
            header%columns=1_header_ip
            header%size='8BYTE00'//char(0)
            header%type='REAL000'//char(0)
            header%lines=int(dimtmp, header_ip)
            header%object=object
            header%resultson=resultson
            if (present(nsubd)) then
                header%nsubd=int(nsubd, header_ip)
            else
                header%nsubd=1_header_ip
            end if
            call fill_header(header, divi, tag1, tag2, ittim, cutim, opt)
            call memory_alloca(mpio_memor,'bufrcv',vacal,bufrcv,dimtmp)
        end if
        call par_broadcast_header(header)
        call start_timer(header%object(1:5), 'WRITE', 'MPIAL', 1_ip, 1_ip, header%file_size)
        if (associated(buf)) then
            call PAR_GATHERV(buf, bufrcv, int(dim,4), vec, wherein_code)
        else
            nullify(bufsend)
            call memory_alloca(mpio_memor,'bufsend',vacal,bufsend,1_ip)
            call PAR_GATHERV(bufsend, bufrcv, 0_4, vec, wherein_code)
            call memory_deallo(mpio_memor,'bufsend',vacal,bufsend)
        end if
        if (IMASTER) then
            call SEQ_FILE_OPEN_WRITE(fh, filename)
            call SEQ_FILE_WRITE_HEADER(fh, header)
            call SEQ_FILE_WRITE(fh, bufrcv, dimtmp)
            call SEQ_FILE_CLOSE(fh)
            call memory_deallo(mpio_memor,'bufrcv',vacal,bufrcv)
        end if
        call memory_deallo(mpio_memor,'vec',vacal,vec)
        call end_timer()
    end subroutine

    subroutine PARSEQ_FILE_WRITE_ALL_REAL_M(buf, filename, object, resultson, dim1, dim2, divi, tag1, tag2, ittim, cutim, opt, nsubd)
        real(rp), pointer,  intent(in)                 :: buf(:,:)
        character(*),          intent(in)              :: filename
        character(8),          intent(in)              :: object
        character(8),          intent(in)              :: resultson
        integer(ip),           intent(in)              :: dim1, dim2
        integer(ip),           intent(in)              :: divi
        integer(ip),           intent(in)              :: tag1
        integer(ip),           intent(in)              :: tag2
        integer(ip),           intent(in), optional    :: ittim
        real(rp),              intent(in), optional    :: cutim
        type(mpio_header_options), intent(in), optional:: opt
        integer(ip),           intent(in), optional    :: nsubd
        integer(ip)                                    :: fh
        type(mpio_header)                              :: header
        integer(ip)                                    :: dimtmp
        integer(4)                                     :: dimtmp4
        real(rp), pointer                              :: bufsend(:,:)
        real(rp), pointer                              :: bufrcv(:,:)
        integer(4), pointer                            :: vec(:)
        character(100), PARAMETER :: vacal = "PARSEQ_FILE_WRITE_ALL_REAL_M"
        dimtmp=dim2
        if (IMASTER) then
           dimtmp=0
        end if
        dimtmp4 = int(dimtmp,4)
        nullify(vec)
        nullify(bufrcv)
        nullify(bufsend)

        call memory_alloca(mpio_memor,'vec',vacal,vec,PAR_CODE_SIZE)
        call PAR_GATHER(dimtmp4, vec, wherein_code)
        call PAR_BARRIER()
        call PAR_BROADCAST(vec, wherein_code)
        call PAR_BARRIER()
        vec=vec*int(dim1,4)
        call PAR_SUM(dimtmp, wherein_code)
        if (IMASTER) then
            header%columns   = int(dim1, header_ip)
            header%size      = '8BYTE00'//char(0)
            header%type      = 'REAL000'//char(0)
            header%lines     = int(dimtmp, header_ip)
            header%object    = object
            header%resultson = resultson
            if (present(nsubd)) then
                header%nsubd=int(nsubd, header_ip)
            else
                header%nsubd=1_header_ip
            end if
            call fill_header(header, divi, tag1, tag2, ittim, cutim, opt)
            call memory_alloca(mpio_memor,'bufrcv',vacal,bufrcv,dim1,dimtmp)
        end if
        call par_broadcast_header(header)
        call start_timer(header%object(1:5), 'WRITE', 'MPIAL', 1_ip, 1_ip, header%file_size)
        if (associated(buf)) then
            call PAR_GATHERV(buf, bufrcv, int(dim2*dim1, 4), vec, wherein_code)
        else
            call memory_alloca(mpio_memor,'bufsend',vacal,bufsend,1_ip,1_ip)
            call PAR_GATHERV(bufsend, bufrcv, 0_4, vec, wherein_code)
            call memory_deallo(mpio_memor,'bufsend',vacal,bufsend)
         end if
        if (IMASTER) then
            call SEQ_FILE_OPEN_WRITE(fh, filename)
            call SEQ_FILE_WRITE_HEADER(fh, header)
            call SEQ_FILE_WRITE(fh, bufrcv, int(header%columns, ip),dimtmp)
            call SEQ_FILE_CLOSE(fh)
            call memory_deallo(mpio_memor,'bufrcv',vacal,bufrcv)
        end if
        call memory_deallo(mpio_memor,'vec',vacal,vec)
        call end_timer()
    end subroutine

    subroutine PARSEQ_FILE_READ_ALL_INT_V(buf, filename, dim, header)
        integer(ip), pointer, intent(inout)            :: buf(:)
        character(*),          intent(in)              :: filename
        integer(ip),               intent(in)          :: dim
        integer(ip)                                    :: fh
        integer(ip)                                    :: lines
        type(mpio_header),     intent(inout), optional :: header
        type(mpio_header)                              :: h
        character(100), PARAMETER :: vacal = "PARSEQ_FILE_READ_ALL_INT_V"
        integer(4), pointer                            :: vec(:)
        integer(ip), pointer                           :: bufsend(:)
        logical                                        :: convert4to8=.false.
        integer(4), pointer                            :: buf4(:)
        nullify(vec)
        nullify(buf4)
        lines=dim
        if (IMASTER) then
            lines=0
        end if
        call memory_alloca(mpio_memor,'vec',vacal,vec,PAR_CODE_SIZE)
        call PAR_GATHER(int(lines,4), vec, wherein_code)
        call PAR_BROADCAST(vec, wherein_code)
        call PAR_SUM(lines, wherein_code)
        if (IMASTER .or. ISEQUEN) then
            nullify(bufsend)
            call memory_alloca(mpio_memor,'buf',vacal,bufsend,lines)
            if (.not.present(header)) then
                call SEQ_FILE_READ_HEADER(filename, h)
            else
                h=header
            end if
            if (h%lines/=int(lines, header_ip)) then
                call runend ("File line number and vector size are inconsistent")
            end if
            if (h%columns/=1_header_ip) then
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
            call start_timer(h%object(1:5), 'READ', 'MPIAL', 1_ip, 1_ip, h%file_size)
            call SEQ_FILE_OPEN_READ(fh, filename)
            call SEQ_FILE_SET_VIEW(fh, header_size)
            if (convert4to8) then
                call memory_alloca(mpio_memor,'buf4',vacal,buf4,int(lines,4))
                call SEQ_FILE_READ(fh, buf4, lines)
                do i=1,lines
                    bufsend(i)=int(buf4(i),ip)
                end do
                call memory_deallo(mpio_memor,'buf4',vacal,buf4)
            else
                call SEQ_FILE_READ(fh, bufsend, lines)
            end if
            call SEQ_FILE_CLOSE(fh)
        end if
        call par_broadcast_header(h)
        call PAR_SCATTERV(bufsend, buf, vec, vec(PAR_MY_CODE_RANK+1), wherein=wherein_code)
        if (IMASTER) then
            call memory_deallo(mpio_memor,'bufsend',vacal,bufsend)
        end if
        call memory_deallo(mpio_memor,'vec',vacal,vec)
        call end_timer()
    end subroutine

    subroutine PARSEQ_FILE_READ_ALL_INT_M(buf, filename, dim1, dim2, header)
        integer(ip), pointer,  intent(inout)           :: buf(:,:)
        character(*),          intent(in)              :: filename
        integer(ip),           intent(in)              :: dim1, dim2
        integer(ip)                                    :: fh
        integer(ip)                                    :: lines
        type(mpio_header),     intent(inout), optional :: header
        type(mpio_header)                              :: h
        character(100), PARAMETER :: vacal = "SEQ_FILE_READ_ALL_INT_M"
        integer(4), pointer                            :: vec(:)
        integer(ip), pointer                           :: bufsend(:,:)
        integer(ip), pointer                           :: temp(:,:)
        logical                                        :: convert4to8=.false.
        integer(4), pointer                            :: buf4(:,:)
        nullify(temp)
        nullify(buf4)
        lines=dim2
        if (IMASTER) then
            lines=0
        end if
        nullify(vec)
        call memory_alloca(mpio_memor,'vec',vacal,vec,PAR_CODE_SIZE)
        call PAR_GATHER(int(lines,4), vec, wherein_code)
        call PAR_BROADCAST(vec, wherein_code)
        vec=vec*int(dim1, 4)
        call PAR_SUM(lines, wherein_code)
        if (IMASTER .or. ISEQUEN) then
            nullify(bufsend)
            call memory_alloca(mpio_memor,'buf',vacal,bufsend,dim1,lines)
            if (.not.present(header)) then
                call SEQ_FILE_READ_HEADER(filename, h)
            else
                h=header
            end if
            if (h%lines/=int(lines, header_ip)) then
                call runend ("File line number and vector size are inconsistent")
            end if
            if (h%type(1:7)/='INTEG00') then
                call runend ("Data contained in the file are not reals")
            end if
            if (h%item_size/=ip) then
                if (h%item_size==4.and.ip==8) then
                    convert4to8=.true.
                else
                    call runend ("File contains integer "//intost(int(h%item_size, ip))//" but Alya is compiled with integer "//intost(int(ip, ip)))
                end if
            end if
            call start_timer(h%object(1:5), 'READ', 'MPIAL', 1_ip, 1_ip, h%file_size)
            call SEQ_FILE_OPEN_READ(fh, filename)
            call SEQ_FILE_SET_VIEW(fh, header_size)
            if (h%columns==int(dim1, header_ip)) then
                if (convert4to8) then
                    call memory_alloca(mpio_memor,'buf4',vacal,buf4, int(h%columns,4), int(lines, 4))
                    call SEQ_FILE_READ(fh, buf4, int(h%columns, ip), int(lines, ip))
                    do j=1,h%columns
                        do i=1,lines
                            bufsend(j,i)=int(buf4(j,i),ip)
                        end do
                    end do
                    call memory_deallo(mpio_memor,'buf4',vacal,buf4)
                else
                    call SEQ_FILE_READ(fh, bufsend, int(h%columns, ip), int(lines, ip))
                end if
            else if (h%columns>int(dim1, header_ip)) then
                if (convert4to8) then
                    call memory_alloca(mpio_memor,'buf4',vacal,buf4, int(h%columns,4), int(lines, 4))
                    call SEQ_FILE_READ(fh, buf4, int(h%columns, ip), int(lines, ip))
                    do j=1,h%columns
                        do i=1,lines
                            buf(j,i)=int(buf4(j,i),ip)
                        end do
                    end do
                    call memory_deallo(mpio_memor,'buf4',vacal,buf4)
                else
                    allocate(temp(h%columns, lines))
                    call SEQ_FILE_READ(fh, temp, int(h%columns, ip), int(lines, ip))
                    buf(1:dim1, 1:lines)=temp(1:dim1, 1:lines)
                    deallocate(temp)
                end if
            else
                call runend ("File column number and dimension 1 are inconsistent")
            end if
            call SEQ_FILE_CLOSE(fh)
        end if
        call par_broadcast_header(h)
        call PAR_SCATTERV(bufsend, buf, vec, vec(PAR_MY_CODE_RANK+1), wherein_code)
        if (IMASTER) then
            call memory_deallo(mpio_memor,'bufsend',vacal,bufsend)
        end if
        call memory_deallo(mpio_memor,'vec',vacal,vec)
        call end_timer()
    end subroutine

    subroutine PARSEQ_FILE_READ_ALL_REAL_V(buf, filename, dim, header)
        real(rp), pointer, intent(inout)                :: buf(:)
        character(*),          intent(in)              :: filename
        integer(ip),               intent(in)              :: dim
        integer(ip)                                        :: fh
        integer(ip)                                        :: lines
        type(mpio_header),     intent(inout), optional :: header
        type(mpio_header)                              :: h
        character(100), PARAMETER :: vacal = "SEQ_FILE_READ_ALL_REAL_V"
        integer(4), pointer                            :: vec(:)
        real(rp), pointer                              :: bufsend(:)
        lines=dim
        if (IMASTER) then
            lines=0
        end if
        nullify(vec)
        call memory_alloca(mpio_memor,'vec',vacal,vec,PAR_CODE_SIZE)
        call PAR_GATHER(int(lines,4), vec, wherein_code)
        call PAR_BROADCAST(vec, wherein_code)
        call PAR_SUM(lines, wherein_code)
        if (IMASTER .or. ISEQUEN) then
            nullify(bufsend)
            call memory_alloca(mpio_memor,'buf',vacal,bufsend,lines)
            if (.not.present(header)) then
                call SEQ_FILE_READ_HEADER(filename, h)
            else
                h=header
            end if
            if (h%lines/=int(lines, header_ip)) then
                call runend ("File line number and vector size are inconsistent")
            end if
            if (h%columns/=1_header_ip) then
                call runend ("Trying to fill a vector from a multi column file")
            end if
            if (h%type(1:7)/='REAL000') then
                call runend ("Data contained in the file are not reals")
            end if
            if (h%item_size/=8) then
                call runend ("Data contained in the file have the wrong size")
            end if
            call start_timer(h%object(1:5), 'READ', 'MPIAL', 1_ip, 1_ip, h%file_size)
            call SEQ_FILE_OPEN_READ(fh, filename)
            call SEQ_FILE_SET_VIEW(fh, header_size)
            call SEQ_FILE_READ(fh, bufsend, lines)
            call SEQ_FILE_CLOSE(fh)
        end if
        call par_broadcast_header(h)
        call PAR_SCATTERV(bufsend, buf, vec, vec(PAR_MY_CODE_RANK+1), wherein_code)
        if (IMASTER) then
            call memory_deallo(mpio_memor,'bufsend',vacal,bufsend)
        end if
        call memory_deallo(mpio_memor,'vec',vacal,vec)
        call end_timer()
    end subroutine

    subroutine PARSEQ_FILE_READ_ALL_REAL_M(buf, filename, dim1, dim2, header)
        real(rp), pointer,      intent(inout)          :: buf(:,:)
        character(*),          intent(in)              :: filename
        integer(ip),               intent(in)          :: dim1, dim2
        integer(ip)                                    :: fh
        integer(ip)                                    :: lines
        type(mpio_header),     intent(inout), optional :: header
        type(mpio_header)                              :: h
        character(100), PARAMETER :: vacal = "SEQ_FILE_READ_ALL_REAL_M"
        integer(4), pointer                            :: vec(:)
        real(rp), pointer                              :: bufsend(:,:)
        real(rp), pointer                              :: temp(:,:)
        nullify(temp)
        lines=dim2
        if (IMASTER) then
            lines=0
        end if
        nullify(vec)
        call memory_alloca(mpio_memor,'vec',vacal,vec,PAR_CODE_SIZE)
        call PAR_GATHER(int(lines, 4), vec, wherein_code)
        call PAR_BROADCAST(vec, wherein_code)
        vec=vec*int(dim1, 4)
        call PAR_SUM(lines, wherein_code)
        if (IMASTER .or. ISEQUEN) then
            nullify(bufsend)
            call memory_alloca(mpio_memor,'buf',vacal,bufsend,dim1,lines)
            if (.not.present(header)) then
                call SEQ_FILE_READ_HEADER(filename, h)
            else
                h=header
            end if
            if (h%lines/=int(lines, header_ip)) then
                call runend ("File line number and vector size are inconsistent")
            end if
            if (h%type(1:7)/='REAL000') then
                call runend ("Data contained in the file are not reals")
            end if
            if (h%item_size/=8) then
                call runend ("Data contained in the file are not 8 bytes")
            end if
            call start_timer(h%object(1:5), 'READ', 'MPIAL', 1_ip, 1_ip, h%file_size)
            call SEQ_FILE_OPEN_READ(fh, filename)
            call SEQ_FILE_SET_VIEW(fh, header_size)
            if (h%columns==int(dim1, header_ip)) then
                call SEQ_FILE_READ(fh, bufsend, int(h%columns, ip), int(lines, ip))
            else if (h%columns>int(dim1, header_ip)) then
                allocate(temp(h%columns, lines))
                call SEQ_FILE_READ(fh, temp, int(h%columns, ip), int(lines, ip))
                buf(1:dim1, 1:lines)=temp(1:dim1, 1:lines)
                deallocate(temp)
            else
                call runend ("File column number and dimension 1 are inconsistent")
            end if
            call SEQ_FILE_CLOSE(fh)
        end if
        call par_broadcast_header(h)
        call PAR_SCATTERV(bufsend, buf, vec, vec(PAR_MY_CODE_RANK+1), wherein=wherein_code)
        if (IMASTER) then
            call memory_deallo(mpio_memor,'bufsend',vacal,bufsend)
        end if
        call memory_deallo(mpio_memor,'vec',vacal,vec)
        call end_timer()
    end subroutine

end module
