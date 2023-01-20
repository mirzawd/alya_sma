!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup mpio
!> @{
!> @file    mod_mpio_seq_io.f90
!> @author  Damien Dosimont
!> @date    27/09/2017
!> @brief   Sequential I/O interface
!> @details This module offers a high level interface to I/O sequential operations.
!>          It is responsible for managing MPI-IO file format metadata and file operations.
!> @}
!-----------------------------------------------------------------------

module mod_mpio_seq_io

    use def_master,             only : IIOSLAVE, IMASTER, INOTMASTER, IIOSLAVE, ISEQUEN
    use def_master,             only : namda    
    use def_master,             only : intost, retost
    use def_master,             only : optional_argument
    use def_kintyp,             only : ip, rp, lg
    use mod_memory,             only : memory_alloca, memory_deallo
    use mod_mpio_seq_iowrapper
    use def_mpio
    use mod_mpio_files
    use mod_mpio_seq_log

    implicit none

    integer(ip)                             :: i, j

    private


    interface SEQ_FILE_WRITE_ALL
        module procedure                        SEQ_FILE_WRITE_ALL_INT_V,&
                                                SEQ_FILE_WRITE_ALL_INT_M,&
                                                SEQ_FILE_WRITE_ALL_REAL_V,&
                                                SEQ_FILE_WRITE_ALL_REAL_M
    end interface

    interface SEQ_FILE_READ_ALL
        module procedure                        SEQ_FILE_READ_ALL_INT_V,&
                                                SEQ_FILE_READ_ALL_INT_M,&
                                                SEQ_FILE_READ_ALL_REAL_V,&
                                                SEQ_FILE_READ_ALL_REAL_M
    end interface

    public                                  ::  CHECK_HEADER, FILL_HEADER, FILE_OPTION_ENABLED, SEQ_FILE_READ_ALL, SEQ_FILE_WRITE_ALL, SEQ_FILE_READ_HEADER, SEQ_FILE_WRITE_HEADER


    contains

    subroutine CHECK_HEADER(header)
        type(mpio_header), intent(inout)               :: header
        if (header%magic_number/=header_magic_number) then
            call runend("Magic number is not the expected value (27093). Value detected: "//trim(intost(int(header%magic_number, ip))))
        end if
        if (header%format/=header_format) then
            call runend("File format is not MPIAL! You may want to use the ASCII to MPIO converter available here: https://github.com/dosimont/alya-mpio-tools")
        end if
        if (header%version/=header_version) then
            call runend("The file format version is not compatible!")
        end if
        if (header%dimension(1:7)=='SCALA00' .and. header%columns/=1) then
            call runend("Dimension is defined as scalar but the number of columns is not 1")
        end if
        if (header%dimension(1:7)=='VECTO00' .and. header%columns==1) then
            call runend("Dimension is defined as vector but the number of columns is 1")
        end if
        if (header%dimension(1:7)/='SCALA00'.and.header%dimension(1:7)/='VECTO00') then
            call runend("Unknown dimension!")
        endif
        if (header%par(1:7)/='SEQUE00'.and.header%par(1:7)/='PARAL00') then
            call runend("Is not sequential nor parallel!")
        end if
        if (header%resultson(1:7)=='NELEM00') then
        else if (header%resultson(1:7)=="NPOIN00") then
        else if (header%resultson(1:7)=="NBOUN00") then
        else if (header%resultson(1:7)=="WHATE00") then
        else
            call runend("Unable to detect what the results are based on")
        end if
        if (header%filter(1:7)=='FILTE00') then
        else if (header%filter(1:7)=='NOFIL00') then
        else
            call runend("Unknown filter!")
        endif
        if (header%type(1:7)=='REAL000') then
        else if (header%type(1:7)=='INTEG00') then
        else
            call runend("Unknown type!")
        end if
        if (header%size(1:7)=='4BYTE00') then
        else if (header%size(1:7)=='8BYTE00') then
        else
            call runend("Wrong size!")
        end if
        if (header%sorting(1:7)/='ASCEN00') then
            call runend("Ascending sorting is mandatory!")
        endif
        if (header%id(1:7)/='NOID000') then
            call runend("Id must be absent!")
        endif
    end subroutine

    subroutine FILL_HEADER(header, divi, tag1, tag2, ittim, cutim, opt)
        type(mpio_header), intent(inout)               :: header
        integer(ip),           intent(in)              :: divi
        integer(ip),           intent(in)              :: tag1
        integer(ip),           intent(in)              :: tag2
        integer(ip),            intent(in), optional   :: ittim
        real(rp),               intent(in), optional   :: cutim
        type(mpio_header_options), intent(in), optional:: opt
        if (header%nsubd==1) then
            header%par='SEQUE00'//char(0)
        else
            header%par='PARAL00'//char(0)
        end if
        if (present(ittim)) then
            header%ittim=int(ittim,header_ip)
        else
            header%ittim=0_header_ip
        end if
        if (present(cutim)) then
            header%time=real(cutim, header_rp)
        else
            header%time=0_header_ip
        end if
        if (header%columns>1) then
            header%dimension='VECTO00'//char(0)
        elseif (header%columns==1) then
            header%dimension='SCALA00'//char(0)
        end if
        header%item_size=0_ip
        if (header%size(1:7)=='8BYTE00') then
            header%item_size=8_ip
        elseif (header%size(1:7)=='4BYTE00') then
            header%item_size=4_ip
        end if
        if (present(opt)) then
            header%options=opt
        end if
        header%divi=int(divi,header_ip)
        header%tag1=int(tag1,header_ip)
        header%tag2=int(tag2,header_ip)
        header%file_size=header_size+(header%lines*header%columns*header%item_size)
    end subroutine

    subroutine SEQ_FILE_READ_HEADER(filename, header)
        character(*)                                        ::  filename
        type(mpio_header), intent(inout)                    ::  header
        integer(ip)                                         ::  fh
        if (IMASTER .or. ISEQUEN) then
            call SEQ_FILE_OPEN_READ(fh, filename)
            call SEQ_FILE_READ(fh, header%magic_number)
            call SEQ_FILE_READ(fh, header%format)
            call SEQ_FILE_READ(fh, header%version)
            call SEQ_FILE_READ(fh, header%object)
            call SEQ_FILE_READ(fh, header%dimension)
            call SEQ_FILE_READ(fh, header%resultson)
            call SEQ_FILE_READ(fh, header%type)
            call SEQ_FILE_READ(fh, header%size)
            call SEQ_FILE_READ(fh, header%par)
            call SEQ_FILE_READ(fh, header%filter)
            call SEQ_FILE_READ(fh, header%sorting)
            call SEQ_FILE_READ(fh, header%id)
            call SEQ_FILE_READ(fh, header%align_chars)
            call SEQ_FILE_READ(fh, header%columns)
            call SEQ_FILE_READ(fh, header%lines)
            call SEQ_FILE_READ(fh, header%ittim)
            call SEQ_FILE_READ(fh, header%nsubd)
            call SEQ_FILE_READ(fh, header%divi)
            call SEQ_FILE_READ(fh, header%tag1)
            call SEQ_FILE_READ(fh, header%tag2)
            call SEQ_FILE_READ(fh, header%time)
            call SEQ_FILE_READ(fh, header%align_chars)
            call SEQ_FILE_READ(fh, header%options%opt, option_size)
            header%item_size=0_ip
            if (header%size(1:7)=='8BYTE00') then
                header%item_size=8_ip
            elseif (header%size(1:7)=='4BYTE00') then
                header%item_size=4_ip
            end if
            header%file_size=header_size+(header%lines*header%columns*header%item_size)
            call check_header(header)
            call SEQ_FILE_CLOSE(fh)
        end if
    end subroutine

    function FILE_OPTION_ENABLED(header, opt) result(enabled)
        type(mpio_header), intent(in)   ::  header
        character(7), intent(in)        ::  opt
        logical                         ::  enabled
        enabled=.false.
        do i=1, option_size
            if (opt==header%options%opt(i)(1:7)) then
                enabled=.true.
                return
            end if
        end do
    end function

    subroutine SEQ_FILE_WRITE_HEADER(fh, header,FORCE_WRITE,rank4)
        integer(ip),                 intent(in) :: fh
        type(mpio_header),           intent(in) :: header
        logical(lg),       optional, intent(in) :: FORCE_WRITE
        integer(4),        optional, intent(in) :: rank4
        
        if(  IMASTER .or. &
             ISEQUEN .or. &
             optional_argument(.false.,FORCE_WRITE) .or. &
             optional_argument(1_4,rank4) == 0_4 ) then
           
            call SEQ_FILE_WRITE(fh, header % magic_number)
            call SEQ_FILE_WRITE(fh, header % format)
            call SEQ_FILE_WRITE(fh, header % version)
            call SEQ_FILE_WRITE(fh, header % object)
            call SEQ_FILE_WRITE(fh, header % dimension)
            call SEQ_FILE_WRITE(fh, header % resultson)
            call SEQ_FILE_WRITE(fh, header % type)
            call SEQ_FILE_WRITE(fh, header % size)
            call SEQ_FILE_WRITE(fh, header % par)
            call SEQ_FILE_WRITE(fh, header % filter)
            call SEQ_FILE_WRITE(fh, header % sorting)
            call SEQ_FILE_WRITE(fh, header % id)
            call SEQ_FILE_WRITE(fh, header % align_chars)
            call SEQ_FILE_WRITE(fh, header % columns)
            call SEQ_FILE_WRITE(fh, header % lines)
            call SEQ_FILE_WRITE(fh, header % ittim)
            call SEQ_FILE_WRITE(fh, header % nsubd)
            call SEQ_FILE_WRITE(fh, header % divi)
            call SEQ_FILE_WRITE(fh, header % tag1)
            call SEQ_FILE_WRITE(fh, header % tag2)
            call SEQ_FILE_WRITE(fh, header % time)
            call SEQ_FILE_WRITE(fh, header % align_chars)
            call SEQ_FILE_WRITE(fh, header % options % opt, option_size)
            
         end if
         
    end subroutine

    subroutine SEQ_FILE_WRITE_ALL_INT_V(buf, filename, object, resultson, dim, divi, tag1, tag2, ittim, cutim, opt, nsubd)
        integer(ip), pointer,   intent(in)             :: buf(:)
        character(*),           intent(in)             :: filename
        character(8),           intent(in)             :: object
        character(8),           intent(in)             :: resultson
        integer(ip),            intent(in)             :: dim
        integer(ip),           intent(in)              :: divi
        integer(ip),           intent(in)              :: tag1
        integer(ip),           intent(in)              :: tag2
        integer(ip),            intent(in), optional   :: ittim
        real(rp),               intent(in), optional   :: cutim
        type(mpio_header_options), intent(in), optional:: opt
        integer(ip),            intent(in), optional   :: nsubd
        integer(ip)                                    :: fh
        type(mpio_header)                              :: header
        if (IMASTER .or. ISEQUEN) then
            header%columns=1_header_ip
            if (ip==4) then
                header%size='4BYTE00'//char(0)
            else if (ip==8) then
                header%size='8BYTE00'//char(0)
            end if
            header%type='INTEG00'//char(0)
            header%lines=int(dim,header_ip)
            header%object=object
            header%resultson=resultson
            if (present(nsubd)) then
                header%nsubd=int(nsubd,header_ip)
            else
                header%nsubd=1_header_ip
            end if
            call fill_header(header, divi, tag1, tag2, ittim, cutim, opt)
            call start_timer(header%object(1:5), 'WRITE', 'MPIAL', 1_ip, 1_ip, header%file_size)
            call SEQ_FILE_OPEN_WRITE(fh, filename)
            call SEQ_FILE_WRITE_HEADER(fh, header)
            call SEQ_FILE_SET_VIEW(fh, header_size)
            call SEQ_FILE_WRITE(fh, buf, dim)
            call SEQ_FILE_CLOSE(fh)
            call end_timer()
        end if
    end subroutine

    subroutine SEQ_FILE_WRITE_ALL_INT_M(buf, filename, object, resultson, dim1, dim2, divi, tag1, tag2, ittim, cutim, opt, nsubd)
        integer(ip), pointer,   intent(in)             :: buf(:,:)
        character(*),          intent(in)              :: filename
        character(8),          intent(in)              :: object
        character(8),          intent(in)              :: resultson
        integer(ip),           intent(in)              :: dim1, dim2
        integer(ip),           intent(in)              :: divi
        integer(ip),           intent(in)              :: tag1
        integer(ip),           intent(in)              :: tag2
        integer(ip),           intent(in), optional    :: ittim
        real(rp),               intent(in), optional   :: cutim
        type(mpio_header_options), intent(in), optional:: opt
        integer(ip),            intent(in), optional   :: nsubd
        integer(ip)                                    :: fh
        type(mpio_header)                              :: header
        if (IMASTER .or. ISEQUEN) then
            header%columns=int(dim1,header_ip)
            if (ip==4) then
                header%size='4BYTE00'//char(0)
            else if (ip==8) then
                header%size='8BYTE00'//char(0)
            end if
            header%type='INTEG00'//char(0)
            header%lines=int(dim2,header_ip)
            header%object=object
            header%resultson=resultson
            if (present(nsubd)) then
                header%nsubd=int(nsubd,header_ip)
            else
                header%nsubd=1_header_ip
            end if
            call fill_header(header, divi, tag1, tag2, ittim, cutim, opt)
            call start_timer(header%object(1:5), 'WRITE', 'MPIAL', 1_ip, 1_ip, header%file_size)
            call SEQ_FILE_OPEN_WRITE(fh, filename)
            call SEQ_FILE_WRITE_HEADER(fh, header)
            call SEQ_FILE_SET_VIEW(fh, header_size)
            call SEQ_FILE_WRITE(fh, buf, dim1, dim2)
            call SEQ_FILE_CLOSE(fh)
            call end_timer()
        end if
    end subroutine

    subroutine SEQ_FILE_WRITE_ALL_REAL_V(buf, filename, object, resultson, dim, divi, tag1, tag2, ittim, cutim, opt, nsubd)
        real(rp), pointer,     intent(in)                :: buf(:)
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
        if (IMASTER .or. ISEQUEN) then
            header%columns=1_header_ip
            header%size='8BYTE00'//char(0)
            header%type='REAL000'//char(0)
            header%lines=int(dim,header_ip)
            header%object=object
            header%resultson=resultson
            if (present(nsubd)) then
                header%nsubd=int(nsubd,header_ip)
            else
                header%nsubd=1_header_ip
            end if
            call fill_header(header, divi, tag1, tag2, ittim, cutim, opt)
            call start_timer(header%object(1:5), 'WRITE', 'MPIAL', 1_ip, 1_ip, header%file_size)
            call SEQ_FILE_OPEN_WRITE(fh, filename)
            call SEQ_FILE_WRITE_HEADER(fh, header)
            call SEQ_FILE_SET_VIEW(fh, header_size)
            call SEQ_FILE_WRITE(fh, buf, dim)
            call SEQ_FILE_CLOSE(fh)
            call end_timer()
        end if
    end subroutine

    subroutine SEQ_FILE_WRITE_ALL_REAL_M(buf, filename, object, resultson, dim1, dim2, divi, tag1, tag2, ittim, cutim, opt, nsubd)
        real(rp), pointer,     intent(in)             :: buf(:,:)
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
        integer(ip),            intent(in), optional   :: nsubd
        integer(ip)                                    :: fh
        type(mpio_header)                              :: header
        if (IMASTER .or. ISEQUEN) then
            header%columns=int(dim1,header_ip)
            header%size='8BYTE00'//char(0)
            header%type='REAL000'//char(0)
            header%lines=int(dim2,header_ip)
            header%object=object
            header%resultson=resultson
            if (present(nsubd)) then
                header%nsubd=int(nsubd,header_ip)
            else
                header%nsubd=1_header_ip
            end if
            call fill_header(header, divi, tag1, tag2, ittim, cutim, opt)
            call start_timer(header%object(1:5), 'WRITE', 'MPIAL', 1_ip, 1_ip, header%file_size)
            call SEQ_FILE_OPEN_WRITE(fh, filename)
            call SEQ_FILE_WRITE_HEADER(fh, header)
            call SEQ_FILE_SET_VIEW(fh, header_size)
            call SEQ_FILE_WRITE(fh, buf, dim1, dim2)
            call SEQ_FILE_CLOSE(fh)
            call end_timer()
        end if
    end subroutine

    subroutine SEQ_FILE_READ_ALL_INT_V(buf, filename, dim, header)
        integer(ip), pointer,  intent(inout)           :: buf(:)
        character(*),          intent(in)              :: filename
        integer(ip),           intent(in)              :: dim
        integer(ip)                                    :: fh
        integer(ip)                                    :: lines
        type(mpio_header),     intent(inout), optional :: header
        type(mpio_header)                              :: h
        logical                                        :: convert4to8=.false.
        integer(4), pointer                            :: buf4(:)
        character(100), PARAMETER :: vacal = "SEQ_FILE_READ_ALL_INT_V"
        nullify(buf4)
        if (IMASTER .or. ISEQUEN) then
            if( present(header) ) then
               h=header
            else
               call SEQ_FILE_READ_HEADER(filename, h)               
            end if
            lines=dim
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
                call memory_alloca(mpio_memor,'buf4',vacal,buf4,int(dim,4))
                call SEQ_FILE_READ(fh, buf4, dim)
                do i=1,dim
                    buf(i)=int(buf4(i),ip)
                end do
                call memory_deallo(mpio_memor,'buf4',vacal,buf4)
            else
                call SEQ_FILE_READ(fh, buf, dim)
            end if
            call SEQ_FILE_CLOSE(fh)
            call end_timer()
        end if
    end subroutine

    subroutine SEQ_FILE_READ_ALL_INT_M(buf, filename, dim1, dim2, header)
        integer(ip), pointer,  intent(inout)           :: buf(:,:)
        character(*),          intent(in)              :: filename
        integer(ip),           intent(in)              :: dim1, dim2
        integer(ip)                                    :: fh
        integer(ip)                                    :: lines
        type(mpio_header),     intent(inout), optional :: header
        type(mpio_header)                              :: h
        integer(ip), pointer                           :: temp(:,:)
        logical                                        :: convert4to8=.false.
        integer(4), pointer                            :: buf4(:,:)
        character(100), PARAMETER :: vacal = "SEQ_FILE_READ_ALL_INT_M"
        nullify(temp)
        nullify(buf4)
        if (IMASTER .or. ISEQUEN) then
            if( present(header) ) then
               h=header
            else
               call SEQ_FILE_READ_HEADER(filename, h)              
            end if
            lines=dim2
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
                    call memory_alloca(mpio_memor,'buf4',vacal,buf4, int(h%columns,4), int(dim2, 4))
                    call SEQ_FILE_READ(fh, buf4, int(h%columns, ip), int(dim2, ip))
                    do j=1,h%columns
                        do i=1,dim2
                            buf(j,i)=int(buf4(j,i),ip)
                        end do
                    end do
                    call memory_deallo(mpio_memor,'buf4',vacal,buf4)
                else
                    call SEQ_FILE_READ(fh, buf, int(h%columns, ip), int(dim2, ip))
                end if
            else if (h%columns>int(dim1, header_ip)) then
                if (convert4to8) then
                    call memory_alloca(mpio_memor,'buf4',vacal,buf4, int(h%columns,4), int(dim2, 4))
                    call SEQ_FILE_READ(fh, buf4, int(h%columns, ip), int(dim2, ip))
                    do j=1,dim1
                        do i=1,dim2
                            buf(j,i)=int(buf4(j,i),ip)
                        end do
                    end do
                    call memory_deallo(mpio_memor,'buf4',vacal,buf4)
                else
                    allocate(temp(h%columns, dim2))
                    call SEQ_FILE_READ(fh, temp, int(h%columns, ip), int(dim2, ip))
                    buf(1:dim1, 1:dim2)=temp(1:dim1, 1:dim2)
                    deallocate(temp)
                end if
            else
                call runend ("File column number and dimension 1 are inconsistent")
            end if
            call SEQ_FILE_CLOSE(fh)
            call end_timer()
        end if
    end subroutine

    subroutine SEQ_FILE_READ_ALL_REAL_V(buf, filename, dim, header)
        real(rp), pointer,     intent(inout)           :: buf(:)
        character(*),          intent(in)              :: filename
        integer(ip),           intent(in)              :: dim
        integer(ip)                                    :: fh
        integer(ip)                                    :: lines
        type(mpio_header),     intent(inout), optional :: header
        type(mpio_header)                              :: h
        character(100), PARAMETER :: vacal = "SEQ_FILE_READ_ALL_REAL_V"
        if (IMASTER .or. ISEQUEN) then
            if( present(header) ) then
               h=header
            else
               call SEQ_FILE_READ_HEADER(filename, h)               
            end if
            lines=dim
            if (h%lines/=int(lines, header_ip)) then
                call runend ("File line number and vector size are inconsistent")
            end if
            if (h%columns/=1_header_ip) then
                call runend ("Trying to fill a vector from a multi column file")
            end if
            if (h%type(1:7)/='REAL000') then
                call runend ("Data contained in the file are not reals")
            end if
            if (h%item_size/=rp) then
                call runend ("Data contained in the file have not the good size")
            end if
            call start_timer(h%object(1:5), 'READ', 'MPIAL', 1_ip, 1_ip, h%file_size)
            call SEQ_FILE_OPEN_READ(fh, filename)
            call SEQ_FILE_SET_VIEW(fh, header_size)
            call SEQ_FILE_READ(fh, buf, dim)
            call SEQ_FILE_CLOSE(fh)
            call end_timer()
        end if
    end subroutine

    subroutine SEQ_FILE_READ_ALL_REAL_M(buf, filename, dim1, dim2, header)
        real(rp), pointer,      intent(inout)          :: buf(:,:)
        character(*),          intent(in)              :: filename
        integer(ip),           intent(in)              :: dim1, dim2
        integer(ip)                                    :: fh
        integer(ip)                                    :: lines
        type(mpio_header),     intent(inout), optional :: header
        type(mpio_header)                              :: h
        real(rp), pointer                              :: temp(:,:)
        character(100), PARAMETER :: vacal = "SEQ_FILE_READ_ALL_REAL_M"
        nullify(temp)
        if (IMASTER .or. ISEQUEN) then
            if( present(header)) then
               h=header
            else
               call SEQ_FILE_READ_HEADER(filename, h)               
            end if
            lines=dim2
            if (h%lines/=int(lines, header_ip)) then
                call runend ("File line number and vector size are inconsistent")
            end if
            if (h%type(1:7)/='REAL000') then
                call runend ("Data contained in the file are not reals")
            end if
            if (h%item_size/=rp) then
                call runend ("Data contained in the file have not the good size")
            end if
            call start_timer(h%object(1:5), 'READ', 'MPIAL', 1_ip, 1_ip, h%file_size)
            call SEQ_FILE_OPEN_READ(fh, filename)
            call SEQ_FILE_SET_VIEW(fh, header_size)
            if (h%columns==int(dim1, header_ip)) then
                call SEQ_FILE_READ(fh, buf, int(h%columns, ip), int(dim2, ip))
            else if (h%columns>int(dim1, header_ip)) then
                allocate(temp(h%columns, dim2))
                call SEQ_FILE_READ(fh, temp, int(h%columns, ip), int(dim2, ip))
                buf(1:dim1, 1:dim2)=temp(1:dim1, 1:dim2)
                deallocate(temp)
            else
                call runend ("File column number and dimension 1 are inconsistent")
            end if
            call SEQ_FILE_CLOSE(fh)
            call end_timer()
        end if
    end subroutine

end module
