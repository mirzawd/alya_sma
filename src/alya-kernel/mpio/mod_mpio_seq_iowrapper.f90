!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup mpio
!> @{
!> @file    mod_mpio_seq_iowrapper.f90
!> @author  Damien Dosimont
!> @date    27/09/2017
!> @brief   IO Wrapper
!> @details This module is a wrapper of classic IO functions.
!>          It makes type management transparent to the user.
!> @}
!-----------------------------------------------------------------------

module mod_mpio_seq_iowrapper

    use def_kintyp,                     only : ip,rp,lg
    use def_master,                     only : IPARSLAVE, IMASTER, INOTMASTER, ISEQUEN

    implicit none

    private

    integer(4), parameter                       :: SEEK_SET = 0, SEEK_CUR = 1, SEEK_END = 2
    integer(ip)                                 :: u=195878
!    integer(4)                                  :: ierr

    integer(ip)                                 ::  i

    integer(ip), dimension(1:1000)              :: hm(1000)
    integer(ip), dimension(1:1000)              :: rw(1000)

    interface SEQ_FILE_READ
        module procedure                        SEQ_FILE_READ_INT4,&
                                                SEQ_FILE_READ_INT4_V,&
                                                SEQ_FILE_READ_INT4_M,&
                                                SEQ_FILE_READ_INT8,&
                                                SEQ_FILE_READ_INT8_V,&
                                                SEQ_FILE_READ_INT8_M,&
                                                SEQ_FILE_READ_REAL8,&
                                                SEQ_FILE_READ_REAL8_V,&
                                                SEQ_FILE_READ_REAL8_M,&
                                                SEQ_FILE_READ_CHAR8,&
                                                SEQ_FILE_READ_CHAR8_P
    end interface

    interface SEQ_FILE_WRITE
        module procedure                        SEQ_FILE_WRITE_INT4,&
                                                SEQ_FILE_WRITE_INT4_V,&
                                                SEQ_FILE_WRITE_INT4_M,&
                                                SEQ_FILE_WRITE_INT8,&
                                                SEQ_FILE_WRITE_INT8_V,&
                                                SEQ_FILE_WRITE_INT8_M,&
                                                SEQ_FILE_WRITE_REAL8,&
                                                SEQ_FILE_WRITE_REAL8_V,&
                                                SEQ_FILE_WRITE_REAL8_M,&
                                                SEQ_FILE_WRITE_CHAR8,&
                                                SEQ_FILE_WRITE_CHAR8_P

    end interface

    public                                  ::  SEQ_FILE_OPEN_READ,&
                                                SEQ_FILE_OPEN_WRITE,&
                                                SEQ_FILE_OPEN_WRITE_APPEND,&
                                                SEQ_INFO_CREATE,&
                                                SEQ_INFO_FREE,&
                                                SEQ_INFO_SET,&
                                                SEQ_FILE_CLOSE,&
                                                SEQ_FILE_SET_VIEW,&
                                                SEQ_FILE_SEEK_SET,&
                                                SEQ_FILE_SET_SIZE,&
                                                SEQ_FILE_GET_SIZE,&
                                                SEQ_FILE_READ,&
                                                SEQ_FILE_WRITE



    contains

    subroutine SEQ_INFO_CREATE()
    end subroutine

    subroutine SEQ_INFO_SET(field, value)
        character(*),          intent(in)              :: field, value
    end subroutine

    subroutine SEQ_INFO_FREE()
    end subroutine

    subroutine SEQ_FILE_OPEN_READ(fh, filename)
        integer(ip),           intent(inout)           :: fh
        character(*),          intent(in)              :: filename
        integer(ip)                                    :: stat
        u=u+1
        fh=u
        open(unit=fh, file=filename, iostat=stat, action="read", form='unformatted', access='stream')
        do i=1,100
          if (rw(i)==0) then
            rw(i)=1
            hm(i)=fh
            exit
          end if
        end do
    end subroutine

    subroutine SEQ_FILE_OPEN_WRITE(fh, filename)
        integer(ip),           intent(inout)           :: fh
        character(*),          intent(in)              :: filename
        integer(ip)                                    :: stat
        u=u+1
        fh=u
        open(unit=fh, iostat=stat, file=filename, status='old')
        if (stat == 0) close (fh, status='delete')
        open(fh, file=filename, status="new", action="write", form='unformatted', access='stream')
        do i=1,100
          if (rw(i)==0) then
            rw(i)=2
            hm(i)=fh
            exit
          end if
        end do
    end subroutine

    subroutine SEQ_FILE_OPEN_WRITE_APPEND(fh, filename)
        integer(ip),           intent(inout)           :: fh
        character(*),          intent(in)              :: filename
        u=u+1
        fh=u
        open(fh, file=filename, status="old", action="write", position="append", form='unformatted', access='stream')
        do i=1,100
          if (rw(i)==0) then
            rw(i)=3
            hm(i)=fh
            exit
          end if
        end do
    end subroutine

    subroutine SEQ_FILE_CLOSE(fh)
        integer(ip),           intent(inout)           :: fh
        close(fh)
        do i=1,100
          if (hm(i)==fh) then
            rw(i)=0
          end if
        end do
    end subroutine

    subroutine SEQ_FILE_SET_VIEW(fh, offset)
        integer(ip),           intent(inout)           :: fh
        integer(ip),           intent(in)              :: offset
        call SEQ_FILE_SEEK_SET(fh, offset)
    end subroutine

    subroutine SEQ_FILE_READ_INT4(fh, buf)
        integer(ip),           intent(inout)           :: fh
        integer(4),            intent(inout)           :: buf
        read(fh) buf
    end subroutine

    subroutine SEQ_FILE_READ_INT4_V(fh, buf, bsize)
        integer(ip),           intent(inout)           :: fh
        integer(4), pointer,   intent(inout)           :: buf(:)
        integer(ip),           intent(in)              :: bsize
        read(fh) buf(1:bsize)
    end subroutine

    subroutine SEQ_FILE_READ_INT4_M(fh, buf, bsize1, bsize2)
        integer(ip),           intent(in)              :: fh
        integer(4), pointer,   intent(inout)           :: buf(:,:)
        integer(ip),           intent(in)              :: bsize1, bsize2
        read(fh) buf(1:bsize1, 1:bsize2)
    end subroutine

    subroutine SEQ_FILE_READ_INT8(fh, buf)
        integer(ip),           intent(inout)           :: fh
        integer(8),            intent(inout)           :: buf
        read(fh) buf
    end subroutine

    subroutine SEQ_FILE_READ_INT8_V(fh, buf, bsize)
        integer(ip),           intent(inout)           :: fh
        integer(8), pointer,   intent(inout)           :: buf(:)
        integer(ip),           intent(in)              :: bsize
        read(fh) buf(1:bsize)
    end subroutine

    subroutine SEQ_FILE_READ_INT8_M(fh, buf, bsize1, bsize2)
        integer(ip),           intent(in)              :: fh
        integer(8), pointer,   intent(inout)           :: buf(:,:)
        integer(ip),           intent(in)              :: bsize1, bsize2
        read(fh) buf(1:bsize1, 1:bsize2)
    end subroutine

    subroutine SEQ_FILE_READ_REAL8(fh, buf)
        integer(ip),           intent(in)              :: fh
        real(8),               intent(inout)           :: buf
        read(fh) buf
    end subroutine

    subroutine SEQ_FILE_READ_REAL8_V(fh, buf, bsize)
        integer(ip),           intent(in)              :: fh
        real(8), pointer,      intent(inout)           :: buf(:)
        integer(ip),           intent(in)              :: bsize
        read(fh) buf(1:bsize)
    end subroutine

    subroutine SEQ_FILE_READ_REAL8_M(fh, buf, bsize1, bsize2)
        integer(ip),           intent(in)              :: fh
        real(8), pointer,      intent(inout)           :: buf(:,:)
        integer(ip),           intent(in)              :: bsize1, bsize2
        read(fh) buf(1:bsize1, 1:bsize2)
    end subroutine

    subroutine SEQ_FILE_READ_CHAR8(fh, buf)
        integer(ip),           intent(in)              :: fh
        character(8),          intent(inout)           :: buf
        read(fh) buf
    end subroutine

    subroutine SEQ_FILE_READ_CHAR8_P(fh, buf, bsize)
        integer(ip),           intent(in)              :: fh
        character(8), dimension(*), intent(inout)      :: buf
        integer(ip),           intent(in)              :: bsize
        read(fh) buf(1:bsize)
    end subroutine

    subroutine SEQ_FILE_WRITE_INT4(fh, buf)
        integer(ip),           intent(in)              :: fh
        integer(4),            intent(in)              :: buf
        write(fh) buf
    end subroutine

    subroutine SEQ_FILE_WRITE_INT4_V(fh, buf, bsize)
        integer(ip),           intent(in)              :: fh
        integer(4), pointer,   intent(in)              :: buf(:)
        integer(ip),           intent(in)              :: bsize
        if (associated(buf)) write(fh) buf(1:bsize)
    end subroutine

    subroutine SEQ_FILE_WRITE_INT4_M(fh, buf, bsize1, bsize2)
        integer(ip),           intent(in)              :: fh
        integer(4), pointer,   intent(in)              :: buf(:,:)
        integer(ip),           intent(in)              :: bsize1, bsize2
        if (associated(buf)) write(fh) buf(1:bsize1, 1:bsize2)
    end subroutine

    subroutine SEQ_FILE_WRITE_INT8(fh, buf)
        integer(ip),           intent(in)              :: fh
        integer(8),            intent(in)              :: buf
        write(fh) buf
    end subroutine

    subroutine SEQ_FILE_WRITE_INT8_V(fh, buf, bsize)
        integer(ip),           intent(in)              :: fh
        integer(8), pointer,   intent(in)              :: buf(:)
        integer(ip),           intent(in)              :: bsize
        if (associated(buf)) write(fh) buf(1:bsize)
    end subroutine

    subroutine SEQ_FILE_WRITE_INT8_M(fh, buf, bsize1, bsize2)
        integer(ip),           intent(in)              :: fh
        integer(8), pointer,   intent(in)              :: buf(:,:)
        integer(ip),           intent(in)              :: bsize1, bsize2
        if (associated(buf)) write(fh) buf(1:bsize1, 1:bsize2)
    end subroutine

    subroutine SEQ_FILE_WRITE_REAL8(fh, buf)
        integer(ip),           intent(in)              :: fh
        real(8),               intent(in)              :: buf
        write(fh) buf
    end subroutine

    subroutine SEQ_FILE_WRITE_REAL8_V(fh, buf, bsize)
        integer(ip),           intent(in)              :: fh
        real(8), pointer,      intent(in)              :: buf(:)
        integer(ip),           intent(in)              :: bsize
        if (associated(buf)) write(fh) buf(1:bsize)
    end subroutine

    subroutine SEQ_FILE_WRITE_REAL8_M(fh, buf, bsize1, bsize2)
        integer(ip),           intent(in)              :: fh
        real(8), pointer,      intent(in)              :: buf(:,:)
        integer(ip),           intent(in)              :: bsize1, bsize2
        if (associated(buf)) write(fh) buf(1:bsize1, 1:bsize2)
    end subroutine

    subroutine SEQ_FILE_WRITE_CHAR8(fh, buf)
        integer(ip),           intent(in)              :: fh
        character(8),          intent(in)              :: buf
        write(fh) buf
    end subroutine

    subroutine SEQ_FILE_WRITE_CHAR8_P(fh, buf, bsize)
        integer(ip),           intent(in)              :: fh
        character(8), dimension(*), intent(in)         :: buf
        integer(ip),           intent(in)              :: bsize
        write(fh) buf(1:bsize)
    end subroutine

    subroutine SEQ_FILE_SEEK_SET(fh, offset)
        integer(ip),           intent(in)              :: fh
        integer(ip),           intent(in)              :: offset
        integer(4)                                     :: dummy(offset/4)
        !TODO find a solution to make this call compatible with the standard
        !CALL FSEEK(fh, offset, SEEK_SET, ierr)
        if (hm(i)==fh) then
          if (rw(i)==1) then
            rewind(fh)
            read(fh) dummy(1:offset/4)
          end if
        end if
    end subroutine

    subroutine SEQ_FILE_SET_SIZE(filename, bsize)
        character(*),          intent(in)              :: filename
        integer(ip),           intent(in)              :: bsize
        call runend("SEQ_FILE_SET_SIZE not implemented")
    end subroutine

    subroutine SEQ_FILE_GET_SIZE(filename, bsize)
        character(*),          intent(in)              :: filename
        integer(ip),           intent(inout)           :: bsize
        INQUIRE(FILE=filename, SIZE=bsize)
    end subroutine

end module
