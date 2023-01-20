!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup mpio
!> @{
!> @file    mod_mpio_par_async_io.f90
!> @author  Damien Dosimont
!> @date    27/09/2017
!> @brief   Asynchronous IO manager
!> @details This module provides asynchronous IO management
!> @}
!-----------------------------------------------------------------------

module mod_mpio_par_async_io

    use def_master,             only : IIOSLAVE, IMASTER, INOTMASTER, IIOSLAVE
    use def_master,             only : namda
    use def_master,             only : intost, retost
    use mod_memory,             only : memory_alloca, memory_deallo
    use mod_parall,             only : PAR_MY_CODE_RANK
    use def_kintyp,             only : ip, rp, lg
    use def_mpio,               only : mpio_memor
    use mod_mpio_par_mpiwrapper
    use mod_messages,          only : livinf
    use def_mpi
#include "def_mpi.inc"

    implicit none

    Type asyncfile_int_v
       MY_MPI_FILE                                      :: fh
       integer(ip), pointer                             :: buf(:)
       Type(asyncfile_int_v), pointer                   :: next
    end Type asyncfile_int_v
    
    Type asyncfile_int_m
       MY_MPI_FILE                                      :: fh
       integer(ip), pointer                             :: buf(:,:)
       Type(asyncfile_int_m), pointer                   :: next
    end Type asyncfile_int_m
    
    Type asyncfile_real_v
       MY_MPI_FILE                                      :: fh
       real(rp), pointer                                :: buf(:)
       Type(asyncfile_real_v), pointer                  :: next
    end Type asyncfile_real_v
    
    Type asyncfile_real_m
       MY_MPI_FILE                                      :: fh
       real(rp), pointer                                :: buf(:,:)
       Type(asyncfile_real_m), pointer                  :: next
    end Type asyncfile_real_m
    
    interface PAR_ASYNCHRONOUS_WRITE_PUSH
        module procedure                        PAR_ASYNCHRONOUS_WRITE_PUSH_INT_V,&
                                                PAR_ASYNCHRONOUS_WRITE_PUSH_INT_M,&
                                                PAR_ASYNCHRONOUS_WRITE_PUSH_REAL_V,&
                                                PAR_ASYNCHRONOUS_WRITE_PUSH_REAL_M
    end interface

    interface PAR_ASYNCHRONOUS_WRITE_POP
        module procedure                        PAR_ASYNCHRONOUS_WRITE_POP_INT_V,&
                                                PAR_ASYNCHRONOUS_WRITE_POP_INT_M,&
                                                PAR_ASYNCHRONOUS_WRITE_POP_REAL_V,&
                                                PAR_ASYNCHRONOUS_WRITE_POP_REAL_M
    end interface

    interface ADD_ELEMENT
        module procedure                        ADD_ELEMENT_INT_V,&
                                                ADD_ELEMENT_INT_M,&
                                                ADD_ELEMENT_REAL_V,&
                                                ADD_ELEMENT_REAL_M
    end interface

    type(asyncfile_int_v), pointer                     :: asyncfile_int_v_stack
    type(asyncfile_int_m), pointer                     :: asyncfile_int_m_stack
    type(asyncfile_real_v), pointer                    :: asyncfile_real_v_stack
    type(asyncfile_real_m), pointer                    :: asyncfile_real_m_stack

    integer(ip)                                        :: count_int_v=0
    integer(ip)                                        :: count_int_m=0
    integer(ip)                                        :: count_real_v=0
    integer(ip)                                        :: count_real_m=0

    integer(ip)                                        :: previous_ittim=0

    public                                  ::  PAR_ASYNCHRONOUS_WRITE_PUSH, PAR_ASYNCHRONOUS_WRITE_POP_ALL


    contains

    subroutine PAR_ASYNCHRONOUS_WRITE_PUSH_INT_V(fh, buf, size)
        MY_MPI_FILE   ,        intent(in)              :: fh
        integer(ip), pointer,  intent(in)              :: buf(:)
        integer(ip),           intent(in)              :: size
        integer(ip), pointer                           :: cpbuf(:)
        nullify(cpbuf)
#ifndef MPI_OFF
        if (size==0 .or. .not.associated(buf)) then
            call memory_alloca(mpio_memor,'cpbuf','PAR_ASYNCHRONOUS_WRITE_PUSH_INT_V',  cpbuf, 1_ip)
            call PAR_FILE_WRITE_ALL_BEGIN(fh, cpbuf, size)
        else
            call memory_alloca(mpio_memor,'cpbuf','PAR_ASYNCHRONOUS_WRITE_PUSH_INT_V',  cpbuf, size)
            cpbuf=buf(1:size)
            call PAR_FILE_WRITE_ALL_BEGIN(fh, cpbuf, size)
        end if
        call ADD_ELEMENT(fh, cpbuf)
#endif
    end subroutine

    subroutine PAR_ASYNCHRONOUS_WRITE_PUSH_INT_M(fh, buf, dim1, dim2)
        MY_MPI_FILE   ,        intent(in)              :: fh
        integer(ip), pointer,  intent(in)              :: buf(:,:)
        integer(ip),           intent(in)              :: dim1
        integer(ip),           intent(in)              :: dim2
        integer(ip), pointer                           :: cpbuf(:,:)
        nullify(cpbuf)
#ifndef MPI_OFF
        if (dim1*dim2==0 .or. .not.associated(buf)) then
            call memory_alloca(mpio_memor,'cpbuf','PAR_ASYNCHRONOUS_WRITE_PUSH_INT_M',  cpbuf, 1_ip, 1_ip)
            call PAR_FILE_WRITE_ALL_BEGIN(fh, cpbuf, 0_ip)
        else
            call memory_alloca(mpio_memor,'cpbuf','PAR_ASYNCHRONOUS_WRITE_PUSH_INT_M',  cpbuf, dim1, dim2)
            cpbuf=buf(1:dim1, 1:dim2)
            call PAR_FILE_WRITE_ALL_BEGIN(fh, cpbuf, dim1*dim2)
        end if
        call ADD_ELEMENT(fh, cpbuf)
#endif
    end subroutine


    subroutine PAR_ASYNCHRONOUS_WRITE_PUSH_REAL_V(fh, buf, size)
        MY_MPI_FILE   ,        intent(in)              :: fh
        real(rp), pointer,     intent(in)             :: buf(:)
        integer(ip),           intent(in)              :: size
        real(rp), pointer                              :: cpbuf(:)
        nullify(cpbuf)
#ifndef MPI_OFF
        if (size==0 .or. .not.associated(buf)) then
            call memory_alloca(mpio_memor,'cpbuf','PAR_ASYNCHRONOUS_WRITE_PUSH_REAL_V',  cpbuf, 1_ip)
            call PAR_FILE_WRITE_ALL_BEGIN(fh, cpbuf, 0_ip)
        else
            call memory_alloca(mpio_memor,'cpbuf','PAR_ASYNCHRONOUS_WRITE_PUSH_REAL_V',  cpbuf, size)
            cpbuf=buf(1:size)
            call PAR_FILE_WRITE_ALL_BEGIN(fh, cpbuf, size)
        end if
        call ADD_ELEMENT(fh, cpbuf)
#endif
    end subroutine

    subroutine PAR_ASYNCHRONOUS_WRITE_PUSH_REAL_M(fh, buf, dim1, dim2)
        MY_MPI_FILE   ,          intent(in)              :: fh
        real(rp), pointer,  intent(in)                 :: buf(:,:)
        integer(ip),           intent(in)              :: dim1
        integer(ip),           intent(in)              :: dim2
        real(rp), pointer                              :: cpbuf(:,:)
        nullify(cpbuf)
#ifndef MPI_OFF
        if (dim1*dim2==0 .or. .not.associated(buf)) then
            call memory_alloca(mpio_memor,'CPBUF','PAR_ASYNCHRONOUS_WRITE_PUSH_REAL_M',  cpbuf, 1_ip, 1_ip)
            call PAR_FILE_WRITE_ALL_BEGIN(fh, cpbuf, 0_ip)
        else
            call memory_alloca(mpio_memor,'CPBUF','PAR_ASYNCHRONOUS_WRITE_PUSH_REAL_M',  cpbuf, dim1, dim2)
            cpbuf=buf(1:dim1, 1:dim2)
            call PAR_FILE_WRITE_ALL_BEGIN(fh, cpbuf, dim1*dim2)
        end if
        call ADD_ELEMENT(fh, cpbuf)
#endif
    end subroutine

    subroutine PAR_ASYNCHRONOUS_WRITE_POP_INT_V(element)
        type(asyncfile_int_v), pointer                :: element
#ifndef MPI_OFF
        call PAR_FILE_WRITE_ALL_END(element%fh, element%buf)
        call PAR_FILE_CLOSE(element%fh)
#endif
    end subroutine

    subroutine PAR_ASYNCHRONOUS_WRITE_POP_INT_M(element)
        type(asyncfile_int_m), pointer                :: element
#ifndef MPI_OFF
        call PAR_FILE_WRITE_ALL_END(element%fh, element%buf)
        call PAR_FILE_CLOSE(element%fh)
#endif
    end subroutine

    subroutine PAR_ASYNCHRONOUS_WRITE_POP_REAL_V(element)
        type(asyncfile_real_v), pointer                :: element
#ifndef MPI_OFF
        call PAR_FILE_WRITE_ALL_END(element%fh, element%buf)
        call PAR_FILE_CLOSE(element%fh)
#endif
    end subroutine

    subroutine PAR_ASYNCHRONOUS_WRITE_POP_REAL_M(element)
        type(asyncfile_real_m), pointer                :: element
#ifndef MPI_OFF
        call PAR_FILE_WRITE_ALL_END(element%fh, element%buf)
        call PAR_FILE_CLOSE(element%fh)
#endif
    end subroutine

    subroutine ADD_ELEMENT_INT_V(fh, buf)
        MY_MPI_FILE   ,            intent(in)              :: fh
        integer(ip), pointer,   intent(in)              :: buf(:)
        type(asyncfile_int_v), pointer                :: current
        type(asyncfile_int_v), pointer                :: pre
#ifndef MPI_OFF
        if (count_int_v==0) then
            nullify(asyncfile_int_v_stack)
            allocate(asyncfile_int_v_stack)
            asyncfile_int_v_stack%fh=fh
            asyncfile_int_v_stack%buf=>buf
            nullify(asyncfile_int_v_stack%next)
        else
            nullify(current)
            allocate(current)
            current%fh=fh
            current%buf=>buf
            nullify(current%next)
            pre=>asyncfile_int_v_stack
            do while (associated(pre%next))
                pre=>pre%next
            end do
            pre%next=>current
        end if
        count_int_v=count_int_v+1

#endif
    end subroutine

    subroutine ADD_ELEMENT_INT_M(fh, buf)
        MY_MPI_FILE   ,           intent(in)              :: fh
        integer(ip), pointer,  intent(in)              :: buf(:,:)
        type(asyncfile_int_m), pointer                :: current
        type(asyncfile_int_m), pointer                :: pre
#ifndef MPI_OFF
        if (count_int_m==0) then
            nullify(asyncfile_int_m_stack)
            allocate(asyncfile_int_m_stack)
            asyncfile_int_m_stack%fh=fh
            asyncfile_int_m_stack%buf=>buf
            nullify(asyncfile_int_m_stack%next)
        else
            nullify(current)
            allocate(current)
            current%fh=fh
            current%buf=>buf
            nullify(current%next)
            pre=>asyncfile_int_m_stack
            do while (associated(pre%next))
                pre=>pre%next
            end do
            pre%next=>current
        end if
        count_int_m=count_int_m+1
#endif
    end subroutine

    subroutine ADD_ELEMENT_REAL_V(fh, buf)
        MY_MPI_FILE   ,        intent(in)               :: fh
        real(rp), pointer,     intent(in)               :: buf(:)
        type(asyncfile_real_v), pointer                :: current
        type(asyncfile_real_v), pointer                :: pre
#ifndef MPI_OFF
        if (count_real_v==0) then
            nullify(asyncfile_real_v_stack)
            allocate(asyncfile_real_v_stack)
            asyncfile_real_v_stack%fh=fh
            asyncfile_real_v_stack%buf=>buf
            nullify(asyncfile_real_v_stack%next)
        else
            nullify(current)
            allocate(current)
            current%fh=fh
            current%buf=>buf
            nullify(current%next)
            pre=>asyncfile_real_v_stack
            do while (associated(pre%next))
                pre=>pre%next
            end do
            pre%next=>current
        end if
        count_real_v=count_real_v+1
#endif
    end subroutine

    subroutine ADD_ELEMENT_REAL_M(fh, buf)
        MY_MPI_FILE   ,       intent(in)               :: fh
        real(rp), pointer,  intent(in)                 :: buf(:,:)
        type(asyncfile_real_m), pointer                :: current
        type(asyncfile_real_m), pointer                :: pre
        
#ifndef MPI_OFF
        if (count_real_m==0) then
            nullify(asyncfile_real_m_stack)
            allocate(asyncfile_real_m_stack)
            asyncfile_real_m_stack%fh=fh
            asyncfile_real_m_stack%buf=>buf
            nullify(asyncfile_real_m_stack%next)
        else
            nullify(current)
            allocate(current)
            current%fh=fh
            current%buf=>buf
            nullify(current%next)
            pre=>asyncfile_real_m_stack
            do while (associated(pre%next))
                pre=>pre%next
            end do
            pre%next=>current
        end if
        count_real_m=count_real_m+1
#endif
    end subroutine

    subroutine RM_LAST_ELEMENT_INT_V()
        type(asyncfile_int_v), pointer                              :: element
        type(asyncfile_int_v), pointer                              :: pre
#ifndef MPI_OFF
        count_int_v=count_int_v-1
        element=>asyncfile_int_v_stack
        pre=>element
        do while(associated(element%next))
            pre=>element
            element=>element%next
        end do
        call memory_deallo(mpio_memor,'cpbuf','RM_LAST_ELEMENT_INT_V', element%buf)
        if (associated(pre%next)) then
            deallocate(pre%next)
        else
            deallocate(pre)
        end if
#endif
    end subroutine

    subroutine RM_LAST_ELEMENT_INT_M()
        type(asyncfile_int_m), pointer                              :: element
        type(asyncfile_int_m), pointer                              :: pre
#ifndef MPI_OFF
        count_int_m=count_int_m-1
        element=>asyncfile_int_m_stack
        pre=>element
        do while(associated(element%next))
            pre=>element
            element=>element%next
        end do
        call memory_deallo(mpio_memor,'cpbuf','RM_LAST_ELEMENT_INT_M', element%buf)
        if (associated(pre%next)) then
            deallocate(pre%next)
        else
            deallocate(pre)
        end if
#endif
    end subroutine

    subroutine RM_LAST_ELEMENT_REAL_V()
        type(asyncfile_real_v), pointer                              :: element
        type(asyncfile_real_v), pointer                              :: pre
#ifndef MPI_OFF
        count_real_v=count_real_v-1
        element=>asyncfile_real_v_stack
        pre=>element
        do while(associated(element%next))
            pre=>element
            element=>element%next
        end do
        call memory_deallo(mpio_memor,'cpbuf','RM_LAST_ELEMENT_REAL_V', element%buf)
        if (associated(pre%next)) then
            deallocate(pre%next)
        else
            deallocate(pre)
        end if
#endif
    end subroutine

    subroutine RM_LAST_ELEMENT_REAL_M()
        type(asyncfile_real_m), pointer                              :: element
        type(asyncfile_real_m), pointer                              :: pre
#ifndef MPI_OFF
        element=>asyncfile_real_m_stack
        count_real_m=count_real_m-1
        pre=>element
        do while(associated(element%next))
            pre=>element
            element=>element%next
        end do
        call memory_deallo(mpio_memor,'cpbuf','RM_LAST_ELEMENT_REAL_M', element%buf)
        if (associated(pre%next)) then
            deallocate(pre%next)
        else
            deallocate(pre)
        end if
#endif
    end subroutine

    subroutine PAR_ASYNCHRONOUS_WRITE_POP_ALL(threshold, ittim)
        integer(ip),               intent(in)                        :: threshold
        integer(ip),               intent(in), optional              :: ittim
        type(asyncfile_int_v), pointer                              :: element_int_v
        type(asyncfile_int_m), pointer                              :: element_int_m
        type(asyncfile_real_v), pointer                             :: element_real_v
        type(asyncfile_real_m), pointer                             :: element_real_m

#ifndef MPI_OFF
        if (threshold<0) then
            if (present(ittim)) then
                if (previous_ittim==ittim) then
                    return
                else
                    previous_ittim=ittim
                end if
            end if
        end if
        if (threshold<(count_int_v+count_int_m+count_real_v+count_real_m)) then
            call livinf(0_ip,'FLUSH ASYNCHRONOUS MPI-IO WRITE BUFFER',0_ip)
            do while(count_int_v>0)
                element_int_v=>asyncfile_int_v_stack
                do while (associated(element_int_v%next))
                    element_int_v=>element_int_v%next
                end do
                call PAR_ASYNCHRONOUS_WRITE_POP(element_int_v)
                call RM_LAST_ELEMENT_int_V()
            end do
            do while(count_int_m>0)
                element_int_m=>asyncfile_int_m_stack
                do while (associated(element_int_m%next))
                    element_int_m=>element_int_m%next
                end do
                call PAR_ASYNCHRONOUS_WRITE_POP(element_int_m)
                call RM_LAST_ELEMENT_int_M()
            end do
            do while(count_real_v>0)
                element_real_v=>asyncfile_real_v_stack
                do while (associated(element_real_v%next))
                    element_real_v=>element_real_v%next
                end do
                call PAR_ASYNCHRONOUS_WRITE_POP(element_real_v)
                call RM_LAST_ELEMENT_real_V()
            end do
            do while(count_real_m>0)
                element_real_m=>asyncfile_real_m_stack
                do while (associated(element_real_m%next))
                    element_real_m=>element_real_m%next
                end do
                call PAR_ASYNCHRONOUS_WRITE_POP(element_real_m)
                call RM_LAST_ELEMENT_real_M()
            end do
        end if
#endif
    end subroutine

end module
