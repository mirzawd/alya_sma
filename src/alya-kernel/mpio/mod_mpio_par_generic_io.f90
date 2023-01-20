!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup mpio
!> @{
!> @file    mod_mpio_par_generic_io.f90
!> @author  Damien Dosimont
!> @date    27/09/2017
!> @brief   MPI-IO generic
!> @details This module is a bridge that redirects to the sequential or the parallel
!>          versions of the MPI I/O operations
!> @}
!-----------------------------------------------------------------------

module mod_mpio_par_generic_io

    use def_master
    use def_mpio
    use def_domain,                     only : meshe
    use def_kermod,                     only : ndivi
    use mod_mpio_config,                only : mpio_config
    use mod_parall,                     only : PAR_CODE_SIZE
    use mod_communications,             only : PAR_BROADCAST, PAR_BARRIER, PAR_SUM
    use mod_memory,                     only : memory_alloca, memory_deallo
    use mod_mpio_par_io
    use mod_mpio_par_hybrid_io

    implicit none

    private

    character(150)                              :: wherein_code="IN MY CODE"

    interface MPIO_PAR_READ
        module procedure                        MPIO_PAR_READ_INT_V,&
                                                MPIO_PAR_READ_INT_M,&
                                                MPIO_PAR_READ_REAL_V,&
                                                MPIO_PAR_READ_REAL_M
    end interface

    public :: MPIO_PAR_READ

    contains

    subroutine MPIO_PAR_READ_INT_V(buf, filename, dim)
      integer(ip), pointer, intent(inout)         :: buf(:)
      character(*), intent(in)                    :: filename
      type(mpio_header)                           :: header
      integer(ip), intent(inout)                  :: dim
      logical(lg)                                 :: too_small
      call MPIO_COMPUTE_DIMENSIONS(filename, header, dim)
      too_small=is_too_small(dim, ip)
      call memory_alloca(mpio_memor,'buf', " MPIO_PAR_READ_INT_V", buf, max(dim,1_ip))
      if (too_small) then
        call PARSEQ_FILE_READ_ALL(buf, filename, dim, header)
      else
        call PAR_FILE_READ_ALL(buf, filename, dim, wherein_code, header)
      end if
    end subroutine

    subroutine MPIO_PAR_READ_INT_M(buf, filename, lines, columns)
      integer(ip), pointer, intent(inout)         :: buf(:,:)
      character(*), intent(in)                    :: filename
      type(mpio_header)                           :: header
      integer(ip), intent(inout)                  :: lines
      integer(ip), intent(inout)                  :: columns
      logical(lg)                                 :: too_small
      call MPIO_COMPUTE_DIMENSIONS(filename, header, lines, columns)
      too_small=is_too_small(lines*columns, ip)
      call memory_alloca(mpio_memor,'buf', " MPIO_PAR_READ_INT_M", buf, max(columns,1_ip), max(lines,1_ip))
      if (too_small) then
        call PARSEQ_FILE_READ_ALL(buf, filename, columns, lines, header)
      else
        call PAR_FILE_READ_ALL(buf, filename, columns, lines, wherein_code, header)
      end if
    end subroutine

    subroutine MPIO_PAR_READ_REAL_V(buf, filename, dim)
      real(rp), pointer, intent(inout)            :: buf(:)
      character(*), intent(in)                    :: filename
      type(mpio_header)                           :: header
      integer(ip), intent(inout)                  :: dim
      logical(lg)                                 :: too_small
      call MPIO_COMPUTE_DIMENSIONS(filename, header, dim)
      too_small=is_too_small(dim, rp)
      call memory_alloca(mpio_memor,'buf', " MPIO_PAR_READ_REAL_V", buf, max(dim,1_ip))
      if (too_small) then
        call PARSEQ_FILE_READ_ALL(buf, filename, dim, header)
      else
        call PAR_FILE_READ_ALL(buf, filename, dim, wherein_code, header)
      end if
    end subroutine

    subroutine MPIO_PAR_READ_REAL_M(buf, filename, lines, columns)
      real(rp), pointer, intent(inout)            :: buf(:,:)
      character(*),      intent(in)               :: filename
      type(mpio_header)                           :: header
      integer(ip)                                 :: lines
      integer(ip)                                 :: columns
      logical(lg)                                 :: too_small
      call MPIO_COMPUTE_DIMENSIONS(filename, header, lines, columns)
      too_small=is_too_small(lines*columns, rp)
      call memory_alloca(mpio_memor,'buf', " MPIO_PAR_READ_REAL_M", buf, max(columns,1_ip), max(lines,1_ip))
      if (too_small) then
        call PARSEQ_FILE_READ_ALL(buf, filename, columns, lines, header)
      else
        call PAR_FILE_READ_ALL(buf, filename, columns, lines, wherein_code, header)
      end if
    end subroutine

    subroutine MPIO_COMPUTE_DIMENSIONS(filename, header, lines, columns)
      character(*),      intent(in)               :: filename
      type(mpio_header), intent(out)              :: header
      integer(ip),       intent(out)              :: lines
      integer(ip),optional,intent(out)            :: columns
      character(8)                                :: resultson
      integer(ip)                                 :: imesh
      call PAR_FILE_READ_HEADER(filename, header, wherein_code)
      resultson=header%resultson
      columns=header%columns
      imesh = ndivi
      if( resultson(1:5) == 'NPOIN' ) then
          lines=meshe(imesh) % npoin
      else if( resultson(1:5) == 'NELEM' ) then
          lines=meshe(imesh) % nelem
      else if( resultson(1:5) == 'NBOUN' ) then
          lines=meshe(imesh) % nboun
      else if( resultson(1:5) == 'WHATE' ) then
          call runend("The result type "//resultson(1:5)//" can not be processed by the MPI-IO module...")
      else
          call runend("The result type "//resultson(1:5)//" can not be processed by the MPI-IO module...")
      end if
    end subroutine MPIO_COMPUTE_DIMENSIONS

    function IS_TOO_SMALL(dims, type_size) result(too_small)
      integer(ip), intent(in)    ::  dims
      logical(lg)                ::  too_small
      integer(ip)                ::  dimstmp
      integer(ip)                ::  avgsize
      integer, intent(in)        ::  type_size
      if (mpio_config%output%par_min_block_size .le. 0_ip) then
         too_small=.false.
      else
         dimstmp=dims
         if (IMASTER) then
            dimstmp=0
         end if
         call PAR_SUM(dimstmp, wherein_code)
         avgsize=int((real(type_size,rp)*real(dimstmp,rp))/((1024.0_rp*1024.0_rp*(real(PAR_CODE_SIZE,rp)-1.0_rp))), ip)
         call PAR_BROADCAST(avgsize, wherein=wherein_code)
         if (avgsize < mpio_config%output%par_min_block_size) then
            too_small=.true.
         else
            too_small=.false.
         end if
      end if
      call PAR_BROADCAST(too_small, wherein=wherein_code)
    end function is_too_small


end module
