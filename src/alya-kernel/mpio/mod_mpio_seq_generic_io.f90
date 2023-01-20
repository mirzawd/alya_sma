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

module mod_mpio_seq_generic_io

    use def_master
    use def_mpio
    use def_domain,                     only : meshe
    use def_kermod,                     only : ndivi
    use mod_parall,                     only : PAR_CODE_SIZE
    use mod_communications,             only : PAR_BROADCAST, PAR_BARRIER, PAR_SUM
    use mod_memory,                     only : memory_alloca, memory_deallo
    use mod_mpio_seq_io

    implicit none

    private

    interface MPIO_SEQ_READ
        module procedure                        MPIO_SEQ_READ_INT_V,&
                                                MPIO_SEQ_READ_INT_M,&
                                                MPIO_SEQ_READ_REAL_V,&
                                                MPIO_SEQ_READ_REAL_M
    end interface

    public :: MPIO_SEQ_READ

    contains

    subroutine MPIO_SEQ_READ_INT_V(buf, filename, dim)
      integer(ip), pointer, intent(inout)         :: buf(:)
      character(*), intent(in)                    :: filename
      type(mpio_header)                           :: header
      integer(ip), intent(inout)                  :: dim
      call MPIO_COMPUTE_DIMENSIONS(filename, header, dim)
      call memory_alloca(mpio_memor,'buf', " MPIO_SEQ_READ_INT_V", buf, dim)
      call SEQ_FILE_READ_ALL(buf, filename, dim, header)
    end subroutine

    subroutine MPIO_SEQ_READ_INT_M(buf, filename, lines, columns)
      integer(ip), pointer, intent(inout)         :: buf(:,:)
      character(*), intent(in)                    :: filename
      type(mpio_header)                           :: header
      integer(ip), intent(inout)                  :: lines
      integer(ip), intent(inout)                  :: columns
      call MPIO_COMPUTE_DIMENSIONS(filename, header, lines, columns)
      call memory_alloca(mpio_memor,'buf', " MPIO_SEQ_READ_INT_M", buf, columns, lines)
      call SEQ_FILE_READ_ALL(buf, filename, columns, lines, header)
    end subroutine

    subroutine MPIO_SEQ_READ_REAL_V(buf, filename, dim)
      real(rp), pointer, intent(inout)            :: buf(:)
      character(*), intent(in)                    :: filename
      type(mpio_header)                           :: header
      integer(ip), intent(inout)                  :: dim
      call MPIO_COMPUTE_DIMENSIONS(filename, header, dim)
      call memory_alloca(mpio_memor,'buf', " MPIO_SEQ_READ_REAL_V", buf, dim)
      call SEQ_FILE_READ_ALL(buf, filename, dim, header)
    end subroutine

    subroutine MPIO_SEQ_READ_REAL_M(buf, filename, lines, columns)
      real(rp), pointer, intent(inout)            :: buf(:,:)
      character(*),      intent(in)               :: filename
      type(mpio_header)                           :: header
      integer(ip)                                 :: lines
      integer(ip)                                 :: columns
      call MPIO_COMPUTE_DIMENSIONS(filename, header, lines, columns)
      call memory_alloca(mpio_memor,'buf', " MPIO_SEQ_READ_REAL_M", buf, columns, lines)
      call SEQ_FILE_READ_ALL(buf, filename, columns, lines, header)
    end subroutine

    subroutine MPIO_COMPUTE_DIMENSIONS(filename, header, lines, columns)
      character(*),      intent(in)               :: filename
      type(mpio_header), intent(out)              :: header
      integer(ip),       intent(out)              :: lines
      integer(ip),optional,intent(out)            :: columns
      character(8)                                :: resultson
      integer(ip)                                 :: imesh
      call SEQ_FILE_READ_HEADER(filename, header)
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

end module
