!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!

module def_mpio

    use def_kintyp, only: ip, rp, lg, r1p
#ifndef MPI_OFF
    use def_mpi
#include "def_mpi.inc"
#endif

    !-----------------------------------------------------------------------------------------------------------------------------
    !                              Definitions
    !-----------------------------------------------------------------------------------------------------------------------------

    character(9), parameter                 :: mpio_ext = '.mpio.bin'           ! MPIO file format extension

    character(9), parameter                 :: post_ext = '.post'               ! Post format extension

    character(6), parameter                 :: coord_ext = '-COORD', &
                                               ltype_ext = '-LTYPE', &
                                               lnods_ext = '-LNODS', &
                                               ltypb_ext = '-LTYPB', &
                                               lelbo_ext = '-LELBO', &
                                               lnodb_ext = '-LNODB', &
                                               leinv_ext = '-LEINV', &
                                               lninv_ext = '-LNINV', &
                                               lbinv_ext = '-LBINV', &
                                               lnoch_ext = '-LNOCH', &
                                               lelch_ext = '-LELCH', &
                                               lboch_ext = '-LBOCH', &
                                               lmate_ext = '-LMATE', &
                                               lmast_ext = '-LMAST', &
                                               lesub_ext = '-LESUB', &
                                               codbo_ext = '-CODBO', &
                                               codno_ext = '-CODNO', &
                                               leset_ext = '-LESET', &
                                               lbset_ext = '-LBSET', &
                                               lnset_ext = '-LNSET', &
                                               field_ext = '-XFIEL'

    integer(ip), parameter                  :: header_size = 256                   ! Header size of MPI-IO file

    integer(ip), parameter                  :: option_size = 10                    ! Option size

    integer(ip), parameter                  :: value_count = 1                     ! Count of an integer or real value

    integer(ip), parameter                  :: string_count = 8                    ! Count of a string

    integer, parameter                      :: header_ip = 8

    integer, parameter                      :: header_rp = 8

    integer(8), parameter                   :: header_magic_number = 27093

    character(3), parameter                 :: c5_f = '00'//char(0)
    character(1), parameter                 :: c7_f = char(0)

    character(8), parameter                 :: header_format = 'MPIAL'//c5_f, &
                                               header_version = 'V0004'//c5_f, &
                                               header_align_chars = '00000'//c5_f, &
                                               header_option_init = 'NONE0'//c5_f, &
                                               header_no_filter = 'NOFIL'//c5_f, &
                                               header_asc_sorting = 'ASCEN'//c5_f, &
                                               header_no_id = 'NOID0'//c5_f
#ifndef MPI_OFF
    integer                                 ::  ierr                                ! error flag
    MY_MPI_INFO                             ::  info                                ! info flag
#endif

    type mpio_header_options
        character(string_count), dimension(option_size) :: opt = (/header_option_init, &
                                                                   header_option_init, &
                                                                   header_option_init, &
                                                                   header_option_init, &
                                                                   header_option_init, &
                                                                   header_option_init, &
                                                                   header_option_init, &
                                                                   header_option_init, &
                                                                   header_option_init, &
                                                                   header_option_init/)
    end type mpio_header_options

    type mpio_header
        ! fields contained in the binary file header
        integer(header_ip)                          :: magic_number = header_magic_number
        character(string_count)                     :: format = header_format
        character(string_count)                     :: version = header_version
        character(string_count)                     :: object
        character(string_count)                     :: dimension
        character(string_count)                     :: resultson
        character(string_count)                     :: type
        character(string_count)                     :: size
        character(string_count)                     :: par
        character(string_count)                     :: filter = header_no_filter
        character(string_count)                     :: sorting = header_asc_sorting
        character(string_count)                     :: id = header_no_id
        character(string_count)                     :: align_chars = header_align_chars
        integer(header_ip)                          :: columns
        integer(header_ip)                          :: lines
        integer(header_ip)                          :: ittim
        integer(header_ip)                          :: nsubd
        integer(header_ip)                          :: divi
        integer(header_ip)                          :: tag1
        integer(header_ip)                          :: tag2
        real(header_rp)                             :: time
        type(mpio_header_options)                   :: options
        ! fields that are not exported
        integer(ip)                                 :: item_size
        integer(8)                                  :: file_size
    end type mpio_header

    integer(8)                                      :: mpio_memor(2)

    public                                          :: mpio_ext, &
                                                       coord_ext, &
                                                       ltype_ext, &
                                                       lnods_ext, &
                                                       ltypb_ext, &
                                                       lelbo_ext, &
                                                       lnodb_ext, &
                                                       leinv_ext, &
                                                       lninv_ext, &
                                                       lbinv_ext, &
                                                       lnoch_ext, &
                                                       lelch_ext, &
                                                       lboch_ext, &
                                                       lmate_ext, &
                                                       codbo_ext, &
                                                       codno_ext, &
                                                       leset_ext, &
                                                       lbset_ext, &
                                                       lnset_ext, &
                                                       field_ext, &
                                                       header_size, &
                                                       option_size, &
                                                       value_count, &
                                                       string_count, &
                                                       header_magic_number, &
                                                       header_format, &
                                                       header_version, &
                                                       header_align_chars, &
                                                       mpio_header, &
                                                       mpio_header_options, &
                                                       c5_f, &
                                                       c7_f, &
                                                       header_ip, &
                                                       header_rp, &
                                                       mpio_memor

end module def_mpio
