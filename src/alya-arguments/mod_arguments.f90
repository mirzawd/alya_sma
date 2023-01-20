!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!

!-----------------------------------------------------------------------
!> @addtogroup Master
!> @{
!> @file    mod_arguments.f90
!> @author  houzeaux
!> @date    2022-11-27
!> @brief   Get options
!> @details Get Alya options by pasring command line
!-----------------------------------------------------------------------

module mod_arguments

    use def_kintyp_basic, only: ip
    use mod_strings, only: integer_to_string
    implicit none
    private

    type type_options
        character(1)   :: shortopt = ''   ! Short text for the option
        character(30)  :: longopt = ''   ! Full text for the option
        integer(ip)    :: num_arguments = 0_ip ! Number of arguments of the option
        character(50)  :: description = ''   ! Description of the option
    end type type_options

    public :: type_options
    public :: arguments_help
    public :: arguments_number
    public :: arguments_get
    public :: arguments_find

contains

    !-----------------------------------------------------------------------
    !>
    !> @author  houzeaux
    !> @date    2020-05-11
    !> @brief   Help message
    !> @details Write help message
    !>
    !-----------------------------------------------------------------------

    subroutine arguments_help(options)

        type(type_options), intent(in) :: options(:)
        character(30)                  :: wlong
        character(50)                  :: wdescr
        character(20)                  :: wopt(size(options))
        integer(ip)                    :: ii, kk, mm, nn

        mm = min(maxval(len_trim(options(:)%longopt)) + 1, len(wlong))
        do ii = 1, size(options)
            wopt(ii) = ''
            do kk = 1, options(ii)%num_arguments
                wopt(ii) = trim(wopt(ii))//' arg'//integer_to_string(kk)
            end do
        end do
        nn = min(maxval(len_trim(wopt)), len(wopt))

        write (6, 1)
        do ii = 1, size(options)
            wlong = adjustl('--'//adjustl(trim(options(ii)%longopt)))
            wdescr = adjustl(options(ii)%description)
            write (6, '(a,a,a,a)') &
                '    -'//trim(options(ii)%shortopt), &
                '  '//adjustl(wlong(1:mm)), &
                '  '//wopt(ii) (1:nn), &
                '  '//wdescr
        end do

        write (6, 3)
        stop

1       format('   ', /, &
               '   ALYA USAGE:', /, &
               '   ', /, &
               '   Alya.x [options] case ',/)
3       format('   ', /, &
               '   ', /, &
               '   Runs Alya for problem case. ', /, &
               '   ', /, &
               '   The following I/O files are located/created in current directory (mod is any activated module extension)', /, &
               '   * means optional:', /, &
               '   ', /, &
               '   (I)    case.dat:                     run data', /, &
               '   (I)    case.ker.dat:                 kernel data', /, &
               '   (I)    case.dom.dat:                 mesh data', /, &
               '   (I*)   case.cou.dat:                 coupling data', /, &
               '   (I)    case.mod.dat:                 module data', /, &
               '   ', /, &
               '   (O)    case.log:                     run log', /, &
               '   (O)    case.cvg:                     global convergencesof modules', /, &
               '   (O)    case.ker.log:                 kernel log', /, &
               '   (O)    case-memory.res               memory tracking', /, &
               '   (O)    case-performance.csv          timings and performance metrics', /, &
               '   (O*)   case.liv:                     live info', /, &
               '   (O*)   case-partition.par.post.msh   partition mesh in GiD format', /, &
               '   (O*)   case-partition.par.post.res   partition results in GiD format', /, &
               '   (O*)   case-VAR-X.post.alyabin       postprocess of variable VAR at time step X in alya format', /, &
               '   (O*)   case-VAR-X.post.mpio.bin      postprocess of variable VAR at time step X in mpio format', /, &
               '   (O)    case-VAR.mod.sol              solver information for variable VAR', /, &
               '   (O*)   case-VAR.mod.cso              solver convergence for variable VAR', /, &
               '   (O)    case.mod.cvg:                 module convergence', /, &
               '   (O)    case.mod.rst:                 module restart', /, &
               '   (O*)   case.mod.wit:                 module witness point results', /, &
               '   (O*)   case-VAR-X-Y.nsi.rst          restart file for VAR at time step X, component Y, in alya format', /, &
               '   (O*)   case-VAR-X-Y.nsi.rst.mpio.bin restart file for VAR at time step X, component Y, in mpio format', /, &
               '   (O*)   case-node.mod.set             module node set results', /, &
               '   (O*)   case-element.mod.set          module element set results', /, &
               '   (O*)   case-boundary.mod.set         module boundary set results', /, &
               '   ')

    end subroutine arguments_help

    !-----------------------------------------------------------------------
    !>
    !> @author  houzeaux
    !> @date    2022-10-27
    !> @brief   Number of arguments
    !> @details Number of arguments
    !>
    !-----------------------------------------------------------------------

    integer(ip) function arguments_number()

        arguments_number = command_argument_count()

    end function arguments_number

    !-----------------------------------------------------------------------
    !>
    !> @author  houzeaux
    !> @date    2022-10-27
    !> @brief   Get arguments
    !> @details Get arguments
    !>
    !-----------------------------------------------------------------------

    function arguments_get(iopt) result(opt)

        integer(ip), intent(in)  :: iopt
        character(100)                :: opt

        call get_command_argument(int(iopt, 4), opt)

    end function arguments_get

    !-----------------------------------------------------------------------
    !>
    !> @author  houzeaux
    !> @date    2022-10-27
    !> @brief   Find option
    !> @details Find the option in options
    !>
    !-----------------------------------------------------------------------

    function arguments_find(options, opt) result(iopt)

        type(type_options), intent(in) :: options(:)
        character(*), intent(in) :: opt
        integer(ip)                    :: iopt
        integer(ip)                    :: ii, kk

        iopt = 0
        kk = 0

        if (len_trim(opt) >= 2) then

            do ii = 1, size(options)

                if (opt(1:2) == '--') then
                    kk = 1
                    if (len_trim(opt) >= 3) then
                        if (trim(opt(3:)) == trim(options(ii)%longopt)) then
                            iopt = ii
                            return
                        end if
                    end if
                else if (opt(1:1) == '-') then
                    kk = 1
                    if (trim(opt(2:)) == trim(options(ii)%shortopt)) then
                        iopt = ii
                        return
                    end if
                end if

            end do

        end if

        if (kk == 1) iopt = -1

    end function arguments_find

end module mod_arguments
!> @}
