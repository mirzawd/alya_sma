!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!

!-----------------------------------------------------------------------
!> @addtogroup Maths
!> @{
!> @file    mod_random.f90
!> @author  houzeaux
!> @date    2020-05-11
!> @brief   Random number generator
!> @details Random number generator
!>
!-----------------------------------------------------------------------

module mod_random

  use def_kintyp_basic, only: rp, ip, lg
  use def_mpi
#include "def_mpi.inc"
  use mod_communications,    only : PAR_BROADCAST
  use def_master,            only : ISEQUEN, IMASTER, IPARALL, INOTSLAVE
  use mod_parall,            only : PAR_MY_CODE_RANK, PAR_CODE_SIZE
  use mod_optional_argument, only : optional_argument
  use mod_iofile_basic,      only : iofile_available_unit
  use mod_iofile,            only : iofile_open_unit
  use mod_iofile,            only : iofile_close_unit
  implicit none

  logical(lg)        :: opt_random_from_file ! Random numbers are read from file
  logical(lg)        :: opt_broadcast_seed   ! Seed should be broadcast
  integer(ip)        :: nrand                ! Number of random numbers in file
  integer(ip)        :: irand                ! Current random number
  integer(ip)        :: urand=0              ! Random file unit
  character(LEN=256) :: my_iomsg
  
  public :: random_generate_number
  public :: random_initialization
  public :: random_end

contains

  !-----------------------------------------------------------------------
  !>
  !> @author  calmet
  !> @date    2022-10-26
  !> @brief   Seed generator
  !> @details When a unique seed is used, only MPI rank=0 generates the
  !>          seed and broadcast it to antoher ranks and used the same
  !>          random numbers. If not each rank generate its own seed.
  !>
  !>          1. To have a fixed seed for reproducibility over several runs:
  !>             call random_initialization(.true.)
  !>          2. To have a unique seed over the MPIs:
  !>             call random_initialization(.false.,.true.)
  !>          3. To have a different seed over all MPIs:
  !>             call random_initialization(.false.,.false.)
  !>
  !>          Last mode is useless because a fixed seed does not need
  !>          to be broadcasted
  !>
  !-----------------------------------------------------------------------

  subroutine random_initialization(broadcast_seed, PAR_COMM_IN, random_from_file)

    logical(lg), intent(in), optional  :: broadcast_seed
    MY_MPI_COMM, intent(in), optional  :: PAR_COMM_IN
    logical(lg), intent(in), optional  :: random_from_file
    integer(4)                         :: k,j
    integer(ip)                        :: i,clock
    integer(4),  pointer               :: seed(:)
    integer                            :: values(8)
    real(rp)                           :: r,s

    if( urand == 0 ) then
       !
       ! Fixed seed: should always generate the same number sequence
       !
       opt_random_from_file = optional_argument(.false., random_from_file)
       !
       ! Broadcast seed: decide if seed is unique for all MPIs
       !
       opt_broadcast_seed = optional_argument(.false., broadcast_seed)

       if (opt_random_from_file) then
          !
          ! Fixed seed: random numbers are read from file
          !
          if( urand == 0 ) then
             irand = 0
             urand = iofile_available_unit()
             call iofile_open_unit(urand,'./random.txt','RANDOM NUMBER LIST','old')
             read(UNIT=urand,IOSTAT=i,IOMSG=my_iomsg,FMT='(i10)') nrand
             if( i /= 0 ) call runend('MOD_RANDOM: '//trim(my_iomsg)) !FIRST NUMBER IN RANDOM FILE SHOULD BE THE NUMBER OF RANDOM NUMBERS TO READ')
          end if

       else
          !
          ! Generate seed vicious
          !
          call random_seed(size=k)
          allocate (seed(k))
          seed(:) = 1_4
          !
          ! Factor depending on the MPI rank if required
          !
          r = 1.0_rp
          if (.not. opt_broadcast_seed) then
             if (PAR_CODE_SIZE /= 0_ip) then
                r = real(max(PAR_MY_CODE_RANK+1_ip,1_ip),rp)/real(PAR_CODE_SIZE,rp)
             end if
          end if
          !
          ! Seed
          !
          call system_clock(COUNT=clock)
          seed = int(clock,4) + 37_4 * (/ (j-1, j = 1,k) /)
          if( any(seed==0_4) ) then
             call date_and_time(VALUES=values)
             seed = int(sum(values(1:8)*1000),4) + 37_4 * (/ (j-1, j = 1,k) /)             
          end if          
          call random_seed(put=seed)
          call random_number(s)
          seed = int(real(seed,rp)*s*r,4)
          !
          ! Master broadcast the seed if required
          !
          if (opt_broadcast_seed) then
             if (present(PAR_COMM_IN)) then
                call PAR_BROADCAST(seed, PAR_COMM_IN=PAR_COMM_IN)
             else
                call PAR_BROADCAST(seed)
             end if
          end if
          !
          ! Put seed
          !
          call random_seed(put=seed)
          deallocate (seed)
       end if

    end if

  end subroutine random_initialization

  !-----------------------------------------------------------------------
  !>
  !> @author  calmet, houzeaux
  !> @date    2022-10-26
  !> @brief   Random number generator
  !>          This subroutine must be called in the same way the
  !>          initialization was carried out. This is to check possible
  !>          incoherencies between both.
  !
  !-----------------------------------------------------------------------

  function random_generate_number(broadcast_seed) result(r)

    logical(lg), intent(in), optional :: broadcast_seed 
    real(rp)                          :: r
    integer(ip)                       :: ierr
    
    if ( opt_random_from_file ) then
       
       irand = irand + 1
       read(urand,*,iostat=ierr) r
       if( mod(irand,nrand) == 0 .or. ierr /= 0 ) then
          rewind(urand)
          read(urand,*)
       end if
       
    else if( opt_broadcast_seed .neqv. optional_argument(.false., broadcast_seed) ) then
       
       if( opt_broadcast_seed ) then
          call runend("MOD_RANDOM: SEED WAS INITIALIZED AS UNIQUE AMONG MPIS AND DOES NOT CORRESPOND TO THE GENERATION")
       else
          call runend("MOD_RANDOM: SEED WAS INITIALIZED AS DIFFERENT AMONG MPIS AND DOES NOT CORRESPOND TO THE GENERATION")
       end if
       
    else
       
       call random_number(r)
       !block
       !  use def_master
       !  write(90+kfl_paral,*) r
       !  flush(90+kfl_paral)
       !end block 

    end if

  end function random_generate_number
  
  !-----------------------------------------------------------------------
  !>
  !> @author  houzeaux
  !> @date    2022-10-26
  !> @brief   End a random number generator campaign
  !
  !-----------------------------------------------------------------------

  subroutine random_end()

  end subroutine random_end
  
end module mod_random
!> @}
