!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Parall
!> @{
!> @file    mod_parall_destructor.f90
!> @author  houzeaux
!> @date    2020-03-06
!> @brief   Communicator destruction
!> @details Communicator destruction
!-----------------------------------------------------------------------

module mod_parall_destructor

  use def_kintyp,                  only : ip,lg
  use def_master,                  only : current_code
  use mod_parall,                  only : PAR_COMM_MY_CODE_ARRAY
  use mod_parall,                  only : PAR_COMM_COLOR 
  use mod_parall,                  only : mcolo
  use mod_parall,                  only : par_code_zone_subd_to_color
  use mod_communications,          only : PAR_COMM_FREE
  use mod_communications,          only : PAR_COMM_NULL
  use mod_par_color_communicators, only : par_color_communicators_deallocate

  implicit none

  private

  public :: parall_destructor

contains
  !-----------------------------------------------------------------------
  !>
  !> @author  houzeaux
  !> @date    2019-06-11
  !> @brief   Destroy parallelization
  !> @details Destroy parallelization
  !>
  !-----------------------------------------------------------------------

  subroutine parall_destructor(COMM_MY_CODE_ARRAY)

    logical(lg), optional,   intent(in) :: COMM_MY_CODE_ARRAY
    integer(ip)                         :: icolo,jcolo
    integer(ip)                         :: color_current_code
    integer(ip)                         :: color_world
    logical(lg)                         :: if_comm_my_code_array

    if( present(COMM_MY_CODE_ARRAY) ) then
       if_comm_my_code_array = COMM_MY_CODE_ARRAY
    else
       if_comm_my_code_array = .true.
    end if
    !
    ! Global communicator
    !
    if( if_comm_my_code_array ) &
         call PAR_COMM_MY_CODE_ARRAY(1) % deallo(COMM_NAME='COMMD',INITIALIZE=.true.)
    !
    ! Communicators: free only half has it is symmetric...
    !
    color_current_code = par_code_zone_subd_to_color(current_code,0_ip,0_ip)
    color_world        = par_code_zone_subd_to_color(0_ip,0_ip,0_ip)

    do icolo = 0,mcolo
       do jcolo = icolo,mcolo
          if( ( icolo == color_world .and. jcolo == color_world ) .or. ( icolo == color_current_code .and. jcolo == color_current_code ) ) then
             continue
          else
             call PAR_COMM_FREE(PAR_COMM_COLOR(icolo,jcolo)) ! Free communicator
             PAR_COMM_COLOR(icolo,jcolo) = PAR_COMM_NULL     ! Put MPU_COMM_NLL                    
             PAR_COMM_COLOR(jcolo,icolo) = PAR_COMM_NULL     ! Symmetric communicator already freed
          end if
       end do
    end do
    !
    ! Deallocate colo communicators
    !
    call par_color_communicators_deallocate()

  end subroutine parall_destructor

end module mod_parall_destructor
!> @}
