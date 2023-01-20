!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Domain
!> @{
!> @file    mod_groups.f90
!> @author  houzeaux
!> @date    2020-03-19
!> @brief   Groups
!> @details Tools to manage groups
!-----------------------------------------------------------------------

module mod_groups

  use def_master
  use def_domain
  use mod_htable
  use mod_messages,       only : messages_live
  use mod_communications, only : PAR_MIN
  use mod_communications, only : PAR_SUM
  implicit none
  private

  public :: groups_check_prescription
  
contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-03-19
  !> @brief   Check prescription
  !> @details Check if at least on dof has been prescribed
  !> 
  !-----------------------------------------------------------------------

  subroutine groups_check_prescription(limpo)
    
    integer(ip), pointer, intent(inout) :: limpo(:)
    integer(ip)                         :: num_nodes
    integer(ip)                         :: ipoin
    integer(ip)                         :: ipoin_min
    integer(ip)                         :: ipoin_fix
    !
    ! NUM_NODES: number of prescribed nodes
    !
   num_nodes = 0
    if( associated(limpo) ) then
       do ipoin = 1,npoin
          if( solve_sol(1) % limpo(ipoin) > 0 ) then
             num_nodes = num_nodes + 1
          end if
       end do
    end if
    call PAR_SUM(num_nodes)
    !
    ! No node is prescribed: prescibe at least one, possibly on boundary
    !
    ipoin_min = huge(1_ip)
    if( num_nodes == 0 ) then
       loop_ipoin: do ipoin = 1,npoin
          if( lpoty(ipoin) /= 0 ) then !.and. lmast(ipoin) == 0 ) then
             ipoin_min = lninv_loc(ipoin)
             exit loop_ipoin
          end if
       end do loop_ipoin
       call PAR_MIN(ipoin_min)
       ipoin_fix = 0
       loop_ipoin2: do ipoin = 1,npoin
          if( lninv_loc(ipoin) == ipoin_min ) then
             ipoin_fix = ipoin
             exit loop_ipoin2
          end if
       end do loop_ipoin2
       if( ipoin_fix /= 0 .and. associated(limpo) ) then
          limpo(ipoin_fix) = 1
       end if
       call messages_live('GROUPS: ONE NODE HAS BEEN AUTOMATICALLY PRESCRIBED: '//intost(ipoin_min))
    end if
    
  end subroutine groups_check_prescription

end module mod_groups
!> @}
