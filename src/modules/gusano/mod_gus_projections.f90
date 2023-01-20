!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Gusano
!> @{
!> @file    mod_gus_projections.f90
!> @author  houzeaux
!> @date    2020-10-22
!> @brief   Projections
!> @details Projections for stabilization
!-----------------------------------------------------------------------

module mod_gus_projections

  use def_kintyp_basic,                  only : ip
  use def_master,                        only : ITASK_ENDINN
  use def_domain,                        only : npoin
  use def_domain,                        only : vmass
  use mod_communications_point_to_point, only : PAR_INTERFACE_NODE_EXCHANGE
  use def_gusano
  
  implicit none
  private

  public :: gus_projections_initialization
  public :: gus_projections_solve
  public :: gus_projections_updunk

contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-10-22
  !> @brief   Initialize projection
  !> @details Initialize projections
  !> 
  !-----------------------------------------------------------------------

  subroutine gus_projections_initialization()

    integer(ip) :: ipoin
    
    if( associated(projm_gus) ) then
       do ipoin = 1,npoin
          projm_gus(ipoin,1) = 0.0_rp
       end do
    end if
    if( associated(projc_gus) ) then
       do ipoin = 1,npoin
          projc_gus(ipoin,1) = 0.0_rp
       end do
    end if

  end subroutine gus_projections_initialization
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-10-22
  !> @brief   Solve projection
  !> @details Solve projections
  !> 
  !-----------------------------------------------------------------------

  subroutine gus_projections_solve()

    integer(ip)         :: ipoin
    real(rp),   pointer :: projm1(:)
    
    if( associated(projm_gus) ) then
       projm1 => projm_gus(:,1)
       call PAR_INTERFACE_NODE_EXCHANGE(projm1,'SUM')
       do ipoin = 1,npoin
          projm_gus(ipoin,1) = projm_gus(ipoin,1) / vmass(ipoin)
       end do
    end if
    if( associated(projc_gus) ) then
       projm1 => projc_gus(:,1)
       call PAR_INTERFACE_NODE_EXCHANGE(projm1,'SUM')
       do ipoin = 1,npoin
          projc_gus(ipoin,1) = projc_gus(ipoin,1) / vmass(ipoin)
       end do
    end if
    
  end subroutine gus_projections_solve

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-10-22
  !> @brief   Update projection
  !> @details Update projections
  !> 
  !-----------------------------------------------------------------------

  subroutine gus_projections_updunk(itask)

    integer(ip), intent(in) :: itask
    integer(ip)             :: ipoin

    select case ( itask )
       
    case ( ITASK_ENDINN )
       
       do ipoin = 1,npoin
          projm_gus(ipoin,2) = projm_gus(ipoin,1)
          projc_gus(ipoin,2) = projc_gus(ipoin,1)
       end do
       
    end select
    
  end subroutine gus_projections_updunk
  
end module mod_gus_projections
!> @}

