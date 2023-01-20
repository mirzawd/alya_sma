!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!
!> @defgroup Permutation_Toolbox
!> Tollbox for permutation of arrays
!> @{
!> @name    ToolBox for permutation
!> @file    mod_permut.f90
!> @author  Guillaume Houzeaux, Damien Dosimont
!> @brief   ToolBox for output
!> @details ToolBox for output, mainly for debugging
!> @{
!
!-----------------------------------------------------------------------

module mod_permut

  use def_kintyp, only : ip,rp,lg
  implicit none

  private

  interface permut
     module procedure &
          permut_rp1,&
          permut_rp2,&
          permut_ip1,&
          permut_ip2
  end interface

  public :: permut, permut_rp1, permut_rp2, permut_ip1, permut_ip2

contains

  subroutine permut_rp1(pleng,permu,bridge_in,bridge_out)

    integer(ip), intent(in)  :: pleng
    integer(ip), intent(in)  :: permu(*)
    real(rp),    intent(in)  :: bridge_in(*)
    real(rp),    intent(out) :: bridge_out(*)
    integer(ip)              :: ii,kk

    do ii = 1,pleng
       kk = permu(ii)
       bridge_out(ii) = bridge_in(kk)
    end do
  end subroutine

  subroutine permut_rp2(pdime,pleng,permu,bridge_in,bridge_out)

    integer(ip), intent(in)  :: pdime
    integer(ip), intent(in)  :: pleng
    integer(ip), intent(in)  :: permu(*)
    real(rp),    intent(in)  :: bridge_in(pdime,*)
    real(rp),    intent(out) :: bridge_out(pdime,*)
    integer(ip)              :: ii,kk

    do ii = 1,pleng
       kk = permu(ii)
       bridge_out(1:pdime,ii) = bridge_in(1:pdime,kk)
    end do
  end subroutine

  subroutine permut_ip1(pleng,permu,bridge_in,bridge_out)

    integer(ip), intent(in)  :: pleng
    integer(ip), intent(in)  :: permu(*)
    integer(ip), intent(in)  :: bridge_in(*)
    integer(ip), intent(out) :: bridge_out(*)
    integer(ip)              :: ii,kk

    do ii = 1,pleng
       kk = permu(ii)
       bridge_out(ii) = bridge_in(kk)
    end do
  end subroutine

  subroutine permut_ip2(pdime,pleng,permu,bridge_in,bridge_out)

    integer(ip), intent(in)  :: pdime
    integer(ip), intent(in)  :: pleng
    integer(ip), intent(in)  :: permu(*)
    integer(ip), intent(in)  :: bridge_in(pdime,*)
    integer(ip), intent(out) :: bridge_out(pdime,*)
    integer(ip)              :: ii,kk

    do ii = 1,pleng
       kk = permu(ii)
       bridge_out(1:pdime,ii) = bridge_in(1:pdime,kk)
    end do
  end subroutine

end module mod_permut
!> @}
