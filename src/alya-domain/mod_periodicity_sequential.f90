!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!
!> @addtogroup Domain
!> @{
!> @name    ToolBox for periodicity
!> @file    mod_periodicity.f90
!> @author  Guillaume Houzeaux
!> @brief   ToolBox periodicity
!> @details ToolBox periodicity
!> @{
!
!-----------------------------------------------------------------------

module mod_periodicity_sequential

  use def_kintyp_basic,          only : ip,rp,lg,i1p
  use def_master,                only : ISEQUEN
  implicit none

  public :: periodicity_sequential
  
contains


  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-12-14
  !> @brief   Impose periodicity
  !> @details Impose periodicity
  !> 
  !-----------------------------------------------------------------------

  subroutine periodicity_sequential(ndofn,xx)

    use def_domain
    integer(ip), intent(in)    :: ndofn
    real(rp),    intent(inout) :: xx(ndofn,*)
    integer(ip)                :: ipoin,jpoin,idime

    if( ISEQUEN ) then
       !
       ! Master = master + slaves
       !
       do ipoin = 1,npoin
          jpoin = lmast(ipoin)
          if(jpoin>0) then
             do idime = 1,ndofn
                xx(idime,jpoin) = xx(idime,jpoin) + xx(idime,ipoin) 
             end do
          end if
       end do
       !
       ! Slaves = master
       !
       do ipoin = 1,npoin
          jpoin = lmast(ipoin)
          if(jpoin>0) then
             do idime = 1,ndofn
                xx(idime,ipoin) = xx(idime,jpoin) 
             end do
          end if
       end do

    end if

  end subroutine periodicity_sequential

end module mod_periodicity_sequential
!> @}
