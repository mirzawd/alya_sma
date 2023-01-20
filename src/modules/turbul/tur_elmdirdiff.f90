!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tur_elmdirdiff(pnode,lnods,elresdiff)
  !-----------------------------------------------------------------------
  !****f* Turbul/tur_elmdirdiff
  ! NAME
  !   tur_elmdirdiff
  ! DESCRIPTION
  ! This routine prescribes the boundary conditions for the 
  ! turbulence equations to the elresdiff.
  ! USES
  ! USED BY
  !    tur_elmop1
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp
  use def_turbul, only       :  kfl_fixno_tur,&
       &                        iunkn_tur,kfl_algor_tur
  implicit none
  integer(ip), intent(in)    :: pnode
  integer(ip), intent(in)    :: lnods(pnode)
  real(rp),    intent(inout) :: elresdiff(pnode)
  integer(ip)                :: inode,ipoin

  if( kfl_algor_tur == 1 ) then

        !
        ! Fluid element
        !
        do inode = 1,pnode
           ipoin = lnods(inode)
           if(  kfl_fixno_tur(1,ipoin,iunkn_tur) > 0 ) then
              elresdiff(inode) = 0.0_rp
           end if
        end do


  end if

end subroutine tur_elmdirdiff
