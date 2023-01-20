!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



 subroutine tem_elmext(pnode,elmat,elrhs)
  !-----------------------------------------------------------------------
  !****f* Nastin/tem_elmext
  ! NAME 
  !    tem_elmext
  ! DESCRIPTION
  !    MOdify element matrix when extension elements are used
  ! USES
  ! USED BY
  !    tem_elmext
  !***
  !----------------------------------------------------------------------- 
  use def_kintyp, only       :  ip,rp
  implicit none
  integer(ip), intent(in)    :: pnode
  real(rp),    intent(out)   :: elmat(pnode,pnode)
  real(rp),    intent(out)   :: elrhs(pnode)
  integer(ip)                :: inode,jnode
  !
  ! A, b
  !
  do inode = 2,pnode
     do jnode = 1,pnode
        elmat(inode,jnode) = 0.0_rp
     end do
     elrhs(inode) = 0.0_rp
  end do
  
end subroutine tem_elmext
