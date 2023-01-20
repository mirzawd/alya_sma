!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine shafga(xieta,pdime,ptopo,pgaus,pnode,shaga,ierro)
  !-----------------------------------------------------------------------
  !****f* Domain/shafga
  ! NAME
  !    cshder
  ! DESCRIPTION
  !     This routine evaluates shape functions associated to gauss points
  !     for linear and quadratic isoparametric elements   
  ! OUTPUT
  !    SHAGA(mgaus,mnode,nelty)
  ! USED BY
  !    cshder
  ! SOURCE
  !-----------------------------------------------------------------------
  use def_kintyp
  use mod_extrapolation, only : shaga1
  use mod_extrapolation, only : shaga2
  use mod_extrapolation, only : shaga3
  implicit none
  integer(ip), intent(in)    :: pdime,ptopo,pgaus,pnode
  integer(ip), intent(inout) :: ierro
  real(rp),    intent(inout) :: xieta(pdime,pnode)
  real(rp),    intent(out)   :: shaga(pgaus,pnode) 
  integer(ip)                :: inode

  ierro = 0

  do inode = 1,pnode

     select case( pdime )
        
     case( 1_ip )
        call shaga1(xieta(1,inode),pgaus,shaga,ierro)
        
     case( 2_ip )
        call shaga2([xieta(1,inode),xieta(pdime,inode)],ptopo,&
             pgaus,shaga(1,inode),ierro)
        
     case( 3_ip )
        call shaga3([xieta(1,inode),xieta(2,inode),xieta(3,inode)],ptopo,&
             pgaus,shaga(1,inode),ierro)
        
     end select

  end do

end subroutine shafga
!-----------------------------------------------------------------------
! NOTES
! 
! Shaga is used to extrapolate from the Gauss points to the nodes.
!
! X--*------*--X
! |            |  o Element Gauss point:  igaus
! |   o    o   |  * Boundary Gauss point: igaub
! |            |  x Boundary nodes:       inodb
! |   o    o   |
! |            |
! +------------+
!
! A typical use is the following:
! f(igaub)=        \Sum_{inodb=1}^{nnodb} shape_inodb(igaub) f_inodb
! f(inodb)=f_inodb=\Sum_{igaus=1}^{pgaus} shaga_igaus(inode) f_igaus
! Therefore
! f(igaub)=        \Sum_{inodb=1}^{nnodb} \Sum_{igaus=1}^{pgaus} 
!                  (shape_inodb(igaub) shaga_igaus(inode) f_igaus)
! If f are the Cartesian derivatives known on the Gauss points igaus,
! it enables therefore to obtain them on the boundary Gauss points
! igaub.
!
!***
!-----------------------------------------------------------------------
