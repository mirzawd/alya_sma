!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine elmtyp()
  !-----------------------------------------------------------------------
  !****f* Domain/elmtyp
  ! NAME
  !    cderda
  ! DESCRIPTION
  !    This routine defines the different types of elements:
  !
  !    NNODE ... number fo nodes
  !
  !                 1D           2D           3D
  !                 ----------------------------------
  !    LTOPO ... -1 Lines          -
  !               0    -    Quadrilateral  Hexahedra
  !               1    -    Triangle       Tetrahedra
  !               2    -           -       Pentahedra
  !               3    -           -       Pyramid
  !
  !    LLAPL ... 1 if Hessian matrix should be considered
  !              0 otherwise
  !
  ! USED BY
  !    Domain
  !***
  !-----------------------------------------------------------------------

  use def_kintyp_basic, only : ip
  use def_domain,       only : nelty,lenex
  use def_domain,       only : ldime,lorde
  use def_domain,       only : ltopo,nnode
  use def_domain,       only : llapl
  use def_domain,       only : mface
  use def_elmgeo,       only : element_type

  implicit none
  integer(ip) :: pelty,inode,pnode

  !---------------------------------------------------------------------
  !
  ! Copy element databse to global variable
  !
  !---------------------------------------------------------------------
  
  do pelty = 1,size(element_type)
     ldime(pelty) = element_type(pelty) % dimensions
     lorde(pelty) = element_type(pelty) % order
     ltopo(pelty) = element_type(pelty) % topology
     nnode(pelty) = element_type(pelty) % number_nodes
     if( element_type(pelty) % hessian ) then
        llapl(pelty) = 1
     else
        llapl(pelty) = 0
     end if
  end do
  
  !---------------------------------------------------------------------
  !
  ! Derived parameters
  !
  !---------------------------------------------------------------------
  !
  ! Next element node
  !
  do pelty = 1,nelty
     pnode = nnode(pelty)
     if( pnode > 0 ) then
        do inode = 1,pnode-1
           lenex(inode,pelty) = inode + 1
        end do
        lenex(pnode,pelty) = 1
     end if
  end do
  !
  ! Mirror NNODE
  !
  do pelty = 1,nelty
     nnode(-pelty) = nnode(pelty)
  end do
  mface = maxval(element_type(:) % number_faces)
  
end subroutine elmtyp
