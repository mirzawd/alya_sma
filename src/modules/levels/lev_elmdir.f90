!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine lev_elmdir(pnode,lnods,elmat,elrhs)
  !------------------------------------------------------------------------
  !****f* Levels/lev_elmdir
  ! NAME 
  !    lev_elmdir
  ! DESCRIPTION
  !
  ! Impose Dirichlet conditions
  ! 
  ! USED BY
  !    lev_matrix
  !------------------------------------------------------------------------

  use def_kintyp, only       :  ip,rp
  use def_levels, only       :  bvess_lev,kfl_fixno_lev
  implicit none
  integer(ip), intent(in)    :: pnode
  integer(ip), intent(in)    :: lnods(pnode)
  real(rp),    intent(inout) :: elmat(pnode,pnode),elrhs(pnode)
  real(rp)                   :: adiag
  integer(ip)                :: inode,ipoin,jnode

  ! Dirichlet condition
  do inode = 1,pnode
     ipoin = lnods(inode)    
     if(    kfl_fixno_lev(1,ipoin) == 1 .or. & 
          & kfl_fixno_lev(1,ipoin) == 8 ) then
        adiag = elmat(inode,inode)
        do jnode = 1,pnode
           elmat(inode,jnode) = 0.0_rp
           elrhs(jnode)       = elrhs(jnode) - elmat(jnode,inode) * bvess_lev(1,ipoin,1)
           elmat(jnode,inode) = 0.0_rp
        end do
        elmat(inode,inode) = adiag
        elrhs(inode) = adiag * bvess_lev(1,ipoin,1)
     end if
  end do

end subroutine lev_elmdir
