!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine lev_elmdi2(pnode,lnods,elmat,elrhs)
  !------------------------------------------------------------------------
  !****f* Levels/lev_elmdi2
  ! NAME 
  !    lev_elmdi2
  ! DESCRIPTION
  !
  ! Impose Dirichlet conditions for redistanciation equation
  ! 
  ! USED BY
  !    lev_elmope
  !------------------------------------------------------------------------

  use def_kintyp, only       :  ip,rp
  use def_master
  use def_levels, only       :  icupt_lev
  implicit none
  integer(ip), intent(in)    :: pnode
  integer(ip), intent(in)    :: lnods(pnode)
  real(rp),    intent(inout) :: elmat(pnode,pnode),elrhs(pnode)
  real(rp)                   :: adiag
  integer(ip)                :: inode,ipoin,jnode

  ! Dirichlet condition
  do inode = 1,pnode
     ipoin = lnods(inode)    
     if( icupt_lev(ipoin) > 0_ip ) then
        adiag = elmat(inode,inode)
        do jnode = 1,pnode
           elmat(inode,jnode) = 0.0_rp
           elrhs(jnode)       = elrhs(jnode) - elmat(jnode,inode) * fleve(ipoin,1)
           elmat(jnode,inode) = 0.0_rp
        end do
        elmat(inode,inode) = adiag
        elrhs(inode) = adiag * fleve(ipoin,1)
     end if
  end do

end subroutine lev_elmdi2
