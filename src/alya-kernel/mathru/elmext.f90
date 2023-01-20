!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine elmext(&
     itask,ndofn,pnode,elauu,elaup,elapu,elapp,elmap,&
     elrbu,elrbp,elunk)
  !-----------------------------------------------------------------------
  !****f* Nastin/elmext
  ! NAME 
  !    elmext
  ! DESCRIPTION
  !    Modify element matrix when extension elements are used
  !    Only equation of the first node should be assembled as
  !    it corresponds to its extension test function
  !
  ! USES
  ! USED BY
  !    elmext
  !***
  !----------------------------------------------------------------------- 
  use def_kintyp, only       :  ip,rp
  use def_domain, only       :  ndime
  implicit none
  integer(ip), intent(in)    :: itask,ndofn,pnode
  real(rp),    intent(out)   :: elauu(pnode*ndime,pnode*ndime)
  real(rp),    intent(out)   :: elaup(pnode*ndime,pnode)
  real(rp),    intent(out)   :: elapu(pnode,pnode*ndime)
  real(rp),    intent(out)   :: elapp(pnode,pnode)
  real(rp),    intent(out)   :: elmap(pnode*ndofn,pnode*ndofn)
  real(rp),    intent(out)   :: elrbu(pnode*ndime)
  real(rp),    intent(out)   :: elrbp(pnode*ndofn)
  real(rp),    intent(in)    :: elunk(pnode*ndofn)
  integer(ip)                :: inode,jnode,idofn,jdofn

  if( itask == 1 ) then
     !
     ! Auu, Aup, Apu, App, bu, bp
     !
     do idofn = ndime+1,ndime*pnode
        do jdofn = 1,ndime*pnode
           elauu(idofn,jdofn) = 0.0_rp
        end do
        do jnode = 1,pnode
           elaup(idofn,jnode) = 0.0_rp
        end do
        elrbu(idofn) = 0.0_rp
     end do
     do inode = 2,pnode
        do jdofn = 1,ndime*pnode
           elapu(inode,jdofn) = 0.0_rp
        end do
        do jnode = 1,pnode
           elapp(inode,jnode) = 0.0_rp
        end do
        elrbp(inode) = 0.0_rp
     end do

  else if( itask == 2 ) then
     !
     ! Q
     !
     return
     do inode = 2,pnode
        do jnode = 1,pnode
           elmap(inode,jnode) = 0.0_rp
        end do
     end do

  else if( itask == 3 ) then
     !
     ! b
     !
     do inode = ndofn+1,pnode*ndofn
        elrbp(inode) = 0.0_rp
     end do

  else if( itask == 4 ) then
     !
     ! A and b
     !
     do inode = ndofn+1,ndofn*pnode
        do jnode = 1,pnode*ndofn
           elmap(inode,jnode) = 0.0_rp
        end do
     end do
     do inode = ndofn+1,pnode*ndofn
        elrbp(inode) = 0.0_rp
     end do

  else if( itask == 5 ) then
     !
     ! b-Ax
     !
     do inode = ndofn+1,ndofn*pnode
        do jnode = 1,pnode*ndofn
           elmap(inode,jnode) = 0.0_rp
        end do
     end do
     do inode = ndofn+1,pnode*ndofn
        elrbp(inode) = 0.0_rp
     end do
     do inode = ndofn+1,ndofn*pnode
        elrbp(1) = elrbp(1) - elmap(1,inode) * elunk(inode)
        elmap(1,inode) = 0.0_rp
     end do

  end if

end subroutine elmext
