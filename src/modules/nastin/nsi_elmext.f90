!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_elmext(&
     itask,pnode,elauu,elaup,elapu,elapp,elmap,elrbu,elrbp)
  !-----------------------------------------------------------------------
  !****f* Nastin/nsi_elmext
  ! NAME 
  !    nsi_elmext
  ! DESCRIPTION
  !    Modify element matrix when extension elements are used
  !    Only equation of the first node should be assembled as
  !    it corresponds to its extension test function
  !
  ! USES
  ! USED BY
  !    nsi_elmext
  !***
  !----------------------------------------------------------------------- 
  use def_kintyp, only       :  ip,rp
  use def_domain, only       :  ndime
  implicit none
  integer(ip), intent(in)    :: itask,pnode
  real(rp),    intent(out)   :: elauu(pnode*ndime,pnode*ndime)
  real(rp),    intent(out)   :: elaup(pnode*ndime,pnode)
  real(rp),    intent(out)   :: elapu(pnode,pnode*ndime)
  real(rp),    intent(out)   :: elapp(pnode,pnode)
  real(rp),    intent(out)   :: elmap(pnode,pnode)
  real(rp),    intent(out)   :: elrbu(pnode*ndime)
  real(rp),    intent(out)   :: elrbp(pnode)
  integer(ip)                :: inode,jnode,idofn,jdofn,iwhat

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
     iwhat = 1

     if( iwhat == 1 ) then
        !
        ! Full matrix
        !
        return

     else if( iwhat == 2 ) then
        !
        ! All zero
        !
        do inode = 1,pnode
           do jnode = 1,pnode
              elmap(inode,jnode) = 0.0_rp
           end do
        end do

     else  if( iwhat == 3 ) then
        !
        ! Diagonal
        !
        do inode = 1,pnode
           do jnode = 1,pnode
              if( inode /= jnode ) elmap(inode,jnode) = 0.0_rp
           end do
        end do

     else if( iwhat == 4 ) then
        !
        ! Symmetrized
        !
        do inode = 2,pnode
           do jnode = 2,pnode
              elmap(inode,jnode) = 0.0_rp
           end do
        end do
        do inode = 1,pnode
           do jnode = 1,pnode
              if( abs(elmap(inode,jnode)-elmap(jnode,inode)) > 1.0e-8_rp ) print*,'merde'
           end do
        end do


     end if
     
  end if

end subroutine nsi_elmext
