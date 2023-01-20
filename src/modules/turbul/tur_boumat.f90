!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tur_boumat(&
     pnode,pnodb,lboel,xmmat,xmrhs,gbsha,&
     gbsur,elmat,elrhs)
  !------------------------------------------------------------------------
  !****f* turbul/tur_boumat
  ! NAME 
  !    tur_boumat
  ! DESCRIPTION
  !    Assemble boundary contribution
  ! USES
  ! USED BY
  !    tur_bouope
  !***
  !------------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp
 
  implicit none
  integer(ip), intent(in)    :: pnode,pnodb
  integer(ip), intent(in)    :: lboel(pnodb)
  real(rp),    intent(in)    :: xmmat,xmrhs,gbsha(pnodb),gbsur
  real(rp),    intent(inout) :: elmat(pnode,pnode),elrhs(pnode)
  
  ! local variables
  integer(ip)                :: inodb,jnodb,inode,jnode
  real(rp)                   :: xmuit

  do inodb=1,pnodb  
     inode=lboel(inodb)
     elrhs(inode)=elrhs(inode)+gbsha(inodb)*xmrhs*gbsur
     xmuit=xmmat*gbsha(inodb)*gbsur
     do jnodb=1,pnodb
        jnode=lboel(jnodb)
        elmat(inode,jnode)=elmat(inode,jnode)&
             +xmuit*gbsha(jnodb)
     end do
  end do
 
end subroutine tur_boumat
