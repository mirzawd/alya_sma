!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_elmibm(&
     pnode,pevat,gpeps,gbsur,shaib,elmat,elrhs,gpcib,gptau,gbcod)
  !------------------------------------------------------------------------
  !****f* Temper/nsi_elmibm
  ! NAME 
  !    nsi_elmibm
  ! DESCRIPTION
  ! USES
  ! USED BY
  !***
  !------------------------------------------------------------------------
  use def_kintyp, only     :  ip,rp
  use def_domain, only     :  ndime
  implicit none
  integer(ip), intent(in)  :: pnode,pevat
  real(rp),    intent(in)  :: gpeps,gbsur,shaib(pnode)
  real(rp),    intent(out) :: elmat(pevat,pevat)
  real(rp),    intent(out) :: elrhs(pevat)
  integer(ip)              :: inode,jnode,ievat,jevat,idime
  real(rp)                 :: gpcib(ndime,pnode),gptau,gbcod(ndime)
  !
  ! Initialize
  !
  do ievat=1,pevat
     elrhs(ievat)=0.0_rp
     do jevat=1,pevat
        elmat(jevat,ievat)=0.0_rp
     end do
  end do
  !
  ! Fill in
  !
  do inode=1,pnode
     ievat=(inode-1) * ndime
     do idime=1,ndime
        ievat=ievat+1
        elrhs(ievat)=elrhs(ievat)+gpeps*gbsur*shaib(inode)*0.0_rp
        do jnode=1,pnode
           jevat=(jnode-1)*ndime+idime
           elmat(ievat,jevat)=elmat(ievat,jevat)&
                +gpeps*gbsur*shaib(inode)*shaib(jnode)
        end do
     end do
  end do

  return

  do inode=1,pnode
     ievat=inode*(ndime+1)
     do jnode=1,pnode
        do idime=1,ndime
           jevat=(jnode-1)*(ndime+1)+idime
           elmat(ievat,jevat)=elmat(ievat,jevat)&
                +gpeps*gbsur*gpcib(idime,inode)*shaib(jnode)
        end do
     end do
  end do

  return
  !
  !
  !

  !if( sqrt(gbcod(1)**2+gbcod(2)**2)>=0.5_rp) then
     do inode=1,pnode
        ievat=inode*(ndime+1)
        do jnode=1,pnode
           do idime=1,ndime
              jevat=(jnode-1)*(ndime+1)+idime
              elmat(ievat,jevat)=elmat(ievat,jevat)&
                   +gptau*gpeps*gbsur*gpcib(idime,inode)*shaib(jnode)
           end do
        end do
     end do
  !end if

end subroutine nsi_elmibm
