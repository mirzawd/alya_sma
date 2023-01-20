!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine assrhs(ndofn,pnode,lnods,elrhs,rhsid)
  !------------------------------------------------------------------------
  !****f* mathru/assrhs
  ! NAME 
  !    assrhs
  ! DESCRIPTION
  !    Assembly of the RHS
  ! USES
  ! USED BY
  !    *_elmope
  !***
  !------------------------------------------------------------------------
  use def_kintyp, only        :  ip,rp 
  use def_domain, only        :  npoin
  !use omp_lib
  implicit none
  integer(ip),  intent(in)    :: ndofn,pnode
  integer(ip),  intent(in)    :: lnods(pnode)
  real(rp),     intent(in)    :: elrhs(*)
  real(rp),     intent(inout) :: rhsid(*)
  integer(ip)                 :: inode,ipoin,idofl,idofg,idofn,ndof2

  if(ndofn==1) then
     !
     ! 1 DOF
     !
     do inode=1,pnode

        ipoin=lnods(inode)

        !$OMP ATOMIC
        rhsid(ipoin)=rhsid(ipoin)+elrhs(inode)

     end do


  else if(ndofn==2) then
     !
     ! 2 DOF's
     !
     do inode=1,pnode           

        ipoin=lnods(inode)
        idofg=2*ipoin-1
        idofl=2*inode-1

        !$OMP ATOMIC
        rhsid(idofg)=rhsid(idofg)+elrhs(idofl)
        !$OMP ATOMIC
        rhsid(idofg+1)=rhsid(idofg+1)+elrhs(idofl+1)

     end do

  else if(ndofn>2) then
     !
     ! >2 DOF's
     !
     do inode=1,pnode           

        ipoin=lnods(inode)
        idofg=(ipoin-1)*ndofn
        idofl=(inode-1)*ndofn

        do idofn=1,ndofn
           idofg=idofg+1
           idofl=idofl+1
           !$OMP ATOMIC
           rhsid(idofg)=rhsid(idofg)+elrhs(idofl)
        end do

     end do

  else if(ndofn<0) then
     !
     ! >2 DOF's
     !
     ndof2=abs(ndofn)
     do inode=1,pnode           

        ipoin=lnods(inode)
        idofg=(ipoin-1)*ndof2
        idofl=(inode-1)*ndof2

        do idofn=1,ndof2
           idofg=(idofn-1)*npoin+ipoin
           idofl=idofl+1
           !$OMP ATOMIC
           rhsid(idofg)=rhsid(idofg)+elrhs(idofl)
        end do

     end do

  end if

end subroutine assrhs
