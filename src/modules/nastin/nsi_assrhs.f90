!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_assrhs(&
     itask,ndofn,pnode,pevat,ndime,kfl_algso,lnods,&
     elvel,elrhs,elmat,rhsid)
  !-----------------------------------------------------------------------
  !****f* nastin/nsi_assrhs
  ! NAME 
  !    nsi_assrhs
  ! DESCRIPTION
  !    Assembly of the RHS
  ! USES
  ! USED BY
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only        :  ip,rp
  use def_nastin, only        :  NSI_SCHUR_COMPLEMENT,&
       &                         ndbgs_nsi,NSI_FRACTIONAL_STEP
  use def_kermod, only        :  kfl_adj_prob
  implicit none
  integer(ip),  intent(in)    :: itask,ndofn,pnode,pevat,ndime
  integer(ip),  intent(in)    :: kfl_algso
  integer(ip),  intent(in)    :: lnods(pnode)
  real(rp),     intent(in)    :: elrhs(pevat),elmat(pevat,pevat)
  real(rp),     intent(in)    :: elvel(ndime,pnode)
  real(rp),     intent(inout) :: rhsid(*)
  integer(ip)                 :: inode,ipoin,idofg,idofl,idofn,idime
  integer(ip)                 :: jnode,jdofl,jdime

  if (kfl_adj_prob == 0) then 
    
  if( kfl_algso == -3 .and. itask /= 3 ) then
     !
     ! Matrix-free Richardson
     !
     do inode=1,pnode
        ipoin = lnods(inode)
        idofg = (ipoin-1)*ndofn
        idofl = (inode-1)*ndofn
        do idime=1,ndime
           idofg=idofg+1
           idofl=idofl+1
           do jnode=1,pnode
              jdofl = (jnode-1)*ndofn
              do jdime=1,ndime
                 jdofl=jdofl+1
                 rhsid(idofg)=rhsid(idofg)&
                      -elmat(idofl,jdofl)*elvel(jdime,jnode)
              end do
           end do
           rhsid(idofg)=rhsid(idofg)+elrhs(idofl)
        end do
     end do
     
  else if( NSI_SCHUR_COMPLEMENT .or. NSI_FRACTIONAL_STEP ) then
     !
     ! Schur complement system
     ! RHS(U): rhsid(      1,1:npoin)
     ! RHS(V): rhsid(      2,1:npoin)
     ! RHS(W): rhsid(      3,1:npoin)
     ! RHS(P): rhsid(      1,ndbgs_nsi+1:ndbgs_nsi+npoin)
     !
     do inode = 1,pnode
        ipoin = lnods(inode)
        idofg = (ipoin-1) * ndime
        idofl = (inode-1) * ndofn
        do idime = 1,ndime
           idofg = idofg+1
           idofl = idofl+1
           !$OMP ATOMIC
           rhsid(idofg) = rhsid(idofg) + elrhs(idofl)
        end do
        idofg = ndbgs_nsi + ipoin
        idofl = idofl + 1
        !$OMP ATOMIC
        rhsid(idofg) = rhsid(idofg) + elrhs(idofl)
     end do
    
  else
     !
     ! General case:
     ! RHS(U): rhsid(      1,1:npoin)
     ! RHS(V): rhsid(      2,1:npoin)
     ! RHS(W): rhsid(      3,1:npoin)
     ! RHS(P): rhsid(ndime+1,1:npoin)
     !
     if(ndofn==1) then
        do inode=1,pnode
           ipoin=lnods(inode)
           !$OMP ATOMIC
           rhsid(ipoin)=rhsid(ipoin)+elrhs(inode)
        end do
     else if(ndofn==2) then
        do inode=1,pnode
           ipoin=lnods(inode)
           idofg=(ipoin-1)*ndofn
           idofl=(inode-1)*ndofn
           idofg=idofg+1
           idofl=idofl+1
           !$OMP ATOMIC
           rhsid(idofg)=rhsid(idofg)+elrhs(idofl)
           idofg=idofg+1
           idofl=idofl+1
           !$OMP ATOMIC
           rhsid(idofg)=rhsid(idofg)+elrhs(idofl)
        end do
     else if(ndofn==3) then
        do inode=1,pnode
           ipoin=lnods(inode)
           idofg=(ipoin-1)*ndofn
           idofl=(inode-1)*ndofn
           idofg=idofg+1
           idofl=idofl+1
           !$OMP ATOMIC
           rhsid(idofg)=rhsid(idofg)+elrhs(idofl)
           idofg=idofg+1
           idofl=idofl+1
           !$OMP ATOMIC
           rhsid(idofg)=rhsid(idofg)+elrhs(idofl)
           idofg=idofg+1
           idofl=idofl+1
           !$OMP ATOMIC
           rhsid(idofg)=rhsid(idofg)+elrhs(idofl)
        end do
     else
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
     end if

  end if
  
  elseif (kfl_adj_prob == 1) then
  
  if( ( NSI_SCHUR_COMPLEMENT .or. NSI_FRACTIONAL_STEP ) .and. itask == 1) then
     !
     ! Schur complement system
     ! RHS(U): rhsid(      1,1:npoin)
     ! RHS(V): rhsid(      2,1:npoin)
     ! RHS(W): rhsid(      3,1:npoin)
     ! RHS(P): rhsid(      1,ndbgs_nsi+1:ndbgs_nsi+npoin)
     !
     do inode = 1,pnode
        ipoin = lnods(inode)
        idofg = (ipoin-1) * ndime
        idofl = (inode-1) * ndofn
        do idime = 1,ndime
           idofg = idofg+1
           idofl = idofl+1
           !$OMP ATOMIC
           rhsid(idofg) = rhsid(idofg) - elrhs(idofl) ! rhsid  = rhsid - dF/dU
        end do
        idofg = ndbgs_nsi + ipoin
        idofl = idofl + 1
        !$OMP ATOMIC
        rhsid(idofg) = rhsid(idofg) - elrhs(idofl) ! rhsid  = rhsid - dF/dP
     end do
    
  endif

  if( ( NSI_SCHUR_COMPLEMENT .or. NSI_FRACTIONAL_STEP ) .and. itask == 2) then
     !
     ! Schur complement system
     ! RHS(U): rhsid(      1,1:npoin)
     ! RHS(V): rhsid(      2,1:npoin)
     ! RHS(W): rhsid(      3,1:npoin)
     ! RHS(P): rhsid(      1,ndbgs_nsi+1:ndbgs_nsi+npoin)
     !
     do inode = 1,pnode
        ipoin = lnods(inode)
        idofg = (ipoin-1) * ndime
        idofl = (inode-1) * ndime
        do idime = 1,ndime
           idofg = idofg+1
           idofl = idofl+1
           !$OMP ATOMIC
           rhsid(idofg) = rhsid(idofg) + elrhs(idofl) ! dcost_dx_nsi  = dcost_dx_nsi + dF/dX
        end do
     end do

  endif
  
  
  endif

end subroutine nsi_assrhs
