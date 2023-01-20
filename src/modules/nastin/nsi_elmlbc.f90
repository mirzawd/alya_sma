!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_elmlbc(&
     pnode,pgaus,pevat,gpsha,gpcar,gpvol,elmat,elrhs)
  !----------------------------------------------------------------------
  !****f* Nastin/nsi_elmlbc
  ! NAME 
  !    nsi_elmlap
  ! DESCRIPTION
  !    Compute the Laplacian matrix
  ! USES
  ! USED BY
  !***
  !----------------------------------------------------------------------
  use def_kintyp, only     :  ip,rp 
  use def_domain, only     :  ndime,mnode
  implicit none
  integer(ip), intent(in)  :: pnode,pgaus,pevat
  real(rp),    intent(in)  :: gpsha(pnode,pgaus)
  real(rp),    intent(in)  :: gpcar(ndime,mnode,pgaus)
  real(rp),    intent(in)  :: gpvol(pgaus)
  real(rp),    intent(out) :: elmat(pevat,pevat),elrhs(pevat)
  integer(ip)              :: inode,jnode,kdime,igaus
  integer(ip)              :: idime,idofn,jdofn
  real(rp)                 :: fact1,fact2,fact3
  !
  ! Initialization
  !
  do idofn=1,pevat
     elrhs(idofn)=0.0_rp
     do jdofn=1,pevat
        elmat(jdofn,idofn)=0.0_rp
     end do
  end do
  !
  ! RHS
  !
  do igaus=1,pgaus
     fact1=gpvol(igaus)
     do inode=1,pnode
        idofn=(inode-1)*ndime
        do idime=1,ndime
           idofn=idofn+1
           elrhs(idofn)=elrhs(idofn)+fact1*gpsha(inode,igaus)
        end do
     end do
  end do
  !
  ! Laplacian matrix
  !
  call runend('NSI_ELMLBCS: NOT CODED')
  do igaus=1,pgaus
     do inode=1,pnode
        do jnode=inode+1,pnode
           fact1=0.0_rp
           do kdime=1,ndime
              fact1=fact1+gpcar(kdime,inode,igaus)&
                   &     *gpcar(kdime,jnode,igaus)
           end do
           !!!fact1=fact1*gpvol(igaus)*visco_nsi(1,1)
           do idime=1,ndime
              idofn=(inode-1)*ndime+idime
              jdofn=(jnode-1)*ndime+idime
              elmat(idofn,jdofn)=elmat(idofn,jdofn)+fact1
              elmat(jdofn,idofn)=elmat(jdofn,idofn)+fact1
           end do
        end do
        fact1=0.0_rp
        do kdime=1,ndime
           fact1=fact1+gpcar(kdime,inode,igaus)&
                &     *gpcar(kdime,inode,igaus)
        end do
        !!!fact1=fact1*gpvol(igaus)*visco_nsi(1,1)
        do idime=1,ndime
           idofn=(inode-1)*ndime+idime
           elmat(idofn,idofn)=elmat(idofn,idofn)+fact1
        end do
     end do
  end do
  !
  ! Penalty method 
  !
  !
  ! Velocity: div(u)*div(v)*mu*alpha, alpha>>1
  !
  if(ndime==2) then
     do igaus=1,pgaus
        !!!fact1=visco_nsi(1,1)*10000.0_rp*gpvol(igaus)
        do inode=1,pnode
           idofn=(inode-1)*ndime
           do idime=1,ndime
              idofn=idofn+1
              fact2=fact1*gpcar(idime,inode,igaus)
              do jnode=1,pnode              
                 jdofn=(jnode-1)*ndime   
                 jdofn=jdofn+1
                 fact3=fact2*gpcar(1,jnode,igaus)
                 elmat(idofn,jdofn)=elmat(idofn,jdofn)+fact3
                 jdofn=jdofn+1
                 fact3=fact2*gpcar(2,jnode,igaus)
                 elmat(idofn,jdofn)=elmat(idofn,jdofn)+fact3
              end do
           end do
        end do
     end do
  else
     do igaus=1,pgaus
        !!!fact1=visco_nsi(1,1)*10000.0_rp*gpvol(igaus)
        do inode=1,pnode
           idofn=(inode-1)*ndime
           do idime=1,ndime
              idofn=idofn+1
              fact2=fact1*gpcar(idime,inode,igaus)
              do jnode=1,pnode              
                 jdofn=(jnode-1)*ndime   
                 jdofn=jdofn+1
                 fact3=fact2*gpcar(1,jnode,igaus)
                 elmat(idofn,jdofn)=elmat(idofn,jdofn)+fact3
                 jdofn=jdofn+1
                 fact3=fact2*gpcar(2,jnode,igaus)
                 elmat(idofn,jdofn)=elmat(idofn,jdofn)+fact3
                 jdofn=jdofn+1
                 fact3=fact2*gpcar(3,jnode,igaus)
                 elmat(idofn,jdofn)=elmat(idofn,jdofn)+fact3
              end do
           end do
        end do
     end do
  end if

end subroutine nsi_elmlbc
