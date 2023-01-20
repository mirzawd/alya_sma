!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tur_elmadj(&
     pnode,plapl,pgaus,gpsha,gpcar,gphes,gprea,gpden,&
     gpdif,gpvel,gpgrd,gpadj)
  !-----------------------------------------------------------------------
  !****f* Turbul/tur_elmadj
  ! NAME
  !   tur_elmadj
  ! DESCRIPTION
  !    Compute the adjoint operator
  !    -[rho*a+grad(k)].grad(v) -k*Lapl(v) + r*v
  ! OUTPUT 
  !    GPADG ... Adjoint operator at Gauss point
  ! USES
  ! USED BY
  !    tur_elmope
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only     :  ip,rp
  use def_domain, only     :  ndime,ntens,mnode
  use def_turbul, only     :  kfl_advec_tur,kfl_algor_tur
  implicit none
  integer(ip), intent(in)  :: pnode,plapl,pgaus
  real(rp),    intent(in)  :: gpsha(pnode,pgaus)
  real(rp),    intent(in)  :: gpcar(ndime,mnode,pgaus)
  real(rp),    intent(in)  :: gphes(ntens,mnode,pgaus)
  real(rp),    intent(in)  :: gpdif(kfl_algor_tur,pgaus)
  real(rp),    intent(in)  :: gpgrd(kfl_algor_tur,ndime,pgaus)
  real(rp),    intent(in)  :: gprea(kfl_algor_tur,kfl_algor_tur,pgaus)
  real(rp),    intent(in)  :: gpden(pgaus)
  real(rp),    intent(in)  :: gpvel(ndime,pgaus)
  real(rp),    intent(out) :: gpadj(kfl_algor_tur,pnode,pgaus)
  integer(ip)              :: inode,idime,igaus
  real(rp)                 :: dummr

  if(kfl_algor_tur==1) then

     if(kfl_advec_tur/=0) then
        !
        ! r*v - [rho*a+grad(k)].grad(v)
        !
        do igaus=1,pgaus
           do inode=1,pnode
              gpadj(1,inode,igaus)=gprea(1,1,igaus)*gpsha(inode,igaus)
              do idime=1,ndime
                 gpadj(1,inode,igaus)=gpadj(1,inode,igaus)&
                      -(gpden(igaus)*gpvel(idime,igaus)+gpgrd(1,idime,igaus))&
                      *gpcar(idime,inode,igaus)
              end do
           end do
        end do

     else

        do igaus=1,pgaus
           do inode=1,pnode
              gpadj(1,inode,igaus)=gprea(1,1,igaus)*gpsha(inode,igaus)
              do idime=1,ndime
                 gpadj(1,inode,igaus)=gpadj(1,inode,igaus)&
                      -gpgrd(1,idime,igaus)*gpcar(idime,inode,igaus)
              end do
           end do
        end do

     end if

     if(plapl==1) then
        !
        ! -k*Lapl(v)
        !
        do igaus=1,pgaus
           do inode=1,pnode
              gpadj(1,inode,igaus)=gpadj(1,inode,igaus)&
                   -gpdif(1,igaus)*sum(gphes(1:ndime,inode,igaus))
           end do
        end do
     end if

  else

     if(kfl_advec_tur/=0) then
        !
        ! r*v - [rho*a+grad(k)].grad(v)
        !
        do igaus=1,pgaus
           do inode=1,pnode
              gpadj(1,inode,igaus)=gprea(1,1,igaus)*gpsha(inode,igaus)
              gpadj(2,inode,igaus)=gprea(2,2,igaus)*gpsha(inode,igaus)
              do idime=1,ndime
                 dummr=gpden(igaus)*gpvel(idime,igaus)
                 gpadj(1,inode,igaus)=gpadj(1,inode,igaus)&
                      -(dummr+gpgrd(1,idime,igaus))*gpcar(idime,inode,igaus)
                 gpadj(2,inode,igaus)=gpadj(2,inode,igaus)&
                      -(dummr+gpgrd(2,idime,igaus))*gpcar(idime,inode,igaus)
              end do
           end do
        end do

     else

        do igaus=1,pgaus
           do inode=1,pnode
              gpadj(1,inode,igaus)=gprea(1,1,igaus)*gpsha(inode,igaus)
              gpadj(2,inode,igaus)=gprea(2,2,igaus)*gpsha(inode,igaus)
              do idime=1,ndime
                 gpadj(1,inode,igaus)=gpadj(1,inode,igaus)&
                      -gpgrd(1,idime,igaus)*gpcar(idime,inode,igaus)
                 gpadj(2,inode,igaus)=gpadj(2,inode,igaus)&
                      -gpgrd(2,idime,igaus)*gpcar(idime,inode,igaus)
              end do
           end do
        end do

     end if

     if(plapl==1) then
        !
        ! -k*Lapl(v)
        !
        do igaus=1,pgaus
           do inode=1,pnode
              dummr=sum(gphes(1:ndime,inode,igaus))
              gpadj(1,inode,igaus)=gpadj(1,inode,igaus)&
                   -gpdif(1,igaus)*dummr
              gpadj(2,inode,igaus)=gpadj(2,inode,igaus)&
                   -gpdif(2,igaus)*dummr
           end do
        end do
     end if

  end if

end subroutine tur_elmadj
