!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine matreorder(ff,xx,amatr,ia,ja,ndof,nbrows,korder,amatr_ndof,ja_ndof,ia_ndof,faux,xaux)
  !
  ! generates amatr_ndof,ia_ndof,ja_ndof
  ! it is a matrix obtained as if each direction of the velocity is a separate unknown.
  ! not using submatrices of size ndof*ndof as is used in the input matrix
  ! if ndof == 1 it just copies it
  !
  use def_kintyp
  implicit none

  integer(ip), intent(in)    :: ndof,nbrows  
  integer(ip), intent(in)    :: ia(nbrows+1)
  real(rp),    intent(in)    :: amatr(ndof,ndof,ia(nbrows+1)-1)
  integer(ip), intent(in)    :: ja(ia(nbrows+1)-1)
  real(rp),    intent(in)    :: xx(ndof*nbrows)
  real(rp),    intent(in)    :: ff(ndof*nbrows)
  integer(ip), intent(in)    :: korder

  real(rp),    intent(out)   :: amatr_ndof(ndof*ndof*ia(nbrows+1)-1)
  integer(ip), intent(out)   :: ja_ndof(ndof*ndof*ia(nbrows+1)-1)
  integer(ip), intent(out)   :: ia_ndof(ndof*nbrows+1)
  real(rp),    intent(out)   :: faux(ndof*nbrows)
  real(rp),    intent(out)   :: xaux(ndof*nbrows)
  
  integer(ip)          :: iz,jj,ii_ndof,iz_ndof,j,idof,irow
  integer(ip)          :: kk1,kk2,nz,neqn,ii,i,nz_save  
  
  neqn = nbrows*ndof
  nz  = ia(nbrows+1)-1
  nz_save = nz
  
  if( ndof > 1 ) then

     ii_ndof = 0
     iz_ndof = 0

     ! korder = 1 ->  1x,1y,2x,2y......npx,npy
     ! 2 ->  1x,2x,...npx,1y,2y,...npy
     if(korder==1) then 
        do ii = 1,nbrows
           do kk2 = 1,ndof
              do iz = ia(ii),ia(ii+1)-1
                 jj = ja(iz)
                 do kk1 = 1,ndof
                    iz_ndof = iz_ndof + 1
                    amatr_ndof(iz_ndof) = amatr(kk1,kk2,iz)
                    ja_ndof(iz_ndof)    = (jj-1)*ndof + kk1
                 end do
              end do
           end do
           nz = ia(ii+1)-ia(ii)
           do kk1 = 1,ndof
              ii_ndof = ii_ndof + 1
              ia_ndof(ii_ndof) = nz*ndof
           end do
        end do

        kk1        = ia_ndof(1)
        ia_ndof(1) = 1 
        do ii = 2,neqn+1
           kk2         = ia_ndof(ii)
           ia_ndof(ii) = ia_ndof(ii-1) + kk1
           kk1         = kk2
        end do
        faux=ff
        xaux=xx

     else if (korder==2) then
        !
        ! Check of diag =0
        !
        if (1==2) then
           do iz=1,nz_save
              if (abs(amatr(1,2,iz))>1.0e-12) print*,'12,iz,amatr(1,2,iz)',iz,amatr(1,2,iz)
              if (abs(amatr(1,3,iz))>1.0e-12) print*,'13,iz,amatr(1,3,iz)',iz,amatr(1,3,iz)
              if (abs(amatr(2,3,iz))>1.0e-12) print*,'23,iz,amatr(2,3,iz)',iz,amatr(2,3,iz)
              if (abs(amatr(2,1,iz))>1.0e-12) print*,'21,iz,amatr(2,1,iz)',iz,amatr(2,1,iz)
              if (abs(amatr(3,1,iz))>1.0e-12) print*,'31,iz,amatr(3,1,iz)',iz,amatr(3,1,iz)
              if (abs(amatr(3,2,iz))>1.0e-12) print*,'32,iz,amatr(3,2,iz)',iz,amatr(3,2,iz)
           end do
        end if
        do kk2 = 1,ndof
           do ii = 1,nbrows
              do kk1 = 1,ndof
                 do iz = ia(ii),ia(ii+1)-1
                    jj = ja(iz)
                    iz_ndof = iz_ndof + 1
                    amatr_ndof(iz_ndof) = amatr(kk1,kk2,iz)
                    ja_ndof(iz_ndof)    = jj + (kk1-1) *nbrows
                 end do
              end do  !kk1 = 1,ndof
              
              nz = ia(ii+1)-ia(ii)
              ii_ndof = ii_ndof + 1
              ia_ndof(ii_ndof) = nz*ndof
              
           end do  !ii = 1,nbrows
        end do  !kk2 = 1,ndof

        kk1        = ia_ndof(1)
        ia_ndof(1) = 1
        do ii = 2,neqn+1
           kk2         = ia_ndof(ii)
           ia_ndof(ii) = ia_ndof(ii-1) + kk1
           kk1         = kk2
        end do
        do idof=1,ndof
           do irow=1,nbrows
              i=idof+(irow-1)*ndof
              j=irow+(idof-1)*nbrows
              faux(j)=ff(i)
              xaux(j)=xx(i)
!              print*,'i,j,ff(i),faux(j),xx(i),xaux(j)',i,j,ff(i),faux(j),xx(i),xaux(j)
           end do
        end do
     end if


  else if ( ndof == 1 ) then   ! just copy 

     amatr_ndof(1:nz_save) = amatr(1,1,1:nz_save)
     ja_ndof(1:nz_save)    = ja(1:nz_save)
     ia_ndof(1:nbrows+1)   = ia(1:nbrows+1)
     faux(1:nbrows)       = ff(1:nbrows)
     xaux(1:nbrows)       = xx(1:nbrows)

  end if

 end subroutine matreorder
