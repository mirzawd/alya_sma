!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_matndof(ff,xx,amatr,ia,ja,ndof,nbrows,ipass_out)
  !
  ! outputs reordered matrices
  !
  use def_kintyp
  implicit none

  integer(ip), intent(in)    :: ndof,nbrows
  integer(ip), intent(in)    :: ia(nbrows+1)
  real(rp),    intent(in)    :: amatr(ndof,ndof,ia(nbrows+1)-1)
  integer(ip), intent(in)    :: ja(ia(nbrows+1)-1)

  real(rp),    intent(inout)    :: xx(ndof*nbrows)
  real(rp),    intent(inout)    :: ff(ndof*nbrows)

  integer(ip),save           :: ipass = 0_ip
  integer(ip),intent(out)    :: ipass_out

  real(rp),allocatable       :: amatr_ndof(:),faux(:),xaux(:)
  integer(ip),allocatable    :: ja_ndof(:)
  integer(ip),allocatable    :: ia_ndof(:)

  integer(ip)          :: nz,neqn,i,nz_save
  integer(ip)          :: kfl_file_formated
  character(50)        :: char_aux

  ipass = ipass+1_ip
  ipass_out = ipass
 
  neqn = nbrows*ndof
  nz  = ia(nbrows+1)-1
  nz_save = nz

  allocate(amatr_ndof(ndof*ndof*nz))
  allocate(ja_ndof(ndof*ndof*nz))
  allocate(ia_ndof(ndof*nbrows+1))
  allocate(faux(ndof*nbrows))
  allocate(xaux(ndof*nbrows))

  call matreorder(ff,xx,amatr,ia,ja,ndof,nbrows,1_ip,&  ! 1_ip means 1x,1y,2x,2y
       amatr_ndof,ja_ndof,ia_ndof,faux,xaux)

  kfl_file_formated=1
  
  if (kfl_file_formated==1) then    ! formated
     
     ! output matrix
     write(char_aux,'(i1)')ipass
     char_aux='matrix'//trim(char_aux)//'.txt'
     print*,char_aux
     open(777,file=char_aux)
     nz=nz_save
     write(777,'(i15,1x,i15)')ndof*nbrows,nz*ndof*ndof
     ! ja=c_dom,matrix
     do i=1,nz*ndof*ndof
        write(777,'(i15,1x,e17.10)')ja_ndof(i),amatr_ndof(i)
     end do
     ! rr,zz
     do i=1,ndof*nbrows
        write(777,'(e17.10,1x,e17.10)')faux(i),xaux(i)
     end do
     ! ia=r_dom
     do i=1,ndof*nbrows+1
        write(777,'(i15)')ia_ndof(i)
     end do
     close(777)
     
  else
     
     ! output matrix
     write(char_aux,'(i1)')ipass
     char_aux='matrix'//trim(char_aux)//'.txt'
     print*,char_aux
     open(777,file=char_aux,form='unformatted')
     nz=nz_save
     write(777)ndof*nbrows,nz*ndof*ndof
     ! ja=c_dom,matrix
     write(777)ja_ndof
     write(777)amatr_ndof
     ! rr,zz
     write(777)faux
     write(777)xaux
     ! ia=r_dom
     write(777)ia_ndof
     close(777)

  end if

  deallocate(amatr_ndof)
  deallocate(ja_ndof)
  deallocate(ia_ndof)
  deallocate(faux)
  deallocate(xaux)

 end subroutine nsi_matndof
