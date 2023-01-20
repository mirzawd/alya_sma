!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tur_elmdif(&
     pnode,gpvis,gpmut,gpcar,eledd,gpdif,gpgrd)
  !------------------------------------------------------------------------
  !****f* Turbul/tur_elmdif
  ! NAME
  !   tur_elmdif
  ! DESCRIPTION
  !    Compute diffusion and diffusion gradient
  ! OUTPUT 
  !    GPDIF .......... k 
  !    GPGRD(NDIME) ... grad(k) coefficient  ! it is calculated in tur_elmco2, this should be cleaned
  ! USES
  ! USED BY
  !***
  !------------------------------------------------------------------------
  use def_kintyp, only     :  ip,rp
  use def_domain, only     :  ndime
  use def_turbul, only     :  iunkn_tur,kfl_algor_tur,&
       &                      param_tur
  implicit none
  integer(ip), intent(in)  :: pnode
  real(rp),    intent(in)  :: eledd(pnode)
  real(rp),    intent(in)  :: gpvis,gpmut
  real(rp),    intent(out) :: gpdif(kfl_algor_tur)
  real(rp),    intent(out) :: gpgrd(kfl_algor_tur,ndime)
  real(rp),    intent(in)  :: gpcar(ndime,pnode)
!  integer(ip)              :: idime,inode
  !
  ! Diffusion GPDIF = mu + mu_t/sig
  !
  if(kfl_algor_tur==1) then
     gpdif(1) = gpvis + gpmut/param_tur(iunkn_tur)
  else
     gpdif(1) = gpvis + gpmut/param_tur(1)
     gpdif(2) = gpvis + gpmut/param_tur(2)
  end if
  !
  ! Diffusion gradient GPGRD = grad(mu_t)
  !  
  !do idime = 1,ndime
  !   gpgrd(1,idime) = 0.0_rp
  !end do
  !do inode = 1,pnode
  !   do idime = 1,ndime
  !      gpgrd(1,idime) = gpgrd(1,idime)+gpcar(idime,inode)*eledd(inode)
  !   end do
  !end do
  !
  ! Diffusion gradient GPGRD = grad(mu_t)/sig
  !  
  !if(kfl_algor_tur==1) then
  !   !
  !   ! Decoupled
  !   !
  !   do idime=1,ndime 
  !      gpgrd(1,idime)=gpgrd(1,idime)/param_tur(iunkn_tur)
  !   end do
  !else
  !   !
  !   ! Coupled
  !   !
  !   do idime=1,ndime 
  !      gpgrd(2,idime)=gpgrd(1,idime)/param_tur(2)
  !      gpgrd(1,idime)=gpgrd(1,idime)/param_tur(1)
  !   end do     
  !end if

end subroutine tur_elmdif
