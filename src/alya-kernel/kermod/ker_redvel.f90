!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Kermod
!> @{
!> @file    ker_redvel.f90
!> @author  Hadrien Calmet
!> @date    05/07/2012
!> @brief   Compute reduced velocity
!> @details Compute reduced velocity
!> @}
!-----------------------------------------------------------------------
subroutine ker_redvel(vered)

  use def_parame
  use def_master
  use def_domain
  use mod_gradie
  use mod_memory


  implicit none
  real(rp),   intent(in) :: vered(ndime,*)
  integer(ip)            :: ipoin
  real(rp),   pointer    :: gradu(:,:,:)

  if( INOTMASTER ) then


     nullify(gradu)
     call memory_alloca(mem_modul(1:2,modul),'GRADU','redvel',gradu,ndime,ndime,npoin)
     call gradie(veloc(:,:,1),gradu)

     call runend('KER_REDVEL: THIS SHOULD NO LONGER BE USED')
     points: do ipoin = 1,npoin

        call hadri3(gradu(1,1,ipoin),veloc(1,ipoin,1),veloc(1,ipoin,3),ndime,vered(1,ipoin))

     end do points

     call memory_deallo(mem_modul(1:2,modul),'GRADU','redvel',gradu)

  end if

end subroutine ker_redvel

subroutine hadri3(gvelo,velo,veli,ndime,vered)

  use def_kintyp, only : ip,rp,lg
  use def_kermod

  implicit none

  integer(ip),intent(in)  :: ndime
  real(rp),   intent(in)  :: gvelo(ndime,ndime)
  real(rp),   intent(in)  :: velo(ndime)
  real(rp),   intent(in)  :: veli(ndime)
  real(rp),   intent(out) :: vered(ndime)
  real(rp)                :: accel(ndime)
  integer(ip)             :: ierr,jdime,idime,i
  real(rp)                :: jaco(ndime,ndime)
  real(rp)                :: autovec(ndime,ndime)
  real(rp)                :: lambda_r(ndime)
  real(rp)                :: lambda_i(ndime)
  real(rp)                :: fv1(ndime)
  real(rp)                :: fv2(ndime),m,vr(ndime),n(ndime),velon
  logical(lg)             :: matzz

  i = 0
  ierr = 0 ! TODO: it was left uninitialized

  do idime = 1,ndime
     lambda_r(idime) = 0.0_rp
     lambda_i(idime) = 0.0_rp
     fv1(idime)      = 0.0_rp
     fv2(idime)      = 0.0_rp
     vr(idime)       = 0.0_rp
     n(idime)        = 0.0_rp
     do jdime = 1,ndime
        jaco(idime,jdime)    = gvelo(idime,jdime)
        autovec(idime,jdime) = 0.0_rp
     end do
  end do


  matzz = .true.

  !call hadri4(&
  !     ndime,ndime,jaco,lambda_r,lambda_i,matzz,&
  !     autovec,fv1,ierr)

  if(      lambda_i(1) == 0.0_rp .and. lambda_i(2) == -lambda_i(3) ) then
     i = 1
  else if( lambda_i(2) == 0.0_rp .and. lambda_i(1) == -lambda_i(3) ) then
     i = 2
  else if( lambda_i(3) == 0.0_rp .and. lambda_i(1) == -lambda_i(2) ) then
     i = 3
  end if

  if( i /= 0 ) then
     !
     ! calcul of Vr reduced velocity
     !
     ! normalized real eigenvector n
     !
     m = sqrt(&
          autovec(i,1)*autovec(i,1) + &
          autovec(i,2)*autovec(i,2) + &
          autovec(i,3)*autovec(i,3) )
     !
     !
     !
     n(1) = autovec(i,1) / m
     n(2) = autovec(i,2) / m
     n(3) = autovec(i,3) / m
     !
     ! calcul del scalar product velo*n
     !
     if( kfl_vortx_thres == 0 .or. kfl_vortx_thres == 1) then
        !
        ! Eigen method
        !
        velon = (velo(1)*n(1)+velo(2)*n(2)+velo(3)*n(3))
        !
        vered(1) = velo(1) - velon * n(1)
        vered(2) = velo(2) - velon * n(2)
        vered(3) = velo(3) - velon * n(3)

     elseif ( kfl_vortx_thres == 2) then
        !
        ! Parallel Method
        !
        do idime=1,ndime
           accel(idime) = (gvelo(1,idime)*velo(1)+gvelo(2,idime)*velo(2)+gvelo(3,idime)*velo(3))+(velo(idime)-veli(idime))
        enddo
        !
        velon = (accel(1)*n(1)+accel(2)*n(2)+accel(3)*n(3))
        !
        vered(1) = accel(1) - velon * n(1)
        vered(2) = accel(2) - velon * n(2)
        vered(3) = accel(3) - velon * n(3)
        !
     end if
  else
     vered(1) = 0.0_rp
     vered(2) = 0.0_rp
     vered(3) = 0.0_rp
  end if

  if (ierr/=0) then
     write(*,*)' todo mal '
     stop
  end if

end subroutine hadri3
