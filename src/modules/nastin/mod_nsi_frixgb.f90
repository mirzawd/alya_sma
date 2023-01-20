!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Nastin
!> @{
!> @file    mod_nsi_frixgb.f90
!> @author  Radhakrishnan Sarath
!> @brief   Functions for machine learning
!> @details Functions to couple as well as generate input to the model
!-----------------------------------------------------------------------

module mod_nsi_frixgb
  use def_kermod
  use def_domain
  use def_nastin
  use mod_ker_proper 
  use def_master
  use iso_c_binding                   ! to interface with C
  use def_kintyp,                     only : ip,rp
  use mod_bouder
  implicit none
  real(rp), allocatable     :: tveno_all(:,:), local_re(:,:), dim_u(:,:), tvscd_all(:,:)
!  real(rp), allocatable     :: tveno_all(:,:)
  real(rp), allocatable     :: ustars(:,:) 
  real(rp), allocatable     :: delta_ml(:), dim_y(:), delsc_ml(:)!, local_re(:), dim_u(:) 
  real(rp)                  :: ustarb
#ifdef ML_C_BINDING
  interface
    subroutine frixgb(lcre,y,u,uplus,n) bind (c, name='frixgb')
!    subroutine frixgb(lcre,y,u,uplus,n,kfl) bind (c, name='frixgb')
    use, intrinsic :: iso_c_binding
    implicit none

!    integer(c_int), value :: kfl                   ! kfl_paral for debugging
    integer(c_int), value :: n                     ! array size
    real(c_double), dimension(n) :: lcre, y, u, uplus

    end subroutine frixgb
  end interface
  
#else
		public :: frixgb
#endif


 ! interface
 !   subroutine frixgb(y,u,nu,ustar,n) bind (c, name='frixgb')
 !   use, intrinsic :: iso_c_binding
 !   implicit none
 !
 !   integer(c_int), value :: n                     ! array size
 !   real(c_double), dimension(n) :: y, u, ustar
 !   real(c_double), value :: nu                    ! viscosity
!!    real(c_double), value :: y, u, nu             
!!    type(c_ptr), value :: ustar                   ! I used this when I was sending a single value
 !
 !   end subroutine frixgb
 ! end interface

  contains
  

#ifndef ML_C_BINDING
   subroutine frixgb(lcre,y,u,uplus,n)
!   subroutine frixgb(lcre,y,u,uplus,n,kfl)
    use, intrinsic :: iso_c_binding
    implicit none

!    integer(c_int), intent(inout) :: kfl
    integer(c_int), intent(in)    :: n                   
    real(c_double), intent(inout) :: lcre(n), y(n), u(n), uplus(n)

  end subroutine frixgb
#endif

!    subroutine nsi_ml_ustar_all(vec_size,list_boundaries,pnodb,pgaub)
    subroutine nsi_ml_ustar_all(vec_size,list_boundaries,vikin,pnodb,pgaub)
  implicit none
!    !
!    ! Input and output variables
!    !
  integer(ip), intent(in) :: vec_size                                      ! Vector_size, may be i can just use VECTOR_SIZE
  integer(ip), intent(in), dimension(vec_size) :: list_boundaries          ! List of boundaries
  integer(ip), intent(in) :: pnodb                                         ! Number of nodes
  integer(ip), intent(in) :: pgaub                                         ! Number of Gauss points
  !    real(rp),    intent(in) :: gbsha(pnodb)
  real(rp)                :: bovfi(ndime,pnodb)
  real(rp)                :: gbsha(pnodb)
  real(rp)                :: baloc(ndime,ndime)                            
  real(rp)                :: bocod(ndime,mnodb)
  real(rp)                :: eucta
  real(rp)                :: velex(3),vewal(3),velfi(3),velsc(3),vewml(3) 
  real(rp)                :: tvelo(3),tvefi(3),tvedi(3),tvesc(3),tvscd(3)
  real(rp)                :: vikin
  integer(ip)             :: idime,jdime
  integer(ip)             :: igaub,inodb,ipoin
  integer(ip)             :: pblty
  integer(ip)             :: iboun,ielem,kboun
  
  ! There is an ambiguity on ndimb
  do kboun = 1, vec_size
  iboun = list_boundaries(kboun)
  if (iboun > 0_ip) then
    pblty = ltypb(iboun)
    ielem = lelbo(iboun)
    do inodb = 1,pnodb
      ipoin = lnodb(inodb,iboun)
      if( kfl_coupl(ID_NASTIN,ID_ALEFOR) /= 0 ) then
        bovfi(1:ndime,inodb) = velom(1:ndime,ipoin)
      else
        bovfi(1:ndime,inodb) = 0.0_rp
      end if
    end do
    do igaub = 1,pgaub
      call bouder(pnodb,ndime,ndimb,elmar(pblty)%deriv(:,:,igaub),bocod,baloc,eucta)
      if ( kfl_boexc_ker(kfl_codbo(iboun)) == 1 ) then
      velex(1:ndime) = velel_ker(1:ndime,lexlo_ker(igaub,iboun))  ! exchange loctn velocity
      velsc(1:ndime) = velml_ker(1:ndime,lexml_ker(igaub,iboun))  ! exchange scale velocity
      end if
      do idime = 1,ndime
        vewal(idime) = velex(idime)
        vewml(idime) = velsc(idime)
        velfi(idime) = 0.0_rp
      end do
      do inodb = 1,pnodb
        do idime = 1,ndime            
          gbsha = elmar(pblty)%shape(1,igaub) 
          velfi(idime) = velfi(idime) &
     + gbsha(inodb) * bovfi(idime,inodb)
        end do
      end do
      do idime = 1,ndime     
        tvelo(idime) = vewal(idime)
        tvesc(idime) = vewml(idime)
        tvefi(idime) = velfi(idime)
        do jdime = 1,ndime
          tvelo(idime) = tvelo(idime)   &
     - baloc(idime,ndime) &
     * baloc(jdime,ndime) * vewal(jdime)
          tvesc(idime) = tvesc(idime)   &
     - baloc(idime,ndime) &
     * baloc(jdime,ndime) * vewml(jdime)
          tvefi(idime) = tvefi(idime)   &
     - baloc(idime,ndime) &
     * baloc(jdime,ndime) * velfi(jdime)
        end do
        tvedi(idime) = tvelo(idime) - tvefi(idime)
        tvscd(idime) = tvesc(idime) - tvefi(idime)
!        if (isnan(tvedi(idime))) then
!          write(2000+kfl_paral,*)"tvelo, tvefi, idime",tvelo,tvefi,idime
!        end if
      end do
      call vecnor(tvedi,ndime,tveno_all(igaub,iboun),2_ip)
      call vecnor(tvscd,ndime,tvscd_all(igaub,iboun),2_ip)
      if (kfl_delta == 0_ip) then
         delta_ml((iboun-1)*mgaub+igaub) = delta_dom     ! constant wall distance
         delsc_ml((iboun-1)*mgaub+igaub) = delta_sc      ! constant wall distance
      else if (kfl_delta == 1_ip) then
         delta_ml((iboun-1)*mgaub+igaub) = ywalb(iboun)  ! variable wall distance
         delsc_ml((iboun-1)*mgaub+igaub) = yscab(iboun)  ! to be modified
      end if
      if (tvscd_all(igaub,iboun) .lt. 1e-10) then 
         dim_u(igaub,iboun) = tveno_all(igaub,iboun) / 1.0
      else
         dim_u(igaub,iboun) = tveno_all(igaub,iboun) / tvscd_all(igaub,iboun)
      end if
      dim_y((iboun-1)*mgaub+igaub) = delta_ml((iboun-1)*mgaub+igaub) / delsc_ml((iboun-1)*mgaub+igaub)
      local_re(igaub,iboun)        = tveno_all(igaub,iboun)*delta_ml((iboun-1)*mgaub+igaub)/vikin
    end do
    end if
  end do
  end subroutine nsi_ml_ustar_all
end module mod_nsi_frixgb
