!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Kernel
!> @{
!> @file    mod_ker_nsw_viscf90
!> @author  Herbert Owen
!> @brief   Obtains viscosity for no slip wall law 
!> @details - 
!>          - 
!> @} 
!-----------------------------------------------------------------------
module mod_ker_nsw_visc
#include "def_vector_size.inc"
  use def_kintyp,                   only : ip,rp
  use def_master,                   only : ittim,kfl_paral    ! solucion trucha
  use def_kermod,                   only : kfl_nswel_ker,avwei_ker,avta1_nsw_ker
  use def_kermod,                   only : lnsw_exch,kfl_rough,rough_dom,kfl_ustar
  use def_domain,                   only : lnods,ltype,nnode,rough,lpoty,mnode
  use mod_frivel,                   only : frivel

  implicit none
  private
  public :: ker_nsw_visc

  contains   

subroutine ker_nsw_visc(ndime,pnode,pgaus,list_elements,elavv,gpcar,elnnsw,gpvis,gpden,gpmut,elibopo,elywal,elcod,gpvis_nsw)   ! elcod just for debugging REMOVE

  integer(ip), intent(in)    :: ndime,pnode,pgaus
  integer(ip), intent(in)    :: list_elements(VECTOR_SIZE)
  real(rp),    intent(in)    :: elavv(VECTOR_SIZE,ndime,pnode)
  real(rp),    intent(in)    :: gpcar(VECTOR_SIZE,ndime,mnode,pgaus)
  real(rp),    intent(in)    :: elnnsw(VECTOR_SIZE,ndime)
  real(rp),    intent(in)    :: gpvis(VECTOR_SIZE,pgaus)
  real(rp),    intent(in)    :: gpden(VECTOR_SIZE,pgaus)
  real(rp),    intent(in)    :: gpmut(VECTOR_SIZE,pgaus)
  real(rp),    intent(in)    :: elibopo(VECTOR_SIZE,pnode)
  real(rp),    intent(in)    :: elywal(VECTOR_SIZE)
  real(rp),    intent(in)    :: elcod(VECTOR_SIZE,ndime,pnode)
  real(rp),    intent(out)   :: gpvis_nsw(VECTOR_SIZE,pgaus)


  integer(ip)      :: kount(VECTOR_SIZE)
  real(rp)         :: gpgnavv(VECTOR_SIZE,ndime,pgaus)  ! Gauss Point Gradient in Normal dir of AVerage Velocity 
  real(rp)         :: elgnavv(VECTOR_SIZE,ndime)        ! ELement Gradient in Normal dir of AVerage Velocity 
  real(rp)         :: elgnavvt(VECTOR_SIZE)             ! ELement Gradient in Normal dir of AVerage Tangent Velocity 

  real(rp)         :: auxvi(VECTOR_SIZE)
  real(rp)         :: tveno_aux(VECTOR_SIZE)
  real(rp)         :: auxde(VECTOR_SIZE)
  real(rp)         :: auxmut(VECTOR_SIZE)
  real(rp)         :: avelavv(VECTOR_SIZE,ndime)        ! AVerage ELement AVerage Velocity
  real(rp)         :: avta1_aux(VECTOR_SIZE,ndime)
  real(rp)         :: auxi(VECTOR_SIZE)
  real(rp)         :: velfr(VECTOR_SIZE)

  real(rp)         :: avtan_fric_grad_based(VECTOR_SIZE)    ! just the magnitude
  real(rp)         :: av_mu_mut(VECTOR_SIZE)


  
  real(rp)         :: rough_aux,vikin,tveno,kount1
  integer(ip)      :: inode,ivect,idime,jdime,igaus,ielem,ipoin,ibopo
  integer(ip),parameter      :: imethod = 4   !now the default will be method 4 that is teh one taht is working best


#ifdef OPENACC
#define DEF_VECT ivect
#else
#define DEF_VECT 1:VECTOR_SIZE
#endif

#ifdef OPENACC
    do ivect = 1,VECTOR_SIZE
#endif 
       !----------------------------------------------------------------------
       !
       ! Gauss point values
       !
       !----------------------------------------------------------------------
       !
       ! GPGNAVV = dj uav_i nj =  dj N_I_i nj  Uav_I  ! gauss point gradient in normal direction of the average velocity 
       !
       gpgnavv(DEF_VECT,:,:)   = 0.0_rp

       do igaus = 1,pgaus
          do inode = 1,pnode
             do idime = 1,ndime
                do jdime = 1,ndime
                   gpgnavv(DEF_VECT,idime,igaus) = gpgnavv(DEF_VECT,idime,igaus) + elavv(DEF_VECT,idime,inode) * &
                        gpcar(DEF_VECT,jdime,inode,igaus) * elnnsw(DEF_VECT,jdime)
                end do
             end do
          end do
       end do
       !
       ! Obtain the average over gauss points of the normal gradient of the average velocity for all gauss points
       ! I could merge this loop with the upper one and obtain elgnavv directly
       ! Also average density and viscosity
       !
       elgnavv  = 0.0_rp
       auxvi = 0.0_rp
       auxde = 0.0_rp
       auxmut = 0.0_rp
       do igaus = 1,pgaus
          do idime = 1,ndime
!             elgnavvn(DEF_VECT) = elgnavvn(DEF_VECT) + gpgnavv(DEF_VECT,idime,igaus) * elnnsw(DEF_VECT,idime)
             elgnavv(DEF_VECT,idime)  = elgnavv(DEF_VECT,idime)  + gpgnavv(DEF_VECT,idime,igaus)
          end do
          auxvi(DEF_VECT)  = auxvi(DEF_VECT)  + gpvis(DEF_VECT,igaus)
          auxde(DEF_VECT)  = auxde(DEF_VECT)  + gpden(DEF_VECT,igaus)
          auxmut(DEF_VECT) = auxmut(DEF_VECT) + gpmut(DEF_VECT,igaus)
       end do
       do idime = 1,ndime
          elgnavv(DEF_VECT,idime) = elgnavv(DEF_VECT,idime) / real(pgaus,rp)
       end do
       auxvi(DEF_VECT)  = auxvi(DEF_VECT)  / real(pgaus,rp)
       auxde(DEF_VECT)  = auxde(DEF_VECT)  / real(pgaus,rp)
       auxmut(DEF_VECT) = auxmut(DEF_VECT) / real(pgaus,rp)
       !
       ! Substract normal component to keep only tangential one
       !
       !
       auxi(DEF_VECT) = 0.0_rp
       do idime=1,ndime
          auxi(DEF_VECT) = auxi(DEF_VECT) + elgnavv(DEF_VECT,idime) * elnnsw(DEF_VECT,idime)
       end do
       do idime=1,ndime
          elgnavv(DEF_VECT,idime) = elgnavv(DEF_VECT,idime) - auxi(DEF_VECT) * elnnsw(DEF_VECT,idime)
       end do
       !
       ! obtain modulus of tangential component
       !
       auxi(DEF_VECT) = 0.0_rp
       do idime=1,ndime
          auxi(DEF_VECT) = auxi(DEF_VECT) + elgnavv(DEF_VECT,idime) * elgnavv(DEF_VECT,idime)
       end do
       elgnavvt(DEF_VECT) = sqrt (auxi(DEF_VECT))




       if (.false.) then   ! this part will no longer be used
          !
          ! Obtain the average velocity for non boundary nodes
          !
          avelavv(DEF_VECT,:) = 0.0_rp
          kount = 0_ip
          do inode = 1,pnode
             !
             ! I do not need an if ibopo > 0 that would not be suitable in the vector case
             ! Instead I usar elibopo elibopo(DEF_VECT,inode)  that is =1.0 if ibopo>0 and 0.0 else
             !
             do idime=1,ndime
                avelavv(DEF_VECT,idime) = avelavv(DEF_VECT,idime) + elavv(DEF_VECT,idime,inode) * ( 1.0_rp - elibopo(DEF_VECT,inode) ) 
             end do
             kount(DEF_VECT) = kount(DEF_VECT) + nint(elibopo(DEF_VECT,inode))
          end do

          do idime=1,ndime
             avelavv(DEF_VECT,idime) = avelavv(DEF_VECT,idime) / (real(kount,rp)+1.0e-30_rp)   ! OJO trucheada para safar de que estoy haciendo cuentas en elementso interiores al pedo
          end do                                                                               ! esta relacionado con que no se como poner un if paar que haga algunos elem y otrso no
          !
          ! Substract normal component to keep only tangential one
          !
          !
          auxi(DEF_VECT) = 0.0_rp
          do idime=1,ndime
             auxi(DEF_VECT) = auxi(DEF_VECT) + avelavv(DEF_VECT,idime) * elnnsw(DEF_VECT,idime)
          end do
          do idime=1,ndime
             avelavv(DEF_VECT,idime) = avelavv(DEF_VECT,idime) - auxi(DEF_VECT) * elnnsw(DEF_VECT,idime)
          end do
       else  ! NEW OPTION
          !This comes directly in lnsw_exch(ielem)%velav(1:ndime)  - obtained in nsi_wallav -- ojo se hace al final elpaso de tiempo en restrat habra que volver a calcularlo - ya lo corregi
          ! for the moment it is not vectorized - not sure how easy/practical it would be 
          do ivect = 1,VECTOR_SIZE
             ielem = list_elements(ivect)
             if (ielem/=0_ip) then
                avelavv(ivect,1:ndime) = lnsw_exch(ielem)%velav(1:ndime)
             else
                avelavv(ivect,1:ndime) = 0.0_rp
             end if
          end do
          !
          ! Substract normal component to keep only tangential one
          !
          !
          auxi(DEF_VECT) = 0.0_rp
          do idime=1,ndime
             auxi(DEF_VECT) = auxi(DEF_VECT) + avelavv(DEF_VECT,idime) * elnnsw(DEF_VECT,idime)
          end do
          do idime=1,ndime
             avelavv(DEF_VECT,idime) = avelavv(DEF_VECT,idime) - auxi(DEF_VECT) * elnnsw(DEF_VECT,idime)
          end do
       end if

#ifdef OPENACC
    end do
#endif
    !
    ! Obtain tange from wall_law.
    ! Compute U*: VELFR   this part is not vectorized for the moment  I would need a vector frivel , not difficult
    ! also Time average of (mu+mut) d u_t / dn  - to be used later
    !
    do ivect = 1,VECTOR_SIZE
       ielem = list_elements(ivect)
       if( ielem > 0 ) then
          if(kfl_nswel_ker(ielem) >0) then
             vikin = auxvi(ivect) / auxde(ivect)                           ! nu
             call vecnor(avelavv(ivect,1:ndime),ndime,tveno,2_ip)          ! |u_tan-u_fix_tan|
             if(tveno > 1.0e-20_rp) then
                tveno_aux(ivect) = tveno    ! I save it for later use
             else
                tveno_aux(ivect) = 1.0_rp
             end if

             if( kfl_rough == 0 ) then 
                !
                ! Constant roughness
                !
                rough_aux = rough_dom

             else if( kfl_rough > 0 ) then
                rough_aux = 0.0_rp
                kount1 = 0_ip
                do inode = 1,pnode
                   ipoin =  lnods(inode,ielem)
                   ibopo = lpoty(ipoin)
                   if (ibopo /= 0) then
                      kount1 = kount1 + 1_ip
                      rough_aux = rough_aux + rough(ipoin)
                   end if
                end do
                rough_aux = rough_aux / kount1
!             else
!                rough_aux = 1.0_rp   ! this value will not be used    - but debug might complain if no value is given
             end if
             
             call frivel(kfl_ustar,elywal(ivect),rough_aux,tveno,vikin,velfr(ivect))      ! U*
             !
             ! Time average of (mu+mut) d u_t / dn  - to be used later
             !
             avta1_aux(ivect,:) = avta1_nsw_ker(:,kfl_nswel_ker(ielem))

          else   ! Element is not boundary element
             velfr(ivect) = 0.0_rp
             tveno_aux(ivect) = 1.0_rp   ! just some value so that it does not give problems
             avta1_aux(ivect,:) = 0.0_rp 
          end if
       else   ! Element number is null
          velfr(ivect) = 0.0_rp
          tveno_aux(ivect) = 1.0_rp      ! just some value so that it does not give problems
          avta1_aux(ivect,:) = 0.0_rp
       end if
    end do

#ifdef OPENACC
    do ivect = 1,VECTOR_SIZE
#endif

       
       !
       ! Tangential force that comes from the friction velocity - gradient based because it has already been transformed using fact
       !
       avtan_fric_grad_based(DEF_VECT) =   velfr(DEF_VECT) * velfr(DEF_VECT) * auxde(DEF_VECT)
       !
       ! divide by (du/dy)   
       !
       gpvis_nsw(DEF_VECT,1) = ( avtan_fric_grad_based(DEF_VECT) / (abs(elgnavvt(DEF_VECT))+1.0e-30_rp) )
       !
       ! Here I am going to sustract average mu+mut --- av_mu_mut
       ! to obtain it I first obtain modulus of avta1_aux(DEF_VECT,idime)
       !
       av_mu_mut(DEF_VECT) = avta1_aux(DEF_VECT,1) * avta1_aux(DEF_VECT,1) +  avta1_aux(DEF_VECT,2) * avta1_aux(DEF_VECT,2)
       if (ndime==3) av_mu_mut(DEF_VECT) = av_mu_mut(DEF_VECT) +  avta1_aux(DEF_VECT,3) * avta1_aux(DEF_VECT,3)
       av_mu_mut(DEF_VECT) = sqrt(av_mu_mut(DEF_VECT)) / (abs(elgnavvt(DEF_VECT))+1.0e-30_rp)
       !
       ! Now I substract -- moreover I take into account that the total traction can not be smaller than the one coming from mu+mut
       ! ESTE ES EL PASO QUE STAB FALTANDO (el max) - aca se supone que ambas tienne la misma direccion pero creo que es algo lógico
       ! pero se podría repensar
       !
       gpvis_nsw(DEF_VECT,1) = gpvis_nsw(DEF_VECT,1) - av_mu_mut(DEF_VECT)
       !
       ! max between  - gpmut(DEF_VECT,igaus) and the calculated values that for the moment is in igaus =1  
       !
       do igaus=2,pgaus   
          gpvis_nsw(DEF_VECT,igaus) = max ( - gpmut(DEF_VECT,igaus)  , gpvis_nsw(DEF_VECT,1) )
       end do
       gpvis_nsw(DEF_VECT,1) = max ( - gpmut(DEF_VECT,1)  , gpvis_nsw(DEF_VECT,1) )
       
#ifdef OPENACC
    end do
#endif

  end subroutine ker_nsw_visc

end module mod_ker_nsw_visc
