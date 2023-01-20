!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tur_frivel()
  !-----------------------------------------------------------------------
  !****f* Turbul/tur_frivel
  ! NAME 
  !    tur_frivel
  ! DESCRIPTION
  !    This routine computes the friction velocity on the wall
  ! USES
  !
  ! USED BY
  !    tur_begite
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_kermod
  use mod_ker_proper
  use def_domain
  use def_turbul
  use mod_gradie
  use mod_frivel, only : frivel
  implicit none
  integer(ip) :: ipoin,ibopo,dummi
  real(rp)    :: n1,n2,n3,s11,s12,s22,s13,s23,s33,stnn
  real(rp)    :: sn1,sn2,sn3,stx,sty,stz,nu,tveno,twall,rho(1),mu(1)

  if( kfl_ustar_tur/=0 .and. INOTMASTER ) then
     !
     ! Compute velocity strain rates GEVEC = grad(U)
     !
     call memgen(zero,ntens,npoin)
     call gradie(advec(1:ndime,1:npoin,1),gevec)

     if(kfl_ustar_tur/=0) then
        !
        ! Compute friction velocity USTAR_TUR
        !     
        do ipoin=1,npoin
           ibopo=lpoty(ipoin)
           if(ibopo>=1) then
              call ker_proper('DENSI','IPOIN',ipoin,dummi,rho(1:))
              call ker_proper('VISCO','IPOIN',ipoin,dummi,mu(1:))
              if(kfl_fixno_tur(1,ipoin,1)==3) then
                 !
                 ! USTAR_TUR using wall law
                 !
                 nu = mu(1)/rho(1)
                 call vecnor(advec(:,ipoin,1),ndime,tveno,2_ip)
                 if( kfl_rough > 0 ) rough_dom = rough(ipoin) 
                 if (kfl_delta == 1 ) then
                    call frivel(kfl_ustar,ywalp(ibopo),rough_dom,tveno,nu,ustar_tur(ibopo))
                 else
                    call frivel(kfl_ustar,delta_tur,rough_dom,tveno,nu,ustar_tur(ibopo))
                 end if
              else
                 !
                 ! USTAR_TUR = sqrt(tau_wall/rho)
                 !
                 n1  = exnor(1,1,ibopo)
                 n2  = exnor(2,1,ibopo)
                 s11 = mu(1)*gevec(1,ipoin)
                 s22 = mu(1)*gevec(2,ipoin)
                 s12 = mu(1)*gevec(3,ipoin)
                 sn1 = s11*n1+s12*n2
                 sn2 = s12*n1+s22*n2
                 if(ndime==3) then
                    n3    = exnor(ndime,1,ibopo) 
                    s33   = mu(1)*gevec(4,ipoin)
                    s13   = mu(1)*gevec(5,ipoin)
                    s23   = mu(1)*gevec(6,ipoin)
                    sn1   = sn1+s13*n3
                    sn2   = sn2+s23*n3
                    sn3   = s13*n1+s23*n2+s33*n3
                    stnn  = sn1*n1-sn2*n2-sn3*n3
                    stx   = sn1-stnn*n1
                    sty   = sn2-stnn*n2
                    stz   = sn3-stnn*n3
                    twall = sqrt(stx*stx+sty*sty+stz*stz)
                 else
                    stnn  = sn1*n1+sn2*n2
                    stx   = sn1-stnn*n1
                    sty   = sn2-stnn*n2
                    twall = sqrt(stx*stx+sty*sty)
                 end if
                 ustar_tur(ibopo)=sqrt(twall/rho(1))
              end if
           end if
        end do
     end if
     !
     ! Deallocate memory
     !
     call memgen(two,ntens,npoin)

  end if

end subroutine tur_frivel
!------------------------------------------------------------------------ 
!
! - Velocity strain rates gravv(ntens,npoin)
!   gradv_nsi(1,ipoin)=mu*(du/dx+du/dx)     
!   gradv_nsi(2,ipoin)=mu*(dv/dy+dv/dy)     
!   gradv_nsi(3,ipoin)=mu*(du/dy+dv/dx)     
!   gradv_nsi(4,ipoin)=mu*(dw/dz+dw/dz)     
!   gradv_nsi(5,ipoin)=mu*(du/dz+dw/dz)     
!   gradv_nsi(6,ipoin)=mu*(dv/dz+dw/dy)
! 
! In the case of the law of the wall, on walls:
! 
! (sig.n).t = rho*(U*)^2
! 
! Otherwise:      
!  
! Let sij=1/2(ui,j+uj,i), the traction is given by      
!         +-                                               -+
!         | -p*n1 + 2*mu*s11*n1 + 2*mu*s12*n2 + 2*mu*s13*n3 |
! sig.n = | -p*n2 + 2*mu*s12*n1 + 2*mu*s22*n2 + 2*mu*s23*n3 |
!         | -p*n3 + 2*mu*s13*n1 + 2*mu*s23*n2 + 2*mu*s33*n3 |
!         +-                                               -+
!       = (sn1,sn2,sn3)
! 
! The tangential component of the traction is computed as:
! (sig.n).t = sig.n-(n.sig.n)n
!           = (sn1,sn2,sn3)-(sn1*n1+sn2*n2+sn3*n3)(n1,n2,n3)
!             +-                                       -+
!             | sn1 - sn1*n1*n1 - sn2*n2*n1 - sn3*n3*n1 |
!           = | sn2 - sn1*n1*n2 - sn2*n2*n2 - sn3*n3*n2 |
!             | sn3 - sn1*n1*n3 - sn2*n2*n3 - sn3*n3*n3 |
!             +-                                       -+
! 
! NB: the pressure does not intervene in the tangential stress
!     and grave(j,i)=mu*(ui,j+uj,i) so that sni=grave(j,i)*nj      
!***
!-----------------------------------------------------------------------
