!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_outtan(kfl_varia_outtan_nsi)
  !------------------------------------------------------------------------
  !****f* Nastin/nsi_outtan
  ! NAME 
  !    nsi_outtan
  ! DESCRIPTION
  !    This routine computes the tangential component of the traction -  comes in gevec
  ! USES
  ! USED BY
  !    nsi_outvar
  !    nsi_averag
  ! output
  !    gevec
  !------------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_nastin
  use mod_nsi_schur_operations, only : nsi_solini
  use mod_memory
  use def_kermod
  use mod_gradie
  use mod_ker_proper,     only : ker_proper
  use mod_communications, only : PAR_INTERFACE_NODE_EXCHANGE
  use mod_ker_nsw_visc2,  only : ker_nod_nsw_visc_0
  use mod_projec,         only : projec_elements_to_nodes
  use mod_frivel,         only : frivel

  implicit none
  integer(ip),intent(in)  :: kfl_varia_outtan_nsi

  integer(ip)             :: ipoin,ibopo,idime,dummi
  real(rp)                :: n1,n2,n3,s11,s12,s22,s13,s23,s33,stnn,us2
  real(rp)                :: sn1,sn2,sn3,stx,sty,stz,rho,nu,u,ustar,yplus
  real(rp),    pointer    :: visco_nsw(:)
  real(rp),    pointer    :: visco_tmp(:)
  real(rp),    pointer    :: densi_tmp(:)

  nullify(visco_nsw)
  nullify(visco_tmp)
  nullify(densi_tmp)
  ! 
  ! Velocity gradients GRADV_NSI
  !
  if( INOTEMPTY ) then
     call memory_alloca(mem_modul(1:2,modul),'GRADV_NSI','nsi_outtan',gradv_nsi,ntens,npoin)
     call memory_alloca(mem_modul(1:2,modul),'VISCO_TMP','nsi_outtan',visco_tmp,npoin)
     call gradie(veloc(1:ndime,1:npoin,1),gradv_nsi)
     visco_tmp(1:npoin) = prope_nsi(2,1:npoin)  !viscosity
     if( kfl_noslw_ker /= 0 ) then
        call ker_nod_nsw_visc_0
        call memory_alloca(mem_modul(1:2,modul),'VISCO_NSW','nsi_outtan',visco_nsw,npoin)
        visco_nsw = 0.0_rp
        call projec_elements_to_nodes(el_nsw_visc,visco_nsw)
        visco_tmp = visco_tmp + visco_nsw
        !
        ! Turbulent viscosity was missing !!!! -- Oriol says it should be close to 0 with a good model
        ! for the moment I only put it for the nos slip wall law case 
        !
        visco_nsw =0.0_rp
        call ker_proper('TURBU','NPOIN',dummi,dummi,visco_nsw)
        !     I reuse visco_nsw   but now it is the turbulent viscosity
        !     call ker_proper('TURBU','PGAUS',dummi,list_elements_p,gpmut,pnode,pgaus,gpsha,gpcar)     ! mut  mod_nsi_element_operations
        !     call ker_proper('TURBU','PGAUS',dummi,ielem,gpmut,pnode,pgaus,gpsha,gpcar)  ! nsi_elmope_omp
        call memory_alloca(mem_modul(1:2,modul),'DENSI_TMP','nsi_outtan',densi_tmp,npoin)
        call ker_proper('DENSI','NPOIN',dummi,dummi,densi_tmp)
        visco_nsw = densi_tmp * visco_nsw
        call memory_deallo(mem_modul(1:2,modul),'DENSI_TMP','nsi_outtan',densi_tmp)
        visco_tmp = visco_tmp + visco_nsw
        call memory_deallo(mem_modul(1:2,modul),'VISCO_NSW','nsi_outtan',visco_nsw)
     end if
     do ipoin = 1,npoin   !mu grad (u) * (viscosity+deltavisc)
        gradv_nsi(1:ntens,ipoin) = gradv_nsi(1:ntens,ipoin) * visco_tmp(ipoin) 
     end do
     call memory_deallo(mem_modul(1:2,modul),'VISCO_TMP','nsi_outtan',visco_tmp)     
  end if
  !
  ! Obtains notra_nsi and massb_nsi for exchange location (but not no slip wall law) or two layer
  ! in the no slip wall law notra_nsi has already been obtained in nsi_nsw_averages from vafor_nsi
  ! if it enters here for the noslip wall case, it steps over the value the correct value that has been obtained before 
  !
  if ( (kfl_waexl_ker == 1_ip.and.kfl_noslw_ker == 0_ip) .or. kfl_twola_ker == 1_ip ) then 
     call nsi_solini(3_ip)     ! it does nothing 
     call nsi_bouope_all(1_ip) ! calculates traction on nodes (notra_nsi), and lumped mass matrix (massb_nsi)
     if( INOTMASTER) then
        call PAR_INTERFACE_NODE_EXCHANGE(notra_nsi,'SUM','IN MY CODE')
        call PAR_INTERFACE_NODE_EXCHANGE(massb_nsi,'SUM','IN MY CODE')
     end if
  end if
  !-------------------------------------------------------------------
  !
  ! Tangential traction on all nodes. Set to 0 on interior nodes
  !
  !-------------------------------------------------------------------

  do ipoin = 1,npoin
     if (kfl_wlawf_nsi(ipoin) >= 1_ip) then    ! wall law but not - no slip wall law
        !
        ! Wall law
        !
        if ( (kfl_waexl_ker == 1_ip .and. kfl_noslw_ker == 0_ip) .or. kfl_twola_ker == 1_ip ) then !exchange location(but not no slip wall law - treated below)  and I believe all methods
           ibopo = lpoty(ipoin)
           if( ibopo >= 1 ) then
              notra_nsi(1:ndime,ipoin) = notra_nsi(1:ndime,ipoin) / massb_nsi(ipoin)
              do idime=1,ndime
                 gevec(idime,ipoin) = notra_nsi(idime,ipoin)
              end do
           else
              do idime = 1,ndime
                 gevec(idime,ipoin) = 0.0_rp
              end do
           end if
        else           ! no exchange location
           ibopo = lpoty(ipoin)
           if( ibopo >= 1 ) then
              rho = prope_nsi(1,ipoin)
              nu  = prope_nsi(2,ipoin)/rho
              u   = 0.0_rp
              if ( (kfl_wlaav_ker == 1_ip) .and. (ittim .gt. 0_ip) ) then
                 do idime=1,ndime
                    u = u + avupo_ker(idime,ipoin)*avupo_ker(idime,ipoin)
                 end do
              else
                 do idime=1,ndime
                    u = u + veloc(idime,ipoin,1)*veloc(idime,ipoin,1)
                 end do
              end if
              u = sqrt(u)
              if( kfl_rough >  0 ) rough_dom = rough(ipoin) 
              if( kfl_delta == 1 ) then
                 call frivel(kfl_ustar,ywalp(ibopo),rough_dom,u,nu,ustar)
                 if(u/=0.0_rp) then
                    us2 = rho*ustar*ustar/u
                 else
                    us2 = 0.0_rp
                 end if
              else
                 call frivel(kfl_ustar,delta_nsi,rough_dom,u,nu,ustar)
                 yplus = delta_nsi * u / nu
                 if( yplus < 5.0_rp  .and. kfl_ustar == 0 ) then
                    us2 = rho * nu/ delta_nsi
                 else if( kfl_ustar == 1 ) then
                    ustar = 0.41_rp  / log((delta_nsi+rough_dom) / rough_dom)
                    us2 = rho * u * ustar*ustar
                 else if( kfl_ustar == 2 ) then ! ABL 2
                    !                       velfr = 0.41_rp / log((delta_nsi+rough_dom) / rough_dom)
                    !                       velf2 = kinen**0.5 *cmu_st**0.25 ! frveloc in terms of kinetic turbulent energy 

                    !                       fact1 = gbden * velfr * velf2
                 else
                    us2 = rho * ustar *ustar / u
                 end if
              end if
              if(u/=0.0_rp) then
                 do idime=1,ndime
                    gevec(idime,ipoin) = -us2*veloc(idime,ipoin,1)
                 end do
              end if
           else  ! ibopo =0
              do idime = 1,ndime
                 gevec(idime,ipoin) = 0.0_rp
              end do
           end if
        end if
     else
        !
        ! Up to the wall (no wall law & no slip wall law ) --Actually the no slip wall law could be moved upwards - with the only subtelty that (/ massb_nsi(ipoin)) must not be applied
        !
        if ( kfl_varia_outtan_nsi ==1 .and. kfl_noslw_ker == 1_ip ) then ! forces from variational reaction at nodes
           gevec(:,ipoin) = notra_nsi(:,ipoin)
        else  ! ! default (origina)l behaviour  - from velocity gradients
           ibopo = lpoty(ipoin)
           if( ibopo >= 1 ) then
              n1  = exnor(1,1,ibopo)
              n2  = exnor(2,1,ibopo)
              s11 = gradv_nsi(1,ipoin)
              s22 = gradv_nsi(2,ipoin)
              s12 = gradv_nsi(3,ipoin)
              sn1 = s11*n1+s12*n2
              sn2 = s12*n1+s22*n2
              if( ndime == 3 ) then
                 n3   = exnor(3,1,ibopo) 
                 s33  = gradv_nsi(4,ipoin)
                 s13  = gradv_nsi(5,ipoin)
                 s23  = gradv_nsi(6,ipoin)
                 sn1  = sn1+s13*n3
                 sn2  = sn2+s23*n3
                 sn3  = s13*n1+s23*n2+s33*n3
                 stnn = sn1*n1+sn2*n2+sn3*n3
                 stx  = sn1-stnn*n1
                 sty  = sn2-stnn*n2
                 stz  = sn3-stnn*n3
                 gevec(1,ipoin)=stx         
                 gevec(2,ipoin)=sty       
                 gevec(3,ipoin)=stz      
              else
                 stnn = sn1*n1+sn2*n2
                 stx  = sn1-stnn*n1
                 sty  = sn2-stnn*n2
                 gevec(1,ipoin)=stx         
                 gevec(2,ipoin)=sty       
              end if
           else
              do idime = 1,ndime
                 gevec(idime,ipoin) = 0.0_rp
              end do
           end if
        end if
     end if
  end do
  !
  ! Deallocate memory
  !
  call memory_deallo(mem_modul(1:2,modul),'GRADV_NSI','nsi_outtan',gradv_nsi)

end subroutine nsi_outtan

!------------------------------------------------------------------------ 
!
! - Velocity strain rates gravb(ntens,npoin)
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
