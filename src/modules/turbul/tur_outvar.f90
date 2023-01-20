!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tur_outvar(ivari,imesh)
  !-----------------------------------------------------------------------
  !****f* Turbul/tur_outvar
  ! NAME
  !   tur_outvar
  ! DESCRIPTION
  !    Compute and output turbulent variables
  ! USES
  !    postpr
  ! USED BY
  !    tur_output
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_turbul 
  use def_kermod
  use mod_ker_proper,         only : ker_proper
  use mod_ker_regularization, only : regul_k, regul_e, kfl_regularization
  use mod_outvar,             only : outvar 
  use mod_arrays,             only : arrays_name
  implicit none
  integer(ip),  intent(in) :: ivari
  integer(ip),  intent(in) :: imesh   !< Mesh to postprocess on
  integer(ip)              :: ipoin,ibopo,idime,icont,iline,jpoin,dummi
  real(rp)                 :: xfact,rho(1),mu(1),auxi,nu,dummr

  select case ( arrays_name(ivari) )  

  case ( 'KEY  ' )
     !
     ! k: turbulent kinetic energy
     !
     if(TUR_K_XU_CHIEN.or.TUR_FAMILY_K_EPS.or.TUR_FAMILY_K_OMEGA.or.TUR_TKE_SGS) then
        if( INOTMASTER ) then
           call memgen(zero,npoin,zero)
           do ipoin = 1,npoin
              gesca(ipoin) = untur(1,ipoin,1)
           end do
           if (kfl_logva==1) gesca(1:npoin) = exp(gesca(1:npoin))
           if (kfl_regularization ) then
              do ipoin =1, npoin
                 gesca(ipoin) = regul_k(untur(1,ipoin,1))
              end do
           end if
        end if
     end if

  case ( 'EPSIL' )
     !
     ! eps: turbulent dissipation
     !
     if(TUR_K_XU_CHIEN) then
        if( INOTMASTER ) then
           call memgen(zero,npoin,zero)
           call tur_secvar(5_ip,gesca)
        end if

     else if(TUR_FAMILY_K_EPS) then
        if(TUR_K_EPS_CHIEN .and. INOTMASTER ) then
           call memgen(zero,npoin,zero)
           call tur_secvar(8_ip,gesca)  ! eps = eps0 + eps'
        else
           if( INOTMASTER ) then
              call memgen(zero,npoin,zero)
              do ipoin = 1,npoin
                 gesca(ipoin) = untur(2,ipoin,1)
              end do
              if (kfl_logva==1) gesca(1:npoin) = exp(gesca(1:npoin))
              if (kfl_regularization ) then
                 do ipoin =1, npoin
                    gesca(ipoin) = regul_e(untur(2,ipoin,1))
                 end do
              end if
           end if
        end if

     else if(TUR_FAMILY_K_OMEGA) then
        if( INOTMASTER ) then
           call memgen(zero,npoin,zero)
           call tur_secvar(1_ip,gesca)
        end if

     end if

  case ( 'OMEGA' )
     !
     ! w: turbulent frequency
     !
     if(TUR_FAMILY_K_OMEGA) then
        if( INOTMASTER ) then
           call memgen(zero,npoin,zero)
           do ipoin = 1,npoin
              gesca(ipoin) = untur(2,ipoin,1)
           end do
        end if
     else if(TUR_FAMILY_K_EPS) then
        if( INOTMASTER ) then
           call memgen(zero,npoin,zero)
           call tur_secvar(2_ip,gesca)
        end if
     end if

  case ( 'NUTIL' )
     !
     ! nu tilde (Spalart-Allmaras model)
     ! 
     if(TUR_SPALART_ALLMARAS) then
        if( INOTMASTER ) gesca => untur(1,1:npoin,1) 
     end if

  case ( 'LENGT' )
     ! 
     ! l: turbulent length scale
     !
     if(TUR_K_XU_CHIEN) then
        if( INOTMASTER ) then
           call memgen(zero,npoin,zero)
           call tur_secvar(6_ip,gesca)
        end if

     else if(TUR_FAMILY_K_EPS) then
        if( INOTMASTER ) then
           call memgen(zero,npoin,zero)
           call tur_secvar(3_ip,gesca)
        end if

     else if(TUR_FAMILY_K_OMEGA) then
        if( INOTMASTER ) then
           call memgen(zero,npoin,zero)
           call tur_secvar(4_ip,gesca)
        end if

     end if

  case ( 'TURVI' )
     !
     ! mu_t: turbulent viscosity
     !
     if( INOTMASTER ) gesca  => turmu

  case ( 'YPLUS' )
     ! 
     ! y+: dimensionless distance to the wall
     !  
     if( INOTMASTER ) then
        call memgen(zero,npoin,zero)
        call tur_secvar(9_ip,gesca)
     end if

  case ( 'USTAR' )
     ! 
     ! U*: friction velocity
     !  
     if( INOTMASTER ) then
        call memgen(zero,npoin,zero)
        do ipoin=1,npoin
           ibopo=lpoty(ipoin)
           if(ibopo>=1) then
              gesca(ipoin)=ustar_tur(ibopo)
           else
              gesca(ipoin)=0.0_rp
           end if
        end do
     end if

  case ( 'YVSTA' )
     ! 
     ! yv*: y*sqrt(vv)/nu
     !  
     if( INOTMASTER ) then
        call memgen(zero,npoin,zero)
        if(TUR_K_XU_CHIEN) then
           call tur_secvar(7_ip,gesca)
        end if
     end if

  case ( 'YSTAR' )
     ! 
     ! y*: y*sqrt(k)/nu
     !  
     if( INOTMASTER ) then
        call memgen(zero,npoin,zero)
        do ipoin=1,npoin
           call ker_proper('DENSI','IPOIN',ipoin,dummi,rho(1:))
           call ker_proper('VISCO','IPOIN',ipoin,dummi,mu(1:))
           gesca(ipoin) = walld(ipoin)*sqrt(untur(1,ipoin,1))*rho(1)/mu(1)
        end do
     end if

  case ( 'RESI1' )
     ! 
     ! UNTUR(1) difference
     !  
     if( INOTMASTER ) then
        call memgen(zero,npoin,zero)
        do ipoin=1,npoin
           gesca(ipoin)=abs(untur(1,ipoin,1)-unold_tur(1,ipoin))
        end do
     end if

  case ( 'RESI2' )
     ! 
     ! UNTUR(2) difference
     !  
     if( INOTMASTER ) then
        call memgen(zero,npoin,zero)
        do ipoin=1,npoin
           gesca(ipoin)=abs(untur(nturb_tur,ipoin,1)-unold_tur(nturb_tur,ipoin))
        end do
     end if

  case ( 'RESIT' )
     ! 
     ! TURMU difference
     !  
     if( INOTMASTER ) then
        call memgen(zero,npoin,zero)
        do ipoin=1,npoin
           gesca(ipoin)=abs(turmu(ipoin)-unold_tur(nturb_tur+1,ipoin))
        end do
     end if

  case ( 'PHI  ' )
     !
     ! PHI or V2
     !
     if( TUR_K_EPS_V2_F .or. TUR_K_EPS_PHI_F ) then
        gesca => untur(3,1:npoin,1) 
     end if

  case ( 'F    ' )
     !
     ! F
     !
     if( TUR_K_EPS_V2_F .or. TUR_K_EPS_PHI_F ) then
        gesca => untur(4,1:npoin,1) 
     end if

  case ( 'GRSQK' )
     !
     ! F
     !
     if( TUR_K_EPS_LAUNDER_SHARMA ) then
        gesca => grsqk_tur
     end if

  case ( 'INTEN' )
     ! 
     ! Energy of turbulence I = sqrt(2/3*k)/u
     !  
     if( .not. TUR_SPALART_ALLMARAS ) then
        if( INOTMASTER ) then
           call memgen(zero,npoin,zero)
           do ipoin=1,npoin
              xfact=0.0_rp
              do idime=1,ndime
                 xfact=xfact+advec(idime,ipoin,1)*advec(idime,ipoin,1)
              end do
              xfact=sqrt(xfact)
              if(xfact==0.0_rp) then
                 gesca(ipoin)=0.0_rp
              else
                 gesca(ipoin)=sqrt(untur(1,ipoin,1)*2.0_rp/3.0_rp)/xfact
              end if
           end do
        end if
     end if

  case ( 'PROJE' )
     !
     ! L2 projection
     !
     if( kfl_ortho_tur >=1 ) then
        if( INOTMASTER ) gesca => unpro_tur(1,1:npoin)
     end if

  case ( 'LIMIT' )
     !
     ! Limiter
     !
     if( kfl_limit_tur /= 0 ) then
        if( INOTMASTER ) then
           do ipoin = 1,npoin
              rhsid(ipoin) = 0.0_rp
           end do
           iunkn_tur = 1
           call tur_elmop2(5_ip)
           call rhsmod(1_ip,rhsid)
           do ipoin = 1,npoin
              rhsid(ipoin) = rhsid(ipoin) / vmass(ipoin)
           end do
           gesca => rhsid
        end if
     end if

  case ( 'RATIO' )
     ! 
     ! Viscosity ratio: mut/mu
     !  
     if( INOTMASTER ) then
        call memgen(zero,npoin,zero)
        do ipoin=1,npoin
           call ker_proper('DENSI','IPOIN',ipoin,dummi,rho(1:))
           call ker_proper('VISCO','IPOIN',ipoin,dummi,mu(1:))
           gesca(ipoin)=turmu(ipoin)/mu(1)
        end do
     end if

  case ( 'UNKN1' )
     !
     ! Turbulent variable 1
     !
     if( INOTMASTER ) then
        call memgen(zero,npoin,zero)
        do ipoin=1,npoin
           gesca(ipoin)=untur(1,ipoin,1)
        end do
     end if

  case ( 'UNKN2' )
     !
     ! Turbulent variable 2
     !
     if( nturb_tur >= 2 ) then
        if( INOTMASTER ) then
           call memgen(zero,npoin,zero)
           do ipoin=1,npoin
              gesca(ipoin)=untur(2,ipoin,1)
           end do
        end if
     end if

  case ( 'UNKN3' )
     !
     ! Turbulent variable 3
     !
     if( nturb_tur >= 3 ) then
        if( INOTMASTER ) then
           call memgen(zero,npoin,zero)
           do ipoin=1,npoin
              gesca(ipoin)=untur(3,ipoin,1)
           end do
        end if
     end if

  case ( 'UNKN4' )
     !
     ! Turbulent variable 4
     !
     if( nturb_tur >= 4 ) then
        if( INOTMASTER ) then
           call memgen(zero,npoin,zero)
           do ipoin=1,npoin
              gesca(ipoin)=untur(4,ipoin,1)
           end do
        end if
     end if

  case ( 'LINEL' )
     !
     ! LINEL: Linelets of preconditioner CG
     !
     if( INOTMASTER ) then
        icont=0
        do ipoin=1,npoin
           rhsid(ipoin)=0.0_rp
        end do
        do iline=1,solve(5)%nline
           icont=icont+1
           do ipoin=solve(5)%lline(iline),solve(5)%lline(iline+1)-1
              jpoin=solve(5)%lrenup(ipoin)
              rhsid(jpoin)=real(icont,rp)
           end do
        end do
        gesca => rhsid
     end if

  case ( 'GROUP' )
     !
     ! GROUP: GROUPS FOR DEFLATED CG
     !
     if( INOTMASTER ) then
        do ipoin=1,npoin
           rhsid(ipoin)=real(solve(5)%lgrou(ipoin),rp)
        end do
        gesca => rhsid
     end if

  case ( 'FDDES' )
     !
     ! FDDES - blending factor used for DDES
     !
     if( kfl_ddesm_tur /= 0 ) then
        if( INOTMASTER ) then
           call rhsmod(1_ip,fddes_tur)
           do ipoin = 1,npoin
              fddes_tur(ipoin) = fddes_tur(ipoin) / vmass(ipoin)
           end do
           gesca => fddes_tur
        end if
     end if

  case ( 'GDDES' )
     !
     ! GDDES - blending function used for DDES - fd
     !
     if( kfl_ddesm_tur /= 0 ) then
        if( INOTMASTER ) then
           call rhsmod(1_ip,gddes_tur)
           do ipoin = 1,npoin
              gddes_tur(ipoin) = gddes_tur(ipoin) / vmass(ipoin)
           end do
           gesca => gddes_tur
        end if
     end if

  case ( 'SSTF1' )
     !
     ! SSTF1 - blending function used for SST
     !
     if( TUR_SST_K_OMEGA ) then
        if( INOTMASTER ) then
           call rhsmod(1_ip,sstf1_tur)
           do ipoin = 1,npoin
              sstf1_tur(ipoin) = sstf1_tur(ipoin) / vmass(ipoin)
           end do
           gesca => sstf1_tur
        end if
     end if

  case ( 'SSTF2' )
     !
     ! SSTF2 - blending function used for SST
     !
     if( TUR_SST_K_OMEGA ) then
        if( INOTMASTER ) gesca  => sstf2_tur
     end if

  case ( 'SASSO' )
     !
     ! SASSO - SAS source term for SST
     !
     if( TUR_SST_K_OMEGA ) then
        if( INOTMASTER ) then
           call rhsmod(1_ip,sasso_tur)
           do ipoin = 1,npoin
              sasso_tur(ipoin) = sasso_tur(ipoin) / vmass(ipoin)
           end do
           gesca => sasso_tur
        end if
     end if

  case ( 'AVKEY' )
     !
     !  AVKEY
     !
     if (cutim > avtim_tur) then
        auxi = cutim - avtim_tur
        if ( INOTMASTER ) then
           do ipoin= 1,npoin
              unkno(ipoin) = avkey_tur(ipoin) / auxi
              avkey_tur(ipoin) = 0.0_rp
           end do
           gesca => unkno
        endif
     else
        return
     endif

  case ( 'AVOME' )
     !
     !  AVOME
     !
     if (cutim > avtim_tur) then
        auxi = cutim - avtim_tur
        if ( INOTMASTER ) then
           do ipoin= 1,npoin
              unkno(ipoin)     = avome_tur(ipoin) / auxi
              avome_tur(ipoin) = 0.0_rp
           end do
           gesca => unkno
        endif
     else
        return
     endif

  case ( 'AVTVI' )
     !
     !  AVTVI
     !
     if (cutim > avtim_tur) then
        auxi = cutim - avtim_tur
        if ( INOTMASTER ) then
           do ipoin= 1,npoin
              unkno(ipoin)     = avtvi_tur(ipoin) / auxi
              avtvi_tur(ipoin) = 0.0_rp
           end do
           gesca => unkno
        endif
     else
        return
     endif

  case ( 'VORTI' )
     !
     !  VORTI
     !
     gesca => vorti_tur

  case ( 'CACAS' )
     !
     ! CACAS
     !
     if( INOTMASTER ) then
        call ker_proper('DENSI','IPOIN',ipoin,dummi,rho(1:))
        call ker_proper('VISCO','IPOIN',ipoin,dummi,mu(1:))
        do ipoin=1,npoin
           nu = mu(1)/rho(1) 
           call tur_nut2nd(3_ip,ipoin,nu,dummr,untur(1,ipoin,1),rhsid(ipoin))
        end do
        gesca => rhsid
     end if

  case ( 'ERROR' )
     !
     ! Manufactured solution
     !
     !if( INOTMASTER ) then
     !   call memgen(zero,npoin,zero)
     !   call tur_manufactured_nodal_error(gesca)
     !end if

  case ( 'FIXT1' )
     !
     ! KFL_FIXNO_TUR 1st variable
     !
     if( INOTMASTER ) then
        call memgen(zero,npoin,zero)
        do ipoin = 1,npoin
           gesca(ipoin) = real(kfl_fixno_tur(1,ipoin,1),rp)
        end do
     end if

  case ( 'FIXT2' )
     !
     ! KFL_FIXNO_TUR 2nd variable
     !
     if( INOTMASTER ) then
        call memgen(zero,npoin,zero)
        do ipoin = 1,npoin
           gesca(ipoin) = real(kfl_fixno_tur(1,ipoin,2),rp)
        end do
     end if

  case ( 'SENSM' )
     !
     ! SENSM: mesh sensitivities
     !
     if( INOTMASTER ) then
        call memgen(zero,ndime,npoin)
        do ipoin = 1,npoin
           do idime = 1,ndime
              gevec(idime,ipoin) = 0.0_rp
           end do
        end do
        do ipoin = 1,npoin
           do idime = 1,ndime
              gevec(idime,ipoin) = sens_mesh(idime,ipoin) 
           end do
        end do
     end if
     
  case ( 'MAXLE' ) 
     !
     ! Max_mixing_length
     !
     if( INOTMASTER ) then
        call memgen(zero,npoin,zero)
        do ipoin = 1,npoin
           gesca(ipoin) = tur_max_mixlen(ipoin)       
        end do
     end if
  end select
  !
  ! Postprocess
  !
  call outvar(&
       ivari,&
       ittim,cutim,postp(1) % wopos(:,ivari),MESH_ID=imesh)

end subroutine tur_outvar
