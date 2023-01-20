!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tur_addarr
  !----------------------------------------------------------------------
  !****f* Turbul/tur_addarr
  ! NAME 
  !    tur_addarr
  ! DESCRIPTION
  !    This routine computes additional arrays needed by the turbulence
  !    models. 
  !
  !    WALLD_TUR: the distance to the wall is needed for the following
  !    models and to perform the following tasks:
  !    --------------------------------------------
  !    Model                  #      bc  mut matrix
  !    --------------------------------------------
  !    Spalart-Allmaras       1      x    x    x
  !    Chien k-eps           13           x    x
  !    Lam-Bremhorst k-eps   14           x    x
  !    2-layer Rodi k-eps    15      x    x
  !    2-layer Xu-Chen k-eps 16      x    x
  !    Wilcox k-w            20      x
  !    --------------------------------------------
  !
  !    USTAR_TUR: the friction velocity
  !    --------------------------------------------
  !    Model                  #      bc  mut matrix
  !    --------------------------------------------
  !    Chien k-eps           13           x    x
  !    Two-layer Rodi k-eps  15      x    x    x
  !    Wilcox k-w            20      x              (only if rough wall)
  !    Bredberg k-w          21      x              (only if rough wall)
  !    --------------------------------------------
  !   
  !    LWNEI_TUR: the shortest node on the wall is needed for the following 
  !    models to compute y+:
  !    --------------------------------------------
  !    Model                  #      bc  mut matrix
  !    --------------------------------------------
  !    Chien k-eps           13                x
  !    2-layer Rodi k-eps    15      x         x     
  !    2-layer Xu-Chen k-eps 16      x         x     
  ! 
  ! OUTPUT
  !    LWNEI_TUR ... Nearest node on the wall
  !    WALLD_TUR ... Distance to nearest node on the wall
  ! USES
  ! USED BY
  !    tur_turnon
  !***
  !----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_kermod
  use def_domain
  use def_turbul
  use mod_communications, only : PAR_MAX
  use mod_arrays,         only : arrays_number
  implicit none
  integer(ip) :: ipoin,dummi
!  real(rp)    :: dista,dimin

  kfl_walld_tur = 0
  kfl_lwnei_tur = 0
  kfl_ustar_tur = 0
  kfl_grve2_tur = 0
  kfl_grsqk_tur = 0
  kfl_grono_tur = 0
  kfl_greps_tur = 0
  kfl_grk12_tur = 0
  kfl_grphi_tur = 0
  kfl_vorti_tur = 0
  kfl_avvel_tur = 0
  kfl_adapt_tur = 0
  kfl_fixn6_tur = 0
  kfl_fixn8_tur = 0

  !-------------------------------------------------------------------
  !
  ! Averaged velocity needed and adaptive condition?
  ! 
  !-------------------------------------------------------------------

  if( INOTMASTER ) then
     ipoin=0                                         ! Nodes with code 6    
     do while(ipoin<npoin)
        ipoin=ipoin+1
        if(lpoty(ipoin)/=0) then
           if(kfl_fixno_tur(1,ipoin,1)==6) then
              kfl_fixn6_tur=1
           end if
        end if
     end do
     ipoin=0                                         ! Nodes with code 8   
     do while(ipoin<npoin)
        ipoin=ipoin+1
        if(lpoty(ipoin)/=0) then
           if(kfl_fixno_tur(1,ipoin,1)==8) then
              kfl_fixn8_tur=1
           end if
        end if
     end do
  end if
  call PAR_MAX(kfl_fixn6_tur)
  call PAR_MAX(kfl_fixn8_tur)

  !-------------------------------------------------------------------
  !
  ! Flags
  ! 
  !-------------------------------------------------------------------
  !
  ! u* = friction velocity
  !
  if(  TUR_K_EPS_CHIEN         .or.& 
       TUR_K_EPS_NAGANO_TAGAWA .or.&
       TUR_TWO_LAYER_RODI      .or.&
       TUR_TWO_LAYER_XU_CHEN) then
     kfl_ustar_tur=2
  else if((kfl_rough==0.and.(rough_dom<-zetur.or.rough_dom>zetur)).and.(&
       TUR_K_OMEGA             .or.&
       TUR_K_OMEGA_BREDBERG    .or.&
       TUR_SST_K_OMEGA)        .or.&
       postp(1) % npp_stepi (arrays_number('YPLUS'),0)/=0    .or.&
       postp(1) % npp_stepi (arrays_number('USTAR'),0)/=0) then
     kfl_ustar_tur=1
  end if
  call tur_memarr(5_ip)
  !                                               
  ! Nearsest wall node
  !
  if(  TUR_K_XU_CHIEN          .or.&
       TUR_K_EPS_CHIEN         .or.&
       TUR_K_EPS_NAGANO_TAGAWA .or.&
       TUR_TWO_LAYER_RODI      .or.&
       TUR_TWO_LAYER_XU_CHEN   ) then
     kfl_lwnei_tur=1
  end if
  !
  ! y = Wall distance
  !
  if(  TUR_SPALART_ALLMARAS               .or.&
       TUR_K_XU_CHIEN                     .or.&
       TUR_K_EPS_NAGANO_TAGAWA            .or.&
       TUR_K_EPS_CHIEN                    .or.&
       TUR_K_EPS_LAM_BREMHORST            .or.& 
       TUR_K_EPS_JAW_HWANG                .or.&
       TUR_TWO_LAYER_RODI                 .or.&
       TUR_TWO_LAYER_XU_CHEN              .or.&
       TUR_K_OMEGA                        .or.&
       TUR_SST_K_OMEGA                    .or.&
       TUR_TKE_SGS                        .or.&
       (TUR_K_EPS_PHI_F.and.inits_tur/=0) .or.&
       postp(1) % npp_stepi (arrays_number('YPLUS'),0)/=0               .or.&
       postp(1) % npp_stepi (arrays_number('YSTAR'),0)/=0) then
     kfl_walld_tur=1
  end if
  if( (kfl_fixn6_tur==1.or.kfl_fixn8_tur==1)&
       .and.(kfl_infl1_tur==2.or.kfl_infl2_tur==2&
       .or.(kfl_infl2_tur==1.and.turle_tur<0.0_rp))) then
     kfl_walld_tur=1
  end if
  if( (kfl_fixn6_tur==1.or.kfl_fixn8_tur==1)&
       .and.(kfl_infl1_tur==4.or.kfl_infl2_tur==4)) then
     kfl_walld_tur=1
  end if
  !
  ! grad(phi)
  !
  if( TUR_K_EPS_PHI_F ) then
     kfl_grphi_tur=1
     call tur_memarr(11_ip)
  end if
  !
  ! Vorticity: Omega
  !
  if( TUR_SST_K_OMEGA ) then
     kfl_vorti_tur=1
  end if
  !
  ! grad(u^2), grad(sqrt(k))^2
  !
  if(  TUR_K_EPS_LAUNDER_SHARMA ) then                 
     kfl_grve2_tur=1
     kfl_grsqk_tur=1
     call tur_memarr( 9_ip)
     call tur_memarr(10_ip)
  end if
  !
  ! Average velocity + adaptive b.c.
  !
  if(kfl_fixn6_tur==1) kfl_avvel_tur=1            ! Average velocity needed
  if(kfl_fixn8_tur==1) kfl_avvel_tur=1            ! Average velocity needed
  if(kfl_fixn8_tur==1) kfl_adapt_tur=1            ! Adaptive b.c.

  !-------------------------------------------------------------------
  !
  ! Distance to the wall + nearest wall node: WALLD_TUR and LWNEI_TUR
  ! 
  !-------------------------------------------------------------------

  if(kfl_walld_tur/=0.and.(kfl_walgo_tur==0.or.kfl_lwnei_tur==1)) then

     call runend('TUR_ADARR: THIS HAS NOT BEEN CODED IN PARALLEL')
     if (kfl_delta == 1) call runend('Tur_addarr: not ready for kfl_delta == 1 , variable wall law') ! needs to be adapted similarly to tur_walgen
     call runend('TUR_ADARR: THIS HAS NOT BEEN CODED IN PARALLEL nor in serial!!!!!')   ! I commented this lines to eliminate walld_tur 
!     call tur_memarr(1_ip)
!     call tur_memarr(2_ip) 
!     !
!     ! Define wall nodes: they are the nodes located on wall boundaries. A wall 
!     ! boundary is a boundary having at least one wall node
!     !
!     if(nmate>1) then
!        do ipoin=1,npoin
!           if(  kfl_fixno_tur(1,ipoin,1)==3.or.&
!                kfl_fixno_tur(1,ipoin,1)==4.or.&
!                kfl_fixno_tur(1,ipoin,1)==7) then
!              lwnei_tur(ipoin)=1
!           end if
!        end do
!        do ipoin=1,npoin
!           if(lwnei_tur(ipoin)==0) then
!              dimin=1.0e9
!              do jpoin=1,npoin
!                 if(lwnei_tur(jpoin)==1) then
!                    dista=0.0_rp
!                    do idime=1,ndime
!                       dista=dista+(coord(idime,ipoin)-coord(idime,jpoin))&
!                            &     *(coord(idime,ipoin)-coord(idime,jpoin))
!                    end do
!                    if(dista<dimin) then
!                       dimin=dista
!                    end if
!                 end if
!              end do
!              walld_tur(ipoin)=sqrt(dimin)+delta_tur
!           else
!              walld_tur(ipoin)=delta_tur
!           end if
!        end do
!     else
!        do ipoin=1,npoin
!           if(  kfl_fixno_tur(1,ipoin,1)==3.or.&
!                kfl_fixno_tur(1,ipoin,1)==4.or.&
!                kfl_fixno_tur(1,ipoin,1)==7) then
!              lwnei_tur(ipoin)=lpoty(ipoin)
!           end if
!        end do
!        do ipoin=1,npoin
!           if(lwnei_tur(ipoin)==0) then
!              dimin=1.0e9
!              do jpoin=1,npoin
!                 if(lwnei_tur(jpoin)>0) then
!                    dista=0.0_rp
!                    do idime=1,ndime
!                       dista=dista+(coord(idime,ipoin)-coord(idime,jpoin))&
!                            &     *(coord(idime,ipoin)-coord(idime,jpoin))
!                    end do
!                    if(dista<dimin) then
!                       dimin=dista
!                       lwnei_tur(ipoin)=-lpoty(jpoin)
!                    end if
!                 end if
!              end do
!              walld_tur(ipoin)=sqrt(dimin)+delta_tur
!           else
!              walld_tur(ipoin)=delta_tur
!           end if
!        end do
!     end if
!     !
!     ! Deallocate memory
!     !
!     if(kfl_lwnei_tur==0) then
!        call tur_memarr(6_ip)    
!     else
!        do ipoin=1,npoin
!           lwnei_tur(ipoin)=abs(lwnei_tur(ipoin))
!        end do
!     end if
  end if

  !-------------------------------------------------------------------
  !
  ! Put flag and allocate memory for projected functions (tur_updbcs)
  ! 
  !-------------------------------------------------------------------

  if( INOTMASTER ) then
     ipoin=0
     do while(ipoin<npoin)
        ipoin=ipoin+1
        if(kfl_fixno_tur(1,ipoin,1)==4) then
           if(TUR_FAMILY_K_OMEGA) then
              kfl_grono_tur = 1                                  ! Smooth 1/y**2
              if(kfl_wallw_tur==1) then
                 kfl_greps_tur = 1                               ! Smooth eps^{1/4}*y
                 kfl_grk12_tur = 1                               ! Smooth sqrt(k)/y
              end if
           else if(TUR_FAMILY_K_EPS.and.&
                .not.(TUR_K_EPS_CHIEN.and.TUR_K_XU_CHIEN)) then
              kfl_grk12_tur = 1                                  ! Smooth sqrt(k)/y
           end if
           ipoin = npoin
        end if
     end do
  end if
  call PAR_MAX(kfl_grono_tur)
  call PAR_MAX(kfl_greps_tur)
  call PAR_MAX(kfl_grk12_tur)
  call tur_memarr(8_ip)

  !-------------------------------------------------------------------
  !
  ! Check errors
  ! 
  !-------------------------------------------------------------------

  dummi = 0
  if( INOTMASTER ) then
     do ipoin=1,npoin
        if( (kfl_fixno_tur(1,ipoin,1)==3.or.kfl_fixno_tur(1,ipoin,1)==4).and.&
             lpoty(ipoin)==0.and.nmate==1) then
           dummi = 1
        end if
     end do
  end if
  call PAR_MAX(dummi)
  if( dummi /= 0 ) call runend('TUR_ADDARR: CANNOT IMPOSE WALL BC ON INTERIOR NODES')

end subroutine tur_addarr
