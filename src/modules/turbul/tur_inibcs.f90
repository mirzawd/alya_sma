!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tur_inibcs()
  !-----------------------------------------------------------------------
  !****f* Turbul/tur_inibcs
  ! NAME 
  !    tur_inibcs
  ! DESCRIPTION
  !    This routine reads TURBUL boundary conditions.
  !
  !    The different codes for KFL_FIXNO_TUR(1,IPOIN) are:
  !    = 0 ... Free or initial
  !    = 1 ... Dirichlet: fixed value
  !    = 2 ... Nothing
  !    = 3 ... Wall law 
  !    = 4 ... Wall
  !    = 5 ... Prescribe eps knowing the turbulent length scale
  !    = 6 ... Inflow condition
  !    = 7 ... Adaptive depending on the angle between normal to the boundary and velocity. 
  !    = 8 ... Adaptive inflow condition
  ! USES
  ! USED BY
  !    tur_turnon
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_inpout
  use def_master
  use def_kermod
  use def_domain
  use def_turbul
  use mod_ker_regularization, only : inv_regul_k, inv_regul_e, kfl_regularization
  use mod_communications, only : PAR_MAX
  use mod_run_config
  implicit none
  integer(ip) :: ipoin,iturb,icode,ncode,kfl_value

  if( INOTMASTER ) then
     !
     ! Allocate memory for KFL_FIXNO_TUR and BVESS_TUR
     !
     call tur_membcs(1_ip)
     !-------------------------------------------------------------
     !
     ! Geometrical codes
     !
     !-------------------------------------------------------------
     if( kfl_geome == 1 ) then
        do iturb = 1,nturb_tur
           iunkn_tur  =  iturb
           tncod      => tgcod_tur(iturb:)
           kfl_value  =  0
           param(1:5) =  0.0_rp
           icode = 10 ; call tur_fixgeo( icode, 6_ip, kfl_value, param(2) )
           icode = 11 ; call tur_fixgeo( icode, 6_ip, kfl_value, param(2) )
           icode = 20 ; call tur_fixgeo( icode, 0_ip, kfl_value, param(2) )
           if (iturb ==1.and.kfl_ustar==2) then
              icode = 30 ;  call tur_fixgeo( icode, 0_ip, kfl_value, param(2) )
           else
              icode = 30 ;  call tur_fixgeo( icode, 3_ip, kfl_value, param(2) )
           end if
           
           icode = 31 ; call tur_fixgeo( icode, 3_ip, kfl_value, param(2) )
           icode = 40 ; call tur_fixgeo( icode, 0_ip, kfl_value, param(2) )
           icode = 41 ; call tur_fixgeo( icode, 0_ip, kfl_value, param(2) )
           icode = 50 ; call tur_fixgeo( icode, 4_ip, kfl_value, param(2) )
           icode = 60 ; call tur_fixgeo( icode, 0_ip, kfl_value, param(2) )
           do ncode = 1,tncod(1) % ncode
              icode     = tncod(1) % l(ncode) % lcode(1)
              kfl_value = tncod(1) % l(ncode) % kfl_value
              if( icode == 10 ) call tur_fixgeo( icode, 6_ip, kfl_value, param(2) )
              if( icode == 11 ) call tur_fixgeo( icode, 6_ip, kfl_value, param(2) )
              if( icode == 20 ) call tur_fixgeo( icode, 0_ip, kfl_value, param(2) )
              if (iturb ==1.and.kfl_ustar==2) then
                 if( icode == 30 ) call tur_fixgeo( icode, 0_ip, kfl_value, param(2) )
              else
                 if( icode == 30 ) call tur_fixgeo( icode, 3_ip, kfl_value, param(2) )
              end if
              if( icode == 31 ) call tur_fixgeo( icode, 3_ip, kfl_value, param(2) )
              if( icode == 40 ) call tur_fixgeo( icode, 0_ip, kfl_value, param(2) )
              if( icode == 41 ) call tur_fixgeo( icode, 0_ip, kfl_value, param(2) )
              if( icode == 50 ) call tur_fixgeo( icode, 4_ip, kfl_value, param(2) )
              if( icode == 60 ) call tur_fixgeo( icode, 0_ip, kfl_value, param(2) )
           end do
        end do
     end if

     !-------------------------------------------------------------
     !
     ! Node codes (I have just put it after geometrical as in nsi_initur)
     ! so that one can modify by hand - I am not sure if this is the only thing needed
     ! I will test it
     !
     !-------------------------------------------------------------

     if( kfl_icodn > 0 ) then
        ifloc     =  0           
        ifbop     =  0
        do iturb = 1,nturb_tur
           iffun     =  1
           kfl_funno => kfl_funno_tur(:,iturb)
           kfl_funtn => kfl_funtn_tur(:,iturb)
           kfl_fixno => kfl_fixno_tur(:,:,iturb)
           bvess     => bvess_tur(:,:,iturb)
           tncod     => tncod_tur(iturb:)
           call reacod(10_ip)
           iffun     =  0
        end do
     end if

     !-------------------------------------------------------------
     !
     ! Boundary codes
     !
     !-------------------------------------------------------------

     if( kfl_icodb > 0 ) then

        kfl_fixbo => kfl_fixbo_tur
        do iturb = 1, nturb_tur
           bvnat     => bvnat_tur(:,:,iturb)
           tbcod     => tbcod_tur(iturb:)
           call reacod(20_ip)
        end do        
     end if
     !-------------------------------------------------------------------
     !
     ! Use user boundary conditions 
     !
     !-------------------------------------------------------------------

     if( kfl_usrbc_tur /= 0 ) then
        call tur_usrbcs(kfl_usrbc_tur)
     end if

     !-------------------------------------------------------------------
     !
     ! Initial values
     !
     !-------------------------------------------------------------------

     if( kfl_initi_tur == 1 ) then
        do ipoin=1,npoin
           if(kfl_fixno_tur(1,ipoin,1)<=0) &
                bvess_tur(1,ipoin,1:nturb_tur) = xinit_tur(1:nturb_tur)
        end do
     end if
     !
     ! If DELTA_DOM is precribed, use DELTA_DOM instead of DELTA_TUR
     !
     if(delta_dom>zetur.or.delta_dom<-zetur) then
        delta_tur=delta_dom
     end if
     !
     ! If distance to the wall is negative, transform wall law to wall
     !
     if( (delta_tur < -zetur) .or. (kfl_delta == -1) ) then

        do ipoin=1,npoin
           do iturb=1,nturb_tur
              if(kfl_fixno_tur(1,ipoin,iturb)==3) then
                 kfl_fixno_tur(1,ipoin,iturb) = 4
                 bvess_tur(1,ipoin,iturb)     = 0.0_rp
              end if
           end do
        end do
     end if
     delta_tur=max(delta_tur,0.0_rp)
     if( kfl_delta == -1) delta_tur = 0.0_rp
     !
     ! For two-layer models, take off wall condition for epsilon
     !
     if(kfl_model_tur==15.or.kfl_model_tur==16) then
        do ipoin=1,npoin
           if(kfl_fixno_tur(1,ipoin,2)==4.or.kfl_fixno_tur(1,ipoin,2)==3) then
              kfl_fixno_tur(1,ipoin,2)=0
           else if(kfl_fixno_tur(1,ipoin,2)==1) then
           end if
        end do
     end if
     !
     ! Customer: CFDWind1.0 with base flow
     !

     if( run_config%customer == cust_cfwind0 ) then
        do ipoin = 1,npoin
           do iunkn_tur = 1,nturb_tur
              if( kfl_fixno_tur(1,ipoin,iunkn_tur) == 6 ) kfl_fixno_tur(1,ipoin,iunkn_tur) = 1
              bvess_tur(1,ipoin,iunkn_tur) = 0.0_rp
           end do          
        end do
     end if
     ! ABL2, wall law fixing grad k =0 over wall
     if (kfl_ustar==2) then
        iunkn_tur =1
        do ipoin =1, npoin
           if( kfl_fixno_tur(1,ipoin,iunkn_tur) == 3 )   kfl_fixno_tur(1,ipoin,iunkn_tur)=0
        end do
     end if

     if (kfl_regularization) then

           do ipoin = 1,npoin
              !  if (abs(bvess_tur(1, ipoin, iunkn_tur)).gt.1.0e-10_rp) &
              bvess_tur(1,ipoin,1) =  inv_regul_k(bvess_tur(1, ipoin, 1))
              bvess_tur(1,ipoin,2) =  inv_regul_e(bvess_tur(1, ipoin, 2))
              
           end do
     else if (kfl_logva==1) then
        do iunkn_tur = 1,nturb_tur
           do ipoin = 1,npoin
              if (abs(bvess_tur(1, ipoin, iunkn_tur)).gt.1.0e-10_rp) &
                   bvess_tur(1,ipoin,iunkn_tur) = log( bvess_tur(1, ipoin, iunkn_tur))
           end do
        end do
     end if
  end if
  !
  ! exists fixity 7 ?  - I had to take it out of if master so that master enters into parari 
  !
  kfl_exist_fixi7_tur = 0
  do iturb = 1,nturb_tur
     do ipoin = 1,npoin
        if ( kfl_fixno_tur(1,ipoin,iturb) == 7 ) kfl_exist_fixi7_tur = 1
     end do
  end do
  call PAR_MAX(kfl_exist_fixi7_tur)
  
end subroutine tur_inibcs
