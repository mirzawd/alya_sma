!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_geobcs()
  !------------------------------------------------------------------------
  !****f* Nastin/nsi_geobcs
  ! NAME 
  !    nsi_geobcs
  ! DESCRIPTION
  !    This routine interprets geometrical codes:
  !    reabcs(): KFL_GEOBO(NBOUN)
  !    =>
  !    geonor(): KFL_GEONO(NPOIN), SKCOS(NDIME,NDIME,NBOPO)
  !    =>
  !    nsi_geobcs(): KFL_FIXNO_NSI, KFL_FIXRS
  !                     
  !                  KFL_GEONO     FIXNO   FIXRS
  !    Prescribed ......... 10 :   111       0
  !                         11 :   111       0
  !    Freestream ......... 20 :   000       0 
  !    Wall law ........... 30 :   100      -3
  !                         31 :   101      -3
  !    Symmetry / Splip ... 40 :   100      -3
  !                         41 :   101      -3
  !    No slip / Wall ..... 50 :   111       0
  !    Free / Outflow ..... 60 :   000       0
  !
  !    FIXRS = -3 means that local basis uses SKCOS(:,:,IBOPO)
  !  
  ! USES
  !    nsi_fixgeo
  ! USED BY
  !    nsi_inibcs
  !***
  !------------------------------------------------------------------------
  use def_parame
  use def_inpout
  use def_master
  use def_kermod
  use def_domain
  use def_nastin 
  use mod_memchk
  implicit none
  integer(ip)  :: iboun,icode,ncode
  integer(ip)  :: kfl_value

  if( INOTMASTER ) then

     !-------------------------------------------------------------
     !
     ! Geometrical codes
     !
     !-------------------------------------------------------------

     if( kfl_geome == 1 ) then
        tncod => tgcod_nsi(1:)        
        !tncod => momod(ID_NASTIN) % tgcod(1:)
        param(1:4) = 0.0_rp
        kfl_value  = 0 
        icode = 10 ; call nsi_fixgeo( icode, 1_ip,1_ip,1_ip,  0_ip, kfl_value, param(2) )
        icode = 11 ; call nsi_fixgeo( icode, 4_ip,4_ip,4_ip,  0_ip, kfl_value, param(2) )
        icode = 20 ; call nsi_fixgeo( icode, 0_ip,0_ip,0_ip,  0_ip, kfl_value, param(2) )
        if( kfl_nopen_nsi == 1 ) then
           icode = 30 ; call nsi_fixgeo( icode, 0_ip,0_ip,0_ip,  0_ip, kfl_value, param(2) )
           icode = 31 ; call nsi_fixgeo( icode, 0_ip,0_ip,0_ip,  0_ip, kfl_value, param(2) )
        else
           icode = 30 ; call nsi_fixgeo( icode, 1_ip,0_ip,0_ip, -3_ip, kfl_value, param(2) )
           icode = 31 ; call nsi_fixgeo( icode, 1_ip,0_ip,1_ip, -3_ip, kfl_value, param(2) )
        end if
        icode = 40 ; call nsi_fixgeo( icode, 1_ip,0_ip,0_ip, -3_ip, kfl_value, param(2) )
        icode = 41 ; call nsi_fixgeo( icode, 1_ip,0_ip,1_ip, -3_ip, kfl_value, param(2) )
        icode = 50 ; call nsi_fixgeo( icode, 1_ip,1_ip,1_ip,  0_ip, kfl_value, param(2) )
        icode = 60 ; call nsi_fixgeo( icode, 0_ip,0_ip,0_ip,  0_ip, kfl_value, param(2) )
        do ncode = 1,tncod(1) % ncode
           icode     = tncod(1) % l(ncode) % lcode(1)
           kfl_value = tncod(1) % l(ncode) % kfl_value
           if( icode == 10 ) call nsi_fixgeo( icode, 1_ip,1_ip,1_ip,  0_ip, kfl_value, tncod(1) % l(ncode) % bvess )
           if( icode == 11 ) call nsi_fixgeo( icode, 1_ip,1_ip,1_ip,  0_ip, kfl_value, tncod(1) % l(ncode) % bvess )
           if( icode == 20 ) call nsi_fixgeo( icode, 0_ip,0_ip,0_ip,  0_ip, kfl_value, tncod(1) % l(ncode) % bvess )
           if( kfl_nopen_nsi == 1 ) then
              if( icode == 30 ) call nsi_fixgeo( icode, 0_ip,0_ip,0_ip,  0_ip, kfl_value, tncod(1) % l(ncode) % bvess )
              if( icode == 31 ) call nsi_fixgeo( icode, 0_ip,0_ip,0_ip,  0_ip, kfl_value, tncod(1) % l(ncode) % bvess )
           else
              if( icode == 30 ) call nsi_fixgeo( icode, 1_ip,0_ip,0_ip, -3_ip, kfl_value, tncod(1) % l(ncode) % bvess )
              if( icode == 31 ) call nsi_fixgeo( icode, 1_ip,0_ip,1_ip, -3_ip, kfl_value, tncod(1) % l(ncode) % bvess )
           end if
           if( icode == 40 ) call nsi_fixgeo( icode, 1_ip,0_ip,0_ip, -3_ip, kfl_value, tncod(1) % l(ncode) % bvess )
           if( icode == 41 ) call nsi_fixgeo( icode, 1_ip,0_ip,1_ip, -3_ip, kfl_value, tncod(1) % l(ncode) % bvess )
           if( icode == 50 ) call nsi_fixgeo( icode, 1_ip,1_ip,1_ip,  0_ip, kfl_value, tncod(1) % l(ncode) % bvess )
           if( icode == 60 ) call nsi_fixgeo( icode, 0_ip,0_ip,0_ip,  0_ip, kfl_value, tncod(1) % l(ncode) % bvess )
           if( ( icode == 20 ) .and. kfl_initi_nsi == 0 ) then
              kfl_initi_nsi      = 1
              velin_nsi          = 0
              velin_nsi(1:ndime) = tncod(1) % l(ncode) % bvess(1:ndime)
           end if
        end do
        do iboun = 1,nboun
           if( kfl_geobo(iboun) >= 30 .and. kfl_geobo(iboun) <= 39 ) then
              if( kfl_nopen_nsi == 1 ) then
                 kfl_fixbo_nsi(iboun) = 18
              else
                 kfl_fixbo_nsi(iboun) = 3
              endif
           else if( kfl_geobo(iboun) >= 60 .and. kfl_geobo(iboun) <= 69 ) then
              if( kfl_colev_nsi /= 0 ) then
                 kfl_fixbo_nsi(iboun) = 12
              end if
           end if
        end do
     end if

  end if

end subroutine nsi_geobcs
