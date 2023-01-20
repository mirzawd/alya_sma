!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine lev_inibcs()
  !-----------------------------------------------------------------------
  !****f* wavequ/lev_inibcs
  ! NAME 
  !    lev_inibcs
  ! DESCRIPTION
  !    This routine reads the boundary conditions
  ! USES
  !    lev_membc
  ! USED BY
  !    lev_turnon
  !***
  !-----------------------------------------------------------------------
  use def_parame, only       :  ip,rp
  use def_inpout
  use def_master, only       :  INOTMASTER
  use def_master, only       :  mem_modul
  use def_master, only       :  modul
  use def_domain
  use def_levels
  use mod_opebcs
  use mod_memory
  implicit none
  integer(ip) :: icode,ncode

  if( INOTMASTER ) then

     !-------------------------------------------------------------
     !
     ! Allocate memory
     !
     !-------------------------------------------------------------

     call memory_alloca(mem_modul(1:2,modul),'KFL_FIXNO_LEV','lev_inibcs',kfl_fixno_lev,1_ip,npoin)
     call memory_alloca(mem_modul(1:2,modul),'BVESS_LEV'    ,'lev_inibcs',bvess_lev,1_ip,npoin,1_ip)
     if ( kfl_conbc_lev == 0 ) then
        call memory_alloca(mem_modul(1:2,modul),'KFL_FUNNO_LEV','lev_inibcs',kfl_funno_lev,npoin)
     end if
     
     !-------------------------------------------------------------
     !
     ! Geometrical codes
     !
     !-------------------------------------------------------------

     if( kfl_geome == 1 ) then

        tncod => tgcod_lev   
        param(1:2) = 0.0_rp
        icode = 10 ; call lev_fixgeo( icode, 8_ip, param(2) )
        do ncode = 1,tncod(1) % ncode
           icode = tncod(1) % l(ncode) % lcode(1)
           if( icode == 10 ) call lev_fixgeo( icode, 8_ip, tncod(1) % l(ncode) % bvess(1) )
        end do

     end if

     !-------------------------------------------------------------
     !
     ! Node codes
     !
     !-------------------------------------------------------------

     if( kfl_icodn > 0 ) then

        if( kfl_conbc_lev == 0 ) then
           iffun      =  1
           kfl_funno  => kfl_funno_lev
        else
           iffun      =  0
        end if
        kfl_fixno => kfl_fixno_lev
        bvess     => bvess_lev(:,:,1)
        tncod     => tncod_lev
        call reacod(10_ip)

     end if

  end if

end subroutine lev_inibcs

subroutine lev_fixgeo(icode,ifixx,bvalu)
  !------------------------------------------------------------------------
  !****f* Levels/lev_fixgeo
  ! NAME 
  !    lev_fixgeo
  ! DESCRIPTION
  ! USES
  ! USED BY
  !    lev_inibcs
  !***
  !------------------------------------------------------------------------
  use def_parame
  use def_inpout
  use def_master
  use def_domain
  use def_levels
  implicit none
  integer(ip), intent(in) :: icode,ifixx
  real(rp),    intent(in) :: bvalu
  integer(ip)             :: ipoin,ibopo

  do ipoin = 1,npoin
     if( lpoty(ipoin) > 0 ) then
        ibopo = lpoty(ipoin)
        if( kfl_geono(ibopo) == icode ) then
           ibopo                  = lpoty(ipoin)
           kfl_fixno_lev(1,ipoin) = ifixx
           bvess_lev(1,ipoin,1)   = bvalu
        end if
     end if
  end do

end subroutine lev_fixgeo
