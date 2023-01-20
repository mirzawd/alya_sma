!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tur_solsgs()
  !-----------------------------------------------------------------------
  !****f* Nastin/tur_solsgs
  ! NAME 
  !    tur_solsgs
  ! DESCRIPTION
  !    This routine solves the SGS equation
  ! USES
  ! USED BY
  !    tur_endite
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_domain
  use def_master
  use def_turbul
  use mod_memory
  implicit none
  integer(ip) :: ipoin
  real(rp)    :: time1,time2

  if( kfl_ortho_tur >= 1 ) then
     !
     ! Initialization
     !
     call cputim(time1)
     !
     ! Residual projections
     !
     if( INOTMASTER ) then

        call memory_alloca(mem_modul(1:2,modul),'TUPR2_TUR','tur_solsgs',tupr2_tur,npoin)         !
        if ( kfl_ortho_tur==2) then ! Split oss
           call memory_alloca(mem_modul(1:2,modul),'TUPRR_TUR','tur_solsgs',tuprr_tur,npoin)     
        end if

        ! Update SGS
        !
        call tur_elmop2(4_ip)
        !
        ! Residual projections
        !
        call rhsmod(1_ip,tupr2_tur)        
        do ipoin = 1,npoin
           unpro_tur(iunkn_tur,ipoin) = tupr2_tur(ipoin) / vmass(ipoin)
        end do
        call memory_deallo(mem_modul(1:2,modul),'TUPR2_TEM','tur_solsgs',tupr2_tur) 
        
        if ( kfl_ortho_tur==2) then ! Split oss
           call rhsmod(1_ip,tuprr_tur)        
           do ipoin = 1,npoin
              unprr_tur(iunkn_tur,ipoin) = tuprr_tur(ipoin) / vmass(ipoin)
           end do
        call memory_deallo(mem_modul(1:2,modul),'TUPRR_TUR','tur_solsgs',tuprr_tur) 
        end if
     end if

     call cputim(time2)
     !cputi_tur(3) = cputi_tur(3) + time2 - time1

  end if
  if (kfl_shock_tur /=0) then
     !
     ! Initialization
     !
     call cputim(time1)
     !
     ! Residual projections
     !
     if( INOTMASTER ) then
        call memory_alloca(mem_modul(1:2,modul),'TUPGR_TUR','tur_solsgs',tupgr_tur,ndime, npoin)  
!        call memory_alloca(mem_modul(1:2,modul),'UNTUR','tur_memall',untur,nturb_tur,npoin,ncomp_tur)

        ! Update gradient projection 
        !
        call tur_elmop2(4_ip)
        !
        ! Residual projections
        !
        call rhsmod(ndime,tupgr_tur)        
        do ipoin = 1,npoin
           unpgr_tur(iunkn_tur,1:ndime,ipoin) = tupgr_tur(1:ndime,ipoin) / vmass(ipoin)
        end do
        call memory_deallo(mem_modul(1:2,modul),'TUPGR_TUR','tur_solsgs',tupgr_tur) 
        
        
     end if

     call cputim(time2)
  end if

end subroutine tur_solsgs
