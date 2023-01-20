!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_solsgs(itask)
  !-----------------------------------------------------------------------
  !****f* Nastin/nsi_solsgs
  ! NAME 
  !    nsi_solsgs
  ! DESCRIPTION
  !
  !    This routine solves the SGS equation and the projection
  !
  !    Global    Elem. Gauss  | Full Oss       Split
  !    -----------------------+-----------------------------------------
  !    VEPRO_NSI ELVEP GPVEP  | -tau1'*R(u)    tau1'*rho*(a.grad)u
  !    PRPRO_NSI ELPRP GPPRP  | tau2'*div(u)   tau2'*div(u)  
  !    GRPRO_NSI ELGRP GPGRP  |    -           tau1'*( grad(p) - rho*f )
  !
  !    Using open rule, the projection of a linear function is NOT the
  !    exact linear function near the boundary
  !
  ! USES
  ! USED BY
  !    nsi_endite
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_domain
  use def_master
  use def_nastin
  use mod_memory
  use mod_nsi_elmope_all  
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: ipoin,imeth
  real(rp)                :: dummr

  if( kfl_sgsti_nsi /= 0 .or. kfl_sgsco_nsi /= 0 .or. kfl_stabi_nsi == NSI_OSS .or. kfl_stabi_nsi == NSI_SPLIT_OSS ) then

     select case ( itask ) 

     case ( 1_ip )
        !
        ! Initialization
        !
        if( kfl_sgsco_nsi == 1 ) itsta_nsi = 0
        resgs_nsi(1) = 0.0_rp
        resgs_nsi(2) = 0.0_rp
        rmsgs_nsi    = 0.0_rp
        if( kfl_sgsco_nsi == 1 .or. kfl_sgsti_nsi == 1 ) then
           itsta_nsi = 0
           resis_nsi = 0.0_rp
           resis_nsi = 0.0_rp
        end if

     case ( 2_ip )
        !
        ! Residual projections
        !
        call times(3) % ini() 

        if( INOTMASTER .and. kfl_stabi_nsi > 0 ) then

           call memory_alloca(mem_modul(1:2,modul),'VEPR2_NSI','nsi_solsgs',vepr2_nsi,ndime,npoin)
           call memory_alloca(mem_modul(1:2,modul),'PRPR2_NSI','nsi_solsgs',prpr2_nsi,npoin)
           if( kfl_stabi_nsi == 2 ) then
              call memory_alloca(mem_modul(1:2,modul),'GRPR2_NSI','nsi_solsgs',grpr2_nsi,ndime,npoin)
           end if

        end if
        !
        ! Update SGS and projections
        !
        imeth = 1
        if( INOTMASTER .and. kfl_sgscp_nsi == 0 ) then
           if( imeth == 2 ) call opeclo(1_ip)
           call nsi_elmope_all(4_ip) 
           if( imeth == 2 ) call opeclo(2_ip)
        end if
        !
        ! Residual projections
        !
        if( INOTEMPTY .and. kfl_stabi_nsi > 0 ) then

           call rhsmod(ndime,vepr2_nsi)
           call rhsmod(1_ip, prpr2_nsi)

           !
           ! This repetition of lines is not nice
           ! 2 options. a) use auxilially vmas_aux o b) put the lines inside 'if imeth' in a separate subroutine and call it. 
           !

           if( imeth == 1 ) then

              do ipoin = 1,npoin
                 dummr                    = 1.0_rp / vmass(ipoin)
                 vepro_nsi(1:ndime,ipoin) = vepr2_nsi(1:ndime,ipoin) * dummr 
                 prpro_nsi(ipoin)         = prpr2_nsi(ipoin) * dummr 
              end do
              if( kfl_stabi_nsi == 2 ) then
                 call rhsmod(ndime,grpr2_nsi)
                 do ipoin = 1,npoin
                    dummr                    = 1.0_rp / vmass(ipoin)
                    grpro_nsi(1:ndime,ipoin) = grpr2_nsi(1:ndime,ipoin) * dummr 
                 end do
              end if

           else

              do ipoin = 1,npoin
                 dummr                    = 1.0_rp / vmasc(ipoin)
                 vepro_nsi(1:ndime,ipoin) = vepr2_nsi(1:ndime,ipoin) * dummr 
                 prpro_nsi(ipoin)         = prpr2_nsi(ipoin) * dummr 
              end do
              if( kfl_stabi_nsi == 2 ) then
                 call rhsmod(ndime,grpro_nsi)
                 do ipoin = 1,npoin
                    dummr                    = 1.0_rp / vmasc(ipoin)
                    grpro_nsi(1:ndime,ipoin) = grpr2_nsi(1:ndime,ipoin) * dummr 
                 end do
              end if

           end if

           call memory_deallo(mem_modul(1:2,modul),'VEPR2_NSI','nsi_solsgs',vepr2_nsi)
           call memory_deallo(mem_modul(1:2,modul),'PRPR2_NSI','nsi_solsgs',prpr2_nsi)
           if( kfl_stabi_nsi == 2 ) then
              call memory_deallo(mem_modul(1:2,modul),'GRPR2_NSI','nsi_solsgs',grpr2_nsi)
           end if

        end if
        !
        ! Output SGS convergence
        !
        call nsi_cvgsgs()

        call times(3) % add() 

     end select

  end if

end subroutine nsi_solsgs
