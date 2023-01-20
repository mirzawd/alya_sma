!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine ale_inibcs()
  !-----------------------------------------------------------------------
  !****f* Alefor/ale_inibcs
  ! NAME
  !    ale_inibcs
  ! DESCRIPTION
  !    This routine applied boundary conditions
  ! OUTPUT 
  ! USES
  ! USED BY
  !    ale_turnon
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_inpout
  use def_master
  use def_domain
  use def_alefor
  use mod_opebcs
  use mod_memory,         only : memory_alloca_min
  use mod_communications, only : PAR_INTERFACE_NODE_EXCHANGE
  use mod_alefor,         only : alefor_memory_allocate
  implicit none
  integer(ip)  :: ipoin,iboun,inodb,idime

  if( INOTMASTER ) then

     !-------------------------------------------------------------
     !
     ! Allocate memory
     !
     !-------------------------------------------------------------

     call alefor_memory_allocate('BOUNDARY CONDITIONS') 

     !-------------------------------------------------------------
     !
     ! Node codes
     !
     !-------------------------------------------------------------

     if( kfl_icodn > 0 ) then
        iffun     =  1
        kfl_funno => kfl_funno_ale
        kfl_funtn => kfl_funtn_ale
        ifloc     =  1
        ifbop     =  0
        kfl_fixrs => kfl_fixrs_ale
        kfl_fixno => kfl_fixno_ale
        bvess     => bvess_ale(:,:,1)
        tncod     => tncod_ale
        call reacod(10_ip)

     end if

     !-------------------------------------------------------------
     !
     ! Boundary codes
     !
     !-------------------------------------------------------------

     if( kfl_icodb > 0 ) then

        !kfl_fixbo => kfl_fixbo_ale
        !tbcod     => tbcod_ale
        !call reacod(20_ip)

     end if

     !------------------------------------------------------------------
     !
     ! Non-constant boundary conditions
     !
     !------------------------------------------------------------------

     do ipoin = 1,npoin
        bvess_ale(:,ipoin,2) = bvess_ale(:,ipoin,1)
     end do

     !-------------------------------------------------------------
     !
     ! Boundary to node fixity and function
     !
     !-------------------------------------------------------------

     call memgen(1_ip,npoin,0_ip)
     do iboun = 1,nboun
        if( kfl_fixbo_ale(iboun) == 1 ) then
           do inodb = 1,nnode(ltypb(iboun))
              ipoin = lnodb(inodb,iboun)
              kfl_fixno_ale(1,ipoin) = 1
              kfl_funno_ale(ipoin)   = kfl_funbo_ale(iboun)
           end do
        end if
     end do
     call PAR_INTERFACE_NODE_EXCHANGE(gisca,'SUM')
     do ipoin = 1,npoin
        if( gisca(ipoin) >= 1 ) then
           do idime = 1,ndime
              kfl_fixno_ale(idime,ipoin) = 1
           end do
        end if
     end do
     call memgen(3_ip,npoin,0_ip)

  else

     call memory_alloca_min(mem_modul(1:2,modul),'BVESS_ALE','ale_inibcs',bvess_ale) 

  end if

end subroutine ale_inibcs
