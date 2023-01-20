!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tem_membcs(itask)
  !-----------------------------------------------------------------------
  !****f* Temper/tem_membcs
  ! NAME
  !    tem_membcs
  ! DESCRIPTION
  !    Allocate memory for the physical problem
  ! OUTPUT 
  ! USES
  ! USED BY
  !    tem_reaphy
  !***
  !-----------------------------------------------------------------------
  use def_master
  use def_domain
  use def_temper
  use mod_memory, only : memory_alloca
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: ipoin,icomp

  select case ( itask )

  case (1_ip )
     !
     ! Fixity and boundary values
     !
     call memory_alloca(mem_modul(1:2,modul),'KFL_FIXNO_TEM','tem_membcs',kfl_fixno_tem,1_ip,npoin)
     call memory_alloca(mem_modul(1:2,modul),'KFL_FIXBO_TEM','tem_membcs',kfl_fixbo_tem,nboun)
     do ipoin=1,npoin
        kfl_fixno_tem(1,ipoin)=-1
     end do
     if(kfl_conbc_tem==0) then
        icomp = 2_ip
     else
        icomp = 1_ip
     end if
     call memory_alloca(mem_modul(1:2,modul),'BVESS_TEM','tem_membcs',bvess_tem,1_ip,npoin,icomp)
     call memory_alloca(mem_modul(1:2,modul),'BVNAT_TEM','tem_membcs',bvnat_tem,npnat_tem,nboun,icomp)
     if( kfl_regim_tem == 4 ) then 
        call memory_alloca(mem_modul(1:2,modul),'BVTEM_TEM','tem_membcs',bvtem_tem,1_ip,npoin,icomp)
     end if
     
  case ( 2_ip )

     !
     ! Non-constant b.c.'s : Functions
     !
     call memory_alloca(mem_modul(1:2,modul),'KFL_FUNNO_TEM','tem_membcs',kfl_funno_tem,npoin)
     call memory_alloca(mem_modul(1:2,modul),'KFL_FUNBO_TEM','tem_membcs',kfl_funbo_tem,nboun)
     call memory_alloca(mem_modul(1:2,modul),'KFL_FUNTN_TEM','tem_membcs',kfl_funtn_tem,npoin)
     call memory_alloca(mem_modul(1:2,modul),'KFL_FUNTB_TEM','tem_membcs',kfl_funtb_tem,nboun)

  case ( 4_ip )

     call memory_alloca(mem_modul(1:2,modul),'KFL_FUNTY_TEM','tem_membcs',kfl_funty_tem,10_ip)
     call memory_alloca(mem_modul(1:2,modul),'FUNPA_TEM','tem_membcs',funpa_tem,6_ip,10_ip)

  end select

end subroutine tem_membcs
