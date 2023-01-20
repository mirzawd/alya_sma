!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tem_outhfl()
  !------------------------------------------------------------------------
  !****f* Temper/tem_outhfl
  ! NAME 
  !    tem_outhfl
  ! DESCRIPTION
  !    This routine computes the heat flux
  ! USES
  ! USED BY
  !    tem_output
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_temper
  use mod_postpr
  use mod_gradie
  use mod_memory,      only : memory_alloca, memory_deallo
  implicit none
  integer(ip)             :: ipoin,ibopo,idime
  real(rp), pointer       :: gradt(:,:)
  !
  ! Allocate memory
  ! 
  nullify(gradt)
  call memory_alloca(mem_modul(1:2,modul),'GRADT','tem_outhfl',gradt,ndime,npoin)
  !
  ! Compute temperature gradients
  !
 
  call tem_heatfl(gradt)

  !
  ! Compute heat flux
  !
  do ipoin=1,npoin
     ibopo=lpoty(ipoin)
     if(ibopo>=1) then
        do idime=1,ndime
           gesca(ipoin)=gesca(ipoin)&
                +gradt(idime,ipoin)*exnor(idime,1,ibopo)
        end do
     else
        gesca(ipoin)=0.0_rp
     end if
  end do
  !
  ! Deallocate memory
  !
  call memory_deallo(mem_modul(1:2,modul),'GRADT','tem_outhfl',gradt)

end subroutine tem_outhfl
 
