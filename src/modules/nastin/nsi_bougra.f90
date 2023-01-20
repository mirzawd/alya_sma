!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> 
!> @author  david oks 
!> @date    2021-08-12
!> @brief   Compute number of Neumann boundary elements per node
!> @details Compute number of boundary elements per node, considering
!>          only Neumann boundary condition elements
!> 
!-----------------------------------------------------------------------

subroutine nsi_bougra()

  use def_kintyp,         only :  ip, rp,lg
  use def_nastin,         only :  kfl_fixbo_nsi, NSI_FRACTIONAL_STEP
  use def_nastin,         only :  bvnat_nsi
  use def_nastin,         only :  bpess_nsi,kfl_fixpp_nsi, bvnat_nsi
  use def_nastin,         only :  kfl_exist_fib02_nsi, kfl_fixpr_nsi
  use def_nastin,         only :  kfl_exist_fib20_nsi
  use def_master,         only :  modul, mem_modul, gesca
  use def_domain,         only :  npoin, nboun, ltypb, lnodb, meshe
  use def_domain,         only :  elmar
  use def_kermod,         only :  ndivi
  use def_elmgeo,         only :  element_type
  use mod_memory_basic,   only :  memory_alloca
  use mod_memory_basic,   only :  memory_deallo
  use mod_communications, only :  PAR_MAX
  use mod_communications, only :  PAR_INTERFACE_NODE_EXCHANGE
  use mod_projec,         only :  projec_boundaries_to_nodes
  implicit none

  integer(ip)          :: inodb, ipoin, iboun
  real(rp),    pointer :: pp(:)
  logical(lg), pointer :: bmask(:)

  nullify(bmask)
  !
  ! If pressure is prescribe:
  ! - Fix Laplacian un mod_nsi_multi_step_fs by puttin rc = 0
  ! - Fix pressure in nsi_updbcs using bpess_nsi and fixpp_nsi
  ! - Modify fixpr_nsi following fixpp_nsi
  !
  if( ( kfl_exist_fib02_nsi == 1 .or. kfl_exist_fib20_nsi == 1 ) .and. NSI_FRACTIONAL_STEP ) then

     call memory_alloca( mem_modul(1:2,modul),'BMASK','nsi_bougra',bmask,nboun)

     call memgen(0_ip,nboun,0_ip)
     do iboun = 1,nboun
        gesca(iboun) = bvnat_nsi(1,iboun,1)
        if( kfl_fixbo_nsi(iboun) == 2 .or. kfl_fixbo_nsi(iboun) == 20 ) then
           bmask(iboun) = .true.
           do inodb = 1,element_type(abs(ltypb(iboun))) % number_nodes
              ipoin = lnodb(inodb,iboun)
              kfl_fixpp_nsi(1,ipoin) = 1
           end do
        end if
     end do
     call PAR_INTERFACE_NODE_EXCHANGE(kfl_fixpp_nsi,'MAX')
     do ipoin = 1,npoin
        kfl_fixpr_nsi(1,ipoin) = max(kfl_fixpr_nsi(1,ipoin),kfl_fixpp_nsi(1,ipoin))
     end do
     pp => bpess_nsi(1,:,1) 
     call projec_boundaries_to_nodes(gesca,meshe(ndivi),elmar,pp,boundary_mask=bmask)

     call memgen(2_ip,nboun,0_ip)
     call memory_deallo( mem_modul(1:2,modul),'BMASK','nsi_bougra',bmask)

  end if
  
end subroutine nsi_bougra 
