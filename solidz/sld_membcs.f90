!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine sld_membcs(itask)
!-----------------------------------------------------------------------
!****f* Solidz/sld_membcs
! NAME
!    sld_memcbs
! DESCRIPTION
!    This routine allocates memory for the boundary conditions arrays
! USES
!    ecoute
! USED BY
!    sld_reabcs
!***
!-----------------------------------------------------------------------
  use def_parame
  use def_inpout
  use def_master
  use def_domain
  use mod_memchk
  use def_solidz
  use mod_memory,    only : memory_alloca
  use mod_memory,    only : memory_alloca_min
  use def_coupli,    only : kfl_immer

  implicit none
  integer(ip), intent(in) :: itask
  integer(4)              :: istat
  integer(ip)             :: ifunc,dummi

  if (itask == 1_ip) then
     !
     ! Mandatory arrays
     !
     call memory_alloca( mem_modul(1:2,modul),'KFL_FIXNO_SLD','sld_membcs',kfl_fixno_sld,ndime,npoin     ) ! Fixity on nodes
     call memory_alloca( mem_modul(1:2,modul),'KFL_FIXBO_SLD','sld_membcs',kfl_fixbo_sld,nboun           ) ! Fixity on boundaries
     call memory_alloca( mem_modul(1:2,modul),'KFL_FIXRS_SLD','sld_membcs',kfl_fixrs_sld,npoin           ) ! Local axes

     if( kfl_conbc_sld == 1_ip ) then ! constant bc
        call memory_alloca( mem_modul(1:2,modul),'BVESS_SLD'    ,'sld_membcs',bvess_sld,ndime,npoin,1_ip ) ! Prescribed displacement
        call memory_alloca( mem_modul(1:2,modul),'BVNAT_SLD'    ,'sld_membcs',bvnat_sld,ndime,nboun,1_ip ) ! Natural condition values
     else
        call memory_alloca( mem_modul(1:2,modul),'KFL_FUNNO_SLD','sld_membcs',kfl_funno_sld,npoin        )
        call memory_alloca( mem_modul(1:2,modul),'KFL_FUNBO_SLD','sld_membcs',kfl_funbo_sld,nboun        )
        call memory_alloca( mem_modul(1:2,modul),'KFL_FUNTN_SLD','sld_membcs',kfl_funtn_sld,npoin        )
        call memory_alloca( mem_modul(1:2,modul),'KFL_FUNTB_SLD','sld_membcs',kfl_funtb_sld,nboun        )
        call memory_alloca( mem_modul(1:2,modul),'BVESS_SLD'    ,'sld_membcs',bvess_sld,ndime,npoin,2_ip )
        call memory_alloca( mem_modul(1:2,modul),'BVNAT_SLD'    ,'sld_membcs',bvnat_sld,ndime,nboun,3_ip )
     end if
     !
     ! Rotation matrix for local axes / PDN-contact
     !
     call memory_alloca( mem_modul(1:2,modul),'JACROT_DU_DQ_SLD','sld_membcs',jacrot_du_dq_sld,ndime,ndime,npoin ) ! Rot. matrix
     call memory_alloca( mem_modul(1:2,modul),'JACROT_DQ_DU_SLD','sld_membcs',jacrot_dq_du_sld,ndime,ndime,npoin ) ! Rot. matrix (transpose)
     !
     ! Stent arrays
     !
     if (kfl_conta_stent > 0_ip) then 
        call memory_alloca( mem_modul(1:2,modul),'kfl_contn_stent','sld_membcs',kfl_contn_stent, npoin        )
     else  
        call memory_alloca( mem_modul(1:2,modul),'kfl_contn_stent','sld_membcs',kfl_contn_stent, npoin       ) ! Prescribed displacement
     end if
     !
     if (kfl_immer) then
        call memory_alloca( mem_modul(1:2,modul),'KFL_IMMER_SLD','sld_membcs',kfl_immer_sld,ndime,npoin     ) ! Fixity on nodes for immersed boundary
     end if

  else if (itask == 3_ip) then

     allocate(tload_sld(10),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'TLOAD_SLD','sld_reabcs',tload_sld)

  else if ((itask > 10_ip).and.(itask < 21_ip)) then
     ifunc = itask - 10
     allocate(tload_sld(ifunc)%a(ndime+1,mtloa_sld(ifunc)),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'TLOAD_SLD','sld_reabcs',tload_sld(ifunc)%a)

  else if( itask ==  30 ) then

     allocate(crkpo_sld(ndime,ncrak_sld),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'CRKPO_SLD','sld_reabcs',crkpo_sld)
     allocate(crkno_sld(ndime,ncrak_sld),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'CRKNO_SLD','sld_reabcs',crkno_sld)

  else if( itask == -30 ) then

     call memchk(two,istat,mem_modul(1:2,modul),'CRKPO_SLD','sld_membcs',crkpo_sld)
     deallocate(crkpo_sld,stat=istat)
     if( istat /= 0 ) call memerr(two,'CRKPO_SLD','sld_membcs',0_ip)
     call memchk(two,istat,mem_modul(1:2,modul),'CRKNO_SLD','sld_membcs',crkno_sld)
     deallocate(crkno_sld,stat=istat)
     if( istat /= 0 ) call memerr(two,'CRKNO_SLD','sld_membcs',0_ip)

  else if( itask ==  31 ) then

     if( ndime == 2 ) then
        dummi = 2
     else
        dummi = 4
     end if
     allocate(crkco_sld(ndime,dummi,ncrak_sld),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'CRKCO_SLD','sld_reabcs',crkco_sld)

  else if( itask == -31 ) then

     call memchk(two,istat,mem_modul(1:2,modul),'CRKCO_SLD','sld_membcs',crkco_sld)
     deallocate(crkco_sld,stat=istat)
     if( istat /= 0 ) call memerr(two,'CRKCO_SLD','sld_membcs',0_ip)

  end if

end subroutine sld_membcs

