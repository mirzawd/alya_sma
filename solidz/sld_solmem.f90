!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    sld_solmem.f90
!> @date    14/06/2019
!> @author  hozeaux
!> @brief   Allocate memory
!> @details Allocate memory for solver and variables that should be
!>          reallocated and not redistributed in case of repartitioning
!> @}
!-----------------------------------------------------------------------

subroutine sld_solmem()

  use def_kintyp,               only : ip
  use def_master,               only : INOTMASTER
  use def_master,               only : ID_EXMEDI
  use def_master,               only : TIME_N
  use def_master,               only : solve
  use def_master,               only : press
  use def_master,               only : velom
  use def_master,               only : modul, mem_modul
  use def_master,               only : gdepo,gdeinv,gdedet
  use def_master,               only : coupling
  use def_master,               only : kfl_modul
  use def_master,               only : donna_gp
  use def_domain,               only : ndime,npoin,nelem,ngaus
  use def_domain,               only : ltype
  use def_solver,               only : solve_sol
  use mod_memory,               only : memory_alloca
  use mod_biofibers,            only : kfl_biofibers
  use mod_output_postprocess,   only : output_postprocess_check_variable_postprocess
  use mod_output_postprocess,   only : output_postprocess_check_variable_witness
  use mod_output_postprocess,   only : output_postprocess_check_variable_node_sets
  use mod_output_postprocess,   only : output_postprocess_check_variable_element_sets
  use mod_eccoupling,           only : kfl_exmsld_ecc
  use mod_sld_rbo,              only : sld_rbo_solmem
  use mod_sld_atm,              only : kfl_therm_sld
  use mod_sld_atm,              only : sld_atm_allocate_memmory
  use def_solidz

  implicit none

  integer(ip) :: ielem,pgaus
  integer(ip) :: istat

  if( kfl_rigid_sld == 0 ) then
     
     !----------------------------------------------------------------------
     !
     ! Solver
     !
     !----------------------------------------------------------------------
     !
     ! Memory
     !
     solve_sol => solve(1:2)
     call soldef(4_ip)
     !
     ! Boundary conditions
     !
     solve_sol => solve(1:)
     solve(1) % bvess     => bvess_sld(:,:,1)
     solve(1) % bvnat     => bvnat_sld(:,:,1)
     solve(1) % kfl_fixno => kfl_fixno_sld
     !
     ! Dirichlet value is 0 (we solve increments)
     !
     if( solve(1) % kfl_iffix /= 0 ) solve(1) % kfl_iffix = 2

     !----------------------------------------------------------------------
     !
     ! Global variables
     !
     !----------------------------------------------------------------------
     !
     ! Force vectors
     !
     call memory_alloca(mem_modul(1:2,modul),'FINTE_SLD','sld_solmem',finte_sld,ndime*npoin)
     call memory_alloca(mem_modul(1:2,modul),'FINTT_SLD','sld_solmem',fintt_sld,ndime,npoin,ncomp_sld)
     call memory_alloca(mem_modul(1:2,modul),'FEXTE_SLD','sld_solmem',fexte_sld,ndime*npoin)
     call memory_alloca(mem_modul(1:2,modul),'FEXTT_SLD','sld_solmem',fextt_sld,ndime,npoin,ncomp_sld)
     call memory_alloca(mem_modul(1:2,modul),'MACCE_SLD','sld_solmem',macce_sld,ndime*npoin)
     call memory_alloca(mem_modul(1:2,modul),'FRXID_SLD','sld_solmem',frxid_sld,ndime*npoin)
     !
     ! Mass matrix
     !
     call memory_alloca(mem_modul(1:2,modul),'VMASS_SLD',   'sld_solmem',vmass_sld,npoin)
     !
     ! Variables solution methods
     !
     call memory_alloca(mem_modul(1:2,modul),'DUNKN_SLD',   'sld_solmem',dunkn_sld,   ndime*npoin)
     call memory_alloca(mem_modul(1:2,modul),'UNKNOTMP_SLD','sld_solmem',unknotmp_sld,ndime*npoin)
     call memory_alloca(mem_modul(1:2,modul),'VELOCTMP_SLD','sld_solmem',veloctmp_sld,ndime,npoin)
     call memory_alloca(mem_modul(1:2,modul),'DDISP_SLD',   'sld_solmem',ddisp_sld,   ndime,npoin,ncomp_sld)
     !
     ! Element characteristic length
     !
     call memory_alloca(mem_modul(1:2,modul),'CELEN_SLD','sld_solmem',celen_sld,nelem)
     if( kfl_damag_sld == 1 ) then
        call memory_alloca(mem_modul(1:2,modul),'SRPRO_SLD','sld_solmem',srpro_sld,8_ip,nelem)
     end if
     !
     ! Global gradient deformation operator at nodes
     !
     if( kfl_gdepo /= 0 ) then
        if( .not. kfl_biofibers )then
           call memory_alloca(mem_modul(1:2,modul),'GDEPO',  'sld_solmem',gdepo,  ndime,ndime,npoin)
        endif
        call memory_alloca(mem_modul(1:2,modul),'GDEINV', 'sld_solmem',gdeinv, ndime,ndime,npoin)
        call memory_alloca(mem_modul(1:2,modul),'GDEDET', 'sld_solmem',gdedet, npoin)
     end if
     !
     ! Energies
     !
     call memory_alloca(mem_modul(1:2,modul),'ALLWK_SLD','sld_solmem',allwk_sld,ncomp_sld)
     call memory_alloca(mem_modul(1:2,modul),'ALLIE_SLD','sld_solmem',allie_sld,ncomp_sld)
     call memory_alloca(mem_modul(1:2,modul),'ALLKE_SLD','sld_solmem',allke_sld,ncomp_sld)
     call memory_alloca(mem_modul(1:2,modul),'ETOTA_SLD','sld_solmem',etota_sld,ncomp_sld)
     !
     ! Material axes (element level)
     !
     if( kfl_fiber_sld > 4_ip ) then
        !
        ! Material axes
        !
        call memory_alloca(mem_modul(1:2,modul),'AXIS1_SLD','sld_memphy',axis1_sld,ndime,nelem)
        call memory_alloca(mem_modul(1:2,modul),'AXIS2_SLD','sld_memphy',axis2_sld,ndime,nelem)
        if( ndime == 3_ip ) then
           call memory_alloca(mem_modul(1:2,modul),'AXIS3_SLD','sld_memphy',axis3_sld,ndime,nelem)
        end if
     end if
     
     !----------------------------------------------------------------------
     !
     ! Post-process
     !
     !----------------------------------------------------------------------
     !
     ! Cauchy stress
     !
     call memory_alloca(mem_modul(1:2,modul),'CAUST_SLD','sld_solmem',caust_sld,nvoig_sld,npoin)
     !
     ! Green streen
     !
     call memory_alloca(mem_modul(1:2,modul),'GREEN_SLD','sld_solmem',green_sld,nvoig_sld,npoin)
     call memory_alloca(mem_modul(1:2,modul),'EPSEL_SLD','sld_solmem',epsel_sld,nelem)
     do ielem = 1,nelem
        pgaus = ngaus(abs(ltype(ielem)))
        call memory_alloca(mem_modul(1:2,modul),'EPSEL_SLD % A','sld_solmem',epsel_sld(ielem)%a,ndime,ndime,pgaus)
     end do
     !
     ! Von Mises stress
     !
     if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='SEQVM') .or. &
         output_postprocess_check_variable_witness(    VARIABLE_NAME='SEQVM') .or. &
         output_postprocess_check_variable_node_sets(  VARIABLE_NAME='SEQVM') ) then
        call memory_alloca(mem_modul(1:2,modul),'SEQVM_SLD','sld_solmem',seqvm_sld,npoin)
     end if
     !
     ! Logarithmic strain
     !
     if( output_postprocess_check_variable_postprocess( VARIABLE_NAME='LNEPS') .or. &
         output_postprocess_check_variable_witness(     VARIABLE_NAME='LEPSI',NCHAR=3_ip) .or. &
         output_postprocess_check_variable_node_sets(   VARIABLE_NAME='LEPSI',NCHAR=3_ip) .or. &
         output_postprocess_check_variable_element_sets(VARIABLE_NAME='LEPRT') ) then
        call memory_alloca(mem_modul(1:2,modul),'LEPSI_SLD','sld_solmem',lepsi_sld,nvoig_sld,npoin)
        call memory_alloca(mem_modul(1:2,modul),'LEPSE_SLD','sld_solmem',lepse_sld,nelem)
        do ielem = 1,nelem
           pgaus = ngaus(abs(ltype(ielem)))
           call memory_alloca(mem_modul(1:2,modul),'LEPSE_SLD % A','sld_solmem',lepse_sld(ielem)%a,ndime,ndime,pgaus)
        end do
     end if
     !
     ! Strain for Bio fibers in longitudinal direction
     !
     if( output_postprocess_check_variable_element_sets(VARIABLE_NAME='EBFIL') ) then
        call memory_alloca(mem_modul(1:2,modul),'EBFIL_SLD','sld_solmem',ebfil_sld,nelem)
     end if
     !
     ! Stress for Bio fibers in longitudinal direction
     !
     if( output_postprocess_check_variable_element_sets(VARIABLE_NAME='SBFIL') ) then
        call memory_alloca(mem_modul(1:2,modul),'EBFIL_SLD','sld_solmem',sbfil_sld,nelem)
     end if
     !
     ! Material orientations (element level)
     !
     if( kfl_fiber_sld > 3_ip ) then
        !
        ! Rotation matrix for material axes
        !
        call memory_alloca(mem_modul(1:2,modul),'RMATE_SLD','sld_memphy',rmate_sld,nelem)
        do ielem = 1,nelem
           pgaus = ngaus(abs(ltype(ielem)))
           call memory_alloca(mem_modul(1:2,modul),'RMATE_SLD % A','sld_memphy',rmate_sld(ielem)%a,ndime,ndime,pgaus)
        end do
     end if
     !
     !Elastic strain energy
     !
     call memory_alloca(mem_modul(1:2,modul),'ENDEN_SLD','sld_solmem',enden_sld,nelem)
     do ielem = 1,nelem
        pgaus = ngaus(abs(ltype(ielem)))
        call memory_alloca(mem_modul(1:2,modul),'ENDEN_SLD % A','sld_memphy',enden_sld(ielem)%a, pgaus)
     end do
     !
     !Water content for sm400
     !
     call memory_alloca(mem_modul(1:2,modul),'WATER_SLD','sld_solmem',water_sld,nelem)
     do ielem = 1,nelem
        pgaus = ngaus(abs(ltype(ielem)))
        call memory_alloca(mem_modul(1:2,modul),'WATER_SLD % A','sld_memphy',water_sld(ielem)%a, pgaus)
     end do
     !
     !Donnan osmosis for sm400
     !
     call memory_alloca(mem_modul(1:2,modul),'DONNA_GP','sld_solmem',donna_gp,nelem)
     do ielem = 1,nelem
        pgaus = ngaus(abs(ltype(ielem)))
        call memory_alloca(mem_modul(1:2,modul),'DONNA_GP % A','sld_memphy',donna_gp(ielem)%a, pgaus)
     end do

     !----------------------------------------------------------------------
     !
     ! Couplings
     !
     !----------------------------------------------------------------------

     call memory_alloca(mem_modul(1:2,modul),'PRESS',   'sld_solmem',press,npoin,ncomp_sld)
     !
     ! Coupling with NASTIN
     !
     if( coupling('SOLIDZ','NASTIN') >= 1 ) then
        ! the mesh velocity is the solid deformation velocity
        ! velom => veloc_sld(:,:,1)
        call memory_alloca(mem_modul(1:2,modul),'DDISM_SLD','sld_solmem',ddism_sld,ndime,npoin)
     end if
     !
     ! Coupling with EXMEDI
     !
     if( kfl_exmsld_ecc )then
        if( kfl_gdepo == 0_ip )then
           call runend('SLD_SOLMEM: ROTATION ON SHOULD BE INCLUDED IN SOLIDZ-EXMEDI SIMULATIONS ')
        endif
     else
        if( kfl_modul(ID_EXMEDI) == 1_ip )then
           call runend('SLD_SOLMEM: COUPLING WITH EXMEDI NOW IN THE KER.DAT FILE!! CHECK DOCS.')
        endif
     end if
     !
     ! Coupling with TEMPER
     !
     if( kfl_therm_sld ) then
        call sld_atm_allocate_memmory()
     end if

     !----------------------------------------------------------------------
     !
     ! Others: (Requires revision/deletion)
     !
     !----------------------------------------------------------------------

     if( output_postprocess_check_variable_witness(     VARIABLE_NAME='SIGFI') ) then 
        call memory_alloca(mem_modul(1:2,modul),'FIBDE_SLD','sld_solmem',fibde_sld,ndime,npoin)
     end if
     if( output_postprocess_check_variable_element_sets(VARIABLE_NAME='GRLST',NCHAR=3_ip) ) then
        call memory_alloca(mem_modul(1:2,modul),'GRLST_SLD','sld_solmem',grlst_sld,nelem,ndime*ndime)
     end if
     if( output_postprocess_check_variable_postprocess( VARIABLE_NAME='EPRIN') .or. &
         output_postprocess_check_variable_node_sets(   VARIABLE_NAME='EPRIN') .or. &
         output_postprocess_check_variable_postprocess( VARIABLE_NAME='SIGEI') ) then
        call memory_alloca(mem_modul(1:2,modul),'EPRIN_SLD','sld_solmem',eprin_sld,npoin)
        call memory_alloca(mem_modul(1:2,modul),'SIGEI_SLD','sld_solmem',sigei_sld,npoin)
        call memory_alloca(mem_modul(1:2,modul),'ROLOC_SLD','sld_solmem',roloc_sld,ndime,ndime,npoin)
     end if
     if( output_postprocess_check_variable_postprocess( VARIABLE_NAME='NSIGN') ) then
        call memory_alloca(mem_modul(1:2,modul),'CAUNN_SLD','sld_solmem',caunn_sld,npoin)
     end if
     if( INOTMASTER ) then
        call memory_alloca(mem_modul(1:2,modul),'GPGDI_SLD','sld_solmem',gpgdi_sld,nelem)
        call memory_alloca(mem_modul(1:2,modul),'GPPIO_SLD','sld_solmem',gppio_sld,nelem)
        call memory_alloca(mem_modul(1:2,modul),'DEDEF_SLD','sld_solmem',dedef_sld,nelem)
        do ielem = 1,nelem
           pgaus = ngaus(abs(ltype(ielem)))
           call memory_alloca(mem_modul(1:2,modul),'GPGDI_SLD % A','sld_solmem',gpgdi_sld(ielem)%a,ndime,ndime,pgaus)
           call memory_alloca(mem_modul(1:2,modul),'GPPIO_SLD % A','sld_solmem',gppio_sld(ielem)%a,ndime,ndime,pgaus)
           call memory_alloca(mem_modul(1:2,modul),'DEDEF_SLD % A','sld_solmem',dedef_sld(ielem)%a,pgaus)
        end do
     endif
     if( kfl_plast_sld == 1 ) then
        allocate(gpsl0_sld(ndime,ndime,30),stat=istat); gpsl0_sld=0.0_rp
        allocate(rstr0_sld(ndime,ndime,30),stat=istat); rstr0_sld=0.0_rp
     end if
     !
     ! Assigns the last veloc_sld to velom. It will be used by exmedi when coupled with exm, to compute the ALE mesh velocity terms
     !     
     velom => veloc_sld(:,:,TIME_N)
     
  else

     !-------------------------------------------------------------------
     !
     ! Rigid body
     !
     !-------------------------------------------------------------------

     call sld_rbo_solmem()

  end if

end subroutine sld_solmem
