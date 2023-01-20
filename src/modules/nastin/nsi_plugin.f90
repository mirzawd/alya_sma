!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup Nastin 
!> @{
!> @file    nsi_plugin.f90
!> @date    14/10/2014
!> @author  Guillaume Houzeaux
!> @brief   Plugin for coupling
!> @details Plugin for coupling
!>          1. Allocate a minimum memory so that all ranks can enter 
!>             COU_INTERPOLATE_NODAL_VALUES without blowing up
!>             (see nsi_membcs.f90 as an example)
!> @}
!------------------------------------------------------------------------
subroutine nsi_plugin(icoup)

  use def_kintyp,                only :  ip,rp
  use def_master,                only :  momod
  use def_master,                only :  modul,mem_modul
  use def_master,                only :  vbset
  use def_domain,                only :  ndime
  use def_domain,                only :  lbsec
  use def_domain,                only :  lbset
  use def_coupli,                only :  coupling_type
  use def_coupli,                only :  UNKNOWN
  use def_coupli,                only :  RESIDUAL
  use def_coupli,                only :  SCALAR
  use def_coupli,                only :  scala_cou
  use mod_couplings,             only :  COU_INTERPOLATE_NODAL_VALUES
  use mod_couplings,             only :  cou_residual_scalar
  use def_solver,                only :  solve_sol
  use mod_matrix,                only :  matrix_initialize
  use mod_maths_arrays,          only :  maths_findloc
  use mod_memory,                only :  memory_alloca
  use mod_memory,                only :  memory_deallo
  use mod_communications_global, only :  PAR_BROADCAST,PAR_SUM
  !
  ! Possible variables => 
  ! 
  use def_master,                only : veloc
  use def_master,                only : tempe
  use def_master,                only : press,INOTMASTER
  use def_nastin,                only : bvess_nsi
  use def_nastin,                only : bpess_nsi
  use def_kermod,                only : kfl_twola_ker
  use def_nastin,                only : btrac_nsi, tracr_nsi, tluav_nsi
  use def_master,                only : momentum_sink
  use def_master,                only : ITER_K
  use def_domain,                only : npoin,nboun
  use def_nastin,                only : kfl_fixno_nsi,vafor_nsi
  use def_nastin,                only : NSI_FRACTIONAL_STEP
  use def_nastin,                only : fsifo_nsi
  use def_nastin,                only : vefix
  use def_nastin_aux,            only : bvnat_nsi
  use mod_communications,        only : PAR_MAX
  use def_coupli,                only : kfl_efect
  use mod_nsi_efect
  use mod_projec
  use mod_gradie

  implicit none
  integer(ip) :: ipoin, idime, iboun

  !
  ! <= end coupling variables
  !
  integer(ip), intent(in) :: icoup    !< Coupling number
  integer(ip)             :: iset,ibset
  character(5)            :: variable
  real(rp),    pointer    :: dumm2(:,:)
  real(rp),    pointer    :: force_nsi(:,:)
  real(rp),    pointer    :: veloc_tmp(:,:)
  real(rp)                :: press_gus

  nullify(dumm2)
  nullify(force_nsi)
  nullify(veloc_tmp)

  variable = coupling_type(icoup) % variable 
  !
  ! Velocity 
  ! 
  if( variable == 'VELOC' .or. variable == 'UNKNO' ) then  
     if (kfl_twola_ker == 0_ip) then
        !
        ! Usual velocity interpolation
        !
        call COU_INTERPOLATE_NODAL_VALUES(icoup,ndime,bvess_nsi,veloc,kfl_fixno_nsi)
     else if (kfl_efect) then
        !
        ! Interpolates velocity from solid to fluid fringe nodes in Embedded Finite Element
        ! Coupling Technique (EFECT) 
        !
        call COU_INTERPOLATE_NODAL_VALUES(icoup,ndime,vefix,dumm2)
     else                      
        !
        ! Two-layer wall model
        !
        call COU_INTERPOLATE_NODAL_VALUES(icoup,ndime,bvess_nsi,tluav_nsi,kfl_fixno_nsi)
     end if
  end if
  !
  ! Pressure
  ! 
  if( variable == 'PRESS' .or. variable == 'UNKNO' ) then   
     if (coupling_type(icoup) % what /= SCALAR  )then  !!OJO para no confundir con gusi
      call COU_INTERPOLATE_NODAL_VALUES(icoup,1_ip,bpess_nsi,press) 
     end if
  end if
  !
  ! Temperature
  !
  if( variable == 'TEMPE' ) then   
     call COU_INTERPOLATE_NODAL_VALUES(icoup,1_ip,tempe,dumm2)      
  end if
  !
  ! Momentum residual
  !
  if( variable == 'MOMEN' .or. variable == 'RESID' ) then
     !
     if( kfl_efect ) then
        !
        ! Embedded Finite Element Coupling Technique 
        !
        call memory_alloca(mem_modul(1:2,modul),'FORCE_NSI','nsi_plugin',force_nsi,ndime,max(1_ip,npoin))
        call nsi_algebraic_reaction_force(icoup, force_nsi) ! TODO: Check that I'm not using old LHS/RHS to compute forces
        if( INOTMASTER ) fsifo_nsi(1:ndime,1:npoin) = force_nsi(1:ndime,1:npoin)
        call COU_INTERPOLATE_NODAL_VALUES(icoup,ndime,dumm2,force_nsi)
     else
        !
        ! Fractional-step ALE FSI
        !
        call matrix_initialize(momod(modul) % solve(1) % block_array(1) % bvnat)
        if( NSI_FRACTIONAL_STEP ) then
           if( INOTMASTER ) then
              do ipoin = 1,npoin
                 solve_sol(1) % reaction(1:ndime,ipoin) = vafor_nsi(1:ndime,ipoin)
              end do
           end if
        end if
        call COU_INTERPOLATE_NODAL_VALUES(icoup,ndime,momod(modul) % solve(1) % block_array(1) % bvnat,momod(modul) % solve(1) % reaction)
     end if
     !
  end if
  !
  ! Continuity residual
  !
  if( variable == 'CONTI' .or. variable == 'RESID' ) then
     call matrix_initialize(momod(modul) % solve(1) % block_array(2) % bvnat)
     call COU_INTERPOLATE_NODAL_VALUES(icoup,1_ip,momod(modul) % solve(1) % block_array(2) % bvnat,momod(modul) % solve(2) % reaction)
  end if
  !
  ! Traction for two-layer wall modelling (RANS/LES coupling)
  !
  if( variable == 'TRACT' .or. variable == 'UNKNO' ) then
     call COU_INTERPOLATE_NODAL_VALUES(icoup,ndime,btrac_nsi,tracr_nsi,kfl_fixno_nsi)
  end if
  !
  ! Momentum sink
  ! 
  if( variable == 'MOMSK' ) then
     if( .not. associated(momentum_sink) ) then 
        call memory_alloca(mem_modul(1:2,modul),'MOMENTUM_SINK','nsi_plugin',momentum_sink,ndime,max(1_ip,npoin))
     endif
     call COU_INTERPOLATE_NODAL_VALUES(icoup,ndime,momentum_sink)
  end if
  !
  ! Immersed boundary: solid node advection
  !
  if( variable == 'ADVEC' ) then
     if( .not. associated(veloc_tmp) ) then 
        call memory_alloca(mem_modul(1:2,modul),'VELOC_TMP','nsi_plugin',veloc_tmp,ndime,max(1_ip,npoin))
     endif
     do ipoin = 1,npoin
        do idime = 1,ndime
           veloc_tmp(idime,ipoin) = veloc(idime,ipoin,ITER_K)
        end do
     end do
     call COU_INTERPOLATE_NODAL_VALUES(icoup,ndime,dumm2,veloc_tmp)
  end if
  !
  ! Immersed boundary: Residual momentum sink
  !
  if ( variable == 'REMSK' ) then
     !
     ! Extrapolate solid body force and indicator field to fluid
     !
     do ipoin = 1,npoin
        do idime = 1,ndime
           fsifo_nsi(idime,ipoin) = 0.0_rp
        end do
     end do
     call COU_INTERPOLATE_NODAL_VALUES(icoup,ndime,fsifo_nsi,dumm2)
     !
     ! Copy FSI forces to BVNAT
     !
     if ( INOTMASTER ) then
        do ipoin = 1,npoin
           do idime = 1,ndime
              solve_sol(1) % bvnat(idime,ipoin) = fsifo_nsi(idime,ipoin)
           end do
        end do
     end if
     !
  end if
  !
  ! Flow rate 
  ! 
  if( variable == 'FLOWR' .and. coupling_type(icoup) % what == SCALAR ) then
     iset  = coupling_type(icoup) % where_number_source
     ibset = maths_findloc(lbsec,iset)

     call nsi_bouset(lbsec(ibset),ibset)
     call PAR_SUM(vbset(2,ibset))

     call PAR_BROADCAST(vbset(2,ibset),'IN THE WORLD')
     scala_cou(icoup,1) = vbset(2,ibset)
     call cou_residual_scalar(coupling_type(icoup))
  end if
  !
  ! Pressure
  ! 
  if( variable == 'PRESS' .and. coupling_type(icoup) % what == SCALAR ) then
     press_gus = -huge(1.0_rp)
     iset  = coupling_type(icoup) % where_number_source
     ibset = maths_findloc(lbsec,iset)
     call PAR_MAX(press_gus,'IN THE WORLD')
     !if(kfl_paral==0) print*,'NSI RECV PRESS: ',icoup,press_gus
     do iboun = 1,nboun
        if(ibset == lbset(iboun))then
           bvnat_nsi(1,iboun,1) = press_gus
        end if
     end do
  end if

  if(associated(dumm2))      deallocate(dumm2)
  if(associated(force_nsi))  call memory_deallo(mem_modul(1:2,modul),'FORCE_NSI','nsi_plugin',force_nsi)
  if(associated(veloc_tmp))  call memory_deallo(mem_modul(1:2,modul),'VELOC_TMP','nsi_plugin',veloc_tmp)

end subroutine nsi_plugin
