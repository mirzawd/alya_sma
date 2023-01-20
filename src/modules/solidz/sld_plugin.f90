!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!----------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    sld_plugin.f90
!> @author  J.C. Cajas
!> @date    17/04/2014
!> @brief   Receive force from Nastin and send displacement to Alefor
!> @details Receive force from Nastin and send displacement to Alefor
!> @        usind the coupling structures and functions.
!> @} 
!----------------------------------------------------------------------
subroutine sld_plugin(icoup)
  !
  ! Obligatory variables 
  !
  use def_coupli,         only :  coupling_type
  use def_domain,         only :  npoin, ndime
  use def_master,         only :  solve_sol,modul,mem_modul
  use def_kintyp,         only :  ip,rp,lg
  use def_master,         only :  INOTMASTER
  use def_master,         only :  TIME_N, ITER_K
  use mod_couplings,      only :  COU_INTERPOLATE_NODAL_VALUES
  use mod_couplings,      only :  COU_SET_FIXITY_ON_TARGET
  use mod_memory,         only :  memory_deallo
  use mod_memory,         only :  memory_alloca,memory_alloca_min
  use mod_matrix,         only :  matrix_initialize
  use mod_bouder
  !
  ! Possible variables 
  !
  use def_master,        only :  displ
  use def_master,        only :  ID_NASTIN
  use def_master,        only :  gdepo
  use def_master,        only :  dtime 
  use def_solidz,        only :  bvess_sld
  use def_solidz,        only :  kfl_gdepo
  use def_solidz,        only :  kfl_fixno_sld
  use def_solidz,        only :  kfl_minco_sld
  use def_solidz,        only :  kfl_immer_sld
  use def_solidz,        only :  veloc_sld
  use def_solidz,        only :  fintt_sld
  use def_solidz,        only :  fextt_sld
  use mod_eccoupling,    only :  eccou_do_plugin
  use def_coupli,        only :  kfl_efect 
  use mod_sld_efect

  implicit none 
  real(rp),    pointer    :: xvalu(:,:)
  real(rp),    pointer    :: svalu(:,:)
  real(rp)                :: foref  ! for the coupling with nastin
  real(rp),    pointer    :: force_sld(:,:) 
  real(rp),    pointer    :: stres(:,:)
  real(rp),    pointer    :: tract(:,:) 
  real(rp),    pointer    :: rdata(:,:)

  integer(ip)             :: idime, jdime
  integer(ip)             :: ipoin, kpoin
  integer(ip), intent(in) :: icoup
  character(5)            :: variable
  integer(ip), save       :: ipass = 0_ip

  nullify(xvalu) 
  nullify(svalu)
  nullify(force_sld)
  nullify(stres)
  nullify(tract)
  nullify(rdata)

  variable = coupling_type(icoup) % variable

  if( variable == 'DISPL' ) then   
     !
     ! Displacement
     !
     if( ipass == 0_ip ) then
        ipass = 1_ip
        !call COU_SET_FIXITY_ON_TARGET('DISPL',modul,kfl_fixno_sld)
     end if

     call memory_alloca(mem_modul(1:2,modul),'SVALU','sld_plugin',svalu,ndime,max(1_ip,npoin))
     call memory_alloca(mem_modul(1:2,modul),'XVALU','sld_plugin',xvalu,ndime,max(1_ip,npoin))
     if( INOTMASTER ) then
        svalu = displ(:,:,1)
     else
        svalu = 0.0_rp
     endif
    
     call COU_INTERPOLATE_NODAL_VALUES(icoup,ndime,xvalu,svalu)

  else if( variable == 'VELOC' ) then
     !
     ! Velocity, used in Embedded Finite Element Coupling Technique (EFECT)
     !
     call memory_alloca(mem_modul(1:2,modul),'SVALU','sld_plugin',svalu,ndime,max(1_ip,npoin))
     call memory_alloca(mem_modul(1:2,modul),'XVALU','sld_plugin',xvalu,ndime,max(1_ip,npoin))
     if( INOTMASTER ) then
        svalu = displ(:,:,ITER_K)
     else
        svalu = 0.0_rp
     endif
     !
     ! Interpolate solid velocities to fluid domain
     !
     call COU_INTERPOLATE_NODAL_VALUES(icoup,ndime,xvalu,svalu)

  elseif( variable == 'ECCOU' )then

     call eccou_do_plugin(icoup)

  else if( variable == 'RESID' .or. variable == 'MOMEN' .or. variable == 'FWALL' ) then
     !
     ! Residual
     !
     if( kfl_efect ) then
        !
        ! Embedded Finite Element Coupling Techinque (EFECT) for FSI
        !
        if(kfl_minco_sld .eq. 0_ip) kfl_minco_sld = icoup
        if (icoup .eq. kfl_minco_sld) call matrix_initialize(solve_sol(1) % bvnat)
        call memory_alloca(mem_modul(1:2,modul),'TRACT','sld_plugin',tract,ndime,max(1_ip,npoin))
        ! Set fixity on solid nodes
        if( ipass == 0_ip ) then
           ipass = 1_ip
           call COU_SET_FIXITY_ON_TARGET('RESID',modul,kfl_fixno_sld,'FREE FIXITY')
        end if
        ! Interpolate fluid algebraic reaction forces to solid boundary
        call COU_INTERPOLATE_NODAL_VALUES(icoup,ndime,tract,svalu)
        ! Copy tractions to bvnat
        if( INOTMASTER ) then
           do ipoin = 1,npoin
              do idime = 1,ndime
                 if( kfl_fixno_sld(idime,ipoin) /= 1 ) then
                    solve_sol(1) % bvnat(idime,ipoin) = tract(idime,ipoin)
                 end if
              end do
           end do
        end if
     else if( kfl_gdepo == 1_ip ) then
        !
        ! Push forward for the force coming from nastin
        !
        if(kfl_minco_sld.eq.0_ip) kfl_minco_sld=icoup
        if (icoup .eq. kfl_minco_sld) call matrix_initialize(solve_sol(1) % bvnat)
        call COU_INTERPOLATE_NODAL_VALUES(icoup,ndime,solve_sol(1) % bvnat,solve_sol(1) % reaction)
        if( coupling_type(icoup) % module_source == ID_NASTIN )then
           if( INOTMASTER ) then
              do ipoin = 1,npoin
                 do idime = 1,ndime
                    foref = 0.0_rp
                    if( kfl_fixno_sld(idime,ipoin) /= 1 ) then
                       do jdime = 1,ndime
                          foref = foref + gdepo(idime,jdime,ipoin) * solve_sol(1) % bvnat(jdime,ipoin) 
                       end do
                       solve_sol(1) % bvnat(idime,ipoin) = foref
                    end if
                 end do
              end do
           end if
        end if
     else
        call matrix_initialize(solve_sol(1) % bvnat)
        call COU_INTERPOLATE_NODAL_VALUES(icoup,ndime,solve_sol(1) % bvnat,solve_sol(1) % reaction)
        if( coupling_type(icoup) % module_source == ID_NASTIN )then
           if( INOTMASTER ) then
              do kpoin = 1, coupling_type(icoup) % wet % npoin_wet
                 ipoin = coupling_type(icoup) % wet % lpoin_wet(kpoin)
                 do idime = 1,ndime
                    if( kfl_fixno_sld(idime,ipoin) == 1_ip ) then
                       solve_sol(1) % bvnat(idime,ipoin) = 0.0_rp
                    end if
                 end do
              end do
           end if
        end if
     end if


  else if( variable == 'ALEFO' ) then
     !
     ! Coupling with alefor
     !
     call memory_alloca(mem_modul(1:2,modul),'SVALU','sld_plugin',svalu,ndime,max(1_ip,npoin))
     svalu = 0.0_rp
     
     if( INOTMASTER ) then
        if( coupling_type(icoup) % frequ_send > 1_ip .or. coupling_type(icoup) % frequ_recv > 1_ip ) then
           !
           ! Subcycling
           !
           do kpoin = 1,coupling_type(icoup) % geome % npoin_source
              ipoin = coupling_type(icoup) % geome % lpoin_source(kpoin)
              svalu(1:ndime,ipoin) = displ(1:ndime,ipoin,1_ip) - coupling_type(icoup) % values_frequ(1:ndime,kpoin,1_ip)
              coupling_type(icoup) % values_frequ(1:ndime,kpoin,2_ip) = displ(1:ndime,ipoin,1_ip)
           end do
        else
           !
           ! If the exchanges are on each time step, the displacement in calculated with displ only
           !
           do kpoin = 1,coupling_type(icoup) % geome % npoin_source
              ipoin = coupling_type(icoup) % geome % lpoin_source(kpoin)
              svalu(1:ndime,ipoin) = displ(1:ndime,ipoin,1_ip) - displ(1:ndime,ipoin,3_ip)
           end do
        end if
     end if

     call COU_INTERPOLATE_NODAL_VALUES(icoup,ndime,bvess_sld,svalu)

  else if( variable == 'ADVEC' ) then
     !
     ! Immersed boundary: advection of solid nodes
     !
     call memory_alloca(mem_modul(1:2,modul),'XVALU','sld_plugin',xvalu,ndime,max(1_ip,npoin))
     !
     ! Set fixity on solid nodes
     !
     if( ipass == 0_ip ) then
        ipass = 1_ip
        ! Save fixities of target before coupling modifies it
        do ipoin = 1,npoin
           do idime = 1,ndime
              if (kfl_fixno_sld(idime,ipoin) == 0_ip) then
                 kfl_immer_sld(idime,ipoin) = 1_ip
              else
                 kfl_immer_sld(idime,ipoin) = 0_ip
              end if
           end do
        end do
        ! Modify fixity of target
        call COU_SET_FIXITY_ON_TARGET('ADVEC',modul,kfl_fixno_sld)
     end if
     !
     ! Interpolate velocity from NASTIN to SOLID boundary
     !
     call COU_INTERPOLATE_NODAL_VALUES(icoup,ndime,xvalu,svalu)
     !
     ! Integrate displacements from interpolated velocities
     !
     if ( INOTMASTER ) then
        do ipoin=1,npoin
           do idime=1,ndime
              !if ( kfl_immer_sld(idime,ipoin) == 1_ip) then
                 !
                 ! First order Forward Euler
                 !
                 veloc_sld(idime,ipoin,ITER_K) = xvalu(idime,ipoin)
                 bvess_sld(idime,ipoin,ITER_K) = displ(idime,ipoin,TIME_N) + dtime*xvalu(idime,ipoin)
              !end if
           end do
        end do
     end if

  else if( variable == 'REMSK' ) then
     !
     ! Immersed boundary: Residual momentum sink
     !
     call memory_alloca(mem_modul(1:2,modul),'SVALU','sld_plugin',xvalu,ndime,max(1_ip,npoin))
     !
     ! TODO: Ignore forces from real dirichlet nodes
     !
     do ipoin = 1,npoin
        do idime = 1,ndime
           xvalu(idime,ipoin) = fextt_sld(idime,ipoin,ITER_K) - fintt_sld(idime,ipoin,ITER_K)
        end do
     end do
     !
     ! Spread solid forces to fluid
     !
     call COU_INTERPOLATE_NODAL_VALUES(icoup,ndime,svalu,xvalu)

  end if

  if( associated(xvalu) )     call memory_deallo(mem_modul(1:2,modul),'XVALU','sld_plugin',xvalu)
  if( associated(svalu) )     call memory_deallo(mem_modul(1:2,modul),'SVALU','sld_plugin',svalu)
  if( associated(stres) )     call memory_deallo(mem_modul(1:2,modul),'STRES','sld_plugin',stres)
  if( associated(tract) )     call memory_deallo(mem_modul(1:2,modul),'TRACT','sld_plugin',tract)
  if( associated(rdata) )     call memory_deallo(mem_modul(1:2,modul),'RDATA','sld_plugin',rdata)
  if( associated(force_sld) ) call memory_deallo(mem_modul(1:2,modul),'FORCE_SLD','sld_plugin',force_sld)

end subroutine sld_plugin
