!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine ker_memall()
  !-----------------------------------------------------------------------
  !****f* Kermod/ker_memall
  ! NAME 
  !    ker_memall
  ! DESCRIPTION
  !    This routine allocates memory for the arrays needed to solve the
  !    NS equations
  ! USES
  !    memchk
  !    mediso
  ! USED BY
  !    ker_turnon
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_inpout
  use def_master
  use def_kermod
  use def_domain
  use def_solver
  use mod_memory
  use mod_maths,          only : maths_equalize_arrays
  use mod_chktyp,         only : check_type
  use mod_communications, only : PAR_SUM
  use mod_ker_subdomain,  only : ker_subdomain_motion_exists
  use mod_ker_arrays,     only : ker_arrays
  use mod_eccoupling,     only : kfl_exmsld_ecc, eccou_assign_fields, eccou_allocate_memory
  implicit none
  !
  ! Primary arrays
  !
  call ker_arrays('ALLOCATE')
  
  if( INOTMASTER ) then
     !
     ! WALLD: wall distance
     ! 
     if( kfl_walld /= 0 ) then
        call memory_alloca(mem_modul(1:2,modul),'WALLD'    ,'ker_memall',walld,    npoin)
        if (kfl_walld == 2 .or. kfl_walld == 3) then 
          call memory_alloca(mem_modul(1:2,modul),'WALLO'    ,'ker_memall',wallo,    npoin)
          call memory_alloca(mem_modul(1:2,modul),'WALLCOOR'    ,'ker_memall',wallcoor,  ndime,  npoin)
        endif
        call memory_alloca(mem_modul(1:2,modul),'UWALL_KER','ker_memall',uwall_ker,npoin)
        if( kfl_delta == 1 ) then
           call memory_alloca(mem_modul(1:2,modul),'UWAL2_KER','ker_memall',uwal2_ker,npoin)
        end if
     end if
     !
     ! WALLN: wall normal
     ! 
     if( kfl_walln /= 0 ) then
        call memory_alloca(mem_modul(1:2,modul),'WALLN'    ,'ker_memall',walln,ndime,npoin)
     end if
     !
     ! ROUGH: Allocate roughness
     !
     if( kfl_rough == 0 ) then 
        call memory_alloca(mem_modul(1:2,modul),'ROUGH','ker_memall',rough,npoin)
     end if
     !
     ! CANHE: Allocate canopy heigh 
     !
     if( kfl_canhe == 0 ) then 
        call memory_alloca(mem_modul(1:2,modul),'CANHE','ker_memall',canhe,npoin)
     end if
     !
     ! Support boundary
     !
     if( kfl_suppo == 1 ) then
        call memory_alloca(mem_modul(1:2,modul),'DISPL_KER','ker_memall',displ_ker,ndime,npoin)
     end if
     !
     ! adjoint optimization variables
     !
     if( kfl_adj_prob == 1 .and. kfl_dvar_type /= 5 .and. kfl_dvar_type /= 6) then
        call memory_alloca(mem_modul(1:2,modul),'SENS','ker_memall',sens,kfl_ndvars_opt)
     elseif( kfl_adj_prob == 1 .and. kfl_dvar_type == 5) then
        call memory_alloca(mem_modul(1:2,modul),'SENS_MESH','ker_memall',sens_mesh,ndime,npoin)
     elseif( kfl_adj_prob == 1 .and. kfl_dvar_type == 6) then
        call memory_alloca(mem_modul(1:2,modul),'SENS_MESH','ker_memall',sens_mesh,ndime,1_ip)
     end if
     !
     ! Mesh velocity
     !
     if( kfl_difun /= 0 ) then
        call memory_alloca(mem_modul(1:2,modul),'VELOM','ker_memall',velom,ndime,npoin)       
     end if
     !
     ! Force based on residual for rigid body motion
     !
     if( kfl_forca_res == 1_ip ) then
        call memory_alloca(mem_modul(1:2,modul),'FORCA','ker_memall',forca,ndime,npoin,3_ip)
     end if
     !
     ! DISPM and VELOM
     !
     if( ker_subdomain_motion_exists() ) then
        call memory_alloca(mem_modul(1:2,modul),'KFL_FIXBO_WALLN_KER','ker_memory' , dispm , ndime , npoin , 3_ip )
        call memory_alloca(mem_modul(1:2,modul),'KFL_FIXBO_WALLN_KER','ker_memory' , velom , ndime , npoin )
     end if
     !
     ! no slip wall
     !
     if( kfl_noslw_ker /= 0 ) then
        call memory_alloca(mem_modul(1:2,modul),'EL_NSW_VISC'    ,'ker_memall',el_nsw_visc    ,nelem)
     end if

  else
     !
     ! Wall distance
     !
     if( kfl_walld /= 0 ) then
        call memory_alloca(mem_modul(1:2,modul),'WALLD'    ,'ker_memall',walld,    1_ip)
        if (kfl_walld == 2 .or. kfl_walld == 3) then
          call memory_alloca(mem_modul(1:2,modul),'WALLO'    ,'ker_memall',wallo,    1_ip)
          call memory_alloca(mem_modul(1:2,modul),'WALLCOOR'    ,'ker_memall',wallcoor,  ndime,  1_ip)
        endif  
        call memory_alloca(mem_modul(1:2,modul),'UWALL_KER','ker_memall',uwall_ker,1_ip)
        if( kfl_delta == 1 ) then
           call memory_alloca(mem_modul(1:2,modul),'UWAL2_KER','ker_memall',uwal2_ker,1_ip)
        end if
     end if
     !
     ! ROUGH: Allocate roughness
     !
     if( kfl_rough == 0 ) then 
        call memory_alloca(mem_modul(1:2,modul),'ROUGH','ker_memall',rough,1_ip)
     end if
     !
     ! CANHE: Allocate canopy heigh 
     !
     if( kfl_canhe == 0 ) then 
        call memory_alloca(mem_modul(1:2,modul),'CANHE','ker_memall',canhe,1_ip)
     end if
     !
     ! Support boundary
     !
     if( kfl_suppo == 1 ) then
        call memory_alloca(mem_modul(1:2,modul),'DISPL_KER','ker_memall',displ_ker,1_ip,1_ip)
     end if
     !
     ! Fields prescribed by master
     !
     if( kfl_vefun /= 0 ) then
        call memory_alloca(mem_modul(1:2,modul),'ADVEC','ker_memall',advec,1_ip,1_ip,3_ip)
     else
        advec => veloc
     end if
     if( kfl_tefun /= 0 ) then
        call memory_alloca(mem_modul(1:2,modul),'TEMPE','ker_memall',tempe,1_ip,3_ip)
     end if
     if( kfl_cofun /= 0 ) then
        call memory_alloca(mem_modul(1:2,modul),'CONCE','ker_memall',conce,1_ip,1_ip,3_ip)
     end if
     if( kfl_arfun /= 0 ) then
        call memory_alloca(mem_modul(1:2,modul),'AREAS','ker_memall',areas,1_ip,3_ip)
     end if
     !
     ! Force based on residual for rigid body motion
     !
     if( kfl_forca_res == 1_ip  )then
        call memory_alloca_min(mem_modul(1:2,modul),'FORCA','ker_memall',forca)
     end if
     !
     ! adjoint optimization variables
     !
     if( kfl_adj_prob == 1 .and. kfl_dvar_type /= 5 .and. kfl_dvar_type /= 6) then
        call memory_alloca(mem_modul(1:2,modul),'SENS','ker_memall',sens,1_ip)
     elseif( kfl_adj_prob == 1 .and. kfl_dvar_type == 5) then
        call memory_alloca(mem_modul(1:2,modul),'SENS_MESH','ker_memall',sens_mesh,ndime,1_ip)
     elseif( kfl_adj_prob == 1 .and. kfl_dvar_type == 6) then
        call memory_alloca(mem_modul(1:2,modul),'SENS_MESH','ker_memall',sens_mesh,ndime,1_ip)
     end if
     !
     ! no slip wall
     !
     if( kfl_noslw_ker /= 0 ) then
        call memory_alloca(mem_modul(1:2,modul),'EL_NSW_VISC'    ,'ker_memall',el_nsw_visc    ,1_ip)   ! guess his is the correect for master        
     end if

  end if

  !
  ! Allocate memmory for the ECCOUPLING
  !
  if( kfl_exmsld_ecc )then
     call eccou_assign_fields()
     call eccou_allocate_memory(1_ip)
  endif

end subroutine ker_memall
