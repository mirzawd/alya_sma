!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tem_iniunk()
  !-----------------------------------------------------------------------
  !****f* Temper/tem_iniunk
  ! NAME 
  !    tem_iniunk
  ! DESCRIPTION
  !    This routine sets up the initial condition for the temperature.
  ! USED BY
  !    tem_begste
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_temper
  use def_kermod
  use mod_ker_proper 
  use mod_ker_space_time_function
  use mod_chktyp,         only : check_type
  use mod_communications, only : PAR_MAX, PAR_MIN
  use mod_physics,        only : physics_T_2_HCp
  use mod_ker_tendencies, only : ker_tendencies_initialize_th, kfl_tendencies_ker
  use mod_ker_updpro,     only : ker_updpro

  implicit none
  integer(ip) :: ipoin,ifunc
  integer(ip) :: kfl_advec_old,kfl_timei_old, kfl_value
  real(rp)    :: dtinv_tmp

  if( kfl_rstar == 0 ) then 

     !-------------------------------------------------------------------
     !
     ! Initial conditions
     !
     !-------------------------------------------------------------------

     if( kfl_inico_tem < - 100 ) then 
        !
        ! Take initial condition from a value function field
        ! 
        kfl_value = -  kfl_inico_tem - 100
        if( INOTEMPTY ) call check_type(xfiel,kfl_value,1_ip,nunkn_tem) ! Check if value function exist
        do ipoin = 1,nunkn_tem
           therm(ipoin,1) = xfiel(kfl_value) % a(1,ipoin,1)
        end do

     else if( kfl_inico_tem < 0 ) then 
        ! 
        ! Space time function
        !
        ifunc = -kfl_inico_tem
        do ipoin = 1,nunkn_tem
           call ker_space_time_function(&
                ifunc,coord(1:ndime,ipoin),cutim,therm(ipoin,1))
        end do

     else if( kfl_inico_tem == 1 ) then
        !
        ! Constant initial
        !
        do ipoin = 1,nunkn_tem
           therm(ipoin,1) = initial_tem
        end do

     else if( kfl_inico_tem == 3 ) then
        !
        ! Solve diffusion problem
        !        
        kfl_advec_old = kfl_advec_tem
        kfl_timei_old = kfl_timei_tem
        dtinv_tmp     = dtinv_tem
        kfl_advec_tem = 0
        kfl_timei_tem = 0
        dtinv_tem     = 0.0_rp
        call tem_updunk(ITASK_BEGITE)
        if( kfl_prope /= 0 ) call ker_updpro() ! Force update of kernel properties
        call tem_solite()
        call tem_updunk(ITASK_ENDITE)
        kfl_inidi_tem = 0
        kfl_advec_tem = kfl_advec_old 
        kfl_timei_tem = kfl_timei_old 
        dtinv_tem     = dtinv_tmp
        
     else if (kfl_tendencies_ker) then
        !  when no initialized, initialized from WRF Th in tendencies file
         
          if (INOTEMPTY) call ker_tendencies_initialize_th(therm, walld)
     else

        do ipoin = 1,nunkn_tem              
           therm(ipoin,1) = bvess_tem(1,ipoin,1)  
        end do

     end if

     !-------------------------------------------------------------------
     !
     ! Impose boundary conditions
     !
     !-------------------------------------------------------------------
     !
     ! Smooth Dirichlet boundary conditions
     !
     call tem_smobcs() 
     !
     ! Apply Dirichlet boundary conditions on (:,1)
     !
     call tem_updbcs(ITASK_INIUNK)

     !-------------------------------------------------------------------
     !
     ! Interpolation from coarse to fine mesh
     !
     !-------------------------------------------------------------------

     if( kfl_meshi_tem /= 0_ip ) call tem_coarfine(1_ip)

     call tem_coupli(ITASK_INIUNK)
     
     !-------------------------------------------------------------------
     !
     ! Copy (:,3) <- (:,1)
     !
     !-------------------------------------------------------------------

     call tem_updunk(ITASK_INIUNK)

  end if
  !
  ! Compute temperature from enthalpy if needed (:,1), even in the case of restart.
  ! This is needed to compute properties
  !
  call tem_temperature_from_enthalpy()

end subroutine tem_iniunk
