!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Partis
!> @{
!> @file    mod_pts_thermodynamic.f90
!> @author  houzeaux
!> @date    2018-09-21
!> @brief   Module for particle thermodynamic
!> @details Contains all subroutine to compute particles thermodynamic
!>          - Specific enthalpy: Sensible Heat, it is the quantity of heat
!>            contained in 1 kg of water according to the selected temperature.
!>          
!-----------------------------------------------------------------------

module mod_pts_thermodynamic

  use def_partis
  use def_master
  use def_partis,      only : PTS_PARTICLE_EVAPORATED
  use def_partis,      only : PTS_PARTICLE_EXISTS
  use def_parame,      only : pi,oneo3
  use mod_ker_proper,  only : ker_proper
  use mod_ker_updpro,  only : ker_updpro
  use mod_physics,     only : physics_sphere_diameter
  use mod_physics,     only : physics_sphere_mass
  use mod_physics,     only : physics_air_water_vapor_diffusion_coefficient
  use mod_physics,     only : universal_gas_constant

  use mod_physics,     only : physics_set_liquid_temperature
  use mod_physics,     only : physics_sphere_Nusselt
  use mod_physics,     only : physics_sphere_Sherwood
  use mod_physics,     only : physics_mole_2_mass 
  use mod_physics,     only : physics_mass_2_w 
  use mod_physics,     only : physics_T_2_HCp
  use mod_physics,     only : physics_H_2_TCp
  use mod_physics,     only : liquid_state 
   
  use mod_interp_tab,  only : fw_lookup,tab_interp 

  use mod_pts_particle, only: pts_particle_diameter
  
  implicit none

  real(rp), parameter :: w_fluid   = 28.9647e-3_rp
  
  private

  public :: pts_thermodynamic_properties
  public :: pts_thermodynamic_transport

contains


  !
  ! Get mean material properties for flamelet model
  !
  subroutine pts_lookup_mean_properties(itype, Tp, liq, xvap, Therm_seen, conce_seen, &
                                       confl, sphfl, denfl, visfl, Dvg_m, y_surf, y_seen, T_seen,w_nonf_seen)
        integer(ip),             intent(in)  :: itype 
        real(rp),                intent(in)  :: Tp 
        type(liquid_state),      intent(inout)  :: liq
        real(rp),                intent(in)  :: xvap 
        real(rp),                intent(in)  :: Therm_seen 
        real(rp),                intent(in)  :: conce_seen(:) 
        real(rp),                intent(out) :: confl 
        real(rp),                intent(out) :: sphfl
        real(rp),                intent(out) :: denfl
        real(rp),                intent(out) :: visfl
        real(rp),                intent(out) :: Dvg_m
        real(rp),                intent(out) :: y_surf
        real(rp),                intent(out) :: y_seen
        real(rp),                intent(out) :: T_seen
        real(rp),                intent(out) :: w_nonf_seen
     
        real(rp)    :: w_seen, &
                       T_mean, dT_mean, w_mean, w_surf, &
                       con, dcon_dT, vis, dvis_dT, Dvg, dDvg_dT, &
                       cpcoef(6,2), H_mean, error
                       
        real(rp)    :: seen_retva(parttyp(itype) % table_fw % main_table % nvar)
        real(rp)    :: mean_retva(parttyp(itype) % table_fw % main_table % nvar)
        real(rp)    :: seen_spr_retva(parttyp(itype) % spr_tab_fw % main_table % nvar)
        real(rp)    :: mean_spr_retva(parttyp(itype) % spr_tab_fw % main_table % nvar)
        real(rp)    :: control(parttyp(itype) % table_fw % main_table % ndim)
        real(rp)    :: mean_control(parttyp(itype) % table_fw % main_table % ndim)
        real(rp)    :: tab_scale_control(parttyp(itype) % table_fw % main_table % ndim)
        integer(ip) :: ind(parttyp(itype) % table_fw % main_table % ndim) 

        integer(ip) :: jj, ind_entha, idimt

        !
        ! Lookup seen state 
        !        
        control = 0.0_rp
        ind     = 1_ip
        ind_entha = 0
        do idimt = 1, parttyp(itype) % table_fw % main_table % ndim
           if (parttyp(itype) % table_fw % kfl_chm_control(idimt) > 0) then
              !
              ! >0: one of the conces
              !
              control(idimt) = conce_seen(parttyp(itype) % table_fw % kfl_chm_control(idimt))
           else
              if (parttyp(itype) % table_fw % kfl_chm_control(idimt) == -1) then
                 !
                 ! -1: enthalpy
                 !
                 control(idimt) = Therm_seen
                 ind_entha      = idimt
              elseif (parttyp(itype) % table_fw % kfl_chm_control(idimt) == -2) then
                 !
                 ! -2: scalar dissipation rate
                 !
                 call runend('pts_lookup_mean_properties: droplet model in not compatible with unsteady flamelet progress variable method.')
              endif
           endif
        enddo
        call fw_lookup( control, tab_scale_control, parttyp(itype) % table_fw, seen_retva, ind )
        call tab_interp( parttyp(itype) % spr_tab_fw % main_table, tab_scale_control, seen_spr_retva, ind )

        !
        ! Get seen liquid mass fraction, and seen molar mass
        !
        y_seen = seen_spr_retva(parttyp(itype) % kfl_Yfuel_spr_index)
        w_seen = seen_retva(parttyp(itype) % kfl_W_tab_index) 

        !
        ! First is low temperature, second is high temperature in table
        !
        cpcoef(1:6,1) = seen_retva(parttyp(itype)%kfl_cpCoefLT_tab_index:parttyp(itype)%kfl_cpCoefLT_end_tab_index)
        cpcoef(1:6,2) = seen_retva(parttyp(itype)%kfl_cpCoefHT_tab_index:parttyp(itype)%kfl_cpCoefHT_end_tab_index)

        !
        ! This is temperature from enthalpy and tabulated NASA polynomials
        ! Use tabulated temperature as initial guess
        ! sphfl is calculated here in case weight_seen == 1
        ! 
        T_seen = seen_spr_retva(parttyp(itype) % kfl_T_spr_index)
        call physics_H_2_TCp(Therm_seen, cpcoef, T_seen, sphfl) 
        

        !
        ! Molar mass of seen mixture without the fuel
        !
        if (y_seen < 1.0_rp) then
            w_nonf_seen = (1.0_rp - y_seen) / (1.0_rp / w_seen - y_seen / liq % W )
        else
            w_nonf_seen = 28.0e-3_rp
        endif
        
        !
        ! Mass fraction of mixing the fuel with this mixture
        !
        y_surf = physics_mole_2_mass(xvap, liq % W, w_nonf_seen)
        w_surf = physics_mass_2_w(y_surf, liq % W, w_nonf_seen) 

        !
        ! Mean control variables, for the surface, only consider the fuel as
        ! mixture fraction (basically taking the mixing line state)
        ! The other control variables are 0 in this mixing step
        !
        mean_control = 0.0_rp
        do idimt = 1, parttyp(itype) % table_fw % main_table % ndim
           if (parttyp(itype) % table_fw % kfl_chm_control(idimt) > 0) then
              !
              ! >0: one of the conces: take seen portion
              !
              mean_control(idimt) = parttyp(itype) % weight_seen * conce_seen(parttyp(itype) % table_fw % kfl_chm_control(idimt))
           else
              if (parttyp(itype) % table_fw % kfl_chm_control(idimt) == -1) then
                 !
                 ! -1: enthalpy
                 !
                 mean_control(idimt) = Therm_seen
              endif
           endif
           !
           ! Add surface portion for mixture fraction
           !
           if (parttyp(itype) % table_fw % main_table % coords(idimt) % name == 'ZMEAN' .or.&
               parttyp(itype) % table_fw % main_table % coords(idimt) % name == 'Z    ') then
               mean_control(idimt) = mean_control(idimt) + ( 1.0_rp - parttyp(itype) % weight_seen ) * y_surf
           endif

        enddo
        
        T_mean = parttyp(itype) % weight_seen * T_seen +  ( 1.0_rp - parttyp(itype) % weight_seen ) * Tp 
        
        
       
        if (parttyp(itype) % weight_seen == 1.0_rp) then
             mean_retva     = seen_retva
             mean_spr_retva = seen_spr_retva
        else
            if (parttyp(itype) % table_fw % kfl_needs_enthalpy /= 0_ip) then
                !
                ! Solve Newton-Raphson for heat loss
                !
                error = 1.0e15_rp
                jj = 0_ip
                do while( (abs(error) > 1.0e-8_rp) .and. (jj<10_ip)  )
                    !
                    ! Look up cp coefficients and tabulated temperature
                    !
                    call fw_lookup( mean_control, tab_scale_control, parttyp(itype) % table_fw, mean_retva, ind )
                    
                    !
                    ! Specific heat at mean temperature: Plug into CP coeffs
                    !
                    cpcoef(1:6,1) = mean_retva(parttyp(itype)%kfl_cpCoefLT_tab_index:parttyp(itype)%kfl_cpCoefLT_end_tab_index)
                    cpcoef(1:6,2) = mean_retva(parttyp(itype)%kfl_cpCoefHT_tab_index:parttyp(itype)%kfl_cpCoefHT_end_tab_index)
                    call physics_T_2_HCp(T_mean, cpcoef, H_mean, sphfl)

                    !
                    ! Update error
                    !
                    error = (H_mean - mean_control(ind_entha) )/abs(H_mean)
                    jj = jj + 1_ip 
                    
                    !
                    ! Update enthalpy in control variable
                    !
                    mean_control(ind_entha) =  H_mean

                    !
                    ! Check if we hit limits of table
                    !
                    if ( abs(tab_scale_control(ind_entha)-0.0_rp) < 1.0e-5_rp &
                                                      .and. error < 0.0_rp  ) then
                        error = 0.0_rp
                    endif
                    
                    if ( abs(tab_scale_control(ind_entha)-1.0_rp) < 1.0e-5_rp &
                                                      .and. error > 0.0_rp  ) then
                        error = 0.0_rp
                    endif 

                enddo
                   
            else
                !
                ! Pass through for adiabatic tables
                !
                call fw_lookup( mean_control, tab_scale_control, parttyp(itype) % table_fw, mean_retva, ind )

                !
                ! Specific heat at mean temperature: Plug into CP coeffs
                !
                cpcoef(1:6,1) = mean_retva(parttyp(itype)%kfl_cpCoefLT_tab_index:parttyp(itype)%kfl_cpCoefLT_end_tab_index)
                cpcoef(1:6,2) = mean_retva(parttyp(itype)%kfl_cpCoefHT_tab_index:parttyp(itype)%kfl_cpCoefHT_end_tab_index)
                call physics_T_2_HCp(T_mean, cpcoef, H_mean, sphfl)
                 
            endif
        endif

        call tab_interp( parttyp(itype) % spr_tab_fw % main_table, tab_scale_control, mean_spr_retva, ind )
            
        dT_mean = T_mean - mean_spr_retva(parttyp(itype) % kfl_T_spr_index)
        
        !
        ! Correct properties for temperature difference
        ! between T_mean and the tabulated temperature
        !
        con = mean_retva(parttyp(itype) % kfl_k_tab_index)
        vis = mean_retva(parttyp(itype) % kfl_mu_tab_index)
        Dvg = mean_spr_retva(parttyp(itype) % kfl_Dfuel_spr_index)
        dcon_dT = mean_spr_retva(parttyp(itype) % kfl_dkdT_spr_index    )
        dvis_dT = mean_spr_retva(parttyp(itype) % kfl_dmudT_spr_index   )
        dDvg_dT = mean_spr_retva(parttyp(itype) % kfl_dDfueldT_spr_index)

        confl   = con + dT_mean * dcon_dT
        visfl   = vis + dT_mean * dvis_dT
        Dvg_m   = Dvg + dT_mean * dDvg_dT
 
        !
        ! Density from ideal gas law
        ! 
        w_mean  = physics_mass_2_w(parttyp(itype) % weight_seen, w_seen, w_surf)
        denfl   = liq % P * w_mean / (universal_gas_constant * T_mean)

  end subroutine pts_lookup_mean_properties


  !
  ! Get gas properties as mean properties  
  !
  subroutine pts_basic_mean_properties(itype, Tp, liq, xvap, T_seen, ielem, pnode, shapf, &
                                       confl, sphfl, denfl, visfl, Dvg_m, y_surf, w_nonf_seen )
        integer(ip),             intent(in)  :: itype 
        real(rp),                intent(in)  :: Tp 
        type(liquid_state),      intent(inout)  :: liq
        real(rp),                intent(in)  :: xvap 
        real(rp),                intent(in)  :: T_seen 
        integer(ip),             intent(in)  :: ielem 
        integer(ip),             intent(in)  :: pnode 
        real(rp),                intent(in)  :: shapf(pnode)
        real(rp),                intent(out) :: confl 
        real(rp),                intent(out) :: sphfl
        real(rp),                intent(out) :: denfl
        real(rp),                intent(out) :: visfl
        real(rp),                intent(out) :: Dvg_m
        real(rp),                intent(out) :: y_surf
        real(rp),                intent(out) :: w_nonf_seen

        integer(ip) :: dumm0
        
        !
        ! Calculate surface vapor mass fraction
        !
        w_nonf_seen = w_fluid
        y_surf  = physics_mole_2_mass(xvap, liq % W, w_nonf_seen)
        
        !
        ! Calculate binary diffusion coefficient of vapor in the fluid
        !
        Dvg_m   = physics_air_water_vapor_diffusion_coefficient(1_ip,T_seen)

        !
        ! Thermal conductivity, and specific heat from kernel functions
        !
        call ker_proper('CONDU','IGAUS',dumm0,ielem,confl,pnode,1_ip,shapf) 
        call ker_proper('SPHEA','IGAUS',dumm0,ielem,sphfl,pnode,1_ip,shapf)
        call ker_proper('DENSI','IGAUS',dumm0,ielem,denfl,pnode,1_ip,shapf)       
        call ker_proper('VISCO','IGAUS',dumm0,ielem,visfl,pnode,1_ip,shapf)       

  end subroutine pts_basic_mean_properties


  !
  ! Get mean material properties with seen quantities
  !
  subroutine pts_mean_properties_with_seen(itype, Tp, liq, Therm_seen, ielem, pnode, shapf, &
                                           confl, sphfl, denfl, visfl, Dvg_m, Pr_m, Sc_m, LK_m, y_surf, y_seen, T_seen, &
                                           xvap,w_nonf_seen,conce_seen)
        integer(ip),             intent(in)  :: itype
        real(rp),                intent(in)  :: Tp 
        type(liquid_state),      intent(inout)  :: liq
        real(rp),                intent(in)  :: Therm_seen 
        integer(ip),             intent(in)  :: ielem 
        integer(ip),             intent(in)  :: pnode 
        real(rp),                intent(in)  :: shapf(pnode)
        real(rp),                intent(out) :: confl 
        real(rp),                intent(out) :: sphfl
        real(rp),                intent(out) :: denfl
        real(rp),                intent(out) :: visfl
        real(rp),                intent(out) :: Dvg_m
        real(rp),                intent(out) :: Pr_m
        real(rp),                intent(out) :: Sc_m
        real(rp),                intent(out) :: LK_m
        real(rp),                intent(out) :: y_surf
        real(rp),                intent(out) :: y_seen
        real(rp),                intent(out) :: T_seen
        real(rp),                intent(out) :: xvap
        real(rp),                intent(out) :: w_nonf_seen
        real(rp),                intent(in)  :: conce_seen(:) 

        real(rp) :: alfa_eps

        !
        ! Update liquid properties
        !
        call physics_set_liquid_temperature( liq, Tp) 
    
        !
        ! Calculate surface vapor mole fraction
        !
        xvap    = liq % psat / liq % P

        !
        ! Calculate liquid properties, mean material properties, and droplet surface mass fraction
        ! 
        select case(parttyp(itype) % kfl_therm)
            case(1)
                !
                !  Mean properties are the seen gas properties
                !
                call pts_basic_mean_properties(itype, Tp, liq, xvap, Therm_seen, ielem, pnode, shapf, &
                                               confl, sphfl, denfl, visfl, Dvg_m, y_surf, w_nonf_seen )
                y_seen = conce_seen(1)
                T_seen = Therm_seen
            case(2)
                !
                !  Mean properties according to "1/3" law using tabulated values
                !
                call pts_lookup_mean_properties(itype, Tp, liq, xvap, Therm_seen, conce_seen, &
                                                  confl, sphfl, denfl, visfl, Dvg_m, y_surf, y_seen, T_seen, w_nonf_seen)
        end select
        
        !
        ! Prandtl and Schmidt numbers
        !
        Pr_m = sphfl * visfl / confl
        Sc_m = visfl / ( denfl * Dvg_m )

        !
        ! Knudsen layer thickness
        !
        alfa_eps = 1.0_rp ! molecular accommodation coefficient
        LK_m    = visfl * sqrt( 2.0_rp * pi * Tp * universal_gas_constant / liq % W ) / &
                  (alfa_eps * Sc_m * liq % P)


  end subroutine pts_mean_properties_with_seen

  !
  ! Get mean material properties at shape function
  !
  subroutine pts_thermodynamic_properties(itype, Tp, ielem, pnode, lnods, shapf, &
                                          confl, sphfl, denfl, visfl, Dvg_m, Pr_m, Sc_m, LK_m, y_surf, &
                                          y_seen, Therm_seen, T_seen, xvap, w_nonf_seen, conce_seen)
        integer(ip),             intent(in)  :: itype
        real(rp),                intent(in)  :: Tp 
        integer(ip),             intent(in)  :: ielem 
        integer(ip),             intent(in)  :: pnode 
        integer(ip),             intent(in)  :: lnods(pnode) 
        real(rp),                intent(in)  :: shapf(pnode)
        real(rp),                intent(out) :: confl 
        real(rp),                intent(out) :: sphfl
        real(rp),                intent(out) :: denfl
        real(rp),                intent(out) :: visfl
        real(rp),                intent(out) :: Dvg_m
        real(rp),                intent(out) :: Pr_m
        real(rp),                intent(out) :: Sc_m
        real(rp),                intent(out) :: LK_m
        real(rp),                intent(out) :: y_surf
        real(rp),                intent(out) :: y_seen
        real(rp),                intent(out) :: Therm_seen 
        real(rp),                intent(out) :: T_seen 
        real(rp),                intent(out) :: xvap 
        real(rp),                intent(out) :: w_nonf_seen 
        real(rp),                intent(out) :: conce_seen(:) 
        
        integer(ip) :: inode, ipoin
        real(rp) :: conce_seen_loc(nclas_pts)


        Therm_seen = 0.0_rp
        conce_seen_loc = 0.0_rp
        conce_seen = 0.0_rp
        do inode = 1,pnode
            ipoin = lnods(inode)
            Therm_seen  = Therm_seen  + shapf(inode) * therm(ipoin,1)
            if (nclas_pts == 1) then
                conce_seen_loc  = conce_seen_loc  + shapf(inode) * conce(ipoin,1,1)
            else
                conce_seen_loc  = conce_seen_loc  + shapf(inode) * conce(ipoin,1:nclas_pts,1)
            endif
        end do


        call pts_mean_properties_with_seen(itype, Tp, parttyp(itype) % liq, Therm_seen, ielem, pnode, shapf, &
                                           confl, sphfl, denfl, visfl, Dvg_m, Pr_m, Sc_m, LK_m, y_surf, y_seen, T_seen,  &
                                           xvap, w_nonf_seen, conce_seen_loc)
        if (nclas_pts == 1) then
            conce_seen = conce_seen_loc(1_ip)
        else
            conce_seen = conce_seen_loc(1_ip:min(nclas_pts, size(conce_seen, KIND=ip)))
        endif

  end subroutine pts_thermodynamic_properties





  !
  ! Mass transfer potential
  !
  subroutine pts_mass_transfer_poptential(kfl_mass_pot, Yv_surf, Yv_fluid_k, potential, BM)
    integer(ip), intent(in)    :: kfl_mass_pot
    real(rp),    intent(in)    :: Yv_surf
    real(rp),    intent(in)    :: Yv_fluid_k 
    real(rp),    intent(out)   :: potential
    real(rp),    intent(out)   :: BM

    !
    ! Mass transfer number
    !
    select case(abs(kfl_mass_pot))
    case(1,2)
        !
        ! Spalding number
        !
        BM = max(-1.0_rp+1.0e-6_rp, min(1.0e6_rp, (Yv_surf - Yv_fluid_k) / ( 1.0_rp - min(1.0_rp-1.0e-12_rp, Yv_surf))))
    case(3)
        !
        ! Dilute limit
        !
        BM = (Yv_surf - Yv_fluid_k)
    end select

    !
    ! Potential
    !
    select case(abs(kfl_mass_pot))
    case(1)
        !
        ! Quasi-steady 
        !
        potential = log( 1.0_rp + BM )
    case(2,3)
        !
        ! Mass analogy
        !
        potential = BM
    end select
  end subroutine pts_mass_transfer_poptential



  real(rp) function birds_func(beta)
    real(rp), intent(in) :: beta
!    real(rp) :: denom
    !
    ! birds_func = beta / (exp(beta) - 1)
    !
    if (abs(beta) < 1.0e-10_rp) then
        !
        ! At 0: denom ~= 1 + beta + beta**2 / 2 - 1 = beta * ( 1 + 0.5 beta )
        !
        birds_func = 1.0_rp / (1.0_rp + 0.5_rp * beta)
    else
        birds_func = beta / (exp(beta) - 1.0_rp)
    endif
  end function 

  real(rp) function inv_abramzon_thickness_ratio(Bin)
    real(rp), intent(in) :: Bin
    real(rp)             :: B
    !
    ! abramzon_thickness_ratio = ((1+B)^0.7)/B * ln(1+B)
    ! inv_abramzon_thickness_ratio = B/((1+B)^0.7 * ln(1+B))
    !
    B = min(20.0_rp,Bin)
    if (abs(B) < 1.0e-10_rp) then
        inv_abramzon_thickness_ratio = 1.0_rp - 0.2_rp * B 
    else
        inv_abramzon_thickness_ratio = B / ( (1.0_rp+B)**0.7_rp * log(1.0_rp+B) )
    endif
  end function 


  !
  ! Residual evaluation
  !
  subroutine pts_thermodynamic_residual(itype, mass_k, tempe_k, & ! variabbles
                                        u_slip, ielem, pnode, lnods, shapf, & ! external effects
                                        resid, conv_heat_flux, Therm_fluid_k, &
                                        T_fluid_k, Yv_fluid_k, conce_seen, &
                                        interpFluidProps, dt_k, mass_doom, writeProperties, BM) ! residual
    !>          dm                                                          
    !>          -- =  - pi * diam * Dvg_m * Sh_m * Potential             
    !>          dt
    !>         
    !>          Potential:  
    !>             
    !>                   /      Yv_surf - Yv_fluid_k   \
    !>             1) log| 1 + ----------------------- |   Quasi-steady                 
    !>                   \         1    -  Yv_surf     /                    
    !>             
    !>             
    !>                 Yv_surf - Yv_fluid_k                       
    !>             2) -----------------------              Mass analogy
    !>                    1    -  Yv_surf                        
    !>             
    !>             
    !>             3) Yv_surf - Yv_fluid_k                 Dilute mass analogy
    !>                                                                      
    !>                                                                      
    !>               dT                                                 dm   
    !>          m Cp -- = pi * diam * cond_m * Nu_m ( T_fluid_k - T ) + -- Lv
    !>               dt                                                 dt   
    integer(ip), intent(in)    :: itype
    real(rp),    intent(in)    :: mass_k
    real(rp),    intent(in)    :: tempe_k
    real(rp),    intent(in)    :: u_slip
    integer(ip), intent(in)    :: ielem 
    integer(ip), intent(in)    :: pnode 
    integer(ip), intent(in)    :: lnods(pnode) 
    real(rp),    intent(in)    :: shapf(pnode)
    real(rp),    intent(out)   :: resid(2)
    real(rp),    intent(out)   :: conv_heat_flux
    real(rp),    intent(inout) :: Therm_fluid_k 
    real(rp),    intent(inout) :: T_fluid_k 
    real(rp),    intent(inout) :: Yv_fluid_k 
    real(rp),    intent(inout) :: conce_seen(nclas_pts)
    logical(lg), intent(in)    :: interpFluidProps
    real(rp),    intent(in)    :: dt_k
    real(rp),    intent(in)    :: mass_doom
    logical(lg), optional, intent(in)  :: writeProperties
    real(rp),    optional, intent(out) :: BM
   
    real(rp)    :: diame, denfl, visfl, Re_m
    real(rp)    :: Dvg_m, Yv_surf, Sh_m, potential, BM_loc
    real(rp)    :: confl, sphfl, Nu_m
    real(rp)    :: LK_m, Pr_m, Sc_m, cp_vap, dummr, T_mean, beta
    real(rp)    :: mass_flux,xvap,xvap_neq,Yv_surf_neq,w_nonf_seen, Sh0, Nu0
    real(rp)    :: invFM, invFT, BT, Phi, PhiStar, Nu_old 
    integer(ip) :: ias
 
    !
    ! Exit for non-positive droplet size
    ! 
    if (mass_k <= 0.0_rp) then
        resid = 0.0_rp
        return
    endif

    !
    ! Particle properties, and representative gas phase properties
    !
    if (interpFluidProps) then
        call pts_thermodynamic_properties(itype, tempe_k, ielem, pnode, lnods, shapf, &
                                          confl, sphfl, denfl, visfl, Dvg_m, Pr_m, Sc_m, LK_m, Yv_surf, Yv_fluid_k, &
                                          Therm_fluid_k, T_fluid_k, xvap, w_nonf_seen, conce_seen) 
    else
        call pts_mean_properties_with_seen(itype, tempe_k, parttyp(itype) % liq, Therm_fluid_k, ielem, pnode, shapf, &
                                           confl, sphfl, denfl, visfl, Dvg_m, Pr_m, Sc_m, LK_m, Yv_surf, Yv_fluid_k, T_fluid_k,  &
                                           xvap, w_nonf_seen, conce_seen)
    endif

    diame = physics_sphere_diameter(mass_k,parttyp(itype) % liq % rho) 
  
    !
    ! Dimensionless numbers
    !
    Re_m      = denfl * diame * u_slip / visfl
    Sh0       = physics_sphere_Sherwood(1_ip,Re_m,Dvg_m,Schmidt=Sc_m)
    Nu0       = physics_sphere_Nusselt(1_ip,Re_m,confl,Prandtl=Pr_m)
    Sh_m      = Sh0
    Nu_m      = Nu0
 
    !
    ! Eqilibrium mass flux 
    !
    call pts_mass_transfer_poptential(parttyp(itype) % kfl_mass_pot, Yv_surf, Yv_fluid_k, potential, BM_loc)

    !
    ! Special corrections
    !
    if ((parttyp(itype) % kfl_mass_pot < 0) .or. parttyp(itype) % kfl_heattr_corr /= 0) then
        !
        ! Evaporation rate factor (beta)
        ! 
        T_mean = parttyp(itype) % weight_seen * T_fluid_k +  ( 1.0_rp - parttyp(itype) % weight_seen ) * tempe_k 
        call physics_T_2_HCp(T_mean, parttyp(itype) % cpcoef_v_chm, dummr, cp_vap)
        Phi  = cp_vap * Sh0 * Pr_m / (sphfl * Nu0 * Sc_m)
        beta = Phi * log(1.0_rp+BM_loc)

        !
        ! Non-eqilibrium mass flux 
        !
        if (parttyp(itype) % kfl_mass_pot < 0) then
            xvap_neq  = max(0.0_rp,min(1.0_rp,xvap - (2.0_rp*LK_m/diame) * beta))
            Yv_surf_neq  = physics_mole_2_mass(xvap_neq, parttyp(itype) % liq % W, w_nonf_seen)
            call pts_mass_transfer_poptential(parttyp(itype) % kfl_mass_pot, Yv_surf_neq, Yv_fluid_k, potential, BM_loc) 
            ! 
            ! Update beta with new mass transfer number
            !
            beta = Phi * log(1.0_rp+BM_loc)
        endif

        !
        ! Heat flux corrections
        !
        select case(abs(parttyp(itype) % kfl_heattr_corr))
        case(1)
            !
            ! Bird's correction
            !
            Nu_m = Nu_m * birds_func(beta)
        case(2)
            !
            ! Abramzon & Sirignano model
            !
            invFM   = inv_abramzon_thickness_ratio(BM_loc)
            Sh_m    = 2.0_rp + invFM * ( Sh0-2.0_rp )
            ias     = 0_ip
            Nu_old  = 1e10_rp
            !
            ! Iterate
            !
            do while( (abs(Nu_old-Nu_m)/Nu_m > 1e-7_rp) .and. (ias <= 20) )
                Nu_old  = Nu_m
                PhiStar = cp_vap * Sh_m * Pr_m / (sphfl * Nu_m * Sc_m)
                BT      = (1.0_rp + BM_loc)**PhiStar - 1.0_rp
                invFT   = inv_abramzon_thickness_ratio(BT)
                Nu_m    = 2.0_rp + invFT * ( Nu0-2.0_rp )
                ias     = ias +1
            enddo

            !
            ! Add Stefan flow effect
            !
            Nu_m = Nu_m * birds_func(beta)

        end select

    endif
    !
    ! Culate mass flux
    !
    mass_flux = -1.0_rp * pi * diame * denfl * Dvg_m * Sh_m * potential

    !
    ! Limit mass flux to maximum flux that brings particle to mass_doom in dt_k time
    !
    if (mass_k + mass_flux * dt_k < mass_doom) mass_flux = (mass_doom-mass_k)/dt_k


    conv_heat_flux = pi * diame * confl * Nu_m * (T_fluid_k - tempe_k)
    resid(1) =  mass_flux
    resid(2) = ( conv_heat_flux &
                 + parttyp(itype) % liq % Lv * mass_flux )  / (mass_k * parttyp(itype) % liq % cp )
 

    !
    ! Output Spalding mass transfer number
    ! 
    if (present(BM)) BM=BM_loc


    !
    ! Output for debugging
    ! 
    if (present(writeProperties)) then
        if (writeProperties) then
            write(666,'(30E18.10)') cutim, tempe_k, mass_k, diame, Yv_surf, &
                                    Yv_fluid_k, T_fluid_k, sphfl, confl, visfl,&
                                    denfl, Dvg_m, Re_m, Nu_m, Sh_m
        endif
    endif

  end subroutine pts_thermodynamic_residual




  !
  ! Numerical evaluation of jacobian
  !
  subroutine pts_thermodynamic_jacobian(itype, mass_k, tempe_k, & ! variabbles
                                        u_slip, ielem, pnode, lnods, shapf, & ! external effects
                                        resid,jacob,conv_heat_flux, dt_k, mass_doom, T_fluid_k, Yv_fluid_k,BM)
    !        /   dRm     dRm  \
    !        |   ---     ---  |                                                                
    !        |   dm      dT   |                                                                
    !    J = |                |                                                                
    !        |   dRT     dRT  |                                                                
    !        |   ---     ---  |                                                                
    !        \   dm      dT   /                                                                

    integer(ip),        intent(in)  :: itype
    real(rp),           intent(in)  :: mass_k
    real(rp),           intent(in)  :: tempe_k
    real(rp),           intent(in)  :: u_slip
    integer(ip),        intent(in)  :: ielem 
    integer(ip),        intent(in)  :: pnode 
    integer(ip),        intent(in)  :: lnods(pnode) 
    real(rp),           intent(in)  :: shapf(pnode)
    real(rp),           intent(out) :: resid(2)
    real(rp),           intent(out) :: jacob(2,2)
    real(rp),           intent(out) :: conv_heat_flux
    real(rp),           intent(in)  :: dt_k
    real(rp),           intent(in)  :: mass_doom
    real(rp), optional, intent(out) :: T_fluid_k 
    real(rp), optional, intent(out) :: Yv_fluid_k 
    real(rp), optional, intent(out) :: BM 
   
    real(rp)    :: resim(2), resiT(2), drdy(2), dT, dm, Therm_fluid_k, T_fluid_k_loc, Yv_fluid_k_loc,&
                   conce_seen(nclas_pts)

    !
    ! Residual at m+dm, T
    !
    dm = max(10.0_rp*zeror, 0.001_rp * mass_k)
    call pts_thermodynamic_residual(itype, mass_k + dm, tempe_k, & 
                                    u_slip, ielem, pnode, lnods, shapf, & 
                                    resim, conv_heat_flux, Therm_fluid_k, &
                                    T_fluid_k_loc, Yv_fluid_k_loc, conce_seen, &
                                    .true., dt_k, mass_doom) 
    !
    ! Residual at m, T+dT
    !
    dT = -0.0001_rp * tempe_k
    call pts_thermodynamic_residual(itype, mass_k, tempe_k + dT, & 
                                    u_slip, ielem, pnode, lnods, shapf, & 
                                    resiT, conv_heat_flux, Therm_fluid_k, &
                                    T_fluid_k_loc, Yv_fluid_k_loc, conce_seen, &
                                    .false., dt_k, mass_doom) 
    !
    ! Residual at m,T
    !
    call pts_thermodynamic_residual(itype, mass_k, tempe_k, & 
                                    u_slip, ielem, pnode, lnods, shapf, & 
                                    resid, conv_heat_flux, Therm_fluid_k, &
                                    T_fluid_k_loc, Yv_fluid_k_loc, conce_seen, &
                                    .false., dt_k, mass_doom, BM=BM) 

    drdy = (resim-resid)/dm
    jacob(1,1) = drdy(1)
    jacob(2,1) = drdy(2)

    drdy = (resiT-resid)/dT
    jacob(1,2) = drdy(1)
    jacob(2,2) = drdy(2) 

    if (present(T_fluid_k))  T_fluid_k  = T_fluid_k_loc
    if (present(Yv_fluid_k)) Yv_fluid_k = Yv_fluid_k_loc
            
  end subroutine pts_thermodynamic_jacobian





  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-03-18
  !> @brief   Mass and temperature transport    
  !> @details Mass and temperature transport of droplets
  !>
  !>          w_fluid .............. molecular weight of the fluid   [ kg/mol ] 
  !>          T_fluid_k ............ fluid temperature               [ K ]
  !>          Yv_fluid_k ........... fluid vapor mass fraction       [ kg/kg ]
  !>          V .................... cell volume                     [ m^3 ]
  !>
  !>          w .................... molecular weight of the droplet [ kg/mol ]
  !>          A .................... sphere surface area             [ m2 ] 
  !>          diam ................. droplet diameter                [ m ] 
  !>          Cp ................... heat capacity                   [ J/(kg.K) ]
  !>          Lv ................... heat of vaporization            [ J/kg ]
  !>          h_heat ............... heat transfer coefficient       [ W/(m^2.K) ]
  !>          h_mass ............... mass transfer coefficient       [ m/s ]
  !>          qr ................... droplet radiation absorption    [ W ]
  !>          Yv_surf .............. vapor mass fraction on droplet surface [ kg/kg ]
  !>          particle % mass_k .... droplet mass                    [ kg ]
  !>          particle % tempe_k ... droplet temperature             [ K ]
  !>
  !>          The easy way to write these equations: 
  !>          dm
  !>          -- = A * h_mass * rho_f * (Yv_fluid_k - Yl)
  !>          dt
  !>
  !>               dT
  !>          m Cp -- = A h_heat(T_fluid_k-T) 
  !>               dt  
  !>                  + A h_mass * rho_f * (Yv_fluid_k - Yl) * Lv + qr
  !>                                                                      
  !>                                                                      
  !>          The common way to write these equations:                    
  !>          dm                                                          
  !>          -- =  - pi * diam * Dvg_m * Sh_m * Potential             
  !>          dt
  !>         
  !>          Potential:  
  !>             a) Yv_surf - Yv_fluid_k  
  !>                                                                      
  !>                   /      Yv_surf - Yv_fluid_k   \
  !>             b) log| 1 + ----------------------- |                    
  !>                   \         1    -  Yv_surf     /                    
  !>                                                                      
  !>               dT                                                 dm   
  !>          m Cp -- = pi * diam * cond_m * Nu_m ( T_fluid_k - T ) + -- Lv
  !>               dt                                                 dt   
  !>                                                              
  !-----------------------------------------------------------------------

  subroutine pts_thermodynamic_transport( &
       ilagr,particle,Re,ielem,pnode,lnods,shapf,dt_k, mass_kp1, tempe_kp1, conv_heat_flux)

    integer(ip), intent(in)    :: ilagr
    type(latyp), intent(inout) :: particle          !< Particle
    real(rp),    intent(in)    :: Re
    integer(ip), intent(in)    :: ielem
    integer(ip), intent(in)    :: pnode
    integer(ip), intent(in)    :: lnods(pnode)
    real(rp),    intent(in)    :: shapf(pnode)
    real(rp),    intent(in)    :: dt_k
    real(rp),    intent(out)   :: mass_kp1 
    real(rp),    intent(out)   :: tempe_kp1
    real(rp),    intent(out)   :: conv_heat_flux

    integer(ip)                :: itype
    integer(ip)                :: kfl_doom
    real(rp)                   :: mass_doom, dmini
    real(rp)                   :: T_fluid_k
    real(rp)                   :: Therm_fluid_k
    real(rp)                   :: diame
    real(rp)                   :: tempe_k,mass_k, cp_old
    real(rp)                   :: resid(2),jacob(2,2),aa(2,2), bb(2), invaa(2,2), detaa
    real(rp)                   :: conce_seen(nclas_pts), Yv_fluid_k, u_slip, urel(3)
    real(rp)                   :: BM 

    
    itype       = particle % itype
    urel        = particle % v_fluid_k - particle % veloc
    u_slip      = sqrt(dot_product(urel,urel))

    !
    ! Initial guess
    !
    tempe_k   = particle % tempe_k 
    mass_k    = particle % mass_k
    call physics_set_liquid_temperature( parttyp(itype) % liq , tempe_k)
    cp_old    = parttyp(itype) % liq % cp


    ! 
    ! Limit final diameter 
    !
    select case(abs(parttyp(itype) % kfl_dmini)) 
    case(0,1)
        !
        ! Use a constant diameter to eliminate or keep droplets
        !
        dmini     = parttyp(itype) % param_dmini 
        mass_doom = physics_sphere_mass(dmini, parttyp(itype) % liq % rho)
    case(2)
        !
        ! Use a secified fraction of the initial mass to eliminate or keep static droplets
        !
        mass_doom = parttyp(itype) % param_dmini * particle % mass_0
        dmini     =  physics_sphere_diameter(mass_doom, parttyp(itype) % liq % rho)
    case(3)
        !
        ! Use a secified fraction of the initial diameter to eliminate or keep static droplets
        !
        dmini     = parttyp(itype) % param_dmini * particle % diam_0
        mass_doom = physics_sphere_mass(dmini, parttyp(itype) % liq % rho)
    end select

    if ( parttyp(itype) % kfl_dmini < 0_ip ) then
        !
        ! Keep droplet at given size
        !
        kfl_doom  = PTS_PARTICLE_EXISTS
    else
        !
        ! Eliminate droplet
        !
        kfl_doom  = PTS_PARTICLE_EVAPORATED
        mass_doom = 0.0_rp
    endif



    

    !
    ! Time integration
    !
    select case(kfl_thermo_timsch_pts)
        case(0)
            !
            ! Explicit:
            ! 
            call pts_thermodynamic_residual(itype, mass_k, tempe_k,   & 
                                u_slip, ielem, pnode, lnods, shapf,   & 
                                resid, conv_heat_flux, Therm_fluid_k, &
                                T_fluid_k, Yv_fluid_k, conce_seen,    &
                                .true., dt_k, mass_doom, BM=BM) 
            mass_kp1  = mass_k  + dt_k * resid(1)
            tempe_kp1 = tempe_k + dt_k * resid(2)

        case(1)
            !
            ! Implicit  aa u^k+1 = r(u^k) + aa * u^k = bb
            !            with aa = I/dt_k - J
            !
            call pts_thermodynamic_jacobian(itype, mass_k, tempe_k, & 
                                u_slip, ielem, pnode, lnods, shapf, & 
                                resid, jacob, conv_heat_flux, dt_k, &
                                mass_doom, T_fluid_k, Yv_fluid_k, BM=BM) 

            aa      = -1.0_rp * jacob
            aa(1,1) = aa(1,1) + 1.0_rp/dt_k 
            aa(2,2) = aa(2,2) + 1.0_rp/dt_k 

            bb      = resid
            bb(1)   = bb(1) + mass_k * aa(1,1) + tempe_k * aa(1,2)
            bb(2)   = bb(2) + mass_k * aa(2,1) + tempe_k * aa(2,2)
            
            detaa   = aa(1,1)*aa(2,2) - aa(1,2)*aa(2,1)

            !
            ! Invert Jacobian, and use it for the timestep
            !
            invaa   = 1.0_rp / detaa
            invaa(1,1) = invaa(1,1) * aa(2,2)
            invaa(2,2) = invaa(2,2) * aa(1,1)
            invaa(1,2) = invaa(1,2) * (-1.0_rp) * aa(1,2)
            invaa(2,1) = invaa(2,1) * (-1.0_rp) * aa(2,1)

            mass_kp1  = invaa(1,1) * bb(1) + invaa(1,2) * bb(2)
            tempe_kp1 = invaa(2,1) * bb(1) + invaa(2,2) * bb(2)

    end select 
   
    !
    ! Postprocessing fluid variables
    !
    particle % Temp_fluid_k = T_fluid_k
    particle % Yvap_fluid_k = Yv_fluid_k
    particle % BM           = BM

    !
    ! Clip
    !
    mass_kp1   = max(0.0_rp, mass_kp1)
    tempe_kp1  = max(200.0_rp,min(parttyp(itype) % liq % Tsat, tempe_kp1 ))

    ! 
    ! Limit based on diameter 
    !
    call physics_set_liquid_temperature( parttyp(itype) % liq , tempe_kp1)
    diame = physics_sphere_diameter(mass_kp1, parttyp(itype) % liq % rho) 

    !
    ! Eliminate or keep static
    !
    if (diame <= dmini) then
        mass_kp1  = mass_doom
        tempe_kp1 = particle % tempe_k
        particle % kfl_exist = kfl_doom 
        particle % BM = 0.0_rp
    endif
   


  end subroutine pts_thermodynamic_transport
  

end module mod_pts_thermodynamic
  !> @}
