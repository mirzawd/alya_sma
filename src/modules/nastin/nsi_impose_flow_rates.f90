!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup NastinInput
!> @{
!> @file    nsi_impose_flow_rates.f90
!> @date    09/11/2022
!> @author  Guillaume Houzeaux
!> @brief   Impose flow rates
!> @details Impose flow rates
!> @} 
!-----------------------------------------------------------------------

subroutine nsi_impose_flow_rates()

  use def_parame,                  only :  pi
  use def_elmtyp
  use def_master
  use def_kermod
  use def_domain
  use def_nastin
  use mod_ker_space_time_function
  use mod_communications,          only : PAR_INTERFACE_NODE_EXCHANGE
  use mod_communications,          only : PAR_MIN
  use mod_couplings,               only : COU_PUT_VALUE_ON_TARGET
  use mod_projec,                  only : projec_mass_conservation
  use mod_memory,                  only : memory_alloca
  use mod_memory,                  only : memory_alloca_min
  use mod_memory,                  only : memory_deallo
  use mod_memory,                  only : memory_size
  use mod_pump,                    only : pump_curve_flow
  use mod_ker_functions,           only : ker_functions
  use mod_messages,                only : messages_live
  use mod_local_basis,             only : local_basis_global_to_local
  use mod_local_basis,             only : local_basis_local_to_global
  use mod_output,                  only : output_mesh_gid_format
  use mod_func,                    only : func_ptr, func_initialization, func_dotn
  use mod_integrals,               only : field_arrays, integrals_boundary
  use mod_std
  implicit none
  integer(ip)             :: ipoin,idime,iboun,dummi,itotv,izdom,jpoin
  integer(ip)             :: inodb,inode,ifunc,ibopo,jdime,kpoin,ii
  integer(ip)             :: iflow,ierro,itype,incode, itable, ibtim
  integer(ip)             :: i_left,i_right,ind_left,ind_right
  real(rp)                :: venew,veold,vefun(3),xxmin(3)
  real(rp)                :: rot_matrix(ndime,ndime),udotn
  real(rp)                :: bvess_new(3),xnorm,bvnew,bvold
  real(rp)                :: alpha,center_gravity(3),p_in,dummr,Q
  real(rp)                :: flow_rate_stf_value,flow_rate_tf_value, flow_rate_pf_value
  real(rp), save          :: flow_rate_pf_value_old=0.0_rp
  real(rp)                :: mean_pressure_in, mean_pressure_out, delta_P
  real(rp)                :: prfun,flowr_current,w
  logical(lg), pointer    :: lboun_nsi(:)
  logical(lg), pointer    :: lpoin_nsi(:)
  real(rp),    pointer    :: bvess_loc(:,:)
  type(func_ptr)          :: array_func(1,1)
  type(field_arrays)      :: array_fields(1)
  real(rp)                :: integrals(1,1)
  integer(ip), pointer    :: legro(:)
  real(rp)                :: flowr_imposed
  real(rp),    external   :: funcre

  nullify(lboun_nsi)
  nullify(lpoin_nsi)
  nullify(bvess_loc)
  nullify(legro)


  if( nflow_rates > 0 ) then

     call memory_alloca(mem_modul(1:2,modul),'LBOUN_NSI','nsi_updbcs',lboun_nsi,max(nboun,1_ip))
     call memory_alloca(mem_modul(1:2,modul),'LPOIN_NSI','nsi_updbcs',lpoin_nsi,max(npoin,1_ip))
     call memory_alloca(mem_modul(1:2,modul),'BVESS_LOC','nsi_updbcs',bvess_loc,ndime,npoin)

     do iflow = 1,nflow_rates

        flow_rate_stf_value = 1.0_rp
        flow_rate_tf_value  = 1.0_rp
        flow_rate_pf_value  = 1.0_rp

        if( INOTEMPTY ) then
           !
           ! LBOUN_NSI: Boundaries to consider
           !
           do iboun = 1,nboun
              if( kfl_codbo(iboun) == kfl_flow_rate_codes_nsi(iflow) ) then
                 lboun_nsi(iboun) = .true.
              else
                 lboun_nsi(iboun) = .false.
              end if
           end do
           !
           ! LPOIN_NSI: Nodes with zero velocity will be considered as fixed Dirichlet 
           !
           do ipoin = 1,npoin
              if( kfl_funno_nsi(ipoin) == 0 ) then 
                 bvess_loc(:,ipoin) = bvess_nsi(:,ipoin,2)
              else
                 bvess_loc(:,ipoin) = bvess_nsi(:,ipoin,1)
              end if
              if(  any(kfl_codno(:,ipoin) == kfl_flow_rate_codes_nsi(iflow)) ) then
                 if( sqrt(dot_product(bvess_loc(:,ipoin),bvess_loc(:,ipoin))) > zeror ) then
                    lpoin_nsi(ipoin) = .true.
                 else
                    lpoin_nsi(ipoin) = .false.
                 end if
              else
                 lpoin_nsi(ipoin) = .false.
              end if

           end do

           !
           ! Variable flow rate with space_time function
           !
           if( kfl_flow_rate_stfn_nsi(iflow) > 0 ) then
              ifunc = kfl_flow_rate_stfn_nsi(iflow)
              ipoin = 1
              call ker_space_time_function(&
                   ifunc,coord(1,ipoin),coord(2,ipoin),coord(ndime,ipoin),cutim,vefun(1:ndime))
              flow_rate_stf_value = vefun(1)
           end if

           !
           ! Variable flow rate with time function
           !
           if( kfl_flow_rate_tfn_nsi(iflow) > 0 ) then 
              ifunc = kfl_flow_rate_tfn_nsi(iflow)
              vefun(1) = funcre( &
                   time_function(ifunc) % parameters,    &
                   time_function(ifunc) % npara,         &
                   time_function(ifunc) % kfl_type,      &
                   cutim)
              flow_rate_tf_value = vefun(1)
           end if

           !
           ! Variable flow rate with pump HQ function
           !
           if( kfl_flow_rate_pfn_nsi(iflow) > 0 ) then
              ifunc = kfl_flow_rate_pfn_nsi(iflow)
              incode =  kfl_flow_rate_set_nsi(iflow)

              ! compute mean pressure on set
              mean_pressure_in  = momod(ID_NASTIN) % postp(1) % vbset(1,incode)
              mean_pressure_out = flow_rate_press_nsi(iflow)


              if(flow_rate_values_nsi(iflow) .le. 0.0_rp)then                   
                 ! Blow condition
                 delta_P = mean_pressure_in + mean_pressure_out
              elseif(flow_rate_values_nsi(iflow) .gt. 0.0_rp)then                  
                 ! Suck condition
                 delta_P = mean_pressure_out - mean_pressure_in
              endif

              call  pump_curve_flow( ifunc, delta_P, Q)


              ! Optional time relaxation to avoid oscillations
              if(flow_rate_relax_nsi(iflow).eq.1.0_rp .or. flow_rate_pf_value_old.eq.0.0_rp) then
                 w=1.0_rp
              else
                 w=flow_rate_relax_nsi(iflow)
              endif

              flow_rate_pf_value     = Q*w + flow_rate_pf_value_old*(1.0_rp-w)
              flow_rate_pf_value_old = flow_rate_pf_value
           end if

           call local_basis_local_to_global(kfl_fixrs_nsi,bvess_loc) ! Local to global                    
        end if
        !
        ! Flowr_Current conservation constraint
        !        
        if( kfl_flow_rate_alg_nsi(iflow) == 0 ) then
           call projec_mass_conservation(&
                bvess_loc,lboun_nsi,lpoin_nsi,'LOCAL MASS',&
                flow_rate_values_nsi(iflow)*flow_rate_tf_value*flow_rate_stf_value*flow_rate_pf_value,&
                'IN MY CODE WITHOUT MASTER',ERROR=ierro)
        else

           ierro = 0
           call memory_alloca(mem_modul(1:2,modul),'legro','mod_metric',legro,nboun)          
           if( nboun > 0 ) where(lboun_nsi) legro = 1

           call func_initialization(array_func)
           array_func(1,1) % f => func_dotn
           if( INOTEMPTY ) array_fields(1) % a => bvess_loc

           call integrals_boundary(array_fields,array_func,legro,integrals)
           call memory_deallo(mem_modul(1:2,modul),'legro','mod_metric',legro)

           flowr_current = integrals(1,1)
           flowr_imposed = flow_rate_values_nsi(iflow)*flow_rate_tf_value*flow_rate_stf_value*flow_rate_pf_value

           if( abs(flowr_current) > zeror ) then 
              do ipoin = 1,npoin
                 if( lpoin_nsi(ipoin) ) then
                    bvess_loc(:,ipoin) = bvess_loc(:,ipoin) * flowr_imposed / flowr_current
                 end if
              end do
           else
              do ipoin = 1,npoin
                 if( lpoin_nsi(ipoin) ) then
                    bvess_loc(:,ipoin) = 0.0_rp
                 end if
              end do              
           end if
           
        end if
        !
        ! Update solution
        !
        if( INOTEMPTY ) then
           if( ierro /= 0 ) then 
              call output_mesh_gid_format(meshe(ndivi),'FLOW_RATE_PROBLEM_'//trim(intost(iflow)),lun_outpu,lpoin_nsi)                    
              call messages_live( 'Error imposing flowrate on boundary with code '//trim(intost(kfl_flow_rate_codes_nsi(iflow))), 'WARNING' )
              if( associated(lboun_nsi) ) then 
                 if ( count(lboun_nsi,KIND=ip)==0_ip ) &
                      call messages_live( "The boundary "//trim(intost(kfl_flow_rate_codes_nsi(iflow)))&
                      //" has 0 faces. It is likely the boundary does not exist on the mesh", 'WARNING')
              end if
              if( associated(lpoin_nsi) ) then
                 if ( count(lpoin_nsi,KIND=ip)==0_ip ) &
                      call messages_live( "The boundary "//trim(intost(kfl_flow_rate_codes_nsi(iflow)))&
                      //" has 0 nodes. It is likely the boundary does not exist on the mesh", 'WARNING')
              end if

              write(lun_outpu,*) 'ERROR IN FLOW RATE'
              write(lun_outpu,*) 'iflow=',iflow, ' code=', kfl_flow_rate_codes_nsi(iflow),' flowrate=',&
                   flow_rate_values_nsi(iflow), ' flow_rate_stf_value=',flow_rate_stf_value
              if( associated(lboun_nsi) ) write(lun_outpu,*) 'LBOUN= ',count(lboun_nsi,KIND=ip)
              if( associated(lpoin_nsi) ) write(lun_outpu,*) 'LPOIN= ',count(lpoin_nsi,KIND=ip)
           end if

           call local_basis_global_to_local(kfl_fixrs_nsi,bvess_loc) ! Global to local

           do ipoin = 1,npoin
              if( any(kfl_codno(:,ipoin) == kfl_flow_rate_codes_nsi(iflow)) ) then
                 bvess_nsi(:,ipoin,1) = bvess_loc(:,ipoin)
              end if
           end do

        end if

        if( ierro /= 0 ) call runend('NSI_UPDBCS: INCORRECT FLOWRATE VALUE OR FLOWRATE IMPOSED ON NONEXISTENT BOUNDARY CODE')

     end do

     call memory_deallo(mem_modul(1:2,modul),'BVESS_LOC','nsi_updbcs',bvess_loc)
     call memory_deallo(mem_modul(1:2,modul),'LBOUN_NSI','nsi_updbcs',lboun_nsi)
     call memory_deallo(mem_modul(1:2,modul),'LPOIN_NSI','nsi_updbcs',lpoin_nsi)             

  end if

end subroutine nsi_impose_flow_rates
