!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Nastin
!> @{
!> @file    mod_nsi_arrays.f90
!> @author  houzeaux
!> @date    2019-11-16
!> @brief   Nastin arrays
!> @details Nastin arrays
!-----------------------------------------------------------------------

module mod_nsi_arrays

  use def_master
  use def_domain 
  use def_nastin
  use def_kermod
  use mod_output_postprocess,  only : output_postprocess_check_variable_postprocess
  use mod_arrays,              only : arrays
  use mod_memory,              only : memory_alloca
  use mod_arrays,              only : arrays_number

  
  implicit none

  private

  public :: nsi_arrays
  
contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-11-16
  !> @brief   Nastin arrays
  !> @details Do what you have to do with nastin arrays
  !> 
  !-----------------------------------------------------------------------
  
  subroutine nsi_arrays(wtask)

    character(len=*), intent(in) :: wtask
    integer(ip)                  :: ncomp_loc,ielem
    integer(ip)                  :: pgaus,ncsgs,pelty

    if (NSI_FRACTIONAL_STEP.and.kfl_stabi_nsi == NSI_ASGS) then
       ncomp_loc = ncomp_nsi+1
    else
       ncomp_loc = ncomp_nsi
    end if
    !
    ! Unknowns
    !
    call arrays(arrays_number('VELOC'),trim(wtask),veloc,ndime,npoin,ncomp_loc)
    call arrays(arrays_number('PRESS'),trim(wtask),press,npoin,ncomp_nsi)
    !
    ! DENSI: Compressible regime
    !
    if( kfl_regim_nsi == 1 ) then
       call arrays(arrays_number('DENSI'),trim(wtask),densi,npoin,1_ip)
    else if( kfl_regim_nsi == 2 ) then
       call arrays(arrays_number('DENSI'),trim(wtask),densi,npoin,ncomp_nsi)
    end if
    !
    ! VEOLD_NSI: Allocate memory for old (last iteration) velocity
    !
    if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='RESID') ) then
       kfl_resid_nsi=1
       call arrays(arrays_number('VEOLD'),trim(wtask),veold_nsi,ndime,npoin)        
    else
       kfl_resid_nsi=0
    end if
    !
    ! VEPRO_NSI, PRPRO_NSI, GRPRO_NSI: Projections for orthogonal SGS
    ! 
    if( kfl_stabi_nsi > 0 ) then
       call arrays(arrays_number('VEPRO'),trim(wtask),vepro_nsi,ndime,npoin)        
       call arrays(arrays_number('PRPRO'),trim(wtask),prpro_nsi,npoin)  
       if( kfl_stabi_nsi == 2 ) then
          call arrays(arrays_number('GRPRO'),trim(wtask),grpro_nsi,ndime,npoin)
       end if
    end if
    !
    ! Immersed boundary method
    !
    if(      kfl_immer_nsi == 1 ) then
       call arrays(arrays_number('LAGRA'),trim(wtask),lagra_nsi,ndime,npoin,1_ip)
    else if( kfl_immer_nsi == 2 ) then
       call arrays(arrays_number('TAUIB'),trim(wtask),tauib_nsi,ndime,ndime,npoin)
    end if
    !
    ! VESGS: Subgrid scale velocity
    !
    if( kfl_sgsco_nsi == 1 .or. kfl_sgsti_nsi == 1 ) then
       call arrays(arrays_number('VESGS'),trim(wtask),vesgs,nelem)
       if( trim(wtask) == 'ALLOCATE' ) then
          ncsgs = min(2_ip,2_ip*kfl_sgsti_nsi+kfl_sgsco_nsi)
          do ielem = 1,nelem
             pelty = abs(ltype(ielem))
             pgaus = ngaus(pelty)
             call memory_alloca(mem_modul(1:2,modul),'VESGS % A','nsi_memall',vesgs(ielem)%a,ndime,pgaus,ncsgs)
          end do
       end if
    end if
    !
    ! DUNKN_NSI, DUNKP_NSI: Delta velocity and pressure for Aitken relaxation strategy
    !
    if(kfl_relax_nsi==2) then
       call arrays(arrays_number('DUNKN'),trim(wtask),dunkn_nsi,ndime,npoin)
    end if
    if(kfl_relap_nsi==2) then
       call arrays(arrays_number('DUNKP'),trim(wtask),dunkp_nsi,npoin)
    end if
    !
    ! MASS_RHO_NSI
    !
    if( kfl_predi_nsi == 7 .or. kfl_predi_nsi == 8 .or. kfl_predi_nsi == 9 ) then
       call arrays(arrays_number('MASRH'),trim(wtask),mass_rho_nsi,npoin,ncomp_nsi)
    end if
    !
    ! arrays for no slip wall law - later they can be used for other cases too
    ! I have intialized tehm to 0 here but ther better option might be possible
    !
    if( kfl_noslw_ker /= 0                                                   .or. &
        output_postprocess_check_variable_postprocess(VARIABLE_NAME='VAFOR') .or. &
        (NSI_FRACTIONAL_STEP .and. kfl_tisch_nsi == 4 )                           ) then
       call arrays(arrays_number('VAFOR'),trim(wtask),vafor_nsi,ndime,npoin)       
    end if
    if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVVAF') .or. kfl_noslw_ker /= 0 )  then
       call arrays(arrays_number('AVVAF'),trim(wtask),avvaf_nsi,ndime,npoin)       
    end if
    if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVNTR') .or. kfl_noslw_ker /= 0 )  then
       call arrays(arrays_number('AVNTR'),trim(wtask),avntr_nsi,ndime,npoin)       
    end if
    if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVGTR') .or. kfl_noslw_ker /= 0 )  then
       call arrays(arrays_number('AVGTR'),trim(wtask),avgtr_nsi,ndime,npoin)       
    end if
    !
    ! Coupling with SOLIDZ
    !
    if( coupling('SOLIDZ','NASTIN') >= 1 .or. coupling('NASTIN','IMMBOU') >= 1 ) then
       call arrays(arrays_number('FORCF'),trim(wtask),forcf,ndime,npoin)       
    end if
    !
    ! Traction on boundary nodes for alternative AVTAN postprocess
    !
    if(    kfl_noslw_ker /= 0                                                   .or. &
         & output_postprocess_check_variable_postprocess(VARIABLE_NAME='NOTRA') .or. &
         & output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVTAN') .or. &
         & output_postprocess_check_variable_postprocess(VARIABLE_NAME='ENTAN') .or. &
         & output_postprocess_check_variable_postprocess(VARIABLE_NAME='TANGE') .or. &
         & output_postprocess_check_variable_postprocess(VARIABLE_NAME='FTANG')      )  then
       call arrays(arrays_number('NOTRA'),trim(wtask),notra_nsi,ndime,npoin)       
    end if
    !
    ! AVVEL_NSI: average velocity (postprocess)
    !
    if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVVEL') ) then
       call arrays(arrays_number('AVVEL'),trim(wtask),avvel_nsi,ndime,npoin)       
    end if
    !
    ! AVPRE_NSI: average pressure (postprocess)
    !
    if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVPRE') ) then
       call arrays(arrays_number('AVPRE'),trim(wtask),avpre_nsi,npoin)       
    end if
    !
    ! AVVE2_NSI: average velocity**2 (postprocess)
    !
    if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVVE2') ) then
       call arrays(arrays_number('AVVE2'),trim(wtask),avve2_nsi,ndime,npoin)       
    end if
    !
    ! AVVXY_NSI: average vx*vy (postprocess)
    !
    if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVVXY') ) then
       call arrays(arrays_number('AVVXY'),trim(wtask),avvxy_nsi,ndime,npoin)       
    end if
    !
    ! AVPR2_NSI: average pressure (postprocess)
    !
    if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVPR2') ) then
       call arrays(arrays_number('AVPR2'),trim(wtask),avpr2_nsi,npoin)       
    end if
    !
    ! AVTAN_NSI: average TANGE (postprocess)
    !
    if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVTAN') ) then
       call arrays(arrays_number('AVTAN'),trim(wtask),avtan_nsi,ndime,npoin)       
    end if
    !
    ! AVMUT_NSI: Average turbulent viscosity (postprocess)
    !
    if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVMUT') ) then
       call arrays(arrays_number('AVMUT'),trim(wtask),avmut_nsi,npoin)       
    end if
    !
    ! AVSTX_NSI: Average stress mu_t * grad(u)
    !
    if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVSTX') ) then
       call arrays(arrays_number('AVSTX'),trim(wtask),avstx_nsi,ndime,npoin)       
    end if
    !
    ! AVSTY_NSI: Average stress mu_t * grad(v)
    !
    if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVSTY') ) then
       call arrays(arrays_number('AVSTY'),trim(wtask),avsty_nsi,ndime,npoin)       
    end if
    !
    ! AVSTZ_NSI: Average stress mu_t * grad(w)
    !
    if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVSTZ') ) then
       call arrays(arrays_number('AVSTZ'),trim(wtask),avstz_nsi,ndime,npoin)       
    end if
    !
    ! AVMOS_NSI: average momentum source from spray
    !
    if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVMOS') ) then
       call arrays(arrays_number('AVMOS'),trim(wtask),avmos_nsi,ndime,npoin)       
    end if
    !
    ! AVMFL_NSI: average mass flux
    !
    if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVMFL') ) then
       call arrays(arrays_number('AVMFL'),trim(wtask),av_mass_flux_nsi,ndime,npoin)       
    end if
    !
    ! AVRUU_NSI: average momentum flux, diagonal terms
    !
    if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVRUU') ) then
       call arrays(arrays_number('AVRUU'),trim(wtask),av_mom_flux_diag_nsi,ndime,npoin)       
    end if
    !
    ! AVRUV_NSI: average momentum flux, off-diagonal terms
    !
    if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVRUV') ) then
       call arrays(arrays_number('AVRUV'),trim(wtask),av_mom_flux_off_nsi,ndime,npoin)       
    end if
    !
    ! VMAXP_NSI: Maximum nodewise velocity
    !
    if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='VMAXP') .or. &
        output_postprocess_check_variable_postprocess(VARIABLE_NAME='PINDE') ) then
       call arrays(arrays_number('VMAXP'),trim(wtask),vmaxp_nsi,ndime,npoin)       
        !vmaxp_nsi=0.0_rp
    end if
    !
    ! VMINP_NSI: Minumum nodewise velocity
    !
    if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='VMINP') .or. &
        output_postprocess_check_variable_postprocess(VARIABLE_NAME='PINDE') ) then
       call arrays(arrays_number('VMINP'),trim(wtask),vminp_nsi,ndime,npoin)       
       !vminp_nsi=huge(1.0_rp)
    end if
    !
    ! VAVEP_NSI: Average nodewise velocity
    !
    if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='VAVEP') .or. &
        output_postprocess_check_variable_postprocess(VARIABLE_NAME='PINDE') ) then
       call arrays(arrays_number('VAVEP'),trim(wtask),vavep_nsi,ndime,npoin)       
    end if
    !
    ! PINDE_NSI: Nodewise pulsatility index
    !
    if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='PINDE') ) then
       call arrays(arrays_number('PINDE'),trim(wtask),pinde_nsi,1_ip,npoin)       
    end if
    !
    ! ENVEL_NSI: ensemble velocity (postprocess)
    !
    if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='ENVEL') ) then
       call arrays(arrays_number('ENVEL'),trim(wtask),envel_nsi,ndime,npoin)       
    end if
    !
    ! ENPRE_NSI: ensemble pressure (postprocess)
    !
    if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='ENPRE') ) then
       call arrays(arrays_number('ENPRE'),trim(wtask),enpre_nsi,npoin)       
    end if
    !
    ! ENVE2_NSI: ensemble velocity**2 (postprocess)
    !
    if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='ENVE2') ) then
       call arrays(arrays_number('ENVE2'),trim(wtask),enve2_nsi,ndime,npoin)       
    end if
    !
    ! ENVXY_NSI: ensemble vx*vy (postprocess)
    !
    if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='ENVXY') ) then
       call arrays(arrays_number('ENVXY'),trim(wtask),envxy_nsi,ndime,npoin)       
    end if
    !
    ! ENPRE_NSI: ensemble pressure (postprocess)
    !
    if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='ENPR2') ) then
       call arrays(arrays_number('ENPR2'),trim(wtask),enpr2_nsi,npoin)       
    end if
    !
    ! ENTAN_NSI: ensemble TANGE (postprocess)
    !
    if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='ENTAN') ) then
       call arrays(arrays_number('ENTAN'),trim(wtask),entan_nsi,ndime,npoin)       
    end if
    !
    ! ENMUT_NSI: Ensemble turbulent viscosity (postprocess)
    !
    if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='ENMUT') ) then
       call arrays(arrays_number('ENMUT'),trim(wtask),enmut_nsi,npoin)       
    end if
    !
    ! ENSTX_NSI: Ensemble stress mu_t * grad(u)
    !
    if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='ENSTX') ) then
       call arrays(arrays_number('ENSTX'),trim(wtask),enstx_nsi,ndime,npoin)       
    end if
    !
    ! ENSTY_NSI: Ensemble stress mu_t * grad(v)
    !
    if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='ENSTY') ) then
       call arrays(arrays_number('ENSTY'),trim(wtask),ensty_nsi,ndime,npoin)       
    end if
    !
    ! ENSTZ_NSI: Ensemble stress mu_t * grad(w)
    !
    if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='ENSTZ') ) then
       call arrays(arrays_number('ENSTZ'),trim(wtask),enstz_nsi,ndime,npoin)       
    end if
    !
    ! Bubble
    !
    if( kfl_bubbl_nsi /= 0 ) then
       call arrays(arrays_number('BUBBL'),trim(wtask),bubble_nsi    ,nelem)            
       call arrays(arrays_number('BUAQQ'),trim(wtask),bubble_aqq_nsi,nelem)            
       call arrays(arrays_number('BUAQU'),trim(wtask),bubble_aqu_nsi,mnode*ndime,nelem)
       call arrays(arrays_number('BUAQP'),trim(wtask),bubble_aqp_nsi,mnode,nelem)      
       call arrays(arrays_number('BUBQ '),trim(wtask),bubble_bq_nsi ,nelem)                    
    end if
    !
    ! RANS/LES two-layer model
    !
    if( kfl_twola_ker > 0 ) then
       call arrays(arrays_number('BTRAC'),trim(wtask),btrac_nsi,ndime,npoin)                    
       call arrays(arrays_number('TRACR'),trim(wtask),tracr_nsi,ndime,npoin)  
       call arrays(arrays_number('TLUAV'),trim(wtask),tluav_nsi,ndime,npoin) 
    end if

  end subroutine nsi_arrays
   
end module mod_nsi_arrays
!> @}
