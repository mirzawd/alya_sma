!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!Sets vconc(1,ipoin,1) -- calcium
!     vconc(2,ipoin,1) -- reference potential specific to Fitzhugh-Nagumo model
module mod_exm_fitzhugh_nagumo
    use def_master
    use def_domain
    use def_exmedi
    use mod_messages,   only : messages_live

    implicit none

    logical(lg), protected :: EXM_FITZHUGH_NAGUMO_MODEL = .FALSE.         ! use to check if the model is used for atleast one material
    public :: EXM_FITZHUGH_NAGUMO_MODEL

    private

    integer(ip), parameter :: nparams = 13_ip                                  ! number of model paramters per material
    integer(ip), parameter :: &
        FHN_PAR_A           = 1_ip, &
        FHN_PAR_B           = 2_ip, &
        FHN_PAR_C1          = 3_ip, &
        FHN_PAR_C2          = 4_ip, &
        FHN_PAR_D           = 5_ip, &
        FHN_PAR_TS          = 6_ip, &
        FHN_PAR_VMIN        = 7_ip, &
        FHN_PAR_VMAX        = 8_ip, &
        FHN_PAR_TAU         = 9_ip, &
        FHN_PAR_CAM         = 10_ip, & 
        FHN_PAR_CA0         = 11_ip, & 
        FHN_PAR_CA_DELAY    = 12_ip, & 
        FHN_PAR_CA_TRESHOLD = 13_ip 
        
    integer(ip), parameter :: vconc_nvars = 4_ip  !number of variables below (used from vconc)
    integer(ip), parameter :: &
        CALCIUM_VCONCID = 1_ip, &
        REF_POTENTIAL_VCONCID = 2_ip, &                     ! which element of vconc to use for the reference potential
        ACT_TIME_VCONCID = 3_ip, &
        U_PREV_VCONCID = 4_ip
    
    real(rp), pointer :: params_fhn_exm(:,:)   ! Physical model parameters, for now only Fitzhugh-Nagumo model is using it


    public :: exm_fitzhugh_nagumo_ReadDat
    public :: exm_fitzhugh_nagumo_GetNumberOfVariables
    public :: exm_fitzhugh_nagumo_NodalIonicCurrents
    public :: exm_fitzhugh_nagumo_SendData
    public :: exm_fitzhugh_nagumo_InitVariables
    public :: exm_fitzhugh_nagumo_Initialize
    public :: exm_fitzhugh_nagumo_WriteInitialValuesToLog
    public :: exm_elmagToFitzhugh_nagumo
    public :: exm_elmagFromFitzhugh_nagumo
    public :: exm_fitzhugh_nagumo_GetDT
    public :: exm_fitzhugh_nagumo_GetCa0
    
contains

function exm_fitzhugh_nagumo_GetDT(imate)
    implicit none
    integer(ip), intent(in)  :: imate
    real(rp) :: exm_fitzhugh_nagumo_GetDT

    exm_fitzhugh_nagumo_GetDT = dtime * 1000.0_rp * 1.0_rp / params_fhn_exm( FHN_PAR_TS,  imate )  
end function

function exm_elmagToFitzhugh_nagumo(elmag, imate)
    implicit none
    real(rp), intent(in)     :: elmag
    integer(ip), intent(in)  :: imate
    real(rp)    :: exm_elmagToFitzhugh_nagumo
    real(rp)    :: Vref_min, Vref_max
    
    Vref_min = params_fhn_exm(FHN_PAR_VMIN, imate) ! Min ref potential
    Vref_max = params_fhn_exm(FHN_PAR_VMAX, imate) ! Max ref potential

    exm_elmagToFitzhugh_nagumo = (elmag - Vref_min) / (Vref_max - Vref_min)
end function

function exm_elmagFromFitzhugh_nagumo(elmag_fhn, imate)
    implicit none
    real(rp), intent(in)     :: elmag_fhn
    integer(ip), intent(in)  :: imate
    real(rp)    :: exm_elmagFromFitzhugh_nagumo
    real(rp)    :: Vref_min, Vref_max
    
    Vref_min = params_fhn_exm(FHN_PAR_VMIN, imate) ! Min ref potential
    Vref_max = params_fhn_exm(FHN_PAR_VMAX, imate) ! Max ref potential

    exm_elmagFromFitzhugh_nagumo = elmag_fhn * (Vref_max - Vref_min) + Vref_min
end function




subroutine exm_fitzhugh_nagumo_WriteInitialValuesToLog()
    use mod_outfor
    implicit none


    integer(ip) :: imate

    if ( EXM_FITZHUGH_NAGUMO_MODEL ) then
        call outfor(25_ip,momod(modul) % lun_outpu,'FITZHUGH-NAGUMO PARAMETERS')
        write(momod(modul)%lun_outpu,*) NEW_LINE('a')

        do imate=1,nmate
            write(momod(modul)%lun_outpu,*) '   Material '//trim(intost(imate))
            write(momod(modul)%lun_outpu,*) '      a            = ', params_fhn_exm(FHN_PAR_A,           imate)
            write(momod(modul)%lun_outpu,*) '      b            = ', params_fhn_exm(FHN_PAR_B,           imate) 
            write(momod(modul)%lun_outpu,*) '      c1           = ', params_fhn_exm(FHN_PAR_C1,          imate) 
            write(momod(modul)%lun_outpu,*) '      c2           = ', params_fhn_exm(FHN_PAR_C2,          imate) 
            write(momod(modul)%lun_outpu,*) '      d            = ', params_fhn_exm(FHN_PAR_D,           imate) 
            write(momod(modul)%lun_outpu,*) '      time_scale   = ', params_fhn_exm(FHN_PAR_TS,          imate) 
            write(momod(modul)%lun_outpu,*) '      Vmin         = ', params_fhn_exm(FHN_PAR_VMIN,        imate) 
            write(momod(modul)%lun_outpu,*) '      Vmax         = ', params_fhn_exm(FHN_PAR_VMAX,        imate) 
            write(momod(modul)%lun_outpu,*) '      Tau          = ', params_fhn_exm(FHN_PAR_TAU,         imate) 
            write(momod(modul)%lun_outpu,*) '      cam          = ', params_fhn_exm(FHN_PAR_CAM,         imate) / sms_conversion_currents_exm
            write(momod(modul)%lun_outpu,*) '      Ca0          = ', params_fhn_exm(FHN_PAR_CA0,         imate) / sms_conversion_currents_exm
            write(momod(modul)%lun_outpu,*) '      Ca_delay     = ', params_fhn_exm(FHN_PAR_CA_DELAY,    imate) 
            write(momod(modul)%lun_outpu,*) '      Ca_threshold = ', params_fhn_exm(FHN_PAR_CA_TRESHOLD, imate)    
        end do

        write(momod(modul)%lun_outpu,*) NEW_LINE('a')
    end if

end subroutine exm_fitzhugh_nagumo_WriteInitialValuesToLog


subroutine exm_fitzhugh_nagumo_Initialize(ipoin, imate)
    implicit none
    integer(ip), intent(in) :: ipoin, imate

    elmag(ipoin,:) = params_fhn_exm(FHN_PAR_VMIN, imate)
    vconc(U_PREV_VCONCID, ipoin, : ) = 0.0_rp
    vconc(ACT_TIME_VCONCID, ipoin, :) = -1.0_rp !calcium activation hsa not happened yet
    vconc(CALCIUM_VCONCID, ipoin, :) = exm_fitzhugh_nagumo_GetCa0( imate )
    vconc(REF_POTENTIAL_VCONCID, ipoin, :) = 0.0_rp
end subroutine exm_fitzhugh_nagumo_Initialize


subroutine exm_fitzhugh_nagumo_InitVariables()
    use mod_memory, only : memory_alloca
    implicit none
    integer(ip) :: i,j

    call memory_alloca(mem_modul(1:2,modul),'params_fhn_exm', 'exm_fitzhugh_nagumo', params_fhn_exm, nparams,  nmate)

    do i=1,size(params_fhn_exm,1,kind=ip)
        do j=1,size(params_fhn_exm,2,kind=ip)
            params_fhn_exm(i,j) = 0.0_rp
        end do
    end do

    params_fhn_exm(FHN_PAR_A,    :) = 0.13_rp  ! original values of paper: a=0.13;
    params_fhn_exm(FHN_PAR_B,    :) = 0.013_rp   !0.0035_rp ! original values of paper: b=0.013;
    params_fhn_exm(FHN_PAR_C1,   :) = 0.26_rp  ! original values of paper: c1=0.26;
    params_fhn_exm(FHN_PAR_C2,   :) = 0.1_rp   ! original values of paper: c2=0.1;
    params_fhn_exm(FHN_PAR_D,    :) = 1.0_rp    ! original values of paper: d=1.0;
    params_fhn_exm(FHN_PAR_TS,   :) = 1.0_rp    ! "cycle length" to scale APD
    params_fhn_exm(FHN_PAR_VMIN, :) = -86.0_rp  ! Min ref potential
    params_fhn_exm(FHN_PAR_VMAX, :) =  35.0_rp  ! Max ref potential
    params_fhn_exm(FHN_PAR_TAU, :)  =  0.06_rp  ! % for calcium calculation sec
    params_fhn_exm(FHN_PAR_CAM, :)  =  0.001_rp  ! % for calcium calculation,microM, units of the paper
    params_fhn_exm(FHN_PAR_CA0, :)  =  0.00001_rp  ! % for calcium calculation,microM, units of the paper
    params_fhn_exm(FHN_PAR_CA_DELAY, :)  =  0.008_rp  ! % for calcium calculation,s. Calcium onset delay from activation
    params_fhn_exm(FHN_PAR_CA_TRESHOLD, :)  =  0.2_rp  ! % for calcium calculation. Threshold of u for deciding if to release calcium

end subroutine exm_fitzhugh_nagumo_InitVariables

function exm_fitzhugh_nagumo_GetCa0(imate)
    implicit none
    integer(ip) :: imate
    real(rp)    :: exm_fitzhugh_nagumo_GetCa0

    exm_fitzhugh_nagumo_GetCa0 = params_fhn_exm(FHN_PAR_CA0, imate)

end function exm_fitzhugh_nagumo_GetCa0

subroutine exm_fitzhugh_nagumo_SendData()
    use      mod_exchange,            only : exchange_add, exchange_end, exchange_init

    implicit none

    call exchange_init()
    call exchange_add(params_fhn_exm)
    call exchange_add(EXM_FITZHUGH_NAGUMO_MODEL)
    call exchange_end()

end subroutine exm_fitzhugh_nagumo_SendData




subroutine exm_fitzhugh_nagumo_ReadDat(imate)
    use mod_ecoute,              only :  ecoute
    use def_inpout

    implicit none
    integer(ip), intent(in) :: imate 

    EXM_FITZHUGH_NAGUMO_MODEL = .TRUE.

    !default initial values            
    !    ALPHA= 0.1_rp   $ old alpha
    !    BETA = 0.01_rp  $ old epsilon
    !    DELTA= 0.5_rp   $ old gamma_membrane
    !    C1FHN= 10.0_rp  $ old c_membrane_const
    !    C2FHN= 0.5_rp   $ no estaba

    do while(words(1)/='ENDCE')

        if(words(1)=='CMCON') then               ! Model constant Cm
            xmccmmate_exm(imate)=param(1)                          
        end if

        if(words(1)==     'A    ') then               ! Model constant 
            params_fhn_exm( FHN_PAR_A, imate)=param(1)
        else if(words(1)=='B    ') then               ! Model constant 
            params_fhn_exm( FHN_PAR_B, imate)=param(1)
        else if(words(1)=='C1   ') then               ! Model constant 
            params_fhn_exm( FHN_PAR_C1, imate)=param(1)
        else if(words(1)=='C2   ') then               ! Model constant 
            params_fhn_exm( FHN_PAR_C2, imate)=param(1)
        else if(words(1)=='D    ') then               ! Model constant 
            params_fhn_exm( FHN_PAR_D, imate)=param(1)                            
        else if(words(1)=='TIMES') then               ! Model constant 
            params_fhn_exm( FHN_PAR_TS, imate)=param(1)                             
        else if(words(1)=='VMIN ') then
            params_fhn_exm( FHN_PAR_VMIN, imate)=param(1)                             
        else if(words(1)=='VMAX ') then
            params_fhn_exm( FHN_PAR_VMAX, imate)=param(1)                             
        else if(words(1)=='TAU  ') then
            params_fhn_exm( FHN_PAR_TAU, imate)=param(1)                             
        else if(words(1)=='CAM  ') then
            params_fhn_exm( FHN_PAR_CAM, imate)=param(1) * sms_conversion_currents_exm                            
        else if(words(1)=='CA0  ') then
            params_fhn_exm( FHN_PAR_CA0, imate)=param(1) * sms_conversion_currents_exm
        else if(words(1)=='CADEL') then
            params_fhn_exm( FHN_PAR_CA_DELAY, imate)=param(1)                             
        else if(words(1)=='CATHR') then
            params_fhn_exm( FHN_PAR_CA_TRESHOLD, imate)=param(1)                             
        end if    

        call ecoute('exm_reaphy')
    end do
end subroutine exm_fitzhugh_nagumo_ReadDat


subroutine exm_fitzhugh_nagumo_GetNumberOfVariables(imate, nvars_vauxi, nvars_vconc)
    implicit none

    integer(ip), intent(in)  :: imate
    integer(ip), intent(out) :: nvars_vauxi, nvars_vconc

    nvars_vauxi = 0_ip  !  number of variables used for activation/inactivation 
    nvars_vconc = vconc_nvars  !  number of variables used for concentrations: Ca, ref. potential

end subroutine exm_fitzhugh_nagumo_GetNumberOfVariables



!------------------------------------------------------------------------
  !> @addtogroup ExmediIonicurrents
  !> @{
  !> @file    exm_mdfhn_ionicurrents.f90
  !> @date    17/10/2016
  !> @author  Jazmin A-S
  !> @brief   Computes the modified Fitzhugh-Nagumo EP model
  !> @details Implementation from Rogers & McCulloch 1994
  !> @}
  !
  ! Model solved in this subroutine:
  ! (see that the last term contains u, this is a modification suggested in R&M 1994)
  !                   Iion= C_1*u(u-A)(1-u)-C_2*v*u 
  !                   dv/dt=  B*(u-Dv)   
  ! 
  ! This subroutine computes Iion, which can be then added to the original equation
  ! 
  !                   du/dt= diffusion terms + (Istim+Iion)/Cm
  ! Where:
  !       u           => voltage (elmag(:,ITER_K))
  !       v           => recovery variable (V_fhn_exm)
  !       afhn=0.13_rp  => original values of paper: a=0.13;    old alpha
  !       bfhn=0.0035_rp => original values of paper: b=0.013;  old epsilon
  !       c1fhn=0.26_rp  => original values of paper: c1=0.26;  old cmemb
  !       c2fhn=0.1_rp   => original values of paper: c2=0.1;   
  !       dfhn=1.0_rp    => original values of paper: d=1.0;    old gamma
  !------------------------------------------------------------------------
subroutine exm_fitzhugh_nagumo_NodalIonicCurrents(ipoin,xioni,dioni,cai)
    use      def_master
    use      def_domain
    use      def_elmtyp
    use      def_exmedi

  
    implicit none
    integer(ip), intent(in) :: ipoin
    real(rp), intent(out)   :: xioni
    real(rp), intent(out)   :: dioni
    real(rp), intent(out)   :: cai
    integer(ip)             :: imate    
  
    real(rp)   :: a, & ! A
         b, & ! B
         c1, & ! C1
         c2, & ! C2
         d,  & !  D
!         dv, & !  D
         Vref_max, Vref_min, dt, u_old, u_new, v_old, v_new, thau, cam, ca0, t, t_activated, ca_delay, u_prev
    
  
    imate= nodemat(ipoin)
  
    a     = params_fhn_exm( FHN_PAR_A,     imate)           
    b     = params_fhn_exm( FHN_PAR_B,     imate)
    c1    = params_fhn_exm( FHN_PAR_C1,    imate)           
    c2    = params_fhn_exm( FHN_PAR_C2,    imate)           
    d     = params_fhn_exm( FHN_PAR_D,     imate)         
    Vref_min = params_fhn_exm( FHN_PAR_VMIN,  imate)                   
    Vref_max = params_fhn_exm( FHN_PAR_VMAX,  imate)    

    dt =  exm_fitzhugh_nagumo_GetDT(imate) 
    
      
    u_old = (elmag(ipoin,TIME_N  ) - Vref_min) / (Vref_max - Vref_min)
    v_old = vconc(REF_POTENTIAL_VCONCID, ipoin, 2_ip)


    ! solve the refhn time equation ! dv/dt = b*(u - d*v)  
    v_new = v_old + dt * b * ( u_old - d * v_old )  
    vconc(REF_POTENTIAL_VCONCID, ipoin, 1_ip) = v_new
  
    dioni = 0.0_rp
    u_new = ( c1 * u_old * ( u_old-a )*(1.0_rp - u_old) - c2*u_old*v_old )
    xioni = -u_new*(Vref_max-Vref_min) ! the minus is needed to cancel out the "-" in mod_exm_ionicurrents
      
    
    !---------------------------------------
    !
    !  Decide if activation happened
    !
    ! If u is > threshold and it is an upstroke, start releasing calcium
    u_prev = vconc(U_PREV_VCONCID, ipoin, 1_ip )
    if ( (u_prev < params_fhn_exm( FHN_PAR_CA_TRESHOLD, imate) ) .AND. &
        (u_old >= params_fhn_exm( FHN_PAR_CA_TRESHOLD, imate) ) .AND. &
        (u_old > u_prev) ) then 
        !print *,'activating'
        vconc(ACT_TIME_VCONCID, ipoin, ITER_K) = cutim
    end if
    !----------------------------------------
    !
    !   calculate calcium
    !
    thau = params_fhn_exm( FHN_PAR_TAU, imate)
    cam = params_fhn_exm( FHN_PAR_CAM, imate)
    ca0 = params_fhn_exm( FHN_PAR_CA0, imate)
    ca_delay = params_fhn_exm( FHN_PAR_CA_DELAY, imate)
    t_activated = vconc(ACT_TIME_VCONCID, ipoin, ITER_K)

    if (t_activated>=0.0_rp) then !initially vconc(ACT_TIME_VCONCID, ipoin, ITER_K) is -1.0_rp
        t = cutim - t_activated 

        if ( t-ca_delay >= 0.0_rp ) then
            cai = ( ca0 + (cam-ca0) * ( (t-ca_delay)/thau )*exp(1.0_rp-((t-ca_delay)/thau)) ) ! convert to exmedi units
        else
            cai = ca0
        end if
    else 
        cai = ca0
    end if
    vconc(CALCIUM_VCONCID,ipoin,1_ip) = cai

    !if (ipoin==1) then
    !    write(2000,*) cutim, ipoin, u_old, v_old, cai, elmag(ipoin,1)
    !end if

    vconc(U_PREV_VCONCID, ipoin, 1_ip ) = u_old

  end subroutine exm_fitzhugh_nagumo_NodalIonicCurrents


end module mod_exm_fitzhugh_nagumo
