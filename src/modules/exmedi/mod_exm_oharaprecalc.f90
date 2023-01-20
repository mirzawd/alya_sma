!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



module mod_exm_oharaprecalc

   use def_kintyp, only : ip, rp, lg
   use def_master, only : intost
   use mod_eccoupling, only : EXMSLD_CELL_OHARA, EXMSLD_CELL_OHARA_INAPA

   implicit none

   
   private


   integer(ip), public, parameter :: &
      OHR_M_INF             =1,  OHR_H_INF       =2,  OHR_H_CAMK_INF=3,  OHR_M_L_INF    =4,  &
      OHR_H_L_INF           =5,  OHR_H_L_CAMK_INF=6,  OHR_A_INF     =7,  OHR_INV_TAU_A_1=8,  &
      OHR_INV_TAU_A_2       =9,  OHR_I_INF       =10, OHR_A_I_FAST  =11, OHR_A_CAMK_INF =12, &
      OHR_DELTA_CAMK_RECOVER=13, OHR_D_INF       =14, OHR_F_INF     =15, OHR_A_F_CA_FAST=16, &
      OHR_X_R_INF           =17, OHR_A_XR_FAST   =18, OHR_R_KR_1    =19, OHR_R_KR_2     =20, &
      OHR_X_S1_INF          =21, OHR_X_KB        =22, OHR_DELTA_EPI =23
   integer(ip), public, parameter :: EXM_OHR_INF_EQUATIONS = 23

   integer(ip), public, parameter :: &
      OHR_TAU_M      =1,  OHR_TAU_H_FAST =2,  OHR_TAU_H_SLOW        =3,  OHR_TAU_J        =4, &
      OHR_TAU_I_FAST =5,  OHR_TAU_I_SLOW =6,  OHR_DELTA_CAMK_DEVELOP=7,  OHR_TAU_D        =8, &
      OHR_TAU_F_FAST =9,  OHR_TAU_F_SLOW =10, OHR_TAU_F_CA_FAST     =11, OHR_TAU_F_CA_SLOW=12, &
      OHR_TAU_XR_FAST=13, OHR_TAU_XR_SLOW=14, OHR_TAU_X_S1          =15, OHR_TAU_X_S2     =16, &
      OHR_TAU_X_K1   =17
   integer(ip), public, parameter :: EXM_OHR_TAU_EQUATIONS = 17




   TYPE, public :: ohara_precalc
      private
      logical(lg)                     :: modified
      integer(ip)                     :: use_safe_exp !1 - safe exponent, 0 - normal exponent
      real(rp)                        :: v 
      real(rp), dimension(:), pointer :: inf_storage => null()
      real(rp), dimension(:), pointer :: tau_storage => null()
      integer(ip)                     :: model_type
      real(rp)                        :: tol  !tolerance to decide if recalculation is needed (only is |v_old-v_new|>=tol)

   contains
      private
      procedure, public   :: calculate        => ohara_precalc_calculate
      procedure           :: calculate_infs   => ohara_precalc_calculate_infs
      procedure           :: calculate_taus   => ohara_precalc_calculate_taus
      procedure, public   :: init             => ohara_precalc_init
      procedure, public   :: calculate_gates  => ohara_precalc_calculate_gates
      procedure, public   :: set_voltage      => ohara_precalc_set_voltage
      procedure, public   :: get_infs         => ohara_precalc_get_infs
      procedure, public   :: get_taus         => ohara_precalc_get_taus
      procedure, public   :: get_inf          => ohara_precalc_get_inf
      procedure, public   :: get_tau          => ohara_precalc_get_tau
      procedure, public   :: get_use_safe_exp => ohara_precalc_get_use_safe_exp
      procedure, public   :: ismodified       => ohara_precalc_modified
      final               :: ohara_precalc_destructor
   END TYPE ohara_precalc

   !public :: exm_ohara_infs
   !public :: exm_ohara_taus

contains

#include "exm_safe_exp.f90.inc"


function ohara_precalc_get_use_safe_exp(this) Result(r)
    implicit none
    class (ohara_precalc) :: this
    integer(ip)           :: r
 
    r = this % use_safe_exp
 end function ohara_precalc_get_use_safe_exp
 




function ohara_precalc_modified(this) Result(r)
   implicit none
   class (ohara_precalc) :: this
   logical(lg)             :: r

   r = this % modified
end function ohara_precalc_modified




subroutine ohara_precalc_get_infs(this, res)
   implicit none
   class (ohara_precalc) :: this
   real(rp), intent(out) :: res(EXM_OHR_INF_EQUATIONS)
   integer(ip)           :: i

   if (this % modified) call this % calculate()
   do i=1,size(this % inf_storage, 1, ip)
        res(i) = this % inf_storage(i)
   end do
end subroutine ohara_precalc_get_infs


subroutine ohara_precalc_get_taus(this, res)
   implicit none
   class (ohara_precalc) :: this
   real(rp), intent(out) :: res(EXM_OHR_TAU_EQUATIONS)
   integer(ip)           :: i

   if (this % modified) call this % calculate()
   do i=1,size(this % tau_storage, 1, ip)
       res(i) = this % tau_storage(i)
   end do
end subroutine ohara_precalc_get_taus


function ohara_precalc_get_inf(this, i) Result(r)
   implicit none
   class (ohara_precalc) :: this
   real(rp)               :: r
   integer(ip)            :: i

   if (this % modified) call this % calculate()
   if ( i>EXM_OHR_INF_EQUATIONS .or. i<1_ip ) then
      call runend("ohara_precalc_get_inf index "//trim(intost(i))//" out of range")
   end if
   r = this % inf_storage(i)
end function ohara_precalc_get_inf


function ohara_precalc_get_tau(this, i) Result(r)
   implicit none
   class (ohara_precalc) :: this
   real(rp)               :: r
   integer(ip)            :: i

   if (this % modified) call this % calculate()
   if ( i>EXM_OHR_TAU_EQUATIONS .or. i<1_ip ) then
      call runend("ohara_precalc_get_tau index "//trim(intost(i))//" out of range")
   end if
   r = this % tau_storage(i)
end function ohara_precalc_get_tau


subroutine ohara_precalc_set_voltage(this, v)
   implicit none
   class (ohara_precalc) :: this
   real(rp)               :: v

   if ( abs(this%v - v) >= this%tol ) then
      this % v = v
      this % modified = .TRUE.
   end if
end subroutine ohara_precalc_set_voltage

subroutine ohara_precalc_init(this, model_type, use_safe_exp) 
   implicit none
   class (ohara_precalc)   :: this
   integer(ip), intent(in) :: model_type
   logical(lg), intent(in) :: use_safe_exp
   integer(4)              :: stat
   character(100)          :: errmsg
  
   !initialize variables directly
   this%modified     = .TRUE.
   this%v            = 0.0_rp
   this%tol          = epsilon(1.0_rp)
   this%model_type   = model_type
   if ( use_safe_exp ) then
        this%use_safe_exp = 1_ip
   else 
        this%use_safe_exp = 0_ip
   end if

   if (.not. any(model_type == (/EXMSLD_CELL_OHARA, EXMSLD_CELL_OHARA_INAPA/)) ) call runend("ohara_precalc_init: unknown model type passed")

   if( associated(this%inf_storage) ) then
        deallocate( this%inf_storage ) 
        nullify   ( this%inf_storage )
   end if
   if( associated(this%tau_storage) ) then 
        deallocate( this%tau_storage ) 
        nullify   ( this%tau_storage )
   end if

   allocate( this%inf_storage( EXM_OHR_INF_EQUATIONS ), STAT=stat, ERRMSG=errmsg ) 
   if( stat .ne. 0 ) call runend( "ohara_precalc_init: allocation of inf_storage failed, #elements="//trim(intost(EXM_OHR_INF_EQUATIONS))//", error="//errmsg )
   allocate( this%tau_storage( EXM_OHR_TAU_EQUATIONS ), STAT=stat, ERRMSG=errmsg ) 
   if( stat .ne. 0 ) call runend( "ohara_precalc_init: allocation of tau_storage failed, #elements="//trim(intost(EXM_OHR_TAU_EQUATIONS))//", error="//errmsg )

   this%inf_storage(:) = 0.0_rp
   this%tau_storage(:) = 0.0_rp

end subroutine ohara_precalc_init

subroutine ohara_precalc_destructor(this)
   implicit none
   type(ohara_precalc) :: this

   if( associated(this%inf_storage) ) then
        deallocate( this%inf_storage ) 
        nullify   ( this%inf_storage )
   end if
   if( associated(this%tau_storage) ) then 
        deallocate( this%tau_storage ) 
        nullify   ( this%tau_storage )
   end if
end subroutine ohara_precalc_destructor

subroutine ohara_precalc_calculate(this)
   implicit none
   class (ohara_precalc) :: this
 
   if( this % modified ) then
      call this % calculate_infs()
      call this % calculate_taus()

      this % modified = .FALSE.
   end if

end subroutine ohara_precalc_calculate


!subroutine exm_ohara_infs(v, res,ipoin,kfl_paral)
!cellmodel = EXMSLD_CELL_OHARA | EXMSLD_CELL_OHARA_INAPA 
subroutine ohara_precalc_calculate_infs(this)
   implicit none
   class (ohara_precalc) :: this

   ! Computes for each i in [1..neqs]:
   ! res(i) = a(i) + b(i) / ( 1 + exp(c(i) ( v + d(i) ) ) )
   !----------------------------------------
   !   Original O'Hara rudy model
   !--------------------------------------
   real(rp), target  :: eq_inf_coefs_ohara(EXM_OHR_INF_EQUATIONS,4) = reshape( (/ &
   0.0_rp,  1.0_rp,             -1.0_rp / 9.871_rp,    39.57_rp,   & ! M_INF
   0.0_rp,  1.0_rp,              1.0_rp / 6.086_rp,     82.9_rp,   & ! H_INF !Passini et al. 6.22_rp,78.5_rp, (original: 6.086_rp,    82.9_rp,)
   0.0_rp,  1.0_rp,              1.0_rp / 6.086_rp,     89.1_rp,   & ! H_CAMK_INF !Passini et al. 6.22_rp,84.7_rp, (original: 6.086_rp,    89.1_rp,) 
   0.0_rp,  1.0_rp,             -1.0_rp / 5.264_rp,    42.85_rp,   & ! M_L_INF
   0.0_rp,  1.0_rp,              1.0_rp / 7.488_rp,    87.61_rp,   & ! H_L_INF
   0.0_rp,  1.0_rp,              1.0_rp / 7.488_rp,    93.81_rp,   & ! H_L_CAMK_INF
   0.0_rp,  1.0_rp,             -1.0_rp / 14.82_rp,   -14.34_rp,   & ! A_INF
   0.0_rp,  1.0_rp / 1.2089_rp, -1.0_rp / 29.3814_rp, -18.4099_rp, & ! INV_TAU_A_1
   0.0_rp,  3.5_rp,              1.0_rp / 29.3814_rp,  100.0_rp,   & ! INV_TAU_A_2
   0.0_rp,  1.0_rp,              1.0_rp / 5.711_rp,    43.94_rp,   & ! I_INF
   0.0_rp,  1.0_rp,              1.0_rp / 151.2_rp,   -213.6_rp,   & ! A_I_FAST 
   0.0_rp,  1.0_rp,             -1.0_rp / 14.82_rp,   -24.34_rp,   & ! A_CAMK_INF
   1.0_rp, -0.5_rp,              1.0_rp / 20.0_rp,     70.0_rp,    & ! DELTA_CAMK_RECOVER
   0.0_rp,  1.0_rp,             -1.0_rp / 4.230_rp,    3.940_rp,   & ! D_INF
   0.0_rp,  1.0_rp,              1.0_rp / 3.696_rp,    19.58_rp,   & ! F_INF
   0.3_rp,  0.6_rp,              1.0_rp / 10.0_rp,    -10.0_rp,    & ! A_F_CA_FAST
   0.0_rp,  1.0_rp,             -1.0_rp / 6.789_rp,    8.337_rp,   & ! X_R_INF
   0.0_rp,  1.0_rp,              1.0_rp / 38.21_rp,    54.81_rp,   & ! A_XR_FAST
   0.0_rp,  1.0_rp,              1.0_rp / 75.0_rp,     55.0_rp,    & ! R_KR_1
   0.0_rp,  1.0_rp,              1.0_rp / 30.0_rp,    -10.0_rp,    & ! R_KR_2
   0.0_rp,  1.0_rp,             -1.0_rp / 8.932_rp,    11.60_rp,   & ! X_S1_INF
   0.0_rp,  1.0_rp,             -1.0_rp / 18.34_rp,   -14.48_rp,   & ! X_KB
   1.0_rp, -0.95_rp,             1.0_rp / 5.0_rp,      70.0_rp     & ! DELTA_EPI
   /), (/ EXM_OHR_INF_EQUATIONS,4_ip /), ORDER= (/ 2_ip,1_ip /) )

   !----------------------------------------
   !   O'Hara Passini model
   !--------------------------------------
   real(rp), parameter :: & 
       inf1=6.22_rp, &   !Hinf gate coef 3
       inf2=78.5_rp, &   !Hinf gate coef 4
       inf3=6.22_rp, &   !Hsp gate coef 3
       inf4=84.7_rp   !Hsp gate coef 4

   real(rp), target  :: eq_inf_coefs_passini(EXM_OHR_INF_EQUATIONS,4) = reshape( (/ &
   0.0_rp,  1.0_rp,             -1.0_rp / 9.871_rp,    39.57_rp,   & ! M_INF
   0.0_rp,  1.0_rp,              1.0_rp / inf1,        inf2,       & ! H_INF !Passini et al. 6.22_rp,78.5_rp, (original: 6.086_rp,    82.9_rp,)
   0.0_rp,  1.0_rp,              1.0_rp / inf3,        inf4,       & ! H_CAMK_INF !Passini et al. 6.22_rp,84.7_rp, (original: 6.086_rp,    89.1_rp,) 
   0.0_rp,  1.0_rp,             -1.0_rp / 5.264_rp,    42.85_rp,   & ! M_L_INF
   0.0_rp,  1.0_rp,              1.0_rp / 7.488_rp,    87.61_rp,   & ! H_L_INF
   0.0_rp,  1.0_rp,              1.0_rp / 7.488_rp,    93.81_rp,   & ! H_L_CAMK_INF
   0.0_rp,  1.0_rp,             -1.0_rp / 14.82_rp,   -14.34_rp,   & ! A_INF
   0.0_rp,  1.0_rp / 1.2089_rp, -1.0_rp / 29.3814_rp, -18.4099_rp, & ! INV_TAU_A_1
   0.0_rp,  3.5_rp,              1.0_rp / 29.3814_rp,  100.0_rp,   & ! INV_TAU_A_2
   0.0_rp,  1.0_rp,              1.0_rp / 5.711_rp,    43.94_rp,   & ! I_INF
   0.0_rp,  1.0_rp,              1.0_rp / 151.2_rp,   -213.6_rp,   & ! A_I_FAST
   0.0_rp,  1.0_rp,             -1.0_rp / 14.82_rp,   -24.34_rp,   & ! A_CAMK_INF
   1.0_rp, -0.5_rp,              1.0_rp / 20.0_rp,     70.0_rp,    & ! DELTA_CAMK_RECOVER
   0.0_rp,  1.0_rp,             -1.0_rp / 4.230_rp,    3.940_rp,   & ! D_INF
   0.0_rp,  1.0_rp,              1.0_rp / 3.696_rp,    19.58_rp,   & ! F_INF
   0.3_rp,  0.6_rp,              1.0_rp / 10.0_rp,    -10.0_rp,    & ! A_F_CA_FAST
   0.0_rp,  1.0_rp,             -1.0_rp / 6.789_rp,    8.337_rp,   & ! X_R_INF
   0.0_rp,  1.0_rp,              1.0_rp / 38.21_rp,    54.81_rp,   & ! A_XR_FAST
   0.0_rp,  1.0_rp,              1.0_rp / 75.0_rp,     55.0_rp,    & ! R_KR_1
   0.0_rp,  1.0_rp,              1.0_rp / 30.0_rp,    -10.0_rp,    & ! R_KR_2
   0.0_rp,  1.0_rp,             -1.0_rp / 8.932_rp,    11.60_rp,   & ! X_S1_INF
   0.0_rp,  1.0_rp,             -1.0_rp / 18.34_rp,   -14.48_rp,   & ! X_KB
   1.0_rp, -0.95_rp,             1.0_rp / 5.0_rp,      70.0_rp     & ! DELTA_EPI
   /), (/ EXM_OHR_INF_EQUATIONS,4_ip /), ORDER= (/ 2_ip,1_ip /) )   

   real(rp), pointer, dimension (:,:) :: eq_inf_coefs => null()
   real(rp)                           :: exps(EXM_OHR_INF_EQUATIONS)
   integer(ip) :: ieq

    select case ( this % model_type )
       case (EXMSLD_CELL_OHARA) !default
          eq_inf_coefs => eq_inf_coefs_ohara
       case (EXMSLD_CELL_OHARA_INAPA) 
          eq_inf_coefs => eq_inf_coefs_passini
       case default
          call runend("exm_ohara_infs(): undefined cell type - "//trim(intost(this % model_type)))
    end select

    do ieq = 1,EXM_OHR_INF_EQUATIONS
        exps(ieq) = safe_exp_cond( eq_inf_coefs(ieq,3) * ( this % v + eq_inf_coefs(ieq,4) ), this%use_safe_exp )
    end do

    do ieq = 1,EXM_OHR_INF_EQUATIONS
        this % inf_storage(ieq) = eq_inf_coefs(ieq,1) + eq_inf_coefs(ieq,2) / ( 1.0_rp + exps(ieq) )
    end do

    contains

#include "exm_safe_exp.f90.inc"

end subroutine ohara_precalc_calculate_infs


 !cellmodel = EXMSLD_CELL_OHARA | EXMSLD_CELL_OHARA_INAPA 
subroutine ohara_precalc_calculate_taus(this)
   implicit none
   class (ohara_precalc) :: this

   ! Computes for each i in [1..neqs]:
   ! res(i) = a(i) + b(i) / ( c(i) exp( d(i) (v + e(i)) ) + f(i) exp( g(i) (v + h(i)) ))

   !----------------------------------------
   !   Original O'Hara rudy model
   !--------------------------------------
   real(rp), target  :: eq_tau_coefs_ohara(EXM_OHR_TAU_EQUATIONS,8) = reshape( (/ &
   0.0_rp,    1.0_rp,    6.765_rp,      1.0_rp / 34.77_rp,  11.64_rp, 8.552_rp,         -1.0_rp / 5.955_rp,   77.42_rp,  & ! TAU_M
   0.0_rp,    1.0_rp,    0.00001432_rp,-1.0_rp / 6.285_rp,  1.196_rp, 6.149_rp,          1.0_rp / 20.27_rp,   0.5096_rp,      & ! TAU_H_FAST !Ohara Original
   0.0_rp,    1.0_rp,    0.009794_rp,  -1.0_rp / 28.05_rp,  17.95_rp, 0.3343_rp,         1.0_rp / 56.66_rp,   5.730_rp,  & ! TAU_H_SLOW
   2.038_rp,  1.0_rp,    0.02136_rp,   -1.0_rp / 8.281_rp,  100.6_rp, 0.3052_rp,         1.0_rp / 38.45_rp,   0.9941_rp,     & ! TAU_J  !Ohara original
   4.562_rp,  1.0_rp,    0.3933_rp,    -1.0_rp / 100.0_rp,  100.0_rp, 0.08004_rp,        1.0_rp / 16.59_rp,   50.0_rp,   & ! TAU_I_FAST
   23.62_rp,  1.0_rp,    0.001416_rp,  -1.0_rp / 59.05_rp,  96.52_rp, 0.00000001780_rp,  1.0_rp / 8.079_rp,   114.1_rp,  & ! TAU_I_SLOW
   1.354_rp,  0.0001_rp, 1.0_rp,        1.0_rp / 15.89_rp, -167.4_rp, 1.0_rp,           -1.0_rp / 0.2154_rp, -12.23_rp,  & ! DELTA_CAMK_DEVELOP
   0.6_rp,    1.0_rp,    1.0_rp,       -0.05_rp,            6.0_rp,   1.0_rp,            0.09_rp,             14.0_rp,   & ! TAU_D
   7.0_rp,    1.0_rp,    0.0045_rp,    -1.0_rp / 10.0_rp,   20.0_rp,  0.0045_rp,         1.0_rp / 10.0_rp,    20.0_rp,   & ! TAU_F_FAST
   1000.0_rp, 1.0_rp,    0.000035_rp,  -1.0_rp / 4.0_rp,    5.0_rp,   0.000035_rp,       1.0_rp / 6.0_rp,     5.0_rp,    & ! TAU_F_SLOW
   7.0_rp,    1.0_rp,    0.04_rp,      -1.0_rp / 7.0_rp,   -4.0_rp,   0.04_rp,           1.0_rp / 7.0_rp,    -4.0_rp,    & ! TAU_F_CA_FAST
   100.0_rp,  1.0_rp,    0.00012_rp,   -1.0_rp / 3.0_rp,    0.0_rp,   0.00012_rp,        1.0_rp / 7.0_rp,     0.0_rp,    & ! TAU_F_CA_SLOW
   12.98_rp,  1.0_rp,    0.3652_rp,     1.0_rp / 3.869_rp, -31.66_rp, 0.00004123_rp,    -1.0_rp / 20.38_rp,  -47.78_rp,  & ! TAU_XR_FAST
   1.865_rp,  1.0_rp,    0.06629_rp,    1.0_rp / 7.355_rp, -34.70_rp, 0.00001128_rp,    -1.0_rp / 25.94_rp,  -29.74_rp,  & ! TAU_XR_SLOW
   817.3_rp,  1.0_rp,    0.0002326_rp,  1.0_rp / 17.80_rp,  48.28_rp, 0.001292_rp,      -1.0_rp / 230.0_rp,   210.0_rp,  & ! TAU_X_S1
   0.0_rp,    1.0_rp,    0.01_rp,       1.0_rp / 20.0_rp,  -50.0_rp,  0.0193_rp,        -1.0_rp / 31.0_rp,    66.54_rp,  & ! TAU_X_S2
   0.0_rp,    122.2_rp,  1.0_rp,       -1.0_rp / 20.36_rp,  127.2_rp, 1.0_rp,            1.0_rp / 69.33_rp,   236.8_rp   & ! TAU_X_K1
   /), (/ EXM_OHR_TAU_EQUATIONS, 8_ip /), ORDER = (/ 2_ip, 1_ip /) )


   !----------------------------------------
   !   O'Hara Passini model
   !--------------------------------------

   real(rp), parameter :: & 
    tau1=0.000003686_rp, & !Hf tau coef 3
    tau2=7.8579_rp, &    !Hf tau coef 4
    tau3=3.8875_rp, &    !Hf tau coef 5
    tau4=16.0_rp, &      !Hf tau coef 6
    tau5=9.1843_rp, &    !Hf tau coef 7
    tau6=-0.4963_rp, &    !Hf tau coef 8           
    tau7=4.859_rp, &     !J tau coef 1
    tau8=0.8628_rp, &    !J tau coef 3
    tau9=7.6005_rp, &    !J tau coef 4
    tau10=116.7258_rp, & !J tau coef 5
    tau11=1.1096_rp, &   !J tau coef 6
    tau12=9.0358_rp, &   !J tau coef 7
    tau13=6.2719_rp   !J tau coef 8   
   
   real(rp), target  :: eq_tau_coefs_passini(EXM_OHR_TAU_EQUATIONS,8) = reshape( (/ &
   0.0_rp,    1.0_rp,    6.765_rp,      1.0_rp / 34.77_rp,  11.64_rp, 8.552_rp,         -1.0_rp / 5.955_rp,   77.42_rp,  & ! TAU_M
   0.0_rp,    1.0_rp,    tau1,         -1.0_rp / tau2,      tau3,     tau4,              1.0_rp / tau5,       tau6,      & ! TAU_H_FAST !Ohara Original
   0.0_rp,    1.0_rp,    0.009794_rp,  -1.0_rp / 28.05_rp,  17.95_rp, 0.3343_rp,         1.0_rp / 56.66_rp,   5.730_rp,  & ! TAU_H_SLOW
   tau7,      1.0_rp,    tau8,         -1.0_rp / tau9,      tau10,    tau11,             1.0_rp / tau12,      tau13,     & ! TAU_J  !Ohara original
   4.562_rp,  1.0_rp,    0.3933_rp,    -1.0_rp / 100.0_rp,  100.0_rp, 0.08004_rp,        1.0_rp / 16.59_rp,   50.0_rp,   & ! TAU_I_FAST
   23.62_rp,  1.0_rp,    0.001416_rp,  -1.0_rp / 59.05_rp,  96.52_rp, 0.00000001780_rp,  1.0_rp / 8.079_rp,   114.1_rp,  & ! TAU_I_SLOW
   1.354_rp,  0.0001_rp, 1.0_rp,        1.0_rp / 15.89_rp, -167.4_rp, 1.0_rp,           -1.0_rp / 0.2154_rp, -12.23_rp,  & ! DELTA_CAMK_DEVELOP
   0.6_rp,    1.0_rp,    1.0_rp,       -0.05_rp,            6.0_rp,   1.0_rp,            0.09_rp,             14.0_rp,   & ! TAU_D
   7.0_rp,    1.0_rp,    0.0045_rp,    -1.0_rp / 10.0_rp,   20.0_rp,  0.0045_rp,         1.0_rp / 10.0_rp,    20.0_rp,   & ! TAU_F_FAST
   1000.0_rp, 1.0_rp,    0.000035_rp,  -1.0_rp / 4.0_rp,    5.0_rp,   0.000035_rp,       1.0_rp / 6.0_rp,     5.0_rp,    & ! TAU_F_SLOW
   7.0_rp,    1.0_rp,    0.04_rp,      -1.0_rp / 7.0_rp,   -4.0_rp,   0.04_rp,           1.0_rp / 7.0_rp,    -4.0_rp,    & ! TAU_F_CA_FAST
   100.0_rp,  1.0_rp,    0.00012_rp,   -1.0_rp / 3.0_rp,    0.0_rp,   0.00012_rp,        1.0_rp / 7.0_rp,     0.0_rp,    & ! TAU_F_CA_SLOW
   12.98_rp,  1.0_rp,    0.3652_rp,     1.0_rp / 3.869_rp, -31.66_rp, 0.00004123_rp,    -1.0_rp / 20.38_rp,  -47.78_rp,  & ! TAU_XR_FAST
   1.865_rp,  1.0_rp,    0.06629_rp,    1.0_rp / 7.355_rp, -34.70_rp, 0.00001128_rp,    -1.0_rp / 25.94_rp,  -29.74_rp,  & ! TAU_XR_SLOW
   817.3_rp,  1.0_rp,    0.0002326_rp,  1.0_rp / 17.80_rp,  48.28_rp, 0.001292_rp,      -1.0_rp / 230.0_rp,   210.0_rp,  & ! TAU_X_S1
   0.0_rp,    1.0_rp,    0.01_rp,       1.0_rp / 20.0_rp,  -50.0_rp,  0.0193_rp,        -1.0_rp / 31.0_rp,    66.54_rp,  & ! TAU_X_S2
   0.0_rp,    122.2_rp,  1.0_rp,       -1.0_rp / 20.36_rp,  127.2_rp, 1.0_rp,            1.0_rp / 69.33_rp,   236.8_rp   & ! TAU_X_K1
   /), (/ EXM_OHR_TAU_EQUATIONS, 8_ip /), ORDER = (/ 2_ip, 1_ip /) )

   real(rp), pointer, dimension (:,:) :: eq_tau_coefs => null()
   integer(ip) :: ieq
   real(rp)                           :: exps1(EXM_OHR_TAU_EQUATIONS), exps2(EXM_OHR_TAU_EQUATIONS)
   
    select case ( this % model_type )
       case (EXMSLD_CELL_OHARA) !default
          eq_tau_coefs => eq_tau_coefs_ohara
       case (EXMSLD_CELL_OHARA_INAPA) 
          eq_tau_coefs => eq_tau_coefs_passini
       case default
          call runend("exm_ohara_infs(): undefined cell type - "//trim(intost(this % model_type)))
    end select

    do ieq = 1,EXM_OHR_TAU_EQUATIONS
        exps1(ieq) = safe_exp_cond( eq_tau_coefs(ieq,4) * ( this%v + eq_tau_coefs(ieq,5)), this%use_safe_exp ) 
        exps2(ieq) = safe_exp_cond( eq_tau_coefs(ieq,7) * ( this%v + eq_tau_coefs(ieq,8)), this%use_safe_exp ) 
    end do


    do ieq = 1,EXM_OHR_TAU_EQUATIONS
        this % tau_storage ( ieq ) = eq_tau_coefs(ieq,1) &
           + eq_tau_coefs(ieq,2) / ( eq_tau_coefs(ieq,3) * exps1(ieq) + eq_tau_coefs(ieq,6) * exps2(ieq) )        
    end do

    contains
#include "exm_safe_exp.f90.inc"
end subroutine ohara_precalc_calculate_taus




 !cellmodel = EXMSLD_CELL_OHARA | EXMSLD_CELL_OHARA_INAPA 
subroutine ohara_precalc_calculate_gates(this, gates_in, gates_out, cell_type, ko, ca_ss, dtimeEP, tauHL, kmn, k2n)
   use def_kintyp, only : ip,rp
   use def_exmedi
   use mod_exm_cellmodel
   implicit none
   class (ohara_precalc) :: this
   real(rp), dimension(:), intent(in)  :: gates_in
   real(rp), dimension(:), intent(out) :: gates_out


   integer(ip), intent(in) :: cell_type !< node
   real(rp), intent(in)  :: ko !< ko
   real(rp), intent(in)  :: ca_ss !< ca_ss,  vcolo(6), ca current in the subspace ss
   real(rp), intent(in)  :: dtimeEP !< dtimeEP
   real(rp), intent(in)  :: tauHL ! tauHL
   real(rp), intent(in)  :: kmn != 0.002_rp
   real(rp), intent(in)  :: k2n != 1000.0_rp

   integer(ip) :: iauxi
   real(rp) :: infs(nauxi_exm),taus(nauxi_exm)
   real(rp) :: anca,dnca,km2n,vaux1,vaux4, vaux2
   real(rp) :: exps(29)

   if( this % modified ) call this % calculate()

   infs(1) = this % inf_storage(OHR_M_INF)
   taus(1) = this % tau_storage(OHR_TAU_M)

   infs(2) = this % inf_storage(OHR_H_INF)
   taus(2) = this % tau_storage(OHR_TAU_H_FAST)

   infs(3) = this % inf_storage(OHR_H_INF)
   taus(3) = this % tau_storage(OHR_TAU_H_SLOW)

   infs(4) = this % inf_storage(OHR_H_INF)  ! j_inf = h_inf
   taus(4) = this % tau_storage(OHR_TAU_J)

   infs(5) = this % inf_storage(OHR_H_CAMK_INF)
   taus(5) = 3.0_rp * this % tau_storage(OHR_TAU_H_SLOW)

   infs(6) = this % inf_storage(OHR_H_INF)  ! j_camk_inf = j_inf = h_inf
   taus(6) = 1.46_rp * this % tau_storage(OHR_TAU_J)

   infs(7) = this % inf_storage(OHR_M_L_INF)
   taus(7) = this % tau_storage(OHR_TAU_M)  ! tau_m_l = tau_m

   infs(8) = this % inf_storage(OHR_H_L_INF)
   taus(8) = tauHL

   infs(9) = this % inf_storage(OHR_H_L_CAMK_INF)
   taus(9) = 3.0_rp * 200.0_rp

   ! NOTE: OHR_INV_TAU_A_1/2 are infs, not taus!
   infs(10) = this % inf_storage(OHR_A_INF)
   taus(10) = 1.0515_rp / (this % inf_storage(OHR_INV_TAU_A_1) + this % inf_storage(OHR_INV_TAU_A_2))

   infs(11) = this % inf_storage(OHR_I_INF)
   taus(11) = this % tau_storage(OHR_TAU_I_FAST)
   if(cell_type == EXM_CELLTYPE_EPI) then
      taus(11) = taus(11) * this % inf_storage(OHR_DELTA_EPI)
   end if

   infs(12) = this % inf_storage(OHR_I_INF)
   taus(12) = this % tau_storage(OHR_TAU_I_SLOW)
   if(cell_type == EXM_CELLTYPE_EPI) then
      taus(12) = taus(12) * this % inf_storage(OHR_DELTA_EPI)
   end if

   infs(13) = this % inf_storage(OHR_A_CAMK_INF)
   taus(13) = taus(10) ! tau_a_camk_inf = tau_a

   infs(14) = this % inf_storage(OHR_I_INF) ! i_camk_inf = i_inf
   taus(14) = this % tau_storage(OHR_DELTA_CAMK_DEVELOP) * this % inf_storage(OHR_DELTA_CAMK_RECOVER) * taus(11)
   !print *,"this % tau_storage(OHR_DELTA_CAMK_DEVELOP)=",this % tau_storage(OHR_DELTA_CAMK_DEVELOP)
   !print *,"this % inf_storage(OHR_DELTA_CAMK_RECOVER)=",this % inf_storage(OHR_DELTA_CAMK_RECOVER) 
   !print *,"taus(11)=", taus(11)
   !print *,"PRE: ", infs(14), gates_in(14), taus(14)


   infs(15) = this % inf_storage(OHR_I_INF) ! i_camk_inf = i_inf
   taus(15) = this % tau_storage(OHR_DELTA_CAMK_DEVELOP) * this % inf_storage(OHR_DELTA_CAMK_RECOVER) * taus(12)

   infs(16) = this % inf_storage(OHR_D_INF)
   taus(16) = this % tau_storage(OHR_TAU_D)

   infs(17) = this % inf_storage(OHR_F_INF)
   taus(17) = this % tau_storage(OHR_TAU_F_FAST)

   infs(18) = this % inf_storage(OHR_F_INF)
   taus(18) = this % tau_storage(OHR_TAU_F_SLOW)

   infs(19) = this % inf_storage(OHR_F_INF) ! f_ca_inf = f_inf
   taus(19) = this % tau_storage(OHR_TAU_F_CA_FAST)

   infs(20) = this % inf_storage(OHR_F_INF) ! f_ca_inf = f_inf
   taus(20) = this % tau_storage(OHR_TAU_F_CA_SLOW)

   infs(21) = this % inf_storage(OHR_F_INF) ! j_ca_inf = f_ca_inf = f_inf
   taus(21) = 75.0_rp

   infs(23) = this % inf_storage(OHR_F_INF) ! f_camk_inf = f_inf
   taus(23) = 2.5_rp * this % tau_storage(OHR_TAU_F_FAST)

   infs(24) = this % inf_storage(OHR_F_INF) ! f_ca_camk_inf = f_inf
   taus(24) = 2.5_rp * this % tau_storage(OHR_TAU_F_CA_FAST)

   infs(25) = this % inf_storage(OHR_X_R_INF)
   taus(25) = this % tau_storage(OHR_TAU_XR_FAST)

   infs(26) = this % inf_storage(OHR_X_R_INF)
   taus(26) = this % tau_storage(OHR_TAU_XR_SLOW)

   infs(27) = this % inf_storage(OHR_X_S1_INF)
   taus(27) = this % tau_storage(OHR_TAU_X_S1)

   infs(28) = this % inf_storage(OHR_X_S1_INF)
   taus(28) = this % tau_storage(OHR_TAU_X_S2)

   infs(29) = 1.0_rp / (1.0_rp + safe_exp_cond(-(this % v + 2.5538_rp * ko + 144.59_rp) / (1.5692_rp * ko + 3.8115_rp), this%use_safe_exp))
   taus(29) = this % tau_storage(OHR_TAU_X_K1)

   do iauxi=1,21
      exps(iauxi) = safe_exp_cond(-dtimeEP/taus(iauxi),this%use_safe_exp)
   end do
   do iauxi=23,29
      exps(iauxi) = safe_exp_cond(-dtimeEP/taus(iauxi),this%use_safe_exp)
   end do


   do iauxi=1,21
      gates_out(iauxi) = infs(iauxi) - (infs(iauxi) - gates_in(iauxi)) * safe_exp_cond(-dtimeEP/taus(iauxi),this%use_safe_exp)
   end do
   do iauxi=23,29
      gates_out(iauxi) = infs(iauxi) - (infs(iauxi) - gates_in(iauxi)) * safe_exp_cond(-dtimeEP/taus(iauxi),this%use_safe_exp)
   end do
        
   !!!!! calculate icana and icak
   km2n = gates_out(21)
   vaux1 = kmn / ca_ss
   vaux4 = (1.0_rp+vaux1) *(1.0_rp+vaux1) * (1.0_rp+vaux1) * (1.0_rp+vaux1)
   vaux2 = (k2n/km2n) + vaux4
   anca  = 1.0_rp/vaux2
   dnca  = (anca*k2n/km2n)
   gates_out(22) = dnca - (dnca - gates_in(22)) * safe_exp_cond(-dtimeEP*km2n, this%use_safe_exp)

   contains
#include "exm_safe_exp.f90.inc"

end subroutine ohara_precalc_calculate_gates



   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
   !!
   !!
   !!  CLASS ENDS HERE
   !!
   !!
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11

!old!
!old!
!old!   !subroutine exm_ohara_infs(v, res,ipoin,kfl_paral)
!old!   !cellmodel = EXMSLD_CELL_OHARA | EXMSLD_CELL_OHARA_INAPA 
!old!   subroutine exm_ohara_infs(v, res, cellmodel)
!old!
!old!     ! Computes for each i in [1..neqs]:
!old!     ! res(i) = a(i) + b(i) / ( 1 + exp(c(i) ( v + d(i) ) ) )
!old!     real(rp), intent(in)    :: v
!old!     real(rp), intent(out)   :: res(EXM_OHR_INF_EQUATIONS)
!old!     integer(ip), intent(in) :: cellmodel
!old!
!old!     !----------------------------------------
!old!     !   Original O'Hara rudy model
!old!     !--------------------------------------
!old!     real(rp), target  :: eq_inf_coefs_ohara(EXM_OHR_INF_EQUATIONS,4) = reshape( (/ &
!old!     0.0_rp,  1.0_rp,             -1.0_rp / 9.871_rp,    39.57_rp,   & ! M_INF
!old!     0.0_rp,  1.0_rp,              1.0_rp / 6.086_rp,     82.9_rp,   & ! H_INF !Passini et al. 6.22_rp,78.5_rp, (original: 6.086_rp,    82.9_rp,)
!old!     0.0_rp,  1.0_rp,              1.0_rp / 6.086_rp,     89.1_rp,   & ! H_CAMK_INF !Passini et al. 6.22_rp,84.7_rp, (original: 6.086_rp,    89.1_rp,) 
!old!     0.0_rp,  1.0_rp,             -1.0_rp / 5.264_rp,    42.85_rp,   & ! M_L_INF
!old!     0.0_rp,  1.0_rp,              1.0_rp / 7.488_rp,    87.61_rp,   & ! H_L_INF
!old!     0.0_rp,  1.0_rp,              1.0_rp / 7.488_rp,    93.81_rp,   & ! H_L_CAMK_INF
!old!     0.0_rp,  1.0_rp,             -1.0_rp / 14.82_rp,   -14.34_rp,   & ! A_INF
!old!     0.0_rp,  1.0_rp / 1.2089_rp, -1.0_rp / 29.3814_rp, -18.4099_rp, & ! INV_TAU_A_1
!old!     0.0_rp,  3.5_rp,              1.0_rp / 29.3814_rp,  100.0_rp,   & ! INV_TAU_A_2
!old!     0.0_rp,  1.0_rp,              1.0_rp / 5.711_rp,    43.94_rp,   & ! I_INF
!old!     0.0_rp,  1.0_rp,              1.0_rp / 151.2_rp,   -213.6_rp,   & ! A_I_FAST
!old!     0.0_rp,  1.0_rp,             -1.0_rp / 14.82_rp,   -24.34_rp,   & ! A_CAMK_INF
!old!     1.0_rp, -0.5_rp,              1.0_rp / 20.0_rp,     70.0_rp,    & ! DELTA_CAMK_RECOVER
!old!     0.0_rp,  1.0_rp,             -1.0_rp / 4.230_rp,    3.940_rp,   & ! D_INF
!old!     0.0_rp,  1.0_rp,              1.0_rp / 3.696_rp,    19.58_rp,   & ! F_INF
!old!     0.3_rp,  0.6_rp,              1.0_rp / 10.0_rp,    -10.0_rp,    & ! A_F_CA_FAST
!old!     0.0_rp,  1.0_rp,             -1.0_rp / 6.789_rp,    8.337_rp,   & ! X_R_INF
!old!     0.0_rp,  1.0_rp,              1.0_rp / 38.21_rp,    54.81_rp,   & ! A_XR_FAST
!old!     0.0_rp,  1.0_rp,              1.0_rp / 75.0_rp,     55.0_rp,    & ! R_KR_1
!old!     0.0_rp,  1.0_rp,              1.0_rp / 30.0_rp,    -10.0_rp,    & ! R_KR_2
!old!     0.0_rp,  1.0_rp,             -1.0_rp / 8.932_rp,    11.60_rp,   & ! X_S1_INF
!old!     0.0_rp,  1.0_rp,             -1.0_rp / 18.34_rp,   -14.48_rp,   & ! X_KB
!old!     1.0_rp, -0.95_rp,             1.0_rp / 5.0_rp,      70.0_rp     & ! DELTA_EPI
!old!     /), (/ EXM_OHR_INF_EQUATIONS,4_ip /), ORDER= (/ 2_ip,1_ip /) )
!old!
!old!     !----------------------------------------
!old!     !   O'Hara Passini model
!old!     !--------------------------------------
!old!     real(rp), parameter :: & 
!old!         inf1=6.22_rp, &   !Hinf gate coef 3
!old!         inf2=78.5_rp, &   !Hinf gate coef 4
!old!         inf3=6.22_rp, &   !Hsp gate coef 3
!old!         inf4=84.7_rp   !Hsp gate coef 4
!old!
!old!     real(rp), target  :: eq_inf_coefs_passini(EXM_OHR_INF_EQUATIONS,4) = reshape( (/ &
!old!     0.0_rp,  1.0_rp,             -1.0_rp / 9.871_rp,    39.57_rp,   & ! M_INF
!old!     0.0_rp,  1.0_rp,              1.0_rp / inf1,        inf2,       & ! H_INF !Passini et al. 6.22_rp,78.5_rp, (original: 6.086_rp,    82.9_rp,)
!old!     0.0_rp,  1.0_rp,              1.0_rp / inf3,        inf4,       & ! H_CAMK_INF !Passini et al. 6.22_rp,84.7_rp, (original: 6.086_rp,    89.1_rp,) 
!old!     0.0_rp,  1.0_rp,             -1.0_rp / 5.264_rp,    42.85_rp,   & ! M_L_INF
!old!     0.0_rp,  1.0_rp,              1.0_rp / 7.488_rp,    87.61_rp,   & ! H_L_INF
!old!     0.0_rp,  1.0_rp,              1.0_rp / 7.488_rp,    93.81_rp,   & ! H_L_CAMK_INF
!old!     0.0_rp,  1.0_rp,             -1.0_rp / 14.82_rp,   -14.34_rp,   & ! A_INF
!old!     0.0_rp,  1.0_rp / 1.2089_rp, -1.0_rp / 29.3814_rp, -18.4099_rp, & ! INV_TAU_A_1
!old!     0.0_rp,  3.5_rp,              1.0_rp / 29.3814_rp,  100.0_rp,   & ! INV_TAU_A_2
!old!     0.0_rp,  1.0_rp,              1.0_rp / 5.711_rp,    43.94_rp,   & ! I_INF
!old!     0.0_rp,  1.0_rp,              1.0_rp / 151.2_rp,   -213.6_rp,   & ! A_I_FAST
!old!     0.0_rp,  1.0_rp,             -1.0_rp / 14.82_rp,   -24.34_rp,   & ! A_CAMK_INF
!old!     1.0_rp, -0.5_rp,              1.0_rp / 20.0_rp,     70.0_rp,    & ! DELTA_CAMK_RECOVER
!old!     0.0_rp,  1.0_rp,             -1.0_rp / 4.230_rp,    3.940_rp,   & ! D_INF
!old!     0.0_rp,  1.0_rp,              1.0_rp / 3.696_rp,    19.58_rp,   & ! F_INF
!old!     0.3_rp,  0.6_rp,              1.0_rp / 10.0_rp,    -10.0_rp,    & ! A_F_CA_FAST
!old!     0.0_rp,  1.0_rp,             -1.0_rp / 6.789_rp,    8.337_rp,   & ! X_R_INF
!old!     0.0_rp,  1.0_rp,              1.0_rp / 38.21_rp,    54.81_rp,   & ! A_XR_FAST
!old!     0.0_rp,  1.0_rp,              1.0_rp / 75.0_rp,     55.0_rp,    & ! R_KR_1
!old!     0.0_rp,  1.0_rp,              1.0_rp / 30.0_rp,    -10.0_rp,    & ! R_KR_2
!old!     0.0_rp,  1.0_rp,             -1.0_rp / 8.932_rp,    11.60_rp,   & ! X_S1_INF
!old!     0.0_rp,  1.0_rp,             -1.0_rp / 18.34_rp,   -14.48_rp,   & ! X_KB
!old!     1.0_rp, -0.95_rp,             1.0_rp / 5.0_rp,      70.0_rp     & ! DELTA_EPI
!old!     /), (/ EXM_OHR_INF_EQUATIONS,4_ip /), ORDER= (/ 2_ip,1_ip /) )   
!old!
!old!     real(rp), pointer, dimension (:,:) :: eq_inf_coefs
!old!     integer(ip) :: ieq
!old!
!old!      select case (cellmodel)
!old!         case (EXMSLD_CELL_OHARA) !default
!old!            eq_inf_coefs => eq_inf_coefs_ohara
!old!         case (EXMSLD_CELL_OHARA_INAPA) 
!old!            eq_inf_coefs => eq_inf_coefs_passini
!old!         case default
!old!            call runend("exm_ohara_infs(): undefined cell type - "//trim(intost(cellmodel)))
!old!      end select
!old!
!old!     do ieq = 1,EXM_OHR_INF_EQUATIONS
!old!        res(ieq) = eq_inf_coefs(ieq,1) &
!old!                 + eq_inf_coefs(ieq,2) / ( 1.0_rp + safe_exp( eq_inf_coefs(ieq,3) * ( v + eq_inf_coefs(ieq,4) ) ) )
!old!     end do
!old!
!old!  end subroutine exm_ohara_infs
!old!
!old!
!old!   !cellmodel = EXMSLD_CELL_OHARA | EXMSLD_CELL_OHARA_INAPA 
!old!  subroutine exm_ohara_taus(v, res, cellmodel)
!old!
!old!     ! Computes for each i in [1..neqs]:
!old!     ! res(i) = a(i) + b(i) / ( c(i) exp( d(i) (v + e(i)) ) + f(i) exp( g(i) (v + h(i)) ))
!old!     real(rp), intent(in)    :: v
!old!     real(rp), intent(out)   :: res(EXM_OHR_TAU_EQUATIONS)
!old!     integer(ip), intent(in) :: cellmodel
!old!
!old!     !----------------------------------------
!old!     !   Original O'Hara rudy model
!old!     !--------------------------------------
!old!     real(rp), target  :: eq_tau_coefs_ohara(EXM_OHR_TAU_EQUATIONS,8) = reshape( (/ &
!old!     0.0_rp,    1.0_rp,    6.765_rp,      1.0_rp / 34.77_rp,  11.64_rp, 8.552_rp,         -1.0_rp / 5.955_rp,   77.42_rp,  & ! TAU_M
!old!     0.0_rp,    1.0_rp,    1.432E-05_rp, -1.0_rp / 6.285_rp,  1.196_rp, 6.149_rp,          1.0_rp / 20.27_rp,   0.5096_rp,      & ! TAU_H_FAST !Ohara Original
!old!     0.0_rp,    1.0_rp,    0.009794_rp,  -1.0_rp / 28.05_rp,  17.95_rp, 0.3343_rp,         1.0_rp / 56.66_rp,   5.730_rp,  & ! TAU_H_SLOW
!old!     2.038_rp,  1.0_rp,    0.02136_rp,   -1.0_rp / 8.281_rp,  100.6_rp, 0.3052_rp,         1.0_rp / 38.45_rp,   0.9941_rp,     & ! TAU_J  !Ohara original
!old!     4.562_rp,  1.0_rp,    0.3933_rp,    -1.0_rp / 100.0_rp,  100.0_rp, 0.08004_rp,        1.0_rp / 16.59_rp,   50.0_rp,   & ! TAU_I_FAST
!old!     23.62_rp,  1.0_rp,    0.001416_rp,  -1.0_rp / 59.05_rp,  96.52_rp, 0.00000001780_rp,  1.0_rp / 8.079_rp,   114.1_rp,  & ! TAU_I_SLOW
!old!     1.354_rp,  0.0001_rp, 1.0_rp,        1.0_rp / 15.89_rp, -167.4_rp, 1.0_rp,           -1.0_rp / 0.2154_rp, -12.23_rp,  & ! DELTA_CAMK_DEVELOP
!old!     0.6_rp,    1.0_rp,    1.0_rp,       -0.05_rp,            6.0_rp,   1.0_rp,            0.09_rp,             14.0_rp,   & ! TAU_D
!old!     7.0_rp,    1.0_rp,    0.0045_rp,    -1.0_rp / 10.0_rp,   20.0_rp,  0.0045_rp,         1.0_rp / 10.0_rp,    20.0_rp,   & ! TAU_F_FAST
!old!     1000.0_rp, 1.0_rp,    0.000035_rp,  -1.0_rp / 4.0_rp,    5.0_rp,   0.000035_rp,       1.0_rp / 6.0_rp,     5.0_rp,    & ! TAU_F_SLOW
!old!     7.0_rp,    1.0_rp,    0.04_rp,      -1.0_rp / 7.0_rp,   -4.0_rp,   0.04_rp,           1.0_rp / 7.0_rp,    -4.0_rp,    & ! TAU_F_CA_FAST
!old!     100.0_rp,  1.0_rp,    0.00012_rp,   -1.0_rp / 3.0_rp,    0.0_rp,   0.00012_rp,        1.0_rp / 7.0_rp,     0.0_rp,    & ! TAU_F_CA_SLOW
!old!     12.98_rp,  1.0_rp,    0.3652_rp,     1.0_rp / 3.869_rp, -31.66_rp, 4.123E-5_rp,      -1.0_rp / 20.38_rp,  -47.78_rp,  & ! TAU_XR_FAST
!old!     1.865_rp,  1.0_rp,    0.06629_rp,    1.0_rp / 7.355_rp, -34.70_rp, 1.128E-5_rp,      -1.0_rp / 25.94_rp,  -29.74_rp,  & ! TAU_XR_SLOW
!old!     817.3_rp,  1.0_rp,    2.326E-4_rp,   1.0_rp / 17.80_rp,  48.28_rp, 0.001292_rp,      -1.0_rp / 230.0_rp,   210.0_rp,  & ! TAU_X_S1
!old!     0.0_rp,    1.0_rp,    0.01_rp,       1.0_rp / 20.0_rp,  -50.0_rp,  0.0193_rp,        -1.0_rp / 31.0_rp,    66.54_rp,  & ! TAU_X_S2
!old!     0.0_rp,    122.2_rp,  1.0_rp,       -1.0_rp / 20.36_rp,  127.2_rp, 1.0_rp,            1.0_rp / 69.33_rp,   236.8_rp   & ! TAU_X_K1
!old!     /), (/ EXM_OHR_TAU_EQUATIONS, 8_ip /), ORDER = (/ 2_ip, 1_ip /) )
!old!
!old!
!old!     !----------------------------------------
!old!     !   O'Hara Passini model
!old!     !--------------------------------------
!old!
!old!     real(rp), parameter :: & 
!old!      tau1=3.686E-06_rp, & !Hf tau coef 3
!old!      tau2=7.8579_rp, &    !Hf tau coef 4
!old!      tau3=3.8875_rp, &    !Hf tau coef 5
!old!      tau4=16.0_rp, &      !Hf tau coef 6
!old!      tau5=9.1843_rp, &    !Hf tau coef 7
!old!      tau6=-0.4963_rp, &    !Hf tau coef 8           
!old!      tau7=4.859_rp, &     !J tau coef 1
!old!      tau8=0.8628_rp, &    !J tau coef 3
!old!      tau9=7.6005_rp, &    !J tau coef 4
!old!      tau10=116.7258_rp, & !J tau coef 5
!old!      tau11=1.1096_rp, &   !J tau coef 6
!old!      tau12=9.0358_rp, &   !J tau coef 7
!old!      tau13=6.2719_rp   !J tau coef 8   
!old!     
!old!     real(rp), target  :: eq_tau_coefs_passini(EXM_OHR_TAU_EQUATIONS,8) = reshape( (/ &
!old!     0.0_rp,    1.0_rp,    6.765_rp,      1.0_rp / 34.77_rp,  11.64_rp, 8.552_rp,         -1.0_rp / 5.955_rp,   77.42_rp,  & ! TAU_M
!old!     0.0_rp,    1.0_rp,    tau1,         -1.0_rp / tau2,      tau3,     tau4,              1.0_rp / tau5,       tau6,      & ! TAU_H_FAST !Ohara Original
!old!     0.0_rp,    1.0_rp,    0.009794_rp,  -1.0_rp / 28.05_rp,  17.95_rp, 0.3343_rp,         1.0_rp / 56.66_rp,   5.730_rp,  & ! TAU_H_SLOW
!old!     tau7,      1.0_rp,    tau8,         -1.0_rp / tau9,      tau10,    tau11,             1.0_rp / tau12,      tau13,     & ! TAU_J  !Ohara original
!old!     4.562_rp,  1.0_rp,    0.3933_rp,    -1.0_rp / 100.0_rp,  100.0_rp, 0.08004_rp,        1.0_rp / 16.59_rp,   50.0_rp,   & ! TAU_I_FAST
!old!     23.62_rp,  1.0_rp,    0.001416_rp,  -1.0_rp / 59.05_rp,  96.52_rp, 0.00000001780_rp,  1.0_rp / 8.079_rp,   114.1_rp,  & ! TAU_I_SLOW
!old!     1.354_rp,  0.0001_rp, 1.0_rp,        1.0_rp / 15.89_rp, -167.4_rp, 1.0_rp,           -1.0_rp / 0.2154_rp, -12.23_rp,  & ! DELTA_CAMK_DEVELOP
!old!     0.6_rp,    1.0_rp,    1.0_rp,       -0.05_rp,            6.0_rp,   1.0_rp,            0.09_rp,             14.0_rp,   & ! TAU_D
!old!     7.0_rp,    1.0_rp,    0.0045_rp,    -1.0_rp / 10.0_rp,   20.0_rp,  0.0045_rp,         1.0_rp / 10.0_rp,    20.0_rp,   & ! TAU_F_FAST
!old!     1000.0_rp, 1.0_rp,    0.000035_rp,  -1.0_rp / 4.0_rp,    5.0_rp,   0.000035_rp,       1.0_rp / 6.0_rp,     5.0_rp,    & ! TAU_F_SLOW
!old!     7.0_rp,    1.0_rp,    0.04_rp,      -1.0_rp / 7.0_rp,   -4.0_rp,   0.04_rp,           1.0_rp / 7.0_rp,    -4.0_rp,    & ! TAU_F_CA_FAST
!old!     100.0_rp,  1.0_rp,    0.00012_rp,   -1.0_rp / 3.0_rp,    0.0_rp,   0.00012_rp,        1.0_rp / 7.0_rp,     0.0_rp,    & ! TAU_F_CA_SLOW
!old!     12.98_rp,  1.0_rp,    0.3652_rp,     1.0_rp / 3.869_rp, -31.66_rp, 4.123E-5_rp,      -1.0_rp / 20.38_rp,  -47.78_rp,  & ! TAU_XR_FAST
!old!     1.865_rp,  1.0_rp,    0.06629_rp,    1.0_rp / 7.355_rp, -34.70_rp, 1.128E-5_rp,      -1.0_rp / 25.94_rp,  -29.74_rp,  & ! TAU_XR_SLOW
!old!     817.3_rp,  1.0_rp,    2.326E-4_rp,   1.0_rp / 17.80_rp,  48.28_rp, 0.001292_rp,      -1.0_rp / 230.0_rp,   210.0_rp,  & ! TAU_X_S1
!old!     0.0_rp,    1.0_rp,    0.01_rp,       1.0_rp / 20.0_rp,  -50.0_rp,  0.0193_rp,        -1.0_rp / 31.0_rp,    66.54_rp,  & ! TAU_X_S2
!old!     0.0_rp,    122.2_rp,  1.0_rp,       -1.0_rp / 20.36_rp,  127.2_rp, 1.0_rp,            1.0_rp / 69.33_rp,   236.8_rp   & ! TAU_X_K1
!old!     /), (/ EXM_OHR_TAU_EQUATIONS, 8_ip /), ORDER = (/ 2_ip, 1_ip /) )
!old!
!old!     real(rp), pointer, dimension (:,:) :: eq_tau_coefs
!old!     integer(ip) :: ieq
!old!
!old!      select case (cellmodel)
!old!         case (EXMSLD_CELL_OHARA) !default
!old!            eq_tau_coefs => eq_tau_coefs_ohara
!old!         case (EXMSLD_CELL_OHARA_INAPA) 
!old!            eq_tau_coefs => eq_tau_coefs_passini
!old!         case default
!old!            call runend("exm_ohara_infs(): undefined cell type - "//trim(intost(cellmodel)))
!old!      end select
!old!
!old!
!old!     do ieq = 1,EXM_OHR_TAU_EQUATIONS
!old!        res(ieq) = eq_tau_coefs(ieq,1) &
!old!             + eq_tau_coefs(ieq,2) / ( eq_tau_coefs(ieq,3) * safe_exp( eq_tau_coefs(ieq,4) * (v + eq_tau_coefs(ieq,5)) ) +&
!old!             eq_tau_coefs(ieq,6) * safe_exp( eq_tau_coefs(ieq,7) * (v + eq_tau_coefs(ieq,8)) ) )        
!old!     end do
!old!  end subroutine exm_ohara_taus
!old!  


  

end module mod_exm_oharaprecalc

