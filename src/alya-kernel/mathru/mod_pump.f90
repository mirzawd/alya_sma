!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!
!> @defgroup Module for generic windkessel functions
!> Subroutine to implement and use Windkessel functions
!> the module allows to use different modules as input
!> and output
!> @{
!> @name    Module for generic windkessel functions
!> @file    mod_pump.f90
!> @author  Beatriz Eguzkitza
!> @brief    Module for generic pump curve functions
!> @details  Module for generic pump_curve functions
!
!-----------------------------------------------------------------------

module mod_pump
  use def_kintyp,         only : lg,ip,rp
  use def_master,         only : ID_NASTIN
  use def_master,         only : INOTMASTER
  use def_master,         only : mem_modul,modul
  use def_kermod,         only : number_pump_curve
  use def_kermod,         only : pump_curve
  use mod_memory,         only : memory_alloca
  use mod_memory,         only : memory_deallo
  use mod_parall,         only : PAR_MY_CODE_RANK
  use mod_communications, only : PAR_SUM
  use mod_communications, only : PAR_BARRIER

  implicit none



  public  :: pump_curve_flow
  public  :: pump_curve_par_exchange

  contains



   !-----------------------------------------------------------------------
   !>
   !> @author  
   !> @date    2019-12-17
   !> @brief   Distribute values of the windkessel system among the slaves
   !> @details Distribute values of the windkessel system among the slaves
   !>
   !>
   !-----------------------------------------------------------------------
   subroutine pump_curve_par_exchange()
     implicit none
!     integer(ip)                       :: ifunc,i


  !   do ifunc = 1,number_windk_systems
  !     call PAR_SUM( windk_systems(ifunc) % sysid     ,"IN MY CODE", "INCLUDE MASTER")
       !wdks_model already exchanged
       !nparam     already exchanged
       !params     already exchanged
  !     call PAR_SUM( windk_systems(ifunc) % ID_IN     ,"IN MY CODE", "INCLUDE MASTER") 
  !     call PAR_SUM( windk_systems(ifunc) % ID_OUT    ,"IN MY CODE", "INCLUDE MASTER")   
  !     call PAR_SUM( windk_systems(ifunc) % tag_in    ,"IN MY CODE", "INCLUDE MASTER") 
  !     call PAR_SUM( windk_systems(ifunc) % tag_out   ,"IN MY CODE", "INCLUDE MASTER") 
       !ndxs        already exchanged
  !     call PAR_SUM( windk_systems(ifunc) % iflow_nsi ,"IN MY CODE", "INCLUDE MASTER") 
  !   end do


   end subroutine pump_curve_par_exchange


   !-----------------------------------------------------------------------
   !>
   !> @author  Beatriz Eguzkitza
   !> @date    2020-07-02
   !> @brief   
   !> @details 
   !>
   !>
   !>
   !>
   !-----------------------------------------------------------------------
   subroutine pump_curve_flow(ifunc,delta_P,Q)
        use def_master,        only: cutim
        use def_kermod,        only: time_function, number_time_function
        implicit none
        integer(ip), intent(in) :: ifunc
        real(rp), intent(in)    :: delta_P
        real(rp), intent(out)   :: Q
        real(rp), external      :: funcre
        real(rp)                :: a, b, c !Parameters of the parabole
        real(rp)                :: s, H0rpm, kpp, kpn !Parameters of the HVAD
        real(rp)                :: Hpp, dpaux
        real(rp)                :: Q_calc, Q_max, Q_min, delta_P_min, delta_P_max !Maximum and minimum values of the domain and image
        integer(ip)             :: auxi, fun_code

        if (pump_curve(ifunc) %model .eq. 1_ip ) then
          ! PARABOLIC_CURVE: Q= param(1) * Delta_p^2 +param(2)*Delta_p + param(3) 
          ! inside a DeltP range: param(4), param(5)
          !
          ! Q = a * AP^2 + b * AP + c
          !
          ! delta_P =  timev
          ! param(1) = a
          ! param(2) = b
          ! param(3) = c
          ! range for deltaP = param(4), param(5)
          a = pump_curve(ifunc) %params(1)  
          b = pump_curve(ifunc) %params(2)
          c = pump_curve(ifunc) %params(3)
          delta_P_min = pump_curve(ifunc) % params(4) 
          delta_P_max = pump_curve(ifunc) % params(5) 
    
          Q_max = a*delta_p_max*delta_P_max + b*delta_P_max + c
          Q_min = a*delta_p_min*delta_P_min + b*delta_P_min + c

          Q_calc= a*delta_P*delta_P + b*delta_P + c
        elseif (pump_curve(ifunc) %model .eq. 2_ip ) then
          !HVAD pump curve.
          ! ORIG:
          !       Hp= Hpp(s)-K_{pp}*Q^2   / Q>=0
          !       Hp= Hpp(s)+K_{pn}*Q^2   / Q<0
          !
          !       Hpp=(S[RPM]/H0_RPM)^2
          !
          ! MODIFIED:
          !      Q=sqrt((Hpp-Hp)/Kpp)  / suck
          !      Q=sqrt((Hp-Hpp)/Kpn)  / blow
          s = pump_curve(ifunc) %params(1)  
          H0rpm = pump_curve(ifunc) %params(2)
          Kpp   = pump_curve(ifunc) %params(3)
          Kpn   = pump_curve(ifunc) %params(4)
          delta_P_min = pump_curve(ifunc) % params(5) 
          delta_P_max = pump_curve(ifunc) % params(6) 

          ! Time variable speed
          if( pump_curve(ifunc) % vhvad .ne. '     ') then
            do auxi = 1,number_time_function
             if( trim(pump_curve(ifunc) % vhvad ) == trim(time_function(auxi) % name) ) then
                fun_code = auxi
                exit
             end if
            end do
            s = funcre( time_function(fun_code) % parameters, &
                        time_function(fun_code) % npara,      &
                        time_function(fun_code) % kfl_type,   &
                        cutim)

          endif

          !call runend ('die')


          Hpp=(s/H0rpm)**2.0_rp
          delta_P_min = 0.0_rp
          delta_P_max = Hpp
          Q_min=sqrt(Hpp/Kpp) ! Flow for minimum delta_P
          Q_max=0.0_rp        ! Flow for maximum delta_p


          dpaux=delta_P
          if(delta_P .le. delta_P_min) then
                dpaux=delta_P_min
          elseif((delta_P .gt. delta_P_min) .and. (delta_P .lt. delta_P_max)) then
                dpaux=delta_P
          elseif(delta_P .ge. delta_P_max) then
                dpaux=delta_P_max
          endif

          Q_calc=sqrt((Hpp-dpaux)/Kpp)
          
        else
          WRITE(6,*) ifunc, pump_curve(ifunc) %model
          call runend('mathru/mod_pump: no pump model detected')
        endif


        ! Clipping flow to imposed limits
        if (delta_P > delta_P_min .and. delta_P < delta_P_max )then
         Q = Q_calc
        else if (delta_P < delta_P_min ) then
         Q = Q_min
        else if (delta_P > delta_P_max ) then
         Q = Q_max
        end if
        
      end subroutine pump_curve_flow



end module mod_pump

