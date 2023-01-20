!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Alefor
!> @{
!> @file    mod_ale_rigid_body_vtx.f90
!> @author  J.C. Cajas
!> @date    5/04/2017
!> @brief   Rigid body equations solver with RK4 method
!> @details Rigid body equations solver with RK4 method
!> @} 
!-----------------------------------------------------------------------
module mod_ale_rigid_body_rk4

  use def_kintyp, only : ip, rp
  use def_parame, only : pi
  use def_master, only : INOTMASTER, cutim
  use mod_parall, only : PAR_MY_CODE_RANK
  use mod_memory, only : memory_alloca
  use mod_memory, only : memory_deallo

  implicit none

  private


  public :: ale_rigid_body_rk4_solution
  public :: ale_rigid_body_rk4_ang_trans1
  public :: ale_rigid_body_rk4_rot_angs1
  public :: ale_rigid_body_rk4_rot_angs2
  
contains 

  !-----------------------------------------------------------------------
  !
  ! Implements : Solution of the equations of motion for a Rigid body with 
  !              a fixed point for the vortex-bladeless case, four
  !              equations of motion are solved
  !
  !         dot{theta}   = v_theta
  !
  !         dot{phi}     = v_phi
  !
  !         dot{v}_theta = frac{I-I_Z}{I} v_phi^2 sin (theta) cos (theta) + frac{mg-F_z}{I} h sin (theta) +
  !                          + frac{F_x cos(phi)+F_y sin(phi)}{I} h cos(theta) - frac{k theta}{I} -
  !                          - frac{alpha_theta v_theta}{I}
  !
  !         dot{v}_phi   = - frac{k phi}{I sin^2 (theta) + I_Z cos^2 (theta)}+
  !                          + frac{F_x sin (phi)+F_y cos (phi)}{I sin^2 (theta) + I_Z cos^2 (theta)} h -
  !                          - frac{alpha_phi v_phi}{I sin^2(theta) + I_Z cos^2 (theta)}.
  !
  !
  !    theta and phi are the Euler angles. These equations are obtained using the Lagrangian
  !    formulation and can be obtained from an analysis similar to the symmetrical top
  !-----------------------------------------------------------------------
  subroutine ale_rigid_body_rk4_solution(theta, phi, vel_theta, vel_phi, dtrb, Irb, Frb, hcm, krb, difrb, mrb, grb, option)

    implicit none

    real(rp), intent(inout) :: theta                ! generalized coordinate 
    real(rp), intent(inout) :: phi                  ! generalized coordinate
    real(rp), intent(inout) :: vel_theta            ! generalized velocity
    real(rp), intent(inout) :: vel_phi              ! generalized velocity
    
    real(rp), intent(in)    :: dtrb, krb, difrb     ! time increment, elasticity constant, dissipation constant
    real(rp), intent(in)    :: hcm                  ! distance between the center of mass of the r.b. and the fixed point
    real(rp), intent(in)    :: mrb, grb             ! mass of the R.B.
    real(rp), intent(in)    :: Irb(3_ip), Frb(3_ip) ! inertia tensor(moving frame), forces on the body(inertial frame)

    character(5), intent(in):: option

    
    real(rp)                :: epsilon              ! amplitude perturbation parameter
    real(rp)                :: tau_c                ! control parameter for the lenght of the perturbation
    real(rp)                :: tim_in               ! initial time for the perturbation
    real(rp)                :: omega_p              ! perturbation frequency
    
    real(rp)                :: Frbp(3_ip)           ! perturbed forces on the body(inertial frame)
    
    integer(ip) :: kk, ii

    if( INOTMASTER ) then
       
       !-----------------------------------------------------------------
       !
       ! Allocate memory if necessary
       !
       !-----------------------------------------------------------------
       
       call ale_rigid_body_rk4_alloca()
       
       !
       ! Solve the equations of motion
       !
       
       do kk = 1_ip, 1_ip
          !
          ! generalized coordinates
          !
          call ale_rigid_body_rk4_gen_coor(theta, theta, vel_theta, dtrb)
          call ale_rigid_body_rk4_gen_coor(phi,   phi,   vel_phi,   dtrb)
          !
          ! generalized velocities
          !
          if ( option == 'TOPMO') then
             !
             ! Non-spinning top model
             !
             call ale_rigid_body_rk4_gen_vel (vel_phi,   vel_phi,   phi,   theta, vel_theta, dtrb, &
                  & Irb, Frb, hcm, krb, difrb, mrb, grb, 'phi'  )
             call ale_rigid_body_rk4_gen_vel (vel_theta, vel_theta, theta, phi,   vel_phi,   dtrb, &
                  & Irb, Frb, hcm, krb, difrb, mrb, grb, 'the')
          else if( option == 'SPRIN' )then
             !
             ! Linear spring model
             !
             call ale_rigid_body_rk4_gen_vel (vel_theta, vel_theta,   theta,   theta, vel_theta, dtrb, &
                  & Irb, Frb, hcm, krb, difrb, mrb, grb, 'vex'  )
             call ale_rigid_body_rk4_gen_vel (vel_phi, vel_phi, phi, phi,   vel_phi,   dtrb, &
                  & Irb, Frb, hcm, krb, difrb, mrb, grb, 'vey') 

!!$             call ale_rigid_body_rk4_gen_vel (vel_phi,   vel_phi,   phi,   theta, vel_theta, dtrb, &
!!$                  & Irb, Frb, hcm, krb, difrb, mrb, grb, 'vex'  )
!!$             call ale_rigid_body_rk4_gen_vel (vel_theta, vel_theta, theta, phi,   vel_phi,   dtrb, &
!!$                  & Irb, Frb, hcm, krb, difrb, mrb, grb, 'vey')
          else if( option == 'PERTU' )then
             !
             ! Modify the aerodynamic forces to perturb the base flow
             !
             tim_in = 2e-3_rp
             if( cutim > tim_in )then
                epsilon = 4e-1_rp
             else
                epsilon = 0_rp
             end if
             tau_c   = 0.40_rp
             omega_p = 188.88_rp
             do ii = 1_ip, 3_ip
                
                Frbp(ii) = Frb(ii) * ( 1_rp - epsilon * dexp(-(cutim-tim_in)/tau_c) * dsin(omega_p*cutim) )

             end do
             !
             ! Perturbation of non-spinning top model
             !
             call ale_rigid_body_rk4_gen_vel (vel_phi,   vel_phi,   phi,   theta, vel_theta, dtrb, &
                  & Irb, Frbp, hcm, krb, difrb, mrb, grb, 'phi'  )
             call ale_rigid_body_rk4_gen_vel (vel_theta, vel_theta, theta, phi,   vel_phi,   dtrb, &
                  & Irb, Frbp, hcm, krb, difrb, mrb, grb, 'the')
             
          end if

          
       end do

       
    end if   ! Inotmaster
       
  end subroutine ale_rigid_body_rk4_solution

  
  !-----------------------------------------------------------------------
  !
  !> @brief   Allocate
  !> @details Allocate
  !
  !-----------------------------------------------------------------------

  subroutine ale_rigid_body_rk4_alloca()

    implicit none
    

    if( INOTMASTER ) then

    end if   ! Inotmaster

  end subroutine ale_rigid_body_rk4_alloca
  
  !-----------------------------------------------------------------
  !
  ! Runge-Kutta4 algorithm for constant function (the generalized coordinates theta and phi):
  ! dot{theta}   = v_theta
  ! dot{phi}     = v_phi
  !
  !-----------------------------------------------------------------
  subroutine ale_rigid_body_rk4_gen_coor(gen_coor_aux,gen_coor,gen_vel_aux,dtk4)

    implicit none
    
    real(rp), intent(out) :: gen_coor_aux
    real(rp), intent(in)  :: dtk4, gen_vel_aux, gen_coor
    real(rp)              :: k1, k2, k3, k4
    
    if( INOTMASTER ) then
       
       k1 = gen_vel_aux
       k2 = gen_vel_aux
       k3 = gen_vel_aux
       k4 = gen_vel_aux
       
       gen_coor_aux = gen_coor + dtk4 * ( k1 + 2_rp * k2 + 2_rp * k3 + k4 )/6_rp
       
    end if  ! Inotmaster
    
  end subroutine ale_rigid_body_rk4_gen_coor

  !-----------------------------------------------------------------
  !
  ! Runge-Kutta4 algorithm for the generalized velocities (v_theta, v_phi)
  !
  ! dot{v}_theta = frac{I-I_Z}{I} v_phi^2 sin (theta) cos (theta) + frac{mg-F_z}{I} h sin (theta) +
  !              + frac{F_x cos(phi)+F_y sin(phi)}{I} h cos(theta) - frac{k theta}{I} -
  !              - frac{alpha_theta v_theta}{I}
  !
  ! dot{v}_phi   = - frac{k phi}{I sin^2 (theta) + I_Z cos^2 (theta)}+
  !              + frac{F_x sin (phi)+F_y cos (phi)}{I sin^2 (theta) + I_Z cos^2 (theta)} h -
  !              - frac{alpha_phi v_phi}{I sin^2(theta) + I_Z cos^2 (theta)}
  !
  !-----------------------------------------------------------------
  subroutine ale_rigid_body_rk4_gen_vel(gen_vel1_aux,gen_vel1,gen_coor1_aux,gen_coor2_aux,gen_vel2_aux,dtk4,Irb,Frb,hcm,krb,difrb,mrb,grb,variable)

    implicit none
    
    real(rp), intent(out)     :: gen_vel1_aux
    real(rp), intent(in)      :: gen_vel1,gen_coor1_aux,gen_coor2_aux,gen_vel2_aux,dtk4,hcm,krb,difrb,mrb,grb
    real(rp), intent(in)      :: Frb(3_ip), Irb(3_ip)
    character(3), intent(in)  :: variable
    
    real(rp)                  :: k1, k2, k3, k4

    if( INOTMASTER ) then
       
       select case(trim(variable))
          
       case('the')
          
          k1 = f_vtheta_dot(Irb, Frb, hcm, krb, mrb, grb, difrb, gen_coor1_aux, gen_coor2_aux, gen_vel1                     , gen_vel2_aux)
          k2 = f_vtheta_dot(Irb, Frb, hcm, krb, mrb, grb, difrb, gen_coor1_aux, gen_coor2_aux, gen_vel1 + dtk4 * k1 * 0.5_rp, gen_vel2_aux)
          k3 = f_vtheta_dot(Irb, Frb, hcm, krb, mrb, grb, difrb, gen_coor1_aux, gen_coor2_aux, gen_vel1 + dtk4 * k2 * 0.5_rp, gen_vel2_aux)
          k4 = f_vtheta_dot(Irb, Frb, hcm, krb, mrb, grb, difrb, gen_coor1_aux, gen_coor2_aux, gen_vel1 + dtk4 * k3         , gen_vel2_aux)
          
       case('phi')
          
          k1 =   f_vphi_dot(Irb, Frb, hcm, krb, mrb, grb, difrb, gen_coor1_aux, gen_coor2_aux, gen_vel1                     , gen_vel2_aux)
          k2 =   f_vphi_dot(Irb, Frb, hcm, krb, mrb, grb, difrb, gen_coor1_aux, gen_coor2_aux, gen_vel1 + dtk4 * k1 * 0.5_rp, gen_vel2_aux)
          k3 =   f_vphi_dot(Irb, Frb, hcm, krb, mrb, grb, difrb, gen_coor1_aux, gen_coor2_aux, gen_vel1 + dtk4 * k2 * 0.5_rp, gen_vel2_aux)
          k4 =   f_vphi_dot(Irb, Frb, hcm, krb, mrb, grb, difrb, gen_coor1_aux, gen_coor2_aux, gen_vel1 + dtk4 * k3         , gen_vel2_aux)
       case('alp')
          
          k1 =   f_valpha_dot(Irb, Frb, hcm, krb, mrb, grb, difrb, gen_coor1_aux, gen_coor2_aux, gen_vel1                     , gen_vel2_aux)
          k2 =   f_valpha_dot(Irb, Frb, hcm, krb, mrb, grb, difrb, gen_coor1_aux, gen_coor2_aux, gen_vel1 + dtk4 * k1 * 0.5_rp, gen_vel2_aux)
          k3 =   f_valpha_dot(Irb, Frb, hcm, krb, mrb, grb, difrb, gen_coor1_aux, gen_coor2_aux, gen_vel1 + dtk4 * k2 * 0.5_rp, gen_vel2_aux)
          k4 =   f_valpha_dot(Irb, Frb, hcm, krb, mrb, grb, difrb, gen_coor1_aux, gen_coor2_aux, gen_vel1 + dtk4 * k3         , gen_vel2_aux)
       
        case('bet')
          
          k1 =   f_vbeta_dot(Irb, Frb, hcm, krb, mrb, grb, difrb, gen_coor1_aux, gen_coor2_aux, gen_vel1                     , gen_vel2_aux)
          k2 =   f_vbeta_dot(Irb, Frb, hcm, krb, mrb, grb, difrb, gen_coor1_aux, gen_coor2_aux, gen_vel1 + dtk4 * k1 * 0.5_rp, gen_vel2_aux)
          k3 =   f_vbeta_dot(Irb, Frb, hcm, krb, mrb, grb, difrb, gen_coor1_aux, gen_coor2_aux, gen_vel1 + dtk4 * k2 * 0.5_rp, gen_vel2_aux)
          k4 =   f_vbeta_dot(Irb, Frb, hcm, krb, mrb, grb, difrb, gen_coor1_aux, gen_coor2_aux, gen_vel1 + dtk4 * k3         , gen_vel2_aux)

       case('vex')

          k1 =   f_spring(Frb(1_ip), krb, mrb, difrb, gen_coor1_aux, gen_vel1                     )
          k2 =   f_spring(Frb(1_ip), krb, mrb, difrb, gen_coor1_aux, gen_vel1 + dtk4 * k1 * 0.5_rp)
          k3 =   f_spring(Frb(1_ip), krb, mrb, difrb, gen_coor1_aux, gen_vel1 + dtk4 * k2 * 0.5_rp)
          k4 =   f_spring(Frb(1_ip), krb, mrb, difrb, gen_coor1_aux, gen_vel1 + dtk4 * k3         )

       case('vey')

          k1 =   f_spring(Frb(2_ip), krb, mrb, difrb, gen_coor1_aux, gen_vel1                     )
          k2 =   f_spring(Frb(2_ip), krb, mrb, difrb, gen_coor1_aux, gen_vel1 + dtk4 * k1 * 0.5_rp)
          k3 =   f_spring(Frb(2_ip), krb, mrb, difrb, gen_coor1_aux, gen_vel1 + dtk4 * k2 * 0.5_rp)
          k4 =   f_spring(Frb(2_ip), krb, mrb, difrb, gen_coor1_aux, gen_vel1 + dtk4 * k3         )

          
       end select
       
       gen_vel1_aux = gen_vel1 + dtk4 * (k1 + 2_rp * k2 + 2_rp * k3 + k4)/6_rp
       
    end if  ! Inotmaster
    
  end subroutine ale_rigid_body_rk4_gen_vel
  
  !-----------------------------------------------------------------
  !
  !  Tranformation of angles to find 'equivalent' rotations, takes
  !  R_z(ang2)R_y(ang1) k^ = R_x(out_ang2)R_y(out_ang1) k^ using
  !  the relations
  !
  !  out_ang1 = arcsin( cos(ang2)sin(ang1) )
  !  out_ang2 = arctan(-sin(ang2)tan(ang1) )
  !
  !-----------------------------------------------------------------
  subroutine ale_rigid_body_rk4_ang_trans1(ang1,ang2,out_ang1,out_ang2)
    
    implicit none
    
    real(rp), intent(in)      :: ang1,ang2
    real(rp), intent(out)     :: out_ang1, out_ang2
    
    out_ang1 = dasin( dcos(ang2) * dsin(ang1) )
    out_ang2 = datan(-dsin(ang2) * dtan(ang1) )
    
  end subroutine ale_rigid_body_rk4_ang_trans1
  
  !-----------------------------------------------------------------
  !
  !  Calculate rotation matrix composition
  !  first rotate ang_1 around y axis, then ang2 around x axis
  !
  !  Rot_mat = R_x(ang2) R_y(ang1)
  !
  !-----------------------------------------------------------------
  subroutine ale_rigid_body_rk4_rot_angs1(rot_mat,ang1,ang2)

    implicit none
    
    real(kind=8), intent(in)      :: ang1,ang2
    real(kind=8), intent(out)     :: rot_mat(3,3)
    
    rot_mat(1,1) = dcos(ang1)
    rot_mat(2,2) = dcos(ang2)
    rot_mat(3,3) = dcos(ang1) * dcos(ang2)
    
    rot_mat(1,2) = 0.d0
    rot_mat(2,1) =-dsin(ang1) * dsin(ang2)

    rot_mat(1,3) =-dsin(ang1)
    rot_mat(3,1) = dsin(ang1) * dcos(ang2)
    
    rot_mat(2,3) =-dcos(ang1) * dsin(ang2)
    rot_mat(3,2) = dsin(ang2)
    
  end subroutine ale_rigid_body_rk4_rot_angs1
 
!!$  !-----------------------------------------------------------------
!!$  !
!!$  !  Calculate rotation matrix composition
!!$  !  first rotate ang_1 around y axis, then ang2 around x axis
!!$  !
!!$  !  Rot_mat = R_x(ang2) R_y(ang1)
!!$  !
!!$  !-----------------------------------------------------------------
!!$  subroutine ale_rigid_body_rk4_rot_angs1(rot_mat,ang1,ang2)
!!$
!!$    implicit none
!!$    
!!$    real(rp), intent(in)      :: ang1,ang2
!!$    real(rp), intent(out)     :: rot_mat(3_ip,3_ip)
!!$    
!!$    rot_mat(1_ip,1_ip) = dcos(ang1)
!!$    rot_mat(2_ip,2_ip) = dcos(ang2)
!!$    rot_mat(3_ip,3_ip) = dcos(ang2)*dcos(ang1)
!!$    
!!$    rot_mat(1_ip,2_ip) = 0_rp
!!$    rot_mat(2_ip,1_ip) = dsin(ang2) * dsin(ang1)
!!$    
!!$    rot_mat(1_ip,3_ip) = dsin(ang1)
!!$    rot_mat(3_ip,1_ip) =-dcos(ang2)*dsin(ang1)
!!$    
!!$    rot_mat(2_ip,3_ip) =-dsin(ang2)*dcos(ang1)
!!$    rot_mat(3_ip,2_ip) = dsin(ang2)
!!$    
!!$  end subroutine ale_rigid_body_rk4_rot_angs1

  !-----------------------------------------------------------------
  !
  !  Calculate rotation matrix composition
  !  first rotate ang_1 around x axis, then ang2 around y axis
  !
  !  Rot_mat = R_y(ang2) R_x(ang1)
  !
  !-----------------------------------------------------------------
  subroutine ale_rigid_body_rk4_rot_angs2(rot_mat,ang1,ang2)

    implicit none
    
    real(rp), intent(in)      :: ang1,ang2
    real(rp), intent(out)     :: rot_mat(3_ip,3_ip)
    
    rot_mat(1_ip,1_ip) = dcos(ang2)
    rot_mat(2_ip,2_ip) = dcos(ang1)
    rot_mat(3_ip,3_ip) = dcos(ang2)*dcos(ang1)
    
    rot_mat(1_ip,2_ip) = dsin(ang2) * dsin(ang1)
    rot_mat(2_ip,1_ip) = 0_rp
    
    rot_mat(1_ip,3_ip) = dcos(ang1)*dsin(ang2)
    rot_mat(3_ip,1_ip) =-dsin(ang2)
    
    rot_mat(2_ip,3_ip) =-dsin(ang1)
    rot_mat(3_ip,2_ip) = dsin(ang1)*dcos(ang2)
    
  end subroutine ale_rigid_body_rk4_rot_angs2
   !-----------------------------------------------------------------
  !
  ! Functions for the generalized velocities
  !
  !-----------------------------------------------------------------
  !-----------------------------------------------------------------
  ! Generalized velocity in theta direction
  ! I\dot{\vtheta} = I \vphi^2\sin \theta \cos \theta +
  !                  mgh \cos \theta \cos \phi +
  !                  kh^2(\sin \theta \cos \theta \cos^2 \phi) +
  !                  (F_z \cos\theta-F_x\sin\theta)\,h -
  !                  \alpha_\theta h^2 \vtheta
  !
  real(kind=8) function  f_vtheta_dot(Irbo, Frbo, hcmo, krbo, mrbo, grbo, difrbo, gen_coor1o, gen_coor2o, gen_vel1o, gen_vel2o) result(eval)

    implicit none

    real(kind=8), intent(in)      :: gen_vel1o,gen_coor1o,gen_coor2o,gen_vel2o,hcmo,krbo,difrbo,mrbo,grbo
    real(kind=8), intent(in)      :: Frbo(3), Irbo(3)

    eval =  ( Irbo(1) * gen_vel2o * gen_vel2o * dsin(gen_coor1o) * dcos(gen_coor1o) + &
         & mrbo * grbo * hcmo * dcos(gen_coor1o) * dcos(gen_coor2o) + &
         & krbo * hcmo * hcmo * dsin(gen_coor1o) * dcos(gen_coor1o) * dcos(gen_coor2o) * dcos(gen_coor2o) + &
         & Frbo(3) * dcos(gen_coor1o) * hcmo - Frbo(1) * dsin(gen_coor1o) * hcmo - &
         & difrbo * hcmo * hcmo * gen_vel1o ) &
         & / Irbo(1)

  end function f_vtheta_dot

  !-----------------------------------------------------------------  
  ! Generalized velocity in phi direction
  ! \dot{\vphi} = - ( 2I\sin\theta\cos\theta\,\vtheta\,\vphi -
  !                 mgh \sin \theta \sin \phi +
  !                 kh^2(\sin\phi\cos\phi\sin^2\theta) -
  !                 (F_z\sin\phi+F_y\cos\phi)\,h +
  !                 \alpha_\phi h^2 \vphi ) / I\sin^2\theta \,
  !
  real(kind=8) function  f_vphi_dot(Irbo, Frbo, hcmo, krbo, mrbo, grbo, difrbo, gen_coor1o, gen_coor2o, gen_vel1o, gen_vel2o) result(eval)

    implicit none

    real(kind=8), intent(in)      :: gen_vel1o,gen_coor1o,gen_coor2o,gen_vel2o,hcmo,krbo,difrbo,mrbo,grbo
    real(kind=8), intent(in)      :: Frbo(3), Irbo(3)   

    eval = ( -mrbo * grbo * hcmo * dsin(gen_coor1o) * dsin( gen_coor2o ) - &
         & 2.d0 * Irbo(1) * dsin(gen_coor2o) * dcos(gen_coor2o) * gen_vel1o * gen_vel2o - &
         & Frbo(3) * dsin(gen_coor1o) * hcmo - Frbo(2) * dcos(gen_coor1o) * hcmo - &
         & krbo * hcmo * hcmo * dsin(gen_coor1o) * dcos(gen_coor1o) * dsin(gen_coor2o) * dsin(gen_coor2o) - &
         & difrbo * hcmo * hcmo * gen_vel1o ) &
         &/ (Irbo(1) * dsin(gen_coor2o) * dsin(gen_coor2o) )

  end function f_vphi_dot
  !-----------------------------------------------------------------  
  !
  !  f_valpha_dot = (-mgh cos(alpha) cos(beta) - k * h^2 (alpha-pi/2) - F_x sen(alpha) h + F_z cos(alpha) h) / I
  !
  real(rp) function  f_valpha_dot(Irbo, Frbo, hcmo, krbo, mrbo, grbo, difrbo, gen_coor1o, gen_coor2o, gen_vel1o, gen_vel2o) result(eval)

    implicit none

    real(rp), intent(in)      :: gen_vel1o,gen_coor1o,gen_coor2o,gen_vel2o,hcmo,krbo,difrbo,mrbo,grbo
    real(rp), intent(in)      :: Frbo(3_ip), Irbo(3_ip)
    real(kind=8)              :: disc
    ! real(kind=8)                  :: eval

    disc = 0.d0
    if( dabs(gen_coor1o-pi*0.5) > 0.0379 )disc=0.d0

    eval = (  (Irbo(1_ip) ) * gen_vel2o * gen_vel2o * dsin(gen_coor1o) * dcos(gen_coor1o) - mrbo * grbo * hcmo * dcos(gen_coor1o) * dsin(gen_coor2o) - &
         &krbo * hcmo * hcmo * ( gen_coor1o - pi*0.5_rp ) - Frbo(1_ip) * dsin(gen_coor1o) * hcmo + Frbo(3_ip) * dcos(gen_coor1o) * hcmo - &
         &difrbo * hcmo * hcmo * gen_vel1o  - sign(1.d0,gen_coor1o - pi*0.5)*256.d0*disc*dabs( dabs( gen_coor1o - pi * 0.5 )- 0.0379 )**0.725 ) / (Irbo(1_ip) )


  end function f_valpha_dot

  !-----------------------------------------------------------------  
  !
  !  f_beta_dot = ( mgh sen(alpha) sen(beta) - k * h^2 (beta-pi/2)  + F_y cos(beta) h + F_z sen(beta) h) / I
  !
  real(rp) function  f_vbeta_dot(Irbo, Frbo, hcmo, krbo, mrbo, grbo, difrbo, gen_coor1o, gen_coor2o, gen_vel1o, gen_vel2o) result(eval)

    implicit none

    real(rp), intent(in)      :: gen_vel1o,gen_coor1o,gen_coor2o,gen_vel2o,hcmo,krbo,difrbo,mrbo,grbo
    real(rp), intent(in)      :: Frbo(3_ip), Irbo(3_ip)
    real(kind=8)              :: disc
    ! real(kind=8)                  :: eval

    disc = 0.d0
    if( dabs(gen_coor1o-pi*0.5) > 0.0379 )disc=0.d0

    eval = (- mrbo * grbo * hcmo * dcos(gen_coor1o) * dsin(gen_coor2o) - krbo * hcmo * hcmo * ( gen_coor1o - pi*0.5_rp ) - &
         & Frbo(2_ip) * dsin(gen_coor1o) * hcmo + Frbo(3_ip) * dcos(gen_coor1o) * hcmo - difrbo * hcmo * hcmo * gen_vel1o  - &
         & sign(1.d0,gen_coor1o - pi*0.5)*256.d0*disc*dabs( dabs( gen_coor1o - pi * 0.5 )- 0.0379 )**0.725 ) / (Irbo(1_ip) * dsin(gen_coor2o) * dsin(gen_coor2o))

  end function f_vbeta_dot

  !-----------------------------------------------------------------
  !
  ! f_spring = -kx - c \frac{dx}{dt} + Fx
  !
  real(rp) function  f_spring(Frbo, krbo, mrbo, crbo, gen_coor1o, gen_vel1o) result(eval)

    implicit none
  
    real(rp), intent(in)        :: gen_coor1o, gen_vel1o, krbo, mrbo, crbo 
    real(rp), intent(in)        :: Frbo 
  
    eval = - (krbo * gen_coor1o + crbo * gen_vel1o - Frbo) / mrbo

  end function 

end module mod_ale_rigid_body_rk4
