!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine ale_sorbjc(dt)
  !-----------------------------------------------------------------------
  !****f* Alefor/ale_sorbhh
  ! NAME
  !    ale_sorbjc
  ! DESCRIPTION
  !    This routines solves the rigid body equations in generalized coordinates
  ! USED BY
  !    ale_solrbo
  !***
  !-----------------------------------------------------------------------
  use def_master
  use def_domain
  use def_alefor
  use def_parame,              only : pi
  use def_kermod,              only : grnor 
  use mod_ale_rigid_body_rk4,  only : ale_rigid_body_rk4_solution,  ale_rigid_body_rk4_ang_trans1
  use mod_ale_rigid_body_rk4,  only : ale_rigid_body_rk4_rot_angs1, ale_rigid_body_rk4_rot_angs2
  use mod_communications,      only : PAR_MAX

  implicit none

  real(rp), intent(in) :: dt
  integer(ip)          :: iimbo,idime,ipoib
  real(rp)             :: mass, ITrb(3_ip)
  real(rp)             :: xcoor(3_ip),rotco(3_ip)
  real(rp)             :: force_rb(3_ip) 
  real(rp)             :: rot_ang1,rot_ang2 
  real(rp)             :: hcm, krb, difrb
  real(rp), pointer    :: vpfor(:,:)
  real(rp), pointer    :: veloa(:,:), posia(:,:)            
  real(rp), pointer    :: rotac(:,:), IT(:)

  do iimbo = 1,nrbod

     mass  =  rbbou(iimbo) % massa       ! Mass of the rigid body
     IT    => rbbou(iimbo) % momin       ! Inertia tensor (here we use the IT of the moving frame)
     posia => rbbou(iimbo) % posia       ! Angular position of the body, here they are the euler angles
     veloa => rbbou(iimbo) % veloa       ! Angular velocity of the body related to the euler angles
     rotac => rbbou(iimbo) % rotac       ! Rotation matrix
     vpfor => rbbou(iimbo) % vpfor       ! viscous and pressure forces

     krb   =  sprin_ale(1_ip)
     difrb =  sprin_ale(2_ip)
     hcm   =  sprin_ale(3_ip)

     do idime = 1_ip, 3_ip

        ITrb(idime)     = IT(idime)
        force_rb(idime) = vpfor(idime,1_ip)

     end do

     !----------------------------------------------------------------------
     !     
     ! Solve the equations of motion in generalized coordinates
     !
     !----------------------------------------------------------------------

     if( kfl_topmo_ale == 1_ip ) then       

        if( kfl_pertu_ale == 0_ip ) then
           !
           ! Rigid body motion with non-spinning top model
           !                    
           call ale_rigid_body_rk4_solution(posia(1_ip,1_ip), posia(2_ip,1_ip), veloa(1_ip,1_ip), veloa(2_ip,1_ip),&
                & dt, ITrb, force_rb, hcm, krb, difrb, mass, grnor, 'TOPMO')
           
        else if( kfl_pertu_ale == 1_ip ) then
           !
           ! Perturbation of the rigid body motion with non-spinning top model
           ! 
           call ale_rigid_body_rk4_solution(posia(1_ip,1_ip), posia(2_ip,1_ip), veloa(1_ip,1_ip), veloa(2_ip,1_ip),&
                & dt, ITrb, force_rb, hcm, krb, difrb, mass, grnor, 'PERTU')             
        end if

     else if( kfl_sprin_ale == 1_ip )then
        !
        ! Rigid body motion with linear spring model
        ! 
        call ale_rigid_body_rk4_solution(posia(1_ip,1_ip), posia(2_ip,1_ip), veloa(1_ip,1_ip), veloa(2_ip,1_ip),&
             & dt, ITrb, force_rb, hcm, krb, difrb, mass, grnor, 'SPRIN')
     end if

     !----------------------------------------------------------------------
     !     
     ! Exchange the result in order to have exactly the same value in all
     ! the MPIs
     !
     !----------------------------------------------------------------------

     do idime = 1_ip, ndime
        call PAR_MAX( posia(idime,1_ip) )
        call PAR_MAX( veloa(idime,1_ip) )
        call PAR_MAX( vpfor(idime,1_ip) )
        call PAR_MAX( vpfor(idime,2_ip) )
        call PAR_MAX( vpfor(idime,3_ip) )
        call PAR_MAX( vpfor(idime,4_ip) )
     end do

     !----------------------------------------------------------------------
     !
     ! Postprocess of the solution ( if necessary )
     !
     !----------------------------------------------------------------------

     if( kfl_topmo_ale == 1_ip ) then
        
        !-------------------------------------------------------------------
        !
        ! Rigid body motion with non-spinning top model
        ! Transform the angles and rotate the rigid body.
        !
        ! IMPORTANT:  The rotation matrix is calculated taking the x axis
        ! as the original axis of the rigid body. For example if the body
        ! is aligned with the z axis, it is necessary to calculate
        ! rot_ang1 = theta - pi / 2
        !-------------------------------------------------------------------
        rot_ang1 = posia(1_ip,1_ip)-0.5_rp*pi !   0.5_rp*pi - posia(1_ip,1_ip)
        rot_ang2 = posia(2_ip,1_ip) !-( 0.5_rp*pi - posia(2_ip,1_ip) )
        !--------------------------------------------------------------------
        !
        ! Calculate the rotation matrix        
        ! 
        !--------------------------------------------------------------------
        call ale_rigid_body_rk4_rot_angs1(rotac,rot_ang1,rot_ang2)
        !--------------------------------------------------------------------
        !
        ! Update coordinates
        !
        !--------------------------------------------------------------------
        do ipoib = 1,rbbou(iimbo) % npoib

           do idime = 1,ndime
              xcoor(idime) = rbbou(iimbo) % cooin(idime,ipoib)
           end do

           call mbvab0(rotco,rotac,xcoor,3_ip,3_ip)     
           do idime = 1,ndime
              rbbou(iimbo) % cooib(idime,ipoib) = rotco(idime)
           end do

        end do

     else if( kfl_sprin_ale == 1_ip )then
        !-------------------------------------------------------------------
        !
        ! Rigid body motion with linear spring model
        ! Translate the rigid body
        !-------------------------------------------------------------------

        !-------------------------------------------------------------------
        !
        ! Update coordinates
        !
        !-------------------------------------------------------------------
        do ipoib = 1,rbbou(iimbo) % npoib

           do idime = 1,ndime

              rbbou(iimbo) % cooib(idime,ipoib) = rbbou(iimbo) % cooin(idime,ipoib) + posia(idime, 1_ip)

           end do

        end do
     end if

  end do   ! nrbod


end subroutine ale_sorbjc

