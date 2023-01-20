!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine ale_forces()
  !-----------------------------------------------------------------------
  !****f* Alefor/ale_forces
  ! NAME
  !    ale_forces
  ! DESCRIPTION
  !    This routines obtains forces for rigid bodies ! similar to ibm_forces
  ! USED BY
  !    domain
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_kermod
  use def_domain
  use def_alefor
  implicit none
  integer(ip)       :: iimbo,idime,indau
  real(rp)          :: grafo
  real(rp), pointer :: force(:,:),vpfor(:,:),  xi(:,:), vi(:,:)      ! Linear motion
  real(rp), pointer :: torqu(:,:),vptor(:,:)     ! Angular motion

  real(rp)          :: xxxxx(3),F1(3),T1(3)

  if ( kfl_crist_ale == 1 ) then 
     indau = 3_ip   ! lo de cristobal hay que revisarlo un poco más
  else
     indau = 1_ip
  end if
  !----------------------------------------------------------------------
  !
  ! Initialize
  ! 
  !----------------------------------------------------------------------

  do iimbo = 1,nrbod
     force  => rbbou(iimbo) % force
     torqu  => rbbou(iimbo) % torqu

     do idime = 1,3
        force(idime,indau) = 0.0_rp
        torqu(idime,indau) = 0.0_rp
     end do
  end do

  !----------------------------------------------------------------------
  !
  ! Add force and torque due to coupling with other modules
  ! 
  !----------------------------------------------------------------------


  if (  coupling('ALEFOR','NASTIN') >= 1_ip  ) then

    do iimbo = 1,nrbod
        force   => rbbou(iimbo) % force
        vpfor  => rbbou(iimbo) % vpfor
        torqu   => rbbou(iimbo) % torqu
        vptor  => rbbou(iimbo) % vptor

        do idime = 1,3
           force(idime,indau) = force(idime,indau) + vpfor(idime,1) 
           torqu(idime,indau) = 0.0_rp ! CHANGE !!!  torqu(idime,indau) + vptor(idime,1) 

        end do

     end do

  end if

  !----------------------------------------------------------------------
  !
  ! Gravity and Buoyancy
  ! 
  !----------------------------------------------------------------------

  grafo = real( kfl_grafo_ale , rp )  ! Gravity  force = 1.0

  do iimbo = 1,nrbod
     force  => rbbou(iimbo) % force
     xi => rbbou(iimbo) % posil
     vi => rbbou(iimbo) % velol
     do idime = 1,ndime
        force(idime,indau) = force(idime,indau) + grnor * gravi(idime) &
             &     * rbbou(iimbo) % volum &
             &     * ( grafo * rbbou(iimbo) % densi )

!          force(idime,indau) = force(idime,indau) - 1.0e3 * xi(idime,indau) &
!                               & - 0.0_rp * vi(idime,indau)

     end do

  end do

  !----------------------------------------------------------------------
  !
  ! Catamaran force
  ! 
  !----------------------------------------------------------------------

  if( kfl_catfo_ale /= 0 ) then ! when I move it to ale put it more general and erase this
     
     force   => rbbou(1) % force
     torqu   => rbbou(1) % torqu
     F1(1) = 775.0_rp   ! without - because we have x and y axis oposite to CHE
     F1(2) = 1600.0_rp  ! without - because we have x and y axis oposite to CHE         
     F1(3) = 200.0_rp
     xxxxx(1) =  0.0_rp   ! x -X_cog
     xxxxx(2) =  0.0_rp
     xxxxx(3) =  7.0_rp

     call vecpro(xxxxx,F1,T1,3_ip)                  ! T1 = (X) x F1  (pressure) 

     do idime = 1,3
        force(idime,indau) = force(idime,indau) + F1(idime)
        torqu(idime,indau) = torqu(idime,indau) + T1(idime)
     end do

     force(3,:) = 0.0_rp !!!!CHANGE
     
  else if( kfl_sprin_ale /= 0 )then
     !----------------------------------------------------------------------
     !
     ! Damped spring forces
     !
     !----------------------------------------------------------------------
     do iimbo = 1,nrbod
        force  => rbbou(iimbo) % force
        xi     => rbbou(iimbo) % posil
        vi     => rbbou(iimbo) % velol
        do idime = 1,ndime        
           force(idime,indau) = sprin_ale(3_ip)*force(idime,indau) - sprin_ale(1_ip) *xi(idime,indau) &
                & - sprin_ale(2_ip) * vi(idime,indau)    
        end do
     end do

  end if

end subroutine ale_forces
