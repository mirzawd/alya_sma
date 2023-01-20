!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine ale_sorbcr(dt)
  !-----------------------------------------------------------------------
  !****f* Alefor/ale_sorbcr
  ! NAME
  !    ale_sorbcr
  ! DESCRIPTION
  !    This routines solves the Euler and Newton equations for rigid bodies using what cristobal
  !    was using in ibm_eulers and ibm_newmak
  !    La idea es luego eliminar esta sub
  ! USED BY
  !    domain
  !***
  !-----------------------------------------------------------------------
  use def_master
  use def_domain
  use def_alefor

  implicit none
  real(rp), intent(in) :: dt
  integer(ip)          :: iimbo,idime
  real(rp)             :: onvma,auxir
  real(rp)             :: Faver(3),Taver(3)

  real(rp), pointer    :: force(:,:),accel(:,:),velol(:,:)        ! Linear motion
  real(rp), pointer    :: torqu(:,:),accea(:,:),IT(:),rotac(:,:)  ! Angular motion


! lo de ibm_newmak
  integer(ip)           :: ipoib
  real(rp)              :: numer
  real(rp)              :: qtem(4),xcoor(3),rotco(3)
  real(rp), pointer     :: ai(:,:),vi(:,:),xi(:,:),zi(:,:),wi(:,:),si(:,:),qi(:,:),Ri(:,:)
  real(rp)              :: gamma_ale,beta_ale

  gamma_ale = 1.5_rp
  beta_ale = 1.0_rp
  !----------------------------------------------------------------------
  !
  ! Solver Euler equations
  ! 
  !----------------------------------------------------------------------

  do iimbo = 1,nrbod

     onvma = 1.0_rp / rbbou(iimbo)%massa

     force =>  rbbou(iimbo) % force
     accel =>  rbbou(iimbo) % accel
     velol =>  rbbou(iimbo) % velol


     torqu =>  rbbou(iimbo) % torqu
     accea =>  rbbou(iimbo) % accea

     
     if (kfl_nforc_ale == 1) then
        do idime = 1,3
           Faver(idime) = force(idime,3)
           Taver(idime) = torqu(idime,3)
        end do
     elseif (kfl_nforc_ale == 2) then
        do idime = 1,3
           Faver(idime) = 0.5_rp*(force(idime,3) + force(idime,4))
           Taver(idime) = 0.5_rp*(torqu(idime,3) + torqu(idime,4))
        end do
     end if
     
     !----------------------------------------------------------------------
     !     
     ! Compute linear quantities
     !
     !----------------------------------------------------------------------     
     do idime = 1,ndime
        if ( kfl_ralei_ale == 0 ) then          
           accel(idime,1) =              &
                xline_ale(idime) * Faver(idime) * onvma
        else

           if ( ittim <= nstra_ale ) then 
              auxir = 1.0_rp
           else if ( ittim >= nenra_ale ) then 
              auxir = 0.0_rp
           else
              auxir = dble ( nenra_ale - ittim ) / dble ( nenra_ale - nstra_ale )
           end if

           accel(idime,1) =              &
                xline_ale(idime) * ( Faver(idime) * onvma - &
                ( auxir * ralei_ale * velol(idime,3) ) )
        end if

     end do
     !----------------------------------------------------------------------
     !     
     ! Compute angular quantities
     !
     !----------------------------------------------------------------------
     rotac => rbbou(iimbo) % rotac
     IT => rbbou(iimbo)%momin
     if ( ndime == 2 ) then
        !
        ! Compute angular quantities in 2D
        !
        accea(3,1) = xrota_ale(3) * Taver(3) / IT(1)        
     else if ( ndime == 3 ) then
        call ale_fixpo1(dt,iimbo,0.000001_rp,100_ip)
     end if
  end do
  !
  ! Aca mismo meto de de ibm_newmak
  !
  !-----------------------------------------------------------------------
  !****f* ibm_newmak/ibm_newmak
  ! NAME
  !    ibm_newmak
  ! DESCRIPTION
  !    This routines detect a collision between particles A PRIORI
  !                             gamma   beta
  !     Fox Goodwin             1 / 2   1 / 12    conditionnaly stable
  !     Linear acceleration     1 / 2   1 / 6     conditionnaly stable
  !     Average acceleration    1 / 2   1 / 4     inconditionnaly stable
  !     Nosotros                0.75    0.390625  diffusive
  !     Nosotros                1.0     0.5625    super diffusive
  !     External                1 / 2   0         explicit 
  !
  !     Relation between beta and gamma: beta = 0.25*(gamma+1/2)^2
  !
  ! USED BY
  !***
  !-----------------------------------------------------------------------



  !----------------------------------------------------------------------
  !
  ! Update the particle data
  ! 
  !----------------------------------------------------------------------
  do iimbo = 1,nrbod  
     xi => rbbou(iimbo) % posil
     vi => rbbou(iimbo) % velol
     ai => rbbou(iimbo) % accel
     si => rbbou(iimbo) % posia
     wi => rbbou(iimbo) % veloa
     zi => rbbou(iimbo) % accea
     qi => rbbou(iimbo) % quate
     Ri => rbbou(iimbo) % rotac

     do idime = 1,3        
        xi(idime,1) = xi(idime,3) + xline_ale(idime) * (dt*vi(idime,3) + 0.5_rp*dt*dt*ai(idime,3) + &
             beta_ale*dt*dt*( ai(idime,1) - ai(idime,3)))
        vi(idime,1) = vi(idime,3) + xline_ale(idime) * (dt*ai(idime,3) + gamma_ale*dt*(ai(idime,1) - ai(idime,3)))

        si(idime,1) = si(idime,3) + xrota_ale(idime) * (dt*wi(idime,3) + 0.5_rp*dt*dt*zi(idime,3) + &
             beta_ale*dt*dt*( zi(idime,1) - zi(idime,3)))
        wi(idime,1) = wi(idime,3) + xrota_ale(idime) * (dt*zi(idime,3) + gamma_ale*dt*(zi(idime,1) - zi(idime,3)))

     end do

     if (ndime == 2) then
        Ri(1,1) =  COS(si(3,1))
        Ri(1,2) = -SIN(si(3,1))
        Ri(1,3) =  0.0_rp
        Ri(2,1) =  SIN(si(3,1))
        Ri(2,2) =  COS(si(3,1))
        Ri(2,3) =  0.0_rp
        Ri(3,1) =  0.0_rp
        Ri(3,2) =  0.0_rp
        Ri(3,3) =  1.0_rp
     else
        !
        ! 2. Add this result, multiplied by 0.5, with the quarterion from the previous time step
        !
        numer=0.0_rp

        qtem(1) = qi(1,2) + 0.5_rp * dt *( -wi(1,1)*qi(2,1) - wi(2,1)*qi(3,1) - wi(3,1)*qi(4,1))
        qtem(2) = qi(2,2) + 0.5_rp * dt *(  qi(1,1)*wi(1,1) + wi(2,1)*qi(4,1) - wi(3,1)*qi(3,1))
        qtem(3) = qi(3,2) + 0.5_rp * dt *(  qi(1,1)*wi(2,1) + wi(3,1)*qi(2,1) - wi(1,1)*qi(4,1))
        qtem(4) = qi(4,2) + 0.5_rp * dt *(  qi(1,1)*wi(3,1) + wi(1,1)*qi(3,1) - wi(2,1)*qi(2,1))

        numer   = qtem(1)*qtem(1) + qtem(2)*qtem(2) + qtem(3)*qtem(3) + qtem(4)*qtem(4)
        !
        ! 3. Normalize actual quaternion
        !
        if (numer/=0.0_rp) then
           numer = sqrt(numer)
           qtem(1) = qtem(1) / numer
           qtem(2) = qtem(2) / numer
           qtem(3) = qtem(3) / numer
           qtem(4) = qtem(4) / numer
        end if
        do idime=1,ndime+1
           qi(idime,1)= qtem(idime)
        end do
        Ri(1,1)= 1_rp - 2_rp*qi(3,1)*qi(3,1) - 2_rp*qi(4,1)*qi(4,1)
        Ri(2,2)= 1_rp - 2_rp*qi(2,1)*qi(2,1) - 2_rp*qi(4,1)*qi(4,1)
        Ri(3,3)= 1_rp - 2_rp*qi(2,1)*qi(2,1) - 2_rp*qi(3,1)*qi(3,1)
        Ri(1,2)=        2_rp*qi(2,1)*qi(3,1) - 2_rp*qi(1,1)*qi(4,1)
        Ri(2,1)=        2_rp*qi(2,1)*qi(3,1) + 2_rp*qi(1,1)*qi(4,1)
        Ri(1,3)=        2_rp*qi(2,1)*qi(4,1) + 2_rp*qi(1,1)*qi(3,1)
        Ri(3,1)=        2_rp*qi(2,1)*qi(4,1) - 2_rp*qi(1,1)*qi(3,1)
        Ri(2,3)=        2_rp*qi(3,1)*qi(4,1) - 2_rp*qi(1,1)*qi(2,1)
        Ri(3,2)=        2_rp*qi(3,1)*qi(4,1) + 2_rp*qi(1,1)*qi(2,1)
     end if
     !
     ! Update coordinates
     !          
     do ipoib = 1,rbbou(iimbo) % npoib
        do idime = 1,ndime
           xcoor(idime) = rbbou(iimbo) % cooin(idime,ipoib)
        end do
        call mbvab0(rotco,Ri,xcoor,3_ip,3_ip)     
        do idime = 1,ndime
           rbbou(iimbo) % cooib(idime,ipoib) = rotco(idime) + xi(idime,1) 
        end do
     end do

  end do

end subroutine ale_sorbcr
