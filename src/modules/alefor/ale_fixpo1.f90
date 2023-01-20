!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine ale_fixpo1(dt,iimbo,tol,maxit)
  !-----------------------------------------------------------------------
  !****f* Alefor/ale_fixpo1
  ! NAME
  !    ale_fixpo1
  ! DESCRIPTION
  !    This routines calculate the angular quantities for a rigid solid
  !    with fixed point iterative method. Full coupled algorithm.   
  ! USED BY
  !    ale_sorbcr
  !***
  !-----------------------------------------------------------------------
  use def_master
  use def_domain
  use def_alefor
  implicit none
  real(rp),    intent(in) :: dt
  integer(ip), intent(in) :: iimbo,maxit
  real(rp),    intent(in) :: tol
  integer(ip)             :: idime,itera
  real(rp)                :: deter,numer,denom,auxir
  real(rp)                :: Taver(3)
  real(rp),    pointer    :: torqu(:,:),accea(:,:),veloa(:,:),IT(:),rotac(:,:),quate(:,:) ! Angular motion
  real(rp)                :: Io(3,3),In1(3,3),Iinv(3,3),wi(3),qtem(4),vtem(3),v2tem(3),Ttem(3,3),error
  real(rp)                :: gamma_ale

  gamma_ale = 1.5_rp

  torqu => rbbou(iimbo)%torqu
  IT => rbbou(iimbo)%momin
  accea  => rbbou(iimbo)%accea
  veloa  => rbbou(iimbo)%veloa
  rotac  => rbbou(iimbo)%rotac
  quate => rbbou(iimbo)%quate

  Taver(3) = 0.0_rp
  do idime = 1,3
     Taver(idime) = torqu(idime,3)
  end do
  ! Initialitation
  Io(1,1)=IT(1); Io(1,2)=IT(4); Io(1,3)=IT(5)
  Io(2,1)=IT(4); Io(2,2)=IT(2); Io(2,3)=IT(6)
  Io(3,1)=IT(5); Io(3,2)=IT(6); Io(3,3)=IT(3)

  wi(1)=veloa(1,1); wi(2)=veloa(2,1) ;wi(3)=veloa(3,1)

  vtem(1)=0.0_rp; vtem(2)=0.0_rp; vtem(3)=0.0_rp

  error = 1.0_rp
  itera = 0_ip

  ! Iterate to find the angular quantities
  do while (error > tol .and. itera < maxit)  
  
     ! Compute the actual quaternion 
     ! 1. Compute the product the actual quaternion with the actual angular velocity
     qtem(1)= -wi(1)*quate(2,1) - wi(2)*quate(3,1) - wi(3)*quate(4,1)
     qtem(2)=  quate(1,1)*wi(1) + wi(2)*quate(4,1) - wi(3)*quate(3,1)
     qtem(3)=  quate(1,1)*wi(2) + wi(3)*quate(2,1) - wi(1)*quate(4,1)
     qtem(4)=  quate(1,1)*wi(3) + wi(1)*quate(3,1) - wi(2)*quate(2,1)
     ! 2. Add this result, multiplied by 0.5, with the quarterion from the previous time step
     numer = 0.0_rp
     do idime=1,ndime+1
        quate(idime,1) = quate(idime,3) + 0.5_rp * dt * qtem(idime)              
        numer = numer + quate(idime,1) * quate(idime,1)
     end do     

     ! 3. Normalize actual quaternion
     if ( numer > zeror) then
        do idime=1,ndime+1
           quate(idime,1) = quate(idime,1) / sqrt(numer)
        end do
     end if

     ! Obtain rotation matrix from actual quaternion
     rotac(1,1)= 1_rp - 2_rp*quate(3,1)*quate(3,1) - 2_rp*quate(4,1)*quate(4,1)
     rotac(2,2)= 1_rp - 2_rp*quate(2,1)*quate(2,1) - 2_rp*quate(4,1)*quate(4,1)
     rotac(3,3)= 1_rp - 2_rp*quate(2,1)*quate(2,1) - 2_rp*quate(3,1)*quate(3,1)
     rotac(1,2)=        2_rp*quate(2,1)*quate(3,1) - 2_rp*quate(1,1)*quate(4,1)
     rotac(2,1)=        2_rp*quate(2,1)*quate(3,1) + 2_rp*quate(1,1)*quate(4,1)
     rotac(1,3)=        2_rp*quate(2,1)*quate(4,1) + 2_rp*quate(1,1)*quate(3,1)
     rotac(3,1)=        2_rp*quate(2,1)*quate(4,1) - 2_rp*quate(1,1)*quate(3,1)
     rotac(2,3)=        2_rp*quate(3,1)*quate(4,1) - 2_rp*quate(1,1)*quate(2,1)
     rotac(3,2)=        2_rp*quate(3,1)*quate(4,1) + 2_rp*quate(1,1)*quate(2,1)

     ! Compute the actual inertia tensor
     call mbmab0(Ttem,rotac,Io,3,3,3)
     call mbmabt(In1,Ttem,rotac,3,3,3)

     ! Compute the actual inverse inertia tensor
     call invmtx(In1,Iinv,deter,3)

     ! Compute the actual angular acceleration
     call mbvab0(vtem,In1,wi,3,3)
     call vecpro(wi,vTem,v2tem,3)

     do idime = 1,ndime
        vtem(idime) = Taver(idime) - v2Tem(idime)                      
     end do

     ! Compute the actual angular quantities           
     do idime = 1,ndime 
        if ( kfl_ralei_ale == 0 ) then          
           accea(idime,1) = xrota_ale(idime) * ( Iinv(idime,1)*vTem(1) + Iinv(idime,2)*vTem(2) + Iinv(idime,3)*vTem(3) )
        else
           
           if ( ittim <= nstra_ale ) then 
              auxir = 1.0_rp
           else if ( ittim >= nenra_ale ) then 
              auxir = 0.0_rp
           else
              auxir = dble ( nenra_ale - ittim ) / dble ( nenra_ale - nstra_ale )
           end if

           accea(idime,1) = xrota_ale(idime) * ( Iinv(idime,1)*vTem(1) + Iinv(idime,2)*vTem(2) + Iinv(idime,3)*vTem(3) &
                - auxir * ralei_ale * veloa(idime,3) )

        end if

        if( ittim == 1 ) accea(idime,3) = accea(idime,1)
        veloa(idime,1) = veloa(idime,3) &
             + xrota_ale(idime) * ( dt*accea(idime,3) + gamma_ale*dt*(accea(idime,1) - accea(idime,3))   )   

     end do

     ! Compute the angular velocity error
     numer=0.0_rp
     denom=0.0_rp
     do idime = 1,ndime              
        numer = numer + (veloa(idime,1) - wi(idime))*(veloa(idime,1) - wi(idime))
        denom = denom + veloa(idime,1)*veloa(idime,1)
     end do

     if (denom < zeror) then
        error = 0.0_rp
     else
        error = sqrt(numer)/sqrt(denom)
     end if

     ! Compute the actual angular velocity iteration
     wi(1)=veloa(1,1); wi(2)=veloa(2,1) ;wi(3)=veloa(3,1) 

     vtem(1)=0_rp; vtem(2)=0_rp; vtem(3)=0_rp     
     itera=itera+1_ip

  end do

end subroutine ale_fixpo1

