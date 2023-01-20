!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine ale_sorbhh(dt)
  !-----------------------------------------------------------------------
  !****f* Alefor/ale_sorbhh
  ! NAME
  !    ale_sorbhh
  ! DESCRIPTION
  !    This routines solves the rigid body equations
  ! USED BY
  !    ale_solrbo
  !***
  !-----------------------------------------------------------------------
  use def_master
  use def_domain
  use def_alefor

  implicit none
  real(rp), intent(in) :: dt
  real(rp)             :: tol,numer,deter,denom

  integer(ip)          :: iimbo,idime,maxit,ipoib,itera
  real(rp)             :: onvma,auxir,error
  real(rp)             :: xcoor(3),rotco(3)

  real(rp), pointer    :: force(:,:),accel(:,:),velol(:,:),posil(:,:)             ! Linear motion
  real(rp), pointer    :: torqu(:,:),accea(:,:),veloa(:,:),posia(:,:)             ! Angular motion
  real(rp), pointer    :: quate(:,:),q_dot(:,:),rotac(:,:),IT(:)
  real(rp)             :: Io(3,3),In1(3,3),Iinv(3,3),wiii(3),qtem(4),vtem(3),v2tem(3),Ttem(3,3),rotpr(3,3)
  real(rp)             :: i_aux(3),j_aux(3),r_aux,eulan(3,2),v_aux(3),carax(3,3)
  real(rp)             :: C1,C2

  if ( kfl_disor_ale == 1_ip ) then ! order for the calculation of the displacements
     C1 = 1.0_rp; C2 = 0.0_rp
  else if( kfl_disor_ale == 2_ip ) then
     C1 = 0.5_rp; C2 = 0.5_rp
  else
     call runend ('ALE_SORBHH: incorrect value for kfl_disor_ale')
  end if
  !----------------------------------------------------------------------
  !
  ! Solver Euler equations
  ! 
  !----------------------------------------------------------------------

  do iimbo = 1,nrbod

     onvma = 1.0_rp / rbbou(iimbo)%massa
     force =>  rbbou(iimbo) % force
     torqu =>  rbbou(iimbo) % torqu
     IT => rbbou(iimbo)%momin
     posil => rbbou(iimbo) % posil
     velol => rbbou(iimbo) % velol
     accel => rbbou(iimbo) % accel
     posia => rbbou(iimbo) % posia
     veloa => rbbou(iimbo) % veloa
     accea => rbbou(iimbo) % accea
     quate => rbbou(iimbo) % quate
     q_dot => rbbou(iimbo) % q_dot
     rotac => rbbou(iimbo) % rotac

     !----------------------------------------------------------------------
     !     
     ! Compute linear quantities - accel & velol
     !
     !----------------------------------------------------------------------
     do idime = 1,ndime
        if ( kfl_ralei_ale == 0 ) then          
           accel(idime,1) =              &
                xline_ale(idime) * force(idime,1) * onvma
        else
           if ( ittim <= nstra_ale ) then 
              auxir = 1.0_rp
           else if ( ittim >= nenra_ale ) then 
              auxir = 0.0_rp
           else
              auxir = dble ( nenra_ale - ittim ) / dble ( nenra_ale - nstra_ale )
           end if

           accel(idime,1) = xline_ale(idime) * ( force(idime,1) * onvma - &
                ( auxir * ralei_ale * velol(idime,3) ) )
        end if
        velol(idime,1) = velol(idime,3) + xline_ale(idime) * dt*accel(idime,1)     
     end do
     !----------------------------------------------------------------------
     !     
     ! Compute angular quantities
     !
     !----------------------------------------------------------------------
     if ( ndime == 2 ) then
        !
        ! Compute angular quantities in 2D
        !
        accea(3,1) = xrota_ale(3) * torqu(3,1) / IT(1)  
        veloa(3,1) = veloa(3,3) + xrota_ale(3) * dt*accea(3,1) 

     else if ( ndime == 3 ) then
        tol = 0.000001_rp
        maxit = 100_ip

        ! Initialitation
        Io(1,1)=IT(1); Io(1,2)=IT(4); Io(1,3)=IT(5)
        Io(2,1)=IT(4); Io(2,2)=IT(2); Io(2,3)=IT(6)
        Io(3,1)=IT(5); Io(3,2)=IT(6); Io(3,3)=IT(3)

        wiii(1) = veloa(1,1); wiii(2) = veloa(2,1) ;wiii(3) = veloa(3,1)

        vtem(1) = 0.0_rp; vtem(2) = 0.0_rp; vtem(3) = 0.0_rp ! Actually I do not see much use in initializing this
                                                             ! since it is calcualed inside do while before it is used

        error = 1.0_rp
        itera = 0_ip

        ! Iterate to find the angular quantities
        do while ( error > tol .and. itera < maxit )  

           ! Compute the actual quaternion 
           ! 1. Compute the product the actual quaternion with the actual angular velocity
           qtem(1) =  0.5_rp * ( -wiii(1)*quate(2,1) - wiii(2)*quate(3,1) - wiii(3)*quate(4,1) )
           qtem(2) =  0.5_rp * ( quate(1,1)*wiii(1) + wiii(2)*quate(4,1) - wiii(3)*quate(3,1) )
           qtem(3) =  0.5_rp * ( quate(1,1)*wiii(2) + wiii(3)*quate(2,1) - wiii(1)*quate(4,1) )
           qtem(4) =  0.5_rp * ( quate(1,1)*wiii(3) + wiii(1)*quate(3,1) - wiii(2)*quate(2,1) )
           ! 2. Add this result, multiplied by 0.5, with the quarterion from the previous time step
           numer = 0.0_rp
           do idime=1,ndime+1
              quate(idime,1) = quate(idime,3) + dt * ( C1 * qtem(idime) + C2 * q_dot(idime,3) )              
              numer = numer + quate(idime,1) * quate(idime,1)
           end do

           ! 3. Normalize actual quaternion
           if ( numer > zeror) then
              do idime=1,ndime+1
                 quate(idime,1) = quate(idime,1) / sqrt(numer)
              end do
           end if

           ! Obtain rotation matrix from actual quaternion
           rotac(1,1) = 1_rp - 2_rp*quate(3,1)*quate(3,1) - 2_rp*quate(4,1)*quate(4,1)
           rotac(2,2) = 1_rp - 2_rp*quate(2,1)*quate(2,1) - 2_rp*quate(4,1)*quate(4,1)
           rotac(3,3) = 1_rp - 2_rp*quate(2,1)*quate(2,1) - 2_rp*quate(3,1)*quate(3,1)
           rotac(1,2) =        2_rp*quate(2,1)*quate(3,1) - 2_rp*quate(1,1)*quate(4,1)
           rotac(2,1) =        2_rp*quate(2,1)*quate(3,1) + 2_rp*quate(1,1)*quate(4,1)
           rotac(1,3) =        2_rp*quate(2,1)*quate(4,1) + 2_rp*quate(1,1)*quate(3,1)
           rotac(3,1) =        2_rp*quate(2,1)*quate(4,1) - 2_rp*quate(1,1)*quate(3,1)
           rotac(2,3) =        2_rp*quate(3,1)*quate(4,1) - 2_rp*quate(1,1)*quate(2,1)
           rotac(3,2) =        2_rp*quate(3,1)*quate(4,1) + 2_rp*quate(1,1)*quate(2,1)

           ! Compute the actual inertia tensor
           call mbmab0(Ttem,rotac,Io,3,3,3)
           call mbmabt(In1,Ttem,rotac,3,3,3)

           ! Compute the actual inverse inertia tensor
           call invmtx(In1,Iinv,deter,3)

           ! Compute the actual angular acceleration
           call mbvab0(vtem,In1,wiii,3,3)
           call vecpro(wiii,vtem,v2tem,3)

           do idime = 1,ndime
              vtem(idime) = torqu(idime,1) - v2Tem(idime)                      
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
              veloa(idime,1) = veloa(idime,3) + xrota_ale(idime) * dt * accea(idime,1)   
           end do

           ! Compute the angular velocity error
           numer=0.0_rp
           denom=0.0_rp
           do idime = 1,ndime              
              numer = numer + (veloa(idime,1) - wiii(idime))*(veloa(idime,1) - wiii(idime))
              denom = denom + veloa(idime,1) * veloa(idime,1)
           end do

           if (denom < zeror) then
              error = 0.0_rp
           else
              error = sqrt(numer)/sqrt(denom)
           end if

           ! Actualice the actual angular velocity iteration
           wiii(1) = veloa(1,1); wiii(2) = veloa(2,1) ; wiii(3) = veloa(3,1) 
           vtem(1) = 0_rp; vtem(2) = 0_rp; vtem(3) = 0_rp    ! nor is there any use in this 
           itera = itera + 1_ip

        end do

     end if
     !----------------------------------------------------------------------
     !
     ! Update the particle data
     ! 
     !----------------------------------------------------------------------
     do idime = 1,3        
        posil(idime,1) =  posil(idime,3) + dt * ( C1 * velol(idime,1) + C2 * velol(idime,3) )              
        if (ndime == 2) posia(idime,1) =  posia(idime,3) + dt * ( C1 * veloa(idime,1) + C2 * veloa(idime,3) ) ! 3d case from R             
     end do

     if (ndime == 2) then
        rotac(1,1) =  COS(posia(3,1))
        rotac(1,2) = -SIN(posia(3,1))
        rotac(1,3) =  0.0_rp
        rotac(2,1) =  SIN(posia(3,1))
        rotac(2,2) =  COS(posia(3,1))
        rotac(2,3) =  0.0_rp
        rotac(3,1) =  0.0_rp
        rotac(3,2) =  0.0_rp
        rotac(3,3) =  1.0_rp
     end if
     !
     ! Update coordinates
     !          
     do ipoib = 1,rbbou(iimbo) % npoib
        do idime = 1,ndime
           xcoor(idime) = rbbou(iimbo) % cooin(idime,ipoib)
        end do
        call mbvab0(rotco,rotac,xcoor,3_ip,3_ip)     
        do idime = 1,ndime
           rbbou(iimbo) % cooib(idime,ipoib) = rotco(idime) + posil(idime,1) 
        end do
     end do

     !
     ! obtain angle from R
     ! 
     if (ndime == 3) then
        ! obtain R_n  from q_n
        rotpr(1,1) = 1_rp - 2_rp*quate(3,3)*quate(3,3) - 2_rp*quate(4,3)*quate(4,3)
        rotpr(2,2) = 1_rp - 2_rp*quate(2,3)*quate(2,3) - 2_rp*quate(4,3)*quate(4,3)
        rotpr(3,3) = 1_rp - 2_rp*quate(2,3)*quate(2,3) - 2_rp*quate(3,3)*quate(3,3)
        rotpr(1,2) =        2_rp*quate(2,3)*quate(3,3) - 2_rp*quate(1,3)*quate(4,3)
        rotpr(2,1) =        2_rp*quate(2,3)*quate(3,3) + 2_rp*quate(1,3)*quate(4,3)
        rotpr(1,3) =        2_rp*quate(2,3)*quate(4,3) + 2_rp*quate(1,3)*quate(3,3)
        rotpr(3,1) =        2_rp*quate(2,3)*quate(4,3) - 2_rp*quate(1,3)*quate(3,3)
        rotpr(2,3) =        2_rp*quate(3,3)*quate(4,3) - 2_rp*quate(1,3)*quate(2,3)
        rotpr(3,2) =        2_rp*quate(3,3)*quate(4,3) + 2_rp*quate(1,3)*quate(2,3)
        !
        ! obtain relative rotations
        ! 
        ! i_n+1 is given by rotac(:,1), j_n+1 is given by rotac(:,2) ....
        ! i_n is given by rotpr(:,1)), j_n is given by rotpr(:,2) ....
        ! now using ecs 2.29 - 2.33 Rodrigo Azcueta's thesis
        ! Obtain i_aux
        r_aux = dot_product(rotac(:,1),rotpr(:,3))
        i_aux = r_aux * rotpr(:,3)
        i_aux = rotac(:,1) - i_aux
        r_aux = sqrt ( dot_product(i_aux,i_aux))
        i_aux = i_aux / r_aux
        ! Obtain j_aux
        call vecpro(rotpr(:,3),i_aux,j_aux,3_ip)
        !
        call vecpro(rotpr(:,1),i_aux,v_aux,3_ip)
        r_aux = dot_product(v_aux,rotpr(:,3))
        eulan(1,1) = asin(r_aux)
        !
        call vecpro(i_aux,rotac(:,1),v_aux,3_ip)
        r_aux = dot_product(v_aux,j_aux)
        eulan(2,2) = asin(r_aux)
        !
        call vecpro(j_aux,rotac(:,2),v_aux,3_ip)
        r_aux = dot_product(v_aux,rotac(:,1))
        eulan(3,2) = asin(r_aux)

        !
        ! obtain absolute rotations
        ! 
        ! i_n+1 is given by rotac(:,1), j_n+1 is given by rotac(:,2) ....
        ! now instead of i_n,  j_n , k_n  we use I=(1,0,0), J=(0,1,0), K=(0,0,1), 
        ! as explained in Rodrigo Azcueta's thesis
        ! Obtain i_aux
        carax = 0.0_rp
        do idime = 1,3
           carax(idime,idime) = 1.0_rp
        end do
        
        r_aux = dot_product(rotac(:,1),carax(:,3))
        i_aux = r_aux * carax(:,3)
        i_aux = rotac(:,1) - i_aux
        r_aux = sqrt ( dot_product(i_aux,i_aux))
        i_aux = i_aux / r_aux
        ! Obtain j_aux
        call vecpro(carax(:,3),i_aux,j_aux,3_ip)
        !
        call vecpro(carax(:,1),i_aux,v_aux,3_ip)
        r_aux = dot_product(v_aux,carax(:,3))
        eulan(1,2) = asin(r_aux)
        !
        call vecpro(i_aux,rotac(:,1),v_aux,3_ip)
        r_aux = dot_product(v_aux,j_aux)
        eulan(2,2) = asin(r_aux)
        !
        call vecpro(j_aux,rotac(:,2),v_aux,3_ip)
        r_aux = dot_product(v_aux,rotac(:,1))
        eulan(3,2) = asin(r_aux)

        if (kfl_paral==1) print*, 'ale_sorbhh:eulan_relat',eulan(:,1)
        if (kfl_paral==1) print*, 'ale_sorbhh:eulan_absol',eulan(:,2)
        posia(:,1) = eulan(:,2)

     end if
  end do

end subroutine ale_sorbhh
