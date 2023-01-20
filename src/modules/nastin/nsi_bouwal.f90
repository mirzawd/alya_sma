!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_bouwal(&
     itask,pevat,pnodb,ndofn,iboun,lboel,gbsha,bovel,bovfi,tract,gbvis,&
     gbden,baloc,velfr,wmatr,rough,kinen,velex,tveav,igaub,lelbo)

  !-----------------------------------------------------------------------
  !****f* Nastin/nsi_bouwal
  ! NAME 
  !    nsi_bouwal
  ! DESCRIPTION
  !    This routine computes the surface traction for the NS equations at
  !    a given integration point of a boundary IBOUN received by argument
  !    due to the use of a turbulent wall law. The algorithm is:
  !    - Compute the tangential velocity u at the integration point x.
  !      In fact, u is not necessarily tangential at x. For example:
  !      -> u    -> u      
  !      o---x---o---x---o
  !                      |\ u
  !                      | 
  !    - Given the distance to the wall y, compute U*
  !    - Compute y+=yU*/nu
  !      if(y+>5) then
  !        t=-rho*U*^2*(u_tan-u_fix_tan)/|u_tan-u_fix_tan|
  !      else
  !        u+=y+ => U*^2=u*nu/y so that
  !        t=-mu*u/y
  !      end if
  ! USES
  !    vecnor
  !    frivel
  ! USED BY
  !    nsi_bouope
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only      :  ip,rp
  use mod_run_config
  use def_kermod, only      :  kfl_delta,usref,kfl_ustar,cmu_st
  use def_kermod, only      :  dexlo_ker,kfl_waexl_ker,kfl_waexl_imp_ker, kfl_mlwm_ker

  use def_kermod, only      :  shape_waexl,lexlo_ker,kfl_wlaav_ker,kfl_aveme_ker
  use def_domain, only      :  ndime,ywalb,ltype,nnode
  use def_nastin, only      :  zensi,delta_nsi
  use mod_ker_regularization, only : regul_k, regul_e, kfl_regularization
  use mod_nsi_frixgb, only  : nsi_ml_ustar_all,ustarb     ! machine learning
  use mod_frivel,     only  : frivel

  implicit none
  integer(ip), intent(in)  :: itask
  integer(ip), intent(in)  :: pevat
  integer(ip), intent(in)  :: pnodb
  integer(ip), intent(in)  :: ndofn
  integer(ip), intent(in)  :: iboun
  integer(ip), intent(in)  :: lboel(pnodb)
  real(rp),    intent(in)  :: gbsha(pnodb)
  real(rp),    intent(in)  :: bovel(ndime,pnodb)
  real(rp),    intent(in)  :: bovfi(ndime,pnodb)
  real(rp),    intent(inout) :: tract(ndime)
  real(rp),    intent(in)  :: gbvis
  real(rp),    intent(in)  :: gbden
  real(rp),    intent(in)  :: baloc(ndime,ndime)
  real(rp),    intent(out) :: velfr
  real(rp),    intent(inout) :: wmatr(pevat,pevat)

  real(rp),    intent(in)  :: rough
  real(rp),    intent(in)  :: kinen
  real(rp),    intent(in)  :: velex(ndime)   ! velocity at the exchange location - only needed for  kfl_waexl_nsi /= 0
  real(rp),    intent(in)  :: tveav(ndime)   ! Average velocity for wall law  - kfl_wlaav_ker = 1
  integer(ip), intent(in)  :: igaub
  integer(ip), intent(in)  :: lelbo

  integer(ip)              :: ldime,idime,inodb,idofn
  integer(ip)              :: jnode,pelty,pnode
  integer(ip)              :: jevab,ievab,jdime,jnodb
  integer(ip)              :: imeth,ielem,indei
  real(rp)                 :: vewal(3),fact1,fact2,velfi(3)
  real(rp)                 :: vikin,yplus                      !  nu, U*, y+, y
  real(rp)                 :: tveno,tven2                      ! |u_tan-u_fix_tan|
  real(rp)                 :: tvelo(3)                         !  u_tan
  real(rp)                 :: tvefi(3)                         !  u_fix_tan  - it comes from bvess_nsi put in global system
  real(rp)                 :: tvedi(3)                         !  u_tan - u_fix_tan
  real(rp)                 :: velsh(ndime,pnodb*ndime)
  real(rp)                 :: vels2(ndime,ndime,pevat)
  real(rp)                 :: delta_aux
  !
  ! Base solution
  !
  real(rp)                ::  ub(ndime)
  real(rp)                ::  velf2,umbra,avvfr
  !
  !
  !
  real(rp)                ::  h1,h2,hs, reguk
!  real(rp)                ::  u1,u2,u3
  !  average method 1) instantaneous direction  2) average direction   3) all average
                                         ! = 1) tau = -C mod(avev)*veloc,
                                         ! = 2) tau = -C mod(veloc)*avev
                                         ! = 3) tau = -C mod(avev)avev
  

  if ( kfl_waexl_ker == 0_ip ) then  ! normal behaviour
     if( kfl_delta == 1 ) then
        delta_aux = ywalb(iboun)                                  ! variable wall distance
     else
        delta_aux = delta_nsi                                     ! fixed wall distance
     end if
  else      ! using exchange location
     if( kfl_delta == 1 ) then
        delta_aux = ywalb(iboun)                                  ! variable wall distance
     else
        delta_aux = dexlo_ker                                     ! fixed wall distance
     end if
  end if
  
  imeth = 2
  if ( kfl_waexl_ker == 1_ip.and.kfl_waexl_imp_ker==0 )  imeth = 1     ! exchange location is done explicitly (default)

  umbra = 1.0e-6_rp
   
  if( ( delta_nsi > zensi ) .or. ( kfl_delta == 1 ) ) then

     if ( kfl_waexl_ker == 0_ip ) then  ! normal behaviour

        do idime = 1,ndime                                        ! Velocity
           vewal(idime) = 0.0_rp
           velfi(idime) = 0.0_rp                                  ! velocity of the wall, due to moving domain
        end do
        do inodb = 1,pnodb
           do idime = 1,ndime            
              vewal(idime) = vewal(idime) &
                   + gbsha(inodb) * bovel(idime,inodb)
              velfi(idime) = velfi(idime) &
                   + gbsha(inodb) * bovfi(idime,inodb)
           end do
        end do

     else  ! using velocity from exchange location
        do idime = 1,ndime
           vewal(idime) = velex(idime)                            ! Velocity - this comes as an input to this subroutine
           velfi(idime) = 0.0_rp                                  ! velocity of the wall, due to moving domain
        end do
        do inodb = 1,pnodb
           do idime = 1,ndime            
              velfi(idime) = velfi(idime) &
                   + gbsha(inodb) * bovfi(idime,inodb)
           end do
        end do

     end if

     if( imeth == 1 ) then

        !-----------------------------------------------------------------
        !
        ! Explicit treatment: Compute u_tan and |u_tan|
        !
        !-----------------------------------------------------------------
        !
        ! Tangent velocity (TVELO), tangent component of the prescribed velocity (TVEFI) and TVENO = |TVELO-TVEFI| 
        !
        do idime = 1,ndime     
           tvelo(idime) = vewal(idime)
           tvefi(idime) = velfi(idime)
           do ldime = 1,ndime
              tvelo(idime) = tvelo(idime)   &
                   - baloc(idime,ndime) &
                   * baloc(ldime,ndime) * vewal(ldime)
              tvefi(idime) = tvefi(idime)   &
                   - baloc(idime,ndime) &
                   * baloc(ldime,ndime) * velfi(ldime)
           end do
           tvedi(idime) = tvelo(idime) - tvefi(idime)
        end do
        if (kfl_aveme_ker==3) tvedi(1:ndime) = tveav(1:ndime)
        !
        ! Compute U*: VELFR
        !
        vikin = gbvis/gbden                                    ! nu
        call vecnor(tvedi,ndime,tveno,2_ip)                    ! |u_tan-u_fix_tan|
		! Coding here for machine learning interfacing
		if (kfl_mlwm_ker == 1_ip) then
            ustarb = tveno/velfr
			velfr  = ustarb 
		else  
        	call frivel(kfl_ustar,delta_aux,rough,tveno,vikin,velfr)         ! U*
		end if
        if (kfl_wlaav_ker == 1_ip.and.kfl_aveme_ker/=2) then  ! if averaged velocity for wall law
           call vecnor(tveav,ndime,tven2,2_ip)               
           call frivel(kfl_ustar,delta_aux,rough,tven2,vikin,avvfr)   
        else  
           avvfr = velfr
           tven2 = tveno
        end if
        !
        ! Compute prescribed traction
        !
        if( itask == 1 ) then 
           if (kfl_mlwm_ker == 0_ip) then
              yplus = delta_aux*velfr/vikin
              if( yplus < 5.0_rp .and. kfl_ustar == 0 ) then
                 !
                 ! t = - mu/y * u
                 !
                 fact1 = gbden * vikin / delta_aux       
         
              else if( kfl_ustar == 1 ) then
                 velfr = 0.41_rp / log( (delta_aux+rough) / rough )
                 fact1 = gbden * tven2 * velfr * velfr    
              else if( kfl_ustar == 2 ) then
                 velfr = 0.41_rp / log( (delta_aux+rough) / rough )
                 fact1 = gbden * velfr * ( kinen**0.5_rp *cmu_st**0.25_rp)           
              else
                 !
                 ! t = - rho*U*^2/|u_tan-u_fix_tan| * (u_tan-u_fix_tan)
                 !
                 fact1 = gbden * avvfr * velfr / tveno   ! avvfr = velfr if time-averaging is not used
              end if
           else
              fact1 = gbden * avvfr * velfr / tveno   ! avvfr = velfr if time-averaging is not used
           end if
           if (kfl_aveme_ker==2.or.kfl_aveme_ker==3) then
              do idime = 1,ndime
                 tract(idime) = tract(idime) - fact1 * tveav(idime) 
              end do
           else
              do idime = 1,ndime
                 tract(idime) = tract(idime) - fact1 * tvedi(idime) 
              end do
           end if

           if( run_config%customer == cust_cfwind0 ) then
              !
              ! t = - rho*U*^2/|u| *u 
              !
              call vecnor(ub,ndime,tveno,2_ip)
              fact1 = gbden * usref * usref / tveno            
              do idime = 1,ndime
                 tract(idime) = tract(idime) + fact1 * ub(idime) 
              end do
           end if

        end if

     else

        !-----------------------------------------------------------------
        !
        ! Implicit treatment
        !
        !-----------------------------------------------------------------

        velsh = 0.0_rp
        do idime = 1,ndime                                     ! Tangent veloc.
           tvelo(idime) = vewal(idime)
           tvefi(idime) = velfi(idime)
           do inodb = 1,pnodb
              idofn = (inodb-1)*ndime+idime
              velsh(idime,idofn) = velsh(idime,idofn) &
                   + gbsha(inodb)
           end do
           do ldime = 1,ndime
              tvelo(idime) = tvelo(idime) &
                       - baloc(idime,ndime)   &
                       * baloc(ldime,ndime) * vewal(ldime)
              tvefi(idime) = tvefi(idime)   &
                       - baloc(idime,ndime) &
                       * baloc(ldime,ndime) * velfi(ldime)
              do inodb = 1,pnodb
                 idofn = (inodb-1)*ndime+ldime
                 velsh(idime,idofn) = velsh(idime,idofn) &
                      - baloc(idime,ndime)               &
                      * baloc(ldime,ndime)*gbsha(inodb)
              end do
           end do
           tvedi(idime) = tvelo(idime) - tvefi(idime)
        end do        

        !
        ! Compute U*
        !
        vikin = gbvis/gbden                                    ! nu
        call vecnor(tvedi,ndime,tveno,2_ip)                    ! |u_tan-u_fix_tan|
        
		! Coding here for machine learning interfacing
        if (kfl_mlwm_ker == 1) then
!            if (tveno .gt. 1e-10) tveno = 1e-6
!			velfr = velfr
			velfr = velfr/tveno
		else  
        	call frivel(kfl_ustar,delta_aux,rough,tveno,vikin,velfr)         ! U*
		end if
        if (kfl_wlaav_ker == 1_ip) then
           call vecnor(tveav,ndime,tven2,2_ip)               
           call frivel(kfl_ustar,delta_aux,rough,tven2,vikin,avvfr)  
        else
           avvfr = velfr
           tven2 = tveno
        end if

        if( itask == 1 ) then
           !
           ! Compute prescribed traction
           !
           yplus = delta_aux * velfr / vikin
           if( yplus < 5.0_rp  .and. kfl_ustar == 0 ) then
              fact1 = gbden * vikin / delta_aux 
           else if( kfl_ustar == 1 ) then 
              velfr = 0.41_rp  / log((delta_aux+rough) / rough)
              fact1 = gbden * tven2 * velfr * velfr    
           else if( kfl_ustar == 2 ) then ! ABL 2
              velfr = 0.41_rp / log((delta_aux+rough) / rough)
              
              if (kfl_regularization) then
                 reguk = regul_k(kinen)
              else
                 reguk = kinen
              end if
              velf2 = reguk**0.5_rp *cmu_st**0.25_rp ! frveloc in terms of kinetic turbulent energy 
!              if (velf2.gt.1.0e-1) then
                 fact1 = gbden * velfr * velf2
!              else
!                 fact1 = gbden * tveno * velfr * velfr
!              end if
           else
              fact1 = gbden * avvfr * velfr / tveno   ! avvfr = velfr if time-averaging is not used     
           end if

           do idime = 1,ndime
              tract(idime) = tract(idime) + fact1 * tvefi(idime) 
           end do
           if (kfl_waexl_ker==0) then
              do inodb = 1,pnodb
                 fact2 = fact1 * gbsha(inodb)
                 ievab = (lboel(inodb)-1)*ndofn
                 do idime = 1,ndime
                    ievab = ievab+1
                    do jnodb = 1,pnodb
                       jevab=(lboel(jnodb)-1)*ndofn
                       do jdime = 1,ndime                          
                          jevab = jevab+1
                          wmatr(ievab,jevab) = wmatr(ievab,jevab) &
                               + fact2 * velsh(idime,(jnodb-1)*ndime+jdime)
                       end do
                    end do
                 end do
              end do
           else  ! exchange location implicit
              ielem = lelbo
              pelty = ltype(ielem)
              pnode = nnode(pelty)
              !  Utan,j = Uk,a *shape(a)*(delta (j,k)-nj*nk)
              !  Utan,j = Uk,a* vels2(a,j,k)        
             
              vels2 = 0.0_rp
              indei = lexlo_ker(igaub,iboun)
              do idime =1, ndime          
                 vels2(idime,idime,1:pnode) = shape_waexl(1:pnode, indei)
                 do jdime =1, ndime
                    vels2(idime,jdime,1:pnode) = vels2(idime, jdime, 1:pnode) - baloc(idime, ndime)*baloc(jdime, ndime)*shape_waexl(1:pnode, indei)
                 end do
              end do
              do jnode =1, pnode
                 jevab=(jnode-1)*ndofn                   
                 do jdime = 1,ndime                          
                    jevab = jevab+1
                    do inodb = 1,pnodb
                       fact2 = fact1 * gbsha(inodb)
                       ievab = (lboel(inodb)-1)*ndofn
                       do idime = 1,ndime
                          ievab = ievab+1
                          wmatr(ievab,jevab) = wmatr(ievab,jevab) &
                          !     + fact2 * velsh(idime,(jnodb-1)*ndime+jdime)
                               + fact2 * vels2(idime,jdime,jnode)                               
                       end do
                    end do
                 end do                
              end do          
           end if
           
        end if

     end if

  end if
  !
  ! Second order correction: OJO QUITAR
  !
  h1 = 0.213767_rp
  h2 = 0.255452_rp
  hs = h2*(h1+h2)
 ! u1 = veloc(1,81,1)
 ! u2 = veloc(1,80,1)
 ! u3 = veloc(1,78,1)
 ! tract(1) = tract(1) + gbvis/hs* ( h2*u1 - (h1+h2)*u2 + h1*u3 ) * baloc(ndime,ndime)

end subroutine nsi_bouwal
