!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



function funcre(param,npara,ifuge,timev)

  !------------------------------------------------------------------------
  !
  ! This function yields a parabolic or periodic evolution
  !
  !------------------------------------------------------------------------

  use def_kintyp, only     :  ip,rp
  use def_master, only     :  funin,funou
  use def_domain, only     :  ndime
  use def_parame, only     :  pi
  implicit none
  real(rp)                 :: funcre,funcre_prim,funcre_mod
  integer(ip), intent(in)  :: npara,ifuge
  real(rp),    intent(in)  :: param(npara),timev
  integer(ip)              :: ipara,idime
  real(rp)                 :: timea,timeb,funca,funcb,zerom,timec
  real(rp)                 :: timei,timef
  real(rp)                 :: coor0(3),betaa,radiu,alpha,deltx,deltz
  real(rp)                 :: slmix1,slmix2,whmix1,whmix2
  real(rp)                 :: f,g,H1,T1,H1_prim,T1_prim,timev_prim,H2
  real(rp)                 :: a1,b1,c1,a2,b2,c2,amp

  funcre = 0.0_rp
  funcre_prim = 0.0_rp
  zerom  = epsilon(1.0_rp)
  if( ifuge == 0 ) then
     !
     ! No time dependence 
     !
     funcre = 1.0_rp

  else if( ifuge == 1 ) then
     !
     ! Parabolic evolution
     !
     if( param(1)-zerom <= timev .and. timev <= param(2)+zerom) then 
        funcre = param(3)*timev*timev+param(4)*timev+param(5)
     else if ( timev > param(2)+zerom ) then
        timea  = param(2)
        funcre = param(3)*timea*timea+param(4)*timea+param(5)
     else if ( timev < param(1)-zerom ) then
        timea  = param(1)
        funcre = param(3)*timea*timea+param(4)*timea+param(5)
     end if

  else if( ifuge == 2 ) then
     !
     ! Periodic evolution
     !
     if( param(1)-zerom <= timev .and. timev <= param(2)+zerom) then 
        funcre = param(3)*cos(param(4)*timev+param(5))+param(6)
     else if ( timev > param(2)+zerom ) then
        timea  = param(2)
        funcre = param(3)*cos(param(4)*timea+param(5))+param(6)
     else if ( timev < param(1)-zerom ) then
        timea  = param(1)
        funcre = param(3)*cos(param(4)*timea+param(5))+param(6)
     end if

  else if( ifuge == 3 .or. ifuge == -3 ) then
     !
     ! Discrete evolution
     !
     timei = param(1)
     timef = param((npara/2-1)*2+1)

     if( timev <= timei ) then
        funcre = param(2)
     else
        if( timev >= timef ) then            ! Look for the time inside the period
           timec = timev
           do while( timec > timef )
              timec = timec - (timef-timei)
           end do
        else
           timec = timev
        end if
        ipara = 0
        do while( ipara < npara/2 )
           ipara = ipara+1
           
           if(timec <param((ipara-1)*2+1)) then
              timea  = param((ipara-1)*2+1)
              funca  = param((ipara-1)*2+2)
              timeb  = param((ipara-2)*2+1)
              funcb  = param((ipara-2)*2+2)
              funcre = (funcb-funca)/(timeb-timea)*(timec-timea)+funca
              ipara  = npara/2
           end if
        end do
     end if

  else if( ifuge == 4 ) then
     !
     ! Special function to change boundary values
     !
     funcre = param(2)

  else if( ifuge == 5 ) then
     !
     ! Marek Prymon's function
     !
     if( timev < param(1) ) then
        funcre = param(2)
     else if( timev < param(3) ) then
        funcre = param(4) - param(5) * timev
     else
        funcre = param(6)
     end if

  else if( ifuge == 6 ) then
     !
     ! Translation
     !
     funcre = 1.0_rp
     do idime = 1,ndime
        funou(idime) = param(idime) 
     end do

  else if( ifuge == 7 ) then
     !
     ! Rotation
     !
     funcre   = 1.0_rp
     coor0(1) = param(1)
     coor0(2) = param(2)
     coor0(3) = param(3)
     betaa    = param(4)*pi/180.0_rp
     deltx    = funin(1)-coor0(1)
     deltz    = funin(3)-coor0(3)
     radiu    = sqrt(  (deltx*deltx) + (deltz*deltz) )     
     alpha    = atan2( deltz , deltx )
     
     funou(1) =  ( radiu * cos(alpha+betaa) ) + coor0(1)
     funou(2) =  funin(2)
     funou(3) =  ( radiu * sin(alpha+betaa) ) + coor0(3) 
     do idime = 1,ndime
        funou(idime) = funou(idime) - funin(idime)
     end do
     
  else if( ifuge == 8 ) then
     !
     ! CYCLE function for cyclic breathing (sinusoidal time variation)
     ! Q=Q_peak*sin(2*pi*f*t)
     ! 
     !
     funcre=param(1)*sin(2.0_rp*pi*param(2)*timev)
     !
     ! CYCLE function for the nose 'CYCLE'
     !  
     ! periodic evolution with 6 ordre
     !
     !funcre = -0.000000582_rp&
     !     +0.000058387_rp*cos(2.579_rp*timev)+0.000598270_rp*sin(2.579_rp*timev)&
     !     -0.000004015_rp*cos(5.158_rp*timev)-0.000006084_rp*sin(5.158_rp*timev)&
     !     -0.000043093_rp*cos(7.737_rp*timev)+0.000111284_rp*sin(7.737_rp*timev)&
     !     -0.00000087_rp*cos(10.316_rp*timev)-0.000006765_rp*sin(10.316_rp*timev)&
     !     -0.000014982_rp*cos(12.895_rp*timev)+0.000028346_rp*sin(12.895_rp*timev)&  
     !     +0.000000056_rp*cos(15.474_rp*timev)-0.000005522_rp*sin(15.474_rp*timev)&
     !     -0.000002663_rp*cos(18.053_rp*timev)+0.000003572_rp*sin(18.053_rp*timev) 

  else if( ifuge == 9 ) then
     !
     ! SNIFF function for the nose 'SNIF2'
     !  
     ! Parabolic evolution 10th order
     !
     if( param(1)-zerom <= timev .and. timev <= param(2)+zerom) then
        funcre = param(3)*timev**10_rp + param(4)*timev**9_rp + param(5)*timev**8_rp+&
             param(6)*timev**7_rp + param(7)*timev**6_rp + param(8)*timev**5_rp +&
             param(9)*timev**4_rp + param(10)*timev**3_rp + param(11)*timev**2_rp +&
             param(12)*timev + param(13)
        
     else if ( timev > param(2)+zerom ) then
        timea  = param(2)
        funcre = param(3)*timev**10_rp + param(4)*timev**9_rp + param(5)*timev**8_rp+&
             param(6)*timev**7_rp + param(7)*timev**6_rp + param(8)*timev**5_rp +&
             param(9)*timev**4_rp + param(10)*timev**3_rp + param(11)*timev**2_rp +&
             param(12)*timev + param(13)
        
     else if ( timev < param(1)-zerom ) then
        timea  = param(1)
        funcre = param(3)*timev**10_rp + param(4)*timev**9_rp + param(5)*timev**8_rp+&
             param(6)*timev**7_rp + param(7)*timev**6_rp + param(8)*timev**5_rp +&
             param(9)*timev**4_rp + param(10)*timev**3_rp + param(11)*timev**2_rp +&
             param(12)*timev + param(13)
     end if
 
 else if( ifuge == 10 ) then
     !
     ! EAWAV evolution
     !
     ! function: (A/2)*(1+cos((2*pi/(t1-t0))*(t-tc)))
     ! Ewave: 1) t0e = 0 2) t1e = 450 3) Ae = 55 4) 2*pi/(t1e-t0e) 5) tce = 220
     ! Awave: 6) t0a = 500 7) t1a = 700 8) Aa = 39 9) 2*pi/(t1a-t0a) 10) tca = 600
     !
     if( param(1) <= timev .and. timev <= param(2)) then 
        funcre = (param(3)/2_rp)*(1_rp+cos(param(4)*(timev-param(5))))
     else if ( param(2) < timev .and. timev <= param(6) ) then
        funcre = 0_rp*timev
     else if ( param(6) < timev .and. timev <= param(7) ) then
        funcre = (param(8)/2_rp)*(1_rp+cos(param(9)*(timev-param(10))))
     end if

 else if( ifuge == 11 ) then
    !
    ! WINDKESSEL pressure function
    !
    ! P_T(t)=qt* (Rs+R*(1-e^(-t/(R*C))))
    !

    ! param(1) = flow rate

    funcre = param(1) * (param(2)+param(3)*(1.0_rp - exp(-timev/param(3)/param(4))))

 else if( ifuge == 12 ) then
    !
    ! KIAO function for cyclic breathing (sinusoidal time variation) KIAO INTHAVONG
    ! inspiration  lps=a1*np.sin(b1*time) a1=0.5229; b1=1.899; ojo (lps) then mass flow rate = lps /1000*density
    ! exhalation  mfr =(a1*np.sin(b1*time) + a2*np.sin(b2*time) + a3*np.sin(b3*time)) ojo (mfr)
    !
    ! parameter for mixing
    slmix1 = 10.0_rp
    whmix1 = 1.65_rp
    slmix2 = 10.0_rp
    whmix2 = 4.0_rp
    !
    funcre_mod = mod(timev,4.0_rp)
    !
    ! f(x) inhalation; g(x)exhalation
    !
    f      = (0.5229_rp * sin(1.899_rp *funcre_mod))*0.01164_rp
    g      = -(0.3674_rp * sin(1.429_rp *(funcre_mod-1.65_rp))+0.0351_rp& 
         *sin(6.223_rp*(funcre_mod-1.65_rp))+0.08071_rp*sin(3.572_rp*(funcre_mod-1.65_rp)))*0.01164_rp
    !
    ! mixing function for meeting point (1.65s and 4s) T1 for the first and T2 both
    !
    H1     = 0.5_rp*(tanh(slmix1*(funcre_mod-whmix1))+1.0_rp)
    T1     = g*H1+(1.0_rp-H1)*f
    !
    ! new timev and dependance
    !
    timev_prim  = funcre_mod - 4_rp
    H1_prim     = 0.5_rp*(tanh(slmix1*(timev_prim-whmix1))+1.0_rp)
    T1_prim     = g*H1_prim+(1.0_rp-H1)*f
    H2          = 0.5_rp*(tanh(slmix2*(funcre_mod-whmix2))+1.0_rp)
    !
    funcre = T1_prim*H2+(1.0_rp-H2)*T1*param(1) 
    !
 else if( ifuge == 13 ) then
    !
    ! MERCEDEZ-BENZ SIMULATION 
    !
    ! during tracheo 2 second jet of exhalation (or cough?)
    !
    amp=0.0003353_rp
    !
    if( 0.0_rp <= timev .and. timev <= 0.2_rp) then
       funcre = (amp/0.2_rp)*timev *param(1)
    else if ( 0.2_rp < timev .and. timev <= 1.8_rp ) then
       funcre = amp *param(1)
    else if ( 1.8_rp < timev .and. timev <= 2.0_rp ) then
       funcre = (-(amp/0.2_rp)*timev + (10.0_rp*amp))*param(1)
    end if
    !
 else if( ifuge == 14 ) then
    !    
    ! Cough time function (Gupta 2009)
    !
    !param(1) = CPFR(L/s) cough peak flow rate
    !param(2) = CEV(L)    cough expiration volume
    !param(3) = PVT(s)    peak velocity time
    !param(4) = factor corrector 
    !
    a1=1.68
    b1=3.338
    c1=0.428
    !
    a2=(param(2)/(param(3)*param(1)))-a1
    b2=((-2.158*param(2))/(param(3)*param(1)))+10.457
    c2=1.8/(b2-1)
    !
    ! new timev (dimensionless)
    !
    timev_prim  = timev/param(3)
    !
    !
    if ( 0.0_rp - zerom <= timev_prim .and. timev_prim <= 1.2_rp) then
       funcre_prim = (a1 * timev_prim ** (b1-1.0_rp) * exp(-timev_prim/c1))/(gamma(b1) * c1 ** b1)
       funcre= funcre_prim * param(1) * param(4)
    else if ( 1.2_rp - zerom <= timev_prim .and. timev_prim <= 9.375_rp) then
       funcre_prim = (a1 * timev_prim ** (b1-1.0_rp) * exp(-timev_prim/c1))/(gamma(b1) * c1 ** b1) + &
            (a2 * (timev_prim-1.2_rp) ** (b2-1.0_rp) * exp(-(timev_prim-1.2_rp)/c2))/(gamma(b2) * c2 ** b2)
       funcre=funcre_prim * param(1) * param(4)
    else if ( timev_prim >= 9.375_rp) then ! End of the cough (9.375=>0.9s) nulify the flow  
       funcre=0.0_rp 
    endif
 else if( ifuge == 15 ) then
    !                                 amp|  /\
    ! Sawtooth time function             | /  \
    !                                    |/____\___
    !param(1) = h half of sawtooth       0   h  e
    !param(2) = e end of sawtooth
    !param(3) = amp amplitude of sawtooth            Y1=(amp/h)*t
    !param(4) = alpha factor corrector               Y2=(amp/(h-e))*t+(e*amp)/(e-h)   
    !
    !
    if( 0.0_rp <= timev .and. timev < param(1)) then
       funcre = (param(3)/param(1))*timev *param(4)
    else if ( param(1) <= timev .and. timev <= param(2) ) then
       funcre = ((param(3)/(param(1)-param(2)))*timev +(param(2)*param(3)/(param(2)-param(1)))) * param(4)
    else if (timev > param(2) ) then
       funcre = 0.0_rp
    end if 

 end if
 
end function funcre

