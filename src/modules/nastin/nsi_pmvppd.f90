!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_pmvppd(CLO,MET,WME,TA,TR,VEL,RH,PMV,PPD,ISTAT)
  !-----------------------------------------------------------------------
  !****f* Nastin/nsi_pmvppd
  ! NAME 
  !    nsi_pmvppd
  ! DESCRIPTION
  !    
  !    This routine computes:
  !    Predicted mean vote (PMV) and Predicted percentage dissatisfied (PPD)
  !    according to ISO 7730.
  !    according to ISO 7730.
  !
  !    Prof Ole Fanger proposed that the condition for thermal comfort 
  !    is that the skin temperature and sweat secretion lies within 
  !    narrow limits. Fanger obtained data from climate chamber 
  !    experiments, in which sweat rate and skin temperature were measured 
  !    on people who considered themselves comfortable at various metabolic 
  !    rates. Fanger proposed that optimal conditions for thermal comfort were 
  !    expressed by the regression line of skin temperature and sweat rate on 
  !    metabolic rate in data from these experiments. In this way an 
  !    expression for optimal thermal comfort can be deduced from the metabolic 
  !    rate, clothing insulation and environmental conditions.
  !
  !    The final equation for optimal thermal comfort is fairly complex and 
  !    need not concern us here. Fanger has solved the equations by 
  !    computer and presented the results in the form of diagrams from which 
  !    optimal comfort conditions can be read given a knowledge of metabolic 
  !    rate and clothing insulation.
  !  
  !    If CLO, MET and WME are positive, the value is in Clo and Met units.
  !    Otherwise, they are in K/W and W/m^2.
  !    IF RH is positive, the value is in %. Otherwise, RH is the pressure
  !    given in Pa.
  !
  !    CLO ... Clothing factor       [0,2]     Clo   [0,0.31]  K/W
  !    MET ... Metabolism            [0.8,4.0] Met   [46,232] W/m^2 
  !    WME ... Wet metabolism        [0.8,4.0] Met   [46,232] W/m^2 
  !    TA .... Ambient temperature   [10,30]   C
  !    TR .... Radiative temperature [10,40]   C
  !    VEL ... Air velocity          [0,1]     m/s
  !    RH .... Relative Humidity     [30,70]   %     [0,2700] Pa  
  !
  ! USES
  ! USED BY
  !    Nastin
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only    :  ip,rp
  implicit none
  real(rp),    intent(in)  :: CLO,MET,WME,TA,TR,VEL,RH
  real(rp),    intent(out) :: PMV,PPD
  integer(ip), intent(out) :: ISTAT
  integer(ip)              :: NMAX,N
  real(rp)                 :: PA,P2,P3,ICL,M,W,MW,HCF
  real(rp)                 :: FCL,XN,HC,XD,TCL,HL1,HL2
  real(rp)                 :: HL3,HL4,HL5,HL6,TS,EPS,TSK,FCLTA,ICLFCL

  ISTAT = 0
  NMAX  = 150
  EPS   = 0.000015_rp
  !
  ! Conversions
  !
  if(RH<0.0_rp) then
     PA = -RH
  else
     PA =  RH*10.0_rp*exp(16.6536_rp-4030.183_rp/(TA+235.0_rp)) ! Pa
  end if
  if(MET<0.0_rp) then
     M = -MET
  else
     M =  MET*58.15_rp
  end if
  if(WME<0.0_rp) then
     W = -WME
  else
     W =  WME*58.15_rp
  end if
  iF(CLO<0.0_rp) then
     ICL = -CLO
  else
     ICL = 0.155_rp*CLO
  end if

  MW  = M-W
  if (ICL < 0.078_rp ) then
     FCL = 1.0_rp  + 1.29_rp*ICL
  else
     FCL = 1.05_rp + 0.645_rp*ICL
  end if
  HCF = 12.1_rp*sqrt(VEL)                             ! Forced convection
  TSK = 35.7_rp-0.0275_rp*MW                          ! Skin temprature
  P2  = 3.96_rp*1.0d-8*FCL   
  P3  = (TR+273.0_rp)**4
  !
  ! TCL: Surface temperature of clothing
  !
  N      = 0
  XD     = 1.0_rp
  TCL    = TA+(35.5_rp-TA)/(3.5_rp*(6.45_rp*ICL+0.1_rp))    ! Initial guess
  FCLTA  = FCL*TA
  ICLFCL = FCL*ICL
  do while( N < NMAX .AND. XD>EPS )
     N   = N+1
     HC  = max(2.38_rp*(TCL-TA)**0.25_rp,HCF)         ! Heat transfer coef: max(natural,forced convection)
     XN  = (TSK-ICL*(P2*((TCL+273.0_rp)**4-P3)-FCLTA*HC))&
          /(1.0_rp+ICLFCL*HC)
     XD  = abs(XN-TCL)
     if(XN/=0.0_rp) XD = XD/abs(XN)
     TCL = XN
  end do
  !
  ! PMV
  !
  if( N == NMAX ) then
     ISTAT = 1
  else 
     HL1 = 3.05_rp*1.0d-3*(5733.0_rp-6.99_rp*MW-PA)
     HL2 = 0.42_rp*(MW-58.15_rp)
     HL3 = 1.7_rp*1.0d-5*M*(5867.0_rp-PA)
     HL4 = 0.0014_rp*M*(34.0_rp-TA)
     HL5 = P2*((TCL+273.0_rp)**4-P3)
     HL6 = FCL*HC*(TCL-TA)
     TS  = 0.303_rp*exp(-0.036_rp*M)+0.028_rp         ! Thermal sensation transfer coef. 
     PMV = TS*(MW-HL1-HL2-HL3-HL4-HL5-HL6)
  end if
  !
  ! PPD
  !
  if( ISTAT == 0 ) then
     PPD = 100.0_rp - 95.0_rp*exp(-(0.03353_rp*PMV**4+0.2179_rp*PMV**2))
  else
     PPD = 999.0_rp
  end if

end subroutine nsi_pmvppd

