!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine sld_cohesive_law_901(itask,ielem,iboun,igaub,delta,ctant,cnorm,tcbar,yybar)

  !-----------------------------------------------------------------------
  !****f* Solidz/sld_cohesive_law_901
  ! NAME
  !    sld_cohesive_law_901
  ! DESCRIPTION
  !    Compute the traction and tangent stiffness from Smith--Ferante cohesive law
  !            combined with contact penalty model
  !
  !    DELTA ... Crack opening (normal and tangential) ........... delta
  !    TCBAR ... Cohesive traction (normal and tangential) ....... tc_bar
  !    YYBAR ... Cohesive tangent stiffness ...................... dtc_bar/ddelta
  !    CNORM ... Unit vector normal to the crack plane ........... n
  !    CTANT ... Unit vector perpendicular to the crack plane .... s
  ! USES
  ! USED BY
  !    sld_cohpre
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp
  use def_domain, only       :  ndime,nmate
  use def_solidz, only       :  dcmax_sld,dceff_sld,lmate_sld,parch_sld


  implicit none
  integer(ip), intent(in)    :: igaub,itask,ielem,iboun
  real(rp),    intent(in)    :: ctant(ndime),cnorm(ndime),delta(2)
  real(rp),    intent(out)   :: tcbar(2),yybar(2,2)
  integer(ip)                :: pmate,icomp,jcomp
  real(rp)                   :: tceff,dceff,tcrit,dcinf,tcmax,dtcdd
  real(rp)                   :: dcdot,dcmax,lambd,dummr,dcrit
  
  logical debugging
  
  debugging = .false. ! (ielem == 81) .and. (igaub == 2) 
  
  pmate = 1
  if( nmate > 1 ) then
    pmate = lmate_sld(ielem)
  end if
  tcrit = parch_sld(1,pmate)
  dcrit = parch_sld(2,pmate)
  dcinf = parch_sld(3,pmate)
  lambd = parch_sld(4,pmate)
 
  !
  ! Effective scalar opening (mode I and mode II/III)
  !
  dceff = sqrt(delta(1)**2 + lambd**2*delta(2)**2)
  dceff_sld(ielem,iboun,igaub,1) = dceff
  dcmax = dcmax_sld(ielem,iboun,igaub)
  dcdot = dceff - dceff_sld(ielem,iboun,igaub,2)
  if (debugging) write(*,*)' '

  !
  ! Effective cohesive traction based on Smith--Ferrante law
  !
  if (dceff>=dcmax .and. dcdot>=0.0_rp .and. dceff<=dcinf) then
     !
     ! Case 1: Cohesive zone in loading
     !
     if (debugging) print*,'loading'
     dummr = dceff/dcrit
     tceff = exp(1.0_rp)*tcrit*dummr*exp(-1.0_rp*dummr)
     dtcdd = exp(1.0_rp)*tcrit/dcrit*exp(-1.0_rp*dummr)*(1.0_rp - dummr)
    
  else if ((dceff<dcmax .or. dcdot<0.0_rp) .and. dceff<=dcinf .and. dcmax<=dcinf) then
     !
     ! Case 2: Cohesive zone in unloading/reloading
     !
     if (debugging) then
        if (dcdot<0.0_rp) print*,'unloading'
        if (dcdot>0.0_rp) print*,'reloading'
     end if
     tcmax = tcrit * (1.0_rp - dcmax/dcinf )
     tceff = tcmax / dcmax * dceff
     dtcdd = tcmax / dcmax
      
  else
     !
     ! Case 3: Cohesive zone is fully-broken
     !
     if (debugging) print*,'fully-broken'
     tceff = 0.0_rp
     dtcdd = 0.0_rp
          
  end if
    
  !
  ! Traction vector in local coordinate system
  !
  tcbar(1) = tceff / dceff * delta(1)
  tcbar(2) = tceff / dceff * delta(2) * lambd**2
  !
  ! Tangent stiffness in local coordinate system (if implicit)
  !
  if (itask == 2) then
     yybar(1,1) = tceff/dceff + (dtcdd-tceff/dceff)*(delta(1)/dceff)**2
     yybar(1,2) = (dtcdd-tceff/dceff)*(delta(1)/dceff)*(delta(2)/dceff)*lambd**2
     yybar(2,1) = (dtcdd-tceff/dceff)*(delta(1)/dceff)*(delta(2)/dceff)*lambd**2
     yybar(2,2) = tceff/dceff*lambd**2 + (dtcdd-tceff/dceff)*((delta(2)/dceff)**2)*lambd**4
  end if

  if (debugging) then
     write(*,*)'ielem = ',ielem,' ** iboun = ',iboun,' ** igaub = ',igaub
     write(*,'(1x,a8,1x,e12.5)')'dceff = ',dceff
     write(*,'(1x,a8,1x,e12.5)')'dcmax = ',dcmax
     write(*,'(1x,a8,1x,e12.5)')'dcdot = ',dcdot
     write(*,'(1x,a8,1x,e12.5)')'tceff = ',tceff
     write(*,*)'tcbar =     ** yybar =                  ** delta = '
     do icomp = 1,2
        write(*,'(1x,e12.5,2x,2(1x,e12.5),2x,1x,e12.5)')tcbar(icomp),(yybar(icomp,jcomp),jcomp =1,2),delta(icomp)
     end do
  end if
  
100 format (5(E16.8,','))
end subroutine sld_cohesive_law_901

