!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine sld_contac(itask,ielem,iboun,igaub,delta,ctant,cnorm,tcbar,yybar)

  !-----------------------------------------------------------------------
  !****f* Solidz/sld_contac
  ! NAME
  !    sld_contac
  ! DESCRIPTION
  !    Compute the traction and tangent stiffness from contact penalty model
  !
  !    DELTA ... Crack opening (normal and tangential) ........... delta
  !    TCBAR ... Contact traction (normal and tangential) ........ tc_bar
  !    YYBAR ... Contact tangent stiffness ....................... dtc_bar/ddelta
  !    CNORM ... Unit vector normal to the crack plane ........... n
  !    CTANT ... Unit vector perpendicular to the crack plane .... s
  ! USES
  ! USED BY
  !    sld_cohpre
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp,lg
  use def_domain, only       :  ndime,nmate
  use def_solidz, only       :  dslip_sld,lmate_sld,parcf_sld,dfric_sld


  implicit none
  integer(ip), intent(in)    :: igaub,itask,ielem,iboun
  real(rp),    intent(in)    :: ctant(ndime),cnorm(ndime),delta(2)
  real(rp),    intent(out)   :: tcbar(2),yybar(2,2)
  integer(ip)                :: pmate,icomp,jcomp
  real(rp)                   :: tnorm,dnorm,dtang,dtdot,ttang,fslip
  real(rp)                   :: frico,nopen,tapen,ttold,ttria
  
  logical(lg)                ::debugging
  
  debugging = .false. !(ielem == 23) .and. (igaub == 2) 
 
  pmate = 1
  if( nmate > 1 ) then
    pmate = lmate_sld(ielem)
  end if
  nopen = parcf_sld(1,pmate)
  tapen = parcf_sld(2,pmate)
  frico = parcf_sld(3,pmate)

  dnorm = delta(1)
  dtang = delta(2)
  dtdot = dtang - dslip_sld(ielem,iboun,igaub,2)
  ttold = dfric_sld(ielem,iboun,igaub,2)
  if (debugging) write(*,*)' '
  if (debugging) print*,'compression'
  if (debugging) write(*,'(1x,a8,1x,e12.5)')'dtang = ',dtang
  if (debugging) write(*,'(1x,a8,1x,e12.5)')'dtdot = ',dtdot

  !
  !  Compute normal traction
  !
  tnorm = nopen * dnorm

  !
  !  Compute tangential traction
  ! 
  ttria = ttold + tapen * dtdot
  fslip = ttria + tnorm * frico
  if (debugging) write(*,'(1x,a8,1x,e12.5)')'ttria = ',ttria
  if (debugging) write(*,'(1x,a8,1x,e12.5)')'fslip = ',fslip
  if (fslip <= 0.0_rp) then
     if (debugging) print*,'stick'
     ttang = ttria
  else
     if (debugging) print*,'slide'
     ttang = -1.0_rp * tnorm * frico
  end if
  dfric_sld(ielem,iboun,igaub,1) = ttang
  
  !
  ! Traction vector in local coordinate system
  !
  tcbar(1) = tnorm
  tcbar(2) = ttang

  !
  ! Tangent stiffness in local coordinate system (if implicit)
  !
  if (itask == 2) then
     yybar(1,1) = nopen
     yybar(1,2) = 0.0_rp     
     if (fslip <= 0.0_rp) then
        yybar(2,1) = 0.0_rp
        yybar(2,2) = tapen
     else
        yybar(2,1) = -1.0_rp * nopen * frico
        yybar(2,2) = 0.0_rp
     end if
  end if

  if (debugging) then
     write(*,*)'ielem = ',ielem,' ** iboun = ',iboun,' ** igaub = ',igaub
     write(*,'(1x,a8,1x,e12.5)')'dtdot = ',dtdot
     write(*,'(1x,a8,1x,e12.5)')'ttold = ',ttold
     write(*,*)'tcbar =     ** yybar =                  ** delta = '
     do icomp = 1,2
        write(*,'(1x,e12.5,2x,2(1x,e12.5),2x,1x,e12.5)')tcbar(icomp),(yybar(icomp,jcomp),jcomp =1,2),delta(icomp)
     end do
  end if
  
100 format (5(E16.8,','))
end subroutine sld_contac

