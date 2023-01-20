!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tem_exasol(itask,kfl_exacs_tem,ndime,&
     gpcod,exadv,exdif,exrea,extem,exgrt) 

!-----------------------------------------------------------------------
!    
! Exact solution for the temperature equation
!
!-----------------------------------------------------------------------
  use def_kintyp
  use def_parame
  implicit none
  integer(ip), intent(in)  :: itask,kfl_exacs_tem,ndime
  real(rp),    intent(in)  :: gpcod(ndime)
  real(rp),    intent(out) :: exadv,exdif,exrea,extem,exgrt(ndime)
  real(rp)                 :: x,y
  
  x=gpcod(1)
  y=gpcod(2)

    
  select case(kfl_exacs_tem)

  case(1)
     if(itask==0) then
        exadv    = 0.0_rp
        exdif    =   exp(x*y)*sin(pi*y)*(2.0_rp*y*pi*cos(pi*x)+(y*y-pi*pi)*sin(pi*x))&
             &     + exp(x*y)*sin(pi*x)*(2.0_rp*x*pi*cos(pi*y)+(x*x-pi*pi)*sin(pi*y))
        exrea    = 0.0_rp
     else if(itask==1) then
        extem    = sin(pi*x)*sin(pi*y)*exp(x*y)
        exgrt(1) = sin(pi*y)*exp(x*y)*(pi*cos(pi*x)+y*sin(pi*x))
        exgrt(2) = sin(pi*x)*exp(x*y)*(pi*cos(pi*y)+x*sin(pi*y))
     end if

  case(2)
     if(itask==0) then
        exadv    = 0.0_rp
        exdif    = 0.0_rp
        exrea    = 0.0_rp
     else if(itask==1) then
        extem    = x+y
        exgrt(1) = 1.0_rp
        exgrt(2) = 1.0_rp
     end if

  case(3)
     if(itask==0) then
        exadv    = -y+0.60_rp+5.0_rp*(x-0.60_rp)
        exdif    = 0.0_rp
        exrea    = x+5.0_rp*y
     else if(itask==1) then
        extem    = x+5.0_rp*y
        exgrt(1) = 1.0_rp
        exgrt(2) = 5.0_rp
     end if

  end select
  
end subroutine tem_exasol
  
