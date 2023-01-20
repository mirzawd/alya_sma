!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine coleve(&
     wtask,igaui,igauf,pgaus,gpfle,densi,densa,&
     visco,visca,thicl,hleng,gpden,gpvis)
  !-----------------------------------------------------------------------
  !****f* mathru/coleve
  ! NAME 
  !    coleve
  ! DESCRIPTION
  !
  ! USES
  !
  ! USED BY
  !    modules
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp  
  use def_parame, only       :  pi
  implicit none
  character(2), intent(in)   :: wtask
  integer(ip),  intent(in)   :: igaui,igauf,pgaus
  real(rp),     intent(in)   :: gpfle(pgaus)
  real(rp),     intent(in)   :: densi
  real(rp),     intent(in)   :: densa
  real(rp),     intent(in)   :: visco
  real(rp),     intent(in)   :: visca
  real(rp),     intent(in)   :: thicl
  real(rp),     intent(in)   :: hleng
  real(rp),     intent(out)  :: gpden(pgaus)
  real(rp),     intent(out)  :: gpvis(pgaus)
  integer(ip)                :: igaus
  real(rp)                   :: phi,eps,f

  if( thicl > 0 ) then
     eps = thicl
  else
     eps = -thicl*hleng
  end if

  if( wtask == '11' ) then
     !
     ! Density and viscosity
     !
     do igaus = igaui,igauf
        phi = gpfle(igaus)
        if( phi < -eps ) then
           gpden(igaus) = densa
           gpvis(igaus) = visca
        else if( phi > eps ) then
           gpden(igaus) = densi
           gpvis(igaus) = visco
        else 
           f = 0.5_rp*(1.0_rp+phi/eps+sin(pi*phi/eps)/pi)
           !f = 0.5_rp*phi/eps + 0.5_rp
           gpden(igaus) = densa + (densi-densa)*f
           gpvis(igaus) = visca + (visco-visca)*f
        end if
     end do

  else

     if( wtask(1:1) == '1' ) then
        !
        ! Density
        !
        do igaus = igaui,igauf
           phi = gpfle(igaus)
           if( phi < -eps ) then
              gpden(igaus) = densa
           else if( phi > eps ) then
              gpden(igaus) = densi
           else 
              f = 0.5_rp*(1.0_rp+phi/eps+sin(pi*phi/eps)/pi)
              !f = 0.5_rp*phi/eps + 0.5_rp
              gpden(igaus) = densa + (densi-densa)*f
           end if
        end do
     end if

     if( wtask(2:2) == '1' ) then
        !
        ! Viscosity
        !
        do igaus = igaui,igauf
           phi = gpfle(igaus)
           if( phi < -eps ) then
              gpvis(igaus) = visca
           else if( phi > eps ) then
              gpvis(igaus) = visco
           else 
              f = 0.5_rp*(1.0_rp+phi/eps+sin(pi*phi/eps)/pi)
              !f = 0.5_rp*phi/eps + 0.5_rp
              gpvis(igaus) = visca + (visco-visca)*f
           end if
        end do
     end if
     
  end if

end subroutine coleve
