!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tur_deslsc(&
     ndime,gptur,gpden,gpvis,gpmut,&
     gpwal,gpgrv,hleng,fddes,gddes,betas)
  !------------------------------------------------------------------------
  !****f* Turbul/tur_sstkom
  ! NAME
  !   tur_deslsc
  ! DESCRIPTION
  !    Computes the (I)DDES length scale
  ! OUTPUT 
  !    FDDES
  !    GDDES    DES zone
  ! USES
  ! USED BY
  !    tur_sstkom
  !***
  !------------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp
  use def_turbul, only       :  nturb_tur, kfl_ddesm_tur, cddes_tur
  implicit none
  integer(ip), intent(in)    :: ndime
  real(rp),    intent(in)    :: gptur(nturb_tur)
  real(rp),    intent(out)   :: fddes, gddes
  real(rp),    intent(in)    :: gpden,gpvis,gpmut
  real(rp),    intent(in)    :: gpwal, betas
  real(rp),    intent(in)    :: gpgrv(ndime,ndime),hleng(ndime)
  integer(ip)                :: idime
  real(rp)                   :: gpkin,gpome, kap

  real(rp)                   :: a,rd,fd,sum_gpgrv,L_rans,L_les, L_ddes
  real(rp)                   :: ft, fl, fe2, fe1, fe, fb, fdt, fdtil, rdt, rdl, cl, ct, alpha, cw
  real(rp)                   :: nu
  integer(ip)                :: jdime


  ! turbulence variables
  gpkin = max(0.0_rp,gptur(1))
  gpome = gptur(2)
  kap   = 0.41_rp

  ! output values
  fddes = 0.0_rp
  gddes = 0.0_rp

  a = 0.0_rp
  L_rans = 0.0_rp
  L_les = 0.0_rp
  L_ddes = 0.0_rp

  ft = 0.0_rp
  fl = 0.0_rp
  fe2 = 0.0_rp
  fe1 = 0.0_rp
  fe = 0.0_rp
  fb = 0.0_rp
  fdt = 0.0_rp
  fdtil = 0.0_rp
  rdt = 0.0_rp
  rdl = 0.0_rp
  cl = 5.0_rp
  ct = 1.87_rp
  rd= 0.0_rp
  fd= 0.0_rp
  alpha = 0.0_rp
  cw = 0.15_rp

  ! warning nut = mut = gpmut(IN for this subroutine)
  nu    = gpvis/gpden
  !
  ! Wall distance GPWAL(dw==gpwal (in of this subroutine))
  !
  ! Velocity gradient GPGRV(j,i)(in of this subroutine)) = grad(a) = da(i)/dx(j)
  !
  !---sqrt(UijUij)
  !
  sum_gpgrv = 0.0_rp
  do idime = 1,ndime
     do jdime = 1,ndime
        sum_gpgrv = sum_gpgrv + gpgrv(jdime,idime)*gpgrv(jdime,idime)
     end do
  end do
  sum_gpgrv=sqrt(sum_gpgrv)
  sum_gpgrv = max(sum_gpgrv,10.0_rp**(-10_ip))
  !------------------------------------------------


  ! calculate blending functions
  alpha = 0.25_rp - gpwal / hleng(1)

  a = kap*kap*gpwal*gpwal
  rdl = nu/(a*sum_gpgrv)
  rdt = gpmut/(a*sum_gpgrv)
  rd = rdl + rdt

  ft = tanh((ct*ct*rdt)**3_ip)
  fl = tanh((cl*cl*rdl)**10_ip)

  fd = 1.0_rp - tanh((8.0_rp*rd)**3_ip)
  fdt = 1.0_rp - tanh((8.0_rp*rdt)**3_ip)

  fe2 = 1.0_rp - max(ft,fl)
  if(alpha >= 0.0_rp) then
    fe1 = 2.0_rp * exp(-11.09_rp*alpha*alpha)
  else
    fe1 = 2.0_rp * exp(-9.0_rp*alpha*alpha)
  endif
  fe = max((fe1 - 1.0_rp),0.0_rp) * fe2
  fb = min(2.0_rp * exp(-9.0_rp*alpha*alpha),1.0_rp)
  fdtil = max(1.0_rp - fdt, fb)

  ! calculate basic length scales

  L_rans=sqrt (gpkin) / (betas * gpome)

  if (kfl_ddesm_tur == 1) then

    if (ndime==2) then
        L_les= ( hleng(1) * hleng(2) )** 0.5_rp
    else if(ndime==3) then
        L_les= ( hleng(1) * hleng(2) * hleng(3) )** (1.0_rp / 3.0_rp)
    end if

  else if (kfl_ddesm_tur == 2) then

    if (ndime==2) then
        L_les= min(max(cw*gpwal,cw*hleng(1),hleng(2)),hleng(1))
    else if(ndime==3) then
        L_les= min(max(cw*gpwal,cw*hleng(1),hleng(3)),hleng(1))
    end if

  else
    call runend("not supported")
  endif
  L_les = L_les * cddes_tur

  ! calculate (I)DDES length scale

  !------------------------------------------------

  if (kfl_ddesm_tur == 1) then
    L_ddes = L_rans - fd * max(0.0_rp,L_rans - L_les)
  else if (kfl_ddesm_tur == 2) then
    L_ddes = fdtil * (1.0_rp + fe) * L_rans + (1.0_rp - fdtil) * L_les
  else
    call runend("not supported")
  endif

  fddes= L_rans / max(L_ddes,1.0e-10_rp)

  ! for postprocessing only
  if (kfl_ddesm_tur == 1) then
    if(L_rans > L_les) then
        gddes= fd
    else
        gddes = 0.0_rp
    end if
  else if (kfl_ddesm_tur == 2) then
        gddes = (1_rp-fdtil)  !TODO better expresion!
  else
    call runend("not supported")
  endif


end subroutine tur_deslsc
