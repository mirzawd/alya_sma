!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine mag_inimat()

  use def_domain, only: meshe, ndime
  use mod_communications, only: PAR_MAX  
  use def_kermod, only: ndivi

  use def_magnet

  implicit none

  integer(ip) :: &
    ielem

  !##############################################################
  maxmat_mag = 0_ip
  do ielem = 1, meshe(ndivi) % nelem
    maxmat_mag = max(meshe(ndivi) % lmate(ielem), maxmat_mag)
  end do
  call PAR_MAX(maxmat_mag)
  !##############################################################

  !##############################################################
  allocate(constrlist_mag(maxmat_mag))
  allocate(selfList_mag(maxmat_mag))
  !##############################################################

  !##############################################################
  allocate(Ec0_mag(maxmat_mag), nc0_mag(maxmat_mag), Jc0_mag(maxmat_mag), mur_mag(maxmat_mag), rho_mag(maxmat_mag))
  allocate(B0k_mag(maxmat_mag))
  allocate(resistOpt_mag(maxmat_mag), scalinOpt_mag(maxmat_mag))

  resistOpt_mag = -1_ip
  scalinOpt_mag = -1_ip
  !##############################################################

  !##############################################################
  allocate(kfl_resiso_mag(maxmat_mag), kfl_Jcriso_mag(maxmat_mag))
  allocate(kfl_Ecriso_mag(maxmat_mag), kfl_ncriso_mag(maxmat_mag))
  allocate(kfl_Jc0iso_mag(maxmat_mag), kfl_Bc0iso_mag(maxmat_mag), kfl_Tc0iso_mag(maxmat_mag))
  allocate(kfl_muriso_mag(maxmat_mag), kfl_rhoiso_mag(maxmat_mag))

  kfl_resiso_mag = .false.
  kfl_Jcriso_mag = .false.
  kfl_Ecriso_mag = .false.
  kfl_ncriso_mag = .false.
  kfl_Jc0iso_mag = .false.
  kfl_Bc0iso_mag = .false. 
  kfl_Tc0iso_mag = .false.
  kfl_muriso_mag = .false.
  kfl_rhoiso_mag = .false.

  allocate(resmat_mag(ncomp_mag,maxmat_mag), Jcrmat_mag(ncomp_mag,maxmat_mag))
  allocate(Ecrmat_mag(ncomp_mag,maxmat_mag), ncrmat_mag(ncomp_mag,maxmat_mag), Jc0mat_mag(ncomp_mag,maxmat_mag), murmat_mag(ncomp_mag,maxmat_mag), &
           rhomat_mag(ncomp_mag,maxmat_mag), Bc0mat_mag(ncomp_mag,maxmat_mag), Tc0mat_mag(ncomp_mag,maxmat_mag))

  resmat_mag = -1_ip
  Jcrmat_mag = -1_ip
  Ecrmat_mag = 0.0_rp
  ncrmat_mag = 0.0_rp
  Jc0mat_mag = 0.0_rp
  murmat_mag = 0.0_rp
  rhomat_mag = 0.0_rp
  Bc0mat_mag = 0.0_rp
  Tc0mat_mag = 0.0_rp
  !##############################################################

  !##############################################################
  allocate(magnen_mag(maxmat_mag), joulen_mag(maxmat_mag))

  magnen_mag = 0.0_rp
  joulen_mag = 0.0_rp
  !##############################################################

  !##############################################################
  allocate(magtiz_mag(3_ip, maxmat_mag))

  magtiz_mag = 0.0_rp

  allocate(cursum_mag(3_ip, maxmat_mag))

  cursum_mag = 0.0_rp

  allocate(magsum_mag(3_ip, maxmat_mag))

  magsum_mag = 0.0_rp

  allocate(volume_mag(maxmat_mag))

  volume_mag = 0.0_rp
  !##############################################################

  !##############################################################
  allocate(momori_mag(ndime, maxmat_mag))

  momori_mag = 0.0_rp
  !##############################################################

end subroutine mag_inimat
