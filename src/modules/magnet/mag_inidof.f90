!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine mag_inidof()

  use def_domain, only: lexis, iesta_dom, iesto_dom
  use def_magnet
  use mod_mag_nedele, only: mag_edgdof
  use mod_mag_lagran, only: mag_noddof

  implicit none

  integer(ip) :: &
    pelty

  maxdof_mag = 0_ip
  mxndof_mag = 0_ip
  do pelty = iesta_dom, iesto_dom
    if( lexis(pelty) /= 0 ) then
      maxdof_mag = max(maxdof_mag, mag_edgdof(pelty))
      mxndof_mag = max(mxndof_mag, mag_noddof(pelty))
    end if
  end do

  allocate(belem(maxdof_mag), rnelem(maxdof_mag), felem(maxdof_mag))

  allocate(edgind_mag(maxdof_mag), nodind_mag(maxdof_mag))
  allocate(sigval_mag(maxdof_mag), lenval_mag(maxdof_mag))

  allocate(elstif(maxdof_mag, maxdof_mag))
  allocate(elstdf(maxdof_mag, maxdof_mag))
  allocate(elmass(maxdof_mag, maxdof_mag))

  allocate(elsrc(maxdof_mag))

  allocate(elheat(mxndof_mag))

  allocate(Aelem(maxdof_mag, maxdof_mag))
  allocate(Anelem(maxdof_mag, maxdof_mag))

!    allocate(nodi(ndime, maxdof_mag), nodq(ndime, maxgau_mag), bq(ndime, maxdof_mag))

end subroutine mag_inidof
