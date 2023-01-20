!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine mag_inisel()

  use def_master, only: ISEQUEN
  use def_domain, only: ndime, meshe
  use def_magnet
  use def_kermod, only: ndivi
  use mod_memory, only: memory_size
  use mod_communications, only: PAR_ALLGATHERV

  implicit none

  integer(ip) :: iedge, ipoin, jpoin, iedgeg

  real(rp) :: n1(ndime), n2(ndime), t0(ndime), c0(ndime)
  !
  ! Communicate boundary edge centroids and tangents for:
  ! - Dirichlet Boundary Conditions
  ! - Self-field calculations
  !
  do iedge = 1_ip, memory_size(diredg_mag)
    !
    iedgeg = diredg_mag(iedge)
    !
    ipoin = meshe(ndivi) % edge_to_node(1, iedgeg)
    jpoin = meshe(ndivi) % edge_to_node(2, iedgeg)
    !
    n1(1:ndime) = meshe(ndivi) % coord(1:ndime, ipoin)
    n2(1:ndime) = meshe(ndivi) % coord(1:ndime, jpoin)
    !
    ! Edge centroid
    !
    c0 = 0.5_rp * (n1 + n2)
    !
    ! Edge tangent
    !
    t0 = n2 - n1
    t0 = t0 / sqrt( dot_product(t0, t0) ) * edgsig_mag(iedgeg)
    !
    ! Local edge boundary data
    !
    locboucen_mag(1:ndime, iedge) = c0
    locboutan_mag(1:ndime, iedge) = t0
    !
  end do
  !
  ! Generate global edge boundary data
  !
  if (ISEQUEN) then
    globoucen_mag = locboucen_mag
    globoutan_mag = locboutan_mag
  else
    call PAR_ALLGATHERV(locboucen_mag, globoucen_mag, proc_nbedg, 'IN MY CODE') 
    call PAR_ALLGATHERV(locboutan_mag, globoutan_mag, proc_nbedg, 'IN MY CODE') 
  end if
  !
  proc_nbedg = proc_nbedg / ndime
  !
end subroutine mag_inisel
