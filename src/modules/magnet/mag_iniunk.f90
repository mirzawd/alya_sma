!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine mag_iniunk()
  !-----------------------------------------------------------------------
  !****f* Magnet/mag_iniunk
  ! NAME
  !    mag_iniunk
  ! DECRIPTION
  !    This routine sets up the initial condition
  ! USED BY
  !    Magnet    
  !***
  !-----------------------------------------------------------------------

  use def_master, only: INOTMASTER, timei
  use def_magnet
  use def_domain, only: meshe, ndime
  use def_kermod, only: ndivi
!  use mod_communications, only: PAR_SUM
  use mod_mag_inpdat
  
  implicit none

  integer(ip) :: &
    ipoin,    &
    jpoin,    &
    iedge

  real(rp) :: &
    n1(ndime),    &
    n2(ndime),    &
    c0(ndime),    &
    t0(ndime)

  magnen_mag = 0.0_rp; joulen_mag = 0.0_rp
  !
  ! Initial condition: field
  !
  if (INOTMASTER) then
    !
    do iedge = 1, meshe(ndivi) % nedge
      !
      ! Edge ends
      !
      ipoin = meshe(ndivi) % edge_to_node(1, iedge)
      jpoin = meshe(ndivi) % edge_to_node(2, iedge)
      !
      ! Coordinates of edge ends
      !
      n1 = meshe(ndivi) % coord(1:ndime, ipoin)
      n2 = meshe(ndivi) % coord(1:ndime, jpoin)
      !
      ! Coordinates of edge centroid
      !
      c0 = 0.5_rp * (n1 + n2)
      !
      ! Edge tangent
      !
      t0 = n2 - n1
      t0 = t0 / sqrt( dot_product(t0, t0) ) * edgsig_mag(iedge)
      !
      ! Add Dirichlet values to solution vector
      !
      He_mag(iedge) = dot_product(mag_inifie(c0), t0)
      !
    end do
    !
    ! Energy calculations
    !
    Hp_mag = He_mag
    Hpav_mag => He_mag
    a0_mag = 1.0_rp
    !
    call mag_elmope()
    !
  end if
  !
  call mag_wridat(0_ip, timei)
  !
!!  do pmate = 1, maxmat_mag
!!    call PAR_SUM(magnen_mag(pmate))
!!    call PAR_SUM(joulen_mag(pmate))
!!  end do
!  !
!!  if (INOTSLAVE) then
!!    open(unit=ioun1_mag, file='magnet.nrj')
!!    write(ioun1_mag, 80, advance='no') timei
!!    do pmate = 1, maxmat_mag-1
!!      write(ioun1_mag, 80, advance='no') magnen_mag(pmate)
!!      write(ioun1_mag, 80, advance='no') joulen_mag(pmate)
!!    end do
!!    write(ioun1_mag, 80, advance='no') magnen_mag(maxmat_mag)
!!    write(ioun1_mag, 80) joulen_mag(maxmat_mag)
!!    80 format (e15.8, e15.8, e15.8)
!!    close(ioun1_mag)
!!  end if
!!
!!  70 format (i6)
!!  80 format (e15.8)
  !
end subroutine mag_iniunk
