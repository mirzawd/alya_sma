!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine mag_wridat(pos, time)

  use def_master, only: INOTSLAVE !, timei
  use def_magnet
  use mod_communications, only: PAR_SUM
  use def_domain, only: ndime

  implicit none

  integer(ip), intent(in) :: pos

  real(rp), intent(in) :: time

  integer(ip) :: pmate, idime

  !
  ! Sum magnetic energy, heat dissipation and magnetization for each material
  !
  do pmate = 1, maxmat_mag
    if (kfl_nrj_mag) then
      call PAR_SUM(magnen_mag(pmate))
      call PAR_SUM(joulen_mag(pmate))
    end if

    if (kfl_vlm_mag) call PAR_SUM(volume_mag(pmate))

    do idime = 1, ndime
      if (kfl_mtz_mag) call PAR_SUM(magtiz_mag(idime, pmate))
      if (kfl_crn_mag) call PAR_SUM(cursum_mag(idime, pmate))
      if (kfl_mtz_mag) call PAR_SUM(magsum_mag(idime, pmate))
    end do
  end do
  !
  if (INOTSLAVE) then
    !
    ! Open magnet.nrj file
    ! Open magnet.mtz file
    ! Open magnet.vlm file
    ! Open magnet.crn file
    !
    if (pos == 0_ip) then
      if (kfl_nrj_mag) open(unit=ioun1_mag, file='magnet.nrj')
      if (kfl_mtz_mag) open(unit=ioun2_mag, file='magnet.mtz')
      if (kfl_vlm_mag) open(unit=ioun3_mag, file='magnet.vlm')
      if (kfl_crn_mag) open(unit=ioun4_mag, file='magnet.crn')
    elseif (pos == 1_ip) then
      if (kfl_nrj_mag) open(unit=ioun1_mag, file='magnet.nrj', position='append')
      if (kfl_mtz_mag) open(unit=ioun2_mag, file='magnet.mtz', position='append')
      if (kfl_vlm_mag) open(unit=ioun3_mag, file='magnet.vlm', position='append')
      if (kfl_crn_mag) open(unit=ioun4_mag, file='magnet.crn', position='append')
    end if
    !
    ! Write data to file
    !
    if (kfl_nrj_mag) write(ioun1_mag, 80, advance='no') time
    if (kfl_vlm_mag) write(ioun3_mag, 80, advance='no') time
    do pmate = 1_ip, maxmat_mag - 1_ip
      if (kfl_nrj_mag) then
        write(ioun1_mag, 80, advance='no') magnen_mag(pmate)
        write(ioun1_mag, 80, advance='no') joulen_mag(pmate)
      end if
      if (kfl_vlm_mag) write(ioun3_mag, 80, advance='no') volume_mag(pmate)
    end do
    if (kfl_nrj_mag) then
      write(ioun1_mag, 80, advance='no') magnen_mag(maxmat_mag)
      write(ioun1_mag, 80) joulen_mag(maxmat_mag)
    end if
    if (kfl_vlm_mag) write(ioun3_mag, 80) volume_mag(maxmat_mag)
    !
    if (kfl_mtz_mag) write(ioun2_mag, 90, advance='no') time
    if (kfl_crn_mag) write(ioun4_mag, 90, advance='no') time
    do pmate = 1_ip, maxmat_mag - 1_ip
      do idime = 1_ip, 3_ip
        if (kfl_mtz_mag) write(ioun2_mag, 90, advance='no') magsum_mag(idime, pmate)
        if (kfl_crn_mag) write(ioun4_mag, 90, advance='no') cursum_mag(idime, pmate)
      end do
    end do
    if (kfl_mtz_mag) then
      write(ioun2_mag, 90, advance='no') magsum_mag(1_ip, maxmat_mag)
      write(ioun2_mag, 90, advance='no') magsum_mag(2_ip, maxmat_mag)
      write(ioun2_mag, 90) magsum_mag(3_ip, maxmat_mag)
    end if
    !
    if (kfl_crn_mag) then
      write(ioun4_mag, 90, advance='no') cursum_mag(1_ip, maxmat_mag)
      write(ioun4_mag, 90, advance='no') cursum_mag(2_ip, maxmat_mag)
      write(ioun4_mag, 90) cursum_mag(3_ip, maxmat_mag)
    end if
    !
    ! Close files
    !
    if (kfl_nrj_mag) close(ioun1_mag)
    if (kfl_mtz_mag) close(ioun2_mag)
    if (kfl_vlm_mag) close(ioun3_mag)
    if (kfl_crn_mag) close(ioun4_mag)
    !
  end if

  80 format (e18.8)
  90 format (e18.8)

end subroutine mag_wridat
