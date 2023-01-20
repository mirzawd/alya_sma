!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



module mod_mag_lagran

  use def_kintyp, only: ip, rp

  implicit none

contains

  !##############################################################
  function mag_noddof(pelty) result(numdof)

    implicit none

    integer(ip), intent(in) :: pelty

    integer(ip) :: numdof

    select case (pelty)

    case (1_ip)
      numdof = 0_ip

    case (10_ip)
      ! TRI03
      numdof = 3_ip

    case (12_ip)
      ! QUA04
      numdof = 4_ip

    case (30_ip)
      ! TET04
      numdof = 4_ip

    case (37_ip)
      ! HEX08
      numdof = 8_ip

    case default
      numdof = 0_ip
      call runend("mag_noddof: error")

    end select

  end function mag_noddof
  !##############################################################


  !##############################################################
  function mag_lagref(x, eletyp) result(basis)

    implicit none

    real(rp), intent(in) ::             x(:,:) ! Point coordinates in reference domain

    integer(ip), intent(in) ::          eletyp ! Element type

    real(rp), allocatable ::            basis(:,:)

    integer(ip) ::                      ipoin

    select case (eletyp)

    case (10_ip)
      !
      ! TRI03
      !
      allocate(basis(3, size(x, 2)))
      do ipoin = 1, size(x, 2)
        basis(:, ipoin) = [ 1 - x(1, ipoin) - x(2, ipoin), &
                            x(1, ipoin), &
                            x(2, ipoin) ]
      end do

    case (12_ip)
      !
      ! QUA04
      !
      allocate(basis(4, size(x, 2)))
      do ipoin = 1, size(x, 2)
        basis(:, ipoin) = 0.25_rp * [ (1 - x(1, ipoin)) * (1 - x(2, ipoin)), &
                                      (1 + x(1, ipoin)) * (1 - x(2, ipoin)), &
                                      (1 + x(1, ipoin)) * (1 + x(2, ipoin)), &
                                      (1 - x(1, ipoin)) * (1 + x(2, ipoin)) ]
      end do

    case (30_ip)
      !
      ! TET04
      !
      allocate(basis(4, size(x, 2)))
      do ipoin = 1, size(x, 2)
        basis(:, ipoin) = [ 1 - x(1, ipoin) - x(2, ipoin) - x(3,ipoin), &
                            x(1, ipoin), &
                            x(2, ipoin), &
                            x(3, ipoin) ]
      end do

    case (37_ip)
      !
      ! HEX08
      !
      allocate(basis(8, size(x, 2)))
      do ipoin = 1, size(x, 2)
        basis(:, ipoin) = 0.125_rp * [ (1 - x(1, ipoin)) * (1 - x(2, ipoin)) * (1 - x(3, ipoin)), &
                                       (1 + x(1, ipoin)) * (1 - x(2, ipoin)) * (1 - x(3, ipoin)), &
                                       (1 + x(1, ipoin)) * (1 + x(2, ipoin)) * (1 - x(3, ipoin)), &
                                       (1 - x(1, ipoin)) * (1 + x(2, ipoin)) * (1 - x(3, ipoin)), &
                                       (1 - x(1, ipoin)) * (1 - x(2, ipoin)) * (1 + x(3, ipoin)), &
                                       (1 + x(1, ipoin)) * (1 - x(2, ipoin)) * (1 + x(3, ipoin)), &
                                       (1 + x(1, ipoin)) * (1 + x(2, ipoin)) * (1 + x(3, ipoin)), &
                                       (1 - x(1, ipoin)) * (1 + x(2, ipoin)) * (1 + x(3, ipoin)) ]
      end do

    case default
      call runend("mag_lagref: unknown element type")

    end select

  end function mag_lagref
  !##############################################################


  !##############################################################
  function mag_lagdif(x, eletyp) result(basis)

    implicit none

    real(rp), intent(in) ::             x(:) ! Point coordinates in reference domain

    integer(ip), intent(in) ::          eletyp ! Element type

    real(rp), allocatable ::            basis(:,:)

    real(rp) :: m1 = -1.0_rp, p1 = 1.0_rp, z0 = 0.0_rp

    select case (eletyp)

    case (10_ip)
      !
      ! TRI03
      !
      allocate(basis(3, 2))
      basis = reshape( [m1, p1, z0, m1, z0, p1], [3, 2] )
!      basis = reshape( [-1.0_rp, 1.0_rp, 0.0_rp, -1.0_rp, 0.0_rp, 1.0_rp], [3, 2] )

    case (12_ip)
      !
      ! QUA04
      !
      allocate(basis(4, 2))
      basis = reshape( [x(2) - 1.0_rp, 1.0_rp - x(2), 1.0_rp + x(2), -1.0_rp - x(2), &
                        x(1) - 1.0_rp, -1.0_rp - x(1), 1.0_rp + x(1), 1.0_rp - x(1)], &
                       [4, 2] ) * 0.25_rp

    case (30_ip)
      !
      ! TET04
      !
      allocate(basis(4, 3))
      basis = reshape( [-1.0_rp, 1.0_rp, 0.0_rp, 0.0_rp, &
                        -1.0_rp, 0.0_rp, 1.0_rp, 0.0_rp, &
                        -1.0_rp, 0.0_rp, 0.0_rp, 1.0_rp ], &
                       [4, 3] )

    case (37_ip)
      allocate(basis(8, 3))
      basis = reshape( [-(p1 - x(2))*(p1 - x(3)), (p1 - x(2))*(p1 - x(3)), (p1 + x(2))*(p1 - x(3)), -(p1 + x(2))*(p1 - x(3)), &
                        -(p1 - x(2))*(p1 + x(3)), (p1 - x(2))*(p1 + x(3)), (p1 + x(2))*(p1 + x(3)), -(p1 + x(2))*(p1 + x(3)), &
                        -(p1 - x(1))*(p1 - x(3)), -(p1 + x(1))*(p1 - x(3)), (p1 + x(1))*(p1 - x(3)), (p1 - x(1))*(p1 - x(3)), &
                        -(p1 - x(1))*(p1 + x(3)), -(p1 + x(1))*(p1 + x(3)), (p1 + x(1))*(p1 + x(3)), (p1 - x(1))*(p1 + x(3)), &
                        -(p1 - x(1))*(p1 - x(2)), -(p1 + x(1))*(p1 - x(2)), -(p1 + x(1))*(p1 + x(2)), -(p1 - x(1))*(p1 + x(2)), &
                        (p1 - x(1))*(p1 - x(2)), (p1 + x(1))*(p1 - x(2)), (p1 + x(1))*(p1 + x(2)), (p1 - x(1))*(p1 + x(2))], &
                        [8, 3] ) * 0.125_rp

    case default
      call runend("mag_lagdif: unknown element type")

    end select

  end function mag_lagdif
  !##############################################################

end module mod_mag_lagran
