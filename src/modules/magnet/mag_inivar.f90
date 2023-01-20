!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine mag_inivar(itask)
  !-----------------------------------------------------------------------
  !****f* magnet/mag_inivar
  ! NAME
  !    mag_inivar
  ! DESCRIPTION
  !    This routine initializes some variables
  !    Master does not run this routine!
  ! USES
  ! USED BY
  !    mag_turnon
  !***
  !-----------------------------------------------------------------------

  use def_magnet

  implicit none

  integer(ip), intent(in) :: itask

  select case(itask)

  case(0_ip)
    !
    ! Initialize post-processing
    !
    call mag_inipos()
    !
    ! Initialize solver
    !
    call mag_inisol()
    !
    ! Nullify pointers
    !
    call mag_nulptr()
    !
    ! Initialize materials
    !
    call mag_inimat()
    !
  case(1_ip)
    !
    ! Check axisymmetric option
    !
    call mag_chkaxs()
    !
  case(2_ip)
    !
    ! Initialize quadrature rules
    !
    call mag_iniqua()
    !
    ! Initialize edge data
    !
    call mag_iniedg()
    !
    ! Maximum number of degrees of freedom
    !
    call mag_inidof()
    !
    ! Initialize boundary data for SF calculations
    !
    call mag_inisel()
    !
    ! Initialize interpolation data
    !
    call mag_inintp()
    !
  case(3_ip)
    !
    if (maxdof_mag <= 0_ip) call runend('mag_inivar: maxdof_mag should be greater than zero')
    if (maxgau_mag <= 0_ip) call runend('mag_inivar: maxgau_mag should be greater than zero')
    if (maxmat_mag <= 0_ip) call runend('mag_inivar: maxmat_mag should be greater than zero')
    !
    ! Initialize time loop
    !
    call mag_initim()
    !
  end select
  !
end subroutine mag_inivar
