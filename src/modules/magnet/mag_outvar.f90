!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine mag_outvar(ivari,imesh)
  
  !------------------------------------------------------------------------
  !****f* Master/mag_output
  ! NAME
  !    mag_output
  ! DESCRIPTION
  !    Output a postprocess variable
  ! USES
  !    postpr
  !    memgen
  ! USED BY
  !    mag_output
  !***
  !------------------------------------------------------------------------

  use def_master
  use def_domain, only : ndime
  use def_magnet, only : Hn_mag,  Hc_mag,  Jcz_mag, Jnz_mag
  use def_magnet, only : Bn_mag,  Bc_mag,  Jc_mag,  Jn_mag
  use def_magnet, only : Fc_mag,  Fn_mag,  Hgp_mag, Jgp_mag
  use def_magnet, only : Bgp_mag, Fgp_mag, Gc_mag,  Gn_mag
  use def_magnet, only : Ggp_mag, postev_mag
  use mod_projec, only : projec_elements_to_nodes
  use mod_outvar, only : outvar

  implicit none

  integer(ip), intent(in) :: ivari   !< Variable number
  integer(ip), intent(in) :: imesh   !< Mesh to postprocess on
  !
  ! Define postprocess variable
  !  
  select case ( ivari )

  case (0_ip)
    !
    ! Do nothing
    !
    return

  case (1_ip)
    !
    ! MAGNE: H @ nodes
    !
    if (INOTMASTER) then
      if (.not. postev_mag) then
        call projec_elements_to_nodes(Hc_mag, Hn_mag)
      else
        call projec_elements_to_nodes(ndime, Hgp_mag, Hn_mag)
      end if
      !
      gevec => Hn_mag(:,:)
    end if

  case (2_ip)
    !
    ! CURCZ: Jz @ elements
    !
    if (INOTMASTER) gesca => Jcz_mag

  case (3_ip)
    !
    ! CURNZ: Jz @ nodes 
    !
    if (INOTMASTER) then
      if (.not. postev_mag) then
        call projec_elements_to_nodes(Jcz_mag, Jnz_mag)
      else
        call projec_elements_to_nodes(Jgp_mag, Jnz_mag)
      end if
      !
      gesca => Jnz_mag
    end if

  case (4_ip)
    !
    ! MAGCE: H @ elements
    !
    if (INOTMASTER) gevec => Hc_mag(:,:)

  case (5_ip)
    !
    ! FLUXN: B @ nodes
    !
    if (INOTMASTER) then
      if (.not. postev_mag) then
        call projec_elements_to_nodes(Bc_mag, Bn_mag)
      else
        call projec_elements_to_nodes(ndime, Bgp_mag, Bn_mag)
      end if
      !
      gevec => Bn_mag(:,:)
    end if

  case (6_ip)
    !
    ! FLUXC: B @ elements
    !
    if (INOTMASTER) gevec => Bc_mag(:,:)

  case (7_ip)
    !
    ! CURRC: J @ elements
    !
    if (INOTMASTER) gevec => Jc_mag(:,:)

  case (8_ip)
    !
    ! CURRN: J @ nodes
    !
    if (INOTMASTER) then
      if (.not. postev_mag) then
        call projec_elements_to_nodes(Jc_mag, Jn_mag)
      else
        call projec_elements_to_nodes(ndime, Jgp_mag, Jn_mag)
      end if
      !
      gevec => Jn_mag(:,:)
    end if

  case (9_ip)
    !
    ! FORCN: F @ nodes
    !
    if (INOTMASTER) then
      if (.not. postev_mag) then
        call projec_elements_to_nodes(Fc_mag, Fn_mag)
      else
        call projec_elements_to_nodes(ndime, Fgp_mag, Fn_mag)
      end if
      !
      gevec => Fn_mag(:,:)
    end if

  case (10_ip)
    !
    ! FORCE: F @ elements
    !
    if (INOTMASTER) gevec => Fc_mag(:,:)

  case (11_ip)
    !
    ! JOULN: G @ nodes
    !
    if (INOTMASTER) then
      if (.not. postev_mag) then
        call projec_elements_to_nodes(Gc_mag, Gn_mag)
      else
        call projec_elements_to_nodes(Ggp_mag, Gn_mag)
      end if
      !
      gesca => Gn_mag
    end if

  case (12_ip)
    !
    ! JOULE: G @ elements
    !
    if (INOTMASTER) gesca => Gc_mag

  end select

  call outvar(ivari, ittim, cutim, postp(1) % wopos(:,ivari), MESH_ID=imesh)

end subroutine mag_outvar
