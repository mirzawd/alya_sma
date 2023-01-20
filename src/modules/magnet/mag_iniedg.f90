!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine mag_iniedg()

  use def_master, only: INOTSLAVE, INOTMASTER
  use mod_communications, only: PAR_SUM
  use def_domain, only: ndime, meshe
  use def_magnet
  use def_kermod, only: ndivi
  use mod_memory, only: memory_size
  use mod_mag_elmgeo

  implicit none

  integer(ip) :: &
    ielem,    &
    iedge,    &
    jedge,    &
    ipoin,    &
    jpoin,    &
    iedgeg,    &
    idime,    &
    pelty,    &
    nedgl,    &
    ipoin_loc,    &
    jpoin_loc,    &
    dofsum,    &
    ki,    &
    kj

  real(rp), dimension(ndime) :: nA, nB, nC, nD, nE, nF, nG, nH

  integer(ip) :: &
    locedgnod(2_ip, meshe(ndivi) % medge)
  !
  ! Check 'domain/edge_data_structures.f90' 'domain/mod_graphs.f90' for further details
  !
  ! Total number of Degrees of Freedom (DoF)
  ! nedg1 - Interior edges
  ! nedg2 - First own boundary edge
  ! nedg3 - Last own boundary edge
  ! nedge - Interior + boundary edges
  !
  dofsum = meshe(ndivi) % nedg3
  call PAR_SUM(dofsum)
  if (dofsum == 0_ip) dofsum = meshe(ndivi) % nedge
  if (INOTSLAVE) write(*,*) "DOFs = ", dofsum
  !
  !##############################################################################
  !
  ! Additional element data needed
  !
  do ielem = 1_ip, meshe(ndivi) % nelem
    !
    ! Element type
    ! TRI03: 10
    ! QUA04: 12
    ! TET04: 30
    ! HEX08: 37
    !
    pelty = meshe(ndivi) % ltype(ielem)
    !
    ! Connectivity Element-Edge, element area/volume
    !
    select case (pelty)

    case (10_ip)
      !
      ! TRI03
      !
      nedgl = 3_ip
      !
      if (struct_mag) then
        eleedg_mag(1:3, ielem) = meshe(ndivi) % ledgs([2, 3, 1], ielem)
        ! Local edge-to-node
        ! Edge 1: 2 -> 3
        ! Edge 2: 3 -> 1
        ! Edge 3: 1 -> 2
        locedgnod(1:2, 1:nedgl) = reshape([2, 3, 3, 1, 1, 2], [2, 3])
      else
        eleedg_mag(1:3, ielem) = meshe(ndivi) % ledgs(1:3, ielem)
        ! Local edge-to-node
        ! Edge 1: 1 -> 2
        ! Edge 2: 2 -> 3
        ! Edge 3: 3 -> 1
        locedgnod(1:2, 1:nedgl) = reshape([1, 2, 2, 3, 3, 1], [2, 3])
      end if
      !
      nA = meshe(ndivi) % coord(:, meshe(ndivi) % lnods(1, ielem))
      nB = meshe(ndivi) % coord(:, meshe(ndivi) % lnods(2, ielem))
      nC = meshe(ndivi) % coord(:, meshe(ndivi) % lnods(3, ielem))
      !
      elevol_mag(ielem) = mag_trivol(nA, nB, nC)
      !
    case (12_ip)
      !
      ! QUA04
      !
      nedgl = 4_ip
      !
      if (struct_mag) then
        eleedg_mag(1:4, ielem) = meshe(ndivi) % ledgs(:, ielem)
        ! Local edge-to-node
        locedgnod(1:2, 1:nedgl) = reshape([1, 2, 2, 3, 3, 4, 4, 1], [2, 4])
      else
        eleedg_mag(1:4, ielem) = meshe(ndivi) % ledgs([1, 3, 4, 2], ielem)
        ! Local edge-to-node
        ! Edge 1: 1 -> 2
        ! Edge 2: 4 -> 3
        ! Edge 3: 1 -> 4
        ! Edge 4: 2 -> 3
        locedgnod(1:2, 1:nedgl) = reshape([1, 2, 4, 3, 1, 4, 2, 3], [2, 4])
      end if
      !
      nA = meshe(ndivi) % coord(:, meshe(ndivi) % lnods(1, ielem))
      nB = meshe(ndivi) % coord(:, meshe(ndivi) % lnods(2, ielem))
      nC = meshe(ndivi) % coord(:, meshe(ndivi) % lnods(3, ielem))
      nD = meshe(ndivi) % coord(:, meshe(ndivi) % lnods(4, ielem))
      !
      elevol_mag(ielem) = mag_quavol(nA, nB, nC, nD)
      !
    case (30_ip)
      !
      ! TET04
      !
      nedgl = 6_ip
      !
      eleedg_mag(1:nedgl, ielem) = meshe(ndivi) % ledgs([1, 4, 2, 3, 5, 6], ielem)
      ! Local edge-to-node
      locedgnod(1:2, 1:nedgl) = reshape([1, 2, 2, 3, 3, 1, 4, 1, 4, 2, 4, 3], [2, 6])
      !
      nA = meshe(ndivi) % coord(:, meshe(ndivi) % lnods(1, ielem))
      nB = meshe(ndivi) % coord(:, meshe(ndivi) % lnods(2, ielem))
      nC = meshe(ndivi) % coord(:, meshe(ndivi) % lnods(3, ielem))
      nD = meshe(ndivi) % coord(:, meshe(ndivi) % lnods(4, ielem))
      !
      elevol_mag(ielem) = mag_tetvol(nA, nB, nC, nD)
      !
    case (37_ip)
      !
      ! HEX08
      !
      nedgl = 12_ip
      !
      eleedg_mag(1:nedgl, ielem) = meshe(ndivi) % ledgs([1, 6, 9, 12, 2, 4, 10, 11, 3, 5, 8, 7], ielem)
      ! Local edge-to-node
      locedgnod(1:2, 1:nedgl) = reshape([1,2,4,3,5,6,8,7,1,4,2,3,5,8,6,7,1,5,2,6,4,8,3,7], [2, 12])
      !
      nA = meshe(ndivi) % coord(:, meshe(ndivi) % lnods(1, ielem))
      nB = meshe(ndivi) % coord(:, meshe(ndivi) % lnods(2, ielem))
      nC = meshe(ndivi) % coord(:, meshe(ndivi) % lnods(3, ielem))
      nD = meshe(ndivi) % coord(:, meshe(ndivi) % lnods(4, ielem))
      nE = meshe(ndivi) % coord(:, meshe(ndivi) % lnods(5, ielem))
      nF = meshe(ndivi) % coord(:, meshe(ndivi) % lnods(6, ielem))
      nG = meshe(ndivi) % coord(:, meshe(ndivi) % lnods(7, ielem))
      nH = meshe(ndivi) % coord(:, meshe(ndivi) % lnods(8, ielem))
      !
      elevol_mag(ielem) = mag_hexvol(nA, nB, nC, nD, nE, nF, nG, nH)
      !
    case default
      !
      nedgl = 0_ip
      call runend("mag_iniedg: Unknown element type")
      !
    end select
    !
    ! Parallel solution
    ! Run over element edges
    !
    do iedge = 1_ip, nedgl
      !
      ! Start point: ipoin ----> End point: jpoin
      !
      ! Local nodes
      !
      ipoin_loc = locedgnod(1, iedge)
      jpoin_loc = locedgnod(2, iedge)
      !
      ! Global nodes
      !
      ipoin = meshe(ndivi) % lnods(ipoin_loc, ielem)
      jpoin = meshe(ndivi) % lnods(jpoin_loc, ielem)
      !
      if (meshe(ndivi) % coord(1, jpoin) > meshe(ndivi) % coord(1, ipoin)) then
        elesig_mag(iedge, ielem) = 1.0_rp
      elseif (meshe(ndivi) % coord(1, jpoin) < meshe(ndivi) % coord(1, ipoin)) then
        elesig_mag(iedge, ielem) = -1.0_rp
      else
        if (meshe(ndivi) % coord(2, jpoin) > meshe(ndivi) % coord(2, ipoin)) then
          elesig_mag(iedge, ielem) = 1.0_rp
        elseif (meshe(ndivi) % coord(2, jpoin) < meshe(ndivi) % coord(2, ipoin)) then
          elesig_mag(iedge, ielem) = -1.0_rp
        else
          if (ndime == 2_ip) call runend("mag_inivar: Different nodes cannot have the same coordinates")

          if (meshe(ndivi) % coord(3, jpoin) > meshe(ndivi) % coord(3, ipoin)) then
            elesig_mag(iedge, ielem) = 1.0_rp
          elseif (meshe(ndivi) % coord(3, jpoin) < meshe(ndivi) % coord(3, ipoin)) then
            elesig_mag(iedge, ielem) = -1.0_rp
          else
            call runend("mag_inivar: Different nodes cannot have the same coordinates")
          end if
        end if
      end if
    end do
    !
  end do
  !##############################################################################
  !
  !
  !##############################################################################
  ! Edge Lengths and Global Signs
  !
  do iedgeg = 1, meshe(ndivi) % nedge
    !
    ! Start point: ipoin ----> End point: jpoin
    !
    ipoin = meshe(ndivi) % edge_to_node(1, iedgeg)
    jpoin = meshe(ndivi) % edge_to_node(2, iedgeg)
    !
    ! Edge lengths
    !
    edglen_mag(iedgeg) = 0.0_rp
    do idime = 1, ndime
      edglen_mag(iedgeg) = edglen_mag(iedgeg) + ( meshe(ndivi) % coord(idime, ipoin) - meshe(ndivi) % coord(idime, jpoin) )**2
    end do
    edglen_mag(iedgeg) = sqrt( edglen_mag(iedgeg) )
    !
    ! Global edge signs
    !
    if (meshe(ndivi) % coord(1, jpoin) > meshe(ndivi) % coord(1, ipoin)) then
      edgsig_mag(iedgeg) = 1.0_rp
    elseif (meshe(ndivi) % coord(1, jpoin) < meshe(ndivi) % coord(1, ipoin)) then
      edgsig_mag(iedgeg) = -1.0_rp
    else
      if (meshe(ndivi) % coord(2, jpoin) > meshe(ndivi) % coord(2, ipoin)) then
        edgsig_mag(iedgeg) = 1.0_rp
      elseif (meshe(ndivi) % coord(2, jpoin) < meshe(ndivi) % coord(2, ipoin)) then
        edgsig_mag(iedgeg) = -1.0_rp
      else
        if (ndime == 2_ip) call runend("mag_inivar: Different nodes cannot have the same coordinates")

        if (meshe(ndivi) % coord(3, jpoin) > meshe(ndivi) % coord(3, ipoin)) then
          edgsig_mag(iedgeg) = 1.0_rp
        elseif (meshe(ndivi) % coord(3, jpoin) < meshe(ndivi) % coord(3, ipoin)) then
          edgsig_mag(iedgeg) = -1.0_rp
        else
          call runend("mag_inivar: Different nodes cannot have the same coordinates")
        end if
      end if
    end if
    !
    if (edgsig_mag(iedgeg) == 1.0_rp) then
      edgnod_mag(1, iedgeg) = ipoin
      edgnod_mag(2, iedgeg) = jpoin
    else
      edgnod_mag(1, iedgeg) = jpoin
      edgnod_mag(2, iedgeg) = ipoin
    end if
    !
  end do
  !##############################################################################
  !
  !
  !##############################################################################
  if (INOTMASTER) then
    !
    allocate(node_to_nedge_mag(meshe(ndivi) % npoin))
    allocate(r_node_to_edge_mag(meshe(ndivi) % npoin + 1_ip))
    allocate(c_node_to_edge_mag(meshe(ndivi) % nedge * 2_ip))
    allocate(i_node_to_edge_mag(meshe(ndivi) % npoin + 1_ip))
    allocate(edge_to_nedge_mag(meshe(ndivi) % nedge))
    allocate(r_edge_to_edge_mag(meshe(ndivi) % nedge + 1_ip))
    allocate(i_edge_to_edge_mag(meshe(ndivi) % nedge + 1_ip))
    !
    !
    ! Get number of edges per node:
    ! node_to_nedge_mag has size npoin
    !
    node_to_nedge_mag = 0_ip
    do iedge = 1, meshe(ndivi) % nedge
      ipoin = meshe(ndivi) % edge_to_node(1, iedge)
      jpoin = meshe(ndivi) % edge_to_node(2, iedge)

      node_to_nedge_mag(ipoin) = node_to_nedge_mag(ipoin) + 1_ip
      node_to_nedge_mag(jpoin) = node_to_nedge_mag(jpoin) + 1_ip
    end do
    !
    !
    ! Get edges per node:
    ! r_node_to_node_mag has size npoin+1
    ! c_node_to_node_mag has size 2*nedge
    ! i_node_to_node_mag has size npoin+1
    !
    r_node_to_edge_mag(1) = 1_ip
    do ipoin = 1, meshe(ndivi) % npoin
      r_node_to_edge_mag(ipoin + 1) = r_node_to_edge_mag(ipoin) + node_to_nedge_mag(ipoin)
    end do
    i_node_to_edge_mag = r_node_to_edge_mag
    !
    do iedge = 1, meshe(ndivi) % nedge
      ipoin = meshe(ndivi) % edge_to_node(1, iedge)
      jpoin = meshe(ndivi) % edge_to_node(2, iedge)

      ki = i_node_to_edge_mag(ipoin)
      c_node_to_edge_mag(ki) = iedge
      i_node_to_edge_mag(ipoin) = ki + 1_ip

      kj = i_node_to_edge_mag(jpoin)
      c_node_to_edge_mag(kj) = iedge
      i_node_to_edge_mag(jpoin) = kj + 1_ip
    end do
    !
    !
    ! Get number of edges per edge
    ! edge_to_nedge has size nedge
    !
    ne2e_mag = 0
    do iedge = 1, meshe(ndivi) % nedge
      ipoin = meshe(ndivi) % edge_to_node(1, iedge)
      jpoin = meshe(ndivi) % edge_to_node(2, iedge)

      edge_to_nedge_mag(iedge) = r_node_to_edge_mag(ipoin + 1) - r_node_to_edge_mag(ipoin) - 1_ip
      edge_to_nedge_mag(iedge) = edge_to_nedge_mag(iedge) + r_node_to_edge_mag(jpoin + 1) - r_node_to_edge_mag(jpoin) - 1_ip

      ne2e_mag = ne2e_mag + edge_to_nedge_mag(iedge)
    end do
    !
    allocate(c_edge_to_edge_mag(ne2e_mag))
    !
    !
    ! Get edges per edge:
    ! r_edge_to_edge_mag has size nedge+1
    ! c_edge_to_edge_mag has size ne2e_mag
    ! i_edge_to_edge_mag has size nedge+1
    !
    r_edge_to_edge_mag(1) = 1_ip
    do iedge = 1, meshe(ndivi) % nedge
      r_edge_to_edge_mag(iedge + 1) = r_edge_to_edge_mag(iedge) + edge_to_nedge_mag(iedge)
    end do
    i_edge_to_edge_mag = r_edge_to_edge_mag
    !
    do iedge = 1, meshe(ndivi) % nedge
      ipoin = meshe(ndivi) % edge_to_node(1, iedge)
      jpoin = meshe(ndivi) % edge_to_node(2, iedge)

      do jedge = r_node_to_edge_mag(ipoin), r_node_to_edge_mag(ipoin+1) - 1_ip
        if (c_node_to_edge_mag(jedge) /= iedge) then
          c_edge_to_edge_mag(i_edge_to_edge_mag(iedge)) = c_node_to_edge_mag(jedge)
          i_edge_to_edge_mag(iedge) = i_edge_to_edge_mag(iedge) + 1_ip
        end if
      end do

      do jedge = r_node_to_edge_mag(jpoin), r_node_to_edge_mag(jpoin+1) - 1_ip
        if (c_node_to_edge_mag(jedge) /= iedge) then
          c_edge_to_edge_mag(i_edge_to_edge_mag(iedge)) = c_node_to_edge_mag(jedge)
          i_edge_to_edge_mag(iedge) = i_edge_to_edge_mag(iedge) + 1_ip
        end if
      end do
    end do
    !
  end if
  !
end subroutine mag_iniedg
