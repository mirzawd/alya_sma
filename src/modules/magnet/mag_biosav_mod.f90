!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine mag_biosav_mod()
  !-----------------------------------------------------------------------
  !****f* magnet/mag_elmope.f90
  ! NAME
  !    mag_biosav_mod
  ! DESCRIPTION
  !    This routine
  !    1. Computes the elementary matrix and RHS for every element in the mesh
  !    2. Assembles the elementary equations into the system equations
  ! USES
  ! USED BY
  !    mag_matrix
  !***
  !-----------------------------------------------------------------------
  use def_parame, only: pi
  use def_domain
  use def_magnet
  use def_kermod, only: ndivi
  use mod_memory, only: memory_size
  use mod_mag_nedele, only: mag_edgdof, mag_nedhex, mag_nedtet, mag_nedqua, mag_nedtri2
  use mod_mag_lagran
  use mod_mag_linalg, only: mag_matvec, mag_vecvec
  
  implicit none

  integer(ip) :: &
    eledof,    &
    elegau,    &
    pelty,    &
    pmate,    &
    ielem,    &
    ipoin,    &
    iquad,    &
    idime

  real(rp) :: &
    nodi(1:ndime, 1:maxdof_mag),    &
    nodq(1:ndime, 1:maxgau_mag),    &
    wq(1:maxgau_mag),    &
    r2,    &
    dx(1:ndime),    &
    bq(1:ndime, maxdof_mag),    &
    curl2(maxdof_mag),    &
    curl3(ndime, maxdof_mag),    &
    Jq(1:ndime, 1:maxgau_mag),    &
    detJ(1:maxgau_mag),    &
    dist,    &
    currDens(1:maxgau_mag)
  !
  Hbs_mag(:,:) = 0.0_rp
  !
  do ielem = 1_ip, nelem
    !
    pmate = lmate(ielem)
    !
    if (.not. selfList_mag(pmate)) cycle
    !
    !---------------------------------------------------------------------------------
    ! Element Data
    !---------------------------------------------------------------------------------
    !
    ! Element type
    ! TRI03: 10
    ! QUA04: 12
    ! TET04: 30
    ! HEX08: 37
    !
    !--------------------------------------------------------------------------------
    pelty = meshe(ndivi) % ltype(ielem)
    !
    ! Number of degrees of freedom in current element
    !
    eledof = mag_edgdof(pelty)
    !
    ! Element nodes
    !
    nodind_mag(1:eledof) = meshe(ndivi) % lnods(1:eledof, ielem)
    !
    ! Element edges
    !
    edgind_mag(1:eledof) = eleedg_mag(1:eledof, ielem)
    !
    ! Edge signs
    !
    sigval_mag(1:eledof) = elesig_mag(1:eledof, ielem)
    !
    ! Edge lengths
    !
    lenval_mag(1:eledof) = edglen_mag(edgind_mag(1:eledof)) * sigval_mag(1:eledof)
    !
    ! Coordinates of element nodes
    !
    nodi(1:ndime, 1:eledof) = coord(1:ndime, nodind_mag(1:eledof))
    ! 
    if (pelty == 10) then
      !
      ! TRI03
      !
      elegau = quadTri % nq
      wq(1:elegau) = quadTri % wq
      !
      ! Coordinates of gauss nodes in the reference element
      !
      nodq(1:ndime, 1:elegau) = matmul(nodi(1:ndime, 1:eledof), mag_lagref(quadTri % xq, pelty))
      !
    elseif (pelty == 12) then
      !
      ! QUA04
      !
      elegau = quadQua % nq
      wq(1:elegau) = quadQua % wq
      !
      ! Coordinates of Gauss quadrature nodes
      !
      nodq(1:ndime, 1:elegau) = matmul(nodi(1:ndime, 1:eledof), mag_lagref(quadQua % xq, pelty))
      !
    elseif (pelty == 30_ip) then
      !
      ! TET04
      !
      elegau = quadTet % nq
      wq(1:elegau) = quadTet % wq
      !
      ! Coordinates of Gauss quadrature nodes
      !
      nodq(1:ndime, 1:elegau) = matmul(nodi(1:ndime, 1:eledof), mag_lagref(quadTet % xq, pelty))
      !
    elseif (pelty == 37_ip) then
      !
      ! HEX08
      !
      elegau = quadHex % nq
      wq(1:elegau) = quadHex % wq
      !
      ! Coordinates of Gauss quadrature nodes
      !
      nodq(1:ndime, 1:elegau) = matmul(nodi(1:ndime, 1:eledof), mag_lagref(quadHex % xq, pelty))
      !
    else
      call runend("mag_elmope: Unknown element type")
    end if
    !
    do iquad = 1, elegau
      !
      if (ndime == 2_ip) then
        !
        if (pelty == 10_ip) then
          !
          call mag_nedtri2(quadTri % xq(1:ndime, iquad), nodi, bq, curl2, detJ(iquad))
          !
        elseif (pelty == 12_ip) then
          !
          call mag_nedqua(quadQua % xq(1:ndime, iquad), nodi, bq, curl2, detJ(iquad))
          !
        end if
        !
        currDens(iquad) = dot_product(curl2(1:eledof), He_mag(edgind_mag(1:eledof)) * sigval_mag(1:eledof))
        !
      elseif (ndime == 3_ip) then
        !
        if (pelty == 30_ip) then
          !
          call mag_nedtet(quadTet % xq(1:ndime, iquad), nodi, bq, curl3, detJ(iquad))
          !
        elseif (pelty == 37_ip) then
          !
          call mag_nedhex(quadHex % xq(1:ndime, iquad), nodi, bq, curl3, detJ(iquad))
          !
        end if
        !
        ! Current density at quadrature node
        !
        Jq(1:ndime, iquad) = mag_matvec(curl3(1:ndime, 1:eledof), He_mag(edgind_mag(1:eledof)) * sigval_mag(1:eledof))
        !
      else
        call runend("mag_biosav_mod: number of dimensions is not correct")
      end if
      !
    end do    
    !
    do ipoin = 1, memory_size(globoucen_mag, 2_ip)
      !
      do iquad = 1, elegau
        !
        ! Vector from quadrature node to point
        !
        r2 = 0.0_rp
        do idime = 1, ndime
          dx(idime) = globoucen_mag(idime, ipoin) - nodq(idime, iquad)
          r2 = r2 + dx(idime)**2
        end do
        if (r2 <= 1.0e-14_rp) cycle
        !
        ! Compute field
        !
        if (ndime == 2_ip) then
          !
          Hbs_mag(1:2, ipoin) = Hbs_mag(1:2, ipoin) + wq(iquad) * currDens(iquad) * abs(detJ(iquad)) / (2.0_rp * pi * r2) * [ -dx(2), dx(1) ]
          !
        elseif (ndime == 3_ip) then
          !
          dist = sqrt(r2)
          !
          Hbs_mag(1:ndime, ipoin) = Hbs_mag(1:ndime, ipoin) + wq(iquad) / (4.0_rp * pi * dist**3) * mag_vecvec(Jq(1:ndime, iquad), dx) * abs(detJ(iquad))
          !
        else
          !
          call runend("mag_biosav_mod: number of dimensions is not correct")
          !
        end if
        !
      end do
      !
    end do
    !
  end do
  !
end subroutine mag_biosav_mod
