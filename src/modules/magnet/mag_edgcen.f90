!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine mag_edgcen()
  !-----------------------------------------------------------------------
  !****f*  magnet/mag_edgcen
  ! NAME
  !    mag_edgcen
  ! DESCRIPTION
  !    This routine computes the solution field components at the element centroids
  ! USES
  ! USED BY
  !    mag_
  !***
  !-----------------------------------------------------------------------

  use def_domain, only: meshe, ndime
  use def_magnet
  use def_kermod, only: ndivi
  use mod_mag_nedele
  use mod_mag_linalg, only: mag_matvec, mag_vecvec
  use mod_mag_lagran
  use mod_mag_matpro

  implicit none

  integer(ip) :: &
    ielem,    &
    eledof,    &
    pelty,    &
    pmate

  real(rp) :: &
    nodi(1:ndime, 1:maxdof_mag),    &
    basis(1:ndime, 1:maxdof_mag),    &
    curl2(1:maxdof_mag),    &
    curl3(1:3, 1:maxdof_mag),    &
    centro(1:ndime),    &
!    nodq(ndime, maxgau_mag),    &
    elevol,    &
    detJ,    &
    tMUR(3,3),    &
    tRES(3,3),    &
    Jc(3),    &
    Hc(1:ndime),    &
    Bc(1:ndime)

  do ielem = 1, meshe(ndivi) % nelem
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
    pelty = meshe(ndivi) % ltype(ielem)
    !
    ! Material
    !
    pmate = meshe(ndivi) % lmate(ielem)
    !
    ! Number of degrees of freedom in current element
    !
    eledof = mag_edgdof(pelty)
    !
    ! Element nodes: indices
    !
    nodind_mag(1:eledof) = meshe(ndivi) % lnods(1:eledof, ielem)
    !
    ! Element edges: indices
    !
    edgind_mag(1:eledof) = eleedg_mag(1:eledof, ielem)
    !
    ! Edge signs
    !
    sigval_mag(1:eledof) = elesig_mag(1:eledof, ielem)
    !
    ! Vertices coordinates
    !
    nodi(1:ndime, 1:eledof) = meshe(ndivi) % coord(1:ndime, nodind_mag(1:eledof))
    !
    ! Centroid coordinates
    !
    centro = sum(nodi, 2) / real(size(nodi, 2), rp)
    !
    ! Basis functions
    !
    select case (pelty)
    case (10_ip)
      !
      ! TRI03
      !
      if (struct_mag) then
        !
        ! Basis
        !
        basis(1:2, 1:eledof) = mag_nedtri(centro, nodi)
        !
      else
        call mag_nedtri2([1.0_rp, 1.0_rp] / 3.0_rp, nodi, basis, curl2, detJ)
      end if
      !
      tMUR(1:2, 1:2) = reshape( [ murmat_mag(1, pmate), 0.0_rp, &
                                  0.0_rp, murmat_mag(2, pmate) ], &
                                          [2, 2] &
                                      )
      !
    case (12_ip)
      !
      ! QUA04
      !
      if (struct_mag) then
        basis(1:2, 1:eledof) = mag_nedrce(nodi)
      else
        call mag_nedqua([0.0_rp, 0.0_rp], nodi, basis, curl2, detJ)  
      end if
      !
      tMUR(1:2, 1:2) = reshape( [ murmat_mag(1, pmate), 0.0_rp, &
                                0.0_rp, murmat_mag(2, pmate) ], &
                                [2, 2] &
                              )
      !
    case (30_ip)
      !
      ! TET04
      !
      call mag_nedtet([0.0_rp, 0.0_rp, 0.0_rp], nodi, basis, curl3, detJ)
      !
      tMUR(1:ndime, 1:ndime) = reshape( [ murmat_mag(1, pmate), 0.0_rp, 0.0_rp, &
                                          0.0_rp, murmat_mag(2, pmate), 0.0_rp, &
                                          0.0_rp, 0.0_rp, murmat_mag(3, pmate) ], &
                                        [3, 3] &
                                      )
      !
    case (37_ip)
      !
      ! HEX08
      !
      call mag_nedhex([0.0_rp, 0.0_rp, 0.0_rp], nodi, basis, curl3, detJ)
      !
      tMUR(1:ndime, 1:ndime) = reshape( [ murmat_mag(1, pmate), 0.0_rp, 0.0_rp, &
                                          0.0_rp, murmat_mag(2, pmate), 0.0_rp, &
                                          0.0_rp, 0.0_rp, murmat_mag(3, pmate) ], &
                                        [3, 3] &
                                      )
      !
    end select
    !
    Hc = matmul( basis(1:ndime, 1:eledof), He_mag(edgind_mag(1:eledof)) * sigval_mag(1:eledof) )
    Bc = mu0_mag * matmul(tMUR(1:ndime, 1:ndime), Hc(1:ndime))

    Hc_mag(1:ndime, ielem) = Hc
    Bc_mag(1:ndime, ielem) = Bc
    !
    !--------------------------------------------------------------------------------
    ! Elementary Current Density
    !--------------------------------------------------------------------------------
    if (ndime == 2_ip .and. struct_mag) then
      !
      ! Element is TRI03 or QUA04 in structured mesh
      ! 
      elevol = elevol_mag(ielem)
      !
      ! Edge lengths
      !
      lenval_mag(1:eledof) = edglen_mag(edgind_mag(1:eledof)) * sigval_mag(1:eledof)
      !
      Jcz_mag(ielem) = dot_product(lenval_mag(1:eledof), He_mag(edgind_mag(1:eledof))) / elevol
      !
      Fc_mag(1:2, ielem) = Jcz_mag(ielem) * [ -Bc(2), Bc(1) ]
      !
      Jc = [0.0_rp, 0.0_rp, Jcz_mag(ielem)]
      !
      Gc_mag(ielem) = mag_resist(pmate, centro, Jc, Hc, 0.d0, .false., 3_ip) * Jcz_mag(ielem) * Jcz_mag(ielem)
      !
    elseif (ndime == 2_ip) then
      !
      ! Element is not TRI03 nor QUA04 in structured mesh and mesh is 2D
      !
      Jcz_mag(ielem) = dot_product(curl2(1:eledof), He_mag(edgind_mag(1:eledof)) * sigval_mag(1:eledof) )
      !
      Fc_mag(1:2, ielem) = Jcz_mag(ielem) * [ -Bc(2), Bc(1) ]
      !
      Jc = [0.0_rp, 0.0_rp, Jcz_mag(ielem)]
      !
      Gc_mag(ielem) = mag_resist(pmate, centro, Jc, Hc, 0.d0, .false., 3_ip) * Jcz_mag(ielem) * Jcz_mag(ielem)
      !
    else
      !
      ! 3D mesh
      !
      Jc = mag_matvec(curl3(1:ndime, 1:eledof), He_mag(edgind_mag(1:eledof)) * sigval_mag(1:eledof))
      !
      Jc_mag(1:ndime, ielem) = Jc
      !
      Fc_mag(1:ndime, ielem) = mag_vecvec(Jc, Bc)
      !
      tRES = reshape( [ mag_resist(pmate, centro, Jc, Hc, 0.d0, .false., 1_ip), 0.0_rp, 0.0_rp, &
               0.0_rp, mag_resist(pmate, centro, Jc, Hc, 0.d0, .false., 2_ip), 0.0_rp, &
               0.0_rp, 0.0_rp, mag_resist(pmate, centro, Jc, Hc, 0.d0, .false., 3_ip) ], &
               [3, 3] &
               )
      !
      Gc_mag(ielem) = dot_product(Jc, matmul(tRES(1:ndime, 1:ndime), Jc))
      !
    end if
    !--------------------------------------------------------------------
    !
  end do
  !
end subroutine mag_edgcen
