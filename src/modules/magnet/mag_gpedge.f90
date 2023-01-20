!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine mag_gpedge()
  !-----------------------------------------------------------------------
  !****f*  magnet/mag_edgcen
  ! NAME
  !    mag_edgcen
  ! DESCRIPTION
  !    This routine computes the solution field components at every element quadrature nodes
  ! USES
  ! USED BY
  !    mag_
  !***
  !-----------------------------------------------------------------------

  use def_domain, only: meshe, ndime, ngaus, elmar
  use def_magnet
  use def_kermod, only: ndivi
  use mod_mag_lagran
  use mod_mag_nedele
  use mod_mag_linalg, only: mag_matvec, mag_vecvec
  use mod_mag_matpro

  implicit none

  integer(ip) :: &
    ielem,    &
    igaus,    &
    eledof,    &
    pelty,    &
    pmate,    &
    elegau

  real(rp) :: &
    nodi(1:ndime, 1:maxdof_mag),    &
    basis(1:ndime, 1:maxdof_mag),    &
    curl2(1:maxdof_mag),    &
    curl3(1:3, 1:maxdof_mag),    &
    detJ,    &
    nodq(ndime, meshe(ndivi) % mgaus),    &
    tMUR(ndime, ndime),    &
    tRES(ndime, ndime),    &
    Hgp(1:ndime),    &
    Bgp(1:ndime),    &
    Jgp(3),    &
    Jzgp,    &
    rc

  do ielem = 1, meshe(ndivi) % nelem
    !
    !---------------------------------------------------------------------------------
    ! Element Data
    !---------------------------------------------------------------------------------
    !
    ! Element type (pelty):
    ! TRI03: 10
    ! QUA04: 12
    ! TET04: 30
    ! HEX08: 37
    !
    pelty = meshe(ndivi) % ltype(ielem)
    !
    ! Number of Gauss nodes
    !
    elegau = ngaus(pelty)
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
    ! Coordinates of gauss nodes in the reference element
    !
    nodq(1:ndime, 1:elegau) = matmul(nodi(1:ndime, 1:eledof), mag_lagref(elmar(pelty) % posgp(1:ndime, 1:elegau), pelty))
    !
    ! Radial coordinate for axisymmetric analysis
    !
    rc = 0.0_rp
    !
    ! Basis functions
    !
    select case (pelty)
    
    case (10_ip)
      !
      ! TRI03
      !
!      tMUR(1:2, 1:2) = reshape( [ murmat_mag(1, pmate), 0.0_rp, 0.0_rp, murmat_mag(2, pmate) ], [2, 2] )
      !
      do igaus = 1, elegau
        !
        call mag_nedtri2(nodq(1:ndime, igaus), nodi, basis, curl2, detJ)
        !
        Hgp = mag_matvec(basis(1:ndime, 1:eledof), He_mag(edgind_mag(1:eledof)) * sigval_mag(1:eledof))
        Jzgp = dot_product(curl2(1:eledof), He_mag(edgind_mag(1:eledof)) * sigval_mag(1:eledof))
        !
        Jgp = [0.0_rp, 0.0_rp, Jzgp]
        !
        tMUR(1:2, 1:2) = reshape( [ mag_permea(pmate, nodq(1:ndime, igaus), Hgp, .false., 1_ip), 0.0_rp, &
                                    0.0_rp, mag_permea(pmate, nodq(1:ndime, igaus), Hgp, .false., 2_ip) &
                                  ], [2, 2] )
        !
        Bgp = mu0_mag * matmul(tMUR(1:ndime, 1:ndime), Hgp)
        !
        Hgp_mag(ielem) % a(1:ndime, igaus, 1) = Hgp
        Jgp_mag(ielem) % a(1, igaus, 1) = Jzgp
        Bgp_mag(ielem) % a(1:ndime, igaus, 1) = Bgp
        Fgp_mag(ielem) % a(1:2, igaus, 1) = Jzgp * [ -Bgp(2), Bgp(1) ]
        Ggp_mag(ielem) % a(1, igaus, 1) = mag_resist(pmate, nodq(1:ndime, igaus), Jgp, Hgp, 0.d0, .false., 3_ip) * Jzgp * Jzgp
        rc = rc + nodq(1, igaus)
        !
      end do

    case (12_ip)
      !
      ! QUA04
      !
!      tMUR(1:2, 1:2) = reshape( [ murmat_mag(1, pmate), 0.0_rp, 0.0_rp, murmat_mag(2, pmate) ], [2, 2] )
      !
      do igaus = 1, elegau
        !
        call mag_nedqua(nodq(1:ndime, igaus), nodi, basis, curl2, detJ)
        !
        Hgp = mag_matvec(basis(1:ndime, 1:eledof), He_mag(edgind_mag(1:eledof)) * sigval_mag(1:eledof))
        Jzgp = dot_product(curl2(1:eledof), He_mag(edgind_mag(1:eledof)) * sigval_mag(1:eledof))
        !
        Jgp = [0.0_rp, 0.0_rp, Jzgp]
        !
        tMUR(1:2, 1:2) = reshape( [ mag_permea(pmate, nodq(1:ndime, igaus), Hgp, .false., 1_ip), 0.0_rp, &
                                    0.0_rp, mag_permea(pmate, nodq(1:ndime, igaus), Hgp, .false., 2_ip) &
                                  ], [2, 2] )
        !
        Bgp = mu0_mag * matmul(tMUR(1:ndime, 1:ndime), Hgp)
        !
        Hgp_mag(ielem) % a(1:ndime, igaus, 1) = Hgp
        Jgp_mag(ielem) % a(1, igaus, 1) = Jzgp
        Bgp_mag(ielem) % a(1:ndime, igaus, 1) = Bgp
        Fgp_mag(ielem) % a(1:2, igaus, 1) = Jzgp * [ -Bgp(2), Bgp(1) ]
        Ggp_mag(ielem) % a(1, igaus, 1) = mag_resist(pmate, nodq(1:ndime, igaus), Jgp, Hgp, 0.d0, .false., 3_ip) * Jzgp * Jzgp
        rc = rc + nodq(1, igaus)
        !
      end do

    case (30_ip)
      !
      ! TET04
      !
!      tMUR(1:ndime, 1:ndime) = reshape( [ murmat_mag(1, pmate), 0.0_rp, 0.0_rp, &
!                                          0.0_rp, murmat_mag(2, pmate), 0.0_rp, &
!                                          0.0_rp, 0.0_rp, murmat_mag(3, pmate) ], &
!                                        [ndime, ndime] &
!                                      )
      !
      do igaus = 1, elegau
        !
        call mag_nedtet(nodq(1:ndime, igaus), nodi, basis, curl3, detJ)
        !
        Hgp = mag_matvec(basis(1:ndime, 1:eledof), He_mag(edgind_mag(1:eledof)) * sigval_mag(1:eledof))
        Jgp = mag_matvec(curl3(1:ndime, 1:eledof), He_mag(edgind_mag(1:eledof)) * sigval_mag(1:eledof))
        !
        tMUR(1:ndime, 1:ndime) = reshape( [ mag_permea(pmate, nodq(1:ndime, igaus), Hgp, .false., 1_ip), 0.0_rp, 0.0_rp, &
                                            0.0_rp, mag_permea(pmate, nodq(1:ndime, igaus), Hgp, .false., 2_ip), 0.0_rp, &
                                            0.0_rp, 0.0_rp, mag_permea(pmate, nodq(1:ndime, igaus), Hgp, .false., 3_ip) &
                                          ], [ndime, ndime] &
                                      )
        !
        Bgp = mu0_mag * matmul(tMUR(1:ndime, 1:ndime), Hgp)
        !
        Hgp_mag(ielem) % a(1:ndime, igaus, 1) = Hgp
        Jgp_mag(ielem) % a(1:ndime, igaus, 1) = Jgp
        Bgp_mag(ielem) % a(1:ndime, igaus, 1) = Bgp
        Fgp_mag(ielem) % a(1:ndime, igaus, 1) = mag_vecvec(Jgp, Bgp)
        !
        tRES = reshape( [ mag_resist(pmate, nodq(1:ndime, igaus), Jgp, Hgp, 0.d0, .false., 1_ip), 0.0_rp, 0.0_rp, &
                0.0_rp, mag_resist(pmate, nodq(1:ndime, igaus), Jgp, Hgp, 0.d0, .false., 2_ip), 0.0_rp, &
                0.0_rp, 0.0_rp, mag_resist(pmate, nodq(1:ndime, igaus), Jgp, Hgp, 0.d0, .false., 3_ip) ], &
                [3, 3] &
                )
        !
        Ggp_mag(ielem) % a(1, igaus, 1) = dot_product(Jgp, matmul(tRES(1:ndime, 1:ndime), Jgp))
        !
      end do
      !
    case (37_ip)
      !
      ! HEX08
      !
!      tMUR(1:ndime, 1:ndime) = reshape( [ murmat_mag(1, pmate), 0.0_rp, 0.0_rp, &
!                                          0.0_rp, murmat_mag(2, pmate), 0.0_rp, &
!                                          0.0_rp, 0.0_rp, murmat_mag(3, pmate) ], &
!                                        [ndime, ndime] &
!                                      )
      !
      do igaus = 1, elegau
        !
        call mag_nedhex(nodq(1:ndime, igaus), nodi, basis, curl3, detJ)
        !
        Hgp = mag_matvec(basis(1:ndime, 1:eledof), He_mag(edgind_mag(1:eledof)) * sigval_mag(1:eledof))
        Jgp = mag_matvec(curl3(1:ndime, 1:eledof), He_mag(edgind_mag(1:eledof)) * sigval_mag(1:eledof))
        !
        tMUR(1:ndime, 1:ndime) = reshape( [ mag_permea(pmate, nodq(1:ndime, igaus), Hgp, .false., 1_ip), 0.0_rp, 0.0_rp, &
                                            0.0_rp, mag_permea(pmate, nodq(1:ndime, igaus), Hgp, .false., 2_ip), 0.0_rp, &
                                            0.0_rp, 0.0_rp, mag_permea(pmate, nodq(1:ndime, igaus), Hgp, .false., 3_ip) &
                                          ], [ndime, ndime] &
                                      )
        !
        Bgp = mu0_mag * matmul(tMUR(1:ndime, 1:ndime), Hgp)
        !
        Hgp_mag(ielem) % a(1:ndime, igaus, 1) = Hgp
        Jgp_mag(ielem) % a(1:ndime, igaus, 1) = Jgp
        Bgp_mag(ielem) % a(1:ndime, igaus, 1) = Bgp
        Fgp_mag(ielem) % a(1:ndime, igaus, 1) = mag_vecvec(Jgp, Bgp)
        !
        tRES = reshape( [ mag_resist(pmate, nodq(1:ndime, igaus), Jgp, Hgp, 0.d0, .false., 1_ip), 0.0_rp, 0.0_rp, &
                0.0_rp, mag_resist(pmate, nodq(1:ndime, igaus), Jgp, Hgp, 0.d0, .false., 2_ip), 0.0_rp, &
                0.0_rp, 0.0_rp, mag_resist(pmate, nodq(1:ndime, igaus), Jgp, Hgp, 0.d0, .false., 3_ip) ], &
                [3, 3] &
                )
        !
        Ggp_mag(ielem) % a(1, igaus, 1) = dot_product(Jgp, matmul(tRES(1:ndime, 1:ndime), Jgp))
        !
      end do
      !
    end select
    !
  end do
  !
end subroutine mag_gpedge
