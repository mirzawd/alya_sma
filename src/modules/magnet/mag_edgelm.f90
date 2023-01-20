!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine mag_edgelm()
  !-----------------------------------------------------------------------
  !****f*  magnet/mag_edgelm
  ! NAME
  !    mag_edgelm
  ! DESCRIPTION
  !    This routine computes the solution field components by averaging over the elements
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
    elegau,    &
    pelty,    &
    pmate,    &
    qgauss

  real(rp) :: &
    nodi(1:ndime, 1:maxdof_mag),    &
    basis(1:ndime, 1:maxdof_mag),    &
    curl2(1:maxdof_mag),    &
    curl3(1:3, 1:maxdof_mag),    &
    nodq(ndime, maxgau_mag),    &
    detJ,    &
    gpvol,    &
    tMUR(3,3),    &
    tRES(3,3),    &
    Hgp(1:ndime),    &
    Jgp(3),    &
    Bgp(1:ndime),    &
    Fgp(1:ndime),    &
    Jzgp,    &
    Ggp,    &
    res,    &
    rc

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
    Hc_mag(1:ndime, ielem) = 0.0_rp
    Jc_mag(1:ndime, ielem) = 0.0_rp
    Bc_mag(1:ndime, ielem) = 0.0_rp
    Fc_mag(1:ndime, ielem) = 0.0_rp
    Jcz_mag(ielem) = 0.0_rp
    Gc_mag(ielem) = 0.0_rp
    gpvol = 0.0_rp
    rc = 0.0_rp
    !
    ! Basis functions
    !
    select case (pelty)
    case (10_ip)
      !
      ! TRI03
      !
      elegau = size(quadTri % wq)
      !
      ! Coordinates of gauss nodes in the reference element
      !
      nodq(1:ndime, 1:elegau) = matmul(nodi(1:ndime, 1:eledof), mag_lagref(quadTri % xq, pelty))
      !
      do qgauss = 1_ip, elegau
        !
        call mag_nedtri2(nodq(1:ndime, qgauss), nodi, basis, curl2, detJ)
        !
        Hgp(1:ndime) = mag_matvec(basis(1:ndime, 1:eledof), He_mag(edgind_mag(1:eledof)) * sigval_mag(1:eledof))
        !
        Jzgp = dot_product(curl2(1:eledof), He_mag(edgind_mag(1:eledof)) * sigval_mag(1:eledof))
        Jgp = [0.0_rp, 0.0_rp, Jzgp]
        !
        res = mag_resist(pmate, nodq(1:ndime, qgauss), Jgp, Hgp, 0.d0, .false., 3_ip)
        !
        tMUR(1:2, 1:2) = reshape( [ mag_permea(pmate, nodq(1:ndime, qgauss), Hgp, .false., 1_ip), 0.0_rp, &
                                            0.0_rp, mag_permea(pmate, nodq(1:ndime, qgauss), Hgp, .false., 2_ip) &
                                          ], [ 2_ip, 2_ip ] )
        !
        Bgp(1:ndime) = mu0_mag * mag_matvec(tMUR(1:ndime, 1:ndime), Hgp(1:ndime))
        !
        Fgp(1:2) = Jzgp * [ -Bgp(2), Bgp(1) ]
        !
        Hc_mag(1:ndime, ielem) = Hc_mag(1:ndime, ielem) + quadTri % wq(qgauss) * abs(detJ) * Hgp(1:ndime)
        !
        Jcz_mag(ielem) = Jcz_mag(ielem) + quadTri % wq(qgauss) * abs(detJ) * Jzgp
        !
        Bc_mag(1:ndime, ielem) = Bc_mag(1:ndime, ielem) + quadTri % wq(qgauss) * abs(detJ) * Bgp(1:ndime)
        !
        Fc_mag(1:ndime, ielem) = Fc_mag(1:ndime, ielem) + quadTri % wq(qgauss) * abs(detJ) * Fgp(1:ndime)
        !
        Gc_mag(ielem) = Gc_mag(ielem) + quadTri % wq(qgauss) * abs(detJ) * res * Jzgp * Jzgp
        !
        gpvol = gpvol + quadTri % wq(qgauss) * abs(detJ)
        !
        rc = rc + nodq(1, qgauss)
        !
      end do
      !
!      tMUR(1:ndime, 1:ndime) = reshape( [ murmat_mag(1, pmate), 0.0_rp, 0.0_rp, &
!                                          0.0_rp, murmat_mag(2, pmate), 0.0_rp, &
!                                          0.0_rp, 0.0_rp, murmat_mag(3, pmate) ], &
!                                        [3, 3] &
!                                      )
      !
    case (12_ip)
      !
      ! QUA04
      !
      elegau = size(quadQua % wq)
      !
      ! Coordinates of gauss nodes in the reference element
      !
      nodq(1:ndime, 1:elegau) = matmul(nodi(1:ndime, 1:eledof), mag_lagref(quadQua % xq, pelty))
      !
      do qgauss = 1_ip, elegau
        !
        call mag_nedqua(nodq(1:ndime, qgauss), nodi, basis, curl2, detJ)
        !
        Hgp(1:ndime) = mag_matvec(basis(1:ndime, 1:eledof), He_mag(edgind_mag(1:eledof)) * sigval_mag(1:eledof))
        !
        Jzgp = dot_product(curl2(1:eledof), He_mag(edgind_mag(1:eledof)) * sigval_mag(1:eledof))
        Jgp = [0.0_rp, 0.0_rp, Jzgp]
        !
        res = mag_resist(pmate, nodq(1:ndime, qgauss), Jgp, Hgp, 0.d0, .false., 3_ip)
        !
        tMUR(1:2, 1:2) = reshape( [ mag_permea(pmate, nodq(1:ndime, qgauss), Hgp, .false., 1_ip), 0.0_rp, &
                                            0.0_rp, mag_permea(pmate, nodq(1:ndime, qgauss), Hgp, .false., 2_ip) &
                                          ], [ 2_ip, 2_ip ] )
        !
        Bgp(1:ndime) = mu0_mag * mag_matvec(tMUR(1:ndime, 1:ndime), Hgp(1:ndime))
        !
        Fgp(1:2) = Jzgp * [ -Bgp(2), Bgp(1) ]
        !
        Hc_mag(1:ndime, ielem) = Hc_mag(1:ndime, ielem) + quadQua % wq(qgauss) * abs(detJ) * Hgp(1:ndime)
        !
        Jcz_mag(ielem) = Jcz_mag(ielem) + quadQua % wq(qgauss) * abs(detJ) * Jzgp
        !
        Bc_mag(1:ndime, ielem) = Bc_mag(1:ndime, ielem) + quadQua % wq(qgauss) * abs(detJ) * Bgp(1:ndime)
        !
        Fc_mag(1:ndime, ielem) = Fc_mag(1:ndime, ielem) + quadQua % wq(qgauss) * abs(detJ) * Fgp(1:ndime)
        !
        Gc_mag(ielem) = Gc_mag(ielem) + quadQua % wq(qgauss) * abs(detJ) * res * Jzgp * Jzgp
        !
        gpvol = gpvol + quadQua % wq(qgauss) * abs(detJ)
        !
        rc = rc + nodq(1, qgauss)
        !
      end do
      !
!      tMUR(1:ndime, 1:ndime) = reshape( [ murmat_mag(1, pmate), 0.0_rp, 0.0_rp, &
!                                          0.0_rp, murmat_mag(2, pmate), 0.0_rp, &
!                                          0.0_rp, 0.0_rp, murmat_mag(3, pmate) ], &
!                                        [3, 3] &
!                                      )
      !
    case (30_ip)
      !
      ! TET04
      !
      elegau = size(quadTet % wq)
      !
      ! Coordinates of gauss nodes in the reference element
      !
      nodq(1:ndime, 1:elegau) = matmul(nodi(1:ndime, 1:eledof), mag_lagref(quadTet % xq, pelty))
      !
      do qgauss = 1, elegau
        !
        call mag_nedtet(nodq(1:ndime, qgauss), nodi, basis, curl3, detJ)
        !
        Hgp(1:ndime) = mag_matvec(basis(1:ndime, 1:eledof), He_mag(edgind_mag(1:eledof)) * sigval_mag(1:eledof))
        !
        Jgp(1:ndime) = mag_matvec(curl3(1:ndime, 1:eledof), He_mag(edgind_mag(1:eledof)) * sigval_mag(1:eledof))
        !
        tRES = reshape( [ mag_resist(pmate, nodq(1:ndime, qgauss), Jgp, Hgp, 0.d0, .false., 1_ip), 0.0_rp, 0.0_rp, &
                          0.0_rp, mag_resist(pmate, nodq(1:ndime, qgauss), Jgp, Hgp, 0.d0, .false., 2_ip), 0.0_rp, &
                          0.0_rp, 0.0_rp, mag_resist(pmate, nodq(1:ndime, qgauss), Jgp, Hgp, 0.d0, .false., 3_ip) &
                        ], [3, 3] )
        !      
        tMUR = reshape( [ mag_permea(pmate, nodq(1:ndime, qgauss), Hgp, .false., 1_ip), 0.0_rp, 0.0_rp, &
                          0.0_rp, mag_permea(pmate, nodq(1:ndime, qgauss), Hgp, .false., 2_ip), 0.0_rp, &
                          0.0_rp, 0.0_rp, mag_permea(pmate, nodq(1:ndime, qgauss), Hgp, .false., 3_ip) &
                        ], [3, 3] )
        !
        Bgp = mu0_mag * mag_matvec(tMUR, Hgp)
        !
        Ggp = dot_product(Jgp, mag_matvec(tRES, Jgp)) 
        !
        Fgp(1:ndime) = mag_vecvec(Jgp(1:ndime), Bgp(1:ndime))
        !
        Hc_mag(1:ndime, ielem) = Hc_mag(1:ndime, ielem) + quadTet % wq(qgauss) * abs(detJ) * Hgp(1:ndime)
        !
        Jc_mag(1:ndime, ielem) = Jc_mag(1:ndime, ielem) + quadTet % wq(qgauss) * abs(detJ) * Jgp(1:ndime)
        !
        Bc_mag(1:ndime, ielem) = Bc_mag(1:ndime, ielem) + quadTet % wq(qgauss) * abs(detJ) * Bgp(1:ndime)
        !
        Fc_mag(1:ndime, ielem) = Fc_mag(1:ndime, ielem) + quadTet % wq(qgauss) * abs(detJ) * Fgp(1:ndime)
        !
        Gc_mag(ielem) = Gc_mag(ielem) + quadTet % wq(qgauss) * abs(detJ) * Ggp
        !
        gpvol = gpvol + quadTet % wq(qgauss) * abs(detJ)
        !
      end do
      !
!      tMUR(1:ndime, 1:ndime) = reshape( [ murmat_mag(1, pmate), 0.0_rp, 0.0_rp, &
!                                          0.0_rp, murmat_mag(2, pmate), 0.0_rp, &
!                                          0.0_rp, 0.0_rp, murmat_mag(3, pmate) ], &
!                                        [3, 3] &
!                                      )
      !
    case (37_ip)
      !
      ! HEX08
      !
      elegau = size(quadHex % wq)
      !
      ! Coordinates of gauss nodes in the reference element
      !
      nodq(1:ndime, 1:elegau) = matmul(nodi(1:ndime, 1:eledof), mag_lagref(quadHex % xq, pelty))
      !
      do qgauss = 1, elegau
        !
        call mag_nedhex(nodq(1:ndime, qgauss), nodi, basis, curl3, detJ)
        !
        Hgp(1:ndime) = mag_matvec(basis(1:ndime, 1:eledof), He_mag(edgind_mag(1:eledof)) * sigval_mag(1:eledof))
        !
        Jgp(1:ndime) = mag_matvec(curl3(1:ndime, 1:eledof), He_mag(edgind_mag(1:eledof)) * sigval_mag(1:eledof))
        !
        tRES = reshape( [ mag_resist(pmate, nodq(1:ndime, qgauss), Jgp, Hgp, 0.d0, .false., 1_ip), 0.0_rp, 0.0_rp, &
                          0.0_rp, mag_resist(pmate, nodq(1:ndime, qgauss), Jgp, Hgp, 0.d0, .false., 2_ip), 0.0_rp, &
                          0.0_rp, 0.0_rp, mag_resist(pmate, nodq(1:ndime, qgauss), Jgp, Hgp, 0.d0, .false., 3_ip) &
                        ], [3, 3] )
        !
        tMUR = reshape( [ mag_permea(pmate, nodq(1:ndime, qgauss), Hgp, .false., 1_ip), 0.0_rp, 0.0_rp, &
                          0.0_rp, mag_permea(pmate, nodq(1:ndime, qgauss), Hgp, .false., 2_ip), 0.0_rp, &
                          0.0_rp, 0.0_rp, mag_permea(pmate, nodq(1:ndime, qgauss), Hgp, .false., 3_ip) &
                          ], [3, 3] )
        !
        Bgp(1:ndime) = mu0_mag * mag_matvec(tMUR(1:ndime, 1:ndime), Hgp(1:ndime))
        !
        Fgp(1:ndime) = mag_vecvec(Jgp(1:ndime), Bgp(1:ndime))
        !
        Ggp = dot_product(Jgp, mag_matvec(tRES, Jgp))
        !
        Hc_mag(1:ndime, ielem) = Hc_mag(1:ndime, ielem) + quadHex % wq(qgauss) * abs(detJ) * Hgp(1:ndime)
        !
        Jc_mag(1:ndime, ielem) = Jc_mag(1:ndime, ielem) + quadHex % wq(qgauss) * abs(detJ) * Jgp(1:ndime)
        !
        Fc_mag(1:ndime, ielem) = Fc_mag(1:ndime, ielem) + quadHex % wq(qgauss) * abs(detJ) * Fgp(1:ndime)
        !
        Gc_mag(ielem) = Gc_mag(ielem) + quadHex % wq(qgauss) * abs(detJ) * Ggp
        !
        gpvol = gpvol + quadHex % wq(qgauss) * abs(detJ)
        !
      end do
      !
!      tMUR(1:ndime, 1:ndime) = reshape( [ murmat_mag(1, pmate), 0.0_rp, 0.0_rp, &
!                                          0.0_rp, murmat_mag(2, pmate), 0.0_rp, &
!                                          0.0_rp, 0.0_rp, murmat_mag(3, pmate) ], &
!                                        [3, 3] &
!                                      )
      !
    end select
    !
    Hc_mag(1:ndime, ielem) = Hc_mag(1:ndime, ielem) / gpvol
    Jc_mag(1:ndime, ielem) = Jc_mag(1:ndime, ielem) / gpvol
    Bc_mag(1:ndime, ielem) = Bc_mag(1:ndime, ielem) / gpvol
    Fc_mag(1:ndime, ielem) = Fc_mag(1:ndime, ielem) / gpvol
    Jcz_mag(ielem) = Jcz_mag(ielem) / gpvol
    Gc_mag(ielem) = Gc_mag(ielem) / gpvol
!    !
!    Bc_mag(1:ndime, ielem) = mu0_mag * matmul(tMUR(1:ndime, 1:ndime), Hc_mag(1:ndime, ielem))
!    Fc_mag(1:ndime, ielem) = mu0_mag * matmul(tMUR(1:ndime, 1:ndime), Fc_mag(1:ndime, ielem))
    !
  end do
  !
end subroutine mag_edgelm
