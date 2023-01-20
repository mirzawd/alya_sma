!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine mag_elmope()
  !-----------------------------------------------------------------------
  !****f* magnet/mag_elmope.f90
  ! NAME
  !    mag_elmope
  ! DESCRIPTION
  !    This routine
  !    1. Computes the elementary matrix and RHS for every element in the mesh
  !    2. Assembles the elementary equations into the system equations
  ! USES
  ! USED BY
  !    mag_matrix
  !***
  !-----------------------------------------------------------------------
  use def_master, only: solve_sol, amatr, rhsid
  use def_domain, only: meshe, ndime
!  use def_master, only: solve
  use def_magnet
  use def_kermod, only: ndivi
  use mod_mag_matpro
  use mod_mag_linalg, only: mag_matvec, mag_vecvec, mag_vecmod
  use mod_mag_nedele
  use mod_mag_lagran
  use mod_mag_inpdat, only: mag_source

  implicit none

  integer(ip) :: &
       ielem,    &
       iedge,    &
       jedge,    &
       inode,    &
       qgauss,   &
       k1, k2, j1

  real(rp) :: &
       elevol,    &
       currDens,    &
       piv,    &
       Hq(1:ndime),    &
       Jq(3),    &
       resq,    &
       res,    &
       resDiff

  integer(ip), pointer :: &
       r_edg(:),    &
       c_edg(:)

  integer(ip) :: &
       eledof,    &
       elegau,    &
       pelty,    &
       pmate,    &
       noddof

  real(rp) :: &
!       elstif(maxdof_mag, maxdof_mag),    &
!       elstdf(maxdof_mag, maxdof_mag),    &
!       elmass(maxdof_mag, maxdof_mag),    &
!       Aelem(maxdof_mag, maxdof_mag),    &
!       Anelem(maxdof_mag, maxdof_mag),    &
!       belem(maxdof_mag),    &
!       rnelem(maxdof_mag),    &
!       elsrc(maxdof_mag),    &
       nodi(ndime, maxdof_mag),    &
       nodq(ndime, maxgau_mag),    &
       source(ndime),    &
       bq(ndime, maxdof_mag),    &
       pq(mxndof_mag, maxgau_mag),    &
       curl2(maxdof_mag),    &
       curl3(ndime, maxdof_mag),    &
       detJ
  !
  ! Material Properties Tensors
  !
  real(rp) :: &
       tMUR(3,3),    &
       tRES(3,3),    &
       tRDF(3,3)

  logical :: flag
  !
  ! Magnetic Energy, Joule dissipation
  !
  magnen_mag = 0.0_rp; joulen_mag = 0.0_rp
  !
  ! Magnetization
  !
  magtiz_mag = 0.0_rp
  cursum_mag = 0.0_rp
  magsum_mag = 0.0_rp
  volume_mag = 0.0_rp
  !
  ! Crank-Nicholson
  !
  Ha_mag = theta_mag * He_mag + (1-theta_mag) * Hp_mag
  !
  ! Initialize matrix and RHS vector before assembly
  !
  amatr = 0.0_rp
  rhsid = 0.0_rp
  bhsid = 0.0_rp
  fhsid = 0.0_rp
  hhsid = 0.0_rp
  !
  ! CSR indices
  !
  r_edg => solve_sol(1) % ia
  c_edg => solve_sol(1) % ja
  ! 
  ! Start building elementary matrices
  !
  do ielem = 1_ip, meshe(ndivi) % nelem
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
     pmate = meshe(ndivi) % lmate(ielem)
     !
     ! Number of degrees of freedom in current element
     !
     eledof = mag_edgdof(pelty)
     !
     noddof = mag_noddof(pelty)
     !
     ! Element area
     !
     elevol = elevol_mag(ielem)
     if (elevol <= 0) call runend("mag_elmope: negative volume")
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
     !--------------------------------------------------------------------------------
     !
     ! Initialize local matrices: mass, stiffness, stiffness derivative
     !
     elmass(1:eledof, 1:eledof) = 0.0_rp
     elstif(1:eledof, 1:eledof) = 0.0_rp
     elstdf(1:eledof, 1:eledof) = 0.0_rp
     elsrc(1:eledof) = 0.0_rp
     elheat(1:noddof) = 0.0_rp
     !
     if (pelty == 10_ip) then
        !
        ! TRI03
        !
        elegau = size(quadTri % wq)
        !
        ! Coordinates of element vertices
        !
        nodi(1:ndime, 1:eledof) = meshe(ndivi) % coord(1:ndime, nodind_mag(1:eledof))
        !
        ! Coordinates of gauss nodes in the reference element
        !
        pq(1:noddof, 1:elegau) = mag_lagref(quadTri % xq, pelty)
        !
        nodq(1:ndime, 1:elegau) = matmul(nodi(1:ndime, 1:eledof), mag_lagref(quadTri % xq, pelty))
        !
        if (struct_mag) then
           ! Elementary Current Density
           !
           currDens = dot_product(lenval_mag(1:eledof), Ha_mag(edgind_mag(1:eledof))) / elevol
           Jq = [ 0.0_rp, 0.0_rp, currDens]
           !
           ! Gauss quadrature
           !
           do qgauss = 1_ip, elegau
              !
              ! Vector functions evaluated at current quadrature node
              !
              bq(1:2, 1:eledof) = mag_nedtri(nodq(1:ndime, qgauss), nodi(1:ndime, 1:eledof))
              !
              ! Magnetic field evaluated at current quadrature node
              ! 
              Hq = mag_matvec( bq(1:ndime, 1:eledof), Ha_mag(edgind_mag(1:eledof)) * sigval_mag(1:eledof) )
              !
              ! Resistivity and its derivative with respect to H
              !
              res = mag_resist(pmate, nodq(1:ndime, qgauss), Jq, Hq, 0.d0, kfl_axsym_mag, 3_ip)
              resDiff = mag_resdif(pmate, nodq(1:ndime, qgauss), Jq, Hq, 0.d0, kfl_axsym_mag, 3_ip)
              !
              tMUR(1:2, 1:2) = reshape( [ mag_permea(pmate, nodq(1:ndime, qgauss), Hq, kfl_axsym_mag, 1_ip), 0.0_rp, &
                   0.0_rp, mag_permea(pmate, nodq(1:ndime, qgauss), Hq, kfl_axsym_mag, 2_ip) ], &
                   [2, 2] &
                   )
              !
              ! Source term
              !
              source = mag_source(myClockTime_mag + theta_mag * dt_mag, nodq(1:ndime, qgauss))
              !
              ! Elementary matrices at current quadrature node
              !
              do iedge = 1, eledof
                 do jedge = 1, eledof
                    !
                    ! Mass matrix
                    !
                    elmass(iedge, jedge) = elmass(iedge, jedge) + quadTri % wq(qgauss) * &
                         dot_product(bq(1:ndime, iedge), matmul(tMUR(1:ndime, 1:ndime), bq(1:ndime, jedge))) * &
                         sigval_mag(iedge) * sigval_mag(jedge)
                    !
                    ! Stiffness matrix
                    !
                    elstif(iedge, jedge) = elstif(iedge, jedge) + quadTri % wq(qgauss) * &
                         res * lenval_mag(iedge) * lenval_mag(jedge)
                    !
                    ! Stiffness derivative matrix
                    !
                    elstdf(iedge, jedge) = elstdf(iedge, jedge) + quadTri % wq(qgauss) * &
                         resDiff * lenval_mag(iedge) * lenval_mag(jedge)
                 end do
                 !
                 ! Source vector
                 !
                 elsrc(iedge) = elsrc(iedge) + quadTri % wq(qgauss) * dot_product(source(1:ndime), bq(1:ndime, iedge)) * sigval_mag(iedge)
                 !  
              end do
              !
              do inode = 1, noddof
                !
                elheat(inode) = elheat(inode) + quadTri % wq(qgauss) * res * currDens * currDens * pq(inode, qgauss)
                !
              end do
              !
           end do
           !
           ! Elementary Mass Matrix
           !
           elmass(1:eledof, 1:eledof) = elmass(1:eledof, 1:eledof) * 2.0_rp * elevol * mu0_mag !* mur_mag(pmate)
           !
           ! Elementary Stiffness Matrix
           !
           elstif(1:eledof, 1:eledof) = elstif(1:eledof, 1:eledof) * 2.0_rp / elevol
           !
           ! Elementary Stiffness Derivative Matrix
           !
           elstdf(1:eledof, 1:eledof) = elstdf(1:eledof, 1:eledof) * abs(currDens) * 2.0_rp / elevol
           !
           ! Elementary Source Vector
           !
           elsrc(1:eledof) = elsrc(1:eledof) * 2.0_rp * elevol * mu0_mag * mur_mag(pmate)
           !
           ! Elementary Heat Vector
           !
           elheat(1:noddof) = elheat(1:noddof) * 2.0_rp * elevol
           !
        else
           !
           do qgauss = 1, elegau
              !
              ! Vector basis functions at quadrature node and their curl
              !
              call mag_nedtri2(quadTri % xq(1:ndime, qgauss), nodi(1:ndime, 1:eledof), bq, curl2, detJ)
              !
              ! Magnetic field at quadrature node
              !
              Hq = mag_matvec( bq(1:ndime, 1:eledof), Ha_mag(edgind_mag(1:eledof)) * sigval_mag(1:eledof))
              !
              ! Current density at quadrature node
              !
              currDens = dot_product(curl2(1:eledof), Ha_mag(edgind_mag(1:eledof)) * sigval_mag(1:eledof))
              Jq = [ 0.0_rp, 0.0_rp, currDens]
              !
              magtiz_mag(1:2, pmate) = magtiz_mag(1:2, pmate) + 0.5_rp * quadTri % wq(qgauss) * abs(detJ) * &
                currDens * [nodq(2, qgauss) - momori_mag(2, pmate), - nodq(1, qgauss) + momori_mag(1, pmate)]
              !
              cursum_mag(3, pmate) = cursum_mag(3, pmate) + quadTri % wq(qgauss) * abs(detJ) * currDens
              !
              magsum_mag(1:2, pmate) = magsum_mag(1:2, pmate) + quadTri % wq(qgauss) * abs(detJ) * Hq(1:2)
              ! 
              volume_mag(pmate) = volume_mag(pmate) + quadTri % wq(qgauss) * abs(detJ)
              !
              ! Resistivity and its derivative with respect to H
              !
              res = mag_resist(pmate, nodq(1:ndime, qgauss), Jq, Hq, 0.d0, kfl_axsym_mag, 3_ip)
              resDiff = mag_resdif(pmate, nodq(1:ndime, qgauss), Jq, Hq, 0.d0, kfl_axsym_mag, 3_ip)
              !
              tMUR(1:2, 1:2) = reshape( [ mag_permea(pmate, nodq(1:ndime, qgauss), Hq, kfl_axsym_mag, 1_ip), 0.0_rp, &
                   0.0_rp, mag_permea(pmate, nodq(1:ndime, qgauss), Hq, kfl_axsym_mag, 2_ip) ], &
                   [2, 2] &
                   )
              !
              ! Source term
              !
              source = mag_source(myClockTime_mag + theta_mag * dt_mag, nodq(1:ndime, qgauss))
              !
              ! Elementary matrices at current quadrature node
              !
              do iedge = 1, eledof
                 do jedge = 1, eledof
                    !
                    ! Mass matrix
                    !
                    elmass(iedge, jedge) = elmass(iedge, jedge) + quadTri % wq(qgauss) * mu0_mag * & !mur_mag(pmate)
                         dot_product(bq(1:ndime, iedge), matmul(tMUR(1:ndime, 1:ndime), bq(1:ndime, jedge))) * &
                         sigval_mag(iedge) * sigval_mag(jedge) * abs(detJ)
                    !
                    ! Stiffness matrix
                    !
                    elstif(iedge, jedge) = elstif(iedge, jedge) + quadTri % wq(qgauss) * &
                         res * curl2(iedge) * curl2(jedge) * sigval_mag(iedge) * sigval_mag(jedge) * abs(detJ)
                    !
                    ! Stiffness derivative matrix
                    !
                    elstdf(iedge, jedge) = elstdf(iedge, jedge) + quadTri % wq(qgauss) * abs(currDens) * &
                         resDiff * curl2(iedge) * curl2(jedge) * sigval_mag(iedge) * sigval_mag(jedge) * abs(detJ)
                    !
                 end do
                 !
                 ! Source vector
                 !
                 elsrc(iedge) = elsrc(iedge) + quadTri % wq(qgauss) * dot_product(source(1:ndime), bq(1:ndime, iedge)) * &
                      sigval_mag(iedge) * abs(detJ) * mu0_mag * mur_mag(pmate)
                 !
              end do
              !
              do inode = 1, noddof
                !
                elheat(inode) = elheat(inode) + quadTri % wq(qgauss) * abs(detJ) * res * currDens * currDens * pq(inode, qgauss)
                !
              end do
              !
           end do
           !
        end if
        !
     elseif (pelty == 12_ip) then
        !
        ! QUA04
        !
        elegau = size(quadQua % wq)
        !
        ! Coordinates of element vertices
        !
        nodi(1:ndime, 1:eledof) = meshe(ndivi) % coord(1:ndime, nodind_mag(1:eledof))
        !
        ! Coordinates of Gauss quadrature nodes
        !
        pq(1:noddof, 1:elegau) = mag_lagref(quadQua % xq, pelty)
        !
        nodq(1:ndime, 1:elegau) = matmul(nodi(1:ndime, 1:eledof), mag_lagref(quadQua % xq, pelty))
        !
        if (struct_mag) then
           !
           ! Elementary Current Density
           !
           currDens = dot_product(lenval_mag(1:eledof), Ha_mag(edgind_mag(1:eledof))) / elevol
           Jq = [0.0_rp, 0.0_rp, currDens]
           !
           res = 0.0_rp; resDiff = 0.0_rp
           !
           do qgauss = 1, elegau
              !
              ! Vector basis functions
              !
              bq(1:2, 1:eledof) = mag_nedrec(nodq(1:ndime, qgauss), nodi)
              !
              ! Magnetic field at quadrature node
              !
              Hq = mag_matvec( bq(1:ndime, 1:eledof), Ha_mag(edgind_mag(1:eledof)) * sigval_mag(1:eledof)) 
              !
              ! Resistivity and its derivative with respect to H
              !
              resq = quadQua % wq(qgauss) * mag_resist(pmate, nodq(1:ndime, qgauss), Jq, Hq, 0.d0, kfl_axsym_mag, 3_ip)
              res = res + resq ! quadQua % wq(qgauss) * mag_resist(pmate, nodq(1:ndime, qgauss), Jq, Hq, 0.d0, kfl_axsym_mag, 3_ip)
              resDiff = resDiff + quadQua % wq(qgauss) * mag_resdif(pmate, nodq(1:ndime, qgauss), Jq, Hq, 0.d0, kfl_axsym_mag, 3_ip)
              !
              !          tMUR(1:ndime, 1:ndime) = reshape( [ mag_permea(pmate, nodq(1:ndime, qgauss), Hq, kfl_axsym_mag, 1_ip), 0.0_rp, &
              !                                              0.0_rp, mag_permea(pmate, nodq(1:ndime, qgauss), Hq, kfl_axsym_mag, 2_ip) ], &
              !                                            [ndime, ndime] &
              !                                          )
              !
              ! Source term
              !
              source = mag_source(myClockTime_mag + theta_mag * dt_mag, nodq(1:ndime, qgauss))
              !
              do iedge = 1, eledof
                 !
                 ! Source vector
                 !
                 elsrc(iedge) = elsrc(iedge) + quadQua % wq(qgauss) * dot_product(source(1:ndime), bq(1:ndime, iedge)) * sigval_mag(iedge)
                 !
              end do
              !
              do inode = 1, noddof
                !
                elheat(inode) = elheat(inode) + quadQua % wq(qgauss) * resq * currDens * currDens * pq(inode, qgauss)
                !
              end do
              !
           end do
           !
           res = res / 4.0_rp; resDiff = resDiff / 4.0_rp
           !
           !
           ! Elementary Mass matrix
           !
           elmass(1:eledof, 1:eledof) = reshape( [ 2.0_rp, 0.0_rp, -1.0_rp, 0.0_rp, &
                0.0_rp, 2.0_rp, 0.0_rp, -1.0_rp, &
                -1.0_rp, 0.0_rp, 2.0_rp, 0.0_rp, &
                0.0_rp, -1.0_rp, 0.0_rp, 2.0_rp ], [ 4, 4 ] ) / 6.0_rp
           !
           do iedge = 1, eledof
              do jedge = 1, eledof
                 !
                 ! Mass matrix
                 !
                 elmass(iedge, jedge) = elmass(iedge, jedge) * sigval_mag(iedge) * sigval_mag(jedge)
                 !
                 ! Stiffness matrix
                 !
                 elstif(iedge, jedge) = res * lenval_mag(iedge) * lenval_mag(jedge)
                 !
                 ! Stiffness Derivative matrix
                 !
                 elstdf(iedge, jedge) = resDiff * lenval_mag(iedge) * lenval_mag(jedge)
                 !
              end do
           end do
           !
           ! Elementary Mass Matrix
           !
           elmass(1:eledof, 1:eledof) = elmass(1:eledof, 1:eledof) * elevol * mu0_mag * mur_mag(pmate)
           !
           ! Elementary Stiffness Matrix
           !
           elstif(1:eledof, 1:eledof) = elstif(1:eledof, 1:eledof) / elevol
           !
           ! Elementary Stiffness Matrix Derivative
           !
           elstdf(1:eledof, 1:eledof) = elstdf(1:eledof, 1:eledof) / elevol * abs(currDens)
           !
           ! Elementary Source Vector
           !
           elsrc(1:eledof) = elsrc(1:eledof) * elevol / 4.0_rp * mu0_mag * mur_mag(pmate)
           !
           ! Elementary Heat Vector
           !
           elheat(1:noddof) = elheat(1:noddof) * elevol / 4.0_rp
           !
        else
           !
           do qgauss = 1, elegau
              !                                          
              ! Vector basis functions at quadrature node and their curl
              !
              call mag_nedqua(quadQua % xq(1:ndime, qgauss), nodi, bq, curl2, detJ)
              !
              ! Magnetic field at quadrature node
              !
              Hq = mag_matvec(bq(1:ndime, 1:eledof), Ha_mag(edgind_mag(1:eledof)) * sigval_mag(1:eledof))
              !
              ! Current density at quadrature node
              !
              currDens = dot_product(curl2(1:eledof), Ha_mag(edgind_mag(1:eledof)) * sigval_mag(1:eledof))
              Jq = [ 0.0_rp, 0.0_rp, currDens]
              !
              magtiz_mag(1:2, pmate) = magtiz_mag(1:2, pmate) + 0.5_rp * quadQua % wq(qgauss) * abs(detJ) * &
                currDens * [nodq(2, qgauss) - momori_mag(2, pmate), - nodq(1, qgauss) + momori_mag(1, pmate)]
              !
              cursum_mag(3, pmate) = cursum_mag(3, pmate) + quadQua % wq(qgauss) * abs(detJ) * currDens
              !
              magsum_mag(1:2, pmate) = magsum_mag(1:2, pmate) + quadQua % wq(qgauss) * abs(detJ) * Hq(1:2)
              !
              volume_mag(pmate) = volume_mag(pmate) + quadQua % wq(qgauss) * abs(detJ)
              !
              ! Resistivity and its derivative with respect to H
              !
              res = mag_resist(pmate, nodq(1:ndime, qgauss), Jq, Hq, 0.d0, kfl_axsym_mag, 3_ip)
              resDiff = mag_resdif(pmate, nodq(1:ndime, qgauss), Jq, Hq, 0.d0, kfl_axsym_mag, 3_ip)
              !
              tMUR(1:2, 1:2) = reshape( [ mag_permea(pmate, nodq(1:ndime, qgauss), Hq, kfl_axsym_mag, 1_ip), 0.0_rp, &
                   0.0_rp, mag_permea(pmate, nodq(1:ndime, qgauss), Hq, kfl_axsym_mag, 2_ip) ], &
                   [2, 2] &
                   )
              !
              ! Source term
              !
              source = mag_source(myClockTime_mag + theta_mag * dt_mag, nodq(1:ndime, qgauss))
              !
              ! Elementary matrices at current quadrature node
              !
              do iedge = 1, eledof
                 do jedge = 1, eledof
                    !
                    ! Mass matrix
                    !
                    elmass(iedge, jedge) = elmass(iedge, jedge) + quadQua % wq(qgauss) * mu0_mag * &
                         dot_product(bq(1:ndime, iedge), matmul(tMUR(1:ndime, 1:ndime), bq(1:ndime, jedge))) * &
                         sigval_mag(iedge) * sigval_mag(jedge) * abs(detJ)
                    !
                    ! Stiffness matrix
                    !
                    elstif(iedge, jedge) = elstif(iedge, jedge) + quadQua % wq(qgauss) * &
                         res * curl2(iedge) * curl2(jedge) * sigval_mag(iedge) * sigval_mag(jedge) * abs(detJ)
                    !
                    ! Stiffness derivative matrix
                    !
                    elstdf(iedge, jedge) = elstdf(iedge, jedge) + quadQua % wq(qgauss) * abs(currDens) * &
                         resDiff * curl2(iedge) * curl2(jedge) * sigval_mag(iedge) * sigval_mag(jedge) * abs(detJ)
                    !
                 end do
                 !
                 ! Source vector
                 !
                 elsrc(iedge) = elsrc(iedge) + quadQua % wq(qgauss) * dot_product(source(1:ndime), bq(1:ndime, iedge)) * &
                      sigval_mag(iedge) * abs(detJ) * mu0_mag * mur_mag(pmate)
                 !
              end do
              !
              do inode = 1, noddof
                !
                elheat(inode) = elheat(inode) + quadQua % wq(qgauss) * abs(detJ) * res * currDens * currDens * pq(inode, qgauss)
                !
              end do
              !
           end do
           !
        end if
        !
     elseif (pelty == 30_ip) then
        !
        ! TET04
        !
        elegau = size(quadTet % wq)
        !
        ! Coordinates of element vertices
        !
        nodi(1:ndime, 1:eledof) = meshe(ndivi) % coord(1:ndime, nodind_mag(1:eledof))
        !
        ! Coordinates of gauss nodes in the reference element
        !
        pq(1:noddof, 1:elegau) = mag_lagref(quadTet % xq, pelty)
        !
        nodq(1:ndime, 1:elegau) = matmul(nodi(1:ndime, 1:eledof), mag_lagref(quadTet % xq, pelty))
        !
        do qgauss = 1, elegau
           !
           ! Vector basis functions at quadrature node and their curl
           !
           call mag_nedtet(quadTet % xq(1:ndime, qgauss), nodi, bq, curl3, detJ)
           !
           ! Magnetic field at quadrature node
           !
           Hq = mag_matvec( bq(1:ndime, 1:eledof), Ha_mag(edgind_mag(1:eledof)) * sigval_mag(1:eledof))
           !
           ! Current density at quadrature node
           !
           Jq = mag_matvec(curl3(1:ndime, 1:eledof), Ha_mag(edgind_mag(1:eledof)) * sigval_mag(1:eledof))
           !
           ! Magnetization
           !
           magtiz_mag(1:ndime, pmate) = magtiz_mag(1:ndime, pmate) + 0.5_rp * quadTet % wq(qgauss) * abs(detJ) * &
                mag_vecvec(nodq(1:ndime, qgauss) - momori_mag(1:ndime, pmate), Jq)
           !
           cursum_mag(1:ndime, pmate) = cursum_mag(1:ndime, pmate) + quadTet % wq(qgauss) * abs(detJ) * Jq
           !
           magsum_mag(1:ndime, pmate) = magsum_mag(1:ndime, pmate) + quadTet % wq(qgauss) * abs(detJ) * Hq(1:ndime)
           !
           volume_mag(pmate) = volume_mag(pmate) + quadTet % wq(qgauss) * abs(detJ)
           !
           ! Resistivity and its derivative at quadrature node
           !
           tRES = reshape( [ mag_resist(pmate, nodq(1:ndime, qgauss), Jq, Hq, 0.d0, .false., 1_ip), 0.0_rp, 0.0_rp, &
                0.0_rp, mag_resist(pmate, nodq(1:ndime, qgauss), Jq, Hq, 0.d0, .false., 2_ip), 0.0_rp, &
                0.0_rp, 0.0_rp, mag_resist(pmate, nodq(1:ndime, qgauss), Jq, Hq, 0.d0, .false., 3_ip) ], &
                [3, 3] &
                )
           tRDF = reshape( [ mag_resdif(pmate, nodq(1:ndime, qgauss), Jq, Hq, 0.d0, .false., 1_ip), 0.0_rp, 0.0_rp, &
                0.0_rp, mag_resdif(pmate, nodq(1:ndime, qgauss), Jq, Hq, 0.d0, .false., 2_ip), 0.0_rp, &
                0.0_rp, 0.0_rp, mag_resdif(pmate, nodq(1:ndime, qgauss), Jq, Hq, 0.d0, .false., 3_ip) ], &
                [3, 3] &
                )
           tMUR = reshape( [ mag_permea(pmate, nodq(1:ndime, qgauss), Hq, .false., 1_ip), 0.0_rp, 0.0_rp, &
                0.0_rp, mag_permea(pmate, nodq(1:ndime, qgauss), Hq, .false., 2_ip), 0.0_rp, &
                0.0_rp, 0.0_rp, mag_permea(pmate, nodq(1:ndime, qgauss), Hq, .false., 3_ip) ], &
                [3, 3] &
                )
           !
           ! Source term
           !
           source = mag_source(myClockTime_mag + theta_mag * dt_mag, nodq(1:ndime, qgauss))
           !
           ! Elementary matrices at current quadrature node
           !
           do iedge = 1, eledof
              do jedge = 1, eledof
                 !
                 ! Mass matrix at quadrature node
                 !
                 !            elmass(iedge, jedge) = elmass(iedge, jedge) + quadTet % wq(qgauss) * mu0_mag * mur_mag(pmate) * &
                 !                dot_product(bq(1:ndime, iedge), bq(1:ndime, jedge)) * sigval_mag(iedge) * sigval_mag(jedge) * abs(detJ)

                 elmass(iedge, jedge) = elmass(iedge, jedge) + quadTet % wq(qgauss) * mu0_mag * &
                      dot_product(bq(1:ndime, iedge), matmul(tMUR(1:ndime, 1:ndime), bq(1:ndime, jedge))) * &
                      sigval_mag(iedge) * sigval_mag(jedge) * abs(detJ)
                 !
                 ! Stiffness matrix at quadrature node
                 !
                 !            elstif(iedge, jedge) = elstif(iedge, jedge) + quadTet % wq(qgauss) * &
                 !                res * dot_product(curl3(1:ndime, iedge), curl3(1:ndime, jedge)) * sigval_mag(iedge) * sigval_mag(jedge) * abs(detJ)

                 elstif(iedge, jedge) = elstif(iedge, jedge) + quadTet % wq(qgauss) * &
                      dot_product(curl3(1:ndime, iedge), matmul(tRES(1:ndime, 1:ndime), curl3(1:ndime, jedge))) * &
                      sigval_mag(iedge) * sigval_mag(jedge) * abs(detJ)
                 !
                 ! Stiffness derivative matrix
                 !
                 if (mag_vecmod(Jq) > SCL_MAG * ZER_MAG) then
                    !              elstdf(iedge, jedge) = elstdf(iedge, jedge) + quadTet % wq(qgauss) * resDiff * abs(detJ) / mag_vecmod(Jq) * &
                    !                   dot_product(Jq, curl3(1:ndime, iedge)) * sigval_mag(iedge) * &
                    !                   dot_product(Jq, curl3(1:ndime, jedge)) * sigval_mag(jedge)

                    elstdf(iedge, jedge) = elstdf(iedge, jedge) + quadTet % wq(qgauss) * abs(detJ) / mag_vecmod(Jq) * &
                         dot_product(matmul(tRDF(1:ndime, 1:ndime), Jq), curl3(1:ndime, iedge)) * sigval_mag(iedge) * &
                         dot_product(Jq, curl3(1:ndime, jedge)) * sigval_mag(jedge)
                 else
                    elstdf(iedge, jedge) = 0.0_rp
                 end if
                 !
              end do
           end do
           !
           do inode = 1, noddof
             !
             elheat(inode) = elheat(inode) + quadTet % wq(qgauss) * abs(detJ) * dot_product(Jq(1:ndime), matmul(tRES(1:ndime, 1:ndime), Jq(1:ndime))) * pq(inode, qgauss)
             !
           end do
           !
        end do
        !
     elseif (pelty == 37_ip) then
        !
        ! HEX08
        !
        elegau = size(quadHex % wq)
        !
        ! Coordinates of element vertices
        !
        nodi(1:ndime, 1:eledof) = meshe(ndivi) % coord(1:ndime, nodind_mag(1:eledof))
        !
        ! Coordinates of gauss nodes in the reference element
        !
        pq(1:noddof, 1:elegau) = mag_lagref(quadHex % xq, pelty)
        !
        nodq(1:ndime, 1:elegau) = matmul(nodi(1:ndime, 1:eledof), mag_lagref(quadHex % xq, pelty))
        !
        do qgauss = 1, elegau
           !
           ! Vector basis functions at quadrature node and their curl
           !
           call mag_nedhex(quadHex % xq(1:ndime, qgauss), nodi, bq, curl3, detJ)
           ! 
           ! Magnetic field at quadrature node
           !
           Hq = mag_matvec( bq(1:ndime, 1:eledof), Ha_mag(edgind_mag(1:eledof)) * sigval_mag(1:eledof))
           !
           ! Current density at quadrature node
           !
           Jq = mag_matvec(curl3(1:ndime, 1:eledof), Ha_mag(edgind_mag(1:eledof)) * sigval_mag(1:eledof))
           !
           ! Magnetization
           !
           magtiz_mag(1:ndime, pmate) = magtiz_mag(1:ndime, pmate) + 0.5_rp * quadHex % wq(qgauss) * abs(detJ) * &
                mag_vecvec(nodq(1:ndime, qgauss) - momori_mag(1:ndime, pmate), Jq)
           !
           cursum_mag(1:ndime, pmate) = cursum_mag(1:ndime, pmate) + quadHex % wq(qgauss) * abs(detJ) * Jq
           !
           magsum_mag(1:ndime, pmate) = magsum_mag(1:ndime, pmate) + quadHex % wq(qgauss) * abs(detJ) * Hq(1:ndime)
           !
           volume_mag(pmate) = volume_mag(pmate) +  quadHex % wq(qgauss) * abs(detJ)
           !
           ! Resistivity and its derivative at quadrature node
           !
           tRES = reshape( [ mag_resist(pmate, nodq(1:ndime, qgauss), Jq, Hq, 0.d0, .false., 1_ip), 0.0_rp, 0.0_rp, & 
                0.0_rp, mag_resist(pmate, nodq(1:ndime, qgauss), Jq, Hq, 0.d0, .false., 2_ip), 0.0_rp, &
                0.0_rp, 0.0_rp, mag_resist(pmate, nodq(1:ndime, qgauss), Jq, Hq, 0.d0, .false., 3_ip) ], &
                [3, 3] &
                )
           tRDF = reshape( [ mag_resdif(pmate, nodq(1:ndime, qgauss), Jq, Hq, 0.d0, .false., 1_ip), 0.0_rp, 0.0_rp, & 
                0.0_rp, mag_resdif(pmate, nodq(1:ndime, qgauss), Jq, Hq, 0.d0, .false., 2_ip), 0.0_rp, &
                0.0_rp, 0.0_rp, mag_resdif(pmate, nodq(1:ndime, qgauss), Jq, Hq, 0.d0, .false., 3_ip) ], &
                [3, 3] &
                )
           tMUR = reshape( [ mag_permea(pmate, nodq(1:ndime, qgauss), Hq, .false., 1_ip), 0.0_rp, 0.0_rp, &
                0.0_rp, mag_permea(pmate, nodq(1:ndime, qgauss), Hq, .false., 2_ip), 0.0_rp, &
                0.0_rp, 0.0_rp, mag_permea(pmate, nodq(1:ndime, qgauss), Hq, .false., 3_ip) ], &
                [3, 3] &
                )
           !
           ! Source term
           !
           source = mag_source(myClockTime_mag + theta_mag * dt_mag, nodq(1:ndime, qgauss))
           !
           ! Elementary matrices at current quadrature node
           !
           do iedge = 1, eledof
              do jedge = 1, eledof
                 !
                 ! Mass matrix at quadrature node
                 !
                 !            elmass(iedge, jedge) = elmass(iedge, jedge) + quadHex % wq(qgauss) * mu0_mag * mur_mag(pmate) * &
                 !                dot_product(bq(1:ndime, iedge), bq(1:ndime, jedge)) * sigval_mag(iedge) * sigval_mag(jedge) * abs(detJ)

                 elmass(iedge, jedge) = elmass(iedge, jedge) + quadHex % wq(qgauss) * mu0_mag * &
                      dot_product(bq(1:ndime, iedge), matmul(tMUR(1:ndime, 1:ndime), bq(1:ndime, jedge))) * &
                      sigval_mag(iedge) * sigval_mag(jedge) * abs(detJ)
                 !
                 ! Stiffness matrix at quadrature node
                 !
                 !            elstif(iedge, jedge) = elstif(iedge, jedge) + quadHex % wq(qgauss) * &
                 !                res * dot_product(curl3(1:ndime, iedge), curl3(1:ndime, jedge)) * sigval_mag(iedge) * sigval_mag(jedge) * abs(detJ)

                 elstif(iedge, jedge) = elstif(iedge, jedge) + quadHex % wq(qgauss) * &
                      dot_product(curl3(1:ndime, iedge), matmul(tRES(1:ndime, 1:ndime), curl3(1:ndime, jedge))) * &
                      sigval_mag(iedge) * sigval_mag(jedge) * abs(detJ)
                 !
                 ! Stiffness derivative matrix
                 !
                 if (mag_vecmod(Jq) > SCL_MAG * ZER_MAG) then
                    !              elstdf(iedge, jedge) = elstdf(iedge, jedge) + quadHex % wq(qgauss) * resDiff * abs(detJ) / mag_vecmod(Jq) * &
                    !                   dot_product(Jq, curl3(1:ndime, iedge)) * sigval_mag(iedge) * &
                    !                   dot_product(Jq, curl3(1:ndime, jedge)) * sigval_mag(jedge)

                    elstdf(iedge, jedge) = elstdf(iedge, jedge) + quadHex % wq(qgauss) * abs(detJ) / mag_vecmod(Jq) * &
                         dot_product(matmul(tRDF(1:ndime, 1:ndime), Jq), curl3(1:ndime, iedge)) * sigval_mag(iedge) * &
                         dot_product(Jq, curl3(1:ndime, jedge)) * sigval_mag(jedge)
                 else
                    elstdf(iedge, jedge) = 0.0_rp
                 end if
                 !
              end do
           end do
           !
           do inode = 1, noddof
             !
             elheat(inode) = elheat(inode) + quadHex % wq(qgauss) * abs(detJ) * dot_product(Jq(1:ndime), matmul(tRDF(1:ndime, 1:ndime), Jq(1:ndime))) * pq(inode, qgauss)
             !
           end do
           !
        end do
        !
     else
        call runend("mag_elmope: Unknown element type")
     end if
     !
     !---------------------------------------------------------
     ! Original System A(x) x = b(x)
     !---------------------------------------------------------
     !
     ! Elementary Matrix for non-linear problem
     ! A = M + dt * theta * K
     !
     !    Aelem(1:eledof, 1:eledof) = elmass(1:eledof, 1:eledof) + dt_mag * theta_mag * elstif(1:eledof, 1:eledof)
     Aelem(1:eledof, 1:eledof) = a0_mag * elmass(1:eledof, 1:eledof) + dt_mag * theta_mag * elstif(1:eledof, 1:eledof)
     !
     ! Elementary RHS for non-linear problem
     !
     ! b = (M - dt * (1 - theta) * K) * Hp
     ! b = M * Hpav - dt * (1 -theta) * K * Hp
     !
     !    belem(1:eledof) = dt_mag * elsrc(1:eledof) + mag_matvec(elmass(1:eledof, 1:eledof), Hp_mag(edgind_mag(1:eledof))) - &
     !                      dt_mag * (1 - theta_mag) * mag_matvec(elstif(1:eledof, 1:eledof), Hp_mag(edgind_mag(1:eledof)))
     belem(1:eledof) = dt_mag * elsrc(1:eledof) + mag_matvec(elmass(1:eledof, 1:eledof), Hpav_mag(edgind_mag(1:eledof))) - &
          dt_mag * (1 - theta_mag) * mag_matvec(elstif(1:eledof, 1:eledof), Hp_mag(edgind_mag(1:eledof)))
     !
     !---------------------------------------------------------
     ! System from Newton method J(x^(k)) dx = r(x^(k))
     !---------------------------------------------------------
     !
     ! Elementary residual and matrix for Newton's method
     ! r = b - A * H
     !
     felem(1:eledof) = mag_matvec(Aelem(1:eledof, 1:eledof), He_mag(edgind_mag(1:eledof)))
     rnelem(1:eledof) = belem(1:eledof) - felem(1:eledof)
!     rnelem(1:eledof) = belem(1:eledof) - mag_matvec(Aelem(1:eledof, 1:eledof), He_mag(edgind_mag(1:eledof)))
     !
     ! Jacobian
     !
     Anelem(1:eledof, 1:eledof) = Aelem(1:eledof, 1:eledof) + dt_mag * theta_mag * elstdf(1:eledof, 1:eledof)
     !
     !-------------------------------------------------------------------------------------
     ! Dirichlet Boundary Conditions
     !-------------------------------------------------------------------------------------
     !
     do iedge = 1_ip, eledof
        if (edgdir_mag(edgind_mag(iedge))) then
           !-----------------------------------------------------
           ! Original system: Dirichlet B.C.
           !-----------------------------------------------------
!
           piv = Aelem(iedge, iedge)
           !
           ! Modify RHS
           !
           ! bi' = bi - Kid * Hd
           !
           !**DIRBC**
           belem(1:eledof) = belem(1:eledof) - Aelem(1:eledof, iedge) * bvess_mag(1,edgind_mag(iedge),1)
           belem(iedge) = piv * bvess_mag(1,edgind_mag(iedge),1)

!           belem(1:eledof) = belem(1:eledof) - Aelem(1:eledof, iedge) * He_mag(edgind_mag(iedge))
!           belem(iedge) = piv * He_mag(edgind_mag(iedge))

           !**DIRBC**
           !
           ! Modify matrix
           !
           Aelem(iedge, 1:eledof) = 0.0_rp
           Aelem(1:eledof, iedge) = 0.0_rp
           Aelem(iedge, iedge) = piv
!
           !
           !-----------------------------------------------------
           ! Newton's method: Dirichlet B.C.
           !-----------------------------------------------------
           !
           !**DIRBC**
           ! ALGEBRAIC_SOLVER
           ! OPTION: FIXITY - kfl_iffix = 1 - Dirichlet B.C. imposed by solver
           ! OPTION: ZERO_FIXITY kfl_iffix = 2 - Null Dirichlet B.C imposed by solver
           !**DIRBC**
           ! Keep commented for the moment as RHS will always be zero at boundary edges
!           if( solve(1) % kfl_iffix == 0 ) then 
              !
              ! Modify RHS
              !
              rnelem(iedge) = 0.0_rp
              !
              ! Modify matrix
              !
              piv                     = Anelem(iedge, iedge)
              Anelem(iedge, 1:eledof) = 0.0_rp
              Anelem(1:eledof, iedge) = 0.0_rp
              Anelem(iedge, iedge)    = piv
!           end if
           !-----------------------------------------------------
        end if
     end do
     !-------------------------------------------------------------------------------------
     !
     ! SPD tests
     ! GGU: Commented lines due to a segmentation fault in line 877
     !if (dot_product(rnelem(1:eledof), mag_matvec(anelem(1:eledof, 1:eledof), rnelem(1:eledof))) < 0.0_rp &
     !     .and. dot_product(rnelem(1:eledof), rnelem(1:eledof)) > 0.0_rp) then
     !   call runend("mag_elmope: Matrix is not PD")
     !end if
     !
     !    if (mag_vecmod(mag_matvec(Anelem(1:eledof, 1:eledof), rnelem(1:eledof)) - & 
     !    mag_matvec(transpose(Anelem(1:eledof, 1:eledof)), rnelem(1:eledof))) >= 1.0e-15_rp) then
     !      call runend("mag_elmope: Matrix is not symmetric")
     !    end if
     !
     !-------------------------------------------------------------------------------------
     ! Global Matrix Assembly
     !-------------------------------------------------------------------------------------
     !
     !call solver_assemble_element_matrix(&
     !      solve_sol,1_ip,size(ee),size(ee),ielem,ee,Anelem,amatr)
     do iedge = 1_ip, eledof
        k1 = r_edg(edgind_mag(iedge))
        k2 = r_edg(edgind_mag(iedge) + 1_ip) - 1_ip
        do jedge = 1, eledof
           flag = .false.
           do j1 = k1, k2
              if (edgind_mag(jedge) == c_edg(j1)) then
                 amatr(j1) = amatr(j1) + Anelem(iedge, jedge)
                 flag = .true.
                 exit
              end if
           end do
           if (.not. flag) then
              call runend("mag_elmope: place is missing!")
           end if
        end do
        !
        rhsid(edgind_mag(iedge)) = rhsid(edgind_mag(iedge)) + rnelem(iedge)
        bhsid(edgind_mag(iedge)) = bhsid(edgind_mag(iedge)) + belem(iedge)
        fhsid(edgind_mag(iedge)) = fhsid(edgind_mag(iedge)) + felem(iedge)
        !
     end do
     !
     do inode = 1_ip, noddof
       !
       hhsid(nodind_mag(inode)) = hhsid(nodind_mag(inode)) + elheat(inode)
       !
     end do
     !
     ! Energy calculation
     !
     magnen_mag(pmate) = magnen_mag(pmate) + &
          dot_product(Ha_mag(edgind_mag(1:eledof)), mag_matvec(elmass(1:eledof, 1:eledof), Ha_mag(edgind_mag(1:eledof)))) / 2.0_rp
     !
     joulen_mag(pmate) = joulen_mag(pmate) + &
          dot_product(Ha_mag(edgind_mag(1:eledof)), mag_matvec(elstif(1:eledof, 1:eledof), Ha_mag(edgind_mag(1:eledof))))
     !
  end do
  !
end subroutine mag_elmope
