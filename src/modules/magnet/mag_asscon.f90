!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine mag_asscon()
  !-----------------------------------------------------------------------
  !****f* magnet/mag_asscon.f90
  ! NAME
  !    mag_asscon
  ! DESCRIPTION
  !    This routine
  !    1. Computes the elementary matrix and RHS for every element in the mesh
  !    2. Assembles the elementary equations into the system equations
  ! USES
  ! USED BY
  !    mag_begste
  !    mag_matrix
  !***
  !-----------------------------------------------------------------------
  use def_domain, only: meshe, ndime
  use def_magnet
  use def_kermod, only: ndivi
  use mod_mag_nedele, only: mag_edgdof, mag_nedtri2, mag_nedqua

  implicit none

  integer(ip) :: &
       ielem,    &
       eledof,    &
       pelty,    &
       pmate,    &
       iconstr, idof, qgauss

  real(rp) :: &
       celem(maxdof_mag),    &
       nodi(ndime, maxdof_mag),    &
       bq(ndime, maxdof_mag),    &
       curl2(maxdof_mag),    &
       detJ

  !
  ! Initialize constraint global vector
  !
  do iconstr = 1, constr_total
     constr_mag(iconstr) % cnstr = 0.0_rp
  end do
  !
  do ielem = 1, meshe(ndivi) % nelem
     !
     !---------------------------------------------------------------------------------
     ! Element Data
     !---------------------------------------------------------------------------------
     !
     ! Element type
     ! TRI03: 10
     ! QUA04: 12
     !
     pelty = meshe(ndivi) % ltype(ielem)
     !
     pmate = meshe(ndivi) % lmate(ielem)
     !
     if (constrlist_mag(pmate) <= 0) cycle
     iconstr = constrlist_mag(pmate)
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
     if (struct_mag) then
        !
        ! Edge lengths
        !
        lenval_mag(1:eledof) = edglen_mag(edgind_mag(1:eledof))
        !
        celem(1:eledof) = lenval_mag(1:eledof) * sigval_mag(1:eledof)
        !
     else
        !
        ! Coordinates of element vertices
        !
        nodi(1:ndime, 1:eledof) = meshe(ndivi) % coord(1:ndime, nodind_mag(1:eledof))
        !
        celem(1:eledof) = 0.0_rp
        !
        select case(pelty)
           !
        case (10_ip)
           !
           do qgauss = 1, size(quadTri % wq)
              !
              call mag_nedtri2(quadTri % xq(1:ndime, qgauss), nodi(1:ndime, 1:eledof), bq, curl2, detJ)
              !
              do idof = 1, eledof
                 celem(idof) = celem(idof) + quadTri % wq(qgauss) * curl2(idof) * sigval_mag(idof) * abs(detJ)
              end do
              !
           end do
           !
        case (12_ip)
           !
           do qgauss = 1, size(quadQua % wq)
              !
              call mag_nedqua(quadQua % xq(1:ndime, qgauss), nodi(1:ndime, 1:eledof), bq, curl2, detJ)
              !
              do idof = 1, eledof
                 celem(idof) = celem(idof) + quadQua % wq(qgauss) * curl2(idof) * sigval_mag(idof) * abs(detJ)
              end do
              !
           end do
           !
        case default
           call runend('mag_asscon: unknown element type')
        end select
        !
     end if
     !
     ! Assembly in global vector
     !
     do idof = 1, eledof
        constr_mag(iconstr) % cnstr(edgind_mag(idof)) = constr_mag(iconstr) % cnstr(edgind_mag(idof)) + celem(idof)
     end do

     !
  end do
  !
end subroutine mag_asscon
