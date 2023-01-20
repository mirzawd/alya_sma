!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine mag_turnof()
  !-----------------------------------------------------------------------
  !****f* Magnet/mag_turnof
  ! NAME
  !    mag_turnof
  ! DESCRIPTION
  !    This routine closes the run for the H-formulation of Maxwell's 
  !    equations.
  ! USES
  ! USED BY
  !    Magnet
  !***
  !-----------------------------------------------------------------------

  use def_magnet
  use def_master, only: INOTMASTER
  
  implicit none

  if (INOTMASTER) then
    deallocate(node_to_nedge_mag)
    deallocate(r_node_to_edge_mag)
    deallocate(c_node_to_edge_mag)
    deallocate(i_node_to_edge_mag)
    deallocate(edge_to_nedge_mag)
    deallocate(r_edge_to_edge_mag)
    deallocate(i_edge_to_edge_mag)
    deallocate(c_edge_to_edge_mag)
  end if

  deallocate(belem)
  deallocate(rnelem)
  deallocate(felem)

  deallocate(elstif)
  deallocate(elstdf)
  deallocate(elmass)

  deallocate(elsrc)

  deallocate(Aelem)
  deallocate(Anelem)

  deallocate(edgind_mag)
  deallocate(nodind_mag)
  deallocate(sigval_mag)
  deallocate(lenval_mag)

  deallocate(magnen_mag)
  deallocate(joulen_mag)
  deallocate(magtiz_mag)
  deallocate(cursum_mag)
  deallocate(magsum_mag)
  deallocate(volume_mag)
  deallocate(momori_mag)

  deallocate(Ec0_mag)
  deallocate(nc0_mag)
  deallocate(Jc0_mag)
  deallocate(mur_mag)
  deallocate(rho_mag)
  deallocate(B0k_mag)

  deallocate(refres_mag)
  deallocate(residu_mag)
  deallocate(refval_mag)
  deallocate(delval_mag)

end subroutine mag_turnof
