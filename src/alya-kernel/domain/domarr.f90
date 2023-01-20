!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine domarr(itask)
  !-----------------------------------------------------------------------
  !****f* domain/domarr
  ! NAME
  !    domarr
  ! DESCRIPTION
  !    This routines computes some arrays that depend the mesh.
  !    If the mesh changes, all these arrays should be recomputed.
  !
  !    VMASS ... Lumped mass matrix
  !    VMASC ... Close rule mass matrix
  !    EXNOR ... Exterior normals
  !    YWALP ... Physical nodal wall distance (includes distance to wall)
  !    YWALB ... Physical boundary wall distance (includes distance to wall)
  !    WALLD ... Mesh nodal wall distance 
  !    SKCOS ... Geometrical local basis
  !    COORD ... Fringe coordinates
  !
  !    ITASK = 1 ... Allocate and compute all variables
  !          = 2 ... Recompute all variables except LPOTY
  !
  !    The order of creation of geometrical arrays is the following
  !    because some solver variables must be allocated:
  !   
  !    From domain():
  !    1. domarr(1)    computes: VMASS, VMASC, EXNOR, LPOTY, YWALP, YWALD
  !    2. opebcs()     computes: KGL_GEONO, SKCOS
  !    From Iniunk(): 
  !    3. ker_iniunk() computes: WALLD
  !    
  !    All these variables are then updated together when mesh is changing.
  !    LPOTY is not recomputed as boundary nodes are assumed not to change.
  !
  ! USED BY
  !    Turnon
  !***
  !-----------------------------------------------------------------------
  use def_master
  use def_kermod
  use def_domain
  use mod_mass_matrix,       only : mass_matrix_open_lumped
  use mod_mass_matrix,       only : mass_matrix_close
  use mod_mass_matrix,       only : mass_matrix_diagonal_lumped
  use mod_mass_matrix,       only : mass_matrix_consistent
  use mod_mass_matrix,       only : mass_matrix_approx_inv_mass
  use mod_communications,    only : PAR_GHOST_NODE_EXCHANGE
  use mod_exterior_normal,   only : extnor
  use mod_local_basis,       only : local_basis
  use mod_par_bin_structure, only : par_bin_structure
  use mod_par_bin_structure, only : par_bin_structure_deallocate
  use mod_elsest,            only : elsest_deallocate
  use mod_messages,          only : messages_live
  implicit none
  integer(ip), intent(in) :: itask

  if( itask == 3 ) then
     !
     ! Mesh arrays that should be computed in the Alya world
     !
     ! This was a good idea but this does not work if we have a coupling between
     ! different instances of Alya and one of the codes does several time steps
     ! before the coupling takes place. The bin structure should thus be recomputed
     ! at a point where codes can be synchronized
     !
     !call par_bin_structure_deallocate()
     !call par_bin_structure()

  else
     !
     ! Fringe coordinates
     !
     if( ISLAVE .and. itask == 2 ) call PAR_GHOST_NODE_EXCHANGE(coord,'SUBSTITUTE','IN MY CODE')
     !
     ! VMASS: Diagonal mass matrix using lumping 
     ! 
     call mass_matrix_open_lumped()
     !
     ! VMASC: Diagonal mass matrix using close rule
     !
     call mass_matrix_close()
     !
     ! DMASS: Diagonaly scaled lumped matrix
     !
     call mass_matrix_diagonal_lumped()
     !
     ! Consistent mass matrix
     !
     call mass_matrix_consistent(CONSISTENT_MASS=.true.)
     !
     ! Approximate inverse
     !
     call mass_matrix_approx_inv_mass()
     !
     ! Check projections
     !
     !call chkmas()  
     !
     ! EXNOR: External normal
     !
     call extnor(itask)
     !
     ! COSYS(:) % LOCAL_BASIS
     !
     call local_basis()  
     !
     ! YWALP and YWALB: Physical wall distance on boundary nodes
     !
     call waldis(itask)

     if( itask == 2 ) then
        !
        ! Elsest
        !
        call elsest_deallocate(ielse)
        call elsini()
        !
        ! KFL_GEONO, SKCOS: Geometrical boundary conditions
        !
        call geonor(itask)
        !
        ! WALLD, WALLN: Mesh wall distance and normal. Computed by Kermod: must go through
        !        main subroutine (it uses a solver)
        !
        call ker_walgen(2_ip)
        call ker_walnor(2_ip)
        !
        ! Element bin 
        !
        call elebin()
     end if

  end if

end subroutine domarr
