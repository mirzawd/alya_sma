!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Wtiness
!> @{
!> @file    def_kintyp_witness_mesh.f90
!> @author  houzeaux
!> @date    2020-10-27
!> @brief   Witness mesh
!> @details Class defintion of witness mesh
!-----------------------------------------------------------------------

module def_kintyp_witness_mesh
  
  use def_kintyp_basic,       only : ip
  use def_kintyp_mesh,        only : mesh_type_basic
  use def_kintyp_postprocess, only : witness_geo
  use def_kintyp_basic,       only : typ_interp
  use mod_memory,             only : memory_alloca
  use mod_memory,             only : memory_deallo
  implicit none
  private
  
  type typ_witness_mesh
     type(mesh_type_basic)  :: mesh                           ! Witness mesh
     type(witness_geo)      :: geom                           ! Witness mesh geometry
     type(typ_interp)       :: inte                           ! Witness interpolation
     character(200)         :: name                           ! Witness mesh name
     integer(8)             :: memor(2)                       ! Memory counter
     integer(ip), pointer   :: perm(:)                        ! Witness permutation
   contains
     procedure,   pass      :: init    => init_witness_mesh   ! Initialization     
     procedure,   pass      :: deallo  => deallo_witness_mesh ! Deallocate     
  end type typ_witness_mesh

  public :: typ_witness_mesh

contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-10-27
  !> @brief   Initialization
  !> @details Initialization of class WITNESS_MESH
  !> 
  !-----------------------------------------------------------------------

  subroutine init_witness_mesh(self)
    
    class(typ_witness_mesh) :: self
    
    call self % mesh % init()
    call self % geom % init()
    call self % inte % init()
    self % name  = 'NULL'
    self % memor = 0_8
    nullify(self % perm)
    
  end subroutine init_witness_mesh
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-10-27
  !> @brief   Initialization
  !> @details Initialization of class WITNESS_MESH
  !> 
  !-----------------------------------------------------------------------

  subroutine deallo_witness_mesh(self)
    
    class(typ_witness_mesh) :: self
    
    call self % mesh % deallo()
    call self % geom % deallo()
    call self % inte % deallo()
    call memory_deallo(self % memor,'self % perm','def_kintyp_witness_mesh',self % perm)

    nullify(self % perm)
    
  end subroutine deallo_witness_mesh
  
end module def_kintyp_witness_mesh
!> @}
