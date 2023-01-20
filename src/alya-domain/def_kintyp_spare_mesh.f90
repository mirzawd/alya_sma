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
!> @details Class defintion of spare mesh
!-----------------------------------------------------------------------

module def_kintyp_spare_mesh
  
  use def_kintyp_basic,       only : ip,i1p,r2p,rp
  use def_kintyp_mesh_basic,  only : mesh_type_basic
  use def_kintyp_mesh,        only : mesh_type
  use mod_memory,             only : memory_alloca
  use mod_memory,             only : memory_deallo
  implicit none
  private
  
  type typ_spare_mesh
     type(mesh_type_basic)           :: mesh                ! Spare mesh
     type(i1p),             pointer  :: eleme(:)
     type(r2p),             pointer  :: shapf(:)
     real(rp),              pointer  :: dista(:)
     integer(8)                      :: memor(2)
   contains
     procedure,   pass      :: init    => init_spare_mesh   ! Initialization     
     procedure,   pass      :: deallo  => deallo_spare_mesh ! Deallocate     
  end type typ_spare_mesh
 
  public :: typ_spare_mesh

contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-10-27
  !> @brief   Initialization
  !> @details Initialization of class SPARE_MESH
  !> 
  !-----------------------------------------------------------------------

  subroutine init_spare_mesh(self)
    
    class(typ_spare_mesh) :: self
    
    call self % mesh % init()
    nullify(self % eleme)
    nullify(self % shapf)
    nullify(self % dista)
    self % memor = 0_8
    
  end subroutine init_spare_mesh
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-10-27
  !> @brief   Initialization
  !> @details Initialization of class SPARE_MESH
  !> 
  !-----------------------------------------------------------------------

  subroutine deallo_spare_mesh(self)
    
    class(typ_spare_mesh) :: self
    
    call self % mesh % deallo()
    call memory_deallo(self % memor,'self % eleme','def_kintyp_spare_mesh',self % eleme)
    call memory_deallo(self % memor,'self % shapf','def_kintyp_spare_mesh',self % shapf)
    call memory_deallo(self % memor,'self % dista','def_kintyp_spare_mesh',self % dista)

  end subroutine deallo_spare_mesh
     
end module def_kintyp_spare_mesh
!> @}
