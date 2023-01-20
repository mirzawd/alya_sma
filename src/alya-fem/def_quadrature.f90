!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Domain
!> @{
!> @file    def_mesh_type.f90
!> @author  houzeaux
!> @date    2020-04-24
!> @brief   Mesh integration 
!> @details Mesh integration strategy
!>   
!-----------------------------------------------------------------------

module def_quadrature

  use def_kintyp_basic,           only : ip,rp
  use mod_memory_basic,           only : memory_alloca
  use mod_memory_basic,           only : memory_deallo
  !
  ! Gauss-Legendre quadrature rules
  !
  use mod_quadrature_gauss,       only : gauss_poi
  use mod_quadrature_gauss,       only : gauss_bar
  use mod_quadrature_gauss,       only : gauss_qua
  use mod_quadrature_gauss,       only : gauss_tri
  use mod_quadrature_gauss,       only : gauss_pen
  use mod_quadrature_gauss,       only : gauss_pyr
  !
  ! Closed quadrature rules
  !
  use mod_quadrature_closed,      only : close_poi
  use mod_quadrature_closed,      only : close_bar
  use mod_quadrature_closed,      only : close_qua
  use mod_quadrature_closed,      only : close_tri
  use mod_quadrature_closed,      only : close_pen
  use mod_quadrature_closed,      only : close_pyr
  !
  ! Trapezoidal rules
  !
  use mod_quadrature_trapezoidal, only : trape_bar
  use mod_quadrature_trapezoidal, only : trape_qua
  !
  ! Chebyshev quadrature (closed rule)
  !
  use mod_quadrature_chebyshev,   only : cheby_qua

  implicit none
  private

  integer(ip), parameter  :: GAUSS_LEGENDRE_RULE = 0
  integer(ip), parameter  :: TRAPEZOIDAL_RULE    = 2
  integer(ip), parameter  :: CLOSED_RULE         = 1
  integer(ip), parameter  :: CHEBYSHEV_RULE      = 3

  type :: quadrature
     integer(ip)          :: ndime            ! Space dimension
     integer(ip)          :: ngaus            ! Number of Gauss points
     integer(ip)          :: type             ! Type of rule
     integer(ip)          :: topo             ! Element topology
     real(rp),    pointer :: weigp(:)         ! Weights
     real(rp),    pointer :: posgp(:,:)       ! Integration points positions
     integer(8)           :: memor(2)         ! Memory counter
   contains
     procedure,      pass :: init             ! Initialization
     procedure,      pass :: alloca           ! Allocate     
     procedure,      pass :: deallo           ! Deallocate  
     procedure,      pass :: set              ! Set integration rule   
  end type quadrature

  public :: set_quadrature
  public :: quadrature
  public :: GAUSS_LEGENDRE_RULE   
  public :: TRAPEZOIDAL_RULE      
  public :: CLOSED_RULE           
  public :: CHEBYSHEV_RULE           

contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-01-26
  !> @brief   Set integration rule
  !> @details Set integration rule
  !> 
  !-----------------------------------------------------------------------

  subroutine set(self,ndime,ngaus,name,topo,ierr)
    
    class(quadrature),            intent(inout) :: self
    integer(ip),        optional, intent(in)    :: ndime
    integer(ip),        optional, intent(in)    :: ngaus
    integer(ip),        optional, intent(in)    :: topo
    character(LEN=*),   optional, intent(in)    :: name
    integer(ip),        optional, intent(out)   :: ierr
    integer(ip)                                 :: ierro

    if( present(ndime) ) self % ndime = ndime
    if( present(ngaus) ) self % ngaus = ngaus
    if( present(topo ) ) self % topo  = topo
    if( present(name ) ) then
       select case ( name )
       case ( 'GAUSS-LEGENDRE'   ) ; self % type = GAUSS_LEGENDRE_RULE
       case ( 'CLOSE' , 'CLOSED' ) ; self % type = CLOSED_RULE
       case ( 'TRAPEZOIDAL'      ) ; self % type = TRAPEZOIDAL_RULE
       end select
    end if
    !
    ! Allocate rule
    !
    call self % alloca()
    ierro = 1

    call set_quadrature(self % ndime,self % ngaus,self % topo,self % type,self % posgp,self % weigp,ierro)

    if( present(ierr) ) then
       ierr = ierro
    else
       if( ierro /= 0 ) stop
    end if
    
  end subroutine set
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-01-26
  !> @brief   Initialization
  !> @details Initialization
  !> 
  !-----------------------------------------------------------------------

  subroutine init(self)

    class(quadrature), intent(inout) :: self

    self % ndime = 0 
    self % ngaus = 0 
    self % type  = 0
    self % topo  = 0 
    self % memor = 0_8 
    nullify(self % weigp)
    nullify(self % posgp)

  end subroutine init

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-01-26
  !> @brief   Allocate
  !> @details Allocate
  !> 
  !-----------------------------------------------------------------------

  subroutine alloca(self)

    class(quadrature), intent(inout) :: self

    call memory_alloca(self % memor,'SELF % WEIGP','alloca',self % weigp,self % ngaus)
    call memory_alloca(self % memor,'SELF % POSGP','alloca',self % posgp,max(1_ip,self % ndime),self % ngaus)

  end subroutine alloca

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-01-26
  !> @brief   Allocate
  !> @details Allocate
  !> 
  !-----------------------------------------------------------------------

  subroutine deallo(self)

    class(quadrature), intent(inout) :: self

    call memory_deallo(self % memor,'SELF % WEIGP','alloca',self % weigp)
    call memory_deallo(self % memor,'SELF % POSGP','alloca',self % posgp)

  end subroutine deallo

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-01-26
  !> @brief   Set quadrature rule
  !> @details Set quadrature rule
  !> 
  !-----------------------------------------------------------------------

  subroutine set_quadrature(ndime,ngaus,topo,type,posgp,weigp,ierro)
    
    integer(ip),        intent(in)    :: ndime
    integer(ip),        intent(in)    :: ngaus
    integer(ip),        intent(in)    :: topo
    integer(ip),        intent(in)    :: type
    real(rp),           intent(out)   :: posgp(max(1_ip,ndime),ngaus)
    real(rp),           intent(out)   :: weigp(ngaus)
    integer(ip),        intent(out)   :: ierro

    select case ( type )
       
    case ( GAUSS_LEGENDRE_RULE )
       !
       ! This choice of quadrature weights wi and quadrature nodes xi
       ! is the unique choice that allows the quadrature rule to integrate
       ! degree 2*ngaus − 1 polynomials exactly.
       !
       select case (  topo ) 
       case(-2_ip) ; call gauss_poi(        ngaus,        weigp, ierro) ! Point
       case(-1_ip) ; call gauss_bar(        ngaus, posgp, weigp, ierro) ! Bar         
       case( 0_ip) ; call gauss_qua( ndime, ngaus, posgp, weigp, ierro) ! Quad/Hexa   
       case( 1_ip) ; call gauss_tri( ndime, ngaus, posgp, weigp, ierro) ! Tria/Tetra    
       case( 2_ip) ; call gauss_pen( ndime, ngaus, posgp, weigp, ierro) ! Prism(Penta)  
       case( 3_ip) ; call gauss_pyr( ndime, ngaus, posgp, weigp, ierro) ! Pyramid      
       end select
       
    case ( CLOSED_RULE )
       !
       ! Uses a close rule where Gauss-points are on the nodes
       !
       select case (  topo ) 
       case(-2_ip) ; call close_poi(        ngaus,        weigp, ierro) ! Point
       case(-1_ip) ; call close_bar(        ngaus, posgp, weigp, ierro) ! Bar         
       case( 0_ip) ; call close_qua( ndime, ngaus, posgp, weigp, ierro) ! Quad/Hexa   
       case( 1_ip) ; call close_tri( ndime, ngaus, posgp, weigp, ierro) ! Tria/Tetra    
       case( 2_ip) ; call close_pen( ndime, ngaus, posgp, weigp, ierro) ! Prism(Penta)  
       case( 3_ip) ; call close_pyr( ndime, ngaus, posgp, weigp, ierro) ! Prism(Penta)  
       end select
       
    case ( TRAPEZOIDAL_RULE )
       !
       ! Uses a close rule where Gauss-points are on the nodes
       !
       select case (  topo ) 
       case(-1_ip) ; call trape_bar(        ngaus, posgp, weigp, ierro) ! Bar  
       case( 0_ip) ; call trape_qua( ndime, ngaus, posgp, weigp, ierro) ! Quad/Hexa 
       end select
       
    case ( CHEBYSHEV_RULE  )
       !
       ! Uses a close rule where Gauss-points are on the nodes
       !
       select case ( topo ) 
       case( 0_ip) ; call cheby_qua( ndime, ngaus, posgp, weigp, ierro) ! Quad/Hexa 
       end select

    end select

  end subroutine set_quadrature
  
end module def_quadrature
!> @}
