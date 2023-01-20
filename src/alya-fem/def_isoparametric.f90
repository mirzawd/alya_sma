!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Kinds_and_types
!> @{
!> @file    def_isoparametric.f90
!> @author  houzeaux
!> @date    2020-04-04
!> @brief   Iso-parametric definition
!> @details Shape funcion, derivative and Hessian
!-----------------------------------------------------------------------

module def_isoparametric

  use def_kintyp_basic, only : ip,rp
  use mod_memory_basic, only : memory_alloca
  use mod_memory_basic, only : memory_deallo
  use mod_lagrange,     only : shape0
  use mod_lagrange,     only : shape1
  use mod_lagrange,     only : shape2
  use mod_lagrange,     only : shape3
  use mod_chebyshev,    only : cheby2
  use mod_chebyshev,    only : cheby3
  use def_quadrature,   only : quadrature
  
  integer(ip), parameter :: LAGRANGE_INTERPOLATION  = 0
  integer(ip), parameter :: CHEBYSHEV_INTERPOLATION = 1

  type isoparametric
     integer(ip)          :: ndime
     integer(ip)          :: ntens
     integer(ip)          :: nnode
     integer(ip)          :: ngaus
     integer(ip)          :: inter            ! Interpolation
     real(rp),    pointer :: shape(:,:)       ! pnode,pgaus
     real(rp),    pointer :: deriv(:,:,:)     ! ndime,pnode,pgaus
     real(rp),    pointer :: heslo(:,:,:)     ! ntens,pnode,pgaus
     integer(8)           :: memor(2)         ! Memory counter
   contains
     procedure,   pass    :: init             ! Initialize
     procedure,   pass    :: alloca           ! Allocate
     procedure,   pass    :: deallo           ! Deallocate
     procedure,   pass    :: set              ! Set 
  end type isoparametric

  private

  public :: isoparametric
  public :: set_isoparametric
  public :: shape0
  public :: shape1
  public :: shape2
  public :: shape3
  public :: cheby2
  public :: cheby3
  public :: LAGRANGE_INTERPOLATION
  public :: CHEBYSHEV_INTERPOLATION

contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-02-12
  !> @brief   Initialize
  !> @details Initialize
  !> 
  !-----------------------------------------------------------------------

  subroutine init(self)
    
    class(isoparametric), intent(inout) :: self

    self % ndime = 0
    self % nnode = 0
    self % ngaus = 0
    self % ntens = 0
    self % memor = 0_8
    self % inter = LAGRANGE_INTERPOLATION
    
    nullify(self % shape)
    nullify(self % deriv)
    nullify(self % heslo)

  end subroutine init

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-01-26
  !> @brief   Allocate
  !> @details Allocate and set iso parametric functions
  !>          SHAPE Ni(g), DERIV and HESLO of node i at gauss point g.
  !>
  !>          The derivatives DERIV(:,i,g) components are:
  !>
  !>          DERIV(1,i,g) = d N_i(g) / ds
  !>          DERIV(2,i,g) = d N_i(g) / dt
  !>          DERIV(3,i,g) = d N_i(g) / dz
  !>
  !>          The Hessian HESLO(:,i,g) components are::
  !>
  !>          In 2D:
  !>
  !>          HESLO(1,i,g) = d^2 N_i(g) / ds^2
  !>          HESLO(2,i,g) = d^2 N_i(g) / dt^2
  !>          HESLO(3,i,g) = d^2 N_i(g) / dsdt
  !>
  !>          In 3D:
  !>
  !>          HESLO(1,i,g) = d^2 N_i(g) / ds^2
  !>          HESLO(2,i,g) = d^2 N_i(g) / dt^2
  !>          HESLO(3,i,g) = d^2 N_i(g) / dz^2
  !>          HESLO(4,i,g) = d^2 N_i(g) / dsdt
  !>          HESLO(5,i,g) = d^2 N_i(g) / dsdz
  !>          HESLO(6,i,g) = d^2 N_i(g) / dtdz
  !> 
  !-----------------------------------------------------------------------

  subroutine set(self,nnode,quad)

    class(isoparametric),  intent(inout) :: self
    integer(ip),           intent(in)    :: nnode
    type(quadrature),      intent(in)    :: quad
    integer(ip)                          :: igaus
    
    self % nnode = nnode
    self % ndime = quad % ndime
    self % ngaus = quad % ngaus
    self % ntens = max(1_ip,3 * self % ndime - 3)

    call self % alloca()

    do igaus = 1,self % ngaus
       call set_isoparametric(self % ndime,self % nnode,self % inter,&
            quad % posgp(:,igaus),self % shape(:,igaus),&
            self % deriv(:,:,igaus),self % heslo(:,:,igaus))
    end do
   
  end subroutine set

  pure subroutine set_isoparametric(ndime,nnode,inter,posgp,shapf,deriv,heslo,ierro)

    integer(ip),           intent(in)    :: ndime
    integer(ip),           intent(in)    :: nnode
    integer(ip),           intent(in)    :: inter
    real(rp),              intent(in)    :: posgp(ndime)
    real(rp),              intent(out)   :: shapf(nnode)
    real(rp),    optional, intent(out)   :: deriv(ndime,nnode)
    real(rp),    optional, intent(out)   :: heslo(max(1_ip,3_ip*ndime-3_ip),nnode)
    integer(ip), optional, intent(inout) :: ierro

    select case ( inter )

    case ( LAGRANGE_INTERPOLATION ) 

       select case ( ndime )
       case ( 0_ip ) ; call shape0(              nnode,shapf,            ierro)
       case ( 1_ip ) ; call shape1(  posgp(  1), nnode,shapf,deriv,heslo,ierro)
       case ( 2_ip ) ; call shape2([ posgp(1:2)],nnode,shapf,deriv,heslo,ierro)
       case ( 3_ip ) ; call shape3([ posgp(1:3)],nnode,shapf,deriv,heslo,ierro)
       end select
       
    case ( CHEBYSHEV_INTERPOLATION )
       
       select case ( ndime )
       case ( 2_ip ) ; call cheby2([ posgp(1:2)],nnode,shapf,deriv,heslo,ierro)
       case ( 3_ip ) ; call cheby3([ posgp(1:3)],nnode,shapf,deriv,heslo,ierro)
       end select
       
    end select
    
  end subroutine set_isoparametric

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2022-01-26
  !> @brief   Allocate
  !> @details Allocate
  !> 
  !-----------------------------------------------------------------------

  subroutine alloca(self)

    class(isoparametric), intent(inout) :: self

    call memory_alloca(self % memor,'SELF % SHAPE','alloca',self % shape,                       self % nnode,self % ngaus)
    call memory_alloca(self % memor,'SELF % DERIV','alloca',self % deriv,max(1_ip,self % ndime),self % nnode,self % ngaus)
    call memory_alloca(self % memor,'SELF % HESLO','alloca',self % heslo,max(1_ip,self % ntens),self % nnode,self % ngaus)

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

    class(isoparametric), intent(inout) :: self

    call memory_deallo(self % memor,'SELF % SHAPE','deallo',self % shape)
    call memory_deallo(self % memor,'SELF % DERIV','deallo',self % deriv)
    call memory_deallo(self % memor,'SELF % HESLO','deallo',self % heslo)

  end subroutine deallo

end module def_isoparametric
!> @}
