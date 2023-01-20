!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Kermod
!> @{
!> @file    def_ker_proper.f90
!> @author  houzeaux
!> @date    2020-11-23
!> @brief   Properties
!> @details Definitions of properties
!-----------------------------------------------------------------------

module def_ker_proper

  use def_kintyp_basic, only : ip,rp,lg,r1p,r2p,r3p
  use def_domain,       only : ndime,mnode
  implicit none
  public

  integer(ip),   parameter       :: mlaws_ker            = 20
  integer(ip),   parameter       :: mresp_ker            = 5
  integer(ip),   parameter       :: mlapa_ker            = 9  ! Max # parameters
  integer(ip),   parameter       :: NUMBER_OF_PROPERTIES = 12
  integer(ip),   parameter       :: SCALAR_PROPERTY      = 0
  integer(ip),   parameter       :: VECTOR_PROPERTY      = 1
  integer(ip),   parameter       :: MATRIX_PROPERTY      = 2
  integer(ip),   parameter       :: UPDATE_PROPERTIES    = 1
  real(rp),      parameter       :: eps                  = epsilon(1.0_rp)
  !
  ! Class for gather operations
  !
  type elgat_typ
     logical(lg)         :: kfl_elcod
     logical(lg)         :: kfl_elvel
     logical(lg)         :: kfl_elwal
     logical(lg)         :: kfl_eltem
     logical(lg)         :: kfl_eltur
     logical(lg)         :: kfl_elcon
     real(rp),   pointer :: elcod(:,:,:) ! coord
     real(rp),   pointer :: elvel(:,:,:) ! veloc
     real(rp),   pointer :: elwal(:,:)   ! walld
     real(rp),   pointer :: eltem(:,:)   ! tempe
     real(rp),   pointer :: eltur(:,:,:) ! untur
     real(rp),   pointer :: elcon(:,:,:) ! conce
   contains
     procedure,  pass    :: init     
     procedure,  pass    :: alloca     
     procedure,  pass    :: deallo   
  end type elgat_typ
  !
  ! Interface for initialization
  ! Decide which variables should be gathered
  !
  abstract interface
     subroutine xxx_init(elgat)
       
       import elgat_typ
       
       type(elgat_typ), intent(inout) :: elgat
       
     end subroutine xxx_init
  end interface
  !
  ! Interface for computing the property
  !
  abstract interface
     subroutine xxx_operations(&
          VECTOR_DIM,list_elements_loc,pelty,pnode,pgaus,gpsha,gpcar,&
          hleng,self,param,gppro,grpro)
       
       import ip,rp,elgat_typ
       import ndime
       import mnode
       import mlapa_ker
       
       integer(ip),     intent(in)    :: VECTOR_DIM 
       integer(ip),     intent(in)    :: list_elements_loc(VECTOR_DIM) 
       integer(ip),     intent(in)    :: pelty
       integer(ip),     intent(in)    :: pnode
       integer(ip),     intent(in)    :: pgaus
       real(rp),        intent(in)    :: gpsha(pnode,pgaus)
       real(rp),        intent(in)    :: gpcar(VECTOR_DIM,ndime,mnode,pgaus)
       real(rp),        intent(in)    :: hleng(VECTOR_DIM,ndime)
       type(elgat_typ), intent(in)    :: self
       real(rp),        intent(in)    :: param(mlapa_ker)
       real(rp),        intent(inout) :: gppro(VECTOR_DIM,pgaus)
       real(rp),        intent(inout) :: grpro(VECTOR_DIM,ndime,pgaus)
       
     end subroutine xxx_operations
  end interface
        
  real(rp)                       :: cpu_prope(2)

  type typ_laws
     integer(ip)                 :: kfl_gradi
     integer(ip)                 :: kfl_deriv
     integer(ip)                 :: kfl_grder
     integer(ip)                 :: kfl_deriv_tur
     integer(ip)                 :: kfl_deriv_vel
     character(5)                :: wname
     integer(ip)                 :: lresp(mresp_ker)
     character(5)                :: where
  end type typ_laws

  type typ_valpr_ker
     integer(ip)                 :: kfl_exist
     integer(ip)                 :: kfl_nedsm
     integer(ip)                 :: kfl_type      ! Scalar,matrix,vector
     integer(ip)                 :: dim           ! Dimension of properties
     character(5)                :: name          ! Name of property
     integer(ip),       pointer  :: comp(:)       ! Component to be used when updating on-the-fly
     integer(ip),       pointer  :: on_the_fly(:) ! Compute property on the fly
     character(5),      pointer  :: wlaws(:)
     integer(ip),       pointer  :: ilaws(:)
     real(rp),          pointer  :: rlaws(:,:)
     type(r1p),         pointer  :: value_ielem(:)
     type(r1p),         pointer  :: value_iboun(:)
     type(r2p),         pointer  :: grval_ielem(:)
     type(r2p),         pointer  :: drval_ielem(:)
     type(r3p),         pointer  :: drval_tur_ielem(:)
     type(r3p),         pointer  :: drval_vel_ielem(:)
     type(r3p),         pointer  :: gdval_ielem(:)
     real(rp),          pointer  :: value_ipoin(:)
     real(rp),          pointer  :: grval_ipoin(:,:)
     real(rp),          pointer  :: value_const(:,:)
     real(rp),          pointer  :: value_default(:)
     type(typ_laws)              :: llaws(mlaws_ker)
     integer(ip),       pointer  :: update(:,:)
     character(len=5),  pointer  :: time_function(:)
  end type typ_valpr_ker
  !
  ! Just a list of the properties to improve readability
  !
  type lis_prope_ker   
     type(typ_valpr_ker), pointer :: prop
  end type lis_prope_ker
  
contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-11-02
  !> @brief   Elgat init
  !> @details Initialization of the class
  !> 
  !-----------------------------------------------------------------------
  
  pure subroutine init(self)
    
    class(elgat_typ), intent(inout) :: self

    self % kfl_elcod = .true.
    self % kfl_elvel = .false.
    self % kfl_elwal = .false.
    self % kfl_eltem = .false.
    self % kfl_eltur = .false.
    self % kfl_elcon = .false.
    nullify(self % elcod)
    nullify(self % elvel)
    nullify(self % elwal)
    nullify(self % eltem)
    nullify(self % eltur)
    nullify(self % elcon)
    
  end subroutine init

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-11-02
  !> @brief   Elgat allocate
  !> @details Allocatation of gather arrys 
  !> 
  !-----------------------------------------------------------------------
  
  pure subroutine alloca(self,VECTOR_DIM,pnode)
    
    class(elgat_typ), intent(inout) :: self
    integer(ip),      intent(in)    :: VECTOR_DIM
    integer(ip),      intent(in)    :: pnode
    
    allocate(self % elcod(VECTOR_DIM,ndime,pnode))
    
    if( self % kfl_elvel ) allocate(self % elvel(VECTOR_DIM,ndime,pnode))
    if( self % kfl_elwal ) allocate(self % elwal(VECTOR_DIM,pnode))
    if( self % kfl_eltem ) allocate(self % eltem(VECTOR_DIM,pnode))
    if( self % kfl_eltur ) allocate(self % eltur(VECTOR_DIM,2,pnode))
    if( self % kfl_elcon ) allocate(self % elcon(VECTOR_DIM,pnode,10))
    
  end subroutine alloca

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-11-02
  !> @brief   Elgat deallocate
  !> @details Deallocatation of gather arrys 
  !> 
  !-----------------------------------------------------------------------
  
  pure subroutine deallo(self)
    
    class(elgat_typ), intent(inout) :: self
    
    deallocate(self % elcod)
    
    if( self % kfl_elvel ) deallocate(self % elvel)
    if( self % kfl_elwal ) deallocate(self % elwal)
    if( self % kfl_eltem ) deallocate(self % eltem)
    if( self % kfl_eltur ) deallocate(self % eltur)
    if( self % kfl_elcon ) deallocate(self % elcon)

  end subroutine deallo
  
end module def_ker_proper
!> @}
