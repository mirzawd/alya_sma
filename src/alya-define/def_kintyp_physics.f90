!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Kinds_and_types
!> @{
!> @file    def_kintyp_physics.g90
!> @author  houzeaux
!> @date    2020-04-04
!> @brief   Functions
!> @details Communications
!-----------------------------------------------------------------------

module def_kintyp_physics

  use def_kintyp_basic
  use def_kintyp_domain, only : netyp
  
  !----------------------------------------------------------------------
  !
  ! atomos y especies pointers
  !
  !----------------------------------------------------------------------

  type atomo
     character(50)        :: wprob         ! archivo atomos
     integer(ip)          :: natoms        ! cantidad de atomos
     integer(ip)          :: nespec        ! cantidad especies
     integer(ip), pointer :: espe(:)       ! numero de especie
     integer(ip), pointer :: Tiespe(:)     ! tipo interno de especie
     integer(ip), pointer :: spin(:)       ! numero de spin por atomo
     integer(ip), pointer :: tipoPP(:)     ! tipo de PP
     real(rp),    pointer :: coorx(:)      ! pos x en order por atomo
     real(rp),    pointer :: coory(:)      ! pos y en order por atomo
     real(rp),    pointer :: coorz(:)      ! pos z en order por atomo
  end type atomo

  type especie
     character(50)        :: wprob1        ! name especie
     character(50)        :: wprob2        ! archivo ppseudo potencial
     integer(ip)          :: atnumber      ! numero atomico
     integer(ip)          :: ncp           ! nuemro cuantico ppal
     integer(ip)          :: ncl           ! nuemro cuantico l
     integer(ip)          :: ncm           ! nuemro cuantico m
     integer(ip)          :: nspin         ! nuemro spin
     integer(ip)          :: nlmax         ! orbital maximo
     integer(ip)          :: nlocal        ! orbital local
     real(rp)             :: valencia      ! electrones de valencia
     integer(ip)          :: nrad          ! numero de puntos que entrega el PP
     real(rp),allocatable :: radio(:)      ! coordenada radial del PP esferico
     real(rp),allocatable :: ppseu(:,:)    ! matriz dim nrad x nlmax con el PP
     real(rp),allocatable :: ppphi(:,:)    ! matriz dim nrad x nlmax con la PPPhi
     real(rp),allocatable :: nocupa(:)     ! ocupacion bien detalladita!
  end type especie

  !
  ! Chemical Species with their properties
  !
  type typ_speci
     character(5) :: name
     real(rp)     :: visco(2)
     integer(ip)  :: lawvi
     real(rp)     :: weigh
     real(rp)     :: densi(4)
     real(rp)     :: entha(2)
     real(rp)     :: prand
     real(rp)     :: lewis
     real(rp)     :: trang(5)
     real(rp)     :: cpcoe(10,4)
     real(rp)     :: activ(2)
  end type typ_speci

  !------------------------------------------------------------------------
  !
  ! IMMBOU module
  !
  !------------------------------------------------------------------------

  type ibtyp

     integer(ip)          :: npoib         ! Number of boundary nodes
     integer(ip)          :: nboib         ! Number of boundary elements
     integer(ip)          :: npoin         ! Number of nodes
     integer(ip)          :: nelem         ! Number of elements
     integer(ip)          :: kfl_insid     ! =1 if particle inside
     integer(ip)          :: kfl_typeb     ! =0 if embedded particle. -1 for volume type. >0 for body fitted (set number)
     integer(ip)          :: kfl_coupl     ! =0 if fringe. =1 of volume.
     integer(ip)          :: kfl_model     ! =0 no specific model

     real(rp),    pointer :: cooin(:,:)    ! Initial coordinates
     real(rp),    pointer :: cooib(:,:)    ! Coordinates
     real(rp),    pointer :: cooi2(:,:)    ! Coordinates
     integer(ip), pointer :: lnoib(:,:)    ! Boundaries
     integer(ip), pointer :: ltyib(:)      ! Boundary types
     integer(ip), pointer :: lninv(:)      ! Permutation for nodes
     integer(ip), pointer :: lbinv(:)      ! Permutation for boundaries

     real(rp),    pointer :: coord(:,:)    ! Coordinates
     integer(ip), pointer :: lnods(:,:)    ! Element connectivity
     integer(ip), pointer :: ltype(:)      ! Element types

     real(rp)             :: bobox(3,2)    ! Particle boundaring box
     real(rp), pointer    :: fabox(:,:,:)  ! Face boundaring boxes
     real(rp)             :: massa         ! Mass
     real(rp)             :: densi         ! Density
     real(rp)             :: volum         ! Volume
     real(rp)             :: momin(6)      ! Momentum of inertia tensor
     real(rp)             :: posgr(3)      ! Initial position of center of gravity
     !
     ! Linear motion
     !
     real(rp)             :: posil(3,2)    ! Linear position
     real(rp)             :: velol(3,2)    ! Linear velocity
     real(rp)             :: accel(3,3)    ! Linear accelaration
     real(rp)             :: force(3,2)    ! Force
     real(rp)             :: pforce(3)     ! Pressure Force
     real(rp)             :: vforce(3)     ! Viscous Force
     !
     ! Angular motion
     !
     real(rp)             :: posia(3,2)    ! Angular position
     real(rp)             :: veloa(3,2)    ! Angular velocity
     real(rp)             :: accea(3,3)    ! Angular accelaration
     real(rp)             :: rotac(3,3)    ! Rotation matrix
     real(rp)             :: torqu(3,2)    ! Torque
     real(rp)             :: ptorqu(3)     ! Pressure Torque
     real(rp)             :: vtorqu(3)     ! Viscous Torque
     real(rp)             :: quate(4,2)    ! quaternion (an alternative representation of the rotation matrix)

     real(rp),    pointer :: sabox(:,:,:)  ! Tree structure of face bounding boxes
     integer(ip), pointer :: blink(:)      ! Index of faces for tree structure
     integer(ip), pointer :: struc(:)      ! Structure use by tree search
     real(rp),    pointer :: ldist(:)      ! Structure use by tree search

     integer(ip)          :: npofb         ! Number of boundaries build with fringes nodes (fringe boundary)
     integer(ip)          :: nbofb         ! Number of boundary nodes build with fringes nodes (fringe boundary)

     integer(ip), pointer :: lnofb(:,:)    ! List of boundaries build with fringes nodes (fringe boundary)
     integer(ip), pointer :: ltyfb(:)      ! List of boundary types build with fringes nodes (fringe boundary)
     integer(ip), pointer :: lelfb(:)      ! List of father elements for the boundaries build with fringes nodes (fringe boundary)

     real(rp)             :: maxdi         ! Maximum distance in the particle measure from his center of gravity
     real(rp)             :: cotim         ! Estimated time of collision

     type(netyp), pointer :: lnele(:)      ! List of elements with a common node

     integer(ip)          :: cucon         ! Current contact
     integer(ip), pointer :: icont(:)      ! Index of contacts
     real   (rp), pointer :: ncont(:,:)    ! Normal contact
     real   (rp), pointer :: rcont(:,:)    ! Vector from the mass center to the contact point

     integer(ip)          :: npari         ! Number of int. model parameters
     integer(ip)          :: nparr         ! Number of real model parameters
     integer(ip), pointer :: ipara(:)      ! Int. parameters
     real(rp),    pointer :: rpara(:)      ! Real parameters
     !
     ! Conservation mass matrix variables
     !
     real   (rp)          :: mass1
     real   (rp)          :: mass2
  end type ibtyp

  !------------------------------------------------------------------------
  !
  ! Rigid body module  - similar to ibtyp but reduced - for rigid body inside ale -
  ! I do it this way in case I have more than 1 rigid body and also for consistency with cristo
  !
  !------------------------------------------------------------------------

  type rbtyp

     integer(ip)          :: npoib         ! Number of boundary nodes
     integer(ip)          :: nboib         ! Number of boundary elements

     integer(ip)          :: nrbse         ! Number of sets the rigid body is formed by
     integer(ip)          :: lrbse(10)     ! List of sets the rigid body is formed by

     real(rp),    pointer :: cooin(:,:)    ! Initial coordinates
     real(rp),    pointer :: cooib(:,:)    ! Coordinates
     integer(ip), pointer :: lnoib(:,:)    ! Boundaries
     integer(ip), pointer :: ltyib(:)      ! Boundary types
     integer(ip), pointer :: lninv(:)      ! Permutation for nodes
     integer(ip), pointer :: lbinv(:)      ! Permutation for boundaries

     real(rp)             :: massa         ! Mass
     real(rp)             :: densi         ! Density
     real(rp)             :: volum         ! Volume
     real(rp)             :: momin(6)      ! Momentum of inertia tensor
     real(rp)             :: posgr(3)      ! Initial position of center of gravity
     !
     ! Linear motion
     !
     real(rp)             :: posil(3,4)    ! Linear position
     real(rp)             :: velol(3,4)    ! Linear velocity
     real(rp)             :: accel(3,4)    ! Linear accelaration
     real(rp)             :: force(3,4)    ! Force
     real(rp)             :: vpfor(3,4)    ! Viscous + pressure Force
     real(rp)             :: pforce(3)     ! Pressure Force
     real(rp)             :: vforce(3)     ! Viscous Force
     !
     ! Angular motion
     !
     real(rp)             :: posia(3,4)    ! Angular position
     real(rp)             :: veloa(3,4)    ! Angular velocity
     real(rp)             :: accea(3,4)    ! Angular accelaration
     real(rp)             :: rotac(3,4)    ! Rotation matrix
     real(rp)             :: torqu(3,4)    ! Torque
     real(rp)             :: vptor(3,4)    ! Viscous + pressure Torque
     real(rp)             :: ptorqu(3)     ! Pressure Torque
     real(rp)             :: vtorqu(3)     ! Viscous Torque
     real(rp)             :: quate(4,4)    ! quaternion (an alternative representation of the rotation matrix)
     real(rp)             :: q_dot(4,4)    ! time derivative of the quaternion
     !
     ! Postprocess
     !
     real(rp)             :: pp_pf(3,10)   ! Pressure Force for post process
     real(rp)             :: pp_vf(3,10)   ! Viscous Force for post process
     real(rp)             :: pp_pt(3,10)   ! Pressure Torque for post process
     real(rp)             :: pp_vt(3,10)   ! Viscous Torque for post process

  end type rbtyp

  type ibint
     integer(ip), pointer :: lnode(:)
     real(rp),    pointer :: shapl(:)
     integer(ip)          :: limit
  end type ibint
   
end module def_kintyp_physics
!> @}
