!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Maths
!> @{
!> @file    def_search_method.f90
!> @author  houzeaux
!> @date    2020-10-02
!> @brief   Projection method
!> @details 
!>
!-----------------------------------------------------------------------

module def_projection_method

  use def_kintyp_basic,                   only : ip,rp,lg,i1p,i2p,r1p,r2p,r3p
  use def_kintyp_mesh_basic,              only : mesh_type_basic
  use def_kintyp_mesh_basic,              only : bmsh_type_basic
  use def_kintyp_comm,                    only : comm_data_par_basic
  use def_mat_coo,                        only : mat_coo
  use mod_memory_basic,                   only : memory_alloca
  use mod_memory_basic,                   only : memory_deallo
  use mod_memory_basic,                   only : memory_size
  use mod_memory_basic,                   only : memory_resize
  use mod_memory_tools,                   only : memory_counter_ini
  use mod_memory_tools,                   only : memory_counter_end
  use mod_optional_argument,              only : optional_argument
  use mod_communications_global,          only : PAR_SUM
  use mod_comm_basic,                     only : par_interface_exchange
  use def_interpolation_method,           only : interpolation
  use mod_elmgeo,                         only : element_type
  use mod_elmgeo,                         only : elmgeo_natural_coordinates
  use mod_elmgeo,                         only : elmgeo_cartesian_derivatives
  use mod_elmgeo,                         only : elmgeo_natural_coordinates_on_boundaries
  use mod_elmgeo,                         only : elmgeo_projection_on_a_face
  use mod_bouder,                         only : bouder
  use def_kintyp_domain,                  only : elm
  use mod_elmgeo,                         only : elmgeo_jacobian
  use mod_elmgeo,                         only : elmgeo_jacobian_matrix
  use mod_elmgeo,                         only : elmgeo_jacobian_boundary
  use mod_communications,    only : PAR_BARRIER
  use mod_std

  implicit none
  private

  character(21), parameter :: vacal = 'def_projection_method'
  integer(ip),   parameter :: ELEMENT_PROJECTION  =  0
  integer(ip),   parameter :: BOUNDARY_PROJECTION =  1

  
  type :: projection
     character(LEN=:),      allocatable    :: name                           ! Name
     integer(ip)                           :: method                         ! Projection method
     type(interpolation),       pointer    :: interp                         ! Interpolation method on source 
     type(elm),                 pointer    :: elmar(:)                       ! Finite element interpolation
     type(bmsh_type_basic),     pointer    :: mesh_source                    ! Source mesh
     type(bmsh_type_basic),     pointer    :: mesh_target                    ! Target mesh
     type(comm_data_par_basic), pointer    :: comm                           ! Communication arrays
     type(mat_coo)                         :: matrix                         ! Interpolation matrix on target side
     type(mat_coo)                         :: mass_matrix                    ! Mass matrix on target side
     integer(ip),               pointer    :: perz(:)                        ! Permutation for matrix
     real(rp),                  pointer    :: xg(:,:)                        ! Gauss point coordinates
     integer(ip),               pointer    :: perm(:)                        ! Node permutation
     integer(ip),               pointer    :: invp(:)                        ! Inverse permutation
     integer(ip)                           :: nn                             ! Number of unknowns
     integer(ip)                           :: ng                             ! Number of Gauss points
     integer(ip)                           :: nz                             ! Number of non-zero
     integer(8)                            :: memor(2)                       ! Memory counter
   contains
     procedure,                pass        :: init                           ! Initialization
     procedure,                pass        :: deallo                         ! Deallocate
     procedure,                pass        :: input     
     procedure,                pass        :: integration_points     
     procedure,                pass        :: target_matrix 
     procedure,                pass        :: preprocess
     procedure,                pass        :: values_11
     generic                               :: values => &
          &                                   values_11
  end type projection

  public :: projection
  public :: ELEMENT_PROJECTION
  public :: BOUNDARY_PROJECTION
  
contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-06-16
  !> @brief   Initialization
  !> @details Initialization
  !> 
  !-----------------------------------------------------------------------

  subroutine init(self)

    class(projection), intent(inout) :: self

    call self % matrix % init()
    call self % mass_matrix % init()
    nullify(self % perz)
    nullify(self % perm)
    nullify(self % invp)
    nullify(self % xg)
    nullify(self % elmar)
    nullify(self % interp)
    nullify(self % mesh_source)
    nullify(self % mesh_target)
    nullify(self % comm)

    self % method = BOUNDARY_PROJECTION
    self % nn     = 0
    self % ng     = 0
    self % nz     = 0
    self % memor  = 0_8

  end subroutine init

  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-06-16
  !> @brief   Input
  !> @details Input data
  !> 
  !-----------------------------------------------------------------------

  subroutine input(     &  
       self,            &
       method,          &
       interp,          &
       elmar,           &
       mesh_source,     &
       mesh_target,     &
       comm             )

    class(projection),                            intent(inout) :: self
    integer(ip)                                                 :: method
    type(interpolation),                  target, intent(in)    :: interp
    type(elm),                            target, intent(in)    :: elmar(:)          !< Finite element interpolation
    type(bmsh_type_basic),                target, intent(in)    :: mesh_source       !< Source mesh
    type(bmsh_type_basic),                target, intent(in)    :: mesh_target       !< Target mesh
    class(comm_data_par_basic), optional, target, intent(in)    :: comm              !< Communication
   
    self % method      =  method
    self % interp      => interp
    self % elmar       => elmar
    self % mesh_source => mesh_source
    self % mesh_target => mesh_target
    if( present(comm) ) self % comm => comm
  end subroutine input

  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-06-16
  !> @brief   Preprocess
  !> @details Preprocess
  !>          XX ... 
  !> 
  !-----------------------------------------------------------------------

  subroutine preprocess(self,mask,MEMORY_COUNTER)

    class(projection),                      intent(inout) :: self              !< Projection
    logical(lg),         optional, pointer, intent(in)    :: mask(:)           !< Mask for boundaries
    integer(8),          optional,          intent(inout) :: MEMORY_COUNTER(2) !< Memory counters
    integer(8)                                            :: memor_loc(2)
    
    
    call memory_counter_ini(memor_loc,self % memor,MEMORY_COUNTER)
    !
    ! Check errors
    !
    if( .not. associated(self % elmar)  )      call runend('DEF_PROJECTION_METHOD: ELEMENT INTEGRATION DATA STRUCTURE ELMAR IS NOR DEFINED')
    if( .not. associated(self % interp) )      call runend('DEF_PROJECTION_METHOD: INTERPOLATIOBN STRATEGY ON SOURCE IS NOT DEFINED')
    if( .not. associated(self % mesh_source) ) call runend('DEF_PROJECTION_METHOD: SOURCE MESH IS NOT DEFINED')
    if( .not. associated(self % mesh_target) ) call runend('DEF_PROJECTION_METHOD: TARGET MESH IS NOT DEFINED')
    !
    ! Find integration points
    !
    call self % integration_points(self % mesh_target,self % elmar,mask,MEMORY_COUNTER=memor_loc)
    !
    ! Preprocess interpolation
    
    
    
    !
    call self % interp % preprocess(self % xg,self % mesh_source)
    !
    ! Compute projection matrix on target
    !
    call self % target_matrix(self % mesh_target,self % elmar,mask,MEMORY_COUNTER=memor_loc) 
    !
    ! Deallocate temporary arrays
    !
    call memory_deallo(memor_loc,'PERZ',vacal,self % perz)
    call memory_deallo(memor_loc,'XG'  ,vacal,self % xg)
    call memory_deallo(memor_loc,'PERM',vacal,self % perm)

    call memory_counter_end(memor_loc,self % memor,MEMORY_COUNTER)

  end subroutine preprocess

  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-06-16
  !> @brief   Deallocate
  !> @details Deallocate
  !> 
  !-----------------------------------------------------------------------

  subroutine deallo(self,MEMORY_COUNTER)

    class(projection),             intent(inout) :: self              !< Projection
    integer(8),          optional, intent(inout) :: MEMORY_COUNTER(2) !< Memory counters
    integer(8)                                   :: memor_loc(2)

    call memory_counter_ini(memor_loc,self % memor,MEMORY_COUNTER)

    call self % matrix % deallo(MEMORY_COUNTER=memor_loc)
    call self % mass_matrix % deallo(MEMORY_COUNTER=memor_loc)
    call memory_deallo(memor_loc,'SELF % INVP',vacal,self % invp)

    call memory_counter_end(memor_loc,self % memor,MEMORY_COUNTER)

  end subroutine deallo

  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-06-16
  !> @brief   Preprocess
  !> @details Preprocess
  !>          XX ... 
  !> 
  !-----------------------------------------------------------------------

  subroutine target_matrix(self,mesh,elmar,mask,lumped,MEMORY_COUNTER)

    class(projection),                      intent(inout) :: self              !< Projection
    class(mesh_type_basic),        target,  intent(in)    :: mesh              !< target mesh
    type(elm),                              intent(in)    :: elmar(:)          !< Finite element interpolation
    integer(8),          optional,          intent(inout) :: MEMORY_COUNTER(2) !< Memory counters
    logical(lg),         optional, pointer, intent(in)    :: mask(:)           !< Mask for boundaries
    logical(lg),         optional,          intent(in)    :: lumped            !< lumped the matrix or not
    integer(ip)                                           :: ipoin,iboun
    integer(ip)                                           :: igaub,ii
    integer(ip)                                           :: iz,nb
    integer(ip)                                           :: inodb,pblty
    integer(ip)                                           :: nn,pnodb
    integer(ip)                                           :: ig,nz,nd!,my_rank
    real(rp)                                              :: gbsur
    real(rp)                                              :: bocod(mesh % ndime,mesh % mnode)
    integer(8)                                            :: memor_loc(2)
    logical(lg)                                           :: if_mask
    real(rp),         pointer                             :: mass(:)
    class(mesh_type_basic),        pointer                :: mesh_target
    procedure(elmgeo_jacobian),    pointer                :: jacobian
!    integer(ip)                                           :: ipoin_i,ipoin_j,ii_i,ii_j,iz_i,iz_j,inodb_i,inodb_j
!    integer(ip)                                           :: inc
    logical(lg)                                           :: glumped


    call memory_counter_ini(memor_loc,self % memor,MEMORY_COUNTER)

    glumped =.true.
    if( present(lumped)) then 
      glumped = lumped
    end if
    !
    ! Mesh
    !
    if(     self % method == BOUNDARY_PROJECTION ) then
       jacobian => elmgeo_jacobian_boundary
       select type ( mesh )
       class is ( bmsh_type_basic ) ; mesh_target => mesh
       class is ( mesh_type_basic ) ; mesh_target => mesh % boundary
       end select
     else if( self % method == ELEMENT_PROJECTION ) then
       jacobian => elmgeo_jacobian_matrix
       select type ( mesh )
       class is ( bmsh_type_basic ) ; call runend('A VOLUME MESH SHOULD BE GIVEN AS INPUT')
       class is ( mesh_type_basic ) ; mesh_target => mesh 
       end select       
    end if
    !
    ! Initialization
    !

    nullify(mass)    
    if_mask = .true.
    nd      = mesh_target % ndime
    nb      = nd - 1
    nn      = self % nn
    nz      = self % nz
    call memory_alloca(memor_loc,'MASS',vacal,mass,mesh_target % npoin)
    !
    ! Compute matrix
    !
    call self % matrix % init()
    self % matrix % ndof1 = 1
    self % matrix % ndof2 = 1
    self % matrix % nrows = nn
    self % matrix % nz    = nz
    call self % matrix % alloca(MEMORY_COUNTER=memor_loc)

    call self % mass_matrix % init()
    self % mass_matrix % ndof1 = 1
    self % mass_matrix % ndof2 = 1
    self % mass_matrix % nrows = self % nn
    self % mass_matrix % nz    = self % nn * self % nn
    call self % mass_matrix % alloca(MEMORY_COUNTER=memor_loc)
    
    !
    ! Compute Ri = \int rr Ni dx
    !
    ig = 0
    iz = 0
    !my_rank =int(self%comm % RANK4,ip) 
    do iboun = 1,mesh_target % nelem

       if( present(mask) ) if_mask = mask(iboun)
       if( if_mask ) then
          pblty = mesh_target % ltype(iboun)
          pnodb = element_type(pblty) % number_nodes  

          do inodb = 1,pnodb
             ipoin                = mesh_target % lnods(inodb,iboun)
             bocod(1:nd,inodb) = mesh_target % coord(1:nd,ipoin)
          end do
          do igaub = 1,elmar(pblty) % pgaus
             call jacobian(nd,pnodb,bocod,elmar(pblty) % deriv(:,:,igaub),gbsur) 
             !call bouder(&
             !     pnodb,nd,nb,elmar(pblty)%deriv(:,:,igaub),&    ! Cartesian derivative
             !     bocod,baloc,gbsur)                             ! and Jacobian   
             gbsur = elmar(pblty) % weigp(igaub) * gbsur           
             ig    = ig + 1           
             do inodb = 1,pnodb
                ipoin                      = mesh_target % lnods(inodb,iboun)                  
                ii                         = self % perm(ipoin)
                iz                         = self % perz(ii)
                mass(ipoin)                = mass(ipoin)     + gbsur * elmar(pblty) % shape(inodb,igaub)      
                self % perz(ii)            = self % perz(ii) + 1
                self % matrix % xA(iz)     = ii
                self % matrix % yA(iz)     = ig
                self % matrix % vA(1,1,iz) = self % matrix % vA(1,1,iz) + gbsur * elmar(pblty) % shape(inodb,igaub) 
             end do    
          end do
       end if
    end do
    
    
    
    
    !inc = 0
    !do iboun = 1,mesh_target % nelem
    !  if( present(mask) ) if_mask = mask(iboun)
    !  if( if_mask ) then
    !    pblty = mesh_target % ltype(iboun)
    !    pnodb = element_type(pblty) % number_nodes  
    !    do inodb = 1,pnodb
    !      ipoin             = mesh_target % lnods(inodb,iboun)
    !      bocod(1:nd,inodb) = mesh_target % coord(1:nd,ipoin)
    !    end do
    !    do inodb_i = 1,pnodb
    !      do inodb_j = 1,pnodb!!!!!!

    !        ipoin_i                    = mesh_target % lnods(inodb_i,iboun) 
    !        ipoin_j                    = mesh_target % lnods(inodb_j,iboun)                 
    !        ii_i                       = self % perm(ipoin_i)
    !        ii_j                       = self % perm(ipoin_j)
    !        inc = (ii_i-1) * self%nn + ii_j
    !        self % mass_matrix % xA(inc)= ii_i
    !        self % mass_matrix % yA(inc)= ii_j               
    !        do igaub = 1,elmar(pblty) % pgaus
    !          call jacobian(nd,pnodb,bocod,elmar(pblty) % deriv(:,:,igaub),gbsur)    
    !          gbsur = elmar(pblty) % weigp(igaub) * gbsur                
    !          self % mass_matrix % vA(1,1,inc) = self % mass_matrix % vA(1,1,inc) + &
    !            gbsur * elmar(pblty) % shape(inodb_i,igaub) * elmar(pblty) % shape(inodb_j,igaub)             
    !         end do
    !       end do
    !     end do
    !   end if
    !end do

    
    
    
    
    
    
    
    !
    ! Parallel exchange
    !
    if( associated(self % comm) ) call par_interface_exchange(mass,self % comm)
    !
    ! Scale by lumped mass
    !


    if(glumped) then

      do iz = 1,nz
         ii                         = self % matrix % xA(iz)
         ipoin                      = self % invp(ii)
         self % matrix % vA(1,1,iz) = self % matrix % vA(1,1,iz) / mass(ipoin)
      end do

    end if

    !block
    !  use def_master
    !  do ii = 1,nn
    !     ipoin = self % invp(ii)
    !     print*,'mass=',lninv_loc(ipoin),mass(ipoin)
    !  end do
    !end block

    call memory_deallo(memor_loc,'MASS',vacal,mass)
    call memory_counter_end(memor_loc,self % memor,MEMORY_COUNTER)
  end subroutine target_matrix

  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-06-16
  !> @brief   Preprocess
  !> @details Preprocess
  !>          SELF % XG ... 
  !> 
  !-----------------------------------------------------------------------

  subroutine integration_points(self,mesh,elmar,mask,MEMORY_COUNTER)

    class(projection),                      intent(inout) :: self              !< Projection
    class(mesh_type_basic),        target,  intent(in)    :: mesh              !< target mesh
    type(elm),                              intent(in)    :: elmar(:)          !< Finite element interpolation
    logical(lg),         optional, pointer, intent(in)    :: mask(:)           !< Mask for boundaries
    integer(8),          optional,          intent(inout) :: MEMORY_COUNTER(2) !< Memory counter
    integer(ip)                                           :: ipoin,iboun
    integer(ip)                                           :: ng,igaub,ii,nb
    integer(ip)                                           :: inodb,pblty
    integer(ip)                                           :: nn,pnodb
    integer(ip)                                           :: ig,nz,nd
    integer(ip)                                           :: kk,ll
    real(rp)                                              :: gbsur
    real(rp)                                              :: bocod(mesh % ndime,mesh % mnode)
    integer(8)                                            :: memor_loc(2)
    logical(lg)                                           :: if_mask
    class(mesh_type_basic),        pointer                :: mesh_target
    procedure(elmgeo_jacobian),    pointer                :: jacobian
    !
    ! Permutation to identify nodes involved on the target surface
    !
    call memory_counter_ini(memor_loc,self % memor,MEMORY_COUNTER)
    !
    ! Mesh
    !
    if(self % method == BOUNDARY_PROJECTION ) then
       jacobian => elmgeo_jacobian_boundary
       select type ( mesh )
       class is ( bmsh_type_basic ) ; mesh_target => mesh
       class is ( mesh_type_basic ) ; mesh_target => mesh % boundary
       end select
    else if( self % method == ELEMENT_PROJECTION ) then
       jacobian => elmgeo_jacobian_matrix
       select type ( mesh )
       class is ( bmsh_type_basic ) ; call runend('A VOLUME MESH SHOULD BE GIVEN AS INPUT')
       class is ( mesh_type_basic ) ; mesh_target => mesh 
       end select       
    end if  
    if_mask = .true.
    if( present(mask) ) then
       if( associated(mask) ) then
          if( memory_size(mask) < mesh_target % nelem ) call runend('DEF_PROJECTION_METHOD: MASK SIZE IS TOO SMALL')
       end if
    end if
    call memory_alloca(memor_loc,'SELF % PERM',vacal,self % perm,mesh_target % npoin)
    nd = mesh_target % ndime
    nb = nd - 1
    ng = 0
    nn = 0

    do iboun = 1,mesh_target % nelem
       if( present(mask) ) if_mask = mask(iboun)
       if( if_mask ) then
          pblty = mesh_target % ltype(iboun)
          ng    = ng + elmar(pblty) % pgaus
          do inodb = 1,element_type(pblty) % number_nodes
             ipoin = mesh_target % lnods(inodb,iboun)
             if( self % perm(ipoin) == 0 ) then
                nn                 = nn + 1
                self % perm(ipoin) = nn
             end if
          end do
       end if
    end do
    call memory_alloca(memor_loc,'SELF % INVP',vacal,self % invp,nn)
    do ipoin = 1,memory_size(self % perm)
       ii = self % perm(ipoin)
       if( ii > 0 ) self % invp(ii) = ipoin

    end do
    !
    ! Compute matrix permutation and Gauss points coordinates
    !
    call memory_alloca(memor_loc,'SELF % XG'  ,vacal,self % xg,nd,ng)
    call memory_alloca(memor_loc,'SELF % PERZ',vacal,self % perz,nn+1_ip)
    ig = 0
    nz = 0
    do iboun = 1,mesh_target % nelem
       if( present(mask) ) if_mask = mask(iboun)

       if( if_mask ) then
          pblty = mesh_target % ltype(iboun)
          pnodb = element_type(pblty) % number_nodes
          do inodb = 1,pnodb
             ipoin             = mesh_target % lnods(inodb,iboun)
             bocod(1:nd,inodb) = mesh_target % coord(1:nd,ipoin)
          end do
          do igaub = 1,elmar(pblty) % pgaus
             call jacobian(nd,pnodb,bocod,elmar(pblty) % deriv(:,:,igaub),gbsur) 

             gbsur = elmar(pblty) % weigp(igaub) * gbsur
             ig    = ig + 1
             do inodb = 1,pnodb
                ipoin              = mesh_target % lnods(inodb,iboun)
                ii                 = self % perm(ipoin)
                nz                 = nz                 + 1
                self % perz(ii)    = self % perz(ii)    + 1
                self % xg(1:nd,ig) = self % xg(1:nd,ig) + bocod(1:nd,inodb) * elmar(pblty) % shape(inodb,igaub)  
             end do
          end do
       end if
    end do

    !
    ! Convert permutation array in matrix into a graph
    !
    kk             = self % perz(1)
    self % perz(1) = 1 
    do ii = 2,nn+1
       ll              = self % perz(ii)
       self % perz(ii) = self % perz(ii-1) + kk
       kk              = ll
    end do
    !
    ! Dimensions
    !
    self % nn = nn
    self % ng = ng
    self % nz = nz

    call memory_counter_end(memor_loc,self % memor,MEMORY_COUNTER)
  end subroutine integration_points

  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-01-18
  !> @brief   Compute value
  !> @details Compute value
  !> 
  !-----------------------------------------------------------------------

  subroutine values_11(self,xx_in,xx_out,POINT,INITIALIZATION,MEMORY_COUNTER)

    class(projection),                      intent(inout) :: self
    real(rp),                      pointer, intent(in)    :: xx_in(:)
    real(rp),                      pointer, intent(inout) :: xx_out(:)
    integer(ip),          optional,         intent(in)    :: POINT
    logical(lg),          optional,         intent(in)    :: INITIALIZATION
    integer(8),           optional,         intent(inout) :: MEMORY_COUNTER(2) !< Memory counters
    integer(ip)                                           :: ii,ipoin!,my_rank
    integer(8)                                            :: memor_loc(2)
    real(rp),                      pointer                :: xx_ng(:)
    real(rp),                      pointer                :: xx_nn(:)



    call memory_counter_ini(memor_loc,self % memor,MEMORY_COUNTER)

    nullify(xx_ng,xx_nn)
    call memory_alloca(memor_loc,'XX_NG',vacal,xx_ng,self % ng)
    call memory_alloca(memor_loc,'XX_NN',vacal,xx_nn,self % nn)
    !
    ! Initialization
    !
    if( optional_argument(.true.,INITIALIZATION) ) then
       do ii = 1,self % nn
          ipoin     = self % invp(ii)
          xx_nn(ii) = xx_out(ipoin) 
       end do
    end if
    !
    ! Source nodes (global) => Target Gauss-points
    !
    call self % interp % values(xx_in,xx_ng,POINT,INITIALIZATION,MEMORY_COUNTER=memor_loc) 
    !my_rank =int(self%comm % RANK4,ip) 
    !
    ! Target Gauss points => Target nodes (local)
    !
    call self % matrix % mV(xx_ng,xx_nn,INITIALIZATION=INITIALIZATION)
    !
    ! Permute Target nodes (local) => Target nodes (global)
    !
    do ii = 1,self % nn
       ipoin         = self % invp(ii)
       xx_out(ipoin) = xx_nn(ii)
    end do
    !
    ! Parallel exchange
    !
    if( associated(self % comm) ) call par_interface_exchange(xx_out,self % comm)
    call memory_deallo(memor_loc,'XX_NN',vacal,xx_nn)
    call memory_deallo(memor_loc,'XX_NG',vacal,xx_ng)
    call memory_counter_end(memor_loc,self % memor,MEMORY_COUNTER)

  end subroutine values_11

end module def_projection_method

