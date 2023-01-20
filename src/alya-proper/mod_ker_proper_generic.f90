!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> 
!> @author  houzeaux
!> @date    2020-11-19
!> @brief   Generic subroutines for properties
!> @details Generic module for element loop to compute properties.
!>          Only valid for properties defined on IELEM.
!>
!>          To add the gather of a global variable:
!>
!>          1. Add the element value ELVAR
!>          2. Add the flag KFL_ELVAR
!>          3. Fill in the class procedures: init, alloca, deallo
!>          4. Add the gather in ker_proper_gather_generic: ELVAR <= VARIA
!>          5. Use it in your property subroutine as you wish
!>
!-----------------------------------------------------------------------

module mod_ker_proper_generic

#include "def_vector_size.inc"
  use def_vectorization
  use def_vectorization
  use def_kintyp_basic,      only : ip,rp,lg
  use def_master,            only : kfl_paral
  use def_master,            only : veloc
  use def_master,            only : conce
  use def_master,            only : tempe
  use def_master,            only : untur 
  use def_domain,            only : walld
  use def_domain,            only : lnods
  use def_domain,            only : ltype
  use def_domain,            only : lmate
  use def_domain,            only : lpoty
  use def_domain,            only : ndime
  use def_domain,            only : coord
  use def_domain,            only : ngaus
  use def_domain,            only : elmar
  use def_domain,            only : mnode
  use def_kermod,            only : mlapa_ker
  use def_kermod,            only : typ_valpr_ker
  use mod_elmgeo,            only : element_type 
  use mod_elmgeo_vector,     only : elmgeo_cartesian_derivatives_jacobian
  use mod_elmgeo_vector,     only : elmgeo_element_characteristic_length
  use mod_parall,            only : num_subd_norace_par
  use mod_parall,            only : num_pack_norace_par
  use mod_parall,            only : list_elements_norace_par
  use mod_parall,            only : typ_list_elements_par
  use mod_parall,            only : par_omp_nelem_chunk
  use mod_optional_argument, only : optional_argument
  use def_ker_proper,        only : elgat_typ
  use def_ker_proper,        only : xxx_init
  use def_ker_proper,        only : xxx_operations

  implicit none
  private
  !
  ! Interface for gradient calculation
  !
  interface ker_proper_gradient
     module procedure ker_proper_gradient_scalar,      &
          &           ker_proper_gradient_scalar_pgaus,&
          &           ker_proper_gradient_vector,      &
          &           ker_proper_gradient_vector_pgaus
  end interface ker_proper_gradient
  
  interface ker_proper_value
     module procedure ker_proper_value_scalar,         &
          &           ker_proper_value_scalar_pgaus,   &
          &           ker_proper_value_vector,         &
          &           ker_proper_value_vector_pgaus
  end interface ker_proper_value
  
  public :: ker_proper_element_operations ! Main element loop
  public :: ker_proper_gradient           ! Gradient from element to Gauss point
  public :: ker_proper_value              ! Value from element to Gauss point
  public :: ker_proper_product            ! Product of the components of a vector
  public :: ker_proper_heaviside          ! Vectorized Heaviside
  public :: elgat_typ                     ! Element gather
  public :: ker_element_operations_value  ! Assembly
 
  !public :: ker_proper_product_nonpure    ! Product of the components of a vector without pure funct
  

#ifdef OPENACCHHH
#define DEF_VECT ivect
#define DEF_VECT_IN 1
#else
#define DEF_VECT 1:VECTOR_SIZE
#define DEF_VECT_IN 1:VECTOR_SIZE
#endif
  
contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-11-02
  !> @brief   Element assembly
  !> @details Element assembly for material PMATE
  !> 
  !-----------------------------------------------------------------------
  
  subroutine ker_proper_element_operations(prope_ker,pmate,subru_init,subru_operations,ON_BOUNDARIES)

    integer(ip),                           intent(in) :: pmate
    type(typ_valpr_ker),         pointer,  intent(in) :: prope_ker
    procedure(xxx_init)                               :: subru_init
    procedure(xxx_operations)                         :: subru_operations
    logical(lg),                 optional, intent(in) :: ON_BOUNDARIES
    integer(ip)                                       :: isubd
    integer(ip)                                       :: ipack,ielem,pnode,pgaus
    integer(ip)                                       :: num_subd,plaws,pelty
    integer(ip),                 pointer              :: num_pack(:)
    type(typ_list_elements_par), pointer              :: list_elements(:)    
    logical(lg)                                       :: if_on_boundaries
    type(elgat_typ)                                   :: elgat
#ifdef ALYA_OMPSS
    integer(ip)                                       :: num_neigh
#endif
    !
    ! Some checks
    !
    if_on_boundaries = optional_argument(.false.,ON_BOUNDARIES)
    plaws = prope_ker % ilaws(pmate)   
    if( prope_ker % llaws(plaws) % where /= 'IELEM' ) then
       call runend('MOD_KER_PROPER_TURBUL: WRONG PROPERTY TYPE')
    end if
    !
    ! Initialization
    !
    call elgat % init()    ! allocatesand initializes gather structures
    call subru_init(elgat) ! activates gather flags for your property
    !
    ! No race condition for element properties
    !
    num_subd      =  num_subd_norace_par
    num_pack      => num_pack_norace_par
    list_elements => list_elements_norace_par
    
#ifndef OPENACCHHH
#endif

    do isubd = 1,num_subd

#ifndef OPENACCHHH
#ifdef ALYA_OMPSS
       num_neigh = size(ompss_domains(isubd) % neighbours,KIND=ip)

       !-----------------------------------------------------------------------------
       !$OMP TASK         COMMUTATIVE(                                              &
       !$OMP              [ompss_domains(ompss_domains(isubd) % neighbours(jsubd)), &
       !$OMP              jsubd = 1,num_neigh] ) PRIORITY(num_neigh)                &
       !$OMP FIRSTPRIVATE ( num_neigh,jsubd,isubd )                                 &
       !$OMP SHARED       ( ompss_domains )                                         &
       !-----------------------------------------------------------------------------
#else
       !-----------------------------------------------------------------------------
       !$OMP PARALLEL DO                                                            &
       !$OMP SCHEDULE     ( DYNAMIC , par_omp_nelem_chunk )                         &
       !$OMP SHARED       ( isubd,par_omp_nelem_chunk )                             &
       !-----------------------------------------------------------------------------
#endif
       !-----------------------------------------------------------------------------
       !$OMP DEFAULT      ( SHARED )                                                &
       !$OMP FIRSTPRIVATE ( elgat )                                                 &
       !$OMP PRIVATE      ( ipack,pelty,pnode,pgaus,ielem,ivect )                         &
       !$OMP SHARED       ( list_elements,num_pack,element_type,ngaus,ltype,lmate,  &
       !$OMP                prope_ker,pmate                                         )              
       !-----------------------------------------------------------------------------
#endif

       do ipack = 1,num_pack(isubd)

          ielem = list_elements(isubd) % packs(ipack) % l(1)  ! Select first element
          pelty = ltype(ielem)
          pnode = element_type(pelty) % number_nodes          ! Number of nodes
          pgaus = ngaus(pelty)                                ! Number of Gauss points
          
          if( ltype(ielem) > 0 .and. lmate(ielem) == pmate ) & 
               call ker_element_operations_assembly(&
               size(list_elements(isubd) % packs(ipack) % l,kind=ip),&
               pelty,pnode,pgaus,pmate,list_elements(isubd) % packs(ipack) % l,&
               elgat,prope_ker,subru_operations)

       end do

#ifndef OPENACCHHH
#ifdef ALYA_OMPSS
       !$OMP END TASK
#else
       !$OMP END PARALLEL DO
#endif
#endif

    end do

#ifndef OPENACCHHH

#ifdef ALYA_OMPSS
    !$OMP  TASKWAIT
#endif
#else
    !$acc wait 
#endif
    
  end subroutine ker_proper_element_operations

  !-----------------------------------------------------------------------
  !> 
  !> @author  guillaume
  !> @date    2021-11-12
  !> @brief   Compute and assemble a property
  !> @details Compute and assemble a property:
  !>          Gather, element loop and scatter
  !> 
  !-----------------------------------------------------------------------
  
  subroutine ker_element_operations_assembly(&
       VECTOR_DIM,pelty,pnode,pgaus,pmate,list_elements,elgat,prope_ker,&
       subru_operations)

    integer(ip),                  intent(in)    :: VECTOR_DIM                !< Vector size 
    integer(ip),                  intent(in)    :: pelty                     !< Element type
    integer(ip),                  intent(in)    :: pnode                     !< Number of nodes
    integer(ip),                  intent(in)    :: pgaus                     !< Number of Gauss points
    integer(ip),                  intent(in)    :: pmate                     !< Material number
    integer(ip),                  intent(in)    :: list_elements(VECTOR_DIM) !< List of elements
    type(elgat_typ),              intent(inout) :: elgat
    type(typ_valpr_ker),          intent(in)    :: prope_ker
    procedure(xxx_operations)                   :: subru_operations
    integer(ip)                                 :: ivect
    
    integer(ip)                                 :: ielem,igaus     
    integer(ip)                                 :: plaws,ielem_last
    integer(ip)                                 :: lnods_loc(VECTOR_DIM,pnode)         ! Local connectivity
    real(rp)                                    :: gpdet(VECTOR_DIM,pgaus)             ! |J|
    real(rp)                                    :: xjaci(VECTOR_DIM,ndime,ndime,pgaus) ! J^-1
    real(rp)                                    :: gpcar(VECTOR_DIM,ndime,mnode,pgaus) ! grad(N)
    real(rp)                                    :: hleng(VECTOR_DIM,ndime)             ! h
    real(rp)                                    :: gppro(VECTOR_DIM,pgaus)             ! phi
    real(rp)                                    :: grpro(VECTOR_DIM,ndime,pgaus)       ! grad(phi)
    integer(ip)                                 :: list_elements_loc(VECTOR_DIM)       !< List of elements    
    !
    ! Local LNODS_LOC and LIST_ELEMENTS_LOC
    !
    grpro(:,:,:) = 0.0_rp

    do ivect = 1,VECTOR_DIM                      
       ielem = list_elements(ivect)
       if( ielem /= 0 ) then
          lnods_loc(ivect,1:pnode) = lnods(1:pnode,ielem)
          list_elements_loc(ivect) = ielem
          ielem_last               = ielem
       else
          lnods_loc(ivect,1:pnode) = lnods(1:pnode,ielem_last)
          list_elements_loc(ivect) = ielem_last
       end if
    end do
    !
    ! Gather variables
    !
    call elgat % alloca(VECTOR_DIM,pnode)

    call ker_proper_gather_generic(&
         VECTOR_DIM,pnode,list_elements,lnods_loc,elgat,prope_ker % comp(pmate))


#ifdef OPENACCHHH
    !$acc enter data copyin(elgat, elgat%elcod, elmar(pelty), &
    !$acc                   elmar(pelty)%deriv, elmar(pelty)%dercg)   

    !$acc enter data create(xjaci, gpcar, gpdet, hleng) 
#endif
    !
    ! GPSHA, GPCAR, HLENG
    !
    call elmgeo_cartesian_derivatives_jacobian(&
         VECTOR_DIM,ndime,mnode,pnode,pgaus,elgat % elcod,&
         elmar(pelty) % deriv,xjaci,gpcar,gpdet)   


    call elmgeo_element_characteristic_length(&
         VECTOR_DIM,ndime,pnode,elmar(pelty)%dercg,elgat % elcod,&
         hleng,element_type(pelty) % natural_length)


#ifdef OPENACCHHH
    !$acc exit data delete (elgat, elgat%elcod, elmar(pelty), &
    !$acc                   elmar(pelty)%deriv, elmar(pelty)%dercg)    
#endif
    !
    ! Output: gauss point values of gppro
    !
    call subru_operations(&
         VECTOR_DIM,list_elements_loc,pelty,pnode,pgaus,&
         elmar(pelty) % shape,gpcar,hleng,elgat,&
         prope_ker % rlaws(:,pmate),gppro,grpro)

#ifdef OPENACCHHH
    !$acc exit data delete (xjaci, gpcar, gpdet, hleng)
#endif


    !
    ! PROPE_KER: Scatter
    !
    do ivect = 1,VECTOR_DIM
       ielem = list_elements(ivect)
       if( ielem /= 0 ) then
          do igaus = 1,pgaus
             prope_ker % value_ielem(ielem) % a(igaus) = gppro(ivect,igaus)
          end do
       end if
    end do
    !do ivect = 1,VECTOR_DIM
    !   ielem = list_elements(ivect)
    !   if( ielem /= 0 ) then
    !      do igaus = 1,pgaus
    !         prope_ker % array_ielem(1:nelem,1:pgaus) = gppro(ivect,igaus)
    !      end do
    !   end if
    !end do
    plaws = prope_ker % ilaws(pmate)
    if( prope_ker % llaws(plaws) % kfl_gradi == 1 ) then
       do ivect = 1,VECTOR_DIM
          ielem = list_elements(ivect)
          if( ielem /= 0 ) then
             do igaus = 1,pgaus
                prope_ker % grval_ielem(ielem) % a(1:ndime,igaus) = grpro(ivect,1:ndime,igaus)
             end do
          end if
       end do
    end if



    call elgat % deallo()
    
  end subroutine ker_element_operations_assembly

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-11-02
  !> @brief   Element value
  !> @details Compute value on list of elements
  !> 
  !-----------------------------------------------------------------------

  subroutine ker_element_operations_value(&
       VECTOR_DIM,pelty,pnode,pgaus,pmate,list_elements,&
       lnods_loc,gpcar,prope_ker,subru_init,subru_operations,&
       gppro)

    integer(ip),                         intent(in) :: VECTOR_DIM                           !< Vector size 
    integer(ip),                         intent(in) :: pelty                                !< Element type
    integer(ip),                         intent(in) :: pnode                                !< Number of nodes
    integer(ip),                         intent(in) :: pgaus                                !< Number of Gauss points
    integer(ip),                         intent(in) :: pmate                                !< Material number
    integer(ip),                         intent(in) :: list_elements(VECTOR_DIM)            !< List of elements
    integer(ip),                         intent(in) :: lnods_loc(VECTOR_DIM,pnode)          !< Connectivity
    real(rp),                            intent(in) :: gpcar(VECTOR_DIM,ndime,mnode,pgaus)  !< grad(N)
    type(typ_valpr_ker),                 intent(in) :: prope_ker
    procedure(xxx_operations)                       :: subru_operations
    procedure(xxx_init)                             :: subru_init
    real(rp),                           intent(out) :: gppro(VECTOR_DIM,pgaus)              !< phi
    
    real(rp)                                        :: hleng(VECTOR_DIM,ndime)              ! h
    real(rp)                                        :: grpro(VECTOR_DIM,ndime,pgaus)        ! grad(phi)
    type(elgat_typ)                                 :: elgat
    ! 
    ! Initialization
    !
    call elgat % init()
    call subru_init(elgat)
    !
    ! Gather variables
    !
    call elgat % alloca(VECTOR_DIM,pnode)

    call ker_proper_gather_generic(&
         VECTOR_DIM,pnode,list_elements,lnods_loc,elgat,&
         prope_ker % comp(pmate))

#ifdef OPENACCHHH
    !$acc enter data copyin(elgat, elgat%elcod, elmar(pelty), &
    !$acc                   elmar(pelty)%deriv, elmar(pelty)%dercg)   
    !$acc enter data create(hleng) 
#endif
    !
    ! HLENG
    !
    call elmgeo_element_characteristic_length(&
         VECTOR_DIM,ndime,pnode,elmar(pelty)%dercg,elgat % elcod,&
         hleng,element_type(pelty) % natural_length)
    

#ifdef OPENACCHHH
    !$acc exit data delete (elgat, elgat%elcod, elmar(pelty), &
    !$acc                   elmar(pelty)%deriv, elmar(pelty)%dercg)    
#endif
    !
    ! Output: gauss point values of gppro
    !
    call subru_operations(&
         VECTOR_DIM,list_elements,pelty,pnode,pgaus,&
         elmar(pelty) % shape,gpcar,hleng,elgat,&
         prope_ker % rlaws(:,pmate),gppro,grpro)

#ifdef OPENACCHHH
    !$acc exit data delete (hleng)
#endif

    call elgat % deallo()
    
  end subroutine ker_element_operations_value

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-11-20
  !> @brief   Gather
  !> @details Generic gather
  !> 
  !-----------------------------------------------------------------------
  
  pure subroutine ker_proper_gather_generic(VECTOR_DIM,pnode,list_elements,lnods_loc,elgat,icomp)
    integer(ip),     intent(in)    :: VECTOR_DIM
    integer(ip),     intent(in)    :: pnode
    integer(ip),     intent(in)    :: list_elements(VECTOR_DIM) !< List of elements
    integer(ip),     intent(in)    :: lnods_loc(VECTOR_DIM,pnode)
    type(elgat_typ), intent(inout) :: elgat
    integer(ip),     intent(in)    :: icomp
    integer(ip)                    :: ielem,inode,ipoin
    integer(ip)                    :: ivect
    
    do ivect = 1,VECTOR_DIM                      
       ielem = list_elements(ivect)
       do inode = 1,pnode
          ipoin                              = lnods_loc(ivect,inode)
          elgat % elcod(ivect,1:ndime,inode) = coord(1:ndime,ipoin)
       end do
    end do
    if( elgat % kfl_elvel ) then
       do ivect = 1,VECTOR_DIM                      
          ielem = list_elements(ivect)
          do inode = 1,pnode
             ipoin                              = lnods_loc(ivect,inode)
             elgat % elvel(ivect,1:ndime,inode) = veloc(1:ndime,ipoin,icomp)
          end do
       end do
    end if
    if( elgat % kfl_elwal ) then
       do ivect = 1,VECTOR_DIM                      
          ielem = list_elements(ivect)
          do inode = 1,pnode
             ipoin                      = lnods_loc(ivect,inode)
             elgat % elwal(ivect,inode) = walld(ipoin)
          end do
       end do
    end if
    if( elgat % kfl_eltem ) then
       if (associated(tempe)) then
          do ivect = 1,VECTOR_DIM                      
             ielem = list_elements(ivect)
             do inode = 1,pnode
                ipoin                      = lnods_loc(ivect,inode)
                elgat % eltem(ivect,inode) = tempe(ipoin,1)
             end do
          end do
       else
          elgat % eltem = 300.0_rp
       end if
    end if
    if( elgat % kfl_eltur ) then
       do ivect = 1,VECTOR_DIM                      
          ielem = list_elements(ivect)
          do inode = 1,pnode
             ipoin                      = lnods_loc(ivect,inode)
             elgat % eltur(ivect,1,inode) = untur(1,ipoin,1) !tke
          end do
       end do
    end if
    
  end subroutine ker_proper_gather_generic

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-11-20
  !> @brief   Gradient 
  !> @details Gradient of scalar and vectors at a specific Gauss point
  !> 
  !-----------------------------------------------------------------------
  
!@  pure function ker_proper_gradient_scalar(VECTOR_DIM,pnode,gpcar,elval) result(grval)
!@
!@    integer(ip), intent(in) :: VECTOR_DIM
!@    integer(ip), intent(in) :: pnode
!@    real(rp),    intent(in) :: gpcar(VECTOR_DIM,ndime,pnode)
!@    real(rp),    intent(in) :: elval(VECTOR_DIM,pnode)
!@    real(rp)                :: grval(VECTOR_DIM,ndime)
!@    integer(ip)             :: inode,idime
!@    
!@    grval(DEF_VECT,:) = 0.0_rp
!@    do inode = 1,pnode
!@       do idime = 1,ndime
!@          grval(DEF_VECT,idime) = grval(DEF_VECT,idime) &
!@               + elval(DEF_VECT,inode) * gpcar(DEF_VECT,idime,inode)
!@       end do
!@    end do
!@        
!@  end function ker_proper_gradient_scalar
  
!@  pure function ker_proper_gradient_scalar_pgaus(VECTOR_DIM,pnode,pgaus,gpcar,elval) result(grval)
!@
!@    integer(ip), intent(in) :: VECTOR_DIM
!@    integer(ip), intent(in) :: pnode
!@    integer(ip), intent(in) :: pgaus
!@    real(rp),    intent(in) :: gpcar(VECTOR_DIM,ndime,pnode,pgaus)
!@    real(rp),    intent(in) :: elval(VECTOR_DIM,pnode)
!@    real(rp)                :: grval(VECTOR_DIM,ndime,pgaus)
!@    integer(ip)             :: inode,idime,igaus
!@
!@    do igaus = 1,pgaus
!@       grval(DEF_VECT,:,igaus) = 0.0_rp
!@       do inode = 1,pnode
!@          do idime = 1,ndime
!@             grval(DEF_VECT,idime,igaus) = grval(DEF_VECT,idime,igaus) &
!@                  + elval(DEF_VECT,inode) * gpcar(DEF_VECT,idime,inode,igaus)
!@          end do
!@       end do
!@    end do
!@    
!@  end function ker_proper_gradient_scalar_pgaus
!@  
!@  pure function ker_proper_gradient_vector(VECTOR_DIM,pnode,gpcar,elval) result(grval)
!@
!@    integer(ip), intent(in) :: VECTOR_DIM
!@    integer(ip), intent(in) :: pnode
!@    real(rp),    intent(in) :: gpcar(VECTOR_DIM,ndime,pnode)
!@    real(rp),    intent(in) :: elval(VECTOR_DIM,ndime,pnode)
!@    real(rp)                :: grval(VECTOR_DIM,ndime,ndime)
!@    integer(ip)             :: inode,idime,jdime
!@    
!@    grval(DEF_VECT,:,:) = 0.0_rp
!@    do inode = 1,pnode
!@       do idime = 1,ndime
!@          do jdime = 1,ndime
!@             grval(DEF_VECT,jdime,idime) = grval(DEF_VECT,jdime,idime) &
!@                  + elval(DEF_VECT,idime,inode) * gpcar(DEF_VECT,jdime,inode)
!@          end do
!@       end do
!@    end do
!@        
!@  end function ker_proper_gradient_vector
!@  
!@  pure function ker_proper_gradient_vector_pgaus(VECTOR_DIM,pnode,pgaus,gpcar,elval) result(grval)
!@
!@    integer(ip), intent(in) :: VECTOR_DIM
!@    integer(ip), intent(in) :: pnode
!@    integer(ip), intent(in) :: pgaus
!@    real(rp),    intent(in) :: gpcar(VECTOR_DIM,ndime,pnode,pgaus)
!@    real(rp),    intent(in) :: elval(VECTOR_DIM,ndime,pnode)
!@    real(rp)                :: grval(VECTOR_DIM,ndime,ndime,pgaus)
!@    integer(ip)             :: inode,idime,jdime,igaus
!@
!@    do igaus = 1,pgaus
!@       grval(DEF_VECT,:,:,igaus) = 0.0_rp
!@       do inode = 1,pnode
!@          do idime = 1,ndime
!@             do jdime = 1,ndime
!@                grval(DEF_VECT,jdime,idime,igaus) = grval(DEF_VECT,jdime,idime,igaus) &
!@                     + elval(DEF_VECT,idime,inode) * gpcar(DEF_VECT,jdime,inode,igaus)
!@             end do
!@          end do
!@       end do
!@    end do
!@    
!@  end function ker_proper_gradient_vector_pgaus

  pure subroutine ker_proper_gradient_scalar(VECTOR_DIM,pnode,gpcar,elval,grval)
#if defined(OPENACCHHH) && defined(NDIMEPAR)
   !$acc routine seq
#endif

    integer(ip), intent(in)    :: VECTOR_DIM
    integer(ip), intent(in)    :: pnode
    real(rp),    intent(in)    :: gpcar(VECTOR_DIM,ndime,pnode)
    real(rp),    intent(in)    :: elval(VECTOR_DIM,pnode)
    real(rp),    intent(inout) :: grval(VECTOR_DIM,ndime)
    integer(ip)                :: inode,idime
    
    grval(DEF_VECT_IN,:) = 0.0_rp
    do inode = 1,pnode
       do idime = 1,ndime
          grval(DEF_VECT_IN,idime) = grval(DEF_VECT_IN,idime) &
               + elval(DEF_VECT_IN,inode) * gpcar(DEF_VECT_IN,idime,inode)
       end do
    end do
        
  end subroutine ker_proper_gradient_scalar
 
  pure subroutine ker_proper_gradient_scalar_pgaus(VECTOR_DIM,gnode,pnode,pgaus,gpcar,elval,grval)
#if defined(OPENACCHHH) && defined(NDIMEPAR)
   !$acc routine seq
#endif

    integer(ip), intent(in)    :: VECTOR_DIM
    integer(ip), intent(in)    :: gnode
    integer(ip), intent(in)    :: pnode
    integer(ip), intent(in)    :: pgaus
    real(rp),    intent(in)    :: gpcar(VECTOR_DIM,ndime,gnode,pgaus)
    real(rp),    intent(in)    :: elval(VECTOR_DIM,pnode)
    real(rp),    intent(inout) :: grval(VECTOR_DIM,ndime,pgaus)
    integer(ip)                :: inode,idime,igaus

    do igaus = 1,pgaus
       grval(DEF_VECT_IN,:,igaus) = 0.0_rp
       do inode = 1,pnode
          do idime = 1,ndime
             grval(DEF_VECT_IN,idime,igaus) = grval(DEF_VECT_IN,idime,igaus) &
                  + elval(DEF_VECT_IN,inode) * gpcar(DEF_VECT_IN,idime,inode,igaus)
          end do
       end do
    end do
    
  end subroutine ker_proper_gradient_scalar_pgaus
 
  pure subroutine ker_proper_gradient_vector(VECTOR_DIM,pnode,gpcar,elval,grval)
#if defined(OPENACCHHH) && defined(NDIMEPAR)
   !$acc routine(ker_proper_gradient_vector) seq
#endif

    integer(ip), value, intent(in)    :: VECTOR_DIM
    integer(ip), value, intent(in)    :: pnode
    real(rp),    intent(in)    :: gpcar(VECTOR_DIM,ndime,mnode)
    real(rp),    intent(in)    :: elval(VECTOR_DIM,ndime,pnode)
    real(rp),    intent(inout) :: grval(VECTOR_DIM,ndime,ndime)
    integer(ip)                :: inode,idime,jdime
  
    grval(DEF_VECT_IN,:,:) = 0.0_rp
    do inode = 1,pnode
       do idime = 1,ndime
          do jdime = 1,ndime
             grval(DEF_VECT_IN,jdime,idime) =  grval(DEF_VECT_IN,jdime,idime) &
                  + elval(DEF_VECT_IN,idime,inode) * gpcar(DEF_VECT_IN,jdime,inode)
          end do
       end do
    end do
         
  end subroutine ker_proper_gradient_vector
 
  pure subroutine ker_proper_gradient_vector_pgaus(VECTOR_DIM,gnode,pnode,pgaus,gpcar,elval,grval)
#if defined(OPENACCHHH) && defined(NDIMEPAR)
   !$acc routine seq
#endif

    integer(ip), intent(in)    :: VECTOR_DIM
    integer(ip), intent(in)    :: gnode
    integer(ip), intent(in)    :: pnode
    integer(ip), intent(in)    :: pgaus
    real(rp),    intent(in)    :: gpcar(VECTOR_DIM,ndime,gnode,pgaus)
    real(rp),    intent(in)    :: elval(VECTOR_DIM,ndime,pnode)
    real(rp),    intent(inout) :: grval(VECTOR_DIM,ndime,ndime,pgaus)
    integer(ip)                :: inode,idime,jdime,igaus

    do igaus = 1,pgaus
       grval(DEF_VECT_IN,:,:,igaus) = 0.0_rp
       do inode = 1,pnode
          do idime = 1,ndime
             do jdime = 1,ndime
                grval(DEF_VECT_IN,jdime,idime,igaus) = grval(DEF_VECT_IN,jdime,idime,igaus) &
                     + elval(DEF_VECT_IN,idime,inode) * gpcar(DEF_VECT_IN,jdime,inode,igaus)
             end do
          end do
       end do
    end do
    
  end subroutine ker_proper_gradient_vector_pgaus
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-11-20
  !> @brief   Value 
  !> @details Value of scalar and vectors at a specific Gauss point
  !> 
  !-----------------------------------------------------------------------
  
  pure function ker_proper_value_scalar(VECTOR_DIM,pnode,gpsha,elval) result(gpval)

    integer(ip), intent(in) :: VECTOR_DIM
    integer(ip), intent(in) :: pnode
    real(rp),    intent(in) :: gpsha(pnode)
    real(rp),    intent(in) :: elval(VECTOR_DIM,pnode)
    real(rp)                :: gpval(VECTOR_DIM)
    integer(ip)             :: inode
    
    gpval(DEF_VECT) = 0.0_rp
    do inode = 1,pnode
       gpval(DEF_VECT) = gpval(DEF_VECT) &
            + elval(DEF_VECT,inode) * gpsha(inode)
    end do
        
  end function ker_proper_value_scalar

  pure function ker_proper_value_scalar_pgaus(VECTOR_DIM,pnode,pgaus,gpsha,elval) result(gpval)

    integer(ip), intent(in) :: VECTOR_DIM
    integer(ip), intent(in) :: pnode
    integer(ip), intent(in) :: pgaus
    real(rp),    intent(in) :: gpsha(pnode,pgaus)
    real(rp),    intent(in) :: elval(VECTOR_DIM,pnode)
    real(rp)                :: gpval(VECTOR_DIM,pgaus)
    integer(ip)             :: inode,igaus

    do igaus = 1,pgaus
       gpval(DEF_VECT,igaus) = 0.0_rp
       do inode = 1,pnode
          gpval(DEF_VECT,igaus) = gpval(DEF_VECT,igaus) &
               + elval(DEF_VECT,inode) * gpsha(inode,igaus)
       end do
    end do
    
  end function ker_proper_value_scalar_pgaus

  pure function ker_proper_value_vector(VECTOR_DIM,pnode,gpsha,elval) result(gpval)

    integer(ip), intent(in) :: VECTOR_DIM
    integer(ip), intent(in) :: pnode
    real(rp),    intent(in) :: gpsha(pnode)
    real(rp),    intent(in) :: elval(VECTOR_DIM,ndime,pnode)
    real(rp)                :: gpval(VECTOR_DIM,ndime)
    integer(ip)             :: inode,idime
    
    gpval(DEF_VECT,:) = 0.0_rp
    do inode = 1,pnode
       do idime = 1,ndime
          gpval(DEF_VECT,idime) = gpval(DEF_VECT,idime) &
               + elval(DEF_VECT,idime,inode) * gpsha(inode)
       end do
    end do
        
  end function ker_proper_value_vector

  pure function ker_proper_value_vector_pgaus(VECTOR_DIM,pnode,pgaus,gpsha,elval) result(gpval)

    integer(ip), intent(in) :: VECTOR_DIM
    integer(ip), intent(in) :: pnode
    integer(ip), intent(in) :: pgaus
    real(rp),    intent(in) :: gpsha(pnode,pgaus)
    real(rp),    intent(in) :: elval(VECTOR_DIM,ndime,pnode)
    real(rp)                :: gpval(VECTOR_DIM,ndime,pgaus)
    integer(ip)             :: inode,idime,igaus

    do igaus = 1,pgaus
       gpval(DEF_VECT,:,igaus) = 0.0_rp
       do inode = 1,pnode
          do idime = 1,ndime
             gpval(DEF_VECT,idime,igaus) = gpval(DEF_VECT,idime,igaus) &
                  + elval(DEF_VECT,idime,inode) * gpsha(inode,igaus)
          end do
       end do
    end do
    
  end function ker_proper_value_vector_pgaus

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-11-20
  !> @brief   Product
  !> @details Compute the product of the components of a vector
  !> 
  !-----------------------------------------------------------------------
!@@   subroutine ker_proper_product(VECTOR_DIM,nn,xx,xx_prod) 
!@@#ifdef OPENACCHHH
!@@   !$acc routine seq
!@@#endif    
!@@    integer(ip), intent(in)    :: VECTOR_DIM
!@@    integer(ip), intent(in)    :: nn
!@@    real(rp),    intent(in)    :: xx(VECTOR_DIM,nn)
!@@    real(rp),    intent(inout) :: xx_prod(VECTOR_DIM)
!@@    integer(ip)                :: ii
!@@    xx_prod(DEF_VECT_IN) = 1.0_rp
!@@    do ii = 1,nn
!@@       xx_prod(DEF_VECT_IN) = xx_prod(DEF_VECT_IN) * xx(DEF_VECT_IN,ii)
!@@    end do
!@@  end subroutine ker_proper_product

  pure function ker_proper_product(VECTOR_DIM,nn,ivect,xx) result(xx_prod) 
#ifdef OPENACCHHH
   !$acc routine seq
#endif    

    integer(ip), intent(in)  :: VECTOR_DIM
    integer(ip), intent(in)  :: nn
    integer(ip), intent(in)  :: ivect
    real(rp),    intent(in)  :: xx(VECTOR_DIM,nn)
    real(rp)                 :: xx_prod(VECTOR_DIM)
    integer(ip)              :: ii

    xx_prod(DEF_VECT) = 1.0_rp
    do ii = 1,nn
       xx_prod(DEF_VECT) = xx_prod(DEF_VECT) * xx(DEF_VECT,ii)
    end do
 
!!    xx_prod(:) = 1.0_rp
!!    do ii = 1,nn
!!       xx_prod(:) = xx_prod(:) * xx(:,ii)
!!    end do
    
  end function ker_proper_product


  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-11-20
  !> @brief   Heavside
  !> @details Vectorized Heaviside H(XX-XX_0)
  !> 
  !-----------------------------------------------------------------------
  
!!  pure function ker_proper_heaviside(VECTOR_DIM,ivect,xx,xx_0) result(H)
!!#ifdef OPENACCHHH
!!   !$acc routine seq
!!#endif    
!!   
!!    integer(ip), intent(in)  :: VECTOR_DIM
!!    integer(ip), intent(in)  :: ivect
!!    real(rp),    intent(in)  :: xx(VECTOR_DIM)
!!    real(rp),    intent(in)  :: xx_0
!!    real(rp)                 :: H(VECTOR_DIM)
!!
!!    H(DEF_VECT) = 0.5_rp*(1.0_rp+sign(1.0_rp,xx(DEF_VECT)-xx_0))
!!
!!  end function ker_proper_heaviside
  
  pure subroutine ker_proper_heaviside(VECTOR_DIM,xx,xx_0,H)
#ifdef OPENACCHHH
   !$acc routine seq
#endif    
    integer(ip), intent(in)    :: VECTOR_DIM
    real(rp),    intent(in)    :: xx(VECTOR_DIM)
    real(rp),    intent(in)    :: xx_0
    real(rp),    intent(inout) :: H(VECTOR_DIM)

    H(DEF_VECT_IN) = 0.5_rp*(1.0_rp+sign(1.0_rp,xx(DEF_VECT_IN)-xx_0))

  end subroutine ker_proper_heaviside
 
end module mod_ker_proper_generic
!> @}
