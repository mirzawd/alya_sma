!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup NastinMatrixAssembly
!> @{
!> @file    mod_nsi_element_operations_fast.f90  - started from Alya 6954.
!> @author  Guillaume Houzeaux & Herbert Owen
!> @brief   Compute Navier-Stokes momentum residual r(u) for the explicit and
!>          semi-implicit schemes. In the case of semi-implicit scheme, the
!>          viscous term is assembled.
!>
!>          r(u)   = f - rho*(u.grad)u - div[(mu+mut)*2*eps(+u) - sig*u - SIG*u
!>
!>          where:
!>
!>          - f ..... -rho*beta*g*(T-Tr)
!>          - sig ... scalar permeability
!>          - SIG ... tensor permeability
!>          - rho ... density
!>          - mu .... viscosity
!>          - mut ... turbulent viscosity
!>
!>          For convection: non-conservative and 2 EMACS schemes
!>          For viscous term: divergence, Laplacian and full schemes
!>
!> @}
!------------------------------------------------------------------------

module mod_nsi_element_operations_fast

#include "def_vector_size.inc"
#ifdef _OPENACC
  use openacc
#endif
  use def_kintyp_basic,               only : ip,rp                                     ! in_const
  use def_domain,                     only : ndime                                     !i n_const_scalar  ! actually paarmeter if ndimepar
  use def_master,                     only : kfl_paral
  use mod_getpro,                     only : getpro_val                                ! Properties
  use def_master,                     only : rhsid                                     ! out       ! real(rp), pointer     :: rhsid(:)
  use def_master,                     only : times                                     ! out       ! timings
  use def_master,                     only : veloc                                     ! in_var    ! real(rp), pointer     :: veloc(:,:,:)
  use def_master,                     only : velom                                     ! in_var    ! real(rp), pointer     :: velom(:,:,:)
  use def_master,                     only : tempe                                     ! in_var    ! real(rp), pointer     :: tempe(:,:)
  use def_domain,                     only : coord                                     ! in_const  ! real(rp), pointer     :: coord(:,:)
  use def_domain,                     only : ltype                                     ! in_const  ! integer(ip), pointer  :: ltype(:) -- type of element
  use def_domain,                     only : lnods                                     ! in_const  ! integer(ip), pointer  :: lnods(:,:)
  use def_domain,                     only : lmate                                     ! in_const  ! integer(ip), pointer  :: lmate(:)
  use def_domain,                     only : mnode                                     ! in_const scalar
  use def_domain,                     only : elmar
  use def_domain,                     only : kfl_savda
  use def_master,                     only : solve            
  use def_domain,                     only : elmda_gpvol
  use def_domain,                     only : elmda_gpcar
  use def_domain,                     only : lorde  
  use def_elmtyp,                     only : TET04
  use def_elmgeo,                     only : element_type

  use def_nastin,                     only : rhsid_gravb                               ! out       ! real(rp), pointer     :: rhsid_gravb(:) 
  use def_nastin,                     only : dtinv_nsi                                 ! in_var   scalar
  use def_nastin,                     only : dt_rho_nsi                                ! out  ! real(rp), pointer     :: dt_rho_nsi(:)
  use def_nastin,                     only : mass_rho_nsi                              ! out  ! real(rp), pointer     :: mu_rho_nsi(:)
  use def_nastin,                     only : grnor_nsi                                 ! in_const_scalar
  use def_nastin,                     only : gravi_nsi                                 ! in_const  ! real(rp)         ::gravi_nsi(3)
  use def_nastin,                     only : kfl_stabi_nsi

  use def_nastin,                     only : kfl_cotem_nsi,bougr_nsi
  use def_nastin,                     only : boube_nsi,gravb_nsi
  use def_nastin,                     only : kfl_convection_type_nsi
  use def_nastin,                     only : NSI_CONVECTION_EMAC2
  use def_nastin,                     only : NSI_CONVECTION_NON_CONSERVATIVE
  use def_nastin,                     only : kfl_force_nsi 
  use def_nastin,                     only : lforc_material_nsi
  use def_nastin,                     only : xforc_material_nsi
  use def_nastin,                     only : fvins_nsi
  use def_nastin,                     only : kfl_anipo_nsi
  use def_nastin,                     only : NSI_ASSEMBLY_VISCOUS
  use def_nastin,                     only : NSI_SUPG
  use def_nastin,                     only : visco_nsi
  use def_nastin,                     only : ndbgs_nsi
  use def_nastin,                     only : kfl_fsgrb_nsi

  use def_nastin,                     only : prdivcor_nsi,kfl_stab_corio_nsi
  use mod_nsi_corio_stab,             only : nsi_corio_stab_element_operations
  use mod_nsi_corio_stab,             only : element_length_corio
  use mod_nsi_corio_stab,             only : nsi_element_stabilization_corio
  !
  ! Solver
  !
  use mod_solver,                     only : solver_assemble_element_matrix_vector
  use def_nastin,                     only : NSI_SOLVER_VISCOUS_TERM
  !
  ! include trucho
  !
  use def_kermod,                     only : kfl_nswel_ker
  use def_kermod,                     only : kfl_noslw_ker
  use def_parame,                     only : pi
  use def_nastin,                     only : boutr_nsi
  use def_nastin,                     only : kfl_dampi_nsi
  use def_nastin,                     only : top_r_damp_nsi
  use def_nastin,                     only : bot_r_damp_nsi
  use def_nastin,                     only : val_r_damp_nsi 
  use def_nastin,                     only : v_geo_damp_nsi 
  use def_nastin,                     only : k_dim_damp_nsi 

  implicit none

  private
  public :: nsi_element_operations_fast

contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-05-06
  !> @brief   Fast assembly of NS
  !> @details Bridge to enable inlining
  !> 
  !-----------------------------------------------------------------------

  subroutine nsi_element_operations_fast(&
       VECTOR_DIM,pnode,pgaus,list_elements,time1)

    integer(ip), intent(in)    :: VECTOR_DIM                   !< Number of nodes
    integer(ip), intent(in)    :: pnode                        !< Number of nodes
    integer(ip), intent(in)    :: pgaus                        !< Number of Gauss points
    integer(ip), intent(in)    :: list_elements(VECTOR_DIM)    !< List of elements
    real(rp),    intent(inout) :: time1(10)                    !< Timings

    if(      pnode == 4 .and. pgaus == 4 ) then
       !DIR$ FORCEINLINE  
       call nsi_element_operations_fast_all(&
            VECTOR_DIM,4_ip,4_ip,list_elements,time1)
    else if( pnode == 6 .and. pgaus == 6 ) then
       !DIR$ FORCEINLINE  
       call nsi_element_operations_fast_all(&
            VECTOR_DIM,6_ip,6_ip,list_elements,time1)
    else 
       call nsi_element_operations_fast_all(&
            VECTOR_DIM,pnode,pgaus,list_elements,time1)
    end if

  end subroutine nsi_element_operations_fast

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-05-06
  !> @brief   Fast assembly of NS
  !> @details Assembly of Navier-Stokes for CPU and GPU, with
  !>          restricted options
  !> 
  !-----------------------------------------------------------------------

#if defined PNODE_VALUE && defined PGAUS_VALUE 
  subroutine nsi_element_operations_fast_all(&
       VECTOR_DIM,qnode,qgaus,list_elements,time1)
#else
    subroutine nsi_element_operations_fast_all(&
         VECTOR_DIM,pnode,pgaus,list_elements,time1)
#endif
      !
      ! Input and output variables
      !
      integer(ip), intent(in)          :: VECTOR_DIM                                       !< Number of nodes
#if defined PNODE_VALUE && defined PGAUS_VALUE 
      integer(ip), intent(in)          :: qnode                        !< Number of nodes
      integer(ip), intent(in)          :: qgaus                        !< Number of Gauss points
      integer(ip), parameter           :: pnode = PNODE_VALUE
      integer(ip), parameter           :: pgaus = PGAUS_VALUE 
#else
      integer(ip), intent(in)          :: pnode                        !< Number of nodes
      integer(ip), intent(in)          :: pgaus                        !< Number of Gauss points
#endif
      integer(ip), intent(in)          :: list_elements(VECTOR_DIM)                        !< List of elements

      real(rp),    intent(inout)       :: time1(10)                                        ! Timings
      !
      ! Element matrices and vectors (stiffness)
      !
      real(rp)                         :: elauu(VECTOR_SIZE,pnode*ndime,pnode*ndime)        ! Auu
      real(rp)                         :: elrbu(VECTOR_SIZE,ndime,pnode)                    ! bu
      real(rp)                         :: elrbu2(VECTOR_SIZE,ndime,pnode)                   ! bu - Just gravity and Boussinesq
      real(rp)                         :: elrbp(VECTOR_SIZE,pnode)                          ! bp - for gravity in fractional step
      real(rp)                         :: eldtrho(VECTOR_SIZE,pnode)                        ! Projection of rho/dt
      real(rp)                         :: elmurho(VECTOR_SIZE,pnode)                        ! Projection of mu/rho
      !
      ! Gather
      !
      integer(ip)                      :: elwvi(VECTOR_SIZE)                                ! Wall viscosity element
      real(rp)                         :: elvel(VECTOR_SIZE,ndime,pnode)                    ! u
      real(rp)                         :: elvem(VECTOR_SIZE,ndime,pnode)                    ! um 
      real(rp)                         :: elcod(VECTOR_SIZE,ndime,pnode)                    ! x
      real(rp)                         :: elprdivcor(VECTOR_SIZE,pnode)                     ! Projection of the divergence of the coriolis term
      real(rp)                         :: eltem(VECTOR_SIZE,pnode)                          ! Temperature
      !
      ! Indices and dimensions
      !
      integer(ip)                      :: ielem,inode,ivect
      integer(ip)                      :: pelty,j,k,porde
      integer(ip)                      :: ipoin,igaus,ielem_last
      integer(ip)                      :: lnods_loc(VECTOR_SIZE,pnode)
      integer(ip)                      :: list_elements_p(VECTOR_SIZE)                      ! List of elements (always positive)
      !
      ! Gauss point values
      !
      real(rp)                         :: gpsha(pnode,pgaus)                                ! N
      real(rp)                         :: gpcar(VECTOR_SIZE,ndime,mnode,pgaus)              ! dN/dxi
      real(rp)                         :: gpvol(VECTOR_SIZE,pgaus)                          ! w*|J|, |J|
      real(rp)                         :: gpvis(VECTOR_SIZE,pgaus)                          ! Viscosity
      real(rp)                         :: gpwvi(VECTOR_SIZE,pgaus)                          ! Viscosity  for no slip wall
      real(rp)                         :: gpmut(VECTOR_SIZE,pgaus)                          ! mut
      real(rp)                         :: gpden(VECTOR_SIZE,pgaus)                          ! Porosity
      real(rp)                         :: gpapo(VECTOR_SIZE,ndime,ndime,pgaus)              ! Anisotropic porosity
      real(rp)                         :: gppor(VECTOR_SIZE,pgaus)                          ! Porosity

      real(rp)                         :: gptau(VECTOR_SIZE)                                ! Tau
      real(rp)                         :: gpper(VECTOR_SIZE,pnode,pgaus)                    ! Perturbation
      real(rp)                         :: normu(VECTOR_SIZE)                                ! |u|
      real(rp)                         :: chale(VECTOR_SIZE)                                ! h
      real(rp)                         :: hleng(VECTOR_SIZE,ndime)                          ! h(ndime)
      real(rp)                         :: hnatu

      real(rp)                         :: gptem(VECTOR_SIZE)                                ! Temperature
      real(rp)                         :: gpadv(VECTOR_SIZE,ndime,pgaus)                    ! u+u'
      real(rp)                         :: gprhs(VECTOR_SIZE,ndime,pgaus)                    ! RHS
      real(rp)                         :: gpvel(VECTOR_SIZE,ndime,pgaus)                    ! u
      real(rp)                         :: gpgve(VECTOR_SIZE,ndime,ndime,pgaus)              ! grad(u)
      real(rp)                         :: xjacm(VECTOR_SIZE,ndime,ndime)
      real(rp)                         :: xjaci(VECTOR_SIZE,ndime,ndime)
      real(rp)                         :: gpdet(VECTOR_SIZE,pgaus)

      real(rp)                         :: ugradu (VECTOR_SIZE,ndime)
      real(rp)                         :: ugradut(VECTOR_SIZE,ndime)
      real(rp)                         :: gradu2 (VECTOR_SIZE,ndime)
      real(rp)                         :: divu   (VECTOR_SIZE)

      real(rp)                         :: gpcod(1:VECTOR_SIZE,3)                            ! Axes motion

      real(rp)                         :: T_dtrho(VECTOR_SIZE) 
      real(rp)                         :: d_dtrho(VECTOR_SIZE) 
      real(rp)                         :: T_murho(VECTOR_SIZE) 
      real(rp)                         :: d_murho(VECTOR_SIZE) 
      real(rp)                         :: bulkv,conve
      real(rp)                         :: timea,timeb
      integer(ip)                      :: idime,pmate,iauxi,pgau1
      integer(ip)                      :: jdime,jnode,idofv,jdofv
      !
      ! Internal
      !
#ifdef OPENACCHHH
#define DEF_VECT ivect
#else
#define DEF_VECT 1:VECTOR_SIZE
#endif

#ifdef OPENACCHHH

#define FACT1X     fact1
#define FACT2X     fact2
#define DTINV_LOCX dtinv_loc
#define T1X        t1
#define T2X        t2
#define T3X        t3
#define DENOMX     denom

#else

#define FACT1X     fact1(1:VECTOR_SIZE)
#define FACT2X     fact2(1:VECTOR_SIZE)
#define DTINV_LOCX dtinv_loc(1:VECTOR_SIZE)
#define T1X        t1(1:VECTOR_SIZE)
#define T2X        t2(1:VECTOR_SIZE)
#define T3X        t3(1:VECTOR_SIZE)
#define DENOMX     denom(1:VECTOR_SIZE)

#endif

      real(rp)    :: FACT1X
      real(rp)    :: FACT2X
      real(rp)    :: DTINV_LOCX
      real(rp)    :: T1X
      real(rp)    :: T2X
      real(rp)    :: T3X
      real(rp)    :: DENOMX

      logical     :: using_velom  ! workaround to compile with nvfortran >= 21.9

      !--------------------------------------------------------------------
      !
      ! Gather: global to local
      !
      !--------------------------------------------------------------------

#ifndef OPENACCHHH
      call cputim(timea)
      call times(9) % ini()
#endif

      ielem      = list_elements(1)
      ielem_last = ielem
      pelty      = abs(ltype(ielem))
      porde      = lorde(pelty)
      pmate      = lmate(ielem)
      DTINV_LOCX = dtinv_nsi
      !
      ! Bulk velocity
      !
      if( fvins_nsi < 0.5_rp ) then 
         bulkv   = 0.0_rp
      else
         bulkv   = (1.0_rp-(fvins_nsi-1.0_rp)*2.0_rp/3.0_rp)
      end if
      !
      ! Convection
      !
      if( kfl_convection_type_nsi == NSI_CONVECTION_NON_CONSERVATIVE ) then
         conve = 0.0_rp
      else
         conve = 1.0_rp
      end if
      !
      ! Only compute once Cartesian derivatives and then copy for TET04
      !
      if( pelty == TET04 ) then
         pgau1 = 1
      else
         pgau1 = pgaus
      end if
      using_velom = .false.
      !
      ! LIST_ELEMENTS_P is always > 0. Same for LNODS_LOC
      ! Example:
      ! LIST_ELEMENTS(:)   = [ 2,3,56,0,0]
      ! LIST_ELEMENTS_P(:) = [ 2,3,56,56,56]
      !
      !$acc enter data create( gpsha, list_elements_p, lnods_loc) &
      !$acc copyin(list_elements, elmar, elmar(pelty)%shape)

      !$acc parallel loop gang vector default(present)
      do ivect = 1,VECTOR_SIZE                  
         ielem = abs(list_elements(ivect))
         if( ielem /= 0 ) then
            lnods_loc(ivect,1:pnode) = lnods(1:pnode,ielem)
            list_elements_p(ivect)   = list_elements(ivect)
            ielem_last               = ielem
         else
            list_elements_p(ivect)   = ielem_last
            lnods_loc(ivect,1:pnode) = lnods(1:pnode,ielem_last)
         end if
         gpsha(1:pnode,1:pgaus) = elmar(pelty) % shape(1:pnode,1:pgaus)
      end do
      !$acc end parallel loop
      !$acc update self ( gpsha, list_elements_p,lnods_loc)


      call cputim(timeb)
      time1(3) = time1(3) + timeb - timea
      call cputim(timea)

      !--------------------------------------------------------------------
      !
      ! Shape function and derivatives   !ver si gpvis_nsw tine que entrar aca
      ! added avelavv  becaise it was giving me error at run time not sure if it is correct
      ! tambien avupo_ker,elavv
      ! Aunque estoy corriendo un caso sin no slip walllaw y esas solo aparcen dentro de if las quiere igual
      !
      !--------------------------------------------------------------------
      if (associated(velom)) then
         if (.false.) elvem(1,1,1) = velom(1,1) ! workaround to compile with nvfortran
         using_velom = .true.
         !$acc enter data create(elvem)
      end if

      !$acc enter data create( elvel        , gpvol   , elauu   , elwvi   , &   
      !$acc                    eldtrho      , elcod   , gpadv   , gpvel   , &
      !$acc                    gpgve        , gpcar   , gpdet   , xjacm   , &
      !$acc                    xjaci        , T_dtrho , d_dtrho , T_murho , &
      !$acc                    d_murho      , gprhs   , elmurho , elrbu   , &
      !$acc                    elvem        , divu    , ugradu  , ugradut , &
      !$acc                    gradu2       , gpper   , gptau             ) &
      !$acc copyin(            elmar              , gpsha                 , &
      !$acc                    xforc_material_nsi , lforc_material_nsi    , &
      !$acc                    kfl_nswel_ker      , elmar(pelty)%deriv    , &
      !$acc                    elmar(pelty)%weigp , gravi_nsi             , &
      !$acc                    visco_nsi                                  )
      !
      ! This do starts here both in the openacc version and in the nonopenacc
      ! for the non opencc it ends 30  lines later,
      ! in the openacc case it covers all the subroutine.
      ! Similarly the scatter in the nonopenacc case needs a do ivect
      !    
      !$acc parallel loop gang vector default(present)
      !
      do ivect = 1,VECTOR_SIZE                    
         !
         ! Transient
         !
         do inode = 1,pnode
            ipoin                      = lnods_loc(ivect,inode)
            elvel(ivect,1:ndime,inode) = veloc(1:ndime,ipoin,1)
            elcod(ivect,1:ndime,inode) = coord(1:ndime,ipoin)
         end do
         if( using_velom ) then
            do inode = 1,pnode
               ipoin                      = lnods_loc(ivect,inode)
               elvem(ivect,1:ndime,inode) = velom(1:ndime,ipoin)
            end do
         end if
         !
         ! Coriolis stabilisation
         !
         if( kfl_stab_corio_nsi > 0_ip ) then
            do inode = 1,pnode
               ipoin                   = lnods_loc(ivect,inode)
               elprdivcor(ivect,inode) = prdivcor_nsi(ipoin)
            end do
         else
            elprdivcor(ivect,:) = 0.0_rp
         end if
         !
         ! Temperature coupling
         !
         if( kfl_cotem_nsi == 1) then
            do inode = 1,pnode
               ipoin              = lnods_loc(ivect,inode)
               eltem(ivect,inode) = tempe(ipoin,1)
            end do
         end if
         !
         ! no slip wall law - I could add flag so that i does it only on those elements where it is needed kfl_nswel_ker(ielem)
         !
         if ( kfl_noslw_ker /= 0_ip ) then
            ielem        = list_elements_p(ivect)
            elwvi(ivect) = kfl_nswel_ker(ielem)
         else          
            elwvi(ivect) = 0_ip
         end if

#ifndef OPENACCHHH
      end do

      call cputim(timeb)
      time1(1) = time1(1) + timeb - timea
      call times(9) % add()
      call times(10) % ini()
      call cputim(timea)
#endif

      !--------------------------------------------------------------------
      !
      ! Element Cartesian derivatives and Jacobian: GPCAR, GPVOL
      !
      !--------------------------------------------------------------------

      if( kfl_savda == 2 ) then
#ifndef OPENACCHHH
         do ivect = 1,VECTOR_SIZE     
#endif
            ielem = abs(list_elements(ivect))
            if( ielem > 0 ) then
               gpvol(ivect,1:pgaus)                 = elmda_gpvol(1:pgaus,ielem)
               gpcar(ivect,1:ndime,1:mnode,1:pgaus) = elmda_gpcar(1:ndime,1:mnode,1:pgaus,ielem)
            else
               gpvol(ivect,1:pgaus)                 = 0.0_rp
               gpcar(ivect,1:ndime,1:mnode,1:pgaus) = 0.0_rp
            end if
#ifndef OPENACCHHH
         end do
#endif
      else 

         if(      ndime == 2 ) then

            do igaus = 1,pgau1
               xjacm(DEF_VECT,1,1)    =  0.0_rp
               xjacm(DEF_VECT,1,2)    =  0.0_rp
               xjacm(DEF_VECT,2,1)    =  0.0_rp
               xjacm(DEF_VECT,2,2)    =  0.0_rp
               do k = 1,pnode
                  xjacm(DEF_VECT,1,1) =  xjacm(DEF_VECT,1,1) + elcod(DEF_VECT,1,k) * elmar(pelty) % deriv(1,k,igaus)
                  xjacm(DEF_VECT,1,2) =  xjacm(DEF_VECT,1,2) + elcod(DEF_VECT,1,k) * elmar(pelty) % deriv(2,k,igaus)
                  xjacm(DEF_VECT,2,1) =  xjacm(DEF_VECT,2,1) + elcod(DEF_VECT,2,k) * elmar(pelty) % deriv(1,k,igaus)
                  xjacm(DEF_VECT,2,2) =  xjacm(DEF_VECT,2,2) + elcod(DEF_VECT,2,k) * elmar(pelty) % deriv(2,k,igaus)
               end do
               gpdet(DEF_VECT,igaus)  =  xjacm(DEF_VECT,1,1) * xjacm(DEF_VECT,2,2) - xjacm(DEF_VECT,2,1) * xjacm(DEF_VECT,1,2)
               denom                  =  1.0_rp / gpdet(DEF_VECT,igaus)
               xjaci(DEF_VECT,1,1)    =  xjacm(DEF_VECT,2,2) * denom
               xjaci(DEF_VECT,2,2)    =  xjacm(DEF_VECT,1,1) * denom
               xjaci(DEF_VECT,2,1)    = -xjacm(DEF_VECT,2,1) * denom
               xjaci(DEF_VECT,1,2)    = -xjacm(DEF_VECT,1,2) * denom
               do j = 1, pnode
                  gpcar(DEF_VECT,1,j,igaus) =   xjaci(DEF_VECT,1,1) * elmar(pelty) % deriv(1,j,igaus) &
                       &                      + xjaci(DEF_VECT,2,1) * elmar(pelty) % deriv(2,j,igaus)

                  gpcar(DEF_VECT,2,j,igaus) =   xjaci(DEF_VECT,1,2) * elmar(pelty) % deriv(1,j,igaus) &
                       &                      + xjaci(DEF_VECT,2,2) * elmar(pelty) % deriv(2,j,igaus)
               end do
               gpvol(DEF_VECT,igaus) = elmar(pelty) % weigp(igaus) * gpdet(DEF_VECT,igaus)
            end do

         else if( ndime == 3 ) then

            do igaus = 1,pgau1
               xjacm(DEF_VECT,1:3,1:3)   = 0.0_rp
               do k = 1,pnode 
                  xjacm(DEF_VECT,1,1)    =  xjacm(DEF_VECT,1,1) + elcod(DEF_VECT,1,k) * elmar(pelty) % deriv(1,k,igaus)
                  xjacm(DEF_VECT,1,2)    =  xjacm(DEF_VECT,1,2) + elcod(DEF_VECT,1,k) * elmar(pelty) % deriv(2,k,igaus)
                  xjacm(DEF_VECT,1,3)    =  xjacm(DEF_VECT,1,3) + elcod(DEF_VECT,1,k) * elmar(pelty) % deriv(3,k,igaus)
                  xjacm(DEF_VECT,2,1)    =  xjacm(DEF_VECT,2,1) + elcod(DEF_VECT,2,k) * elmar(pelty) % deriv(1,k,igaus)
                  xjacm(DEF_VECT,2,2)    =  xjacm(DEF_VECT,2,2) + elcod(DEF_VECT,2,k) * elmar(pelty) % deriv(2,k,igaus)
                  xjacm(DEF_VECT,2,3)    =  xjacm(DEF_VECT,2,3) + elcod(DEF_VECT,2,k) * elmar(pelty) % deriv(3,k,igaus)
                  xjacm(DEF_VECT,3,1)    =  xjacm(DEF_VECT,3,1) + elcod(DEF_VECT,3,k) * elmar(pelty) % deriv(1,k,igaus)
                  xjacm(DEF_VECT,3,2)    =  xjacm(DEF_VECT,3,2) + elcod(DEF_VECT,3,k) * elmar(pelty) % deriv(2,k,igaus)
                  xjacm(DEF_VECT,3,3)    =  xjacm(DEF_VECT,3,3) + elcod(DEF_VECT,3,k) * elmar(pelty) % deriv(3,k,igaus)
               end do
               t1                        =  xjacm(DEF_VECT,2,2) * xjacm(DEF_VECT,3,3) - xjacm(DEF_VECT,3,2) * xjacm(DEF_VECT,2,3)
               t2                        = -xjacm(DEF_VECT,2,1) * xjacm(DEF_VECT,3,3) + xjacm(DEF_VECT,3,1) * xjacm(DEF_VECT,2,3)
               t3                        =  xjacm(DEF_VECT,2,1) * xjacm(DEF_VECT,3,2) - xjacm(DEF_VECT,3,1) * xjacm(DEF_VECT,2,2)
               gpdet(DEF_VECT,igaus)     =  xjacm(DEF_VECT,1,1) * t1 + xjacm(DEF_VECT,1,2) * t2 + xjacm(DEF_VECT,1,3) * t3
               denom                     =  1.0_rp / gpdet(DEF_VECT,igaus)
               xjaci(DEF_VECT,1,1)       =  t1 * denom
               xjaci(DEF_VECT,2,1)       =  t2 * denom
               xjaci(DEF_VECT,3,1)       =  t3 * denom
               xjaci(DEF_VECT,2,2)       =  ( xjacm(DEF_VECT,1,1) * xjacm(DEF_VECT,3,3) - xjacm(DEF_VECT,3,1) * xjacm(DEF_VECT,1,3)) * denom
               xjaci(DEF_VECT,3,2)       =  (-xjacm(DEF_VECT,1,1) * xjacm(DEF_VECT,3,2) + xjacm(DEF_VECT,1,2) * xjacm(DEF_VECT,3,1)) * denom
               xjaci(DEF_VECT,3,3)       =  ( xjacm(DEF_VECT,1,1) * xjacm(DEF_VECT,2,2) - xjacm(DEF_VECT,2,1) * xjacm(DEF_VECT,1,2)) * denom
               xjaci(DEF_VECT,1,2)       =  (-xjacm(DEF_VECT,1,2) * xjacm(DEF_VECT,3,3) + xjacm(DEF_VECT,3,2) * xjacm(DEF_VECT,1,3)) * denom
               xjaci(DEF_VECT,1,3)       =  ( xjacm(DEF_VECT,1,2) * xjacm(DEF_VECT,2,3) - xjacm(DEF_VECT,2,2) * xjacm(DEF_VECT,1,3)) * denom
               xjaci(DEF_VECT,2,3)       =  (-xjacm(DEF_VECT,1,1) * xjacm(DEF_VECT,2,3) + xjacm(DEF_VECT,2,1) * xjacm(DEF_VECT,1,3)) * denom

               do j = 1, pnode
                  gpcar(DEF_VECT,1,j,igaus) =   xjaci(DEF_VECT,1,1) * elmar(pelty) % deriv(1,j,igaus) &
                       &                      + xjaci(DEF_VECT,2,1) * elmar(pelty) % deriv(2,j,igaus) &
                       &                      + xjaci(DEF_VECT,3,1) * elmar(pelty) % deriv(3,j,igaus)

                  gpcar(DEF_VECT,2,j,igaus) =   xjaci(DEF_VECT,1,2) * elmar(pelty) % deriv(1,j,igaus) &
                       &                      + xjaci(DEF_VECT,2,2) * elmar(pelty) % deriv(2,j,igaus) &
                       &                      + xjaci(DEF_VECT,3,2) * elmar(pelty) % deriv(3,j,igaus)

                  gpcar(DEF_VECT,3,j,igaus) =   xjaci(DEF_VECT,1,3) * elmar(pelty) % deriv(1,j,igaus) &
                       &                      + xjaci(DEF_VECT,2,3) * elmar(pelty) % deriv(2,j,igaus) &
                       &                      + xjaci(DEF_VECT,3,3) * elmar(pelty) % deriv(3,j,igaus)
               end do
               gpvol(DEF_VECT,igaus) = elmar(pelty) % weigp(igaus) * gpdet(DEF_VECT,igaus)
            end do
         end if
         do igaus = pgau1+1,pgaus
            gpvol(DEF_VECT,igaus)                 = gpvol(DEF_VECT,1) 
            gpcar(DEF_VECT,1:ndime,1:pnode,igaus) = gpcar(DEF_VECT,1:ndime,1:pnode,1) 
         end do

      end if

#ifdef OPENACCHHH
   end do
#endif

   !--------------------------------------------------------------------
   !
   ! Properties 
   !
   !--------------------------------------------------------------------

   !$acc enter data create( gpden, gpvis, gpmut, gpwvi, gppor, gpapo) 

   call getpro_val('DENSI',gpden,pelty,pnode,pgaus,pmate,list_elements_p,lnods_loc,gpcar)
   call getpro_val('VISCO',gpvis,pelty,pnode,pgaus,pmate,list_elements_p,lnods_loc,gpcar)
   call getpro_val('TURBU',gpmut,pelty,pnode,pgaus,pmate,list_elements_p,lnods_loc,gpcar)
   call getpro_val('WALLV',gpwvi,pelty,pnode,pgaus,pmate,list_elements_p,lnods_loc,gpcar)
   call getpro_val('POROS',gppor,pelty,pnode,pgaus,pmate,list_elements_p,lnods_loc,gpcar)
   if( kfl_anipo_nsi == 1 ) &
        call getpro_val('ANIPO',gpapo,pelty,pnode,pgaus,pmate,list_elements_p,lnods_loc,gpcar) 

   !$acc update device( gpmut , gpwvi )

#ifdef OPENACCHHH
   !$acc parallel loop gang vector default(present)
   do ivect = 1,VECTOR_SIZE
#endif

      gpvis(DEF_VECT,:) = gpvis(DEF_VECT,:) + gpden(DEF_VECT,:) * gpmut(DEF_VECT,:) ! mu <= mu + mut

      if( kfl_noslw_ker /= 0 ) then
#ifndef OPENACCHHH
         do ivect = 1,VECTOR_SIZE
            if( elwvi(ivect) /= 0 ) then
               gpvis(ivect,:) = gpwvi(ivect,:)                                      ! Wall viscosity
            end if
         end do
#else
         if( elwvi(DEF_VECT) /= 0 ) then
            gpvis(DEF_VECT,:) = gpwvi(DEF_VECT,:)                                   ! Wall viscosity
         end if
#endif
      end if

      !----------------------------------------------------------------------
      !
      ! Gauss point values & eldtrho, elmurho
      !
      !----------------------------------------------------------------------

      gpgve(DEF_VECT,:,:,:) = 0.0_rp
      gprhs(DEF_VECT,:,:)   = 0.0_rp
      gpvel(DEF_VECT,:,:)   = 0.0_rp
      gptau(DEF_VECT)       = 0.0_rp
      gpper(DEF_VECT,:,:)   = 0.0_rp
      eldtrho(DEF_VECT,:)   = 0.0_rp
      elmurho(DEF_VECT,:)   = 0.0_rp
      !
      ! Gauss points values: u, grad(u), a=u-um, f
      !
      do igaus = 1,pgaus
         FACT2X =  gpden(DEF_VECT,igaus) * grnor_nsi
         do idime = 1,ndime
            do inode = 1,pnode
               gpvel(DEF_VECT,idime,igaus) = gpvel(DEF_VECT,idime,igaus) + elvel(DEF_VECT,idime,inode) * gpsha(inode,igaus)
               do jdime = 1,ndime
                  gpgve(DEF_VECT,jdime,idime,igaus) = gpgve(DEF_VECT,jdime,idime,igaus) &
                       + gpcar(DEF_VECT,jdime,inode,igaus) * elvel(DEF_VECT,idime,inode)
               end do
            end do
            gpadv(DEF_VECT,idime,igaus) = gpvel(DEF_VECT,idime,igaus) 
            gprhs(DEF_VECT,idime,igaus) = gprhs(DEF_VECT,idime,igaus) + FACT2X * gravi_nsi(idime)
         end do
      end do
      !
      ! Mesh velocity um
      !
      if( using_velom ) then
         do igaus = 1,pgaus
            do inode = 1,pnode
               do idime = 1,ndime
                  gpadv(DEF_VECT,idime,igaus)  = gpadv(DEF_VECT,idime,igaus) &
                       - elvem(DEF_VECT,idime,inode) * gpsha(inode,igaus)
               end do
            end do
         end do
      end if
      !
      ! Material force
      !
      if( kfl_force_nsi == 1 ) then
         if( lforc_material_nsi(pmate) == 2 ) then
            do igaus = 1,pgaus
               do idime = 1,ndime
                  gprhs(DEF_VECT,idime,igaus) = gprhs(DEF_VECT,idime,igaus) +  xforc_material_nsi(idime,pmate)
               end do
            end do
         end if
      end if


#ifndef OPENACCHHH  
      !
      ! Boussinesq coupling: -rho*beta*g*(T-Tr)
      !
      if( kfl_cotem_nsi == 1) then
         FACT1X = bougr_nsi * boube_nsi
         do igaus = 1,pgaus
            gptem(DEF_VECT) = 0.0_rp
            do inode = 1,pnode
               gptem(DEF_VECT) = gptem(DEF_VECT) + eltem(DEF_VECT,inode) * gpsha(inode,igaus)
            end do
            FACT2X = gpden(DEF_VECT,igaus) * FACT1X * ( gptem(DEF_VECT) - boutr_nsi )
            do idime = 1,ndime
               gprhs(DEF_VECT,idime,igaus) = gprhs(DEF_VECT,idime,igaus) - FACT2X * gravb_nsi(idime)
            end do
         end do
      end if
#endif
      !
      ! Tau and Tim
      !
      if( porde == 1 ) then

         do igaus = 1,pgaus
            FACT2X = gpvol(DEF_VECT,igaus) / ( gpden(DEF_VECT,igaus) * DTINV_LOCX )
            ! This was mu_rho FACT1X = gpvol(DEF_VECT,igaus) * gpvis(DEF_VECT,igaus) / gpden(DEF_VECT,igaus) 
            FACT1X = gpvol(DEF_VECT,igaus) * gpden(DEF_VECT,igaus) 
            do inode = 1,pnode
               eldtrho(DEF_VECT,inode) = eldtrho(DEF_VECT,inode) + gpsha(inode,igaus) * FACT2X
               elmurho(DEF_VECT,inode) = elmurho(DEF_VECT,inode) + gpsha(inode,igaus) * FACT1X
            end do
         end do

      else

         T_dtrho(DEF_VECT) = 0.0_rp
         d_dtrho(DEF_VECT) = 0.0_rp
         T_murho(DEF_VECT) = 0.0_rp
         d_murho(DEF_VECT) = 0.0_rp
         do igaus = 1,pgaus
            FACT2X                 = gpvol(DEF_VECT,igaus) / ( gpden(DEF_VECT,igaus) * DTINV_LOCX )
            !! AB fractional step fix FACT1X                 = gpvol(DEF_VECT,igaus) * gpvis(DEF_VECT,igaus) / gpden(DEF_VECT,igaus)
            FACT1X = gpvol(DEF_VECT,igaus) * gpden(DEF_VECT,igaus) 


            T_dtrho(DEF_VECT) = T_dtrho(DEF_VECT) + FACT2X
            T_murho(DEF_VECT) = T_murho(DEF_VECT) + FACT1X
            do inode = 1,pnode
               eldtrho(DEF_VECT,inode) = eldtrho(DEF_VECT,inode) + gpsha(inode,igaus)**2 * FACT2X
               elmurho(DEF_VECT,inode) = elmurho(DEF_VECT,inode) + gpsha(inode,igaus)**2 * FACT1X
            end do
         end do

         do inode = 1,pnode
            d_dtrho(DEF_VECT) = d_dtrho(DEF_VECT) + eldtrho(DEF_VECT,inode)
            d_murho(DEF_VECT) = d_murho(DEF_VECT) + elmurho(DEF_VECT,inode)
         end do

         do inode = 1,pnode
            eldtrho(DEF_VECT,inode) = eldtrho(DEF_VECT,inode) * T_dtrho(DEF_VECT) / d_dtrho(DEF_VECT)
            elmurho(DEF_VECT,inode) = elmurho(DEF_VECT,inode) * T_murho(DEF_VECT) / d_murho(DEF_VECT)
         end do

      end if

      !--------------------------------------------------------------------
      !
      ! Stabilization perturbation: tau * (u-um).grad v
      !
      !--------------------------------------------------------------------

      if( kfl_stabi_nsi == NSI_SUPG ) then
         hnatu = element_type(pelty) % natural_length 
         do igaus = 1,pgaus
            !
            ! h
            !
            if( ndime == 2 ) then
               hleng(DEF_VECT,1) = hnatu / sqrt(  xjaci(DEF_VECT,1,1) * xjaci(DEF_VECT,1,1) &
                    &                           + xjaci(DEF_VECT,1,2) * xjaci(DEF_VECT,1,2) )
               hleng(DEF_VECT,2) = hnatu / sqrt(  xjaci(DEF_VECT,2,1) * xjaci(DEF_VECT,2,1) &
                    &                           + xjaci(DEF_VECT,2,2) * xjaci(DEF_VECT,2,2) )
               chale(DEF_VECT)   = min(hleng(DEF_VECT,1),hleng(DEF_VECT,2))/real(porde,rp)
            else
               hleng(DEF_VECT,1) = hnatu / sqrt(  xjaci(DEF_VECT,1,1) * xjaci(DEF_VECT,1,1) &
                    &                           + xjaci(DEF_VECT,1,2) * xjaci(DEF_VECT,1,2) &
                    &                           + xjaci(DEF_VECT,1,3) * xjaci(DEF_VECT,1,3) )
               hleng(DEF_VECT,2) = hnatu / sqrt(  xjaci(DEF_VECT,2,1) * xjaci(DEF_VECT,2,1) &
                    &                           + xjaci(DEF_VECT,2,2) * xjaci(DEF_VECT,2,2) &
                    &                           + xjaci(DEF_VECT,2,3) * xjaci(DEF_VECT,2,3) )
               hleng(DEF_VECT,3) = hnatu / sqrt(  xjaci(DEF_VECT,3,1) * xjaci(DEF_VECT,3,1) &
                    &                           + xjaci(DEF_VECT,3,2) * xjaci(DEF_VECT,3,2) &
                    &                           + xjaci(DEF_VECT,3,3) * xjaci(DEF_VECT,3,3) )               
               chale(DEF_VECT)   = min(hleng(DEF_VECT,1),hleng(DEF_VECT,2),hleng(DEF_VECT,3))/real(porde,rp)               
            end if
            !
            ! |u|
            !
            normu(DEF_VECT) = 0.0_rp
            do idime = 1,ndime
               normu(DEF_VECT) = normu(DEF_VECT) + gpadv(DEF_VECT,idime,igaus)*gpadv(DEF_VECT,idime,igaus)
            end do
            normu(DEF_VECT) = sqrt(normu(DEF_VECT))
            !
            ! tau * (u-um).grad v
            !
            gptau(DEF_VECT) = 1.0_rp / ( 2.0_rp*normu(DEF_VECT)/chale(DEF_VECT) + 4.0_rp * gpvis(DEF_VECT,igaus)/(gpden(DEF_VECT,igaus)*chale(DEF_VECT)**2))
            do inode = 1,pnode
               do idime = 1,ndime 
                  gpper(DEF_VECT,inode,igaus) =  gpper(DEF_VECT,inode,igaus) + gpadv(DEF_VECT,idime,igaus) * gpcar(DEF_VECT,idime,inode,igaus) 
               end do
               gpper(DEF_VECT,inode,igaus) = gptau(DEF_VECT) * gpper(DEF_VECT,inode,igaus)
            end do
         end do
      end if

      !----------------------------------------------------------------------
      !
      ! Element matrices 
      !
      !----------------------------------------------------------------------

      elrbu(DEF_VECT,:,:) = 0.0_rp

#ifndef OPENACCHHH
      if(kfl_fsgrb_nsi==1) elrbu2(DEF_VECT,:,:)  = 0.0_rp
      if(kfl_fsgrb_nsi==1) elrbp(DEF_VECT,:)     = 0.0_rp
#endif

      do igaus = 1,pgaus
         
         ugradu(DEF_VECT,:)  = 0.0_rp
         ugradut(DEF_VECT,:) = 0.0_rp
         gradu2(DEF_VECT,:)  = 0.0_rp
         divu(DEF_VECT)      = 0.0_rp
         !
         ! (a.grad) u, a.grad(u)^t, div u, a=u-um
         !
         do idime = 1,ndime
            do jdime = 1,ndime
               ugradu (DEF_VECT,idime) = ugradu(DEF_VECT,idime) + gpadv(DEF_VECT,jdime,igaus) * gpgve(DEF_VECT,jdime,idime,igaus)
            end do
            divu(DEF_VECT) = divu(DEF_VECT) + gpgve(DEF_VECT,idime,idime,igaus)
         end do

         !--------------------------------------------------------------------
         !
         ! Convective term (using EMAC) and RHS. The EMAC scheme is:
         !
         ! 2*u.eps(u) + (div u) u - 0.5 * rho * grad(u**2)
         !
         ! EMAC2:
         ! rho (u.grad) u + rho * u.grad(u)^t + rho * (div u) u - 0.5 * rho * grad(u**2)
         !
         ! EMAC:
         ! If  0.5 * rho * grad(u**2) is substituted by rho * u.grad(u)^t, the the EMAC
         ! reads:
         ! rho (u.grad) u + rho * (div u) u 
         !
         !--------------------------------------------------------------------

         if( kfl_convection_type_nsi == NSI_CONVECTION_EMAC2 ) then

            if( using_velom ) then
               do inode = 1,pnode
                  do idime = 1,ndime
                     do jdime = 1,ndime
                        gradu2(DEF_VECT,idime) = gradu2(DEF_VECT,idime) &
                             + gpcar(DEF_VECT,idime,inode,igaus) * (elvel(DEF_VECT,jdime,inode) - elvem(DEF_VECT,jdime,inode))**2
                        ugradut(DEF_VECT,idime) = ugradut(DEF_VECT,idime) + gpadv(DEF_VECT,jdime,igaus) * gpgve(DEF_VECT,idime,jdime,igaus)
                     end do
                  end do
               end do

            else

               do inode = 1,pnode
                  do idime = 1,ndime
                     do jdime = 1,ndime
                        gradu2(DEF_VECT,idime) = gradu2(DEF_VECT,idime) &
                             + gpcar(DEF_VECT,idime,inode,igaus) * elvel(DEF_VECT,jdime,inode) * elvel(DEF_VECT,jdime,inode)
                        ugradut(DEF_VECT,idime) = ugradut(DEF_VECT,idime) + gpadv(DEF_VECT,jdime,igaus) * gpgve(DEF_VECT,idime,jdime,igaus)
                     end do
                  end do
               end do
            end if

            do inode = 1,pnode 
               FACT1X = (gpsha(inode,igaus) + gpper(DEF_VECT,inode,igaus))* gpvol(DEF_VECT,igaus)
               do idime = 1,ndime
                  !
                  ! rho * (a.grad) u + rho * a.grad(u)^t + rho * (div u) a + sig*u
                  !
                  elrbu(DEF_VECT,idime,inode) = elrbu(DEF_VECT,idime,inode) - gpden(DEF_VECT,igaus) * FACT1X * (&
                       ugradu (DEF_VECT,idime) +  ugradut(DEF_VECT,idime) + divu(DEF_VECT) * gpadv(DEF_VECT,idime,igaus) + &
                       gppor(DEF_VECT,igaus) * gpvel(DEF_VECT,idime,igaus) )                    
                  ! 
                  ! - 0.5 * grad(a**2)
                  ! 
                  elrbu(DEF_VECT,idime,inode) = elrbu(DEF_VECT,idime,inode) + 0.5_rp * gpden(DEF_VECT,igaus) * FACT1X * gradu2(DEF_VECT,idime) 
                  !
                  ! bu = ( f , v ) 
                  !
                  elrbu(DEF_VECT,idime,inode) = elrbu(DEF_VECT,idime,inode) + FACT1X * gprhs(DEF_VECT,idime,igaus)
               end do
            end do

         else

            do inode = 1,pnode 
               FACT1X = (gpsha(inode,igaus) + gpper(DEF_VECT,inode,igaus)) * gpvol(DEF_VECT,igaus)
               do idime = 1,ndime
                  !
                  ! rho * (a.grad) u + rho * (div u) a 
                  !
                  elrbu(DEF_VECT,idime,inode) = elrbu(DEF_VECT,idime,inode) - gpden(DEF_VECT,igaus) * FACT1X * (&
                       ugradu (DEF_VECT,idime) + conve * divu(DEF_VECT) * gpadv(DEF_VECT,idime,igaus) + &
                       gppor(DEF_VECT,igaus) * gpvel(DEF_VECT,idime,igaus) )                
                  !
                  ! bu = ( f , v ) 
                  !
                  elrbu(DEF_VECT,idime,inode)  = elrbu(DEF_VECT,idime,inode) + FACT1X * gprhs(DEF_VECT,idime,igaus)
               end do
            end do

         end if

      end do

      !--------------------------------------------------------------------
      !
      ! Anisotropic porosity
      !
      !--------------------------------------------------------------------

      if( kfl_anipo_nsi == 1 ) then

         do igaus = 1,pgaus
            do inode = 1,pnode
               FACT1X = (gpsha(inode,igaus) + gpper(DEF_VECT,inode,igaus)) * gpvol(DEF_VECT,igaus)               
               do idime = 1,ndime
                  do jdime = 1,ndime
                     elrbu(DEF_VECT,idime,inode) = elrbu(DEF_VECT,idime,inode) - FACT1X * gpapo(DEF_VECT,jdime,idime,igaus) * gpvel(DEF_VECT,jdime,igaus)
                  end do
               end do
            end do
         end do

      end if

#ifndef OPENACCHHH
      if (kfl_dampi_nsi == 1_ip) then

         do igaus = 1,pgaus

            gpcod = 0.0_rp
            !
            ! Here I could just define the y coordinate
            !
            do idime = 1,ndime
               do inode = 1,pnode
                  gpcod(DEF_VECT,idime) = gpcod(DEF_VECT,idime) + gpsha(inode,igaus) * elcod(DEF_VECT,idime,inode)
               end do
            end do

            do inode = 1,pnode
               FACT1X = gpvol(DEF_VECT,igaus) * gpsha(inode,igaus)  ! ( f , v )
               do idime = 1,ndime

                  elrbu(DEF_VECT,idime,inode)  = elrbu(DEF_VECT,idime,inode)  -  FACT1X  *  &
                       gpden(DEF_VECT,igaus) * ( gpvel(DEF_VECT,idime,igaus) - v_geo_damp_nsi(idime)) * val_r_damp_nsi  * &
                       0.5_rp * ( 1.0_rp - cos ( pi * max ( (gpcod(DEF_VECT,k_dim_damp_nsi ) - bot_r_damp_nsi ) / &
                       (top_r_damp_nsi - bot_r_damp_nsi ) , 0.0_rp) ) )        
               end do
            end do
         end do
      end if


      if (kfl_fsgrb_nsi == 1_ip) then

         do igaus = 1,pgaus
            do inode = 1,pnode

               FACT1X = gpvol(DEF_VECT,igaus) * gpsha(inode,igaus)  ! ( f , v )

               do idime = 1,ndime 
                  elrbu2(DEF_VECT,idime,inode) = elrbu2(DEF_VECT,idime,inode) + FACT1X * gprhs(DEF_VECT,idime,igaus)
               end do
               !
               ! bp = - dt ( f , grad(q) )  - correct treatment of gravity with fractional step - Orthogonal projection of grad(p)-f 
               !
               FACT2X = gpvol(DEF_VECT,igaus) / dtinv_nsi
               do idime = 1,ndime
                  elrbp(DEF_VECT,inode) = elrbp(DEF_VECT,inode) - FACT2X * gpcar(DEF_VECT,idime,inode,igaus) * gprhs(DEF_VECT,idime,igaus) 
               end do

            end do
         end do

      end if
#endif

      !--------------------------------------------------------------------
      !
      ! Viscous term
      !
      !--------------------------------------------------------------------

      if( NSI_ASSEMBLY_VISCOUS ) then
         
         elauu(DEF_VECT,:,:) = 0.0_rp
         do igaus = 1,pgaus
            FACT1X = gpvis(DEF_VECT,igaus) * gpvol(DEF_VECT,igaus) 
            do inode = 1,pnode
               do idime = 1,ndime
                  idofv = (inode-1)*ndime+idime
                  do jnode = 1,pnode
                     jdofv           = (jnode-1)*ndime+idime
                     elauu(DEF_VECT,idofv,jdofv) = elauu(DEF_VECT,idofv,jdofv) &
                          + FACT1X * gpcar(DEF_VECT,idime,inode,igaus) * gpcar(DEF_VECT,idime,jnode,igaus)
                  end do
               end do
            end do
         end do

         if( fvins_nsi > 0.9_rp ) then

            do igaus = 1,pgaus
               do inode = 1,pnode
                  do idime = 1,ndime
                     idofv = (inode-1)*ndime + idime
                     do jnode = 1,pnode
                        FACT1X = bulkv * gpvis(DEF_VECT,igaus) * gpvol(DEF_VECT,igaus) * gpcar(DEF_VECT,idime,jnode,igaus)
                        do jdime = 1,ndime
                           jdofv                       = (jnode-1)*ndime + jdime
                           elauu(DEF_VECT,idofv,jdofv) = elauu(DEF_VECT,idofv,jdofv) + FACT1X * gpcar(DEF_VECT,jdime,inode,igaus)
                        end do
                     end do
                  end do
               end do
            end do

         end if

      else

         do igaus = 1,pgaus
            FACT1X = gpvis(DEF_VECT,igaus) * gpvol(DEF_VECT,igaus) 
            do inode = 1,pnode 
               do jdime = 1,ndime
                  do idime = 1,ndime
                     elrbu(DEF_VECT,idime,inode) = elrbu(DEF_VECT,idime,inode) &
                          - FACT1X * gpcar(DEF_VECT,jdime,inode,igaus) *       &
                          &       ( gpgve(DEF_VECT,jdime,idime,igaus)          &
                          + bulkv * gpgve(DEF_VECT,idime,jdime,igaus) )
                  end do
               end do
            end do
         end do
      end if

      !--------------------------------------------------------------------
      !
      ! Assembly in global system
      !
      !--------------------------------------------------------------------        
      !
      ! Scatter element matrix to global one 
      ! 
#ifndef OPENACCHHH
      call cputim(timeb)
      time1(7) = time1(7) + timeb - timea
      call times(10) % add()
      call times(11) % ini()
      call cputim(timea)
      do ivect = 1,VECTOR_SIZE                        
#endif 
         ielem = list_elements(ivect)
         if ( ielem > 0 ) then
            do inode = 1,pnode
               ipoin = lnods_loc(ivect,inode)
               do idime = 1,ndime
                  iauxi = idime + (ipoin-1) * ndime 
                  !$acc atomic update
#ifdef NO_COLORING
                  !$OMP ATOMIC
#endif
                  rhsid(iauxi) = rhsid(iauxi) + elrbu(ivect,idime,inode)

#ifndef OPENACCHHH
                  if(kfl_fsgrb_nsi==1) rhsid_gravb(iauxi) = rhsid_gravb(iauxi) + elrbu2(ivect,idime,inode)
#endif

               end do


#ifndef OPENACCHHH
               if(kfl_fsgrb_nsi==1) then 
                  iauxi = ipoin + ndbgs_nsi 
                  !$acc atomic update
#ifdef NO_COLORING
                  !$OMP ATOMIC
#endif
                  rhsid(iauxi) = rhsid(iauxi) + elrbp(ivect,inode)
               end if
#endif


               !$acc atomic update
#ifdef NO_COLORING
               !$OMP ATOMIC
#endif
               dt_rho_nsi(ipoin) = dt_rho_nsi(ipoin) + eldtrho(ivect,inode)

               !$acc atomic update
#ifdef NO_COLORING
               !$OMP ATOMIC
#endif
               mass_rho_nsi(ipoin,1) = mass_rho_nsi(ipoin,1) + elmurho(ivect,inode)

            end do
         end if
      end do
      !$acc end parallel loop

      !
      ! Assemble matrix: This operation is performed on the GPU as well
      !
      if( NSI_ASSEMBLY_VISCOUS ) then
         call solver_assemble_element_matrix_vector(&
              solve(NSI_SOLVER_VISCOUS_TERM:),ndime,pnode,pnode*ndime,&
              list_elements,lnods_loc,elauu,visco_nsi)
      end if


      !$acc update self(visco_nsi)


      if (using_velom) then
         !$acc exit data delete(elvem)
      end if

      !$acc exit data delete( lnods_loc , gpvol, elauu ,                    &   
      !$acc              eldtrho , elvel     ,  gpadv , gpvel ,             &
      !$acc              gpgve   , elcod     ,  gpcar ,                     &
      !$acc              gpdet   ,                                          &
      !$acc              xjacm   , xjaci     , elwvi  ,                     &
      !$acc              T_dtrho , d_dtrho   , T_murho, d_murho,            &
      !$acc              gprhs   , elmurho   , elrbu  , list_elements,      &
      !$acc              divu, ugradu , ugradut , gradu2,                   &
      !$acc              xforc_material_nsi  , lforc_material_nsi,          &
      !$acc              gpsha , list_elements_p,                           &
      !$acc              gpden ,                                            &
      !$acc              gpvis ,                                            &
      !$acc              gpmut ,                                            &
      !$acc              gpwvi ,                                            &
      !$acc              elmar,                                             &
      !$acc              elmar(pelty)%shape ,                               &
      !$acc              elmar(pelty)%deriv ,                               &
      !$acc              elmar(pelty)%weigp, visco_nsi, elvem)

#ifndef OPENACCHHH
      call cputim(timeb)
      time1(8) = time1(8) + timeb - timea                       
      call times(11) % add()
#endif

    end subroutine nsi_element_operations_fast_all

  end module mod_nsi_element_operations_fast
  !> @}
