!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!





!- 78   gpcar sin pgaus;  gpvol xmile gpgve tampoco  xnutu, Aalph, Bbeta  G__ij

!------------------------------------------------------------------------
!> @addtogroup NastinMatrixAssembly
!> @{
!> @file    mod_nsi_element_operations_fast.f90  - started from Alya 6954.
!> @author  Guillaume Houzeaux & Herbert Owen
!> @brief   Navier-Stokes system element assembly only of Auu part and send it to RHS
!>          The matrices are created at the elemental level but then they are multiplid by u_n and sent to RHS (explicit only).
!>          Moreover, boundary conditions will be applied after calling to the solver as is usually done in the explicit case. MISSING hhh
!>          For the moment I leave it as it is in mod_nsi_element_operations - when guillaume uploads dirichlet algorithmic I will change it
!> @details Elemental operations. Now there are no tasks.
!>
!>          It is only valid for LES (no RANS), no thermal coupling, emac , divergence form , force just gravity , no stablization, 3d
!>          We only have calls to:
!>                 1) element_shape_function_derivatives_jacobian - could be eliminated if I leave only tetras
!>                 2) cputim
!>                 3) ker_proper - could also be eliminated
!>
!>          OTHER SIMPLICATIONS:
!>
!>         1) When I inlined nsi_rhodt_rhotau_nu_vector   I only included the part for porde==1  (linear elements)   !! added in nsi_outerr
!>                    I am still missing add an error in nsi_outerr  if there are more than linear elements
!>                    I did this to avoid an if and to avid analysing whta is done for high order elements.
!>
!>         2) local time step does not work --  I set dtinv_loc(1:VECTOR_SIZE) = dtinv_nsi added
!>
!>         3) This only works for kfl_matdi_nsi ==2 (NSI_DIRICHLET_ALGORITHM) -added
!>
!>         4) added runened if nor diverg and emac
!>
!>
!>
!>          \verbatim
!>
!>          1 ........ Element calculations and assembly of global system:
!>                     b <= b^(e) - A^(e)*u^(e): RHS ................ RHSID
!>
!>          \endverbatim
!>
!>          CORRESPONDANCE OLD TO NEW SUBROUTINES: hhh delete
!>          --------------------------------------
!>
!>          nsi_elmma4         <= nsi_element_assembly_split_oss
!>          nsi_elmga3         <= nsi_element_operations_gather
!>          nsi_elmlen         <= elmgeo_element_characteristic_length
!>          elmchl             <= elmgeo_element_length
!>          elmca2             <= element_shape_function_derivatives_jacobian 
!>          nsi_elmtss         <= nsi_element_time_step
!>          nsi_elmres         <= nsi_element_residual
!>          nsi_updsgs         <= nsi_element_subgrid_scale
!>          nsi_elmsgs         <= nsi_element_stabilization
!>          nsi_elmort         <= nsi_assembly_projections
!>          nsi_elmope_omp     <= nsi_element_operations
!>          nsi_elmmat         <= nsi_element_assembly_asgs_oss
!>                                nsi_element_assembly_asgs_oss_old
!>          nsi_elmdi3         <= nsi_element_dirichlet
!>                                nsi_element_schur
!>          nsi_assemble_schur <= nsi_assembly_schur_method
!>          nsi_elmext         <= nsi_element_extension
!>          nsi_elmexa         <= nsi_element_manufactured_solution
!>          nsi_elmexf         <= nsi_element_external_force
!>
!>
!> @}
!------------------------------------------------------------------------

module mod_nsi_element_operations_hh78

!hhh  use def_kintyp_basic,               only : ip,rp         !in_const
!hhh#ifndef VECTOR_SIZE
!hhh  use def_master,                     only : VECTOR_SIZE   !in_const
!hhh#endif
!hhh  use def_domain,                     only : ndime         !in_const_scalar  ! actually paarmeter if ndimepar
!hhh#ifndef SUPER_FAST
!hhh  use mod_element_integration,        only : element_shape_function_derivatives_jacobian
!hhh#endif
!hhh
!hhh#ifdef _OPENACC
!hhh  use openacc
!hhh#endif
!hhh
!hhh  use def_master,                     only : kfl_paral
!hhh  
!hhh  implicit none
!hhh  
!hhh  private
!hhh  public :: nsi_element_operations_hh78
!hhh
!hhhcontains
!hhh
!hhh  !-----------------------------------------------------------------------
!hhh  !> 
!hhh  !> @author  houzeaux
!hhh  !> @date    2020-05-06
!hhh  !> @brief   Fast assembly of NS
!hhh  !> @details Assembly of Navier-Stokes for CPU and GPU, with
!hhh  !>          restricted options
!hhh  !> 
!hhh  !-----------------------------------------------------------------------
!hhh
!hhh#if defined PNODE_VALUE && defined PGAUS_VALUE 
!hhh  subroutine nsi_element_operations_hh78(&
!hhh       VECTOR_DIM,qnode,qgaus,list_elements,time1)
!hhh#else
!hhh  subroutine nsi_element_operations_hh78(&
!hhh       VECTOR_DIM,pnode,pgaus,list_elements,time1)
!hhh#endif
!hhh!@@@    !$acc routine(vecnor,frivel) seq
!hhh    use def_nastin,         only :  kfl_surte_nsi
!hhh    use def_nastin,         only :  fvnoa_nsi
!hhh    use def_kintyp,                     only : ip,rp                                     ! in_const_scalar
!hhh    use def_master,                     only : rhsid                                     ! out       ! real(rp), pointer     :: rhsid(:)
!hhh    use def_nastin,                     only : rhsid_gravb                               ! out       ! real(rp), pointer     :: rhsid_gravb(:) 
!hhh    use def_master,                     only : veloc                                     ! in_var    ! real(rp), pointer     :: veloc(:,:,:)
!hhh    use def_master,                     only : tempe                                     ! in_var    ! real(rp), pointer     :: tempe(:,:)
!hhh    use def_domain,                     only : coord                                     ! in_const  ! real(rp), pointer     :: coord(:,:)
!hhh    use def_domain,                     only : ltype                                     ! in_const  ! integer(ip), pointer  :: ltype(:) -- type of element
!hhh    use def_domain,                     only : lnods                                     ! in_const  ! integer(ip), pointer  :: lnods(:,:)
!hhh    use def_domain,                     only : lmate                                     ! in_const  ! integer(ip), pointer  :: lmate(:)
!hhh    use def_domain,                     only : mnode                                     ! in_const scalar
!hhh    use def_domain,                     only : elmar
!hhh    use def_domain,                     only : lorde  
!hhh    use def_nastin,                     only : dtinv_nsi                                 ! in_var   scalar
!hhh    use def_nastin,                     only : grnor_nsi                                 ! in_const_scalar
!hhh    use def_nastin,                     only : gravi_nsi                                 ! in_const  ! real(rp)         ::gravi_nsi(3)
!hhh    use def_nastin,                     only : densi_aux,visco_aux
!hhh
!hhh    use def_nastin,                     only : kfl_cotem_nsi,bougr_nsi,boube_nsi,gravb_nsi
!hhh
!hhh    use def_nastin,                     only : kfl_force_nsi 
!hhh    use def_nastin,                     only : lforc_material_nsi
!hhh    use def_nastin,                     only : xforc_material_nsi
!hhh    use def_nastin,                     only : fvins_nsi
!hhh    use def_nastin,                     only : NSI_ASSEMBLY_VISCOUS
!hhh    use def_nastin,                     only : ndbgs_nsi
!hhh    use def_nastin,                     only : kfl_fsgrb_nsi
!hhh    use mod_nsi_assembly_global_system, only : nsi_assembly_semi_implicit_method
!hhh    use mod_ker_proper,                 only : ker_proper
!hhh
!hhh    use def_nastin,                     only : corio_nsi,facca_nsi,fvela_nsi,frotc_nsi,centr_nsi,faccl_nsi
!hhh    use def_nastin,                     only : prdivcor_nsi,kfl_stab_corio_nsi
!hhh    use mod_nsi_corio_stab,             only : nsi_corio_stab_element_operations
!hhh    use mod_nsi_corio_stab,             only : element_length_corio
!hhh    use mod_nsi_corio_stab,             only : nsi_element_stabilization_corio
!hhh    !
!hhh    ! Solver
!hhh    !
!hhh    use mod_solver,                     only : solver_assemble_element_matrix_vector
!hhh    use def_nastin,                     only : NSI_SOLVER_VISCOUS_TERM
!hhh    !
!hhh    ! include trucho
!hhh    !
!hhh    use def_kermod,                     only : kfl_nswel_ker
!hhh    use def_kermod,                     only : kfl_noslw_ker
!hhh    use def_parame,                     only : pi
!hhh    !
!hhh    ! TEMP GABLS1
!hhh    !
!hhh#ifdef GABLS1TREF
!hhh    use def_master,                     only : cutim
!hhh#else
!hhh    use def_nastin,                     only : boutr_nsi
!hhh#endif
!hhh
!hhh    use def_nastin,                     only : kfl_dampi_nsi
!hhh    use def_nastin,                     only : top_r_damp_nsi
!hhh    use def_nastin,                     only : top_v_damp_nsi
!hhh    use def_nastin,                     only : bot_r_damp_nsi
!hhh    use def_nastin,                     only : bot_v_damp_nsi
!hhh    use def_nastin,                     only : val_r_damp_nsi 
!hhh    use def_nastin,                     only : mul_v_damp_nsi 
!hhh    use def_nastin,                     only : v_geo_damp_nsi 
!hhh    use def_nastin,                     only : k_dim_damp_nsi 
!hhh    !
!hhh    !
!hhh    implicit none
!hhh    !
!hhh    ! Input and output variables
!hhh    !
!hhh    integer(ip), intent(in)          :: VECTOR_DIM                                       !< Number of nodes
!hhh#if defined PNODE_VALUE && defined PGAUS_VALUE 
!hhh    integer(ip), intent(in)          :: qnode                        !< Number of nodes
!hhh    integer(ip), intent(in)          :: qgaus                        !< Number of Gauss points
!hhh    integer(ip), parameter           :: pnode = PNODE_VALUE
!hhh    integer(ip), parameter           :: pgaus = PGAUS_VALUE 
!hhh#else
!hhh    integer(ip), intent(in)          :: pnode                        !< Number of nodes
!hhh    integer(ip), intent(in)          :: pgaus                        !< Number of Gauss points
!hhh#endif
!hhh    integer(ip), intent(in)          :: list_elements(VECTOR_DIM)                        !< List of elements
!hhh
!hhh    real(rp),    intent(inout)       :: time1(10)                                        ! Timings
!hhh    !
!hhh    ! Element matrices and vectors (stiffness and preconditioner)
!hhh    !
!hhh    real(rp)                         :: elrbu(VECTOR_DIM,ndime,pnode)                    ! bu
!hhh    real(rp)                         :: elrbu2(VECTOR_DIM,ndime,pnode)                   ! bu - Just gravity and Boussinesq
!hhh    real(rp)                         :: elrbp(VECTOR_DIM,pnode)                          ! bp - for gravity in fractional step
!hhh    !
!hhh    ! Gather
!hhh    !
!hhh    integer(ip)                      :: elwvi(VECTOR_DIM)                                ! Wall viscosity element
!hhh    real(rp)                         :: elvel(VECTOR_DIM,ndime,pnode)                    ! u
!hhh    real(rp)                         :: elcod(VECTOR_DIM,ndime,pnode)                    ! x
!hhh    real(rp)                         :: elprdivcor(VECTOR_DIM,pnode)                    ! Projection of the divergence of the coriolis term
!hhh    !
!hhh    ! Indices and dimensions
!hhh    !
!hhh    integer(ip)                      :: ielem,inode,ivect
!hhh    integer(ip)                      :: pelty,porde
!hhh    integer(ip)                      :: ipoin,igaus
!hhh    integer(ip)                      :: lnods_loc(VECTOR_DIM,pnode)
!hhh    integer(ip)                      :: list_elements_p(VECTOR_DIM)                      ! List of elements (always positive)
!hhh    !
!hhh    ! Gauss point values
!hhh    !
!hhh    real(rp)                         :: gpsha(VECTOR_DIM,pnode,pgaus)                    ! N
!hhh    real(rp)                         :: gpcar(VECTOR_DIM,ndime,pnode)                    ! dN/dxi
!hhh    real(rp)                         :: gpvol(VECTOR_DIM)                          ! w*|J|, |J|
!hhh    real(rp)                         :: gpvis(VECTOR_DIM,pgaus)                          ! Viscosity
!hhh    real(rp)                         :: gpwvi(VECTOR_DIM,pgaus)                          ! Viscosity  for no slip wall
!hhh    real(rp)                         :: gpmut(VECTOR_DIM,pgaus)                          ! mut
!hhh    real(rp)                         :: gpden(VECTOR_DIM,pgaus)                          ! Density
!hhh    real(rp)                         :: gptem(VECTOR_DIM,pgaus)                          ! Temperature
!hhh    real(rp)                         :: gpadv(VECTOR_DIM,ndime,pgaus)                    ! u+u'
!hhh    real(rp)                         :: gprhs(VECTOR_DIM,ndime,pgaus)                    ! RHS
!hhh    real(rp)                         :: gpvel(VECTOR_DIM,ndime,pgaus)                    ! u
!hhh    real(rp)                         :: gpgve(VECTOR_DIM,ndime,ndime)              ! grad(u)
!hhh    real(rp)                         :: xjaci(VECTOR_DIM,ndime,ndime)
!hhh    real(rp)                         :: gpdet(VECTOR_DIM)
!hhh
!hhh    real(rp)                         :: gpcod(1:VECTOR_SIZE,3)                   ! Axes motion
!hhh    ! Beware gpcod is being used several times and it recalculated each time. For the moment I leave it like this even if it is not efficient because it can be either the coordinate
!hhh    ! or the coordinate minus teh center of rotation. Needs to be tidier. 
!hhh    real(rp)                         :: alpha(1:VECTOR_SIZE,3)                   ! Axes motion
!hhh    real(rp)                         :: dummr(1:VECTOR_SIZE,3)                   ! Axes motion
!hhh    real(rp)                         :: centf(1:VECTOR_SIZE,3)                   ! Axes motion
!hhh
!hhh    real(rp)                         :: bulkv
!hhh    real(rp)                         :: timea,timeb
!hhh    integer(ip)                      :: idime,pmate
!hhh    integer(ip)                      :: jdime
!hhh    integer(ip)                      :: iauxi,dummi
!hhh    !
!hhh    ! TEMP GABLS1
!hhh    !
!hhh#ifdef GABLS1TREF
!hhh    real(rp)    :: boutr_aux(VECTOR_SIZE)             ! reference temperature
!hhh    real(rp)    :: twall,ttop
!hhh#endif
!hhh    !
!hhh    ! Possibly needed or coriolis stab
!hhh    !
!hhh    real(rp)         :: chale(VECTOR_SIZE,2)
!hhh    real(rp)         :: gpstcor(VECTOR_SIZE,pgaus)
!hhh    real(rp)         :: gppor(VECTOR_DIM,pgaus)               ! porosity - for the moment left to 0 - but nsi_element_stabilization_corio is erady to receive it
!hhh    !
!hhh    ! For Vreman
!hhh    !
!hhh    real(rp),parameter        :: hnatu = 1.0   ! would need to correct for elem /=1!!!!
!hhh    integer(ip)               :: kdime
!hhh    real(rp)                  :: xnutu(VECTOR_DIM),xmile(VECTOR_DIM)
!hhh    real(rp)                  :: G__ij(VECTOR_DIM,3,3)
!hhh    real(rp)                  :: Aalph(VECTOR_DIM),Bbeta(VECTOR_DIM)
!hhh    real(rp),parameter        :: const = 0.1_rp   ! seria mejor poner el valor leido
!hhh    real(rp),parameter        :: zeror = epsilon(1.0_rp) ! podria usarlo de def_master 
!hhh
!hhh! P1
!hhh
!hhh    real(rp),parameter                   :: weigp= 0.041666666666666664D+00
!hhh    real(rp),parameter                   :: epsilgeo_div = epsilon(1.0_rp) !epsilgeo usad to avoid divisions by zero
!hhh    real(rp),parameter                   :: sha_aux(4,4) = reshape((/ 0.585410196624969D+00 ,  0.138196601125011D+00 , &
!hhh         0.138196601125011D+00 ,  0.138196601125011D+00 ,  0.138196601125010D+00 ,  &
!hhh         0.585410196624969D+00 ,  0.138196601125011D+00 ,  0.138196601125011D+00 ,  &
!hhh         0.138196601125010D+00 ,  0.138196601125011D+00 ,  0.585410196624969D+00 , &
!hhh         0.138196601125011D+00 ,  0.138196601125010D+00 ,&
!hhh         0.138196601125011D+00 ,  0.138196601125011D+00 ,  0.585410196624969D+00/), (/4,4/))
!hhh
!hhh
!hhh
!hhh    !
!hhh    ! Internal
!hhh    !
!hhh
!hhh#ifdef OPENACCHHH
!hhh#define DEF_VECT ivect
!hhh#else
!hhh#ifdef VECTOR_SIZE_VARIABLE
!hhh#define DEF_VECT 1:VECTOR_DIM
!hhh#else
!hhh#define DEF_VECT 1:VECTOR_SIZE
!hhh#endif
!hhh#endif
!hhh
!hhh#ifdef OPENACCHHH
!hhh
!hhh#define FACT1X     fact1
!hhh#define FACT2X     fact2
!hhh#define DTINV_LOCX dtinv_loc
!hhh#define T1X        t1
!hhh#define T2X        t2
!hhh#define T3X        t3
!hhh#define DENOMX     denom
!hhh
!hhh#else
!hhh
!hhh#ifdef VECTOR_SIZE_VARIABLE
!hhh#define FACT1X     fact1(1:VECTOR_DIM)
!hhh#define FACT2X     fact2(1:VECTOR_DIM)
!hhh#define DTINV_LOCX dtinv_loc(1:VECTOR_DIM)
!hhh#define T1X        t1(1:VECTOR_DIM)
!hhh#define T2X        t2(1:VECTOR_DIM)
!hhh#define T3X        t3(1:VECTOR_DIM)
!hhh#define DENOMX     denom(1:VECTOR_DIM)
!hhh#else
!hhh#define FACT1X     fact1(1:VECTOR_SIZE)
!hhh#define FACT2X     fact2(1:VECTOR_SIZE)
!hhh#define DTINV_LOCX dtinv_loc(1:VECTOR_SIZE)
!hhh#define T1X        t1(1:VECTOR_SIZE)
!hhh#define T2X        t2(1:VECTOR_SIZE)
!hhh#define T3X        t3(1:VECTOR_SIZE)
!hhh#define DENOMX     denom(1:VECTOR_SIZE)
!hhh#endif
!hhh
!hhh#endif
!hhh
!hhh    real(rp)    :: FACT1X
!hhh    real(rp)    :: FACT2X
!hhh    real(rp)    :: DTINV_LOCX
!hhh    real(rp)    :: T1X
!hhh    real(rp)    :: T2X
!hhh    real(rp)    :: T3X
!hhh    real(rp)    :: DENOMX
!hhh
!hhh    real(rp)         :: eltem(VECTOR_SIZE,pnode)
!hhh
!hhh    !--------------------------------------------------------------------
!hhh    !
!hhh    ! Gather: global to local
!hhh    !
!hhh    !--------------------------------------------------------------------
!hhh
!hhh
!hhh
!hhh#if defined PNODE_VALUE && defined PGAUS_VALUE 
!hhh    if(mnode/=pnode) call runend('mnode/=pnode')   ! from 76 onwards I oly accept mnode==pnode==4
!hhh#endif
!hhh    call cputim(timea)
!hhh
!hhh    ielem = list_elements(1)
!hhh    pelty = abs(ltype(ielem))
!hhh    porde = lorde(pelty)
!hhh    pmate = lmate(ielem)
!hhh
!hhh   !$acc enter data create( gpsha, list_elements_p) &
!hhh   !$acc copyin(list_elements, elmar, elmar(pelty)%shape)
!hhh
!hhh    !$acc parallel loop gang vector default(present)
!hhh    do ivect = 1,VECTOR_DIM                      
!hhh       ielem = abs(list_elements(ivect))
!hhh       if( ielem /= 0 ) then
!hhh          list_elements_p(ivect)   = list_elements(ivect)
!hhh       else
!hhh          list_elements_p(ivect)   = list_elements(1)
!hhh       end if
!hhh
!hhh        gpsha(ivect,1:pnode,1:pgaus) = elmar(pelty) % shape(1:pnode,1:pgaus)
!hhh    end do
!hhh    !$acc end parallel loop
!hhh
!hhh    !$acc update self ( gpsha, list_elements_p)
!hhh
!hhh!@    do ivect = 1,VECTOR_DIM      
!hhh!@       ielem = list_elements_p(ivect)
!hhh!@       gpsha(ivect,1:pnode,1:pgaus) = elmar(pelty) % shape(1:pnode,1:pgaus)
!hhh!@    end do
!hhh
!hhh
!hhh    !$acc enter data create( gpden, gpvis, gpmut, gpwvi) 
!hhh 
!hhh    call ker_proper('DENSI','PGAUS',dummi,list_elements_p,gpden,pnode,pgaus,gpsha,VECTOR_DIM=VECTOR_DIM) ! rho
!hhh    call ker_proper('VISCO','PGAUS',dummi,list_elements_p,gpvis,pnode,pgaus,gpsha,VECTOR_DIM=VECTOR_DIM) ! mu
!hhh    
!hhh!   call ker_proper('TURBU','PGAUS',dummi,list_elements_p,gpmut,pnode,pgaus,gpsha,VECTOR_DIM=VECTOR_DIM) ! mut
!hhh!   call ker_proper('WALLV','PGAUS',dummi,list_elements_p,gpwvi,pnode,pgaus,gpsha,VECTOR_DIM=VECTOR_DIM) ! mu_w
!hhh
!hhh    !$acc update device(gpmut, gpwvi )
!hhh
!hhh
!hhh    if(kfl_surte_nsi == 2_ip) then
!hhh        call runend('kfl_surte_nsi == 2_ip not ready need revision')
!hhh#ifdef OPENACCHHH
!hhh        !$acc parallel loop gang vector default(present)
!hhh        do ivect = 1,VECTOR_DIM
!hhh#endif
!hhh           do igaus = 1,pgaus
!hhh              gpvis(DEF_VECT, igaus) = gpvis(DEF_VECT, igaus)/gpden(DEF_VECT, igaus)
!hhh              gpden(DEF_VECT, igaus) = 1.0_rp
!hhh           end do
!hhh#ifdef OPENACCHHH
!hhh        end do
!hhh        !$acc end parallel loop
!hhh#endif
!hhh    endif
!hhh
!hhh
!hhh    DTINV_LOCX = dtinv_nsi
!hhh    bulkv      = (1.0_rp-(fvins_nsi-1.0_rp)*2.0_rp/3.0_rp)
!hhh
!hhh    call cputim(timeb)
!hhh    time1(3) = time1(3) + timeb - timea
!hhh    call cputim(timea)
!hhh
!hhh
!hhh    !--------------------------------------------------------------------
!hhh    !
!hhh    ! Shape function and derivatives   !ver si gpvis_nsw tine que entrar aca
!hhh    ! added avelavv  becaise it was giving me error at run time not sure if it is correct
!hhh    ! tambien avupo_ker,elavv
!hhh    ! Aunque estoy corriendo un caso sin no slip walllaw y esas solo aparcen dentro de if las quiere igual
!hhh    !
!hhh    !--------------------------------------------------------------------
!hhh
!hhh    !$acc enter data create(  lnods_loc , gpvol     ,                     &   
!hhh    !$acc              elvel     ,  gpadv , gpvel ,                       &
!hhh    !$acc              gpgve   , elcod     ,  gpcar ,                     &
!hhh    !$acc              gpdet   ,                                          &
!hhh    !$acc              xjaci     , elwvi     ,                  &
!hhh    !$acc              gprhs   , elrbu                         )          &
!hhh    !$acc copyin(                                               &
!hhh    !$acc              xforc_material_nsi,lforc_material_nsi,             &
!hhh    !$acc              kfl_nswel_ker,                                     &
!hhh    !$acc              elmar(pelty)%deriv ,                               &
!hhh    !$acc              xmile,xnutu, G__ij   ,alpha  , Bbeta,              &
!hhh    !$acc              elmar(pelty)%weigp, gravi_nsi                      )
!hhh        !
!hhh    ! This do starts here both in the openacc version and in the nonopenacc
!hhh    ! for the non opencc it ends 30  lines later,
!hhh    ! in the openacc case it covers all the subroutine.
!hhh    ! Similarly the scatter in the nonopenacc case needs a do ivect
!hhh    !    
!hhh    !$acc parallel loop gang vector default(present)
!hhh    !
!hhh    do ivect = 1,VECTOR_DIM                      
!hhh       ielem = abs(list_elements(ivect))
!hhh       if( ielem /= 0 ) then
!hhh          lnods_loc(ivect,1:pnode) = lnods(1:pnode,ielem)
!hhh          ielem                    = list_elements(ivect)
!hhh       else
!hhh          lnods_loc(ivect,1:pnode) = lnods(1:pnode,list_elements(1))
!hhh          ielem                    = list_elements(1)
!hhh       end if
!hhh       !
!hhh       ! Transient
!hhh       !
!hhh       do inode = 1,pnode
!hhh          ipoin = lnods_loc(ivect,inode)
!hhh          elvel(ivect,1:ndime,inode) = veloc(1:ndime,ipoin,1)
!hhh          elcod(ivect,1:ndime,inode) = coord(1:ndime,ipoin)
!hhh       end do
!hhh       !
!hhh       ! Coriolis stabilisation
!hhh       !
!hhh       if( kfl_stab_corio_nsi > 0_ip ) then
!hhh          do inode = 1,pnode
!hhh             ipoin = lnods_loc(ivect,inode)
!hhh             elprdivcor(ivect,inode) = prdivcor_nsi(ipoin)
!hhh          end do
!hhh       else
!hhh          elprdivcor(ivect,:) = 0.0_rp
!hhh       end if 
!hhh       !
!hhh       ! Temperature coupling
!hhh       !
!hhh       if( kfl_cotem_nsi == 1) then
!hhh          do inode = 1,pnode
!hhh             ipoin = lnods_loc(ivect,inode)
!hhh             eltem(ivect,inode) = tempe(ipoin,1)
!hhh          end do
!hhh       end if
!hhh       !
!hhh       ! no slip wall law - I could add flag so that i does it only on those elements where it is needed kfl_nswel_ker(ielem)
!hhh       !
!hhh       if ( kfl_noslw_ker /= 0_ip ) then
!hhh          elwvi(ivect) = kfl_nswel_ker(ielem)
!hhh       else          
!hhh          elwvi(ivect) = 0.0_rp
!hhh       end if
!hhh       
!hhh#ifndef OPENACCHHH
!hhh    end do
!hhh    !
!hhh    ! Obtain chale
!hhh    !
!hhh    if( kfl_stab_corio_nsi > 0_ip )  call element_length_corio(ndime,pnode,pelty,porde,elcod,chale)
!hhh
!hhh    call cputim(timeb)
!hhh    time1(1) = time1(1) + timeb - timea
!hhh    call cputim(timea)
!hhh#endif
!hhh    !
!hhh    ! Why not trying this which can be vectorized easily
!hhh    !
!hhh    !do inode = 1,pnode
!hhh    !   do ivect = 1,VECTOR_DIM                      
!hhh    !      ipoin = lnods_loc(ivect,inode)
!hhh    !      elvel(ivect,1:ndime,inode) = veloc(1:ndime,ipoin,1)
!hhh    !      elcod(ivect,1:ndime,inode) = coord(1:ndime,ipoin)          
!hhh    !   end do
!hhh    !end do
!hhh
!hhh
!hhh    !--------------------------------------------------------------------
!hhh    !
!hhh    ! Element Cartesian derivatives and Jacobian: GPCAR, GPVOL
!hhh    !
!hhh    !--------------------------------------------------------------------
!hhh
!hhh
!hhh    do igaus=1,pgaus
!hhh       do inode=1,pnode
!hhh          gpsha(DEF_VECT,inode,igaus) = sha_aux(inode,igaus)
!hhh       end do
!hhh    end do
!hhh    !
!hhh    ! GPCAR (from elmgeo_cartesian_derivatives_jacobian_vector), and GPVOL
!hhh    !
!hhh    !
!hhh    ! 3D P1 element
!hhh    !
!hhh    gpcar(DEF_VECT,1,1) =  elcod(DEF_VECT,1,2)   - elcod(DEF_VECT,1,1)
!hhh    gpcar(DEF_VECT,1,2) =  elcod(DEF_VECT,1,3)   - elcod(DEF_VECT,1,1)
!hhh    gpcar(DEF_VECT,1,3) =  elcod(DEF_VECT,1,4)   - elcod(DEF_VECT,1,1)
!hhh    gpcar(DEF_VECT,2,1) =  elcod(DEF_VECT,2,2)   - elcod(DEF_VECT,2,1)
!hhh    gpcar(DEF_VECT,2,2) =  elcod(DEF_VECT,2,3)   - elcod(DEF_VECT,2,1)
!hhh    gpcar(DEF_VECT,2,3) =  elcod(DEF_VECT,2,4)   - elcod(DEF_VECT,2,1)
!hhh    gpcar(DEF_VECT,3,1) =  elcod(DEF_VECT,3,2)   - elcod(DEF_VECT,3,1)
!hhh    gpcar(DEF_VECT,3,2) =  elcod(DEF_VECT,3,3)   - elcod(DEF_VECT,3,1)
!hhh    gpcar(DEF_VECT,3,3) =  elcod(DEF_VECT,3,4)   - elcod(DEF_VECT,3,1)
!hhh    T1X          =  gpcar(DEF_VECT,2,2) * gpcar(DEF_VECT,3,3) - gpcar(DEF_VECT,3,2) * gpcar(DEF_VECT,2,3)
!hhh    T2X          = -gpcar(DEF_VECT,2,1) * gpcar(DEF_VECT,3,3) + gpcar(DEF_VECT,3,1) * gpcar(DEF_VECT,2,3)
!hhh    T3X          =  gpcar(DEF_VECT,2,1) * gpcar(DEF_VECT,3,2) - gpcar(DEF_VECT,3,1) * gpcar(DEF_VECT,2,2)
!hhh    gpdet(DEF_VECT)     =  gpcar(DEF_VECT,1,1) * T1X + gpcar(DEF_VECT,1,2) * T2X + gpcar(DEF_VECT,1,3) * T3X
!hhh
!hhh
!hhh    DENOMX       =  1.0_rp / (sign(1.0_rp,gpdet(DEF_VECT))*max(abs(gpdet(DEF_VECT)),epsilgeo_div))
!hhh
!hhh    xjaci(DEF_VECT,1,1) =  T1X * DENOMX
!hhh    xjaci(DEF_VECT,2,1) =  T2X * DENOMX
!hhh    xjaci(DEF_VECT,3,1) =  T3X * DENOMX
!hhh    xjaci(DEF_VECT,2,2) = ( gpcar(DEF_VECT,1,1) * gpcar(DEF_VECT,3,3) - gpcar(DEF_VECT,3,1) * gpcar(DEF_VECT,1,3)) * DENOMX
!hhh    xjaci(DEF_VECT,3,2) = (-gpcar(DEF_VECT,1,1) * gpcar(DEF_VECT,3,2) + gpcar(DEF_VECT,1,2) * gpcar(DEF_VECT,3,1)) * DENOMX
!hhh    xjaci(DEF_VECT,3,3) = ( gpcar(DEF_VECT,1,1) * gpcar(DEF_VECT,2,2) - gpcar(DEF_VECT,2,1) * gpcar(DEF_VECT,1,2)) * DENOMX
!hhh    xjaci(DEF_VECT,1,2) = (-gpcar(DEF_VECT,1,2) * gpcar(DEF_VECT,3,3) + gpcar(DEF_VECT,3,2) * gpcar(DEF_VECT,1,3)) * DENOMX
!hhh    xjaci(DEF_VECT,1,3) = ( gpcar(DEF_VECT,1,2) * gpcar(DEF_VECT,2,3) - gpcar(DEF_VECT,2,2) * gpcar(DEF_VECT,1,3)) * DENOMX
!hhh    xjaci(DEF_VECT,2,3) = (-gpcar(DEF_VECT,1,1) * gpcar(DEF_VECT,2,3) + gpcar(DEF_VECT,2,1) * gpcar(DEF_VECT,1,3)) * DENOMX
!hhh
!hhh    gpcar(DEF_VECT,1,1) = -xjaci(DEF_VECT,1,1) - xjaci(DEF_VECT,2,1) - xjaci(DEF_VECT,3,1)
!hhh    gpcar(DEF_VECT,1,2) =  xjaci(DEF_VECT,1,1)
!hhh    gpcar(DEF_VECT,1,3) =  xjaci(DEF_VECT,2,1)
!hhh    gpcar(DEF_VECT,1,4) =  xjaci(DEF_VECT,3,1)
!hhh    gpcar(DEF_VECT,2,1) = -xjaci(DEF_VECT,1,2) - xjaci(DEF_VECT,2,2) - xjaci(DEF_VECT,3,2)
!hhh    gpcar(DEF_VECT,2,2) =  xjaci(DEF_VECT,1,2)
!hhh    gpcar(DEF_VECT,2,3) =  xjaci(DEF_VECT,2,2)
!hhh    gpcar(DEF_VECT,2,4) =  xjaci(DEF_VECT,3,2)
!hhh    gpcar(DEF_VECT,3,1) = -xjaci(DEF_VECT,1,3) - xjaci(DEF_VECT,2,3) - xjaci(DEF_VECT,3,3)
!hhh    gpcar(DEF_VECT,3,2) =  xjaci(DEF_VECT,1,3)
!hhh    gpcar(DEF_VECT,3,3) =  xjaci(DEF_VECT,2,3)
!hhh    gpcar(DEF_VECT,3,4) =  xjaci(DEF_VECT,3,3)
!hhh
!hhh    gpvol(DEF_VECT) = weigp * gpdet(DEF_VECT)
!hhh!!!!!!!!!!!!!!!!!!!
!hhh    xmile(DEF_VECT) = ((  hnatu/sqrt(xjaci(DEF_VECT,1,1) * xjaci(DEF_VECT,1,1) &
!hhh          &     + xjaci(DEF_VECT,1,2) * xjaci(DEF_VECT,1,2)           &
!hhh          &     + xjaci(DEF_VECT,1,3) * xjaci(DEF_VECT,1,3)) ) *      &
!hhh          &     (hnatu/sqrt(xjaci(DEF_VECT,2,1) * xjaci(DEF_VECT,2,1) &
!hhh          &     + xjaci(DEF_VECT,2,2) * xjaci(DEF_VECT,2,2)           &
!hhh          &     + xjaci(DEF_VECT,2,3) * xjaci(DEF_VECT,2,3))) *       &
!hhh          &     (hnatu/sqrt(xjaci(DEF_VECT,3,1) * xjaci(DEF_VECT,3,1) &
!hhh          &     + xjaci(DEF_VECT,3,2) * xjaci(DEF_VECT,3,2)           &
!hhh          &     + xjaci(DEF_VECT,3,3) * xjaci(DEF_VECT,3,3))) )** 0.3333333_rp
!hhh
!hhh!!!!!!!!!!!!!!!!!!!
!hhh
!hhh    !--------------------------------------------------------------------
!hhh    !
!hhh    ! Properties 
!hhh    !
!hhh    !--------------------------------------------------------------------
!hhh
!hhh!   gpmut(DEF_VECT,:) = gpden(DEF_VECT,:) * gpmut(DEF_VECT,:)           ! mut   --- need to pt it later
!hhh    
!hhh    if(kfl_dampi_nsi /= 1_ip) then   ! normal behaviour
!hhh!      gpvis(DEF_VECT,:) = gpvis(DEF_VECT,:) + gpmut(DEF_VECT,:)  ! Effective viscosity <= mu+mut     !later
!hhh#ifndef OPENACCHHH
!hhh    else                           ! with damping
!hhh       call runend('damping not ready')
!hhh       do igaus=1,pgaus
!hhh          gpcod = 0.0_rp
!hhh          !
!hhh          ! Here I could just define the y coordinate
!hhh          !
!hhh          do idime = 1,ndime
!hhh             do inode = 1,pnode
!hhh                gpcod(DEF_VECT,idime) = gpcod(DEF_VECT,idime) + gpsha(DEF_VECT,inode,igaus) * elcod(DEF_VECT,idime,inode)
!hhh             end do
!hhh          end do
!hhh
!hhh          gpvis(DEF_VECT,igaus) = gpvis(DEF_VECT,igaus) + gpmut(DEF_VECT,igaus) + gpvis(DEF_VECT,igaus) * mul_v_damp_nsi  *      &
!hhh               0.5_rp * ( 1.0_rp - cos ( pi * max ( (gpcod(DEF_VECT,2) - bot_v_damp_nsi ) / (top_v_damp_nsi - bot_v_damp_nsi )  &
!hhh               ,0.0_rp) ) )    
!hhh       end do
!hhh#endif
!hhh    endif
!hhh
!hhh    if( kfl_noslw_ker /= 0 ) then
!hhh#ifndef OPENACCHHH
!hhh       call runend('no slip wall not ready')
!hhh       do ivect = 1,VECTOR_DIM
!hhh          if( elwvi(ivect) /= 0 ) then
!hhh             gpvis(ivect,:) = gpwvi(ivect,:)                            ! Wall viscosity
!hhh          end if
!hhh       end do
!hhh#else
!hhh       if( elwvi(DEF_VECT) /= 0 ) then
!hhh           gpvis(DEF_VECT,:) = gpwvi(DEF_VECT,:)                            ! Wall viscosity
!hhh       end if
!hhh#endif
!hhh    end if
!hhh    
!hhh    !----------------------------------------------------------------------
!hhh    !
!hhh    ! Gauss point values
!hhh    !
!hhh    !----------------------------------------------------------------------
!hhh
!hhh    gpgve(DEF_VECT,:,:) = 0.0_rp
!hhh    gprhs(DEF_VECT,:,:)   = 0.0_rp
!hhh    gpvel(DEF_VECT,:,:)   = 0.0_rp
!hhh    gptem(DEF_VECT,:)     = 0.0_rp
!hhh
!hhh    do idime = 1,ndime
!hhh       do inode = 1,pnode
!hhh          do jdime = 1,ndime
!hhh             gpgve(DEF_VECT,jdime,idime) = gpgve(DEF_VECT,jdime,idime) &
!hhh                  + gpcar(DEF_VECT,jdime,inode) * elvel(DEF_VECT,idime,inode)
!hhh          end do
!hhh       end do
!hhh    end do
!hhh    do igaus = 1,pgaus
!hhh       FACT2X =  gpden(DEF_VECT,igaus) * grnor_nsi
!hhh       do idime = 1,ndime
!hhh          do inode = 1,pnode
!hhh             gpvel(DEF_VECT,idime,igaus) = gpvel(DEF_VECT,idime,igaus) + elvel(DEF_VECT,idime,inode) * gpsha(DEF_VECT,inode,igaus)
!hhh          end do
!hhh          gpadv(DEF_VECT,idime,igaus) = gpvel(DEF_VECT,idime,igaus)
!hhh          gprhs(DEF_VECT,idime,igaus) = gprhs(DEF_VECT,idime,igaus) + FACT2X * gravi_nsi(idime)
!hhh       end do
!hhh    end do
!hhh    !
!hhh    ! Material force
!hhh    !
!hhh    if( kfl_force_nsi == 1 ) then
!hhh       if( lforc_material_nsi(pmate) == 2 ) then
!hhh          do igaus = 1,pgaus
!hhh             do idime = 1,ndime
!hhh                gprhs(DEF_VECT,idime,igaus) = gprhs(DEF_VECT,idime,igaus) +  xforc_material_nsi(idime,pmate)
!hhh             end do
!hhh          end do
!hhh       end if
!hhh    end if
!hhh
!hhh
!hhh!!!!!!!!!!!!!!!!!!!   VREMAN  - I need to add igaus & DEF_VEC   - beware I also need to obtain xmile a priori
!hhh! en ker_turbul.f90   pense que se usba DEF_VECT pero no supongo que caca va a ser todo mucho mas rapido 
!hhh    G__ij = 0.0_rp ! G = g^T*g
!hhh    do idime = 1_ip,3_ip
!hhh        do jdime = 1_ip,3_ip
!hhh           do kdime = 1_ip,3_ip
!hhh              G__ij(DEF_VECT,idime,jdime) = G__ij(DEF_VECT,idime,jdime) + (gpgve(DEF_VECT,idime,kdime) &
!hhh                   *gpgve(DEF_VECT,jdime,kdime)*xmile(DEF_VECT)*xmile(DEF_VECT))
!hhh           end do
!hhh        end do
!hhh    end do
!hhh
!hhh!    I guess I do not need to add DEF_VECT  and igauss everywhere it shoul be don automaticaly if Aalph,Bbeta & xmile are of size(DEF_VECT,pgaus)
!hhh
!hhh    Aalph = (G__ij(:,1_ip,1_ip) + G__ij(:,2_ip,2_ip) + G__ij(:,3_ip,3_ip))/(xmile*xmile)
!hhh    Bbeta =  G__ij(:,1_ip,1_ip)*G__ij(:,2_ip,2_ip) + G__ij(:,2_ip,2_ip)*G__ij(:,3_ip,3_ip) + G__ij(:,3_ip,3_ip)*G__ij(:,1_ip,1_ip) &
!hhh             - G__ij(:,1_ip,2_ip)*G__ij(:,1_ip,2_ip) - G__ij(:,2_ip,3_ip)*G__ij(:,2_ip,3_ip) - G__ij(:,1_ip,3_ip)*G__ij(:,1_ip,3_ip)
!hhh
!hhh    xnutu = const
!hhh    where ( Aalph > 10.0_rp*zeror )                        ! Avoid divide by zero
!hhh       xnutu = xnutu  * sqrt (max(( Bbeta ) / ( Aalph ), zeror))
!hhh    elsewhere ( Aalph <= 10.0_rp*zeror )
!hhh       xnutu = 0.0_rp
!hhh    endwhere
!hhh    do igaus=1,pgaus
!hhh       gpvis(:,igaus) = gpvis(:,igaus) + gpden(:,igaus) * xnutu   ! creo  correcto  : gpvis(VECTOR_DIM,pgaus), xnutu(VECTOR_DIM,pgaus)
!hhh    end do
!hhh!!!!!!!
!hhh
!hhh#ifndef OPENACCHHH  
!hhh    !
!hhh    ! Boussinesq coupling: -rho*beta*g*(T-Tr)
!hhh    !
!hhh    if( kfl_cotem_nsi == 1) then
!hhh      FACT1X = bougr_nsi * boube_nsi
!hhh       do igaus = 1,pgaus
!hhh          do inode = 1,pnode
!hhh             gptem(DEF_VECT,igaus) = gptem(DEF_VECT,igaus) + eltem(DEF_VECT,inode) * gpsha(DEF_VECT,inode,igaus)
!hhh          end do
!hhh       end do
!hhh
!hhh       do igaus = 1,pgaus
!hhh#ifdef GABLS1TREF
!hhh          !
!hhh          ! Here I could just define the y coordinate
!hhh          !
!hhh          gpcod = 0.0_rp
!hhh          do idime = 1,ndime
!hhh             do inode = 1,pnode
!hhh                gpcod(DEF_VECT,idime) = gpcod(DEF_VECT,idime) + gpsha(DEF_VECT,inode,igaus) * elcod(DEF_VECT,idime,inode)
!hhh             end do
!hhh          end do
!hhh
!hhh          twall     = 265.0_rp - 6.9444e-05 * cutim
!hhh          ttop      = 268.0_rp
!hhh          boutr_aux = twall + (ttop-twall) * gpcod(DEF_VECT,2) / 400.0_rp
!hhh          FACT2X = gpden(DEF_VECT,igaus) * FACT1X * ( gptem(DEF_VECT,igaus) - boutr_aux )
!hhh#else
!hhh          FACT2X = gpden(DEF_VECT,igaus) * FACT1X * ( gptem(DEF_VECT,igaus) - boutr_nsi )
!hhh#endif
!hhh          do idime = 1,ndime
!hhh             gprhs(DEF_VECT,idime,igaus) = gprhs(DEF_VECT,idime,igaus) - FACT2X * gravb_nsi(idime)
!hhh          end do
!hhh       end do
!hhh    end if
!hhh    !----------------------------------------------------------------------
!hhh    !
!hhh    ! GPRHS
!hhh    !
!hhh    ! Rotation term:            f = f - rho * ( w x w x r + dw/dt x r )
!hhh    ! Linear acceleration term: f = f - rho * a
!hhh    !
!hhh    !----------------------------------------------------------------------
!hhh    ! Me lo traje directo de la version no fast -
!hhh    ! Tratare de meter los menos cambiso possibles
!hhh    ! paar eliminar el tema de ifs usar el truco de Gabriel Stafelbach de force inline -- tambien paar lo de boussi
!hhh    ! OJO notar que esto tambien va a parar a gprhs  con las fuerzas de garvedad etc  y va a mi tartamiento paar FS con GRAVEDAD - en pricipio creo esta ok
!hhh    if( corio_nsi > 1.0e-12_rp ) then
!hhh
!hhh       do igaus = 1,pgaus
!hhh          !
!hhh          ! Gauss point coordinates wrt center of rotation
!hhh          !
!hhh          gpcod = 0.0_rp
!hhh          do idime = 1,ndime
!hhh             do inode = 1,pnode
!hhh                gpcod(DEF_VECT,idime) = gpcod(DEF_VECT,idime) + gpsha(DEF_VECT,inode,igaus) * elcod(DEF_VECT,idime,inode)
!hhh             end do
!hhh             gpcod(DEF_VECT,idime) = gpcod(DEF_VECT,idime) - frotc_nsi(idime)
!hhh          end do
!hhh          !
!hhh          ! Angular acceleration dw/dt x r
!hhh          !
!hhh          if( ndime == 2 ) then
!hhh             alpha(DEF_VECT,1) =-facca_nsi(3) * gpcod(DEF_VECT,2)
!hhh             alpha(DEF_VECT,2) = facca_nsi(3) * gpcod(DEF_VECT,1)
!hhh          else if( ndime == 3 ) then
!hhh             alpha(DEF_VECT,1) = facca_nsi(2) * gpcod(DEF_VECT,3) - facca_nsi(3) * gpcod(DEF_VECT,2)
!hhh             alpha(DEF_VECT,2) = facca_nsi(3) * gpcod(DEF_VECT,1) - facca_nsi(1) * gpcod(DEF_VECT,3)
!hhh             alpha(DEF_VECT,3) = facca_nsi(1) * gpcod(DEF_VECT,2) - facca_nsi(2) * gpcod(DEF_VECT,1)
!hhh          end if
!hhh          !
!hhh          ! Centrifugal force w x (w x r)
!hhh          !
!hhh          if( frotc_nsi(1) < 1.0e10_rp ) then
!hhh             if( ndime == 2 ) then
!hhh                dummr(DEF_VECT,1) =-fvela_nsi(3) * gpcod(DEF_VECT,2)
!hhh                dummr(DEF_VECT,2) = fvela_nsi(3) * gpcod(DEF_VECT,1)
!hhh                centf(DEF_VECT,1) =-fvela_nsi(3) * dummr(DEF_VECT,2)
!hhh                centf(DEF_VECT,2) = fvela_nsi(3) * dummr(DEF_VECT,1)
!hhh             else if( ndime ==3 ) then
!hhh                dummr(DEF_VECT,1) = fvela_nsi(2) * gpcod(DEF_VECT,3) - fvela_nsi(3) * gpcod(DEF_VECT,2)
!hhh                dummr(DEF_VECT,2) = fvela_nsi(3) * gpcod(DEF_VECT,1) - fvela_nsi(1) * gpcod(DEF_VECT,3)
!hhh                dummr(DEF_VECT,3) = fvela_nsi(1) * gpcod(DEF_VECT,2) - fvela_nsi(2) * gpcod(DEF_VECT,1)
!hhh                centf(DEF_VECT,1) = fvela_nsi(2) * dummr(DEF_VECT,3) - fvela_nsi(3) * dummr(DEF_VECT,2)
!hhh                centf(DEF_VECT,2) = fvela_nsi(3) * dummr(DEF_VECT,1) - fvela_nsi(1) * dummr(DEF_VECT,3)
!hhh                centf(DEF_VECT,3) = fvela_nsi(1) * dummr(DEF_VECT,2) - fvela_nsi(2) * dummr(DEF_VECT,1)
!hhh             end if
!hhh             centf(DEF_VECT,1:3) = centf(DEF_VECT,1:3) * centr_nsi
!hhh          else
!hhh             centf = 0.0_rp
!hhh          end if
!hhh          !
!hhh          ! Total force: rho * [ - w x (w x r) - dw/dt x r - a ]
!hhh          !
!hhh          do idime = 1,ndime
!hhh             gprhs(DEF_VECT,idime,igaus) = gprhs(DEF_VECT,idime,igaus) - gpden(DEF_VECT,igaus) &
!hhh                  &               * (   centf(DEF_VECT,idime)                &     ! w x (w x r)
!hhh                  &                   + alpha(DEF_VECT,idime)                &     ! dw/dt x r
!hhh                  &                   + faccl_nsi(idime)                     )     ! a
!hhh          end do
!hhh       end do
!hhh
!hhh    end if
!hhh#endif
!hhh
!hhh    !
!hhh    ! Tau and Tim
!hhh    !
!hhh#ifndef OPENACCHHH
!hhh    !
!hhh    ! Obtain tau for coriolis stabilization -- gpstcor
!hhh    !
!hhh    if( kfl_stab_corio_nsi > 0_ip) then  ! tal ves meter este if dentro
!hhh       call runend('corio stab not ready')
!hhh       gppor = 0.0_rp
!hhh       call nsi_element_stabilization_corio(&
!hhh            pgaus,pnode,chale,gpadv,gpvis,gpden,  &
!hhh            gppor,gpstcor)
!hhh    else
!hhh       gpstcor = 0.0_rp
!hhh    end if
!hhh#endif
!hhh    !----------------------------------------------------------------------
!hhh    !
!hhh    ! Element matrices 
!hhh    !
!hhh    !----------------------------------------------------------------------
!hhh
!hhh    elrbu(DEF_VECT,:,:)   = 0.0_rp
!hhh
!hhh#ifndef OPENACCHHH
!hhh    if(kfl_fsgrb_nsi==1) elrbu2(DEF_VECT,:,:)  = 0.0_rp
!hhh    if(kfl_fsgrb_nsi==1) elrbp(DEF_VECT,:)     = 0.0_rp
!hhh#endif
!hhh
!hhh    !--------------------------------------------------------------------
!hhh    !
!hhh    ! Convective,viscous term and RHS  -- todo tidy up
!hhh    !
!hhh    !--------------------------------------------------------------------
!hhh
!hhh         do inode = 1, pnode
!hhh            do idime = 1, ndime
!hhh               elrbu(DEF_VECT, idime, inode) = 0.0
!hhh               do igaus = 1, pgaus
!hhh                elrbu(DEF_VECT, idime, inode) = elrbu(DEF_VECT, idime, inode) + gpvol(DEF_VECT)*gpsha(DEF_VECT,inode, igaus)*gprhs(DEF_VECT, idime, igaus)
!hhh                  do jdime = 1, ndime
!hhh                     elrbu(DEF_VECT, idime, inode) = elrbu(DEF_VECT, idime, inode) &
!hhh                        - densi_aux * gpvol(DEF_VECT) * gpsha(DEF_VECT,inode,igaus) * gpvel(DEF_VECT,jdime,igaus) * gpgve(DEF_VECT,jdime,idime)     &     ! non-conservative - (u.grad) u
!hhh                        - gpvis(DEF_VECT,igaus) * gpvol(DEF_VECT) * gpcar(DEF_VECT,jdime,inode) * (gpgve(DEF_VECT,jdime,idime) + gpgve(DEF_VECT,idime,jdime))         ! viscous fvins_nsi==1 leads to bulkv=1.0
!hhh                  end do
!hhh               end do
!hhh            end do
!hhh         end do
!hhh!
!hhh!    rho * (div u) u 
!hhh!
!hhh         do inode = 1, pnode
!hhh            do idime = 1, ndime
!hhh               do igaus = 1, pgaus
!hhh                  elrbu(DEF_VECT, idime, inode) = elrbu(DEF_VECT, idime, inode) &
!hhh                      - densi_aux * gpvol(DEF_VECT) * gpsha(DEF_VECT,inode,igaus) * gpvel(DEF_VECT,idime,igaus) *  &
!hhh                        (gpgve(DEF_VECT,1,1) + gpgve(DEF_VECT,2,2) + gpgve(DEF_VECT,3,3) )  ! Addit oriol
!hhh               end do
!hhh            end do
!hhh         end do
!hhh
!hhh
!hhh
!hhh
!hhh
!hhh
!hhh#ifndef OPENACCHHH
!hhh    if (kfl_dampi_nsi == 1_ip) then
!hhh
!hhh        do igaus = 1,pgaus
!hhh
!hhh            gpcod = 0.0_rp
!hhh            !
!hhh            ! Here I could just define the y coordinate
!hhh            !
!hhh            do idime = 1,ndime
!hhh               do inode = 1,pnode
!hhh                  gpcod(DEF_VECT,idime) = gpcod(DEF_VECT,idime) + gpsha(DEF_VECT,inode,igaus) * elcod(DEF_VECT,idime,inode)
!hhh               end do
!hhh            end do
!hhh
!hhh            do inode = 1,pnode
!hhh                FACT1X = gpvol(DEF_VECT) * gpsha(DEF_VECT,inode,igaus)  ! ( f , v )
!hhh                do idime = 1,ndime
!hhh
!hhh                    elrbu(DEF_VECT,idime,inode)  = elrbu(DEF_VECT,idime,inode)  -  FACT1X  *  &
!hhh                    gpden(DEF_VECT,igaus) * ( gpvel(DEF_VECT,idime,igaus) - v_geo_damp_nsi(idime)) * val_r_damp_nsi  * &
!hhh                    0.5_rp * ( 1.0_rp - cos ( pi * max ( (gpcod(DEF_VECT,k_dim_damp_nsi ) - bot_r_damp_nsi ) / &
!hhh                        (top_r_damp_nsi - bot_r_damp_nsi ) , 0.0_rp) ) )        
!hhh                end do
!hhh            end do
!hhh        end do
!hhh    end if
!hhh
!hhh
!hhh    if (kfl_fsgrb_nsi == 1_ip) then
!hhh
!hhh        do igaus = 1,pgaus
!hhh            do inode = 1,pnode
!hhh    
!hhh              FACT1X = gpvol(DEF_VECT) * gpsha(DEF_VECT,inode,igaus)  ! ( f , v )
!hhh    
!hhh                 do idime = 1,ndime 
!hhh                    elrbu2(DEF_VECT,idime,inode) = elrbu2(DEF_VECT,idime,inode) + FACT1X * gprhs(DEF_VECT,idime,igaus)
!hhh                 end do
!hhh          !
!hhh          ! bp = - dt ( f , grad(q) )  - correct treatment of gravity with fractional step - Orthogonal projection of grad(p)-f 
!hhh          !
!hhh   
!hhh                 FACT2X = gpvol(DEF_VECT) / dtinv_nsi
!hhh                 do idime = 1,ndime
!hhh                    elrbp(DEF_VECT,inode) = elrbp(DEF_VECT,inode) - FACT2X * gpcar(DEF_VECT,idime,inode) * gprhs(DEF_VECT,idime,igaus) 
!hhh                 end do
!hhh    
!hhh            end do
!hhh        end do
!hhh
!hhh    end if
!hhh
!hhh#endif
!hhh
!hhh    !--------------------------------------------------------------------
!hhh    !
!hhh    ! Viscous term
!hhh    !
!hhh    !--------------------------------------------------------------------
!hhh
!hhh#ifndef OPENACCHHH
!hhh
!hhh    if( kfl_stab_corio_nsi > 0_ip ) then  ! i could perhaps unify this better
!hhh       !
!hhh       ! Coriolis stabilization --  ( tau_cor * div(omega X u) , div(omega X v) )   -- also galerking
!hhh       !
!hhh       call runend('corio_stab not ready')
!hhh!      call nsi_corio_stab_element_operations(pnode,pgaus,gpvol,gpden,gpsha,gpcar,gpstcor,elprdivcor,elrbu,elauu)
!hhh    else
!hhh       !
!hhh       ! just Coriolis galerkin term - adapted from whta I had in nsi_corio_stab_element_operations -is170
!hhh       !
!hhh       if(abs(fvnoa_nsi)>1.0e-6) call runend('corio_stab not ready')
!hhh!      call nsi_corio_element_operationshh78(pnode,pgaus,gpvol,gpden,gpsha,elauu)
!hhh    end if    
!hhh#endif
!hhh
!hhh    !--------------------------------------------------------------------
!hhh    !
!hhh    ! Assembly in global syste,
!hhh    !
!hhh    !--------------------------------------------------------------------        
!hhh    !
!hhh    ! Scatter element matrix to global one 
!hhh    ! 
!hhh#ifndef OPENACCHHH
!hhh    call cputim(timeb)
!hhh    time1(7) = time1(7) + timeb - timea
!hhh    call cputim(timea)
!hhh    do ivect = 1,VECTOR_DIM                        
!hhh#endif 
!hhh       ielem = list_elements(ivect)
!hhh       if ( ielem > 0 ) then
!hhh          do inode = 1,pnode
!hhh             ipoin = lnods_loc(ivect,inode)
!hhh             do idime = 1,ndime
!hhh                iauxi = idime + (ipoin-1) * ndime 
!hhh                !$acc atomic update
!hhh#ifdef NO_COLORING
!hhh                !$OMP ATOMIC
!hhh#endif
!hhh                rhsid(iauxi) = rhsid(iauxi) + elrbu(ivect,idime,inode)
!hhh 
!hhh#ifndef OPENACCHHH
!hhh                if(kfl_fsgrb_nsi==1) rhsid_gravb(iauxi) = rhsid_gravb(iauxi) + elrbu2(ivect,idime,inode)
!hhh#endif
!hhh
!hhh             end do
!hhh
!hhh
!hhh#ifndef OPENACCHHH
!hhh             if(kfl_fsgrb_nsi==1) then 
!hhh                iauxi = ipoin + ndbgs_nsi 
!hhh                !$acc atomic update
!hhh#ifdef NO_COLORING
!hhh                !$OMP ATOMIC
!hhh#endif
!hhh                rhsid(iauxi) = rhsid(iauxi) + elrbp(ivect,inode)
!hhh             end if
!hhh#endif
!hhh
!hhh          end do
!hhh       end if
!hhh    end do
!hhh    !$acc end parallel loop
!hhh
!hhh    !
!hhh    ! Assemble matrix: This operation is performed on the GPU as well
!hhh    !
!hhh    if( NSI_ASSEMBLY_VISCOUS ) then
!hhh
!hhh!      call solver_assemble_element_matrix_vector(&
!hhh!           solve(NSI_SOLVER_VISCOUS_TERM:),ndime,pnode,pnode*ndime,&
!hhh!           list_elements,lnods_loc,elauu,visco_nsi)
!hhh    end if
!hhh
!hhh
!hhh
!hhh    !$acc exit data delete( lnods_loc , gpvol ,                    &
!hhh    !$acc              elvel     ,  gpadv , gpvel ,                       &
!hhh    !$acc              gpgve   , elcod     ,  gpcar ,                     &
!hhh    !$acc              gpdet   ,                                          &
!hhh    !$acc              xjaci     , elwvi  ,                     &
!hhh    !$acc              gprhs   , elrbu  , list_elements,                  &
!hhh    !$acc              xforc_material_nsi  , lforc_material_nsi,          &
!hhh    !$acc              gpsha , list_elements_p,                           &
!hhh    !$acc              gpden ,                                            &
!hhh    !$acc              gpvis ,                                            &
!hhh    !$acc              gpmut ,                                            &
!hhh    !$acc              gpwvi ,                                            &
!hhh    !$acc              elmar,                                             &
!hhh    !$acc              elmar(pelty)%shape ,                               &
!hhh    !$acc              elmar(pelty)%deriv ,                               &
!hhh    !$acc              elmar(pelty)%weigp                      )
!hhh
!hhh#ifndef OPENACCHHH
!hhh    call cputim(timeb)
!hhh    time1(8) = time1(8) + timeb - timea                       
!hhh#endif
!hhh
!hhh  end subroutine nsi_element_operations_hh78

end module mod_nsi_element_operations_hh78

