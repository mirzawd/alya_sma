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

module mod_nsi_element_operations_fast_dev

#include "def_vector_size.inc"
  use def_kintyp_basic,               only : ip,rp         !in_const
  use def_domain,                     only : ndime         !in_const_scalar  ! actually paarmeter if ndimepar
#ifndef SUPER_FAST
  use mod_element_integration,        only : element_shape_function_derivatives_jacobian
#endif

#ifdef _OPENACC
  use openacc
#endif

  use def_master,                     only : kfl_paral
  
  implicit none
  
  private
  public :: nsi_element_operations_fast_dev

contains

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
  subroutine nsi_element_operations_fast_dev(&
       VECTOR_DIM,qnode,qgaus,list_elements,time1)
#else
  subroutine nsi_element_operations_fast_dev(&
       VECTOR_DIM,pnode,pgaus,list_elements,time1)
#endif
!@@@    !$acc routine(vecnor,frivel) seq
    use def_nastin,         only :  kfl_surte_nsi
    use def_kintyp,                     only : ip,rp                                     ! in_const_scalar
    use def_master,                     only : rhsid                                     ! out       ! real(rp), pointer     :: rhsid(:)
    use def_nastin,                     only : rhsid_gravb                               ! out       ! real(rp), pointer     :: rhsid_gravb(:) 
    use def_master,                     only : veloc                                     ! in_var    ! real(rp), pointer     :: veloc(:,:,:)
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
    use def_nastin,                     only : dtinv_nsi                                 ! in_var   scalar
    use def_nastin,                     only : dt_rho_nsi                                ! out  ! real(rp), pointer     :: dt_rho_nsi(:)
    use def_nastin,                     only : mass_rho_nsi                              ! out  ! real(rp), pointer     :: mu_rho_nsi(:)
    use def_nastin,                     only : grnor_nsi                                 ! in_const_scalar
    use def_nastin,                     only : gravi_nsi                                 ! in_const  ! real(rp)         ::gravi_nsi(3)

    use def_nastin,                     only : kfl_cotem_nsi,bougr_nsi,boube_nsi,gravb_nsi

    use def_nastin,                     only : kfl_force_nsi 
    use def_nastin,                     only : lforc_material_nsi
    use def_nastin,                     only : xforc_material_nsi
    use def_nastin,                     only : fvins_nsi
    use def_nastin,                     only : NSI_ASSEMBLY_VISCOUS
    use def_nastin,                     only : visco_nsi
    use def_nastin,                     only : ndbgs_nsi
    use def_nastin,                     only : kfl_fsgrb_nsi
    use mod_nsi_assembly_global_system, only : nsi_assembly_semi_implicit_method
    use mod_ker_proper,                 only : ker_proper

    use def_nastin,                     only : corio_nsi,facca_nsi,fvela_nsi,frotc_nsi,centr_nsi,faccl_nsi
    use def_nastin,                     only : prdivcor_nsi,kfl_stab_corio_nsi, kfl_anipo_nsi
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
    use def_kermod,                     only : kfl_dampi_ker
    use def_parame,                     only : pi
    !
    ! TEMP GABLS1
    !
#ifdef GABLS1TREF
    use def_master,                     only : cutim
#else
    use def_nastin,                     only : boutr_nsi
#endif
    !
    ! AB: densi
    !
!   use def_master, only     :  wmean_gp,tempe_gp,prthe
!   use def_kermod, only     :  gasco
    use mod_ker_tendencies,  only : kfl_tendencies_ker, get_tendencies_u, get_tendencies_uwrf 
    use def_master,          only : cutim
    use def_domain,          only : walld
!   use def_nastin,          only : ntabr_nsi, fvela_nsi
    use mod_ker_discrete_function,   only : ker_discrete_function


    use def_nastin,                     only : kfl_dampi_nsi
    use def_nastin,                     only : top_r_damp_nsi
    use def_nastin,                     only : top_v_damp_nsi
    use def_nastin,                     only : bot_r_damp_nsi
    use def_nastin,                     only : bot_v_damp_nsi
    use def_nastin,                     only : val_r_damp_nsi 
    use def_nastin,                     only : mul_v_damp_nsi 
    use def_nastin,                     only : v_geo_damp_nsi 
    use def_nastin,                     only : k_dim_damp_nsi 
    use def_nastin,                     only : kfl_vegeo_time_nsi
    !
    !
    implicit none
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
    ! Element matrices and vectors (stiffness and preconditioner)
    !
    real(rp)                         :: elauu(VECTOR_DIM,pnode*ndime,pnode*ndime)        ! Auu
    real(rp)                         :: elrbu(VECTOR_DIM,ndime,pnode)                    ! bu
    real(rp)                         :: elrbu2(VECTOR_DIM,ndime,pnode)                   ! bu - Just gravity and Boussinesq
    real(rp)                         :: elrbp(VECTOR_DIM,pnode)                          ! bp - for gravity in fractional step
    real(rp)                         :: eldtrho(VECTOR_DIM,pnode)                        ! Projection of rho/dt
    real(rp)                         :: elmurho(VECTOR_DIM,pnode)                        ! Projection of mu/rho
    !
    ! Gather
    !
    integer(ip)                      :: elwvi(VECTOR_DIM)                                ! Wall viscosity element
    real(rp)                         :: elvel(VECTOR_DIM,ndime,pnode)                    ! u
    real(rp)                         :: elcod(VECTOR_DIM,ndime,pnode)                    ! x
    real(rp)                         :: elprdivcor(VECTOR_DIM,pnode)                    ! Projection of the divergence of the coriolis term
    !
    ! Indices and dimensions
    !
    integer(ip)                      :: ielem,inode,ivect
    integer(ip)                      :: pelty,j,k,porde
    integer(ip)                      :: ipoin,igaus
    integer(ip)                      :: lnods_loc(VECTOR_DIM,pnode)
    integer(ip)                      :: list_elements_p(VECTOR_DIM)                      ! List of elements (always positive)
    !
    ! Gauss point values
    !
    real(rp)                         :: gpsha(VECTOR_DIM,pnode,pgaus)                    ! N
    real(rp)                         :: gpcar(VECTOR_DIM,ndime,mnode,pgaus)              ! dN/dxi
    real(rp)                         :: gpvol(VECTOR_DIM,pgaus)                          ! w*|J|, |J|
    real(rp)                         :: gpvis(VECTOR_DIM,pgaus)                          ! Viscosity
    real(rp)                         :: gpwvi(VECTOR_DIM,pgaus)                          ! Viscosity  for no slip wall
    real(rp)                         :: gpmut(VECTOR_DIM,pgaus)                          ! mut
    real(rp)                         :: gpden(VECTOR_DIM,pgaus)                          ! Density
    real(rp)                         :: gptem(VECTOR_DIM,pgaus)                          ! Temperature
    real(rp)                         :: gpadv(VECTOR_DIM,ndime,pgaus)                    ! u+u'
    real(rp)                         :: gprhs(VECTOR_DIM,ndime,pgaus)                    ! RHS
    real(rp)                         :: gpvel(VECTOR_DIM,ndime,pgaus)                    ! u
    real(rp)                         :: gpgve(VECTOR_DIM,ndime,ndime,pgaus)              ! grad(u)
    real(rp)                         :: xjacm(VECTOR_DIM,ndime,ndime)
    real(rp)                         :: xjaci(VECTOR_DIM,ndime,ndime)
    real(rp)                         :: gpdet(VECTOR_DIM,pgaus)
    real(rp)                         :: wgrgr(VECTOR_DIM,pnode,pnode,pgaus)
    real(rp)                         :: agrau(VECTOR_DIM,pnode,pgaus)

    real(rp)                         :: gpcod(1:VECTOR_SIZE,3)                   ! Axes motion
    ! Beware gpcod is being used several times and it recalculated each time. For the moment I leave it like this even if it is not efficient because it can be either the coordinate
    ! or the coordinate minus teh center of rotation. Needs to be tidier. 
    real(rp)                         :: alpha(1:VECTOR_SIZE,3)                   ! Axes motion
    real(rp)                         :: dummr(1:VECTOR_SIZE,3)                   ! Axes motion
    real(rp)                         :: centf(1:VECTOR_SIZE,3)                   ! Axes motion

    real(rp)                         :: T_dtrho(VECTOR_DIM) 
    real(rp)                         :: d_dtrho(VECTOR_DIM) 
    real(rp)                         :: T_murho(VECTOR_DIM) 
    real(rp)                         :: d_murho(VECTOR_DIM) 
    real(rp)                         :: bulkv
    real(rp)                         :: timea,timeb
    integer(ip)                      :: idime,pmate
    integer(ip)                      :: jdime,jnode,idofv,jdofv
    integer(ip)                      :: ievat,jevat,iauxi,dummi
    !
    ! TEMP GABLS1
    !
#ifdef GABLS1TREF
    real(rp)    :: boutr_aux(VECTOR_SIZE)             ! reference temperature
    real(rp)    :: twall,ttop
#endif
    real(rp)    :: fcori
    real(rp)    :: gphei, elhei(pnode)
    real(rp)    :: gpten_geo(2), gpten_adv(2), gpvel_ref(2), ugeos(2)
    !
    ! Possibly needed or coriolis stab
    !
    real(rp)         :: chale(VECTOR_SIZE,2)
    real(rp)         :: gpstcor(VECTOR_SIZE,pgaus)
    real(rp)         :: gppor(VECTOR_DIM,pgaus)               ! porosity - for the moment left to 0 - but nsi_element_stabilization_corio is erady to receive it
    !
    ! Internal
    !
#ifdef OPENACCHHH
#define DEF_VECT ivect
#else
#define DEF_VECT 1:VECTOR_SIZE
#endif

#ifdef OPENACCHHH

#define FACT0X     fact0
#define FACT1X     fact1
#define FACT2X     fact2
#define FACT4X     fact4
!#define FACT5X     fact5
#define FACT6X     fact6
#define DTINV_LOCX dtinv_loc
#define T1X        t1
#define T2X        t2
#define T3X        t3
#define DENOMX     denom

#else

#define FACT0X     fact0(1:VECTOR_SIZE)
#define FACT1X     fact1(1:VECTOR_SIZE)
#define FACT2X     fact2(1:VECTOR_SIZE)
#define FACT4X     fact4(1:VECTOR_SIZE)
!#define FACT5X     fact5(1:VECTOR_SIZE)
#define FACT6X     fact6(1:VECTOR_SIZE)
#define DTINV_LOCX dtinv_loc(1:VECTOR_SIZE)
#define T1X        t1(1:VECTOR_SIZE)
#define T2X        t2(1:VECTOR_SIZE)
#define T3X        t3(1:VECTOR_SIZE)
#define DENOMX     denom(1:VECTOR_SIZE)

#endif

    real(rp)    :: FACT0X
    real(rp)    :: FACT1X
    real(rp)    :: FACT2X
    real(rp)    :: FACT4X
!    real(rp)    :: FACT5X
    real(rp)    :: FACT6X
    real(rp)    :: DTINV_LOCX
    real(rp)    :: T1X
    real(rp)    :: T2X
    real(rp)    :: T3X
    real(rp)    :: DENOMX

    real(rp)         :: eltem(VECTOR_SIZE,pnode)

    !--------------------------------------------------------------------
    !
    ! Gather: global to local
    !
    !--------------------------------------------------------------------

    call cputim(timea)

    ielem = list_elements(1)
    pelty = abs(ltype(ielem))
    porde = lorde(pelty)
    pmate = lmate(ielem)

   !$acc enter data create( gpsha, list_elements_p) &
   !$acc copyin(list_elements, elmar, elmar(pelty)%shape)

    !$acc parallel loop gang vector default(present)
    do ivect = 1,VECTOR_DIM                      
       ielem = abs(list_elements(ivect))
       if( ielem /= 0 ) then
          list_elements_p(ivect)   = list_elements(ivect)
       else
          list_elements_p(ivect)   = list_elements(1)
       end if

        gpsha(ivect,1:pnode,1:pgaus) = elmar(pelty) % shape(1:pnode,1:pgaus)
    end do
    !$acc end parallel loop

    !$acc update self ( gpsha, list_elements_p)

!@    do ivect = 1,VECTOR_DIM      
!@       ielem = list_elements_p(ivect)
!@       gpsha(ivect,1:pnode,1:pgaus) = elmar(pelty) % shape(1:pnode,1:pgaus)
!@    end do


    !$acc enter data create( gpden, gpvis, gpmut, gpwvi) 
 
    call ker_proper('DENSI','PGAUS',dummi,list_elements_p,gpden,pnode,pgaus,gpsha,VECTOR_DIM=VECTOR_DIM) ! rho
    call ker_proper('VISCO','PGAUS',dummi,list_elements_p,gpvis,pnode,pgaus,gpsha,VECTOR_DIM=VECTOR_DIM) ! mu
    call ker_proper('POROS','PGAUS',dummi,list_elements_p,gppor,pnode,pgaus,gpsha,VECTOR_DIM=VECTOR_DIM) ! porosity
    
    call ker_proper('TURBU','PGAUS',dummi,list_elements_p,gpmut,pnode,pgaus,gpsha,VECTOR_DIM=VECTOR_DIM) ! mut
    call ker_proper('WALLV','PGAUS',dummi,list_elements_p,gpwvi,pnode,pgaus,gpsha,VECTOR_DIM=VECTOR_DIM) ! mu_w


    !$acc update device(gpmut, gpwvi )


    if(kfl_surte_nsi == 2_ip) then
#ifdef OPENACCHHH
        !$acc parallel loop gang vector default(present)
        do ivect = 1,VECTOR_DIM
#endif
           do igaus = 1,pgaus
              gpvis(DEF_VECT, igaus) = gpvis(DEF_VECT, igaus)/gpden(DEF_VECT, igaus)
              gpden(DEF_VECT, igaus) = 1.0_rp
           end do
#ifdef OPENACCHHH
        end do
        !$acc end parallel loop
#endif
    endif


    DTINV_LOCX = dtinv_nsi
    bulkv      = (1.0_rp-(fvins_nsi-1.0_rp)*2.0_rp/3.0_rp)

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

    !$acc enter data create( agrau   , lnods_loc , gpvol     , elauu ,    &   
    !$acc              eldtrho , elvel     ,  gpadv , gpvel ,             &
    !$acc              gpgve   , elcod     ,  gpcar ,                     &
    !$acc              wgrgr   , gpdet     ,                              &
    !$acc              xjacm   , xjaci     , elwvi     ,                  &
    !$acc              T_dtrho , d_dtrho   , T_murho, d_murho,            &
    !$acc              gprhs   , elmurho   , elrbu             )          &
    !$acc copyin(      visco_nsi,                                         &
    !$acc              xforc_material_nsi,lforc_material_nsi,             &
    !$acc              kfl_nswel_ker,                                     &
    !$acc              elmar(pelty)%deriv ,                               &
    !$acc              elmar(pelty)%weigp, gravi_nsi                      )
        !
    ! This do starts here both in the openacc version and in the nonopenacc
    ! for the non opencc it ends 30  lines later,
    ! in the openacc case it covers all the subroutine.
    ! Similarly the scatter in the nonopenacc case needs a do ivect
    !    
    !$acc parallel loop gang vector default(present)
    !
    do ivect = 1,VECTOR_DIM                      
       ielem = abs(list_elements(ivect))
       if( ielem /= 0 ) then
          lnods_loc(ivect,1:pnode) = lnods(1:pnode,ielem)
          ielem                    = list_elements(ivect)
       else
          lnods_loc(ivect,1:pnode) = lnods(1:pnode,list_elements(1))
          ielem                    = list_elements(1)
       end if
       !
       ! Transient
       !
       do inode = 1,pnode
          ipoin = lnods_loc(ivect,inode)
          elvel(ivect,1:ndime,inode) = veloc(1:ndime,ipoin,1)
          elcod(ivect,1:ndime,inode) = coord(1:ndime,ipoin)
       end do
       !
       ! Coriolis stabilisation
       !
       if( kfl_stab_corio_nsi > 0_ip ) then
          do inode = 1,pnode
             ipoin = lnods_loc(ivect,inode)
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
             ipoin = lnods_loc(ivect,inode)
             eltem(ivect,inode) = tempe(ipoin,1)
          end do
       end if
       !
       ! no slip wall law - I could add flag so that i does it only on those elements where it is needed kfl_nswel_ker(ielem)
       !
       if ( kfl_noslw_ker /= 0_ip ) then
          elwvi(ivect) = kfl_nswel_ker(ielem)
       else          
          elwvi(ivect) = 0.0_rp
       end if
       
#ifndef OPENACCHHH
    end do
    !
    ! Obtain chale
    !
    if( kfl_stab_corio_nsi > 0_ip )  call element_length_corio(ndime,pnode,pelty,porde,elcod,chale)

    call cputim(timeb)
    time1(1) = time1(1) + timeb - timea
    call cputim(timea)
#endif
    !
    ! Why not trying this which can be vectorized easily
    !
    !do inode = 1,pnode
    !   do ivect = 1,VECTOR_DIM                      
    !      ipoin = lnods_loc(ivect,inode)
    !      elvel(ivect,1:ndime,inode) = veloc(1:ndime,ipoin,1)
    !      elcod(ivect,1:ndime,inode) = coord(1:ndime,ipoin)          
    !   end do
    !end do

    !--------------------------------------------------------------------
    !
    ! Element Cartesian derivatives and Jacobian: GPCAR, GPVOL
    !
    !--------------------------------------------------------------------

    if( kfl_savda == 2 ) then
#ifndef OPENACCHHH
       do ivect = 1,VECTOR_DIM      
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

          do igaus = 1,pgaus
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

          do igaus = 1,pgaus
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


    end if

    !--------------------------------------------------------------------
    !
    ! Properties 
    !
    !--------------------------------------------------------------------

    gpmut(DEF_VECT,:) = gpden(DEF_VECT,:) * gpmut(DEF_VECT,:)           ! mut
    
    if(kfl_dampi_nsi /= 1_ip) then   ! normal behaviour
       gpvis(DEF_VECT,:) = gpvis(DEF_VECT,:) + gpmut(DEF_VECT,:)  ! Effective viscosity <= mu+mut
#ifndef OPENACCHHH
    else                           ! with damping
       do igaus=1,pgaus
          gpcod = 0.0_rp
          !
          ! Here I could just define the y coordinate
          !
          do idime = 1,ndime
             do inode = 1,pnode
                gpcod(DEF_VECT,idime) = gpcod(DEF_VECT,idime) + gpsha(DEF_VECT,inode,igaus) * elcod(DEF_VECT,idime,inode)
             end do
          end do

          gpvis(DEF_VECT,igaus) = gpvis(DEF_VECT,igaus) + gpmut(DEF_VECT,igaus) + gpvis(DEF_VECT,igaus) * mul_v_damp_nsi  *      &
               0.5_rp * ( 1.0_rp - cos ( pi * max ( (gpcod(DEF_VECT,2) - bot_v_damp_nsi ) / (top_v_damp_nsi - bot_v_damp_nsi )  &
               ,0.0_rp) ) )    
       end do
#endif
    endif

    if( kfl_noslw_ker /= 0 ) then
#ifndef OPENACCHHH

       do ivect = 1,VECTOR_DIM
          if( elwvi(ivect) /= 0 ) then
             gpvis(ivect,:) = gpwvi(ivect,:)                            ! Wall viscosity
          end if
       end do
#else
       if( elwvi(DEF_VECT) /= 0 ) then
           gpvis(DEF_VECT,:) = gpwvi(DEF_VECT,:)                            ! Wall viscosity
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
    gptem(DEF_VECT,:)     = 0.0_rp
    eldtrho(DEF_VECT,:)   = 0.0_rp
    elmurho(DEF_VECT,:)   = 0.0_rp

    do igaus = 1,pgaus
       FACT2X =  gpden(DEF_VECT,igaus) * grnor_nsi
       do idime = 1,ndime
          do inode = 1,pnode
             gpvel(DEF_VECT,idime,igaus) = gpvel(DEF_VECT,idime,igaus) + elvel(DEF_VECT,idime,inode) * gpsha(DEF_VECT,inode,igaus)
             do jdime = 1,ndime
                gpgve(DEF_VECT,jdime,idime,igaus) = gpgve(DEF_VECT,jdime,idime,igaus) &
                     + gpcar(DEF_VECT,jdime,inode,igaus) * elvel(DEF_VECT,idime,inode)
             end do
          end do
          gpadv(DEF_VECT,idime,igaus) = gpvel(DEF_VECT,idime,igaus)
          gprhs(DEF_VECT,idime,igaus) = gprhs(DEF_VECT,idime,igaus) + FACT2X * gravi_nsi(idime) &
                                                 - gppor(DEF_VECT,igaus) * gpvel(DEF_VECT,idime,igaus)
       end do
    end do
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
    !
    ! Tendecies
    !
#ifndef OPENACCHHH  
    if( kfl_tendencies_ker ) then
          fcori = 2.0_rp*sqrt(fvela_nsi(1)*fvela_nsi(1) + &
          fvela_nsi(2)*fvela_nsi(2) + &
          fvela_nsi(3)*fvela_nsi(3))
       loop_ivect1:    do ivect = 1,VECTOR_SIZE
          ielem = list_elements(ivect)
          if( ielem > 0 ) then

             ! interpolate tendencies to present height (over terrain?)
             do inode = 1,pnode
                elhei(inode) = walld(lnods_loc(ivect,inode))  !elcod(ivect, 3, 1:pnode)  ! by the moment height is z coordinate
             end do
             !             end do
             do igaus = 1,pgaus
                gphei =0.0_rp
                do inode =1, pnode
                   gphei = gphei + gpsha(ivect,inode, igaus)*elhei(inode)
                end do

                call get_tendencies_u(gphei,cutim, gpten_geo, gpten_adv) ! get geostrophic and advection forces
                gprhs(ivect,1, igaus) =  gprhs(ivect,1, igaus) + gpden(ivect,igaus)*fcori*(-gpten_geo(2) + gpten_adv(1))
                gprhs(ivect,2, igaus) =  gprhs(ivect,2, igaus) + gpden(ivect,igaus)*fcori*( gpten_geo(1) + gpten_adv(2))
                !
                ! RHS DAMPING damp*vel  (puts vel= WRF velocity) (not possible to differentiate between damping and forest)
                ! 
                if (kfl_dampi_ker(pmate).gt.0.and. kfl_anipo_nsi == 0) then
                   call get_tendencies_uwrf(gphei, cutim, gpvel_ref) 
                   gprhs(ivect,1:2, igaus) =  gprhs(ivect,1:2, igaus) + gppor(ivect, igaus)*gpvel_ref(1:2)
                end if
             end do

          end if
       end do loop_ivect1
    else if (kfl_vegeo_time_nsi.gt.0) then  !Time dependent pressure gradient (acts like a time dependent and uniform tendency) 
       fcori = 2.0_rp*sqrt(fvela_nsi(1)*fvela_nsi(1) + &
            fvela_nsi(2)*fvela_nsi(2) + &
            fvela_nsi(3)*fvela_nsi(3))
       loop_ivect2:    do ivect = 1,VECTOR_SIZE
          ielem = list_elements(ivect)
          if( ielem > 0 ) then
!            Time interpolation of geostrophic velocity      
             call ker_discrete_function(kfl_vegeo_time_nsi,cutim,ugeos)
             gprhs(ivect,1, 1:pgaus) =  gprhs(ivect,1, 1:pgaus) -  gpden(ivect,1:pgaus)*fcori*ugeos(2)
             gprhs(ivect,2, 1:pgaus) =  gprhs(ivect,2, 1:pgaus) +  gpden(ivect,1:pgaus)*fcori*ugeos(1)
             ! add RHS damping
             if (kfl_dampi_ker(pmate).gt.0.and. kfl_anipo_nsi == 0) then
                gprhs(ivect,1, 1:pgaus) =  gprhs(ivect,1, 1:pgaus) + gppor(ivect, 1:pgaus)*ugeos(1)
                gprhs(ivect,2, 1:pgaus) =  gprhs(ivect,2, 1:pgaus) + gppor(ivect, 1:pgaus)*ugeos(2)
             end if

          end if
       end do loop_ivect2
   
    end if

    !
    ! Boussinesq coupling: -rho*beta*g*(T-Tr)
    !
    if( kfl_cotem_nsi == 1) then
      FACT1X = bougr_nsi * boube_nsi
       do igaus = 1,pgaus
          do inode = 1,pnode
             gptem(DEF_VECT,igaus) = gptem(DEF_VECT,igaus) + eltem(DEF_VECT,inode) * gpsha(DEF_VECT,inode,igaus)
          end do
       end do

       do igaus = 1,pgaus
#ifdef GABLS1TREF
          !
          ! Here I could just define the y coordinate
          !
          gpcod = 0.0_rp
          do idime = 1,ndime
             do inode = 1,pnode
                gpcod(DEF_VECT,idime) = gpcod(DEF_VECT,idime) + gpsha(DEF_VECT,inode,igaus) * elcod(DEF_VECT,idime,inode)
             end do
          end do

          twall     = 265.0_rp - 6.9444e-05 * cutim
          ttop      = 268.0_rp
          boutr_aux = twall + (ttop-twall) * gpcod(DEF_VECT,2) / 400.0_rp
          FACT2X = gpden(DEF_VECT,igaus) * FACT1X * ( gptem(DEF_VECT,igaus) - boutr_aux )
#else
          FACT2X = gpden(DEF_VECT,igaus) * FACT1X * ( gptem(DEF_VECT,igaus) - boutr_nsi )
#endif
          do idime = 1,ndime
             gprhs(DEF_VECT,idime,igaus) = gprhs(DEF_VECT,idime,igaus) - FACT2X * gravb_nsi(idime)
          end do
       end do
    end if
    !----------------------------------------------------------------------
    !
    ! GPRHS
    !
    ! Rotation term:            f = f - rho * ( w x w x r + dw/dt x r )
    ! Linear acceleration term: f = f - rho * a
    !
    !----------------------------------------------------------------------
    ! Me lo traje directo de la version no fast -
    ! Tratare de meter los menos cambiso possibles
    ! paar eliminar el tema de ifs usar el truco de Gabriel Stafelbach de force inline -- tambien paar lo de boussi
    ! OJO notar que esto tambien va a parar a gprhs  con las fuerzas de garvedad etc  y va a mi tartamiento paar FS con GRAVEDAD - en pricipio creo esta ok
    if( corio_nsi > 1.0e-12_rp ) then

       do igaus = 1,pgaus
          !
          ! Gauss point coordinates wrt center of rotation
          !
          gpcod = 0.0_rp
          do idime = 1,ndime
             do inode = 1,pnode
                gpcod(DEF_VECT,idime) = gpcod(DEF_VECT,idime) + gpsha(DEF_VECT,inode,igaus) * elcod(DEF_VECT,idime,inode)
             end do
             gpcod(DEF_VECT,idime) = gpcod(DEF_VECT,idime) - frotc_nsi(idime)
          end do
          !
          ! Angular acceleration dw/dt x r
          !
          if( ndime == 2 ) then
             alpha(DEF_VECT,1) =-facca_nsi(3) * gpcod(DEF_VECT,2)
             alpha(DEF_VECT,2) = facca_nsi(3) * gpcod(DEF_VECT,1)
          else if( ndime == 3 ) then
             alpha(DEF_VECT,1) = facca_nsi(2) * gpcod(DEF_VECT,3) - facca_nsi(3) * gpcod(DEF_VECT,2)
             alpha(DEF_VECT,2) = facca_nsi(3) * gpcod(DEF_VECT,1) - facca_nsi(1) * gpcod(DEF_VECT,3)
             alpha(DEF_VECT,3) = facca_nsi(1) * gpcod(DEF_VECT,2) - facca_nsi(2) * gpcod(DEF_VECT,1)
          end if
          !
          ! Centrifugal force w x (w x r)
          !
          if( frotc_nsi(1) < 1.0e10_rp ) then
             if( ndime == 2 ) then
                dummr(DEF_VECT,1) =-fvela_nsi(3) * gpcod(DEF_VECT,2)
                dummr(DEF_VECT,2) = fvela_nsi(3) * gpcod(DEF_VECT,1)
                centf(DEF_VECT,1) =-fvela_nsi(3) * dummr(DEF_VECT,2)
                centf(DEF_VECT,2) = fvela_nsi(3) * dummr(DEF_VECT,1)
             else if( ndime ==3 ) then
                dummr(DEF_VECT,1) = fvela_nsi(2) * gpcod(DEF_VECT,3) - fvela_nsi(3) * gpcod(DEF_VECT,2)
                dummr(DEF_VECT,2) = fvela_nsi(3) * gpcod(DEF_VECT,1) - fvela_nsi(1) * gpcod(DEF_VECT,3)
                dummr(DEF_VECT,3) = fvela_nsi(1) * gpcod(DEF_VECT,2) - fvela_nsi(2) * gpcod(DEF_VECT,1)
                centf(DEF_VECT,1) = fvela_nsi(2) * dummr(DEF_VECT,3) - fvela_nsi(3) * dummr(DEF_VECT,2)
                centf(DEF_VECT,2) = fvela_nsi(3) * dummr(DEF_VECT,1) - fvela_nsi(1) * dummr(DEF_VECT,3)
                centf(DEF_VECT,3) = fvela_nsi(1) * dummr(DEF_VECT,2) - fvela_nsi(2) * dummr(DEF_VECT,1)
             end if
             centf(DEF_VECT,1:3) = centf(DEF_VECT,1:3) * centr_nsi
          else
             centf = 0.0_rp
          end if
          !
          ! Total force: rho * [ - w x (w x r) - dw/dt x r - a ]
          !
          do idime = 1,ndime
             gprhs(DEF_VECT,idime,igaus) = gprhs(DEF_VECT,idime,igaus) - gpden(DEF_VECT,igaus) &
                  &               * (   centf(DEF_VECT,idime)                &     ! w x (w x r)
                  &                   + alpha(DEF_VECT,idime)                &     ! dw/dt x r
                  &                   + faccl_nsi(idime)                     )     ! a
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
          ! This was mu_rho FACT4X = gpvol(DEF_VECT,igaus) * gpvis(DEF_VECT,igaus) / gpden(DEF_VECT,igaus) 
          FACT4X = gpvol(DEF_VECT,igaus) * gpden(DEF_VECT,igaus) 
          do inode = 1,pnode
             eldtrho(DEF_VECT,inode) = eldtrho(DEF_VECT,inode) + gpsha(DEF_VECT,inode,igaus) * FACT2X
             elmurho(DEF_VECT,inode) = elmurho(DEF_VECT,inode) + gpsha(DEF_VECT,inode,igaus) * FACT4X
          end do
       end do
    else

       T_dtrho(DEF_VECT) = 0.0_rp
       d_dtrho(DEF_VECT) = 0.0_rp
       T_murho(DEF_VECT) = 0.0_rp
       d_murho(DEF_VECT) = 0.0_rp
       do igaus = 1,pgaus
          FACT2X                 = gpvol(DEF_VECT,igaus) / ( gpden(DEF_VECT,igaus) * DTINV_LOCX )
          !! AB fractional step fix FACT4X                 = gpvol(DEF_VECT,igaus) * gpvis(DEF_VECT,igaus) / gpden(DEF_VECT,igaus)
          FACT4X = gpvol(DEF_VECT,igaus) * gpden(DEF_VECT,igaus) 


          T_dtrho(DEF_VECT) = T_dtrho(DEF_VECT) + FACT2X
          T_murho(DEF_VECT) = T_murho(DEF_VECT) + FACT4X
          do inode = 1,pnode
             eldtrho(DEF_VECT,inode) = eldtrho(DEF_VECT,inode) + gpsha(DEF_VECT,inode,igaus)**2 * FACT2X
             elmurho(DEF_VECT,inode) = elmurho(DEF_VECT,inode) + gpsha(DEF_VECT,inode,igaus)**2 * FACT4X
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
#ifndef OPENACCHHH
    !
    ! Obtain tau for coriolis stabilization -- gpstcor
    !
    if( kfl_stab_corio_nsi > 0_ip) then  ! tal ves meter este if dentro
       gppor = 0.0_rp
       call nsi_element_stabilization_corio(&
            pgaus,pnode,chale,gpadv,gpvis,gpden,  &
            gppor,gpstcor)
    else
       gpstcor = 0.0_rp
    end if
#endif
    !----------------------------------------------------------------------
    !
    ! Element matrices 
    !
    !----------------------------------------------------------------------

    elauu(DEF_VECT,:,:)   = 0.0_rp
    elrbu(DEF_VECT,:,:)   = 0.0_rp
    agrau(DEF_VECT,:,:)   = 0.0_rp
    wgrgr(DEF_VECT,:,:,:) = 0.0_rp

#ifndef OPENACCHHH
    if(kfl_fsgrb_nsi==1) elrbu2(DEF_VECT,:,:)  = 0.0_rp
    if(kfl_fsgrb_nsi==1) elrbp(DEF_VECT,:)     = 0.0_rp
#endif

    !
    ! AGRAU = rho * (a.grad) Ni
    ! WGRGR = grad(Ni) . grad(Nj)
    !
    do igaus = 1,pgaus
       do inode = 1,pnode
          do idime = 1,ndime
             agrau(DEF_VECT,inode,igaus) =  agrau(DEF_VECT,inode,igaus) + &
                  &                         gpadv(DEF_VECT,idime,igaus) * gpcar(DEF_VECT,idime,inode,igaus)
          end do
          agrau(DEF_VECT,inode,igaus) =  gpden(DEF_VECT,igaus) * agrau(DEF_VECT,inode,igaus)
          do jnode = 1,pnode
             do idime = 1,ndime
                wgrgr(DEF_VECT,inode,jnode,igaus) = wgrgr(DEF_VECT,inode,jnode,igaus) + &
                     &                              gpcar(DEF_VECT,idime,inode,igaus)*gpcar(DEF_VECT,idime,jnode,igaus)
             end do
          end do
       end do
    end do

    !--------------------------------------------------------------------
    !
    ! Convective term and RHS
    !
    !--------------------------------------------------------------------

    do igaus = 1,pgaus
       FACT0X = gpden(DEF_VECT,igaus) * gpvol(DEF_VECT,igaus)
       FACT2X = 0.0_rp
       do idime = 1,ndime
          FACT2X = FACT2X + gpgve(DEF_VECT,idime,idime,igaus)
       end do
       FACT2X = FACT2X * FACT0X

       do inode = 1,pnode
          FACT4X = gpsha(DEF_VECT,inode,igaus) * gpvol(DEF_VECT,igaus)
          do idime = 1,ndime
             idofv = (inode-1)*ndime+idime
             do jnode = 1,pnode
                jdofv = (jnode-1)*ndime+idime
                !
                ! rho (u.grad) u
                !
                elauu(DEF_VECT,idofv,jdofv) = elauu(DEF_VECT,idofv,jdofv) + FACT4X * agrau(DEF_VECT,jnode,igaus)
                !
                ! rho * (div u) u
                !
                elauu(DEF_VECT,idofv,jdofv) = elauu(DEF_VECT,idofv,jdofv) &
                     + FACT2X * gpsha(DEF_VECT,inode,igaus) * gpsha(DEF_VECT,jnode,igaus)
                !
                ! rho * u.grad(u)^t
                !
                do jdime = 1,ndime
                   jdofv = (jnode-1)*ndime+jdime
                   elauu(DEF_VECT,idofv,jdofv) = elauu(DEF_VECT,idofv,jdofv) &
                        + FACT0X * gpsha(DEF_VECT,inode,igaus) * gpcar(DEF_VECT,idime,jnode,igaus) * gpvel(DEF_VECT,jdime,igaus)
                end do
             end do
             ! 
             ! 0.5 * grad(u**2)
             ! test oriol
             do jdime = 1,ndime
                elrbu(DEF_VECT,idime,inode) = elrbu(DEF_VECT,idime,inode) &
                     +  FACT0X * gpsha(DEF_VECT,inode,igaus) * gpvel(DEF_VECT,jdime,igaus)*gpgve(DEF_VECT,idime,jdime,igaus)
             end do
          end do
          !
          ! bu = ( f , v ) 
          !
          FACT1X = gpvol(DEF_VECT,igaus) * gpsha(DEF_VECT,inode,igaus)  ! ( f , v )
 
          do idime = 1,ndime 
             elrbu(DEF_VECT,idime,inode)  = elrbu(DEF_VECT,idime,inode)  + FACT1X * gprhs(DEF_VECT,idime,igaus)
          end do


       end do
    end do

#ifndef OPENACCHHH
    if (kfl_dampi_nsi == 1_ip) then

        do igaus = 1,pgaus

            gpcod = 0.0_rp
            !
            ! Here I could just define the y coordinate
            !
            do idime = 1,ndime
               do inode = 1,pnode
                  gpcod(DEF_VECT,idime) = gpcod(DEF_VECT,idime) + gpsha(DEF_VECT,inode,igaus) * elcod(DEF_VECT,idime,inode)
               end do
            end do

            do inode = 1,pnode
                FACT1X = gpvol(DEF_VECT,igaus) * gpsha(DEF_VECT,inode,igaus)  ! ( f , v )
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
    
              FACT1X = gpvol(DEF_VECT,igaus) * gpsha(DEF_VECT,inode,igaus)  ! ( f , v )
    
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
    !
    ! Send only convective term to RHS
    !
    if( NSI_ASSEMBLY_VISCOUS ) then
       do jnode = 1,pnode
          do jdime = 1,ndime
             jevat = (jnode-1)*ndime+jdime
             do inode = 1,pnode
                do idime = 1,ndime
                   ievat = (inode-1)*ndime+idime
                   elrbu(DEF_VECT,idime,inode) = elrbu(DEF_VECT,idime,inode) &
                        - elauu(DEF_VECT,ievat,jevat) * elvel(DEF_VECT,jdime,jnode)
                   elauu(DEF_VECT,ievat,jevat) = 0.0_rp
                end do
             end do
          end do
       end do
    end if

    !--------------------------------------------------------------------
    !
    ! Viscous term
    !
    !--------------------------------------------------------------------

    do igaus = 1,pgaus
       FACT6X = gpvis(DEF_VECT,igaus) * gpvol(DEF_VECT,igaus) 
       do inode = 1,pnode
          do idime = 1,ndime
             idofv = (inode-1)*ndime+idime
             do jnode = 1,pnode
                jdofv           = (jnode-1)*ndime+idime
                elauu(DEF_VECT,idofv,jdofv) = elauu(DEF_VECT,idofv,jdofv) + FACT6X * wgrgr(DEF_VECT,inode,jnode,igaus)
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


#ifndef OPENACCHHH

    if( kfl_stab_corio_nsi > 0_ip ) then  ! i could perhaps unify this better
       !
       ! Coriolis stabilization --  ( tau_cor * div(omega X u) , div(omega X v) )   -- also galerking
       !
       call nsi_corio_stab_element_operations(pnode,pgaus,gpvol,gpden,gpsha,gpcar,gpstcor,elprdivcor,elrbu,elauu)
    else
       !
       ! just Coriolis galerkin term - adapted from whta I had in nsi_corio_stab_element_operations -is170
       !
       call nsi_corio_element_operations(pnode,pgaus,gpvol,gpden,gpsha,elauu)
    end if    
#endif

    !
    ! Send both convective and viscous terms to RHS  also corio if present
    !
    if( .not. NSI_ASSEMBLY_VISCOUS ) then
       do jnode = 1,pnode
          do jdime = 1,ndime
             jevat = (jnode-1)*ndime+jdime
             do inode = 1,pnode
                do idime = 1,ndime
                   ievat = (inode-1)*ndime+idime
                   elrbu(DEF_VECT,idime,inode) = elrbu(DEF_VECT,idime,inode) &
                        - elauu(DEF_VECT,ievat,jevat) * elvel(DEF_VECT,jdime,jnode)
                end do
             end do
          end do
       end do
    end if


    !--------------------------------------------------------------------
    !
    ! Assembly in global syste,
    !
    !--------------------------------------------------------------------        
    !
    ! Scatter element matrix to global one 
    ! 
#ifndef OPENACCHHH
    call cputim(timeb)
    time1(7) = time1(7) + timeb - timea
    call cputim(timea)
    do ivect = 1,VECTOR_DIM                        
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
       !call nsi_assembly_semi_implicit_method(&
       !     pnode,list_elements,lnods_loc,elauu,visco_nsi)!amatr(poauu_nsi:))
    end if


    !$acc update self(visco_nsi)

    !$acc exit data delete( agrau, lnods_loc , gpvol, elauu ,             &   
    !$acc              eldtrho , elvel     ,  gpadv , gpvel ,             &
    !$acc              gpgve   , elcod     ,  gpcar ,                     &
    !$acc              wgrgr   , gpdet     ,                              &
    !$acc              xjacm   , xjaci     , elwvi  ,                     &
    !$acc              T_dtrho , d_dtrho   , T_murho, d_murho,            &
    !$acc              gprhs   , elmurho   , elrbu  , list_elements,      &
    !$acc              xforc_material_nsi  , lforc_material_nsi,          &
    !$acc              gpsha , list_elements_p,                           &
    !$acc              gpden ,                                            &
    !$acc              gpvis ,                                            &
    !$acc              gpmut ,                                            &
    !$acc              gpwvi ,                                            &
    !$acc              elmar,                                             &
    !$acc              elmar(pelty)%shape ,                               &
    !$acc              elmar(pelty)%deriv ,                               &
    !$acc              elmar(pelty)%weigp, visco_nsi                      )

#ifndef OPENACCHHH
    call cputim(timeb)
    time1(8) = time1(8) + timeb - timea                       
#endif

  end subroutine nsi_element_operations_fast_dev


  !-----------------------------------------------------------------------
  !>
  !> @author  Herbert Owen
  !> @date    2020-01-10
  !> @brief   Contribution to elauu -  Coriolis galerkin  
  !> @details Contribution to elauu -  Coriolis galerkin - Adapted from is170_coriolis_stab nsi_corio_stab_element_operations
  !>
  !-----------------------------------------------------------------------

  subroutine nsi_corio_element_operations(pnode,pgaus,gpvol,gpden,gpsha,elauu)
    
    use def_nastin,                     only : fvela_nsi,corio_nsi
    implicit none

    integer(ip),intent(in)                                 :: pnode,pgaus
    real(rp),intent(in)                                    :: gpvol(VECTOR_SIZE,pgaus)
    real(rp),intent(in)                                    :: gpden(VECTOR_SIZE,pgaus)
    real(rp),intent(in)                                    :: gpsha(VECTOR_SIZE,pnode,pgaus)
    real(rp),intent(inout)                                 :: elauu(VECTOR_SIZE,ndime*pnode,ndime*pnode)
    !yor! for memory padding purpose : 4*((pnode+3)/4) will give the closest number >= pnode which is a multiple of 4.
    real(rp), dimension(VECTOR_SIZE,4*((pnode+3)/4),pnode) :: elauu12g,elauu23g,elauu31g   ! I will treat the galerkin part separatelly so that I can take advange that the stab is symmetric
    real(rp), dimension(VECTOR_SIZE,4*((pnode+3)/4),pnode) :: elauu21g,elauu32g,elauu13g   ! that I can take advange that the stabilization is symmetric - rethink if is best for speed
    real(rp)                                               :: fact0(VECTOR_SIZE),fact1(VECTOR_SIZE),fact2(VECTOR_SIZE),fact3(VECTOR_SIZE)
    integer(ip)                                            :: igaus,inode,jnode

    real(rp)                                               :: rmom2(VECTOR_SIZE,ndime,ndime,pnode)
    real(rp)                                               :: factvec1(VECTOR_SIZE,4*((pnode+3)/4))

#define DEF_VECT_IN 1:VECTOR_SIZE

!@#ifndef OPENACCHHH
!@#define DEF_VECT ivect
!@#else
!@#define DEF_VECT 1:VECTOR_SIZE
!@#endif
    !
    ! Coriolis Galerkin term
    !
    if( corio_nsi > 1.0e-12_rp .and. ndime == 3_ip ) then ! 2d not implemented
       elauu12g(DEF_VECT_IN,:,:) = 0.0_rp
       elauu23g(DEF_VECT_IN,:,:) = 0.0_rp
       elauu31g(DEF_VECT_IN,:,:) = 0.0_rp
       elauu21g(DEF_VECT_IN,:,:) = 0.0_rp
       elauu32g(DEF_VECT_IN,:,:) = 0.0_rp
       elauu13g(DEF_VECT_IN,:,:) = 0.0_rp

       do igaus = 1,pgaus
          !
          ! Coriolis force - Galerkin term
          !
          ! 2*rho*(w x u)
          ! x-equation: w x u = wy*uz - wz*uy
          ! y-equation: w x u = wz*ux - wx*uz
          ! z-equation: w x u = wx*uy - wy*ux
          ! Borrowed from nsi_element_operations - here recalculated for each igaus
          !
          rmom2 = 0.0_rp
          fact0(DEF_VECT_IN) = 2.0_rp * gpden(DEF_VECT_IN,igaus)  * gpvol(DEF_VECT_IN,igaus)   ! este 2 tengo que revisarlo !!! en otras partes creo no esta
          fact1(DEF_VECT_IN) = fact0(DEF_VECT_IN) * fvela_nsi(1)
          fact2(DEF_VECT_IN) = fact0(DEF_VECT_IN) * fvela_nsi(2)
          fact3(DEF_VECT_IN) = fact0(DEF_VECT_IN) * fvela_nsi(3)
          do inode = 1,pnode
             rmom2(DEF_VECT_IN,1,2,inode) = - fact3(DEF_VECT_IN) * gpsha(DEF_VECT_IN,inode,igaus)  ! -wz*uy
             rmom2(DEF_VECT_IN,1,3,inode) =   fact2(DEF_VECT_IN) * gpsha(DEF_VECT_IN,inode,igaus)  !  wy*uz
             rmom2(DEF_VECT_IN,2,1,inode) =   fact3(DEF_VECT_IN) * gpsha(DEF_VECT_IN,inode,igaus)  !  wz*ux
             rmom2(DEF_VECT_IN,2,3,inode) = - fact1(DEF_VECT_IN) * gpsha(DEF_VECT_IN,inode,igaus)  ! -wx*uz
             rmom2(DEF_VECT_IN,3,1,inode) = - fact2(DEF_VECT_IN) * gpsha(DEF_VECT_IN,inode,igaus)  ! -wy*ux
             rmom2(DEF_VECT_IN,3,2,inode) =   fact1(DEF_VECT_IN) * gpsha(DEF_VECT_IN,inode,igaus)  !  wx*uy
          end do
          !
          ! Borrowed from nsi_element_assembly
          !
          factvec1(DEF_VECT_IN,1:pnode) = gpsha(DEF_VECT_IN,1:pnode,igaus)
          do jnode = 1,pnode
             do inode = 1,pnode
                elauu21g(DEF_VECT_IN,inode,jnode) = elauu21g(DEF_VECT_IN,inode,jnode) + factvec1(DEF_VECT_IN,inode)*rmom2(DEF_VECT_IN,2,1,jnode)
                elauu31g(DEF_VECT_IN,inode,jnode) = elauu31g(DEF_VECT_IN,inode,jnode) + factvec1(DEF_VECT_IN,inode)*rmom2(DEF_VECT_IN,3,1,jnode)
                elauu12g(DEF_VECT_IN,inode,jnode) = elauu12g(DEF_VECT_IN,inode,jnode) + factvec1(DEF_VECT_IN,inode)*rmom2(DEF_VECT_IN,1,2,jnode)
                elauu32g(DEF_VECT_IN,inode,jnode) = elauu32g(DEF_VECT_IN,inode,jnode) + factvec1(DEF_VECT_IN,inode)*rmom2(DEF_VECT_IN,3,2,jnode)
                elauu13g(DEF_VECT_IN,inode,jnode) = elauu13g(DEF_VECT_IN,inode,jnode) + factvec1(DEF_VECT_IN,inode)*rmom2(DEF_VECT_IN,1,3,jnode)
                elauu23g(DEF_VECT_IN,inode,jnode) = elauu23g(DEF_VECT_IN,inode,jnode) + factvec1(DEF_VECT_IN,inode)*rmom2(DEF_VECT_IN,2,3,jnode)
             end do
          end do
       end do  !igauss

       do jnode = 1,pnode
          do inode = 1,pnode
             elauu(DEF_VECT_IN,3*inode-1,3*jnode-2) = elauu(DEF_VECT_IN,3*inode-1,3*jnode-2) + elauu21g(DEF_VECT_IN,inode,jnode)
             elauu(DEF_VECT_IN,3*inode  ,3*jnode-2) = elauu(DEF_VECT_IN,3*inode  ,3*jnode-2) + elauu31g(DEF_VECT_IN,inode,jnode)
             elauu(DEF_VECT_IN,3*inode-2,3*jnode-1) = elauu(DEF_VECT_IN,3*inode-2,3*jnode-1) + elauu12g(DEF_VECT_IN,inode,jnode)
             elauu(DEF_VECT_IN,3*inode  ,3*jnode-1) = elauu(DEF_VECT_IN,3*inode  ,3*jnode-1) + elauu32g(DEF_VECT_IN,inode,jnode)
             elauu(DEF_VECT_IN,3*inode-2,3*jnode  ) = elauu(DEF_VECT_IN,3*inode-2,3*jnode  ) + elauu13g(DEF_VECT_IN,inode,jnode)
             elauu(DEF_VECT_IN,3*inode-1,3*jnode  ) = elauu(DEF_VECT_IN,3*inode-1,3*jnode  ) + elauu23g(DEF_VECT_IN,inode,jnode)
          end do
       end do
    end if

  end subroutine nsi_corio_element_operations

end module mod_nsi_element_operations_fast_dev
!> @}
