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

module mod_nsi_element_operations_test

#include "def_vector_size.inc"
  use def_kintyp,                     only : ip,rp   !in_const
  use def_domain,                     only : ndime  !in_const_scalar  ! actually paarmeter if ndimepar
#ifndef SUPER_FAST
  use mod_element_integration,        only : element_shape_function_derivatives_jacobian
#endif

#ifdef _OPENACC
  use openacc
#endif

  use def_master,                     only : kfl_paral
  
  implicit none
  private
  
  public :: nsi_element_operations_fast3
  public :: nsi_element_operations_fast6
  public :: nsi_element_operations_fast8



contains

  subroutine nsi_element_operations_fast3(&
       pnode,pgaus,list_elements,time1)

    use def_kintyp,            only : ip,rp  ! in_const_scalar
    use def_master,            only : rhsid  ! out       ! real(rp), pointer     :: rhsid(:)
    use def_master,            only : veloc  ! in_var    ! real(rp), pointer     :: veloc(:,:,:)
    use def_domain,            only : coord  ! in_const  ! real(rp), pointer     :: coord(:,:)
    use def_domain,            only : ltype  ! in_const  ! integer(ip), pointer  :: ltype(:) -- type of element
    use def_domain,            only : lnods  ! in_const  ! integer(ip), pointer  :: lnods(:,:)
    use def_domain,            only : mnode  ! in_const scalar
#ifndef SUPER_FAST
    use def_domain,            only : elmar
    use mod_ker_proper,        only : ker_proper
#endif
    use def_domain,            only : ntens        ! in_const scalar
    use def_nastin,            only : dtinv_nsi    ! in_var   scalar
    use def_nastin,            only : dt_rho_nsi   ! out  ! real(rp), pointer     :: dt_rho_nsi(:)
    use def_nastin,            only : mass_rho_nsi ! out  ! real(rp), pointer     :: mass_rho_nsi(:,:)
    use def_nastin,            only : grnor_nsi    ! in_const_scalar
    use def_nastin,            only : gravi_nsi    ! in_const  ! real(rp)         ::gravi_nsi(3)
    !
    ! For skews this part will be eliminated when guillaume uploads dirichlet algorithmic
#ifdef SUPER_FAST
    use def_kermod, only       :  densi_ker,visco_ker  ! in_const type(typ_valpr_ker), target   :: visco_ker
    ! to avoid problems we usa a real densi_aux = densi_ker % rlaws(1,1)
    ! once it works we could revert this to see if derived data types work fine
#endif

    implicit none
    !
    ! Input and output variables
    !
    integer(ip), intent(in)          :: pnode                        !< Number of nodes
    integer(ip), intent(in)          :: pgaus                        !< Number of Gauss points
    integer(ip), intent(in)          :: list_elements(VECTOR_SIZE)   !< List of elements

    real(rp),    intent(inout)       :: time1(10)                    ! Timings
    !
    ! Element matrices and vectors (stiffness and preconditioner)
    !
    real(rp)    :: elauu(VECTOR_SIZE,pnode*ndime,pnode*ndime)        ! Auu
    real(rp)    :: elrbu(VECTOR_SIZE,ndime,pnode)                    ! bu
    real(rp)    :: eldtrho(VECTOR_SIZE,pnode)                        ! Projection of rho/dt
    real(rp)    :: elmurho(VECTOR_SIZE,pnode)                        ! Projection of mu/rho
    !
    ! Gather
    !
    real(rp)    :: elvel(VECTOR_SIZE,ndime,pnode)                       ! u
    real(rp)    :: elcod(VECTOR_SIZE,ndime,pnode)                         ! x
    !
    ! Indices and dimensions
    !
    integer(ip) :: ielem,inode,ivect
    integer(ip) :: pevat
    integer(ip) :: pelty,plapl
    integer(ip) :: ipoin,igaus
    integer(ip) :: lnods_loc(VECTOR_SIZE,pnode)
    integer(ip) :: list_elements_p(VECTOR_SIZE)                      ! List of elements (always positive)
    !
    ! Gauss point values
    !
    real(rp)    :: gpsha(VECTOR_SIZE,pnode,pgaus)                    ! N
    real(rp)    :: gpder(VECTOR_SIZE,ndime,pnode,pgaus)              ! dN/dsi
    real(rp)    :: gpcar(VECTOR_SIZE,ndime,mnode,pgaus)              ! dN/dxi
    real(rp)    :: gphes(VECTOR_SIZE,ntens,mnode,pgaus)              ! d2N/dxidxj
    real(rp)    :: gpvol(VECTOR_SIZE,pgaus)                          ! w*|J|, |J|
    real(rp)    :: gpvis(VECTOR_SIZE,pgaus)                          ! Viscosity
    real(rp)    :: gpmut(VECTOR_SIZE,pgaus)                          ! mut
    real(rp)    :: gpden(VECTOR_SIZE,pgaus)                          ! Density
    real(rp)    :: gpadv(VECTOR_SIZE,ndime,pgaus)                    ! u+u'
    real(rp)    :: gprhs(VECTOR_SIZE,ndime,pgaus)                    ! RHS
    real(rp)    :: gpvel(VECTOR_SIZE,ndime,pgaus)                    ! u
    real(rp)    :: gpgve(VECTOR_SIZE,ndime,ndime,pgaus)              ! grad(u)
    !
    ! Internal
    !
#ifdef OPENACCHHH

#define FACT0X     fact0
#define FACT1X     fact1
#define FACT2X     fact2
#define FACT4X     fact4
#define FACT5X     fact5
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
#define FACT5X     fact5(1:VECTOR_SIZE)
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
    real(rp)    :: FACT5X
    real(rp)    :: FACT6X
    real(rp)    :: DTINV_LOCX
#ifdef SUPER_FAST
    real(rp)    :: T1X
    real(rp)    :: T2X
    real(rp)    :: T3X
    real(rp)    :: DENOMX
#endif

    real(rp)    :: timea,timeb
    integer(ip) :: idime,jdime,jnode,idofv,jdofv,ievat,jevat,iauxi

    ! Local arrays - previously in element_assembly
    real(rp)    :: wgrgr(VECTOR_SIZE,pnode,pnode,pgaus)
    real(rp)    :: agrau(VECTOR_SIZE,pnode,pgaus)
    !
    ! For shape when SUPER_FAST
    !
#ifdef SUPER_FAST
    real(rp)                             :: weigp(4)
    real(rp)                             :: gpdet(VECTOR_SIZE,pgaus)
    real(rp)                             :: densi_aux
    real(rp)                             :: visco_aux
    real(rp)                             :: xjaci(VECTOR_SIZE,ndime,ndime,pgaus)
    real(rp),parameter                   :: epsilgeo_div = epsilon(1.0_rp) !epsilgeo usad to avoid divisions by zero
    real(rp)                             :: sha_aux(4,4) = reshape((/ 0.585410196624969D+00 ,  0.138196601125011D+00 , &
         0.138196601125011D+00 ,  0.138196601125011D+00 ,  0.138196601125010D+00 ,  &
         0.585410196624969D+00 ,  0.138196601125011D+00 ,  0.138196601125011D+00 ,  &
         0.138196601125010D+00 ,  0.138196601125011D+00 ,  0.585410196624969D+00 , &
         0.138196601125011D+00 ,  0.138196601125010D+00 ,&
         0.138196601125011D+00 ,  0.138196601125011D+00 ,  0.585410196624969D+00/), (/4,4/))
#endif


#ifndef SUPER_FAST
    integer(ip) :: dummi
#endif


#ifdef OPENACCHHH
#define DEF_VECT ivect
#else
#define DEF_VECT 1:VECTOR_SIZE
#endif

    !--------------------------------------------------------------------
    !
    ! Gather: global to local
    !
    !--------------------------------------------------------------------
#ifdef SUPER_FAST
    densi_aux = densi_ker % rlaws(1,1)
    visco_aux = visco_ker % rlaws(1,1)
#endif

    call cputim(timea)
    DTINV_LOCX = dtinv_nsi
    ielem = list_elements(1)
    pelty = abs(ltype(ielem))
    pevat = ndime * pnode
    plapl = 0
#ifdef SUPER_FAST
    weigp(1:4)= 1.0_rp/24.0_rp
#endif

    list_elements_p = list_elements

    !$acc enter data create(agrau , lnods_loc , gpvol , elauu ,      &
    !$acc    eldtrho , elvel , gpmut     , gpadv , gpvel ,      &
    !$acc    gpgve   , elcod , gpsha     , gpcar ,              &
    !$acc    gpden   , xjaci , wgrgr     ,                      &
    !$acc    gpdet   , gprhs , elmurho   , elrbu , gpvis     )  &
    !$acc    copyin( weigp  , list_elements , sha_aux )


#ifndef OPENACCHHH
    call cputim(timea)
#endif

    !
    ! This do starts here both in the openacc version and in the nonopenacc
    ! for the non opencc it ends 30  lines later,
    ! in the openacc case it covers all the subroutine.
    ! Similarly the scatter in the nonopenacc case needs a do ivect
    !
    !$acc parallel loop gang vector default(present) async
    do ivect = 1,VECTOR_SIZE
       ielem = abs(list_elements(ivect))
       if ( ielem /= 0 ) then
          lnods_loc(ivect,1:pnode) = lnods(1:pnode,ielem)
       else
          list_elements_p(ivect)   = list_elements(1)
          lnods_loc(ivect,1:pnode) = 0
       end if

       ielem = list_elements(ivect)
       if ( ielem > 0 ) then
          !
          ! Transient
          !
          do inode = 1,pnode
             ipoin = lnods_loc(ivect,inode)
             do idime = 1,ndime
                elvel(ivect,idime,inode) = veloc(idime,ipoin,1)
                elcod(ivect,idime,inode) = coord(idime,ipoin)
             end do
          end do
       else
          !
          ! Element number is null
          !
          elcod(ivect,:,:)   = 0.0_rp
          elvel(ivect,:,:) = 0.0_rp
       end if

#ifndef OPENACCHHH
    end do
    call cputim(timeb)
    time1(1) = time1(1) + timeb - timea
    call cputim(timea)
#endif


    !--------------------------------------------------------------------
    !
    ! Element shape functions and derivatives and properties
    ! Here I coded the superfast version taht is with tetras inlined and without using ker_proper
    ! The non superfast version is not ready but we should try to recover it at some point at least teh sahpe functions part
    !
    !--------------------------------------------------------------------

#ifdef SUPER_FAST

    do igaus=1,pgaus
       do inode=1,pnode
          gpsha(DEF_VECT,inode,igaus) = sha_aux(inode,igaus)   !elmar(pelty) % shape(inode,igaus)
       end do
    end do
    !
    ! GPCAR (from elmgeo_cartesian_derivatives_jacobian_vector), and GPVOL
    !
    !
    ! 3D P1 element
    !
    gpcar(DEF_VECT,1,1,1) =  elcod(DEF_VECT,1,2)   - elcod(DEF_VECT,1,1)
    gpcar(DEF_VECT,1,2,1) =  elcod(DEF_VECT,1,3)   - elcod(DEF_VECT,1,1)
    gpcar(DEF_VECT,1,3,1) =  elcod(DEF_VECT,1,4)   - elcod(DEF_VECT,1,1)
    gpcar(DEF_VECT,2,1,1) =  elcod(DEF_VECT,2,2)   - elcod(DEF_VECT,2,1)
    gpcar(DEF_VECT,2,2,1) =  elcod(DEF_VECT,2,3)   - elcod(DEF_VECT,2,1)
    gpcar(DEF_VECT,2,3,1) =  elcod(DEF_VECT,2,4)   - elcod(DEF_VECT,2,1)
    gpcar(DEF_VECT,3,1,1) =  elcod(DEF_VECT,3,2)   - elcod(DEF_VECT,3,1)
    gpcar(DEF_VECT,3,2,1) =  elcod(DEF_VECT,3,3)   - elcod(DEF_VECT,3,1)
    gpcar(DEF_VECT,3,3,1) =  elcod(DEF_VECT,3,4)   - elcod(DEF_VECT,3,1)
    T1X          =  gpcar(DEF_VECT,2,2,1) * gpcar(DEF_VECT,3,3,1) - gpcar(DEF_VECT,3,2,1) * gpcar(DEF_VECT,2,3,1)
    T2X          = -gpcar(DEF_VECT,2,1,1) * gpcar(DEF_VECT,3,3,1) + gpcar(DEF_VECT,3,1,1) * gpcar(DEF_VECT,2,3,1)
    T3X          =  gpcar(DEF_VECT,2,1,1) * gpcar(DEF_VECT,3,2,1) - gpcar(DEF_VECT,3,1,1) * gpcar(DEF_VECT,2,2,1)
    gpdet(DEF_VECT,1)     =  gpcar(DEF_VECT,1,1,1) * T1X + gpcar(DEF_VECT,1,2,1) * T2X + gpcar(DEF_VECT,1,3,1) * T3X

    DENOMX       =  1.0_rp / (sign(1.0_rp,gpdet(DEF_VECT,1))*max(abs(gpdet(DEF_VECT,1)),epsilgeo_div))

    xjaci(DEF_VECT,1,1,1) =  T1X * DENOMX
    xjaci(DEF_VECT,2,1,1) =  T2X * DENOMX
    xjaci(DEF_VECT,3,1,1) =  T3X * DENOMX
    xjaci(DEF_VECT,2,2,1) = ( gpcar(DEF_VECT,1,1,1) * gpcar(DEF_VECT,3,3,1) - gpcar(DEF_VECT,3,1,1) * gpcar(DEF_VECT,1,3,1)) * DENOMX
    xjaci(DEF_VECT,3,2,1) = (-gpcar(DEF_VECT,1,1,1) * gpcar(DEF_VECT,3,2,1) + gpcar(DEF_VECT,1,2,1) * gpcar(DEF_VECT,3,1,1)) * DENOMX
    xjaci(DEF_VECT,3,3,1) = ( gpcar(DEF_VECT,1,1,1) * gpcar(DEF_VECT,2,2,1) - gpcar(DEF_VECT,2,1,1) * gpcar(DEF_VECT,1,2,1)) * DENOMX
    xjaci(DEF_VECT,1,2,1) = (-gpcar(DEF_VECT,1,2,1) * gpcar(DEF_VECT,3,3,1) + gpcar(DEF_VECT,3,2,1) * gpcar(DEF_VECT,1,3,1)) * DENOMX
    xjaci(DEF_VECT,1,3,1) = ( gpcar(DEF_VECT,1,2,1) * gpcar(DEF_VECT,2,3,1) - gpcar(DEF_VECT,2,2,1) * gpcar(DEF_VECT,1,3,1)) * DENOMX
    xjaci(DEF_VECT,2,3,1) = (-gpcar(DEF_VECT,1,1,1) * gpcar(DEF_VECT,2,3,1) + gpcar(DEF_VECT,2,1,1) * gpcar(DEF_VECT,1,3,1)) * DENOMX

    gpcar(DEF_VECT,1,1,1) = -xjaci(DEF_VECT,1,1,1) - xjaci(DEF_VECT,2,1,1) - xjaci(DEF_VECT,3,1,1)
    gpcar(DEF_VECT,1,2,1) =  xjaci(DEF_VECT,1,1,1)
    gpcar(DEF_VECT,1,3,1) =  xjaci(DEF_VECT,2,1,1)
    gpcar(DEF_VECT,1,4,1) =  xjaci(DEF_VECT,3,1,1)
    gpcar(DEF_VECT,2,1,1) = -xjaci(DEF_VECT,1,2,1) - xjaci(DEF_VECT,2,2,1) - xjaci(DEF_VECT,3,2,1)
    gpcar(DEF_VECT,2,2,1) =  xjaci(DEF_VECT,1,2,1)
    gpcar(DEF_VECT,2,3,1) =  xjaci(DEF_VECT,2,2,1)
    gpcar(DEF_VECT,2,4,1) =  xjaci(DEF_VECT,3,2,1)
    gpcar(DEF_VECT,3,1,1) = -xjaci(DEF_VECT,1,3,1) - xjaci(DEF_VECT,2,3,1) - xjaci(DEF_VECT,3,3,1)
    gpcar(DEF_VECT,3,2,1) =  xjaci(DEF_VECT,1,3,1)
    gpcar(DEF_VECT,3,3,1) =  xjaci(DEF_VECT,2,3,1)
    gpcar(DEF_VECT,3,4,1) =  xjaci(DEF_VECT,3,3,1)

    do igaus = 2,pgaus
       gpdet(DEF_VECT,igaus) = gpdet(DEF_VECT,1)
       do idime = 1,3
          do jdime = 1,3
             xjaci(DEF_VECT,idime,jdime,igaus) = xjaci(DEF_VECT,idime,jdime,1)
          end do
          do inode = 1,4
             gpcar(DEF_VECT,idime,inode,igaus) = gpcar(DEF_VECT,idime,inode,1)
          end do
       end do
    end do

    do igaus = 1,pgaus
       gpvol(DEF_VECT,igaus) = weigp(igaus) * gpdet(DEF_VECT,igaus)
    end do
    !
    ! No Hessian
    !
    gphes(DEF_VECT,1:ntens,1:pnode,1:pgaus) = 0.0_rp


    !--------------------------------------------------------------------
    !
    ! Properties and local time step
    !
    !--------------------------------------------------------------------

    gpden(DEF_VECT,:) = densi_aux
    gpvis(DEF_VECT,:) = visco_aux
    gpmut(DEF_VECT,:) = 0.0_rp

    gpmut(DEF_VECT,:) = gpden(DEF_VECT,:) * gpmut(DEF_VECT,:)
    gpvis(DEF_VECT,:) = gpvis(DEF_VECT,:) + gpmut(DEF_VECT,:)  ! Effective viscosity <= mu+mut

#else
! NON SUPERFAST part missing

#ifdef OPENACCHHH
    call runend('nsi_element_operations_fast: non superfast version not ready for openacc case')
#endif

    call element_shape_function_derivatives_jacobian(&
         pnode,pgaus,plapl,elmar(pelty) % weigp,elmar(pelty) % shape,&
         elmar(pelty) % deriv,elmar(pelty) % heslo,&
         elcod(DEF_VECT,:,:),gpvol(DEF_VECT,:),gpsha(DEF_VECT,:,:),gpder(DEF_VECT,:,:,:),gpcar(DEF_VECT,:,:,:),gphes(DEF_VECT,:,:,:),&
         list_elements(DEF_VECT))

    !--------------------------------------------------------------------
    !
    ! Properties and local time step
    !
    !--------------------------------------------------------------------

    call ker_proper('DENSI','PGAUS',dummi,list_elements_p(DEF_VECT),gpden(DEF_VECT,:),  &
         pnode,pgaus,gpsha(DEF_VECT,:,:),gpcar(DEF_VECT,:,:,:))     ! rho
    call ker_proper('VISCO','PGAUS',dummi,list_elements_p(DEF_VECT),gpvis(DEF_VECT,:),  &
         pnode,pgaus,gpsha(DEF_VECT,:,:),gpcar(DEF_VECT,:,:,:))     ! mu
    call ker_proper('TURBU','PGAUS',dummi,list_elements_p(DEF_VECT),gpmut(DEF_VECT,:),  &
         pnode,pgaus,gpsha(DEF_VECT,:,:),gpcar(DEF_VECT,:,:,:))     ! mut

    gpmut = gpden * gpmut
    gpvis = gpvis + gpmut  ! Effective viscosity <= mu+mut

#endif
! SUPERFAST

    !----------------------------------------------------------------------
    !
    ! Gauss point values & eldtrho, elmurho
    !
    !----------------------------------------------------------------------
    gpgve(DEF_VECT,:,:,:)   = 0.0_rp
    gprhs(DEF_VECT,:,:)     = 0.0_rp
    gpvel(DEF_VECT,:,:) = 0.0_rp
    eldtrho(DEF_VECT,:) = 0.0_rp
    elmurho(DEF_VECT,:) = 0.0_rp

    do igaus = 1,pgaus
       FACT2X =  gpden(DEF_VECT,igaus)  * grnor_nsi
       do idime = 1,ndime
          do inode = 1,pnode
             gpvel(DEF_VECT,idime,igaus) = gpvel(DEF_VECT,idime,igaus) + elvel(DEF_VECT,idime,inode) * gpsha(DEF_VECT,inode,igaus)

             do jdime = 1,ndime
                gpgve(DEF_VECT,jdime,idime,igaus) = gpgve(DEF_VECT,jdime,idime,igaus) &
                     + gpcar(DEF_VECT,jdime,inode,igaus) * elvel(DEF_VECT,idime,inode)
             end do

          end do
          gpadv(DEF_VECT,idime,igaus) = gpvel(DEF_VECT,idime,igaus)
          gprhs(DEF_VECT,idime,igaus) = gprhs(DEF_VECT,idime,igaus) + FACT2X * gravi_nsi(idime)
       end do
       FACT2X = gpvol(DEF_VECT,igaus) / ( gpden(DEF_VECT,igaus) * DTINV_LOCX )
       !FACT4X = gpvol(DEF_VECT,igaus) * gpvis(DEF_VECT,igaus) / ( gpden(DEF_VECT,igaus)  )
       FACT4X = gpvol(DEF_VECT,igaus) * gpden(DEF_VECT,igaus) 
       do inode = 1,pnode
          eldtrho(DEF_VECT,inode) = eldtrho(DEF_VECT,inode) + gpsha(DEF_VECT,inode,igaus) * FACT2X
          elmurho(DEF_VECT,inode) = elmurho(DEF_VECT,inode) + gpsha(DEF_VECT,inode,igaus) * FACT4X
       end do
    end do

    !----------------------------------------------------------------------
    !
    ! Element matrices
    !
    !----------------------------------------------------------------------

    elauu(DEF_VECT,:,:)   = 0.0_rp
    elrbu(DEF_VECT,:,:)   = 0.0_rp
    agrau(DEF_VECT,:,:)   = 0.0_rp
    wgrgr(DEF_VECT,:,:,:) = 0.0_rp

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

    !
    !  rho*(a.grad)Nj ) Ni + mu * grad(Ni) . grad(Nj)
    !
    do igaus = 1,pgaus

       FACT6X = gpvis(DEF_VECT,igaus)   * gpvol(DEF_VECT,igaus)

       do inode = 1,pnode
          do idime = 1,ndime

             idofv           = (inode-1)*ndime+idime

             do jnode = 1,pnode
                jdofv           = (jnode-1)*ndime+idime
                FACT4X = gpsha(DEF_VECT,inode,igaus) * gpvol(DEF_VECT,igaus)
                FACT5X = FACT4X * ( agrau(DEF_VECT,jnode,igaus) ) &
                     &           + FACT6X *   wgrgr(DEF_VECT,inode,jnode,igaus)
                elauu(DEF_VECT,idofv,jdofv) = elauu(DEF_VECT,idofv,jdofv) + FACT5X
             end do
          end do
       end do
    end do


    do igaus = 1,pgaus

       FACT0X = gpden(DEF_VECT,igaus) * gpvol(DEF_VECT,igaus)
       FACT2X = 0.0_rp
       do idime = 1,ndime
          FACT2X = FACT2X + gpgve(DEF_VECT,idime,idime,igaus)
       end do
       FACT2X = FACT2X * FACT0X

       do inode = 1,pnode
          do idime = 1,ndime
             idofv = (inode-1)*ndime+idime
             do jnode = 1,pnode
                jdofv = (jnode-1)*ndime+idime
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
             ! 0.5*grad(u**2)
             !
             do jdime = 1,ndime
                elrbu(DEF_VECT,idime,inode) = elrbu(DEF_VECT,idime,inode) &
                     +  FACT0X * gpsha(DEF_VECT,inode,igaus) * gpvel(DEF_VECT,jdime,igaus)*gpgve(DEF_VECT,idime,jdime,igaus)
             end do
          end do
          !
          ! ( mu*duj/dxi , dv/dxj ) (only div form)
          !
          do idime = 1,ndime
             idofv = (inode-1)*ndime + idime
             do jnode = 1,pnode
                FACT1X = gpvis(DEF_VECT,igaus) * gpvol(DEF_VECT,igaus) * gpcar(DEF_VECT,idime,jnode,igaus)
                do jdime = 1,ndime
                   jdofv                       = (jnode-1)*ndime + jdime
                   elauu(DEF_VECT,idofv,jdofv) = elauu(DEF_VECT,idofv,jdofv) + FACT1X * gpcar(DEF_VECT,jdime,inode,igaus)
                end do
             end do
          end do
          !
          ! bu = ( f , v )
          !
          FACT1X = gpvol(DEF_VECT,igaus) * gpsha(DEF_VECT,inode,igaus)  ! ( f , v )
          do idime = 1,ndime
             elrbu(DEF_VECT,idime,inode) = elrbu(DEF_VECT,idime,inode) + FACT1X * gprhs(DEF_VECT,idime,igaus)
          end do
       end do
    end do
    !
    ! Send matrix to RHS
    !
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

    !
    ! Scatter element matrix to global one
    !
#ifndef OPENACCHHH
    call cputim(timeb)
    time1(7) = time1(7) + timeb - timea
    call cputim(timea)
    do ivect = 1,VECTOR_SIZE
#endif
       ielem = list_elements(ivect)
       if ( ielem > 0 ) then
          do inode = 1,pnode
             ipoin = lnods_loc(ivect,inode)
             do idime = 1,ndime
                iauxi = idime + (ipoin-1) * ndime
#ifdef NO_COLORING
                !$OMP ATOMIC
#endif
               !$acc atomic update
                rhsid(iauxi) = rhsid(iauxi) + elrbu(ivect,idime,inode)

             end do
             !$acc atomic update
#ifdef NO_COLORING
             !$OMP ATOMIC
# endif
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
    !$acc wait
    !!$acc end data


#ifndef OPENACCHHH
    call cputim(timeb)
    time1(8) = time1(8) + timeb - timea
#endif

  end subroutine nsi_element_operations_fast3

    subroutine nsi_element_operations_fast6(&
       VECTOR_DIM,pnode,pgaus,list_elements,time1)

    use def_kintyp,            only : ip,rp                                              ! in_const_scalar
    use def_master,            only : rhsid                                              ! out       ! real(rp), pointer     :: rhsid(:)
    use def_master,            only : veloc                                              ! in_var    ! real(rp), pointer     :: veloc(:,:,:)
    use def_domain,            only : coord                                              ! in_const  ! real(rp), pointer     :: coord(:,:)
    use def_domain,            only : ltype                                              ! in_const  ! integer(ip), pointer  :: ltype(:) -- type of element
    use def_domain,            only : lnods                                              ! in_const  ! integer(ip), pointer  :: lnods(:,:)
    use def_domain,            only : mnode                                              ! in_const scalar
    use def_domain,            only : elmar
    use def_domain,            only : kfl_savda
    use def_domain,            only : elmda_gpvol
    use def_domain,            only : elmda_gpcar
    use def_nastin,            only : dtinv_nsi                                          ! in_var   scalar
    use def_nastin,            only : dt_rho_nsi                                         ! out  ! real(rp), pointer     :: dt_rho_nsi(:)
    use def_nastin,            only : mass_rho_nsi                                       ! out  ! real(rp), pointer     :: mu_rho_nsi(:)
    use def_nastin,            only : grnor_nsi                                          ! in_const_scalar
    use def_nastin,            only : gravi_nsi                                          ! in_const  ! real(rp)         ::gravi_nsi(3)
    use mod_ker_proper,        only : ker_proper
    use def_kermod,            only : densi_ker 
    use def_kermod,            only : visco_ker 
    use def_kermod,            only : turmu_ker 
    implicit none
    !
    ! Input and output variables
    !
    integer(ip), intent(in)          :: VECTOR_DIM                                       !< Number of nodes
    integer(ip), intent(in)          :: pnode                                            !< Number of nodes
    integer(ip), intent(in)          :: pgaus                                            !< Number of Gauss points
    integer(ip), intent(in)          :: list_elements(VECTOR_DIM)                        !< List of elements

    real(rp),    intent(inout)       :: time1(10)                                        ! Timings
    !
    ! Element matrices and vectors (stiffness and preconditioner)
    !
    real(rp)                         :: elauu(VECTOR_DIM,pnode*ndime,pnode*ndime)        ! Auu
    real(rp)                         :: elrbu(VECTOR_DIM,ndime,pnode)                    ! bu
    real(rp)                         :: eldtrho(VECTOR_DIM,pnode)                        ! Projection of rho/dt
    real(rp)                         :: elmurho(VECTOR_DIM,pnode)                        ! Projection of mu/rho
    !
    ! Gather
    !
    real(rp)                         :: elvel(VECTOR_DIM,ndime,pnode)                    ! u
    real(rp)                         :: elcod(VECTOR_DIM,ndime,pnode)                    ! x
    !
    ! Indices and dimensions
    !
    integer                          :: ielem,inode,ivect
    integer(ip)                      :: pelty,j,k
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
    real(rp)                         :: gpmut(VECTOR_DIM,pgaus)                          ! mut
    real(rp)                         :: gpden(VECTOR_DIM,pgaus)                          ! Density
    real(rp)                         :: gpadv(VECTOR_DIM,ndime,pgaus)                    ! u+u'
    real(rp)                         :: gprhs(VECTOR_DIM,ndime,pgaus)                    ! RHS
    real(rp)                         :: gpvel(VECTOR_DIM,ndime,pgaus)                    ! u
    real(rp)                         :: gpgve(VECTOR_DIM,ndime,ndime,pgaus)              ! grad(u)
    real(rp)                         :: xjacm(VECTOR_DIM,ndime,ndime)
    real(rp)                         :: xjaci(VECTOR_DIM,ndime,ndime)
    real(rp)                         :: gpdet(VECTOR_DIM,pgaus)
    real(rp)                         :: wgrgr(VECTOR_DIM,pnode,pnode,pgaus)
    real(rp)                         :: agrau(VECTOR_DIM,pnode,pgaus)
    real(rp)                         :: timea,timeb
    integer                          :: idime
    integer(ip)                      :: jdime,jnode,idofv,jdofv
    integer(ip)                      :: ievat,jevat,iauxi
!    integer(ip)                      :: dummi
    !
    ! Internal
    !
#ifdef OPENACCHHH

#define FACT0X     fact0
#define FACT1X     fact1
#define FACT2X     fact2
#define FACT4X     fact4
#define FACT5X     fact5
#define FACT6X     fact6
#define DTINV_LOCX dtinv_loc
#define T1X        t1
#define T2X        t2
#define T3X        t3
#define DENOMX     denom

#else

#define FACT0X     fact0(1:VECTOR_DIM)
#define FACT1X     fact1(1:VECTOR_DIM)
#define FACT2X     fact2(1:VECTOR_DIM)
#define FACT4X     fact4(1:VECTOR_DIM)
#define FACT5X     fact5(1:VECTOR_DIM)
#define FACT6X     fact6(1:VECTOR_DIM)
#define DTINV_LOCX dtinv_loc(1:VECTOR_DIM)
#define T1X        t1(1:VECTOR_DIM)
#define T2X        t2(1:VECTOR_DIM)
#define T3X        t3(1:VECTOR_DIM)
#define DENOMX     denom(1:VECTOR_DIM)

#endif

    real(rp)    :: FACT0X
    real(rp)    :: FACT1X
    real(rp)    :: FACT2X
    real(rp)    :: FACT4X
    real(rp)    :: FACT5X
    real(rp)    :: FACT6X
    real(rp)    :: DTINV_LOCX
    real(rp)    :: T1X
    real(rp)    :: T2X
    real(rp)    :: T3X
    real(rp)    :: DENOMX
    
#ifdef OPENACCHHH
#define DEF_VECT ivect
#else
#define DEF_VECT 1:VECTOR_SIZE
#endif
    
    !--------------------------------------------------------------------
    !
    ! Gather: global to local
    !
    !--------------------------------------------------------------------

    call cputim(timea)
    
    ielem = list_elements(1)
    pelty = abs(ltype(ielem))
    
    do ivect = 1,VECTOR_DIM                      
       ielem = abs(list_elements(ivect))
       if( ielem /= 0 ) then
          list_elements_p(ivect)   = list_elements(ivect)
       else
          list_elements_p(ivect)   = list_elements(1)
       end if
    end do
    do ivect = 1,VECTOR_DIM      
       ielem = list_elements_p(ivect)
       gpsha(ivect,1:pnode,1:pgaus) = elmar(pelty) % shape(1:pnode,1:pgaus)
    end do

    !call ker_proper('DENSI','PGAUS',dummi,list_elements_p,gpden,pnode,pgaus,gpsha,VECTOR_DIM=VECTOR_DIM) ! rho
    !call ker_proper('VISCO','PGAUS',dummi,list_elements_p,gpvis,pnode,pgaus,gpsha,VECTOR_DIM=VECTOR_DIM) ! mu
    !call ker_proper('TURBU','PGAUS',dummi,list_elements_p,gpmut,pnode,pgaus,gpsha,VECTOR_DIM=VECTOR_DIM) ! mut
    
    DTINV_LOCX = dtinv_nsi

    !--------------------------------------------------------------------
    !
    ! Shape function and derivatives
    !
    !--------------------------------------------------------------------
    
    !$acc data create( agrau   , lnods_loc , gpvol     , elauu ,          &
    !$acc              eldtrho , elvel     , gpmut     , gpadv , gpvel ,  &
    !$acc              gpgve   , elcod     , gpsha     , gpcar ,          &
    !$acc              gpden   , wgrgr     , gpdet     ,                  &
    !$acc              xjacm   , xjaci     ,                              &
    !$acc              elmar,                                             &
    !$acc              elmar(pelty)%shape,                                &
    !$acc              elmar(pelty)%deriv,                                &
    !$acc              elmar(pelty)%weigp,                                &
    !$acc              gprhs   , elmurho   , elrbu     , gpvis )          &
    !$acc copyin(      list_elements        ,                             &
    !$acc              gpsha ,                                            &
    !$acc              gpden ,                                            &
    !$acc              gpvis ,                                            &
    !$acc              gpmut ,                                            &
    !$acc              elmar,                                             &
    !$acc              elmar(pelty)%shape ,                               &
    !$acc              elmar(pelty)%deriv ,                               &
    !$acc              elmar(pelty)%weigp                                 )

#ifndef OPENACCHHH
    call cputim(timea)
#endif
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
          do idime = 1,ndime
             elvel(ivect,idime,inode) = veloc(idime,ipoin,1)
             elcod(ivect,idime,inode) = coord(idime,ipoin)
          end do
       end do
       !
       ! Properties
       !
       gpden(ivect,1:pgaus) = densi_ker % value_const(1,1)
       gpvis(ivect,1:pgaus) = visco_ker % value_const(1,1)
       gpmut(ivect,1:pgaus) = turmu_ker % value_ielem(ielem) % a(1:pgaus)       

#ifndef OPENACCHHH
    end do

    call cputim(timeb)
    time1(1) = time1(1) + timeb - timea
    call cputim(timea)
#endif

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
    ! Properties and local time step
    !
    !--------------------------------------------------------------------

    gpmut(DEF_VECT,:) = gpden(DEF_VECT,:) * gpmut(DEF_VECT,:)
    gpvis(DEF_VECT,:) = gpvis(DEF_VECT,:) + gpmut(DEF_VECT,:)  ! Effective viscosity <= mu+mut
    
    !----------------------------------------------------------------------
    !
    ! Gauss point values & eldtrho, elmurho
    !
    !----------------------------------------------------------------------

    gpgve(DEF_VECT,:,:,:) = 0.0_rp
    gprhs(DEF_VECT,:,:)   = 0.0_rp
    gpvel(DEF_VECT,:,:)   = 0.0_rp
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
          gprhs(DEF_VECT,idime,igaus) = gprhs(DEF_VECT,idime,igaus) + FACT2X * gravi_nsi(idime)
       end do
       FACT2X = gpvol(DEF_VECT,igaus) / ( gpden(DEF_VECT,igaus) * DTINV_LOCX )
       FACT4X = gpvol(DEF_VECT,igaus) * gpvis(DEF_VECT,igaus) / ( gpden(DEF_VECT,igaus)  )
       do inode = 1,pnode
          eldtrho(DEF_VECT,inode) = eldtrho(DEF_VECT,inode) + gpsha(DEF_VECT,inode,igaus) * FACT2X
          elmurho(DEF_VECT,inode) = elmurho(DEF_VECT,inode) + gpsha(DEF_VECT,inode,igaus) * FACT4X
       end do
    end do

    !----------------------------------------------------------------------
    !
    ! Element matrices 
    !
    !----------------------------------------------------------------------

    elauu(DEF_VECT,:,:)   = 0.0_rp
    elrbu(DEF_VECT,:,:)   = 0.0_rp
    agrau(DEF_VECT,:,:)   = 0.0_rp
    wgrgr(DEF_VECT,:,:,:) = 0.0_rp
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
    !
    !  rho*(a.grad)Nj ) Ni + mu * grad(Ni) . grad(Nj)
    !
    do igaus = 1,pgaus

       FACT6X = gpvis(DEF_VECT,igaus)   * gpvol(DEF_VECT,igaus)

       do inode = 1,pnode
          do idime = 1,ndime

             idofv           = (inode-1)*ndime+idime

             do jnode = 1,pnode
                jdofv           = (jnode-1)*ndime+idime
                FACT4X = gpsha(DEF_VECT,inode,igaus) * gpvol(DEF_VECT,igaus)
                FACT5X = FACT4X * ( agrau(DEF_VECT,jnode,igaus) ) &
                     &           + FACT6X *   wgrgr(DEF_VECT,inode,jnode,igaus)
                elauu(DEF_VECT,idofv,jdofv) = elauu(DEF_VECT,idofv,jdofv) + FACT5X
             end do
          end do
       end do
    end do

    do igaus = 1,pgaus

       FACT0X = gpden(DEF_VECT,igaus) * gpvol(DEF_VECT,igaus)
       FACT2X = 0.0_rp
       do idime = 1,ndime
          FACT2X = FACT2X + gpgve(DEF_VECT,idime,idime,igaus)
       end do
       FACT2X = FACT2X * FACT0X

       do inode = 1,pnode
          do idime = 1,ndime
             idofv = (inode-1)*ndime+idime
             do jnode = 1,pnode
                jdofv = (jnode-1)*ndime+idime
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
             !
             do jdime = 1,ndime
                elrbu(DEF_VECT,idime,inode) = elrbu(DEF_VECT,idime,inode) &
                     +  FACT0X * gpsha(DEF_VECT,inode,igaus) * gpvel(DEF_VECT,jdime,igaus)*gpgve(DEF_VECT,idime,jdime,igaus)
             end do
          end do
          !
          ! ( mu*duj/dxi , dv/dxj ) (only div form)
          !      
          do idime = 1,ndime
             idofv = (inode-1)*ndime + idime
             do jnode = 1,pnode
                FACT1X = gpvis(DEF_VECT,igaus) * gpvol(DEF_VECT,igaus) * gpcar(DEF_VECT,idime,jnode,igaus)
                do jdime = 1,ndime
                   jdofv                       = (jnode-1)*ndime + jdime
                   elauu(DEF_VECT,idofv,jdofv) = elauu(DEF_VECT,idofv,jdofv) + FACT1X * gpcar(DEF_VECT,jdime,inode,igaus)
                end do
             end do
          end do
          !
          ! bu = ( f , v ) 
          !
          FACT1X = gpvol(DEF_VECT,igaus) * gpsha(DEF_VECT,inode,igaus)  ! ( f , v )
          do idime = 1,ndime
             elrbu(DEF_VECT,idime,inode) = elrbu(DEF_VECT,idime,inode) + FACT1X * gprhs(DEF_VECT,idime,igaus) 
          end do
       end do
    end do
    !
    ! Send matrix to RHS
    !
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
#ifdef NO_COLORING
                !$OMP ATOMIC
#endif
                !$acc atomic update
                rhsid(iauxi) = rhsid(iauxi) + elrbu(ivect,idime,inode)

             end do
#ifdef NO_COLORING
             !$OMP ATOMIC
#endif
             !$acc atomic update
             dt_rho_nsi(ipoin) = dt_rho_nsi(ipoin) + eldtrho(ivect,inode)
#ifdef NO_COLORING
             !$OMP ATOMIC
#endif
             !$acc atomic update
             mass_rho_nsi(ipoin,1) = mass_rho_nsi(ipoin,1) + elmurho(ivect,inode)

          end do
       end if
    end do
    !$acc end parallel loop
    !$acc end data


#ifndef OPENACCHHH
    call cputim(timeb)
    time1(8) = time1(8) + timeb - timea                       
#endif

  end subroutine nsi_element_operations_fast6


    subroutine nsi_element_operations_fast8(&
       VECTOR_DIM,pnode,pgaus,list_elements,time1,streamid)

    use def_kintyp,            only : ip,rp                                              ! in_const_scalar
    use def_master,            only : rhsid                                              ! out       ! real(rp), pointer     :: rhsid(:)
    use def_master,            only : veloc                                              ! in_var    ! real(rp), pointer     :: veloc(:,:,:)
    use def_domain,            only : coord                                              ! in_const  ! real(rp), pointer     :: coord(:,:)
    use def_domain,            only : ltype                                              ! in_const  ! integer(ip), pointer  :: ltype(:) -- type of element
    use def_domain,            only : lnods                                              ! in_const  ! integer(ip), pointer  :: lnods(:,:)
    use def_domain,            only : mnode                                              ! in_const scalar
    use def_domain,            only : elmar
    use def_domain,            only : kfl_savda
    use def_domain,            only : elmda_gpvol
    use def_domain,            only : elmda_gpcar
    use def_nastin,            only : dtinv_nsi                                          ! in_var   scalar
    use def_nastin,            only : dt_rho_nsi                                         ! out  ! real(rp), pointer     :: dt_rho_nsi(:)
    use def_nastin,            only : mass_rho_nsi                                       ! out  ! real(rp), pointer     :: mu_rho_nsi(:)
    use def_nastin,            only : grnor_nsi                                          ! in_const_scalar
    use def_nastin,            only : gravi_nsi                                          ! in_const  ! real(rp)         ::gravi_nsi(3)
    use mod_ker_proper,        only : ker_proper
    use def_kermod,            only : densi_ker 
    use def_kermod,            only : visco_ker 
    implicit none
    !
    ! Input and output variables
    !
    integer(ip), intent(in)          :: VECTOR_DIM                                       !< Number of nodes
    integer(ip), intent(in)          :: pnode                                            !< Number of nodes
    integer(ip), intent(in)          :: pgaus                                            !< Number of Gauss points
    integer(ip), intent(in)          :: list_elements(VECTOR_DIM)                        !< List of elements

    real(rp),    intent(inout)       :: time1(10)                                        ! Timings
    integer(ip)                      :: streamid
    !
    ! Element matrices and vectors (stiffness and preconditioner)
    !
    real(rp)                         :: elauu(VECTOR_DIM,pnode*ndime,pnode*ndime)        ! Auu
    real(rp)                         :: elrbu(VECTOR_DIM,ndime,pnode)                    ! bu
    real(rp)                         :: eldtrho(VECTOR_DIM,pnode)                        ! Projection of rho/dt
    real(rp)                         :: elmurho(VECTOR_DIM,pnode)                        ! Projection of mu/rho
    !
    ! Gather
    !
    real(rp)                         :: elvel(VECTOR_DIM,ndime,pnode)                    ! u
    real(rp)                         :: elcod(VECTOR_DIM,ndime,pnode)                    ! x
    !
    ! Indices and dimensions
    !
    integer                          :: ielem,inode,ivect
    integer(ip)                      :: pelty,j,k
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
    real(rp)                         :: gpmut(VECTOR_DIM,pgaus)                          ! mut
    real(rp)                         :: gpden(VECTOR_DIM,pgaus)                          ! Density
    real(rp)                         :: gpadv(VECTOR_DIM,ndime,pgaus)                    ! u+u'
    real(rp)                         :: gprhs(VECTOR_DIM,ndime,pgaus)                    ! RHS
    real(rp)                         :: gpvel(VECTOR_DIM,ndime,pgaus)                    ! u
    real(rp)                         :: gpgve(VECTOR_DIM,ndime,ndime,pgaus)              ! grad(u)
    real(rp)                         :: xjacm(VECTOR_DIM,ndime,ndime)
    real(rp)                         :: xjaci(VECTOR_DIM,ndime,ndime)
    real(rp)                         :: gpdet(VECTOR_DIM,pgaus)
    real(rp)                         :: wgrgr(VECTOR_DIM,pnode,pnode,pgaus)
    real(rp)                         :: agrau(VECTOR_DIM,pnode,pgaus)
    real(rp)                         :: timea,timeb
    integer                          :: idime
    integer(ip)                      :: jdime,jnode,idofv,jdofv
    integer(ip)                      :: ievat,jevat,iauxi,dummi
    real(rp)                         :: scaden,scavis   !
    ! Internal
    !
#ifdef OPENACCHHH

#define FACT0X     fact0
#define FACT1X     fact1
#define FACT2X     fact2
#define FACT4X     fact4
#define FACT5X     fact5
#define FACT6X     fact6
#define DTINV_LOCX dtinv_loc
#define T1X        t1
#define T2X        t2
#define T3X        t3
#define DENOMX     denom

#else

#define FACT0X     fact0(1:VECTOR_DIM)
#define FACT1X     fact1(1:VECTOR_DIM)
#define FACT2X     fact2(1:VECTOR_DIM)
#define FACT4X     fact4(1:VECTOR_DIM)
#define FACT5X     fact5(1:VECTOR_DIM)
#define FACT6X     fact6(1:VECTOR_DIM)
#define DTINV_LOCX dtinv_loc(1:VECTOR_DIM)
#define T1X        t1(1:VECTOR_DIM)
#define T2X        t2(1:VECTOR_DIM)
#define T3X        t3(1:VECTOR_DIM)
#define DENOMX     denom(1:VECTOR_DIM)

#endif

    real(rp)    :: FACT0X
    real(rp)    :: FACT1X
    real(rp)    :: FACT2X
    real(rp)    :: FACT4X
    real(rp)    :: FACT5X
    real(rp)    :: FACT6X
    real(rp)    :: DTINV_LOCX
    real(rp)    :: T1X
    real(rp)    :: T2X
    real(rp)    :: T3X
    real(rp)    :: DENOMX
    
#ifdef OPENACCHHH
#define DEF_VECT ivect
#else
#define DEF_VECT 1:VECTOR_SIZE
#endif
    
    !--------------------------------------------------------------------
    !
    ! Gather: global to local
    !
    !--------------------------------------------------------------------

    call cputim(timea)
    scaden=densi_ker % value_const(1,1)  
    scavis=visco_ker % value_const(1,1)
 
    ielem = list_elements(1)
    pelty = abs(ltype(ielem))
    
    do ivect = 1,VECTOR_DIM                      
       ielem = abs(list_elements(ivect))
       if( ielem /= 0 ) then
          list_elements_p(ivect)   = list_elements(ivect)
       else
          list_elements_p(ivect)   = list_elements(1)
       end if
    end do
    do ivect = 1,VECTOR_DIM      
       ielem = list_elements_p(ivect)
       gpsha(ivect,1:pnode,1:pgaus) = elmar(pelty) % shape(1:pnode,1:pgaus)
    end do

 !   call ker_proper('DENSI','PGAUS',dummi,list_elements_p,gpden,pnode,pgaus,gpsha,VECTOR_DIM=VECTOR_DIM) ! rho
!    call ker_proper('VISCO','PGAUS',dummi,list_elements_p,gpvis,pnode,pgaus,gpsha,VECTOR_DIM=VECTOR_DIM) ! mu
    call ker_proper('TURBU','PGAUS',dummi,list_elements_p,gpmut,pnode,pgaus,gpsha,VECTOR_DIM=VECTOR_DIM) ! mut
    
    DTINV_LOCX = dtinv_nsi
#ifdef OPENACCHHH

    !--------------------------------------------------------------------
    !
    ! Shape function and derivatives
    !
    !--------------------------------------------------------------------
    !$acc enter data create( agrau, lnods_loc , gpvol     , elauu ,          &
    !$acc              eldtrho , elvel     ,  gpadv , gpvel ,  &
    !$acc              gpgve   , elcod     ,  gpcar ,          &
    !$acc              wgrgr     , gpdet     ,gpden, gpvis,           &
    !$acc              xjacm   , xjaci     ,                              &
    !$acc              gprhs   , elmurho   , elrbu)          &
    !$acc copyin(      list_elements        ,coord,                             &
    !$acc              gpsha ,                                            &
    !$acc              gpmut ,                                            &
    !$acc              elmar,                                             &
    !$acc              elmar(pelty)%shape ,                               &
    !$acc              elmar(pelty)%deriv ,                               &
    !$acc              elmar(pelty)%weigp                 ) async(streamid)
#endif   

#ifndef OPENACCHHH
    call cputim(timea)
#endif
    !
    ! This do starts here both in the openacc version and in the nonopenacc
    ! for the non opencc it ends 30  lines later,
    ! in the openacc case it covers all the subroutine.
    ! Similarly the scatter in the nonopenacc case needs a do ivect
    !    
#ifdef OPENACCHHH
    !$acc parallel loop gang vector default(present) &
    !$acc        async(streamid)
    !
#endif
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
          do idime = 1,ndime
             elvel(ivect,idime,inode) = veloc(idime,ipoin,1)
             elcod(ivect,idime,inode) = coord(idime,ipoin)
          end do
       end do
       !
       gpden(ivect,1:pgaus) = scaden 
       gpvis(ivect,1:pgaus) = scavis 
!
#ifndef OPENACCHHH
    end do

    call cputim(timeb)
    time1(1) = time1(1) + timeb - timea
    call cputim(timea)
#endif

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

!!! ----
    !--------------------------------------------------------------------
    !
    ! Properties and local time step
    !
    !--------------------------------------------------------------------

    gpmut(DEF_VECT,:) = gpden(DEF_VECT,:) * gpmut(DEF_VECT,:)
    gpvis(DEF_VECT,:) = gpvis(DEF_VECT,:) + gpmut(DEF_VECT,:)  ! Effective viscosity <= mu+mut
    
    !----------------------------------------------------------------------
    !
    ! Gauss point values & eldtrho, elmurho
    !
    !----------------------------------------------------------------------

    gpgve(DEF_VECT,:,:,:) = 0.0_rp
    gprhs(DEF_VECT,:,:)   = 0.0_rp
    gpvel(DEF_VECT,:,:)   = 0.0_rp
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
          gprhs(DEF_VECT,idime,igaus) = gprhs(DEF_VECT,idime,igaus) + FACT2X * gravi_nsi(idime)
       end do
       FACT2X = gpvol(DEF_VECT,igaus) / ( gpden(DEF_VECT,igaus) * DTINV_LOCX )
       FACT4X = gpvol(DEF_VECT,igaus) * gpvis(DEF_VECT,igaus) / ( gpden(DEF_VECT,igaus)  )
       do inode = 1,pnode
          eldtrho(DEF_VECT,inode) = eldtrho(DEF_VECT,inode) + gpsha(DEF_VECT,inode,igaus) * FACT2X
          elmurho(DEF_VECT,inode) = elmurho(DEF_VECT,inode) + gpsha(DEF_VECT,inode,igaus) * FACT4X
       end do
    end do



    !----------------------------------------------------------------------
    !
    ! Element matrices 
    !
    !----------------------------------------------------------------------

    elauu(DEF_VECT,:,:)   = 0.0_rp
    elrbu(DEF_VECT,:,:)   = 0.0_rp
    agrau(DEF_VECT,:,:)   = 0.0_rp
    wgrgr(DEF_VECT,:,:,:) = 0.0_rp
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

    !
    !  rho*(a.grad)Nj ) Ni + mu * grad(Ni) . grad(Nj)
    !
    do igaus = 1,pgaus

       FACT6X = gpvis(DEF_VECT,igaus)   * gpvol(DEF_VECT,igaus)

       do inode = 1,pnode
          do idime = 1,ndime

             idofv           = (inode-1)*ndime+idime

             do jnode = 1,pnode
                jdofv           = (jnode-1)*ndime+idime
                FACT4X = gpsha(DEF_VECT,inode,igaus) * gpvol(DEF_VECT,igaus)
                FACT5X = FACT4X * ( agrau(DEF_VECT,jnode,igaus) ) &
                     &           + FACT6X *   wgrgr(DEF_VECT,inode,jnode,igaus)
                elauu(DEF_VECT,idofv,jdofv) = elauu(DEF_VECT,idofv,jdofv) + FACT5X
             end do
          end do
       end do
    end do

    do igaus = 1,pgaus

       FACT0X = gpden(DEF_VECT,igaus) * gpvol(DEF_VECT,igaus)
       FACT2X = 0.0_rp
       do idime = 1,ndime
          FACT2X = FACT2X + gpgve(DEF_VECT,idime,idime,igaus)
       end do
       FACT2X = FACT2X * FACT0X

       do inode = 1,pnode
          do idime = 1,ndime
             idofv = (inode-1)*ndime+idime
             do jnode = 1,pnode
                jdofv = (jnode-1)*ndime+idime
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
             !
             do jdime = 1,ndime
                elrbu(DEF_VECT,idime,inode) = elrbu(DEF_VECT,idime,inode) &
                     +  FACT0X * gpsha(DEF_VECT,inode,igaus) * gpvel(DEF_VECT,jdime,igaus)*gpgve(DEF_VECT,idime,jdime,igaus)
             end do
          end do
          !
          ! ( mu*duj/dxi , dv/dxj ) (only div form)
          !      
          do idime = 1,ndime
             idofv = (inode-1)*ndime + idime
             do jnode = 1,pnode
                FACT1X = gpvis(DEF_VECT,igaus) * gpvol(DEF_VECT,igaus) * gpcar(DEF_VECT,idime,jnode,igaus)
                do jdime = 1,ndime
                   jdofv                       = (jnode-1)*ndime + jdime
                   elauu(DEF_VECT,idofv,jdofv) = elauu(DEF_VECT,idofv,jdofv) + FACT1X * gpcar(DEF_VECT,jdime,inode,igaus)
                end do
             end do
          end do
          !
          ! bu = ( f , v ) 
          !
          FACT1X = gpvol(DEF_VECT,igaus) * gpsha(DEF_VECT,inode,igaus)  ! ( f , v )
          do idime = 1,ndime
             elrbu(DEF_VECT,idime,inode) = elrbu(DEF_VECT,idime,inode) + FACT1X * gprhs(DEF_VECT,idime,igaus) 
          end do
       end do
    end do

!!!    !ultima version


    !
    ! Send matrix to RHS
    !
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

#ifdef NO_COLORING
                !$OMP ATOMIC
#endif
 
#ifdef OPENACCHHH
               !$acc atomic update
#endif
                rhsid(iauxi) = rhsid(iauxi) + elrbu(ivect,idime,inode)

             end do
#ifdef NO_COLORING
             !$OMP ATOMIC
#endif
#ifdef OPENACCHHH
            !$acc atomic update
#endif
             dt_rho_nsi(ipoin) = dt_rho_nsi(ipoin) + eldtrho(ivect,inode)
#ifdef NO_COLORING
             !$OMP ATOMIC
#endif
#ifdef OPENACCHHH
           !$acc atomic update
#endif
             mass_rho_nsi(ipoin,1) = mass_rho_nsi(ipoin,1) + elmurho(ivect,inode)

          end do
       end if
    end do
#ifdef OPENACCHHH

    !$acc end parallel loop

    !$acc exit data delete (agrau, lnods_loc , gpvol, elauu, &
    !$acc      eldtrho , elvel     , gpmut     , gpadv , gpvel ,  &
    !$acc      gpgve   , elcod     , gpsha     , gpcar ,          &
    !$acc      gpden   , wgrgr     , gpdet     ,                  &
    !$acc      xjacm   , xjaci     ,                              &
    !$acc      elmar,                                             &
    !$acc      elmar(pelty)%shape,                                &
    !$acc      elmar(pelty)%deriv,                                &
    !$acc      elmar(pelty)%weigp,                                &
    !$acc      gprhs   , elmurho   , elrbu, gpvis,list_elements ) &
    !$acc    async(streamid)

#endif


#ifndef OPENACCHHH
    call cputim(timeb)
    time1(8) = time1(8) + timeb - timea                       
#endif

  end subroutine nsi_element_operations_fast8

end module mod_nsi_element_operations_test
!> @}

