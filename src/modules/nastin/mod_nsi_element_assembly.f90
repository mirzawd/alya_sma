!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



module mod_nsi_element_assembly

  !------------------------------------------------------------------------
  !> @addtogroup NastinMatrixAssembly
  !> @{
  !> @file    nsi_elmmat.f90
  !> @author  Guillaume Houzeaux
  !> @brief   Navier-Stokes system element assembly and other element
  !>          calculations
  !> @details Compute element matrices
  !>
  !>          \verbatim
  !>
  !>          Without enrichement
  !>          -------------------
  !>            +-        +  +- -+     +-  -+
  !>            | Auu Aup |  | u |     | bu |
  !>            |         |  |   |  =  |    |
  !>            | Apu App |  | p |     | bp |
  !>            +-       -+  +- -+     +-  -+
  !>
  !>          With enrichement
  !>          ----------------
  !>            +-            +  +- -+     +-  -+
  !>            | Auu Aup Auq |  | u |     | bu |
  !>            |             |  |   |     |    |
  !>            | Apu App Apq |  | p |  =  | bp |
  !>            |             |  |   |     |    |
  !>            | Aqu Aqp Aqq |  | q |     | bq |
  !>            +-           -+  +- -+     +-  -+
  !>
  !>            q = Aqq^-1 ( bq - Aqu u - Aqp p )
  !>
  !>            Auu <= Auu - Auq Aqq^-1 Aqu
  !>            Aup <= Aup - Auq Aqq^-1 Aqp
  !>            bu  <= bu  - Auq Aqq^-1 bq
  !>
  !>            Apu <= Apu - Apq Aqq^-1 Aqu
  !>            App <= App - Apq Aqq^-1 Aqp
  !>            bp  <= bp  - Apq Aqq^-1 bq
  !>
  !>          \endverbatim
  !>
  !------------------------------------------------------------------------

  use def_kintyp, only :  ip,rp,lg
#include "def_vector_size.inc"
  use def_master, only :  intost
  use def_master, only :  kfl_lumped
  use def_domain, only :  ndime,mnode,ntens
  use def_kermod, only :  kfl_duatss
  use def_kermod, only :  fact_duatss
  use def_nastin, only :  kfl_stabi_nsi
  use def_nastin, only :  fvins_nsi,fcons_nsi,bemol_nsi,kfl_regim_nsi
  use def_nastin, only :  fvela_nsi,kfl_rmom2_nsi,kfl_press_nsi
  use def_nastin, only :  kfl_p1ve2_nsi,kfl_linea_nsi,pabdf_nsi
  use def_nastin, only :  kfl_confi_nsi,nbdfp_nsi,kfl_sgsti_nsi
  use def_nastin, only :  kfl_nota1_nsi,kfl_limit_nsi,kfl_penal_nsi
  use def_nastin, only :  penal_nsi,kfl_immer_nsi
  use def_nastin, only :  kfl_bubbl_nsi
  use def_nastin, only :  kfl_convection_type_nsi
  use def_nastin, only :  NSI_GALERKIN
  use def_nastin, only :  NSI_ASGS
  use def_nastin, only :  NSI_OSS
  use def_nastin, only :  NSI_SPLIT_OSS
  use def_nastin, only :  NSI_ALGEBRAIC_SPLIT_OSS
  use def_nastin, only :  NSI_FRACTIONAL_STEP
  use def_nastin, only :  NSI_CONVECTION_NON_CONSERVATIVE  
  use def_nastin, only :  NSI_CONVECTION_CONSERVATIVE      
  use def_nastin, only :  NSI_CONVECTION_SKEW              
  use def_nastin, only :  NSI_CONVECTION_EMAC
  use def_kermod, only :  kfl_noslw_ker


  implicit none

  real(rp), parameter :: zeror = epsilon(1.0_rp)

  private


  public :: nsi_element_assembly_asgs_oss     ! Assembly using ASGS or OSS
  public :: nsi_element_assembly_split_oss    ! Assembly using split OSS or Galerkin
  public :: nsi_element_assembly_asgs_oss_old ! Old assembly using ASGS or OSS
  public :: nsi_element_system_output         ! Output of element matrix

contains

  !-----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @brief   ASGS and OSS
  !> @details Assembly of Navier Stokes equations using ASGS or OSS
  !>          Variational Multiscale Model.
  !>          Includes subgriscale tracking in time and convection
  !>
  !-----------------------------------------------------------------------

  subroutine nsi_element_assembly_asgs_oss(             &
       pnode,pgaus,gpden,gpvis,gppor,gpgvi,gpsp1,gptt1, &
       gpsp2,gptt2,gpvol,gpsha,gpcar,gplap,gphes,gpadv, &
       gpvep,gprhs,gprhc,rmomu,rcont,elpre,elbub,elauu, &
       elaup,elapp,elapu,elrbu,elrbp,rmom2,gpst1,gpgve, &
       gprh2,gprhs_sgs,elvel,ellum,dtinv_loc,pbubl,     &
       gpsha_bub,gpcar_bub,gppre,elauq,elapq,elaqu,     &
       elaqp,elaqq,elrbq,densi,gplag,gpmix)

    integer(ip), intent(in)    :: pnode,pgaus
    real(rp),    intent(in)    :: gpden(VECTOR_SIZE,pgaus)
    real(rp),    intent(in)    :: gpvis(VECTOR_SIZE,pgaus)
    real(rp),    intent(in)    :: gppor(VECTOR_SIZE,pgaus)
    real(rp),    intent(in)    :: gpgvi(VECTOR_SIZE,ndime,pgaus)
    real(rp),    intent(in)    :: gpsp1(VECTOR_SIZE,pgaus)
    real(rp),    intent(in)    :: gptt1(VECTOR_SIZE,pgaus)
    real(rp),    intent(in)    :: gpsp2(VECTOR_SIZE,pgaus)
    real(rp),    intent(in)    :: gptt2(VECTOR_SIZE,pgaus)
    real(rp),    intent(in)    :: gpvol(VECTOR_SIZE,pgaus)
    real(rp),    intent(in)    :: gpsha(VECTOR_SIZE,pnode,pgaus)
    real(rp),    intent(in)    :: gpcar(VECTOR_SIZE,ndime,mnode,pgaus)
    real(rp),    intent(in)    :: gpadv(VECTOR_SIZE,ndime,pgaus)
    real(rp),    intent(in)    :: gpvep(VECTOR_SIZE,ndime,pgaus)
    real(rp),    intent(in)    :: gplap(VECTOR_SIZE,pnode,pgaus)
    real(rp),    intent(in)    :: gphes(VECTOR_SIZE,ntens,mnode,pgaus)
    real(rp),    intent(in)    :: gprhs(VECTOR_SIZE,ndime,pgaus)
    real(rp),    intent(inout) :: gprhs_sgs(VECTOR_SIZE,ndime,pgaus)
    real(rp),    intent(in)    :: gprhc(VECTOR_SIZE,pgaus)
    real(rp),    intent(in)    :: gprh2(VECTOR_SIZE,pgaus)
    real(rp),    intent(in)    :: rmomu(VECTOR_SIZE,pnode,pgaus)
    real(rp),    intent(in)    :: rcont(VECTOR_SIZE,ndime,pnode,pgaus)
    real(rp),    intent(in)    :: elpre(VECTOR_SIZE,pnode,*)
    real(rp),    intent(in)    :: elbub(VECTOR_SIZE)
    real(rp),    intent(in)    :: densi(VECTOR_SIZE,pgaus,nbdfp_nsi)
    real(rp),    intent(in)    :: gplag(VECTOR_SIZE,ndime,ndime,pgaus)            ! Lagrange multiplier
    real(rp),    intent(in)    :: gpmix(VECTOR_SIZE,pgaus)                        ! Mixing function
    ! Element matrices
    real(rp),    intent(out)   :: elauu(VECTOR_SIZE,pnode*ndime,pnode*ndime)
    real(rp),    intent(out)   :: elaup(VECTOR_SIZE,pnode*ndime,pnode)
    real(rp),    intent(out)   :: elapp(VECTOR_SIZE,pnode,pnode)
    real(rp),    intent(out)   :: elapu(VECTOR_SIZE,pnode,pnode*ndime)
    real(rp),    intent(out)   :: elrbu(VECTOR_SIZE,ndime,pnode)
    real(rp),    intent(out)   :: elrbp(VECTOR_SIZE,pnode)
    real(rp),    intent(in)    :: rmom2(VECTOR_SIZE,ndime,ndime,pnode,pgaus)
    real(rp),    intent(in)    :: gpst1(VECTOR_SIZE,pgaus)
    real(rp),    intent(in)    :: gpgve(VECTOR_SIZE,ndime,ndime,pgaus)
    real(rp),    intent(in)    :: elvel(VECTOR_SIZE,ndime,pnode,*)
    real(rp),    intent(out)   :: ellum(VECTOR_SIZE,pnode)
    real(rp),    intent(in)    :: dtinv_loc(VECTOR_SIZE)
    integer(ip), intent(in)    :: pbubl(VECTOR_SIZE)
    real(rp),    intent(in)    :: gpsha_bub(VECTOR_SIZE,pgaus)                    ! Ne
    real(rp),    intent(in)    :: gpcar_bub(VECTOR_SIZE,ndime,pgaus)              ! dNe/dxi
    real(rp),    intent(in)    :: gppre(VECTOR_SIZE,pgaus)
    ! Enrichement Element matrices
    real(rp),    intent(out)   :: elauq(VECTOR_SIZE,pnode*ndime,1)
    real(rp),    intent(out)   :: elapq(VECTOR_SIZE,pnode,1)
    real(rp),    intent(out)   :: elaqu(VECTOR_SIZE,1,pnode*ndime)
    real(rp),    intent(out)   :: elaqp(VECTOR_SIZE,1,pnode)
    real(rp),    intent(out)   :: elaqq(VECTOR_SIZE,1,1)
    real(rp),    intent(out)   :: elrbq(VECTOR_SIZE,1)
    ! Indices and local variables
    integer(ip)                :: inode,jnode,idofv
    integer(ip)                :: idof1,jdime,itime
#ifdef OPENACC
    integer(ip)                :: ivect
#endif
    integer(ip)                :: igaus,idime,jdofv,kdime
    real(rp)                   :: xvis2,xvisc,one_rp
    ! Vectorized local arrays
    real(rp)                   :: p1ve2(VECTOR_SIZE,ndime,ndime,pnode,pgaus)
    real(rp)                   :: p1vec(VECTOR_SIZE,pnode,pgaus)
    real(rp)                   :: p2vec(VECTOR_SIZE,ndime,pnode,pgaus)
    real(rp)                   :: p2sca(VECTOR_SIZE,pnode,pgaus)
    real(rp)                   :: wgrgr(VECTOR_SIZE,pnode,pnode,pgaus)
    real(rp)                   :: wgrvi(VECTOR_SIZE,pnode,pgaus)
    real(rp)                   :: fact0(VECTOR_SIZE),fact1(VECTOR_SIZE)
    real(rp)                   :: fact2(VECTOR_SIZE),fact3(VECTOR_SIZE)
    real(rp)                   :: fact4(VECTOR_SIZE),fact5(VECTOR_SIZE)
    real(rp)                   :: fact6(VECTOR_SIZE),fact7(VECTOR_SIZE)
    real(rp)                   :: fact8(VECTOR_SIZE),factz(VECTOR_SIZE)
    real(rp)                   :: factx(VECTOR_SIZE),facty(VECTOR_SIZE)
    real(rp)                   :: facx1(VECTOR_SIZE),facy1(VECTOR_SIZE)
    real(rp)                   :: ugraN(VECTOR_SIZE)
    real(rp)                   :: gramugraN(VECTOR_SIZE)
    real(rp)                   :: penal(VECTOR_SIZE)
    real(rp)                   :: gprhh(VECTOR_SIZE,ndime,pgaus)
    real(rp)                   :: taupr(VECTOR_SIZE,pgaus)
    real(rp)                   :: gpveo(VECTOR_SIZE,ndime)
    real(rp)                   :: gpcar1ji(VECTOR_SIZE)
    real(rp)                   :: gpcar2ji(VECTOR_SIZE)
    real(rp)                   :: gpcar3ji(VECTOR_SIZE)
    real(rp)                   :: p2sca_bub(VECTOR_SIZE,pgaus)
    real(rp)                   :: p2vec_bub(VECTOR_SIZE,ndime,pgaus)
    real(rp)                   :: eps_test(VECTOR_SIZE,ndime,ndime)
    real(rp)                   :: delta(ndime,ndime)
    real(rp)                   :: xib
    real(rp)                   :: dtinv_mod(VECTOR_SIZE)
    real(rp)                   :: grpre(VECTOR_SIZE,ndime)
    real(rp)                   :: dudt(VECTOR_SIZE, ndime)
    
    real(rp), dimension(VECTOR_SIZE,4*((pnode+3)/4),pnode) :: elauu11,elauu21,elauu31
    real(rp), dimension(VECTOR_SIZE,4*((pnode+3)/4),pnode) :: elauu12,elauu22,elauu32
    real(rp), dimension(VECTOR_SIZE,4*((pnode+3)/4),pnode) :: elauu13,elauu23,elauu33
    real(rp), dimension(VECTOR_SIZE,4*((pnode+3)/4))       :: factvec1,factvec2,factvec3,factvec4

    !yor! the comments marked with 'yor' explain the recent changes for optimization purpose, within EoCoE project. To be removed later.
    !yor! optimization test : for memory padding purpose : 4*((pnode+3)/4) will give the closest number >= pnode which is a multiple of 4.
    !yor! We have to think of a clean way to write it. The idea could be expanded in the whole code, with the padding size as parameter.

#ifdef OPENACC
#define DEF_VECT ivect
#else
#define DEF_VECT 1:VECTOR_SIZE
#endif

    !----------------------------------------------------------------------
    !
    ! Fractional step: time derivatives are computed using mass matrix
    ! later on, when updating the velocity using AB or RK schemes
    !
    !----------------------------------------------------------------------

    if( NSI_FRACTIONAL_STEP ) then
       dtinv_mod = 0.0_rp
    else
       dtinv_mod = dtinv_loc
    end if
   
    !----------------------------------------------------------------------
    !
    ! Local variables
    !
    !----------------------------------------------------------------------
    !
    ! Viscous term factor
    !
    xvis2 = 0.0_rp
    if( fvins_nsi > 0.9_rp ) then
       xvisc = 1.0_rp
       if( fvins_nsi > 1.9_rp ) xvis2 = 2.0_rp/3.0_rp
    else
       xvisc = 0.0_rp
    end if
    !
    ! Do not stabilize convection, reaction, coriolis
    !
    if( kfl_nota1_nsi == 1 ) then
       one_rp = 0.0_rp
    else
       one_rp = 1.0_rp
    end if

#ifdef OPENACC
    do ivect = 1,VECTOR_SIZE
#endif

       !----------------------------------------------------------------------
       !
       ! Initialization
       !
       !----------------------------------------------------------------------

       elauu(DEF_VECT,:,:)   = 0.0_rp
       elaup(DEF_VECT,:,:)   = 0.0_rp
       elapu(DEF_VECT,:,:)   = 0.0_rp
       elapp(DEF_VECT,:,:)   = 0.0_rp
       elrbu(DEF_VECT,:,:)   = 0.0_rp
       elrbp(DEF_VECT,:)     = 0.0_rp

       elauu11(DEF_VECT,:,:) = 0.0_rp
       elauu21(DEF_VECT,:,:) = 0.0_rp
       elauu12(DEF_VECT,:,:) = 0.0_rp
       elauu22(DEF_VECT,:,:) = 0.0_rp
       elauu31(DEF_VECT,:,:) = 0.0_rp
       elauu32(DEF_VECT,:,:) = 0.0_rp
       elauu13(DEF_VECT,:,:) = 0.0_rp
       elauu23(DEF_VECT,:,:) = 0.0_rp
       elauu33(DEF_VECT,:,:) = 0.0_rp
       !
       ! Incompressible vs low Mach
       !
       if( kfl_regim_nsi == 3 ) then
          taupr(DEF_VECT,1:pgaus) = 1.0_rp
       else
          taupr(DEF_VECT,1:pgaus) = 1.0_rp / max(gpst1(DEF_VECT,1:pgaus),zeror)
       end if

       !----------------------------------------------------------------------
       !
       ! Test functions
       !
       !----------------------------------------------------------------------
       !
       ! P1VEC =  v * [ (tau1'/tau1) - tau1' * sig ] + tau1' * rho * (uc.grad)v
       ! P2SCA =  ( tau2^{-1} * tau2' ) * q
       ! P2VEC =  tau1' * grad(q)  ( and * rho if in Low-Mach, but not yet)
       ! WGRVI =  grad(Ni) . grad(mu) + mu * Delta(Ni)
       ! WGRGR =  grad(Ni) . grad(Nj)
       !
       do igaus = 1,pgaus

          fact1(DEF_VECT) = gpsp1(DEF_VECT,igaus) * one_rp * gpden(DEF_VECT,igaus) * gpvol(DEF_VECT,igaus)
          fact2(DEF_VECT) = (gptt1(DEF_VECT,igaus)-gpsp1(DEF_VECT,igaus) * one_rp * gppor(DEF_VECT,igaus)) * gpvol(DEF_VECT,igaus)
          fact3(DEF_VECT) = gpsp1(DEF_VECT,igaus) * gpvol(DEF_VECT,igaus)
          fact4(DEF_VECT) = gptt2(DEF_VECT,igaus) * gpvol(DEF_VECT,igaus)

          do inode = 1,pnode

             ugraN(DEF_VECT)     = 0.0_rp
             gramugraN(DEF_VECT) = 0.0_rp

             do idime = 1,ndime
                ugraN(DEF_VECT)                   = ugraN(DEF_VECT)     + gpadv(DEF_VECT,idime,igaus) * gpcar(DEF_VECT,idime,inode,igaus)
                gramugraN(DEF_VECT)               = gramugraN(DEF_VECT) + gpgvi(DEF_VECT,idime,igaus) * gpcar(DEF_VECT,idime,inode,igaus)
                p2vec(DEF_VECT,idime,inode,igaus) = fact3(DEF_VECT) * gpcar(DEF_VECT,idime,inode,igaus)
             end do

             p1vec(DEF_VECT,inode,igaus)          = fact2(DEF_VECT) * gpsha(DEF_VECT,inode,igaus) + fact1(DEF_VECT) * ugraN(DEF_VECT)
             p2sca(DEF_VECT,inode,igaus)          = fact4(DEF_VECT) * gpsha(DEF_VECT,inode,igaus)
             wgrvi(DEF_VECT,inode,igaus)          = gramugraN(DEF_VECT) + gpvis(DEF_VECT,igaus) * gplap(DEF_VECT,inode,igaus)

             do jnode = 1,pnode
                wgrgr(DEF_VECT,inode,jnode,igaus) = 0.0_rp
                do idime = 1,ndime
                   wgrgr(DEF_VECT,inode,jnode,igaus) = wgrgr(DEF_VECT,inode,jnode,igaus) + gpcar(DEF_VECT,idime,inode,igaus) * gpcar(DEF_VECT,idime,jnode,igaus)
                end do
             end do

          end do
       end do
       if( NSI_FRACTIONAL_STEP ) p2vec = 0.0_rp 
  
       !
       ! P1VE2: Off-diagonal test function
       !
       if( kfl_p1ve2_nsi /= 0 ) then
          !
          ! tau1'*2*rho*(w x v)
          ! x-equation: v = (v,0,0) => w x v = (    0 , wz*v , -wy*v)
          ! y-equation: v = (0,v,0) => w x v = (-wz*v ,    0 ,  wx*v)
          ! z-equation: v = (0,0,v) => w x v = ( wy*v ,-wx*v ,     0)
          !
          p1ve2 = 0.0_rp
          if( ndime == 2 ) then
             do igaus = 1,pgaus
                fact3(DEF_VECT) = 2.0_rp * gpden(DEF_VECT,igaus) * fvela_nsi(3)
                factz(DEF_VECT) = gpsp1(DEF_VECT,igaus)* fact3(DEF_VECT) * gpvol(DEF_VECT,igaus)
                do inode = 1,pnode
                   fact1(DEF_VECT)                 =   factz(DEF_VECT) * gpsha(DEF_VECT,inode,igaus)
                   p1ve2(DEF_VECT,1,2,inode,igaus) =   fact1(DEF_VECT)
                   p1ve2(DEF_VECT,2,1,inode,igaus) = - fact1(DEF_VECT)
                end do
             end do
          else
             do igaus = 1,pgaus
                if( NSI_FRACTIONAL_STEP ) then
                   fact8(DEF_VECT) = 2.0_rp * gpvol(DEF_VECT,igaus)/dtinv_loc ! tau = dt/rho
                else
                   fact8(DEF_VECT) = 2.0_rp * gpden(DEF_VECT,igaus) * gpsp1(DEF_VECT,igaus)* gpvol(DEF_VECT,igaus)
                end if
                factx(DEF_VECT) = fact8(DEF_VECT)  * fvela_nsi(1)
                facty(DEF_VECT) = fact8(DEF_VECT)  * fvela_nsi(2)
                factz(DEF_VECT) = fact8(DEF_VECT)  * fvela_nsi(3)
                do inode = 1,pnode
                   p1ve2(DEF_VECT,1,2,inode,igaus) =  factz(DEF_VECT) * gpsha(DEF_VECT,inode,igaus)  ! x-equation
                   p1ve2(DEF_VECT,1,3,inode,igaus) = -facty(DEF_VECT) * gpsha(DEF_VECT,inode,igaus)  ! x-equation
                   p1ve2(DEF_VECT,2,1,inode,igaus) = -factz(DEF_VECT) * gpsha(DEF_VECT,inode,igaus)  ! y-equation
                   p1ve2(DEF_VECT,2,3,inode,igaus) =  factx(DEF_VECT) * gpsha(DEF_VECT,inode,igaus)  ! y-equation
                   p1ve2(DEF_VECT,3,1,inode,igaus) =  facty(DEF_VECT) * gpsha(DEF_VECT,inode,igaus)  ! z-equation
                   p1ve2(DEF_VECT,3,2,inode,igaus) = -factx(DEF_VECT) * gpsha(DEF_VECT,inode,igaus)  ! z-equation
                end do
             end do
          end if
       end if

       !----------------------------------------------------------------------
       !
       ! Auu
       !
       !----------------------------------------------------------------------
       !
       !  Laplacian  form: A=0, B=0, eps(u) = 1/2 grad(u)
       !  Divergence form: A=1, B=0, eps(u) = 1/2 ( grad(u) + grad(u)^t )
       !  Complete   form: A=1, B=1, eps(u) = 1/2 ( grad(u) + grad(u)^t ) - 1/3 div(u) I
       !
       !    ( div(u) , tau2' div(v) )             (1)   <=   tau2' * duj/dxj * dvi/dxi
       !  + ( rmomu(DEF_VECT,u) , p1vec(DEF_VECT,v) )               (2)
       !  + ( mu dui/dxj , dvi/dxj )              (3)   <=   ( 2 mu eps(u) , grad(v) ) = mu * ( dui/dxj + duj/dxi - 2/3 div(u) delta_ij ) dvi/dxj
       !  + A * ( mu duj/dxi , dvi/dxj )          (4)        divergence
       !  - 2/3 * B * ( mu (div u) , dvi/dxi )    (5)        complete
       !  + ( mu d^2ui/dxj^2 , vi )               (6)   <= - ( div[-2 mu eps(u)] , v ) = - d/dxj ( -2*mu* ( dui/dxj + duj/dxi - 2/3*mu div(u) delta_ij ) * vi
       !  + ( dmu/dxj dui/dxj , vi )              (7)        laplacian
       !  + A * ( mu  d^2uj/dxidxj , vi )         (8)        divergence
       !  + A * ( dmu/dxj duj/dxi , vi )          (9)        divergence
       !  - 2/3 * B * ( dmu/dxi (div u) , vi )   (10)        complete
       !  - 2/3 * B * ( mu d(div u)/dxi , vi )   (11)        complete
       !
       !  1. The terms ( div(u) , tau2' div(v) ) and - 2/3 * B * ( mu div u, dvi/dxj ) are assembled together as
       !     ( ( tau2' - 2/3 B mu ) div u, dvi/dxi ). Therefore, term (5) is assembled always although it is zero
       !  2. Terms (3), (6) and (7) are assembled just below as they are the same for all momentum equations
       !  3. Terms (4), (8), (9), (10) and (11) are assmbled later on
       !
       !  (3)        x:     mu * ( dux/dx * dv/dx  +  dux/dy * dv/dy  +  dux/dz * dv/dz )
       !             y:     mu * ( duy/dx * dv/dx  +  duy/dy * dv/dy  +  duy/dz * dv/dz )
       !             z:     mu * ( duz/dx * dv/dx  +  duz/dy * dv/dy  +  duz/dz * dv/dz )
       !
       !  (4)        x: A * mu * ( dux/dx * dv/dx  +  duy/dx * dv/dy  +  duz/dx * dv/dz )
       !             y: A * mu * ( dux/dy * dv/dx  +  duy/dy * dv/dy  +  duz/dy * dv/dz )
       !             z: A * mu * ( dux/dz * dv/dx  +  duy/dz * dv/dy  +  duz/dz * dv/dz )
       !
       !  (6)        x:  mu * ( d^2ux/dx^2  +  d^2ux/dy^2  +  d^2ux/dz^2 ) v
       !             y:  mu * ( d^2uy/dx^2  +  d^2uy/dy^2  +  d^2uy/dz^2 ) v
       !             z:  mu * ( d^2uz/dx^2  +  d^2uz/dy^2  +  d^2uz/dz^2 ) v
       !
       !  (7)        x: ( dmu/dx * dux/dx  +  dmu/dy * dux/dy  +  dmu/dz * dux/dz ) v
       !             y: ( dmu/dx * duy/dx  +  dmu/dy * duy/dy  +  dmu/dz * duy/dz ) v
       !             y: ( dmu/dx * duz/dx  +  dmu/dy * duz/dy  +  dmu/dz * duz/dz ) v
       !
       !  (8)        x: A * mu * ( d^2ux/dxdx  +  d^2uy/dxdy  +  d^2uz/dxdz ) v
       !             y: A * mu * ( d^2ux/dxdy  +  d^2uy/dydy  +  d^2uz/dydz ) v
       !             z: A * mu * ( d^2ux/dxdz  +  d^2uy/dzdy  +  d^2uz/dzdz ) v
       !
       !  (9)        x: A * ( dmu/dx * dux/dx  +  dmu/dy * duy/dx  +  dmu/dz * duz/dx ) v
       !             y: A * ( dmu/dx * dux/dy  +  dmu/dy * duy/dy  +  dmu/dz * duz/dy ) v
       !             z: A * ( dmu/dx * dux/dz  +  dmu/dy * duy/dz  +  dmu/dz * duz/dz ) v
       !
       !  (10)       x: -2/3 B * dmu/dx ( dux/dx  +  duy/dy  +  duz/dz ) * v
       !             y: -2/3 B * dmu/dy ( dux/dx  +  duy/dy  +  duz/dz ) * v
       !             z: -2/3 B * dmu/dz ( dux/dx  +  duy/dy  +  duz/dz ) * v
       !
       !  (11)       x: -2/3 B * mu * ( d^2ux/dxdx + d^2uy/dydx + d^2uz/dzdx ) v
       !             y: -2/3 B * mu * ( d^2ux/dxdy + d^2uy/dydy + d^2uz/dzdy ) v
       !             z: -2/3 B * mu * ( d^2ux/dxdz + d^2uy/dydz + d^2uz/dzdz ) v
       !
       !  (1) + (5)  x: ( tau2' - 2/3 B mu ) ( dux/dx  +  duy/dy  +  duz/dz ) * dv/dx
       !             y: ( tau2' - 2/3 B mu ) ( dux/dx  +  duy/dy  +  duz/dz ) * dv/dy
       !             z: ( tau2' - 2/3 B mu ) ( dux/dx  +  duy/dy  +  duz/dz ) * dv/dz
       !
       if( ndime == 2 ) then

          do igaus = 1,pgaus

             fact0(DEF_VECT) = ( gpsp2(DEF_VECT,igaus) - xvis2 * gpvis(DEF_VECT,igaus) ) * gpvol(DEF_VECT,igaus)                                      ! (1) + (5)
             fact6(DEF_VECT) = gpvis(DEF_VECT,igaus) * gpvol(DEF_VECT,igaus)

             do inode = 1,pnode
                fact1(DEF_VECT) = fact0(DEF_VECT) * gpcar(DEF_VECT,1,inode,igaus)
                fact2(DEF_VECT) = fact0(DEF_VECT) * gpcar(DEF_VECT,2,inode,igaus)
                fact5(DEF_VECT) = gpsha(DEF_VECT,inode,igaus) * gpvol(DEF_VECT,igaus)

                do jnode = 1,pnode
                   fact4(DEF_VECT)               =   p1vec(DEF_VECT,inode,igaus) * rmomu(DEF_VECT,jnode,igaus) &                                          ! (2):       ( rmomu(DEF_VECT,u) , p1vec(DEF_VECT,v) )
                        &                          + fact5(DEF_VECT) * wgrvi(DEF_VECT,jnode,igaus)              &                                         ! (3) + (6): ( mu d^2ui/dxj^2 + dmu/dxj dui/dxj, vi )
                        &                          + fact6(DEF_VECT) * wgrgr(DEF_VECT,inode,jnode,igaus)                                                  ! (7):       ( mu dui/dxj , dvi/dxj )

                   elauu11(DEF_VECT,inode,jnode) = elauu11(DEF_VECT,inode,jnode) + fact1(DEF_VECT) * rcont(DEF_VECT,1,jnode,igaus) + fact4(DEF_VECT)      ! Auu_xx
                   elauu21(DEF_VECT,inode,jnode) = elauu21(DEF_VECT,inode,jnode) + fact2(DEF_VECT) * rcont(DEF_VECT,1,jnode,igaus)                        ! Auu_yx
                   elauu12(DEF_VECT,inode,jnode) = elauu12(DEF_VECT,inode,jnode) + fact1(DEF_VECT) * rcont(DEF_VECT,2,jnode,igaus)                        ! Auu_xy
                   elauu22(DEF_VECT,inode,jnode) = elauu22(DEF_VECT,inode,jnode) + fact2(DEF_VECT) * rcont(DEF_VECT,2,jnode,igaus) + fact4(DEF_VECT)      ! Auu_yy
                end do

             end do

          end do

       else

          do igaus = 1,pgaus

             fact0(DEF_VECT) = ( gpsp2(DEF_VECT,igaus) - xvis2 * gpvis(DEF_VECT,igaus) ) * gpvol(DEF_VECT,igaus)                                          ! (1) + (5)
             fact6(DEF_VECT) = gpvis(DEF_VECT,igaus) * gpvol(DEF_VECT,igaus)

             factvec1(DEF_VECT,1:pnode) = gpcar(DEF_VECT,1,1:pnode,igaus)
             factvec2(DEF_VECT,1:pnode) = gpcar(DEF_VECT,2,1:pnode,igaus)
             factvec3(DEF_VECT,1:pnode) = gpcar(DEF_VECT,3,1:pnode,igaus)

             do jnode = 1,pnode
                fact1(DEF_VECT) = fact0(DEF_VECT)*rcont(DEF_VECT,1,jnode,igaus)
                fact2(DEF_VECT) = fact0(DEF_VECT)*rcont(DEF_VECT,2,jnode,igaus)
                fact3(DEF_VECT) = fact0(DEF_VECT)*rcont(DEF_VECT,3,jnode,igaus)
                fact5(DEF_VECT) = wgrvi(DEF_VECT,jnode,igaus)*gpvol(DEF_VECT,igaus)

                do inode = 1,pnode
                   factvec4(DEF_VECT,inode) = p1vec(DEF_VECT,inode,igaus) * rmomu(DEF_VECT,jnode,igaus) + fact5(DEF_VECT) * gpsha(DEF_VECT,inode,igaus) &
                        + fact6(DEF_VECT) * wgrgr(DEF_VECT,inode,jnode,igaus)
                end do

                do inode = 1,pnode
                   elauu11(DEF_VECT,inode,jnode) = elauu11(DEF_VECT,inode,jnode) + fact1(DEF_VECT) * factvec1(DEF_VECT,inode) + factvec4(DEF_VECT,inode)  ! Auu_xx
                   elauu12(DEF_VECT,inode,jnode) = elauu12(DEF_VECT,inode,jnode) + fact2(DEF_VECT) * factvec1(DEF_VECT,inode)                             ! Auu_xy
                   elauu13(DEF_VECT,inode,jnode) = elauu13(DEF_VECT,inode,jnode) + fact3(DEF_VECT) * factvec1(DEF_VECT,inode)                             ! Auu_xz

                   elauu21(DEF_VECT,inode,jnode) = elauu21(DEF_VECT,inode,jnode) + fact1(DEF_VECT) * factvec2(DEF_VECT,inode)                             ! Auu_yx
                   elauu22(DEF_VECT,inode,jnode) = elauu22(DEF_VECT,inode,jnode) + fact2(DEF_VECT) * factvec2(DEF_VECT,inode) + factvec4(DEF_VECT,inode)  ! Auu_yy
                   elauu23(DEF_VECT,inode,jnode) = elauu23(DEF_VECT,inode,jnode) + fact3(DEF_VECT) * factvec2(DEF_VECT,inode)                             ! Auu_yz

                   elauu31(DEF_VECT,inode,jnode) = elauu31(DEF_VECT,inode,jnode) + fact1(DEF_VECT) * factvec3(DEF_VECT,inode)                             ! Auu_zx
                   elauu32(DEF_VECT,inode,jnode) = elauu32(DEF_VECT,inode,jnode) + fact2(DEF_VECT) * factvec3(DEF_VECT,inode)                             ! Auu_zy
                   elauu33(DEF_VECT,inode,jnode) = elauu33(DEF_VECT,inode,jnode) + fact3(DEF_VECT) * factvec3(DEF_VECT,inode) + factvec4(DEF_VECT,inode)  ! Auu_zz
                end do

             end do

          end do

       end if

       if( fvins_nsi > 0.9_rp ) then

          if( ndime == 2 ) then

             do igaus = 1,pgaus

                fact4(DEF_VECT) = gpgvi(DEF_VECT,1,igaus) * gpvol(DEF_VECT,igaus)
                fact5(DEF_VECT) = gpgvi(DEF_VECT,2,igaus) * gpvol(DEF_VECT,igaus)
                fact7(DEF_VECT) = gpvis(DEF_VECT,igaus)   * gpvol(DEF_VECT,igaus)

                do inode = 1,pnode
                   fact0(DEF_VECT) = fact7(DEF_VECT) * gpsha(DEF_VECT,inode,igaus) *(1.0_rp-xvis2)
                   fact1(DEF_VECT) = fact4(DEF_VECT) * gpsha(DEF_VECT,inode,igaus)
                   fact2(DEF_VECT) = fact5(DEF_VECT) * gpsha(DEF_VECT,inode,igaus)
                   facx1(DEF_VECT) = fact1(DEF_VECT) * xvis2
                   facy1(DEF_VECT) = fact2(DEF_VECT) * xvis2
                   do jnode = 1,pnode
                      !
                      !                                         (4)  ( mu duj/dxi , dvi/dxj )
                      elauu11(DEF_VECT,inode,jnode) = elauu11(DEF_VECT,inode,jnode) + fact7(DEF_VECT) * gpcar(DEF_VECT,1,inode,igaus) * gpcar(DEF_VECT,1,jnode,igaus)                   ! Auu_xx
                      elauu12(DEF_VECT,inode,jnode) = elauu12(DEF_VECT,inode,jnode) + fact7(DEF_VECT) * gpcar(DEF_VECT,2,inode,igaus) * gpcar(DEF_VECT,1,jnode,igaus)                   ! Auu_xy
                      elauu21(DEF_VECT,inode,jnode) = elauu21(DEF_VECT,inode,jnode) + fact7(DEF_VECT) * gpcar(DEF_VECT,1,inode,igaus) * gpcar(DEF_VECT,2,jnode,igaus)                   ! Auu_yx
                      elauu22(DEF_VECT,inode,jnode) = elauu22(DEF_VECT,inode,jnode) + fact7(DEF_VECT) * gpcar(DEF_VECT,2,inode,igaus) * gpcar(DEF_VECT,2,jnode,igaus)                   ! Auu_yy
                      !
                      !                                         (10) - 2/3 * B * ( dmu/dxi (div u) , vi )
                      elauu11(DEF_VECT,inode,jnode) = elauu11(DEF_VECT,inode,jnode) - facx1(DEF_VECT) * gpcar(DEF_VECT,1,jnode,igaus)                                                   ! Auu_xx
                      elauu12(DEF_VECT,inode,jnode) = elauu12(DEF_VECT,inode,jnode) - facx1(DEF_VECT) * gpcar(DEF_VECT,2,jnode,igaus)                                                   ! Auu_xy
                      elauu21(DEF_VECT,inode,jnode) = elauu21(DEF_VECT,inode,jnode) - facy1(DEF_VECT) * gpcar(DEF_VECT,1,jnode,igaus)                                                   ! Auu_yx
                      elauu22(DEF_VECT,inode,jnode) = elauu22(DEF_VECT,inode,jnode) - facy1(DEF_VECT) * gpcar(DEF_VECT,2,jnode,igaus)                                                   ! Auu_yy
                      !
                      !                                         (9)  ( dmu/dxj duj/dxi , vi ) + (8)  ( mu d^2uj/dxidxj , vi )
                      !                                               + (11)   - 2/3 * B * ( mu d(div u)/dxi , vi )
                      elauu11(DEF_VECT,inode,jnode) = elauu11(DEF_VECT,inode,jnode) + fact1(DEF_VECT) * gpcar(DEF_VECT,1,jnode,igaus) + fact0(DEF_VECT) * gphes(DEF_VECT,1,jnode,igaus) ! Auu_xx
                      elauu12(DEF_VECT,inode,jnode) = elauu12(DEF_VECT,inode,jnode) + fact2(DEF_VECT) * gpcar(DEF_VECT,1,jnode,igaus) + fact0(DEF_VECT) * gphes(DEF_VECT,3,jnode,igaus) ! Auu_xy
                      elauu21(DEF_VECT,inode,jnode) = elauu21(DEF_VECT,inode,jnode) + fact1(DEF_VECT) * gpcar(DEF_VECT,2,jnode,igaus) + fact0(DEF_VECT) * gphes(DEF_VECT,3,jnode,igaus) ! Auu_yx
                      elauu22(DEF_VECT,inode,jnode) = elauu22(DEF_VECT,inode,jnode) + fact2(DEF_VECT) * gpcar(DEF_VECT,2,jnode,igaus) + fact0(DEF_VECT) * gphes(DEF_VECT,2,jnode,igaus) ! Auu_yy
                      !
                      !

                   end do

                end do
             end do

          else

             do igaus = 1,pgaus

                fact4(DEF_VECT) = gpgvi(DEF_VECT,1,igaus) * gpvol(DEF_VECT,igaus)
                fact5(DEF_VECT) = gpgvi(DEF_VECT,2,igaus) * gpvol(DEF_VECT,igaus)
                fact6(DEF_VECT) = gpgvi(DEF_VECT,3,igaus) * gpvol(DEF_VECT,igaus)
                fact7(DEF_VECT) = gpvis(DEF_VECT,igaus)   * gpvol(DEF_VECT,igaus)
                !          factx(DEF_VECT) = xvis2 * fact4(DEF_VECT)
                !          facty(DEF_VECT) = xvis2 * fact5(DEF_VECT)
                !          factz(DEF_VECT) = xvis2 * fact6(DEF_VECT)

                !yor! [Originally line 451] reverted loops order, refactored the factors, with an effort to vectorize each instruction
                !
                !                                (4)  ( mu duj/dxi , dvi/dxj )
                factvec1(DEF_VECT,1:pnode) = gpcar(DEF_VECT,1,1:pnode,igaus)
                factvec2(DEF_VECT,1:pnode) = gpcar(DEF_VECT,2,1:pnode,igaus)
                factvec3(DEF_VECT,1:pnode) = gpcar(DEF_VECT,3,1:pnode,igaus)

                do jnode = 1,pnode
                   gpcar1ji(DEF_VECT) = gpcar(DEF_VECT,1,jnode,igaus)
                   gpcar2ji(DEF_VECT) = gpcar(DEF_VECT,2,jnode,igaus)
                   gpcar3ji(DEF_VECT) = gpcar(DEF_VECT,3,jnode,igaus)
                   !
                   do inode = 1,pnode
                      elauu11(DEF_VECT,inode,jnode) = elauu11(DEF_VECT,inode,jnode) + fact7(DEF_VECT) * factvec1(DEF_VECT,inode) * gpcar1ji(DEF_VECT)         ! Auu_xx
                      elauu21(DEF_VECT,inode,jnode) = elauu21(DEF_VECT,inode,jnode) + fact7(DEF_VECT) * factvec1(DEF_VECT,inode) * gpcar2ji(DEF_VECT)         ! Auu_yx
                      elauu31(DEF_VECT,inode,jnode) = elauu31(DEF_VECT,inode,jnode) + fact7(DEF_VECT) * factvec1(DEF_VECT,inode) * gpcar3ji(DEF_VECT)         ! Auu_zx

                      elauu12(DEF_VECT,inode,jnode) = elauu12(DEF_VECT,inode,jnode) + fact7(DEF_VECT) * factvec2(DEF_VECT,inode) * gpcar1ji(DEF_VECT)         ! Auu_xy
                      elauu22(DEF_VECT,inode,jnode) = elauu22(DEF_VECT,inode,jnode) + fact7(DEF_VECT) * factvec2(DEF_VECT,inode) * gpcar2ji(DEF_VECT)         ! Auu_yy
                      elauu32(DEF_VECT,inode,jnode) = elauu32(DEF_VECT,inode,jnode) + fact7(DEF_VECT) * factvec2(DEF_VECT,inode) * gpcar3ji(DEF_VECT)         ! Auu_zy

                      elauu13(DEF_VECT,inode,jnode) = elauu13(DEF_VECT,inode,jnode) + fact7(DEF_VECT) * factvec3(DEF_VECT,inode) * gpcar1ji(DEF_VECT)         ! Auu_xz
                      elauu23(DEF_VECT,inode,jnode) = elauu23(DEF_VECT,inode,jnode) + fact7(DEF_VECT) * factvec3(DEF_VECT,inode) * gpcar2ji(DEF_VECT)         ! Auu_yz
                      elauu33(DEF_VECT,inode,jnode) = elauu33(DEF_VECT,inode,jnode) + fact7(DEF_VECT) * factvec3(DEF_VECT,inode) * gpcar3ji(DEF_VECT)         ! Auu_zz
                   end do
                end do
                !
                !                     (10) - 2/3 * B * ( dmu/dxi (div u) , vi )

                do inode = 1,pnode
                   factvec1(DEF_VECT,inode) = fact4(DEF_VECT) * gpsha(DEF_VECT,inode,igaus) * xvis2 !yor! (ex facx1) this copy is for pure esthaetical purpose
                   factvec2(DEF_VECT,inode) = fact5(DEF_VECT) * gpsha(DEF_VECT,inode,igaus) * xvis2 !yor! (ex facy1)
                   factvec3(DEF_VECT,inode) = fact6(DEF_VECT) * gpsha(DEF_VECT,inode,igaus) * xvis2 !yor! (ex facz1)
                end do

                do jnode = 1,pnode
                   gpcar1ji(DEF_VECT) = gpcar(DEF_VECT,1,jnode,igaus)
                   gpcar2ji(DEF_VECT) = gpcar(DEF_VECT,2,jnode,igaus)
                   gpcar3ji(DEF_VECT) = gpcar(DEF_VECT,3,jnode,igaus)

                   do inode = 1,pnode
                      elauu11(DEF_VECT,inode,jnode) = elauu11(DEF_VECT,inode,jnode) - factvec1(DEF_VECT,inode) * gpcar1ji(DEF_VECT)                                ! Auu_xx
                      elauu12(DEF_VECT,inode,jnode) = elauu12(DEF_VECT,inode,jnode) - factvec1(DEF_VECT,inode) * gpcar2ji(DEF_VECT)                                ! Auu_xy
                      elauu13(DEF_VECT,inode,jnode) = elauu13(DEF_VECT,inode,jnode) - factvec1(DEF_VECT,inode) * gpcar3ji(DEF_VECT)                                ! Auu_xz

                      elauu21(DEF_VECT,inode,jnode) = elauu21(DEF_VECT,inode,jnode) - factvec2(DEF_VECT,inode) * gpcar1ji(DEF_VECT)                                ! Auu_yx
                      elauu22(DEF_VECT,inode,jnode) = elauu22(DEF_VECT,inode,jnode) - factvec2(DEF_VECT,inode) * gpcar2ji(DEF_VECT)                                ! Auu_yy
                      elauu23(DEF_VECT,inode,jnode) = elauu23(DEF_VECT,inode,jnode) - factvec2(DEF_VECT,inode) * gpcar3ji(DEF_VECT)                                ! Auu_yz

                      elauu31(DEF_VECT,inode,jnode) = elauu31(DEF_VECT,inode,jnode) - factvec3(DEF_VECT,inode) * gpcar1ji(DEF_VECT)                                ! Auu_zx
                      elauu32(DEF_VECT,inode,jnode) = elauu32(DEF_VECT,inode,jnode) - factvec3(DEF_VECT,inode) * gpcar2ji(DEF_VECT)                                ! Auu_zy
                      elauu33(DEF_VECT,inode,jnode) = elauu33(DEF_VECT,inode,jnode) - factvec3(DEF_VECT,inode) * gpcar3ji(DEF_VECT)                                ! Auu_zz
                   end do
                end do

                !
                !                                         (9)  ( dmu/dxj duj/dxi , vi ) + (8)  ( mu d^2uj/dxidxj , vi )
                !                                                + (11)   - 2/3 * B * ( mu d(div u)/dxi , vi )
                do inode = 1,pnode
                   factvec4(DEF_VECT,inode) = fact7(DEF_VECT) * gpsha(DEF_VECT,inode,igaus)*(1.0_rp-xvis2) !yor!(fact0(DEF_VECT))
                   factvec1(DEF_VECT,inode) = fact4(DEF_VECT) * gpsha(DEF_VECT,inode,igaus) !yor!(ex fact1(DEF_VECT))
                   factvec2(DEF_VECT,inode) = fact5(DEF_VECT) * gpsha(DEF_VECT,inode,igaus) !yor!(ex fact2(DEF_VECT))
                   factvec3(DEF_VECT,inode) = fact6(DEF_VECT) * gpsha(DEF_VECT,inode,igaus) !yor!(ex fact3(DEF_VECT))
                end do

                do jnode = 1,pnode
                   gpcar1ji(DEF_VECT) = gpcar(DEF_VECT,1,jnode,igaus)
                   gpcar2ji(DEF_VECT) = gpcar(DEF_VECT,2,jnode,igaus)
                   gpcar3ji(DEF_VECT) = gpcar(DEF_VECT,3,jnode,igaus)

                   do inode = 1,pnode
                      elauu11(DEF_VECT,inode,jnode) = elauu11(DEF_VECT,inode,jnode) + factvec1(DEF_VECT,inode) * gpcar1ji(DEF_VECT) + factvec4(DEF_VECT,inode) * gphes(DEF_VECT,1,jnode,igaus) ! Auu_xx
                      elauu21(DEF_VECT,inode,jnode) = elauu21(DEF_VECT,inode,jnode) + factvec1(DEF_VECT,inode) * gpcar2ji(DEF_VECT) + factvec4(DEF_VECT,inode) * gphes(DEF_VECT,4,jnode,igaus) ! Auu_yx
                      elauu31(DEF_VECT,inode,jnode) = elauu31(DEF_VECT,inode,jnode) + factvec1(DEF_VECT,inode) * gpcar3ji(DEF_VECT) + factvec4(DEF_VECT,inode) * gphes(DEF_VECT,5,jnode,igaus) ! Auu_zx

                      elauu12(DEF_VECT,inode,jnode) = elauu12(DEF_VECT,inode,jnode) + factvec2(DEF_VECT,inode) * gpcar1ji(DEF_VECT) + factvec4(DEF_VECT,inode) * gphes(DEF_VECT,4,jnode,igaus) ! Auu_xy
                      elauu22(DEF_VECT,inode,jnode) = elauu22(DEF_VECT,inode,jnode) + factvec2(DEF_VECT,inode) * gpcar2ji(DEF_VECT) + factvec4(DEF_VECT,inode) * gphes(DEF_VECT,2,jnode,igaus) ! Auu_yy
                      elauu32(DEF_VECT,inode,jnode) = elauu32(DEF_VECT,inode,jnode) + factvec2(DEF_VECT,inode) * gpcar3ji(DEF_VECT) + factvec4(DEF_VECT,inode) * gphes(DEF_VECT,6,jnode,igaus) ! Auu_zy

                      elauu13(DEF_VECT,inode,jnode) = elauu13(DEF_VECT,inode,jnode) + factvec3(DEF_VECT,inode) * gpcar1ji(DEF_VECT) + factvec4(DEF_VECT,inode) * gphes(DEF_VECT,5,jnode,igaus) ! Auu_xz
                      elauu23(DEF_VECT,inode,jnode) = elauu23(DEF_VECT,inode,jnode) + factvec3(DEF_VECT,inode) * gpcar2ji(DEF_VECT) + factvec4(DEF_VECT,inode) * gphes(DEF_VECT,6,jnode,igaus) ! Auu_yz
                      elauu33(DEF_VECT,inode,jnode) = elauu33(DEF_VECT,inode,jnode) + factvec3(DEF_VECT,inode) * gpcar3ji(DEF_VECT) + factvec4(DEF_VECT,inode) * gphes(DEF_VECT,3,jnode,igaus) ! Auu_zz
                   end do
                end do
             end do
          end if
       end if
       !
       ! Skew symmetric when fcons=0.5
       !
       if( fcons_nsi > 0.1_rp ) then
          if( ndime == 2 ) then
             do igaus = 1,pgaus
                fact0(DEF_VECT) = fcons_nsi * gpden(DEF_VECT,igaus) * gpvol(DEF_VECT,igaus)
                do inode = 1,pnode
                   do jnode = 1,pnode
                      fact1(DEF_VECT) = 0.0_rp
                      fact2(DEF_VECT) = 0.0_rp
                      do idime = 1,2
                         fact1(DEF_VECT) = fact1(DEF_VECT) + gpadv(DEF_VECT,idime,igaus) * gpcar(DEF_VECT,idime,inode,igaus)
                         fact2(DEF_VECT) = fact2(DEF_VECT) + gpadv(DEF_VECT,idime,igaus) * gpcar(DEF_VECT,idime,jnode,igaus)
                      end do
                      fact1(DEF_VECT) = fact1(DEF_VECT) * gpsha(DEF_VECT,jnode,igaus) * fact0(DEF_VECT)
                      fact2(DEF_VECT) = fact2(DEF_VECT) * gpsha(DEF_VECT,inode,igaus) * fact0(DEF_VECT)

                      elauu11(DEF_VECT,inode,jnode) = elauu11(DEF_VECT,inode,jnode) - fact1(DEF_VECT) - fact2(DEF_VECT)                      ! Auu_xx
                      elauu22(DEF_VECT,inode,jnode) = elauu22(DEF_VECT,inode,jnode) - fact1(DEF_VECT) - fact2(DEF_VECT)                      ! Auu_yy
                   end do
                end do
             end do
          else
             do igaus = 1,pgaus
                fact0(DEF_VECT) = fcons_nsi * gpden(DEF_VECT,igaus) * gpvol(DEF_VECT,igaus)
                do inode = 1,pnode
                   do jnode = 1,pnode
                      fact1(DEF_VECT) = 0.0_rp
                      fact2(DEF_VECT) = 0.0_rp
                      do idime = 1,3
                         fact1(DEF_VECT) = fact1(DEF_VECT) + gpadv(DEF_VECT,idime,igaus) * gpcar(DEF_VECT,idime,inode,igaus)
                         fact2(DEF_VECT) = fact2(DEF_VECT) + gpadv(DEF_VECT,idime,igaus) * gpcar(DEF_VECT,idime,jnode,igaus)
                      end do
                      fact1(DEF_VECT) = fact1(DEF_VECT) * gpsha(DEF_VECT,jnode,igaus) * fact0(DEF_VECT)
                      fact2(DEF_VECT) = fact2(DEF_VECT) * gpsha(DEF_VECT,inode,igaus) * fact0(DEF_VECT)

                      elauu11(DEF_VECT,inode,jnode) = elauu11(DEF_VECT,inode,jnode) - fact1(DEF_VECT) - fact2(DEF_VECT)                      ! Auu_xx
                      elauu22(DEF_VECT,inode,jnode) = elauu22(DEF_VECT,inode,jnode) - fact1(DEF_VECT) - fact2(DEF_VECT)                      ! Auu_yy
                      elauu33(DEF_VECT,inode,jnode) = elauu33(DEF_VECT,inode,jnode) - fact1(DEF_VECT) - fact2(DEF_VECT)                      ! Auu_zz
                   end do
                end do
             end do
          end if
          !
          ! Low-Mach with skew symmetric (adding temporal term: d(rho*u)/dt)
          !
          if (kfl_regim_nsi == 3_ip) then
             do itime =2, nbdfp_nsi         ! only rhs with Galerkin needs to be modified
                do igaus =1, pgaus
                   gpveo(DEF_VECT,1:ndime) = 0.0_rp  ! velocity at itime 
                   do inode =1, pnode
                      do idime=1,ndime
                         gpveo(DEF_VECT,idime) = gpveo(DEF_VECT,idime) + elvel(DEF_VECT,idime,inode,itime) * gpsha(DEF_VECT,inode,igaus)
                      end do
                   end do
                   fact0(DEF_VECT) = 0.5_rp * (gpden(DEF_VECT,igaus) - densi(DEF_VECT,igaus, itime))*pabdf_nsi(itime) & 
                                   * dtinv_loc(DEF_VECT) * gpvol(DEF_VECT,igaus)
                   do inode =1, pnode
                      do idime=1,ndime
                         elrbu(DEF_VECT,idime,inode) = elrbu(DEF_VECT,idime, inode) &
                                                     + fact0(DEF_VECT) * gpsha(DEF_VECT,inode,igaus) * gpveo(DEF_VECT,idime)
                      end do
                   end do
                end do
             end do
          end if

       end if

       !----------------------------------------------------------------------
       !
       ! Auu: Off-diagonal terms
       !
       !----------------------------------------------------------------------

       if( kfl_rmom2_nsi /= 0 ) then
          !
          ! p1vec * rmom2(DEF_VECT,
          !
          !yor! [Originally loop on line 577] Reordered loops to solve indirections problems
          if (ndime == 2) then
             do igaus = 1,pgaus
                do jnode = 1,pnode
                   do inode = 1,pnode
                      elauu11(DEF_VECT,inode,jnode) = elauu11(DEF_VECT,inode,jnode) + p1vec(DEF_VECT,inode,igaus)*rmom2(DEF_VECT,1,1,jnode,igaus)
                      elauu21(DEF_VECT,inode,jnode) = elauu21(DEF_VECT,inode,jnode) + p1vec(DEF_VECT,inode,igaus)*rmom2(DEF_VECT,2,1,jnode,igaus)
                      elauu12(DEF_VECT,inode,jnode) = elauu12(DEF_VECT,inode,jnode) + p1vec(DEF_VECT,inode,igaus)*rmom2(DEF_VECT,1,2,jnode,igaus)
                      elauu22(DEF_VECT,inode,jnode) = elauu22(DEF_VECT,inode,jnode) + p1vec(DEF_VECT,inode,igaus)*rmom2(DEF_VECT,2,2,jnode,igaus)
                   end do
                end do
             end do
          else !ndime==3
             do igaus = 1,pgaus
                factvec1(DEF_VECT,1:pnode) = p1vec(DEF_VECT,1:pnode,igaus)
                do jnode = 1,pnode
                   do inode = 1,pnode
                      elauu11(DEF_VECT,inode,jnode) = elauu11(DEF_VECT,inode,jnode) + factvec1(DEF_VECT,inode)*rmom2(DEF_VECT,1,1,jnode,igaus)
                      elauu21(DEF_VECT,inode,jnode) = elauu21(DEF_VECT,inode,jnode) + factvec1(DEF_VECT,inode)*rmom2(DEF_VECT,2,1,jnode,igaus)
                      elauu31(DEF_VECT,inode,jnode) = elauu31(DEF_VECT,inode,jnode) + factvec1(DEF_VECT,inode)*rmom2(DEF_VECT,3,1,jnode,igaus)
                      elauu12(DEF_VECT,inode,jnode) = elauu12(DEF_VECT,inode,jnode) + factvec1(DEF_VECT,inode)*rmom2(DEF_VECT,1,2,jnode,igaus)
                      elauu22(DEF_VECT,inode,jnode) = elauu22(DEF_VECT,inode,jnode) + factvec1(DEF_VECT,inode)*rmom2(DEF_VECT,2,2,jnode,igaus)
                      elauu32(DEF_VECT,inode,jnode) = elauu32(DEF_VECT,inode,jnode) + factvec1(DEF_VECT,inode)*rmom2(DEF_VECT,3,2,jnode,igaus)
                      elauu13(DEF_VECT,inode,jnode) = elauu13(DEF_VECT,inode,jnode) + factvec1(DEF_VECT,inode)*rmom2(DEF_VECT,1,3,jnode,igaus)
                      elauu23(DEF_VECT,inode,jnode) = elauu23(DEF_VECT,inode,jnode) + factvec1(DEF_VECT,inode)*rmom2(DEF_VECT,2,3,jnode,igaus)
                      elauu33(DEF_VECT,inode,jnode) = elauu33(DEF_VECT,inode,jnode) + factvec1(DEF_VECT,inode)*rmom2(DEF_VECT,3,3,jnode,igaus)
                   end do
                end do
             end do
          end if

          if( kfl_p1ve2_nsi /= 0 ) then
             !
             ! p1ve2 * rmomu
             ! p1ve2 * rmom2(DEF_VECT,
             !
             !yor! optim : [originally line 597] reverted loops, refactored
             if(ndime==2)then
                do igaus=1,pgaus
                   do jnode=1,pnode
                      do inode=1,pnode
                         elauu11(DEF_VECT,inode,jnode)=elauu11(DEF_VECT,inode,jnode)+p1ve2(DEF_VECT,1,1,inode,igaus)*rmomu(DEF_VECT,jnode,igaus)
                         elauu21(DEF_VECT,inode,jnode)=elauu21(DEF_VECT,inode,jnode)+p1ve2(DEF_VECT,2,1,inode,igaus)*rmomu(DEF_VECT,jnode,igaus)
                         elauu12(DEF_VECT,inode,jnode)=elauu12(DEF_VECT,inode,jnode)+p1ve2(DEF_VECT,1,2,inode,igaus)*rmomu(DEF_VECT,jnode,igaus)
                         elauu22(DEF_VECT,inode,jnode)=elauu22(DEF_VECT,inode,jnode)+p1ve2(DEF_VECT,2,2,inode,igaus)*rmomu(DEF_VECT,jnode,igaus)
                         do kdime = 1,ndime
                            elauu11(DEF_VECT,inode,jnode)=elauu11(DEF_VECT,inode,jnode)+p1ve2(DEF_VECT,1,kdime,inode,igaus)*rmom2(DEF_VECT,kdime,1,jnode,igaus)
                            elauu21(DEF_VECT,inode,jnode)=elauu21(DEF_VECT,inode,jnode)+p1ve2(DEF_VECT,2,kdime,inode,igaus)*rmom2(DEF_VECT,kdime,1,jnode,igaus)
                            elauu12(DEF_VECT,inode,jnode)=elauu12(DEF_VECT,inode,jnode)+p1ve2(DEF_VECT,1,kdime,inode,igaus)*rmom2(DEF_VECT,kdime,2,jnode,igaus)
                            elauu22(DEF_VECT,inode,jnode)=elauu22(DEF_VECT,inode,jnode)+p1ve2(DEF_VECT,2,kdime,inode,igaus)*rmom2(DEF_VECT,kdime,2,jnode,igaus)
                         end do
                      end do
                   end do
                end do
             else ! ndime==3
                do igaus=1,pgaus
                   do jnode=1,pnode
                      do inode=1,pnode
                         elauu11(DEF_VECT,inode,jnode)=elauu11(DEF_VECT,inode,jnode)+p1ve2(DEF_VECT,1,1,inode,igaus)*rmomu(DEF_VECT,jnode,igaus)
                         elauu21(DEF_VECT,inode,jnode)=elauu21(DEF_VECT,inode,jnode)+p1ve2(DEF_VECT,2,1,inode,igaus)*rmomu(DEF_VECT,jnode,igaus)
                         elauu31(DEF_VECT,inode,jnode)=elauu31(DEF_VECT,inode,jnode)+p1ve2(DEF_VECT,3,1,inode,igaus)*rmomu(DEF_VECT,jnode,igaus)
                         elauu12(DEF_VECT,inode,jnode)=elauu12(DEF_VECT,inode,jnode)+p1ve2(DEF_VECT,1,2,inode,igaus)*rmomu(DEF_VECT,jnode,igaus)
                         elauu22(DEF_VECT,inode,jnode)=elauu22(DEF_VECT,inode,jnode)+p1ve2(DEF_VECT,2,2,inode,igaus)*rmomu(DEF_VECT,jnode,igaus)
                         elauu32(DEF_VECT,inode,jnode)=elauu32(DEF_VECT,inode,jnode)+p1ve2(DEF_VECT,3,2,inode,igaus)*rmomu(DEF_VECT,jnode,igaus)
                         elauu13(DEF_VECT,inode,jnode)=elauu13(DEF_VECT,inode,jnode)+p1ve2(DEF_VECT,1,3,inode,igaus)*rmomu(DEF_VECT,jnode,igaus)
                         elauu23(DEF_VECT,inode,jnode)=elauu23(DEF_VECT,inode,jnode)+p1ve2(DEF_VECT,2,3,inode,igaus)*rmomu(DEF_VECT,jnode,igaus)
                         elauu33(DEF_VECT,inode,jnode)=elauu33(DEF_VECT,inode,jnode)+p1ve2(DEF_VECT,3,3,inode,igaus)*rmomu(DEF_VECT,jnode,igaus)
                         do kdime = 1,ndime
                            elauu11(DEF_VECT,inode,jnode)=elauu11(DEF_VECT,inode,jnode)+p1ve2(DEF_VECT,1,kdime,inode,igaus)*rmom2(DEF_VECT,kdime,1,jnode,igaus)
                            elauu21(DEF_VECT,inode,jnode)=elauu21(DEF_VECT,inode,jnode)+p1ve2(DEF_VECT,2,kdime,inode,igaus)*rmom2(DEF_VECT,kdime,1,jnode,igaus)
                            elauu31(DEF_VECT,inode,jnode)=elauu31(DEF_VECT,inode,jnode)+p1ve2(DEF_VECT,3,kdime,inode,igaus)*rmom2(DEF_VECT,kdime,1,jnode,igaus)
                            elauu12(DEF_VECT,inode,jnode)=elauu12(DEF_VECT,inode,jnode)+p1ve2(DEF_VECT,1,kdime,inode,igaus)*rmom2(DEF_VECT,kdime,2,jnode,igaus)
                            elauu22(DEF_VECT,inode,jnode)=elauu22(DEF_VECT,inode,jnode)+p1ve2(DEF_VECT,2,kdime,inode,igaus)*rmom2(DEF_VECT,kdime,2,jnode,igaus)
                            elauu32(DEF_VECT,inode,jnode)=elauu32(DEF_VECT,inode,jnode)+p1ve2(DEF_VECT,3,kdime,inode,igaus)*rmom2(DEF_VECT,kdime,2,jnode,igaus)
                            elauu13(DEF_VECT,inode,jnode)=elauu13(DEF_VECT,inode,jnode)+p1ve2(DEF_VECT,1,kdime,inode,igaus)*rmom2(DEF_VECT,kdime,3,jnode,igaus)
                            elauu23(DEF_VECT,inode,jnode)=elauu23(DEF_VECT,inode,jnode)+p1ve2(DEF_VECT,2,kdime,inode,igaus)*rmom2(DEF_VECT,kdime,3,jnode,igaus)
                            elauu33(DEF_VECT,inode,jnode)=elauu33(DEF_VECT,inode,jnode)+p1ve2(DEF_VECT,3,kdime,inode,igaus)*rmom2(DEF_VECT,kdime,3,jnode,igaus)
                         end do
                      end do
                   end do
                end do
             end if
          end if
       end if
       !
       ! Bemol:
       ! - b * ( v , rho * (a.grad) u ) - b * ( u , rho * (a.grad) v ) - b ( rho*div(a) u , v )
       !
       if( abs(bemol_nsi).gt.1.0e-9_rp ) then
          if( ndime == 2 ) then
             do igaus = 1,pgaus
                fact0(DEF_VECT) = bemol_nsi * gpden(DEF_VECT,igaus) * gpvol(DEF_VECT,igaus)
                fact5(DEF_VECT) = fact0(DEF_VECT) * ( gpgve(DEF_VECT,1,1,igaus) + gpgve(DEF_VECT,2,2,igaus) )
                do inode = 1,pnode
                   do jnode = 1,pnode
                      !
                      ! fact2(DEF_VECT) = - b * ( v , rho * (a.grad) u )
                      ! fact3(DEF_VECT) = - b * ( u , rho * (a.grad) v )
                      ! fact4(DEF_VECT) = - b * ( rho * div(a) u , v )
                      !
                      fact2(DEF_VECT) = ( gpadv(DEF_VECT,1,igaus) * gpcar(DEF_VECT,1,inode,igaus) + gpadv(DEF_VECT,2,igaus) * gpcar(DEF_VECT,2,inode,igaus) ) * gpsha(DEF_VECT,jnode,igaus)
                      fact3(DEF_VECT) = ( gpadv(DEF_VECT,1,igaus) * gpcar(DEF_VECT,1,jnode,igaus) + gpadv(DEF_VECT,2,igaus) * gpcar(DEF_VECT,2,jnode,igaus) ) * gpsha(DEF_VECT,inode,igaus)
                      fact4(DEF_VECT) = fact5(DEF_VECT) * gpsha(DEF_VECT,inode,igaus) * gpsha(DEF_VECT,jnode,igaus)
                      fact2(DEF_VECT) = fact0(DEF_VECT) * fact2(DEF_VECT)
                      fact3(DEF_VECT) = fact0(DEF_VECT) * fact3(DEF_VECT)

                      elauu11(DEF_VECT,inode,jnode) = elauu11(DEF_VECT,inode,jnode) - fact2(DEF_VECT) - fact3(DEF_VECT) - fact4(DEF_VECT)
                      elauu22(DEF_VECT,inode,jnode) = elauu22(DEF_VECT,inode,jnode) - fact2(DEF_VECT) - fact3(DEF_VECT) - fact4(DEF_VECT)
                   end do
                end do
             end do
          else
             do igaus = 1,pgaus
                fact0(DEF_VECT) = bemol_nsi * gpden(DEF_VECT,igaus) * gpvol(DEF_VECT,igaus)
                fact5(DEF_VECT) = fact0(DEF_VECT) * ( gpgve(DEF_VECT,1,1,igaus) + gpgve(DEF_VECT,2,2,igaus) + gpgve(DEF_VECT,3,3,igaus))
                do inode = 1,pnode
                   do jnode = 1,pnode
                      !
                      ! fact2(DEF_VECT) = - b * ( v , rho * (a.grad) u )
                      ! fact3(DEF_VECT) = - b * ( u , rho * (a.grad) v )
                      ! fact4(DEF_VECT) = - b * ( rho * div(a) u , v )
                      !
                      fact2(DEF_VECT) = ( gpadv(DEF_VECT,1,igaus) * gpcar(DEF_VECT,1,inode,igaus) + gpadv(DEF_VECT,2,igaus) * gpcar(DEF_VECT,2,inode,igaus) &
                           + gpadv(DEF_VECT,3,igaus) * gpcar(DEF_VECT,3,inode,igaus)) * gpsha(DEF_VECT,jnode,igaus)
                      fact3(DEF_VECT) = ( gpadv(DEF_VECT,1,igaus) * gpcar(DEF_VECT,1,jnode,igaus) + gpadv(DEF_VECT,2,igaus) * gpcar(DEF_VECT,2,jnode,igaus) &
                           + gpadv(DEF_VECT,3,igaus) * gpcar(DEF_VECT,3,jnode,igaus)) * gpsha(DEF_VECT,inode,igaus)
                      fact4(DEF_VECT) = fact5(DEF_VECT) * gpsha(DEF_VECT,inode,igaus) * gpsha(DEF_VECT,jnode,igaus)
                      fact2(DEF_VECT) = fact0(DEF_VECT) * fact2(DEF_VECT)
                      fact3(DEF_VECT) = fact0(DEF_VECT) * fact3(DEF_VECT)

                      elauu11(DEF_VECT,inode,jnode) = elauu11(DEF_VECT,inode,jnode) - fact2(DEF_VECT) - fact3(DEF_VECT) - fact4(DEF_VECT)
                      elauu22(DEF_VECT,inode,jnode) = elauu22(DEF_VECT,inode,jnode) - fact2(DEF_VECT) - fact3(DEF_VECT) - fact4(DEF_VECT)
                      elauu33(DEF_VECT,inode,jnode) = elauu33(DEF_VECT,inode,jnode) - fact2(DEF_VECT) - fact3(DEF_VECT) - fact4(DEF_VECT)
                   end do
                end do
             end do
          end if

       end if
       !
       ! Newton-Raphson
       !
       if( kfl_linea_nsi == 2 ) then
          if( ndime == 2 ) then
             do igaus = 1,pgaus
                do inode = 1,pnode
                   fact0(DEF_VECT) = gpvol(DEF_VECT,igaus) * gpden(DEF_VECT,igaus) * gpsha(DEF_VECT,inode,igaus)
                   do jnode = 1,pnode
                      fact1(DEF_VECT) = fact0(DEF_VECT) * gpsha(DEF_VECT,jnode,igaus)
                      ! Auu_1i
                      elauu11(DEF_VECT,inode,jnode) = elauu11(DEF_VECT,inode,jnode) + fact1(DEF_VECT) * gpgve(DEF_VECT,1,1,igaus) ! rho * ux * dux^i-1/dx
                      elauu12(DEF_VECT,inode,jnode) = elauu12(DEF_VECT,inode,jnode) + fact1(DEF_VECT) * gpgve(DEF_VECT,2,1,igaus) ! rho * uy * dux^i-1/dy
                      ! Auu_2i
                      elauu21(DEF_VECT,inode,jnode) = elauu21(DEF_VECT,inode,jnode) + fact1(DEF_VECT) * gpgve(DEF_VECT,1,2,igaus) ! rho * ux * duy^i-1/dx
                      elauu22(DEF_VECT,inode,jnode) = elauu22(DEF_VECT,inode,jnode) + fact1(DEF_VECT) * gpgve(DEF_VECT,2,2,igaus) ! rho * uy * duy^i-1/dy
                   end do
                end do
             end do
          else
             do igaus = 1,pgaus
                do inode = 1,pnode
                   fact0(DEF_VECT) = gpvol(DEF_VECT,igaus) * gpden(DEF_VECT,igaus) * gpsha(DEF_VECT,inode,igaus)
                   do jnode = 1,pnode
                      fact1(DEF_VECT) = fact0(DEF_VECT) * gpsha(DEF_VECT,jnode,igaus)
                      ! Auu_1i
                      elauu11(DEF_VECT,inode,jnode) = elauu11(DEF_VECT,inode,jnode) + fact1(DEF_VECT) * gpgve(DEF_VECT,1,1,igaus) ! rho * ux * dux^i-1/dx
                      elauu12(DEF_VECT,inode,jnode) = elauu12(DEF_VECT,inode,jnode) + fact1(DEF_VECT) * gpgve(DEF_VECT,2,1,igaus) ! rho * uy * dux^i-1/dy
                      elauu13(DEF_VECT,inode,jnode) = elauu13(DEF_VECT,inode,jnode) + fact1(DEF_VECT) * gpgve(DEF_VECT,3,1,igaus) ! rho * uz * dux^i-1/dz
                      ! Auu_2i
                      elauu21(DEF_VECT,inode,jnode) = elauu21(DEF_VECT,inode,jnode) + fact1(DEF_VECT) * gpgve(DEF_VECT,1,2,igaus) ! rho * ux * duy^i-1/dx
                      elauu22(DEF_VECT,inode,jnode) = elauu22(DEF_VECT,inode,jnode) + fact1(DEF_VECT) * gpgve(DEF_VECT,2,2,igaus) ! rho * uy * duy^i-1/dy
                      elauu23(DEF_VECT,inode,jnode) = elauu23(DEF_VECT,inode,jnode) + fact1(DEF_VECT) * gpgve(DEF_VECT,3,2,igaus) ! rho * uz * duy^i-1/dz
                      ! Auu_3i
                      elauu31(DEF_VECT,inode,jnode) = elauu31(DEF_VECT,inode,jnode) + fact1(DEF_VECT) * gpgve(DEF_VECT,1,3,igaus) ! rho * ux * duz^i-1/dx
                      elauu32(DEF_VECT,inode,jnode) = elauu32(DEF_VECT,inode,jnode) + fact1(DEF_VECT) * gpgve(DEF_VECT,2,3,igaus) ! rho * uy * duz^i-1/dy
                      elauu33(DEF_VECT,inode,jnode) = elauu33(DEF_VECT,inode,jnode) + fact1(DEF_VECT) * gpgve(DEF_VECT,3,3,igaus) ! rho * uz * duz^i-1/dz

                   end do
                end do
             end do
          end if
       end if
       !
       ! Lumped evolution matrix (only backward euler)
       !
       if( kfl_lumped == 1 ) then
          !
          ! Remove Galerkin term and add lumped term
          !
          if( ndime == 2 ) then
             call runend('PREGUNTAR A MATIAS QUE LO PROGRAME')
          else
             do igaus =1, pgaus
                gpveo(DEF_VECT,1:3) = 0.0_rp
                do inode = 1,pnode
                   gpveo(DEF_VECT,1) = gpveo(DEF_VECT,1) + elvel(DEF_VECT,1, inode, 2)* gpsha(DEF_VECT,inode, igaus)
                   gpveo(DEF_VECT,2) = gpveo(DEF_VECT,2) + elvel(DEF_VECT,2, inode, 2)* gpsha(DEF_VECT,inode, igaus)
                   gpveo(DEF_VECT,3) = gpveo(DEF_VECT,3) + elvel(DEF_VECT,3, inode, 2)* gpsha(DEF_VECT,inode, igaus)
                end do
                do inode =1, pnode
                   fact0(DEF_VECT) = gpvol(DEF_VECT,igaus) * gpden(DEF_VECT,igaus) * gpsha(DEF_VECT,inode,igaus)*dtinv_mod(DEF_VECT)
                   elauu11(DEF_VECT,inode,inode) =  elauu11(DEF_VECT,inode,inode) + fact0(DEF_VECT)
                   elauu22(DEF_VECT,inode,inode) =  elauu22(DEF_VECT,inode,inode) + fact0(DEF_VECT)
                   elauu33(DEF_VECT,inode,inode) =  elauu33(DEF_VECT,inode,inode) + fact0(DEF_VECT)
                   do idime = 1,ndime
                      elrbu(DEF_VECT,idime,inode)     =  elrbu(DEF_VECT,idime,inode)   - fact0(DEF_VECT)*gpveo(DEF_VECT,idime)
                      elrbu(DEF_VECT,idime,inode)     =  elrbu(DEF_VECT,idime,inode)   + fact0(DEF_VECT)*elvel(DEF_VECT,idime, inode, 2)
                   end do
                   do jnode =1, pnode !yor! for this case, better write a second loop and revert inode and jnode
                      elauu11(DEF_VECT,inode,jnode) =  elauu11(DEF_VECT,inode,jnode) - fact0(DEF_VECT)*gpsha(DEF_VECT,jnode, igaus)
                      elauu22(DEF_VECT,inode,jnode) =  elauu22(DEF_VECT,inode,jnode) - fact0(DEF_VECT)*gpsha(DEF_VECT,jnode, igaus)
                      elauu33(DEF_VECT,inode,jnode) =  elauu33(DEF_VECT,inode,jnode) - fact0(DEF_VECT)*gpsha(DEF_VECT,jnode, igaus)
                   end do
                end do
             end do
          end if
       else if( kfl_lumped == 2 ) then
          !
          ! No time term have been added up to now: add Galerkin term
          !
          do igaus = 1,pgaus
             fact0(DEF_VECT) = gpvol(DEF_VECT,igaus) * gpden(DEF_VECT,igaus) * dtinv_mod(DEF_VECT)
             if( ndime == 2 ) then
                do inode = 1, pnode
                   elauu11(DEF_VECT,inode,inode) = elauu11(DEF_VECT,inode,inode) + fact0(DEF_VECT) * gpsha(DEF_VECT,inode,igaus)
                   elauu22(DEF_VECT,inode,inode) = elauu22(DEF_VECT,inode,inode) + fact0(DEF_VECT) * gpsha(DEF_VECT,inode,igaus)
                end do
             else ! ndime==3
                do inode = 1, pnode
                   elauu11(DEF_VECT,inode,inode) = elauu11(DEF_VECT,inode,inode) + fact0(DEF_VECT) * gpsha(DEF_VECT,inode,igaus)
                   elauu22(DEF_VECT,inode,inode) = elauu22(DEF_VECT,inode,inode) + fact0(DEF_VECT) * gpsha(DEF_VECT,inode,igaus)
                   elauu33(DEF_VECT,inode,inode) = elauu33(DEF_VECT,inode,inode) + fact0(DEF_VECT) * gpsha(DEF_VECT,inode,igaus)
                end do
             end if
             do inode = 1, pnode
                do idime = 1,ndime
                   elrbu(DEF_VECT,idime,inode) = elrbu(DEF_VECT,idime,inode) + fact0(DEF_VECT) * gpsha(DEF_VECT,inode,igaus) * elvel(DEF_VECT,idime,inode,2)
                end do
             end do
          end do
       end if
       !
       ! Dual time step preconditioner
       !
       if( kfl_duatss == 1 ) then
          ellum = 0.0_rp
          do igaus =1, pgaus
             do inode =1, pnode
                fact0(DEF_VECT)       = gpvol(DEF_VECT,igaus) * gpden(DEF_VECT,igaus) * gpsha(DEF_VECT,inode,igaus) * dtinv_loc(DEF_VECT) * real(fact_duatss-1,rp)
                ellum(DEF_VECT,inode) = ellum(DEF_VECT,inode) + fact0(DEF_VECT)
             end do
          end do
       end if

       !yor! elauu final matrix assembly
       if( ndime == 2 ) then
          do jnode = 1,pnode
             do inode = 1,pnode
                elauu(DEF_VECT,2*inode-1,2*jnode-1)=elauu11(DEF_VECT,inode,jnode)
                elauu(DEF_VECT,2*inode  ,2*jnode-1)=elauu21(DEF_VECT,inode,jnode)
                elauu(DEF_VECT,2*inode-1,2*jnode  )=elauu12(DEF_VECT,inode,jnode)
                elauu(DEF_VECT,2*inode  ,2*jnode  )=elauu22(DEF_VECT,inode,jnode)
             end do
          end do
       else ! ndime==3
          do jnode = 1,pnode
             do inode = 1,pnode
                elauu(DEF_VECT,3*inode-2,3*jnode-2)=elauu11(DEF_VECT,inode,jnode)
                elauu(DEF_VECT,3*inode-1,3*jnode-2)=elauu21(DEF_VECT,inode,jnode)
                elauu(DEF_VECT,3*inode  ,3*jnode-2)=elauu31(DEF_VECT,inode,jnode)
                elauu(DEF_VECT,3*inode-2,3*jnode-1)=elauu12(DEF_VECT,inode,jnode)
                elauu(DEF_VECT,3*inode-1,3*jnode-1)=elauu22(DEF_VECT,inode,jnode)
                elauu(DEF_VECT,3*inode  ,3*jnode-1)=elauu32(DEF_VECT,inode,jnode)
                elauu(DEF_VECT,3*inode-2,3*jnode  )=elauu13(DEF_VECT,inode,jnode)
                elauu(DEF_VECT,3*inode-1,3*jnode  )=elauu23(DEF_VECT,inode,jnode)
                elauu(DEF_VECT,3*inode  ,3*jnode  )=elauu33(DEF_VECT,inode,jnode)
             end do
          end do
       end if
       !----------------------------------------------------------------------
       !
       ! Aup
       !
       !----------------------------------------------------------------------
       !
       ! Pressure: - ( p , div(v)  ) + ( grad(p) , p1vec(DEF_VECT,v)-v ) ( if kfl_press_nsi = 1 )
       ! Pressure: + ( grad(p) , p1vec(DEF_VECT,v) )                     ( if kfl_press_nsi = 0 )
       !
       if( kfl_press_nsi == 1 ) then
          do igaus = 1,pgaus
             do jnode = 1,pnode
                fact2(DEF_VECT) = gpvol(DEF_VECT,igaus)*gpsha(DEF_VECT,jnode,igaus)
                do inode = 1,pnode
                   fact1(DEF_VECT) = -gpvol(DEF_VECT,igaus)*gpsha(DEF_VECT,inode,igaus)+p1vec(DEF_VECT,inode,igaus)
                   idof1 = (inode-1)*ndime
                   do idime = 1,ndime
                      idofv = idof1 + idime
                      elaup(DEF_VECT,idofv,jnode) = elaup(DEF_VECT,idofv,jnode) - fact2(DEF_VECT) * gpcar(DEF_VECT,idime,inode,igaus) &
                           &                                  + fact1(DEF_VECT) * gpcar(DEF_VECT,idime,jnode,igaus)
                   end do
                end do
             end do
          end do
       else
          do igaus = 1,pgaus
             do jnode = 1,pnode
                do inode = 1,pnode
                   idof1 = (inode-1)*ndime
                   do idime = 1,ndime
                      idofv = idof1 + idime
                      elaup(DEF_VECT,idofv,jnode) = elaup(DEF_VECT,idofv,jnode) + p1vec(DEF_VECT,inode,igaus) * gpcar(DEF_VECT,idime,jnode,igaus)
                   end do
                end do
             end do
          end do
       end if
       if( kfl_p1ve2_nsi /= 0 ) then 
          if (NSI_FRACTIONAL_STEP) then ! Residual components: 
             !       pressure gradient and temporal derivative term on rhs of momentum
             do igaus = 1,pgaus
                grpre =0.0_rp ! pressure gradient at igaus
                dudt = 0.0_rp ! temporal derivative of velocity at igaus
                do idime =1, ndime
                   do inode =1, pnode
                      grpre(DEF_VECT,idime) =  grpre(DEF_VECT,idime) + elpre(DEF_VECT,inode,1) * gpcar(DEF_VECT,idime,inode,igaus)
                      dudt (DEF_VECT,idime)  = dudt( DEF_VECT,idime) + (elvel(DEF_VECT,1, inode, 2) - elvel(DEF_VECT,1, inode, 3))* gpsha(DEF_VECT,inode, igaus)
                   end do
                   dudt (DEF_VECT,idime)  = dudt (DEF_VECT,idime) *dtinv_loc(DEF_VECT)
                end do               
                do inode = 1,pnode                 
                   do idime = 1,ndime
                      do jdime = 1,ndime
                         elrbu(DEF_VECT,idime, inode) = elrbu(DEF_VECT,idime, inode) - &
                             ( grpre(DEF_VECT,jdime)+ gpden(DEF_VECT,igaus) * dudt(DEF_VECT,jdime))*p1ve2(DEF_VECT,idime,jdime,inode,igaus)
                      end do
                   end do
                   
                end do
             end do
          else
             do igaus = 1,pgaus
                do jnode = 1,pnode
                   do inode = 1,pnode
                      idof1 = (inode-1)*ndime
                      do idime = 1,ndime
                         idofv = idof1 + idime
                         do jdime = 1,ndime
                            elaup(DEF_VECT,idofv,jnode) = elaup(DEF_VECT,idofv,jnode) + &
                                 gpcar(DEF_VECT,jdime,jnode,igaus) * p1ve2(DEF_VECT,idime,jdime,inode,igaus)
                         end do
                      end do
                   end do
                end do
             end do
          end if
       end if

       !----------------------------------------------------------------------
       !
       ! Apu
       !
       !----------------------------------------------------------------------
       !
       ! ( div(u) , (tau2^{-1}*tau2')*q )   --- or ( div (rho u), q) if Low-Mach
       ! + ( rho*(uc.grad)u + 2*rho*(w x u) + sig*u -div[2*mu*eps(u)], tau1' grad(q) )  (only laplacian form of div[2 mu eps(u)] )
       !
       if (kfl_regim_nsi == 3 .and. kfl_confi_nsi == 1) then ! Penalization of pressure, never imposse pressure
          do igaus =1, pgaus
             penal(DEF_VECT) = 1.0e-4_rp*gpden(DEF_VECT,igaus)/gpvis(DEF_VECT,igaus)
             do inode =1, pnode
                fact0(DEF_VECT) = penal(DEF_VECT) * p2sca(DEF_VECT,inode,igaus)
                do jnode =inode +1, pnode
                   fact1(DEF_VECT) = fact0(DEF_VECT)*gpsha(DEF_VECT,jnode,igaus)
                   elapp(DEF_VECT,inode,jnode) = elapp(DEF_VECT,inode,jnode) + fact1(DEF_VECT)
                   elapp(DEF_VECT,jnode,inode) = elapp(DEF_VECT,jnode,inode) + fact1(DEF_VECT)
                end do
                elapp(DEF_VECT,inode,inode) = elapp(DEF_VECT,inode,inode) + fact0(DEF_VECT) * gpsha(DEF_VECT,inode, igaus)
                elrbp(DEF_VECT,inode)       = elrbp(DEF_VECT,inode)       + fact0(DEF_VECT) * gppre(DEF_VECT,igaus)
             end do
          end do
       end if

#ifdef matiaslma

       if (ndime==2) then
          do igaus =1, pgaus
             do inode =1,pnode
                p2vec(DEF_VECT,1,inode,igaus)= p2vec(DEF_VECT,1,inode,igaus)*gpden(DEF_VECT,igaus)
                p2vec(DEF_VECT,2,inode,igaus)= p2vec(DEF_VECT,2,inode,igaus)*gpden(DEF_VECT,igaus)
                p2sca(DEF_VECT,inode, igaus) = p2sca(DEF_VECT,inode, igaus)*gpden(DEF_VECT,igaus)
             end do
             fact0(DEF_VECT) = gpden(DEF_VECT,igaus)* gpvol(DEF_VECT,igaus)
             do inode =1, pnode
                idof1 = 2*inode-1
                idof2 = idof1+1
                fact1(DEF_VECT) = gpsha(DEF_VECT,inode,igaus)* fact0(DEF_VECT)
                do jnode =1, pnode
                   elapu(DEF_VECT,jnode, idof1) = elapu(DEF_VECT,jnode, idof1) -fact1(DEF_VECT)             * gpcar(DEF_VECT,1,jnode,igaus)  &
                        &                                  + rmomu(DEF_VECT,inode,igaus) * p2vec(DEF_VECT,1,jnode,igaus)
                   elapu(DEF_VECT,jnode, idof2) = elapu(DEF_VECT,jnode, idof2) -fact1(DEF_VECT)             * gpcar(DEF_VECT,2,jnode,igaus)  &
                        &                                  + rmomu(DEF_VECT,inode,igaus) * p2vec(DEF_VECT,2,jnode,igaus)
                end do
             end do
          end do
       else !ndime =3
          do igaus =1, pgaus
             do inode =1,pnode
                p2vec(DEF_VECT,1:ndime,inode,igaus) = p2vec(DEF_VECT,1:ndime,inode,igaus)*gpden(DEF_VECT,igaus)
                p2sca(DEF_VECT,inode, igaus)        = p2sca(DEF_VECT,inode, igaus) *gpden(DEF_VECT,igaus)
             end do
             fact0(DEF_VECT) = gpden(DEF_VECT,igaus)* gpvol(DEF_VECT,igaus)
             do inode =1, pnode
                idof1 = 3*inode - 2
                idof2 = idof1+1
                idof3 = idof2+1
                fact1(DEF_VECT) = gpsha(DEF_VECT,inode,igaus)* fact0(DEF_VECT)
                do jnode =1, pnode
                   elapu(DEF_VECT,jnode, idof1) = elapu(DEF_VECT,jnode, idof1) - fact1(DEF_VECT)              * gpcar(DEF_VECT,1,jnode,igaus) &
                        &                                    + rmomu(DEF_VECT,inode,igaus) * p2vec(DEF_VECT,1,jnode,igaus)
                   elapu(DEF_VECT,jnode, idof2) = elapu(DEF_VECT,jnode, idof2) - fact1(DEF_VECT)              * gpcar(DEF_VECT,2,jnode,igaus) &
                        &                                    + rmomu(DEF_VECT,inode,igaus) * p2vec(DEF_VECT,2,jnode,igaus)
                   elapu(DEF_VECT,jnode, idof3) = elapu(DEF_VECT,jnode, idof3) - fact1(DEF_VECT)              * gpcar(DEF_VECT,3,jnode,igaus) &
                        &                                    + rmomu(DEF_VECT,inode,igaus) * p2vec(DEF_VECT,3,jnode,igaus)
                end do
             end do
          end do
       end if

#else
       do igaus = 1,pgaus
          do inode = 1,pnode
             do jnode = 1,pnode
                do idime = 1,ndime
                   idof1 = (inode-1)*ndime+idime
                   elapu(DEF_VECT,jnode,idof1) = elapu(DEF_VECT,jnode,idof1) + rcont(DEF_VECT,idime,inode,igaus) * p2sca(DEF_VECT,jnode,igaus)&
                        &                                  + rmomu(DEF_VECT,inode,igaus)       * p2vec(DEF_VECT,idime,jnode,igaus)
                end do
             end do
          end do
       end do

#endif
       !yor! optim [originally loop on line 985] reordered inner loops
       if( kfl_rmom2_nsi /= 0 ) then
          do igaus = 1,pgaus
             do jnode = 1,pnode
                do jdime = 1,ndime
                   jdofv = (jnode-1)*ndime+jdime
                   do inode = 1,pnode
                      do kdime = 1,ndime
                         elapu(DEF_VECT,inode,jdofv) = elapu(DEF_VECT,inode,jdofv) + rmom2(DEF_VECT,kdime,jdime,jnode,igaus) * p2vec(DEF_VECT,kdime,inode,igaus)
                      end do
                   end do
                end do
             end do
          end do
       end if

       !----------------------------------------------------------------------
       !
       ! App
       !
       !----------------------------------------------------------------------
       !
       ! Pressure: ( grad(p) , tau1' grad(q) )
       !
       do igaus = 1,pgaus
          do inode = 1,pnode
             do jnode = 1,pnode
                do idime = 1,ndime
                   elapp(DEF_VECT,jnode,inode) = elapp(DEF_VECT,jnode,inode) + p2vec(DEF_VECT,idime,jnode,igaus) * gpcar(DEF_VECT,idime,inode,igaus)
                end do
             end do
          end do
       end do
       !call nsi_element_system_output(&
       !     pnode,elauu(1,:,:),elaup(1,:,:),elapp(1,:,:),elapu(1,:,:),elrbu(1,:,:),elrbp(1,:),&
       !     elauq(1,:,:),elapq(1,:,:),elaqu(1,:,:),elaqp(1,:,:),elaqq(1,:,:),elrbq(1,:))
       !stop

       !
       ! Penalization
       !
       do igaus = 1,pgaus
          fact1(DEF_VECT) = penal_nsi * gpvol(DEF_VECT,igaus)
          do inode = 1,pnode
             elapp(DEF_VECT,inode,inode) = elapp(DEF_VECT,inode,inode) + fact1(DEF_VECT) * gpsha(DEF_VECT,inode,igaus)
             elrbp(DEF_VECT,inode)       = elrbp(DEF_VECT,inode)       + fact1(DEF_VECT) * gpsha(DEF_VECT,inode,igaus) * elpre(DEF_VECT,inode,1)
          end do
       end do

       !----------------------------------------------------------------------
       !
       ! bu and bp
       !
       !----------------------------------------------------------------------
       !
       ! bu = ( f , p1vec(DEF_VECT,v) )
       ! bp = ( f , tau1' grad(q) )
       !
       if( ndime == 2 ) then
          do igaus = 1,pgaus
             gprhs_sgs(DEF_VECT,1,igaus) = gprhs_sgs(DEF_VECT,1,igaus) + taupr(DEF_VECT,igaus) * gpvep(DEF_VECT,1,igaus)
             gprhs_sgs(DEF_VECT,2,igaus) = gprhs_sgs(DEF_VECT,2,igaus) + taupr(DEF_VECT,igaus) * gpvep(DEF_VECT,2,igaus)
             gprhh(DEF_VECT,1,igaus)     = gprhs(DEF_VECT,1,igaus)     + gprhs_sgs(DEF_VECT,1,igaus)
             gprhh(DEF_VECT,2,igaus)     = gprhs(DEF_VECT,2,igaus)     + gprhs_sgs(DEF_VECT,2,igaus)
             do inode = 1,pnode
                elrbu(DEF_VECT,1,inode) = elrbu(DEF_VECT,1,inode) + p1vec(DEF_VECT,inode,igaus)   * gprhh(DEF_VECT,1,igaus)
                elrbu(DEF_VECT,2,inode) = elrbu(DEF_VECT,2,inode) + p1vec(DEF_VECT,inode,igaus)   * gprhh(DEF_VECT,2,igaus)
                elrbp(DEF_VECT,inode)   = elrbp(DEF_VECT,inode)   + p2vec(DEF_VECT,1,inode,igaus) * gprhh(DEF_VECT,1,igaus) &
                     &                                            + p2vec(DEF_VECT,2,inode,igaus) * gprhh(DEF_VECT,2,igaus)
             end do
          end do
       else
          do igaus = 1,pgaus
             gprhs_sgs(DEF_VECT,1,igaus) = gprhs_sgs(DEF_VECT,1,igaus) + taupr(DEF_VECT,igaus) * gpvep(DEF_VECT,1,igaus)
             gprhs_sgs(DEF_VECT,2,igaus) = gprhs_sgs(DEF_VECT,2,igaus) + taupr(DEF_VECT,igaus) * gpvep(DEF_VECT,2,igaus)
             gprhs_sgs(DEF_VECT,3,igaus) = gprhs_sgs(DEF_VECT,3,igaus) + taupr(DEF_VECT,igaus) * gpvep(DEF_VECT,3,igaus)
             gprhh(DEF_VECT,1, igaus)    = gprhs(DEF_VECT,1,igaus)     + gprhs_sgs(DEF_VECT,1,igaus)
             gprhh(DEF_VECT,2, igaus)    = gprhs(DEF_VECT,2,igaus)     + gprhs_sgs(DEF_VECT,2,igaus)
             gprhh(DEF_VECT,3, igaus)    = gprhs(DEF_VECT,3,igaus)     + gprhs_sgs(DEF_VECT,3,igaus)
             do inode = 1,pnode
                elrbu(DEF_VECT,1,inode) = elrbu(DEF_VECT,1,inode) + p1vec(DEF_VECT,inode,igaus)   * gprhh(DEF_VECT,1,igaus)
                elrbu(DEF_VECT,2,inode) = elrbu(DEF_VECT,2,inode) + p1vec(DEF_VECT,inode,igaus)   * gprhh(DEF_VECT,2,igaus)
                elrbu(DEF_VECT,3,inode) = elrbu(DEF_VECT,3,inode) + p1vec(DEF_VECT,inode,igaus)   * gprhh(DEF_VECT,3,igaus)
                elrbp(DEF_VECT,inode)   = elrbp(DEF_VECT,inode)   + p2vec(DEF_VECT,1,inode,igaus) * gprhh(DEF_VECT,1,igaus) &
                     &                                            + p2vec(DEF_VECT,2,inode,igaus) * gprhh(DEF_VECT,2,igaus) &
                     &                                            + p2vec(DEF_VECT,3,inode,igaus) * gprhh(DEF_VECT,3,igaus)
             end do
          end do
       end if
       if( kfl_p1ve2_nsi /= 0 ) then
          do igaus = 1,pgaus
             do inode = 1,pnode
                do idime = 1,ndime
                   do kdime = 1,ndime
                      elrbu(DEF_VECT,idime,inode) = elrbu(DEF_VECT,idime,inode) &
                           + p1ve2(DEF_VECT,idime,kdime,inode,igaus) * gprhh(DEF_VECT,kdime,igaus)
                   end do
                end do
             end do
          end do
       end if
       !
       ! low Mach regime: we add to the right hand side the residue of the continuity eq
       !
       do igaus = 1,pgaus
          fact0(DEF_VECT) = gprhc(DEF_VECT,igaus) * gpvol(DEF_VECT,igaus)
          fact1(DEF_VECT) = gpsp2(DEF_VECT,igaus) * gprh2(DEF_VECT,igaus) * gpvol(DEF_VECT,igaus)
          do inode = 1,pnode
             do idime = 1,ndime
                elrbu(DEF_VECT,idime,inode) = elrbu(DEF_VECT,idime,inode) + fact1(DEF_VECT)        * gpcar(DEF_VECT,idime,inode,igaus) ! ( rhs, tau2' * div(v) )
             end do
             elrbp(DEF_VECT,inode) = elrbp(DEF_VECT,inode) + fact0(DEF_VECT) * gpsha(DEF_VECT,inode,igaus)                             ! ( rhs , q)
          end do
       end do
       !
       ! Newton-Raphson
       !
       if( kfl_linea_nsi == 2 ) then
          if( ndime == 2 ) then
             do igaus = 1,pgaus
                fact0(DEF_VECT) =  gpden(DEF_VECT,igaus) * gpvol(DEF_VECT,igaus)
                do inode = 1,pnode
                   fact1(DEF_VECT)         = fact0(DEF_VECT) * gpsha(DEF_VECT,inode,igaus)
                   elrbu(DEF_VECT,1,inode) = elrbu(DEF_VECT,1,inode) + fact1(DEF_VECT) * ( gpadv(DEF_VECT,1,igaus) * gpgve(DEF_VECT,1,1,igaus) + gpadv(DEF_VECT,2,igaus) * gpgve(DEF_VECT,2,1,igaus) )
                   elrbu(DEF_VECT,2,inode) = elrbu(DEF_VECT,2,inode) + fact1(DEF_VECT) * ( gpadv(DEF_VECT,1,igaus) * gpgve(DEF_VECT,1,2,igaus) + gpadv(DEF_VECT,2,igaus) * gpgve(DEF_VECT,2,2,igaus) )
                end do
             end do
          else
             do igaus = 1,pgaus
                fact0(DEF_VECT) =  gpden(DEF_VECT,igaus) * gpvol(DEF_VECT,igaus)
                do inode = 1,pnode
                   fact1(DEF_VECT)         = fact0(DEF_VECT) * gpsha(DEF_VECT,inode,igaus)
                   elrbu(DEF_VECT,1,inode) = elrbu(DEF_VECT,1,inode) &
                        + fact1(DEF_VECT) * ( gpadv(DEF_VECT,1,igaus) * gpgve(DEF_VECT,1,1,igaus) + gpadv(DEF_VECT,2,igaus) * gpgve(DEF_VECT,2,1,igaus) &
                        + gpadv(DEF_VECT,3,igaus) * gpgve(DEF_VECT,3,1,igaus) )
                   elrbu(DEF_VECT,2,inode) = elrbu(DEF_VECT,2,inode) &
                        + fact1(DEF_VECT) * ( gpadv(DEF_VECT,1,igaus) * gpgve(DEF_VECT,1,2,igaus) + gpadv(DEF_VECT,2,igaus) * gpgve(DEF_VECT,2,2,igaus) &
                        + gpadv(DEF_VECT,3,igaus) * gpgve(DEF_VECT,3,2,igaus) )
                   elrbu(DEF_VECT,3,inode) = elrbu(DEF_VECT,3,inode) &
                        + fact1(DEF_VECT) * ( gpadv(DEF_VECT,1,igaus) * gpgve(DEF_VECT,1,3,igaus) + gpadv(DEF_VECT,2,igaus) * gpgve(DEF_VECT,2,3,igaus) &
                        + gpadv(DEF_VECT,3,igaus) * gpgve(DEF_VECT,3,3,igaus) )
                end do
             end do
          end if
       end if
       !
       ! Projection: ( L*(v) , P )
       !
       if( kfl_stabi_nsi == 1 ) then

          do igaus = 1, pgaus
             do inode = 1,pnode
                fact0(DEF_VECT) = gpsha(DEF_VECT,inode,igaus) * gpvol(DEF_VECT,igaus)
                do idime = 1,ndime
                   elrbu(DEF_VECT,idime,inode) = elrbu(DEF_VECT,idime,inode) - gprhs_sgs(DEF_VECT,idime,igaus) * fact0(DEF_VECT)
                end do
             end do
          end do

       end if
       !
       ! Immersed boundary
       !
       if( kfl_immer_nsi /= 0 ) then

          if( kfl_immer_nsi == 1 ) then
             xib = 2.0_rp
          else
             xib = 1.0_rp
          end if
          
          delta = 0.0_rp
          do idime = 1,ndime
             delta(idime,idime) = 1.0_rp
          end do
          
          do igaus = 1,pgaus             
             do inode = 1,pnode
                do kdime = 1,ndime
                   !
                   ! eps^{ij} = 1/2 ( dv_i/dx_j * delta_ik + dv_i/dx_j * delta_jk ) 
                   !               
                   eps_test(DEF_VECT,:,:) = 0.0_rp
                   do jdime = 1,ndime
                      do idime = 1,ndime
                         eps_test(DEF_VECT,idime,jdime) = eps_test(DEF_VECT,idime,jdime) &
                              + 0.5_rp * ( gpcar(DEF_VECT,jdime,inode,igaus) * delta(idime,kdime)    &
                              +            gpcar(DEF_VECT,idime,inode,igaus) * delta(jdime,kdime)    )
                      end do
                   end do
                   do jdime = 1,ndime
                      do idime = 1,ndime
                         elrbu(DEF_VECT,kdime,inode) = elrbu(DEF_VECT,kdime,inode)   &
                              - xib * gpvol(DEF_VECT,igaus) * gplag(DEF_VECT,idime,jdime,igaus) &
                              * eps_test(DEF_VECT,idime,jdime) * gpmix(DEF_VECT,igaus)
                      end do
                   end do
                end do
             end do
          end do
 
       end if
       
#ifdef OPENACC
    end do
#endif

    !--------------------------------------------------------------------
    !
    ! Pressure bubble
    !
    !--------------------------------------------------------------------

    if( maxval(pbubl) == 1 ) then
       !
       ! Initialization
       !
       elauq = 0.0_rp
       elapq = 0.0_rp
       elaqu = 0.0_rp
       elaqp = 0.0_rp
       elaqq = 0.0_rp
       elrbq = 0.0_rp
       !
       ! Test functions
       !
       do igaus = 1,pgaus
          p2sca_bub(DEF_VECT,igaus) = gptt2(DEF_VECT,igaus) * gpvol(DEF_VECT,igaus) * gpsha_bub(DEF_VECT,igaus)
          do idime = 1,ndime
             p2vec_bub(DEF_VECT,idime,igaus) = gpsp1(DEF_VECT,igaus) * gpvol(DEF_VECT,igaus) * gpcar_bub(DEF_VECT,idime,igaus)
          end do
       end do
       !
       ! Momentum equations
       !
       if( kfl_press_nsi == 1 ) then
          do igaus = 1,pgaus
             fact2(DEF_VECT) = gpvol(DEF_VECT,igaus) * gpsha_bub(DEF_VECT,igaus)
             do inode = 1,pnode
                fact1(DEF_VECT) = -gpvol(DEF_VECT,igaus) * gpsha(DEF_VECT,inode,igaus) + p1vec(DEF_VECT,inode,igaus)
                idof1 = (inode-1)*ndime
                do idime = 1,ndime
                   idofv = idof1 + idime
                   elauq(DEF_VECT,idofv,1) = elauq(DEF_VECT,idofv,1) - fact2(DEF_VECT) * gpcar(DEF_VECT,idime,inode,igaus) &
                        &                                            + fact1(DEF_VECT) * gpcar_bub(DEF_VECT,idime,igaus)
                end do
             end do
          end do
       else
          do igaus = 1,pgaus
             do inode = 1,pnode
                idof1 = (inode-1)*ndime
                do idime = 1,ndime
                   idofv = idof1 + idime
                   elauq(DEF_VECT,idofv,1) = elauq(DEF_VECT,idofv,1) + p1vec(DEF_VECT,inode,igaus) * gpcar_bub(DEF_VECT,idime,igaus)
                end do
             end do
          end do
       end if
       if( kfl_p1ve2_nsi /= 0 ) then
          do igaus = 1,pgaus
             do inode = 1,pnode
                idof1 = (inode-1)*ndime
                do idime = 1,ndime
                   idofv = idof1 + idime
                   do jdime = 1,ndime
                      elauq(DEF_VECT,idofv,1) = elauq(DEF_VECT,idofv,1) + &
                           gpcar_bub(DEF_VECT,jdime,igaus) * p1ve2(DEF_VECT,idime,jdime,inode,igaus)
                   end do
                end do
             end do
          end do
       end if
       !
       ! Bubble equation
       !
       do igaus = 1,pgaus
          do inode = 1,pnode
             do idime = 1,ndime
                idof1 = (inode-1)*ndime+idime
                elaqu(DEF_VECT,1,idof1) = elaqu(DEF_VECT,1,idof1) + rcont(DEF_VECT,idime,inode,igaus) * p2sca_bub(DEF_VECT,igaus) &
                     &                                            + rmomu(DEF_VECT,inode,igaus)       * p2vec_bub(DEF_VECT,idime,igaus)
             end do
          end do
       end do
       if( kfl_rmom2_nsi /= 0 ) then
          do igaus = 1,pgaus
             do jnode = 1,pnode
                do jdime = 1,ndime
                   jdofv = (jnode-1)*ndime+jdime
                   do kdime = 1,ndime
                      elaqu(DEF_VECT,1,jdofv) = elaqu(DEF_VECT,1,jdofv) + rmom2(DEF_VECT,kdime,jdime,jnode,igaus) * p2vec_bub(DEF_VECT,kdime,igaus)
                   end do
                end do
             end do
          end do
       end if
       !
       ! Also pressure equation...
       !
       do igaus = 1,pgaus
          do inode = 1,pnode
             do idime = 1,ndime
                elapq(DEF_VECT,inode,1) = elapq(DEF_VECT,inode,1) + p2vec(DEF_VECT,idime,inode,igaus) * gpcar_bub(DEF_VECT,idime,igaus)
                elaqp(DEF_VECT,1,inode) = elaqp(DEF_VECT,1,inode) + p2vec_bub(DEF_VECT,idime,igaus)   * gpcar(DEF_VECT,idime,inode,igaus)
             end do
          end do
          do idime = 1,ndime
             elaqq(DEF_VECT,1,1)             = elaqq(DEF_VECT,1,1)             + p2vec_bub(DEF_VECT,idime,igaus) * gpcar_bub(DEF_VECT,idime,igaus)
             gprhs_sgs(DEF_VECT,idime,igaus) = gprhs_sgs(DEF_VECT,idime,igaus) + taupr(DEF_VECT,igaus) * gpvep(DEF_VECT,idime,igaus)
             gprhh(DEF_VECT,idime,igaus)     = gprhs(DEF_VECT,idime,igaus)     + gprhs_sgs(DEF_VECT,idime,igaus)
             elrbq(DEF_VECT,1)               = elrbq(DEF_VECT,1)               + p2vec_bub(DEF_VECT,idime,igaus) * gprhh(DEF_VECT,idime,igaus)
          end do
       end do
       do igaus = 1,pgaus
          elrbq(DEF_VECT,1) = elrbq(DEF_VECT,1) + gprhc(DEF_VECT,igaus) * p2sca_bub(DEF_VECT,igaus) ! ( rhs , q_bub)
       end do
       !
       ! Penalization
       !
       do igaus = 1,pgaus
          elaqq(DEF_VECT,1,1) = elaqq(DEF_VECT,1,1) + penal_nsi * gpvol(DEF_VECT,igaus) * gpsha_bub(DEF_VECT,igaus)
          elrbq(DEF_VECT,1)   = elrbq(DEF_VECT,1)   + penal_nsi * gpvol(DEF_VECT,igaus) * gpsha_bub(DEF_VECT,igaus) * elbub(DEF_VECT)
       end do

    end if

  end subroutine nsi_element_assembly_asgs_oss

  !-----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @brief   ASGS and OSS
  !> @details Assembly of Navier Stokes equations using the split OSS
  !>          Variational Multiscale Model.
  !>
  !-----------------------------------------------------------------------

  subroutine nsi_element_assembly_split_oss(&
       pnode,pgaus,gpden,gpvis,gpvis_nsw,gppor,gpsp1,gpsp2,gpvol,   &
       gpsha,gpcar,gpadv,gpvep,gpgrp,gprhs,gprhc,gpvel,   &
       gpgve,gpsgs,elvel,elpre,elbub,elauu,elaup,elapp,   &
       elapu,elrbu,elrbp,dtinv_loc,dtsgs,pbubl,           &
       gpsha_bub,gpcar_bub,elauq,elapq,elaqu,elaqp,elaqq, &
       elrbq,densi,elavv,elporrbu)

    integer(ip), intent(in)    :: pnode,pgaus
    real(rp),    intent(in)    :: gpden(VECTOR_SIZE,pgaus)
    real(rp),    intent(in)    :: gpvis(VECTOR_SIZE,pgaus)
    real(rp),    intent(in)    :: gpvis_nsw(VECTOR_SIZE,pgaus)
    real(rp),    intent(in)    :: gppor(VECTOR_SIZE,pgaus)
    real(rp),    intent(in)    :: gpsp1(VECTOR_SIZE,pgaus)
    real(rp),    intent(in)    :: gpsp2(VECTOR_SIZE,pgaus)
    real(rp),    intent(in)    :: gpvol(VECTOR_SIZE,pgaus)
    real(rp),    intent(in)    :: gpsha(VECTOR_SIZE,pnode,pgaus)
    real(rp),    intent(in)    :: gpcar(VECTOR_SIZE,ndime,mnode,pgaus)
    real(rp),    intent(in)    :: gpadv(VECTOR_SIZE,ndime,pgaus)
    real(rp),    intent(inout) :: gpvep(VECTOR_SIZE,ndime,pgaus)
    real(rp),    intent(inout) :: gpgrp(VECTOR_SIZE,ndime,pgaus)
    real(rp),    intent(inout) :: gprhs(VECTOR_SIZE,ndime,pgaus)
    real(rp),    intent(inout) :: gprhc(VECTOR_SIZE,pgaus)
    real(rp),    intent(in)    :: gpvel(VECTOR_SIZE,ndime,pgaus,*)
    real(rp),    intent(in)    :: gpgve(VECTOR_SIZE,ndime,ndime,pgaus)     !< Velocity gradient
    real(rp),    intent(in)    :: gpsgs(VECTOR_SIZE,ndime,pgaus,*)
    real(rp),    intent(in)    :: elvel(VECTOR_SIZE,ndime,pnode,*)
    real(rp),    intent(in)    :: elavv(VECTOR_SIZE,ndime,pnode)
    real(rp),    intent(in)    :: elpre(VECTOR_SIZE,pnode,*)
    real(rp),    intent(in)    :: elbub(VECTOR_SIZE)
    ! Matrices
    real(rp),    intent(out)   :: elauu(VECTOR_SIZE,pnode*ndime,pnode*ndime)
    real(rp),    intent(out)   :: elaup(VECTOR_SIZE,pnode*ndime,pnode)
    real(rp),    intent(out)   :: elapp(VECTOR_SIZE,pnode,pnode)
    real(rp),    intent(out)   :: elapu(VECTOR_SIZE,pnode,pnode*ndime)
    real(rp),    intent(out)   :: elrbu(VECTOR_SIZE,ndime,pnode)
    real(rp),    intent(out)   :: elporrbu(VECTOR_SIZE,ndime,pnode)
    real(rp),    intent(out)   :: elrbp(VECTOR_SIZE,pnode)
    ! Others
    real(rp),    intent(in)    :: dtinv_loc(VECTOR_SIZE)
    real(rp),    intent(in)    :: dtsgs(VECTOR_SIZE)
    integer(ip), intent(in)    :: pbubl(VECTOR_SIZE)
    real(rp),    intent(in)    :: gpsha_bub(VECTOR_SIZE,pgaus)
    real(rp),    intent(in)    :: gpcar_bub(VECTOR_SIZE,ndime,pgaus)
    real(rp),    intent(in)    :: densi(VECTOR_SIZE,pgaus,nbdfp_nsi)
    ! Enrichement Element matrices
    real(rp),    intent(out)   :: elauq(VECTOR_SIZE,pnode*ndime,1)
    real(rp),    intent(out)   :: elapq(VECTOR_SIZE,pnode,1)
    real(rp),    intent(out)   :: elaqu(VECTOR_SIZE,1,pnode*ndime)
    real(rp),    intent(out)   :: elaqp(VECTOR_SIZE,1,pnode)
    real(rp),    intent(out)   :: elaqq(VECTOR_SIZE,1,1)
    real(rp),    intent(out)   :: elrbq(VECTOR_SIZE,1)

    !converting bool to int
    integer(ip)                :: NSI_FRACTIONAL_STEP_int

    if (NSI_FRACTIONAL_STEP) then
        NSI_FRACTIONAL_STEP_int=1
    else
        NSI_FRACTIONAL_STEP_int=0
    end if

#ifndef BOAST
  call nsi_element_assembly_split_oss_default(&
#else
#ifdef VECTOR_SIZE
  call nsi_element_assembly_split_oss_boast(&    ! para boast la voy a estropear  -- veremso que pasa
#else
#error "BOAST: VECTOR_SIZE is necessary"
#endif
#endif
       int(VECTOR_SIZE,ip), kfl_lumped, ndime,&
       mnode, ntens, kfl_stabi_nsi,&
       fvins_nsi, kfl_regim_nsi,&
       kfl_press_nsi,&
       kfl_linea_nsi, pabdf_nsi, nbdfp_nsi,&
       kfl_sgsti_nsi, kfl_nota1_nsi, kfl_limit_nsi,&
       penal_nsi, kfl_convection_type_nsi, &
       NSI_GALERKIN, NSI_ALGEBRAIC_SPLIT_OSS, &
       NSI_FRACTIONAL_STEP_int, &
       NSI_CONVECTION_CONSERVATIVE, NSI_CONVECTION_SKEW, &
       NSI_CONVECTION_EMAC,kfl_noslw_ker, &
       pnode,pgaus,gpden,gpvis,gpvis_nsw,gppor,gpsp1,gpsp2,gpvol,   &
       gpsha,gpcar,gpadv,gpvep,gpgrp,gprhs,gprhc,gpvel,   &
       gpgve,gpsgs,elvel,elpre,elbub,elauu,elaup,elapp,   &
       elapu,elrbu,elrbp,dtinv_loc,dtsgs,pbubl,           &
       gpsha_bub,gpcar_bub,elauq,elapq,elaqu,elaqp,elaqq, &
       elrbq,densi,elavv,elporrbu)       ! ojo voy a joder lo boast!!!!

  end subroutine nsi_element_assembly_split_oss

  !-----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @brief   ASGS and OSS
  !> @details Assembly of Navier Stokes equations using ASGS or OSS
  !>          Variational Multiscale Model
  !>
  !>          BEFORE YACINE INTERVENTION
  !>
  !-----------------------------------------------------------------------

  subroutine nsi_element_assembly_asgs_oss_old(         &
       pnode,pgaus,pevat,gpden,gpvis,gppor,gpgvi,gpsp1, &
       gptt1,gpsp2_v,gptt2,gpvol,gpsha,gpcar,gplap,gphes, &
       gpadv,gpvep,gprhs,gprhc,rmomu,rcont,p1vec,p2vec, &
       p2sca,wgrgr,wgrvi,elauu,elaup,elapp,elapu,elrbu, &
       elrbp,rmom2,p1ve2,gpst1,gpgve,gprh2,gprhs_sgs,   &
       elvel,ellum,dtinv_loc,pbubl,gpsha_bub,gpcar_bub, &
       gppre)

    integer(ip), intent(in)    :: pnode
    integer(ip), intent(in)    :: pgaus
    integer(ip), intent(in)    :: pevat
    real(rp),    intent(in)    :: gpden(VECTOR_SIZE,pgaus)
    real(rp),    intent(in)    :: gpvis(VECTOR_SIZE,pgaus)
    real(rp),    intent(in)    :: gppor(VECTOR_SIZE,pgaus)
    real(rp),    intent(in)    :: gpgvi(VECTOR_SIZE,ndime,pgaus)
    real(rp),    intent(inout) :: gpsp1(VECTOR_SIZE,pgaus)
    real(rp),    intent(in)    :: gptt1(VECTOR_SIZE,pgaus)
    real(rp),    intent(inout) :: gpsp2_v(VECTOR_SIZE,pgaus)
    real(rp),    intent(in)    :: gptt2(VECTOR_SIZE,pgaus)
    real(rp),    intent(in)    :: gpvol(VECTOR_SIZE,pgaus)
    real(rp),    intent(in)    :: gpsha(VECTOR_SIZE,pnode,pgaus)
    real(rp),    intent(in)    :: gpcar(VECTOR_SIZE,ndime,mnode,pgaus)
    real(rp),    intent(in)    :: gpadv(VECTOR_SIZE,ndime,pgaus)
    real(rp),    intent(inout) :: gpvep(VECTOR_SIZE,ndime,pgaus)
    real(rp),    intent(in)    :: gplap(VECTOR_SIZE,pnode,pgaus)
    real(rp),    intent(in)    :: gphes(VECTOR_SIZE,ntens,mnode,pgaus)
    real(rp),    intent(in)    :: gprhs(VECTOR_SIZE,ndime,pgaus)
    real(rp),    intent(inout) :: gprhs_sgs(VECTOR_SIZE,ndime,pgaus)
    real(rp),    intent(in)    :: gprhc(VECTOR_SIZE,pgaus)
    real(rp),    intent(in)    :: gprh2(VECTOR_SIZE,pgaus)
    real(rp),    intent(in)    :: rmomu(VECTOR_SIZE,pnode,pgaus)
    real(rp),    intent(in)    :: rcont(VECTOR_SIZE,ndime,pnode,pgaus)
    real(rp),    intent(out)   :: p1vec(VECTOR_SIZE,pnode,pgaus)
    real(rp),    intent(out)   :: p2vec(VECTOR_SIZE,ndime,pnode,pgaus)
    real(rp),    intent(out)   :: p2sca(VECTOR_SIZE,pnode,pgaus)
    real(rp),    intent(out)   :: wgrgr(VECTOR_SIZE,pnode,pnode,pgaus)
    real(rp),    intent(out)   :: wgrvi(VECTOR_SIZE,pnode,pgaus)
    ! Element matrices
    real(rp),    intent(out)   :: elauu(VECTOR_SIZE,pnode*ndime,pnode*ndime)
    real(rp),    intent(out)   :: elaup(VECTOR_SIZE,pnode*ndime,pnode)
    real(rp),    intent(out)   :: elapp(VECTOR_SIZE,pnode,pnode)
    real(rp),    intent(out)   :: elapu(VECTOR_SIZE,pnode,pnode*ndime)
    real(rp),    intent(out)   :: elrbu(VECTOR_SIZE,ndime,pnode)
    real(rp),    intent(out)   :: elrbp(VECTOR_SIZE,pnode)
    real(rp),    intent(in)    :: rmom2(VECTOR_SIZE,ndime,ndime,pnode,pgaus)
    real(rp),    intent(out)   :: p1ve2(VECTOR_SIZE,ndime,ndime,pnode,pgaus)
    real(rp),    intent(in)    :: gpst1(VECTOR_SIZE,pgaus)
    real(rp),    intent(in)    :: gpgve(VECTOR_SIZE,ndime,ndime,pgaus)
    real(rp),    intent(in)    :: elvel(VECTOR_SIZE,ndime,pnode,*)
    real(rp),    intent(out)   :: ellum(VECTOR_SIZE,pnode)
    real(rp),    intent(in)    :: dtinv_loc(VECTOR_SIZE)
    integer(ip), intent(in)    :: pbubl(VECTOR_SIZE)
    real(rp),    intent(in)    :: gpsha_bub(VECTOR_SIZE,pgaus)                    ! Ne
    real(rp),    intent(in)    :: gpcar_bub(VECTOR_SIZE,ndime,pgaus)              ! dNe/dxi
    real(rp),    intent(in)    :: gppre(VECTOR_SIZE,pgaus)
    ! Enrichement Element matrices
    real(rp)                   :: elauq(VECTOR_SIZE,pnode*ndime,1)
    real(rp)                   :: elapq(VECTOR_SIZE,pnode,1)
    real(rp)                   :: elaqu(VECTOR_SIZE,1,pnode*ndime)
    real(rp)                   :: elaqp(VECTOR_SIZE,1,pnode)
    real(rp)                   :: elaqq(VECTOR_SIZE,1,1)
    real(rp)                   :: elrbq(VECTOR_SIZE,1)

    integer(ip)                :: inode,jnode,idofv,jdime
    integer(ip)                :: idof1,idof2,idof3
#ifdef OPENACC
    integer(ip)                :: ivect
#endif
    integer(ip)                :: jdof1,jdof2,jdof3
    integer(ip)                :: igaus,idime,jdofv,kdime
    real(rp)                   :: fact1(VECTOR_SIZE),fact2(VECTOR_SIZE)
    real(rp)                   :: fact3(VECTOR_SIZE),fact4(VECTOR_SIZE)
    real(rp)                   :: fact5(VECTOR_SIZE),fact6(VECTOR_SIZE)
    real(rp)                   :: fact7(VECTOR_SIZE),fact0(VECTOR_SIZE)
    real(rp)                   :: ugraN(VECTOR_SIZE)
    real(rp)                   :: gramugraN(VECTOR_SIZE)
    real(rp)                   :: factx(VECTOR_SIZE)
    real(rp)                   :: facty(VECTOR_SIZE)
    real(rp)                   :: factz(VECTOR_SIZE)
    real(rp)                   :: gpveo(VECTOR_SIZE,ndime)
    real(rp)                   :: facx1(VECTOR_SIZE)
    real(rp)                   :: facy1(VECTOR_SIZE)
    real(rp)                   :: facz1(VECTOR_SIZE)
    real(rp)                   :: penal(VECTOR_SIZE)
    real(rp)                   :: taupr(VECTOR_SIZE,pgaus)
    real(rp)                   :: gprhh(VECTOR_SIZE,ndime,pgaus)
    real(rp)                   :: p2sca_bub(VECTOR_SIZE,pgaus)
    real(rp)                   :: p2vec_bub(VECTOR_SIZE,ndime,pgaus)
    real(rp)                   :: xvisc,xvis2,one_rp

#ifdef OPENACC
#define DEF_VECT ivect
#else
#define DEF_VECT 1:VECTOR_SIZE
#endif

    !----------------------------------------------------------------------
    !
    ! Initialization
    !
    !----------------------------------------------------------------------

    elrbp = 0.0_rp
    elrbu = 0.0_rp
    elauu = 0.0_rp
    elaup = 0.0_rp
    elapu = 0.0_rp
    elapp = 0.0_rp

    xvis2 = 0.0_rp
    if( fvins_nsi > 0.9_rp ) then
       xvisc = 1.0_rp
       if( fvins_nsi > 1.9_rp ) xvis2 = 2.0_rp / 3.0_rp
    else
       xvisc = 0.0_rp
    end if

    !----------------------------------------------------------------------
    !
    ! Test functions
    !
    !----------------------------------------------------------------------
    !
    ! P1VEC =  v * [ (tau1'/tau1) - tau1' * sig ] + tau1' * rho * (uc.grad)v
    ! P2SCA =  ( tau2^{-1} * tau2' ) * q
    ! P2VEC =  tau1' * grad(q)  ( and * rho if in Low-Mach, but not yet)
    ! WGRVI =  grad(Ni) . grad(mu) + mu * Delta(Ni)
    ! WGRGR =  grad(Ni) . grad(Nj)
    !
    if( kfl_regim_nsi == 3 ) then
       taupr = 1.0_rp
    else
       taupr = 1.0_rp / max(gpst1,zeror)
    end if

    one_rp = 1.0_rp
    if( kfl_nota1_nsi == 1 ) one_rp = 0.0_rp  ! do not stabilize convection, reaction, coriolis

    do igaus = 1,pgaus

       fact1(DEF_VECT) = gpsp1(DEF_VECT,igaus) * one_rp * gpden(DEF_VECT,igaus) * gpvol(DEF_VECT,igaus)
       fact2(DEF_VECT) = (gptt1(DEF_VECT,igaus)-gpsp1(DEF_VECT,igaus) * one_rp * gppor(DEF_VECT,igaus)) * gpvol(DEF_VECT,igaus)
       fact3(DEF_VECT) = gpsp1(DEF_VECT,igaus) * gpvol(DEF_VECT,igaus)
       fact4(DEF_VECT) = gptt2(DEF_VECT,igaus) * gpvol(DEF_VECT,igaus)

       do inode = 1,pnode
          p1vec(DEF_VECT,inode,igaus) = 0.0_rp
          ugraN(DEF_VECT)             = 0.0_rp
          gramugraN(DEF_VECT)         = 0.0_rp

          do idime = 1,ndime
             ugraN(DEF_VECT)                   = ugraN(DEF_VECT)     + gpadv(DEF_VECT,idime,igaus) * gpcar(DEF_VECT,idime,inode,igaus)
             gramugraN(DEF_VECT)               = gramugraN(DEF_VECT) + gpgvi(DEF_VECT,idime,igaus) * gpcar(DEF_VECT,idime,inode,igaus)
             p2vec(DEF_VECT,idime,inode,igaus) = fact3(DEF_VECT) * gpcar(DEF_VECT,idime,inode,igaus)
          end do

          p1vec(DEF_VECT,inode,igaus)          = fact2(DEF_VECT) * gpsha(DEF_VECT,inode,igaus) + fact1(DEF_VECT) * ugraN(DEF_VECT)
          p2sca(DEF_VECT,inode,igaus)          = fact4(DEF_VECT) * gpsha(DEF_VECT,inode,igaus)

          wgrvi(DEF_VECT,inode,igaus)          = gramugraN(DEF_VECT) + gpvis(DEF_VECT,igaus) * gplap(DEF_VECT,inode,igaus)

          do jnode = 1,pnode
             wgrgr(DEF_VECT,inode,jnode,igaus) = 0.0_rp
             do idime = 1,ndime
                wgrgr(DEF_VECT,inode,jnode,igaus) = wgrgr(DEF_VECT,inode,jnode,igaus) + gpcar(DEF_VECT,idime,inode,igaus) * gpcar(DEF_VECT,idime,jnode,igaus)
             end do
          end do
       end do

    end do
    !
    ! P1VE2: Off-diagonal test function
    !
    if( kfl_p1ve2_nsi /= 0 ) then
       p1ve2(DEF_VECT,:,:,:,:) = 0.0_rp
       !
       ! tau1'*2*rho*(w x v)
       ! x-equation: v = (v,0,0) => w x v = (    0 , wz*v , -wy*v)
       ! y-equation: v = (0,v,0) => w x v = (-wz*v ,    0 ,  wx*v)
       ! z-equation: v = (0,0,v) => w x v = ( wy*v ,-wx*v ,     0)
       !
       if( ndime == 2 ) then
          do igaus = 1,pgaus
             fact3(DEF_VECT) = 2.0_rp * gpden(DEF_VECT,igaus) * fvela_nsi(3)
             factz(DEF_VECT) = gpsp1(DEF_VECT,igaus) * fact3(DEF_VECT) * gpvol(DEF_VECT,igaus)
             do inode = 1,pnode
                fact1(DEF_VECT)                 =   factz(DEF_VECT) * gpsha(DEF_VECT,inode,igaus)
                p1ve2(DEF_VECT,1,2,inode,igaus) =   fact1(DEF_VECT)
                p1ve2(DEF_VECT,2,1,inode,igaus) = - fact1(DEF_VECT)
             end do
          end do
       else
          do igaus = 1,pgaus
             fact3(DEF_VECT) = 2.0_rp * gpden(DEF_VECT,igaus) * gpsp1(DEF_VECT,igaus) * gpvol(DEF_VECT,igaus)
             factx(DEF_VECT) = fact3(DEF_VECT)  * fvela_nsi(1)
             facty(DEF_VECT) = fact3(DEF_VECT)  * fvela_nsi(2)
             factz(DEF_VECT) = fact3(DEF_VECT)  * fvela_nsi(3)
             do inode = 1,pnode
                p1ve2(DEF_VECT,1,2,inode,igaus) =  factz(DEF_VECT) * gpsha(DEF_VECT,inode,igaus)  ! x-equation
                p1ve2(DEF_VECT,1,3,inode,igaus) = -facty(DEF_VECT) * gpsha(DEF_VECT,inode,igaus)  ! x-equation
                p1ve2(DEF_VECT,2,1,inode,igaus) = -factz(DEF_VECT) * gpsha(DEF_VECT,inode,igaus)  ! y-equation
                p1ve2(DEF_VECT,2,3,inode,igaus) =  factx(DEF_VECT) * gpsha(DEF_VECT,inode,igaus)  ! y-equation
                p1ve2(DEF_VECT,3,1,inode,igaus) =  facty(DEF_VECT) * gpsha(DEF_VECT,inode,igaus)  ! z-equation
                p1ve2(DEF_VECT,3,2,inode,igaus) = -factx(DEF_VECT) * gpsha(DEF_VECT,inode,igaus)  ! z-equation
             end do
          end do
       end if
    end if

    !----------------------------------------------------------------------
    !
    ! Auu
    !
    !----------------------------------------------------------------------
    !
    !  Laplacian  form: A=0, B=0, eps(u) = 1/2 grad(u)
    !  Divergence form: A=1, B=0, eps(u) = 1/2 ( grad(u) + grad(u)^t )
    !  Complete   form: A=1, B=1, eps(u) = 1/2 ( grad(u) + grad(u)^t ) - 1/3 div(u) I
    !
    !    ( div(u) , tau2' div(v) )             (1)   <=   tau2' * duj/dxj * dvi/dxi
    !  + ( rmomu(u) , p1vec(v) )               (2)
    !  + ( mu dui/dxj , dvi/dxj )              (3)   <=   ( 2 mu eps(u) , grad(v) ) = mu * ( dui/dxj + duj/dxi - 2/3 div(u) delta_ij ) dvi/dxj
    !  + A * ( mu duj/dxi , dvi/dxj )          (4)        divergence
    !  - 2/3 * B * ( mu (div u) , dvi/dxi )    (5)        complete
    !  + ( mu d^2ui/dxj^2 , vi )               (6)   <= - ( div[-2 mu eps(u)] , v ) = - d/dxj ( -2*mu* ( dui/dxj + duj/dxi - 2/3*mu div(u) delta_ij ) * vi
    !  + ( dmu/dxj dui/dxj , vi )              (7)        laplacian
    !  + A * ( mu  d^2uj/dxidxj , vi )         (8)        divergence
    !  + A * ( dmu/dxj duj/dxi , vi )          (9)        divergence
    !  - 2/3 * B * ( dmu/dxi (div u) , vi )   (10)        complete
    !  - 2/3 * B * ( mu d(div u)/dxi , vi )   (11)        complete
    !
    !  1. The terms ( div(u) , tau2' div(v) ) and - 2/3 * B * ( mu div u, dvi/dxj ) are assembled together as
    !     ( ( tau2' - 2/3 B mu ) div u, dvi/dxi ). Therefore, term (5) is assembled always although it is zero
    !  2. Terms (3), (6) and (7) are assembled just below as they are the same for all momentum equations
    !  3. Terms (4), (8), (9), (10) and (11) are assmbled later on
    !
    !  (3)        x:     mu * ( dux/dx * dv/dx  +  dux/dy * dv/dy  +  dux/dz * dv/dz )
    !             y:     mu * ( duy/dx * dv/dx  +  duy/dy * dv/dy  +  duy/dz * dv/dz )
    !             z:     mu * ( duz/dx * dv/dx  +  duz/dy * dv/dy  +  duz/dz * dv/dz )
    !
    !  (4)        x: A * mu * ( dux/dx * dv/dx  +  duy/dx * dv/dy  +  duz/dx * dv/dz )
    !             y: A * mu * ( dux/dy * dv/dx  +  duy/dy * dv/dy  +  duz/dy * dv/dz )
    !             z: A * mu * ( dux/dz * dv/dx  +  duy/dz * dv/dy  +  duz/dz * dv/dz )
    !
    !  (6)        x:  mu * ( d^2ux/dx^2  +  d^2ux/dy^2  +  d^2ux/dz^2 ) v
    !             y:  mu * ( d^2uy/dx^2  +  d^2uy/dy^2  +  d^2uy/dz^2 ) v
    !             z:  mu * ( d^2uz/dx^2  +  d^2uz/dy^2  +  d^2uz/dz^2 ) v
    !
    !  (7)        x: ( dmu/dx * dux/dx  +  dmu/dy * dux/dy  +  dmu/dz * dux/dz ) v
    !             y: ( dmu/dx * duy/dx  +  dmu/dy * duy/dy  +  dmu/dz * duy/dz ) v
    !             y: ( dmu/dx * duz/dx  +  dmu/dy * duz/dy  +  dmu/dz * duz/dz ) v
    !
    !  (8)        x: A * mu * ( d^2ux/dxdx  +  d^2uy/dxdy  +  d^2uz/dxdz ) v
    !             y: A * mu * ( d^2ux/dxdy  +  d^2uy/dydy  +  d^2uz/dydz ) v
    !             z: A * mu * ( d^2ux/dxdz  +  d^2uy/dzdy  +  d^2uz/dzdz ) v
    !
    !  (9)        x: A * ( dmu/dx * dux/dx  +  dmu/dy * duy/dx  +  dmu/dz * duz/dx ) v
    !             y: A * ( dmu/dx * dux/dy  +  dmu/dy * duy/dy  +  dmu/dz * duz/dy ) v
    !             z: A * ( dmu/dx * dux/dz  +  dmu/dy * duy/dz  +  dmu/dz * duz/dz ) v
    !
    !  (10)       x: -2/3 B * dmu/dx ( dux/dx  +  duy/dy  +  duz/dz ) * v
    !             y: -2/3 B * dmu/dy ( dux/dx  +  duy/dy  +  duz/dz ) * v
    !             z: -2/3 B * dmu/dz ( dux/dx  +  duy/dy  +  duz/dz ) * v
    !
    !  (11)       x: -2/3 B * mu * ( d^2ux/dxdx + d^2uy/dydx + d^2uz/dzdx ) v
    !             y: -2/3 B * mu * ( d^2ux/dxdy + d^2uy/dydy + d^2uz/dzdy ) v
    !             z: -2/3 B * mu * ( d^2ux/dxdz + d^2uy/dydz + d^2uz/dzdz ) v
    !
    !  (1) + (5)  x: ( tau2' - 2/3 B mu ) ( dux/dx  +  duy/dy  +  duz/dz ) * dv/dx
    !             y: ( tau2' - 2/3 B mu ) ( dux/dx  +  duy/dy  +  duz/dz ) * dv/dy
    !             z: ( tau2' - 2/3 B mu ) ( dux/dx  +  duy/dy  +  duz/dz ) * dv/dz
    !

    do igaus = 1,pgaus

       fact0(DEF_VECT) = ( gpsp2_v(DEF_VECT,igaus) - xvis2 * gpvis(DEF_VECT,igaus) ) * gpvol(DEF_VECT,igaus)  ! (1) + (5)
       fact6(DEF_VECT) = gpvis(DEF_VECT,igaus) * gpvol(DEF_VECT,igaus)

       do inode = 1,pnode
          do idime = 1,ndime
             idof1 = (inode-1)*ndime + idime
             fact1(DEF_VECT) = fact0(DEF_VECT) * gpcar(DEF_VECT,idime,inode,igaus)
             fact5(DEF_VECT) = gpsha(DEF_VECT,inode,igaus)*gpvol(DEF_VECT,igaus)
             do jnode = 1,pnode
                do jdime = 1,ndime
                   jdof1 = (jnode-1)*ndime + jdime
                   elauu(DEF_VECT,idof1,jdof1) = elauu(DEF_VECT,idof1,jdof1) &
                        + fact1(DEF_VECT) * rcont(DEF_VECT,jdime,jnode,igaus)                               ! Auu_xixj
                end do
                fact4(DEF_VECT) = p1vec(DEF_VECT,inode,igaus)*rmomu(DEF_VECT,jnode,igaus) &                 ! (2):       ( rmomu(u) , p1vec(v) )
                     &  + fact5(DEF_VECT) * wgrvi(DEF_VECT,jnode,igaus)&                                    ! (3) + (6): ( mu*d^2ui/dxk^2 + grad(mu).grad(u), v )
                     &  + fact6(DEF_VECT) * wgrgr(DEF_VECT,inode,jnode,igaus)                               ! (7):       ( mu*dui/dxk , dv/dxk )
                jdof1 =  (jnode-1)*ndime + idime
                elauu(DEF_VECT,idof1,jdof1) = elauu(DEF_VECT,idof1,jdof1) + fact4(DEF_VECT)
             end do
          end do
       end do
    end do

!!$  if( ndime == 2 ) then
!!$
!!$     do igaus = 1,pgaus
!!$
!!$        fact0 = ( gpsp2_v(igaus) - xvis2 * gpvis(igaus) ) * gpvol(igaus)                        ! (1) + (5)
!!$        fact6 = gpvis(igaus) * gpvol(igaus)
!!$
!!$        do inode = 1,pnode
!!$           idof1 = 2*inode-1
!!$           idof2 = idof1+1
!!$           fact1 = fact0 * gpcar(1,inode,igaus)
!!$           fact2 = fact0 * gpcar(2,inode,igaus)
!!$           fact5 = gpsha(inode,igaus) * gpvol(igaus)
!!$
!!$           do jnode = 1,pnode
!!$              jdof1              =   2*jnode-1
!!$              jdof2              =   2*jnode
!!$              fact4              =   p1vec(inode,igaus) * rmomu(jnode,igaus) &                ! (2):       ( rmomu(u) , p1vec(v) )
!!$                   &               + fact5 * wgrvi(jnode,igaus)              &                ! (3) + (6): ( mu d^2ui/dxj^2 + dmu/dxj dui/dxj, vi )
!!$                   &               + fact6 * wgrgr(inode,jnode,igaus)                         ! (7):       ( mu dui/dxj , dvi/dxj )
!!$
!!$              elauu(idof1,jdof1) = elauu(idof1,jdof1) + fact1 * rcont(1,jnode,igaus) + fact4  ! Auu_xx
!!$              elauu(idof2,jdof1) = elauu(idof2,jdof1) + fact2 * rcont(1,jnode,igaus)          ! Auu_yx
!!$              elauu(idof1,jdof2) = elauu(idof1,jdof2) + fact1 * rcont(2,jnode,igaus)          ! Auu_xy
!!$              elauu(idof2,jdof2) = elauu(idof2,jdof2) + fact2 * rcont(2,jnode,igaus) + fact4  ! Auu_yy
!!$
!!$           end do
!!$
!!$        end do
!!$     end do
!!$
!!$  else
!!$
!!$     do igaus = 1,pgaus
!!$
!!$        fact0 = ( gpsp2_v(igaus) - xvis2 * gpvis(igaus) ) * gpvol(igaus)                        ! (1) + (5)
!!$        fact6 = gpvis(igaus) * gpvol(igaus)
!!$
!!$        do inode = 1,pnode
!!$           idof1 = 3*inode-2
!!$           idof2 = idof1+1
!!$           idof3 = idof2+1
!!$
!!$           fact1 = fact0*gpcar(1,inode,igaus)
!!$           fact2 = fact0*gpcar(2,inode,igaus)
!!$           fact3 = fact0*gpcar(3,inode,igaus)
!!$           fact5 = gpsha(inode,igaus)*gpvol(igaus)
!!$
!!$           do jnode = 1,pnode
!!$              jdof1              = 3*jnode-2
!!$              jdof2              = 3*jnode-1
!!$              jdof3              = 3*jnode
!!$
!!$              fact4              =   p1vec(inode,igaus)*rmomu(jnode,igaus) &                  ! (2):       ( rmomu(u) , p1vec(v) )
!!$                   &               + fact5*wgrvi(jnode,igaus)&                                ! (3) + (6): ( mu*d^2ui/dxk^2 + grad(mu).grad(u), v )
!!$                   &               + fact6*wgrgr(inode,jnode,igaus)                           ! (7):       ( mu*dui/dxk , dv/dxk )
!!$
!!$              elauu(idof1,jdof1) = elauu(idof1,jdof1) + fact1 * rcont(1,jnode,igaus) + fact4  ! Auu_xx
!!$              elauu(idof2,jdof1) = elauu(idof2,jdof1) + fact2 * rcont(1,jnode,igaus)          ! Auu_yx
!!$              elauu(idof3,jdof1) = elauu(idof3,jdof1) + fact3 * rcont(1,jnode,igaus)          ! Auu_zx
!!$              elauu(idof1,jdof2) = elauu(idof1,jdof2) + fact1 * rcont(2,jnode,igaus)          ! Auu_xy
!!$              elauu(idof2,jdof2) = elauu(idof2,jdof2) + fact2 * rcont(2,jnode,igaus) + fact4  ! Auu_yy
!!$              elauu(idof3,jdof2) = elauu(idof3,jdof2) + fact3 * rcont(2,jnode,igaus)          ! Auu_zy
!!$              elauu(idof1,jdof3) = elauu(idof1,jdof3) + fact1 * rcont(3,jnode,igaus)          ! Auu_xz
!!$              elauu(idof2,jdof3) = elauu(idof2,jdof3) + fact2 * rcont(3,jnode,igaus)          ! Auu_yz
!!$              elauu(idof3,jdof3) = elauu(idof3,jdof3) + fact3 * rcont(3,jnode,igaus) + fact4  ! Auu_zz
!!$           end do
!!$
!!$        end do
!!$     end do
!!$
!!$  end if

    if( fvins_nsi > 0.9_rp ) then

       if( ndime == 2 ) then

          do igaus = 1,pgaus

             fact4(DEF_VECT) = gpgvi(DEF_VECT,1,igaus) * gpvol(DEF_VECT,igaus)
             fact5(DEF_VECT) = gpgvi(DEF_VECT,2,igaus) * gpvol(DEF_VECT,igaus)
             fact7(DEF_VECT) = gpvis(DEF_VECT,igaus)   * gpvol(DEF_VECT,igaus)
             do inode = 1,pnode
                idof1           = 2*inode-1
                idof2           = 2*inode
                fact0(DEF_VECT) = fact7(DEF_VECT) * gpsha(DEF_VECT,inode,igaus) * (1.0_rp-xvis2)
                fact1(DEF_VECT) = fact4(DEF_VECT) * gpsha(DEF_VECT,inode,igaus)
                fact2(DEF_VECT) = fact5(DEF_VECT) * gpsha(DEF_VECT,inode,igaus)
                facx1(DEF_VECT) = fact1(DEF_VECT) * xvis2
                facy1(DEF_VECT) = fact2(DEF_VECT) * xvis2

                do jnode = 1,pnode
                   jdof1              = 2*jnode-1
                   jdof2              = 2*jnode
                   !
                   !                                         (4)  ( mu duj/dxi , dvi/dxj )
                   !
                   elauu(DEF_VECT,idof1,jdof1) = elauu(DEF_VECT,idof1,jdof1) + fact7(DEF_VECT) * gpcar(DEF_VECT,1,inode,igaus) * gpcar(DEF_VECT,1,jnode,igaus)                    ! Auu_xx
                   elauu(DEF_VECT,idof1,jdof2) = elauu(DEF_VECT,idof1,jdof2) + fact7(DEF_VECT) * gpcar(DEF_VECT,2,inode,igaus) * gpcar(DEF_VECT,1,jnode,igaus)                    ! Auu_xy
                   elauu(DEF_VECT,idof2,jdof1) = elauu(DEF_VECT,idof2,jdof1) + fact7(DEF_VECT) * gpcar(DEF_VECT,1,inode,igaus) * gpcar(DEF_VECT,2,jnode,igaus)                    ! Auu_yx
                   elauu(DEF_VECT,idof2,jdof2) = elauu(DEF_VECT,idof2,jdof2) + fact7(DEF_VECT) * gpcar(DEF_VECT,2,inode,igaus) * gpcar(DEF_VECT,2,jnode,igaus)                    ! Auu_yy
                   !
                   !                                         (10) - 2/3 * B * ( dmu/dxi (div u) , vi )
                   !2
                   elauu(DEF_VECT,idof1,jdof1) = elauu(DEF_VECT,idof1,jdof1) - facx1(DEF_VECT) * gpcar(DEF_VECT,1,jnode,igaus)                                                    ! Auu_xx
                   elauu(DEF_VECT,idof1,jdof2) = elauu(DEF_VECT,idof1,jdof2) - facx1(DEF_VECT) * gpcar(DEF_VECT,2,jnode,igaus)                                                    ! Auu_xy
                   elauu(DEF_VECT,idof2,jdof1) = elauu(DEF_VECT,idof2,jdof1) - facy1(DEF_VECT) * gpcar(DEF_VECT,1,jnode,igaus)                                                    ! Auu_yx
                   elauu(DEF_VECT,idof2,jdof2) = elauu(DEF_VECT,idof2,jdof2) - facy1(DEF_VECT) * gpcar(DEF_VECT,2,jnode,igaus)                                                    ! Auu_yy
                   !
                   !                                         (9)  ( dmu/dxj duj/dxi , vi ) + (8)  ( mu d^2uj/dxidxj , vi )
                   !                                               + (11)   - 2/3 * B * ( mu d(div u)/dxi , vi )
                   !
                   elauu(DEF_VECT,idof1,jdof1) = elauu(DEF_VECT,idof1,jdof1) + fact1(DEF_VECT) * gpcar(DEF_VECT,1,jnode,igaus) + fact0(DEF_VECT) * gphes(DEF_VECT,1,jnode,igaus)  ! Auu_xx
                   elauu(DEF_VECT,idof1,jdof2) = elauu(DEF_VECT,idof1,jdof2) + fact2(DEF_VECT) * gpcar(DEF_VECT,1,jnode,igaus) + fact0(DEF_VECT) * gphes(DEF_VECT,3,jnode,igaus)  ! Auu_xy
                   elauu(DEF_VECT,idof2,jdof1) = elauu(DEF_VECT,idof2,jdof1) + fact1(DEF_VECT) * gpcar(DEF_VECT,2,jnode,igaus) + fact0(DEF_VECT) * gphes(DEF_VECT,3,jnode,igaus)  ! Auu_yx
                   elauu(DEF_VECT,idof2,jdof2) = elauu(DEF_VECT,idof2,jdof2) + fact2(DEF_VECT) * gpcar(DEF_VECT,2,jnode,igaus) + fact0(DEF_VECT) * gphes(DEF_VECT,2,jnode,igaus)  ! Auu_yy

                end do

             end do
          end do

       else

          do igaus = 1,pgaus

             fact4(DEF_VECT) = gpgvi(DEF_VECT,1,igaus) * gpvol(DEF_VECT,igaus)
             fact5(DEF_VECT) = gpgvi(DEF_VECT,2,igaus) * gpvol(DEF_VECT,igaus)
             fact6(DEF_VECT) = gpgvi(DEF_VECT,3,igaus) * gpvol(DEF_VECT,igaus)
             fact7(DEF_VECT) = gpvis(DEF_VECT,igaus)   * gpvol(DEF_VECT,igaus)

             do inode = 1,pnode
                idof1           = 3*inode-2
                idof2           = 3*inode-1
                idof3           = 3*inode
                fact0(DEF_VECT) = fact7(DEF_VECT) * gpsha(DEF_VECT,inode,igaus) * (1.0_rp-xvis2)
                fact1(DEF_VECT) = fact4(DEF_VECT) * gpsha(DEF_VECT,inode,igaus)
                fact2(DEF_VECT) = fact5(DEF_VECT) * gpsha(DEF_VECT,inode,igaus)
                fact3(DEF_VECT) = fact6(DEF_VECT) * gpsha(DEF_VECT,inode,igaus)
                facx1(DEF_VECT) = fact1(DEF_VECT) * xvis2
                facy1(DEF_VECT) = fact2(DEF_VECT) * xvis2
                facz1(DEF_VECT) = fact3(DEF_VECT) * xvis2

                do jnode = 1,pnode
                   jdof1              = 3*jnode-2
                   jdof2              = 3*jnode-1
                   jdof3              = 3*jnode
                   !
                   !                                         (4)  ( mu duj/dxi , dvi/dxj )
                   !
                   elauu(DEF_VECT,idof1,jdof1) = elauu(DEF_VECT,idof1,jdof1) + fact7(DEF_VECT) * gpcar(DEF_VECT,1,inode,igaus) * gpcar(DEF_VECT,1,jnode,igaus)                    ! Auu_xx
                   elauu(DEF_VECT,idof1,jdof2) = elauu(DEF_VECT,idof1,jdof2) + fact7(DEF_VECT) * gpcar(DEF_VECT,2,inode,igaus) * gpcar(DEF_VECT,1,jnode,igaus)                    ! Auu_xy
                   elauu(DEF_VECT,idof1,jdof3) = elauu(DEF_VECT,idof1,jdof3) + fact7(DEF_VECT) * gpcar(DEF_VECT,3,inode,igaus) * gpcar(DEF_VECT,1,jnode,igaus)                    ! Auu_xz

                   elauu(DEF_VECT,idof2,jdof1) = elauu(DEF_VECT,idof2,jdof1) + fact7(DEF_VECT) * gpcar(DEF_VECT,1,inode,igaus) * gpcar(DEF_VECT,2,jnode,igaus)                    ! Auu_yx
                   elauu(DEF_VECT,idof2,jdof2) = elauu(DEF_VECT,idof2,jdof2) + fact7(DEF_VECT) * gpcar(DEF_VECT,2,inode,igaus) * gpcar(DEF_VECT,2,jnode,igaus)                    ! Auu_yy
                   elauu(DEF_VECT,idof2,jdof3) = elauu(DEF_VECT,idof2,jdof3) + fact7(DEF_VECT) * gpcar(DEF_VECT,3,inode,igaus) * gpcar(DEF_VECT,2,jnode,igaus)                    ! Auu_yz

                   elauu(DEF_VECT,idof3,jdof1) = elauu(DEF_VECT,idof3,jdof1) + fact7(DEF_VECT) * gpcar(DEF_VECT,1,inode,igaus) * gpcar(DEF_VECT,3,jnode,igaus)                    ! Auu_zx
                   elauu(DEF_VECT,idof3,jdof2) = elauu(DEF_VECT,idof3,jdof2) + fact7(DEF_VECT) * gpcar(DEF_VECT,2,inode,igaus) * gpcar(DEF_VECT,3,jnode,igaus)                    ! Auu_zy
                   elauu(DEF_VECT,idof3,jdof3) = elauu(DEF_VECT,idof3,jdof3) + fact7(DEF_VECT) * gpcar(DEF_VECT,3,inode,igaus) * gpcar(DEF_VECT,3,jnode,igaus)                    ! Auu_zz
                   !
                   !                                         (10) - 2/3 * B * ( dmu/dxi (div u) , vi )

                   elauu(DEF_VECT,idof1,jdof1) = elauu(DEF_VECT,idof1,jdof1) - facx1(DEF_VECT) * gpcar(DEF_VECT,1,jnode,igaus)                                                    ! Auu_xx
                   elauu(DEF_VECT,idof1,jdof2) = elauu(DEF_VECT,idof1,jdof2) - facx1(DEF_VECT) * gpcar(DEF_VECT,2,jnode,igaus)                                                    ! Auu_xy
                   elauu(DEF_VECT,idof1,jdof3) = elauu(DEF_VECT,idof1,jdof3) - facx1(DEF_VECT) * gpcar(DEF_VECT,3,jnode,igaus)                                                    ! Auu_xz

                   elauu(DEF_VECT,idof2,jdof1) = elauu(DEF_VECT,idof2,jdof1) - facy1(DEF_VECT) * gpcar(DEF_VECT,1,jnode,igaus)                                                    ! Auu_yx
                   elauu(DEF_VECT,idof2,jdof2) = elauu(DEF_VECT,idof2,jdof2) - facy1(DEF_VECT) * gpcar(DEF_VECT,2,jnode,igaus)                                                    ! Auu_yy
                   elauu(DEF_VECT,idof2,jdof3) = elauu(DEF_VECT,idof2,jdof3) - facy1(DEF_VECT) * gpcar(DEF_VECT,3,jnode,igaus)                                                    ! Auu_yz

                   elauu(DEF_VECT,idof3,jdof1) = elauu(DEF_VECT,idof3,jdof1) - facz1(DEF_VECT) * gpcar(DEF_VECT,1,jnode,igaus)                                                    ! Auu_zx
                   elauu(DEF_VECT,idof3,jdof2) = elauu(DEF_VECT,idof3,jdof2) - facz1(DEF_VECT) * gpcar(DEF_VECT,2,jnode,igaus)                                                    ! Auu_zy
                   elauu(DEF_VECT,idof3,jdof3) = elauu(DEF_VECT,idof3,jdof3) - facz1(DEF_VECT) * gpcar(DEF_VECT,3,jnode,igaus)                                                    ! Auu_zz
                   !
                   !                                         (9)  ( dmu/dxj duj/dxi , vi ) + (8)  ( mu d^2uj/dxidxj , vi )
                   !                                                + (11)   - 2/3 * B * ( mu d(div u)/dxi , vi )
                   !
                   elauu(DEF_VECT,idof1,jdof1) = elauu(DEF_VECT,idof1,jdof1) + fact1(DEF_VECT) * gpcar(DEF_VECT,1,jnode,igaus) + fact0(DEF_VECT) * gphes(DEF_VECT,1,jnode,igaus)  ! Auu_xx
                   elauu(DEF_VECT,idof1,jdof2) = elauu(DEF_VECT,idof1,jdof2) + fact2(DEF_VECT) * gpcar(DEF_VECT,1,jnode,igaus) + fact0(DEF_VECT) * gphes(DEF_VECT,4,jnode,igaus)  ! Auu_xy
                   elauu(DEF_VECT,idof1,jdof3) = elauu(DEF_VECT,idof1,jdof3) + fact3(DEF_VECT) * gpcar(DEF_VECT,1,jnode,igaus) + fact0(DEF_VECT) * gphes(DEF_VECT,5,jnode,igaus)  ! Auu_xz

                   elauu(DEF_VECT,idof2,jdof1) = elauu(DEF_VECT,idof2,jdof1) + fact1(DEF_VECT) * gpcar(DEF_VECT,2,jnode,igaus) + fact0(DEF_VECT) * gphes(DEF_VECT,4,jnode,igaus)  ! Auu_yx
                   elauu(DEF_VECT,idof2,jdof2) = elauu(DEF_VECT,idof2,jdof2) + fact2(DEF_VECT) * gpcar(DEF_VECT,2,jnode,igaus) + fact0(DEF_VECT) * gphes(DEF_VECT,2,jnode,igaus)  ! Auu_yy
                   elauu(DEF_VECT,idof2,jdof3) = elauu(DEF_VECT,idof2,jdof3) + fact3(DEF_VECT) * gpcar(DEF_VECT,2,jnode,igaus) + fact0(DEF_VECT) * gphes(DEF_VECT,6,jnode,igaus)  ! Auu_yz

                   elauu(DEF_VECT,idof3,jdof1) = elauu(DEF_VECT,idof3,jdof1) + fact1(DEF_VECT) * gpcar(DEF_VECT,3,jnode,igaus) + fact0(DEF_VECT) * gphes(DEF_VECT,5,jnode,igaus)  ! Auu_zx
                   elauu(DEF_VECT,idof3,jdof2) = elauu(DEF_VECT,idof3,jdof2) + fact2(DEF_VECT) * gpcar(DEF_VECT,3,jnode,igaus) + fact0(DEF_VECT) * gphes(DEF_VECT,6,jnode,igaus)  ! Auu_zy
                   elauu(DEF_VECT,idof3,jdof3) = elauu(DEF_VECT,idof3,jdof3) + fact3(DEF_VECT) * gpcar(DEF_VECT,3,jnode,igaus) + fact0(DEF_VECT) * gphes(DEF_VECT,3,jnode,igaus)  ! Auu_zz

                end do

             end do
          end do

       end if

    end if

    

    if( fcons_nsi > 0.1_rp ) then

     do igaus = 1,pgaus
        fact0(DEF_VECT) = fcons_nsi * gpden(DEF_VECT,igaus) * gpvol(DEF_VECT,igaus)
        do inode = 1,pnode
           do jnode = 1,pnode
              fact1(DEF_VECT) = 0.0_rp
              fact2(DEF_VECT) = 0.0_rp
              do idime = 1,ndime
                 fact1(DEF_VECT) = fact1(DEF_VECT) + gpadv(DEF_VECT,idime,igaus) * gpcar(DEF_VECT,idime,inode,igaus)
                 fact2(DEF_VECT) = fact2(DEF_VECT) + gpadv(DEF_VECT,idime,igaus) * gpcar(DEF_VECT,idime,jnode,igaus)
              end do
              fact1(DEF_VECT) = fact1(DEF_VECT)  * gpsha(DEF_VECT,jnode,igaus) * fact0(DEF_VECT)
              fact2(DEF_VECT) = fact2(DEF_VECT)  * gpsha(DEF_VECT,inode,igaus) * fact0(DEF_VECT)
              do idime = 1,ndime
                 idof1 = (inode-1)*ndime+idime
                 jdof1 = (jnode-1)*ndime+idime
                 elauu(DEF_VECT,idof1,jdof1) = elauu(DEF_VECT,idof1,jdof1) - fact1(DEF_VECT) - fact2(DEF_VECT)
              end do
           end do
        end do
     end do

!!$       if( ndime == 2 ) then
!!$
!!$          do igaus = 1,pgaus
!!$
!!$             fact0(DEF_VECT) = fcons_nsi * gpden(DEF_VECT,igaus) * gpvol(DEF_VECT,igaus)
!!$
!!$             do inode = 1,pnode
!!$                idof1 = 2*inode-1
!!$                idof2 = 2*inode
!!$
!!$                do jnode = 1,pnode
!!$                   jdof1           = 2*jnode-1
!!$                   jdof2           = 2*jnode
!!$
!!$                   fact1(DEF_VECT) = (gpadv(DEF_VECT,1,igaus) * gpcar(DEF_VECT,1,inode,igaus) + gpadv(DEF_VECT,2,igaus) * gpcar(DEF_VECT,2,inode,igaus)) * gpsha(DEF_VECT,jnode,igaus) * fact0(DEF_VECT)
!!$                   fact2(DEF_VECT) = (gpadv(DEF_VECT,1,igaus) * gpcar(DEF_VECT,1,jnode,igaus) + gpadv(DEF_VECT,2,igaus) * gpcar(DEF_VECT,2,jnode,igaus)) * gpsha(DEF_VECT,inode,igaus) * fact0(DEF_VECT)
!!$
!!$                   elauu(DEF_VECT,idof1,jdof1) = elauu(DEF_VECT,idof1,jdof1) - fact1(DEF_VECT) - fact2(DEF_VECT)                     ! Auu_xx
!!$                   elauu(DEF_VECT,idof2,jdof2) = elauu(DEF_VECT,idof2,jdof2) - fact1(DEF_VECT) - fact2(DEF_VECT)                     ! Auu_yy
!!$
!!$                end do
!!$
!!$             end do
!!$
!!$          end do
!!$
!!$       else
!!$
!!$          do igaus = 1,pgaus
!!$
!!$             fact0(DEF_VECT) = fcons_nsi * gpden(DEF_VECT,igaus) * gpvol(DEF_VECT,igaus)
!!$
!!$             do inode = 1,pnode
!!$                idof1 = 3*inode-2
!!$                idof2 = 3*inode-1
!!$                idof3 = 3*inode
!!$
!!$                do jnode = 1,pnode
!!$                   jdof1           = 3*jnode-2
!!$                   jdof2           = 3*jnode-1
!!$                   jdof3           = 3*jnode
!!$

!!$
!!$                   elauu(DEF_VECT,idof1,jdof1) = elauu(DEF_VECT,idof1,jdof1) - fact1(DEF_VECT) - fact2(DEF_VECT)                      ! Auu_xx
!!$                   elauu(DEF_VECT,idof2,jdof2) = elauu(DEF_VECT,idof2,jdof2) - fact1(DEF_VECT) - fact2(DEF_VECT)                     ! Auu_yy
!!$                   elauu(DEF_VECT,idof3,jdof3) = elauu(DEF_VECT,idof3,jdof3) - fact1(DEF_VECT) - fact2(DEF_VECT)                      ! Auu_zz
!!$                end do
!!$
!!$             end do
!!$
!!$          end do
!!$
!!$       end if

    end if

    !----------------------------------------------------------------------
    !
    ! Auu: Off-diagonal terms
    !
    !----------------------------------------------------------------------

    if( kfl_rmom2_nsi /= 0 ) then
       !
       ! p1vec * rmom2
       !
       do igaus = 1,pgaus
          do inode = 1,pnode
             do idime = 1,ndime
                idof1 = (inode-1)*ndime + idime
                do jnode = 1,pnode
                   do jdime = 1,ndime
                      jdof1 = (jnode-1)*ndime + jdime
                      elauu(DEF_VECT,idof1,jdof1) = elauu(DEF_VECT,idof1,jdof1) &
                           + p1vec(DEF_VECT,inode,igaus) * rmom2(DEF_VECT,idime,jdime,jnode,igaus)
                   end do
                end do
             end do
          end do
       end do

       if( kfl_p1ve2_nsi /= 0 ) then
          !
          ! p1ve2 * rmomu
          ! p1ve2 * rmom2
          !
          do igaus = 1,pgaus
             do inode = 1,pnode
                do idime = 1,ndime
                   idofv = (inode-1)*ndime + idime
                   do jnode = 1,pnode
                      do jdime = 1,ndime
                         jdofv = (jnode-1)*ndime + jdime
                         elauu(DEF_VECT,idofv,jdofv) = elauu(DEF_VECT,idofv,jdofv)    &
                              + p1ve2(DEF_VECT,idime,jdime,inode,igaus)     &
                              * rmomu(DEF_VECT,jnode,igaus)
                         do kdime = 1,ndime
                            elauu(DEF_VECT,idofv,jdofv) = elauu(DEF_VECT,idofv,jdofv) &
                                 + p1ve2(DEF_VECT,idime,kdime,inode,igaus)  &
                                 * rmom2(DEF_VECT,kdime,jdime,jnode,igaus)
                         end do
                      end do
                   end do
                end do
             end do
          end do
       end if

    end if
    !
    ! Bemol:
    ! - b * ( v , rho * (a.grad) u ) - b * ( u , rho * (a.grad) v ) - b ( rho*div(a) u , v )
    !
    if( abs(bemol_nsi) > 1.0e-9_rp ) then

       if( ndime == 2 ) then

          do igaus = 1,pgaus

             fact0(DEF_VECT) = bemol_nsi * gpden(DEF_VECT,igaus) * gpvol(DEF_VECT,igaus)
             fact5(DEF_VECT) = fact0(DEF_VECT) * ( gpgve(DEF_VECT,1,1,igaus) + gpgve(DEF_VECT,2,2,igaus) )

             do inode = 1,pnode
                idof1 = 2*inode-1
                idof2 = 2*inode

                do jnode = 1,pnode
                   jdof1 = 2*jnode-1
                   jdof2 = 2*jnode
                   !
                   ! fact2 = - b * ( v , rho * (a.grad) u )
                   ! fact3 = - b * ( u , rho * (a.grad) v )
                   ! fact4 = - b * ( rho * div(a) u , v )
                   !
                   fact2(DEF_VECT) = ( gpadv(DEF_VECT,1,igaus) * gpcar(DEF_VECT,1,inode,igaus) + gpadv(DEF_VECT,2,igaus) * gpcar(DEF_VECT,2,inode,igaus) ) * gpsha(DEF_VECT,jnode,igaus)
                   fact3(DEF_VECT) = ( gpadv(DEF_VECT,1,igaus) * gpcar(DEF_VECT,1,jnode,igaus) + gpadv(DEF_VECT,2,igaus) * gpcar(DEF_VECT,2,jnode,igaus) ) * gpsha(DEF_VECT,inode,igaus)
                   fact4(DEF_VECT) = fact5(DEF_VECT) * gpsha(DEF_VECT,inode,igaus) * gpsha(DEF_VECT,jnode,igaus)
                   fact2(DEF_VECT) = fact0(DEF_VECT) * fact2(DEF_VECT)
                   fact3(DEF_VECT) = fact0(DEF_VECT) * fact3(DEF_VECT)
                   fact6(DEF_VECT) = fact2(DEF_VECT) + fact3(DEF_VECT) + fact4(DEF_VECT)

                   elauu(DEF_VECT,idof1,jdof1) = elauu(DEF_VECT,idof1,jdof1) - fact6(DEF_VECT)
                   elauu(DEF_VECT,idof2,jdof2) = elauu(DEF_VECT,idof2,jdof2) - fact6(DEF_VECT)

                end do

             end do

          end do

       else

          do igaus = 1,pgaus

             fact0(DEF_VECT) = bemol_nsi * gpden(DEF_VECT,igaus) * gpvol(DEF_VECT,igaus)
             fact5(DEF_VECT) = fact0(DEF_VECT) * ( gpgve(DEF_VECT,1,1,igaus) + gpgve(DEF_VECT,2,2,igaus) + gpgve(DEF_VECT,3,3,igaus) )

             do inode = 1,pnode
                idof1 = 3*inode-2
                idof2 = 3*inode-1
                idof3 = 3*inode

                do jnode = 1,pnode
                   jdof1 = 3*jnode-2
                   jdof2 = 3*jnode-1
                   jdof3 = 3*jnode
                   !
                   ! fact2 = - b * ( v , rho * (a.grad) u )
                   ! fact3 = - b * ( u , rho * (a.grad) v )
                   ! fact4 = - b * ( rho * div(a) u , v )
                   !
                   fact2(DEF_VECT) = ( gpadv(DEF_VECT,1,igaus) * gpcar(DEF_VECT,1,inode,igaus) + gpadv(DEF_VECT,2,igaus) * gpcar(DEF_VECT,2,inode,igaus) + &
                        gpadv(DEF_VECT,3,igaus) * gpcar(DEF_VECT,3,inode,igaus)) * gpsha(DEF_VECT,jnode,igaus)
                   fact3(DEF_VECT) = ( gpadv(DEF_VECT,1,igaus) * gpcar(DEF_VECT,1,jnode,igaus) + gpadv(DEF_VECT,2,igaus) * gpcar(DEF_VECT,2,jnode,igaus) + & 
                        gpadv(DEF_VECT,3,igaus) * gpcar(DEF_VECT,3,jnode,igaus)) * gpsha(DEF_VECT,inode,igaus)
                   fact4(DEF_VECT) = fact5(DEF_VECT) * gpsha(DEF_VECT,inode,igaus) * gpsha(DEF_VECT,jnode,igaus)
                   fact2(DEF_VECT) = fact0(DEF_VECT) * fact2(DEF_VECT)
                   fact3(DEF_VECT) = fact0(DEF_VECT) * fact3(DEF_VECT)
                   fact6(DEF_VECT) = fact2(DEF_VECT) + fact3(DEF_VECT) + fact4(DEF_VECT)

                   elauu(DEF_VECT,idof1,jdof1) = elauu(DEF_VECT,idof1,jdof1) - fact6(DEF_VECT)
                   elauu(DEF_VECT,idof2,jdof2) = elauu(DEF_VECT,idof2,jdof2) - fact6(DEF_VECT)
                   elauu(DEF_VECT,idof3,jdof3) = elauu(DEF_VECT,idof3,jdof3) - fact6(DEF_VECT)

                end do

             end do

          end do

       end if

    end if
    !
    ! Newton-Raphson
    !
    if( kfl_linea_nsi == 2 ) then
       if( ndime == 2 ) then
          do igaus = 1,pgaus
             do inode = 1,pnode
                idof1           = 2*inode-1
                idof2           = 2*inode
                fact0(DEF_VECT) = gpvol(DEF_VECT,igaus) * gpden(DEF_VECT,igaus) * gpsha(DEF_VECT,inode,igaus)
                do jnode = 1,pnode
                   jdof1           = 2*jnode-1
                   jdof2           = 2*jnode
                   fact1(DEF_VECT) = fact0(DEF_VECT) * gpsha(DEF_VECT,jnode,igaus)
                   ! Auu_1i
                   elauu(DEF_VECT,idof1,jdof1) = elauu(DEF_VECT,idof1,jdof1) + fact1(DEF_VECT) * gpgve(DEF_VECT,1,1,igaus) ! rho * ux * dux^i-1/dx
                   elauu(DEF_VECT,idof1,jdof2) = elauu(DEF_VECT,idof1,jdof2) + fact1(DEF_VECT) * gpgve(DEF_VECT,2,1,igaus) ! rho * uy * dux^i-1/dy
                   ! Auu_2i
                   elauu(DEF_VECT,idof2,jdof1) = elauu(DEF_VECT,idof2,jdof1) + fact1(DEF_VECT) * gpgve(DEF_VECT,1,2,igaus) ! rho * ux * duy^i-1/dx
                   elauu(DEF_VECT,idof2,jdof2) = elauu(DEF_VECT,idof2,jdof2) + fact1(DEF_VECT) * gpgve(DEF_VECT,2,2,igaus) ! rho * uy * duy^i-1/dy
                end do
             end do
          end do
       else
          do igaus = 1,pgaus
             do inode = 1,pnode
                idof1           = 3*inode-2
                idof2           = 3*inode-1
                idof3           = 3*inode
                fact0(DEF_VECT) = gpvol(DEF_VECT,igaus) * gpden(DEF_VECT,igaus) * gpsha(DEF_VECT,inode,igaus)
                do jnode = 1,pnode
                   jdof1           = 3*jnode-2
                   jdof2           = 3*jnode-1
                   jdof3           = 3*jnode
                   fact1(DEF_VECT) = fact0(DEF_VECT) * gpsha(DEF_VECT,jnode,igaus)
                   ! Auu_1i
                   elauu(DEF_VECT,idof1,jdof1) = elauu(DEF_VECT,idof1,jdof1) + fact1(DEF_VECT) * gpgve(DEF_VECT,1,1,igaus) ! rho * ux * dux^i-1/dx
                   elauu(DEF_VECT,idof1,jdof2) = elauu(DEF_VECT,idof1,jdof2) + fact1(DEF_VECT) * gpgve(DEF_VECT,2,1,igaus) ! rho * uy * dux^i-1/dy
                   elauu(DEF_VECT,idof1,jdof3) = elauu(DEF_VECT,idof1,jdof3) + fact1(DEF_VECT) * gpgve(DEF_VECT,3,1,igaus) ! rho * uz * dux^i-1/dz
                   ! Auu_2i
                   elauu(DEF_VECT,idof2,jdof1) = elauu(DEF_VECT,idof2,jdof1) + fact1(DEF_VECT) * gpgve(DEF_VECT,1,2,igaus) ! rho * ux * duy^i-1/dx
                   elauu(DEF_VECT,idof2,jdof2) = elauu(DEF_VECT,idof2,jdof2) + fact1(DEF_VECT) * gpgve(DEF_VECT,2,2,igaus) ! rho * uy * duy^i-1/dy
                   elauu(DEF_VECT,idof2,jdof3) = elauu(DEF_VECT,idof2,jdof3) + fact1(DEF_VECT) * gpgve(DEF_VECT,3,2,igaus) ! rho * uz * duy^i-1/dz
                   ! Auu_3i
                   elauu(DEF_VECT,idof3,jdof1) = elauu(DEF_VECT,idof3,jdof1) + fact1(DEF_VECT) * gpgve(DEF_VECT,1,3,igaus) ! rho * ux * duz^i-1/dx
                   elauu(DEF_VECT,idof3,jdof2) = elauu(DEF_VECT,idof3,jdof2) + fact1(DEF_VECT) * gpgve(DEF_VECT,2,3,igaus) ! rho * uy * duz^i-1/dy
                   elauu(DEF_VECT,idof3,jdof3) = elauu(DEF_VECT,idof3,jdof3) + fact1(DEF_VECT) * gpgve(DEF_VECT,3,3,igaus) ! rho * uz * duz^i-1/dz

                end do
             end do
          end do
       end if
    end if
    !
    ! Lumped evolution matrix (only backward euler)
    !
    if( kfl_lumped == 1 ) then
       !
       ! Remove Galerkin term and add lumped term
       !
       if( ndime == 2 ) then
          call runend('PREGUNTAR A MATIAS QUE LO PROGRAME')
       else
          do igaus =1, pgaus
             gpveo(DEF_VECT,1:3) = 0.0_rp
             do inode = 1,pnode
                gpveo(DEF_VECT,1) = gpveo(DEF_VECT,1) + elvel(DEF_VECT,1,inode,2) * gpsha(DEF_VECT,inode,igaus)
                gpveo(DEF_VECT,2) = gpveo(DEF_VECT,2) + elvel(DEF_VECT,2,inode,2) * gpsha(DEF_VECT,inode,igaus)
                gpveo(DEF_VECT,3) = gpveo(DEF_VECT,3) + elvel(DEF_VECT,3,inode,2) * gpsha(DEF_VECT,inode,igaus)
             end do
             do inode = 1,pnode
                idof1                       = 3*inode-2
                idof2                       = 3*inode-1
                idof3                       = 3*inode
                fact0(DEF_VECT)             = gpvol(DEF_VECT,igaus) * gpden(DEF_VECT,igaus) * gpsha(DEF_VECT,inode,igaus) * dtinv_loc(DEF_VECT)
                elauu(DEF_VECT,idof1,idof1) = elauu(DEF_VECT,idof1,idof1) + fact0(DEF_VECT)
                elauu(DEF_VECT,idof2,idof2) = elauu(DEF_VECT,idof2,idof2) + fact0(DEF_VECT)
                elauu(DEF_VECT,idof3,idof3) = elauu(DEF_VECT,idof3,idof3) + fact0(DEF_VECT)
                do idime = 1,ndime
                   elrbu(DEF_VECT,idime,inode)   = elrbu(DEF_VECT,idime,inode)   - fact0(DEF_VECT) * gpveo(DEF_VECT,idime)
                   elrbu(DEF_VECT,idime,inode)   = elrbu(DEF_VECT,idime,inode)   + fact0(DEF_VECT) * elvel(DEF_VECT,idime,inode,2)
                end do
                do jnode = 1,pnode
                   jdof1                       = 3*jnode-2
                   jdof2                       = 3*jnode-1
                   jdof3                       = 3*jnode
                   elauu(DEF_VECT,idof1,jdof1) = elauu(DEF_VECT,idof1,jdof1) - fact0(DEF_VECT) * gpsha(DEF_VECT,jnode,igaus)
                   elauu(DEF_VECT,idof2,jdof2) = elauu(DEF_VECT,idof2,jdof2) - fact0(DEF_VECT) * gpsha(DEF_VECT,jnode,igaus)
                   elauu(DEF_VECT,idof3,jdof3) = elauu(DEF_VECT,idof3,jdof3) - fact0(DEF_VECT) * gpsha(DEF_VECT,jnode,igaus)
                end do
             end do
          end do
       end if

    else if( kfl_lumped == 2 ) then
       !
       ! No time term have been added up to now: add Galerkin term
       !
       do igaus = 1,pgaus
          fact0(DEF_VECT) = gpvol(DEF_VECT,igaus) * gpden(DEF_VECT,igaus) * dtinv_loc(DEF_VECT)
          do inode = 1, pnode
             fact1(DEF_VECT) = fact0(DEF_VECT) * gpsha(DEF_VECT,inode,igaus)
             do idime = 1,ndime
                idof1                       = (inode-1) * ndime + idime
                elauu(DEF_VECT,idof1,idof1) = elauu(DEF_VECT,idof1,idof1) + fact1(DEF_VECT)
                elrbu(DEF_VECT,idime,inode) = elrbu(DEF_VECT,idime,inode) + fact1(DEF_VECT) * elvel(DEF_VECT,idime,inode,2)
             end do
          end do
       end do

    end if
    !
    ! Dual time step preconditioner
    !
    if( kfl_duatss == 1 ) then
       ellum(DEF_VECT,:) = 0.0_rp
       do igaus = 1,pgaus
          do inode = 1,pnode
             fact0(DEF_VECT)       = gpvol(DEF_VECT,igaus) * gpden(DEF_VECT,igaus) * gpsha(DEF_VECT,inode,igaus) * dtinv_loc(DEF_VECT) * real(fact_duatss-1,rp)
             ellum(DEF_VECT,inode) = ellum(DEF_VECT,inode) + fact0(DEF_VECT)
          end do
       end do
    end if

    !----------------------------------------------------------------------
    !
    ! Aup
    !
    !----------------------------------------------------------------------
    !
    ! Pressure: - ( p , div(v)  ) + ( grad(p) , p1vec(v)-v ) ( if kfl_press_nsi = 1 )
    ! Pressure: + ( grad(p) , p1vec(v) )                     ( if kfl_press_nsi = 0 )
    !
    if( kfl_press_nsi == 1 ) then
       do igaus = 1,pgaus
          do jnode = 1,pnode
             fact2(DEF_VECT) = gpvol(DEF_VECT,igaus) * gpsha(DEF_VECT,jnode,igaus)
             do inode = 1,pnode
                fact1 = -gpvol(DEF_VECT,igaus) * gpsha(DEF_VECT,inode,igaus) + p1vec(DEF_VECT,inode,igaus)
                idof1 = (inode-1)*ndime
                do idime = 1,ndime
                   idofv = idof1 + idime
                   elaup(DEF_VECT,idofv,jnode) = elaup(DEF_VECT,idofv,jnode) - fact2(DEF_VECT) * gpcar(DEF_VECT,idime,inode,igaus) &
                        &                                  + fact1(DEF_VECT) * gpcar(DEF_VECT,idime,jnode,igaus)
                end do
             end do
          end do
       end do
    else
       do igaus = 1,pgaus
          do jnode = 1,pnode
             do inode = 1,pnode
                idof1 = (inode-1)*ndime
                do idime = 1,ndime
                   idofv = idof1 + idime
                   elaup(DEF_VECT,idofv,jnode) = elaup(DEF_VECT,idofv,jnode) + p1vec(DEF_VECT,inode,igaus) * gpcar(DEF_VECT,idime,jnode,igaus)
                end do
             end do
          end do
       end do
    end if
    if( kfl_p1ve2_nsi /= 0 ) then
       do igaus = 1,pgaus
          do jnode = 1,pnode
             do inode = 1,pnode
                idof1 = (inode-1)*ndime
                do idime = 1,ndime
                   idofv = idof1 + idime
                   do jdime = 1,ndime
                      elaup(DEF_VECT,idofv,jnode) = elaup(DEF_VECT,idofv,jnode) + &
                           gpcar(DEF_VECT,jdime,jnode,igaus) * p1ve2(DEF_VECT,idime,jdime,inode,igaus)
                   end do
                end do
             end do
          end do
       end do
    end if

    !----------------------------------------------------------------------
    !
    ! Apu
    !
    !----------------------------------------------------------------------
    !
    ! ( div(u) , (tau2^{-1}*tau2')*q )   --- or ( div (rho u), q) if Low-Mach
    ! + ( rho*(uc.grad)u + 2*rho*(w x u) + sig*u -div[2*mu*eps(u)], tau1' grad(q) )  (only laplacian form of div[2 mu eps(u)] )
    !
    if( kfl_regim_nsi == 3 .and. kfl_confi_nsi == 1 ) then ! Penalization of pressure, never imposse pressure
       do igaus =1, pgaus
          penal(DEF_VECT)  = 1.0e-4_rp * gpden(DEF_VECT,igaus) / gpvis(DEF_VECT,igaus)
          do inode = 1,pnode
             fact0(DEF_VECT)  = penal(DEF_VECT) * p2sca(DEF_VECT,inode,igaus)
             do jnode = inode+1,pnode
                fact1(DEF_VECT)             = fact0(DEF_VECT) * gpsha(DEF_VECT,jnode,igaus)
                elapp(DEF_VECT,inode,jnode) = elapp(DEF_VECT,inode,jnode) + fact1(DEF_VECT)
                elapp(DEF_VECT,jnode,inode) = elapp(DEF_VECT,jnode,inode) + fact1(DEF_VECT)
             end do
             elapp(DEF_VECT,inode,inode) = elapp(DEF_VECT,inode,inode) + fact0(DEF_VECT) * gpsha(DEF_VECT,inode,igaus)
             elrbp(DEF_VECT,inode)       = elrbp(DEF_VECT,inode)       + fact0(DEF_VECT) * gppre(DEF_VECT,igaus)
          end do
       end do
    end if

#ifdef matiaslma

    if( ndime == 2 ) then
       do igaus = 1,pgaus
          do inode = 1,pnode
             p2vec(DEF_VECT,1,inode,igaus) = p2vec(DEF_VECT,1,inode,igaus) * gpden(DEF_VECT,igaus)
             p2vec(DEF_VECT,2,inode,igaus) = p2vec(DEF_VECT,2,inode,igaus) * gpden(DEF_VECT,igaus)
             p2sca(DEF_VECT,inode,igaus)   = p2sca(DEF_VECT,inode,igaus)   * gpden(DEF_VECT,igaus)
          end do
          fact0(DEF_VECT) = gpden(DEF_VECT,igaus) * gpvol(DEF_VECT,igaus)
          do inode = 1,pnode
             idof1           = 2*inode-1
             idof2           = 2*inode
             fact1(DEF_VECT) = gpsha(DEF_VECT,inode,igaus) * fact0(DEF_VECT)
             do jnode = 1,pnode
                elapu(DEF_VECT,jnode,idof1) = elapu(DEF_VECT,jnode,idof1) - fact1(DEF_VECT)             * gpcar(DEF_VECT,1,jnode,igaus)  &
                     &                                                    + rmomu(DEF_VECT,inode,igaus) * p2vec(DEF_VECT,1,jnode,igaus)
                elapu(DEF_VECT,jnode,idof2) = elapu(DEF_VECT,jnode,idof2) - fact1(DEF_VECT)             * gpcar(DEF_VECT,2,jnode,igaus)  &
                     &                                                    + rmomu(DEF_VECT,inode,igaus) * p2vec(DEF_VECT,2,jnode,igaus)
             end do
          end do
       end do
    else !ndime =3
       do igaus = 1,pgaus
          do inode = 1,pnode
             p2vec(DEF_VECT,1,inode,igaus) = p2vec(DEF_VECT,1,inode,igaus) * gpden(DEF_VECT,igaus)
             p2vec(DEF_VECT,2,inode,igaus) = p2vec(DEF_VECT,2,inode,igaus) * gpden(DEF_VECT,igaus)
             p2vec(DEF_VECT,3,inode,igaus) = p2vec(DEF_VECT,3,inode,igaus) * gpden(DEF_VECT,igaus)
             p2sca(DEF_VECT,inode,igaus)   = p2sca(DEF_VECT,inode,igaus)   * gpden(DEF_VECT,igaus)
          end do
          fact0(DEF_VECT) = gpden(DEF_VECT,igaus)* gpvol(DEF_VECT,igaus)
          do inode = 1, pnode
             idof1           = 3*inode-2
             idof2           = 3*inode-1
             idof3           = 3*inode
             fact1(DEF_VECT) = gpsha(DEF_VECT,inode,igaus)* fact0(DEF_VECT)
             do jnode = 1,pnode
                elapu(DEF_VECT,jnode,idof1) = elapu(DEF_VECT,jnode,idof1) - fact1(DEF_VECT)              * gpcar(DEF_VECT,1,jnode,igaus) &
                     &                                                    + rmomu(DEF_VECT,inode,igaus)  * p2vec(DEF_VECT,1,jnode,igaus)
                elapu(DEF_VECT,jnode,idof2) = elapu(DEF_VECT,jnode,idof2) - fact1(DEF_VECT)              * gpcar(DEF_VECT,2,jnode,igaus) &
                     &                                                    + rmomu(DEF_VECT,inode,igaus)  * p2vec(DEF_VECT,2,jnode,igaus)
                elapu(DEF_VECT,jnode,idof3) = elapu(DEF_VECT,jnode,idof3) - fact1(DEF_VECT)              * gpcar(DEF_VECT,3,jnode,igaus) &
                     &                                                     + rmomu(DEF_VECT,inode,igaus) * p2vec(DEF_VECT,3,jnode,igaus)
             end do
          end do
       end do
    end if

#else

    do igaus = 1,pgaus
       do inode = 1,pnode
          do jnode = 1,pnode
             do idime = 1,ndime
                idof1 = (inode-1)*ndime+idime
                elapu(DEF_VECT,jnode,idof1) = elapu(DEF_VECT,jnode,idof1) + rcont(DEF_VECT,idime,inode,igaus) * p2sca(DEF_VECT,jnode,igaus) &
                     &                                                    + rmomu(DEF_VECT,inode,igaus)       * p2vec(DEF_VECT,idime,jnode,igaus)
             end do
          end do
       end do
    end do

#endif
    !
    !yor! optim [originally loop on line 985] reordered inner loops
    !
    if( kfl_rmom2_nsi /= 0 ) then
       do igaus = 1,pgaus
          do jnode = 1,pnode
             do jdime = 1,ndime
                jdofv = (jnode-1)*ndime+jdime
                do inode = 1,pnode
                   do kdime = 1,ndime
                      elapu(DEF_VECT,inode,jdofv) = elapu(DEF_VECT,inode,jdofv) + rmom2(DEF_VECT,kdime,jdime,jnode,igaus) * p2vec(DEF_VECT,kdime,inode,igaus)
                   end do
                end do
             end do
          end do
       end do
    end if

    !----------------------------------------------------------------------
    !
    ! App
    !
    !----------------------------------------------------------------------
    !
    ! Pressure: ( grad(p) , tau1' grad(q) )
    !
    do igaus = 1,pgaus
       do inode = 1,pnode
          do jnode = 1,pnode
             do idime = 1,ndime
                elapp(DEF_VECT,jnode,inode) = elapp(DEF_VECT,jnode,inode) + p2vec(DEF_VECT,idime,jnode,igaus) * gpcar(DEF_VECT,idime,inode,igaus)
             end do
          end do
       end do
    end do
    !
    ! Penalization
    !
    stop
    ! do igaus = 1,pgaus
    !    fact1(DEF_VECT) = penal_nsi * gpvol(DEF_VECT,igaus)
    !    do inode = 1,pnode
    !       elapp(DEF_VECT,inode,inode) = elapp(DEF_VECT,inode,inode) + fact1(DEF_VECT) * gpsha(DEF_VECT,inode,igaus)
    !       elrbp(DEF_VECT,inode)       = elrbp(DEF_VECT,inode)       + fact1(DEF_VECT) * gpsha(DEF_VECT,inode,igaus) * elpre(DEF_VECT,inode,1)
    !    end do
    ! end do

    !----------------------------------------------------------------------
    !
    ! bu and bp
    !
    !----------------------------------------------------------------------
    !
    ! bu = ( f , p1vec(v) )
    ! bp = ( f , tau1' grad(q) )
    !
    if( ndime == 2 ) then
       do igaus = 1,pgaus
          gprhs_sgs(DEF_VECT,1,igaus) = gprhs_sgs(DEF_VECT,1,igaus) + taupr(DEF_VECT,igaus) * gpvep(DEF_VECT,1,igaus)
          gprhs_sgs(DEF_VECT,2,igaus) = gprhs_sgs(DEF_VECT,2,igaus) + taupr(DEF_VECT,igaus) * gpvep(DEF_VECT,2,igaus)
          gprhh(DEF_VECT,1,igaus)     = gprhs(DEF_VECT,1,igaus)     + gprhs_sgs(DEF_VECT,1,igaus)
          gprhh(DEF_VECT,2,igaus)     = gprhs(DEF_VECT,2,igaus)     + gprhs_sgs(DEF_VECT,2,igaus)
          do inode = 1,pnode
             elrbu(DEF_VECT,1,inode)  = elrbu(DEF_VECT,1,inode)     + p1vec(DEF_VECT,inode,igaus)   * gprhh(DEF_VECT,1,igaus)
             elrbu(DEF_VECT,2,inode)  = elrbu(DEF_VECT,2,inode)     + p1vec(DEF_VECT,inode,igaus)   * gprhh(DEF_VECT,2,igaus)
             elrbp(DEF_VECT,inode)    = elrbp(DEF_VECT,inode)       + p2vec(DEF_VECT,1,inode,igaus) * gprhh(DEF_VECT,1,igaus) &
                  &                                                 + p2vec(DEF_VECT,2,inode,igaus) * gprhh(DEF_VECT,2,igaus)
          end do
       end do
    else
       do igaus = 1,pgaus
          gprhs_sgs(DEF_VECT,1,igaus) = gprhs_sgs(DEF_VECT,1,igaus) + taupr(DEF_VECT,igaus) * gpvep(DEF_VECT,1,igaus)
          gprhs_sgs(DEF_VECT,2,igaus) = gprhs_sgs(DEF_VECT,2,igaus) + taupr(DEF_VECT,igaus) * gpvep(DEF_VECT,2,igaus)
          gprhs_sgs(DEF_VECT,3,igaus) = gprhs_sgs(DEF_VECT,3,igaus) + taupr(DEF_VECT,igaus) * gpvep(DEF_VECT,3,igaus)
          gprhh(DEF_VECT,1,igaus)     = gprhs(DEF_VECT,1,igaus)     + gprhs_sgs(DEF_VECT,1,igaus)
          gprhh(DEF_VECT,2,igaus)     = gprhs(DEF_VECT,2,igaus)     + gprhs_sgs(DEF_VECT,2,igaus)
          gprhh(DEF_VECT,3,igaus)     = gprhs(DEF_VECT,3,igaus)     + gprhs_sgs(DEF_VECT,3,igaus)
          do inode = 1,pnode
             elrbu(DEF_VECT,1,inode)  = elrbu(DEF_VECT,1,inode)     + p1vec(DEF_VECT,inode,igaus)   * gprhh(DEF_VECT,1,igaus)
             elrbu(DEF_VECT,2,inode)  = elrbu(DEF_VECT,2,inode)     + p1vec(DEF_VECT,inode,igaus)   * gprhh(DEF_VECT,2,igaus)
             elrbu(DEF_VECT,3,inode)  = elrbu(DEF_VECT,3,inode)     + p1vec(DEF_VECT,inode,igaus)   * gprhh(DEF_VECT,3,igaus)
             elrbp(DEF_VECT,inode)    = elrbp(DEF_VECT,inode)       + p2vec(DEF_VECT,1,inode,igaus) * gprhh(DEF_VECT,1,igaus) &
                  &                                                 + p2vec(DEF_VECT,2,inode,igaus) * gprhh(DEF_VECT,2,igaus) &
                  &                                                 + p2vec(DEF_VECT,3,inode,igaus) * gprhh(DEF_VECT,3,igaus)
          end do
       end do
    end if
    if( kfl_p1ve2_nsi /= 0 ) then
       do igaus = 1,pgaus
          do inode = 1,pnode
             do idime = 1,ndime
                do kdime = 1,ndime
                   elrbu(DEF_VECT,idime,inode) = elrbu(DEF_VECT,idime,inode) &
                        + p1ve2(DEF_VECT,idime,kdime,inode,igaus) * gprhh(DEF_VECT,kdime,igaus)
                end do
             end do
          end do
       end do
    end if
    !
    ! low Mach regime: we add to the right hand side the residue of the continuity eq
    !
    do igaus = 1,pgaus
       fact1(DEF_VECT) = gpsp2_v(DEF_VECT,igaus) * gprh2(DEF_VECT,igaus) * gpvol(DEF_VECT,igaus)
       do inode = 1,pnode
          do idime = 1,ndime
             elrbu(DEF_VECT,idime,inode) = elrbu(DEF_VECT,idime,inode) + fact1(DEF_VECT) * gpcar(DEF_VECT,idime,inode,igaus)
          end do
          elrbp(DEF_VECT,inode) = elrbp(DEF_VECT,inode) + gprhc(DEF_VECT,igaus) * p2sca(DEF_VECT,inode,igaus)
       end do
    end do
    !
    ! Newton-Raphson
    !
    if( kfl_linea_nsi == 2 ) then
       if( ndime == 2 ) then
          do igaus = 1,pgaus
             fact0(DEF_VECT) = gpden(DEF_VECT,igaus) * gpvol(DEF_VECT,igaus)
             do inode = 1,pnode
                fact1 = fact0(DEF_VECT) * gpsha(DEF_VECT,inode,igaus)
                elrbu(DEF_VECT,1,inode) = elrbu(DEF_VECT,1,inode) + fact1(DEF_VECT) * ( gpadv(DEF_VECT,1,igaus) * gpgve(DEF_VECT,1,1,igaus) + gpadv(DEF_VECT,2,igaus) * gpgve(DEF_VECT,2,1,igaus) )
                elrbu(DEF_VECT,2,inode) = elrbu(DEF_VECT,2,inode) + fact1(DEF_VECT) * ( gpadv(DEF_VECT,1,igaus) * gpgve(DEF_VECT,1,2,igaus) + gpadv(DEF_VECT,2,igaus) * gpgve(DEF_VECT,2,2,igaus) )
             end do
          end do
       else
          do igaus = 1,pgaus
             fact0(DEF_VECT) = gpden(DEF_VECT,igaus) * gpvol(DEF_VECT,igaus)
             do inode = 1,pnode
                fact1(DEF_VECT) = fact0(DEF_VECT) * gpsha(DEF_VECT,inode,igaus)
                elrbu(DEF_VECT,1,inode) = elrbu(DEF_VECT,1,inode) &
                     & + fact1(DEF_VECT) * (   gpadv(DEF_VECT,1,igaus) * gpgve(DEF_VECT,1,1,igaus) &
                     &                       + gpadv(DEF_VECT,2,igaus) * gpgve(DEF_VECT,2,1,igaus) &
                     &                       + gpadv(DEF_VECT,3,igaus) * gpgve(DEF_VECT,3,1,igaus) )
                elrbu(DEF_VECT,2,inode) = elrbu(DEF_VECT,2,inode) &
                     & + fact1(DEF_VECT) * (   gpadv(DEF_VECT,1,igaus) * gpgve(DEF_VECT,1,2,igaus) &
                     &                       + gpadv(DEF_VECT,2,igaus) * gpgve(DEF_VECT,2,2,igaus) &
                     &                       + gpadv(DEF_VECT,3,igaus) * gpgve(DEF_VECT,3,2,igaus) )
                elrbu(DEF_VECT,3,inode) = elrbu(DEF_VECT,3,inode) &
                     & + fact1(DEF_VECT) * (   gpadv(DEF_VECT,1,igaus) * gpgve(DEF_VECT,1,3,igaus) &
                     &                       + gpadv(DEF_VECT,2,igaus) * gpgve(DEF_VECT,2,3,igaus) &
                     &                       + gpadv(DEF_VECT,3,igaus) * gpgve(DEF_VECT,3,3,igaus) )
             end do
          end do
       end if
    end if
    !
    ! Projection: ( L*(v) , P )
    !
    if( kfl_stabi_nsi == 1 ) then
       do igaus = 1, pgaus
          do inode = 1,pnode
             do idime = 1,ndime
                elrbu(DEF_VECT,idime,inode) = elrbu(DEF_VECT,idime,inode) - gprhs_sgs(DEF_VECT,idime,igaus) * gpsha(DEF_VECT,inode,igaus) * gpvol(DEF_VECT,igaus)
             end do
          end do
       end do
    end if

    !--------------------------------------------------------------------
    !
    ! Pressure bubble
    !
    !--------------------------------------------------------------------

    if( maxval(pbubl) == 1 ) then
       !
       ! Initialization
       !
       elauq = 0.0_rp
       elapq = 0.0_rp
       elaqu = 0.0_rp
       elaqp = 0.0_rp
       elaqq = 0.0_rp
       elrbq = 0.0_rp
       !
       ! Test functions
       !
       do igaus = 1,pgaus
          p2sca_bub(DEF_VECT,igaus) = gptt2(DEF_VECT,igaus) * gpvol(DEF_VECT,igaus) * gpsha_bub(DEF_VECT,igaus)
          do idime = 1,ndime
             p2vec_bub(DEF_VECT,idime,igaus) = gpsp1(DEF_VECT,igaus) * gpvol(DEF_VECT,igaus) * gpcar_bub(DEF_VECT,idime,igaus)
          end do
       end do
       !
       ! Momentum equations
       !
       if( kfl_press_nsi == 1 ) then
          do igaus = 1,pgaus
             fact2(DEF_VECT) = gpvol(DEF_VECT,igaus) * gpsha_bub(DEF_VECT,igaus)
             do inode = 1,pnode
                fact1(DEF_VECT) = -gpvol(DEF_VECT,igaus) * gpsha(DEF_VECT,inode,igaus) + p1vec(DEF_VECT,inode,igaus)
                idof1 = (inode-1)*ndime
                do idime = 1,ndime
                   idofv = idof1 + idime
                   elauq(DEF_VECT,idofv,1) = elauq(DEF_VECT,idofv,1) - fact2(DEF_VECT) * gpcar(DEF_VECT,idime,inode,igaus) &
                        &                                            + fact1(DEF_VECT) * gpcar_bub(DEF_VECT,idime,igaus)
                end do
             end do
          end do
       else
          do igaus = 1,pgaus
             do inode = 1,pnode
                idof1 = (inode-1)*ndime
                do idime = 1,ndime
                   idofv = idof1 + idime
                   elauq(DEF_VECT,idofv,1) = elauq(DEF_VECT,idofv,1) + p1vec(DEF_VECT,inode,igaus) * gpcar_bub(DEF_VECT,idime,igaus)
                end do
             end do
          end do
       end if
       if( kfl_p1ve2_nsi /= 0 ) then
          do igaus = 1,pgaus
             do inode = 1,pnode
                idof1 = (inode-1)*ndime
                do idime = 1,ndime
                   idofv = idof1 + idime
                   do jdime = 1,ndime
                      elauq(DEF_VECT,idofv,1) = elauq(DEF_VECT,idofv,1) + &
                           gpcar_bub(DEF_VECT,jdime,igaus) * p1ve2(DEF_VECT,idime,jdime,inode,igaus)
                   end do
                end do
             end do
          end do
       end if
       !
       ! Bubble equation
       !
       do igaus = 1,pgaus
          do inode = 1,pnode
             do idime = 1,ndime
                idof1 = (inode-1)*ndime+idime
                elaqu(DEF_VECT,1,idof1) = elaqu(DEF_VECT,1,idof1) + rcont(DEF_VECT,idime,inode,igaus) * p2sca_bub(DEF_VECT,igaus) &
                     &                                            + rmomu(DEF_VECT,inode,igaus)       * p2vec_bub(DEF_VECT,idime,igaus)
             end do
          end do
       end do
       if( kfl_rmom2_nsi /= 0 ) then
          do igaus = 1,pgaus
             do jnode = 1,pnode
                do jdime = 1,ndime
                   jdofv = (jnode-1)*ndime+jdime
                   do kdime = 1,ndime
                      elaqu(DEF_VECT,1,jdofv) = elaqu(DEF_VECT,1,jdofv) + rmom2(DEF_VECT,kdime,jdime,jnode,igaus) * p2vec_bub(DEF_VECT,kdime,igaus)
                   end do
                end do
             end do
          end do
       end if
       !
       ! Also pressure equation...
       !
       do igaus = 1,pgaus
          do inode = 1,pnode
             do idime = 1,ndime
                elapq(DEF_VECT,inode,1) = elapq(DEF_VECT,inode,1) + p2vec(DEF_VECT,idime,inode,igaus) * gpcar_bub(DEF_VECT,idime,igaus)
                elaqp(DEF_VECT,1,inode) = elaqp(DEF_VECT,1,inode) + p2vec_bub(DEF_VECT,idime,igaus)   * gpcar(DEF_VECT,idime,inode,igaus)
             end do
          end do
          do idime = 1,ndime
             elaqq(DEF_VECT,1,1)             = elaqq(DEF_VECT,1,1)             + p2vec_bub(DEF_VECT,idime,igaus) * gpcar_bub(DEF_VECT,idime,igaus)
             gprhs_sgs(DEF_VECT,idime,igaus) = gprhs_sgs(DEF_VECT,idime,igaus) + taupr(DEF_VECT,igaus) * gpvep(DEF_VECT,idime,igaus)
             gprhh(DEF_VECT,idime,igaus)     = gprhs(DEF_VECT,idime,igaus)     + gprhs_sgs(DEF_VECT,idime,igaus)
             elrbq(DEF_VECT,1)               = elrbq(DEF_VECT,1)               + p2vec_bub(DEF_VECT,idime,igaus) * gprhh(DEF_VECT,idime,igaus)
          end do
       end do
       do igaus = 1,pgaus
          elrbq(DEF_VECT,1) = elrbq(DEF_VECT,1) + gprhc(DEF_VECT,igaus) * p2sca_bub(DEF_VECT,igaus) ! ( rhs , q_bub)
       end do
       !
       ! Penalization
       !
       stop
       !do igaus = 1,pgaus
       !   elaqq(DEF_VECT,1,1) = elaqq(DEF_VECT,1,1) + penal_nsi * gpvol(DEF_VECT,igaus) * gpsha_bub(DEF_VECT,igaus)
       !   elrbq(DEF_VECT,1)   = elrbq(DEF_VECT,1)   + penal_nsi * gpvol(DEF_VECT,igaus) * gpsha_bub(DEF_VECT,igaus) * elbub(DEF_VECT)
       !end do
       !
    end if

  end subroutine nsi_element_assembly_asgs_oss_old

  !-----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @brief   Output Navier-Stokes system
  !> @details Output Navier-Stokes elemental system
  !>
  !-----------------------------------------------------------------------

  subroutine nsi_element_system_output(&
       pnode,elauu,elaup,elapp,elapu,elrbu,elrbp,&
       elauq,elapq,elaqu,elaqp,elaqq,elrbq)

    integer(ip), intent(in)           :: pnode
    ! Element matrices
    real(rp),    intent(in)           :: elauu(pnode*ndime,pnode*ndime)
    real(rp),    intent(in)           :: elaup(pnode*ndime,pnode)
    real(rp),    intent(in)           :: elapp(pnode,pnode)
    real(rp),    intent(in)           :: elapu(pnode,pnode*ndime)
    real(rp),    intent(in)           :: elrbu(ndime,pnode)
    real(rp),    intent(in)           :: elrbp(pnode)
    ! Enrichement Element matrices
    real(rp),    intent(in), optional :: elauq(pnode*ndime,1)
    real(rp),    intent(in), optional :: elapq(pnode,1)
    real(rp),    intent(in), optional :: elaqu(1,pnode*ndime)
    real(rp),    intent(in), optional :: elaqp(1,pnode)
    real(rp),    intent(in), optional :: elaqq(1,1)
    real(rp),    intent(in), optional :: elrbq(1)

    integer(ip),   parameter          :: lreal=9
    character(13)                     :: FMT1
    integer(ip)                       :: inode,idime,idofn,jnode,jdime,jdofn
    character(3)                      :: vanam
    logical(lg)                       :: bubble

    bubble = .false.
    if( kfl_bubbl_nsi /= 0 .and. present(elauq) ) bubble = .true.
    !
    ! Format
    !
    FMT1 = '(1x,e'//trim(intost(lreal))//'.'//trim(intost(lreal-7))
    !
    ! First line
    !
    write(99,'(a)',advance='no') '+-'
    do idofn = 1,pnode*ndime
       do idime = 1,lreal+1
          write(99,'(a)',advance='no') ' '
       end do
    end do
    write(99,'(a)',advance='no') ' '
    do idofn = 1,pnode
       do idime = 1,lreal+1
          write(99,'(a)',advance='no') ' '
       end do
    end do

    if( bubble ) then
       do idime = 1,lreal+4
          write(99,'(a)',advance='no') ' '
       end do
    end if

    write(99,'(a)',advance='no') '-+'
    write(99,'(a)',advance='no') '  +-   -+'

    write(99,'(a)',advance='no') '     +-'
    do idime = 1,lreal
       write(99,'(a)',advance='no') ' '
    end do
    write(99,'(a)',advance='no') '-+'

    write(99,*)
    flush(99)
    
    do inode = 1,pnode
       do idime = 1,ndime
          !
          ! Auu
          !
          write(99,'(a)',advance='no') '|'
          idofn = (inode-1)*ndime+idime
          do jnode = 1,pnode
             do jdime = 1,ndime
                jdofn = (jnode-1)*ndime+jdime
                write(99,FMT1,advance='no') elauu(idofn,jdofn)
             end do
          end do
          !
          ! Aup
          !
          write(99,'(a)',advance='no') ' |'
          do jnode = 1,pnode
             write(99,FMT1,advance='no') elaup(idofn,jnode)
          end do
          !
          ! Auq
          !
          if( bubble ) then
             write(99,'(a)',advance='no') ' ||'
             write(99,FMT1,advance='no') elauq(idofn,1)
          end if
          !
          ! ui
          !
          write(99,'(a)',advance='no') ' |'
          if( idime == 1 ) then
             vanam = 'u'
          else if( idime == 2 ) then
             vanam = 'v'
          else if( idime == 3 ) then
             vanam = 'w'
          end if
          vanam = trim(vanam) // trim(intost(inode))
          write(99,'(a)',advance='no') '  | '//vanam//' |'
          !
          ! bu
          !
          write(99,'(a)',advance='no') '     |'
          write(99,FMT1,advance='no') elrbu(idime,inode)
          write(99,'(a)',advance='no') ' |'

          write(99,*)
       end do
    end do
    !
    ! u -> p block
    !
    write(99,'(a)',advance='no') '| '
    do idofn = 1,pnode*ndime
       do idime = 1,lreal+1
          write(99,'(a)',advance='no') '-'
       end do
    end do
    write(99,'(a)',advance='no') '+'
    do idofn = 1,pnode
       do idime = 1,lreal+1
          write(99,'(a)',advance='no') '-'
       end do
    end do

    if( bubble ) then
       write(99,'(a)',advance='no') '-++'
       do idime = 1,lreal+1
          write(99,'(a)',advance='no') '-'
       end do
    end if

    write(99,'(a)',advance='no') ' |'
    write(99,'(a)',advance='no') '  | --- |  =  | '
    do idime = 1,lreal
       write(99,'(a)',advance='no') '-'
    end do
    write(99,'(a)',advance='no') ' |'

    write(99,*)

    do inode = 1,pnode
       !
       ! Apu
       !
       write(99,'(a)',advance='no') '|'
       do jnode = 1,pnode
          do jdime = 1,ndime
             jdofn = (jnode-1)*ndime+jdime
             write(99,FMT1,advance='no') elapu(inode,jdofn)
          end do
       end do
       !
       ! App
       !
       write(99,'(a)',advance='no') ' |'
       do jnode = 1,pnode
          write(99,FMT1,advance='no') elapp(inode,jnode)
       end do
       !
       ! Apq
       !
       if( bubble ) then
          write(99,'(a)',advance='no') ' ||'
          write(99,FMT1,advance='no') elapq(inode,1)
       end if
       !
       ! pi
       !
       write(99,'(a)',advance='no') ' |'
       vanam = 'p' // trim(intost(inode))
       write(99,'(a)',advance='no') '  | '//vanam//' |'
       !
       ! bp
       !
       write(99,'(a)',advance='no') '     |'
       write(99,FMT1,advance='no') elrbp(inode)
       write(99,'(a)',advance='no') ' |'

       write(99,*)
    end do
    !
    ! p -> q block
    !
    write(99,'(a)',advance='no') '| '
    do idofn = 1,pnode*ndime
       do idime = 1,lreal+1
          write(99,'(a)',advance='no') '='
       end do
    end do
    write(99,'(a)',advance='no') '+'
    do idofn = 1,pnode
       do idime = 1,lreal+1
          write(99,'(a)',advance='no') '='
       end do
    end do

    if( bubble ) then
       write(99,'(a)',advance='no') '=++'
       do idime = 1,lreal+1
          write(99,'(a)',advance='no') '='
       end do
    end if

    write(99,'(a)',advance='no') ' |'
    write(99,'(a)',advance='no') '  | === |  =  | '
    do idime = 1,lreal
       write(99,'(a)',advance='no') '='
    end do
    write(99,'(a)',advance='no') ' |'

    write(99,*)
    !
    ! Aqu
    !
    inode = 1
    write(99,'(a)',advance='no') '|'
    do jnode = 1,pnode
       do jdime = 1,ndime
          jdofn = (jnode-1)*ndime+jdime
          write(99,FMT1,advance='no') elaqu(inode,jdofn)
       end do
    end do
    !
    ! App
    !
    write(99,'(a)',advance='no') ' |'
    do jnode = 1,pnode
       write(99,FMT1,advance='no') elaqp(inode,jnode)
    end do
    !
    ! Apq
    !
    if( bubble ) then
       write(99,'(a)',advance='no') ' ||'
       write(99,FMT1,advance='no') elaqq(inode,1)
    end if
    !
    ! pi
    !
    write(99,'(a)',advance='no') ' |'
    vanam = 'q  '
    write(99,'(a)',advance='no') '  | '//vanam//' |'
    !
    ! bp
    !
    write(99,'(a)',advance='no') '     |'
    write(99,FMT1,advance='no') elrbq(inode)
    write(99,'(a)',advance='no') ' |'

    write(99,*)
    !
    ! Last line
    !
    write(99,'(a)',advance='no') '+-'
    do idofn = 1,pnode*ndime
       do idime = 1,lreal+1
          write(99,'(a)',advance='no') ' '
       end do
    end do

    write(99,'(a)',advance='no') ' '
    do idofn = 1,pnode
       do idime = 1,lreal+1
          write(99,'(a)',advance='no') ' '
       end do
    end do

    if( bubble ) then
       do idime = 1,lreal+4
          write(99,'(a)',advance='no') ' '
       end do
    end if

    write(99,'(a)',advance='no') '-+'
    write(99,'(a)',advance='no') '  +-   -+'

    write(99,'(a)',advance='no') '     +-'
    do idime = 1,lreal
       write(99,'(a)',advance='no') ' '
    end do
    write(99,'(a)',advance='no') '-+'

    write(99,*)
    flush(99)
    
  end subroutine nsi_element_system_output

end module mod_nsi_element_assembly
!> @}
