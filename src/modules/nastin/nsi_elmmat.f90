!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



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
!> @} 
!------------------------------------------------------------------------

subroutine nsi_elmmat(                                &
     pnode,pgaus,pevat,gpden,gpvis,gppor,gpgvi,gpsp1, &
     gptt1,gpsp2,gptt2,gpvol,gpsha,gpcar,gplap,gphes, &
     gpadv,gpvep,gprhs,gprhc,rmomu,rcont,p1vec,p2vec, &
     p2sca,wgrgr,wgrvi,elauu,elaup,elapp,elapu,elrbu, &
     elrbp,rmom2,p1ve2,gpst1,gpgve,gprh2,gprhs_sgs,   &
     elvel,ellum,dtinv_loc,pbubl,gpsha_bub,gpcar_bub, &
     gppre,elauq,elapq,elaqu,elaqp,elaqq,elrbq,elpre, &
     elbub, densi)

  use def_kintyp, only       :  ip,rp
  use def_domain, only       :  ndime,mnode,ntens
  use def_nastin, only       :  kfl_stabi_nsi,nbdfp_nsi,NSI_FRACTIONAL_STEP
  use def_nastin, only       :  fvins_nsi,fcons_nsi,bemol_nsi,kfl_regim_nsi,pabdf_nsi
  use def_nastin, only       :  fvela_nsi,kfl_rmom2_nsi,kfl_press_nsi
  use def_nastin, only       :  kfl_p1ve2_nsi,kfl_linea_nsi
  use def_nastin, only       :  kfl_confi_nsi
  use def_nastin, only       :  kfl_nota1_nsi,penal_nsi
  use def_nastin, only       :  kfl_sgsli_nsi
  use def_nastin, only       :  nbdfp_nsi, pabdf_nsi
  use def_master, only       :  kfl_lumped
  use def_kermod, only       :  kfl_duatss,fact_duatss
  use mod_maths,  only       :  maths_schur_complement
  implicit none

  integer(ip), intent(in)    :: pnode,pgaus,pevat
  real(rp),    intent(in)    :: gpden(pgaus),gpvis(pgaus)
  real(rp),    intent(in)    :: gppor(pgaus),gpgvi(ndime,pgaus)
  real(rp),    intent(inout) :: gpsp1(pgaus),gptt1(pgaus)
  real(rp),    intent(inout) :: gpsp2(pgaus),gptt2(pgaus)
  real(rp),    intent(in)    :: gpvol(pgaus),gpsha(pnode,pgaus)
  real(rp),    intent(in)    :: gpcar(ndime,mnode,pgaus)
  real(rp),    intent(in)    :: gpadv(ndime,pgaus)
  real(rp),    intent(inout) :: gpvep(ndime,pgaus)
  real(rp),    intent(in)    :: gplap(pnode,pgaus)
  real(rp),    intent(in)    :: gphes(ntens,mnode,pgaus)
  real(rp),    intent(in)    :: gprhs(ndime,pgaus)
  real(rp),    intent(inout) :: gprhs_sgs(ndime,pgaus)
  real(rp),    intent(in)    :: gprhc(pgaus)
  real(rp),    intent(in)    :: gprh2(pgaus)
  real(rp),    intent(in)    :: rmomu(pnode,pgaus)
  real(rp),    intent(in)    :: rcont(ndime,pnode,pgaus)
  real(rp),    intent(out)   :: p1vec(pnode,pgaus)
  real(rp),    intent(out)   :: p2vec(ndime,pnode,pgaus)
  real(rp),    intent(out)   :: p2sca(pnode,pgaus)
  real(rp),    intent(out)   :: wgrgr(pnode,pnode,pgaus)
  real(rp),    intent(out)   :: wgrvi(pnode,pgaus) 
  ! Element matrices
  real(rp),    intent(out)   :: elauu(pnode*ndime,pnode*ndime)
  real(rp),    intent(out)   :: elaup(pnode*ndime,pnode)
  real(rp),    intent(out)   :: elapp(pnode,pnode)
  real(rp),    intent(out)   :: elapu(pnode,pnode*ndime)
  real(rp),    intent(out)   :: elrbu(ndime,pnode)
  real(rp),    intent(out)   :: elrbp(pnode)
  real(rp),    intent(in)    :: rmom2(ndime,ndime,pnode,pgaus)
  real(rp),    intent(out)   :: p1ve2(ndime,ndime,pnode,pgaus)
  real(rp),    intent(in)    :: gpst1(pgaus)
  real(rp),    intent(in)    :: gpgve(ndime,ndime,pgaus)
  real(rp),    intent(in)    :: elvel(ndime, pnode, *)
  real(rp),    intent(out)   :: ellum(pnode)
  real(rp),    intent(in)    :: dtinv_loc
  integer(ip), intent(in)    :: pbubl
  real(rp),    intent(in)    :: gpsha_bub(pgaus)                    ! Ne
  real(rp),    intent(in)    :: gpcar_bub(ndime,pgaus)              ! dNe/dxi
  real(rp),    intent(in)    :: gppre(pgaus) 
  ! Enrichement Element matrices
  real(rp),    intent(out)   :: elauq(pnode*ndime,1)
  real(rp),    intent(out)   :: elapq(pnode,1)
  real(rp),    intent(out)   :: elaqu(1,pnode*ndime)
  real(rp),    intent(out)   :: elaqp(1,pnode)
  real(rp),    intent(out)   :: elaqq(1,1)
  real(rp),    intent(out)   :: elrbq(1)
  real(rp),    intent(in)    :: elpre(pnode,*)
  real(rp),    intent(in)    :: elbub
  real(rp),    intent(in)    :: densi(pgaus,nbdfp_nsi)

  integer(ip)                :: inode,jnode,idofv
  integer(ip)                :: idof1,jdime,itime
  integer(ip)                :: igaus,idime,jdofv,kdime
  real(rp)                   :: fact2,fact3,fact4,fact5,fact6
  real(rp)                   :: factx,facty,factz,fact8,fact0
  real(rp)                   :: facx1,facy1,fact1,xvisc
  real(rp)                   :: xvis2,fact7,penal
  real(rp)                   :: gprhh(ndime,pgaus), taupr(pgaus), gpveo(ndime)
  real(rp)                   :: gpcar1ji, gpcar2ji, gpcar3ji
  real(rp)                   :: p2sca_bub(pgaus),p2vec_bub(ndime,pgaus), one_rp                  
  real(rp), dimension(4*((pnode+3)/4),pnode) :: elauu11,elauu21,elauu31,elauu12,elauu22,elauu32,elauu13,elauu23,elauu33
  real(rp), dimension(4*((pnode+3)/4))       :: factvec1,factvec2,factvec3,factvec4
  real(rp)                   :: auxva, sgsli(ndime), gpvel(ndime)
  real(rp)                   :: dtinv_mod
  !yor! the comments marked with 'yor' explain the recent changes for optimization purpose, within EoCoE project. To be removed later.
  !yor! optimization test : for memory padding purpose : 4*((pnode+3)/4) will give the closest number >= pnode which is a multiple of 4.
  !yor! We have to think of a clean way to write it. The idea could be expanded in the whole code, with the padding size as parameter.

  if( NSI_FRACTIONAL_STEP ) then
     dtinv_mod = 0.0_rp
  else
     dtinv_mod = dtinv_loc
  end if  
  
  !----------------------------------------------------------------------
  !
  ! Initialization
  !
  !----------------------------------------------------------------------

  !yor! rewrote initialization
  elrbp = 0.0_rp
  elrbu = 0.0_rp
  elapp = 0.0_rp
  !yor! For performance and vectorization purpose, elauu is built as 4 to 9 separate tables that will be reassembled together at the end
  !yor! Another idea, for same purpose, would be to use a single array : elauu_new(pnode,pnode, ndime, ndime)   
  elauu11 = 0.0_rp
  elauu21 = 0.0_rp
  elauu12 = 0.0_rp
  elauu22 = 0.0_rp
  if( ndime == 3 ) then
     elauu31 = 0.0_rp
     elauu32 = 0.0_rp
     elauu13 = 0.0_rp
     elauu23 = 0.0_rp
     elauu33 = 0.0_rp
  end if

  elaup = 0.0_rp
  elapu = 0.0_rp
  !
  xvis2 = 0.0_rp
  if( fvins_nsi > 0.9_rp ) then
     xvisc = 1.0_rp
     if( fvins_nsi > 1.9_rp ) xvis2 = 2.0_rp/3.0_rp
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
     taupr(1:pgaus) = 1.0_rp
  else
     taupr(1:pgaus) = 1.0_rp / gpst1(1:pgaus)
  end if

  one_rp =1.0_rp
  if (kfl_nota1_nsi==1) one_rp=0.0_rp  ! do not stabilize convection, reaction, coriolis


  do igaus = 1,pgaus
     fact1 = gpsp1(igaus) *one_rp * gpden(igaus)             * gpvol(igaus) 
     fact2 = (gptt1(igaus)-gpsp1(igaus)*one_rp*gppor(igaus)) * gpvol(igaus)
     fact3 = gpsp1(igaus) * gpvol(igaus)
     fact4 = gptt2(igaus) * gpvol(igaus)
     do inode = 1,pnode
        p1vec(inode,igaus)          = fact2 * gpsha(inode,igaus) + fact1 * dot_product(gpadv(1:ndime,igaus),gpcar(1:ndime,inode,igaus))
        p2sca(inode,igaus)          = fact4 * gpsha(inode,igaus)
        p2vec(1:ndime,inode,igaus)  = fact3 * gpcar(1:ndime,inode,igaus)
        wgrvi(inode,igaus)          = dot_product(gpgvi(1:ndime,igaus),gpcar(1:ndime,inode,igaus)) & 
             &                        + gpvis(igaus) * gplap(inode,igaus)
        do jnode = 1,pnode
           wgrgr(inode,jnode,igaus) = dot_product(gpcar(1:ndime,inode,igaus),gpcar(1:ndime,jnode,igaus))
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
     p1ve2(1:ndime,1:ndime,1:pnode,1:pgaus) = 0.0_rp
     if( ndime == 2 ) then
        do igaus = 1,pgaus           
           fact3 = 2.0_rp * gpden(igaus) * fvela_nsi(3)
           factz = gpsp1(igaus)* fact3 * gpvol(igaus)
           do inode = 1,pnode 
              fact1                  =   factz * gpsha(inode,igaus)
              p1ve2(1,2,inode,igaus) =   fact1
              p1ve2(2,1,inode,igaus) = - fact1
           end do
        end do
     else
        do igaus = 1,pgaus   
           if( NSI_FRACTIONAL_STEP ) then
              fact8 = 2.0_rp * gpvol(igaus)/dtinv_loc ! tau = dtime when corio
           else             
              fact8 = 2.0_rp * gpden(igaus) * gpsp1(igaus)* gpvol(igaus)
           end if
           factx = fact8  * fvela_nsi(1)
           facty = fact8  * fvela_nsi(2)
           factz = fact8  * fvela_nsi(3)
           do inode = 1,pnode
              p1ve2(1,2,inode,igaus) =  factz * gpsha(inode,igaus)  ! x-equation
              p1ve2(1,3,inode,igaus) = -facty * gpsha(inode,igaus)  ! x-equation
              p1ve2(2,1,inode,igaus) = -factz * gpsha(inode,igaus)  ! y-equation
              p1ve2(2,3,inode,igaus) =  factx * gpsha(inode,igaus)  ! y-equation
              p1ve2(3,1,inode,igaus) =  facty * gpsha(inode,igaus)  ! z-equation
              p1ve2(3,2,inode,igaus) = -factx * gpsha(inode,igaus)  ! z-equation    
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
  if( ndime == 2 ) then

     do igaus = 1,pgaus

        fact0 = ( gpsp2(igaus) - xvis2 * gpvis(igaus) ) * gpvol(igaus)                        ! (1) + (5)
        fact6 = gpvis(igaus) * gpvol(igaus)
        do inode = 1,pnode
           fact1 = fact0 * gpcar(1,inode,igaus) 
           fact2 = fact0 * gpcar(2,inode,igaus)
           fact5 = gpsha(inode,igaus) * gpvol(igaus)

           do jnode = 1,pnode              
              fact4              =   p1vec(inode,igaus) * rmomu(jnode,igaus) &       ! (2):       ( rmomu(u) , p1vec(v) ) 
                   &               + fact5 * wgrvi(jnode,igaus)              &       ! (3) + (6): ( mu d^2ui/dxj^2 + dmu/dxj dui/dxj, vi )
                   &               + fact6 * wgrgr(inode,jnode,igaus)                ! (7):       ( mu dui/dxj , dvi/dxj ) 

              elauu11(inode,jnode) = elauu11(inode,jnode) + fact1 * rcont(1,jnode,igaus) + fact4  ! Auu_xx
              elauu21(inode,jnode) = elauu21(inode,jnode) + fact2 * rcont(1,jnode,igaus)          ! Auu_yx
              elauu12(inode,jnode) = elauu12(inode,jnode) + fact1 * rcont(2,jnode,igaus)          ! Auu_xy
              elauu22(inode,jnode) = elauu22(inode,jnode) + fact2 * rcont(2,jnode,igaus) + fact4  ! Auu_yy
           end do
        end do
     end do

  else
  
     do igaus = 1,pgaus    
     
        fact0 = ( gpsp2(igaus) - xvis2 * gpvis(igaus) ) * gpvol(igaus)                        ! (1) + (5)
        fact6 = gpvis(igaus) * gpvol(igaus)

        !yor! [originally line 344] loops were reverted and refactored, with a vectorization effort
        factvec1(1:pnode) = gpcar(1,1:pnode,igaus)
        factvec2(1:pnode) = gpcar(2,1:pnode,igaus)
        factvec3(1:pnode) = gpcar(3,1:pnode,igaus)

        do jnode = 1,pnode              
           fact1 = fact0*rcont(1,jnode,igaus) 
           fact2 = fact0*rcont(2,jnode,igaus)
           fact3 = fact0*rcont(3,jnode,igaus)
           fact5 = wgrvi(jnode,igaus)*gpvol(igaus)

           factvec4(1:pnode) = p1vec(1:pnode,igaus) * rmomu(jnode,igaus) + fact5 * gpsha(1:pnode,igaus) &
                + fact6 * wgrgr(1:pnode,jnode,igaus)

           do inode = 1,pnode
              elauu11(inode,jnode) = elauu11(inode,jnode) + fact1 * factvec1(inode) + factvec4(inode)  ! Auu_xx
              elauu12(inode,jnode) = elauu12(inode,jnode) + fact2 * factvec1(inode)                    ! Auu_xy
              elauu13(inode,jnode) = elauu13(inode,jnode) + fact3 * factvec1(inode)                    ! Auu_xz

              elauu21(inode,jnode) = elauu21(inode,jnode) + fact1 * factvec2(inode)                    ! Auu_yx
              elauu22(inode,jnode) = elauu22(inode,jnode) + fact2 * factvec2(inode) + factvec4(inode)  ! Auu_yy
              elauu23(inode,jnode) = elauu23(inode,jnode) + fact3 * factvec2(inode)                    ! Auu_yz

              elauu31(inode,jnode) = elauu31(inode,jnode) + fact1 * factvec3(inode)                    ! Auu_zx
              elauu32(inode,jnode) = elauu32(inode,jnode) + fact2 * factvec3(inode)                    ! Auu_zy
              elauu33(inode,jnode) = elauu33(inode,jnode) + fact3 * factvec3(inode) + factvec4(inode)  ! Auu_zz
           end do
        end do

     end do

  end if

  if( fvins_nsi > 0.9_rp ) then

     if( ndime == 2 ) then

        do igaus = 1,pgaus

           fact4 = gpgvi(1,igaus) * gpvol(igaus)
           fact5 = gpgvi(2,igaus) * gpvol(igaus)
           fact7 = gpvis(igaus)   * gpvol(igaus)
           !           factx = xvis2 * fact4
           !           facty = xvis2 * fact5

           do inode = 1,pnode
              fact0 = fact7 * gpsha(inode,igaus) *(1.0_rp-xvis2)
              fact1 = fact4 * gpsha(inode,igaus)
              fact2 = fact5 * gpsha(inode,igaus)
              facx1 = fact1 * xvis2
              facy1 = fact2 * xvis2 
              !              facx1 = factx * gpsha(inode,igaus)
              !              facy1 = facty * gpsha(inode,igaus)         
              do jnode = 1,pnode              
                 !
                 !                                         (4)  ( mu duj/dxi , dvi/dxj )
                 elauu11(inode,jnode) = elauu11(inode,jnode) + fact7 * gpcar(1,inode,igaus) * gpcar(1,jnode,igaus)         ! Auu_xx
                 elauu12(inode,jnode) = elauu12(inode,jnode) + fact7 * gpcar(2,inode,igaus) * gpcar(1,jnode,igaus)         ! Auu_xy
                 elauu21(inode,jnode) = elauu21(inode,jnode) + fact7 * gpcar(1,inode,igaus) * gpcar(2,jnode,igaus)         ! Auu_yx
                 elauu22(inode,jnode) = elauu22(inode,jnode) + fact7 * gpcar(2,inode,igaus) * gpcar(2,jnode,igaus)         ! Auu_yy
                 !
                 !                                         (10) - 2/3 * B * ( dmu/dxi (div u) , vi )
                 elauu11(inode,jnode) = elauu11(inode,jnode) - facx1 * gpcar(1,jnode,igaus)                                ! Auu_xx
                 elauu12(inode,jnode) = elauu12(inode,jnode) - facx1 * gpcar(2,jnode,igaus)                                ! Auu_xy
                 elauu21(inode,jnode) = elauu21(inode,jnode) - facy1 * gpcar(1,jnode,igaus)                                ! Auu_yx
                 elauu22(inode,jnode) = elauu22(inode,jnode) - facy1 * gpcar(2,jnode,igaus)                                ! Auu_yy
                 !
                 !                                         (9)  ( dmu/dxj duj/dxi , vi ) + (8)  ( mu d^2uj/dxidxj , vi )
                 !                                               + (11)   - 2/3 * B * ( mu d(div u)/dxi , vi )
                 elauu11(inode,jnode) = elauu11(inode,jnode) + fact1 * gpcar(1,jnode,igaus) + fact0 * gphes(1,jnode,igaus) ! Auu_xx
                 elauu12(inode,jnode) = elauu12(inode,jnode) + fact2 * gpcar(1,jnode,igaus) + fact0 * gphes(3,jnode,igaus) ! Auu_xy
                 elauu21(inode,jnode) = elauu21(inode,jnode) + fact1 * gpcar(2,jnode,igaus) + fact0 * gphes(3,jnode,igaus) ! Auu_yx
                 elauu22(inode,jnode) = elauu22(inode,jnode) + fact2 * gpcar(2,jnode,igaus) + fact0 * gphes(2,jnode,igaus) ! Auu_yy
                 !
                 !                                        

              end do

           end do
        end do

     else

        do igaus = 1,pgaus

           fact4 = gpgvi(1,igaus) * gpvol(igaus)
           fact5 = gpgvi(2,igaus) * gpvol(igaus)
           fact6 = gpgvi(3,igaus) * gpvol(igaus)
           fact7 = gpvis(igaus)   * gpvol(igaus)
           !          factx = xvis2 * fact4
           !          facty = xvis2 * fact5
           !          factz = xvis2 * fact6

           !yor! [Originally line 451] reverted loops order, refactored the factors, with an effort to vectorize each instruction 
           !
           !                                (4)  ( mu duj/dxi , dvi/dxj )
           factvec1(1:pnode) = gpcar(1,1:pnode,igaus)
           factvec2(1:pnode) = gpcar(2,1:pnode,igaus)
           factvec3(1:pnode) = gpcar(3,1:pnode,igaus)

           do jnode = 1,pnode
              gpcar1ji = gpcar(1,jnode,igaus)
              gpcar2ji = gpcar(2,jnode,igaus)
              gpcar3ji = gpcar(3,jnode,igaus)
              !
              do inode = 1,pnode
                 elauu11(inode,jnode) = elauu11(inode,jnode) + fact7 * factvec1(inode) * gpcar1ji         ! Auu_xx
                 elauu21(inode,jnode) = elauu21(inode,jnode) + fact7 * factvec1(inode) * gpcar2ji         ! Auu_yx
                 elauu31(inode,jnode) = elauu31(inode,jnode) + fact7 * factvec1(inode) * gpcar3ji         ! Auu_zx

                 elauu12(inode,jnode) = elauu12(inode,jnode) + fact7 * factvec2(inode) * gpcar1ji         ! Auu_xy
                 elauu22(inode,jnode) = elauu22(inode,jnode) + fact7 * factvec2(inode) * gpcar2ji         ! Auu_yy
                 elauu32(inode,jnode) = elauu32(inode,jnode) + fact7 * factvec2(inode) * gpcar3ji         ! Auu_zy

                 elauu13(inode,jnode) = elauu13(inode,jnode) + fact7 * factvec3(inode) * gpcar1ji         ! Auu_xz
                 elauu23(inode,jnode) = elauu23(inode,jnode) + fact7 * factvec3(inode) * gpcar2ji         ! Auu_yz
                 elauu33(inode,jnode) = elauu33(inode,jnode) + fact7 * factvec3(inode) * gpcar3ji         ! Auu_zz
              end do
           end do
           !
           !                     (10) - 2/3 * B * ( dmu/dxi (div u) , vi )


           factvec1(1:pnode) = fact4 * gpsha(1:pnode,igaus) * xvis2 !yor! (ex facx1) this copy is for pure esthaetical purpose
           factvec2(1:pnode) = fact5 * gpsha(1:pnode,igaus) * xvis2 !yor! (ex facy1)
           factvec3(1:pnode) = fact6 * gpsha(1:pnode,igaus) * xvis2 !yor! (ex facz1)  

           do jnode = 1,pnode
              gpcar1ji = gpcar(1,jnode,igaus)
              gpcar2ji = gpcar(2,jnode,igaus)
              gpcar3ji = gpcar(3,jnode,igaus)

              do inode = 1,pnode
                 elauu11(inode,jnode) = elauu11(inode,jnode) - factvec1(inode) * gpcar1ji                                ! Auu_xx 
                 elauu12(inode,jnode) = elauu12(inode,jnode) - factvec1(inode) * gpcar2ji                                ! Auu_xy
                 elauu13(inode,jnode) = elauu13(inode,jnode) - factvec1(inode) * gpcar3ji                                ! Auu_xz

                 elauu21(inode,jnode) = elauu21(inode,jnode) - factvec2(inode) * gpcar1ji                                ! Auu_yx
                 elauu22(inode,jnode) = elauu22(inode,jnode) - factvec2(inode) * gpcar2ji                                ! Auu_yy
                 elauu23(inode,jnode) = elauu23(inode,jnode) - factvec2(inode) * gpcar3ji                                ! Auu_yz

                 elauu31(inode,jnode) = elauu31(inode,jnode) - factvec3(inode) * gpcar1ji                                ! Auu_zx
                 elauu32(inode,jnode) = elauu32(inode,jnode) - factvec3(inode) * gpcar2ji                                ! Auu_zy
                 elauu33(inode,jnode) = elauu33(inode,jnode) - factvec3(inode) * gpcar3ji                                ! Auu_zz
              end do
           end do

           !                                                                   
           !                                         (9)  ( dmu/dxj duj/dxi , vi ) + (8)  ( mu d^2uj/dxidxj , vi )
           !                                                + (11)   - 2/3 * B * ( mu d(div u)/dxi , vi ) 

           factvec4(1:pnode) = fact7 * gpsha(1:pnode,igaus)*(1.0_rp-xvis2) !yor!(fact0)

           factvec1(1:pnode) = fact4 * gpsha(1:pnode,igaus) !yor!(ex fact1)
           factvec2(1:pnode) = fact5 * gpsha(1:pnode,igaus) !yor!(ex fact2)
           factvec3(1:pnode) = fact6 * gpsha(1:pnode,igaus) !yor!(ex fact3)

           do jnode = 1,pnode
              gpcar1ji = gpcar(1,jnode,igaus)
              gpcar2ji = gpcar(2,jnode,igaus)
              gpcar3ji = gpcar(3,jnode,igaus)

              do inode = 1,pnode
                 elauu11(inode,jnode) = elauu11(inode,jnode) + factvec1(inode) * gpcar1ji + factvec4(inode) * gphes(1,jnode,igaus) ! Auu_xx
                 elauu21(inode,jnode) = elauu21(inode,jnode) + factvec1(inode) * gpcar2ji + factvec4(inode) * gphes(4,jnode,igaus) ! Auu_yx
                 elauu31(inode,jnode) = elauu31(inode,jnode) + factvec1(inode) * gpcar3ji + factvec4(inode) * gphes(5,jnode,igaus) ! Auu_zx

                 elauu12(inode,jnode) = elauu12(inode,jnode) + factvec2(inode) * gpcar1ji + factvec4(inode) * gphes(4,jnode,igaus) ! Auu_xy
                 elauu22(inode,jnode) = elauu22(inode,jnode) + factvec2(inode) * gpcar2ji + factvec4(inode) * gphes(2,jnode,igaus) ! Auu_yy
                 elauu32(inode,jnode) = elauu32(inode,jnode) + factvec2(inode) * gpcar3ji + factvec4(inode) * gphes(6,jnode,igaus) ! Auu_zy

                 elauu13(inode,jnode) = elauu13(inode,jnode) + factvec3(inode) * gpcar1ji + factvec4(inode) * gphes(5,jnode,igaus) ! Auu_xz
                 elauu23(inode,jnode) = elauu23(inode,jnode) + factvec3(inode) * gpcar2ji + factvec4(inode) * gphes(6,jnode,igaus) ! Auu_yz
                 elauu33(inode,jnode) = elauu33(inode,jnode) + factvec3(inode) * gpcar3ji + factvec4(inode) * gphes(3,jnode,igaus) ! Auu_zz
              end do
           end do
        end do
     end if
  end if

  if( fcons_nsi > 0.1_rp ) then   ! Skew symmetric when fcons=0.5
     if( ndime == 2 ) then
        do igaus = 1,pgaus
           fact0 = fcons_nsi * gpden(igaus) * gpvol(igaus)
           do inode = 1,pnode
              do jnode = 1,pnode              
                 fact1 = dot_product(gpadv(1:ndime,igaus),gpcar(1:ndime,inode,igaus)) * gpsha(jnode,igaus) * fact0
                 fact2 = dot_product(gpadv(1:ndime,igaus),gpcar(1:ndime,jnode,igaus)) * gpsha(inode,igaus) * fact0

                 elauu11(inode,jnode) = elauu11(inode,jnode) - fact1 - fact2                      ! Auu_xx
                 elauu22(inode,jnode) = elauu22(inode,jnode) - fact1 - fact2                      ! Auu_yy
              end do
           end do
        end do
     else
        do igaus = 1,pgaus
           fact0 = fcons_nsi * gpden(igaus) * gpvol(igaus)
           do inode = 1,pnode
              do jnode = 1,pnode
                 fact1 = dot_product(gpadv(1:ndime,igaus),gpcar(1:ndime,inode,igaus)) * gpsha(jnode,igaus) * fact0
                 fact2 = dot_product(gpadv(1:ndime,igaus),gpcar(1:ndime,jnode,igaus)) * gpsha(inode,igaus) * fact0

                 elauu11(inode,jnode) = elauu11(inode,jnode) - fact1 - fact2                      ! Auu_xx
                 elauu22(inode,jnode) = elauu22(inode,jnode) - fact1 - fact2                      ! Auu_yy
                 elauu33(inode,jnode) = elauu33(inode,jnode) - fact1 - fact2                      ! Auu_zz
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
              gpveo(1:ndime) = 0.0_rp  ! velocity at itime 
              do inode =1, pnode
                 gpveo(1:ndime) = gpveo(1:ndime) + elvel(1:ndime, inode, itime)* gpsha(inode, igaus)
              end do
              fact0 = 0.5_rp*(gpden(igaus) - densi(igaus, itime))*pabdf_nsi(itime)*dtinv_loc*gpvol(igaus)
              do inode =1, pnode
                 elrbu(1:ndime, inode) = elrbu(1:ndime, inode) &
                                       + fact0*gpsha(inode,igaus)*gpveo(1:ndime)
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
     ! p1vec * rmom2
     !
     !yor! [Originally loop on line 577] Reordered loops to solve indirections problems 
     if (ndime == 2) then
        do igaus = 1,pgaus
           do jnode = 1,pnode
              do inode = 1,pnode
                 elauu11(inode,jnode) = elauu11(inode,jnode) + p1vec(inode,igaus)*rmom2(1,1,jnode,igaus)
                 elauu21(inode,jnode) = elauu21(inode,jnode) + p1vec(inode,igaus)*rmom2(2,1,jnode,igaus)
                 elauu12(inode,jnode) = elauu12(inode,jnode) + p1vec(inode,igaus)*rmom2(1,2,jnode,igaus)
                 elauu22(inode,jnode) = elauu22(inode,jnode) + p1vec(inode,igaus)*rmom2(2,2,jnode,igaus)            
              end do
           end do
        end do
     else !ndime==3
        do igaus = 1,pgaus
           factvec1(1:pnode) = p1vec(1:pnode,igaus)
           do jnode = 1,pnode
              do inode = 1,pnode
                 elauu11(inode,jnode) = elauu11(inode,jnode) + factvec1(inode)*rmom2(1,1,jnode,igaus)
                 elauu21(inode,jnode) = elauu21(inode,jnode) + factvec1(inode)*rmom2(2,1,jnode,igaus)
                 elauu31(inode,jnode) = elauu31(inode,jnode) + factvec1(inode)*rmom2(3,1,jnode,igaus)
                 elauu12(inode,jnode) = elauu12(inode,jnode) + factvec1(inode)*rmom2(1,2,jnode,igaus)
                 elauu22(inode,jnode) = elauu22(inode,jnode) + factvec1(inode)*rmom2(2,2,jnode,igaus)
                 elauu32(inode,jnode) = elauu32(inode,jnode) + factvec1(inode)*rmom2(3,2,jnode,igaus)
                 elauu13(inode,jnode) = elauu13(inode,jnode) + factvec1(inode)*rmom2(1,3,jnode,igaus)
                 elauu23(inode,jnode) = elauu23(inode,jnode) + factvec1(inode)*rmom2(2,3,jnode,igaus)
                 elauu33(inode,jnode) = elauu33(inode,jnode) + factvec1(inode)*rmom2(3,3,jnode,igaus)
              end do
           end do
        end do
     end if

     if( kfl_p1ve2_nsi /= 0 ) then
        !
        ! p1ve2 * rmomu
        ! p1ve2 * rmom2
        !
        !yor! optim : [originally line 597] reverted loops, refactored 
        if(ndime==2)then   
           do igaus=1,pgaus
              do jnode=1,pnode
                 do inode=1,pnode
                    elauu11(inode,jnode)=elauu11(inode,jnode)+p1ve2(1,1,inode,igaus)*rmomu(jnode,igaus)
                    elauu21(inode,jnode)=elauu21(inode,jnode)+p1ve2(2,1,inode,igaus)*rmomu(jnode,igaus)
                    elauu12(inode,jnode)=elauu12(inode,jnode)+p1ve2(1,2,inode,igaus)*rmomu(jnode,igaus)
                    elauu22(inode,jnode)=elauu22(inode,jnode)+p1ve2(2,2,inode,igaus)*rmomu(jnode,igaus)
                    do kdime = 1,ndime
                       elauu11(inode,jnode)=elauu11(inode,jnode)+p1ve2(1,kdime,inode,igaus)*rmom2(kdime,1,jnode,igaus)
                       elauu21(inode,jnode)=elauu21(inode,jnode)+p1ve2(2,kdime,inode,igaus)*rmom2(kdime,1,jnode,igaus)
                       elauu12(inode,jnode)=elauu12(inode,jnode)+p1ve2(1,kdime,inode,igaus)*rmom2(kdime,2,jnode,igaus)
                       elauu22(inode,jnode)=elauu22(inode,jnode)+p1ve2(2,kdime,inode,igaus)*rmom2(kdime,2,jnode,igaus)
                    end do
                 end do
              end do
           end do
        else ! ndime==3
           do igaus=1,pgaus
              do jnode=1,pnode
                 do inode=1,pnode
                    elauu11(inode,jnode)=elauu11(inode,jnode)+p1ve2(1,1,inode,igaus)*rmomu(jnode,igaus)
                    elauu21(inode,jnode)=elauu21(inode,jnode)+p1ve2(2,1,inode,igaus)*rmomu(jnode,igaus)
                    elauu31(inode,jnode)=elauu31(inode,jnode)+p1ve2(3,1,inode,igaus)*rmomu(jnode,igaus)
                    elauu12(inode,jnode)=elauu12(inode,jnode)+p1ve2(1,2,inode,igaus)*rmomu(jnode,igaus)
                    elauu22(inode,jnode)=elauu22(inode,jnode)+p1ve2(2,2,inode,igaus)*rmomu(jnode,igaus)
                    elauu32(inode,jnode)=elauu32(inode,jnode)+p1ve2(3,2,inode,igaus)*rmomu(jnode,igaus)
                    elauu13(inode,jnode)=elauu13(inode,jnode)+p1ve2(1,3,inode,igaus)*rmomu(jnode,igaus)
                    elauu23(inode,jnode)=elauu23(inode,jnode)+p1ve2(2,3,inode,igaus)*rmomu(jnode,igaus)
                    elauu33(inode,jnode)=elauu33(inode,jnode)+p1ve2(3,3,inode,igaus)*rmomu(jnode,igaus)
                    do kdime = 1,ndime
                       elauu11(inode,jnode)=elauu11(inode,jnode)+p1ve2(1,kdime,inode,igaus)*rmom2(kdime,1,jnode,igaus)
                       elauu21(inode,jnode)=elauu21(inode,jnode)+p1ve2(2,kdime,inode,igaus)*rmom2(kdime,1,jnode,igaus)
                       elauu31(inode,jnode)=elauu31(inode,jnode)+p1ve2(3,kdime,inode,igaus)*rmom2(kdime,1,jnode,igaus)
                       elauu12(inode,jnode)=elauu12(inode,jnode)+p1ve2(1,kdime,inode,igaus)*rmom2(kdime,2,jnode,igaus)
                       elauu22(inode,jnode)=elauu22(inode,jnode)+p1ve2(2,kdime,inode,igaus)*rmom2(kdime,2,jnode,igaus)
                       elauu32(inode,jnode)=elauu32(inode,jnode)+p1ve2(3,kdime,inode,igaus)*rmom2(kdime,2,jnode,igaus)
                       elauu13(inode,jnode)=elauu13(inode,jnode)+p1ve2(1,kdime,inode,igaus)*rmom2(kdime,3,jnode,igaus)
                       elauu23(inode,jnode)=elauu23(inode,jnode)+p1ve2(2,kdime,inode,igaus)*rmom2(kdime,3,jnode,igaus)
                       elauu33(inode,jnode)=elauu33(inode,jnode)+p1ve2(3,kdime,inode,igaus)*rmom2(kdime,3,jnode,igaus)
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
           fact0 = bemol_nsi * gpden(igaus) * gpvol(igaus)
           fact5 = fact0 * ( gpgve(1,1,igaus) + gpgve(2,2,igaus) )
           do inode = 1,pnode
              do jnode = 1,pnode              
                 !
                 ! fact2 = - b * ( v , rho * (a.grad) u )
                 ! fact3 = - b * ( u , rho * (a.grad) v ) 
                 ! fact4 = - b * ( rho * div(a) u , v ) 
                 !
                 fact2 = ( gpadv(1,igaus) * gpcar(1,inode,igaus) + gpadv(2,igaus) * gpcar(2,inode,igaus) ) * gpsha(jnode,igaus)
                 fact3 = ( gpadv(1,igaus) * gpcar(1,jnode,igaus) + gpadv(2,igaus) * gpcar(2,jnode,igaus) ) * gpsha(inode,igaus)
                 fact4 = fact5 * gpsha(inode,igaus) * gpsha(jnode,igaus) 
                 fact2 = fact0 * fact2
                 fact3 = fact0 * fact3

                 elauu11(inode,jnode) = elauu11(inode,jnode) - fact2 - fact3 - fact4
                 elauu22(inode,jnode) = elauu22(inode,jnode) - fact2 - fact3 - fact4
              end do
           end do
        end do
     else
        do igaus = 1,pgaus
           fact0 = bemol_nsi * gpden(igaus) * gpvol(igaus)
           fact5 = fact0 * ( gpgve(1,1,igaus) + gpgve(2,2,igaus) + gpgve(3,3,igaus))
           do inode = 1,pnode
              do jnode = 1,pnode
                 !
                 ! fact2 = - b * ( v , rho * (a.grad) u )
                 ! fact3 = - b * ( u , rho * (a.grad) v ) 
                 ! fact4 = - b * ( rho * div(a) u , v ) 
                 !
                 fact2 = ( gpadv(1,igaus) * gpcar(1,inode,igaus) + gpadv(2,igaus) * gpcar(2,inode,igaus) &
                      + gpadv(3,igaus) * gpcar(3,inode,igaus)) * gpsha(jnode,igaus)
                 fact3 = ( gpadv(1,igaus) * gpcar(1,jnode,igaus) + gpadv(2,igaus) * gpcar(2,jnode,igaus) &
                      + gpadv(3,igaus) * gpcar(3,jnode,igaus)) * gpsha(inode,igaus)
                 fact4 = fact5 * gpsha(inode,igaus) * gpsha(jnode,igaus)
                 fact2 = fact0 * fact2
                 fact3 = fact0 * fact3

                 elauu11(inode,jnode) = elauu11(inode,jnode) - fact2 - fact3 - fact4
                 elauu22(inode,jnode) = elauu22(inode,jnode) - fact2 - fact3 - fact4
                 elauu33(inode,jnode) = elauu33(inode,jnode) - fact2 - fact3 - fact4
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
              fact0 = gpvol(igaus) * gpden(igaus) * gpsha(inode,igaus)
              do jnode = 1,pnode
                 fact1 = fact0 * gpsha(jnode,igaus)
                 ! Auu_1i
                 elauu11(inode,jnode) = elauu11(inode,jnode) + fact1 * gpgve(1,1,igaus) ! rho * ux * dux^i-1/dx
                 elauu12(inode,jnode) = elauu12(inode,jnode) + fact1 * gpgve(2,1,igaus) ! rho * uy * dux^i-1/dy
                 ! Auu_2i
                 elauu21(inode,jnode) = elauu21(inode,jnode) + fact1 * gpgve(1,2,igaus) ! rho * ux * duy^i-1/dx
                 elauu22(inode,jnode) = elauu22(inode,jnode) + fact1 * gpgve(2,2,igaus) ! rho * uy * duy^i-1/dy
              end do
           end do
        end do
     else
        do igaus = 1,pgaus
           do inode = 1,pnode
              fact0 = gpvol(igaus) * gpden(igaus) * gpsha(inode,igaus)
              do jnode = 1,pnode
                 fact1 = fact0 * gpsha(jnode,igaus)
                 ! Auu_1i
                 elauu11(inode,jnode) = elauu11(inode,jnode) + fact1 * gpgve(1,1,igaus) ! rho * ux * dux^i-1/dx
                 elauu12(inode,jnode) = elauu12(inode,jnode) + fact1 * gpgve(2,1,igaus) ! rho * uy * dux^i-1/dy
                 elauu13(inode,jnode) = elauu13(inode,jnode) + fact1 * gpgve(3,1,igaus) ! rho * uz * dux^i-1/dz
                 ! Auu_2i
                 elauu21(inode,jnode) = elauu21(inode,jnode) + fact1 * gpgve(1,2,igaus) ! rho * ux * duy^i-1/dx
                 elauu22(inode,jnode) = elauu22(inode,jnode) + fact1 * gpgve(2,2,igaus) ! rho * uy * duy^i-1/dy
                 elauu23(inode,jnode) = elauu23(inode,jnode) + fact1 * gpgve(3,2,igaus) ! rho * uz * duy^i-1/dz
                 ! Auu_3i
                 elauu31(inode,jnode) = elauu31(inode,jnode) + fact1 * gpgve(1,3,igaus) ! rho * ux * duz^i-1/dx
                 elauu32(inode,jnode) = elauu32(inode,jnode) + fact1 * gpgve(2,3,igaus) ! rho * uy * duz^i-1/dy
                 elauu33(inode,jnode) = elauu33(inode,jnode) + fact1 * gpgve(3,3,igaus) ! rho * uz * duz^i-1/dz

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
        do igaus =1, pgaus
           do inode =1, pnode
              fact0 = gpvol(igaus) * gpden(igaus) * gpsha(inode,igaus) * dtinv_mod * pabdf_nsi(1)
              elauu11(inode,inode) =  elauu11(inode,inode) + fact0
              elauu22(inode,inode) =  elauu22(inode,inode) + fact0
              do jnode =1, pnode !yor! for this case, better write a second loop and revert inode and jnode
                 elauu11(inode,jnode) =  elauu11(inode,jnode) - fact0*gpsha(jnode, igaus)
                 elauu22(inode,jnode) =  elauu22(inode,jnode) - fact0*gpsha(jnode, igaus)
              end do
           end do
        end do
        do itime =2, nbdfp_nsi ! RHS
           do igaus =1, pgaus
              gpveo(1:2) = 0.0_rp
              do inode =1, pnode
                 gpveo(1:2) = gpveo(1:2) + elvel(1:2, inode, itime)* gpsha(inode, igaus)
              end do
              do inode =1, pnode
                 fact0 = - gpvol(igaus) * gpden(igaus) * gpsha(inode,igaus) * dtinv_mod * pabdf_nsi(itime)
                 elrbu(1:2,inode)   =  elrbu(1:2,inode)   - fact0*gpveo(1:2)
                 elrbu(1:2,inode)   =  elrbu(1:2,inode)   + fact0*elvel(1:2, inode, itime)
              end do
           end do
        end do
     else ! ndime
        do igaus =1, pgaus
           do inode =1, pnode
              fact0 = gpvol(igaus) * gpden(igaus) * gpsha(inode,igaus) * dtinv_mod * pabdf_nsi(1)
              elauu11(inode,inode) =  elauu11(inode,inode) + fact0
              elauu22(inode,inode) =  elauu22(inode,inode) + fact0
              elauu33(inode,inode) =  elauu33(inode,inode) + fact0
              do jnode =1, pnode !yor! for this case, better write a second loop and revert inode and jnode
                 elauu11(inode,jnode) =  elauu11(inode,jnode) - fact0*gpsha(jnode, igaus)
                 elauu22(inode,jnode) =  elauu22(inode,jnode) - fact0*gpsha(jnode, igaus)
                 elauu33(inode,jnode) =  elauu33(inode,jnode) - fact0*gpsha(jnode, igaus)
              end do
           end do
        end do
        do itime =2, nbdfp_nsi ! RHS
           do igaus =1, pgaus
              gpveo(1:3) = 0.0_rp
              do inode =1, pnode
                 gpveo(1:3) = gpveo(1:3) + elvel(1:3, inode, itime)* gpsha(inode, igaus)
              end do
              do inode =1, pnode
                 fact0 = - gpvol(igaus) * gpden(igaus) * gpsha(inode,igaus) * dtinv_mod * pabdf_nsi(itime)
                 elrbu(1:3,inode)   =  elrbu(1:3,inode)   - fact0*gpveo(1:3)
                 elrbu(1:3,inode)   =  elrbu(1:3,inode)   + fact0*elvel(1:3, inode, itime)
              end do
           end do
        end do
     end if
  else if( kfl_lumped == 2 ) then 
     !
     ! No time term have been added up to now: add Galerkin term
     !
     do igaus = 1,pgaus
        fact0 = gpvol(igaus) * gpden(igaus) * dtinv_mod
        if( ndime == 2 ) then
           do inode = 1, pnode
              elauu11(inode,inode) = elauu11(inode,inode) + fact0 * gpsha(inode,igaus)
              elauu22(inode,inode) = elauu22(inode,inode) + fact0 * gpsha(inode,igaus)
           end do
        else ! ndime==3
           do inode = 1, pnode
              elauu11(inode,inode) = elauu11(inode,inode) + fact0 * gpsha(inode,igaus)
              elauu22(inode,inode) = elauu22(inode,inode) + fact0 * gpsha(inode,igaus)
              elauu33(inode,inode) = elauu33(inode,inode) + fact0 * gpsha(inode,igaus)
           end do
        end if
        do inode = 1, pnode
           do idime = 1,ndime   
              elrbu(idime,inode) = elrbu(idime,inode) + fact0 * gpsha(inode,igaus) * elvel(idime,inode,2)
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
           fact0 = gpvol(igaus) * gpden(igaus) * gpsha(inode,igaus) * dtinv_mod * real(fact_duatss-1,rp)
           ellum(inode) = ellum(inode) + fact0 
        end do
     end do
  end if

  !yor! elauu final matrix assembly
  if( ndime == 2 ) then
     do jnode = 1,pnode
        do inode = 1,pnode
           elauu(2*inode-1,2*jnode-1)=elauu11(inode,jnode)
           elauu(2*inode  ,2*jnode-1)=elauu21(inode,jnode)
           elauu(2*inode-1,2*jnode  )=elauu12(inode,jnode)
           elauu(2*inode  ,2*jnode  )=elauu22(inode,jnode)
        end do
     end do
  else ! ndime==3
     do jnode = 1,pnode
        do inode = 1,pnode
           elauu(3*inode-2,3*jnode-2)=elauu11(inode,jnode)
           elauu(3*inode-1,3*jnode-2)=elauu21(inode,jnode)
           elauu(3*inode  ,3*jnode-2)=elauu31(inode,jnode)
           elauu(3*inode-2,3*jnode-1)=elauu12(inode,jnode)
           elauu(3*inode-1,3*jnode-1)=elauu22(inode,jnode)
           elauu(3*inode  ,3*jnode-1)=elauu32(inode,jnode)
           elauu(3*inode-2,3*jnode  )=elauu13(inode,jnode)
           elauu(3*inode-1,3*jnode  )=elauu23(inode,jnode)
           elauu(3*inode  ,3*jnode  )=elauu33(inode,jnode)
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
           fact2 = gpvol(igaus)*gpsha(jnode,igaus)
           do inode = 1,pnode
              fact1 = -gpvol(igaus)*gpsha(inode,igaus)+p1vec(inode,igaus)
              idof1 = (inode-1)*ndime
              do idime = 1,ndime
                 idofv = idof1 + idime
                 elaup(idofv,jnode) = elaup(idofv,jnode) - fact2 * gpcar(idime,inode,igaus) &
                      &                                  + fact1 * gpcar(idime,jnode,igaus)
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
                 elaup(idofv,jnode) = elaup(idofv,jnode) + p1vec(inode,igaus) * gpcar(idime,jnode,igaus)
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
                 elaup(idofv,jnode) = elaup(idofv,jnode) + &
                      dot_product(gpcar(1:ndime,jnode,igaus),p1ve2(idime,1:ndime,inode,igaus))
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
  if (kfl_regim_nsi == 3 .and. kfl_confi_nsi == 1) then ! Penalization of pressure, never imposse pressure
    do igaus =1, pgaus
       penal = 1.0e-4_rp*gpden(igaus)/gpvis(igaus)
       do inode =1, pnode
          fact0 = penal*p2sca(inode,igaus)
          do jnode =inode +1, pnode
             fact1 = fact0*gpsha(jnode,igaus)
             elapp(inode, jnode) = elapp(inode, jnode) + fact1
             elapp(jnode, inode) = elapp(jnode, inode) + fact1
          end do
          elapp (inode, inode) = elapp (inode, inode) + fact0*gpsha(inode, igaus)
          elrbp (inode) = elrbp (inode) + fact0*gppre(igaus)
       end do
    end do
  end if

#ifdef matiaslma
  
  if (ndime==2) then
     do igaus =1, pgaus
        do inode =1,pnode
           p2vec(1,inode,igaus)= p2vec(1,inode,igaus)*gpden(igaus)
           p2vec(2,inode,igaus)= p2vec(2,inode,igaus)*gpden(igaus)
           p2sca(inode, igaus) = p2sca(inode, igaus)*gpden(igaus)
        end do
        fact0 = gpden(igaus)* gpvol(igaus)
        do inode =1, pnode           
           idof1 = 2*inode-1
           idof2 = idof1+1
           fact1 = gpsha(inode,igaus)* fact0
           do jnode =1, pnode
              elapu(jnode, idof1) = elapu(jnode, idof1) -fact1             * gpcar(1,jnode,igaus)  &
                   &                                  + rmomu(inode,igaus) * p2vec(1,jnode,igaus)
              elapu(jnode, idof2) = elapu(jnode, idof2) -fact1             * gpcar(2,jnode,igaus)  &
                   &                                  + rmomu(inode,igaus) * p2vec(2,jnode,igaus)
           end do
        end do
     end do
  else !ndime =3
     do igaus =1, pgaus
        do inode =1,pnode
           p2vec(1:ndime,inode,igaus) = p2vec(1:ndime,inode,igaus)*gpden(igaus)
           p2sca(inode, igaus)        = p2sca(inode, igaus) *gpden(igaus)
        end do
        fact0 = gpden(igaus)* gpvol(igaus)
        do inode =1, pnode
           idof1 = 3*inode - 2 
           idof2 = idof1+1
           idof3 = idof2+1
           fact1 = gpsha(inode,igaus)* fact0
           do jnode =1, pnode 
              elapu(jnode, idof1) = elapu(jnode, idof1) - fact1              * gpcar(1,jnode,igaus) &
                   &                                    + rmomu(inode,igaus) * p2vec(1,jnode,igaus)
              elapu(jnode, idof2) = elapu(jnode, idof2) - fact1              * gpcar(2,jnode,igaus) &
                   &                                    + rmomu(inode,igaus) * p2vec(2,jnode,igaus)
              elapu(jnode, idof3) = elapu(jnode, idof3) - fact1              * gpcar(3,jnode,igaus) &
                   &                                    + rmomu(inode,igaus) * p2vec(3,jnode,igaus)
           end do
        end do
     end do
  end if

#else
  do igaus = 1,pgaus
     if (kfl_sgsli_nsi == 0_ip) then
        gpvel(1:ndime) = 0.0_rp
        do inode =1, pnode
           gpvel(1:ndime) = gpvel(1:ndime) + elvel(1:ndime, inode, 1)* gpsha(inode, igaus)
        end do
        do inode = 1,pnode
           auxva = 0.0_rp
           do idime =1, ndime
              auxva = auxva  + gpden(igaus)*sgsli(idime)*gpcar(idime, inode, igaus)
           end do
           do jnode = 1,pnode
              do idime = 1,ndime
                 idof1 = (inode-1)*ndime+idime
                 elapu(jnode,idof1) = elapu(jnode,idof1) + rcont(idime,inode,igaus)     * p2sca(jnode,igaus)&
                      &                                  + (rmomu(inode,igaus) - auxva) * p2vec(idime,jnode,igaus)
              end do
           end do
        end do
     else
        do inode = 1,pnode
           do jnode = 1,pnode
              do idime = 1,ndime
                 idof1 = (inode-1)*ndime+idime
                 elapu(jnode,idof1) = elapu(jnode,idof1) + rcont(idime,inode,igaus) * p2sca(jnode,igaus)&
                      &                                  + rmomu(inode,igaus)       * p2vec(idime,jnode,igaus)
              end do
           end do
        end do
     end if
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
                    elapu(inode,jdofv) = elapu(inode,jdofv) + rmom2(kdime,jdime,jnode,igaus) * p2vec(kdime,inode,igaus)
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
        do jnode = inode+1,pnode
           fact1 = dot_product(p2vec(1:ndime,jnode,igaus),gpcar(1:ndime,inode,igaus))
           elapp(jnode,inode) = elapp(jnode,inode) + fact1
           elapp(inode,jnode) = elapp(inode,jnode) + fact1
        end do
        fact1 = dot_product(p2vec(1:ndime,inode,igaus),gpcar(1:ndime,inode,igaus))
        elapp(inode,inode) = elapp(inode,inode) + fact1
     end do
  end do
  !
  ! Penalization
  !
  do igaus = 1,pgaus
     fact1 = penal_nsi * gpvol(igaus)
     do inode = 1,pnode
        elapp(inode,inode) = elapp(inode,inode) + fact1 * gpsha(inode,igaus)
        elrbp(inode)       = elrbp(inode)       + fact1 * gpsha(inode,igaus) * elpre(inode,1) 
     end do
  end do

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
        gprhs_sgs(1,igaus) = gprhs_sgs(1,igaus) + taupr(igaus) * gpvep(1,igaus)
        gprhs_sgs(2,igaus) = gprhs_sgs(2,igaus) + taupr(igaus) * gpvep(2,igaus)
        gprhh(1,igaus)     = gprhs(1,igaus)     + gprhs_sgs(1,igaus)
        gprhh(2,igaus)     = gprhs(2,igaus)     + gprhs_sgs(2,igaus)
        do inode = 1,pnode
           elrbu(1,inode) = elrbu(1,inode) + p1vec(inode,igaus)   * gprhh(1,igaus)
           elrbu(2,inode) = elrbu(2,inode) + p1vec(inode,igaus)   * gprhh(2,igaus)
           elrbp(inode)   = elrbp(inode)   + p2vec(1,inode,igaus) * gprhh(1,igaus) &
                &                          + p2vec(2,inode,igaus) * gprhh(2,igaus) 
        end do
     end do
  else
     do igaus = 1,pgaus
        gprhs_sgs(1,igaus) = gprhs_sgs(1,igaus) + taupr(igaus) * gpvep(1,igaus)
        gprhs_sgs(2,igaus) = gprhs_sgs(2,igaus) + taupr(igaus) * gpvep(2,igaus)
        gprhs_sgs(3,igaus) = gprhs_sgs(3,igaus) + taupr(igaus) * gpvep(3,igaus)
        gprhh(1, igaus)    = gprhs(1,igaus)     + gprhs_sgs(1,igaus)
        gprhh(2, igaus)    = gprhs(2,igaus)     + gprhs_sgs(2,igaus)
        gprhh(3, igaus)    = gprhs(3,igaus)     + gprhs_sgs(3,igaus)
        do inode = 1,pnode
           elrbu(1,inode) = elrbu(1,inode) + p1vec(inode,igaus)   * gprhh(1,igaus)
           elrbu(2,inode) = elrbu(2,inode) + p1vec(inode,igaus)   * gprhh(2,igaus)
           elrbu(3,inode) = elrbu(3,inode) + p1vec(inode,igaus)   * gprhh(3,igaus)
           elrbp(inode)   = elrbp(inode)   + p2vec(1,inode,igaus) * gprhh(1,igaus) &
                &                          + p2vec(2,inode,igaus) * gprhh(2,igaus) &
                &                          + p2vec(3,inode,igaus) * gprhh(3,igaus)
        end do
     end do
  end if
  if( kfl_p1ve2_nsi /= 0 ) then
     do igaus = 1,pgaus
        do inode = 1,pnode
           do idime = 1,ndime
              do kdime = 1,ndime
                 elrbu(idime,inode) = elrbu(idime,inode) &
                      + p1ve2(idime,kdime,inode,igaus) * gprhh(kdime,igaus)
              end do
           end do
        end do
     end do
  end if
  !
  ! low Mach regime: we add to the right hand side the residue of the continuity eq
  !
  do igaus = 1,pgaus
     fact0 = gprhc(igaus) * gpvol(igaus) 
     fact1 = gpsp2(igaus) * gprh2(igaus) * gpvol(igaus) 
     do inode = 1,pnode
        elrbu(1:ndime,inode) = elrbu(1:ndime,inode) + fact1 * gpcar(1:ndime,inode,igaus) ! ( rhs, tau2' * div(v) )
        elrbp(inode)         = elrbp(inode)         + fact0 * gpsha(inode,igaus)         ! ( rhs , q)
     end do
  end do
  !
  ! Newton-Raphson
  !
  if( kfl_linea_nsi == 2 ) then
     if( ndime == 2 ) then
        do igaus = 1,pgaus
           fact0 =  gpden(igaus) * gpvol(igaus) 
           do inode = 1,pnode
              fact1 = fact0 * gpsha(inode,igaus)
              elrbu(1,inode) = elrbu(1,inode) + fact1 * ( gpadv(1,igaus) * gpgve(1,1,igaus) + gpadv(2,igaus) * gpgve(2,1,igaus) )
              elrbu(2,inode) = elrbu(2,inode) + fact1 * ( gpadv(1,igaus) * gpgve(1,2,igaus) + gpadv(2,igaus) * gpgve(2,2,igaus) )
           end do
        end do
     else
        do igaus = 1,pgaus
           fact0 =  gpden(igaus) * gpvol(igaus) 
           do inode = 1,pnode
              fact1 = fact0 * gpsha(inode,igaus)
              elrbu(1,inode) = elrbu(1,inode) & 
                   + fact1 * ( gpadv(1,igaus) * gpgve(1,1,igaus) + gpadv(2,igaus) * gpgve(2,1,igaus) &
                   + gpadv(3,igaus) * gpgve(3,1,igaus) )
              elrbu(2,inode) = elrbu(2,inode) &
                   + fact1 * ( gpadv(1,igaus) * gpgve(1,2,igaus) + gpadv(2,igaus) * gpgve(2,2,igaus) &
                   + gpadv(3,igaus) * gpgve(3,2,igaus) )
              elrbu(3,inode) = elrbu(3,inode) &
                   + fact1 * ( gpadv(1,igaus) * gpgve(1,3,igaus) + gpadv(2,igaus) * gpgve(2,3,igaus) &
                   + gpadv(3,igaus) * gpgve(3,3,igaus) )
           end do
        end do
     end if
  end if
  !
  ! Projection: ( L*(v) , P )
  !
  if( kfl_stabi_nsi == 1 ) then

     if( ndime == 2 ) then

        do igaus = 1, pgaus
           do inode = 1,pnode
              elrbu(1,inode) = elrbu(1,inode) - gprhs_sgs(1, igaus)*gpsha(inode, igaus)*gpvol(igaus)
              elrbu(2,inode) = elrbu(2,inode) - gprhs_sgs(2, igaus)*gpsha(inode, igaus)*gpvol(igaus)
           end do
        end do
!!$        do igaus = 1,pgaus 
!!$           fact1 = gpden(igaus) * dtsgs_nsi
!!$           fact2 = gpvol(igaus) / gpst1(igaus)
!!$           factx = fact1 * gpsgs(1,igaus,2) + gpvep(1,igaus) / gpst1(igaus) ! rho/dt * u'^n + tau^-1 P
!!$           facty = fact1 * gpsgs(2,igaus,2) + gpvep(2,igaus) / gpst1(igaus) 
!!$           facx1 = fact2 * gpvep(1,igaus)                                   ! tau^-1 P
!!$           facy1 = fact2 * gpvep(2,igaus)
!!$           do inode = 1,pnode
!!$              elrbu(1,inode) = elrbu(1,inode) + p1vec(inode,igaus)   * factx - gpsha(inode,igaus) * facx1
!!$              elrbu(2,inode) = elrbu(2,inode) + p1vec(inode,igaus)   * facty - gpsha(inode,igaus) * facy1
!!$              elrbp(inode)   = elrbp(inode)   + p2vec(1,inode,igaus) * factx + p2vec(2,inode,igaus) * facty 
!!$           end do
!!$        end do

     else
        do igaus = 1, pgaus
           do inode = 1,pnode
              elrbu(1,inode) = elrbu(1,inode) - gprhs_sgs(1, igaus)*gpsha(inode, igaus)*gpvol(igaus)
              elrbu(2,inode) = elrbu(2,inode) - gprhs_sgs(2, igaus)*gpsha(inode, igaus)*gpvol(igaus)
              elrbu(3,inode) = elrbu(3,inode) - gprhs_sgs(3, igaus)*gpsha(inode, igaus)*gpvol(igaus)
           end do
        end do
!!$        do igaus = 1,pgaus 
!!$           fact1 = gpden(igaus) * dtsgs_nsi
!!$           fact2 = gpvol(igaus) / gpst1(igaus)
!!$           factx = fact1 * gpsgs(1,igaus,2) + gpvep(1,igaus) / gpst1(igaus) ! rho/dt * u'^n + tau^-1 P
!!$           facty = fact1 * gpsgs(2,igaus,2) + gpvep(2,igaus) / gpst1(igaus) 
!!$           factz = fact1 * gpsgs(3,igaus,2) + gpvep(3,igaus) / gpst1(igaus) 
!!$           facx1 = fact2 * gpvep(1,igaus)                                   ! tau^-1 P
!!$           facy1 = fact2 * gpvep(2,igaus)
!!$           facz1 = fact2 * gpvep(3,igaus)
!!$           do inode = 1,pnode
!!$              elrbu(1,inode) = elrbu(1,inode) + p1vec(inode,igaus)   * factx - gpsha(inode,igaus) * facx1
!!$              elrbu(2,inode) = elrbu(2,inode) + p1vec(inode,igaus)   * facty - gpsha(inode,igaus) * facy1
!!$              elrbu(3,inode) = elrbu(3,inode) + p1vec(inode,igaus)   * factz - gpsha(inode,igaus) * facz1
!!$              elrbp(inode)   = elrbp(inode)   + p2vec(1,inode,igaus) * factx &
!!$                   &                          + p2vec(2,inode,igaus) * facty &
!!$                   &                          + p2vec(3,inode,igaus) * factz 
!!$           end do
!!$        end do

     end if
  end if


  if( pbubl == 1 ) then

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
        p2sca_bub(igaus)         = gptt2(igaus) * gpvol(igaus) * gpsha_bub(igaus)
        p2vec_bub(1:ndime,igaus) = gpsp1(igaus) * gpvol(igaus) * gpcar_bub(1:ndime,igaus)
     end do
     !
     ! Momentum equations
     !
     if( kfl_press_nsi == 1 ) then
        do igaus = 1,pgaus
           fact2 = gpvol(igaus)*gpsha_bub(igaus)
           do inode = 1,pnode
              fact1 = -gpvol(igaus)*gpsha(inode,igaus)+p1vec(inode,igaus)
              idof1 = (inode-1)*ndime
              do idime = 1,ndime
                 idofv = idof1 + idime
                 elauq(idofv,1) = elauq(idofv,1) - fact2 * gpcar(idime,inode,igaus) &
                      &                          + fact1 * gpcar_bub(idime,igaus)
              end do
           end do
        end do
     else
        do igaus = 1,pgaus
           do inode = 1,pnode
              idof1 = (inode-1)*ndime
              do idime = 1,ndime
                 idofv = idof1 + idime
                 elauq(idofv,1) = elauq(idofv,1) + p1vec(inode,igaus) * gpcar_bub(idime,igaus)
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
                 elauq(idofv,1) = elauq(idofv,1) + &
                      dot_product(gpcar_bub(1:ndime,igaus),p1ve2(idime,1:ndime,inode,igaus))
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
              elaqu(1,idof1) = elaqu(1,idof1) + rcont(idime,inode,igaus) * p2sca_bub(igaus) &
                   &                          + rmomu(inode,igaus)       * p2vec_bub(idime,igaus)
           end do
        end do
     end do
     if( kfl_rmom2_nsi /= 0 ) then
        do igaus = 1,pgaus
           do jnode = 1,pnode
              do jdime = 1,ndime
                 jdofv = (jnode-1)*ndime+jdime
                 do kdime = 1,ndime
                    elaqu(1,jdofv) = elaqu(1,jdofv) + rmom2(kdime,jdime,jnode,igaus) * p2vec_bub(kdime,igaus)
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
           elapq(inode,1) = elapq(inode,1) + dot_product(p2vec(1:ndime,inode,igaus),gpcar_bub(1:ndime,igaus))
           elaqp(1,inode) = elaqp(1,inode) + dot_product(p2vec_bub(1:ndime,igaus)  ,gpcar(1:ndime,inode,igaus))
        end do
        elaqq(1,1) = elaqq(1,1) + dot_product(p2vec_bub(1:ndime,igaus),gpcar_bub(1:ndime,igaus))
        do idime = 1,ndime
           gprhs_sgs(idime,igaus) = gprhs_sgs(idime,igaus) + taupr(igaus) * gpvep(idime,igaus)
           gprhh(idime,igaus)     = gprhs(idime,igaus)     + gprhs_sgs(idime,igaus) 
           elrbq(1)               = elrbq(1)               + p2vec_bub(idime,igaus) * gprhh(idime,igaus)
        end do
        elrbq(1) = elrbq(1) + gpvol(igaus) * gprhc(igaus) * gpsha_bub(igaus)
     end do
     !
     ! Penalization
     !
     do igaus = 1,pgaus
        elaqq(1,1) = elaqq(1,1) + penal_nsi * gpvol(igaus) * gpsha_bub(igaus) 
        elrbq(1)   = elrbq(1)   + penal_nsi * gpvol(igaus) * gpsha_bub(igaus) * elbub 
     end do
     !
     ! Eliminate bubble
     !
     ! Auu <= Auu - Auq Aqq^-1 Aqu
     ! Aup <= Aup - Auq Aqq^-1 Aqp
     ! bu  <= bu  - Auq Aqq^-1 bq
     !
     ! Apu <= Apu - Apq Aqq^-1 Aqu
     ! App <= App - Apq Aqq^-1 Aqp
     ! bp  <= bp  - Apq Aqq^-1 bq
     ! 
     call maths_schur_complement(pevat,pevat, 1_ip,elauu,elauq,elaqq,elaqu)
     call maths_schur_complement(pevat,pnode, 1_ip,elaup,elauq,elaqq,elaqp)
     call maths_schur_complement(pevat, 1_ip, 1_ip,elrbu,elauq,elaqq,elrbq)
     
     call maths_schur_complement(pnode,pevat, 1_ip,elapu,elapq,elaqq,elaqu)
     call maths_schur_complement(pnode,pnode, 1_ip,elapp,elapq,elaqq,elaqp)
     call maths_schur_complement(pnode, 1_ip, 1_ip,elrbp,elapq,elaqq,elrbq)

  end if

end subroutine nsi_elmmat
