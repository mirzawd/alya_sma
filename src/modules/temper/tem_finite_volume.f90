!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tem_finite_volume

  use def_kintyp
  use def_master
  use def_domain
  use def_kermod
  use def_temper
  use mod_ker_proper, only : ker_proper
  implicit none
  integer(ip) :: ielem,iface,jface,ifacg,jelem
  integer(ip) :: iiz,ijz,jjz,jiz,iboun
  real(rp)    :: Area,phi,value_dirichlet,value_neumann
  real(rp)    :: gpden(2),gpcon(2),gpsph(2)
  real(rp)    :: phi_ij,phi_ji,Fij,Fji,rhocp
  real(rp)    :: val_ci,val_di,val_ui, corr_i
  real(rp)    :: val_cj,val_dj,val_uj, corr_j
  integer(ip)    ::  aci, adi
  integer(ip)    ::  acj, adj
  real(rp)    :: Li,Lj,rf_i,rf_j,fun_i,fun_j
  real(rp)    :: lo_ij,lo_ji,ho_ij,ho_ji
  real(rp)    :: localVelocity(3)
  real(rp)    :: Vol
  real(rp)    :: gradientTempe(1:ndime,nelem) ! no hay una variable global que sea el gradiente de la temp??
  integer(ip), parameter:: fv_grad_method_gauss = 1, fv_grad_method_ls = 2
  integer(ip), parameter:: fv_advec_minmod = 1, fv_advec_vanleer = 2, fv_advec_superbee = 3, fv_advec_upwind = 4, fv_advec_cds = 5
  integer(ip), parameter:: fv_advec_solv_full_implicit = 1,fv_advec_solv_deferred = 2
  integer(ip) :: my_grad = fv_grad_method_gauss, my_advec_scheme = fv_advec_minmod, my_advec_solv = fv_advec_solv_full_implicit

  ! parametros cutres
  localVelocity(1) = 1.0_rp*(cos(0.8726646259971648_rp))
  localVelocity(2) = 1.0_rp*(sin(0.8726646259971648_rp))
  localVelocity(3) = 0.0_rp

  !localVelocity = 0.0_rp
  !localVelocity(1) = 1.0_rp

  ! Gradients evaluation, conservative gauss (consistent with formulation
  ! turbulent flow) or ls (optimal for unstructured
  ! meshes)

  select case( my_grad )

  case ( fv_grad_method_gauss )
     call tem_gradient_gauss(tempe,gradientTempe)
  case (fv_grad_method_ls)
     call tem_gradient_ls(tempe,gradientTempe)
  case default
     call tem_gradient_gauss(tempe,gradientTempe) ! is the one more robust
  end select

  do ifacg = 1,nfacg
     ielem = lfacg(1,ifacg) 
     jelem = lfacg(2,ifacg)
     iface = lfacg(3,ifacg)
     jface = lfacg(4,ifacg)
     iiz   = fv_face_graph(1,ifacg)
     ijz   = fv_face_graph(2,ifacg)
     jjz   = fv_face_graph(3,ifacg)
     jiz   = fv_face_graph(4,ifacg)
     !
     ! Diffusive flux for equation IELEM and JELEM:
     ! 
     ! FOR IELEM: 
     ! A_II   = sum_iface - k * Area_f / phi_IJ
     ! A_IJ   = k   * Area_f / phi_IJ
     ! phi_IJ = | (center_j,center_i) . n_face |
     ! FOR JELEM:
     ! A_JJ   = sum_iface - k * Area_f / phi_JI
     ! A_JI   = k   * Area_f / phi_JI
     ! phi_JI = | (center_i,center_j) . n_face |
     !
     ! Convective flux:
     !
     ! Fij and Fji for now analytical velocity function now diagonal flow benchmark
     ! Convection: TVD scheme from Denner et al., Journal of Computational Physics 298 (2015) 466–479, some parts changed
     !
     if( jelem > 0 ) then
        !
        ! Properties
        !
        call ker_proper('CONDU','COG  ',1_ip,ielem,gpcon)
        call ker_proper('DENSI','COG  ',1_ip,ielem,gpden) ! todo: interpolate at the face with harmonic mean 
        call ker_proper('SPHEA','COG  ',1_ip,ielem,gpsph)
        rhocp = gpden(1) * gpsph(1)
        !
        ! Area and phi 
        !           
        Area   = fv_face_area(ifacg) 
        phi_ij = abs(dot_product(fv_center_vector(1:ndime,ijz),fv_face_normal(1:ndime,ifacg)))
        phi_ji = abs(dot_product(fv_center_vector(1:ndime,jiz),fv_face_normal(1:ndime,ifacg)))
        Fij    = dot_product(localVelocity(1:ndime),fv_face_normal(1:ndime,ifacg)) * Area * rhocp * fv_face_orientation(ijz)
        Fji    = -Fij 

        !
        ! estas distancias se tendran que precomputar y guardar, ahora solo para
        ! ver si tiene sentido el tema
        !
        Li = (    sqrt(dot_product(fv_center_face(1:ndime,iface,ielem), fv_center_face(1:ndime,iface,ielem)) ) + &
             &    sqrt(dot_product(fv_center_face(1:ndime,jface,jelem), fv_center_face(1:ndime,jface,jelem)) )   &
             &) / sqrt(dot_product(fv_center_face(1:ndime,iface,ielem), fv_center_face(1:ndime,iface,ielem))) 
        Lj = (    sqrt(dot_product(fv_center_face(1:ndime,jface,jelem), fv_center_face(1:ndime,jface,jelem)) ) + &
             &    sqrt(dot_product(fv_center_face(1:ndime,iface,ielem), fv_center_face(1:ndime,iface,ielem)) )   &
             &) / sqrt(dot_product(fv_center_face(1:ndime,jface,jelem), fv_center_face(1:ndime,jface,jelem))) 

        if( Fij > 0.0_rp ) then
           val_ci = tempe(ielem,1)
           val_di = tempe(jelem,1)
           val_ui = tempe(jelem,1) - 2.0_rp*dot_product( gradientTempe(1:ndime, ielem) , fv_center_vector(1:ndime,ijz)) 
           corr_i = Li*dot_product( gradientTempe(1:ndime, ielem) , fv_center_face(1:ndime,iface,ielem)) / (val_di - val_ci)

           val_cj = val_ci
           val_dj = val_di
           val_uj = val_ui
           corr_j = corr_i

           aci = iiz
           adi = ijz
           acj = jiz
           adj = jjz
        else 
           val_ci = tempe(jelem,1)
           val_di = tempe(ielem,1)
           val_ui = tempe(ielem,1) - 2.0_rp*dot_product( gradientTempe(1:ndime, jelem) , fv_center_vector(1:ndime,jiz))
           corr_j = Lj*dot_product( gradientTempe(1:ndime, jelem) , fv_center_face(1:ndime,jface,jelem)) / (val_dj - val_cj)

           val_cj = val_ci
           val_dj = val_di
           val_uj = val_ui
           corr_i = corr_j
           acj = jjz
           adj = jiz
           aci = ijz
           adi = iiz
        end if

        rf_i = 0.0_rp
        rf_i = 0.0_rp
        if( abs(val_di-val_ci) > 1e-10_rp ) then
           rf_i = ((val_ci - val_ui) / (val_di - val_ci))! + corr_i
        end if
        if( abs(val_dj-val_cj) > 1e-10_rp ) then 
           rf_j = ((val_cj - val_uj) / (val_dj - val_cj))! + corr_j
        endif

        select case ( my_advec_scheme )

        case ( fv_advec_minmod )        
           !
           ! MINMOD
           !
           fun_i = max(0.0_rp,min(1.0_rp, rf_i))
           fun_j = max(0.0_rp,min(1.0_rp, rf_j))

        case ( fv_advec_vanleer )
           !
           ! VAN LERR
           !
           fun_i = (0.5_rp*Li*rf_i + 0.5_rp*Li*abs(rf_i)) / (Li - 1.0_rp + abs(rf_i))
           fun_j = (0.5_rp*Lj*rf_j + 0.5_rp*Lj*abs(rf_j)) / (Lj - 1.0_rp + abs(rf_j))

        case ( fv_advec_superbee )
           !
           ! SUPERBEE
           !
           fun_i = max(0.0_rp, min(1.0_rp, Li * rf_i), min(rf_i, Li) )
           fun_j = max(0.0_rp, min(1.0_rp, Lj * rf_j), min(rf_j, Lj) )
 
        case ( fv_advec_upwind )
           !
           ! UPWIND
           !
           fun_i = 0.0_rp
           fun_j = 0.0_rp

        case ( fv_advec_cds )
           !
           ! CDS Unstructured
           !
           fun_i = 1.0_rp
           fun_j = 1.0_rp

        case default 
           !
           ! UPWIND as default
           !
           fun_i = 0.0_rp
           fun_j = 0.0_rp

        end select

        if( (my_advec_scheme /= fv_advec_cds) .and. (my_advec_scheme /= fv_advec_upwind) ) then
           fun_i = max(max(0.0_rp, min(rf_i, 1.0_rp)) ,fun_i) !TVD bounding 
           fun_i = min(fun_i, max(0.0_rp, min(rf_i, Li), min(Li*rf_i, 1.0_rp)))
           fun_j = max(max(0.0_rp, min(rf_j, 1.0_rp)) ,fun_j) 
           fun_j = min(fun_j, max(0.0_rp, min(rf_j, Lj), min(Lj*rf_j, 1.0_rp)))
        end if
        
        lo_ij  = Fij * val_ci
        lo_ji  = Fji * val_cj
        ! blended with a central second order unstructured interpolation
        ho_ij  = Fij * ( (1.0_rp - (fun_i / Li) )*val_ci  + (fun_i/Li)*val_di )
        ho_ji  = Fji * ( (1.0_rp - (fun_j / Lj) )*val_cj  + (fun_j/Lj)*val_dj )
        
        !
        ! Assemble diffusion  
        !
        amatr(iiz) = amatr(iiz) + gpcon(1) * Area / phi_ij
        amatr(ijz) = amatr(ijz) - gpcon(1) * Area / phi_ij
        amatr(jjz) = amatr(jjz) + gpcon(1) * Area / phi_ji
        amatr(jiz) = amatr(jiz) - gpcon(1) * Area / phi_ji
        !
        ! Assemble  convection 
        !
        
        select case ( my_advec_solv )

        case ( fv_advec_solv_full_implicit ) !only make sense for upwind and cds, or a blending between them, util we use a NR solver

          amatr(aci)   = amatr(aci)   + Fij * (1.0_rp - (fun_i / Li))
          amatr(adi)   = amatr(adi)   + Fij * (fun_i/Li)                 

          amatr(acj)   = amatr(acj)   + Fji * (1.0_rp - (fun_j / Lj))              
          amatr(adj)   = amatr(adj)   + Fji * (fun_j/Lj)

        case ( fv_advec_solv_deferred )

          amatr(iiz)   = amatr(iiz)   + max(Fij,0.0_rp)              ! Low order
          amatr(ijz)   = amatr(ijz)   + min(Fij,0.0_rp)                 
          rhsid(ielem) = rhsid(ielem) + (lo_ij - ho_ij)      ! Deferred correction

          amatr(jjz)   = amatr(jjz)   + max(Fji,0.0_rp)              ! Low order
          amatr(jiz)   = amatr(jiz)   + min(Fji,0.0_rp)
          rhsid(jelem) = rhsid(jelem) + (lo_ji - ho_ji)      ! Deferred correction

        end select
     end if

  end do
  !
  ! Time derivative
  !
  if( kfl_timei_tem /= 0 ) then
     do ielem = 1,nelem
        call ker_proper('DENSI','COG  ',1_ip,ielem,gpden) 
        call ker_proper('SPHEA','COG  ',1_ip,ielem,gpsph)
        rhocp        = gpden(1) * gpsph(1)
        Vol          = fv_cell_volume(ielem) 
        iiz          = fv_graph_diag(ielem)
        amatr(iiz)   = amatr(iiz)   + dtinv_tem * rhocp * Vol
        rhsid(ielem) = rhsid(ielem) + dtinv_tem * rhocp * Vol * tempe(ielem,3)
     end do
  end if
  !
  ! Boundary conditions
  !
  do iboun = 1,nboun

     ifacg = fv_face_boundary(iboun)
     ielem = lfacg(1,ifacg)

     if( kfl_fixbo_tem(iboun) == 1 ) then
        !
        ! Impose Dirichlet
        !
        ! Take half of distance
        ! phi = | (face,center_i) . n_face |
        !
        call ker_proper('CONDU','COG  ',1_ip,ielem,gpcon)
        call ker_proper('DENSI','COG  ',1_ip,ielem,gpden) 
        call ker_proper('SPHEA','COG  ',1_ip,ielem,gpsph) 

        rhocp           = gpden(1)*gpsph(1)
        value_dirichlet = bvnat_tem(1,iboun,1)
        Area            = fv_face_area(ifacg)
        iface           = lfacg(3,ifacg)
        phi             = dot_product(fv_center_face(1:ndime,iface,ielem),fv_face_normal(1:ndime,ifacg))
        iiz             = fv_face_graph(1,ifacg)
        Fij             = dot_product(localVelocity(1:ndime), fv_face_normal(1:ndime,ifacg)) * Area * rhocp * fv_face_orientation(iiz)

        amatr(iiz)      = amatr(iiz)   + gpcon(1) * Area / abs(phi)

        rhsid(ielem)    = rhsid(ielem) + gpcon(1) * Area / abs(phi) * value_dirichlet
        rhsid(ielem)    = rhsid(ielem) -  Fij * value_dirichlet ! strong bcc

     else if( kfl_fixbo_tem(iboun) == 2 ) then 
        !
        ! Impose Neumann
        !
        call ker_proper('DENSI','COG  ',1_ip,ielem,gpden) 
        call ker_proper('SPHEA','COG  ',1_ip,ielem,gpsph) 

        rhocp         = gpden(1)*gpsph(1)
        value_neumann = bvnat_tem(1,iboun,1)
        Area          = fv_face_area(ifacg)
        iface         = lfacg(3,ifacg)
        iiz           = fv_face_graph(1,ifacg)
        Fij           = dot_product(localVelocity(1:ndime),fv_face_normal(1:ndime,ifacg)) * Area * rhocp * fv_face_orientation(iiz)

        amatr(iiz)    = amatr(iiz)   + Fij 
        rhsid(ielem)  = rhsid(ielem) + Area * value_neumann 

     end if
  end do

end subroutine tem_finite_volume

! evaluate gradients of tempe with a global 2Order skew gauss method not the best in
! tetras but fast to implement, in turbulence consitent with the formulation to
! be utilized, thus, cancelation errors are expected
subroutine tem_gradient_gauss (xvalu,xgrad)

  use def_kintyp
  use def_master
  use def_domain
  use def_kermod
  use def_temper
  implicit none
  real(rp), intent(in)  :: xvalu(nelem)
  real(rp), intent(out) :: xgrad(1:ndime,nelem)
  integer(ip)           :: ielem,iface,ifacg,jelem, jface
  integer(ip)           :: iiz,ijz,jiz
  real(rp)              :: Area
  real(rp)              :: invVol_i, invVol_j,avtem,xfact(ndime)

  
  xgrad = 0.0_rp 

  do ifacg = 1,nfacg
     ielem    = lfacg(1,ifacg) 
     jelem    = lfacg(2,ifacg)
     iface    = lfacg(3,ifacg)
     jface    = lfacg(4,ifacg)
     if( jelem > 0 ) then
        invVol_i                     = 1.0_rp / fv_cell_volume(ielem)  
        invVol_j                     = 1.0_rp / fv_cell_volume(jelem)  
        Area                         = fv_face_area(ifacg) 
        avtem                        = 0.5_rp * ( xvalu(jelem) + xvalu(ielem) )
        ijz                          = fv_face_graph(2,ifacg)
        jiz                          = fv_face_graph(4,ifacg)
        xfact                        = Area * fv_face_normal(1:ndime,ifacg) * fv_face_orientation(ijz)
        xgrad(1:ndime,ielem) = xgrad(1:ndime, ielem) + avtem * invVol_i * xfact
        xgrad(1:ndime,jelem) = xgrad(1:ndime, jelem) - avtem * invVol_j * xfact
     else 
        iiz      = fv_face_graph(1,ifacg)
        invVol_i = 1.0_rp / fv_cell_volume(ielem)  
        Area     = fv_face_area(ifacg)

        xgrad(1:ndime,ielem) = xgrad(1:ndime, ielem) &
                + xvalu(ielem) * ( Area * invVol_i ) * fv_face_normal(1:ndime,ifacg) * fv_face_orientation(iiz)
     end if
  end do

end subroutine tem_gradient_gauss

! evaluate gradients of tempe with LS method based on face connectivity for
! standard numerical methods is enough, for interface problems better utilize vertex
! conectivity
! A first order taylor extrapolation is minimized in order to get the gradient
! at the centroids, tempe at the centroids is supose to be know, look to my phd!
! La parte de la matriz a y su inversa deberian ser un preproceso pero no se
! como ahora mismo... fortran.. :(
! por vertices seria igual pero con ditancias entre centroides que compartenb
! vertices, en general no val la pena

subroutine tem_gradient_ls (xvalu,xgrad)

  use def_kintyp
  use def_master
  use def_domain
  use def_kermod
  use def_temper
  use mod_maths
  implicit none
  real(rp), intent(in)  :: xvalu(nelem)
  real(rp), intent(out) :: xgrad(1:ndime,nelem)
  integer(ip)           :: ielem,ifacg,jelem,idime,jdime, iface
  integer(ip)           :: iiz,ijz,jjz,jiz
  real(rp)              :: a(1:ndime,1:ndime, nelem)
  real(rp)              :: b(1:ndime, nelem)
  real(rp)              :: w

  a     = 0.0_rp
  b     = 0.0_rp
  xgrad = 0.0_rp

  do ifacg = 1,nfacg
     ielem = lfacg(1,ifacg) 
     jelem = lfacg(2,ifacg)
     iiz   = fv_face_graph(1,ifacg)
     ijz   = fv_face_graph(2,ifacg)
     jjz   = fv_face_graph(3,ifacg)
     jiz   = fv_face_graph(4,ifacg)
     iface = lfacg(3,ifacg)

     if( jelem > 0 ) then
         w = 1.0_rp/fv_center_distance(ijz)
         w = w*w
         do idime =1,ndime
            do jdime =1,ndime
               a(idime,jdime,ielem) = a(idime,jdime,ielem) + w * fv_center_vector(jdime,ijz) * fv_center_vector(idime,ijz)
               a(idime,jdime,jelem) = a(idime,jdime,jelem) + w * fv_center_vector(jdime,jiz) * fv_center_vector(idime,jiz)
            end do
         end do
         b(1:ndime,ielem) = b(1:ndime, ielem) + w * (xvalu(jelem) - xvalu(ielem)) * fv_center_vector(1:ndime,ijz)
         b(1:ndime,jelem) = b(1:ndime, jelem) + w * (xvalu(ielem) - xvalu(jelem)) * fv_center_vector(1:ndime,jiz)
     else
       w = 1.0_rp/sqrt(dot_product(fv_center_face(1:ndime,iface,ielem),fv_center_face(1:ndime,iface,ielem)))
       w = w*w
       do idime =1,ndime
          do jdime =1,ndime
             a(idime,jdime,ielem) = a(idime,jdime,ielem) + w * fv_center_face(jdime,iface,ielem) * fv_center_face(idime,iface,ielem)
          end do
       end do
     end if
  end do

  do ielem = 1,nelem
     call maths_invert_matrix(ndime, a(1:ndime,1:ndime,ielem))
     do idime =1,ndime
        do jdime =1,ndime
           xgrad(idime,ielem) = xgrad(idime, ielem) + a(idime, jdime, ielem) *  b(jdime, ielem)
        end do
     end do
  end do

end subroutine tem_gradient_ls
