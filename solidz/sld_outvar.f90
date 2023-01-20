!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    sld_outvar.f90
!> @author  Mariano Vazquez
!> @date    August, 2006
!>          - Subroutine written
!> @brief   Field Output variables for post-process of results
!>
!> @details
!>
!>          \verbatim
!>          Post-processing vectors (gevec):
!>          --------------------------------
!>          Vectors shall be stored in gevec with ndime components
!>
!>          Post-processing scalars at element/point level (gesca):
!>          -------------------------------------------------------
!>          Scalars shall be stored in gesca
!>
!>          Post-processing symmetric tensors (gevec):
!>          ------------------------------------------
!>          Tensors are stored in gevec with a voigt notation. They are
!>          defined as "TENSO" in sld_inivar. Finally, this subroutine
!>          spread the tensor components in scalars, dumped in:
!>          "...XX", "...YY", "...XY", ...
!>
!>          Post-processing scalars at gauss point level (ger3p):
!>          -----------------------------------------------------
!>          Scalars shall be stored in ger3p
!>
!>          \endverbatim
!>
!> @}
!-----------------------------------------------------------------------

subroutine sld_outvar(ivari,imesh)

  use def_kintyp,          only : ip,rp
  use def_parame,          only : zero
  use def_master,          only : coupling
  use def_master,          only : gesca,gevec,ger3p
  use def_master,          only : INOTEMPTY,kfl_paral,leinv_loc,lninv_loc
  use def_master,          only : ITER_K, ITER_K_STATE
  use def_master,          only : nfacg,npoin_type,rhsid,lfacg,displ,therm
  use def_master,          only : postp,solve,cutim,ittim,forcf,solve_sol
  use def_domain,          only : ndime,nelem,npoin,ngaus,nboun,neset,nbset,nnset
  use def_master,          only : gdepo
  use def_master,          only : donna_gp
  use def_domain,          only : nnode
  use def_domain,          only : lmate,lnods,ltype,lnoch,lnodb,ltypb,lpoty,lelch
  use def_domain,          only : leset,lbset,lnset
  use def_domain,          only : vmass,coord
  use def_domain,          only : kfl_icodb, kfl_codbo
  use def_domain,          only : kfl_icodn, kfl_codno
  use def_elmtyp,          only : NODE_CONTACT_SOLID, HEX08, QUA04, PEN06
  use def_solidz,          only : bvess_sld, bvnat_sld
  use def_elmgeo,          only : element_type
  use def_solidz,          only : veloc_sld, accel_sld
  use def_solidz,          only : frxid_sld, fexte_sld
  use def_solidz,          only : lawst_sld,lawch_sld
  use def_solidz,          only : svegm_sld,nsvar_sld
  use def_solidz,          only : kfl_foten_sld
  use def_solidz,          only : caust_sld, cause_sld, caunp_sld, caulo_sld, caunn_sld
  use def_solidz,          only : green_sld, lepsi_sld, eprin_sld
  use def_solidz,          only : epsee_sld, celen_sld
  use def_solidz,          only : dxfem_sld,vxfem_sld,axfem_sld,lnenr_sld
  use def_solidz,          only : kfl_fixno_sld, kfl_fixbo_sld, kfl_fixrs_sld
  use def_solidz,          only : jacrot_du_dq_sld
  use def_solidz,          only : lcrkf_sld,isoch_sld
  use def_solidz,          only : nvoig_sld
  use def_solidz,          only : sigei_sld,rorig_sld
  use def_solidz,          only : kfl_fiber_sld,fibeg_sld,seqvm_sld
  use def_solidz,          only : kfl_cohes_sld, kfl_damag_sld
  use def_solidz,          only : axis1_sld, axis2_sld, axis3_sld, orien_sld, rmate_sld
  use mod_sld_csys,        only : sld_csys_midsurface_element_normal, sld_csys_rotuni
  use def_solidz,          only : kfl_conta_sld, fcont_sld
  use def_solidz,          only : kfl_rigid_sld
  use def_solidz,          only : srpro_sld
  use def_solidz,          only : kfl_donna_sld
  use def_solidz,          only : water_sld, enden_sld
  use mod_communications,  only : PAR_INTERFACE_NODE_EXCHANGE
#if defined COMMDOM && COMMDOM == 2
  use mod_plepp_pdn_contact, only : PLEPP_CPLNG
#endif
  use mod_sld_fe2
  use mod_eccoupling,      only : calcium_ecc, kfl_exmsld_ecc
  use mod_outvar,          only : outvar
  use mod_biofibers,       only : biofib_point_nodal_fibers, kfl_biofibers
  use mod_sld_stress_model_comput, only : SM_voigt_to_tensor_second
  use mod_maths_basic,     only : maths_MULT_VxMxV, maths_MULT_MxV, maths_normalize_vector
  use mod_projec,          only : projec_elements_to_nodes
  use mod_arrays,          only : arrays_name
  use mod_sld_strain,      only : strain_basis

  implicit none

  integer(ip), intent(in) :: ivari   !< checking variable ivari
  integer(ip), intent(in) :: imesh   !< Mesh to postprocess on
  integer(ip)             :: pnode,pgaus,pelty,pblty,pnodb,hnode
  integer(ip)             :: ivoig,icont,ipoin,jdime
  integer(ip)             :: jpoin,idime,iline,inodf,iface,ibopo
  integer(ip)             :: ielem,ifacg,inode,ielty,idofn,igaus,isdva
  integer(ip)             :: iboun,inodb
  real(rp)                :: xauxi(2),anpoi,andis,angul,dummr,foref
  real(rp)                :: u(3),dudx(9),d2udx2(27)
  real(rp)                :: norvec(ndime)
  real(rp)                :: svegm_aux
  real(rp)                :: bvess_aux(ndime)               
  integer(ip)             :: e
  integer(ip)             :: cost
  logical                 :: non_linear, converged
  real(rp)                :: matrix(ndime, ndime)

  !
  ! Define postprocess variable
  !
  select case ( arrays_name(ivari) )

  case( 'DISPL' )
     !
     ! DISPL: Displacement
     !
     if( INOTEMPTY ) gevec => displ(:,:,ITER_K)

  case( 'VELOC' )
     !
     ! VELOC: Velocity
     !
     if( INOTEMPTY ) gevec => veloc_sld(:,:,ITER_K)

  case( 'ACCEL' )
     !
     ! ACCEL: Acceleration
     !
     if( INOTEMPTY ) gevec => accel_sld(:,:,ITER_K)

  case( 'SIGMA' )
     !
     ! SIGMA: Cauchy stress tensor
     !
     if( INOTEMPTY ) then
        kfl_foten_sld = kfl_foten_sld + 1_ip
        if(kfl_foten_sld == 1_ip ) call sld_calculate_pp_tensors
        call memgen(0_ip, nvoig_sld, npoin)
        gevec(1:nvoig_sld,1:npoin) = caust_sld(1:nvoig_sld,1:npoin)
     end if

  case( 'EPSIL' )
     !
     ! EPSIL: Strain tensor (Green)
     !
     if( INOTEMPTY ) then
        kfl_foten_sld = kfl_foten_sld + 1_ip
        if( kfl_foten_sld == 1_ip ) call sld_calculate_pp_tensors
        call memgen(0_ip, nvoig_sld, npoin)
        gevec(1:nvoig_sld,1:npoin) = green_sld(1:nvoig_sld,1:npoin)
     end if

  case( 'LNEPS' )
     !
     ! LNEPS: Logarithmic strain tensor
     !
     if( INOTEMPTY ) then
        kfl_foten_sld = kfl_foten_sld + 1_ip
        if( kfl_foten_sld == 1_ip ) call sld_calculate_pp_tensors
        call memgen(0_ip, nvoig_sld, npoin)
        gevec(1:nvoig_sld,1:npoin) = lepsi_sld(1:nvoig_sld,1:npoin)

     end if

  case( 'SEQVM' )
     !
     ! SEQVM: Von Mises Stress
     !
     if( INOTEMPTY ) then
        kfl_foten_sld = kfl_foten_sld + 1_ip
        if (kfl_foten_sld == 1_ip) call sld_calculate_pp_tensors
        gesca => seqvm_sld
     end if

  case( 'FIBER' )
     !
     ! BFIBL: Bio-fibers longitudinal-fiber direction
     !
     call biofib_point_nodal_fibers( gevec, 'LONGITUDINAL', 'CURRENT' )

  case( 'BVESS' )
     !
     ! BVESS: Essential boundary conditions
     !
     if( INOTEMPTY ) then
        call memgen(0_ip,ndime,npoin)
        do ipoin = 1,npoin
           ibopo = lpoty(ipoin)
           bvess_aux(1:ndime) = bvess_sld(1:ndime,ipoin,ITER_K)
           if( ibopo > 0 ) then
              if( kfl_fixrs_sld(ipoin) /= 0_ip .and. kfl_fixno_sld(1,ipoin) /=3_ip) then
                 ! Local --> Global
                 call sld_csys_rotuni(2_ip,ndime,ipoin,bvess_aux(1:ndime))
              end if
           end if
           gevec(1:ndime,ipoin) = bvess_aux(1:ndime)
        end do
        
     end if

  case( 'FIBEL' )
     !
     ! Fibres on elements
     !
     call sld_elmope(4_ip)
     ger3p => fibeg_sld

  case( 'DXFEM' )
     !
     ! XFEM Displacement
     !
     gevec => dxfem_sld(:,:,1)

  case( 'VXFEM' )
     !
     ! XFEM Velocity
     !
     gevec => vxfem_sld(:,:,1)

  case( 'AXFEM' )
     !
     ! XFEM Acceleration
     !
     gevec => axfem_sld(:,:,1)

  case( 'LINEL' )
     !
     ! LINEL: Linelets of preconditioner CG
     !
     if( INOTEMPTY ) then
        do ipoin = 1,npoin
           rhsid(ipoin) = 0.0_rp
        end do
        icont = 0
        do iline = 1,solve(1)%nline
           icont = icont + 1
           do ipoin = solve(1)%lline(iline),solve(1)%lline(iline+1)-1
              jpoin = solve(1)%lrenup(ipoin)
              rhsid(jpoin) = real(icont,rp)
           end do
        end do
        gesca => rhsid
     end if

  case( 'NSIGN' )
     !
     ! Stress tensor
     !
     if (INOTEMPTY) then
        kfl_foten_sld= kfl_foten_sld + 1
        if (kfl_foten_sld == 1) then   ! fotens computes caust, green and lepsi, so do it only the first kfote
           call sld_calculate_pp_tensors
        end if
     end if

     gesca => caunn_sld

  case( 'XYROT' )
     !
     ! Rotation on plane XY
     !
     if( INOTEMPTY ) then
        call memgen(zero,npoin,zero)
        do ipoin = 1,npoin

           anpoi= 0.0_rp  ! coordinates angle in xy
           if (abs(coord(1,ipoin)) > 0) anpoi= atan(coord(2,ipoin)/coord(1,ipoin))
           andis= 0.0_rp  ! displacement angle in xy
           if (abs(displ(1,ipoin,1)) > 0) anpoi= atan(displ(2,ipoin,1)/displ(1,ipoin,1))
           !!acaaaaaaaaa
           angul= andis - anpoi     ! displacement angle with respect to coord

           gesca(ipoin)= angul

        end do
     end if

  case( 'ZXROT' )
     !
     ! Rotation on plane XZ
     !
     if( INOTEMPTY .and. ndime==3) then
        call memgen(zero,npoin,zero)
        do ipoin = 1,npoin
           xauxi(1)= displ(3,ipoin,1) - rorig_sld(3)
           xauxi(2)= displ(2,ipoin,1) - rorig_sld(2)
           gesca(ipoin)= 0.0_rp
           if (xauxi(1) > 0.0_rp) gesca(ipoin)= atan(xauxi(2)/xauxi(1))
        end do
     end if

  case( 'DCOHE' )
     !
     ! DCOHE: Damage variable for cohesive elements
     !
     if( kfl_cohes_sld == 0_ip ) return
     if( INOTEMPTY ) then
        !
        ! Averaged value from GP (Scalar per element)
        !
        call memgen(0_ip,nelem,0_ip)
        do ielem = 1,nelem
           pgaus = ngaus(ltype(ielem))
           pnode = nnode(ltype(ielem))
           if ( lawch_sld(lmate(ielem)) == 904_ip .or. &
                lawch_sld(lmate(ielem)) == 905_ip  ) then
              pgaus = pnode/2
              gesca(ielem) = sum(svegm_sld(ielem)%a(1,1:pgaus,ITER_K_STATE))/real(pgaus,rp)
           end if
        end do
     end if

  case( 'LNENR' )
     !
     ! LNENR_sld
     !
     if( INOTEMPTY ) then
        call memgen(zero,npoin,zero)
        do ipoin = 1,npoin
           gesca(ipoin) = real(lnenr_sld(ipoin),rp)
        end do
     end if

  case( 'ERROR' )
     !
     ! Error
     !
     if( INOTEMPTY) then
        call memgen(zero,ndime,npoin)
        do ipoin = 1,npoin
           call sld_exacso(1_ip,coord(1,ipoin),u,dudx,d2udx2,dummr,dummr)
           do idime = 1,ndime
              gevec(idime,ipoin) = displ(idime,ipoin,1)-u(idime)
           end do
        end do
     end if

  case( 'CRACK' )
     !
     ! Cracked nodes
     !
     if( INOTEMPTY ) then
        call memgen(zero,npoin,zero)
        do ifacg = 1,nfacg
           if( lcrkf_sld(ifacg) /= 0 ) then
              iface = lfacg(3,ifacg)
              ielem = lfacg(1,ifacg)
              ielty = ltype(ielem)
              do inodf = 1,element_type(ielty) % node_faces(iface)
                 inode = element_type(ielty) % list_faces(inodf,iface)
                 ipoin = lnods(inode,ielem)
                 gesca(ipoin) = 1.0_rp
              end do
           end if
        end do
        call pararr('SLX',NPOIN_TYPE,npoin,gesca)
        do ipoin = 1,npoin
           gesca(ipoin) = min(1.0_rp,gesca(ipoin))
        end do
     end if

  case( 'GROUP' )
     !
     ! GROUP: Groups for solvers
     !
     if( INOTEMPTY ) then
        call memgen(zero,npoin,zero)
        do ipoin = 1,npoin
           gesca(ipoin) = real( solve_sol(1) % lgrou(ipoin) ,rp)
        end do
     end if

  case( 'SIGEI' )
     !
     ! Sigma largest eigenvalue
     !
     if (INOTEMPTY) then
        kfl_foten_sld= kfl_foten_sld + 1
        if (kfl_foten_sld == 1) then   ! fotens computes caust, green and lepsi, so do it only the first kfote
           call sld_calculate_pp_tensors
        end if
     end if

     gesca => sigei_sld

  case( 'SDV1E' )
     !
     ! SDV1E: d1 (sm152)
     !
     if( INOTEMPTY ) then
        !
        ! Averaged value from GP (Scalar per element)
        !
        call memgen(zero,nelem,zero)
        do ielem = 1,nelem
           pgaus = ngaus(ltype(ielem))
           svegm_aux = 0.0_rp
           if ( lawst_sld(lmate(ielem)) == 152_ip ) then
              do igaus = 1,pgaus
                 svegm_aux = svegm_aux + svegm_sld(ielem)%a(4,igaus,2)
              end do
           end if
           gesca(ielem) = svegm_aux/real(pgaus,rp)
        end do
     end if

  case( 'SDV2E' )
     !
     ! SDV2E: dG (sm152)
     !
     if( INOTEMPTY ) then
        !
        ! Averaged value from GP (Scalar per element)
        !
        call memgen(zero,nelem,zero)
        do ielem = 1,nelem
           pgaus = ngaus(ltype(ielem))
           svegm_aux = 0.0_rp
           if ( lawst_sld(lmate(ielem)) == 152_ip ) then
              do igaus = 1,pgaus
                 svegm_aux = svegm_aux + svegm_sld(ielem)%a(5,igaus,2)
              end do
           end if
           gesca(ielem) = svegm_aux/real(pgaus,rp)
        end do
     end if

  case( 'SDV3E' )
     !
     ! SDV3E: dK (sm152)
     !
     if( INOTEMPTY ) then
        !
        ! Averaged value from GP (Scalar per element)
        !
        call memgen(zero,nelem,zero)
        do ielem = 1,nelem
           pgaus = ngaus(ltype(ielem))
           svegm_aux = 0.0_rp
           if ( lawst_sld(lmate(ielem)) == 152_ip ) then
              do igaus = 1,pgaus
                 svegm_aux = svegm_aux + svegm_sld(ielem)%a(6,igaus,2)
              end do
           end if
           gesca(ielem) = svegm_aux/real(pgaus,rp)
        end do
     end if

  case( 'SDV4E' )
     !
     ! SDV4E: d6 (sm152)
     !
     if( INOTEMPTY ) then
        !
        ! Averaged value from GP (Scalar per element)
        !
        call memgen(zero,nelem,zero)
        do ielem = 1,nelem
           pgaus = ngaus(ltype(ielem))
           svegm_aux = 0.0_rp
           if ( lawst_sld(lmate(ielem)) == 152_ip ) then
              do igaus = 1,pgaus
                 svegm_aux = svegm_aux + svegm_sld(ielem)%a(7,igaus,2)
              end do
           end if
           gesca(ielem) = svegm_aux/real(pgaus,rp)
        end do
     end if

  case( 'FIXRS' )
     !
     ! FIXRS: Fixity code for local axes
     !
     if( INOTEMPTY ) then
        call memgen(0_ip,npoin,0_ip)
        do ipoin = 1,npoin
           if ( lpoty(ipoin) > 0 ) then
              gesca(ipoin) = real(kfl_fixrs_sld(ipoin),rp)
           else
              gesca(ipoin) = 0.0_rp
           end if
        end do
     end if

  case( 'SRPRO' )
     !
     ! SRPRO: Strenght reduction properties
     !
     if( kfl_damag_sld == 0_ip ) return
     if( INOTEMPTY ) then
        call memgen(0_ip,8_ip,nelem)
        do ielem = 1,nelem
           if( lawst_sld(lmate(ielem)) == 154_ip ) then
              gevec(1:8,ielem) = srpro_sld(1:8,ielem)
           end if
        end do
     end if

  case( 'SDVAR' )
     !
     ! SDVAR: State dependent variables
     !
     if( kfl_damag_sld == 0_ip ) return
     if( INOTEMPTY ) then
        call memgen(0_ip,nsvar_sld,npoin)
        call projec_elements_to_nodes(nsvar_sld,svegm_sld,gevec)
     end if

  case( 'NOSET' )
     !
     ! NOSET: Node sets
     !
     if( nnset < 1 ) return
     if( INOTEMPTY ) then
        call memgen(0_ip,npoin,0_ip)
        gesca(1:npoin) = real(lnset(1:npoin),rp)
     end if

  case( 'NELEM' )
     !
     ! NELEM: Element number
     !
     if( INOTEMPTY ) then
        call memgen(0_ip,nelem,0_ip)
        gesca(1:nelem) = real(leinv_loc(1:nelem),rp)
     end if
     
  case( 'PELCH' )
     !
     ! PELCH: Element characteristic code
     !
     if( INOTEMPTY ) then
        call memgen(0_ip,nelem,0_ip)
        do ielem = 1,nelem
           gesca(ielem) = real(lelch(ielem),rp)
        end do
     end if
     
  case( 'FIXNO' )
     !
     ! FIXNO: Fixity codes for each DoF
     !
     if( INOTEMPTY ) then
        call memgen(zero,ndime,npoin)
        gevec(1:ndime,1:npoin) = real(kfl_fixno_sld(1:ndime,1:npoin),rp)
     end if
     
  case( 'EXAFO' )
     !
     ! Exact force
     !
     if( INOTEMPTY) then
        call memgen(zero,ndime,npoin)
        do idofn = 1,npoin*ndime
           rhsid(idofn) = 0.0_rp
        end do
        call sld_elmope(9_ip)
        call PAR_INTERFACE_NODE_EXCHANGE(ndime,rhsid,'SUM')
        do ipoin = 1,npoin
           idofn = (ipoin-1)*ndime
           do idime = 1,ndime
              idofn = idofn + 1
              gevec(idime,ipoin) = rhsid(idofn) / vmass(ipoin)
           end do
        end do
     end if

  case( 'FORCF' )
     !
     ! Fluid force
     !
     if( INOTEMPTY ) then
        call memgen(zero,ndime,npoin)
        do ipoin = 1,npoin
           if( lnoch(ipoin) ==  NODE_CONTACT_SOLID ) then
              do idime = 1,ndime
                 if( kfl_fixno_sld(idime,ipoin) /= 1 ) then
                    foref = 0.0_rp
                    do jdime = 1,ndime
                       foref = foref + gdepo(idime,jdime,ipoin) * forcf(jdime,ipoin)
                    end do
                    gevec(idime,ipoin) = foref
                 end if
              end do
           end if
        end do
     end if

  case( 'SEGMA' )
     !
     ! Element-wise Cauchy stress
     !
     if( INOTEMPTY ) then
        do ielem = 1,nelem
           pelty = abs(ltype(ielem))
           pgaus = ngaus(pelty)
           cause_sld(ielem) % a(1:ndime,1:pgaus,1) = 0.0_rp
        end do
        call sld_elmope(10_ip)
        ger3p => cause_sld
     end if

  case( 'FRXID' )
     !
     ! FRXID: Nodal Reactions Forces
     !
     if ( INOTEMPTY ) then
        call memgen(zero,ndime,npoin)
        idofn = 0_ip
        do ipoin = 1,npoin
           do idime = 1,ndime
              idofn = idofn + 1
              gevec(idime,ipoin) = -frxid_sld(idofn) ! opposite sign
           end do
        end do
     end if

  case( 'SIGNP' )
     !
     ! Recovered stresses
     !
     if( INOTEMPTY ) then
        do ielem = 1,nelem
           pelty = abs(ltype(ielem))
           pgaus = ngaus(pelty)
           caunp_sld(ielem) % a(1:ndime,1:pgaus,1) = 0.0_rp
        end do
        call sld_elmope(11_ip)
        ger3p => caunp_sld
     end if

  case( 'SIGLO' )
     !
     ! Element-wise Cauchy stress according to the local coordinate system        ( *AQU* )
     !
     if( INOTEMPTY ) then
        do ielem = 1,nelem
           pelty = abs(ltype(ielem))
           pgaus = ngaus(pelty)
           caulo_sld(ielem) % a(:,:,1) = 0.0_rp
        end do
        call sld_elmope(12_ip)
        ger3p => caulo_sld
     end if

  case( 'AXIS1' )
     !
     ! AXIS1: (Material CSYS)
     !
     if( INOTEMPTY ) then
        call memgen(0_ip, ndime, nelem)
        if( kfl_fiber_sld > 3 ) then
           gevec(1:ndime,1:nelem) = axis1_sld(1:ndime,1:nelem)
        else
           gevec(1:ndime,1:nelem) = 0.0_rp
        end if
     end if

  case( 'AXIS2' )
     !
     ! AXIS2: (Material CSYS)
     !
     if( INOTEMPTY ) then
        call memgen(0_ip, ndime, nelem)
        if( kfl_fiber_sld > 3 ) then
           gevec(1:ndime,1:nelem) = axis2_sld(1:ndime,1:nelem)
        else
           gevec(1:ndime,1:nelem) = 0.0_rp
        end if
     end if

  case( 'AXIS3' )
     !
     ! AXIS3: (Material CSYS)
     !
     if( ndime == 2_ip ) return
     if( INOTEMPTY ) then
        call memgen(0_ip, ndime, nelem)
        if( kfl_fiber_sld > 3 ) then
           gevec(1:ndime,1:nelem) = axis3_sld(1:ndime,1:nelem)
        else
           gevec(1:ndime,1:nelem) = 0.0_rp
        end if
     end if

  case( 'ORIEN' )
     !
     ! ORIEN: Orientation angle
     !
     if( INOTEMPTY ) then
        call memgen(0_ip,nelem,0_ip)
        if ( kfl_fiber_sld > 4 ) then
           gesca(1:nelem) = orien_sld(1:nelem)
        else
           gesca(1:nelem) = 0.0_rp
        end if
     end if

  case( 'SREAC' )
     !
     ! SREAC: Solver reaction
     !
     if( INOTEMPTY ) then
        call memgen(0_ip,ndime,npoin)
        if ( solve_sol(1)%kfl_react == 1_ip ) then
           gevec(1:ndime,1:npoin) = solve_sol(1)%reaction(1:ndime,1:npoin)
        end if
     end if
     
  case( 'PMATE' )
     !
     ! PMATE: Material Numbering
     !
     if( INOTEMPTY ) then
        call memgen(0_ip,nelem,0_ip)
        gesca(1:nelem) = real(lmate(1:nelem),rp)
     end if

  case( 'RMAT1' )
     !
     ! RMAT1: Rotated material axis 1
     !
     if( INOTEMPTY ) then
        call memgen(0_ip,ndime,nelem)
        if( kfl_fiber_sld > 3 ) then
           do ielem = 1,nelem
              pgaus = ngaus(ltype(ielem))
              do idime = 1,ndime
                 gevec(idime,ielem) = sum(rmate_sld(ielem)%a(1,idime,1:pgaus))/real(pgaus,rp)
              end do
           end do
        else
           gevec(1:ndime,1:nelem) = 0.0_rp
        end if
     end if
     
  case( 'RMAT2' )
     !
     ! RMAT2: Rotated material axis 2
     !
     if( INOTEMPTY ) then
        call memgen(0_ip,ndime,nelem)
        if( kfl_fiber_sld > 3 ) then
           do ielem = 1,nelem
              pgaus = ngaus(ltype(ielem))
              do idime = 1,ndime
                 gevec(idime,ielem) = sum(rmate_sld(ielem)%a(2,idime,1:pgaus))/real(pgaus,rp)
              end do
           end do
        else
           gevec(1:ndime,1:nelem) = 0.0_rp
        end if
     end if

  case( 'RMAT3' )
     !
     ! RMAT3: Rotated material axis 3
     !
     if( ndime == 2_ip ) return
     if( INOTEMPTY ) then
        call memgen(0_ip,ndime,nelem)
        if( kfl_fiber_sld > 3 ) then
           do ielem = 1,nelem
              pgaus = ngaus(ltype(ielem))
              do idime = 1,ndime
                 gevec(idime,ielem) = sum(rmate_sld(ielem)%a(3,idime,1:pgaus))/real(pgaus,rp)
              end do
           end do
        else
           gevec(1:ndime,1:nelem) = 0.0_rp
        end if
     end if     

  case ( 'DIISO' )
     ! 'DIISO'
     if( INOTEMPTY ) then
        gesca => isoch_sld(:,1)
     endif

  case ( 'VMISO' )
     ! 'VMISO'
     if( INOTEMPTY ) then
        gesca => isoch_sld(:,2)
     endif

  case( 'BOSET' )
     !
     ! BOSET: Boundary set
     !
     if( nbset < 1 ) return
     if( INOTEMPTY ) then
        call memgen(0_ip,npoin,0_ip)
        gesca(1:npoin) = 0.0_rp
        do iboun = 1,nboun
           do inodb = 1,nnode(abs(ltypb(iboun)))
              ipoin = lnodb(inodb,iboun)
              gesca(ipoin) = max(gesca(ipoin),real(lbset(iboun),rp))
           end do
        end do
        call PAR_INTERFACE_NODE_EXCHANGE(gesca,'MAX','IN MY CODE')
     end if

  case( 'ELSET' )
     !
     ! ELSET: Element set
     !
     if( neset < 1 ) return
     if( INOTEMPTY ) then
        call memgen(0_ip,nelem,0_ip)
        gesca(1:nelem) = real(leset(1:nelem),rp)
     end if

  case( 'EPSIS' )
     !
     ! Infinitessimal strain tensor
     !
     if (INOTEMPTY) then
        !green_sld = 0.0_rp
        call sld_elmope(13_ip)
        call memgen(0_ip, 6_ip, npoin)

        do ipoin = 1,npoin
           do ivoig = 1,nvoig_sld
              gevec(ivoig,ipoin) = green_sld(ivoig,ipoin)
           end do
        end do
     end if

  case( 'PK2  ' )
     !
     ! Infinitesimal stress tensor
     !
     if (INOTEMPTY) then
        !caust_sld = 0.0_rp
        call sld_elmope(14_ip)
        call memgen(0_ip, 6_ip, npoin)

        do ipoin = 1,npoin
           do ivoig = 1,nvoig_sld
              gevec(ivoig,ipoin) = caust_sld(ivoig,ipoin)
           end do
        end do
     end if

  case( 'INTLI' )
     !
     ! INTLI: Interior list PLE++ 
     !
     if( INOTEMPTY ) then
        call memgen(0_ip, npoin, 0_ip)
#if defined COMMDOM && COMMDOM == 2
        if( ittim > 0_ip ) then
           do jpoin = 1,size(PLEPP_CPLNG%interior_list_j)
              gesca(PLEPP_CPLNG%interior_list_j(jpoin)) = 1.0_rp
           end do
        end if
#else
        gesca(1:npoin) = 0.0_rp
#endif
     end if

  case( 'SVEGM' )
     !
     ! SVEGM: All State Dependent Variables (SDVs)
     !
     if ( INOTEMPTY ) ger3p => svegm_sld(:)

  case( 'CELEN' )
     !
     ! CELEN: Characteristic element length
     !
     if ( INOTEMPTY ) then
        call memgen(0_ip,nelem,0_ip)
        gesca(1:nelem) = celen_sld(1:nelem)
     end if

  case( 'FIXBO' )
     !
     ! FIXBO: Fixity code on boundaries
     !
     if( INOTEMPTY ) then
        call memgen(0_ip,npoin,0_ip)
        if( kfl_icodb > 0 ) then
           do iboun = 1,nboun
              pblty = ltypb(iboun)
              pnodb = nnode(pblty)
              do inodb = 1,pnodb
                 ipoin = lnodb(inodb,iboun)
                 gesca(ipoin) = max(real(kfl_fixbo_sld(iboun),rp),real(gesca(ipoin),rp))
              end do
           end do
           call PAR_INTERFACE_NODE_EXCHANGE(gesca,'MAX','IN MY CODE')
        else
           do ipoin = 1,npoin
              gesca(ipoin) = 0.0_rp
           end do
        end if
     end if

  case( 'BOCOD' )
     !
     ! BOCOD: Boundary codes
     !
     if( INOTEMPTY ) then
        call memgen(0_ip, npoin, 0_ip)
        if( kfl_icodb > 0 ) then
           do iboun = 1,nboun
              pblty = ltypb(iboun)
              pnodb = nnode(pblty)
              do inodb = 1,pnodb
                 ipoin = lnodb(inodb,iboun)
                 gesca(ipoin) = max(real(kfl_codbo(iboun),rp), real(gesca(ipoin),rp))
              end do
           end do
           call PAR_INTERFACE_NODE_EXCHANGE(gesca,'MAX','IN MY CODE')
        else
           return
        end if
     end if

  case( 'NOCOD' )
     !
     ! NOCOD: Node codes
     !
     if( INOTEMPTY ) then
        call memgen(0_ip,npoin,0_ip)
        if( kfl_icodn > 0 ) then
           do ipoin = 1,npoin
              gesca(ipoin) = real(kfl_codno(1,ipoin),rp)
           end do
        else
           do ipoin = 1,npoin
              gesca(ipoin) = 0.0_rp
           end do
        end if
     end if

  case( 'ELNOR' )
     !
     ! ELNOR: Normal vector for HEX08 and QUA04 (Interface/contact elements)
     !
     if( INOTEMPTY ) then
        call memgen(0_ip,ndime,nelem)
        do ielem = 1,nelem
           pelty = ltype(ielem)
           pnode = nnode(abs(pelty))
           if( pelty == HEX08 .or. pelty == QUA04 ) then
              call sld_csys_midsurface_element_normal(ielem, ndime, pnode, norvec)
              gevec(1:ndime,ielem) = norvec(1:ndime)
           else
              gevec(1:ndime,ielem) = 0.0_rp
           end if
        end do
     end if

  case( 'FEXTE' )
     !
     ! FEXTE: External nodal forces
     !
     if( INOTEMPTY ) then
        call memgen(0_ip,ndime,npoin)
        do ipoin = 1,npoin
           idofn = (ipoin-1)*ndime
           gevec(1:ndime,ipoin) = fexte_sld(idofn+1:idofn+ndime)
        end do
     end if

  case ( 'WETNO' )
     ! 'WETNO'
     if( INOTEMPTY ) then
        call memgen(zero,npoin,zero)
        call sld_wetno(gesca(1:npoin))
     endif

  case ( 'SBVNA' )
     !
     ! SBVNA: Solver bvnat
     !
     if( INOTEMPTY ) then
        call memgen(0_ip,ndime,npoin)
        if( solve_sol(1)%kfl_bvnat == 1_ip ) then
           gevec(1:ndime,1:npoin) = solve_sol(1)%bvnat(1:ndime,1:npoin)
        end if
        call PAR_INTERFACE_NODE_EXCHANGE(gevec,'SUM')
     end if

  case ( 'PARTI' )
     !
     ! PARTI: MPI Partitions
     !
     if( INOTEMPTY ) then
        call memgen(0_ip,nelem,0_ip)
        gesca(1:nelem) = real(kfl_paral,rp)
     end if

  case ( 'ELEPS' )
     !
     ! elastic strains
     !
     if ( INOTEMPTY ) then
  !     call memgen(zero, npoin, zero)
  !     call projec_elements_to_nodes(svegm_sld,gesca)
        if (kfl_foten_sld == 1_ip) call sld_calculate_pp_tensors
        call memgen(0_ip, 6_ip, npoin)
        gevec(1:nvoig_sld,1:npoin) = epsee_sld(1:nvoig_sld,1:npoin)
     end if

  case( 'ROTM1' )
     !
     ! ROTM1: Rotation matrix (1st vector)
     !
     if( INOTEMPTY ) then
        call memgen(0_ip,ndime,npoin)
        do ipoin = 1,npoin
           ibopo = lpoty(ipoin)
           if( ibopo > 0 ) then
              gevec(1:ndime,ipoin) = jacrot_du_dq_sld(1:ndime,1,ipoin)
           else
              gevec(1:ndime,ipoin) = 0.0_rp
           end if
        end do
     end if

  case( 'ROTM2' )
     !
     ! ROTM2: Rotation vector (2nd vector)
     !
     if( INOTEMPTY ) then
        call memgen(0_ip,ndime,npoin)
        do ipoin = 1,npoin
           ibopo = lpoty(ipoin)
           if( ibopo > 0 ) then
              gevec(1:ndime,ipoin) = jacrot_du_dq_sld(1:ndime,2,ipoin)
           else
              gevec(1:ndime,ipoin) = 0.0_rp
           end if
        end do
     end if

  case( 'ROTM3' )
     !
     ! ROTM3: Rotation vector (3rd vector)
     !
     if( ndime == 2_ip ) return
     if( INOTEMPTY ) then
        call memgen(0_ip,ndime,npoin)
        do ipoin = 1, npoin
           ibopo = lpoty(ipoin)
           if( ibopo > 0 ) then
              gevec(1:ndime,ipoin) = jacrot_du_dq_sld(1:ndime,3,ipoin)
           else
              gevec(1:ndime,ipoin) = 0.0_rp
           end if
        end do
     end if

  case( 'FCONT' )
     !
     ! FCONT: Contact force
     !
     if( INOTEMPTY ) then
        call memgen(0_ip,ndime,npoin)
        if( kfl_conta_sld /= 0_ip .and. kfl_rigid_sld == 0_ip ) then
           do ipoin = 1,npoin
              gevec(1:ndime,ipoin) = fcont_sld(1:ndime,ipoin)
           end do
        else
           gevec(:,:) = 0.0_rp
        end if
     end if

  case( 'DAMAG' )
     !
     ! DAMAG: Damage variables
     ! sm152: d1, dG, dK, d6
     ! sm154: d1, d2, d6
     if( kfl_damag_sld == 0_ip ) return
     if( INOTEMPTY ) then
        call memgen(0_ip,4_ip,nelem)
        do ielem = 1,nelem
           pgaus = ngaus(ltype(ielem))
           if(      lawst_sld(lmate(ielem)) == 152_ip ) then
              do isdva = 1,4
                 gevec(isdva,ielem) = sum(svegm_sld(ielem)%a(isdva+3,1:pgaus,ITER_K_STATE))/real(pgaus,rp)
              end do
           else if( lawst_sld(lmate(ielem)) == 154_ip ) then
              gevec(1,ielem) = sum(svegm_sld(ielem)%a( 5,1:pgaus,ITER_K_STATE))/real(pgaus,rp)
              gevec(2,ielem) = sum(svegm_sld(ielem)%a( 6,1:pgaus,ITER_K_STATE))/real(pgaus,rp)
              gevec(3,ielem) = sum(svegm_sld(ielem)%a(10,1:pgaus,ITER_K_STATE))/real(pgaus,rp)
              gevec(4,ielem) = sum(svegm_sld(ielem)%a(11,1:pgaus,ITER_K_STATE))/real(pgaus,rp)
           else
              gevec(:,ielem) = 0.0_rp
           end if
        end do
     else
        call memgen(0_ip,4_ip,1_ip)
     end if

  case ( 'FFLUI' )
     !
     ! FFLUI: Fluid force (FSI)
     !
     if( INOTEMPTY ) then
        call memgen(0_ip,ndime,npoin)
        if( solve_sol(1)%kfl_bvnat == 1_ip ) then
           do ipoin = 1,npoin
              gevec(1:ndime,ipoin) = solve_sol(1) % bvnat(1:ndime,ipoin)
           end do
           call rhsmod(ndime,gevec)
        end if
     end if
     
  case ( 'MICNL' )

     ! MICNL
     if(INOTEMPTY) then
        call memgen(zero, nelem, zero)
        gesca(1:nelem) = 0.0_rp
        do e = 1, nelem
           non_linear = fe2_is_non_linear(e)
           if (non_linear) then
                   gesca(e) = 1.0_rp
           else
                   gesca(e) = 0.0_rp
           endif
        end do
     endif

  case ( 'MICCO' )

     ! MICCO
     if(INOTEMPTY) then
        call memgen(zero, nelem, zero)
        gesca(1:nelem) = 0.0_rp
        do e = 1, nelem
           cost = fe2_get_cost(e)
           gesca(e) = real(cost,rp)
        end do
     endif

  case ( 'MICCV' )

     ! MICCV
     if(INOTEMPTY) then
        call memgen(zero, nelem, zero)
        gesca(1:nelem) = 1.0_rp
        do e = 1, nelem
           converged = fe2_has_converged(e)
           if (converged) then
                   gesca(e) = 1.0_rp
           else
                   gesca(e) = 0.0_rp
           endif
        end do
     endif

  case ( 'EPRIN' )
     !
     ! Principal stretches
     !
     if (INOTEMPTY) then
        kfl_foten_sld= kfl_foten_sld + 1
        if (kfl_foten_sld == 1) then   ! fotens computes caust, green and lepsi, so do it only the first kfote
           call sld_calculate_pp_tensors
        end if
     end if

     gesca => eprin_sld

  case ( 'CALCI' )
     ! Calcium in solidz
     if(INOTEMPTY) then
        call memgen(0_ip,npoin,0_ip)
        if( kfl_exmsld_ecc )then
          do ipoin=1,npoin
              gesca(ipoin)= calcium_ecc(ipoin)
          enddo
        else
          gesca(:)= 0.0_rp
        endif
     endif
     
  case ( 'EBFIL' )
     ! 
     ! EBFIL: Green strain EPSIL (green_sld) in longitudinal fiber direction, green_sld(1:nvoig_sld,1:npoin)
     !
     if( INOTEMPTY ) then
        if (.not. kfl_biofibers) call runend("SLD_OUTVAR: STRAIN PROJECTION ON FIBERS IS IMPLEMENTED ONLY FOR BIOFIBERS ON NODES")
        call biofib_point_nodal_fibers( gevec, 'LONGITUDINAL', 'CURRENT' )
        ! fotens computes caust, green and lepsi, so do it only the first kfote
        kfl_foten_sld = kfl_foten_sld + 1_ip
        if( kfl_foten_sld == 1_ip ) call sld_calculate_pp_tensors
        call memgen(0_ip,npoin,0_ip)
        do ipoin = 1, npoin
            call SM_voigt_to_tensor_second(ndime, matrix, green_sld(:, ipoin))
            norvec = gevec(1:ndime, ipoin)
            call maths_normalize_vector(ndime, norvec)
            gesca(ipoin) = maths_MULT_VxMxV(norvec, &
                                            ndime, matrix, &
                                            norvec, ndime)
        end do
     endif

  case ( 'SBFIL' )
     !
     ! SBFIL: Stress sigma (caust_sld) in longitudinal fiber direction, caust_sld(1:nvoig_sld,1:npoin)
     !
     if( INOTEMPTY ) then
        if (.not. kfl_biofibers) call runend("SLD_OUTVAR: STRESS PROJECTION ON FIBERS IS IMPLEMENTED ONLY FOR BIOFIBERS ON NODES")
        call biofib_point_nodal_fibers( gevec, 'LONGITUDINAL', 'CURRENT' )
        ! fotens computes caust, green and lepsi, so do it only the first kfote
        kfl_foten_sld = kfl_foten_sld + 1_ip
        if( kfl_foten_sld == 1_ip ) call sld_calculate_pp_tensors
        call memgen(0_ip,npoin,0_ip)
        do ipoin = 1, npoin
            call SM_voigt_to_tensor_second(ndime, matrix, caust_sld(:, ipoin))
            norvec = gevec(1:ndime, ipoin)
            call maths_normalize_vector(ndime, norvec)
            gesca(ipoin) = maths_MULT_VxMxV(norvec, &
                                            ndime, matrix, &
                                            norvec, ndime)
        end do
     endif

  case( 'STACK' )
     !
     ! STACK: Element stacking direction
     !
     if( INOTEMPTY ) then
        call memgen(0_ip, npoin, 0_ip)
        do ielem = 1,nelem
           pelty = ltype(ielem)
           pnode = nnode(pelty)
           if( pelty == HEX08 .or. &
               pelty == QUA04 .or. &
               pelty == PEN06  ) then
              hnode = int(pnode/2,ip)
              do inode = 1,hnode
                 jpoin = lnods(inode+hnode,ielem)
                 gesca(jpoin) = 1.0_rp
              end do
           end if
        end do
        call PAR_INTERFACE_NODE_EXCHANGE(gesca,'MAX','IN MY CODE')
      end if

   case ('EPSUL')
      !
      ! EPSUL: Green strain EPSIL in user supplied long direction
      !
      if (INOTEMPTY) then
         if (.not. strain_basis%on) call runend("SLD_OUTVAR: STRAIN VECTOR FIELD NOT PROVDED FOR EPSUL")
         call postprocess_strain_user_field(strain_basis%lng)
      end if

   case ('EPSUR')
      !
      ! EPSUR: Green strain EPSIL in user supplied radial direction
      !
      if (INOTEMPTY) then
         if (.not. strain_basis%on) call runend("SLD_OUTVAR: STRAIN VECTOR FIELD NOT PROVDED FOR EPSUR")
         call postprocess_strain_user_field(strain_basis%rad)
      end if

   case ('EPSUC')
      !
      ! EPSUC: Green strain EPSIL in user supplied circumferential direction
      !
      if (INOTEMPTY) then
         if (.not. strain_basis%on) call runend("SLD_OUTVAR: STRAIN VECTOR FIELD NOT PROVDED FOR EPSUC")
         call postprocess_strain_user_field(strain_basis%cir)
     end if
   
  case ( 'BVNAT' )
     !
     ! BVNAT: Pressure load on element boundaries
     !
     if( INOTEMPTY ) then
        call memgen(0_ip, npoin, 0_ip)
        do iboun = 1,nboun
           do inodb = 1,nnode(abs(ltypb(iboun)))
              ipoin = lnodb(inodb,iboun)
              do idime = 1,ndime
                 if( kfl_fixno_sld(idime,ipoin) == 0_ip ) then
                    gesca(ipoin) = bvnat_sld(1,iboun,ITER_K)
                 end if
              end do
           end do
        end do
     end if
     
  case ( 'TEMPE' )
     !
     ! TEMPE: Temperature
     !
     if( INOTEMPTY ) then
        call memgen(0_ip, npoin, 0_ip)
        do ipoin = 1,npoin
           gesca(ipoin) = therm(ipoin,ITER_K)
        end do
     end if
     
  case ( 'PRESS' )
     !
     ! PRESS: Pressure load on element faces
     !
     if( INOTEMPTY ) then
        call memgen(0_ip, npoin, 0_ip)
        do iboun = 1,nboun
           do inodb = 1,nnode(abs(ltypb(iboun)))
              ipoin = lnodb(inodb,iboun)
              do idime = 1,ndime
                 if( kfl_fixno_sld(idime,ipoin) == 0_ip ) then
                    gesca(ipoin) = bvnat_sld(1,iboun,ITER_K)
                 end if
              end do
           end do
        end do
     end if

  case ( 'NPOIN' )
     !
     ! NPOIN: Point number
     !
     if( INOTEMPTY ) then
        call memgen(0_ip, npoin, 0_ip)
        gesca(1:npoin) = real(lninv_loc(1:npoin),rp)
     end if

  case( 'PELTY' )
     !
     ! PELTY: Element type
     !
     if( INOTEMPTY ) then
        call memgen(0_ip, nelem, 0_ip)
        do ielem = 1,nelem
           gesca(ielem) = real(ltype(ielem),rp)
        end do
     end if

  case( 'ENDEN' )
     !
     ! ENDEN: Energy density calculated by the material models
     !
     if( INOTEMPTY ) then
        call memgen(0_ip, npoin, 0_ip)
        call smooth(enden_sld, gesca)
     end if

  case( 'WATER' )
     !
     ! WATER: Water content of sm400
     !
     if( INOTEMPTY ) then
        call memgen(0_ip, npoin, 0_ip)
        call smooth(water_sld, gesca)
     end if

  case( 'DONNA' )
     !
     ! DONNA: Donnan osmosis from sm400
     !
     if(kfl_donna_sld == 1_ip) then
        if( INOTEMPTY ) then
           call memgen(0_ip, npoin, 0_ip)
           call smooth(donna_gp, gesca)
        end if
     end if
     
  case default

     return

  end select
  !
  ! Postprocess
  !
  call outvar(&
       ivari,&
       ittim,cutim,postp(1) % wopos(:,ivari),MESH_ID=imesh)

contains
   subroutine postprocess_strain_user_field(ref_basis)
      use def_master, only: gdepo
      implicit NONE
      real(rp), dimension(:, :), intent(in) :: ref_basis
      real(rp)                             :: v(ndime)

      kfl_foten_sld = kfl_foten_sld + 1_ip
      if (kfl_foten_sld == 1_ip) call sld_calculate_pp_tensors

      call memgen(0_ip, npoin, 0_ip)
      do ipoin = 1, npoin
         call SM_voigt_to_tensor_second(ndime, matrix, green_sld(:, ipoin))
         v = maths_MULT_MxV(gdepo(1:ndime, 1:ndime, ipoin), ref_basis(1:ndime, ipoin), ndime)
         call maths_normalize_vector(ndime, v)
         gesca(ipoin) = maths_MULT_VxMxV(v, ndime, matrix, v, ndime)
      end do
   end subroutine postprocess_strain_user_field
  
end subroutine sld_outvar
