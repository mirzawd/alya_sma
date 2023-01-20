!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!
!> @defgroup Properties_Toolbox
!> Toolbox for creating, updating and computing properties
!> @{
!> @name    ToolBox for property creations
!> @file    mod_ker_proper.f90
!> @author  Guillaume Houzeaux
!> @brief   Update properties
!> @details Update properties
!
!------------------------------------------------------------------------

module mod_ker_updpro

  use def_kintyp,                   only : ip,rp,r1p
  use def_parame,                   only : pi,zero,one
  use def_master,                   only : cutim,momod
  use def_domain,                   only : nbset,lbsec
  use def_domain,                   only : xfiel
  use def_domain,                   only : kfl_field
  use def_master,                   only : NPOIN_TYPE
  use def_master,                   only : NELEM_TYPE
  use def_master,                   only : itinn,ittim
  use def_master,                   only : INOTMASTER
  use def_master,                   only : modul,mem_modul
  use def_master,                   only : ITASK_DOITER
  use def_master,                   only : ITASK_ENDSTE
  use def_master,                   only : ITASK_INIUNK
  use def_master,                   only : INOTMASTER
  use def_master,                   only : ID_TEMPER
  use def_master,                   only : ID_PARTIS
  use def_master,                   only : ID_CHEMIC
  use def_master,                   only : visck,sphek,tempe
  use def_master,                   only : visco,densi,condk
  use def_master,                   only : wmean,speci,nspec
  use def_master,                   only : gevec,tesgs,wmean_gp
  use def_master,                   only : densi_gp,conce_forw
  use def_master,                   only : condu_gp,sphec_gp,visco_gp
  use def_master,                   only : iavec
  use def_master,                   only : IMASTER,advec,tempe_gp
  use def_master,                   only : tempe_forw,untur_forw
  use def_master,                   only : fleve,prthe,kfl_coupl
  use def_master,                   only : veloc_forw,untur,conce
  use def_domain,                   only : nmate,nelem,npoin,ltypb
  use def_domain,                   only : ltype,lmate,ngaus,nboun
  use def_domain,                   only : lgaus,ndime,mnode,lnods
  use def_domain,                   only : ndime,nnode,kfl_naxis
  use def_domain,                   only : mgaus,elmar,coord,lelch
  use def_domain,                   only : lnnod,lnodb,nmatn,lelbo
  use def_domain,                   only : lmatn,hnatu,ntens,canhe
  use def_domain,                   only : canla,walld,heiov,rough
  use mod_parall,                   only : par_omp_nelem_chunk
  use mod_ker_polynomial,           only : ker_polynomial_proper
  use mod_ker_polynomial,           only : ker_polynomial_name 
  use mod_memory,                   only : memory_alloca
  use mod_memory,                   only : memory_deallo
  use mod_ker_regularization,       only : regul_k, regul_e, kfl_regularization
  use mod_interpolation,            only : COU_GET_INTERPOLATE_POINTS_VALUES
  use def_ker_proper,               only : UPDATE_PROPERTIES
  use def_ker_proper,               only : NUMBER_OF_PROPERTIES
  use def_ker_proper,               only : typ_valpr_ker
  use def_ker_proper,               only : lis_prope_ker
  use def_ker_proper,               only : mresp_ker,cpu_prope
  use mod_ker_proper_generic,       only : ker_proper_element_operations
  use mod_ker_proper,               only : ker_proper
  use def_kermod,                   only : kfl_reset
  use def_kermod,                   only : densi_ker,visco_ker
  use def_kermod,                   only : time_function
  use def_kermod,                   only : wallcoupling_extr_boun
  use def_kermod,                   only : cmu_st,condu_ker
  use def_kermod,                   only : delta_dom,rough_dom
  use def_kermod,                   only : gasco,kfl_adj_prob
  use def_kermod,                   only : kfl_canhe,kfl_canla
  use def_kermod,                   only : kfl_kxmod_ker,kfl_logva
  use def_kermod,                   only : kfl_prope,kfl_rough
  use def_kermod,                   only : lastm_ker,mixin_ker
  use def_kermod,                   only : number_time_function
  use def_kermod,                   only : is_interior,kfl_noslw_ker
  use def_kermod,                   only : walvi_ker,absor_ker
  use def_kermod,                   only : anipo_ker,poros_ker
  use def_kermod,                   only : dummy_ker,scatt_ker
  use def_kermod,                   only : sphea_ker,turmu_ker
  use def_kermod,                   only : thicl,kfl_nswel_ker
  use def_kermod,                   only : ptb_to_use

#ifndef PROPER_PRIVATE_OFF

  use mod_ker_proper_turmu,         only : ker_proper_turmu_walvi
  use mod_ker_proper_turmu,         only : ker_proper_turmu_smagorinsky_init
  use mod_ker_proper_turmu,         only : ker_proper_turmu_smagorinsky_operations
  use mod_ker_proper_turmu,         only : ker_proper_turmu_sigma_init
  use mod_ker_proper_turmu,         only : ker_proper_turmu_sigma_operations
  use mod_ker_proper_turmu,         only : ker_proper_turmu_wale_init
  use mod_ker_proper_turmu,         only : ker_proper_turmu_wale_operations
  use mod_ker_proper_turmu,         only : ker_proper_turmu_vreman_init
  use mod_ker_proper_turmu,         only : ker_proper_turmu_vreman_operations
  use mod_ker_proper_turmu,         only : ker_proper_turmu_ilsa_init
  use mod_ker_proper_turmu,         only : ker_proper_turmu_ilsa_operations
  use mod_ker_proper_turmu,         only : ker_proper_turmu_tkesgs_init
  use mod_ker_proper_turmu,         only : ker_proper_turmu_tkesgs_operations


#endif
                       
  implicit none
  private

  public :: ker_updpro

contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-11-23
  !> @brief   Property update
  !> @details Property update
  !> 
  !-----------------------------------------------------------------------

  subroutine ker_updpro(itask_xxx)

    use mod_ker_valve_proper,  only: ker_valve_proper_update_porosity
    use mod_communications, only : PAR_MAX
    implicit none

    integer(ip), optional, intent(in) :: itask_xxx
    integer(ip), save                 :: ipass(100) = 0
    integer(ip)                       :: itask, kfl_force
    integer(ip)                       :: imate,pnode,ielem,ipoin,kpoin
    integer(ip)                       :: pgaus,inode,pelty,inodb,pblty
    integer(ip)                       :: iprop,ilaws,pnodb,iboun,pgaub,igaub
    integer(ip)                       :: iresp,kfl_updpr(NUMBER_OF_PROPERTIES),idime,igaus
    integer(ip)                       :: ispec,jspec,jdime,dummi, icomp
    real(rp)                          :: gpsha(mnode,mgaus)
    real(rp)                          :: gpder(ndime,mnode,mgaus)
    real(rp)                          :: gpcar(ndime,mnode,mgaus)
    real(rp)                          :: gphes(ntens,mnode,mgaus)
    real(rp)                          :: gpvol(mgaus)
    real(rp)                          :: dummr(ndime,mgaus)
    real(rp)                          :: elcod(ndime,mnode)
    real(rp)                          :: elvel(ndime,mnode), eltur(2, mnode), elwal(mnode)
    real(rp)                          :: gpcor(ndime,mgaus)
    real(rp)                          :: gp_temperature(mgaus)
    real(rp)                          :: gp_grad_tem(ndime,mgaus), gp_grad_wmean(ndime,mgaus)
    real(rp)                          :: gpvel(ndime,mgaus)
    real(rp)                          :: gptan(ndime,mgaus),gbtan(ndime)
    real(rp)                          :: gpgve(ndime,ndime,mgaus)
    real(rp)                          :: gptur,beta,gpvis_nsw(mgaus)

    real(rp)                          :: ellev(mnode),eltem(mnode),elwmean(mnode),elY(mnode)
    real(rp)                          :: elden(mnode)
    real(rp), pointer                 :: auxde(:)
    real(rp)                          :: prope_wat,prope_air
    real(rp)                          :: hleng(3),tragl(ndime,ndime)
    real(rp)                          :: gpY,rho_l,rho_g
    real(rp)                          :: eps,ooeps,oovpi,f,phi,kap
    real(rp)                          :: densi_phase(2),z0,rho,ustar

    real(rp)                          :: T,T0,mu0,mu,S0,const, heigh, cd, cdlad, auxva , lad
    real(rp)                          :: aux_prop,time1,time2
    real(rp)                          :: lscale,tscale

    real(rp)                          :: phiij,sumphi(nspec),gptem,gbtem,gbwme,gphog(mgaus)
    real(rp)                          :: zmaxi,zmaxf,gpcah(mgaus),gphei(mgaus)
    real(rp)                          :: lamax, zz, n, shaaa, hoveg, gpwmean, gpden(mgaus)
    real(rp)                          :: gpvis(mgaus), gpnut,gpnu,X,cv1,fv1
    real(rp)                          :: F2, kinen, gpwal,gpwal2(mgaus), omega, gpvor, a1, as
    real(rp)                          :: epsil, relat, gbwal, gbden(1), gbvis(1),velmo,tanmo
    real(rp)                          :: dgpkin(mnode),dgpome(mnode), arg2,darg2_dkin(mnode)
    real(rp)                          :: darg2_dome(mnode), dgpnut_du(mnode,ndime)
    real(rp)                          :: gpvor_mat(ndime,ndime)
    real(rp)                          :: dF2_dkin(mnode), dF2_dome(mnode), dgpnut_dkin(mnode)
    real(rp)                          :: dgpnut_dome(mnode), dgpvor_du(mnode,ndime)
    real(rp)                          :: dX_dnut,dfv1_dnut,dgpmut_dnut(mnode)
    real (rp),save                    :: LAD_aver


    integer(ip)                       :: which_time, icoun, ncoun, kdime, iclas_Y
    ! ke models
    real(rp)                          :: A0, As_r , W, phi_r, cmu, divve, simgr(3,3),  Wsqr6 
    real(rp)                          :: uaste,seci4, Cr,f0, f02, alphast, ReT
    real(rp)                          :: regue, reguk, sigmr, dzlay, dcoef, zdamp
    real(rp), save                    :: ztop

    type(typ_valpr_ker), pointer      :: prope_ker
    type(lis_prope_ker), pointer      :: listp_ker(:)

    real(rp),            pointer      :: conce_save(:,:,:)  !!! value to save conce for adjoint
    real(rp),            pointer      :: tempe_save(:,:)    !!! value to save tempe for adjoint
    real(rp),            pointer      :: untur_save(:,:,:)  !!! value to save untur for adjoint
    real(rp),            pointer      :: advec_save(:,:,:)  !!! value to save veloc for adjoint

    !
    ! for VALVE porosity model
    !
    real(rp)                          :: coeff
    !
    ! for mixing length model 
    !
    real(rp)                          :: gbvel(ndime),gbgve(ndime,ndime),eltan(ndime,mnode)
    real(rp)    , pointer             :: tau_wall_boun(:,:)

    real(rp)                          :: gpmve,darcy,forch

    ! Level set / SPRAY

    real(rp)                          :: delta

    ! temporal functions in properties
    integer(ip)                       :: fun_code, ifunc
    real(rp), external                :: funcre


    if( IMASTER .or. kfl_prope == 0 ) return

    call cputim(time1)

    nullify(tau_wall_boun)
    nullify(prope_ker)
    nullify(listp_ker)

    if (present(itask_xxx)) then
       itask=itask_xxx
       kfl_force=0
    else
       itask=-1000
       kfl_force=1
    endif

    if (kfl_reset == 1) then
       which_time=3
    else
       which_time=1
    endif

    !----------------------------------------------------------------------
    !
    ! Assign forward values for adjoint case to temper and conce
    !
    !----------------------------------------------------------------------

    if( kfl_adj_prob == 1_ip) then
       conce_save => conce
       tempe_save => tempe
       untur_save => untur
       advec_save => advec
       conce      => conce_forw
       tempe      => tempe_forw
       untur      => untur_forw
       advec      => veloc_forw
    end if

    allocate(listp_ker(NUMBER_OF_PROPERTIES))

    listp_ker( 1) % prop => densi_ker
    listp_ker( 2) % prop => visco_ker
    listp_ker( 3) % prop => poros_ker
    listp_ker( 4) % prop => condu_ker
    listp_ker( 5) % prop => sphea_ker
    listp_ker( 6) % prop => dummy_ker
    listp_ker( 7) % prop => turmu_ker
    listp_ker( 8) % prop => absor_ker
    listp_ker( 9) % prop => scatt_ker
    listp_ker(10) % prop => mixin_ker
    listp_ker(11) % prop => anipo_ker
    listp_ker(12) % prop => walvi_ker

    !----------------------------------------------------------------------
    !
    ! Constant properties
    !
    !----------------------------------------------------------------------

    do iprop = 1,NUMBER_OF_PROPERTIES
       prope_ker =>  listp_ker(iprop) % prop
       if( prope_ker % kfl_exist == 1 ) then
          do imate = 1,nmate
             if( prope_ker % wlaws(imate) == 'CONST' ) &
                  prope_ker % value_const(1:prope_ker % dim,imate) = prope_ker % rlaws(1:prope_ker % dim,imate)
          end do
       end if
    end do

    !----------------------------------------------------------------------
    !
    ! Check if a property depends on the last module solved: LASTM_KER
    !
    !----------------------------------------------------------------------

    do iprop = 1,NUMBER_OF_PROPERTIES
       prope_ker => listp_ker(iprop)%prop
       kfl_updpr(iprop) = 0
       if( prope_ker % kfl_exist == 1 ) then
          do imate = 1,nmate
             if( prope_ker % wlaws(imate) /= 'CONST' ) then
                ilaws = prope_ker % ilaws(imate)
                do iresp = 1,mresp_ker
                   if( prope_ker % llaws(ilaws) % lresp(iresp) == lastm_ker&
                        .and.(prope_ker % update(1, imate)==itask&
                        .or.prope_ker % update(2, imate)==itask)) then
                      kfl_updpr(iprop) = 1
                   else if ( prope_ker % llaws(ilaws) % lresp(iresp) == -1 ) then
                      kfl_updpr(iprop) = 1
                   end if
                end do
             end if
          end do
          if ( ( kfl_force == 1 ) .or. ( itask == ITASK_INIUNK ) )  then
             kfl_updpr(iprop) = 1
          endif
       end if
    end do

    !----------------------------------------------------------------------
    !
    ! Variable properties
    !
    !----------------------------------------------------------------------

    do iprop = 1,NUMBER_OF_PROPERTIES

       if( kfl_updpr(iprop) == 1 ) then

          prope_ker => listp_ker(iprop) % prop

          if( prope_ker % kfl_exist == 1 ) then

             do imate = 1,nmate
                if( itask == ITASK_INIUNK .and. prope_ker % value_default(imate) > 0.0_rp ) then
                   !----------------------------------------------------------------------
                   !
                   ! Run is starting (INIUNK): impose a default value if it exists (> 0)
                   !
                   !----------------------------------------------------------------------
                   ilaws = prope_ker % ilaws(imate)

                   if( prope_ker % llaws(ilaws) % where == 'IPOIN' ) then

                      do ipoin = 1,npoin
                         prope_ker % value_ipoin(ipoin) = prope_ker % value_default(imate)
                      end do
                      if( prope_ker % llaws(ilaws) % kfl_gradi == 1 ) then
                         do ipoin = 1,npoin
                            prope_ker % grval_ipoin(:,ipoin) = 0.0_rp
                         end do
                      end if

                   else if( prope_ker % llaws(ilaws) % where == 'IELEM' ) then

                      do ielem = 1,nelem
                         prope_ker % value_ielem(ielem) % a = prope_ker % value_default(imate)
                      end do
                      do iboun = 1,nboun
                         prope_ker % value_iboun(iboun) % a = prope_ker % value_default(imate)
                      end do

                      if( prope_ker % llaws(ilaws) % kfl_gradi == 1 ) then
                         do ielem = 1,nelem
                            prope_ker % grval_ielem(ielem) % a = 0.0_rp
                         end do
                      end if

                   else if( prope_ker % wlaws(imate) == 'CONST' ) then

                      continue ! Nothing to do, it was already imposed
                      !prope_ker % value_const(imate) = prope_ker % value_default(imate)

                   end if

                else if( prope_ker % wlaws(imate) == 'IMMER' ) then

                   !----------------------------------------------------------------------
                   !
                   ! Immersed boundary mixing
                   !
                   !----------------------------------------------------------------------

                   if( iprop == 10 ) then

                      prope_wat = prope_ker % rlaws(1,imate)
                      prope_air = prope_ker % rlaws(2,imate)

                      do ielem = 1,nelem
                         pelty = ltype(ielem)
                         if( lmate(ielem) == imate .and. pelty > 0 ) then
                            pgaus = ngaus(pelty)
                            pnode = lnnod(ielem)
                            do inode = 1,pnode
                               ipoin = lnods(inode,ielem)
                               ellev(inode) = fleve(ipoin,which_time)
                               elcod(1:ndime,inode) = coord(1:ndime,ipoin)
                            end do
                            call elmca2(&
                                 pnode,pgaus,0_ip,elmar(pelty)%weigp,elmar(pelty)%shape,&
                                 elmar(pelty)%deriv,elmar(pelty)%heslo,elcod,gpvol,gpsha,&
                                 gpder,gpcar,gphes,ielem)
                            call ker_biflui(&
                                 3_ip,pgaus,pnode,gpsha,gpcar,ellev,prope_wat,prope_air, &
                                 thicl,hleng,prope_ker % value_ielem(ielem) % a,         &
                                 dummr)            
                         end if
                      end do

                      prope_ker % kfl_nedsm = 1

                   else

                      call runend('THIS PROPERTY DOES NOT MAKE SENSE FOR IMMERSED BOUNDARY')

                   end if

                else if( prope_ker % wlaws(imate) == 'BIFLU' ) then

                   !----------------------------------------------------------------------
                   !
                   ! Bifluid model
                   !
                   !----------------------------------------------------------------------

                   prope_wat = prope_ker % rlaws(1,imate)
                   prope_air = prope_ker % rlaws(2,imate)

                   do ielem = 1,nelem
                      pelty = ltype(ielem)
                      if( lmate(ielem) == imate .and. pelty > 0 ) then
                         pgaus = ngaus(pelty)
                         pnode = lnnod(ielem)
                         do inode = 1,pnode
                            ipoin = lnods(inode,ielem)
                            ellev(inode) = fleve(ipoin,which_time)
                            elcod(1:ndime,inode) = coord(1:ndime,ipoin)
                         end do
                         call elmca2(&
                              pnode,pgaus,0_ip,elmar(pelty)%weigp,elmar(pelty)%shape,&
                              elmar(pelty)%deriv,elmar(pelty)%heslo,elcod,gpvol,gpsha,&
                              gpder,gpcar,gphes,ielem)
                         call ker_biflui(&
                              1_ip,pgaus,pnode,gpsha,gpcar,ellev,prope_wat,prope_air, &
                              thicl,hleng,prope_ker % value_ielem(ielem) % a,         &
                              dummr)
                      end if
                   end do

                   do iboun = 1,nboun
                      pblty = abs(ltypb(iboun))
                      pnodb = nnode(pblty)
                      ielem = lelbo(iboun)
                      if( lmate(ielem) == imate ) then
                         pelty = ltype(ielem)
                         if( pelty > 0 ) then
                            pgaub = ngaus(pblty)
                            do inodb = 1,pnodb
                               ipoin = lnodb(inodb,iboun)
                               ellev(inodb) = fleve(ipoin,which_time)
                            end do
                            call ker_biflui(&
                                 2_ip,pgaub,pnodb,elmar(pblty) % shape,gpcar,ellev,prope_wat,prope_air,&
                                 thicl,hleng,prope_ker % value_iboun(iboun) % a,dummr)
                         end if
                      end if
                   end do
                   prope_ker % kfl_nedsm = 1

                else if( prope_ker % wlaws(imate) == 'BIFL2' ) then

                   !----------------------------------------------------------------------
                   !
                   ! Bifluid model - including gradients
                   !
                   !----------------------------------------------------------------------

                   prope_wat = prope_ker % rlaws(1,imate)
                   prope_air = prope_ker % rlaws(2,imate)
                   do ielem = 1,nelem
                      if( lmate(ielem) == imate ) then
                         pelty = ltype(ielem)
                         if( pelty > 0 ) then
                            pgaus = ngaus(pelty)
                            pnode = lnnod(ielem)
                            do inode = 1,pnode
                               ipoin = lnods(inode,ielem)
                               ellev(inode) = fleve(ipoin,which_time)
                               elcod(1:ndime,inode) = coord(1:ndime,ipoin)
                            end do
                            call elmca2(&
                                 pnode,pgaus,0_ip,elmar(pelty)%weigp,elmar(pelty)%shape,&
                                 elmar(pelty)%deriv,elmar(pelty)%heslo,elcod,gpvol,gpsha,&
                                 gpder,gpcar,gphes,ielem)
                            call ker_biflu2(&
                                 1_ip,pgaus,pnode,gpsha,gpcar,ellev,prope_wat,prope_air, &
                                 thicl,hleng,prope_ker % value_ielem(ielem) % a,         &
                                 prope_ker % grval_ielem(ielem) % a)
                         end if
                      end if
                   end do

                   do iboun = 1,nboun
                      pblty = abs(ltypb(iboun))
                      pnodb = nnode(pblty)
                      ielem = lelbo(iboun)
                      if( lmate(ielem) == imate ) then
                         pelty = ltype(ielem)
                         if( pelty > 0 ) then
                            pgaub = ngaus(pblty)
                            do inodb = 1,pnodb
                               ipoin = lnodb(inodb,iboun)
                               ellev(inodb) = fleve(ipoin,which_time)
                            end do
                            call ker_biflu2(&
                                 2_ip,pgaub,pnodb,elmar(pblty) % shape,gpcar,ellev,prope_wat,prope_air,&
                                 thicl,hleng,prope_ker % value_iboun(iboun) % a,dummr)
                         end if
                      end if
                   end do

                   prope_ker % kfl_nedsm = 1

                else if( prope_ker % wlaws(imate) == 'NAHME' ) then

                   !----------------------------------------------------------------------
                   !                               
                   ! Nahme's law: mu = mu_0 exp(-beta*(T-T0))
                   !
                   !----------------------------------------------------------------------

                   mu0  = prope_ker % rlaws(1,imate)
                   beta = prope_ker % rlaws(2,imate)
                   T0   = prope_ker % rlaws(3,imate)
                   do ielem = 1,nelem
                      if( lmate(ielem) == imate ) then
                         pelty = ltype(ielem)
                         if( pelty > 0 ) then
                            pgaus = ngaus(pelty)
                            pnode = lnnod(ielem)
                            do inode = 1,pnode
                               ipoin = lnods(inode,ielem)
                               eltem(inode) = tempe(ipoin,1)
                            end do
                            do igaus = 1,pgaus
                               gp_temperature(1) = 0.0_rp
                               do inode = 1,pnode
                                  gp_temperature(1) = gp_temperature(1) &
                                       + elmar(pelty) % shape(inode,igaus) * eltem(inode)
                               end do
                               prope_ker % value_ielem(ielem) % a(igaus) = mu0 * exp(-beta*(gp_temperature(1)-T0))
                            end do
                         end if
                      end if
                   end do
                   prope_ker % kfl_nedsm = 1

                else if( prope_ker % wlaws(imate) == 'SUTHE' ) then

                   !----------------------------------------------------------------------
                   !
                   !                                        T0 + S
                   ! Sutherland's law: mu/mu_0 = (T/T0)^1.5 --------
                   !                                        T  + S
                   !
                   !----------------------------------------------------------------------

                   mu0   = prope_ker % rlaws(1,imate)
                   T0    = prope_ker % rlaws(2,imate)
                   S0    = prope_ker % rlaws(3,imate)
                   if (S0.lt. epsilon(1.0_rp)) S0  = 110.0_rp
                   if (T0.lt. epsilon(1.0_rp)) T0  = 298.0_rp
                   do kpoin = 1,nmatn(imate)
                      ipoin = lmatn(imate) % l(kpoin)
                      T     = tempe(ipoin,which_time)
                      prope_ker % value_ipoin(ipoin) = mu0 * (T/T0)**1.5_rp * (T0+S0) / (T+S0)
                   end do
                else if( prope_ker % wlaws(imate) == 'CANOP' ) then
                   !----------------------------------------------------------------------
                   !
                   ! Canopy model
                   !
                   !----------------------------------------------------------------------
                   cd    = prope_ker % rlaws(2,imate) ! cd (generally =0.2)
                   lad   = prope_ker % rlaws(3,imate) ! leaf area density
                   zmaxf = prope_ker % rlaws(4,imate) ! zmaxf = zmax/heigh
                   ! if zmaxf is positive , then an analytical LAD will be prescribed
                   !from LALIC and MIHAILOVIC, an Empirical Relation Describing Leaf-Area Density inside the Forest
                   !http://journals.ametsoc.org/doi/pdf/10.1175/1520-0450%282004%29043%3C0641%3AAERDLD%3E2.0.CO%3B2
                   ! only one time 
                   if ( ipass(imate)==0 ) then
                      ! if canopy height field does not exist
                      ! a constant canopy field is used

                      if ( kfl_canhe == 0) then  ! constant canopy height from ker.dat property
                         do ielem = 1,nelem
                            if( lmate(ielem) == imate ) then
                               pelty = ltype(ielem)
                               if( pelty > 0 ) then
                                  pnode = lnnod(ielem)
                                  do inode = 1,pnode
                                     ipoin = lnods(inode,ielem)
                                     canhe(ipoin) = prope_ker % rlaws(1,imate)
                                  end do
                               end if
                            end if
                         end do
                      end if
                      !
                      !  Prescription of vertical varible density function
                      !  through zmaxf, valid for constant LAD
                      if (kfl_canla==0_ip.and.zmaxf.gt.0.01_rp) then 
                         !
                         ! Obtain the integral for a height = 1, later it will be  multiplied by the real height
                         !
                         ncoun =1000 ! integration discretization
                         LAD_aver = 0.0_rp   ! averaged LAD
                         heigh = 1.0_rp
                         zmaxi = zmaxf*heigh

                         do icoun = 1, ncoun -1
                            zz = real(icoun,rp)/real(ncoun,rp)*heigh
                            if (zz.lt.zmaxi) then
                               n = 6.0_rp
                            else
                               n = 0.5_rp
                            end if
                            LAD_aver = LAD_aver + (((heigh -zmaxi) /(heigh - zz))**n )&
                                 *exp(n*(1.0_rp-(heigh -zmaxi) /(heigh - zz)))
                         end do

                         LAD_aver = LAD_aver  + 0.5_rp*(((heigh -zmaxi) /(heigh ))**6.0_rp )&
                              *exp(6.0_rp*(1.0_rp-(heigh -zmaxi) /(heigh )))

                         LAD_aver = LAD_aver  + 0.5_rp*(((heigh -zmaxi) /(1.0e-4_rp))**0.5_rp )&
                              *exp(0.5_rp*(1.0_rp-(heigh -zmaxi) /(1.0e-4_rp)))


                         LAD_aver = LAD_aver * heigh / real(ncoun,rp) ! is the averaged LAD
                      end if
                   end if
                   ipass(imate) = 1


                   !
                   ! Assigns canopy porosity to each gauss point cd*LAD*|veloc|
                   ! Assigns canopy height and height over terrain to each gauss point
                   ! This term is used in nastin and turbul
                   !
                   do ielem = 1,nelem
                      if( lmate(ielem) == imate ) then
                         pelty = ltype(ielem)
                         if( pelty > 0 ) then
                            pgaus = ngaus(pelty)
                            pnode = lnnod(ielem)
                            gpvel = 0.0_rp
                            gphog = 0.0_rp ! height over terrain
                            gpcah = 0.0_rp ! canopy height
                            do inode = 1,pnode
                               ipoin = lnods(inode,ielem)
                               do idime = 1,ndime
                                  elvel(idime, inode) = advec(idime, ipoin,1)
                                  do igaus =1, pgaus
                                     gpvel(idime, igaus) = gpvel(idime,igaus) + &
                                          elvel(idime, inode) * elmar(pelty) % shape(inode,igaus)
                                  end do
                               end do
                               do igaus =1, pgaus
                                  shaaa = elmar(pelty) % shape(inode,igaus)
                                  gpcah(igaus) = gpcah(igaus) + canhe(ipoin)* shaaa
                                  gphog(igaus) = gphog(igaus) + heiov(ipoin)* shaaa
                               end do
                            end do

                            !Assigns canopy porosity to each gauss point cd*LAD*|veloc|

                            do igaus =1, pgaus
                               heigh = gpcah(igaus) ! forest height
                               zmaxi = zmaxf*heigh
                               hoveg = gphog(igaus)
                               !
                               ! only assigns for igaus below forest height,
                               ! when forest height is higher than a lower limit (Iberdrola)
                               !
                               if ((hoveg.lt.heigh)) then !.and.(heigh.gt.h_can_critical)) then 
                                  auxva =0.0_rp
                                  do idime =1, ndime
                                     auxva = auxva + gpvel(idime, igaus)*gpvel(idime, igaus)
                                  end do
                                  !
                                  ! leaf area index
                                  !
                                  if (kfl_canla>0) then !Load LAD from field
                                     lad =0.0_rp
                                     do inode = 1,pnode
                                        lad = lad + elmar(pelty) % shape(inode,igaus)*canla(lnods(inode,ielem))
                                     end do
                                     cdlad = cd*lad
                                  else if (zmaxf.le.0.01_rp) then ! uniform LAD canopy
                                     cdlad = cd*lad
                                  else     ! LAD profile prescribed from vertical function
                                     !non uniform canopy distribution from LALIC and MIHAILOVIC,
                                     !An Empirical Relation Describing Leaf-Area Density inside the Forest
                                     !http://journals.ametsoc.org/doi/pdf/10.1175/1520-0450%282004%29043%3C0641%3AAERDLD%3E2.0.CO%3B2
                                     lamax = lad / LAD_aver
                                     if (hoveg.gt.0.999_rp*heigh ) then
                                        cdlad = 0.0_rp
                                     else if (hoveg.lt.zmaxi) then
                                        cdlad =   cd*lamax * (((heigh -zmaxi) /(heigh -hoveg))**6.0_rp )*&
                                             exp(6.0_rp*(1.0_rp-(heigh -zmaxi) /(heigh -hoveg)))
                                     else
                                        cdlad =   cd*lamax * (((heigh -zmaxi) /(heigh -hoveg))**0.5_rp )*&
                                             exp(0.5_rp*(1.0_rp-(heigh -zmaxi) /(heigh -hoveg)))
                                     end if

                                  end if
                                  prope_ker % value_ielem(ielem) % a(igaus) = densi_ker % value_const(1,imate)&
                                       *cdlad* sqrt(auxva)
                               else  ! above forest height
                                  prope_ker % value_ielem(ielem) % a(igaus) = 0.0_rp
                               end if
                            end do
                         end if
                      end if
                   end do
                else if( prope_ker % wlaws(imate) == 'DAMPI' ) then
                   !----------------------------------------------------------------------
                   !
                   !  Damping porosity term +rho*damping*u
                   !  This term is used to damp inertial  oscilations above the 
                   ! atmosprehic boundary layer
                   !  The damping coeff is imposed smoothly on the top (+Z) of the domain
                   !
                   ! 
                   !----------------------------------------------------------------------
                   dcoef  = prope_ker % rlaws(1_ip,imate) !  coeff
                   dzlay  = prope_ker % rlaws(2_ip,imate) ! thickness of damping layer 
                   ztop   = prope_ker % rlaws(3_ip,imate) ! ztop  (not necessary)
                   if (ipass(imate)==0) then ! damping and forest are mutually excluded
                      ztop = -1.0e25_rp
                      do ipoin=1, npoin
                         ztop = max(coord(3_ip,ipoin),ztop)
                      end do
                      call PAR_MAX(ztop,'IN MY CODE WITHOUT MASTER',INCLUDE_ROOT=.true.)
                   end if
                   ipass(imate) = 1_ip
                   zdamp  = ztop -dzlay                 ! height where damping starts
 
                   !
                   ! Assigns damping to over gauss points 
                   ! This term is used in nastin 
                   !
                   do ielem = 1,nelem
                      if( lmate(ielem) == imate ) then
                         pelty = ltype(ielem)
                         if( pelty > 0_ip ) then
                            !obtains gp height
                            pgaus = ngaus(pelty)
                            pnode = lnnod(ielem)
                            gphei = 0.0_rp
                            do inode = 1_ip,pnode
                               ipoin = lnods(inode,ielem)
                               gphei( 1:pgaus) = gphei(1:pgaus) + &
                                       coord(3,ipoin)* elmar(pelty) % shape(inode,1:pgaus)
                            end do
                            !Assigns damping to each gauss point 
                         
                            if (prope_ker % name(1:5) == 'ANIPO') then !ANISOTROPIC DAMPING 
                               ! it returns a matrix orderer as follows, that enter into NS porous term as
                               !%a (idime, jdime, igaus) *veloc(jdime)* e(idime)
                               !matrix component icomp is: 
                               !icomp=ndime*ndime*(igaus-1) + ndime*(jdime-1) + idime  
                               do igaus =1_ip, pgaus 
                                  
                                  icoun = 9_ip*(igaus-1_ip) 
                                  ! initialization of matrix components idime=1, ndime
                                  prope_ker % value_ielem(ielem) % a(icoun + 1:icoun + 9) = 0.0_rp
                                  if ((gphei(igaus).gt.zdamp)) then ! buffer zone, aply damping
                                     ! damping value
                                     coeff =  densi_ker % value_const(1,imate)&
                                          *0.5_rp*dcoef*(1.0_rp+ cos(pi*(ztop-gphei(igaus))/dzlay) )
                                     !horizontal components 1,1 and 2,2 a= coeff
                                     !vertical   component  3,3  a= 500*coeff (5e-3, 1e-5)                                     
                                     icomp = icoun + 1_ip
                                     prope_ker % value_ielem(ielem) % a(icomp) = coeff
                                     icomp = icoun + ndime +2_ip
                                     prope_ker % value_ielem(ielem) % a(icomp) = coeff
                                     icomp = icoun + 2_ip*ndime +3_ip
                                     prope_ker % value_ielem(ielem) % a(icomp) = coeff*1000.0_rp
                                     !(coeff, 0, 0, 0,coeff,0,0,0,coeff*500)                                  
                                  end if
                               end do
                               
                            else  ! ISOTROPIC DAMPING
                               do igaus =1, pgaus
                               !
                                  if (dzlay.gt.3000.0) then ! damping for inertial oscillations all domain
                                     if (gphei(igaus).gt.1000.0) then 
                                        prope_ker % value_ielem(ielem) % a(igaus) = &
                                             densi_ker % value_const(1,imate)&
                                             *dcoef
                                     else !linear up to 1000m 
                                        prope_ker % value_ielem(ielem) % a(igaus) = &
                                             densi_ker % value_const(1,imate)&
                                             *dcoef*gphei(igaus)/1000.0
                                     end if
                                  else if ((gphei(igaus).gt.zdamp)) then !buffer zone, aply damping
                                     ! damping value
                                     prope_ker % value_ielem(ielem) % a(igaus) = &
                                          densi_ker % value_const(1,imate)&
                                          *0.5_rp*dcoef*(1.0_rp+ cos(pi*(ztop-gphei(igaus))/dzlay) )
                                  else                                      
                                     prope_ker % value_ielem(ielem) % a(igaus) = 0.0_rp
                                  end if
                               end do
                            end if
                            
                         end if
                      end if
                   end do
                   
                else if( prope_ker % wlaws(imate) == 'VALVE' ) then
                   !----------------------------------------------------------------------
                   !
                   ! Valve model
                   !
                   !----------------------------------------------------------------------
                    call ker_valve_proper_update_porosity(prope_ker, imate)

                else if( prope_ker % wlaws(imate) == 'DAFOR' ) then
                   !----------------------------------------------------------------------
                   !
                   ! Darcy-Forchheimer porosity
                   ! S_{i}=-(\mu D+ 1/2 \rho|u_{jj}|F) u_{i}
                   ! beware it only works with contant density and viscosity (visco_ker % value_const(1,imate))
                   !
                   !----------------------------------------------------------------------
                   darcy = prope_ker % rlaws(1_ip,imate)  ! linear term
                   forch = prope_ker % rlaws(2_ip,imate)  ! quadratic term

                   do ielem = 1,nelem
                      if( lmate(ielem) == imate ) then
                         pelty = ltype(ielem)
                         if( pelty > 0 ) then
                            pgaus = ngaus(pelty)
                            pnode = lnnod(ielem)
                            gpvel = 0.0_rp
                            do inode = 1,pnode
                               ipoin = lnods(inode,ielem)

                               do idime = 1,ndime
                                  elvel(idime, inode) = advec(idime, ipoin,1)
                                  do igaus =1, pgaus
                                     gpvel(idime, igaus) = gpvel(idime,igaus) + &
                                          elvel(idime, inode) * elmar(pelty) % shape(inode,igaus)
                                  end do
                               end do
                            end do

                            do igaus =1, pgaus
                               gpmve = sqrt( dot_product(gpvel(1:ndime, igaus), gpvel(1:ndime, igaus)))
                               prope_ker % value_ielem(ielem) % a(igaus) = &
                                    darcy * visco_ker % value_const(1,imate)  &
                                    + 0.5_rp * forch * densi_ker % value_const(1,imate) * gpmve 
                            end do
                         end if
                      end if
                   end do
                else if( prope_ker % wlaws(imate) == 'FUNCT' ) then
                   !----------------------------------------------------------------------
                   !
                   ! Time function variable porosity
                   ! 
                   !
                   !----------------------------------------------------------------------

                   ! Define the value of the property to impose

                   aux_prop=0.0_rp

                   if(prope_ker % time_function(imate) .ne. 'NULL ') then
                      fun_code=0_ip
                      do ifunc = 1,number_time_function
                         if( trim(prope_ker % time_function(imate)) == trim(time_function(ifunc) % name) ) then
                            fun_code = ifunc
                            exit
                         end if
                      end do
                   else
                      call runend('KER_UPDPRO: A TIME FUNCTION MUST BE DEFINED')
                   endif

                   aux_prop = funcre( time_function(fun_code) % parameters, &
                        time_function(fun_code) % npara,      &
                        time_function(fun_code) % kfl_type,   &
                        cutim)


                   do ielem = 1,nelem
                      if( lmate(ielem) == imate ) then
                         prope_ker % value_ielem(ielem) % a = aux_prop
                      endif
                   end do
                   do iboun = 1,nboun
                      ielem = lelbo(iboun)
                      if( lmate(ielem) == imate ) then
                         prope_ker % value_iboun(iboun) % a = aux_prop
                      endif
                   end do


                else if( prope_ker % wlaws(imate) == 'GAUSS' ) then
                   !----------------------------------------------------------------------
                   !
                   ! Gaussian filter to generate turbulent fluctuations by diffusion
                   ! A. Kempf, M. Klein, J. Janicka, Flow Turbul. Combust. (2005)
                   !
                   !       D = C_D * L**2 / T
                   !----------------------------------------------------------------------
                   lscale = prope_ker % rlaws(1,imate)                   ! Target length scale
                   tscale = 1.0_rp / prope_ker % rlaws(2,imate)          ! Temporal integration
                   const  = 0.15915_rp                                   ! C_D = 1 / ( 2*pi)

                   do ielem = 1,nelem
                      if( lmate(ielem) == imate ) then
                         pelty = ltype(ielem)
                         if( pelty > 0 ) then
                            pgaus = ngaus(pelty)
                            pnode = lnnod(ielem)

                            do igaus = 1,pgaus
                               gptur = const * lscale * lscale * tscale              ! D = C_D * L**2 / T
                               prope_ker % value_ielem(ielem) % a(igaus) = gptur

                               do idime = 1,ndime
                                  prope_ker % grval_ielem(ielem) % a(idime,igaus) = 0.0_rp
                               enddo
                            end do
                         end if
                      end if
                   end do

                   prope_ker % kfl_nedsm = 1

                   do iboun = 1,nboun
                      pblty = abs(ltypb(iboun))
                      pnodb = nnode(pblty)
                      ielem = lelbo(iboun)

                      if( lmate(ielem) == imate ) then
                         pelty = ltype(ielem)
                         if( pelty > 0 ) then
                            pgaub = ngaus(pblty)
                            do igaub = 1,pgaub
                               gptur = const * lscale * lscale * tscale              ! D = C_D * L**2 / T
                               prope_ker % value_iboun(iboun) % a(igaub) =  gptur
                            end do
                         end if
                      end if
                   end do

                else if( prope_ker % wlaws(imate) == 'GPSUT' ) then
                   !----------------------------------------------------------------------
                   !
                   !                                                         T0 + S
                   ! Sutherland's law in Gauss points : mu/mu_0 = (T/T0)^1.5 --------
                   !                                                         T  + S
                   !
                   !----------------------------------------------------------------------

                   mu0   = prope_ker % rlaws(1,imate) ! reference value
                   T0    = prope_ker % rlaws(2,imate) ! referentemperature
                   ! Sutherland Temperature S0
                   S0    = prope_ker % rlaws(3,imate)
                   ! if not given, load default value (air)
                   if (S0.lt. epsilon(1.0_rp)) S0    = 110.0_rp

                   do ielem = 1,nelem
                      if( lmate(ielem) == imate ) then
                         pelty = ltype(ielem)
                         if( pelty > 0 ) then
                            pgaus = ngaus(pelty)
                            pnode = lnnod(ielem)
                            do inode = 1,pnode
                               ipoin = lnods(inode,ielem)
                               eltem(inode) = tempe(ipoin,which_time)
                               do idime = 1,ndime
                                  elcod(idime,inode) = coord(idime,ipoin)
                               end do
                            end do
                            call elmcar(&
                                 pnode,pgaus,0_ip,elmar(pelty) % weigp,elmar(pelty) % shape,&
                                 elmar(pelty) % deriv,elmar(pelty) % heslo,elcod,gpvol,gpcar,&
                                 dummr,ielem)
                            gp_temperature=0.0_rp
                            gp_grad_tem=0.0_rp
                            do igaus = 1,pgaus
                               do inode = 1,pnode
                                  ipoin = lnods(inode,ielem)
                                  gp_temperature(igaus) = gp_temperature(igaus) + eltem(inode) * elmar(pelty) % shape(inode,igaus)
                                  do idime = 1,ndime
                                     gp_grad_tem(idime,igaus) = gp_grad_tem(idime,igaus) + eltem(inode) * gpcar(idime,inode,igaus)
                                  end do
                               end do
                               if (gp_temperature(igaus) <= epsilon(1.0_rp) ) then
                                  ! For initialization we default to reference temperature
                                  gp_temperature(igaus) = T0
                                  gp_grad_tem(1:ndime,igaus) = 0.0_rp
                               endif
                               mu =  mu0 * (gp_temperature(igaus)/T0)**1.5_rp * (T0+S0) / (gp_temperature(igaus)+S0)
                               prope_ker % value_ielem(ielem) % a(igaus) = mu
                               do idime = 1,ndime
                                  prope_ker % grval_ielem(ielem) % a(idime,igaus) = gp_grad_tem(idime,igaus) * &
                                       (3.0_rp*S0 + gp_temperature(igaus)) &
                                       / (2.0_rp * gp_temperature(igaus)* (gp_temperature(igaus)+S0))*mu
                               enddo

                               ! derivative of mu (at gauss points) respect to nodal tempe
                               do inode = 1,pnode
                                  prope_ker % drval_ielem(ielem) % a(inode,igaus) = elmar(pelty) % shape(inode,igaus) * &
                                       (3.0_rp*S0 + gp_temperature(igaus)) &
                                       / (2.0_rp * gp_temperature(igaus)* (gp_temperature(igaus)+S0))*mu
                               enddo

                               ! derivative of grad(mu) (at gauss points) respect to nodal tempe
                               do inode = 1,pnode
                                  do idime = 1, ndime
                                     prope_ker % gdval_ielem(ielem) % a(idime, inode,igaus) = gpcar(idime,inode,igaus) * &
                                          (3.0_rp*S0 + gp_temperature(igaus)) &
                                          / (2.0_rp * gp_temperature(igaus)* (gp_temperature(igaus)+S0))*mu + &
                                          gp_grad_tem(idime,igaus) * elmar(pelty) % shape(inode,igaus) * &
                                          ( (3.0_rp*S0 + gp_temperature(igaus))**2 &
                                          / (2.0_rp * gp_temperature(igaus)* (gp_temperature(igaus)+S0))**2*mu + &
                                          ( -2*gp_temperature(igaus)**2 -12*gp_temperature(igaus)*S0 -6*S0**2 ) / &
                                          (4*gp_temperature(igaus)**2*(gp_temperature(igaus)+S0)**2 ) * mu )
                                  enddo
                               enddo

                            end do
                         end if
                      end if
                   end do
                   prope_ker % kfl_nedsm = 1

                   do iboun = 1,nboun
                      pblty = abs(ltypb(iboun))
                      pnodb = nnode(pblty)
                      ielem = lelbo(iboun)
                      if( lmate(ielem) == imate ) then
                         pelty = ltype(ielem)
                         if( pelty > 0 ) then
                            pgaub = ngaus(pblty)
                            do inodb = 1,pnodb
                               eltem(inodb) = tempe(lnodb(inodb,iboun),1)
                            end do
                            do igaub = 1,pgaub
                               gbtem = 0.0_rp
                               do inodb = 1,pnodb
                                  gbtem = gbtem + eltem(inodb) * elmar(pblty) % shape(inodb,igaub)
                               end do
                               if (gbtem <= epsilon(1.0_rp) ) then
                                  ! For initialization we default to reference temperature
                                  gbtem = 300.0_rp
                               endif
                               prope_ker % value_iboun(iboun) % a(igaub) = &
                                    mu0 * (gbtem/T0)**1.5_rp * (T0+S0) / (gbtem+S0)
                            end do
                         end if
                      end if
                   end do
#ifndef PROPER_PRIVATE_OFF
                else if( prope_ker % wlaws(imate) == 'SMAGO' ) then

                   !----------------------------------------------------------------------
                   !
                   !
                   !  SMAGORINSKY MODEL nu_t = (C_s * Vol)^2 * sqrt(2*Sij*Sij)
                   !
                   !
                   !----------------------------------------------------------------------

                   call ker_proper_element_operations(&
                        prope_ker,imate,&
                        ker_proper_turmu_smagorinsky_init,&
                        ker_proper_turmu_smagorinsky_operations)

                else if( prope_ker % wlaws(imate) == 'SIGMA' ) then

                   !----------------------------------------------------------------------
                   !
                   !
                   !  Sigma - Baya et al., 2010, Proceeding of teh Summer
                   !  Program
                   !
                   !----------------------------------------------------------------------

                   call ker_proper_element_operations(&
                        prope_ker,imate,&
                        ker_proper_turmu_sigma_init,&
                        ker_proper_turmu_sigma_operations)

                else if( prope_ker % wlaws(imate) == 'WALE ' ) then

                   !----------------------------------------------------------------------
                   !
                   !  WALE  - Wall-adapting local eddy-viscosity - Nicoud-Ducros 99
                   !
                   !----------------------------------------------------------------------

                   call ker_proper_element_operations(&
                        prope_ker,imate,&
                        ker_proper_turmu_wale_init,&
                        ker_proper_turmu_wale_operations)

                else if( prope_ker % wlaws(imate) == 'VRMAN' ) then

                   !----------------------------------------------------------------------
                   !
                   !  VRMAN - Vreman SGS model static version
                   !
                   !----------------------------------------------------------------------

                   call ker_proper_element_operations(&
                        prope_ker,imate,&
                        ker_proper_turmu_vreman_init,&
                        ker_proper_turmu_vreman_operations)

                else if( prope_ker % wlaws(imate) == 'ILSA ' ) then

                   !---------------------------------------------------------------------------------
                   !
                   !  ILSA
                   !  Dynamic subfilter-scale stress model for large-eddy simulations,
                   !  A. Rouhi, U. Piomelli and B.J. Geurts, PHYSICAL REVIEW FLUIDS 1, 044401 (2016)
                   ! 
                   !---------------------------------------------------------------------------------

                   call ker_proper_element_operations(&
                        prope_ker,imate,&
                        ker_proper_turmu_ilsa_init,&
                        ker_proper_turmu_ilsa_operations)

                else if( prope_ker % wlaws(imate) == 'TKESG' ) then

                   !----------------------------------------------------------------------
                   !
                   !  TKESGS - TKE SGS model, needs to be coupled with turbul, to solve an extra equation for TKE
                   !
                   !----------------------------------------------------------------------

                   call ker_proper_element_operations(&
                        prope_ker,imate,&
                        ker_proper_turmu_tkesgs_init,&
                        ker_proper_turmu_tkesgs_operations)

                else if( prope_ker % wlaws(imate) == 'FIELD' ) then

                   !----------------------------------------------------------------------
                   !
                   !  Property given by a field
                   !
                   !----------------------------------------------------------------------

                   icomp = int(prope_ker % rlaws(1,imate),ip)
                   
                   if( kfl_field(2,icomp) == NPOIN_TYPE ) then
                      
                      do ipoin = 1,npoin
                         prope_ker % value_ipoin(ipoin) = xfiel(icomp) % a(1,ipoin,1)
                      end do
                      
                   else if( kfl_field(2,icomp) == NELEM_TYPE ) then

                      ncoun = ndime*ndime
                      do ielem = 1,nelem
                         pelty = ltype(ielem)
                         do igaus = 1, ngaus(pelty)
                            do icoun = 1,ncoun
                               prope_ker % value_ielem(ielem) % a((igaus-1)*ncoun+icoun) = xfiel(icomp) % a(icoun,ielem,1)
                            end do
                         end do
                      end do
                      
                   else
                      call runend('MOD_KER_UPDPRO: WROING FIELD')
                   end if

                 else if( prope_ker % wlaws(imate) == 'WVSTD' ) then

                   !----------------------------------------------------------------------
                   !
                   !  Wall viscosity model
                   !
                   !----------------------------------------------------------------------

                     if( kfl_noslw_ker /= 0_ip ) then

                        const = prope_ker % rlaws(1,imate)
                        do ielem = 1,nelem
                           if( lmate(ielem) == imate .and. kfl_nswel_ker(ielem) > 0_ip .and. ltype(ielem) > 0 ) then
                              pelty = ltype(ielem)
                              pgaus = ngaus(pelty)
                              pnode = lnnod(ielem)
                              do inode = 1,pnode
                                 ipoin = lnods(inode,ielem)
                                 elcod(1:ndime,inode) = coord(1:ndime,ipoin)
                                 elvel(1:ndime,inode) = advec(1:ndime,ipoin,1)
                              end do
                              call elmcar(&
                                   pnode,pgaus,0_ip,elmar(pelty) % weigp,elmar(pelty) % shape,&
                                   elmar(pelty) % deriv,elmar(pelty) % heslo,elcod,gpvol,gpcar,&
                                   dummr,ielem)

                              call ker_proper('DENSI','PGAUS',igaus,ielem,gpden,pnode,pgaus,elmar(pelty) % shape,gpcar)
                              call ker_proper('VISCO','PGAUS',igaus,ielem,gpvis,pnode,pgaus,elmar(pelty) % shape,gpcar)

                              call ker_proper_turmu_walvi(ielem,pnode,pgaus,gpcar,elvel,gpvis,gpden,gpvis_nsw)
                              prope_ker % value_ielem(ielem) % a(1:pgaus) = gpvis_nsw(1:pgaus)
                           else
                              prope_ker % value_ielem(ielem) % a(:) = 0.0_rp
                           end if
                        end do

                        prope_ker % kfl_nedsm = 1

                        do iboun = 1,nboun
                           prope_ker % value_iboun(iboun) % a(:) =  0.0_rp
                        end do

                     end if
#endif
                else if( prope_ker % wlaws(imate) == 'MIXIN' ) then

                   !----------------------------------------------------------------------
                   !
                   !
                   !  ALGEBRAIC MIXING LENGTH MODEL -- nu_t = (kappa * yplus)^2 * sqrt(2*Sij*Sij) * (1-exp(yplus/A))^2
                   !
                   !
                   !----------------------------------------------------------------------

                   const = prope_ker % rlaws(1,imate)  ! Here comes A+ = 26 
                   !
                   ! 1) Obtain tange - the tangential componenet of the shear stress at the wall - idem postpr
                   !
                   ! TANGE: Tangential stress
                   !
                   if( INOTMASTER ) then
                      call memgen(zero,ndime,npoin)
                      call runend (' can not call nsi_outtan here -- need to find a solution') 
                      !                      call nsi_outtan()              ! generates non zero values in gevec(idime,ipoin)
                   end if
                   !
                   ! 2) Extrapolate it to the interior
                   ! Allocate and obtain tau_wall_boun  -  the shear stress at wall nodes
                   !
                   kpoin = wallcoupling_extr_boun % wet % npoin_wet  ! temporary use of kpoin
                   call memory_alloca(mem_modul(1:2,modul),'tau_wall_boun','mod_ker_proper',tau_wall_boun, ndime, kpoin)

                   call COU_GET_INTERPOLATE_POINTS_VALUES(gevec,tau_wall_boun,wallcoupling_extr_boun)    ! tau_wall_boun(kount)  , tau_wall is the nodal shear stress
                   !
                   ! The extrapolated shear stress at an interior node is obtained - for boundary nodes I already have the correct value in gevec
                   !
                   kpoin =0
                   do ipoin= 1,npoin
                      if(is_interior(ipoin) /= 0_ip) then ! interior
                         kpoin = kpoin + 1
                         gevec(1:ndime,ipoin) = tau_wall_boun( 1:ndime, ptb_to_use(kpoin) )
                      end if
                   end do
                   call memory_deallo(mem_modul(1:2,modul),'tau_wall_boun','mod_ker_proper',tau_wall_boun)   ! ojo tal ves no hacerlo cada vez
                   !
                   ! 3) Prepare everything that is needed for ker_turbul(6_ip
                   !
                   do ielem = 1,nelem
                      if( lmate(ielem) == imate ) then
                         pelty = ltype(ielem)
                         if( pelty > 0 ) then
                            pgaus = ngaus(pelty)
                            pnode = lnnod(ielem)
                            do inode = 1,pnode
                               ipoin = lnods(inode,ielem)
                               do idime = 1,ndime
                                  elcod(idime,inode) = coord(idime,ipoin)
                                  elvel(idime,inode) = advec(idime,ipoin,which_time)
                                  eltan(idime,inode) = gevec(idime,ipoin)
                               end do
                               elwal(inode) = walld(ipoin)
                            end do
                            call elmcar(&
                                 pnode,pgaus,0_ip,elmar(pelty) % weigp,elmar(pelty) % shape,&
                                 elmar(pelty) % deriv,elmar(pelty) % heslo,elcod,gpvol,gpcar,&
                                 dummr,ielem)

                            gpvel = 0.0_rp
                            gptan = 0.0_rp
                            gpgve = 0.0_rp
                            gpwal = 0.0_rp
                            !
                            ! HLENG and TRAGL at center of gravity
                            !
                            call elmlen(&
                                 ndime,pnode,elmar(pelty)%dercg,tragl,elcod,&
                                 hnatu(pelty),hleng)

                            gpwal2(:) = 0.0_rp
                            do igaus = 1,pgaus
                               do inode = 1,pnode
                                  gpwal2(igaus) = gpwal2(igaus) + elmar(pelty) % shape(inode,igaus) * elwal(inode)
                                  do idime = 1,ndime
                                     gpvel(idime,igaus) = gpvel(idime,igaus) + elvel(idime,inode) * elmar(pelty) % shape(inode,igaus)
                                     gptan(idime,igaus) = gptan(idime,igaus) + eltan(idime,inode) * elmar(pelty) % shape(inode,igaus)
                                     do jdime = 1,ndime
                                        gpgve(jdime,idime,igaus) = gpgve(jdime,idime,igaus) + elvel(idime,inode) * gpcar(jdime,inode,igaus)
                                     end do
                                  end do
                               end do

                               velmo = 0.0_rp
                               tanmo = 0.0_rp
                               do idime= 1,ndime
                                  velmo = velmo + gpvel(idime,igaus) * gpvel(idime,igaus)
                                  tanmo = tanmo + gptan(idime,igaus) * gptan(idime,igaus)
                               end do

                               velmo = sqrt(velmo)
                               tanmo = sqrt(tanmo)

                               call ker_proper('DENSI','IGAUS',igaus,ielem,gpden,pnode,pgaus,elmar(pelty) % shape,gpcar)
                               call ker_proper('VISCO','IGAUS',igaus,ielem,gpvis,pnode,pgaus,elmar(pelty) % shape,gpcar)
                               !
                               ! Computing turbulent viscosity nu_t at gauss points
                               ! according to mixing length model: gptur
                               !
                               call ker_turbul(6_ip,gpwal,velmo,const,gpvis(1),gpgve(:,:,igaus),hleng,gptur,gpden(1),tanmo)   

                               prope_ker % value_ielem(ielem) % a(igaus) = gptur

                               do idime = 1,ndime
                                  prope_ker % grval_ielem(ielem) % a(idime,igaus) = 0.0_rp
                               enddo
                            end do
                         end if
                      end if
                   end do

                   prope_ker % kfl_nedsm = 1

                   do iboun = 1,nboun
                      pblty = abs(ltypb(iboun))
                      pnodb = nnode(pblty)
                      ielem = lelbo(iboun)
                      if( lmate(ielem) == imate ) then
                         pelty = ltype(ielem)
                         pnode = nnode(pelty)
                         do inode = 1,pnode
                            ipoin = lnods(inode,ielem)
                            do idime = 1,ndime
                               elcod(idime,inode) = coord(idime,ipoin)
                            end do
                         end do
                         !
                         ! HLENG and TRAGL at center of gravity
                         !
                         call elmlen(&
                              ndime,pnode,elmar(pelty)%dercg,tragl,elcod,&
                              hnatu(pelty),hleng)

                         if( pelty > 0 ) then
                            pgaub = ngaus(pblty)
                            do inodb = 1,pnodb
                               do idime = 1,ndime
                                  elvel(idime,inodb) = advec(idime,lnodb(inodb,iboun),which_time)
                                  eltan(idime,inodb) = gevec(idime,lnodb(inodb,iboun))
                               end do
                               elwal(inodb) = walld(lnodb(inodb,iboun))
                            end do
                            do igaub = 1,pgaub
                               gbvel = 0.0_rp
                               gbtan = 0.0_rp
                               gbgve = 0.0_rp
                               gbwal = 0.0_rp
                               do inodb = 1,pnodb
                                  do idime = 1,ndime
                                     gbvel(idime) = gbvel(idime) + elvel(idime,inodb) * elmar(pblty) % shape(inodb,igaub)
                                     gbtan(idime) = gbvel(idime) + eltan(idime,inodb) * elmar(pblty) % shape(inodb,igaub)
                                     do jdime = 1,ndime
                                        gbgve(jdime,idime) = gbgve(jdime,idime) + elvel(idime,inodb) * gpcar(jdime,inodb,igaub)
                                     end do
                                     gbwal = gbwal + elmar(pblty) % shape(inodb,igaub) * elwal(inodb)
                                  end do
                               end do

                               velmo = 0.0_rp
                               tanmo = 0.0_rp
                               do idime= 1,ndime
                                  velmo = velmo + gbvel(idime) * gbvel(idime)
                                  tanmo = tanmo + gbtan(idime) * gbtan(idime)
                               end do

                               velmo = sqrt(velmo)
                               tanmo = sqrt(tanmo)

                               call ker_proper('DENSI','IGAUB',igaub,iboun,gbden)
                               call ker_proper('VISCO','IGAUB',igaub,iboun,gbvis)
                               !
                               ! Computing turbulent viscosity nu_t at gauss points
                               ! according to mixing length model: gptur
                               !
                               call ker_turbul(6_ip,gbwal,velmo,const,gbvis(1),gbgve(:,:),hleng,gptur,gbden(1),tanmo)

                               prope_ker % value_iboun(iboun) % a(igaub) =  gptur

                            end do
                         end if
                      end if
                   end do

                   if( iavec == 1 ) call memgen(2_ip,one,one)
                   nullify(gevec)

                else if( prope_ker % wlaws(imate) == 'SSTKO' ) then

                   !----------------------------------------------------------------------
                   !
                   !  Turbulent viscosity calculated at gauss points for SST K-Omega model
                   !
                   !----------------------------------------------------------------------

                   do ielem = 1,nelem
                      if( lmate(ielem) == imate ) then
                         pelty = ltype(ielem)
                         if( pelty > 0 ) then
                            pgaus = ngaus(pelty)
                            pnode = lnnod(ielem)
                            do inode = 1,pnode
                               ipoin = lnods(inode,ielem)
                               eltur(1:2,inode) = untur(1:2,ipoin,which_time)
                               elwal(inode) = walld(ipoin)
                               elvel(1:ndime, inode) = advec(1:ndime, ipoin,1)
                               do idime = 1,ndime
                                  elcod(idime,inode) = coord(idime,ipoin)
                               end do
                            end do
                            call elmcar(&
                                 pnode,pgaus,0_ip,elmar(pelty) % weigp,elmar(pelty) % shape,&
                                 elmar(pelty) % deriv,elmar(pelty) % heslo,elcod,gpvol,gpcar,&
                                 dummr,ielem)
                            do igaus = 1,pgaus
                               kinen = 0.0_rp
                               omega = 0.0_rp
                               gpwal = 0.0_rp

                               do idime = 1,ndime
                                  do jdime = 1,ndime
                                     gpgve(jdime,idime, igaus) = 0.0_rp
                                  end do
                               end do
                               do inode = 1,pnode
                                  !gpsha(inode, igaus) =  elmar(pelty) % shape(inode,igaus)
                                  kinen = kinen + elmar(pelty) % shape(inode,igaus) * eltur(1,inode)
                                  omega = omega + elmar(pelty) % shape(inode,igaus) * eltur(2,inode)
                                  gpwal = gpwal + elmar(pelty) % shape(inode,igaus) * elwal(inode)
                                  do idime = 1,ndime
                                     do jdime = 1,ndime
                                        gpgve(jdime,idime,igaus) = gpgve(jdime,idime,igaus) &
                                             + gpcar(jdime,inode,igaus) * elvel(idime,inode)
                                     end do
                                  end do
                               end do
                               if (kinen <= epsilon(1.0_rp) ) kinen = 1.0e10_rp
                               if (omega <= epsilon(1.0_rp) ) omega = 1.0e10_rp
                               !
                               ! ERROR in calculation of gpvor, gpvor = 2 W_ij W_ij where W_ij = sqrt(0.5 (dui_dxj - duj_dxi))
                               !
                               !                                gpvor = 0.0_rp   ! 2 S_ij : S_ij
                               !                                do idime = 1,ndime
                               !                                   do jdime = 1,ndime
                               !                                      gpvor = gpvor + gpgve(idime,jdime,igaus) &    ! 2 S_ij : S_ij
                               !                                           *(gpgve(idime,jdime,igaus)        &
                               !                                           +gpgve(jdime,idime,igaus))
                               !                                   end do
                               !                                end do

                               gpvor_mat = 0.0_rp
                               do idime = 1,ndime
                                  do jdime = 1,ndime
                                     gpvor_mat(idime,jdime) = 0.5_rp*(gpgve(jdime, idime,igaus) - gpgve(idime, jdime,igaus))
                                  end do
                               end do
                               gpvor = 0.0_rp
                               do idime = 1,ndime
                                  do jdime = 1,ndime
                                     gpvor = gpvor + 2.0_rp * gpvor_mat(idime,jdime) * gpvor_mat(idime,jdime)
                                  end do
                               end do

                               call ker_proper('DENSI','IGAUS',igaus,ielem,gpden,pnode,pgaus,elmar(pelty) % shape,gpcar)
                               call ker_proper('VISCO','IGAUS',igaus,ielem,gpvis,pnode,pgaus,elmar(pelty) % shape,gpcar)
                               gpvor = sqrt(gpvor)
                               a1 =0.31_rp
                               gpnut = 0.0_rp
                               if (gpwal.gt.1.0e-8_rp) then
                                  F2 = max( 2.0_rp *sqrt(max(0.0_rp,kinen)) &
                                       /(0.09_rp*omega*gpwal), 500.0_rp * gpvis(1) / (gpden(1)*gpwal*gpwal*omega))
                                  !                                   F2 = max( 2.0_rp *sqrt(max(0.0_rp,kinen)) /(0.09_rp*omega*gpwal), 500.0_rp * gpvis(1) / (gpwal*gpwal*omega)) !!!!!
                                  F2 = tanh( F2 * F2 )
                                  as = max( a1 * omega, gpvor * F2 )
                                  gpnut = a1*kinen/as
                               end if
                               ! check the limits
                               gpnut = max( gpnut, 0.1_rp*gpvis(1)/gpden(1) )
                               gpnut = min( gpnut,1.0_rp * 1.0e20_rp/gpden(1) )
                               prope_ker % value_ielem(ielem) % a(igaus) = gpnut

                               !
                               ! derivative of gpmut (at gauss points) respect to nodal kin and ome
                               ! arg2 = max (c1 , c2)                 F2 = TANH(arg2 * arg2)
                               !
                               ! dF2/dk = 2*arg2*d(arg2)/dk * (1-F2*F2)
                               ! dF2/dw = 2*arg2*d(arg2)/dw * (1-F2*F2)
                               !
                               do inode=1,pnode
                                  dgpkin(inode) = elmar(pelty) % shape(inode,igaus)
                                  dgpome(inode) = elmar(pelty) % shape(inode,igaus)
                               enddo

                               arg2 = max( 2.0_rp *sqrt(max(0.0_rp,kinen)) &
                                    /(0.09_rp*omega*gpwal), 500.0_rp * gpvis(1) / (gpden(1)*gpwal*gpwal*omega))

                               ! calculation of d(arg2)/dk and d(arg2)/dw
                               darg2_dkin = 0.0_rp
                               darg2_dome = 0.0_rp
                               if (kinen > 0.0_rp) then
                                  if (2.0_rp *sqrt(max(0.0_rp,kinen)) &
                                       /(0.09_rp*omega*gpwal) > 500.0_rp * gpvis(1) / (gpden(1)*gpwal*gpwal*omega) ) then
                                     do inode=1,pnode
                                        darg2_dkin(inode) = dgpkin(inode)/(0.09_rp*omega*gpwal*sqrt(max(0.0_rp,kinen)))
                                        darg2_dome(inode) = -2.0_rp*sqrt(max(0.0_rp,kinen))*dgpome(inode)/(0.09_rp*omega*omega*gpwal)
                                     enddo
                                  else
                                     do inode=1,pnode
                                        darg2_dkin(inode) = 0.0_rp
                                        darg2_dome(inode) = -500.0_rp * gpvis(1) * dgpome(inode)/ (gpden(1)*gpwal*gpwal*omega*omega)
                                     enddo
                                  endif
                               else
                                  if (0.0_rp > 500.0_rp * gpvis(1) / (gpden(1)*gpwal*gpwal*omega) ) then
                                     do inode=1,pnode
                                        darg2_dkin(inode) = 0.0_rp
                                        darg2_dome(inode) = 0.0_rp
                                     enddo
                                  else
                                     do inode=1,pnode
                                        darg2_dkin(inode) = 0.0_rp
                                        darg2_dome(inode) = -500.0_rp * gpvis(1) * dgpome(inode)/ (gpden(1)*gpwal*gpwal*omega*omega)
                                     enddo
                                  endif
                               endif

                               ! calculation of dF2/dk and dF2/dw
                               do inode = 1,pnode
                                  dF2_dkin(inode) = 2.0_rp*arg2*darg2_dkin(inode)*(1.0_rp-F2*F2)
                                  dF2_dome(inode) = 2.0_rp*arg2*darg2_dome(inode)*(1.0_rp-F2*F2)
                               enddo

                               ! calculation of d(nu_t)/dk and d(nu_t)/dw
                               dgpnut_dkin = 0.0_rp
                               dgpnut_dome = 0.0_rp
                               if (a1 * omega > gpvor * F2) then
                                  do inode = 1,pnode
                                     dgpnut_dkin(inode) = dgpkin(inode)/omega
                                     dgpnut_dome(inode) = -kinen*dgpome(inode)/(omega*omega)
                                  enddo
                               else
                                  do inode = 1,pnode
                                     dgpnut_dkin(inode) = a1*dgpkin(inode)/(gpvor * F2) - a1*kinen*dF2_dkin(inode)/(gpvor*F2*F2)
                                     dgpnut_dome(inode) = -a1*kinen*dF2_dome(inode)/(gpvor*F2*F2)
                                  enddo
                               endif

                               ! check limits
                               if ( gpden(1)*gpnut > 1.0e20_rp .or. gpden(1)*gpnut < 0.1_rp*gpvis(1) ) then
                                  dgpnut_dkin = 0.0_rp
                                  dgpnut_dome = 0.0_rp
                               endif
                               ! pass d(nu_t)/dk and d(nu_t)/dw to the kernel variables
                               do inode = 1,pnode
                                  prope_ker % drval_tur_ielem(ielem) % a(1,inode,igaus) = gpden(1) * dgpnut_dkin(inode)
                                  prope_ker % drval_tur_ielem(ielem) % a(2,inode,igaus) = gpden(1) * dgpnut_dome(inode)
                               end do

                               !
                               ! derivative of gpmut (at gauss points) respect to nodal veloc
                               ! arg2 = max (c1 , c2)                 F2 = TANH(arg2 * arg2)
                               ! dF2/dU = 2*arg2*d(arg2)/du * (1-F2*F2)
                               ! d(arg2)/du = d(gpvor)/du
                               ! gpvor = sqrt[ (du1/dx2 - du2_dx1)**2 + (du1/dx3 - du3/dx1)**2 + (du3/dx2 - du2/dx3)**2 ]
                               !
                               if (ndime == 2) then
                                  do inode = 1,pnode
                                     dgpvor_du(inode,1) = gpcar(2,inode,igaus)*(gpgve(2,1,igaus)-gpgve(1,2,igaus))/gpvor
                                     dgpvor_du(inode,2) = -gpcar(1,inode,igaus)*(gpgve(2,1,igaus)-gpgve(1,2,igaus))/gpvor
                                  enddo
                               elseif (ndime == 3) then
                                  do inode = 1,pnode
                                     dgpvor_du(inode,1) = ( gpcar(2,inode,igaus)*(gpgve(2,1,igaus)-gpgve(1,2,igaus)) &
                                          + gpcar(3,inode,igaus)*(gpgve(3,1,igaus)-gpgve(1,3,igaus)) ) /gpvor
                                     dgpvor_du(inode,2) = (-gpcar(1,inode,igaus)*(gpgve(2,1,igaus)-gpgve(1,2,igaus)) &
                                          + gpcar(3,inode,igaus)*(gpgve(3,2,igaus)-gpgve(2,3,igaus)) ) /gpvor
                                     dgpvor_du(inode,3) = ( gpcar(1,inode,igaus)*(gpgve(3,1,igaus)-gpgve(1,3,igaus)) &
                                          - gpcar(2,inode,igaus)*(gpgve(3,2,igaus)-gpgve(2,3,igaus)) ) /gpvor
                                  enddo
                               endif

                               ! calculation of d(nu_t)/du
                               dgpnut_du = 0.0_rp
                               if (a1 * omega < gpvor * F2) then
                                  do inode = 1,pnode
                                     do idime = 1,ndime
                                        dgpnut_du(inode,idime) = -a1*kinen*dgpvor_du(inode,idime)/(F2*gpvor*gpvor)
                                     enddo
                                  enddo
                               endif

                               ! check limits
                               if ( gpden(1)*gpnut > 1.0e20_rp .or. gpden(1)*gpnut < 0.1_rp*gpvis(1) ) then
                                  dgpnut_du = 0.0_rp
                               endif

                               ! pass d(nu_t)/du to the kernel variables
                               do inode = 1,pnode
                                  do idime = 1,ndime
                                     prope_ker % drval_vel_ielem(ielem) % a(idime,inode,igaus) = gpden(1) * dgpnut_du(inode,idime)
                                  enddo
                               end do

                            end do ! igaus

                            ! gradient
                            prope_ker % grval_ielem(ielem) % a = 0.0_rp

                         end if
                      end if
                   end do
                else if( prope_ker % wlaws(imate) == 'STDKE' ) then

                   !----------------------------------------------------------------------
                   !
                   !  Turbulent viscosity calculated at gauss points for standard K-Epsilon model
                   !
                   !----------------------------------------------------------------------

                   prope_ker % kfl_nedsm = 1 ! flag to project the property from elements to nodes in postprocess

                   if (kfl_logva==1) then
                      relat =1.0_rp
                   else
                      relat =1.0_rp
                   end if
                   ! Realizable constant depending on cmu value
                   if (kfl_kxmod_ker.ne.0) &
                        A0 = (1.0_rp - 3.0_rp*sqrt(0.5_rp*cmu_st))/cmu_st

                   element_loop:  do ielem = 1,nelem
                      if( lmate(ielem) == imate ) then
                         pelty = ltype(ielem)
                         if( pelty > 0 ) then
                            pgaus = ngaus(pelty)
                            pnode = lnnod(ielem)
                            do inode = 1,pnode
                               ipoin = lnods(inode,ielem)
                               eltur(1:2,inode) = untur(1:2,ipoin,which_time)
                               elvel(1:ndime, inode) = advec(1:ndime, ipoin,1)
                               elcod(1:ndime,inode) = coord(1:ndime,ipoin)    
                            end do
                            gauss_loop:  do igaus = 1,pgaus
                               kinen = 0.0_rp
                               epsil = 0.0_rp
                               do inode = 1,pnode
                                  kinen = kinen + eltur(1,inode) * elmar(pelty)%shape(inode,igaus)
                                  epsil = epsil + eltur(2,inode) * elmar(pelty)%shape(inode,igaus)
                               end do
                               if (kfl_regularization) then
                                  reguk  = regul_k(kinen)
                                  regue  = regul_e(epsil)
                               else
                                  reguk  = kinen
                                  regue  = epsil
                               end if
                               if (kfl_kxmod_ker==0.or.ittim==0) then ! k eps standard
                                  cmu = cmu_st

                               else if (kfl_kxmod_ker==2.or.kfl_kxmod_ker==3) then  ! Realizable or kefp models, modify cmu
                                  call elmca2(&
                                       pnode,pgaus,0_ip,elmar(pelty)%weigp,elmar(pelty)%shape,&
                                       elmar(pelty)%deriv,elmar(pelty)%heslo,elcod,gpvol,gpsha,&
                                       gpder,gpcar,gphes,ielem)

                                  gpgve(1:ndime,1:ndime,igaus) = 0.0_rp
                                  do inode = 1,pnode
                                     do idime = 1,ndime
                                        do jdime = 1,ndime
                                           gpgve(jdime,idime,igaus) = gpgve(jdime,idime,igaus) &
                                                + gpcar(jdime,inode,igaus) * elvel(idime,inode)
                                        end do
                                     end do
                                  end do
                                  divve =0.0_rp       ! divergence of u
                                  do idime = 1,ndime  ! symmetric gradient (strain tensor)
                                     divve = divve + gpgve(idime, idime, igaus)
                                     do jdime = 1,ndime
                                        simgr(idime,jdime) = 0.5_rp*(gpgve(idime, jdime,igaus) &
                                             +  gpgve(jdime, idime,igaus))
                                     end do
                                  end do
                                  seci4 = 0.0_rp
                                  W=0.0_rp
                                  uaste =0.0_rp
                                  do idime = 1,ndime !deviator, necessary to compute W
                                     simgr(idime, idime ) = simgr(idime, idime) - divve/3.0_rp 
                                  end do
                                  do idime = 1,ndime
                                     do jdime = 1,ndime
                                        seci4 = seci4 + simgr(idime, jdime) &        ! S_ij : S_ij
                                             *simgr(idime, jdime)
                                        uaste = uaste +gpgve(idime,jdime,igaus) &    ! D_ij : D_ij
                                             *(gpgve(idime,jdime,igaus))   

                                        do kdime =1, ndime
                                           W = W +  simgr(idime,jdime)* &
                                                simgr(jdime,kdime)* &
                                                simgr(kdime,idime)
                                        end do
                                     end do
                                  end do
                                  uaste = sqrt(uaste -divve*divve /3.0_rp)
                                  if (seci4.gt.epsilon(1.0_rp)) W = W/sqrt(seci4*seci4*seci4)

                                  Wsqr6 = sqrt(6.0_rp)*W     

                                  Wsqr6 = min(Wsqr6,  1.0_rp)  ! corrects limits
                                  Wsqr6 = max(Wsqr6, -1.0_rp)  ! it should not happen

                                  phi_r = 1.0_rp/3.0_rp*acos(Wsqr6)
                                  As_r= sqrt(6.0_rp)*cos(phi_r)

                                  ! calculates turb viscosity
                                  if (kfl_kxmod_ker==2) then ! REALIZABLE MODEL
                                     ! effective cmu
                                     cmu =1.0_rp/(A0+As_r*reguk/regue*uaste)
                                  else if (kfl_kxmod_ker==3) then ! KEFP
                                     Cr  = 4.5_rp
                                     f0  = Cr/(Cr-1.0_rp)
                                     f02 = 4.0_rp*f0*(f0-1.0_rp)
                                     sigmr = cmu_st *(uaste*reguk/regue)*(uaste*reguk/regue)
                                     ! effective cmu
                                     cmu   = cmu_st *2.0_rp*f0 /(1.0_rp+ sqrt(1.0_rp + f02*sigmr))
                                  end if

                                  !                                 Gpmut(igaus)= max(gpden(igaus)*cmu*kinen*kinen/epsil, gpvis(igaus))

                                  !                                  if (kfl_regularization) &
                                  !                                       gpmut(igaus)  = cmu*gpden(igaus)*reguk*reguk/regue

                               end if ! end ke model

                               ! density and viscosity
                               call ker_proper('DENSI','IGAUS',igaus,ielem,gpden,pnode,pgaus,elmar(pelty)%shape)
                               call ker_proper('VISCO','IGAUS',igaus,ielem,gpvis,pnode,pgaus,elmar(pelty)%shape)

                               if (kfl_logva==0) then

                                  if (regue.gt. 10.0*epsilon(1.0_rp)) then
                                     gpnut = cmu *  reguk*reguk/regue
                                  else
                                     gpnut =0.0_rp
                                  end if
                                  if (.not.kfl_regularization) &
                                       gpnut = max(gpnut,gpvis(1)/gpden(1))
                               else
                                  gpnut = cmu_st * exp(2.0_rp*reguk- regue)
                               end if
                               prope_ker % value_ielem(ielem) % a(igaus) = gpnut
                               !  USE RELAXATION
                               !                                    max(relat*gpden(1)*gpnut + (1.0_rp-relat)*  &
                               !                                    prope_ker % value_ielem(ielem) % a(igaus) , gpvis(1))

                            end do gauss_loop
                            ! gradients! convergence problems
                            prope_ker % grval_ielem(ielem) % a = 0.0_rp
                         end if ! pelty > 0
                      end if
                   end do element_loop
                else if( prope_ker % wlaws(imate) == 'KOMEG' ) then
                   !----------------------------------------------------------------------
                   !
                   !  Turbulent viscosity calculated at gauss points for standard K-Omega model
                   !
                   !----------------------------------------------------------------------

                   prope_ker % kfl_nedsm = 1 ! flag to project the property from elements to nodes in postprocess
                   alphast =1.0_rp           ! used for low Reynolds model
                   ielement_loop:  do ielem = 1,nelem
                      if( lmate(ielem) == imate ) then
                         pelty = ltype(ielem)
                         if( pelty > 0 ) then
                            pgaus = ngaus(pelty)
                            pnode = lnnod(ielem)
                            do inode = 1,pnode
                               ipoin = lnods(inode,ielem)
                               eltur(1:2,inode) = untur(1:2,ipoin,which_time)
                               elvel(1:ndime, inode) = advec(1:ndime, ipoin,1)
                               elcod(1:ndime,inode) = coord(1:ndime,ipoin)    
                            end do
                            
                            igauss_loop:  do igaus = 1,pgaus
                               kinen = 0.0_rp
                               omega = 0.0_rp
                               do inode = 1,pnode
                                  kinen = kinen + eltur(1,inode) * elmar(pelty)%shape(inode,igaus)
                                  omega = omega + eltur(2,inode) * elmar(pelty)%shape(inode,igaus)
                               end do
                               call ker_proper('DENSI','IGAUS',igaus,ielem,gpden,pnode,pgaus,elmar(pelty)%shape)
                               call ker_proper('VISCO','IGAUS',igaus,ielem,gpvis,pnode,pgaus,elmar(pelty)%shape)

                               if (kfl_kxmod_ker==3.or.kfl_kxmod_ker==4) then ! Wilcox 2006 model, introduction of production limiter limiting omega by below
                                  call elmca2(&
                                       pnode,pgaus,0_ip,elmar(pelty)%weigp,elmar(pelty)%shape,&
                                       elmar(pelty)%deriv,elmar(pelty)%heslo,elcod,gpvol,gpsha,&
                                       gpder,gpcar,gphes,ielem)

                                  gpgve(1:ndime,1:ndime,igaus) = 0.0_rp
                                  do inode = 1,pnode
                                     do idime = 1,ndime
                                        do jdime = 1,ndime
                                           gpgve(jdime,idime,igaus) = gpgve(jdime,idime,igaus) &
                                                + gpcar(jdime,inode,igaus) * elvel(idime,inode)
                                        end do
                                     end do
                                  end do

                                  divve =0.0_rp       ! divergence of u
                                  do idime = 1,ndime  
                                     divve = divve + gpgve(idime, idime, igaus)
                                  end do
                                  seci4 = 0.0_rp  ! = ! 2 S^D_ij : S_ij =2 S^D_ij : S^D_ij 
                                  do idime = 1,ndime                   
                                     do jdime = 1,ndime         
                                        ! 2 S_ij : S_ij
                                        seci4 = seci4 + gpgve(idime,jdime, igaus) * (gpgve(idime,jdime, igaus) + gpgve(jdime,idime, igaus))
                                     end do
                                  end do
                                  seci4 = seci4 -2.0_rp/3.0_rp*divve*divve
                                  if (kfl_kxmod_ker==4.and.omega.gt. 10.0_rp*epsilon(1.0_rp)) then ! low Reynolds number version
                                     ReT = gpden(1)*kinen/(gpvis(1)*omega)
                                     alphast =  ( 0.0708_rp /3.0_rp + ReT/6.0_rp)/(1.0_rp+ ReT/6.0_rp ) 
                                  end if
                                  ! modification of omega, only to calculate turb viscosity
                                  omega = max(omega, 7.0_rp/8.0_rp*sqrt(alphast*seci4/Cmu_st))

                               end if !Wilcox 2006 model
                          
                               if (omega.gt. 10.0_rp*epsilon(1.0_rp)) then
                                  gpnut = alphast*max(kinen/omega, 0.0_rp)
                               else
                                  gpnut = 0.0_rp
                               end if
                               gpnut = max(gpnut,1.0e-6_rp*gpvis(1)/gpden(1)) ! lower limiter

                               prope_ker % value_ielem(ielem) % a(igaus) = gpnut
                               ! this value is used in momentum equation
                            end do igauss_loop
                            ! gradients! convergence problems
                            prope_ker % grval_ielem(ielem) % a = 0.0_rp
                         end if ! pelty > 0
                      end if
                   end do ielement_loop

                else if( prope_ker % wlaws(imate) == 'SPALA' ) then

                   !----------------------------------------------------------------------
                   !
                   !  Turbulent viscosity calculated at gauss points for Spalart-Allmaras model
                   !
                   !----------------------------------------------------------------------
                   cv1 = 7.1_rp
                   do ielem = 1,nelem
                      if( lmate(ielem) == imate ) then
                         pelty = ltype(ielem)
                         if( pelty > 0 ) then
                            pgaus = ngaus(pelty)
                            pnode = lnnod(ielem)
                            do inode = 1,pnode
                               ipoin = lnods(inode,ielem)
                               eltur(1,inode) = untur(1,ipoin,which_time)
                               !                                eltur(1,inode) = untur_fix(ipoin)
                            end do
                            do igaus = 1,pgaus

                               gpnut = 0.0_rp
                               do inode = 1,pnode
                                  gpsha(inode,igaus) = elmar(pelty)%shape(inode,igaus)
                                  gpnut = gpnut + gpsha(inode, igaus)*eltur(1,inode)
                               end do
                               call ker_proper('DENSI','IGAUS',igaus,ielem,gpden,pnode,pgaus,gpsha)
                               call ker_proper('VISCO','IGAUS',igaus,ielem,gpvis,pnode,pgaus,gpsha)

                               gpnu  = gpvis(1)/gpden(1)
                               X     = gpnut/gpnu
                               fv1 = X*X*X/(X*X*X + cv1*cv1*cv1)
                               gpnut = fv1*gpnut

                               ! pass mu_t to the kernel variables
                               prope_ker % value_ielem(ielem) % a(igaus) = gpnut

                               ! calculation of d(mu_t)/d(nu_t)
                               do inode = 1,pnode
                                  dX_dnut = gpsha(inode,igaus)/gpnu
                                  dfv1_dnut = ( 3.0_rp*X*X*dX_dnut*(X*X*X + cv1*cv1*cv1) &
                                       - 3.0_rp*X*X*dX_dnut*X*X*X )/(X*X*X + cv1*cv1*cv1)**2.0_rp
                                  dgpmut_dnut(inode) = gpden(1)*dfv1_dnut*gpnut + gpden(1)*fv1*gpsha(inode, igaus)
                               enddo

                               ! pass d(mu_t)/d(nu_t) to the kernel variables
                               do inode = 1,pnode
                                  prope_ker % drval_tur_ielem(ielem) % a(1,inode,igaus) = dgpmut_dnut(inode)
                               end do

                            end do

                         end if
                      end if
                   end do


                else if( prope_ker % wlaws(imate) == 'LOWMA' ) then

                   !----------------------------------------------------------------------
                   !
                   ! Low Mach approximation, density = p(th) / R T
                   !
                   !----------------------------------------------------------------------

                   do kpoin = 1,nmatn(imate)
                      ipoin = lmatn(imate) % l(kpoin)

                      if (tempe(ipoin,which_time) /= 0.0_rp ) then
                         prope_ker % value_ipoin(ipoin) = prthe(which_time) /(gasco * tempe(ipoin,which_time))
                      else
                         prope_ker % value_ipoin(ipoin) = 0.0_rp
                      endif
                   end do

                else if( prope_ker % wlaws(imate) == 'LOWMG' ) then

                   !----------------------------------------------------------------------
                   !
                   ! Low Mach approximation, density = p(th) / R T IN GAUSS POINTS
                   !
                   !----------------------------------------------------------------------
                   do ielem = 1,nelem
                      if( lmate(ielem) == imate ) then
                         pelty = ltype(ielem)
                         if( pelty > 0 ) then
                            pgaus = ngaus(pelty)
                            pnode = lnnod(ielem)
                            if (associated(tesgs)) then
                               do igaus = 1,pgaus
                                  gptem = 0.0_rp
                                  do inode = 1,pnode
                                     ipoin = lnods(inode,ielem)
                                     gptem = gptem + tempe(ipoin,which_time) * elmar(pelty) % shape(inode,igaus)
                                  end do
                                  gptem = gptem+ tesgs(ielem)%a(1,igaus,1)
                                  if (gptem <= epsilon(1.0_rp) ) then
                                     ! For initialization we default to reference temperature
                                     gptem = 300.0_rp
                                  endif
                                  prope_ker % value_ielem(ielem) % a(igaus) = prthe(which_time)/gasco/gptem
                               end do
                            else
                               do igaus = 1,pgaus
                                  gptem = 0.0_rp
                                  do inode = 1,pnode
                                     ipoin = lnods(inode,ielem)
                                     gptem = gptem + tempe(ipoin,which_time) * elmar(pelty) % shape(inode,igaus)
                                  end do
                                  if (gptem <= epsilon(1.0_rp) ) then
                                     ! For initialization we default to reference temperature
                                     gptem = 300.0_rp
                                  endif
                                  prope_ker % value_ielem(ielem) % a(igaus) = prthe(which_time)/gasco/gptem
                               end do
                            end if

                         endif
                      end if
                   end do

                   do iboun = 1,nboun
                      pblty = abs(ltypb(iboun))
                      pnodb = nnode(pblty)
                      ielem = lelbo(iboun)
                      if( lmate(ielem) == imate ) then
                         pelty = ltype(ielem)
                         if( pelty > 0 ) then
                            pgaub = ngaus(pblty)
                            do inodb = 1,pnodb
                               eltem(inodb) = tempe(lnodb(inodb,iboun),1)
                            end do
                            do igaub = 1,pgaub
                               gbtem = 0.0_rp
                               do inodb = 1,pnodb
                                  gbtem = gbtem + eltem(inodb) * elmar(pblty) % shape(inodb,igaub)
                               end do
                               if (gbtem <= epsilon(1.0_rp) ) then
                                  ! For initialization we default to reference temperature
                                  gbtem = 300.0_rp
                               endif
                               prope_ker % value_iboun(iboun) % a(igaub) =  prthe(1)/gasco/gbtem
                            end do
                         end if
                      end if
                   end do

                else if( prope_ker % wlaws(imate) == 'KLOWM' .or. prope_ker % wlaws(imate) == 'TLOWM' ) then

                   !----------------------------------------------------------------------
                   !
                   ! Low Mach approximation, density = p(th) W / R T  where W is mean molecular density
                   !
                   ! Low Mach approximation for the syncronized CFI combustion model,
                   ! density is only updated in temper
                   !
                   !----------------------------------------------------------------------

                   do kpoin = 1,nmatn(imate)
                      ipoin = lmatn(imate) % l(kpoin)
                      if (tempe(ipoin,which_time) > 0.0_rp ) then
                         prope_ker % value_ipoin(ipoin) = prthe(which_time) &
                              * wmean(ipoin,which_time) /(gasco * tempe(ipoin,which_time))
                      else
                         prope_ker % value_ipoin(ipoin) = 0.0_rp
                      endif
                   end do

                else if( prope_ker % wlaws(imate) == 'SPRAY' ) then

                   !----------------------------------------------------------------------
                   !
                   ! Liquid gas density at gauss_points for ELSA model
                   !
                   !    rho = Y_L * rho_L + (1-Y_L) * rho_g     
                   !
                   !----------------------------------------------------------------------
                   iclas_Y = nspec - 1
                   rho_l = prope_ker % rlaws(1,imate)
                   rho_g = prope_ker % rlaws(2,imate)

                   do ielem = 1,nelem

                      if( lmate(ielem) == imate ) then
                         pelty = ltype(ielem)

                         if( pelty > 0 ) then
                            pgaus = ngaus(pelty)
                            pnode = lnnod(ielem)
                            do inode = 1,pnode
                               ipoin = lnods(inode,ielem)
                               do idime = 1,ndime
                                  elcod(idime,inode) = coord(idime,ipoin)
                               end do
                            end do

                            call elmcar(&
                                 pnode,pgaus,0_ip,elmar(pelty) % weigp,elmar(pelty) % shape,&
                                 elmar(pelty) % deriv,elmar(pelty) % heslo,elcod,gpvol,gpcar,&
                                 dummr,ielem)
                            call elmlen(&
                                 ndime,pnode,elmar(pelty)%dercg,tragl,elcod,&
                                 hnatu(pelty),hleng)

                            if(ndime == 3) then
                               delta = (hleng(1)*hleng(2)*hleng(3))**(1.0_rp/3.0_rp)
                            else  
                               delta = sqrt(hleng(1)*hleng(2))
                            endif

                            do igaus = 1,pgaus
                               gpY     = 0.0_rp
                               do inode = 1,pnode
                                  ipoin = lnods(inode,ielem)
                                  gpY = gpY + conce(ipoin,iclas_Y,which_time) * elmar(pelty) % shape(inode,igaus)
                               end do

                               !prope_ker % value_ielem(ielem) % a(igaus) = max(rho_g, min(rho_l,gpY * rho_l + (1.0_rp - gpY) * rho_g) )
                               gpY = max(0.0_rp,min(1.0_rp,gpY))
                               prope_ker % value_ielem(ielem) % a(igaus) = gpY * rho_l + (1.0_rp - gpY) * rho_g
                               !H = 0.0_rp
                               !if(2.0_rp*gpY - 1.0_rp > delta) then
                               !   H = 1.0_rp
                               !else if(2.0_rp*gpY-1.0_rp < -delta) then
                               !   H = 0.0_rp
                               !else
                               !   H = (2.0_rp*gpY - 1.0_rp + delta) / (2.0_rp*delta)
                               !end if
                               !prope_ker % value_ielem(ielem) % a(igaus) = rho_l * H + rho_g * (1.0_rp - H)

                            end do !igaus
                         end if

                      end if

                   end do

                   prope_ker % kfl_nedsm = 1
                   elY = 0.0_rp

                   !!DMM do iboun = 1,nboun
                   !!DMM    pblty = abs(ltypb(iboun))
                   !!DMM    pnodb = nnode(pblty)
                   !!DMM    ielem = lboel(pnodb+1,iboun)

                   !!DMM    if( lmate(ielem) == imate ) then
                   !!DMM       pelty = ltype(ielem)
                   !!DMM       if( pelty > 0 ) then
                   !!DMM          pgaub = ngaus(pblty)
                   !!DMM          do inodb = 1,pnodb
                   !!DMM             elY(inodb) = conce(lnodb(inodb,iboun),iclas_Y,which_time)
                   !!DMM          end do
                   !!DMM          do igaub = 1,pgaub
                   !!DMM             gbY = 0.0_rp
                   !!DMM             do inodb = 1,pnodb
                   !!DMM                gbY = gbY + elY(inodb) * elmar(pblty) % shape(inodb,igaub)
                   !!DMM             end do
                   !!DMM             prope_ker % value_iboun(iboun) % a(igaub) = gbY * rho_l + (1.0_rp - gbY) * rho_g
                   !!DMM          end do
                   !!DMM       end if
                   !!DMM    end if

                   !!DMM end do


                else if( prope_ker % wlaws(imate) == 'GKLOW' ) then

                   !----------------------------------------------------------------------
                   !
                   ! Low Mach approximation, density = p(th) W / R T at gauss_points where W is mean molecular density
                   !                         d(density)/dx = p(th) dW/dx / (R T) - p(th) W R/ (R T)^2 dT/dx
                   !
                   !----------------------------------------------------------------------
                   do ielem = 1,nelem
                      if( lmate(ielem) == imate ) then
                         pelty = ltype(ielem)
                         if( pelty > 0 ) then
                            pgaus = ngaus(pelty)
                            pnode = lnnod(ielem)
                            do inode = 1,pnode
                               ipoin = lnods(inode,ielem)
                               do idime = 1,ndime
                                  elcod(idime,inode) = coord(idime,ipoin)
                               end do
                            end do
                            call elmcar(&
                                 pnode,pgaus,0_ip,elmar(pelty) % weigp,elmar(pelty) % shape,&
                                 elmar(pelty) % deriv,elmar(pelty) % heslo,elcod,gpvol,gpcar,&
                                 dummr,ielem)

                            gp_grad_tem   = 0.0_rp
                            gp_grad_wmean = 0.0_rp

                            do igaus = 1,pgaus

                               gptem = 0.0_rp
                               gpwmean = 0.0_rp
                               do inode = 1,pnode
                                  ipoin = lnods(inode,ielem)
                                  gptem = gptem + tempe(ipoin,which_time) * elmar(pelty) % shape(inode,igaus)
                                  if  ( associated(wmean) ) then 
                                     gpwmean = gpwmean + wmean(ipoin,which_time) * elmar(pelty) % shape(inode,igaus)
                                  else
                                     gpwmean = gpwmean + 1.0_rp * elmar(pelty) % shape(inode,igaus)
                                  end if
                                  do idime = 1,ndime
                                     gp_grad_tem(idime,igaus)   = gp_grad_tem(idime,igaus) &
                                          + tempe(ipoin,which_time) * gpcar(idime,inode,igaus)
                                     if ( associated(wmean) ) then 
                                        gp_grad_wmean(idime,igaus) = gp_grad_wmean(idime,igaus) &
                                             + wmean(ipoin,which_time) * gpcar(idime,inode,igaus)
                                     else
                                        gp_grad_wmean(idime,igaus) = 0.0_rp
                                     end if
                                  end do
                               end do

                               if (associated(tesgs)) then
                                  gptem = gptem + tesgs(ielem)%a(1,igaus,1)
                               end if
                               if (gptem <= epsilon(1.0_rp) ) then
                                  ! For initialization we default to reference temperature
                                  gptem = 300.0_rp
                               end if
                               rho = prthe(which_time)*gpwmean/(gasco*gptem)
                               prope_ker % value_ielem(ielem) % a(igaus) = rho

                               do idime = 1,ndime
                                  prope_ker % grval_ielem(ielem) % a(idime,igaus) = &
                                       prthe(which_time) * gp_grad_wmean(idime,igaus)/(gasco*gptem) &
                                       - prthe(which_time)*gpwmean*gp_grad_tem(idime,igaus)/(gasco*gptem**2)
                               end do

                               ! derivative of rho (at gauss points) respect to nodal tempe
                               do inode = 1,pnode
                                  prope_ker % drval_ielem(ielem) % a(inode,igaus) = -elmar(pelty) % shape(inode,igaus) * rho/gptem
                               end do

                               ! derivative of d(rho)/dx (at gauss points) respect to nodal tempe
                               do inode = 1,pnode
                                  do idime = 1, ndime
                                     prope_ker % gdval_ielem(ielem) % a(idime, inode,igaus) = &
                                          - prthe(which_time) * gp_grad_wmean(idime,igaus)&
                                          /(gasco*gptem**2)*elmar(pelty) % shape(inode,igaus) &
                                          - rho/gptem * gpcar(idime,inode,igaus) &
                                          + 2.0_rp*rho/(gptem**2)*gp_grad_tem(idime,igaus) * elmar(pelty) % shape(inode,igaus)
                                  end do
                               end do

                            end do !igaus
                         end if
                      end if
                   end do
                   prope_ker % kfl_nedsm = 1

                   do iboun = 1,nboun
                      pblty = abs(ltypb(iboun))
                      pnodb = nnode(pblty)
                      ielem = lelbo(iboun)
                      if( lmate(ielem) == imate ) then
                         pelty = ltype(ielem)
                         if( pelty > 0 ) then
                            pgaub = ngaus(pblty)
                            do inodb = 1,pnodb
                               eltem(inodb) = tempe(lnodb(inodb,iboun),which_time)
                               if (kfl_coupl(ID_TEMPER,ID_CHEMIC) /= 0 ) then
                                  elwmean(inodb) = wmean(lnodb(inodb,iboun),which_time)
                               else
                                  elwmean(inodb) = 1.0_rp
                               end if
                            end do
                            do igaub = 1,pgaub
                               gbtem = 0.0_rp
                               gbwme = 0.0_rp
                               do inodb = 1,pnodb
                                  gbtem = gbtem + eltem(inodb) * elmar(pblty) % shape(inodb,igaub)
                                  gbwme = gbwme + elwmean(inodb) * elmar(pblty) % shape(inodb,igaub)
                               end do
                               if (gbtem <= epsilon(1.0_rp) ) then
                                  ! For initialization we default to reference
                                  ! temperature
                                  gbtem = 300.0_rp
                               endif
                               prope_ker % value_iboun(iboun) % a(igaub) = prthe(1)*gbwme/(gasco*gbtem)
                            end do
                         end if
                      end if
                   end do

                else if( prope_ker % wlaws(imate) == 'LMGPL' ) then

                   !----------------------------------------------------------------------
                   !
                   !  Low Mach approximation, density = p(th) W / R T at gauss_points where W is mean molecular density
                   !                         d(density)/dx = p(th) dW/dx / (R T) - p(th) W R/ (R T)^2 dT/dx
                   !  The properties are looked up on the Gauss points 
                   !
                   !----------------------------------------------------------------------

                   do ielem = 1,nelem
                      if( lmate(ielem) == imate ) then
                         pelty = ltype(ielem)
                         if( pelty > 0 ) then
                            pgaus = ngaus(pelty)
                            prope_ker % value_ielem(ielem) % a(1:pgaus) =&
                                 prthe(which_time)*wmean_gp(ielem) % a(1:pgaus,which_time,1)&
                                 /(gasco*tempe_gp(ielem) % a(1:pgaus,which_time,1)) 
                         end if
                      end if
                   end do

                   prope_ker % kfl_nedsm = 1



                   nullify( auxde )
                   call memory_alloca(mem_modul(1:2,modul),'AUXDE','ker_proper',auxde,npoin)
                   call ker_proper('DENSI','NPOIN',dummi,dummi,auxde)

                   do iboun = 1,nboun
                      pblty = abs(ltypb(iboun))
                      pnodb = nnode(pblty)
                      ielem = lelbo(iboun)
                      if( lmate(ielem) == imate ) then
                         pelty = ltype(ielem)
                         if( pelty > 0 ) then
                            pgaub = ngaus(pblty)
                            do inodb = 1,pnodb
                               ipoin =  lnodb(inodb,iboun)
                               elden(inodb) = auxde(ipoin)
                            end do
                            do igaub = 1,pgaub
                               gbden(1) = 0.0_rp
                               do inodb = 1,pnodb
                                  gbden(1) = gbden(1) + elden(inodb) * elmar(pblty) % shape(inodb,igaub)
                               end do
                               prope_ker % value_iboun(iboun) % a(igaub) =  gbden(1)
                            end do
                         end if
                      end if
                   end do

                   call memory_deallo(mem_modul(1:2,modul),'AUXDE','ker_proper',auxde)

                else if( prope_ker % wlaws(imate) == 'MIXTU' ) then

                   !----------------------------------------------------------------------
                   !
                   ! Mixture of constant density fluids, 1/density = Sum_k  Y_k / rho_k where Y_k is the mass fraction
                   !
                   !----------------------------------------------------------------------

                   do kpoin = 1,nmatn(imate)
                      ipoin = lmatn(imate) % l(kpoin)
                      prope_ker % value_ipoin(ipoin) = 0.0_rp
                      do ispec = 1,nspec
                         prope_ker % value_ipoin(ipoin) = prope_ker % value_ipoin(ipoin) &
                              + conce(ipoin,ispec,which_time)/speci(ispec)%densi(1)
                      enddo
                      if ( prope_ker % value_ipoin(ipoin) .gt. 0.0_rp) then
                         prope_ker % value_ipoin(ipoin) = 1.0_rp /  prope_ker % value_ipoin(ipoin)
                      else
                         call runend("KERMOD: Density evaluated to <= 0, check species properties")
                      endif
                   end do

                else if( prope_ker % wlaws(imate) == 'BIPHA' ) then

                   !----------------------------------------------------------------------
                   !
                   ! Mixture of constant density fluids in two phases
                   ! 1/density = Sum_k  Y_k / rho_k where Y_k is the mass fraction
                   ! If level set is negative, density corresponds to first phase, positive to second phase
                   !
                   !----------------------------------------------------------------------

                   if( thicl > 0.0_rp ) then
                      eps = thicl
                   else
                      call runend('Really need hleng in Ker_updpro.')
                      !eps = -thicl * hleng(1)
                   end if
                   ooeps = 1.0_rp/eps
                   oovpi = 1.0_rp/pi

                   do kpoin = 1,nmatn(imate)
                      ipoin = lmatn(imate) % l(kpoin)
                      !
                      ! Compute density of phases
                      !
                      densi_phase(1) = 0.0_rp
                      densi_phase(2) = 0.0_rp
                      do ispec = 1,nspec
                         densi_phase(1) = densi_phase(1) + conce(ipoin,ispec,1)/speci(ispec)%densi(1)
                         densi_phase(2) = densi_phase(2) + conce(ipoin,ispec,1)/speci(ispec)%densi(2)
                      enddo
                      if ( densi_phase(1) .gt. 0.0_rp) then
                         densi_phase(1) = 1.0_rp /  densi_phase(1)
                      else
                         call runend("KERMOD: Density in phase 1 evaluated to <= 0, check species properties")
                      endif
                      if ( densi_phase(2) .gt. 0.0_rp) then
                         densi_phase(2) = 1.0_rp /  densi_phase(2)
                      else
                         call runend("KERMOD: Density in phase 2 evaluated to <= 0, check species properties")
                      endif
                      ! Compute density at point
                      phi=fleve(ipoin,which_time)
                      if( phi < -eps ) then
                         ! First phase
                         prope_ker % value_ipoin(ipoin)=densi_phase(1)
                      else if( phi > eps ) then
                         ! Second phase
                         prope_ker % value_ipoin(ipoin)=densi_phase(2)
                      else
                         f = 0.5_rp*(1.0_rp+phi*ooeps+sin(pi*phi*ooeps)*oovpi)
                         prope_ker % value_ipoin(ipoin)=densi_phase(1) + (densi_phase(2)-densi_phase(1))*f
                      endif

                   end do

                else if( prope_ker % wlaws(imate) == 'TBIPH' ) then

                   !----------------------------------------------------------------------
                   !
                   ! Mixture of temperature dependent density fluids in two phases
                   ! 1/density = Sum_k  Y_k / rho_k where Y_k is the mass fraction
                   ! If level set is negative, density corresponds to first phase, positive to second phase
                   ! Temperature dependence is assumed linear, rho = rho_a * T + rho_0
                   !
                   !----------------------------------------------------------------------

                   if( thicl > 0.0_rp ) then
                      eps = thicl
                   else
                      call runend('Really need hleng in Ker_updpro.')
                      !eps = -thicl * hleng(1)
                   end if
                   ooeps = 1.0_rp/eps
                   oovpi = 1.0_rp/pi

                   do kpoin = 1,nmatn(imate)
                      ipoin = lmatn(imate) % l(kpoin)
                      !
                      ! Compute density of phases
                      !
                      densi_phase(1) = 0.0_rp
                      densi_phase(2) = 0.0_rp
                      do ispec = 1,nspec
                         densi_phase(1) = densi_phase(1) + conce(ipoin,ispec,which_time)/ &
                              ( speci(ispec)%densi(1)*tempe(ipoin,which_time)+speci(ispec)%densi(2) )
                         densi_phase(2) = densi_phase(2) + conce(ipoin,ispec,which_time)/ &
                              ( speci(ispec)%densi(3)*tempe(ipoin,which_time)+speci(ispec)%densi(4) )
                      enddo
                      if ( densi_phase(1) .gt. 0.0_rp) then
                         densi_phase(1) = 1.0_rp /  densi_phase(1)
                      else
                         call runend("KERMOD: Density in phase 1 evaluated to <= 0, check species properties")
                      endif
                      if ( densi_phase(2) .gt. 0.0_rp) then
                         densi_phase(2) = 1.0_rp /  densi_phase(2)
                      else
                         call runend("KERMOD: Density in phase 2 evaluated to <= 0, check species properties")
                      endif
                      ! Compute density at point
                      phi=fleve(ipoin,which_time)
                      if( phi < -eps ) then
                         ! First phase
                         prope_ker % value_ipoin(ipoin)=densi_phase(1)
                      else if( phi > eps ) then
                         ! Second phase
                         prope_ker % value_ipoin(ipoin)=densi_phase(2)
                      else
                         f = 0.5_rp*(1.0_rp+phi*ooeps+sin(pi*phi*ooeps)*oovpi)
                         prope_ker % value_ipoin(ipoin)=densi_phase(1) + (densi_phase(2)-densi_phase(1))*f
                      endif

                   end do

                else if( prope_ker % wlaws(imate) == 'DNAST' ) then

                   !----------------------------------------------------------------------
                   !
                   ! Density is controlled by NASTAL module
                   !
                   !----------------------------------------------------------------------

                   do kpoin = 1,nmatn(imate)
                      ipoin = lmatn(imate) % l(kpoin)
                      prope_ker % value_ipoin(ipoin) = densi(ipoin,which_time)
                   end do


                else if( prope_ker % wlaws(imate) == 'DENGP' ) then

                   !----------------------------------------------------------------------
                   !
                   !  Density read from matrix on Gauss points
                   !
                   !----------------------------------------------------------------------

                   do ielem = 1,nelem
                      if( lmate(ielem) == imate ) then
                         pelty = ltype(ielem)
                         if( pelty > 0 ) then
                            pgaus = ngaus(pelty)
                            prope_ker % value_ielem(ielem) % a(1:pgaus) = densi_gp(ielem) % a(1:pgaus,1,1)
                         end if
                      end if
                   end do

                   prope_ker % kfl_nedsm = 1

                else if( prope_ker % wlaws(imate) == 'DENNP' ) then

                   !----------------------------------------------------------------------
                   !
                   !  Density read from matrix on nodes
                   !
                   !----------------------------------------------------------------------

                   do kpoin = 1,nmatn(imate)
                      ipoin = lmatn(imate) % l(kpoin)
                      prope_ker % value_ipoin(ipoin) = densi(ipoin,1)
                   end do


                else if( prope_ker % wlaws(imate) == 'VISNP' ) then

                   !----------------------------------------------------------------------
                   !
                   !  Viscosity read from matrix on nodes
                   !
                   !----------------------------------------------------------------------

                   do kpoin = 1,nmatn(imate)
                      ipoin = lmatn(imate) % l(kpoin)
                      prope_ker % value_ipoin(ipoin) = visco(ipoin,1)
                   end do

                else if( prope_ker % wlaws(imate) == 'MUTAB' .or. prope_ker % wlaws(imate) == 'TMUTA') then

                   !----------------------------------------------------------------------
                   !
                   ! Viscosity is read from table in CHEMIC module (Flamelet combustion model)
                   !
                   !----------------------------------------------------------------------

                   do kpoin = 1,nmatn(imate)
                      ipoin = lmatn(imate) % l(kpoin)
                      prope_ker % value_ipoin(ipoin) = visck(ipoin,1)
                   end do

                else if( prope_ker % wlaws(imate) == 'MUGPT' ) then

                   !----------------------------------------------------------------------
                   !
                   !  Viscosity from CHEMIC, read from table on Gauss points
                   !
                   !----------------------------------------------------------------------

                   do ielem = 1,nelem
                      if( lmate(ielem) == imate ) then
                         pelty = ltype(ielem)
                         if( pelty > 0 ) then
                            pgaus = ngaus(pelty)
                            prope_ker % value_ielem(ielem) % a(1:pgaus) = visco_gp(ielem) % a(1:pgaus,1,1)
                         end if
                      end if
                   end do

                   prope_ker % kfl_nedsm = 1

                else if( prope_ker % wlaws(imate) == 'KGPTA' ) then

                   !----------------------------------------------------------------------
                   !
                   !  Thermal conductivity from CHEMIC, read from table on Gauss points
                   !
                   !----------------------------------------------------------------------

                   do ielem = 1,nelem
                      if( lmate(ielem) == imate ) then
                         pelty = ltype(ielem)
                         if( pelty > 0 ) then
                            pgaus = ngaus(pelty)
                            prope_ker % value_ielem(ielem) % a(1:pgaus) = condu_gp(ielem) % a(1:pgaus,1,1)
                         end if
                      end if
                   end do

                   prope_ker % kfl_nedsm = 1


                else if( prope_ker % wlaws(imate) == 'CPGPO' ) then

                   !----------------------------------------------------------------------
                   !
                   !  Specific heat at Gauss points from Chemic
                   !
                   !----------------------------------------------------------------------

                   do ielem = 1,nelem
                      if( lmate(ielem) == imate ) then
                         pelty = ltype(ielem)
                         if( pelty > 0 ) then
                            pgaus = ngaus(pelty)
                            prope_ker % value_ielem(ielem) % a(1:pgaus) = sphec_gp(ielem) % a(1:pgaus,1,1)
                         end if
                      end if
                   end do

                   prope_ker % kfl_nedsm = 1


                else if( prope_ker % wlaws(imate) == 'CPGPT' ) then

                   !----------------------------------------------------------------------
                   !
                   !  Specific heat from CHEMIC, read from table on Gauss points
                   !
                   !----------------------------------------------------------------------

                   do ielem = 1,nelem
                      if( lmate(ielem) == imate ) then
                         pelty = ltype(ielem)
                         if( pelty > 0 ) then
                            pgaus = ngaus(pelty)
                            prope_ker % value_ielem(ielem) % a(1:pgaus) = sphec_gp(ielem) % a(1:pgaus,1,1)
                         end if
                      end if
                   end do

                   prope_ker % kfl_nedsm = 1

                   !do iboun = 1,nboun
                   !   pblty = abs(ltypb(iboun))
                   !   pnodb = nnode(pblty)
                   !   ielem = lboel(pnodb+1,iboun)
                   !   if( lmate(ielem) == imate ) then
                   !      pelty = ltype(ielem)
                   !      if( pelty > 0 ) then
                   !         pgaub = ngaus(pblty)
                   !         do igaub = 1,pgaub
                   !            prope_ker % value_iboun(iboun) % a(igaub) =  sphec_gp(ielem) % a(igaub)
                   !         end do
                   !      end if
                   !   end if
                   !end do

                else if( prope_ker % wlaws(imate) == 'KTABL' .or. prope_ker % wlaws(imate) == 'TKTAB') then

                   !----------------------------------------------------------------------
                   !
                   ! Heat conductivity is read from table in CHEMIC module (CFI combustion model)
                   !
                   !----------------------------------------------------------------------

                   do kpoin = 1,nmatn(imate)
                      ipoin = lmatn(imate) % l(kpoin)
                      prope_ker % value_ipoin(ipoin) = condk(ipoin,1)
                   end do

                else if( prope_ker % wlaws(imate) == 'KQUAR') then

                   !----------------------------------------------------------------------
                   !
                   ! Heat conductivity for quartz glass is computed from polynomial function
                   ! based on temperature
                   ! Source: O. Sergeev, A. Shashkov, A. Umanskii, Thermophysical properties of
                   !         quartz glass, Journal of engineering physics 43 (1982) 1375–1383
                   !
                   !----------------------------------------------------------------------

                   do kpoin = 1,nmatn(imate)
                      ipoin = lmatn(imate) % l(kpoin)
                      T     = tempe(ipoin,which_time)
                      prope_ker % value_ipoin(ipoin) = 2.14749_rp - 298.76_rp * T**(-1) &
                           + 20720.0_rp * T**(-2) - 540000.0_rp * T**(-3)
                   end do

                else if( prope_ker % wlaws(imate) == 'TKNPO') then

                   !----------------------------------------------------------------------
                   !
                   ! Thermal conductivity is taken at nodal points
                   !
                   !----------------------------------------------------------------------

                   do kpoin = 1,nmatn(imate)
                      ipoin = lmatn(imate) % l(kpoin)
                      prope_ker % value_ipoin(ipoin) = condk(ipoin,1)
                   end do

                else if( prope_ker % wlaws(imate) == 'CPQUA') then

                   !----------------------------------------------------------------------
                   !
                   ! Specific heat for quartz glass is computed from polynomial function
                   ! based on temperature
                   ! Source: O. Sergeev, A. Shashkov, A. Umanskii, Thermophysical properties of
                   !         quartz glass, Journal of engineering physics 43 (1982) 1375–1383
                   !
                   !----------------------------------------------------------------------

                   do kpoin = 1,nmatn(imate)
                      ipoin = lmatn(imate) % l(kpoin)
                      T     = tempe(ipoin,which_time)
                      prope_ker % value_ipoin(ipoin) = 931.3_rp + 0.256_rp * T - 24000000.0_rp * T**(-2)
                   end do

                else if( prope_ker % wlaws(imate) == 'KSTAI') then

                   !----------------------------------------------------------------------
                   !
                   ! Heat conductivity for stainless steel is computed from polynomial
                   ! function based on temperature
                   ! Source: EN1993 1.2 Annex C
                   !
                   !----------------------------------------------------------------------

                   do kpoin = 1,nmatn(imate)
                      ipoin = lmatn(imate) % l(kpoin)
                      T     = tempe(ipoin,which_time) - 273.15_rp
                      prope_ker % value_ipoin(ipoin) = 14.6_rp + 0.0127_rp * T
                   end do

                else if( prope_ker % wlaws(imate) == 'CPSTA') then

                   !----------------------------------------------------------------------
                   !
                   ! Specific heat for stainless steel is computed from polynomial
                   ! function based on temperature
                   ! Source: EN1993 1.2 Annex C
                   !
                   !----------------------------------------------------------------------

                   do kpoin = 1,nmatn(imate)
                      ipoin = lmatn(imate) % l(kpoin)
                      T     = tempe(ipoin,which_time) - 273.15_rp
                      prope_ker % value_ipoin(ipoin) = 450.0_rp + 0.28_rp * T &
                           - 0.000291_rp * T**2_rp + 0.000000134_rp * T**3_rp
                   end do

                else if( prope_ker % wlaws(imate) == 'CPTAB' .or. prope_ker % wlaws(imate) == 'TCPTA') then

                   !----------------------------------------------------------------------
                   !
                   ! Specific heat is read from table in CHEMIC module (CFI combustion model)
                   !
                   !----------------------------------------------------------------------

                   do kpoin = 1,nmatn(imate)
                      ipoin = lmatn(imate) % l(kpoin)
                      prope_ker % value_ipoin(ipoin) = sphek(ipoin,1)
                   end do

                else if( prope_ker % wlaws(imate) == 'CPNPO') then

                   !----------------------------------------------------------------------
                   !
                   ! Specific heat is taken from nodal points
                   !
                   !----------------------------------------------------------------------

                   do kpoin = 1,nmatn(imate)
                      ipoin = lmatn(imate) % l(kpoin)
                      prope_ker % value_ipoin(ipoin) = sphek(ipoin,1)
                   end do

                else if( prope_ker % wlaws(imate) == 'MUMIX' ) then

                   !----------------------------------------------------------------------
                   !
                   ! Viscosity given by a mixture of species (individual viscosities are computed
                   ! inside CHEMIC, here only the mixture is done)
                   !
                   !----------------------------------------------------------------------

                   do kpoin = 1,nmatn(imate)
                      ipoin = lmatn(imate) % l(kpoin)
                      !
                      ! Mixture coefficients see Combustion Physics, Chung K. Law, page 155
                      !
                      ! phi(i,j) = 1/sqrt(8) 1/sqrt(1 + Wi/Wj) [1+ sqtr(mui/muj) (Wj/Wi)^1/4 ]^2
                      ! phi(i,i) = 1.0
                      !
                      ! Took out wmean(ipoin,1) from all terms because it simplifies in the final expression
                      !
                      do ispec = 1,nspec
                         sumphi(ispec) = 0.0_rp
                         do jspec = 1,nspec
                            if (jspec /= ispec) then
                               if (visck(ipoin,jspec) /= 0.0_rp) then
                                  phiij = sqrt(0.125_rp/(1.0_rp+speci(ispec)%weigh/speci(jspec)%weigh)) * &
                                       (1.0_rp+sqrt( visck(ipoin,ispec)&
                                       / visck(ipoin,jspec))*(speci(jspec)%weigh/speci(ispec)%weigh)**0.25_rp)**2_rp
                               else
                                  phiij =  sqrt(0.125_rp/(1.0_rp+speci(ispec)%weigh/speci(jspec)%weigh))
                               endif
                               sumphi(ispec) = sumphi(ispec) + phiij * conce(ipoin,jspec,which_time) / speci(jspec)%weigh
                            endif
                         enddo
                      enddo
                      !
                      prope_ker % value_ipoin(ipoin) = 0.0_rp
                      do ispec = 1,nspec
                         if (conce(ipoin,ispec,which_time) /= 0.0_rp) then
                            prope_ker % value_ipoin(ipoin) = prope_ker % value_ipoin(ipoin) + visck(ipoin,ispec) * &
                                 (conce(ipoin,ispec,which_time) / speci(ispec)%weigh) &
                                 / (conce(ipoin,ispec,which_time) / speci(ispec)%weigh + sumphi(ispec))
                         endif

                      enddo
                   end do

                else if( prope_ker % wlaws(imate) == 'KMIXT' ) then

                   !----------------------------------------------------------------------
                   !
                   ! Heat conductivity given by a mixture of species (individual conductivities are computed
                   ! inside CHEMIC, here only the mixture is done)
                   !
                   !----------------------------------------------------------------------

                   do kpoin = 1,nmatn(imate)
                      ipoin = lmatn(imate) % l(kpoin)
                      !
                      ! Mixture coefficients see Combustion Physics, Chung K. Law, page 155
                      ! Notice that it's ok that viscosities are used...
                      !
                      ! phi(i,j) = 1/sqrt(8) 1/sqrt(1 + Wi/Wj) [1+ sqtr(mui/muj) (Wj/Wi)^1/4 ]^2
                      ! phi(i,i) = 1.0
                      !
                      ! Took out wmean(ipoin,1) from all terms because it simplifies in the final expression
                      !
                      do ispec = 1,nspec
                         sumphi(ispec) = 0.0_rp
                         do jspec = 1,nspec
                            if (jspec /= ispec) then
                               if (visck(ipoin,jspec) /= 0.0_rp) then
                                  phiij = sqrt(0.125_rp/(1.0_rp+speci(ispec)%weigh/speci(jspec)%weigh)) * &
                                       (1.0_rp+sqrt( visck(ipoin,ispec)&
                                       / visck(ipoin,jspec))*(speci(jspec)%weigh/speci(ispec)%weigh)**0.25_rp)**2
                               else
                                  phiij =  sqrt(0.125_rp/(1.0_rp+speci(ispec)%weigh/speci(jspec)%weigh))
                               endif
                               sumphi(ispec) = sumphi(ispec) + phiij * conce(ipoin,jspec,which_time) / speci(jspec)%weigh
                            endif
                         enddo
                      enddo
                      !
                      prope_ker % value_ipoin(ipoin) = 0.0_rp
                      do ispec = 1,nspec
                         if (conce(ipoin,ispec,which_time) /= 0.0_rp) then
                            prope_ker % value_ipoin(ipoin) = prope_ker % value_ipoin(ipoin) + condk(ipoin,ispec) * &
                                 (conce(ipoin,ispec,which_time) / speci(ispec)%weigh) &
                                 / (conce(ipoin,ispec,which_time) / speci(ispec)%weigh + 1.065_rp * sumphi(ispec))
                         endif
                      enddo
                   end do

                else if( prope_ker % wlaws(imate) == 'CPMIX' ) then

                   !----------------------------------------------------------------------
                   !
                   ! Specific heat of a mixture of species (individual speheas are computed
                   ! inside CHEMIC, here only the mixture is done)
                   !
                   !----------------------------------------------------------------------

                   do kpoin = 1,nmatn(imate)
                      ipoin = lmatn(imate) % l(kpoin)
                      !
                      prope_ker % value_ipoin(ipoin) = 0.0_rp
                      do ispec = 1,nspec
                         prope_ker % value_ipoin(ipoin) = prope_ker % value_ipoin(ipoin) &
                              + sphek(ipoin,ispec) * conce(ipoin,ispec,which_time)
                      enddo
                   end do

                else if( prope_ker % wlaws(imate) == 'TEST1' ) then

                   !----------------------------------------------------------------------
                   !
                   ! TEST1
                   !
                   !----------------------------------------------------------------------

                   do ielem = 1,nelem
                      if( lmate(ielem) == imate ) then
                         pelty = ltype(ielem)
                         if( pelty > 0 ) then
                            pgaus = ngaus(pelty)
                            pnode = lnnod(ielem)
                            gpcor = 0.0_rp
                            do igaus = 1,pgaus
                               do inode = 1,pnode
                                  ipoin = lnods(inode,ielem)
                                  do idime = 1,ndime
                                     gpcor(idime,igaus) = gpcor(idime,igaus) + coord(idime,ipoin) * elmar(pelty) % shape(inode,igaus)
                                  end do
                               end do
                               prope_ker % value_ielem(ielem) % a(igaus)   =  2.0_rp * gpcor(1,igaus) - 3.0_rp * gpcor(2,igaus)
                               prope_ker % grval_ielem(ielem) % a(1,igaus) =  2.0_rp
                               prope_ker % grval_ielem(ielem) % a(2,igaus) = -3.0_rp
                            end do
                         end if
                      end if
                   end do
                   prope_ker % kfl_nedsm = 1

                else if( prope_ker % wlaws(imate) == 'TEST2' ) then

                   !----------------------------------------------------------------------
                   !
                   ! TEST2
                   !
                   !----------------------------------------------------------------------

                   do kpoin = 1,nmatn(imate)
                      ipoin = lmatn(imate) % l(kpoin)
                      prope_ker % value_ipoin(ipoin) = 2.0_rp * coord(1,ipoin) - 3.0_rp * coord(2,ipoin)
                   end do

                else if( prope_ker % wlaws(imate) == 'TEST3' ) then

                   !----------------------------------------------------------------------
                   !
                   ! TEST3
                   !
                   !----------------------------------------------------------------------

                   do kpoin = 1,nmatn(imate)
                      ipoin = lmatn(imate) % l(kpoin)
                      prope_ker % value_ipoin(ipoin)   =  2.0_rp * coord(1,ipoin) - 3.0_rp * coord(2,ipoin)
                      prope_ker % grval_ipoin(1,ipoin) =  2.0_rp
                      prope_ker % grval_ipoin(2,ipoin) = -3.0_rp
                   end do

                else if( prope_ker % wlaws(imate) == 'TEST4' ) then

                   !----------------------------------------------------------------------
                   !
                   ! TEST4
                   !
                   !----------------------------------------------------------------------

                   if( ndime == 2 ) then
                      do kpoin = 1,nmatn(imate)
                         ipoin = lmatn(imate) % l(kpoin)
                         prope_ker % value_ipoin(ipoin)   = 1.0_rp + 2.0_rp * coord(1,ipoin) + 3.0_rp * coord(2,ipoin)
                         prope_ker % grval_ipoin(1,ipoin) = 2.0_rp
                         prope_ker % grval_ipoin(2,ipoin) = 3.0_rp
                      end do
                   else
                      do kpoin = 1,nmatn(imate)
                         ipoin = lmatn(imate) % l(kpoin)
                         prope_ker % value_ipoin(ipoin)   = &
                              1.0_rp + 2.0_rp * coord(1,ipoin) + 3.0_rp * coord(2,ipoin) + 4.0_rp * coord(3,ipoin)
                         prope_ker % grval_ipoin(1,ipoin) = 2.0_rp
                         prope_ker % grval_ipoin(2,ipoin) = 3.0_rp
                         prope_ker % grval_ipoin(3,ipoin) = 4.0_rp
                      end do
                   end if

                else if( prope_ker % wlaws(imate) == 'ABL  ' ) then

                   !----------------------------------------------------------------------
                   !
                   ! ABL: mu = rho/Cmu*kap*(z+z0)*U_*
                   !
                   !----------------------------------------------------------------------

                   prope_ker % kfl_nedsm = 1
                   do ielem = 1,nelem
                      if( lmate(ielem) == imate ) then
                         pelty = abs(ltype(ielem))
                         pgaus = ngaus(pelty)
                         pnode = lnnod(ielem)
                         gpcor = 0.0_rp
                         do igaus = 1,pgaus
                            do inode = 1,pnode
                               ipoin = lnods(inode,ielem)
                               do idime = 1,ndime
                                  gpcor(idime,igaus) = gpcor(idime,igaus) &
                                       + coord(idime,ipoin) * elmar(pelty) % shape(inode,igaus)
                               end do
                            end do
                            if( kfl_rough > 0 ) then ! variable roughness
                               z0 = 0.0_rp
                               do inode = 1,pnode
                                  ipoin = lnods(inode,ielem)
                                  z0    = z0 + rough(ipoin) * elmar(pelty) % shape(inode,igaus)
                               end do
                            else ! constant roughness
                               z0 = rough_dom
                            end if
                            rho   = 1.0_rp
                            ustar = 0.23_rp
                            kap   = 0.41_rp
                            prope_ker % value_ielem(ielem) % a(      igaus) =  rho*ustar*kap*(gpcor(ndime,igaus)+z0+delta_dom)
                            prope_ker % grval_ielem(ielem) % a(    1,igaus) =  0.0_rp
                            prope_ker % grval_ielem(ielem) % a(    2,igaus) =  0.0_rp
                            prope_ker % grval_ielem(ielem) % a(ndime,igaus) =  rho*ustar*kap
                         end do
                      end if
                   end do

                else if( prope_ker % wlaws(imate) == ker_polynomial_name ) then
                   call ker_polynomial_proper(prope_ker, imate, which_time)

                end if

             end do

          end if

       end if

    end do

    deallocate(listp_ker)

    !----------------------------------------------------------------------
    !
    ! Assign previous values for adjoint case to temper and conce
    !
    !----------------------------------------------------------------------

    if( kfl_adj_prob == 1_ip ) then
       conce => conce_save
       tempe => tempe_save
       untur => untur_save
       advec => advec_save
       nullify(conce_save)
       nullify(tempe_save)
       nullify(untur_save)
       nullify(advec_save)
    end if

    call cputim(time2)
    cpu_prope(UPDATE_PROPERTIES) = cpu_prope(UPDATE_PROPERTIES) + time2-time1

  end subroutine ker_updpro

end module mod_ker_updpro
!> @}

  
