!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    mod_sld_vect_stress_model_154.f90
!> @author  Adria Quintanas (adria.quintanas@udg.edu)
!> @author  Gerard Guillamet (gerard.guillamet@bsc.es)
!> @date    June, 2022
!> @brief   Damage model for transversally orthotropic materials.
!>          Vectorized version 
!>                   
!> @details 
!>
!>          State Dependent Variables:\n
!>
!>          SDV01 := rFTnew (Damage threshold variables)
!>          SDV02 := rFCnew
!>          SDV03 := rMTnew
!>          SDV04 := rMCnew
!>          SDV05 := d1 (Damage index fiber)  
!>          SDV06 := d2 (Damage index matrix)
!>          SDV07 := d3 (Damage index coupling matrix/fiber compression)
!>          SDV08 := d4 (Damage index matrix shear)
!>          SDV09 := d5 (Damage index fiber tension)
!>          SDV10 := d6 (Damage index matrix shear + coupling with fiber tension)
!>          SDV11 := grePlast (Plastic strain)
!>          SDV12 := greInela (Elastic strain)
!>
!>          SDV13 := lenchar
!>          SDV14 := xT
!>          SDV15 := xT0
!>          SDV16 := xC
!>          SDV17 := xC0
!>          SDV18 := yT
!>          SDV19 := yC
!>          SDV20 := sL
!>
!>          Strength reduction properties:\n
!>
!>          SRP01: lenchar (characteristic element length (max.))
!>          SRP02: xT
!>          SRP03: xT0
!>          SRP04: xC
!>          SRP05: xC0
!>          SRP06: yT
!>          SRP07: yC
!>          SRP08: sL
!>
!>          References:\n
!>
!>          Soto, A. and Gonzalez, E.V. Maimi, P. Mayugo, J.A. Pasquali, P.R. Camanho, P.P.\n
!>          A methodology to simulate low velocity impact and compression after impact in \n
!>          large composite stiffened panels, 2018 \n
!>
!>          Maimi,P., Camanho, P.P., Mayugo, J.A., Davila,C.G. A continuum damage model for \n
!>          composite laminates: Part II - Computational, \n
!>          implementation and validation, 2007 \n
!>
!>          You are invited to use the subroutine for academic research       
!>          purposes only. If you are going to use the subroutine for         
!>          industrial purposes, please notice to the authors.                
!>                                                                      
!>          Please cite the previous papers in your work if you are using the subroutine.
!>          If you have any comment/suggestion you are welcome to send it to the authors.
!>          
!>          To share is to improve.
!>
!> @todo    To do list:\n
!>
!>          - Improve initialization of sm154 in begrun
!>
!> @}
!-----------------------------------------------------------------------

module mod_sld_vect_stress_model_154

#include "def_solidz_vect_dim.inc" 
   use def_kintyp, only  : ip, rp, lg

   implicit none

   integer(ip), parameter :: nsdv154   = 20_ip    !< 12 SDVS + 8 PROPERTIES
   integer(ip), parameter :: nprops154 = 41_ip
   real(rp)               :: d1Max, d2Max, d3Max, d4Max, d5Max, d6Max
   real(rp),    parameter :: zeroapp = 0.001_rp
   real(rp),    parameter :: damMax  = 0.99999_rp !< UdG uses 0.999_rp. This allows to do not dissipate the whole energy of the law. 
   real(rp),    parameter :: minG_i  = 0.001_rp   !< UdG strength reduction limits
   real(rp),    parameter :: minG_j  = 0.0001_rp              

   private

   public :: nsdv154
   public :: sld_stress_model_154_initialisations
   public :: sld_stress_model_154_compute_properties
   public :: sld_stress_model_154_message
   public :: sld_vect_stress_model_154
   public :: sld_vect_stress_model_154_getSDV
   public :: sld_vect_stress_model_154_setSDV

 contains
   
   subroutine sld_stress_model_154_initialisations(&
        ireducedXT,ireducedXTO,ireducedXC,ireducedXCO,ireducedYT,ireducedYC,ireducedSL)

     integer(ip), intent(out) :: ireducedXT
     integer(ip), intent(out) :: ireducedXTO
     integer(ip), intent(out) :: ireducedXC
     integer(ip), intent(out) :: ireducedXCO
     integer(ip), intent(out) :: ireducedYT
     integer(ip), intent(out) :: ireducedYC
     integer(ip), intent(out) :: ireducedSL
     !
     ! Initializations
     ! 
     ireducedXT  = 0_ip
     ireducedXTO = 0_ip
     ireducedXC  = 0_ip
     ireducedXCO = 0_ip
     ireducedYT  = 0_ip
     ireducedYC  = 0_ip
     ireducedSL  = 0_ip
     
   end subroutine sld_stress_model_154_initialisations

   subroutine sld_stress_model_154_compute_properties(ielem,pmate,lench, &
        ireducedXT, ireducedXTO, ireducedXC, ireducedXCO, ireducedYT, ireducedYC, ireducedSL)

     use def_solidz, only : parco_sld, srpro_sld
     
     implicit none
     
     integer(ip), intent(in)    :: ielem
     integer(ip), intent(in)    :: pmate
     real(rp),    intent(in)    :: lench
     integer(ip), intent(inout) :: ireducedXT
     integer(ip), intent(inout) :: ireducedXTO
     integer(ip), intent(inout) :: ireducedXC
     integer(ip), intent(inout) :: ireducedXCO
     integer(ip), intent(inout) :: ireducedYT
     integer(ip), intent(inout) :: ireducedYC
     integer(ip), intent(inout) :: ireducedSL

     real(rp)    :: E11, E22, G12, G23, v12, v21, v23
     real(rp)    :: xT, xTO, xC, xCO, yT, yC, sL, alpha0
     real(rp)    :: xTnew, xTOnew, xCnew, xCOnew, yTnew, yCnew, sLnew
     real(rp)    :: gxT, gxTO, gxC, gXCO, gYT, gYC, gSL, gG
     real(rp)    :: sP, kP
     real(rp)    :: eta
     real(rp)    :: a11, a22, dTemp, b11, b22, dMoist
     real(rp)    :: thick
     real(rp)    :: lench_aux
     real(rp)    :: mugsL, mufxc                                     !< Friction coefficients
     real(rp)    :: yBT, yBC                                         !< Biaxial strenghts
     real(rp)    :: E1c                                              !< Fiber Young Modulus compression
     
     
     call GetProperties( parco_sld(1:nprops154,pmate), E11, E22, v12, v21, v23, G12, G23, &
          xT, xTO, xC, xCO, yT, yC, sL, &
          alpha0, gxT, gxTO, gxC, gxCO, gyT, gyC, gsL, gG, &
          sP, kP, eta, a11, a22, dTemp, b11, b22, dMoist, &
          thick, mugsL, yBT, yBC, mufxc, E1c)
     !
     ! Modify other properties
     !
     lench_aux = lench
     if( thick > 0.0_rp ) then
        lench_aux = sqrt(lench**3/thick)
     end if
     !
     ! Modify strengths (reduce them to ensure correct energy dissipation)
     !
     xTnew  = xT
     xTOnew = xTO
     xCnew  = xC
     xCOnew = xCO
     yTnew  = yT
     yCnew  = yC
     sLnew  = sL
     call reduceStrengths( E11, E22, G12, gxT, gxTO, gxC, gXCO, gyT, gyC, gsL, &
          lench_aux, xTnew, xTOnew, xCnew, xCOnew, yTnew, yCnew, sLnew )

     if( abs(xT-xTnew)   > epsilon(1.0_rp)) ireducedXT = 1_ip
     if( abs(xTO-xTOnew) > epsilon(1.0_rp)) ireducedXTO = 1_ip
     if( abs(xC-xCnew)   > epsilon(1.0_rp)) ireducedXC = 1_ip
     if( abs(xCO-xCOnew) > epsilon(1.0_rp)) ireducedXCO = 1_ip
     if( abs(yT-yTnew)   > epsilon(1.0_rp)) ireducedYT = 1_ip
     if( abs(yC-yCnew)   > epsilon(1.0_rp)) ireducedYC = 1_ip
     if( abs(sL-sLnew)   > epsilon(1.0_rp)) ireducedSL = 1_ip
     !
     ! Save characteristic length and material properties
     !
     srpro_sld(1,ielem) = lench_aux
     srpro_sld(2,ielem) = xTnew
     srpro_sld(3,ielem) = xTOnew
     srpro_sld(4,ielem) = xCnew
     srpro_sld(5,ielem) = xCOnew
     srpro_sld(6,ielem) = yTnew
     srpro_sld(7,ielem) = yCnew
     srpro_sld(8,ielem) = sLnew

   end subroutine sld_stress_model_154_compute_properties
   
   subroutine sld_stress_model_154_message(&
        ireducedXT,ireducedXTO,ireducedXC,ireducedXCO,ireducedYT,ireducedYC,ireducedSL)

     use def_master,         only : INOTEMPTY, INOTSLAVE
     use def_master,         only : retost
     use mod_communications, only : PAR_MIN
     use mod_communications, only : PAR_MAX
     use mod_communications, only : PAR_SUM
     use mod_messages,       only : messages_live
     use def_solidz,         only : srpro_sld

     integer(ip), intent(inout) :: ireducedXT
     integer(ip), intent(inout) :: ireducedXTO
     integer(ip), intent(inout) :: ireducedXC
     integer(ip), intent(inout) :: ireducedXCO
     integer(ip), intent(inout) :: ireducedYT
     integer(ip), intent(inout) :: ireducedYC
     integer(ip), intent(inout) :: ireducedSL
     logical(lg)                :: reducedXT, reducedXTO, reducedXC, reducedXCO, reducedYT, reducedYC, reducedSL
     real(rp)                   :: xTmin, xTOmin, xCmin, xCOmin, yTmin, yCmin, sLmin
     !
     ! Initializations
     !
     reducedXT  = .false.
     reducedXTO = .false.
     reducedXC  = .false.
     reducedXCO = .false.
     reducedYT  = .false.
     reducedYC  = .false.
     reducedSL  = .false.
     xTmin      = huge(1.0_rp)
     xTOmin     = huge(1.0_rp)
     xCmin      = huge(1.0_rp)
     xCOmin     = huge(1.0_rp)
     ytmin      = huge(1.0_rp)
     yCmin      = huge(1.0_rp)
     sLmin      = huge(1.0_rp)
     !
     ! Check if strength reduction is performed in any of the material properties
     !
     call PAR_SUM(ireducedXT, 'IN MY CODE')
     call PAR_SUM(ireducedXTO,'IN MY CODE')
     call PAR_SUM(ireducedXC, 'IN MY CODE')
     call PAR_SUM(ireducedXCO,'IN MY CODE')
     call PAR_SUM(ireducedYT, 'IN MY CODE')
     call PAR_SUM(ireducedYC, 'IN MY CODE')
     call PAR_SUM(ireducedSL, 'IN MY CODE')
     if( INOTSLAVE ) then
        if( ireducedXT  > 0_ip ) reducedXT  = .true.
        if( ireducedXTO > 0_ip ) reducedXTO = .true.
        if( ireducedXC  > 0_ip ) reducedXC  = .true.
        if( ireducedXCO > 0_ip ) reducedXCO = .true.
        if( ireducedYT  > 0_ip ) reducedYT  = .true.
        if( ireducedYC  > 0_ip ) reducedYC  = .true.
        if( ireducedSL  > 0_ip ) reducedSL  = .true.
     end if
     !
     ! Get min. values from whole mesh
     !
     if( INOTEMPTY ) then
        xTmin  = minval(srpro_sld(2,:),MASK=srpro_sld(2,:) > epsilon(1.0_rp))
        xTOmin = minval(srpro_sld(3,:),MASK=srpro_sld(3,:) > epsilon(1.0_rp))
        xCmin  = minval(srpro_sld(4,:),MASK=srpro_sld(4,:) > epsilon(1.0_rp))
        xCOmin = minval(srpro_sld(5,:),MASK=srpro_sld(5,:) > epsilon(1.0_rp))
        yTmin  = minval(srpro_sld(6,:),MASK=srpro_sld(6,:) > epsilon(1.0_rp))
        yCmin  = minval(srpro_sld(7,:),MASK=srpro_sld(7,:) > epsilon(1.0_rp))
        sLmin  = minval(srpro_sld(8,:),MASK=srpro_sld(8,:) > epsilon(1.0_rp))
     end if
     call PAR_MIN(xTmin, 'IN MY CODE')
     call PAR_MIN(xTOmin,'IN MY CODE')
     call PAR_MIN(xCmin, 'IN MY CODE')
     call PAR_MIN(xCOmin,'IN MY CODE')
     call PAR_MIN(yTmin, 'IN MY CODE')
     call PAR_MIN(yCmin, 'IN MY CODE')
     call PAR_MIN(sLmin, 'IN MY CODE')
     !
     ! Live info
     !
     if( INOTSLAVE ) then
        call messages_live('STRENGTH REDUCTION FOR SM154','START SECTION')
        if( reducedXT  ) call messages_live('XT= ' //trim(retost(xTmin, REAL_FORMAT='(F10.1)')),"WARNING")
        if( reducedXTO ) call messages_live('XTO= '//trim(retost(xTOmin,REAL_FORMAT='(F10.1)')),"WARNING")
        if( reducedXC  ) call messages_live('XC= ' //trim(retost(xCmin, REAL_FORMAT='(F10.1)')),"WARNING")
        if( reducedXCO ) call messages_live('XCO= '//trim(retost(xCOmin,REAL_FORMAT='(F10.1)')),"WARNING")
        if( reducedYT  ) call messages_live('YT= ' //trim(retost(yTmin, REAL_FORMAT='(F10.1)')),"WARNING")
        if( reducedYC  ) call messages_live('YC= ' //trim(retost(yCmin, REAL_FORMAT='(F10.1)')),"WARNING")
        if( reducedSL  ) call messages_live('SL= ' //trim(retost(sLmin, REAL_FORMAT='(F10.1)')),"WARNING")
        if( .not. reducedXT  .and. &
            .not. reducedXT  .and. &
            .not. reducedXTO .and. &
            .not. reducedXC  .and. &
            .not. reducedXCO .and. &
            .not. reducedYT  .and. &
            .not. reducedYC  .and. &
            .not. reducedSL ) then
           call messages_live('NO STRENGTHS REDUCTION REQUIRED FOR SM154')
        end if
        call messages_live('STRENGTH REDUCTION FOR SM154','END SECTION')

     end if

   end subroutine sld_stress_model_154_message
   
   subroutine sld_vect_stress_model_154(lelem_loc,pmate,pgaus,elsys,gpgre,gpstr,gpsdv)

      use def_master,         only : dtime
      use def_solidz,         only : parco_sld
      use mod_sld_vect_csm
      use mod_sld_vect_maths
      !!use mod_sld_vect_csm,   only : vcsm_compact_strain_tensor_3D
      !!use mod_sld_vect_csm,   only : vcsm_uncompact_stress_tensor_3D
      !!use mod_sld_vect_csm,   only : vcsm_compact_transformation_matrix_Tg_3D
      !!use mod_sld_vect_maths, only : vmath_MxV
      !!use mod_sld_vect_maths, only : vmath_TXV

      integer(ip), intent(in)           :: lelem_loc(DVS)
      integer(ip), intent(in)           :: pmate
      integer(ip), intent(in)           :: pgaus
      real(rp),    intent(in)           :: elsys(DVS,ZDIME,ZDIME)
      real(rp),    intent(in)           :: gpgre(DVS,ZDIME,ZDIME,ZGAUS)
      real(rp),    intent(out)          :: gpstr(DVS,ZDIME,ZDIME,ZGAUS)
      real(rp),    intent(out)          :: gpsdv(DVS,nsdv154,ZGAUS)
      real(rp)                          :: gpepsXYZ(DVS,ZVOIG)
      real(rp)                          :: gpeps123(DVS,ZVOIG)
      real(rp)                          :: gpsig123(DVS,ZVOIG)
      real(rp)                          :: gpsigXYZ(DVS,ZVOIG)
      real(rp)                          :: gpTg(DVS,ZVOIG,ZVOIG) 
      integer(ip)                       :: ig
      
      ! Get transformation matrix for compacted notation rotations
      gpTg = vcsm_compact_transformation_matrix_Tg_3D(elsys)
      ! Compute the PK2 stress tensor (S)
      do ig = 1,ZGAUS
            gpepsXYZ(:,:) = vcsm_compact_strain_tensor_3D(gpgre(:,:,:,ig))
            gpeps123(:,:) = vmath_MxV(gpTg,gpepsXYZ)
            call evalModelGP(lelem_loc, parco_sld(:,pmate), dtime, gpeps123, gpsig123, gpsdv(:,:,ig) )
            gpsigXYZ(:,:) = vmath_TXV(gpTg,gpsig123)
            gpstr(:,:,:,ig) = vcsm_uncompact_stress_tensor_3D(gpsigXYZ)
      enddo

   end subroutine sld_vect_stress_model_154


   subroutine sld_vect_stress_model_154_getSDV(lelem_loc, pgaus, gpsdv)
     
      use def_solidz, only : svegm_sld, srpro_sld
      use def_master, only : TIME_N_STATE
      
      integer(ip), intent(in)    :: lelem_loc(DVS)
      integer(ip), intent(in)    :: pgaus
      real(rp),    intent(inout) :: gpsdv(:,:,:)
      integer(ip)                :: iv, ig, ielem

      do concurrent( iv = 1:DVS )
         ielem = lelem_loc(iv)
         do ig = 1, ZGAUS
            gpsdv(iv,1:12,ig)       = svegm_sld(ielem)%a(1:12,ig,TIME_N_STATE)
            gpsdv(iv,13:nsdv154,ig) = srpro_sld(1:8,ielem)
         enddo
      enddo

   end subroutine sld_vect_stress_model_154_getSDV


   subroutine sld_vect_stress_model_154_setSDV(lelem_loc, pgaus, gpsdv)
     
      use def_solidz, only : svegm_sld
      use def_master, only : ITER_K_STATE

      integer(ip), intent(in)    :: lelem_loc(DVS)
      integer(ip), intent(in)    :: pgaus
      real(rp),    intent(inout) :: gpsdv(:,:,:)
      integer(ip)                :: iv, ig, ielem

      do concurrent( iv = 1:DVS )
         ielem = lelem_loc(iv)
         do ig = 1, ZGAUS
            svegm_sld(ielem)%a(1:nsdv154,ig,ITER_K_STATE) = gpsdv(iv,1:nsdv154,ig)
         enddo
      enddo

   end subroutine sld_vect_stress_model_154_setSDV


   subroutine evalModelGP(lelem_loc, props, dt, eps123, sig123, sdv )
     
      use def_solidz, only : kfl_insdv_sld
      use def_domain, only : xfiel

      integer(ip),  intent(in)    :: lelem_loc(DVS)
      real(rp),     intent(in)    :: props(nprops154)                         !< Props
      real(rp),     intent(in)    :: dt                                       !< Time increment
      real(rp),     intent(in)    :: eps123(DVS,6)                            !< Deformation tensor
      real(rp),     intent(out)   :: sig123(DVS,6)                            !< Stresses tensor
      real(rp),     intent(inout) :: sdv(DVS,nsdv154)                         !< Stresses tensor
      real(rp)                    :: E11, E22, G12, G23, v12, v21, v23        !< Material properties (- elastic)
      real(rp)                    :: xT, xTO, xC, xCO, yT, yC, sL, alpha0     !<   - Strength
      real(rp)                    :: nT                                       !<
      real(rp)                    :: gxT, gxTO, gxC, gXCO, gYT, gYC, gSL, gG  !<   - Fracture toughness parameters
      real(rp)                    :: sP, kP                                   !<   - Plasticity parameters
      real(rp)                    :: eta                                      !<   - Viscosity parameters
      real(rp)                    :: a11, a22, dTemp, b11, b22, dMoist        !<   - Hygrothermal parameters
      real(rp)                    :: thick                                    !< Out-of-plane thickness
      real(rp)                    :: mugsL, mufxc                             !< Friction coefficients
      real(rp)                    :: yBT, yBC                                 !< Biaxial strenghts
      real(rp)                    :: E1c                                      !< Fiber Young Modulus compression
      logical(lg)                 :: flagVisco, flagHygro                     !< 
      logical(lg), dimension(DVS) :: flagPlast, flagDamg                      !<
      logical(lg), dimension(DVS) :: flagFT, flagFC, flagMT, flagMC           !< Flags 
      real(rp),    dimension(DVS) :: eps11, eps22, eps33, eps23, eps13, eps12 !<
      real(rp),    dimension(DVS) :: sig11, sig22, sig33, sig23, sig13, sig12 !<
      real(rp),    dimension(DVS) :: epsPla, epsIne                           !< Strains
      real(rp),    dimension(DVS) :: epsH11, epsH22, epsH33                   !<
      real(rp),    dimension(DVS) :: epsT11, epsT22, epsT33                   !<
      real(rp),    dimension(DVS) :: phiFT, phiFC, phiMT, phiMC               !< Loading functions
      real(rp),    dimension(DVS) :: rFTold, rFCold, rMTold, rMCold           !< Internal variables
      real(rp),    dimension(DVS) :: rFTnew, rFCnew, rMTnew, rMCnew           !<
      real(rp),    dimension(DVS) :: d1, d2, d3, d4, d5, d6                   !< Damage variables
      real(rp),    dimension(DVS) :: vxT, vxTO, vxC, vxCO, vyT, vyC, vsL      !<
      real(rp),    dimension(DVS) :: vsT, vnL                                 !<
      real(rp),    dimension(DVS) :: velch                                    !<
      integer(ip)                 :: iv, ielem
      
      ! Material properties
      call GetProperties( &
         props, E11, E22, v12, v21, v23, G12, G23, xT, xTO, xC, xCO, yT, yC, sL, &
         alpha0, gxT, gxTO, gxC, gxCO, gyT, gyC, gsL, gG, sP, kP, &
         eta, a11, a22, dTemp, b11, b22, dMoist, &
         thick, mugsL, yBT, yBC, mufxc, E1c)
     
      ! Recover state variables
      rFTold(1:DVS) = sdv(1:DVS,1)
      rFCold(1:DVS) = sdv(1:DVS,2)
      rMTold(1:DVS) = sdv(1:DVS,3)
      rMCold(1:DVS) = sdv(1:DVS,4)
      epsPla(1:DVS) = sdv(1:DVS,11)
      epsIne(1:DVS) = sdv(1:DVS,12)
      ! Recover strengths
      velch(1:DVS) = sdv(1:DVS,13)
      vxT(1:DVS)   = sdv(1:DVS,14)
      vxTO(1:DVS)  = sdv(1:DVS,15)
      vxC(1:DVS)   = sdv(1:DVS,16)
      vxCO(1:DVS)  = sdv(1:DVS,17)
      vyT(1:DVS)   = sdv(1:DVS,18)
      vyC(1:DVS)   = sdv(1:DVS,19)
      vsL(1:DVS)   = sdv(1:DVS,20)

      ! Compute dependent material properties from reduced strengths
      vsT(1:DVS) = vyC(1:DVS)*cos(alpha0)*(sin(alpha0) + cos(alpha0)/tan(2.0_rp*alpha0))
      nT        = -1.0_rp/tan(2.0_rp*alpha0)
      vnL(1:DVS) = -vsL(1:DVS)*cos(2.0_rp*alpha0)/(vyC(1:DVS)*cos(alpha0)**2)

      if( eta > 0.0_rp ) then
         flagVisco = .true.
      else
         flagVisco = .false.
      endif
      if( dTemp > 0.0_rp .or. dMoist > 0.0_rp ) then
         flagHygro = .true.
      else
         flagHygro = .false.
      endif
      if( sP < 0.0_rp ) kP = 0.0_rp
      if( kP > 0.0_rp ) then
         flagPlast(1:DVS) = .true.
      else
         flagPlast(1:DVS) = .false.
      endif

      ! Elastic strain tensor
      ! - get plastic strain
      do concurrent( iv = 1:DVS ) 
          if( rMTold(iv) > 1.0_rp .or. rMCold(iv) > 1.0_rp )then
              flagDamg(iv) = .true.
          else
              flagDamg(iv) = .false. 
          endif
      enddo
      ! Note that hygrothermal effects doesn't influence shear strains
      call GetPlasticStrain12( & 
         G12, sP, kP, &
         eps123(1:DVS,6), epsPla(1:DVS), epsIne(1:DVS), &
         flagDamg(1:DVS), flagPlast(1:DVS) &
          )
      ! - get thermal expansion strain
      epsT11(1:DVS) = a11 * dTemp
      epsT22(1:DVS) = a22 * dTemp
      epsT33(1:DVS) = a22 * dTemp
      ! - get moisture expansion strain
      epsH11(1:DVS) = b11 * dMoist
      epsH22(1:DVS) = b22 * dMoist
      epsH33(1:DVS) = b22 * dMoist
      ! - epsEff = eps_total - eps_temp - epsHyg - epsInes
      eps11(1:DVS) = eps123(1:DVS,1) - epsT11(1:DVS) - epsH11(1:DVS)
      eps22(1:DVS) = eps123(1:DVS,2) - epsT22(1:DVS) - epsH22(1:DVS)
      eps33(1:DVS) = eps123(1:DVS,3) - epsT33(1:DVS) - epsH33(1:DVS)
      eps23(1:DVS) = eps123(1:DVS,4)
      eps13(1:DVS) = eps123(1:DVS,5)
      eps12(1:DVS) = eps123(1:DVS,6) - epsPla(1:DVS) 

      ! Effective stress tensor
      call EvalStresses( &
         E11, E22, G12, G23, v12, v23, &
         0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, & 
         eps11(1:DVS), eps22(1:DVS), eps33(1:DVS), eps23(1:DVS), eps13(1:DVS), eps12(1:DVS), &
         sig11(1:DVS), sig22(1:DVS), sig33(1:DVS), sig23(1:DVS), sig13(1:DVS), sig12(1:DVS) &
         )

      ! Failure criteria / loading functions
      call EvalFailureCriteria( &
         E11, E22, vxT(1:DVS), vxC(1:DVS), vyT(1:DVS), vyC(1:DVS), vnL(1:DVS), nT, gG, vsL(1:DVS), vsT(1:DVS), alpha0, &
         eps11(1:DVS), eps22(1:DVS), sig11(1:DVS), sig22(1:DVS), sig12(1:DVS), &
         phiFT(1:DVS), phiFC(1:DVS), phiMT(1:DVS), phiMC(1:DVS) &
         )

      ! Damage threshold
      call EvalDamageThresholds( &
         flagVisco, eta, dt, &
          phiFT(1:DVS),  phiFC(1:DVS),  phiMT(1:DVS),  phiMC(1:DVS), &
         rFTold(1:DVS), rFCold(1:DVS), rMTold(1:DVS), rMCold(1:DVS), &
         rFTnew(1:DVS), rFCnew(1:DVS), rMTnew(1:DVS), rMCnew(1:DVS), &
         flagFT(1:DVS), flagFC(1:DVS), flagMT(1:DVS), flagMC(1:DVS)  &
         )
            
      ! Damage state
      call EvalDamageState( &
         E11, E22, v12, v21, G12, &
         vxT(1:DVS), vxTO(1:DVS), vxC(1:DVS), vxCO(1:DVS), vyT(1:DVS), vyC(1:DVS), vsL(1:DVS), vsT(1:DVS), vnL(1:DVS), nT, &
         gxT, gxTO, gxC, gxCO, gyT, gyC, gsL, gG, alpha0, kP, & 
         flagPlast(1:DVS), velch(1:DVS), & 
         eps11(1:DVS), eps22(1:DVS), &
         rFTnew(1:DVS), rFCnew(1:DVS), rMTnew(1:DVS), rMCnew(1:DVS), &
         d1(1:DVS), d2(1:DVS), d3(1:DVS), d4(1:DVS), d5(1:DVS), d6(1:DVS) &
         )
      
      ! Initial conditions for damage state
      if( kfl_insdv_sld > 0_ip ) then
         do iv = 1,DVS
            ielem  = lelem_loc(iv)
            d1(iv) = min(d1Max,max(d1(iv),xfiel(kfl_insdv_sld) % a(1,ielem,1)))
            d2(iv) = min(d2Max,max(d2(iv),xfiel(kfl_insdv_sld) % a(2,ielem,1)))
            d3(iv) = min(d3Max,max(d3(iv),xfiel(kfl_insdv_sld) % a(3,ielem,1)))
            d4(iv) = min(d4Max,max(d4(iv),xfiel(kfl_insdv_sld) % a(4,ielem,1)))
            d5(iv) = min(d5Max,max(d5(iv),xfiel(kfl_insdv_sld) % a(5,ielem,1)))
            d6(iv) = min(d6Max,max(d6(iv),xfiel(kfl_insdv_sld) % a(6,ielem,1)))
         end do
      end if
      
      ! Stiffness tensor
      call EvalStresses( &
         E11, E22, G12, G23, v12, v23, &
         d1(1:DVS), d2(1:DVS), d3(1:DVS), d4(1:DVS), d5(1:DVS), d6(1:DVS), &
         eps11(1:DVS), eps22(1:DVS), eps33(1:DVS), eps23(1:DVS), eps13(1:DVS), eps12(1:DVS), &
         sig11(1:DVS), sig22(1:DVS), sig33(1:DVS), sig23(1:DVS), sig13(1:DVS), sig12(1:DVS) &
         )

      ! Set stresses 
      sig123(1:DVS,1) = sig11(1:DVS)
      sig123(1:DVS,2) = sig22(1:DVS)
      sig123(1:DVS,3) = sig33(1:DVS)
      sig123(1:DVS,4) = sig23(1:DVS)
      sig123(1:DVS,5) = sig13(1:DVS)
      sig123(1:DVS,6) = sig12(1:DVS)
         
      ! Set SDV
      sdv(1:DVS,1)  = rFTnew(1:DVS)
      sdv(1:DVS,2)  = rFCnew(1:DVS)
      sdv(1:DVS,3)  = rMTnew(1:DVS)
      sdv(1:DVS,4)  = rMCnew(1:DVS)
      sdv(1:DVS,5)  = d1(1:DVS)
      sdv(1:DVS,6)  = d2(1:DVS)
      sdv(1:DVS,7)  = d3(1:DVS)
      sdv(1:DVS,8)  = d4(1:DVS)
      sdv(1:DVS,9)  = d5(1:DVS)
      sdv(1:DVS,10) = d6(1:DVS)
      sdv(1:DVS,11) = epsPla(1:DVS)
      sdv(1:DVS,12) = epsIne(1:DVS)

   end subroutine evalModelGP
   
   !-----------------------------------------------------------------------
   !> 
   !> @author  A. Quintanas-Corominas 
   !> @date    June, 2022
   !> @brief   Strength reduction based on characteristic element length
   !> @details Strength reduction based on characteristic element length
   !>          
   !-----------------------------------------------------------------------
   
   subroutine reduceStrengths( &
        E11, E22, G12, gxT, gxTO, gxC, gXCO, gyT, gyC, gsL, lenchar, xT, &
        xTO, xC, xCO, yT, yC, sL)
      
      real(rp), intent(in)    :: E11, E22, G12                       !< Elastic
      real(rp), intent(in)    :: gxT, gxTO, gxC, gXCO, gyT, gyC, gsL !< Fracture thoughness
      real(rp), intent(in)    :: lenchar                             !< Characteristic element length
      real(rp), intent(inout) :: xT, xTO, xC, xCO, yT, yC, sL        !< Strengths
      real(rp)                :: minGA, minGB, minGC, lim, d1ch      !< Limits
      
      !-------------------------------------------------------------------
      !
      ! Fibre strengths
      !
      !-------------------------------------------------------------------
      !
      ! Longitudinal - Tension
      !
      minGA = 2.0_rp*(gxT + gxTO)*E11 - lenchar*(xT**2)
      if( minGA < minG_i ) then
         lim = -10.0_rp*E11
         xT  = sqrt(-2.0_rp*(gxT + gxTO)*E11*lim/lenchar/(E11 - lim))
         xTO = xT*gxTO/(gxT + gxTO)
      else
         minGB = 2.0_rp*gXT*E11 - lenchar*xT*(xT - xTO)
         if( minGB < minG_j ) then
            xTO = 1.01_rp*(lenchar*(xT**2) - 2.0_rp*gxT*E11)/(lenchar*xT)
         end if
         d1ch = 2.0_rp*gxT*E11/(xTO*xT*lenchar + 2.0_rp*gxT*E11)
         minGC = 2.0_rp*gXTO*E11*(1.0_rp - d1ch) - lenchar*(xTO**2)
         if( minGC < minG_j ) then
            xTO = 0.95_rp*sqrt(gxTO*E11*(1.0_rp - d1ch)/lenchar)
         end if
      end if
      !
      ! Longitudinal - Compression
      !
      minGA = 2.0_rp*(gxC + gxCO)*E11 - lenchar*(xC**2)
      if( minGA < minG_i ) then
         lim = -10.0_rp*E11
         xC  = sqrt(-2.0_rp*(gxC + gxCO)*E11*lim/lenchar/(E11 - lim))
         xCO = xC*gxCO/(gxC + gxCO)
      else
         minGB = 2.0_rp*gXC*E11 - lenchar*xC*(xC - xCO)
         if( minGB < minG_j ) then
            xCO = 1.01_rp*(lenchar*(xC**2) - 2.0_rp*gxC*E11)/(lenchar*xC)
         end if
         d1ch = 2.0_rp*gxC*E11/(xCO*xC*lenchar + 2.0_rp*gxC*E11)
         minGC = 2.0_rp*gXCO*E11*(1.0_rp - d1ch) - lenchar*(xCO**2)
         if( minGC < minG_j ) then
            xCO = 0.95_rp*sqrt(gxCO*E11*(1.0_rp - d1ch)/lenchar)
         end if
      end if
      
      !-------------------------------------------------------------------
      !
      ! Matrix strengths
      !
      !-------------------------------------------------------------------
      !
      ! Transversal - Tension
      !
      minGA = 2.0_rp*gyT*E22 - lenchar*(yT**2)
      if( minGA < minG_j ) then
         lim = -10.0_rp*E22
         yT = sqrt(-2.0_rp*gyT*E22*lim/lenchar/(E22 - lim))
      endif
      !
      ! Transversal - Compression
      !
      minGA = 2.0_rp*gyC*E22 - lenchar*(yC**2)
      if( minGA < minG_j ) then
         lim = -10.0_rp*E22
         yC = sqrt(-2.0_rp*gyC*E22*lim/lenchar/(E22 - lim))
      endif
      !
      ! In-plane Shear
      !
      minGA = 2.0_rp*gsL*G12 - lenchar*(sL**2)
      if( minGA < minG_j ) then
         lim = -10.0_rp*G12
         sL = sqrt(-2.0_rp*gsL*G12*lim/lenchar/(G12 - lim))
      end if

   end subroutine reduceStrengths

   !-----------------------------------------------------------------------
   !> 
   !> @author  A. Quintanas-Corominas 
   !> @date    June, 2022
   !> @brief   Material properties 
   !> @details Material properties 
   !>          
   !-----------------------------------------------------------------------
   
   subroutine GetProperties( &
      props, E11, E22, v12, v21, v23, G12, G23, &
      xT, xTO, xC, xCO, yT, yC, sL, &
      alpha0, gxT, gxTO, gxC, gxCO, gyT, gyC, gsL, gG, &
      sP, kP, eta, a11, a22, incrTemp, b11, b22, incrMois, &
      thick,mugsL,yBT,yBC,mufxc,E1c)

      real(rp), intent(in)  :: props(nprops154)
      real(rp), intent(out) :: E11, E22, v12, v21, v23, G12, G23       !< Elastic
      real(rp), intent(out) :: xT, xTO, xC, xCO, yT, yC, sL, alpha0    !< Strength
      real(rp), intent(out) :: gxT, gxTO, gxC, gxCO, gyT, gyC, gsL, gG !< Fracture thoughness
      real(rp), intent(out) :: sP, kP                                  !< Plasticity
      real(rp), intent(out) :: eta                                     !< Viscosity
      real(rp), intent(out) :: a11, a22, incrTemp, b11, b22, incrMois  !< Hygrothermal
      real(rp), intent(out) :: thick                                   !< Out-of-plane thickness
      real(rp), intent(out) :: mugsL, mufxc                            !< Friction coefficients
      real(rp), intent(out) :: yBT, yBC                                !< Biaxial strenghts
      real(rp), intent(out) :: E1c                                     !< Fiber Young Modulus compression
    
      ! Get from properties vectors
      E11    = props(1)
      E22    = props(2)
      v12    = props(3)
      v23    = props(4)
      v21    = E22*v12/E11
      G12    = props(5)
      G23    = E22/(1.0_rp + v23)/2.0_rp
      xT     = props(6)
      xTO    = props(7)*xT   ! fxT * xT
      xC     = props(8)
      xCO    = props(9)*xC   ! fxC * xC
      yT     = props(10)
      yC     = props(11)
      sL     = props(12)
      alpha0 = props(13)

      ! Fracture toughness
      gxT  = props(15)*props(14)             ! gxT  = fGT*G1T
      gxTO = (1.0_rp - props(15))*props(14)  ! gxTO = (1 - fG1C)*G1T
      gxC  = props(17)*props(16)             ! gxC  = fGC*G1C
      gxCO = (1.0_rp - props(17))*props(16)  ! gxCO = (1 - fG1C)*G1C
      gyT  = props(18)
      gyC  = props(19)
      gsL  = props(20)
      gG   = gyT/gsL

      ! Plastic properties
      sP = props(21)
      kP = props(22)

      ! Viscosity properties
      eta = props(23)

      ! Hygrothermal elastic properties
      a11 = props(24)
      a22 = props(25)
      incrTemp = props(26)
      b11 = props(27)
      b22 = props(28)
      incrMois = props(29)
      
      ! Max. damage variables
      d1Max = abs(props(30))                                     
      d2Max = abs(props(31))                                 
      d3Max = abs(props(32))
      d4Max = abs(props(33))
      d5Max = abs(props(34))
      d6Max = abs(props(35))

      ! Thickness
      thick = props(36)

      ! Friction GSL and fxC
      mugsL = props(37)
      mufxc = props(40)
      
      ! YBT and YBC
      yBT = props(38)
      yBC = props(39)

      ! E1c
      E1c = props(41)
      
   end subroutine GetProperties

   !-----------------------------------------------------------------------
   !> 
   !> @author  A. Quintanas-Corominas
   !> @date    June, 2022
   !> @brief   Plasticity strain under shear
   !> @details Get plasticicty strain 1-2
   !>
   !-----------------------------------------------------------------------
   
   elemental subroutine GetPlasticStrain12( &
      G12, sP, kP, eps12, epsPlast, epsInela, flagD, flagP)

      real(rp),    intent(in)    :: G12
      real(rp),    intent(in)    :: sP, kP
      real(rp),    intent(in)    :: eps12
      real(rp),    intent(inout) :: epsPlast
      real(rp),    intent(inout) :: epsInela
      logical(lg), intent(in)    :: flagD
      logical(lg), intent(inout) :: flagP
      real(rp)                   :: ftrial, depsInela, epsTrial
      
      ! Plasticity only before damage
      if( flagP .and. .not. flagD )then
         ! Update plasticity
         depsInela = 0.0_rp
         ftrial = abs(eps12 - epsPlast) - sP / G12 - kP * epsInela
         if( ftrial >= 0.0_rp )then
            depsInela = ftrial / (1.0_rp + kP)
            epsInela = epsInela + depsInela
            epsTrial = eps12 - epsPlast
            epsPlast = epsPlast + depsInela * sign(1.0_rp, epsTrial)
            flagP = .true.
         else
            flagP = .false.
         end if
      else
         flagP = .false.
      endif

   end subroutine GetPlasticStrain12

   !-----------------------------------------------------------------------
   !> 
   !> @author  A. Quintanas-Corominas
   !> @date    June, 2022
   !> @brief   Evaluate failure criteria
   !> @details Failure criteria based on LaRC03 and LaRC04
   !> @todo    Fibre compression to be updated with Furtado/Arteiro versions
   !> 
   !-----------------------------------------------------------------------
   
   elemental subroutine EvalFailureCriteria( &
      E11, E22, xT, xC, yT, yC, nL, nT, gG, sL, sT, alpha0,  &
      eps11, eps22, sig11, sig22, sig12, &
      phiFT, phiFC, phiMT, phiMC)

      real(rp), intent(in)  :: E11, E22                       !< Elastic material properties
      real(rp), intent(in)  :: xT, xC, yT, yC, nL, nT, sL, sT !< Strength properties
      real(rp), intent(in)  :: gG                             !< Fracture toughness ratio
      real(rp), intent(in)  :: alpha0                         !< Fracture angle
      real(rp), intent(in)  :: eps11, eps22                   !< Longitudinal strain
      real(rp), intent(in)  :: sig11, sig22, sig12            !< Effective stresses
      real(rp), intent(out) :: phiFT                          !< Fiber tensile
      real(rp), intent(out) :: phiFC                          !< Fiber compression
      real(rp), intent(out) :: phiMT                          !< Matrix tensile
      real(rp), intent(out) :: phiMC                          !< Matrix compression
      real(rp)              :: auxs1                          !< Auxiliary variables
      real(rp)              :: phiC
      real(rp)              :: theta, s22T, s12T, teffT, teffL

      !-------------------------------------------------------------------
      !
      ! Initializations
      !
      !-------------------------------------------------------------------
      
      phiFT = 0.0_rp
      phiFC = 0.0_rp
      phiMT = 0.0_rp
      phiMC = 0.0_rp

      !-------------------------------------------------------------------
      !
      ! Fibre failure
      !
      !-------------------------------------------------------------------
      
      if( sig11 < 0.0_rp ) then
         ! LaRC04 criteria
         auxs1 = sL/xC + nL
         phiC = atan((1.0_rp - sqrt(max(0.0_rp, 1.0_rp - 4.0_rp*auxs1*sL/xC)))/(2.0_rp*auxs1))
         s22T =  sig11*(sin(phiC)**2) + sig22*(cos(phiC)**2) - 2.0_rp*abs(sig12)*sin(phiC)*cos(phiC)
         s12T = (sig22 - sig11)*sin(phiC)*cos(phiC) + abs(sig12)*(cos(phiC)**2 - sin(phiC)**2)
         phiFC = max(0.0_rp, abs(s12T) + nL*s22T)/sL
         phiFC = min(phiFC, abs(eps11)*E11/yC)
         phiFT = 0.0_rp
      else
         ! LaRC03 criteria: max strain
         phiFT = eps11*E11/xT      
         phiFC = 0.0_rp
      end if

      !-------------------------------------------------------------------
      !
      ! Matrix failure
      !
      !-------------------------------------------------------------------
      
      if( sig22 < -zeroapp ) then
         ! LaRC03
         theta = atan(-abs(sig12)/(sig22*sin(alpha0)))
         teffT = sig22*cos(alpha0)*(nT*cos(alpha0)*cos(theta) - sin(alpha0))
         teffT = max(0.0_rp, teffT)
         teffL = cos(alpha0)*(abs(sig12) + nL*cos(alpha0)*sin(theta)*sig22)
         teffL = max(0.0_rp, teffL)
         phiMC = sqrt((teffT/sT)**2 + (teffL/sL)**2)
         phiMC = min(phiMC,-E22*eps22/yT)
         phiMT = abs(sig12) + nL*sig22
         phiMT = max(0.0_rp, phiMT)/sL
      else
         ! LaRC03
         phiMT = (1.0_rp - gG)*(sig22/yT) + gG*(sig22/yT)**2 + (sig12/sL)**2
         if(phiMT > 0.0_rp ) then
            phiMT = sqrt(phiMT)
            phiMC = 0.0_rp
         else
            phiMT = 0.0_rp
            phiMC = 0.0_rp
         end if
        
      end if

   end subroutine EvalFailureCriteria

   !-----------------------------------------------------------------------
   !> 
   !> @author  A. Quintanas-Corominas
   !> @date    June, 2022
   !> @brief   Damage thresholds
   !> @details Update of the damage limits
   !>
   !> @note    Revise Matrix tensile and compression are different respect to Furtado
   !-----------------------------------------------------------------------
   
   elemental subroutine EvalDamageThresholds( &
      flagVisco, eta, dt, &
      phiFT, phiFC, phiMT, phiMC, &
      rFTold, rFCold, rMTold, rMCold, &
      rFTnew, rFCnew, rMTnew, rMCnew, &
      flagFT, flagFC, flagMT, flagMC )

      logical(lg), intent(in)  :: flagVisco                      !< Flag for viscosity
      real(rp),    intent(in)  :: eta, dt                        !< Viscosity parameters
      real(rp),    intent(in)  :: phiFT, phiFC, phiMT, phiMC     !< Loading state variables
      real(rp),    intent(in)  :: rFTold, rFCold, rMTold, rMCold !< Damage threshold variables (old)
      real(rp),    intent(out) :: rFTnew, rFCnew, rMTnew, rMCnew !< Damage threshold variables (new)
      logical(lg), intent(out) :: flagFT, flagFC, flagMT, flagMC !< Flags damage mechanisms

      !-------------------------------------------------------------------
      !
      ! Fibre failure
      !
      !-------------------------------------------------------------------

      ! Fiber - Tensile
      if( phiFT > rFTold ) then
         flagFT = .true.
         if( flagVisco ) then
            rFTnew = rFTold*eta/(eta + dt) + phiFT*dt/(eta + dt)
         else
            rFTnew = phiFT
         end if
      else
         flagFT = .false.
         rFTnew = rFTold
      end if

      ! Fiber - Compresssion
      if( phiFC > rFCold ) then
         flagFC = .true.
         if( flagVisco ) then
            rFCnew = rFCold*eta/(eta + dt) + phiFC*dt/(eta + dt)
         else
            rFCnew = phiFC
         end if
         ! Coupling with the longitudinal tensile damage
         if( rFCnew > rFTnew ) then
            rFTnew = rFCnew
         end if
      else
         flagFC = .false.
         rFCnew = rFCold
      end if
      
      !-------------------------------------------------------------------
      !
      ! Matrix failure
      !
      !-------------------------------------------------------------------

      ! Matrix - Tensile
      if( phiMT>rMTold ) then
         flagMT = .true.
         rMTnew = phiMT
      else
         flagMT = .false.
         rMTnew = rMTold
      end if

      ! Matrix - Compression (transversal damage at 53º)
      if( phiMC>rMCold ) then
         flagMC = .true.
         rMCnew = phiMC
         ! Coupling with the transverse tensile damage
         if( rMCnew > rMTnew ) then
            rMTnew = rMCnew
         end if
      else
         flagMC = .false.
         rMCnew = rMCold
      end if

   end subroutine EvalDamageThresholds

   !-----------------------------------------------------------------------
   !> 
   !> @author  A. Quintanas-Corominas
   !> @date    June, 2022
   !> @brief   Damage thresholds
   !> @details Update of the damage limits
   !> @note    To revise Matrix shear/compression with respect to furtado
   !>
   !-----------------------------------------------------------------------
   
   elemental subroutine EvalDamageState( &
      E11, E22, v12, v21, G12, &
      xT, xTO, xC, xCO, yT, yC, sL, sT, nL, nT, &
      gxT, gxTO, gxC, gxCO, gyT, gyC, gsL, gG, alpha0, kP, &
      flagPlast, lenchar, &
      eps11, eps22, &
      rFT, rFC, rMT, rMC, &
      d1, d2, d3, d4, d5, d6 )

     logical(lg), intent(in)  :: flagPlast                                !< Plasticity flag
     real(rp),    intent(in)  :: E11, E22, v12, v21, G12                  !< Elastic properties
     real(rp),    intent(in)  :: xT, xTO, xC, xCO, yT, yC, sL, sT, nL, nT !< Strength properties
     real(rp),    intent(in)  :: gxT, gxTO, gxC, gxCO, gyT, gyC, gsL, gG  !< Energy properties
     real(rp),    intent(in)  :: alpha0                                   !< Fracture angle
     real(rp),    intent(in)  :: kP                                       !< Plastic property
     real(rp),    intent(in)  :: lenchar                                  !< Characteristic element length
     real(rp),    intent(in)  :: eps11, eps22                             !< Longitudinal and transversal strain
     real(rp),    intent(in)  :: rFT, rFC, rMT, rMC                       !< Damage threshold variables
     real(rp),    intent(out) :: d1, d2, d3, d4, d5, d6                   !< Damage state variables
     real(rp)                 :: dFTch, rFTch, mxT, nxT                   !< Longitudinal direction tensile
     real(rp)                 :: dFCch, rFCch, mxC, nxC                   !< Longtindunal direction compression
     real(rp)                 :: myT, nyT                                 !< Transverse direction tensil
     real(rp)                 :: myC, nyC                                 !< Tranvserse direction compression
     real(rp)                 :: msL, nsL                                 !< Shear direction
     real(rp)                 :: s1, s2, s3, a1, a2, a3, auxs1, invGreEff !< Auxiliar parameters
     real(rp)                 :: AYC, BYC                                 !< Auxiliar parameters
     real(rp)                 :: d1T, d1C, d2T, d2C, d66                  !< Temporal damage variables
     real(rp)                 :: phiC                                     !< Auxilair material properties

     !-------------------------------------------------------------------
     !
     ! Fibre tension
     !
     !-------------------------------------------------------------------
     
     d1T = 0.0_rp
     if( rFT > 1.0_rp ) then
        dFTch = 2.0_rp*gXT*E11/(xT*xTO*lenchar + 2.0_rp*gXT*E11)
        rFTch = xTO/(xT*(1.0_rp - dFTch))
        if (rFT <= rFTch ) then 
           ! First slope
           mxT = lenchar*E11*xT*(xT - xTO)/(lenchar*xT*(xT - xTO) - 2.0_rp*gXT*E11)
           nxT = xT*(1.0_rp - mxT/E11)
        else
           ! Second slope
           mxT = lenchar*E11*(1.0_rp - dFTch)*(xTO**2)/(lenchar*(xTO**2) - &
                 2.0_rp*gXTO*E11*(1.0_rp - dFTch))
           nxT = xTO*(1.0_rp - mxT/(E11*(1.0_rp - dFTch)))
        end if
        !
        ! Damage variable
        d1T = 1.0_rp - mxT/E11 - nxT/(XT*rFT)
        d1T = max(0.0_rp,d1T)
        d1T = min(1.0_rp,d1T)
     end if

     !-------------------------------------------------------------------
     !
     ! Fibre compresion 
     !
     !-------------------------------------------------------------------
     
     d1C = 0.0_rp
     if( rFC > 1.0_rp ) then
        dFCch  = 2.0_rp*gXC*E11/(xC*xCO*lenchar + 2.0_rp*gXC*E11)
        mxC = lenchar*xC*E11*(xC - xCO)/(lenchar*xC*(xC - xCO) - 2.0_rp*gXC*E11)
        nxC = XC*(1.0_rp - mxC/E11)
        auxs1 = abs(sL/xC) + nL
        phiC = atan((1.0_rp - sqrt(max(0.0_rp, 1.0_rp - 4.0_rp*auxs1*sL/xC)))/(2.0_rp*auxs1))
        s1 = E11 - mxC*v12*v21
        s2 = (E11 - mxC)*v21
        a1 = (s2 - s1)*sin(phiC)*cos(phiC) + nL*(s1*(sin(phiC)**2) + s2*(cos(phiC))**2)
        a2 = nxC*v21*((1.0_rp - v12)*sin(phiC)*cos(phiC) + nL*((cos(phiC))**2 + v12*(sin(phiC)**2)))
        rFCch = ((-xCO/E11/(1.0_rp - dFCch))*a1 + a2)/((1.0_rp - v12*v21)*sL)
        if( rFC <= rFCch ) then
           ! First tram
           invGreEff = a1/(E11*(sL*(1.0_rp - v12*v21)*rFC - a2))
        else
           ! Second tram
           mxC = lenchar*E11*(1.0_rp - dFCch)*(xCO**2)/(lenchar*(xCO**2) - 2.0_rp*gXCO*E11*(1.0_rp - dFCch))
           nxC = xCO*(1.0_rp - mxC/(E11*(1.0_rp - dFCch)))
           s1 = E11 - mxC*v12*v21
           s2 = (E11 - mxC)*v21
           a1 = (s2 - s1)*sin(phiC)*cos(phiC) + nL*(s1*(sin(phiC)**2) + s2*(cos(phiC))**2)
           a2 = nxC*v21*((1.0_rp - v12)*sin(phiC)*cos(phiC) + nL*((cos(phiC))**2 + v12*(sin(phiC)**2)))
           invGreEff = a1/(E11*(sL*(1.0_rp - v12*v21)*rFC - a2))
        end if
        !
        ! Damage variable
        d1C = 1.0_rp - mxC/E11 + nxC*invGreEff
        d1C = max(0.0_rp,d1C)
        d1C = min(1.0_rp,d1C)
     end if

     !-------------------------------------------------------------------
     !
     ! Matrix tension (to be revised!)
     !
     !-------------------------------------------------------------------
     
     d2T = 0.0_rp
     if( rMT > 1.0_rp ) then
        myT = lenchar*(yT**2)/(lenchar*(yT**2) - 2.0_rp*gyT*E22)
        nyT = yT*(1.0_rp - myT)/E22
        s1 = E22*(1.0_rp - myT*v12*v21)/yT/(1.0_rp - v12*v21)
        s2 = E22*nyT*v12*v21/yT/(1.0_rp - v12*v21)   
        a1 = gG*(s1**2)
        a2 = s1*(gG - 1.0_rp - 2.0_rp*gG*s2)
        a3 = s2*(s2*gG + 1.0_rp - gG)
        s3 = sqrt(a2**2 - 4.0_rp*a1*a3 + 4.0_rp*a1*rMT**2)
        !
        ! Damage variable
        d2T = 1.0_rp - myT - nyT*(2.0_rp*a1)/(a2 + s3)
        d2T = max(0.0_rp,d2T)
        d2T = min(1.0_rp,d2T)
     end if

     !-------------------------------------------------------------------
     !
     ! Matrix compression
     !
     !-------------------------------------------------------------------
     
     d2C = 0.0_rp
     if( rMC > 1.0_rp ) then
        myC = lenchar*(yC**2)*E22/(lenchar*(yC**2) - 2.0_rp*gyC*E22)
        nyC = yC*(1.0_rp - myC/E22)
        AYC = -nyC*v12*v21
        BYC = -(1.0_rp - v12*v21)/(cos(alpha0)*(sin(alpha0) - nT*cos(alpha0)))
        !
        ! Damage variable
        d2C = 1.0_rp - myC/E22 + nyC*(E22 - myC*v12*v21)/(E22*(AYC + BYC*rMC*sT))
        d2C = max(0.0_rp,d2C)
        d2C = min(1.0_rp,d2C)
     end if
     
     !-------------------------------------------------------------------
     !
     ! Matrix shear 
     !
     !-------------------------------------------------------------------

     d66 = 0.0_rp
     if( rMT > 1.0_rp ) then
        msL = lenchar*(sL**2)*G12/(lenchar*(sL**2) - 2.0_rp*gsL*G12)
        if( flagPlast ) then
           msL = msL*(1.0_rp + kP)/kP
        end if
        nsL = sL*(1.0_rp - msL/G12)
        !
        ! Damage variable
        d66 = 1.0_rp - msL/G12 - nsL/(sL*rMT)
        d66 = max(0.0_rp,d66)
        d66 = min(1.0_rp,d66)
     end if
     
     !-------------------------------------------------------------------
     !
     ! Active damage modes and dependent damage variables
     !
     !-------------------------------------------------------------------

     ! Fibre
     if(      eps11 > (d2C - 1.0_rp)*v21*eps22 ) then
        d1 = d1T
     else if( eps11 < (d2T - 1.0_rp)*v21*eps22 ) then
        d1 = d1C
     else
        d1 = 0.5_rp*(d1T + d1C)
     end if
     ! Matrix
     if( eps22 > (d1 - 1.0_rp)*v12*eps11) then
        d2 = d2T
     else
        d2 = d2C
     end if
     
     d3 = d3Max - (d2Max - d2C)*(d1Max - d1C)
     d4 = d66
     d5 = d1T
     d6 = d6Max - (d6Max - d66)*(d1Max - d1T)
     
     !-------------------------------------------------------------------
     !
     ! User limits for damage variables
     !
     !-------------------------------------------------------------------
     
     d1 = min(d1Max,d1)
     d2 = min(d2Max,d2)
     d3 = min(d3Max,d3)
     d4 = min(d4Max,d4)
     d5 = min(d5Max,d5)
     d6 = min(d6Max,d6)
    
   end subroutine EvalDamageState
   
   !-----------------------------------------------------------------------
   !> 
   !> @author  A. Quintanas-Corominas
   !> @date    June, 2022
   !> @brief   Constitutive tensor and stress evaluation
   !> @details Constitutive tensor and stress evaluation
   !> 
   !-----------------------------------------------------------------------
   
   elemental subroutine EvalStresses( &
      E11, E22, G12, G23, v12, v23, &
      d1, d2, d3, d4, d5, d6, &
      eps11, eps22, eps33, eps23, eps13, eps12, &
      sig11, sig22, sig33, sig23, sig13, sig12)

      real(rp), intent(in)  :: E11, E22, G12, G23, v12, v23             !< Material properties
      real(rp), intent(in)  :: d1, d2, d3, d4, d5, d6                   !< Damage variables
      real(rp), intent(in)  :: eps11, eps22, eps33, eps23, eps13, eps12 !< Material strains
      real(rp), intent(out) :: sig11, sig22, sig33, sig23, sig13, sig12 !< Material stresses
      real(rp)              :: H11, H22, H33, H12, H23, delta           !< Auxiliary variables
      real(rp)              :: C11, C12, C13, C21, C22, C23             !< Tensor values
      real(rp)              :: C31, C32, C33, C44, C55, C66             !< Tensor values 

      !-------------------------------------------------------------------
      !
      ! Constitutive tensor
      !
      !-------------------------------------------------------------------
      
      ! Initializations
      C11 = 0.0_rp
      C22 = 0.0_rp
      C33 = 0.0_rp
      C23 = 0.0_rp
      C13 = 0.0_rp
      C12 = 0.0_rp
      C21 = 0.0_rp
      C31 = 0.0_rp
      C32 = 0.0_rp
      C44 = 0.0_rp
      C55 = 0.0_rp
      C66 = 0.0_rp
      ! Inversion flexibility matrix (3x3)
      if(      d1>=damMax .and. d2<damMax .and. d3<damMax ) then
         H22 = 1.0_rp/E22/(1.0_rp - d2)
         H33 = 1.0_rp/E22/(1.0_rp - d3)
         H23 = -V23/E22
         C22 =  H33/(H33*H22 - H23**2)
         C33 =  H22/(H33*H22 - H23**2)
         C23 = -H23/(H33*H22 - H23**2)
      else if( d1<damMax .and. d2>=damMax .and. d3<damMax ) then
         H11 = 1.0_rp/E11/(1.0_rp - d1)
         H33 = 1.0_rp/E22/(1.0_rp - d3)
         H12 = -V12/E11
         H23 = -V23/E22
         C11 =  H33/(H33*H11 - H12**2)
         C33 =  H11/(H33*H11 - H12**2)
         C13 = -H12/(H33*H11 - H12**2)
      else if( d1<damMax .and. d2<damMax .and. d3>=damMax ) then
         H11 = 1.0_rp/E11/(1.0_rp - d1)
         H22 = 1.0_rp/E22/(1.0_rp - d2)
         H12 = -V12/E11
         H23 = -V23/E22
         C11 =  H22/(H22*H11 - H12**2)
         C22 =  H11/(H22*H11 - H12**2)
         C12 = -H12/(H22*H11 - H12**2)
      else if( d1<damMax .and. d2>=damMax .and. d3>=damMax ) then
         C11 = E11*(1.0_rp - d1)
      else if( d1>=damMax .and. d2<damMax .and. d3>=damMax ) then
         C22 = E22*(1.0_rp - d2)
      else if( d1>=damMax .and. d2>=damMax .and. d3<damMax ) then
         C33 = E22*(1.0_rp - d3)
      else if( d1<damMax .and. d2<damMax .and. d3<damMax   ) then
         H11 = 1.0_rp/E11/(1.0_rp - d1)
         H22 = 1.0_rp/E22/(1.0_rp - d2)
         H33 = 1.0_rp/E22/(1.0_rp - d3)
         H12 = -V12/E11
         H23 = -V23/E22
         delta = H11*(H23**2 - H22*H33) + (H22 + H33 - 2.0_rp*H23)*(H12**2)
         C11 = (H23**2 - H22*H33)/delta
         C12 =  H12*(H33 - H23)/delta
         C13 =  H12*(H22 - H23)/delta
         C22 = (H12**2 - H11*H33)/delta
         C23 = (H11*H23 - H12**2)/delta
         C33 = (H12**2 - H11*H22)/delta
      end if
      C21 = C12
      C31 = C13
      C32 = C23

      ! Inversion of shear contributions
      C44 = (1.0_rp - d4)*G23
      C55 = (1.0_rp - d5)*G12
      C66 = (1.0_rp - d6)*G12
      
      !-------------------------------------------------------------------
      !
      ! Stresses
      !
      !-------------------------------------------------------------------
      
      sig11 = C11 * eps11 + C12 * eps22 + C13 * eps33
      sig22 = C21 * eps11 + C22 * eps22 + C23 * eps33
      sig33 = C31 * eps11 + C32 * eps22 + C33 * eps33
      sig23 = C44 * eps23
      sig13 = C55 * eps13
      sig12 = C66 * eps12

   end subroutine EvalStresses


end module mod_sld_vect_stress_model_154

