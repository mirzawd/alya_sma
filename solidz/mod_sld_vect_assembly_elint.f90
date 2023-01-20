!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    mod_sld_vect_assembly_elint.f90
!> @author  Gerard Guillamet (gerard.guillamet@bsc.es)
!> @author  Adria Quintanas (adria.quintanas@udg.edu)
!> @date    September 2022
!> @brief   Vectorised version for the assembly of interface 
!>          elements
!>
!> @warning This version is optimezed and only works for Explicit
!>          and 3D hexahedron elements. The vectorized cohesive law           
!>          corresponds to Turon et al. 2018.
!>
!> @details
!>          References:\n
!>
!>          A. Turon, E.V. Gonzalez, C. Sarrado, G. Guillamet, P. Maimi. Accurate
!>          simulation of delamination under mixed-mode loading using a cohesive
!>          model with a mode-dependent penalty stiffness. Composite
!>          structures, 2018.\n
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
!>            - To figure out if we can vectorize the cohesive law  
!> @}
!------------------------------------------------------------------

module mod_sld_vect_assembly_elint

#include "def_solidz_vect_dim.inc" 
  use def_kintyp, only : ip, rp, lg
  use def_solidz, only : ndofn_sld

  implicit none
  
  integer(ip), parameter           :: nprop   = 9_ip
  integer(ip), parameter           :: closing = 0_ip                        
  integer(ip), parameter           :: opening = 1_ip                         
  integer(ip), parameter           :: elastic = 0_ip                         
  integer(ip), parameter           :: damage  = 1_ip                          
  integer(ip), parameter           :: broken  = 2_ip                              

  real(rp), parameter              :: tolPe   = 1.0e-12_rp
  
  integer(ip), protected           :: lelem_loc(DVS)
  integer(ip), protected           :: lmate_loc(DVS)
  integer(ip), protected           :: lelch_loc(DVS)
  logical(lg), protected           :: lmask_loc(DVS)
  integer(ip), protected,  pointer :: lnods_loc(:,:)

  public                           :: sld_vect_element_operations_ELINT
  
contains
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  gguillam and aquintana
  !> @date    2022-10-04
  !> @brief   Element operations for interface element 
  !> @details Element operations for interface element (3D HEX08 element)
  !> 
  !-----------------------------------------------------------------------
  
  subroutine sld_vect_element_operations_ELINT(itask,pnode,list_elements)

    use def_elmtyp,      only : ELINT
    use def_master,      only : ITER_K, ITER_K_STATE, TIME_N_STATE
    use def_master,      only : rhsid
    use def_master,      only : displ, dtime
    use def_domain,      only : lnods, lelch, lmate
    use def_domain,      only : coord
    use mod_matrix,      only : matrix_assemble_element_RHS
    use def_solidz,      only : parch_sld
    use def_solidz,      only : svegm_sld 
    use def_solidz,      only : kfl_fixno_sld
    use def_solidz,      only : finte_sld
    use def_solidz,      only : frxid_sld
    use def_solidz,      only : kfl_fixno_sld
    !
    ! Input and output variables
    !
    integer(ip), intent(in)    :: itask
    integer(ip), intent(in)    :: pnode
    integer(ip), intent(in)    :: list_elements(DVS)
    !
    ! Element matrices and vectors
    !
    real(rp)    :: elrhs(DVS,pnode*ndofn_sld)
    real(rp)    :: elfin(DVS,pnode*ndofn_sld)
    real(rp)    :: elfrx(DVS,pnode*ndofn_sld)
    !
    ! Gather
    !
    real(rp)    :: elcod(DVS,3,pnode)
    real(rp)    :: eldis(DVS,3,pnode)
    !
    ! Interface variables
    !
    real(rp)          :: sfcoo(DVS,3,pnode/2)
    real(rp)          :: sfjum(DVS,pnode*ndofn_sld/2)
    real(rp)          :: properties(nprop)
    !
    ! Interface element gauss point variables
    !
    real(rp), pointer :: gpsdvo(:,:,:)
    real(rp), pointer :: gpsdvn(:,:,:)
    real(rp)          :: gprota(DVS,3,3,pnode/2)       
    real(rp)          :: gpshap(DVS,pnode/2,pnode/2)                   
    real(rp)          :: gpderi(DVS,2,pnode/2,pnode/2)               
    real(rp)          :: gpbmat(DVS,3,pnode*ndofn_sld/2,pnode/2)                 
    real(rp)          :: gpjump(DVS,3,pnode/2)                          
    real(rp)          :: gptrac(DVS,3,pnode/2)                                           
    real(rp)          :: gpweig(DVS,pnode/2)                                 
    real(rp)          :: gparea(DVS,pnode/2)
    
    real(rp)          :: t, s
    real(rp)          :: r1(DVS,3), r2(DVS,3), r3(DVS,3), rM(DVS)
    real(rp),  parameter :: gpposi(2,4) = reshape([-1.0_rp, 1.0_rp, 1.0_rp,-1.0_rp, &
                                                   -1.0_rp,-1.0_rp, 1.0_rp, 1.0_rp],[2,4])
    !
    ! Indices and dimensions
    !
    integer(ip) :: idime, jdime, idofn, ivect
    integer(ip) :: ielem, ipoin, inode, jnode, igaus
    integer(ip) :: pmate
    integer(ip) :: sgaus, snode, sdime, pevah, peva1

    !--------------------------------------------------------------------
    !
    ! Initializations
    !
    !--------------------------------------------------------------------

    ielem = list_elements(1)
    pmate = lmate(ielem)

    !--------------------------------------------------------------------
    !
    ! Gather: global to local
    !
    !--------------------------------------------------------------------

    allocate(lnods_loc(DVS,pnode))

    do ivect = 1,DVS
       ielem = abs(list_elements(ivect))
       if( ielem /= 0 ) then
          lelem_loc(ivect)         = list_elements(ivect)
          lnods_loc(ivect,1:pnode) = lnods(1:pnode,ielem)
          lelch_loc(ivect)         = lelch(ielem)
          lmate_loc(ivect)         = lmate(ielem)
          lmask_loc(ivect)         = .true.
       else
          lelem_loc(ivect)         = list_elements(1)
          lnods_loc(ivect,1:pnode) = lnods(1:pnode,1)
          lelch_loc(ivect)         = lelch(1)
          lmate_loc(ivect)         = lmate(1)
          lmask_loc(ivect)         = .false.
       end if
    end do

    !--------------------------------------------------------------------
    !
    ! Get elemental residual (forces) and matrix
    !
    !--------------------------------------------------------------------
    !
    ! Initializations
    !
    elrhs(:,:) = 0.0_rp
    elfin(:,:) = 0.0_rp
    elfrx(:,:) = 0.0_rp
    !
    ! Gather some variables
    !
    do ivect = 1,DVS
       do inode = 1,pnode
          ipoin = lnods_loc(ivect,inode)
          elcod(ivect,1:3,inode) = coord(1:3,ipoin)
          eldis(ivect,1:3,inode) = displ(1:3,ipoin,ITER_K)
       end do
    end do
    !
    ! Read material properties
    !
    properties(1:nprop) = parch_sld(1:nprop,pmate)
    ! 
    ! Surface dimensions
    !
    sgaus = int(pnode/2,ip)           ! Mid-surface gauss points (It has to be equal to half of the number of nodes)
    snode = int(pnode/2,ip)           ! Mid-surface number of nodes
    sdime = 2_ip                      ! Mid-surface number of dimensions
    pevah = int(pnode*ndofn_sld/2,ip) ! Mid-surface number of degree of freedom
    !
    ! Midsurface coordinates and displacement jump
    !
    do inode = 1,snode
       jnode = inode + snode
       sfcoo(1:DVS,1,inode) = 0.5_rp*(elcod(1:DVS,1,inode) + elcod(1:DVS,1,jnode) + eldis(1:DVS,1,inode) + eldis(1:DVS,1,jnode))
       sfcoo(1:DVS,2,inode) = 0.5_rp*(elcod(1:DVS,2,inode) + elcod(1:DVS,2,jnode) + eldis(1:DVS,2,inode) + eldis(1:DVS,2,jnode))
       sfcoo(1:DVS,3,inode) = 0.5_rp*(elcod(1:DVS,3,inode) + elcod(1:DVS,3,jnode) + eldis(1:DVS,3,inode) + eldis(1:DVS,3,jnode))
    end do
    !   - relative displacement / displacement jump
    do inode = 1,snode
       jnode = inode + snode
       do idime = 1,3
          idofn = 3*(inode - 1) + idime
          sfjum(1:DVS,idofn) = eldis(1:DVS,idime,jnode) - eldis(1:DVS,idime,inode)
       end do
    end do

    !--------------------------------------------------------------------
    !
    ! Gauss point calculations
    !
    !--------------------------------------------------------------------
    !
    ! Get GP state old variables
    !
    allocate(gpsdvo(DVS,2,sgaus))
    allocate(gpsdvn(DVS,2,sgaus))
    do concurrent ( ivect = 1:DVS )
       !if( lmask_loc(ivect) ) then
          gpsdvo(ivect,1:2,1:sgaus) = svegm_sld(lelem_loc(ivect))%a(1:2,1:sgaus,TIME_N_STATE)
       !else
          !gpsdvo(ivect,1:2,1:sgaus) = 0.0_rp
       !end if
    end do
    !
    ! Shape functions and derivatives
    !
    gpshap(1:DVS,:,:)   = 0.0_rp
    gpderi(1:DVS,:,:,:) = 0.0_rp
    do igaus = 1,sgaus 
       !
       ! Newton Cotes
       !
       ! GP Weigth
       gpweig(1:DVS,igaus) = 1.0_rp
       t = gpposi(1,igaus)
       s = gpposi(2,igaus)
       !
       ! GP Shape functions (N)
       gpshap(1:DVS,1,igaus) = 0.25_rp*(1.0_rp + s*t - t - s)
       gpshap(1:DVS,2,igaus) = 0.25_rp*(1.0_rp - s*t + t - s)
       gpshap(1:DVS,3,igaus) = 0.25_rp*(1.0_rp + s*t + t + s)
       gpshap(1:DVS,4,igaus) = 0.25_rp*(1.0_rp - s*t - t + s)
       ! GP Derivectativecte Shape function w.r.t. natural coordinates (dN/dn)
       gpderi(1:DVS,1,1,igaus) = 0.25_rp*(-1.0_rp + s)
       gpderi(1:DVS,1,2,igaus) = 0.25_rp*( 1.0_rp - s)
       gpderi(1:DVS,1,3,igaus) = 0.25_rp*( 1.0_rp + s)
       gpderi(1:DVS,1,4,igaus) = 0.25_rp*(-1.0_rp - s)
       !
       gpderi(1:DVS,2,1,igaus) = 0.25_rp*(-1.0_rp + t)
       gpderi(1:DVS,2,2,igaus) = 0.25_rp*(-1.0_rp - t)
       gpderi(1:DVS,2,3,igaus) = 0.25_rp*( 1.0_rp + t)
       gpderi(1:DVS,2,4,igaus) = 0.25_rp*( 1.0_rp - t)
    end do
    !
    ! Rotation matrix
    !
    do igaus = 1,sgaus
       r1(1:DVS,:) = 0.0_rp
       r2(1:DVS,:) = 0.0_rp
       r3(1:DVS,:) = 0.0_rp
       do inode = 1,snode
          do idime = 1,3
             r1(1:DVS,idime) = r1(1:DVS,idime) + gpderi(1:DVS,1,inode,igaus)*sfcoo(1:DVS,idime,inode)
             r2(1:DVS,idime) = r2(1:DVS,idime) + gpderi(1:DVS,2,inode,igaus)*sfcoo(1:DVS,idime,inode)
          end do
       end do
       ! Compute normal direction vector
       r3(1:DVS,1) =  r1(1:DVS,2)*r2(1:DVS,3) - r1(1:DVS,3)*r2(1:DVS,2)
       r3(1:DVS,2) =  r1(1:DVS,3)*r2(1:DVS,1) - r1(1:DVS,1)*r2(1:DVS,3)
       r3(1:DVS,3) =  r1(1:DVS,1)*r2(1:DVS,2) - r1(1:DVS,2)*r2(1:DVS,1)
       ! Rotation matrix: 1s tangential and normal direction
       rM(1:DVS) = sqrt(r1(1:DVS,1)**2 + r1(1:DVS,2)**2 + r1(1:DVS,3)**2)
       !
       do ivect = 1,DVS
          gprota(ivect,:,:,igaus) = 0.0_rp
          gparea(ivect,igaus)     = 0.0_rp
          if( lmask_loc(ivect) ) then
             gprota(ivect,1,1,igaus) = r1(ivect,1)/rM(ivect)
             gprota(ivect,2,1,igaus) = r1(ivect,2)/rM(ivect)
             gprota(ivect,3,1,igaus) = r1(ivect,3)/rM(ivect)
             !
             gparea(ivect,igaus) = sqrt(r3(ivect,1)**2 + r3(ivect,2)**2 + r3(ivect,3)**2)
             gprota(ivect,1,3,igaus) = r3(ivect,1)/gparea(ivect,igaus)
             gprota(ivect,2,3,igaus) = r3(ivect,2)/gparea(ivect,igaus)
             gprota(ivect,3,3,igaus) = r3(ivect,3)/gparea(ivect,igaus)
             ! Rotation matrix: 2nd tangential
             gprota(ivect,1,2,igaus) = &
                  gprota(ivect,2,3,igaus)*gprota(ivect,3,1,igaus) - gprota(ivect,3,3,igaus)*gprota(ivect,2,1,igaus)
             gprota(ivect,2,2,igaus) = &
                  gprota(ivect,3,3,igaus)*gprota(ivect,1,1,igaus) - gprota(ivect,1,3,igaus)*gprota(ivect,3,1,igaus)
             gprota(ivect,3,2,igaus) = &
                  gprota(ivect,1,3,igaus)*gprota(ivect,2,1,igaus) - gprota(ivect,2,3,igaus)*gprota(ivect,1,1,igaus)
          end if
       end do
    end do
    !
    ! B matrix (B)
    !
    gpbmat(1:DVS,:,:,:) = 0.0_rp
    do igaus = 1,sgaus
       do inode = 1,snode
          do jdime = 1,3
             idofn = 3*(inode - 1) + jdime
             do idime = 1,3
                gpbmat(1:DVS,idime,idofn,igaus) = gprota(1:DVS,jdime,idime,igaus)*gpshap(1:DVS,inode,igaus)
             end do
          end do
       end do
    end do
    !
    ! Relative displacements / displacement jump (delta)
    !
    gpjump(1:DVS,:,:) = 0.0_rp
    do igaus = 1,sgaus
       do idime = 1,3
          do idofn = 1,pevah
             gpjump(1:DVS,idime,igaus) = gpjump(1:DVS,idime,igaus) + gpbmat(1:DVS,idime,idofn,igaus)*sfjum(1:DVS,idofn)
          end do
       end do
    end do
    !
    ! Constitutivecte relationship tensors (Traction and dTraction/dJump)
    !
    do ivect = 1,DVS
       if( lmask_loc(ivect) ) then
          do igaus = 1,sgaus
             call sld_vect_coh_current(&
                  dtime,properties(:),gpjump(ivect,:,igaus),gpsdvo(ivect,:,igaus),gpsdvn(ivect,:,igaus),gptrac(ivect,:,igaus))
          end do
       end if
    end do
    !
    ! Residual vector (RHS)
    !
    do igaus = 1,sgaus
       do idofn = 1,pevah
          do idime = 1,3
             elfin(1:DVS,idofn) = elfin(1:DVS,idofn) + &
                  gpbmat(1:DVS,idime,idofn,igaus)*gptrac(1:DVS,idime,igaus)*gpweig(1:DVS,igaus)*gparea(1:DVS,igaus)
          end do
       end do
    end do
    !
    ! Set GP state new variables
    !
    do concurrent ( ivect = 1:DVS )
       !if( lmask_loc(ivect) ) then
          svegm_sld(lelem_loc(ivect))%a(1:2,1:sgaus,ITER_K_STATE) = gpsdvn(ivect,1:2,1:sgaus)
       !end if
    end do
    deallocate(gpsdvo)
    deallocate(gpsdvn)

    !--------------------------------------------------------------------
    !
    ! Assembly RHS
    !
    !--------------------------------------------------------------------
    !
    ! Full RHS vector
    !
    peva1 = pevah + 1
    elrhs(1:DVS,1:pevah)               =  elfin(1:DVS,1:pevah)
    elrhs(1:DVS,peva1:pnode*ndofn_sld) = -elfin(1:DVS,1:pevah)
    !
    ! Save reaction forces at Dirichlet nodes
    !
    do ivect = 1,DVS
       if( lmask_loc(ivect) ) then
          do inode = 1,pnode
             ipoin = lnods_loc(ivect,inode)
             do idime = 1,3
                if( kfl_fixno_sld(idime,ipoin) > 0_ip ) then
                   idofn = ( inode - 1_ip) * ndofn_sld + idime
                   elfrx(ivect,idofn) = elrhs(ivect,idofn) 
                end if
             end do
          end do
       endif
    end do
    !
    ! Assemble to global RHS 
    !
    if( itask == 1_ip ) then
       do ivect = 1,DVS
          if( lmask_loc(ivect) ) then
             call matrix_assemble_element_RHS(ndofn_sld,ndofn_sld,pnode,lnods_loc(ivect,:),elrhs(ivect,:),rhsid )
             call matrix_assemble_element_RHS(ndofn_sld,ndofn_sld,pnode,lnods_loc(ivect,:),elfin(ivect,:),finte_sld)
             call matrix_assemble_element_RHS(ndofn_sld,ndofn_sld,pnode,lnods_loc(ivect,:),elfrx(ivect,:),frxid_sld)
          end if
       end do
    end if

    deallocate(lnods_loc)

  end subroutine sld_vect_element_operations_ELINT

  !-----------------------------------------------------------------------
  !> 
  !> @author  gguillam and aquintan
  !> @date    2022-10-04
  !> @brief   Vectorized version for cohesive law from 2018 (3d)
  !> @details Vectorized version for cohesive law from 2018 (3d)
  !> 
  !-----------------------------------------------------------------------
  
  subroutine sld_vect_coh_current(dt, properties, jumps, stateOld, stateNew, tu)

    implicit none
    external                     :: runend
    
    real(rp),      intent(in)    :: dt                   !< Time increment
    real(rp),      intent(in)    :: properties(:)        !< Cohesive law properties
    real(rp),      intent(in)    :: jumps(:)             !< Displacement jump (local coordinate system) (dime)

    real(rp),      intent(inout) :: stateOld(:)          !< State variables : old (for possible initialize purpouses)
    real(rp),      intent(inout) :: stateNew(:)          !< State variables : 
    real(rp),      intent(out)   :: tu(:)
         
    logical(lg)                  :: flagVisco            ! Viscosity flag
    integer(ip)                  :: kfLoad, kfDamage     ! Keyflags for damage and loading state
    real(rp)                     :: tauI, tauII          ! Properties: strength
    real(rp)                     :: GI, GII              ! Properties: toughness
    real(rp)                     :: eta                  ! Properties: mode-mix
    real(rp)                     :: rho                  ! Viscosity
    real(rp)                     :: dMax                 ! Max. damage variable
    real(rp)                     :: K11, K22, K33, KCONT ! Properties: penalty stiffness
    real(rp)                     :: JI0, JII0, JIC, JIIC ! Properties: pure mode softening (0) and propagation (C) jumps
    real(rp)                     :: J1, J2, J3, J3p, J3n ! Current jumps + effectivecte opening jump
    real(rp)                     :: Ep                   ! Positivecte elastic energy
    real(rp)                     :: B, KB, KS            ! Mixed mode ratio & Mode depenedent penalty stiffness
    real(rp)                     :: h, h0, hc            ! Equivectalent law jumps
    real(rp)                     :: r, d                 ! Damage variables
    real(rp)                     :: denom, numer         ! Auxiliary varaible
    !
    ! Initializations
    !
    ! Init flags
    kfLoad   = opening
    kfDamage = elastic
    !
    ! Material properties
    GI    = properties(1)
    GII   = properties(2)
    tauI  = properties(3)
    tauII = properties(4)
    eta   = properties(5)
    K33   = properties(6)
    KCONT = properties(6) ! Assumed equal to the user-defined
    rho   = properties(7)
    if( rho > 0.0_rp ) then
       flagVisco = .true.
    else
       flagVisco = .false.
    end if
    if( properties(8) > 0.0_rp ) then !(DVS)
       stateOld(1) = 1.0_rp
       stateOld(2) = 1.0_rp
    end if
    dMax  = properties(9)
    ! Defaults
    if( abs(dMax) < epsilon(1.0_rp) ) then
       dMax = 1.0_rp
    end if
    !
    ! Effectivecte opening jumps (DVS)
    J1 = jumps(1) 
    J2 = jumps(2)
    J3 = jumps(3)
    ! Numerical correction (values close to 0)
    ! Positive (J3p) and negativecte (J3n) normal jumps (DVS)
    if( J3 > tolPe ) then
       J3p = J3
       J3n = 0.0_rp
       kfLoad = opening
    else
       J3p = 0.0_rp
       J3n =  min(J3, 0.0_rp)
       kfLoad = closing
    endif
    !
    ! Penalty stiffness (mode-II / shear)
    K11 = K33*(gI/gII)*(tauII/tauI)**2
    K22 = K11         ! 2-d (K22 not exists); 3-d (K22 = K11 assumption)
    !
    ! Pure mode openings
    JI0  = tauI/K33
    JII0 = tauII/K11
    JIC  = 2.0_rp*GI/tauI
    JIIC = 2.0_rp*GII/tauII
    !
    ! Mix-mode ratio (DVS)
    denom = K11*J1**2 + K22*J2**2 + K33*J3p**2
    numer = K11*J1**2 + K22*J2**2
    if( denom /= 0.0_rp ) then
       B = numer/denom
    else
       B = 0.0_rp
    endif
    !
    ! Mode dependent penalty stiffness
    KS = K11
    KB = (1.0_rp - B)*K33 + B*KS
    !
    ! Equivectalent jumps
    ! - softening
    h0 = 0.0_rp
    if( KB > 0.0_rp ) then
       h0 = sqrt((K33*JI0**2 + (KS*JII0**2 - K33*JI0**2)*(B**eta))/KB)
    else
       call runend('MOD_SLD_INTERFACE_ELEMENT: COHESIVECTE TURON HYP A: KB < 0.0')
    end if
    ! - propagation
    hc = 0.0_rp
    if( KB > 0.0_rp .and. h0 > 0.0_rp ) then
       hc = (K33*(JI0*JIC) + (KS*(JII0*JIIC) - K33*(JI0*JIC))*(B**eta))/(KB*h0)
    else
       call runend('MOD_SLD_INTERFACE_ELEMENT: COHESIVECTE TURON HYP A: h0 < 0.0')
    endif
    ! - checking
    if( hc - h0 < 0.0_rp )then
       call runend('MOD_SLD_INTERFACE_ELEMENT: COHESIVECTE TURON HYP A: hc < h0')
    endif
    ! - current jump
    denom = sqrt((K11**2)*(J1**2) + (K22**2)*(J2**2) + (K33**2)*(J3p**2))
    numer = (K11*(J1**2) + K22*(J2**2) + K33*(J3p**2))
    if( denom /= 0.0_rp ) then
       h = numer/denom
    else
       h = 0.0_rp
    endif
    !
    ! Damage state          (DVS)
    ! - current threshold
    if( h < h0 )then
       r = 0.0_rp
       kfDamage = elastic
    else if( h0 <= h .and. h < hc  )then
       ! - damage threshold
       r = (h - h0)/(hc - h0)
       ! - Viscous regularization (Duvaut and Lions)
       if( flagVisco ) then
          r = r*dt/(rho + dt) + stateOld(2)*rho/(rho + dt)
       end if
       kfDamage = damage
    else
       r = 1.0_rp
       kfDamage = broken
    endif
    ! - historic threshold    (DVS)
    if( r < stateOld(2) )then
       r = stateOld(2)
       kfDamage = elastic
    end if
    ! - state                  (DVS)
    denom = r*hc + (1.0_rp - r)*h0
    if( denom > 0.0_rp ) then
       d = min((r*hc)/denom,dMax  )
    else
       d = 0.0_rp
    end if
    !
    ! Tractions (DVS)
    tu(:) = 0.0_rp
    tu(1) = (1.0_rp - d)*K11*jumps(1)  ! K11 = K22 assumption
    tu(2) = (1.0_rp - d)*K11*jumps(2)  ! K11 = K22 assumption
    if( kfLoad == closing )then
       tu(3) = KCONT*J3n
    else
       tu(3) = (1.0_rp - d)*K33*J3p
    end if
    !
    ! Elastic Energy
    Ep = 0.5_rp*( K11*J1**2 + K22*J2**2 + K33*J3p**2 )
    !
    ! Update state variables (DVS)
    !
    stateNew(1) = d
    stateNew(2) = r
 
  end subroutine sld_vect_coh_current
  
end module mod_sld_vect_assembly_elint
