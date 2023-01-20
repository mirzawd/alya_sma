!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Gusano
!> @{
!> @file    mod_gus_element_operations.f90
!> @author  houzeaux
!> @date    2020-10-20
!> @brief   Gusano assembly
!> @details Assembly
!>
!>          dA   dQ
!>          -- + -- = 0
!>          dt   dz 
!>
!>          dQ    d  a*Q^2    A  dp       Q  
!>          -- +  -- ----  + --- -- +  KR - = 0
!>          dt    dz  A      rho dz       A
!>
!>          Or:
!>
!>          dA   dQ
!>          -- + -- = 0
!>          dt   dz 
!>
!>          dQ   2*a   dQ   a*Q^2 dA    A  dp      Q  
!>          -- + --- Q -- - ----  -- + --- -- + KR - = 0
!>          dt    A    dz    A^2  dz   rho dz      A
!>
!>          so that
!>
!>          dQ
!>          -- = 0
!>          dz 
!>              +-                        -+
!>          rho | dQ   2*a   dQ   a*Q^2 dA |   dp         KR     
!>          --- | -- + --- Q -- - ----  -- | + -- + rho * -- Q = 0
!>           A  | dt    A    dz    A^2  dz |   dz         A^2
!>              +-                        -+
!>
!>          Some redifinitions:
!>
!>                       KR 
!>          KQ   = rho * -- = 128 mu / (pi * D^4) for laminar flow
!>                       A^2
!>
!>          rho' = rho / A
!>
!>          c    = 2 * a * u
!>                              a   dA
!>          r    = KQ - rho' * --   -- * Q
!>                             A^2  dz
!>
!>          So the momentum equation is:
!>
!>                 dQ          dQ   dp   
!>          rho' * -- + rho'*c*-- + -- + r*Q = 0, with
!>                 dt          dz   dz
!>
!>          For constant flow rate and areas, the momentum eq. reduces to:
!>
!>          dp
!>          -- + KQ Q = 0
!>          dz
!>
!-----------------------------------------------------------------------

module mod_gus_element_operations

#include "def_vector_size.inc"
  use def_kintyp_basic
  use def_parame
  use def_elmtyp
  use def_master
  use def_domain
  use def_kermod
  use mod_output
  use def_gusano
  use mod_elmgeo,       only : elmgeo_cartesian_derivatives_jacobian
  use mod_matrix,       only : matrix_assemble_element_matrix_to_CSR
  use mod_matrix,       only : matrix_assemble_element_RHS
  use mod_solver,       only : solver_assemble_element_matrix
  use mod_solver,       only : solver_assemble_element_RHS
  use mod_ker_proper,   only : ker_proper
  use mod_physics,      only : physics_resistance_parameter_Q
  use mod_parall,       only : num_subd_par
  use mod_parall,       only : num_pack_par
  use mod_parall,       only : list_elements_par
  use mod_parall,       only : typ_list_elements_par
  use mod_parall,       only : par_omp_nelem_chunk
  use mod_maths_arrays, only : maths_maxloc_nonzero
  use mod_communications_global
  implicit none

  integer(ip) :: mnode_loc
  
contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-11-02
  !> @brief   Element assembly
  !> @details Element assembly
  !> 
  !-----------------------------------------------------------------------

  subroutine gus_element_operations()

    use def_domain
    use def_kermod

#ifdef ALYA_OMPSS
    integer(ip)                          :: num_neigh
#endif
    integer(ip)                          :: isubd
    integer(ip)                          :: ipack,ielem,pnode,pgaus
    integer(ip)                          :: num_subd
    integer(ip),                 pointer :: num_pack(:)
    type(typ_list_elements_par), pointer :: list_elements(:)
    real(rp)                             :: taumin

    num_subd      =  num_subd_par
    num_pack      => num_pack_par
    list_elements => list_elements_par
    taumin = huge(1.0_rp)


#ifndef OPENACCHHH
#endif


    do isubd = 1,num_subd

#ifndef OPENACCHHH
#ifdef ALYA_OMPSS
       num_neigh = size(ompss_domains(isubd) % neighbours,KIND=ip)

       !-----------------------------------------------------------------------------
       !$OMP TASK         COMMUTATIVE(                                              &
       !$OMP              [ompss_domains(ompss_domains(isubd) % neighbours(jsubd)), &
       !$OMP              jsubd = 1,num_neigh] ) PRIORITY(num_neigh)                &
       !$OMP FIRSTPRIVATE ( num_neigh,jsubd,isubd )                                 &
       !$OMP SHARED       ( ompss_domains )                                         &
       !-----------------------------------------------------------------------------
#else
       !-----------------------------------------------------------------------------
       !$OMP PARALLEL DO                                                            &
       !$OMP SCHEDULE     ( DYNAMIC , par_omp_nelem_chunk )                         &
       !$OMP SHARED       ( isubd,par_omp_nelem_chunk )                             &
       !-----------------------------------------------------------------------------
#endif
       !-----------------------------------------------------------------------------
       !$OMP DEFAULT      ( SHARED )                                                &
       !$OMP PRIVATE      ( ipack,pnode,pgaus,ielem )                               &
       !$OMP SHARED       ( list_elements,num_pack,lnnod,lgaus,ltype                )              
       !-----------------------------------------------------------------------------
#endif

       do ipack = 1,num_pack(isubd)

          ielem = list_elements(isubd) % packs(ipack) % l(1)  ! Select first element
          pnode = lnnod(ielem)                                ! Number of nodes
          pgaus = lgaus(ielem)                                ! Number of Gauss points

          if( ltype(ielem) > 0 ) & 
               call gus_element_operations_assembly(&
               size(list_elements(isubd) % packs(ipack) % l,kind=ip),&
               pnode,pgaus,list_elements(isubd) % packs(ipack) % l,taumin)

       end do

#ifndef OPENACCHHH
#ifdef ALYA_OMPSS
       !$OMP END TASK
#else
       !$OMP END PARALLEL DO
#endif
#endif

    end do

#ifndef OPENACCHHH


#ifdef ALYA_OMPSS
    !$OMP  TASKWAIT
#endif
#else
    !$acc wait 
#endif

    call PAR_MIN(taumin)
    !print*,taumin
    
  end subroutine gus_element_operations

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-11-02
  !> @brief   Assembly
  !> @details Assembly
  !> 
  !-----------------------------------------------------------------------
  
  subroutine gus_element_operations_assembly(VECTOR_DIM,pnode,pgaus,list_elements,taumin)

    integer(ip), intent(in) :: VECTOR_DIM                        !< Vector size 
    integer(ip), intent(in) :: pnode                             !< Number of nodes
    integer(ip), intent(in) :: pgaus                             !< Number of Gauss points
    integer(ip), intent(in) :: list_elements(VECTOR_DIM)         !< List of elements
    real(rp),    intent(inout) :: taumin
    integer(ip)             :: ielem,inode,jnode
    integer(ip)             :: pdime,pelty,plapl
    integer(ip)             :: igaus,ipoin1,ipoin2
    integer(ip)             :: ipoin,dummi
    integer(ip)             :: idofv,idofp
    integer(ip)             :: jdofv,jdofp
    
    real(rp)                :: elmat(VECTOR_DIM,2*pnode,2*pnode) ! Matrix
    real(rp)                :: elrhs(VECTOR_DIM,2*pnode)         ! RHS
    real(rp)                :: elrhm(VECTOR_DIM,pnode)           ! Momentum or pressure gradient projection RHS
    real(rp)                :: elrhc(VECTOR_DIM,pnode)           ! Continuity projection RHS
    real(rp)                :: elsch(VECTOR_DIM,pnode,pnode)     ! Schur complement preconditioner
    real(rp)                :: elmat_1(2*pnode,2*pnode)          ! Matrix
    real(rp)                :: elrhs_1(2*pnode)                  ! RHS
    real(rp)                :: elrhm_1(pnode)                    ! Momentum or pressure gradient projection RHS
    real(rp)                :: elrhc_1(pnode)                    ! Continuity projection RHS
    
    real(rp)                :: eltes(VECTOR_DIM,pnode)           ! Momentum test function
    real(rp)                :: elmom(VECTOR_DIM,pnode)
  
    real(rp)                :: elcod(VECTOR_DIM,1,pnode)         ! Element gather coord
    real(rp)                :: elflo(VECTOR_DIM,pnode,2)         ! Element gather vel1d
    real(rp)                :: elpre(VECTOR_DIM,pnode,2)         ! Element gather press
    real(rp)                :: elprm(VECTOR_DIM,pnode)           ! Element gather projm_gus
    real(rp)                :: elprc(VECTOR_DIM,pnode)           ! Element gather projc_gus
    real(rp)                :: elare(VECTOR_DIM,pnode)           ! Element gather are

    real(rp)                :: gpcar(VECTOR_DIM,1,pnode,pgaus)   ! dNi/dz
    real(rp)                :: gpsha(VECTOR_DIM,pnode,pgaus)     ! Ni
    real(rp)                :: gpvol(VECTOR_DIM,pgaus)           ! |J|w
    real(rp)                :: gpdet(VECTOR_DIM,pgaus)           ! |J|
    real(rp)                :: xjaci(VECTOR_DIM,1,1,pgaus)       ! J^-1
    
    real(rp)                :: gpden(VECTOR_DIM,pgaus)           ! rho
    real(rp)                :: gpvis(VECTOR_DIM,pgaus)           ! mu
    real(rp)                :: alpha                             ! alpha
    real(rp)                :: chale(VECTOR_DIM)                 ! h
    real(rp)                :: dAdz(VECTOR_DIM)                  ! dA/dz
    real(rp)                :: A(VECTOR_DIM)                     ! A
    real(rp)                :: D(VECTOR_DIM)                     ! D
    real(rp)                :: gpflo(VECTOR_DIM,2)               ! Q
    real(rp)                :: gpdiv(VECTOR_DIM)                 ! div(Q)
    real(rp)                :: gpadv(VECTOR_DIM)                 ! c
    real(rp)                :: gptau(VECTOR_DIM)                 ! tau
    real(rp)                :: gpta2(VECTOR_DIM)                 ! tau2
    real(rp)                :: gprea(VECTOR_DIM)                 ! r
    real(rp)                :: grpre(VECTOR_DIM)                 ! dp/dz
    real(rp)                :: gpprm(VECTOR_DIM)                 ! Projection momentum
    real(rp)                :: gpprc(VECTOR_DIM)                 ! Projection continuity
    real(rp)                :: gpmom(VECTOR_DIM)                 ! Momentum operator
    real(rp)                :: gprhm(VECTOR_DIM)                 ! RHS of momentum 
    real(rp)                :: gpvel(VECTOR_DIM)                 ! u
    real(rp)                :: fconv
    integer(ip)             :: list_elements_p(VECTOR_DIM)       ! List of elements (always positive)
    integer(ip)             :: pelem,ivect
    
#ifdef OPENACCHHH
#define DEF_VECT ivect
#else
#define DEF_VECT 1:VECTOR_SIZE
#endif

    fconv = real(kfl_conve_gus,rp)
    pelem = maths_maxloc_nonzero(list_elements)
    pelty = ltype(list_elements(1))
    do ivect = 1,pelem 
       list_elements_p(ivect) = list_elements(ivect)
    end do
    do ivect = pelem+1,VECTOR_DIM
       list_elements_p(ivect) = list_elements(pelem)
    end do
    do ivect = 1,VECTOR_DIM      
       ielem = list_elements_p(ivect)
       gpsha(ivect,1:pnode,1:pgaus) = elmar(pelty) % shape(1:pnode,1:pgaus)
    end do
    
    !--------------------------------------------------------------
    !
    ! Gather
    !
    !--------------------------------------------------------------
    
    do ivect = 1,VECTOR_DIM
       ielem             = list_elements_p(ivect)
       ipoin1            = lnods(1,ielem)
       ipoin2            = lnods(2,ielem)
       elcod(ivect,1,1)  = 0.0_rp
       chale(ivect)      = sqrt(dot_product(coord(:,ipoin1)-coord(:,ipoin2),coord(:,ipoin1)-coord(:,ipoin2)))
       elcod(ivect,1,2)  = chale(ivect)
       do inode = 1,pnode
          ipoin                = lnods(inode,ielem)
          elflo(ivect,inode,1) = flowr(ipoin,1)
          elflo(ivect,inode,2) = flowr(ipoin,3)
          elpre(ivect,inode,1) = press(ipoin,1)
          elpre(ivect,inode,2) = press(ipoin,3)
          elprm(ivect,inode)   = projm_gus(ipoin,2)
          elprc(ivect,inode)   = projc_gus(ipoin,2)
          elare(ivect,inode)   = areas(ipoin,1)
       end do       
    end do
       
    !--------------------------------------------------------------
    !
    ! Cartesian derivatives and Jacobian
    !
    !--------------------------------------------------------------
    
    plapl = 0
    pdime = 1
    call elmgeo_cartesian_derivatives_jacobian(pdime,pnode,pnode,pgaus,elcod,elmar(pelty) % deriv,xjaci,gpcar,gpdet)    
    do igaus = 1,pgaus
       gpvol(DEF_VECT,igaus) = elmar(pelty) % weigp(igaus) * gpdet(DEF_VECT,igaus)
    end do

    !--------------------------------------------------------------
    !
    ! Properties
    !
    !--------------------------------------------------------------
    
    call ker_proper('DENSI','PGAUS',dummi,list_elements_p,gpden,pnode,pgaus,gpsha,VECTOR_DIM=VECTOR_DIM) ! rho
    call ker_proper('VISCO','PGAUS',dummi,list_elements_p,gpvis,pnode,pgaus,gpsha,VECTOR_DIM=VECTOR_DIM) ! mu
    !
    ! Momentum correction coefficient
    !
    alpha = 4.0_rp / 3.0_rp
    
    !--------------------------------------------------------------
    !
    ! Gauss point values from element values
    !
    !--------------------------------------------------------------
    
    elmat = 0.0_rp
    elrhs = 0.0_rp
    elrhm = 0.0_rp
    elrhc = 0.0_rp
    elsch = 0.0_rp
    
    do igaus = 1,pgaus
       !
       ! Initialization
       !
       A     = 0.0_rp
       dAdz  = 0.0_rp
       gpflo = 0.0_rp
       gpdiv = 0.0_rp
       grpre = 0.0_rp
       gpprm = 0.0_rp
       gpprc = 0.0_rp
       gpmom = 0.0_rp
       !
       ! Element to gauss point
       !
       do inode = 1,pnode
          A    (DEF_VECT)   = A    (DEF_VECT)   + elmar(pelty) % shape(inode,igaus) * elare(DEF_VECT,inode)            ! A
          dAdz (DEF_VECT)   = dAdz (DEF_VECT)   + gpcar(DEF_VECT,1,inode,igaus)     * elare(DEF_VECT,inode)            ! dA/dz
          gpflo(DEF_VECT,1) = gpflo(DEF_VECT,1) + elmar(pelty) % shape(inode,igaus) * elflo(DEF_VECT,inode,1)          ! Q
          gpflo(DEF_VECT,2) = gpflo(DEF_VECT,2) + elmar(pelty) % shape(inode,igaus) * elflo(DEF_VECT,inode,2)          ! Q^n
          gpdiv(DEF_VECT)   = gpdiv(DEF_VECT)   + gpcar(DEF_VECT,1,inode,igaus)     * elflo(DEF_VECT,inode,1)          ! grad(Q)
          grpre(DEF_VECT)   = grpre(DEF_VECT)   + gpcar(DEF_VECT,1,inode,igaus)     * elpre(DEF_VECT,inode,1)          ! grad(p)
          gpprm(DEF_VECT)   = gpprm(DEF_VECT)   + elmar(pelty) % shape(inode,igaus) * elprm(DEF_VECT,inode)            ! P(tau*momentum) or P(tau*p)
          gpprc(DEF_VECT)   = gpprc(DEF_VECT)   + elmar(pelty) % shape(inode,igaus) * elprc(DEF_VECT,inode)            ! P(tau2*div(Q))
       end do
       gpvel(DEF_VECT) = gpflo(DEF_VECT,1) / A(DEF_VECT)                                                               ! u
       D(DEF_VECT)     = 2.0_rp * sqrt(A(DEF_VECT)/pi)                                                                 ! D
       gpadv(DEF_VECT) = 2.0_rp * fconv * alpha * gpflo(DEF_VECT,1) / A(DEF_VECT)                                      ! c
       gprea(DEF_VECT) = physics_resistance_parameter_Q(&
            &            gpden(DEF_VECT,igaus),gpvel(DEF_VECT),D(DEF_VECT),gpvis(DEF_VECT,igaus),kfl_regim_gus)        ! KR
       gpden(DEF_VECT,igaus) = gpden(DEF_VECT,igaus) / A(DEF_VECT)                                                     ! rho' = rho / A
       gprea(DEF_VECT) = gprea(DEF_VECT) &
            &            - gpden(DEF_VECT,igaus) * alpha / A(DEF_VECT)**2 * gpflo(DEF_VECT,1) * fconv * dAdz(DEF_VECT) ! r
       gptau(DEF_VECT) = 1.0_rp / ( 2.0_rp * gpden(DEF_VECT,igaus) * abs(gpadv(DEF_VECT)) &
            &            / chale(DEF_VECT)  + abs(gprea(DEF_VECT)) )                                                   ! tau = 1/(2*rho'*|c|/h +|r|)
       gpta2(DEF_VECT) = 4.0_rp * chale(DEF_VECT) * chale(DEF_VECT) / gptau(DEF_VECT)                                  ! tau2 = 4 * h^2 / tau
       gprhm(DEF_VECT) = dtinv_gus * gpden(DEF_VECT,igaus) * gpflo(DEF_VECT,2)

       do ivect = 1,VECTOR_DIM
          taumin = min(taumin,gpden(ivect,igaus)*gptau(ivect)/A(ivect))
       end do
       !
       ! Momentum test function and momentum operator
       !
       do inode = 1,pnode
          eltes(DEF_VECT,inode) =  elmar(pelty) % shape(inode,igaus) + gpden(DEF_VECT,igaus) * &        
               &                   gptau(DEF_VECT) * gpadv(DEF_VECT) * gpcar(DEF_VECT,1,inode,igaus) 
          elmom(DEF_VECT,inode) =  ( gprea(DEF_VECT) + gpden(DEF_VECT,igaus) * dtinv_gus ) &
               &                   * elmar(pelty) % shape(inode,igaus) &                                              ! Mi(Ni)
               &                   + gpden(DEF_VECT,igaus) * gpadv(DEF_VECT) * gpcar(DEF_VECT,1,inode,igaus)        
          gpmom(DEF_VECT)       =  gpmom(DEF_VECT) + elmom(DEF_VECT,inode) * elflo(DEF_VECT,inode,1)                  ! M(Q)
       end do
       !
       ! Elemental matrices and RHS
       !
       do jnode = 1,pnode
          jdofv = (jnode-1)*2+1
          jdofp = (jnode-1)*2+2
          do inode = 1,pnode             
             idofv = (inode-1)*2+1
             idofp = (inode-1)*2+2
             !
             ! Auu: (rho * Q/dt + rho * c * dQ/dz + r * Q,v-tau*L^*(v)) + (tau2 * dQdz,dv/dz)
             !
             elmat(DEF_VECT,idofv,jdofv) = elmat(DEF_VECT,idofv,jdofv) &
                  &                        + gpvol(DEF_VECT,igaus) * eltes(DEF_VECT,inode) * elmom(DEF_VECT,jnode) &
                  &                        + gpvol(DEF_VECT,igaus) * gpta2(DEF_VECT)       * gpcar(DEF_VECT,1,inode,igaus) * gpcar(DEF_VECT,1,jnode,igaus)
             !
             ! Aup: - (dv/dz,p) + (tau * c * dQ/dz, dq/dz)
             !
             elmat(DEF_VECT,idofv,jdofp) = elmat(DEF_VECT,idofv,jdofp) &
                  &                        - gpvol(DEF_VECT,igaus) * gpcar(DEF_VECT,1,inode,igaus) * elmar(pelty) % shape(jnode,igaus) &
                  &                        + gpvol(DEF_VECT,igaus) * gptau(DEF_VECT) * gpadv(DEF_VECT) * gpcar(DEF_VECT,1,inode,igaus) * gpcar(DEF_VECT,1,jnode,igaus)
             !
             ! Apu: (dQ/dz,q) + (tau * M(Q),gradq)
             !
             elmat(DEF_VECT,idofp,jdofv) = elmat(DEF_VECT,idofp,jdofv) &
                  &                        + gpvol(DEF_VECT,igaus) * elmar(pelty) % shape(inode,igaus) * gpcar(DEF_VECT,1,jnode,igaus)
             elmat(DEF_VECT,idofp,jdofv) = elmat(DEF_VECT,idofp,jdofv) &
                  &                        + gpvol(DEF_VECT,igaus) * gptau(DEF_VECT) * gpcar(DEF_VECT,1,inode,igaus) * elmom(DEF_VECT,jnode)             
             !
             ! App: (tau * gradp,gradq)
             !
             elmat(DEF_VECT,idofp,jdofp) = elmat(DEF_VECT,idofp,jdofp) &
                  &                        + gpvol(DEF_VECT,igaus) * gptau(DEF_VECT) * gpcar(DEF_VECT,1,inode,igaus) * gpcar(DEF_VECT,1,jnode,igaus)                 
                       
          end do
          elrhs(DEF_VECT,jdofv)   = elrhs(DEF_VECT,jdofv) + gpvol(DEF_VECT,igaus) * eltes(DEF_VECT,jnode) * gprhm(DEF_VECT)                   ! (rho * Q^n / dt,v-tau*L^*(v))
          elrhs(DEF_VECT,jdofp)   = elrhs(DEF_VECT,jdofp) + gpvol(DEF_VECT,igaus) * gpcar(DEF_VECT,1,jnode,igaus) * &
               &                          ( xsoss_gus * gpprm(DEF_VECT) + gptau(DEF_VECT) * gprhm(DEF_VECT) )                                 ! (P + atu * ru,grad(q))
          !
          ! Projections
          !
          elrhm(DEF_VECT,jnode)   = elrhm(DEF_VECT,jnode) + gpvol(DEF_VECT,igaus) * elmar(pelty) % shape(jnode,igaus) * gptau(DEF_VECT) * &
               &                          ( xsoss_gus * grpre(DEF_VECT) + xoss_gus * (gpmom(DEF_VECT)-gprhm(DEF_VECT)))                       ! tau * q * (grad(p)-L(Q)) = Pm
          elrhc(DEF_VECT,jnode)   = elrhc(DEF_VECT,jnode) + gpvol(DEF_VECT,igaus) * elmar(pelty) % shape(jnode,igaus) * gpta2(DEF_VECT) * &
               &                          xoss_gus * gpdiv(DEF_VECT)                                                                          ! tau2 * q * grad(Q) = Pc
       end do
       
       if( kfl_algor_gus == GUS_SCHUR_COMPLEMENT ) then
          do jnode = 1,pnode
             jdofp = (jnode-1)*2+2
             do inode = 1,pnode
                idofp = (inode-1)*2+2
                elsch(DEF_VECT,inode,jnode) = elmat(DEF_VECT,idofp,jdofp) &
                     + gpvol(DEF_VECT,igaus) / ( dtinv_gus * gpden(DEF_VECT,igaus) + 1.0_rp / gptau(DEF_VECT) ) &
                     * gpcar(DEF_VECT,1,inode,igaus) * gpcar(DEF_VECT,1,jnode,igaus)     
             end do
          end do
       end if
       
    end do
    !
    ! Assembly in ALgebraic system
    !
    do ivect = 1,pelem
       ielem        = list_elements(ivect)
       elmat_1(:,:) = elmat(ivect,:,:)
       elrhs_1(:)   = elrhs(ivect,:)
       elrhm_1(:)   = elrhm(ivect,:)
       elrhc_1(:)   = elrhc(ivect,:)
       call solver_assemble_element_matrix(&
            solve(1),2_ip,pnode,2_ip*pnode,ielem,lnods(:,ielem),elmat(ivect,:,:),amatr)
       call solver_assemble_element_RHS(&
            solve(1),2_ip,pnode,2_ip*pnode,lnods(:,ielem),elrhs(ivect,:),rhsid)
       if( kfl_algor_gus == GUS_SCHUR_COMPLEMENT ) then
          call solver_assemble_element_matrix(&
               solve(3),1_ip,pnode,pnode,ielem,lnods(:,ielem),elsch(ivect,:,:),schur_gus)          
       end if
          !call matrix_assemble_element_matrix_to_CSR(&
       !     kfl_element_to_csr,2_ip,pnode,2*pnode,&
       !     ielem,lnods(:,ielem),elmat_1,r_dom,c_dom,amatr,lezdo)
       !call matrix_assemble_element_RHS(&
       !     2_ip,2_ip,pnode,lnods(:,ielem),elrhs_1,rhsid)
       call matrix_assemble_element_RHS(&
            1_ip,1_ip,pnode,lnods(:,ielem),elrhm_1,projm_gus(:,1))
       call matrix_assemble_element_RHS(&
            1_ip,1_ip,pnode,lnods(:,ielem),elrhc_1,projc_gus(:,1))           
    end do

  end subroutine gus_element_operations_assembly
 
end module mod_gus_element_operations
!> @}
